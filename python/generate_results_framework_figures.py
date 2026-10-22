"""Generate simple figures for Docs/Results_framework.md from the full model.

Run from the repository root:
    python3 python/generate_results_framework_figures.py \
        --draft-pdf ~/Downloads/Results_draft.pdf

Reuses the existing finite-rate Python equations with parameters read explicitly
from default_parameters.m (the Python dataclass defaults are historical).
The gene-length calculation reproduces the MATLAB 5-E-point / 10-length
occupancy interpolation and 200-point population quadrature used in the draft.
It is not an equilibrium approximation. Figure labels use the author's intended
assembly/commitment interpretation of the irreversible kc exit.

PNG/SVG figures, CSV data, and provenance go to Results/ResultsFramework/draft/.
The empirical abundance panel is a crop of the supplied PDF. Spacing counts
are approximately digitized from its histogram image, not new measurements. Rendering them requires Poppler's pdftoppm and Pillow.
Use --reuse-data to redraw saved calculations after changing plot styling.
"""
import argparse
import ast
import csv
import gzip
from dataclasses import asdict, fields, replace
import hashlib
import inspect
import json
import operator
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile

os.environ.setdefault('MPLCONFIGDIR', str(Path(tempfile.gettempdir()) / 'cpad-matplotlib'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, FancyArrowPatch
import numpy as np
from PIL import Image
import scipy
from scipy.interpolate import RegularGridInterpolator
from scipy.optimize import brentq
from scipy.sparse.linalg import spsolve

from cpadynamics.full_model import FullModel
from cpadynamics.gene_length import GeneLengthModel
from cpadynamics.observables import cleavage_profile
from cpadynamics.parameters import Parameters
from plot_ebinding_profile import binding_profiles

ROOT = Path(__file__).resolve().parents[1]
TEAL, ORANGE, PURPLE = '#156f8a', '#cf6e17', '#8f438b'
GREEN, INK, GREY = '#33985a', '#253248', '#687382'
COLORS = {1: ORANGE, 3: PURPLE, 5: TEAL}
STRONG_KHOFF = 0.01  # RunResultsAnalyses.m, strong_kHoff()


def read_parameters(path):
    """Read scalar MATLAB assignments using only literal arithmetic (no eval)."""
    values = {}
    allowed = {f.name for f in fields(Parameters)}
    operations = {ast.Add: operator.add, ast.Sub: operator.sub,
                  ast.Mult: operator.mul, ast.Div: operator.truediv}

    def number(node):
        if isinstance(node, ast.Constant) and type(node.value) in (int, float):
            return node.value
        if isinstance(node, ast.Attribute) and isinstance(node.value, ast.Name) and node.value.id == 'P':
            return values[node.attr]
        if isinstance(node, ast.BinOp) and type(node.op) in operations:
            return operations[type(node.op)](number(node.left), number(node.right))
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
            return -number(node.operand)
        raise ValueError(f'Unsupported MATLAB scalar expression: {ast.dump(node)}')

    source = re.sub(r'%[^\n]*', '', path.read_text())
    for name, expression in re.findall(r'\bP\.(\w+)\s*=\s*([^;]+);', source):
        if name in allowed:
            values[name] = float(number(ast.parse(expression.strip(), mode='eval').body))
    if allowed != values.keys():
        raise ValueError(f'Missing MATLAB parameters: {sorted(allowed-values.keys())}')
    ranges = {}
    for name in ('kHon', 'kHoff', 'kEon', 'kEoff', 'kc'):
        match = re.search(r'P\.ranges\.' + name + r'\s*=\s*\[([^]]+)\]', source)
        if not match:
            raise ValueError(f'Missing sensitivity range: {name}')
        ranges[name] = [float(x) for x in match[1].split(',')]
    ranges['kHd'] = [ranges['kHoff'][0]/ranges['kHon'][1], ranges['kHoff'][1]/ranges['kHon'][0]]
    ranges['kEd'] = [ranges['kEoff'][0]/ranges['kEon'][1], ranges['kEoff'][1]/ranges['kEon'][0]]
    return Parameters(**values), ranges


def solve_case(p, m, diagnostics, free_pools=None):
    model = FullModel(p, m)
    if free_pools is None:
        solution = model.solve_steady_state()
    else:
        rf, ef = free_pools
        occupation = spsolve(-(model.a0 + ef*model.ae), model.unit_source)
        solution = model.summarize(np.r_[occupation*p.k_in*rf, ef, rf])
    cdf = cleavage_profile(solution.R, solution.RHE, p)
    avg_e, avg_p = binding_profiles(model, solution)
    residual = abs(cdf.total_outflux-p.k_in*solution.Pol_free)
    if (solution.state.min() < -1e-9 or solution.rhs_max_abs > 1e-6 or residual > 1e-6
            or np.any(np.diff(cdf.cdf) < -1e-10) or cdf.max_exit_cdf > 1+1e-9):
        raise RuntimeError(f'Invalid finite-rate solution: M={m}')
    if free_pools is None and max(abs(solution.e_residual), abs(solution.pol_residual)) > 1e-6:
        raise RuntimeError('Closed-pool conservation failed')
    pre = model.pas-1
    loaded = solution.R_micro[pre, np.asarray(model.r_states)[:, 1] > 0].sum()/solution.R[pre]
    diagnostics.append(dict(M=m, pool_mode='closed' if free_pools is None else 'fixed',
                            rhs_max_abs=solution.rhs_max_abs, flux_residual=float(residual),
                            max_commitment_cdf=cdf.max_exit_cdf))
    profile = np.column_stack(((np.arange(model.pas+1)-model.pas)*p.L_a/1000,
                               avg_e[:model.pas+1], avg_p[:model.pas+1]))
    return cdf, solution, loaded, profile


def gene_length_data(p, diagnostics):
    """Reproduce MATLAB's draft grid rather than silently changing quadrature."""
    p = replace(p, kHoff=STRONG_KHOFF, k_in=p.k_in/10000)
    length_axis = np.geomspace(2500, 200000, 10)
    e_axis = np.linspace(0, p.E_total, 5)
    grid = GeneLengthModel(p, 5, np.floor(length_axis/p.L_a)*p.L_a)
    profiles = [grid.evaluate(ef) for ef in e_axis]
    # Occupancy scales exactly with free Pol II, so the five MATLAB R nodes
    # can be eliminated without changing their linear interpolation result.
    pol = RegularGridInterpolator((e_axis, length_axis), np.array([s.pol_per_flux for s in profiles]))
    e = RegularGridInterpolator((e_axis, length_axis), np.array([s.e_per_flux for s in profiles]))
    lengths = np.geomspace(2500, 200000, 200)
    dx = np.diff(lengths)
    weights = np.r_[dx[0]/2, (dx[:-1]+dx[1:])/2, dx[-1]/2]
    sigma = .68*np.log(10)
    weights *= np.exp(-.5*(np.log(lengths/22000)/sigma)**2)/(lengths*sigma*np.sqrt(2*np.pi))
    mass = weights.sum()
    weights /= mass

    def balance(ef):
        points = np.column_stack((np.full(len(lengths), ef), lengths))
        rf = p.Pol_total/(1+10000*p.k_in*(weights @ pol(points)))
        return ef+10000*p.k_in*rf*(weights @ e(points))-p.E_total, rf

    ef = brentq(lambda x: balance(x)[0], 0, p.E_total, xtol=1e-9)
    e_residual, rf = balance(ef)
    if abs(e_residual) > 1e-6:
        raise RuntimeError('Interpolated genome pool conservation failed')
    nominal = np.geomspace(2500, 200000, 20)
    lattice = np.floor(nominal/p.L_a)*p.L_a
    rows = [nominal, lattice]
    for m in (1, 5):
        conditional = GeneLengthModel(p, m, lattice).evaluate(ef)
        rows.extend([conditional.cad50_bp, conditional.cleavage_cdf[:, -1]])
    data = {'length_curves': np.column_stack(rows)}
    for length in (2500, 5000, 20000, 40000, 200000):
        q = replace(p, PASposition=length, geneLength_bp=length+5000)
        for m in (1, 5):
            cdf, _, _, profile = solve_case(q, m, diagnostics, (rf, ef))
            reference = GeneLengthModel(p, m, [length]).evaluate(ef)
            np.testing.assert_allclose(cdf.cdf, reference.cleavage_cdf[0], atol=1e-10)
            if m == 5 and length in (5000, 20000):
                data[f'length_profile_{length}'] = profile
            data[f'length_check_{m}_{length}'] = np.array([cdf.cad50_bp, cdf.max_exit_cdf])
    info = dict(parameters=asdict(p), E_free=float(ef), Pol_free=float(rf),
                E_conservation_residual=float(e_residual), num_genes=10000,
                population_median_bp=22000, population_log10_sigma=.68,
                population_mass_before_normalization=float(mass),
                grid=dict(E_points=5, L_points=10, L_interval_bp=[2500, 200000], integration_points=200),
                method='Full finite-rate occupancy grid; linear interpolation; truncated lognormal quadrature.',
                comparison='M=1 and M=5 at the same pools solved for M=5, matching the draft; not a genome-wide capacity perturbation.',
                sparse_reference_checks='CDF agreement at 2.5, 5, 20, 40, 200 kb for M=1 and M=5.')
    return data, info


def compute_data(p, ranges, points):
    data, diagnostics, summary = {}, [], []
    print('Computing capacity, preloading and assembly curves...', flush=True)
    for m in range(1, 11):
        cdf, solution, loaded, profile = solve_case(p, m, diagnostics)
        summary.append([m, cdf.cad50_bp, loaded, solution.E_free, solution.Pol_free])
        data[f'profile_M{m}'] = profile
        data[f'cdf_M{m}'] = np.column_stack((cdf.distances_bp, cdf.cdf))
    data['capacity'] = np.array(summary)
    for name in ('kHd', 'kEd', 'kc'):
        data[f'axis_{name}'] = np.geomspace(*ranges[name], points)
    pairs = [('kHd', 'kEd'), ('kc', 'kHd'), ('kc', 'kEd')]
    for m in (1, 5):
        for pair, (xn, yn) in enumerate(pairs):
            print(f'Computing M={m} {xn}/{yn} map...', flush=True)
            cad, reached = np.zeros((points, points)), np.zeros((points, points))
            for iy, y in enumerate(data[f'axis_{yn}']):
                for ix, x in enumerate(data[f'axis_{xn}']):
                    overrides = {}
                    for name, value in ((xn, x), (yn, y)):
                        field, multiplier = {'kHd': ('kHoff', p.kHon),
                                             'kEd': ('kEoff', p.kEon), 'kc': ('kc', 1)}[name]
                        overrides[field] = value*multiplier
                    cdf, _, _, _ = solve_case(replace(p, **overrides), m, diagnostics)
                    cad[iy, ix], reached[iy, ix] = cdf.cad50_bp, cdf.max_exit_cdf
            data[f'map_{m}_{pair}'], data[f'map_maxcdf_{m}_{pair}'] = cad, reached
    print('Computing gene-length profiles at the draft shared pools...', flush=True)
    length_data, length_info = gene_length_data(p, diagnostics)
    data.update(length_data)
    print('Computing APA commitment curves...', flush=True)
    distances = np.arange(0, 5001, 100)
    apa = [distances]
    for kh in (.001, .1, 1, 5, 10):
        cdf, _, _, _ = solve_case(replace(p, kHoff=kh, geneLength_bp=50000), 5, diagnostics)
        apa.append(np.interp(distances, cdf.distances_bp, cdf.cdf))
    data['apa_spacing'] = np.column_stack(apa)
    ratios = np.geomspace(.25, 4, 15)
    response = []
    for ratio in ratios:
        cdf, _, _, _ = solve_case(replace(p, kHoff=STRONG_KHOFF, E_total=p.E_total*ratio), 5, diagnostics)
        response.append([ratio, *np.interp([300, 1500], cdf.distances_bp, cdf.cdf)])
    data['apa_E_supply'] = np.array(response)
    # Compare to labels actually visible in the attached PDF, allowing only
    # rounding of its integer labels; deviations are reported, never fitted.
    comparisons = []
    for label, value, reference in (
            ('CAD50 M1', data['capacity'][0, 1], 783),
            ('CAD50 M5', data['capacity'][4, 1], 331),
            ('Gene length 2.5 kb M5', data['length_check_5_2500'][0], 3460),
            ('Gene length 200 kb M5', data['length_check_5_200000'][0], 288)):
        comparisons.append(dict(label=label, computed=float(value), draft_label=reference,
                                difference=float(value-reference), matches_rounding=bool(abs(value-reference) <= .5)))
    return data, dict(diagnostics=diagnostics, gene_length=length_info, draft_comparisons=comparisons)


def configure_style():
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
                         'axes.titlesize': 11, 'axes.labelsize': 10,
                         'axes.titleweight': 'bold', 'axes.edgecolor': '#555555',
                         'axes.linewidth': .8, 'axes.grid': True, 'grid.alpha': .22,
                         'grid.linewidth': .6, 'lines.linewidth': 2,
                         'legend.fontsize': 9, 'legend.frameon': False,
                         'figure.facecolor': 'white', 'savefig.facecolor': 'white',
                         'svg.fonttype': 'none'})


def save(fig, output, name):
    for suffix in ('png', 'svg'):
        fig.savefig(output/f'{name}.{suffix}', dpi=200, bbox_inches='tight')
    plt.close(fig)


def panel(ax, label, title):
    ax.set_title(f'{label}  {title}', loc='left', pad=12)


def cartoon_axes(ax):
    ax.set(xlim=(0, 10), ylim=(0, 4))
    ax.axis('off')


def arrow(ax, a, b, label=None, color=GREY, both=False):
    ax.add_patch(FancyArrowPatch(a, b, arrowstyle='<->' if both else '-|>',
                                mutation_scale=13, lw=1.2, color=color))
    if label:
        ax.text((a[0]+b[0])/2, (a[1]+b[1])/2+.16, label, ha='center', fontsize=9, color=color)


def polymerase(ax, x, y, occupied=1, capacity=5, available=5):
    ax.add_patch(Ellipse((x, y), .85, .55, color=INK, zorder=3))
    ax.text(x, y, 'Pol II', ha='center', va='center', color='white', fontsize=8, zorder=4)
    # CTD drawn as a bent tether; green dots are sites, large teal dots are E.
    sx = np.linspace(x-1.6, x-.48, capacity)
    sy = y+.48+.12*np.sin(np.linspace(0, np.pi, capacity))
    ax.plot(np.r_[sx, x-.2], np.r_[sy, y+.22], color=GREY, lw=1.5)
    ax.scatter(sx, sy, s=25, color='#cbd0d5', zorder=4)
    if available:
        ax.scatter(sx[-available:], sy[-available:], s=25, color=GREEN, zorder=4)
    if occupied:
        ax.scatter(sx[-occupied:], sy[-occupied:]+.12, s=80, color=TEAL,
                   edgecolors='white', linewidths=.7, zorder=5)


def gene(ax, y, tss=1, pas=7.5, end=9.3, pol=6.5, occupied=1, available=5):
    ax.plot([tss, end], [y, y], color=INK, lw=1.3)
    for x, label in ((tss, 'TSS'), (pas, 'PAS')):
        ax.plot([x, x], [y-.12, y+.13], color=INK)
        ax.text(x, y-.4, label, ha='center', fontsize=9)
    polymerase(ax, pol, y+.08, occupied, available=available)
    arrow(ax, (tss+.25, y+.86), (pol-.25, y+.86))


def plot_cartoons(output):
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.5), constrained_layout=True)
    for ax in axes:
        cartoon_axes(ax)
    panel(axes[0], 'A', 'E loading during elongation')
    gene(axes[0], 1.35, pol=5.6, occupied=2)
    axes[0].text(1, 3.3, 'Ser2P sites', color=GREEN, fontsize=10)
    axes[0].text(1, 2.95, 'E factors', color=TEAL, fontsize=10)
    panel(axes[1], 'B', 'Reversible PAS recognition')
    gene(axes[1], 1.35, pol=8.2, occupied=2)
    axes[1].plot([7.9, 7.2, 6.9, 6.3], [1.55, 2.25, 2.45, 2.35], color=ORANGE, lw=2)
    axes[1].text(4, 3.1, 'Nascent RNA / PAS', color=ORANGE, fontsize=10)
    arrow(axes[1], (6.3, 2.22), (6.8, 1.98), both=True, color=ORANGE)
    panel(axes[2], 'C', 'Irreversible commitment')
    polymerase(axes[2], 3.2, 1.9, occupied=2)
    arrow(axes[2], (4, 1.9), (6.1, 1.9), label=r'$k_c$')
    axes[2].add_patch(Ellipse((7.8, 1.9), 2.6, 1, facecolor='#e1f0ef', edgecolor=TEAL))
    axes[2].text(7.8, 1.9, 'Committed\nCPA complex', ha='center', va='center', fontsize=10)
    axes[2].text(5, .5, 'Schematic; cleavage is not modeled', ha='center', fontsize=9, color=GREY)
    save(fig, output, 'Fig2A_assembly_cartoon')

    fig, axes = plt.subplots(1, 2, figsize=(10, 3.6), constrained_layout=True)
    for ax in axes:
        cartoon_axes(ax)
    panel(axes[0], 'A', 'Local loading benefit')
    gene(axes[0], 1.2, pol=6.7, occupied=3)
    axes[0].text(5, 3.3, 'More sites allow more E to be carried', ha='center', fontsize=10)
    panel(axes[1], 'B', 'Shared factor availability')
    axes[1].add_patch(Ellipse((5, 2.8), 3.1, .9, facecolor='#e1f0ef', edgecolor=TEAL))
    axes[1].text(5, 2.8, 'Shared free E', ha='center', va='center', fontsize=11, color=TEAL)
    for x in (2, 5, 8):
        axes[1].plot([x-.8, x+.8], [.7, .7], color=INK, lw=1.3)
        axes[1].add_patch(Ellipse((x, .85), .7, .45, color=INK))
        axes[1].scatter([x-.25, x], [1.25, 1.35], s=55, color=TEAL)
        arrow(axes[1], (x, 1.55), (5+(x-5)*.3, 2.35), both=True)
    axes[1].text(5, .05, 'E bound on one gene is unavailable to others', ha='center', fontsize=9)
    save(fig, output, 'Fig4_resource_cartoon')

    fig, axes = plt.subplots(1, 2, figsize=(10, 3.1), constrained_layout=True)
    for ax in axes:
        cartoon_axes(ax)
    panel(axes[0], 'A', 'Short gene')
    gene(axes[0], 1.35, tss=4.4, pas=7.8, pol=7.1, occupied=0, available=1)
    axes[0].text(5, 3.35, 'Less opportunity to preload E', ha='center', fontsize=10)
    panel(axes[1], 'B', 'Long gene')
    gene(axes[1], 1.35, tss=.8, pas=7.8, pol=7.1, occupied=3)
    axes[1].text(5, 3.35, 'More opportunity to preload E', ha='center', fontsize=10)
    save(fig, output, 'Fig5A_gene_length_cartoon')

    fig, ax = plt.subplots(figsize=(9, 2.8), constrained_layout=True)
    cartoon_axes(ax)
    ax.plot([.6, 9.5], [1.3, 1.3], color=INK, lw=1.4)
    for x, label in ((2, 'Proximal PAS'), (8.2, 'Distal PAS')):
        ax.plot([x, x], [1.15, 1.55], color=ORANGE, lw=2)
        ax.text(x, .65, label, ha='center', fontsize=11)
    polymerase(ax, 5.6, 1.4, occupied=2)
    arrow(ax, (2, 3.15), (8.2, 3.15), label='Opportunity for proximal commitment')
    ax.text(5.1, 2.6, 'Inter-PAS distance, d', ha='center', fontsize=10)
    ax.text(5, .1, r'Approximation: proximal usage $\approx F_{\mathrm{commit}}(d)$', ha='center', fontsize=10)
    save(fig, output, 'Fig6A_APA_cartoon')


def plot_parameter_maps(data, output, p):
    """Render Figure 2C without changing numerical data or other figures."""
    summary = data['capacity']
    fig, axes = plt.subplots(2, 3, figsize=(11.4, 6.5), constrained_layout=True)
    names = [('kHd', 'kEd'), ('kc', 'kHd'), ('kc', 'kEd')]
    labels = {'kHd': r'$k_{Hd} = k_{Hoff}/k_{Hon}$',
              'kEd': r'$k_{Ed} = k_{Eoff}/k_{Eon}$ (molecules)', 'kc': r'$k_c$ (s$^{-1}$)'}
    baseline = {'kHd': p.kHoff/p.kHon, 'kEd': p.kEoff/p.kEon, 'kc': p.kc}
    for row, m in enumerate((1, 5)):
        for col, (xn, yn) in enumerate(names):
            ax = axes[row, col]
            x, y, z = data[f'axis_{xn}'], data[f'axis_{yn}'], data[f'map_{m}_{col}']
            mesh = ax.contourf(x, y, np.ma.masked_invalid(np.minimum(z, 2000)), levels=np.linspace(0, 2000, 21), cmap='turbo', corner_mask=False)
            # Keep unreached values masked in the data and contours, but fill
            # their regions with exactly the top color used for capped CAD50.
            ax.set_facecolor(mesh.cmap(mesh.norm(mesh.cvalues[-1])))
            for level in (200, 500):
                if np.nanmin(z) < level < np.nanmax(z):
                    contours = ax.contour(x, y, np.ma.masked_invalid(z), levels=[level], colors='#bd2026', linestyles='--', linewidths=1.1)
                    contour_labels = ax.clabel(contours, fmt='%d', fontsize=8)
                    for label in contour_labels:
                        lx, ly = label.get_position()
                        fx = np.log(lx/x[0])/np.log(x[-1]/x[0])
                        fy = np.log(ly/y[0])/np.log(y[-1]/y[0])
                        if not (.07 < fx < .93 and .07 < fy < .93):
                            label.remove()
            ax.plot(baseline[xn], baseline[yn], '*', ms=13, mfc='white', mec=INK)
            ax.set(xscale='log', yscale='log', xlabel=labels[xn], ylabel=labels[yn])
            ax.grid(False)
            ax.set_title(f'M = {m}; baseline {summary[m-1,1]:.0f} bp', fontsize=10)
    fig.colorbar(mesh, ax=axes, fraction=.025, pad=.02, label=r'CAD$_{50}$ (bp; color capped at 2000)')
    fig.supxlabel('Stars: baseline. Dashed contours: 200 and 500 bp. Top color also includes 50% not reached.', fontsize=9)
    save(fig, output, 'Fig2C_parameter_maps')


def plot_results(data, output, p):
    summary = data['capacity']
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8), constrained_layout=True)
    ax = axes[0]
    panel(ax, 'A', 'Binding capacity and commitment')
    ax.plot(summary[:8, 0], summary[:8, 1], '-o', color=TEAL, ms=5)
    ax.set(xlabel='Maximum E-binding capacity (M)', ylabel=r'CAD$_{50}$ (bp)', xticks=range(1, 9), ylim=(100, 850))
    for i in (0, 4):
        ax.annotate(f'{summary[i,1]:.0f} bp', summary[i, :2], xytext=(8, 10), textcoords='offset points', fontsize=9)
    ax = axes[1]
    panel(ax, 'B', 'CPA assembly and commitment')
    for m in (1, 5):
        cdf = data[f'cdf_M{m}']
        ax.plot(cdf[:, 0], 100*cdf[:, 1], color=COLORS[m], label=f'M = {m}', ls='--' if m == 1 else '-')
    ax.axhline(50, color=GREY, ls=':', lw=1)
    ax.set(xlabel='Distance downstream of PAS (bp)', ylabel='CPA commitment (%)', xlim=(0, 1000), ylim=(0, 100))
    ax.legend(loc='lower right')
    save(fig, output, 'Fig2B_capacity_and_commitment')

    plot_parameter_maps(data, output, p)

    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8), constrained_layout=True)
    ax = axes[0]
    panel(ax, 'A', 'Polymerases preloaded with E')
    ax.plot(summary[:, 0], 100*summary[:, 2], '-o', color=TEAL)
    for i in (0, 4, 9):
        ax.annotate(f'{100*summary[i,2]:.0f}%', (summary[i,0], 100*summary[i,2]), xytext=(0, 9), textcoords='offset points', ha='center', fontsize=9)
    ax.set(xlabel='Maximum E-binding capacity (M)', ylabel='Pol II carrying E before PAS (%)', ylim=(0, 90), xticks=range(1, 11))
    ax = axes[1]
    panel(ax, 'B', 'E and Ser2P loading profiles')
    for m in (1, 3, 5):
        profile = data[f'profile_M{m}'][:-1]  # strictly before recognition
        ax.plot(profile[:, 0], profile[:, 1], color=COLORS[m], label=f'M = {m}')
        ax.plot(profile[:, 0], profile[:, 2], '--', color=COLORS[m])
    ax.text(.03, .97, 'Solid: mean E\nDashed: mean Ser2P sites', transform=ax.transAxes, va='top', fontsize=9)
    ax.set(xlabel='Position relative to PAS (kb)', ylabel='Mean count per polymerase', xlim=(-20, 0), ylim=(0, 5), xticks=[-20, -15, -10, -5, 0], xticklabels=['TSS', '-15', '-10', '-5', 'PAS'])
    ax.legend(loc='center left')
    save(fig, output, 'Fig3_preloading')

    length = data['length_curves']
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.7), constrained_layout=True)
    ax = axes[0]
    panel(ax, 'B', 'CAD and gene length')
    ax.plot(length[:, 0]/1000, length[:, 4], '-o', color=TEAL, ms=4)
    ax.axvline(22, color=GREY, ls=':', lw=1)
    ax.text(24, 3000, 'Population\nmedian 22 kb', fontsize=9, color=GREY)
    ax.set(xscale='log', xlabel='TSS-to-PAS distance (kb)', ylabel=r'CAD$_{50}$ (bp)', ylim=(0, 4000), xlim=(2.5, 200))
    ax = axes[1]
    panel(ax, 'C', 'Multiple sites gain more from length')
    for m, column in ((1, 2), (5, 4)):
        ax.plot(length[:, 0]/1000, length[:, column], '-o', color=COLORS[m], label=f'M = {m}', ms=4)
    ax.set(xscale='log', yscale='log', xlabel='TSS-to-PAS distance (kb)', ylabel=r'CAD$_{50}$ (bp)', xlim=(2.5, 200), ylim=(200, 4500))
    ax.text(.97, .97, 'Same shared pools', transform=ax.transAxes, ha='right', va='top', fontsize=9)
    ax.legend(loc='upper right', bbox_to_anchor=(1, .87))
    ax.set_yticks([250, 500, 1000, 2000, 4000], labels=['250', '500', '1000', '2000', '4000'])
    for ax in axes[:2]:
        ax.set_xticks([2.5, 5, 10, 20, 50, 100, 200], labels=['2.5', '5', '10', '20', '50', '100', '200'])
        ax.minorticks_off()
    ax = axes[2]
    panel(ax, 'D', 'Short genes reach PAS with less E')
    for L, style in ((5000, '--'), (20000, '-')):
        profile = data[f'length_profile_{L}'][:-1]
        ax.plot(profile[:, 0], profile[:, 1], ls=style, color=TEAL, label=f'{L//1000} kb gene')
        ax.plot(profile[:, 0], profile[:, 2], ls=style, color=GREEN)
    ax.text(.03, .97, 'Green: Ser2P sites\nBlue: mean E', transform=ax.transAxes, va='top', fontsize=9)
    ax.set(xlabel='Position relative to PAS (kb)', ylabel='Mean count per polymerase', xlim=(-20, 0), ylim=(0, 5))
    ax.legend(loc='center left')
    save(fig, output, 'Fig5B_gene_length')


def plot_human_lengths(data, output, summary_path, gff_path):
    """Plot observed transcript spans; do not refit or rerun the population model.

    MANE summary coordinates are 1-based inclusive, verified against mRNA GFF
    features. Distance between endpoint bases is end-start (introns included).
    Transcript ends are a proxy for the model PAS, not annotated hexamer sites.
    """
    with gzip.open(summary_path, 'rt') as stream:
        rows = list(csv.DictReader(stream, delimiter='\t'))
    selected = [r for r in rows if r['MANE_status'] == 'MANE Select'
                and r['RefSeq_nuc'].startswith('NM_')
                and re.fullmatch(r'NC_0000(?:0[1-9]|1[0-9]|2[0-4])\.\d+', r['GRCh38_chr'])]
    if not selected or len({r['#NCBI_GeneID'] for r in selected}) != len(selected):
        raise ValueError('Expected one MANE Select transcript per retained gene.')
    lookup = {r['RefSeq_nuc']: r for r in selected}
    checked = set()
    with gzip.open(gff_path, 'rt') as stream:
        for line in stream:
            if line.startswith('#'):
                continue
            columns = line.rstrip().split('\t')
            if len(columns) != 9 or columns[2] != 'mRNA':
                continue
            attrs = dict(item.split('=', 1) for item in columns[8].split(';') if '=' in item)
            accession = attrs.get('transcript_id')
            if accession not in lookup:
                continue
            row = lookup[accession]
            chrom_number = int(row['GRCh38_chr'].split('.')[0].split('_')[1])
            chrom = 'chr' + {23: 'X', 24: 'Y'}.get(chrom_number, str(chrom_number))
            if (chrom, row['chr_start'], row['chr_end'], row['chr_strand']) != (columns[0], columns[3], columns[4], columns[6]):
                raise ValueError(f'MANE summary/GFF coordinate mismatch: {accession}')
            checked.add(accession)
    if checked != lookup.keys():
        raise ValueError('Some retained transcript coordinates were not verified in the GFF.')
    distances = np.array([int(r['chr_end'])-int(r['chr_start']) for r in selected])
    if np.any(distances <= 0):
        raise ValueError('Nonpositive transcript span.')
    csv_dir = output/'data'
    csv_dir.mkdir(exist_ok=True)
    with (csv_dir/'human_MANE_TSS_to_end.csv').open('w') as stream:
        writer = csv.writer(stream)
        writer.writerow(['gene_id', 'symbol', 'refseq_transcript', 'chromosome', 'strand',
                         'TSS_1based', 'transcript_3prime_end_1based', 'distance_bp'])
        for row, distance in zip(selected, distances):
            start, end = int(row['chr_start']), int(row['chr_end'])
            tss, end3 = (start, end) if row['chr_strand'] == '+' else (end, start)
            writer.writerow([row['#NCBI_GeneID'], row['symbol'], row['RefSeq_nuc'],
                             row['GRCh38_chr'], row['chr_strand'], tss, end3, distance])
    median = float(np.median(distances))
    bins = np.geomspace(30, 3e6, 46)
    counts, edges = np.histogram(distances, bins=bins)
    if counts.sum() != len(distances):
        raise ValueError('Histogram range must retain every selected transcript.')
    np.savetxt(csv_dir/'human_MANE_length_histogram.csv',
               np.column_stack((edges[:-1], edges[1:], counts)), delimiter=',',
               header='left_bp,right_bp,gene_count')
    fig, axes = plt.subplots(2, 1, figsize=(7.4, 6.2), sharex=True,
                             gridspec_kw={'height_ratios': [1, 1.2]}, constrained_layout=True)
    ax = axes[0]
    panel(ax, 'C', 'Commitment kinetics and human transcript spans')
    curves = data['length_curves']
    ax.plot(curves[:, 0]/1000, curves[:, 4], '-o', ms=3, color=TEAL,
            label=r'Model: M = 5; $k_{Hoff}=0.01$ s$^{-1}$')
    ax.set(ylabel=r'CAD$_{50}$ (bp)', ylim=(0, 4100))
    ax.legend(loc='upper right', fontsize=8)
    ax.text(.02, .94, 'Curve shown only over\nsimulated 2.5–200 kb',
            transform=ax.transAxes, va='top', fontsize=9, color=GREY)
    ax = axes[1]
    ax.bar(edges[:-1]/1000, counts, width=np.diff(edges)/1000, align='edge',
           color='#69a8cf', edgecolor='white', linewidth=.5)
    ax.set_title(f'Human MANE Select v1.5: {len(distances):,} genes', loc='left', fontsize=10)
    ax.set(ylabel='Genes per log-spaced bin', xlabel='Genomic distance (kb)\nModel: TSS to PAS; annotation: TSS to transcript 3′ end',
           xscale='log', xlim=(.03, 3000), ylim=(0, counts.max()*1.23))
    ax.text(.02, .94, f'Observed median: {median/1000:g} kb', transform=ax.transAxes,
            va='top', fontsize=9, color=PURPLE)
    for ax in axes:
        ax.axvspan(2.5, 200, color=TEAL, alpha=.06, zorder=0)
        ax.axvline(median/1000, color=PURPLE, ls='--', lw=1)
        ax.set_xticks([.1, 1, 10, 100, 1000], labels=['0.1', '1', '10', '100', '1,000'])
        ax.minorticks_off()
    save(fig, output, 'Fig5C_human_length_distribution')
    info = dict(release='MANE GRCh38 v1.5',
                source_url='https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.5/',
                documentation='https://www.ncbi.nlm.nih.gov/refseq/MANE/',
                source_files={p.name: dict(sha256=hashlib.sha256(p.read_bytes()).hexdigest(), bytes=p.stat().st_size)
                              for p in (summary_path, gff_path)},
                selection='MANE Select, NM_ transcripts, assembled chromosomes 1–22/X/Y; one transcript per gene; no length filtering.',
                input_rows=len(rows), retained_genes=len(selected), verified_mRNA_coordinates=len(checked),
                median_bp=median, min_bp=int(distances.min()), max_bp=int(distances.max()),
                fraction_below_20kb=float(np.mean(distances < 20000)),
                fraction_below_model_range=float(np.mean(distances < 2500)),
                fraction_above_model_range=float(np.mean(distances > 200000)),
                fraction_within_model_range=float(np.mean((distances >= 2500) & (distances <= 200000))),
                distance_definition='chr_end - chr_start; distance between 1-based endpoint bases, including introns.',
                limitations=['Transcript 3-prime ends proxy PAS positions; hexamer offsets are not resolved.',
                             'Not tissue-specific, expression-weighted, or a distribution over all alternative PASs.',
                             'Protein-coding selection does not guarantee polyadenylation for every retained transcript.',
                             'Model retains its original assumed lognormal population; no empirical refit.'])
    (output/'human_length_provenance.json').write_text(json.dumps(info, indent=2)+'\n')
    return info


def digitize_spacing_histogram(output):
    """Read approximate bar heights from the draft crop, never raw site data.

    Axis calibration is specific to the supplied PDF rendered at 220 dpi:
    x=264 at 100 bp, x=1002 at 1 Mb; y=475 at 0, y=189 at 400.
    Recover only bins needed for the 30--5000 bp display. Retain their original
    logarithmic widths and counts (no truncation renormalization).
    """
    source = output/'GENCODE_spacing_source.png'
    pixels = np.asarray(Image.open(source).convert('RGB'), dtype=float)
    if pixels.shape != (605, 1225, 3):
        raise ValueError('Histogram crop changed; recalibrate axes before digitizing.')
    r, g, b = pixels.transpose(2, 0, 1)
    mask = (r > 75) & (r < 125) & (g > 140) & (g < 190) & (b > 180) & (b < 225)
    mask[:65] = False
    mask[475:] = False
    columns = np.flatnonzero(mask.sum(axis=0) >= 2)
    runs = np.split(columns, np.flatnonzero(np.diff(columns) > 1)+1)[:32]
    centers = np.array([run.mean() for run in runs])
    if len(runs) != 32 or not np.all((np.diff(centers) > 12) & (np.diff(centers) < 16)):
        raise ValueError('Could not identify the expected consecutive histogram bars.')
    step, origin = np.polyfit(np.arange(len(centers)), centers, 1)
    x_edges = origin + step*(np.arange(len(centers)+1)-.5)
    edges = 10**(2+(x_edges-264)*4/(1002-264))
    heights = np.array([(475-np.where(mask[:, run])[0].min())*400/(475-189) for run in runs])
    values = np.column_stack((edges[:-1], edges[1:], heights))
    values = values[values[:, 0] < 5000]
    np.savetxt(output/'data/spacing_histogram_digitized.csv', values, delimiter=',',
               header='left_edge_bp_approx,right_edge_bp_approx,count_approx')
    info = dict(source=str(source), source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                method='Color-based bar digitization with calibrated log-distance and linear-count axes.',
                calibration_pixels=dict(x100bp=264, x1Mbp=1002, y0=475, y400=189),
                status='Approximate image-derived bin counts; raw observations unavailable.',
                display_interval_bp=[30, 5000], renormalized=False)
    return values, info


def plot_apa(data, output, histogram):
    fig, axes = plt.subplots(2, 1, figsize=(7.3, 6.5), sharex=True,
                             gridspec_kw={'height_ratios': [1.5, 1]}, constrained_layout=True)
    ax = axes[0]
    panel(ax, 'B', 'PAS strength and inter-PAS spacing')
    colors = ['#0072bd', '#d95319', '#edb120', '#7e2f8e', '#77ac30']
    spacing = data['apa_spacing']
    # Include the origin in linear interpolation, but never plot zero on log x.
    display_distances = np.unique(np.r_[30., np.geomspace(30, 5000, 300), spacing[1:, 0]])
    for j, kh in enumerate((.001, .1, 1, 5, 10)):
        ax.plot(display_distances, 100*np.interp(display_distances, spacing[:, 0], spacing[:, j+1]),
                color=colors[j], label=rf'$k_{{Hoff}}={kh:g}$')
    ax.set(ylabel='Proximal commitment (%)', ylim=(0, 100))
    ax.legend(loc='upper left', title=r'$k_{Hoff}$ (s$^{-1}$)', ncol=2, fontsize=8)
    ax = axes[1]
    ax.bar(histogram[:, 0], histogram[:, 2], width=histogram[:, 1]-histogram[:, 0],
           align='edge', color='#69a8cf', edgecolor='white', linewidth=.6)
    ax.set_title('GENCODE inter-PAS distances', loc='left', fontsize=10)
    ax.set(xscale='log', xlabel='Inter-PAS distance (bp)', ylabel='Number of segments',
           xlim=(30, 5000), ylim=(0, 600))
    ax.text(.02, .95, 'Approximate counts digitized from draft', transform=ax.transAxes,
            fontsize=8, va='top', color=GREY)
    ax.set_xticks([30, 100, 300, 1000, 1500, 5000], labels=['30', '100', '300', '1,000', '1,500', '5,000'])
    ax.minorticks_off()
    for ax in axes:
        for d in (300, 1500):
            ax.axvline(d, color=GREY, ls=':', lw=.9)
    save(fig, output, 'Fig6B_APA_curves')

    fig, ax = plt.subplots(figsize=(6.4, 3.9), constrained_layout=True)
    panel(ax, 'C', 'Proximal commitment at a strong PAS')
    curve = data['apa_strong_pas']
    ax.plot(curve[:, 0], 100*curve[:, 1], color=TEAL, label=r'$k_{Hoff}=0.01$ s$^{-1}$')
    ax.axvline(300, color=GREY, ls='--', lw=1.2, label='300 bp')
    ax.set(xlabel='Inter-PAS distance (bp)', ylabel='Proximal commitment (%)',
           xlim=(0, 1000), ylim=(0, 100))
    ax.legend(loc='lower right')
    save(fig, output, 'Fig6C_strong_PAS')


def extract_empirical_panels(draft, output):
    renderer = shutil.which('pdftoppm')
    if not renderer:
        raise RuntimeError('Install Poppler (pdftoppm) to reproduce the empirical PDF panels.')
    crops = [('Fig1_abundance', 1, (.12, .108, .925, .455)),
             ('GENCODE_spacing_source', 4, (.085, .42, .74, .67))]
    provenance = []
    with tempfile.TemporaryDirectory(prefix='cpad-draft-panels-') as folder:
        for name, page, box in crops:
            stem = Path(folder)/f'page{page}'
            subprocess.run([renderer, '-f', str(page), '-l', str(page), '-r', '220',
                            '-singlefile', '-png', str(draft), str(stem)], check=True, capture_output=True)
            with Image.open(stem.with_suffix('.png')) as image:
                w, h = image.size
                image.crop((round(box[0]*w), round(box[1]*h), round(box[2]*w), round(box[3]*h))).save(output/f'{name}.png')
            provenance.append(dict(figure=name, source=str(draft), page=page, normalized_crop=box,
                                   status='Original empirical panel cropped from supplied PDF; raw observations unavailable.'))
    return provenance


def signature(p, ranges, points):
    numerical_sources = [ROOT/'default_parameters.m', ROOT/'python/plot_ebinding_profile.py',
                         *sorted((ROOT/'python/cpadynamics').glob('*.py'))]
    return dict(parameters=asdict(p), ranges=ranges, map_points=points,
                source_hashes={str(f.relative_to(ROOT)): hashlib.sha256(f.read_bytes()).hexdigest() for f in numerical_sources},
                builder_hash=hashlib.sha256(''.join(inspect.getsource(f) for f in
                    (read_parameters, solve_case, gene_length_data, compute_data)).encode()).hexdigest(),
                strong_kHoff=STRONG_KHOFF)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--draft-pdf', type=Path, default=Path.home()/'Downloads/Results_draft.pdf')
    results_root = Path(os.environ.get('CPAD_RESULTS_ROOT', ROOT/'Results'))
    parser.add_argument('--output-dir', type=Path, default=results_root/'ResultsFramework/draft')
    parser.add_argument('--map-points', type=int, default=10)
    parser.add_argument('--reuse-data', action='store_true', help='Reuse numerical data only when its parameter/source signature matches')
    parser.add_argument('--mane-summary', type=Path, help='Local MANE v1.5 summary.gz; otherwise use output sources directory if present')
    parser.add_argument('--mane-gff', type=Path, help='Local MANE v1.5 RefSeq GFF.gz for transcript-coordinate verification')
    args = parser.parse_args()
    if not args.draft_pdf.is_file():
        parser.error('Supply the five-page Results_draft.pdf with --draft-pdf.')
    if args.map_points < 2:
        parser.error('--map-points must be at least 2')
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    p, ranges = read_parameters(ROOT/'default_parameters.m')
    key = signature(p, ranges, args.map_points)
    cache, metadata = output/'numerical_data.npz', output/'numerical_metadata.json'
    if args.reuse_data and cache.exists() and metadata.exists() and json.loads(metadata.read_text())['signature'] == key:
        with np.load(cache, allow_pickle=False) as saved:
            data = {name: saved[name] for name in saved.files}
        info = json.loads(metadata.read_text())['info']
        print('Reusing validated numerical data.', flush=True)
    else:
        data, info = compute_data(p, ranges, args.map_points)
    # One additional finite-rate solve; keeps the unchanged expensive grid cache usable.
    strong_diagnostics = []
    strong_p = replace(p, kHoff=STRONG_KHOFF, geneLength_bp=50000)
    strong_cdf, _, _, _ = solve_case(strong_p, 5, strong_diagnostics)
    distances = np.arange(0, 1001, p.L_a)
    data['apa_strong_pas'] = np.column_stack((distances, np.interp(distances, strong_cdf.distances_bp, strong_cdf.cdf)))
    info['strong_pas'] = dict(parameters=asdict(strong_p), M=5, diagnostics=strong_diagnostics)
    np.savez_compressed(cache, **data)
    metadata.write_text(json.dumps(dict(signature=key, info=info), indent=2, allow_nan=False)+'\n')
    csv_dir = output/'data'
    csv_dir.mkdir(exist_ok=True)
    headers = {'capacity': 'M,CAD50_bp,loaded_fraction_before_PAS,E_free,Pol_free',
               'length_curves': 'nominal_length_bp,lattice_length_bp,CAD50_M1,max_cdf_M1,CAD50_M5,max_cdf_M5',
               'apa_strong_pas': 'distance_bp,proximal_commitment',
               'apa_spacing': 'distance_bp,kHoff_0.001,kHoff_0.1,kHoff_1,kHoff_5,kHoff_10',
               'apa_E_supply': 'E_to_baseline,commitment_at_300bp,commitment_at_1500bp'}
    for name, array in data.items():
        header = headers.get(name, 'position_relative_to_PAS_kb,mean_E,mean_Ser2P_sites' if 'profile' in name
                             else 'distance_bp,commitment_cdf' if name.startswith('cdf_')
                             else 'Rows=y axis; columns=x axis; axes supplied separately' if name.startswith('map_')
                             else 'value')
        np.savetxt(csv_dir/f'{name}.csv', array, delimiter=',', header=header, comments='# ')
    configure_style()
    plot_cartoons(output)
    plot_results(data, output, p)
    empirical = extract_empirical_panels(args.draft_pdf.resolve(), output)
    histogram, histogram_info = digitize_spacing_histogram(output)
    plot_apa(data, output, histogram)
    mane_summary = args.mane_summary or output/'sources/MANE.GRCh38.v1.5.summary.txt.gz'
    mane_gff = args.mane_gff or output/'sources/MANE.GRCh38.v1.5.refseq_genomic.gff.gz'
    human_lengths = None
    if args.mane_summary or args.mane_gff or mane_summary.exists() or mane_gff.exists():
        human_lengths = plot_human_lengths(data, output, mane_summary, mane_gff)
    manifest = dict(parameters=asdict(p), ranges=ranges, empirical_panels=empirical,
                    human_lengths=human_lengths,
                    spacing_histogram=histogram_info,
                    source_pdf_sha256=hashlib.sha256(args.draft_pdf.read_bytes()).hexdigest(),
                    script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                    numerical_signature=key, **info,
                    figure_scope='Draft reproduction with cartoons, aligned spacing histogram and strong-PAS commitment detail.',
                    interpretation='kc is irreversible CPA commitment; cleavage/termination kinetics are not resolved.',
                    limitations=['Resources recycle at the irreversible exit, as in the existing solver.',
                                 'APA usage interpretation assumes a deadline at the distal PAS.',
                                 'Resource competition cartoon is conceptual; no new genome-wide demand sweep.',
                                 'Gene length reproduces coarse MATLAB interpolation, not a grid-converged prediction.',
                                 '200/500 bp map contours are reference levels, not an experimentally established CAD50 confidence interval.'],
                    versions=dict(numpy=np.__version__, scipy=scipy.__version__, matplotlib=matplotlib.__version__))
    (output/'manifest.json').write_text(json.dumps(manifest, indent=2, allow_nan=False)+'\n')
    for record in info['draft_comparisons']:
        print(f"{record['label']}: {record['computed']:.3f} bp; draft {record['draft_label']} bp; matches rounding: {record['matches_rounding']}", flush=True)
    print(f'Figures and data: {output}', flush=True)


if __name__ == '__main__':
    main()
