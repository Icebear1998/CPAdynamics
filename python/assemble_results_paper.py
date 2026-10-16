"""Assemble the Results draft: figures and placeholder values from raw data.

Reads the CSV tables written by RunResultsAnalyses.m and writes, under
Results/AssembleResultsPaper/:

  Fig2.png ... Fig5.png, FigS1.png, FigS2.png (300 dpi) and matching .pdf
  values.json   placeholder key -> formatted value, e.g. "R2.cad_M1": "888"
  values.md     the same as a table, with each value's definition

Figure 1 (abundance), the runway cartoon (Fig4A_original.png) and the
inter-PAS histogram (Fig5A_original.png) are the original images in
Docs/Results.pdf, not made here. The plots follow the style of the original MATLAB figures in that file:
boxed axes with a light grid, the teal/orange/purple series of
EBindingNumberVsCad.m, green Ser2P and blue mean-E profiles, turbo CAD maps
with red dashed contours, and red open circles for default values.

Run from the repository root after RunResultsAnalyses:
  python python/assemble_results_paper.py

Set CPAD_RESULTS_ROOT to read and write a different results folder. A
missing input section is reported and its placeholders are left unfilled.
"""
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, LogNorm
from matplotlib.ticker import FuncFormatter

REPO = Path(__file__).resolve().parents[1]
RESULTS = Path(os.environ.get('CPAD_RESULTS_ROOT', REPO / 'Results'))
RAW = RESULTS / 'RunResultsAnalyses'
OUT = RESULTS / 'AssembleResultsPaper'

AVOGADRO = 6.022e23
NUCLEUS_L = 3.74e-13
WINDOW = (200, 600)                       # CAD50 benchmark window (Chao et al. 1999, SV40 early PAS)
CHAO_POINTS = [(435, 25), (550, 50)]      # long antisense: 25% and 50% rescue
CHAO_END_BP = 800                         # rescue complete
STRONG_KHOFF = 0.01

# Original figure palette (EBindingNumberVsCad.m full_colors) and MATLAB lines.
TEAL, ORANGE, PURPLE = '#156f8a', '#cf6e17', '#8f438b'
SER2P_GREEN, AVGE_BLUE = '#1f9d3a', '#1f4fd6'
DEFAULT_RED = '#d62728'
MATLAB = ['#0072bd', '#d95319', '#edb120', '#7e2f8e', '#77ac30', '#4dbeee']
INK, INK2, GREY = '#1a1a1a', '#4d4d4d', '#9a9a9a'
GRID, BAND = '#dedede', '#ececec'
SEQ = LinearSegmentedColormap.from_list('teal_seq', ['#cfe6ee', '#5aa7c0', '#156f8a', '#0b3a4a'])

plt.rcParams.update({
    'font.family': ['Arial', 'DejaVu Sans'], 'font.size': 8, 'axes.titlesize': 8.5,
    'axes.labelsize': 8, 'xtick.labelsize': 7, 'ytick.labelsize': 7, 'legend.fontsize': 7,
    'axes.edgecolor': '#262626', 'axes.linewidth': 0.6, 'axes.labelcolor': INK,
    'text.color': INK, 'xtick.color': INK2, 'ytick.color': INK2,
    'xtick.direction': 'in', 'ytick.direction': 'in', 'xtick.major.size': 3, 'ytick.major.size': 3,
    'xtick.minor.size': 1.5, 'ytick.minor.size': 1.5,
    'axes.grid': True, 'axes.grid.which': 'major', 'grid.color': GRID, 'grid.linewidth': 0.5,
    'lines.linewidth': 1.6, 'lines.markersize': 4.5,
    'legend.frameon': True, 'legend.edgecolor': '#bfbfbf', 'legend.fancybox': False,
    'legend.framealpha': 1.0, 'figure.dpi': 100, 'savefig.dpi': 300,
    'axes.titleweight': 'bold', 'axes.axisbelow': True, 'figure.constrained_layout.use': True,
})

LABELS = {'kc': 'k$_c$', 'kPon_slope': 'Ser2P slope', 'kPoff': 'k$_{Poff}$',
          'E_total': 'E$_{total}$', 'kEoff': 'k$_{Eoff}$', 'kEon': 'k$_{Eon}$',
          'k_e2': 'k$_{e2}$', 'kHon': 'k$_{Hon}$', 'kHoff': 'k$_{Hoff}$',
          'k_e': 'k$_e$', 'Pol_total': 'Pol II$_{total}$'}

values, definitions, missing = {}, {}, []


# --------------------------------------------------------------------------
def put(key, value, definition):
    values[key] = value
    definitions[key] = definition


def read(section, name):
    path = RAW / section / name
    if not path.exists():
        raise FileNotFoundError(path)
    return pd.read_csv(path)


def appendix_range(name):
    """[min, max] of P.ranges.<name>, as exported by RunResultsAnalyses (Abundance/ranges.csv)."""
    r = read('Abundance', 'ranges.csv').set_index('param')
    return float(r.loc[name, 'min']), float(r.loc[name, 'max'])


def fmt_int(x, sig=None):
    if not np.isfinite(x):
        return 'n/a'
    if sig:
        x = float(f'{x:.{sig}g}')
    return f'{x:,.0f}'


def fmt(x):
    """Two significant figures below 10 (0.76, 2.0, 1.5), whole numbers above."""
    if not np.isfinite(x):
        return 'n/a'
    return f'{x:.2g}' if abs(x) < 1 else f'{x:.1f}' if abs(x) < 10 else f'{x:.0f}'


def pct(x):
    return 'n/a' if not np.isfinite(x) else f'{100 * x:.0f}%'


def log_interp(x, xs, ys):
    return float(np.interp(np.log(x), np.log(xs), ys))


def crossing(xs, ys, level):
    """First x where ys crosses level, linear in log x; NaN if none."""
    xs, ys = np.asarray(xs, float), np.asarray(ys, float)
    for i in range(len(xs) - 1):
        a, b = ys[i] - level, ys[i + 1] - level
        if np.isfinite(a) and np.isfinite(b) and a * b <= 0 and a != b:
            t = a / (a - b)
            return float(np.exp(np.log(xs[i]) + t * (np.log(xs[i + 1]) - np.log(xs[i]))))
    return float('nan')


def doubling_gain(frame):
    """Usage at 300 bp vs E/Pol and its change per doubling (grid steps are sqrt 2)."""
    f = frame.sort_values('E_to_Pol')
    r, v = f.E_to_Pol.to_numpy(), f.usage_300.to_numpy() * 100
    return r, v, r[1:-1], v[2:] - v[:-2]


def panel(ax, letter, title=None):
    """Bold panel letter at the top left, centred bold title as in the originals."""
    ax.text(-0.02, 1.02, letter, transform=ax.transAxes, fontsize=11, fontweight='bold',
            ha='right', va='bottom')
    if title:
        ax.set_title(title)


def plain_ticks(axis):
    axis.set_major_formatter(FuncFormatter(lambda v, _: f'{v:g}'))


def sci_k(axis):
    axis.set_major_formatter(FuncFormatter(lambda v, _: f'{v / 1e3:g}k' if v else '0'))


def note(ax, x, y, text, **kw):
    kw.setdefault('fontsize', 6.5)
    kw.setdefault('color', INK2)
    ax.text(x, y, text, **kw)


def window_band(ax, label=True, x=None, where='bottom'):
    ax.axhspan(*WINDOW, color=BAND, zorder=0, lw=0)
    if label and x is not None:
        y, va = (WINDOW[0] + 15, 'bottom') if where == 'bottom' else (WINDOW[1] - 15, 'top')
        note(ax, x, y, f'Chao et al. window\n{WINDOW[0]}–{WINDOW[1]} bp', va=va, ha='right')


def default_marker(ax, x, y, label='Default value'):
    ax.plot(x, y, 'o', ms=7, mfc='none', mec=DEFAULT_RED, mew=1.4, zorder=6, label=label)


def save(fig, name):
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / f'{name}.png', bbox_inches='tight')
    fig.savefig(OUT / f'{name}.pdf', bbox_inches='tight')
    plt.close(fig)


def section(fn):
    try:
        fn()
    except FileNotFoundError as exc:
        missing.append(f'{fn.__name__}: {exc}')
        print(f'[skip] {fn.__name__}: missing {exc}', file=sys.stderr)


def cdf_panel(ax, cdf, cad, Ms, title):
    """% CPA assembly completion vs Poly(A)-antiPoly(A) distance, as Figure 8 of Chao et al."""
    colors = {1: ORANGE, 2: PURPLE, 3: '#5aa7c0', 5: TEAL}
    ax.axhline(50, color=GREY, lw=0.7, ls=':')
    for m in Ms:
        c = cdf[(cdf.M == m) & (cdf.distance_bp <= 1000)]
        ax.plot(c.distance_bp, c.cdf * 100, color=colors.get(m, TEAL), lw=2 if m in (1, 5) else 1.2,
                label=f'M = {m}')
        x50 = cad.loc[m, 'CAD50']
        if x50 <= 1000:
            ax.plot([x50, x50], [0, 50], color=colors.get(m, TEAL), lw=0.7, ls='--')
    for i, (x, y) in enumerate(CHAO_POINTS):
        ax.plot(x, y, 's', ms=5.5, mfc='white', mec=INK, mew=1.1, zorder=6,
                label='Chao et al. (25%, 50%)' if i == 0 else None)
    ax.axvline(CHAO_END_BP, color=INK2, lw=0.8, ls='-.')
    note(ax, CHAO_END_BP - 15, 6, 'rescue\ncomplete', ha='right')
    ax.set_xlim(0, 1000); ax.set_ylim(0, 100)
    ax.set_xlabel('Poly(A)–antiPoly(A) distance (bp)')
    ax.set_ylabel('% CPA assembly completion')
    ax.legend(loc='upper left', fontsize=6.5)
    ax.set_title(title)


def profile_axes(ax, prof, key, keys_colors, rel_kb=True, show_after=True):
    """Mean bound E (solid) and Ser2P sites (dashed) vs position relative to the PAS."""
    for k, color, label in keys_colors:
        p = prof[prof[key] == k]
        pas = p.pas_node.iloc[0]
        q = p if show_after else p[p.node < pas]
        x = (q.position_bp - pas * 100) / 1e3 if rel_kb else q.position_bp / 1e3
        ax.plot(x, q.avg_E, color=color, lw=1.8, label=label)
        ax.plot(x, q.avg_Ser2P, color=color, lw=1.1, ls='--')
    ax.axvline(0, color=INK2, lw=0.8, ls=':')


# --------------------------------------------------------------------------
def abundance():
    copies = read('Abundance', 'copy_numbers.csv')
    pools = read('Abundance', 'pools.csv').iloc[0]
    to_uM = 1e6 / (AVOGADRO * NUCLEUS_L)
    cpsf160 = copies.loc[copies.protein == 'CPSF-160', 'copies'].item()
    put('R1.cpsf160', fmt_int(cpsf160), 'CPSF-160 copies per HeLa cell (Table G1)')
    put('R1.pol', fmt_int(pools.Pol_total), 'P.Pol_total')
    put('R1.ratio', fmt(pools.E_total / pools.Pol_total), 'P.E_total / P.Pol_total')
    put('R1.ratio_range', f'{pools.E_total_min / pools.Pol_total_max:.2f}–'
        f'{pools.E_total_max / pools.Pol_total_min:.0f}',
        'E_total min / Pol_total max to E_total max / Pol_total min (P.ranges)')
    put('R1.E_conc', f'{pools.E_total * to_uM:.2f}', 'E_total in a 374 µm³ nucleus (µM)')
    put('R1.Kd', f'{pools.kEoff / pools.kEon * to_uM:.2f}', 'kEoff/kEon in the same volume (µM)')


# --------------------------------------------------------------------------
def map_panel(ax, m, base, M, xn, yn):
    """One CAD map in the style of Sweep2DCad.m; returns the filled contour set."""
    names = {'kHd': 'k$_{Hd}$ = k$_{Hoff}$ / k$_{Hon}$', 'kEd': 'k$_{Ed}$ = k$_{Eoff}$ / k$_{Eon}$ (molecules)',
             'kc': 'k$_c$ (s$^{-1}$)'}
    d = m[(m.M == M) & (m.x_name == xn) & (m.y_name == yn)]
    xs, ys = np.sort(d.x.unique()), np.sort(d.y.unique())
    Z = d.pivot(index='y', columns='x', values='CAD50').loc[ys, xs].to_numpy()
    Zc = np.where(np.isfinite(Z), np.minimum(Z, 2000), 2000)
    cs = ax.contourf(xs, ys, Zc, levels=np.linspace(0, 2000, 21), cmap='turbo')
    lines = ax.contour(xs, ys, Zc, levels=list(WINDOW), colors=['#ff2a2a'], linewidths=1.3, linestyles='--')
    ax.clabel(lines, fmt=lambda v: f'{v:.0f}', fontsize=6.5, colors='#ff2a2a', inline=True)
    ax.plot(base[xn], base[yn], marker='*', ms=11, mfc='white', mec=INK, mew=0.8)
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.grid(False)
    ax.set_xlabel(names[xn]); ax.set_ylabel(names[yn])
    # Commitment faster than the window is not a failure, so only the upper edge counts.
    within = np.mean(np.nan_to_num(Z, nan=np.inf) <= WINDOW[1])
    ax.set_title(f'M = {M}: ≤ {WINDOW[1]} bp in {within * 100:.0f}% of plane', fontsize=7.5)
    return cs, within


def multivalency():
    """Figure 2: one E site misses the window at the measured E; five sites reach it."""
    weak = read('BaselineWeakPAS', 'cad_vs_M.csv').set_index('M')
    strong = read('BaselineStrongPAS', 'cad_vs_M.csv').set_index('M')
    cdf_w = read('BaselineWeakPAS', 'cdf_by_M.csv')
    esweep = read('BaselineWeakPAS', 'cad_vs_Etotal.csv')
    pools = read('Abundance', 'pools.csv').iloc[0]
    m = read('Maps2D', 'maps.csv')
    base = read('Maps2D', 'baseline_point.csv').iloc[0]

    for k in (1, 2, 3, 5, 10):
        put(f'R2.cad_M{k}', fmt_int(weak.loc[k, 'CAD50']), f'CAD50 (bp), weak PAS, M = {k}')
    for k in (1, 5, 10):
        put(f'R2.cad_strong_M{k}', fmt_int(strong.loc[k, 'CAD50']), f'CAD50 (bp), strong PAS, M = {k}')
    plateau_M = int(weak.index[(weak.CAD50 - weak.CAD50.iloc[-1]) / weak.CAD50.iloc[-1] < 0.01][0])
    put('R2.plateau_M', f'{plateau_M}', 'smallest M within 1% of the M = 10 CAD50, weak PAS')
    for k in (1, 5):
        c = cdf_w[cdf_w.M == k]
        put(f'R2.F550_M{k}', pct(np.interp(550, c.distance_bp, c.cdf)), f'commitment CDF at 550 bp, M = {k}')
        e = esweep[esweep.M == k].sort_values('E_to_Pol')
        put(f'R2.E_needed_M{k}', fmt(crossing(e.E_to_Pol, e.CAD50, WINDOW[1])),
            f'E_total/Pol_total at which CAD50 = {WINDOW[1]} bp (window edge), M = {k}')
        cad = m[m.M == k].CAD50.fillna(np.inf)
        put(f'R2.map_window_M{k}', pct(np.mean(cad <= WINDOW[1])),
            f'share of all map points with CAD50 <= {WINDOW[1]} bp, M = {k}')

    fig = plt.figure(figsize=(7.4, 7.7))
    rows = fig.add_gridspec(2, 1, height_ratios=[1, 1.85])
    top = rows[0].subgridspec(1, 3, width_ratios=[1.1, 1, 1])

    ax = fig.add_subplot(top[0])
    window_band(ax, x=10.3, where='top')
    for frame, color, label in ((weak, TEAL, f'Weak PAS (k$_{{Hoff}}$ = 1 s$^{{-1}}$)'),
                                (strong, ORANGE, f'Strong PAS (k$_{{Hoff}}$ = {STRONG_KHOFF:g} s$^{{-1}}$)')):
        ax.plot(frame.index, frame.CAD50, '-o', color=color, label=label)
    for k in (1, 2, 3, 5, 10):
        note(ax, k + (-0.1 if k == 10 else 0.15), weak.loc[k, 'CAD50'] + 18, f"{weak.loc[k, 'CAD50']:.0f}",
             color=TEAL, ha='right' if k == 10 else 'left')
    for k in (1, 5, 10):
        note(ax, k + (0.15 if k == 1 else 0), strong.loc[k, 'CAD50'] + (15 if k == 1 else -22),
             f"{strong.loc[k, 'CAD50']:.0f}", color=ORANGE, va='bottom' if k == 1 else 'top',
             ha={1: 'left', 5: 'center', 10: 'right'}[k])
    ax.set_xticks(range(1, 11))
    ax.set_xlim(0.5, 10.5); ax.set_ylim(150, 900)
    ax.set_xlabel('Maximum E binding capacity (M)')
    ax.set_ylabel('CAD$_{50}$ (bp)')
    ax.legend(loc='upper right', fontsize=6)
    panel(ax, 'A', 'CAD$_{50}$ vs E-binding capacity')

    ax = fig.add_subplot(top[1])
    cdf_panel(ax, cdf_w, weak, (1, 2, 3, 5), 'Assembly CDF, weak PAS')
    panel(ax, 'B')

    ax = fig.add_subplot(top[2])
    window_band(ax, label=False)
    note(ax, 2.7e4, WINDOW[0] + 15, f'{WINDOW[0]}–{WINDOW[1]} bp', va='bottom')
    for k, color in ((1, ORANGE), (5, TEAL)):
        e = esweep[esweep.M == k].sort_values('E_to_Pol')
        ax.plot(e.E_to_Pol * pools.Pol_total, e.CAD50, '-o', color=color, mfc='white', mew=1.1, ms=3.8,
                label=f'M = {k}')
        default_marker(ax, pools.E_total, e.loc[np.isclose(e.E_to_Pol, 1.0), 'CAD50'].item(),
                       'Measured E' if k == 5 else None)
    e1 = esweep[esweep.M == 1].sort_values('E_to_Pol')
    need = crossing(e1.E_to_Pol, e1.CAD50, WINDOW[1]) * pools.Pol_total
    ax.annotate(f'M = 1 needs {need / pools.E_total:.1f}×\nmeasured E for {WINDOW[1]} bp', xy=(need, WINDOW[1]),
                xytext=(2.2e5, 1500), fontsize=6.3, color=INK2, ha='center',
                arrowprops=dict(arrowstyle='->', color=INK2, lw=0.7))
    ax.set_xscale('log')
    ax.set_xticks([25e3, 1e5, 4e5]); sci_k(ax.xaxis)
    ax.set_xlabel('E$_{total}$ (molecules); Pol II$_{total}$ = 100k')
    ax.set_ylabel('CAD$_{50}$ (bp)')
    ax.legend(loc='upper right', fontsize=6)
    panel(ax, 'C', 'CAD vs E abundance')

    map_grid(fig, rows[1].subgridspec(2, 3), m, base, 'DE', 'R2.map')
    save(fig, 'Fig2')


def map_grid(fig, grid, m, base, letters, key):
    """2 x 3 CAD maps (rows M = 1, M = 5) with a shared colour bar."""
    pairs = [('kHd', 'kEd'), ('kc', 'kHd'), ('kc', 'kEd')]
    map_axes = []
    for r_i, M in enumerate((1, 5)):
        for c_i, (xn, yn) in enumerate(pairs):
            ax = fig.add_subplot(grid[r_i, c_i])
            cs, within = map_panel(ax, m, base, M, xn, yn)
            put(f'{key}_{xn}_{yn}_M{M}', pct(within), f'share of the {xn}-{yn} map with CAD50 <= {WINDOW[1]} bp, M = {M}')
            panel(ax, letters[r_i] + str(c_i + 1))
            map_axes.append(ax)
    cb = fig.colorbar(cs, ax=map_axes, fraction=0.025, pad=0.02, ticks=range(0, 2001, 400))
    cb.set_label('CAD$_{50}$ (bp; display capped at 2000)')


# --------------------------------------------------------------------------
def loading():
    """Figure 3: pre-loading at the PAS, profiles from TSS to PAS, where the E sits."""
    weak = read('BaselineWeakPAS', 'cad_vs_M.csv').set_index('M')
    prof = read('BaselineWeakPAS', 'profiles_by_M.csv')
    for k in (1, 5):
        put(f'R3.loaded_M{k}', pct(weak.loc[k, 'loaded_prePAS']),
            f'share of Pol II with >= 1 E at the node before the PAS, M = {k}')
        put(f'R3.free_M{k}', pct(weak.loc[k, 'free_E_fraction']), f'free E / E_total, M = {k}')
    put('R3.sites_M5', f"{weak.loc[5, 'sites_prePAS']:.1f}", 'mean phosphorylated sites before the PAS, M = 5')

    prof = prof[prof.node < prof.pas_node].copy()          # gene body only: TSS to PAS
    prof['x_kb'] = (prof.position_bp - prof.pas_node * 100) / 1e3
    wide = prof.pivot_table(index='x_kb', columns='M', values='avg_E')
    ok = wide[1] > 1e-3
    ratio = {k: (wide[k] / wide[1])[ok] for k in (3, 5)}
    r5 = ratio[5]
    cross_kb = r5.index[np.argmax(r5.to_numpy() >= 1)]
    early = r5[r5.index <= r5.index.min() + 2]                 # first 2 kb after loading starts
    body = lambda k, lo, hi: prof[(prof.M == k) & (prof.x_kb >= lo) & (prof.x_kb < hi)].avg_E.sum()
    share = {k: body(k, -5, 0) / body(k, -100, 0) for k in (1, 5)}
    put('R3.ratio_TSS_min', f'{early.min():.2f}', 'lowest M5/M1 bound-E ratio in the first 2 kb')
    put('R3.ratio_TSS_max', f'{early.max():.2f}', 'highest M5/M1 bound-E ratio in the first 2 kb')
    put('R3.ratio_PAS', f'{r5.iloc[-1]:.1f}', 'M5/M1 bound-E ratio at the node before the PAS')
    put('R3.cross_kb', f'{-cross_kb:.0f}', 'distance upstream of the PAS (kb) where M5 first exceeds M1')
    for k in (1, 5):
        put(f'R3.share_last5kb_M{k}', pct(share[k]), f'share of gene-body bound E in the last 5 kb, M = {k}')
        put(f'R3.E_PAS_M{k}', f"{prof[prof.M == k].avg_E.iloc[-1]:.2f}", f'mean bound E before the PAS, M = {k}')

    fig, axes = plt.subplots(1, 3, figsize=(7.4, 2.85), gridspec_kw={'width_ratios': [0.95, 1.1, 1.1]})
    ax = axes[0]
    ax.plot(weak.index, weak.loaded_prePAS * 100, '-o', color=TEAL)
    for k in (1, 5, 10):
        note(ax, k, weak.loc[k, 'loaded_prePAS'] * 100 + 2.2, pct(weak.loc[k, 'loaded_prePAS']),
             ha='center', color=TEAL)
    ax.set_xticks(range(1, 11)); ax.set_xlim(0.5, 10.5)
    ax.set_ylim(0, max(60, float(weak.loaded_prePAS.max()) * 100 * 1.15))
    ax.set_xlabel('Maximum E binding capacity (M)')
    ax.set_ylabel('Pol II carrying ≥ 1 E just before PAS (%)')
    panel(ax, 'A', 'Polymerases pre-loaded with E')

    ax = axes[1]
    for k, color in ((1, ORANGE), (3, PURPLE), (5, TEAL)):
        p = prof[prof.M == k]
        ax.plot(p.x_kb, p.avg_E, color=color, lw=1.8, label=f'M = {k}')
        ax.plot(p.x_kb, p.avg_Ser2P, color=color, lw=1.1, ls='--')
    ax.set_xlim(-20, 0)
    ax.set_xticks([-20, -15, -10, -5, 0], ['TSS', '−15', '−10', '−5', 'PAS'])
    ax.set_xlabel('Position relative to PAS (kb)')
    ax.set_ylabel('Mean per polymerase')
    ax.legend(loc='upper left', title='solid: AverageE\ndashed: Ser2P sites', title_fontsize=6)
    panel(ax, 'B', f"Ser2P sites at PAS: {weak.loc[5, 'sites_prePAS']:.1f} of 5 (M = 5)")

    ax = axes[2]
    ax.axhline(1, color=INK2, lw=0.8)
    ax.axhspan(0, 1, color=BAND, lw=0, zorder=0)
    for k, color in ((3, PURPLE), (5, TEAL)):
        ax.plot(ratio[k].index, ratio[k], color=color, lw=1.8, label=f'M = {k} / M = 1')
    ax.axvline(cross_kb, color=GREY, lw=0.7, ls=':')
    top = max(3.0, float(r5.max()) * 1.12)
    note(ax, cross_kb + 0.4, 0.45, f'crosses 1 at\n{-cross_kb:.0f} kb before PAS', va='bottom')
    note(ax, -19.6, top * 0.78, f'near TSS {early.min():.2f}–{early.max():.2f}×:\nfree E {weak.loc[5, "free_E_fraction"] * 100:.0f}% (M = 5)\n'
         f'vs {weak.loc[1, "free_E_fraction"] * 100:.0f}% (M = 1)', va='top')
    note(ax, -0.3, r5.iloc[-1] + 0.06, f'{r5.iloc[-1]:.1f}× at PAS', ha='right', va='bottom', color=TEAL)
    note(ax, -19.6, top * 0.97, f'share of bound E in last 5 kb:\nM = 5 {share[5] * 100:.0f}%, M = 1 {share[1] * 100:.0f}%',
         va='top', color=INK)
    ax.set_xlim(-20, 0); ax.set_ylim(0.4, top)
    ax.set_xticks([-20, -15, -10, -5, 0], ['TSS', '−15', '−10', '−5', 'PAS'])
    ax.set_xlabel('Position relative to PAS (kb)')
    ax.set_ylabel('Bound E relative to M = 1')
    ax.legend(loc='lower right', fontsize=6.3, bbox_to_anchor=(1.0, 0.1))
    panel(ax, 'C', 'Five sites move E toward the PAS')
    save(fig, 'Fig3')


# --------------------------------------------------------------------------
def gene_length():
    """Figure 4B-D (genome-wide CAD vs length, M = 5 vs M = 1, profiles) and Figure S1 (range ends)."""
    g = read('GeneLength', 'genome_cad_vs_length.csv').sort_values('L_bp')
    pools = read('GeneLength', 'pools.csv').iloc[0]
    w = read('GeneLength', 'length_weights.csv')
    d = read('GeneLength', 'designs_vs_length.csv')
    x = read('GeneLength', 'extremes_vs_length.csv')
    sl = read('GeneLength', 'profiles_short_long.csv')

    put('R4.n_genes', fmt_int(pools.num_genes), 'active genes sharing the pools')
    put('R4.cad_2p5kb', fmt_int(g.CAD50.iloc[0]), 'genome-wide CAD50 at 2.5 kb, strong PAS')
    put('R4.cad_5kb', fmt_int(log_interp(5000, g.L_bp, g.CAD50)), 'genome-wide CAD50 at 5 kb')
    put('R4.cad_20kb', fmt_int(log_interp(20000, g.L_bp, g.CAD50)), 'genome-wide CAD50 at 20 kb')
    plateau = g.CAD50.iloc[-1]
    plateau_L = g[g.CAD50 <= plateau * 1.1].L_bp.iloc[0]
    put('R4.cad_plateau', fmt_int(plateau), 'genome-wide CAD50 at 200 kb')
    put('R4.plateau_kb', f'{plateau_L / 1e3:.0f}', 'shortest simulated length within 10% of the 200 kb value (kb)')
    put('R4.free_pol', fmt_int(pools.R_free), 'shared free Pol II')
    put('R4.free_E_frac', pct(pools.E_free / pools.E_total), 'shared free E / E_total')
    curve = {k: v.sort_values('L_bp') for k, v in d.groupby('design')}
    at = lambda k, L: log_interp(L, curve[k].L_bp, curve[k].CAD50)
    for k in ('M5', 'M1'):
        for L, tag in ((2500, '2p5kb'), (5000, '5kb'), (40000, '40kb'), (200000, '200kb')):
            put(f'R4.{k}_{tag}', fmt_int(at(k, L)), f'CAD50 at {tag} at the shared pools, {k}')
    span = {k: at(k, 5000) / at(k, 40000) for k in ('M5', 'M1')}
    for k in span:
        put(f'R4.span_{k}', fmt(span[k]), f'CAD50(5 kb) / CAD50(40 kb), {k}')
    pre = lambda L: sl[(sl.L_bp == L) & (sl.node == sl.pas_node - 1)].iloc[0]
    for L, tag in ((5000, '5kb'), (20000, '20kb')):
        put(f'R4.E_at_PAS_{tag}', f'{pre(L).avg_E:.2f}', f'mean bound E just before the PAS, {tag} gene')
        put(f'R4.sites_at_PAS_{tag}', f'{pre(L).avg_Ser2P:.2f}', f'mean Ser2P sites just before the PAS, {tag} gene')

    fig = plt.figure(figsize=(7.4, 3.0))
    cols = fig.add_gridspec(1, 3, width_ratios=[1.15, 1, 1])
    strip_ax, main_ax = cols[0].subgridspec(2, 1, height_ratios=[0.17, 1], hspace=0.0)
    ax0 = fig.add_subplot(strip_ax)
    ax = fig.add_subplot(main_ax, sharex=ax0)
    ax0.fill_between(w.L_bp / 1e3, w.weight / w.weight.max(), color='#cfcfcf', lw=0)
    ax0.set_yticks([]); ax0.grid(False); ax0.tick_params(labelbottom=False)
    note(ax0, 2.7, 0.5, 'human gene-length distribution', va='center')
    panel(ax0, 'B', 'CAD vs gene length')
    ax.plot(g.L_bp / 1e3, g.CAD50, '-o', color=TEAL, ms=3.2, label='M = 5, 50% threshold')
    ax.axvline(22, color=INK2, lw=0.8, ls=':')
    note(ax, 23.5, g.CAD50.max() * 0.93, 'median\n22 kb')
    ax.axhline(plateau, color=GREY, lw=0.7, ls='--')
    note(ax, 195, plateau + 520, f'plateau ≈ {plateau:.0f} bp (≥ {plateau_L / 1e3:.0f} kb)', ha='right')
    ax.annotate(f'{g.CAD50.iloc[0]:,.0f} bp\nat 2.5 kb', xy=(2.5, g.CAD50.iloc[0]), xytext=(4.0, g.CAD50.iloc[0] * 0.93),
                fontsize=6.3, color=INK2, va='top', arrowprops=dict(arrowstyle='->', color=INK2, lw=0.6))
    ax.set_xscale('log'); ax.set_xlim(2.5, 200); plain_ticks(ax.xaxis)
    ax.set_xlabel('TSS-to-PAS distance (kb)')
    ax.set_ylabel('CAD$_{50}$ (bp)')
    ax.legend(loc='upper right', bbox_to_anchor=(1.0, 0.84), fontsize=6)

    ax = fig.add_subplot(cols[1])
    for k, color, label in (('M5', TEAL, 'M = 5'), ('M1', ORANGE, 'M = 1')):
        c = curve[k]
        ax.plot(c.L_bp / 1e3, c.CAD50, '-o', ms=2.5, color=color, label=label)
        ax.plot([5, 40], [at(k, 5000), at(k, 40000)], 'o', ms=5.5, mfc='white', mec=color, mew=1.1, zorder=6)
    ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlim(2.5, 200)
    ax.set_yticks([300, 500, 1000, 2000, 4000])
    plain_ticks(ax.xaxis); plain_ticks(ax.yaxis)
    ax.set_xlabel('TSS-to-PAS distance (kb)')
    ax.set_ylabel('CAD$_{50}$ (bp)')
    ax.legend(loc='upper right', fontsize=6.3, title='same shared pools', title_fontsize=6)
    note(ax, 18, 330, f'5 → 40 kb (circles):\nM = 5 falls {span["M5"]:.1f}×\nM = 1 falls {span["M1"]:.1f}×',
         color=INK, va='bottom')
    panel(ax, 'C', 'Five sites gain more from length')

    ax = fig.add_subplot(cols[2])
    genes = ((5000, '--', '5 kb gene'), (20000, '-', '20 kb gene'))
    body = sl[sl.node < sl.pas_node]                          # TSS to PAS only
    for L, style, _ in genes:
        p = body[body.L_bp == L]
        xk = (p.position_bp - p.pas_node.iloc[0] * 100) / 1e3
        ax.plot(xk, p.avg_Ser2P, color=SER2P_GREEN, lw=1.8, ls=style)
        ax.plot(xk, p.avg_E, color=AVGE_BLUE, lw=1.8, ls=style)
        ax.plot([xk.iloc[-1]] * 2, [pre(L).avg_Ser2P, pre(L).avg_E], 'o', color=INK, ms=3.5, zorder=6,
                mfc='white' if L == 5000 else INK)
    ax.axvline(-5, color=GREY, lw=0.7, ls=':')
    note(ax, -4.8, body.avg_Ser2P.max() * 0.45, '5 kb gene\nTSS', va='top')
    ymax = body.avg_Ser2P.max() * 1.1
    ax.set_xlim(-20, 0); ax.set_ylim(0, ymax)
    ax.set_xticks([-20, -15, -10, -5, 0], ['−20', '−15', '−10', '−5', 'PAS'])
    ax.set_xlabel('Position relative to PAS (kb)')
    ax.set_ylabel('Mean per polymerase')
    from matplotlib.lines import Line2D
    handles = [Line2D([], [], color=SER2P_GREEN, lw=1.8, label='Ser2P'),
               Line2D([], [], color=AVGE_BLUE, lw=1.8, label='AverageE')] + \
              [Line2D([], [], color=INK2, lw=1.4, ls=style, label=label) for _, style, label in genes]
    ax.legend(handles=handles, loc='upper left', fontsize=6)
    s5, s20 = pre(5000), pre(20000)
    note(ax, -13.2, ymax * 0.27,
         f'At the PAS\n5 kb:  {s5.avg_Ser2P:.2f} sites, {s5.avg_E:.2f} E\n'
         f'20 kb: {s20.avg_Ser2P:.2f} sites, {s20.avg_E:.2f} E', va='top', color=INK,
         bbox=dict(boxstyle='square,pad=0.35', fc='white', ec='#bfbfbf', lw=0.6))
    panel(ax, 'D', f'Short genes arrive with {s20.avg_E / s5.avg_E:.0f}× less E')
    save(fig, 'Fig4')

    # Figure S1: five-site curves at the range ends of the runway parameters.
    params = ['kPon_slope', 'kPoff', 'k_e', 'kc']
    fig, axes = plt.subplots(2, 2, figsize=(5.6, 4.4), sharex=True, sharey=True)
    base = curve['M5']
    for i, p in enumerate(params):
        ax = axes.flat[i]
        ax.plot(base.L_bp / 1e3, base.CAD50, color=TEAL, lw=1.6, label='baseline')
        for lvl, style, color in (('min', ':', ORANGE), ('max', '--', PURPLE)):
            c = x[(x.param == p) & (x.level == lvl)].sort_values('L_bp')
            ax.plot(c.L_bp / 1e3, c.CAD50, color=color, lw=1.3, ls=style,
                    label=f'range {lvl} ({c.value.iloc[0]:.3g})')
        ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlim(2.5, 200)
        plain_ticks(ax.xaxis); plain_ticks(ax.yaxis)
        ax.legend(loc='upper right', fontsize=5.8)
        if i % 2 == 0:
            ax.set_ylabel('CAD$_{50}$ (bp)')
        if i >= 2:
            ax.set_xlabel('TSS-to-PAS distance (kb)')
        panel(ax, 'ABCD'[i], LABELS[p])
    save(fig, 'FigS1')


# --------------------------------------------------------------------------
def apa():
    """Figure 5B-E (usage vs spacing, E abundance and E affinity) and Figure S2 (range ends)."""
    cdf = read('Apa', 'cdf_by_kHoff.csv')
    u = read('Apa', 'usage_vs_Etotal.csv')
    b = read('ApaBinding', 'usage_vs_Etotal_by_kEoff.csv')
    h = read('ApaBinding', 'pas_strength_by_kEoff.csv')
    pools = read('Abundance', 'pools.csv').iloc[0]
    measured = pools.E_total / pools.Pol_total
    base_k = pools.kEoff
    at = lambda xs, ys, x0: float(np.interp(np.log(x0), np.log(xs), ys))

    r, v, mid, g = doubling_gain(u[u.param == 'baseline'])
    put('R5.usage_300', pct(at(r, v, measured) / 100), 'proximal usage at 300 bp, measured E/Pol, strong PAS')
    put('R5.gain_measured', f'{at(mid, g, measured):.0f}', 'points per doubling of E at the measured ratio')
    put('R5.gain_excess', f'{at(mid, g, np.sqrt(4 * 8)):.1f}', 'points per doubling between 4x and 8x E/Pol')
    put('R5.ceiling', pct(v.max() / 100), 'usage at the highest E/Pol tested (16x)')
    put('R5.drop_half', f'{at(r, v, measured) - at(r, v, measured / 2):.0f}', 'usage lost when E is halved (points)')
    put('R5.rise_double', f'{at(r, v, 2 * measured) - at(r, v, measured):.0f}', 'usage gained when E is doubled (points)')
    for k in sorted(cdf.kHoff.unique()):
        c = cdf[cdf.kHoff == k]
        put(f'R5.usage300_kHoff{k:g}', pct(np.interp(300, c.distance_bp, c.cdf)),
            f'proximal usage at 300 bp, kHoff = {k:g}')

    khoff_rng, keoff_rng = appendix_range('kHoff'), appendix_range('kEoff')

    def affinity_table(b, h):
        """Gain at the measured ratio and PAS-strength discrimination for each kEoff design."""
        rows = []
        b = b.dropna(subset=['kEoff'])            # an unreachable match has no kEoff
        for (design, M, k), f in b.groupby(['design', 'M', 'kEoff']):
            rr, vv, mm, gg = doubling_gain(f)
            hh = h[(h.design == design) & (h.M == M) & np.isclose(h.kEoff, k)].set_index('kHoff').usage_300
            f1 = f[np.isclose(f.E_to_Pol, measured)].iloc[0]
            rows.append(dict(design=design, M=M, kEoff=k, gain=at(mm, gg, measured), usage=at(rr, vv, measured),
                             discrim=100 * (hh.loc[STRONG_KHOFF] - hh.loc[1.0]), free=f1.free_E_fraction,
                             CAD50=f1.CAD50, curve=(rr, vv)))
        aff = pd.DataFrame(rows)
        sweep = {M: aff[(aff.design == 'sweep') & (aff.M == M)].sort_values('kEoff') for M in (1, 5)}
        match = aff[aff.design == 'M1_matched']
        return sweep, (match.iloc[0] if len(match) else None)

    sweep, matched = affinity_table(b, h)
    s5, s1 = sweep[5], sweep[1]
    b5 = s5[np.isclose(s5.kEoff, base_k)].iloc[0]
    tight5 = s5.iloc[0]
    tight1 = s1.iloc[0]
    peak5 = s5.loc[s5.gain.idxmax()]
    put('R5.gain_peak_M5', f'{peak5.gain:.0f}', 'highest gain over the kEoff sweep, M = 5')
    put('R5.gain_peak_kEoff', f'{peak5.kEoff:.2g}', 'kEoff of the highest gain, M = 5')
    put('R5.kEoff_tightest', f'{tight5.kEoff:g}', 'tightest kEoff in the sweep')
    put('R5.fold_tightest', f'{base_k / tight5.kEoff:.0f}', 'baseline kEoff / tightest kEoff')
    designs = [('base_M5', b5), ('tight_M5', tight5), ('tight_M1', tight1)]
    if matched is not None:
        designs.append(('matched_M1', matched))
    for tag, row in designs:
        put(f'R5.gain_{tag}', f'{row.gain:.0f}' if row.gain >= 10 else f'{row.gain:.1f}',
            f'points per doubling of E at the measured ratio, {tag}')
        put(f'R5.usage_{tag}', pct(row.usage / 100), f'usage at 300 bp, measured E, {tag}')
        put(f'R5.discrim_{tag}', f'{row.discrim:.0f}', f'usage(kHoff {STRONG_KHOFF:g}) - usage(kHoff 1), points, {tag}')
        put(f'R5.free_{tag}', pct(row.free), f'free E fraction at the measured E, {tag}')
    for tag, row in (('weak_M5', s5.iloc[-1]), ('weak_M1', s1.iloc[-1])):
        put(f'R5.gain_{tag}', f'{row.gain:.0f}' if row.gain >= 10 else f'{row.gain:.1f}',
            f'points per doubling of E at the measured ratio, {tag} (kEoff = {row.kEoff:g})')
        put(f'R5.discrim_{tag}', f'{row.discrim:.0f}', f'strong - weak PAS usage, points, {tag}')
        put(f'R5.CAD_{tag}', fmt_int(row.CAD50), f'CAD50 at the measured E, {tag}')
    put('R5.CAD_base_M5', fmt_int(b5.CAD50), 'CAD50 at the measured E, strong PAS, M = 5')
    put('R5.CAD_tight_M1', fmt_int(tight1.CAD50), f'CAD50 at the measured E, M = 1, kEoff = {tight1.kEoff:g}')
    if matched is not None:
        put('R5.kEoff_matched', f'{matched.kEoff:.3f}', 'kEoff at which one site matches the M = 5 CAD50, strong PAS')
        put('R5.matched_fold', f'{base_k / matched.kEoff:.0f}', 'baseline kEoff / matched kEoff, strong PAS')
    else:
        put('R5.kEoff_matched', 'none', 'no kEoff within 1e5-fold lets one site match the M = 5 CAD50, strong PAS')

    fig, axes = plt.subplots(2, 2, figsize=(7.4, 5.6))
    ax = axes[0, 0]
    for i, k in enumerate(sorted(cdf.kHoff.unique())):
        c = cdf[(cdf.kHoff == k) & (cdf.distance_bp > 0)]
        inside = khoff_rng[0] <= k <= khoff_rng[1]
        ax.plot(c.distance_bp, c.cdf * 100, color=MATLAB[i], lw=1.9 if k == STRONG_KHOFF else 1.3,
                ls='-' if (inside or k < khoff_rng[0]) else '--', label=f'kHoff = {k:g}')
    ax.axvline(300, color=INK2, lw=0.8, ls=':')
    note(ax, 320, 3, '300 bp', color=INK2)
    ax.set_xscale('log'); ax.set_xlim(100, 30000); plain_ticks(ax.xaxis)
    ax.set_ylim(0, 102)
    ax.set_xlabel('Inter-PAS distance (bp, log scale)')
    ax.set_ylabel('Proximal site usage (%)')
    ax.legend(loc='lower right', fontsize=6, title='dashed: above appendix range', title_fontsize=5.5)
    panel(ax, 'B', 'Usage vs inter-PAS distance')

    ax = axes[0, 1]
    ax.axhline(v.max(), color=GREY, lw=0.7, ls='--')
    note(ax, 0.13, v.max() + 2, f'ceiling ≈ {v.max():.0f}% (set by k$_c$)')
    shown = [(s5.iloc[(s5.kEoff - 0.005).abs().argmin()], '#0b3a4a', ':'),
             (s5.iloc[(s5.kEoff - 0.05).abs().argmin()], '#5aa7c0', '--')]
    for row, color, style in shown:
        rr, vv = row.curve
        ax.plot(rr, vv, color=color, lw=1.3, ls=style, label=f'k$_{{Eoff}}$ = {row.kEoff:.2g} (tighter)')
    ax.plot(r, v, '-o', color=TEAL, ms=3.5, label=f'k$_{{Eoff}}$ = {base_k:g} (baseline)')
    default_marker(ax, measured, at(r, v, measured), 'Measured E')
    ax.set_xscale('log', base=2); plain_ticks(ax.xaxis)
    ax.set_ylim(0, 100)
    ax.set_xlabel('E$_{total}$ / Pol II$_{total}$')
    ax.set_ylabel('Proximal usage at 300 bp (%)')
    ax.legend(loc='lower right', fontsize=6)
    panel(ax, 'C', 'Usage at 300 bp vs E abundance')

    def affinity_axes(ax, col, ylabel, title, letter):
        ax.axvspan(*keoff_rng, color='#d9ecf2', lw=0, zorder=0)
        for M, color in ((5, TEAL), (1, ORANGE)):
            ax.plot(sweep[M].kEoff, sweep[M][col], '-o', ms=3.2, color=color, label=f'M = {M}')
        if matched is not None:
            ax.plot(matched.kEoff, matched[col], 'D', ms=5.5, mfc='white', mec=INK, mew=1.1, zorder=6,
                    label=f'M = 1 matched to M = 5 CAD')
        default_marker(ax, base_k, b5[col], 'Baseline')
        ax.set_xscale('log'); plain_ticks(ax.xaxis)
        ax.set_xlabel('k$_{Eoff}$ (s$^{-1}$); tighter binding ←')
        ax.set_ylabel(ylabel)
        panel(ax, letter, title)

    ax = axes[1, 0]
    affinity_axes(ax, 'gain', 'Change per doubling of E (points)', 'Tighter binding dulls the E response', 'D')
    ax.set_ylim(0, max(18, float(sweep[1].gain.max()) * 1.15))
    note(ax, keoff_rng[0] * 1.1, ax.get_ylim()[1] * 0.97, 'measured\nrange', color=TEAL, va='top')
    ax.legend(loc='lower right', fontsize=6)

    ax = axes[1, 1]
    affinity_axes(ax, 'discrim', 'Usage, strong − weak PAS (points)', 'Tighter binding blurs PAS strength', 'E')
    top = float(pd.concat([sweep[M].discrim for M in (1, 5)]).max()) * 1.15
    ax.set_ylim(0, top)
    note(ax, 0.0011, top * 0.96, f'strong PAS: k$_{{Hoff}}$ = {STRONG_KHOFF:g}; weak: k$_{{Hoff}}$ = 1 s$^{{-1}}$', va='top')
    ax.legend(loc='upper left', bbox_to_anchor=(0.0, 0.92), fontsize=6)
    save(fig, 'Fig5')

    # Figure S2: sensitivity of usage to E at each parameter range end.
    case_gain, case_peak, slow_kc = [], [], {}
    fig, ax = plt.subplots(figsize=(3.8, 2.9))
    for (p, lvl), f in u.groupby(['param', 'level']):
        rr, vv, mm, gg = doubling_gain(f)
        if not np.all(np.isfinite(gg)) or p == 'baseline':
            continue
        if p == 'kc' and lvl == 'min':
            slow_kc = {'gain': at(mm, gg, measured), 'usage': at(rr, vv, measured)}
        else:
            case_gain.append(at(mm, gg, measured))
        case_peak.append(mm[np.argmax(gg)])
        ax.plot(mm, gg, color=GREY, lw=0.7, alpha=0.8)
    put('R5.gain_min', f'{min(case_gain):.0f}', 'lowest gain at the measured ratio across range ends, excluding kc min')
    put('R5.gain_max', f'{max(case_gain):.0f}', 'highest gain at the measured ratio across range ends')
    put('R5.peak_min', fmt(min(case_peak)), 'lowest E/Pol of maximal gain across range ends')
    put('R5.peak_max', fmt(max(case_peak)), 'highest E/Pol of maximal gain across range ends')
    if slow_kc:
        put('R5.gain_kc_min', f"{slow_kc['gain']:.0f}", 'gain at the measured ratio, kc at its minimum')
        put('R5.usage_kc_min', pct(slow_kc['usage'] / 100), 'usage at 300 bp, kc at its minimum')
        note(ax, 0.16, slow_kc['gain'] + 0.6, 'k$_c$ at range min')
    ax.plot([], [], color=GREY, lw=0.7, label='9 parameters at range ends')
    ax.plot(mid, g, '-o', color=TEAL, lw=2, ms=3, label='baseline')
    ax.axvline(measured, color=DEFAULT_RED, lw=0.8, ls=':')
    note(ax, measured * 1.08, 1.0, 'measured', color=DEFAULT_RED)
    ax.set_xscale('log', base=2); plain_ticks(ax.xaxis)
    ax.set_ylim(0, 28)
    ax.set_xlabel('E$_{total}$ / Pol II$_{total}$')
    ax.set_ylabel('Change per doubling of E (points)')
    ax.legend(loc='upper right', fontsize=6)
    ax.set_title('Sensitivity of usage to E at range ends')
    save(fig, 'FigS2')


# --------------------------------------------------------------------------
def main():
    for fn in (abundance, multivalency, loading, gene_length, apa):
        section(fn)
    OUT.mkdir(parents=True, exist_ok=True)
    # Outputs of the earlier layout (old Figures 4 and 6, the global-sample figure).
    for stale in ('Fig1B', 'Fig6', 'Fig6A_original'):
        for ext in ('png', 'pdf'):
            (OUT / f'{stale}.{ext}').unlink(missing_ok=True)
    (OUT / 'values.json').write_text(json.dumps(values, indent=1, ensure_ascii=False), encoding='utf-8')
    lines = ['| Placeholder | Value | Definition |', '| --- | --- | --- |']
    lines += [f'| ⟦{k}⟧ | {v} | {definitions[k]} |' for k, v in values.items()]
    if missing:
        lines += ['', 'Missing inputs:', *[f'- {m}' for m in missing]]
    (OUT / 'values.md').write_text('\n'.join(lines) + '\n', encoding='utf-8')
    print(f'{len(values)} values and figures written to {OUT}')
    if missing:
        print('Missing:', *missing, sep='\n  ')


if __name__ == '__main__':
    main()
