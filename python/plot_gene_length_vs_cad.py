"""Run the full finite-rate counterpart of MATLAB's gene-length/CAD pipeline."""
import argparse
from dataclasses import asdict, replace
from datetime import datetime
import hashlib
import json
import os
from pathlib import Path
import time

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy.sparse.linalg import spsolve

from cpadynamics.full_model import FullModel
from cpadynamics.gene_length import GeneLengthModel, length_population
from cpadynamics.observables import cleavage_profile
from cpadynamics.parameters import working_parameters


def validate_sparse_reference(parameters, m, pools, lengths, profiles, downstream):
    """Check representative genes against the complete sparse operator."""
    records = []
    for length in (2500, 20000, 200000):
        index = int(np.flatnonzero(lengths == length)[0])
        p = replace(parameters, PASposition=length, geneLength_bp=length+downstream)
        model = FullModel(p, m)
        ef, rf = pools['E_free'], pools['Pol_free']
        occupation = spsolve(-(model.a0+ef*model.ae), model.unit_source)
        flux = p.k_in*rf
        solution = model.summarize(np.r_[flux*occupation, ef, rf])
        reference = cleavage_profile(solution.R, solution.RHE, p)
        np.testing.assert_allclose(profiles.pol_per_flux[index], occupation.sum(), rtol=1e-9)
        np.testing.assert_allclose(profiles.e_per_flux[index], model.e_counts @ occupation, rtol=1e-9)
        np.testing.assert_allclose(profiles.cleavage_cdf[index], reference.cdf, atol=1e-10)
        np.testing.assert_allclose(profiles.cad50_bp[index], reference.cad50_bp, rtol=1e-9)
        if solution.rhs_max_abs > 1e-6:
            raise RuntimeError('Full sparse kinetic residual failed validation')
        records.append(dict(length_bp=length,
                            cad50_bp=float(reference.cad50_bp) if np.isfinite(reference.cad50_bp) else None,
                            rhs_max_abs=solution.rhs_max_abs,
                            relative_pol_occupancy_error=float(abs(profiles.pol_per_flux[index]/occupation.sum()-1)),
                            relative_e_occupancy_error=float(abs(profiles.e_per_flux[index]/(model.e_counts @ occupation)-1))))
    return records


def run_analysis(output_dir, m=5):
    start = time.perf_counter()
    root = Path(__file__).resolve().parents[1]
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    num_genes, downstream = 10000, 5000
    parameters = replace(working_parameters(), k_in=working_parameters().k_in/num_genes)
    lengths, weights, probability_mass = length_population(spacing=parameters.L_a)
    print(f'Building full-model length solver: M={m}, {len(lengths)} lattice lengths...', flush=True)
    model = GeneLengthModel(parameters, m, lengths, downstream)
    print('Solving shared genome-wide Pol II and E pools...', flush=True)
    pools, profiles = model.solve_shared_pools(weights, num_genes)
    print(json.dumps(pools, indent=2), flush=True)
    print('Checking 2.5, 20 and 200 kb genes with the full sparse solver...', flush=True)
    validation = validate_sparse_reference(parameters, m, pools, lengths, profiles, downstream)
    flux = parameters.k_in*pools['Pol_free']
    values = np.column_stack((lengths, profiles.cad50_bp, profiles.cleavage_cdf[:, -1],
                              flux*profiles.pol_per_flux, flux*profiles.e_per_flux, weights))
    header = 'TSS_to_PAS_bp,CAD50_bp,within_window_cleavage_fraction,Pol_occupied,E_occupied,population_weight'
    np.savetxt(output_dir / 'gene_length_vs_cad.csv', values, delimiter=',', header=header, comments='')
    selected_lengths = np.array([2500, 5000, 10000, 20000, 50000, 100000, 200000])
    selected = np.searchsorted(lengths, selected_lengths)
    np.savetxt(output_dir / 'selected_lengths.csv', values[selected], delimiter=',', header=header, comments='')
    np.savez_compressed(output_dir / 'cleavage_profiles.npz', lengths_bp=lengths,
                        distances_bp=profiles.distances_bp, cleavage_cdf=profiles.cleavage_cdf)

    fig, ax = plt.subplots(figsize=(8.5, 5.5), constrained_layout=True)
    ax.plot(lengths/1000, profiles.cad50_bp, color='#2364aa', linewidth=2.3,
            label=f'Full finite-rate model (M={m})')
    ax.scatter(lengths[selected]/1000, profiles.cad50_bp[selected], color='#2364aa', s=22)
    ax.set(xscale='log', xlabel='TSS-to-PAS distance (kb)', ylabel=r'CAD$_{50}$ (bp)',
           title='Gene length versus cleavage distance', xlim=(2.5, 200))
    ticks = [2.5, 5, 10, 20, 50, 100, 200]
    ax.set_xticks(ticks, labels=[str(t) for t in ticks])
    ax.legend(frameon=False, loc='upper right')
    ax.text(0.98, 0.84, '10,000 genes sharing Pol II and E pools\n5 kb downstream extension',
            transform=ax.transAxes, ha='right', va='top', fontsize=10, color='#555555')
    ax.grid(alpha=0.2)
    ax.spines[['top', 'right']].set_visible(False)
    missing = ~np.isfinite(profiles.cad50_bp)
    if missing.any():
        ax.text(0.02, 0.03, f'50% cleavage not reached for {missing.sum()} lengths; shown as gaps',
                transform=ax.transAxes, fontsize=9)
    for extension in ('png', 'svg'):
        fig.savefig(output_dir / f'gene_length_vs_cad.{extension}', dpi=200)
    plt.close(fig)

    sources = [Path(__file__).resolve(), root/'default_parameters.m',
               root/'GeneLengthGenerateGrid.m', root/'GeneLengthBuildInterpolation.m',
               root/'GeneLengthAnalyze.m', *sorted((root/'python/cpadynamics').glob('*.py'))]
    manifest = dict(
        model_variant='full_kinetics_rapid_EH_disassembly', pool_mode='shared_genome_pools',
        parameters=asdict(parameters), binding_number=m, num_active_genes=num_genes,
        downstream_extension_bp=downstream, pools=pools,
        population=dict(minimum_bp=2500, maximum_bp=200000, median_bp=22000,
                        log10_sigma=0.68, untruncated_mass_in_range=probability_mass,
                        integration='Exact lognormal CDF differences over every 100 bp lattice-length bin; normalized on range'),
        method='Forward block solution of full finite-rate equations, bracketed E root, analytical Pol conservation; no occupancy interpolant',
        cleavage_convention='Includes PAS bin; its CDF is at 100 bp. The 5 kb extension has 51 PAS-associated bins through 5100 bp, matching MATLAB.',
        lattice_length_count=len(lengths), cad_threshold_reached_count=int((~missing).sum()),
        max_unit_flux_residual=profiles.max_flux_residual,
        sparse_reference_validation=validation,
        versions=dict(numpy=np.__version__, scipy=scipy.__version__, matplotlib=matplotlib.__version__),
        elapsed_seconds=time.perf_counter()-start,
        source_sha256={str(path.relative_to(root)): hashlib.sha256(path.read_bytes()).hexdigest() for path in sources},
    )
    (output_dir/'manifest.json').write_text(json.dumps(manifest, indent=2, allow_nan=False)+'\n')
    report = [f'# Full finite-rate gene length versus CAD (M={m})', '',
              f'Free E: {pools["E_free"]:.6f}; free Pol II: {pools["Pol_free"]:.6f}.', '',
              '| TSS-to-PAS (kb) | CAD50 (bp) | Cleaved within window |', '| ---: | ---: | ---: |']
    for row in values[selected]:
        report.append(f'| {row[0]/1000:g} | {row[1]:.2f} | {100*row[2]:.2f}% |')
    report += ['', 'The distribution is normalized over 2.5–200 kb. Its untruncated probability mass in this interval is '
               f'{probability_mass:.6f}. Every lattice-length bin is integrated using its exact lognormal probability.',
               '', 'This solves the same full kinetic equations as MATLAB, but bypasses the MATLAB occupancy interpolation grid. '
               'It therefore need not reproduce coarse-grid interpolation error. MATLAB itself was not executed.',
               '', 'Cleavage includes the PAS bin: the 5 kb extension produces CDF coordinates through 5100 bp, matching MATLAB. '
               'Unreached CAD thresholds remain NaN.',
               '', f'Full sparse reference checks passed at 2.5, 20 and 200 kb. Runtime: {manifest["elapsed_seconds"]:.1f} seconds.']
    (output_dir/'REPORT.md').write_text('\n'.join(report)+'\n')
    print('\n'.join(report), flush=True)
    print(f'Results: {output_dir.resolve()}', flush=True)
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--m', type=int, default=5, help='Maximum bound E/phosphorylation capacity (default 5)')
    parser.add_argument('--output-dir', type=Path)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    results_root = Path(os.environ.get('CPAD_RESULTS_ROOT', root/'SecondVersionResults'))
    output = args.output_dir or results_root/'python_GeneLengthAnalysis'/datetime.now().strftime('%Y%m%d_%H%M%S_%f')
    run_analysis(output, args.m)


if __name__ == '__main__':
    main()
