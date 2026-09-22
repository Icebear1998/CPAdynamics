"""Run the full finite-rate counterpart of CpaMultipleEMain.m.

Saves E/Ser2P and R/RHE steady-state plots plus raw profiles, cleavage CDF,
diagnostics and a manifest. The model is solved once at M=5 with the working
parameters in default_parameters.m.
"""
import argparse
from dataclasses import asdict, replace
from datetime import datetime
import hashlib
import json
import os
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import scipy

from cpadynamics.full_model import FullModel
from cpadynamics.observables import cleavage_profile
from cpadynamics.parameters import working_parameters
from plot_ebinding_profile import binding_profiles


def run_analysis(output_dir, capacity=5, parameters=None):
    """Run one closed-pool full-model steady state and save MATLAB-like output."""
    if capacity < 1:
        raise ValueError('capacity must be a positive integer')
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    parameter_source = ('default_parameters.m, working values checked 2026-09-15'
                        if parameters is None else 'Explicit parameters supplied to run_analysis')
    parameters = working_parameters() if parameters is None else parameters
    model = FullModel(parameters, capacity)
    solution = model.solve_steady_state()
    avg_e, ser2p = binding_profiles(model, solution)
    cleavage = cleavage_profile(solution.R, solution.RHE, parameters)
    positions_bp = (np.arange(model.n) - model.pas) * parameters.L_a
    rhe_full = np.zeros(model.n)
    rhe_full[model.pas:] = solution.RHE
    flux_residual = cleavage.total_outflux - parameters.k_in * solution.Pol_free

    if (not (np.isfinite(avg_e).all() and np.isfinite(ser2p).all()
             and np.isfinite(solution.E_free) and np.isfinite(solution.Pol_free))
            or np.any(avg_e < -1e-9) or np.any(avg_e > ser2p + 1e-9)
            or np.any(ser2p > capacity + 1e-9)
            or max(abs(solution.e_residual), abs(solution.pol_residual),
                   solution.rhs_max_abs, abs(flux_residual)) > 1e-6):
        raise RuntimeError('Steady state failed profile, conservation or flux validation')

    np.savetxt(output_dir / 'full_profiles.csv',
               np.column_stack((positions_bp, solution.R, rhe_full, avg_e, ser2p)),
               delimiter=',', comments='',
               header='PositionRelPAS_bp,R,RHE,Average_E_bound,Average_Ser2P')
    np.savetxt(output_dir / 'cleavage_cdf.csv',
               np.column_stack((cleavage.distances_bp, cleavage.cdf)),
               delimiter=',', comments='', header='DistanceAfterPAS_bp,Cleavage_CDF')
    np.savez_compressed(output_dir / 'microstates.npz',
                        R_micro=solution.R_micro, RHE_micro=solution.RHE_micro,
                        R_states=np.asarray(model.r_states),
                        RHE_states=np.asarray(model.rhe_states), state=solution.state)

    fig, ax = plt.subplots(figsize=(8.8, 5.5), constrained_layout=True)
    ax.plot(positions_bp, ser2p, color='#2ca02c', linewidth=2.5, label='Ser2P')
    ax.plot(positions_bp, avg_e, color='#1f77b4', linewidth=2.5, label='Average E')
    ax.axvline(0, color='#222222', linestyle='--', linewidth=1, label='_nolegend_')
    ax.set(xlabel='Position relative to PAS (bp)', ylabel='Average count per polymerase',
           title='Full finite-rate E binding and phosphorylation')
    ax.legend(frameon=False, loc='best')
    ax.grid(alpha=0.2)
    ax.spines[['top', 'right']].set_visible(False)
    for suffix in ('png', 'svg'):
        fig.savefig(output_dir / f'E_binding_and_Ser2P.{suffix}', dpi=200)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.8, 5.5), constrained_layout=True)
    ax.plot(positions_bp, solution.R, color='#1f77b4', linewidth=2.5, label='R')
    ax.plot(positions_bp, rhe_full, color='#d62728', linewidth=2.5, label='RHE')
    if np.isfinite(cleavage.cad50_bp):
        ax.axvline(cleavage.cad50_bp, color='#d62728', linestyle='--', linewidth=1.5,
                   label=r'CAD$_{50}$')
    ax.axvline(0, color='#222222', linestyle=':', linewidth=1, label='_nolegend_')
    ax.set(xlabel='Position relative to PAS (bp)', ylabel='Polymerase count',
           title='Full finite-rate CPA steady state')
    ax.legend(frameon=False, loc='best')
    ax.grid(alpha=0.2)
    ax.spines[['top', 'right']].set_visible(False)
    for suffix in ('png', 'svg'):
        fig.savefig(output_dir / f'R_and_RHE_steady_state.{suffix}', dpi=200)
    plt.close(fig)

    root = Path(__file__).resolve().parents[1]
    sources = [Path(__file__).resolve(), root / 'CpaMultipleEMain.m',
               root / 'default_parameters.m', root / 'run_full_termination_simulation.m',
               *sorted((root / 'python' / 'cpadynamics').glob('*.py'))]
    diagnostics = dict(
        CAD50_bp=cleavage.cad50_bp,
        within_window_cleavage_fraction=cleavage.max_exit_cdf,
        total_outflux=cleavage.total_outflux,
        E_free=solution.E_free, E_bound=parameters.E_total-solution.E_free,
        Pol_free=solution.Pol_free, Pol_bound=parameters.Pol_total-solution.Pol_free,
        E_conservation_residual=solution.e_residual,
        Pol_conservation_residual=solution.pol_residual,
        rhs_max_abs=solution.rhs_max_abs, flux_residual=flux_residual,
    )
    manifest = dict(
        model_variant='full_kinetics_rapid_EH_disassembly',
        parameter_source=parameter_source,
        parameters=asdict(parameters), EBindingNumber=capacity,
        position_coordinates='Full gene, relative to PAS; PAS node is 0 bp.',
        cleavage_convention='Cleavage in first PAS-associated node represents (0, L_a] bp.',
        diagnostics=diagnostics,
        versions=dict(numpy=np.__version__, scipy=scipy.__version__,
                      matplotlib=matplotlib.__version__),
        source_sha256={str(path.relative_to(root)): hashlib.sha256(path.read_bytes()).hexdigest()
                       for path in sources},
    )
    (output_dir / 'manifest.json').write_text(json.dumps(manifest, indent=2, allow_nan=False) + '\n')
    print(json.dumps(diagnostics, indent=2))
    print(f'Results: {output_dir.resolve()}')
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, help='Directory for figures, data and diagnostics')
    parser.add_argument('--m', type=int, default=5, help='Maximum E/phosphorylation capacity (default: 5)')
    parser.add_argument('--parameters', type=Path,
                        help='JSON object of parameter overrides applied to the working parameter set')
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    results_root = Path(os.environ.get('CPAD_RESULTS_ROOT', root / 'SecondVersionResults'))
    output_dir = args.output_dir or (results_root / 'python_CPA_multipleE_main'
                                     / datetime.now().strftime('%Y%m%d_%H%M%S_%f'))
    parameters = None
    if args.parameters:
        parameters = replace(working_parameters(), **json.loads(args.parameters.read_text()))
    run_analysis(output_dir, args.m, parameters)


if __name__ == '__main__':
    main()
