"""Run the full-model equivalent of PlotEBindingProfile.m and save its figure/data."""
import argparse
from dataclasses import asdict
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
from cpadynamics.parameters import working_parameters as profile_parameters


def binding_profiles(model, solution):
    """Mean E and phosphorylated sites over ALL polymerases at each node.

    RHE contributes from the PAS onward. Empty nodes return NaN, matching
    run_full_termination_simulation.m. State-map columns are (p, e).
    """
    moments = solution.R_micro @ np.asarray(model.r_states)
    moments[model.pas:] += solution.RHE_micro @ np.asarray(model.rhe_states)
    polymerases = solution.R.copy()
    polymerases[model.pas:] += solution.RHE
    averages = np.full_like(moments, np.nan, dtype=float)
    np.divide(moments, polymerases[:, None], out=averages,
              where=polymerases[:, None] > 0)
    return averages[:, 1], averages[:, 0]


def run_profiles(output_dir, capacities=(1, 5, 10)):
    parameters = profile_parameters()
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(9, 6), constrained_layout=True)
    diagnostics = []
    colors = ['#d62728', '#1f77b4', '#2ca02c']
    for index, capacity in enumerate(capacities):
        print(f'Simulating EBindingNumber = {capacity}...', flush=True)
        model = FullModel(parameters, capacity)
        solution = model.solve_steady_state()
        avg_e, avg_p = binding_profiles(model, solution)
        x = (np.arange(model.n) - model.pas) * parameters.L_a / 1000
        flux_residual = (parameters.k_e * solution.R[-1]
                         + parameters.k_e2 * solution.RHE[-1]
                         + parameters.kc * solution.RHE.sum()
                         - parameters.k_in * solution.Pol_free)
        if (not np.isfinite([avg_e, avg_p]).all()
                or np.any(avg_e < -1e-9) or np.any(avg_e > avg_p + 1e-9)
                or np.any(avg_p > capacity + 1e-9)
                or abs(flux_residual) > 1e-6):
            raise RuntimeError(f'Profile or flux validation failed for M={capacity}')
        np.savetxt(output_dir / f'ProfileData_N{capacity}.csv',
                   np.column_stack((x, avg_e, avg_p)), delimiter=',',
                   header='PositionRelPAS_kb,Avg_E_bound,Avg_Ser2P', comments='')
        np.savez_compressed(output_dir / f'microstates_N{capacity}.npz',
                            R_micro=solution.R_micro, RHE_micro=solution.RHE_micro,
                            R_states=model.r_states, RHE_states=model.rhe_states,
                            E_free=solution.E_free, Pol_free=solution.Pol_free)
        diagnostics.append(dict(
            M=int(capacity), E_free=solution.E_free, Pol_free=solution.Pol_free,
            E_conservation_residual=solution.e_residual,
            Pol_conservation_residual=solution.pol_residual,
            rhs_max_abs=solution.rhs_max_abs, flux_residual=float(flux_residual),
            avg_E_at_PAS=float(avg_e[model.pas]),
            avg_Ser2P_at_PAS=float(avg_p[model.pas]),
        ))
        color = colors[index % len(colors)]
        # Plot through the PAS node; retain the full simulation and raw exports.
        plot_nodes = slice(0, model.pas + 1)
        plot_x = np.arange(1, model.pas + 2) * parameters.L_a / 1000
        ax.plot(plot_x, avg_e[plot_nodes], color=color, linewidth=2,
                label=f'Avg E (N={capacity})')
        ax.plot(plot_x, avg_p[plot_nodes], '--', color=color, linewidth=2,
                label=f'Ser2P (N={capacity})')

    pas_kb = (parameters.geometry()[1] + 1) * parameters.L_a / 1000
    ax.axvline(pas_kb, color='black', linestyle='--', linewidth=1.2)
    ax.annotate('PAS', xy=(pas_kb, 0.98), xycoords=ax.get_xaxis_transform(),
                xytext=(-6, 0), textcoords='offset points', ha='right', va='top')
    ax.set(xlabel='Distance from TSS (kb)', xlim=(0, pas_kb),
           xticks=np.linspace(0, pas_kb, 5),
           ylabel='Average bound E / phosphorylated sites per Pol II',
           title='Full finite-rate E binding (solid) and Ser2P (dashed)')
    ax.legend(loc='upper left', frameon=False)
    ax.grid(alpha=0.2)
    ax.spines[['top', 'right']].set_visible(False)
    for extension in ('png', 'svg'):
        fig.savefig(output_dir / f'Average_E_and_Ser2P_Comparison.{extension}', dpi=200)
    plt.close(fig)

    root = Path(__file__).resolve().parents[1]
    sources = [Path(__file__).resolve(), root / 'default_parameters.m',
               root / 'PlotEBindingProfile.m', root / 'run_full_termination_simulation.m',
               *sorted((root / 'python/cpadynamics').glob('*.py'))]
    manifest = dict(
        model_variant='full_kinetics_rapid_EH_disassembly',
        parameters=asdict(parameters), binding_numbers=list(capacities),
        parameter_source='default_parameters.m, working values on 2026-09-14',
        average_convention='Population-weighted over R and RHE at each node',
        plot_window='TSS to PAS, inclusive; distance from TSS in kb',
        raw_profile_coordinates='Full gene; position relative to PAS in kb',
        versions=dict(numpy=np.__version__, scipy=scipy.__version__,
                      matplotlib=matplotlib.__version__),
        diagnostics=diagnostics,
        source_sha256={str(path.relative_to(root)): hashlib.sha256(path.read_bytes()).hexdigest()
                       for path in sources},
    )
    (output_dir / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(json.dumps(diagnostics, indent=2))
    print(f'Results: {output_dir.resolve()}')
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, help='Directory for figure, CSVs and diagnostics')
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    results_root = Path(os.environ.get('CPAD_RESULTS_ROOT', root / 'SecondVersionResults'))
    output = args.output_dir or (results_root / 'python_EBindingProfile'
                                 / datetime.now().strftime('%Y%m%d_%H%M%S_%f'))
    run_profiles(output)


if __name__ == '__main__':
    main()
