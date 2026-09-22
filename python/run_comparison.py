"""Validate the equilibrium port, then save full-model calculations (no plotting)."""
import argparse
import csv
from dataclasses import asdict, replace
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform

import numpy as np
import scipy

from cpadynamics import Parameters, FullModel, cleavage_profile, solve_equilibrium


TARGETS = {1: 1191, 2: 728, 3: 571, 4: 490, 5: 442, 6: 407}


def write_csv(path, rows):
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--m', nargs='+', type=int, choices=range(1, 8), default=[1, 5])
    parser.add_argument('--transient-seconds', type=float, default=0,
                        help='Also integrate from an empty gene for this duration; zero skips integration')
    parser.add_argument('--output', type=Path, help='New output directory; must not already exist')
    args = parser.parse_args()
    if not np.isfinite(args.transient_seconds) or args.transient_seconds < 0:
        parser.error('--transient-seconds must be finite and nonnegative')
    repo = Path(__file__).resolve().parents[1]
    timestamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S_%fZ')
    output = args.output or repo/'SecondVersionResults'/'python_full_model'/timestamp
    output.mkdir(parents=True, exist_ok=False)
    p = Parameters()
    parity_parameters = replace(p, kHon=4, kHoff=2, kc=0.13, kEoff_engaged=0)
    parity = []

    for m, target in TARGETS.items():
        reference = solve_equilibrium(parity_parameters, m)
        profile = cleavage_profile(reference.R, reference.RHE, parity_parameters)
        error = profile.cad50_bp-target
        passed = bool(abs(error) < .5 and abs(reference.e_residual) < 1e-6)
        parity.append(dict(M=m, matlab_reported_bp=target, python_equilibrium_bp=profile.cad50_bp,
                           difference_bp=error, E_free=reference.E_free,
                           E_conservation_residual=reference.e_residual, passed=passed))
        print(f'Parity M={m}: {profile.cad50_bp:.6f} bp vs MATLAB {target} bp: {"PASS" if passed else "FAIL"}', flush=True)
    write_csv(output/'matlab_parity.csv', parity)
    if not all(row['passed'] for row in parity):
        raise RuntimeError('MATLAB parity failed; full-model calculations were not run')

    comparison, dynamics = [], []
    for m in sorted(set(args.m)):
        model = FullModel(p, m)
        full = model.solve_steady_state()
        reference = solve_equilibrium(p, m)
        full_cdf = cleavage_profile(full.R, full.RHE, p)
        ref_cdf = cleavage_profile(reference.R, reference.RHE, p)
        comparison.append(dict(M=m, equilibrium_cad50_bp=ref_cdf.cad50_bp,
                               full_cad50_bp=full_cdf.cad50_bp,
                               percent_change=100*(full_cdf.cad50_bp/ref_cdf.cad50_bp-1),
                               full_E_free=full.E_free, full_Pol_free=full.Pol_free,
                               full_max_exit_cdf=full_cdf.max_exit_cdf,
                               full_E_residual=full.e_residual, full_Pol_residual=full.pol_residual,
                               full_rhs_max_abs=full.rhs_max_abs,
                               full_flux_residual=full_cdf.total_outflux-p.k_in*full.Pol_free,
                               polymerase_variables=model.size))
        write_csv(output/f'cleavage_M{m}.csv', [dict(distance_bp=d, equilibrium_cdf=e, full_cdf=f)
                  for d, e, f in zip(full_cdf.distances_bp, ref_cdf.cdf, full_cdf.cdf)])
        r_counts, h_counts = np.asarray(model.r_states), np.asarray(model.rhe_states)
        mean_r = np.divide(full.R_micro@r_counts, full.R[:, None],
                           out=np.zeros((model.n, 2)), where=full.R[:, None] > 0)
        mean_h = np.divide(full.RHE_micro@h_counts, full.RHE[:, None],
                           out=np.zeros((model.n_pas, 2)), where=full.RHE[:, None] > 0)
        write_csv(output/f'profiles_M{m}.csv', [dict(
            node_1based=i+1, kPon=p.kPon_min+p.kPon_slope*i, R=full.R[i],
            RHE=full.RHE[i-model.pas] if i >= model.pas else 0,
            mean_P_R=mean_r[i, 0], mean_E_R=mean_r[i, 1],
            mean_P_RHE=mean_h[i-model.pas, 0] if i >= model.pas else '',
            mean_E_RHE=mean_h[i-model.pas, 1] if i >= model.pas else '',
            equilibrium_R=reference.R[i],
            equilibrium_RHE=reference.RHE[i-model.pas] if i >= model.pas else 0,
            equilibrium_mean_E=reference.avg_E[i]) for i in range(model.n)])
        np.savez_compressed(output/f'microstates_M{m}.npz', R=full.R_micro, RHE=full.RHE_micro,
                            R_states=r_counts, RHE_states=h_counts, E_free=full.E_free,
                            Pol_free=full.Pol_free, state=full.state)
        print(f'Full M={m}: CAD50={full_cdf.cad50_bp:.6f} bp; RHS residual={full.rhs_max_abs:.3g}', flush=True)

        if args.transient_seconds:
            print(f'Integrating M={m} from empty gene to {args.transient_seconds:g} s...', flush=True)
            times = np.linspace(0, args.transient_seconds, 101)
            result = model.integrate(args.transient_seconds, t_eval=times)
            final = model.summarize(result.y[:, -1])
            # Species-wise scaled error; a large free pool must not hide errors in R/RHE.
            error = float(np.max(abs(final.state-full.state)/(1+abs(full.state))))
            dynamics.append(dict(M=m, time_seconds=args.transient_seconds,
                                 species_scaled_max_error=error, rhs_max_abs=final.rhs_max_abs,
                                 E_residual=final.e_residual, Pol_residual=final.pol_residual,
                                 converged=bool(error < 1e-5 and final.rhs_max_abs < 1e-4),
                                 function_evaluations=result.nfev))
            write_csv(output/f'transient_M{m}.csv', [dict(time_seconds=t,
                      E_free=y[-2], Pol_free=y[-1], R_total=y[:model.r_size].sum(),
                      RHE_total=y[model.r_size:model.size].sum(),
                      E_residual=y[-2]+model.e_counts@y[:model.size]-p.E_total,
                      Pol_residual=y[-1]+y[:model.size].sum()-p.Pol_total)
                      for t, y in zip(result.t, result.y.T)])
            np.savez_compressed(output/f'transient_final_M{m}.npz', state=final.state)
            print(f'  Transient error={error:.3g}, RHS={final.rhs_max_abs:.3g}; converged={dynamics[-1]["converged"]}', flush=True)

    write_csv(output/'comparison.csv', comparison)
    if dynamics:
        write_csv(output/'dynamic_validation.csv', dynamics)
    source_paths = [repo/name for name in ['default_parameters.m', 'build_full_rate_matrices.m',
                    'run_full_termination_simulation.m', 'ode_dynamics_full_multipleE.m',
                    'calculate_full_pas_cleavage_profile.m']]
    source_paths += sorted((repo/'python'/'cpadynamics').glob('*.py'))
    source_paths += [Path(__file__).resolve()]
    manifest = dict(created_utc=timestamp, parameters=asdict(p),
                    parity_parameters=asdict(parity_parameters), M=sorted(set(args.m)),
                    python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__,
                    matlab_execution='Not run: MATLAB/Octave unavailable. Targets supplied by user.',
                    parity_tolerance_bp=.5, transient_seconds=args.transient_seconds,
                    source_sha256={str(path.relative_to(repo)): hashlib.sha256(path.read_bytes()).hexdigest()
                                   for path in source_paths})
    (output/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    lines = ['# Python full-model comparison', '',
             'The equilibrium reference matches all six user-supplied MATLAB CAD50 values within rounding (0.5 bp).',
             'MATLAB was not executed in this environment.', '',
             '| M | Equilibrium CAD50 (bp) | Full CAD50 (bp) | Change |',
             '|---:|---:|---:|---:|']
    lines += [f'| {r["M"]} | {r["equilibrium_cad50_bp"]:.3f} | {r["full_cad50_bp"]:.3f} | {r["percent_change"]:+.1f}% |'
              for r in comparison]
    lines += ['', 'Full-model steady states retain every finite-rate reaction and transport equation.',
              'One E remains protected in RHE; recognition occurs at e*kHon; initiation is R(0,0).',
              'Cleavage starts in the PAS bin (0,100] bp. CDFs include (0,0) and terminal run-off.',
              'The difference includes finite kinetics, E protection in RHE, and occupancy-dependent recognition.',
              'It cannot be attributed solely to removal of timescale separation.', '',
              'Microstate arrays and state maps are stored in NPZ; profiles and CDFs are CSV.',
              'All parameters, versions and source hashes are recorded in manifest.json.']
    if dynamics:
        lines += ['', '## Time integration from an empty gene', '']
        lines += [f'M={r["M"]}: t={r["time_seconds"]:g} s; max species-scaled error '
                  f'{r["species_scaled_max_error"]:.3g}; max RHS {r["rhs_max_abs"]:.3g}; '
                  f'converged={r["converged"]}.' for r in dynamics]
    (output/'REPORT.md').write_text('\n'.join(lines)+'\n')
    print(f'Results: {output}', flush=True)


if __name__ == '__main__':
    main()
