"""Sensitivity to H-engaged E detachment followed by rapid EH disassembly."""
import argparse
from dataclasses import asdict, replace
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform

import numpy as np
import scipy

from cpadynamics import Parameters, FullModel, cleavage_profile, solve_equilibrium
from run_comparison import TARGETS, write_csv


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--m', nargs='+', type=int, choices=range(1, 8), default=list(range(1, 7)))
    parser.add_argument('--relative-off-rates', nargs='+', type=float, default=[0, .1, 1],
                        help='kEoff_engaged / kEoff; illustrative sensitivity values, not fitted rates')
    parser.add_argument('--output', type=Path, help='New output directory; must not exist')
    args = parser.parse_args()
    if any(not np.isfinite(v) or v < 0 for v in args.relative_off_rates):
        parser.error('Relative off-rates must be finite and nonnegative')
    ratios = sorted(set([0.] + args.relative_off_rates))
    capacities = sorted(set(args.m))
    p0 = Parameters()
    repo = Path(__file__).resolve().parents[1]
    stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S_%fZ')
    output = args.output or repo/'SecondVersionResults'/'python_engaged_E_dissociation'/stamp
    output.mkdir(parents=True, exist_ok=False)
    references, reference_rows = {}, []
    # All six supplied MATLAB targets are checked before running the new model.
    for m in sorted(set(TARGETS) | set(capacities)):
        s = solve_equilibrium(p0, m)
        profile = cleavage_profile(s.R, s.RHE, p0)
        references[m] = profile
        target = TARGETS.get(m)
        if target is not None and abs(profile.cad50_bp-target) >= .5:
            raise RuntimeError(f'MATLAB parity failed for M={m}')
        reference_rows.append(dict(M=m, equilibrium_cad50_bp=profile.cad50_bp,
                                   matlab_reported_bp=target if target is not None else '',
                                   E_conservation_residual=s.e_residual))
    write_csv(output/'equilibrium_reference.csv', reference_rows)

    rows = []
    for m in capacities:
        baseline_cad = None
        for index, ratio in enumerate(ratios):
            p = replace(p0, kEoff_engaged=ratio*p0.kEoff)
            model = FullModel(p, m)
            s = model.solve_steady_state()
            profile = cleavage_profile(s.R, s.RHE, p)
            if ratio == 0:
                baseline_cad = profile.cad50_bp
            label = f'M{m}_rate{index}'
            rows.append(dict(M=m, relative_off_rate=ratio, kEoff_engaged=p.kEoff_engaged,
                             cad50_bp=profile.cad50_bp, equilibrium_cad50_bp=references[m].cad50_bp,
                             percent_change_vs_protected=100*(profile.cad50_bp/baseline_cad-1),
                             max_exit_cdf=profile.max_exit_cdf,
                             E_free=s.E_free, Pol_free=s.Pol_free,
                             E_residual=s.e_residual, Pol_residual=s.pol_residual,
                             rhs_max_abs=s.rhs_max_abs,
                             flux_residual=profile.total_outflux-p.k_in*s.Pol_free,
                             engaged_E_loss_flux=p.kEoff_engaged*s.RHE.sum(),
                             data_prefix=label))
            write_csv(output/f'cleavage_{label}.csv', [dict(distance_bp=d, cdf=c)
                      for d, c in zip(profile.distances_bp, profile.cdf)])
            np.savez_compressed(output/f'microstates_{label}.npz', R=s.R_micro, RHE=s.RHE_micro,
                                R_states=np.asarray(model.r_states), RHE_states=np.asarray(model.rhe_states),
                                E_free=s.E_free, Pol_free=s.Pol_free)
            print(f'M={m}, kEoff_engaged={p.kEoff_engaged:g}/s: CAD50={profile.cad50_bp:.3f} bp', flush=True)
    write_csv(output/'sweep.csv', rows)
    sources = sorted((repo/'python'/'cpadynamics').glob('*.py')) + [
        Path(__file__).resolve(), repo/'python'/'run_comparison.py',
        repo/'default_parameters.m', repo/'EBindingNumberVsCad.m']
    manifest = dict(created_utc=stamp, base_parameters=asdict(p0), M=capacities,
                    relative_off_rates=ratios,
                    rate_interpretation='Illustrative sensitivity values, not experimentally calibrated',
                    mechanism='RHE(p,e) -> R(p,e-1) + free E at kEoff_engaged; rapid EH disassembly',
                    python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__,
                    matlab_validation='Six user-supplied rounded CAD values; MATLAB not executed',
                    source_sha256={str(path.relative_to(repo)): hashlib.sha256(path.read_bytes()).hexdigest()
                                   for path in sources})
    (output/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    header = '| M | Equilibrium (bp) | ' + ' | '.join(f'Full: {r*p0.kEoff:g} /s (bp)' for r in ratios) + ' |'
    report = ['# H-engaged E dissociation sensitivity', '',
              'The engaged E can leave the CTD independently of H dissociation. EH then disassembles rapidly.',
              'The net event releases exactly one E, preserves p and polymerase, and returns RHE(p,e) to R(p,e-1).',
              'H-first dissociation still returns RHE(p,e) to R(p,e) at kHoff=2 /s.',
              'Other E factors unbind at (e-1)*kEoff. There is one engaged E, so its off-rate has no e multiplier.',
              '', 'The tested off-rates are illustrative, not measured or fitted.', '', header,
              '|'+'---:|'*(len(ratios)+2)]
    for m in capacities:
        values = [r['cad50_bp'] for r in rows if r['M'] == m]
        report.append(f'| {m} | {references[m].cad50_bp:.3f} | '+' | '.join(f'{v:.3f}' for v in values)+' |')
    report += ['', 'Zero rate reproduces the previously implemented protected-E model.',
               'The original MATLAB-equilibrium model is a reference, not a model containing this new pathway.',
               'No CAD is extrapolated when the threshold is not reached; such results are NaN.',
               'All returned steady states pass full kinetic and pool-conservation residual checks.',
               'Parameters, source hashes and detailed populations are saved with this report.']
    (output/'REPORT.md').write_text('\n'.join(report)+'\n')
    print(f'Results: {output}', flush=True)


if __name__ == '__main__':
    main()
