"""Estimate joint 3D benchmark coverage for the Results framework.

Run from the repository root:
    python3 python/run_framework_monte_carlo.py --samples 10000 --workers 4

Independent log-uniform kHd, kEd and kc draws use the Figure 2C ranges.
Both capacities use the same draws, fixed association rates, and the full
finite-rate closed-pool solver. Unreached CAD50 values remain NaN and fail
the benchmark; solver failures abort rather than becoming failed benchmarks.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor
import csv
from dataclasses import asdict, replace
import hashlib
import json
import os
from pathlib import Path
import time

# Avoid multiplying BLAS threads across independent simulation workers.
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('OMP_NUM_THREADS', '1')
os.environ.setdefault('VECLIB_MAXIMUM_THREADS', '1')

import numpy as np
import scipy
from scipy.stats import binomtest

from generate_results_framework_figures import ROOT, read_parameters, solve_case


AXES = ('kHd', 'kEd', 'kc')


def evaluate_sample(job):
    index, values, parameters, threshold = job
    khd, ked, kc = values
    p = replace(parameters, kHoff=khd*parameters.kHon,
                kEoff=ked*parameters.kEon, kc=kc)
    row = dict(sample=index, kHd=khd, kEd=ked, kc=kc,
               kHoff=p.kHoff, kEoff=p.kEoff)
    for m in (1, 5):
        diagnostics = []
        cdf, solution, _, _ = solve_case(p, m, diagnostics)
        row.update({
            f'CAD50_M{m}_bp': cdf.cad50_bp,
            f'max_exit_cdf_M{m}': cdf.max_exit_cdf,
            f'meets_benchmark_M{m}': int(cdf.cad50_bp <= threshold),
            f'E_free_M{m}': solution.E_free,
            f'Pol_free_M{m}': solution.Pol_free,
            f'E_residual_M{m}': solution.e_residual,
            f'Pol_residual_M{m}': solution.pol_residual,
            f'rhs_max_abs_M{m}': solution.rhs_max_abs,
            f'flux_residual_M{m}': diagnostics[0]['flux_residual'],
        })
    return row


def summarize(rows, m):
    n = len(rows)
    successes = sum(r[f'meets_benchmark_M{m}'] for r in rows)
    interval = binomtest(successes, n).proportion_ci(
        confidence_level=.95, method='wilson')
    return dict(
        samples=n, successes=successes, percent=100*successes/n,
        wilson_95_percent=[100*interval.low, 100*interval.high],
        unreached=sum(not np.isfinite(r[f'CAD50_M{m}_bp']) for r in rows),
        max_E_residual=max(abs(r[f'E_residual_M{m}']) for r in rows),
        max_Pol_residual=max(abs(r[f'Pol_residual_M{m}']) for r in rows),
        max_rhs_residual=max(r[f'rhs_max_abs_M{m}'] for r in rows),
        max_flux_residual=max(r[f'flux_residual_M{m}'] for r in rows),
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--samples', type=int, default=10000)
    parser.add_argument('--seed', type=int, default=20261020)
    parser.add_argument('--workers', type=int, default=4)
    parser.add_argument('--threshold-bp', type=float, default=500)
    root = Path(os.environ.get('CPAD_RESULTS_ROOT', ROOT/'Results'))
    parser.add_argument('--output-dir', type=Path,
                        default=root/'ResultsFramework/draft/MonteCarlo3D')
    args = parser.parse_args()
    if (min(args.samples, args.workers) < 1 or args.seed < 0
            or not np.isfinite(args.threshold_bp) or args.threshold_bp <= 0):
        parser.error('Require positive samples, workers and threshold; nonnegative seed')
    p, ranges = read_parameters(ROOT/'default_parameters.m')
    bounds = np.array([ranges[name] for name in AXES])
    rng = np.random.Generator(np.random.PCG64(args.seed))
    uniform = rng.random((args.samples, 3))
    draws = np.exp(np.log(bounds[:, 0])
                   + uniform*(np.log(bounds[:, 1])-np.log(bounds[:, 0])))
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    sources = [Path(__file__), ROOT/'default_parameters.m',
               ROOT/'python/generate_results_framework_figures.py',
               *sorted((ROOT/'python/cpadynamics').glob('*.py'))]
    metadata = dict(
        status='running', samples=args.samples, seed=args.seed,
        generator='numpy.random.PCG64', axes=list(AXES),
        ranges={name: ranges[name] for name in AXES},
        distribution='Independent uniform draws in log(parameter); identical draws for M=1 and M=5.',
        threshold_bp=args.threshold_bp, parameters=asdict(p),
        association_rates='kHon and kEon fixed at defaults; kHoff=kHd*kHon and kEoff=kEd*kEon.',
        pool_mode='Closed local E and Pol II pools, solved for every draw and capacity.',
        unreached_policy='NaN CAD50 is a failed benchmark; max_exit_cdf is saved.',
        uncertainty='Marginal 95% Wilson intervals quantify Monte Carlo sampling error only.',
        source_sha256={str(f.relative_to(ROOT)): hashlib.sha256(f.read_bytes()).hexdigest()
                       for f in sources},
        versions=dict(numpy=np.__version__, scipy=scipy.__version__),
    )
    metadata_file = output/'metadata.json'
    metadata_file.write_text(json.dumps(metadata, indent=2, allow_nan=False)+'\n')
    jobs = ((i+1, values, p, args.threshold_bp) for i, values in enumerate(draws))
    rows, start = [], time.perf_counter()
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        with (output/'samples.csv').open('w', newline='') as stream:
            writer = None
            for row in pool.map(evaluate_sample, jobs, chunksize=10):
                if writer is None:
                    writer = csv.DictWriter(stream, fieldnames=list(row))
                    writer.writeheader()
                writer.writerow(row)
                rows.append(row)
                if len(rows) % 500 == 0 or len(rows) == args.samples:
                    stream.flush()
                    print(f'{len(rows)}/{args.samples} paired samples; '
                          f'{time.perf_counter()-start:.1f} s', flush=True)
    metadata.update(
        status='complete', elapsed_seconds=time.perf_counter()-start,
        results={f'M{m}': summarize(rows, m) for m in (1, 5)},
        prefix_estimates={
            str(n): {f'M{m}': summarize(rows[:n], m) for m in (1, 5)}
            for n in (1000, 5000) if n < args.samples},
        paired_counts={
            f'M1_{a}_M5_{b}': sum(r['meets_benchmark_M1'] == a
                                  and r['meets_benchmark_M5'] == b for r in rows)
            for a in (0, 1) for b in (0, 1)},
        samples_sha256=hashlib.sha256((output/'samples.csv').read_bytes()).hexdigest(),
    )
    metadata_file.write_text(json.dumps(metadata, indent=2, allow_nan=False)+'\n')
    print(json.dumps(metadata['results'], indent=2), flush=True)
    print(f'Outputs: {output}', flush=True)


if __name__ == '__main__':
    main()
