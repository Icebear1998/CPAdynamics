"""Render saved comparison CSVs; does not rerun the model."""
import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('results', type=Path)
    args = parser.parse_args()
    with (args.results/'comparison.csv').open() as handle:
        summary = list(csv.DictReader(handle))
    fig, axes = plt.subplots(1, len(summary), figsize=(5.2*len(summary), 4.2), squeeze=False,
                             sharey=True, constrained_layout=True)
    for ax, row in zip(axes[0], summary):
        m = int(row['M'])
        data = np.genfromtxt(args.results/f'cleavage_M{m}.csv', delimiter=',', names=True)
        x = data['distance_bp']
        ax.plot(x, 100*data['equilibrium_cdf'], '--', color='#6b7280', linewidth=2,
                label=f'Equilibrium: {float(row["equilibrium_cad50_bp"]):.0f} bp')
        ax.plot(x, 100*data['full_cdf'], color='#176B87', linewidth=2.4,
                label=f'Full R/RHE: {float(row["full_cad50_bp"]):.0f} bp')
        ax.axhline(50, color='#b3b3b3', linewidth=1, linestyle=':')
        ax.set(xlabel='Distance downstream of PAS (bp)', title=f'M = {m}',
               ylim=(0, 100), xlim=(0, min(x[-1], 3000)))
        ax.grid(alpha=.18)
        ax.spines[['top', 'right']].set_visible(False)
        ax.legend(title='CAD₅₀', frameon=False, loc='lower right')
    axes[0, 0].set_ylabel('Cumulative cleavage (% of incoming flux)')
    fig.suptitle('Equilibrium closure versus full kinetics', fontsize=14)
    fig.savefig(args.results/'cleavage_comparison.png', dpi=180)
    fig.savefig(args.results/'cleavage_comparison.svg')
    plt.close(fig)
    print(args.results/'cleavage_comparison.png')


if __name__ == '__main__':
    main()
