"""Python EBindingNumberVsCad plot: compare the saved equilibrium and full sweeps."""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('results', type=Path, help='Directory produced by run_comparison.py')
    args = parser.parse_args()
    data = np.atleast_1d(np.genfromtxt(args.results/'comparison.csv', delimiter=',', names=True))
    data = np.sort(data, order='M')
    fig, ax = plt.subplots(figsize=(8, 5.5), constrained_layout=True)
    series = [('equilibrium_cad50_bp', 'Equilibrium (MATLAB reference)', '#6b7280', '--', -20, 'top'),
              ('full_cad50_bp', 'Full R/RHE kinetics', '#176B87', '-', 16, 'bottom')]
    for column, label, color, style, offset, alignment in series:
        ax.plot(data['M'], data[column], marker='o', linestyle=style, color=color,
                linewidth=2, markersize=7, label=label)
        for m, cad in zip(data['M'], data[column]):
            ax.annotate(f'{cad:.0f}', (m, cad), xytext=(0, offset), textcoords='offset points',
                        ha='center', va=alignment, fontsize=10, color=color)
    ax.set(xlabel='Maximum E binding and phosphorylation capacity (M)',
           ylabel='CAD₅₀ (bp downstream of PAS)', title='E-binding capacity versus cleavage distance',
           xticks=data['M'], xlim=(data['M'].min()-.4, data['M'].max()+.4),
           ylim=(0, 1.18*max(data['full_cad50_bp'].max(), data['equilibrium_cad50_bp'].max())))
    ax.grid(alpha=.2)
    ax.spines[['top', 'right']].set_visible(False)
    ax.legend(frameon=False, loc='upper right')
    for extension in ('png', 'svg'):
        fig.savefig(args.results/f'ebinding_number_vs_cad.{extension}', dpi=180)
    plt.close(fig)
    print(args.results/'ebinding_number_vs_cad.png')


if __name__ == '__main__':
    main()
