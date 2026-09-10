"""Plot saved CAD sensitivity to the H-engaged E off-rate."""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('results', type=Path)
    args = parser.parse_args()
    data = np.atleast_1d(np.genfromtxt(args.results/'sweep.csv', names=True, delimiter=',', dtype=None, encoding='utf-8'))
    reference = np.atleast_1d(np.genfromtxt(args.results/'equilibrium_reference.csv', names=True, delimiter=','))
    reference = reference[np.isin(reference['M'], data['M'])]
    fig, ax = plt.subplots(figsize=(8.4, 5.5), constrained_layout=True)
    ax.plot(reference['M'], reference['equilibrium_cad50_bp'], 'o--', color='#6b7280',
            linewidth=1.8, label='Original equilibrium reference')
    colors = ['#176B87', '#CA6F1E', '#8E4585', '#467342']
    for index, rate in enumerate(np.unique(data['kEoff_engaged'])):
        subset = np.sort(data[data['kEoff_engaged'] == rate], order='M')
        label = f'Full: engaged E off = {rate:g} s⁻¹'
        if rate == 0:
            label += ' (previous model)'
        ax.plot(subset['M'], subset['cad50_bp'], 'o-', color=colors[index % len(colors)],
                linewidth=2, markersize=6, label=label)
    ax.set(xlabel='Maximum E binding and phosphorylation capacity (M)', ylabel='CAD₅₀ (bp)',
           title='Rapid disassembly after engaged E detachment', xticks=np.unique(data['M']),
           ylim=(0, None))
    ax.grid(alpha=.2)
    ax.spines[['top', 'right']].set_visible(False)
    ax.legend(frameon=False, fontsize=9)
    for extension in ('png', 'svg'):
        fig.savefig(args.results/f'engaged_e_off_sweep.{extension}', dpi=180)
    plt.close(fig)
    print(args.results/'engaged_e_off_sweep.png')


if __name__ == '__main__':
    main()
