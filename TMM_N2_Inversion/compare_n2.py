#!/usr/bin/env python3
import argparse
import numpy as np
import matplotlib.pyplot as plt
import os
import sys

def load_csv(path):
    """Safely load CSV, ignoring non-numeric header rows."""
    if not os.path.exists(path):
        print(f"Error: File not found: {path}", file=sys.stderr)
        sys.exit(1)
        
    data = []
    with open(path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                row = [float(x) for x in line.split(',')]
                data.append(row)
            except ValueError:
                pass  # skip header
    return np.array(data)

def setup_plot_style():
    """Apply clean, modern scientific styling."""
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'axes.labelsize': 12,
        'axes.titlesize': 14,
        'legend.fontsize': 10,
        'xtick.labelsize': 10,
        'ytick.labelsize': 10,
        'figure.dpi': 120,
        'figure.figsize': (10, 8),
        'axes.grid': True,
        'grid.alpha': 0.4,
        'grid.linestyle': '--'
    })

def main():
    parser = argparse.ArgumentParser(description="Visually compare two refractive index CSVs.")
    parser.add_argument("file1", help="First CSV, e.g. output_n2.csv")
    parser.add_argument("file2", help="Second CSV, e.g. result_refractive_index.csv")
    args = parser.parse_args()

    data1 = load_csv(args.file1)
    data2 = load_csv(args.file2)

    f1, nr1, ni1 = data1[:, 0], data1[:, 1], data1[:, 2]
    f2, nr2, ni2 = data2[:, 0], data2[:, 1], data2[:, 2]

    setup_plot_style()
    
    # Create a 2x2 grid. 
    # Top row: Re(n) overlay, Im(n) overlay
    # Bottom row: Re(n) difference, Im(n) difference
    fig, axs = plt.subplots(2, 2, sharex='col')
    fig.suptitle(f'Comparison: {os.path.basename(args.file1)} vs {os.path.basename(args.file2)}', 
                 fontweight='bold')

    name1 = os.path.basename(args.file1)
    name2 = os.path.basename(args.file2)

    # Colors
    c1, c2 = '#1f77b4', '#ff7f0e'  # clean distinct colors

    # 1. Real Part (Overlay)
    axs[0, 0].plot(f1, nr1, label=name1, color=c1, lw=2, alpha=0.8)
    axs[0, 0].plot(f2, nr2, label=name2, color=c2, linestyle='--', lw=2, alpha=0.8)
    axs[0, 0].set_ylabel(r'Real part, $\mathrm{Re}(n)$')
    axs[0, 0].set_title('Real Part Overlaid')
    axs[0, 0].legend()

    # 2. Imaginary Part (Overlay)
    axs[0, 1].plot(f1, ni1, label=name1, color=c1, lw=2, alpha=0.8)
    axs[0, 1].plot(f2, ni2, label=name2, color=c2, linestyle='--', lw=2, alpha=0.8)
    axs[0, 1].set_ylabel(r'Imaginary part, $\mathrm{Im}(n)$')
    axs[0, 1].set_title('Imaginary Part Overlaid')
    axs[0, 1].legend()

    # Since frequency grids might not match exactly due to floats, we interpolate
    # for a robust difference computation
    nr2_interp = np.interp(f1, f2, nr2)
    ni2_interp = np.interp(f1, f2, ni2)
    diff_r = nr1 - nr2_interp
    diff_i = ni1 - ni2_interp

    # 3. Real Part (Difference)
    axs[1, 0].plot(f1, diff_r, color='#d62728', lw=1.5)
    axs[1, 0].set_xlabel('Frequency (THz)')
    axs[1, 0].set_ylabel(r'$\Delta \mathrm{Re}(n)$')
    axs[1, 0].set_title('Difference (Real)')
    axs[1, 0].axhline(0, color='black', lw=1, ls=':')

    # 4. Imaginary Part (Difference)
    axs[1, 1].plot(f1, diff_i, color='#d62728', lw=1.5)
    axs[1, 1].set_xlabel('Frequency (THz)')
    axs[1, 1].set_ylabel(r'$\Delta \mathrm{Im}(n)$')
    axs[1, 1].set_title('Difference (Imaginary)')
    axs[1, 1].axhline(0, color='black', lw=1, ls=':')

    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    main()
