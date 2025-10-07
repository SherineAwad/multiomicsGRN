#!/usr/bin/env python3
import argparse
import pyBigWig
import matplotlib.pyplot as plt
import numpy as np

def read_bw_signal(bw_path, chrom, start, end):
    """Read signal values from a bigWig file for the given region."""
    with pyBigWig.open(bw_path) as bw:
        chroms = bw.chroms()
        if chrom not in chroms:
            raise ValueError(f"{chrom} not found in {bw_path}. Available: {list(chroms.keys())[:5]}...")
        end = min(end, chroms[chrom])
        values = np.array(bw.values(chrom, start, end, numpy=True))
    return values, end

def downsample_signal(positions, values, bin_size=1000):
    """Average signal over bins for large regions."""
    if len(positions) <= bin_size:
        return positions, values
    binned_positions = positions[::bin_size]
    binned_values = [np.nanmean(values[i:i+bin_size]) for i in range(0, len(values), bin_size)]
    return binned_positions, np.array(binned_values)

def normalize_signal(values):
    """Normalize signal to 0-100 for plotting."""
    max_val = np.nanmax(values)
    if max_val == 0:
        return values
    return values / max_val * 100

def plot_bw_tracks(bw1, bw2, chrom, start, end, output):
    values1, end1 = read_bw_signal(bw1, chrom, start, end)
    values2, end2 = read_bw_signal(bw2, chrom, start, end)
    end = min(end1, end2)
    positions = np.arange(start, end)

    # Downsample if region is large
    max_points = 20000  # max points to plot for performance
    if len(positions) > max_points:
        bin_size = len(positions) // max_points
        positions, values1 = downsample_signal(positions, values1, bin_size)
        _, values2 = downsample_signal(positions, values2, bin_size)

    # Normalize signals
    values1 = normalize_signal(values1)
    values2 = normalize_signal(values2)

    # Plot
    fig, axes = plt.subplots(2, 1, figsize=(20, 6), sharex=True)

    # Track 1
    axes[0].fill_between(positions, 0, values1, color='steelblue', alpha=1.0, linewidth=0)
    axes[0].set_ylabel(bw1.split('/')[-1].replace('.bw', ''), rotation=0, labelpad=40, fontsize=10)
    axes[0].set_ylim(0, 110)
    axes[0].set_yticks([])
    axes[0].spines['top'].set_visible(False)
    axes[0].spines['right'].set_visible(False)
    axes[0].spines['left'].set_visible(False)

    # Track 2
    axes[1].fill_between(positions, 0, values2, color='darkorange', alpha=1.0, linewidth=0)
    axes[1].set_ylabel(bw2.split('/')[-1].replace('.bw', ''), rotation=0, labelpad=40, fontsize=10)
    axes[1].set_ylim(0, 110)
    axes[1].set_yticks([])
    axes[1].spines['top'].set_visible(False)
    axes[1].spines['right'].set_visible(False)
    axes[1].spines['left'].set_visible(False)

    axes[1].set_xlabel(f"{chrom}:{start}-{end}")
    plt.tight_layout()
    plt.savefig(output, dpi=300)
    print(f"Saved genome-style filled tracks to {output}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot two bigWig tracks as filled genome-style coverage tracks.")
    parser.add_argument("--bw1", required=True, help="Path to first bigWig file")
    parser.add_argument("--bw2", required=True, help="Path to second bigWig file")
    parser.add_argument("--chrom", required=True, help="Chromosome name (e.g. chr1)")
    parser.add_argument("--start", type=int, required=True, help="Start coordinate")
    parser.add_argument("--end", type=int, required=True, help="End coordinate")
    parser.add_argument("--output", required=True, help="Output image file")
    args = parser.parse_args()

    plot_bw_tracks(args.bw1, args.bw2, args.chrom, args.start, args.end, args.output)

