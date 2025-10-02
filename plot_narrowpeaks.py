import argparse
import pandas as pd
import matplotlib.pyplot as plt

def load_narrowpeak(path):
    peaks = pd.read_csv(path, sep="\t", header=None)
    peaks.columns = [
        "chrom", "start", "end", "name", "score",
        "strand", "signal", "pvalue", "qvalue", "peak"
    ]
    return peaks

def plot_peak_scatter(peaks, out):
    # Sort chromosomes
    def chrom_key(c):
        c = c.replace("chr", "")
        return (0, int(c)) if c.isdigit() else (1, c)
    chroms = sorted(peaks["chrom"].unique(), key=chrom_key)

    peaks["length"] = peaks["end"] - peaks["start"]
    chrom_to_y = {chrom: i for i, chrom in enumerate(chroms)}
    peaks["y"] = peaks["chrom"].map(chrom_to_y)
    peaks["mid"] = (peaks["start"] + peaks["end"]) / 2

    plt.figure(figsize=(14, 0.4 * len(chroms)))

    plt.scatter(
        peaks["mid"],
        peaks["y"],
        s=peaks["length"] / 100,   # scale size by length
        c=peaks["signal"],
        cmap="viridis",
        alpha=0.5,
        edgecolor="none"
    )

    plt.yticks(range(len(chroms)), chroms)
    plt.xlabel("Genomic position (bp)")
    plt.ylabel("Chromosome")
    plt.colorbar(label="Signal intensity")
    plt.title("All Peaks Across Chromosomes (dot size = length)")
    plt.tight_layout()
    plt.savefig(out, dpi=300)
    plt.close()

def main():
    parser = argparse.ArgumentParser(description="Scatter plot of all narrowPeak peaks")
    parser.add_argument("narrowpeak", help="Path to .narrowPeak file")
    parser.add_argument("-o", "--output", default="peak_scatter.png", help="Output image")
    args = parser.parse_args()

    peaks = load_narrowpeak(args.narrowpeak)
    plot_peak_scatter(peaks, args.output)
    print(f"Scatter plot saved to {args.output}")

if __name__ == "__main__":
    main()

