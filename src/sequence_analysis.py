"""
Stage 1 - Nucleotide quantification and visualization.

Reads a raw Martian sequence file, cleans it, counts how often each
nucleotide symbol occurs, and renders a frequency bar chart. Because the
Martian alphabet is unknown in advance, the bases are discovered
dynamically from the data rather than assumed.
"""

from pathlib import Path

import matplotlib.pyplot as plt

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = PROJECT_ROOT / "data"
RESULTS_DIR = PROJECT_ROOT / "results"


def read_file_sequence(file_path: str | Path) -> str | None:
    """Read a sequence file and return it as a single cleaned string.

    Newline characters are removed and the 'N' padding at the start and
    end of the sequence is stripped.

    Args:
        file_path: Path to the sequence file.

    Returns:
        The cleaned sequence, or ``None`` if the file does not exist.
    """
    sequence: list[str] = []
    try:
        with open(file_path, "r") as file:
            for line in file:
                line = line.strip()  # remove the newline characters
                if line:
                    sequence.append(line)
    except FileNotFoundError:
        print(f"Error: File not found at {file_path}")
        return None

    return "".join(sequence).strip("N")


def visualize_counts(title: str, counts_dict: dict[str, int], total_bases: int) -> None:
    """Draw a nucleotide-frequency bar chart and save it into ``results/``.

    Args:
        title: Label used for the chart title and output filename.
        counts_dict: Mapping of each nucleotide symbol to its count.
        total_bases: Total number of bases, used to compute percentages.
    """
    # sorted() keeps the bases displayed on the x-axis in ascending order
    bases_to_plot = sorted(b for b in counts_dict if b != "N")
    frequencies = [(counts_dict[base] / total_bases) * 100 for base in bases_to_plot]

    plt.figure(figsize=(10, 6))  # canvas size: width = 10, height = 6
    colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"]
    plt.bar(bases_to_plot, frequencies, color=colors[: len(bases_to_plot)])

    plt.title(f"{title} Nucleotide proportion chart", fontsize=14)
    plt.xlabel("Nucleotide Symbol", fontsize=12)
    plt.ylabel("Frequency (%)", fontsize=12)

    # mark the specific value at the top of each bar
    for i, freq in enumerate(frequencies):
        plt.text(i, freq + 0.5, f"{freq:.2f}%", ha="center", fontsize=10)

    plt.ylim(0, max(frequencies) * 1.2)
    plt.grid(axis="y", linestyle="--")

    RESULTS_DIR.mkdir(exist_ok=True)
    slug = title.lower().replace(" ", "_").replace("(", "").replace(")", "")
    output_path = RESULTS_DIR / f"{slug}_distribution.png"
    plt.savefig(output_path)
    plt.close()
    print(f"Graph saved as: {output_path}")


def analyze_sequence(sequence: str, label: str) -> None:
    """Count nucleotide frequencies in a sequence and visualize them.

    Args:
        sequence: The cleaned nucleotide sequence.
        label: Human-readable label for reporting and the chart.
    """
    if not sequence:
        print(f"Warning: {label} sequence is empty.")
        return

    # count how often each symbol occurs in the sequence
    counts: dict[str, int] = {}
    for base in sequence:
        counts[base] = counts.get(base, 0) + 1

    total_len = len(sequence)
    print(f"\n=== {label} quantitative outcome ===")
    print(f"Total Length: {total_len}")
    for base in sorted(counts):
        count = counts[base]
        print(f"  {base}: {count} ({(count / total_len) * 100:.2f}%)")

    visualize_counts(label, counts, total_len)


if __name__ == "__main__":
    file_paths = {
        "Gene A (DNA)": DATA_DIR / "gene_a.fa",
        "RNA A (mRNA)": DATA_DIR / "rna_a.fa",
    }
    for label, path in file_paths.items():
        seq_data = read_file_sequence(path)
        if seq_data:
            analyze_sequence(seq_data, label)
