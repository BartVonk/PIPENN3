import sys
import csv
import os
import matplotlib.pyplot as plt

def plot_mcc_from_csv(csv_path, output_path=None):
    entries = []
    avg_mccs = []
    counts = []

    # Read CSV file
    with open(csv_path, newline='', encoding='utf-8') as f:
        reader = csv.reader(f)
        for row in reader:
            if len(row) < 3:
                continue  # skip incomplete rows
            name = row[0].strip()
            try:
                mcc = float(row[1])
                count = int(row[2])
            except ValueError:
                continue  # skip rows with invalid numbers

            entries.append(name)
            avg_mccs.append(mcc)
            counts.append(count)

    if not entries:
        print("No valid data found.")
        return

    # Sort by MCC
    sorted_data = sorted(zip(entries, avg_mccs, counts), key=lambda x: x[1])
    entries, avg_mccs, counts = zip(*sorted_data)

    # Determine output filename based on CSV
    if output_path is None:
        base = os.path.splitext(os.path.basename(csv_path))[0]
        output_path = f"Ranking_{base}.png"

    # Create plot
    fig, ax = plt.subplots(figsize=(12, 6))
    bars = ax.bar(entries, avg_mccs, color='dodgerblue')

    ax.set_ylabel("Average MCC")
    ax.set_xlabel("Group")
    ax.set_title("Ranking by Average MCC")
    ax.set_xticks(range(len(entries)))
    ax.set_xticklabels(entries, rotation=45, ha='right')

    # Add padding to the top of y-axis
    y_max = max(avg_mccs) + 0.1
    ax.set_ylim(0, y_max)

    # Annotate bars
    for bar, mcc, count in zip(bars, avg_mccs, counts):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width() / 2, height - 0.05, str(count),
                ha='center', va='top', fontsize=9, color='white', fontweight='bold')
        ax.text(bar.get_x() + bar.get_width() / 2, height + 0.01, f"{mcc:.3f}",
                ha='center', va='bottom', fontsize=8, color='black')

    plt.tight_layout()
    plt.savefig(output_path)
    print(f"Plot saved as {output_path}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python plot_mcc_bar.py <csv_path> [output_image.png] [bar_color]")
    else:
        csv_input = sys.argv[1]
        output_file = sys.argv[2] if len(sys.argv) >= 3 else None
        plot_mcc_from_csv(csv_input, output_file)
