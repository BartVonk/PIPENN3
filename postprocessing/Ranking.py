import pandas as pd
import sys
import os
import csv
import matplotlib.pyplot as plt
from collections import Counter, defaultdict
import argparse

def parse_and_count(column_values, mcc_values, auc_values):
    counter = Counter()
    mcc_accumulator = defaultdict(list)
    auc_accumulator = defaultdict(list)

    for value, mcc, auc in zip(column_values, mcc_values, auc_values):
        if value:
            items = [item.strip() for item in value.split(',')]
            filtered_items = [item for item in items if item and item.lower() != 'unknown']
            for item in filtered_items:
                counter[item] += 1
                try:
                    mcc_float = float(mcc)
                    mcc_accumulator[item].append(mcc_float)
                except ValueError:
                    pass
                try:
                    auc_float = float(auc)
                    auc_accumulator[item].append(auc_float)
                except ValueError:
                    pass
    return counter, mcc_accumulator, auc_accumulator

def group_subcellular_locations(counter, mcc_accumulator, auc_accumulator):
    supergroup_map = {
        'Nucleus': ['nucleus', 'nucleolus', 'nucleoplasm', 'nuclear'],
        'Cytoplasm': ['cytoplasm', 'cytosol', 'ribosome', 'ribonucleoprotein', 'vesicle'],
        'Membrane': ['membrane', 'plasma membrane', 'inner membrane', 'outer membrane'],
        'Mitochondrion': ['mitochondrion', 'mitochondrial'],
        'Chloroplast/Plastid': ['chloroplast', 'plastid', 'thylakoid', 'apicoplast'],
        'Endoplasmic Reticulum': ['endoplasmic reticulum'],
        'Golgi Apparatus': ['golgi'],
        'Extracellular': ['extracellular', 'secreted', 'cell wall', 'exosome'],
        'Lysosome': ['autolysosome', 'lysosome'],
        'Peroxisome': ['peroxisome'],
        'Synapse': ['synapse', 'postsynaptic'],
        'Virus-Related': ['viral', 'virion', 'capsid', 'nucleocapsid'],
        'Other Organelles': ['proteasome', 'complex', 'chromatin', 'chromosome']
    }

    supergroup_counts = Counter()
    supergroup_mcc = defaultdict(list)
    supergroup_auc = defaultdict(list)

    for label, count in counter.items():
        matched = False
        label_lower = label.lower()
        mcc_list = mcc_accumulator.get(label, [])
        auc_list = auc_accumulator.get(label, [])
        for group, keywords in supergroup_map.items():
            if any(keyword in label_lower for keyword in keywords):
                supergroup_counts[group] += count
                supergroup_mcc[group].extend(mcc_list)
                supergroup_auc[group].extend(auc_list)
                matched = True
                break

    return supergroup_counts, supergroup_mcc, supergroup_auc

def plot_rankings(rankings_list, ranking_type):
    if not rankings_list:
        print(f"No data available to plot for {ranking_type}")
        return

    # Sort by MCC (ascending)
    sorted_items = sorted(rankings_list, key=lambda x: x[1])  # x[1] is avg_mcc
    groups = [item[0] for item in sorted_items]
    avg_mccs = [item[1] for item in sorted_items]
    counts = [item[2] for item in sorted_items]

    fig, ax = plt.subplots(figsize=(12, 6))

    bars = ax.bar(groups, avg_mccs, color='dodgerblue')
    ax.set_xlabel('Group')
    ax.set_ylabel('Average MCC')
    ax.set_title(f'Average MCC per {ranking_type}')
    ax.set_xticks(range(len(groups)))
    ax.set_xticklabels(groups, rotation=45, ha='right')

    for bar, mcc, count in zip(bars, avg_mccs, counts):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width() / 2, height - 0.05, str(count),
                ha='center', va='top', fontsize=9, color='white', fontweight='bold')
        ax.text(bar.get_x() + bar.get_width() / 2, height, f"{mcc:.3f}",
                ha='center', va='bottom', fontsize=8, color='black')

    filename = f'Ranking_{ranking_type.replace(" ", "_")}.png'
    plt.tight_layout()
    plt.savefig(filename)
    plt.close()
    print(f"Saved plot to {filename}")


def print_ranking(header, data, all_entries=False):
    print("#" * 30)
    print(f"Ranking top 10 for {header}")
    print(f"{'Entry':<32} {'avg MCC':>10} {'Count':>8}")
    print("#" * 30)

    if all_entries:
        for rank, (name, avg_mcc, count) in enumerate(data, start=1):
            print(f"{rank}. {name:<30} {avg_mcc:.4f}\t{count}")
    else:
        for rank, (name, avg_mcc, count) in enumerate(data[:10], start=1):
            print(f"{rank}. {name:<30} {avg_mcc:.4f}\t{count}")

        #print("\n" + "#" * 30)
        print(f"Bottom 10 for {header}")
        #print("Entry, avg MCC, Count")
        #print("#" * 30)
        for rank, (name, avg_mcc, count) in enumerate(data[-10:], start=len(data) - 9):
            print(f"{rank}. {name:<30} {avg_mcc:.4f}\t{count}")
    print()

def process_files(filepaths, verbose=False):
    combined_df = pd.concat([pd.read_csv(path) for path in filepaths], ignore_index=True)

    if 'MCC' not in combined_df.columns:
        print("Error: MCC column not found in the file.")
        return

    # CLASS RANKING
    class_df = combined_df[combined_df['Class'] != 'Unknown']
    class_avg_mcc = class_df.groupby('Class')['MCC'].mean().to_dict()
    class_counts = class_df['Class'].value_counts().to_dict()
    class_data = [(name, class_avg_mcc[name], class_counts[name]) for name in sorted(class_avg_mcc, key=class_avg_mcc.get, reverse=True)]

    if verbose:
        print_ranking("Class", class_data, all_entries=True)
    plot_rankings(class_data, "Class")

    # MULTI-VALUE CATEGORIES
    for category in ['pfam_ids', 'supfam_ids']:
        if category in combined_df.columns:
            mcc_values = {}
            counts = {}
            for _, row in combined_df.iterrows():
                if pd.notna(row[category]) and row[category] != "Unknown":
                    entries = str(row[category]).split(',')
                    for entry in entries:
                        entry = entry.strip()
                        if entry:
                            if entry in mcc_values:
                                mcc_values[entry].append(row['MCC'])
                                counts[entry] += 1
                            else:
                                mcc_values[entry] = [row['MCC']]
                                counts[entry] = 1

            avg_mcc = {k: sum(v) / len(v) for k, v in mcc_values.items()}
            data_by_mcc = [(k, avg_mcc[k], counts[k]) for k in avg_mcc]
            data_by_count = sorted(data_by_mcc, key=lambda x: x[2], reverse=True)
            data_by_mcc = sorted(data_by_mcc, key=lambda x: x[1], reverse=True)

            if verbose:
                #print_ranking(f"{category} (by MCC)", data_by_mcc)
                print_ranking(f"{category} (by Count)", data_by_count)

    # SUBCELLULAR LOCATION RANKING
    uniprot_values = combined_df.get("Subcellular Location UniProt", [""] * len(combined_df))
    go_values = combined_df.get("Subcellular location GO", [""] * len(combined_df))
    mcc_values = combined_df["MCC"].tolist()
    auc_values = combined_df.get("AUC", [0] * len(combined_df))

    up_counter, up_mcc, up_auc = parse_and_count(uniprot_values, mcc_values, auc_values)
    go_counter, go_mcc, go_auc = parse_and_count(go_values, mcc_values, auc_values)

    up_super_c, up_super_m, up_super_a = group_subcellular_locations(up_counter, up_mcc, up_auc)
    go_super_c, go_super_m, go_super_a = group_subcellular_locations(go_counter, go_mcc, go_auc)

    up_data = [(k, sum(up_super_m[k]) / len(up_super_m[k]), up_super_c[k]) for k in up_super_c]
    go_data = [(k, sum(go_super_m[k]) / len(go_super_m[k]), go_super_c[k]) for k in go_super_c]

    if verbose:
        print_ranking("Uniprot (Supergrouped)", up_data, all_entries=True)
        print_ranking("GO (Supergrouped)", go_data, all_entries=True)

    #plot_rankings(up_data, "Uniprot_Subcellular_Location")
    plot_rankings(go_data, "GO_Subcellular_Location")

def main():
    parser = argparse.ArgumentParser(description="Rank biological data by MCC and counts across multiple files.")
    parser.add_argument("files", nargs='+', help="Input CSV files")
    parser.add_argument("-v", action="store_true", help="Enable verbose output")
    args = parser.parse_args()
    process_files(args.files, verbose=args.v)

if __name__ == "__main__":
    main()
