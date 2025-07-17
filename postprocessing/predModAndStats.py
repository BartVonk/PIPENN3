import csv
import sys
import requests
import sklearn.metrics
import numpy as np
#import pandas as pd

def fetch_uniprot_info(uniprot_id):
    # Calls UniProt REST API to use UniProt ID for taxonomy info
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
    response = requests.get(url)
    
    if response.status_code == 200:
        data = response.json()
        lineage = data.get("organism", {}).get("lineage", [])
        protein = data.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value", "Unknown")
        
        if not protein or protein == "Unknown":
            protein = data.get("proteinDescription", {}).get("submissionNames", [{}])[0].get("fullName", {}).get("value", "Unknown")
        
        # Extract Pfam and SUPFAM IDs
        pfam_ids = []
        supfam_ids = []
        pfam_ids = list({entry["id"] for entry in data.get("uniProtKBCrossReferences", []) if entry.get("database") == "Pfam"}) or "Unknown"
        supfam_ids = list({entry["id"] for entry in data.get("uniProtKBCrossReferences", []) if entry.get("database") == "SUPFAM"}) or "Unknown"
        
        # Convert lists to strings if only one ID exists
        pfam_ids = ", ".join(pfam_ids) if isinstance(pfam_ids, list) else pfam_ids
        supfam_ids = ", ".join(supfam_ids) if isinstance(supfam_ids, list) else supfam_ids

        
        return {
            "root": lineage[0] if len(lineage) > 0 else "Unknown",
            "class": lineage[1] if len(lineage) > 1 else "Unknown",
            "fold": lineage[2] if len(lineage) > 2 else "Unknown",
            "superfamily": lineage[3] if len(lineage) > 3 else "Unknown",
            "family": lineage[4] if len(lineage) > 4 else "Unknown",
            "protein": protein if protein else "Unknown",
            "species": data.get("organism", {}).get("scientificName", "Unknown"),
            "pfam_ids": pfam_ids,
            "supfam_ids": supfam_ids
        }
    
    return {
        "root": "Unknown", "class": "Unknown", "fold": "Unknown", "superfamily": "Unknown", "family": "Unknown", 
        "protein": "Unknown", "species": "Unknown", "pfam_ids": "Unknown", "supfam_ids": "Unknown"
    }

def fetch_ebi_info(uniprot_id):
    url = f"https://www.ebi.ac.uk/proteins/api/proteins/{uniprot_id}.json"
    response = requests.get(url)

    if response.status_code == 200:
        data = response.json()

        # Extract subcellular location from comments
        subcellular_location_uniprot = []
        for comment in data.get("comments", []):
            if comment.get("type") == "SUBCELLULAR_LOCATION":
                for location in comment.get("locations", []):
                    loc_value = location.get("location", {}).get("value")
                    if loc_value:
                        subcellular_location_uniprot.append(loc_value)

        # Extract GO terms starting with "C:" from dbReferences
        subcellular_location_GO = []
        for ref in data.get("dbReferences", []):
            if ref.get("type") == "GO":
                term = ref.get("properties", {}).get("term", "")
                if term.startswith("C:"):
                    subcellular_location_GO.append(term[2:].strip())

        # Format lists to comma-separated strings or default to "Unknown"
        subcellular_location_uniprot = ", ".join(subcellular_location_uniprot) if subcellular_location_uniprot else "Unknown"
        subcellular_location_GO = ", ".join(subcellular_location_GO) if subcellular_location_GO else "Unknown"

        return {
            "Subcellular location UniProt": subcellular_location_uniprot,
            "Subcellular location GO": subcellular_location_GO
        }

    return {
        "Subcellular location UniProt": "Unknown",
        "Subcellular location GO": "Unknown"
    }

def compute_aa_statistics(sequence):
    # AA sequence properties and classification
    # divide AA into sets based on property
    hydrophobic = set("AILMFWVY")
    charged = set("DEKR")
    aromatic = set("FWY")
    sulfur = set("CM")
    
    total = len(sequence)
    # count AA property groups
    hydrophobic_count = sum(1 for aa in sequence if aa in hydrophobic)
    charged_count = sum(1 for aa in sequence if aa in charged)
    aromatic_count = sum(1 for aa in sequence if aa in aromatic)
    sulfur_count = sum(1 for aa in sequence if aa in sulfur)
    
    if total == 0:
        return (0, "Unknown", 0, "Unknown", 0, "Unknown", 0, "Unknown")
    
    # calculate percentage of AA property group with respect to total length
    hydrophobic_percent = round((hydrophobic_count / total) * 100, 3)
    charged_percent = round((charged_count / total) * 100, 3)
    aromatic_percent = round((aromatic_count / total) * 100, 3)
    sulfur_percent = round((sulfur_count / total) * 100, 3)
    
    # simple classification / educated guess based on calculated percentages
    #hydrophobic_class = "Membrane Protein" if hydrophobic_percent > 50 else "Typical Globular Protein" if 40 <= hydrophobic_percent <= 50 else "Disordered Protein"
    #charged_class = "DNA/RNA-binding protein" if charged_percent > 20 else "Typical Enzyme" if 10 <= charged_percent <= 20 else "Membrane Protein"
    #aromatic_class = "Rich" if aromatic_percent > 10 else "Typical" if 5 <= aromatic_percent <= 10 else "Poor"
    #sulfur_class = "High" if sulfur_percent > 3 else "Low"
    
    return (hydrophobic_percent, charged_percent, aromatic_percent, sulfur_percent)

    #return (hydrophobic_percent, hydrophobic_class, charged_percent, charged_class, aromatic_percent, aromatic_class, sulfur_percent, sulfur_class)

def equalCutoff(y_trues, y_preds):
    # cutoff is calculated based on: number of actual positives == number of predicted positives
    numActualPoses = np.count_nonzero(y_trues)
    sortedPreds = np.sort(y_preds)
    cutoffInd = sortedPreds.size - numActualPoses
    cutoff = sortedPreds[cutoffInd]
    return cutoff

def apply_cutoff_conf(input_file, output_file):
    # confusion matrix
    
    with open(input_file, newline='', encoding='utf-8') as infile, open(output_file, 'w', newline='', encoding='utf-8') as outfile:
        reader = csv.reader(infile)
        writer = csv.writer(outfile)
        
        headers = next(reader) # ignore original headers
        # extend original headers to support new information
        headers.extend(["y_dist", "y_dist_mean", "% Hydrophobic", "% Charged", "% Aromatic", "% Sulfur", "TP", "FP", "TN", "FN", "Accuracy", "Precision", "Recall", "F1 Score", "AUC", "MCC", "Specificity", "FPR", "FNR", "FDR", "Root", "Class", "Fold", "Superfamily", "Family", "Protein", "Species", "Subcellular Location UniProt", "Subcellular location GO", "pfam_ids", "supfam_ids"])
        writer.writerow(headers)
        
        processed_rows = []
        for row in reader:
            uniprot_id = row[0] # uniprot ids in col 1
            sequence = row[1].replace(',', '') # sequence in col 2
            aa_stats = compute_aa_statistics(sequence) # fetches amino acid statisics based on sequence
            info = fetch_uniprot_info(uniprot_id) # fetches taxonomy info based on uniprot id
            ebi_info = fetch_ebi_info(uniprot_id) # fetches subcellular location if available

            y_dist = np.abs(np.subtract(np.array(list(map(float, row[2].split(',')))), np.array(list(map(float, row[3].split(','))))))
            y_dist_mean = y_dist
            y_dist = ','.join(map(str, y_dist))
            y_dist_mean = y_dist_mean.mean()

            # making y_trues human readable            
            third_col_values = row[2].split(',')
            int_values = [str(int(float(v))) for v in third_col_values]
            row[2] = ','.join(int_values)
            y_trues = list(map(int, row[2].split(',')))

            # making y_pred human readable by applying a provided cutoff value
            values = row[3].split(',')
            rounded_values = [str(1 if float(v) >= float(equalCutoff(y_trues, values)) else 0) for v in values]
            row[3] = ','.join(rounded_values)
            y_preds = list(map(int, row[3].split(',')))
            
            # calculating confusion matrix
            TP = sum(1 for yt, yp in zip(y_trues, y_preds) if yt == 1 and yp == 1)
            FP = sum(1 for yt, yp in zip(y_trues, y_preds) if yt == 0 and yp == 1)
            TN = sum(1 for yt, yp in zip(y_trues, y_preds) if yt == 0 and yp == 0)
            FN = sum(1 for yt, yp in zip(y_trues, y_preds) if yt == 1 and yp == 0)
            # calculating statistics based on confusion matrix
            Accuracy = (TP + TN) / (TP + FP + TN + FN) if (TP + FP + TN + FN) > 0 else 0
            Precision = TP / (TP + FP) if (TP + FP) > 0 else 0
            Recall = TP / (TP + FN) if (TP + FN) > 0 else 0
            F1_Score = 2 * (Precision * Recall) / (Precision + Recall) if (Precision + Recall) > 0 else 0
            Specificity = TN / (TN + FP) if (TN + FP) > 0 else 0
            FPR = FP / (FP + TN) if (FP + TN) > 0 else 0
            FNR = FN / (FN + TP) if (FN + TP) > 0 else 0
            FDR = FP / (FP + TP) if (FP + TP) > 0 else 0

            AUC = sklearn.metrics.roc_auc_score(y_trues, y_preds, average="weighted")
            MCC = sklearn.metrics.matthews_corrcoef(y_trues, y_preds, sample_weight=None)
            
            row.append(y_dist)
            row.append(y_dist_mean)
            row.extend(aa_stats)
            row.extend([TP, FP, TN, FN, round(Accuracy, 3), round(Precision, 3), round(Recall, 3), round(F1_Score, 3), round(AUC, 3), round(MCC, 3), round(Specificity, 3), round(FPR, 3), round(FNR, 3), round(FDR, 3), info["root"], info["class"], info["fold"], info["superfamily"], info["family"], info["protein"], info["species"], ebi_info["Subcellular location UniProt"], ebi_info["Subcellular location GO"], info["pfam_ids"], info["supfam_ids"]])
            processed_rows.append(row)
        clean_protein_column(processed_rows)
        writer.writerows(processed_rows)

def clean_protein_column(rows):
    # cleaning up taxonomy data, since it sometimes returns (part of) a dictionary
    for row in rows:
        try:
            protein_value = row[-2]
            if isinstance(protein_value, dict):
                row[-2] = protein_value.get("value", "Unknown")
                print(row[-2])
            elif not protein_value or protein_value == "{}":
                row[-2] = "Unknown"
        except IndexError:
            pass  # Ensures we don't break in case of missing values

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python predModAndStats.py input.csv output.csv")
    else:
        apply_cutoff_conf(sys.argv[1], sys.argv[2])

