#!/bin/bash

# ./GenerateData.sh 1 2
# $1 is pisite /all/ location
# $2 is PDBID file
# to start pipeline run the following code from the desired location: . ../../pipeline/GenerateData.sh ../../data/pisite/all/ BioDLIDs_PDB.txt

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <Pisite map location> <PDBID file>"
    return
fi

# Get the directory where this script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Get the current directory (where the script was *called* from)
CALL_DIR="$PWD"

# structure preparation
mkdir -p data

echo 'Assembling Fasta-like file from PiSite data'
# pisite to fasta
#python "$SCRIPT_DIR/scripts/PiSITE_to_FASTA.py" $1 "$CALL_DIR/data/RAW.fasta"

echo 'Cleaning and Filtering data by length'
# fasta to filtered fasta, minimum sequence length is 26
python "$SCRIPT_DIR/scripts/generalFilter.py" "$CALL_DIR/data/RAW.fasta" "$CALL_DIR/data/filtered.fasta" 26 1023
#rm "$CALL_DIR/data/RAW.fasta"

# parse filtered.fasta to filtered.csv
python "$SCRIPT_DIR/scripts/parser.py" "$CALL_DIR/data/filtered.fasta" "$CALL_DIR/data/filtered.csv" -p
# use duplicatesFilter.py, uses pandas, so loading an env that has pandas goes first.
(
	source "/net/sys/pscst000/export/informatica/users/bvo207/pipenn-exp/bart/venv/bin/activate"
	python "$SCRIPT_DIR/scripts/duplicatesFilter.py" "$CALL_DIR/data/filtered.csv" "$CALL_DIR/data/filtered.fasta"
	deactivate
)

#rm "$CALL_DIR/data/filtered.csv"

echo 'Clustering data'
# filtered fasta to clustered fasta

awk 'NR % 3 != 0' "$CALL_DIR/data/filtered.fasta" > "$CALL_DIR/data/filtered_noInterface.fasta"
(
	module use /opt/site-apps/zen2/modules/all
	module load netsurfp
	mmseqs easy-linclust "$CALL_DIR/data/filtered_noInterface.fasta" "$CALL_DIR/data/clusterRes" "$CALL_DIR/data/tmp" --min-seq-id 0.25 --cov-mode 1 -c 0.9
	#rm "$CALL_DIR/data/filtered_noInterface.fasta"
	#rm "$CALL_DIR/data/clusterRes_cluster.tsv"
	#rm "$CALL_DIR/data/clusterRes_rep_seq.fasta"

	echo 'Removing mapped clusters'
	python "$SCRIPT_DIR/scripts/eliminateClusters.py" "$CALL_DIR/data/clusterRes_all_seqs.fasta" $2 "$CALL_DIR/data/filtered_mapped.fasta"
	#rm "$CALL_DIR/data/clusterRes_all_seqs.fasta"

	echo 'Fetching representatives'

	mmseqs easy-linclust "$CALL_DIR/data/filtered_mapped.fasta" "$CALL_DIR/data/clusterRes2" "$CALL_DIR/data/tmp" --min-seq-id 0.25 --cov-mode 1 -c 0.9
	module unload netsurfp
)
#rm "$CALL_DIR/data/filtered_mapped.fasta"
#rm "$CALL_DIR/data/clusterRes2_all_seqs.fasta"
#rm "$CALL_DIR/data/clusterRes2_cluster.tsv"

# create ID lists (Pisite, filterding)
grep '^>' "$CALL_DIR/data/clusterRes2_rep_seq.fasta" | sed 's/^>//' | sed 's/^[[:space:]]*//;s/[[:space:]]*$//' > "$CALL_DIR/data/clusterRes_rep_IDs.txt"
grep '^>' ../PiSite/PiSITE_filtered.fasta | sed 's/^>//' > "$CALL_DIR/data/PiSITE_IDs.txt"
sort "$CALL_DIR/data/PiSITE_IDs.txt" > "$CALL_DIR/data/sorted_PiSITE_IDs.txt"
sort "$CALL_DIR/data/clusterRes_rep_IDs.txt" > "$CALL_DIR/data/sorted_clusterRes_rep_IDs.txt"
# create IDs_no_overlap.txt
comm -23 "$CALL_DIR/data/sorted_PiSITE_IDs.txt" "$CALL_DIR/data/sorted_clusterRes_rep_IDs.txt" > "$CALL_DIR/data/overlapping_ids.txt"
#rm "$CALL_DIR/data/PiSITE_IDs.txt"
#rm "$CALL_DIR/data/clusterRes_rep_IDs.txt"
#rm "$CALL_DIR/data/clusterRes_rep_seq.fasta"
#rm "$CALL_DIR/data/sorted_PiSITE_IDs.txt"
#rm "$CALL_DIR/data/sorted_clusterRes_rep_IDs.txt"

# removing IDs_no_overlap from filterred fasta
echo 'Constructing dataset'
python "$SCRIPT_DIR/scripts/filterByID.py" "$CALL_DIR/data/overlapping_ids.txt" "$CALL_DIR/data/filtered.fasta" "$CALL_DIR/data/prepared_pisite.fasta"
#rm "$CALL_DIR/data/filtered.fasta"
#rm "$CALL_DIR/data/overlapping_ids.txt"

echo 'Parsing to csv'
# clustered fasta to csv
python "$SCRIPT_DIR/scripts/parser.py" "$CALL_DIR/data/prepared_pisite.fasta" "$CALL_DIR/prepared_pisite.csv" -p
#rm "$CALL_DIR/data/prepared_pisite.fasta"

python "$SCRIPT_DIR/scripts/split7030.py" "$CALL_DIR/prepared_pisite.csv"
# rm "$CALL_DIR/prepared_pisite.csv"

python "$SCRIPT_DIR/scripts/updateMapping.py" "prepared_pisite_testing.csv"
# rm -r "$CALL_DIR/data/"
