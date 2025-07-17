import pandas as pd
import json
from urllib.request import urlopen
import sys

def fetch_uniprot(pdb_id, chain_id):
    try:
        url = f"https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pdb_id}"
        response = urlopen(url).read()
        data = json.loads(response.decode('utf-8'))
        mappings = data.get(pdb_id.lower(), {}).get("UniProt", {})
        for uniprot_id, details in mappings.items():
            for mapping in details["mappings"]:
                if mapping["chain_id"] == chain_id:
                    return uniprot_id
    except Exception as e:
        #print(f"Error fetching {pdb_id}{chain_id}: {e}")
        pass
    return None

def main(input_csv):
    df = pd.read_csv(input_csv)

    # Extract PDB and chain from first column
    df['pdb'] = df.iloc[:, 0].astype(str).str[:4]
    df['chain'] = df.iloc[:, 0].astype(str).str[4:]

    # Map to UniProt
    print("Mapping to UniProt...")
    df['mapped_uniprot'] = df.apply(lambda row: fetch_uniprot(row['pdb'], row['chain']), axis=1)

    # Drop rows without mapping
    df = df.dropna(subset=['mapped_uniprot'])

    # Drop duplicates based on UniProt
    df = df.drop_duplicates(subset='mapped_uniprot', keep='first')

    # Replace first column with UniProt IDs
    df.iloc[:, 0] = df['mapped_uniprot']

    # Drop helper columns
    df = df.drop(columns=['pdb', 'chain', 'mapped_uniprot'])

    # Overwrite original file
    df.to_csv(input_csv, index=False)
    #print(f"File successfully overwritten with UniProt-mapped IDs in first column: {input_csv}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python map_pdb_chain_to_uniprot.py <input_csv>")
        sys.exit(1)

    main(sys.argv[1])
