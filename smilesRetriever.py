import re
import pubchempy as pcp
from Bio.KEGG import REST
import pandas as pd
import os
import time

if os.path.exists('kegg_compounds.csv'):
    print('KEGG compounds DataFrame already found.')
    kegg_df = pd.read_csv('kegg_compounds.csv')
    # Retrieve unique codes from DataFrame
    unique_codes = kegg_df['code'].unique().tolist()
else:
    # Fetch the list of KEGG compounds
    compounds = REST.kegg_list("compound").read()
    compounds = compounds.split("\n")

    # Create a list to store the rows for the DataFrame
    rows = []
    unique_codes = []

    # Process each compound entry
    for compound in compounds:
        if '\t' in compound:
            code, name = compound.split('\t')

            if code not in unique_codes:
                unique_codes.append(code)

            # Handle alternative names if present
            if ';' in name:
                alt_names = name.split(';')
                for alt_name in alt_names:
                    # Strip whitespace and append each name to the DataFrame
                    rows.append({'Code': code, 'Name': alt_name.strip()})
            else:
                # Append the single name to the DataFrame
                rows.append({'Code': code, 'Name': name.strip()})

    # Create a DataFrame from the rows list
    df = pd.DataFrame(rows)

    # Save the DataFrame to a CSV file
    df.to_csv('kegg_compounds.csv', index=False)

    print("KEGG compounds saved to 'kegg_compounds.csv'.")

if os.path.exists('smiles.csv'):
    df = pd.read_csv("smiles.csv")
    # Get IDs from 'code' column
    ids = df['code'].tolist()
    # Modify unique codes list
    unique_codes = [code for code in unique_codes if code not in ids]
else:
    df = pd.DataFrame(columns=['code', 'smiles'])

def get_smiles_from_pubchem(pubchem_id):
    # Retrieve the compound using the PubChem ID
    c = pcp.Compound.from_cid(int(pubchem_id))
    # Check if the compound was found
    return c.isomeric_smiles if c else None

count = 0

for compound in unique_codes:
    if count == 12000:
        break

    print(f"{count} out of {len(unique_codes)}")
    comp_info = REST.kegg_get(compound).read()
    pattern = r"PubChem:\s*(\d+)"

    # Search for the pattern in the text
    match = re.search(pattern, comp_info)

    # Check if a match was found
    if match:
        try:
            pubchem_id = match.group(1)  # Retrieve the value (digits)
            smiles = get_smiles_from_pubchem(pubchem_id)

            if smiles:
                print(f"Successfully associated {compound} with PubChem ID {pubchem_id} and SMILES.")
                # Increment count only if SMILES is found
                count += 1
            else:
                print(f"No SMILES found for PubChem ID {pubchem_id}.")
                smiles = None
        except:
            print(f"No SMILES found for PubChem ID {pubchem_id}.")
            smiles = None
    else:
        print(f"PubChem ID not found for: {compound}.")
        smiles = None

    # Create a provisional DataFrame and concatenate
    prov_df = pd.DataFrame({'code': [compound], 'smiles': [smiles]})
    df = pd.concat([df, prov_df], ignore_index=True)

    # Save progress
    df.to_csv('smiles.csv', index=False)

print("Final DataFrame:")
print(df.shape)
