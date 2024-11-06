import pandas as pd

# Load the CSV file into a DataFrame
df = pd.read_csv('smiles.csv')

# Filter rows where the 'smiles' column is empty or NaN
excluded_df = df[df['smiles'].isnull() | (df['smiles'] == '')]

# Open the excluded_compounds.txt file in append mode
with open('excluded_compounds.txt', 'a') as file:
    # Append the 'code' and 'smiles' columns of the excluded rows to the file
    for index, row in excluded_df.iterrows():
        file.write(f"{row['code']}\n")

print(f"Excluded compounds have been added to 'excluded_compounds.txt'.")
