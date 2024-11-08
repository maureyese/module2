import pandas as pd
import torch
import time
from transformers import RobertaTokenizer, RobertaModel
from sklearn.metrics.pairwise import cosine_similarity
from rdkit import DataStructs
from rdkit.Chem import AllChem, MolFromSmiles
from rdkit.Chem import rdFingerprintGenerator, MolFromSmiles

# Load model and tokenizer
def load_model():
    model_name = "seyonec/ChemBERTa-zinc-base-v1"
    tokenizer = RobertaTokenizer.from_pretrained(model_name)
    model = RobertaModel.from_pretrained(model_name)
    return tokenizer, model

tokenizer, model = load_model()

# Define weights for ChemBERTa and Takimoto similarities
chemberta_weight = 0.8
takimoto_weight = 0.2

# Define functions for SMILES to embedding and similarity calculation
def smiles_to_embedding(smiles):
    if not isinstance(smiles, str):  # Ensure SMILES is a valid string
        return None
    inputs = tokenizer(smiles, return_tensors="pt", padding=True, truncation=True, max_length=128)
    with torch.no_grad():
        outputs = model(**inputs)
    embeddings = outputs.last_hidden_state.mean(dim=1)
    return embeddings

def calculate_chemberta_similarity(embedding1, embedding2):
    if embedding1 is None or embedding2 is None:  # Skip similarity if any embedding is missing
        return 0.0
    return cosine_similarity(embedding1.numpy(), embedding2.numpy())[0][0]

# Create a Morgan Fingerprint generator with specified parameters
morgan_gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048, includeChirality=True)

def calculate_takimoto_similarity(smiles1, smiles2):
    mol1 = MolFromSmiles(smiles1)
    mol2 = MolFromSmiles(smiles2)
    if mol1 is None or mol2 is None:  # Skip if either molecule cannot be generated
        return 0.0
    # Generate fingerprints using the new generator
    fp1 = morgan_gen.GetFingerprint(mol1)
    fp2 = morgan_gen.GetFingerprint(mol2)
    return DataStructs.TanimotoSimilarity(fp1, fp2)

def combined_similarity(smiles1, smiles2):
    # Calculate ChemBERTa similarity
    embedding1 = smiles_to_embedding(smiles1)
    embedding2 = smiles_to_embedding(smiles2)
    chemberta_similarity = calculate_chemberta_similarity(embedding1, embedding2)

    # Calculate Takimoto similarity
    takimoto_similarity = calculate_takimoto_similarity(smiles1, smiles2)

    # Combine using weighted average
    return (chemberta_weight * chemberta_similarity) + (takimoto_weight * takimoto_similarity)

# Load data
reactions_df = pd.read_csv('converted_reactions.csv')
smiles_df = pd.read_csv('smiles.csv')

# Create a dictionary for quick SMILES lookup by code
smiles_dict = dict(zip(smiles_df['code'], smiles_df['smiles']))

# Initialize lists to store new columns
substrate_similarities_list = []
product_similarities_list = []

# Total number of rows for progress tracking
total_rows = len(reactions_df)

# Iterate over each reaction with progress display and timer
start_time = time.time()
for idx, row in reactions_df.iterrows():
    # Calculate elapsed time
    elapsed_time = time.time() - start_time
    print(f"\nProcessing row {idx + 1} out of {total_rows} | Elapsed time: {elapsed_time:.2f} seconds")

    substrate_codes = eval(row['substrates_codes'])  # Convert string representation to list
    product_codes = eval(row['products_codes'])
    substrate_names = eval(row['substrates_names'])
    product_names = eval(row['products_names'])

    # Lists to hold substrate and product similarity scores
    substrate_similarity_list = []
    product_similarity_list = []
    substrate_similarity_output = []
    product_similarity_output = []

    # Calculate similarity scores for each substrate-product pair
    for i, sub_code in enumerate(substrate_codes):
        sub_smiles = smiles_dict.get(sub_code)
        substrate_similarities = []

        for j, prod_code in enumerate(product_codes):
            prod_smiles = smiles_dict.get(prod_code)
            if sub_smiles and prod_smiles:  # Only calculate if both SMILES are available
                combined_score = combined_similarity(sub_smiles, prod_smiles)
                substrate_similarities.append(f"{combined_score:.2f}")
            else:
                # Add 0.0 if either SMILES is missing
                substrate_similarities.append("0.0")

        substrate_similarity_list.append([float(score) for score in substrate_similarities])
        substrate_similarity_output.append(f"{substrate_names[i]} {substrate_similarities}")

    # Print final similarity scores for the row
    print("Substrates:", ", ".join(substrate_similarity_output))
    for i, prod_code in enumerate(product_codes):
        product_similarities = [f"{substrate_similarity_list[sub_idx][i]:.2f}" for sub_idx in range(len(substrate_codes))]
        product_similarity_output.append(f"{product_names[i]} {product_similarities}")
    print("Products:", ", ".join(product_similarity_output))

    # Append results for this row
    substrate_similarities_list.append(substrate_similarity_list)
    product_similarities_list.append(product_similarity_list)

# Add new columns to DataFrame
reactions_df['substrate_similarities'] = substrate_similarities_list
reactions_df['product_similarities'] = product_similarities_list

# Save to CSV
reactions_df.to_csv('reaction_similarity_scores.csv', index=False)
print("Similarity calculation completed. Results saved to 'rdkitchembert_similarity_scores.csv'.")