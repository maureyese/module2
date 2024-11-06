import torch
import pandas as pd
from transformers import RobertaTokenizer, RobertaModel
import networkx as nx
import numpy as np
from sklearn.metrics.pairwise import cosine_similarity

# Load the ChemBERTa model and tokenizer
model_name = "seyonec/ChemBERTa-zinc-base-v1"  # Pre-trained ChemBERTa model
tokenizer = RobertaTokenizer.from_pretrained(model_name)  # Use RobertaTokenizer
model = RobertaModel.from_pretrained(model_name)  # Use RobertaModel

def smiles_to_embedding(smiles):
    """Convert SMILES to an embedding using ChemBERTa."""
    # Tokenize SMILES
    inputs = tokenizer(smiles, return_tensors="pt", padding=True, truncation=True, max_length=128)
    with torch.no_grad():
        outputs = model(**inputs)

    # Get the embeddings from the last hidden state
    embeddings = outputs.last_hidden_state.mean(dim=1)  # Average over all tokens
    return embeddings

def calculate_similarity(embedding1, embedding2):
    """Calculate cosine similarity between two embeddings."""
    return cosine_similarity(embedding1.numpy(), embedding2.numpy())[0][0]

# Load the SMILES data
smiles_df = pd.read_csv('smiles.csv')
smiles_dict = dict(zip(smiles_df['code'], smiles_df['smiles']))  # Map codes to SMILES

# Read the excluded compound codes from the 'excluded_compounds.txt' file
with open('excluded_compounds.txt', 'r') as file:
    excluded_codes = set(line.strip() for line in file.readlines())  # Store excluded codes in a set

# Load the Gibbs reactions DataFrame
df = pd.read_csv('gibbs_reactions.csv')

# Create a graph
G = nx.Graph()

# Initialize a set to keep track of processed comparisons
processed_comparisons = set()

def process_reaction(row):
    substrates_codes = eval(row['substrates_codes'])  # Convert string to list
    products_codes = eval(row['products_codes'])  # Convert string to list
    
    # Get the SMILES strings from the codes
    substrates_smiles = [smiles_dict.get(code) for code in substrates_codes]
    products_smiles = [smiles_dict.get(code) for code in products_codes]

    # Remove any excluded compounds from substrates and products
    substrates_smiles = [smiles for code, smiles in zip(substrates_codes, substrates_smiles) if code not in excluded_codes]
    products_smiles = [smiles for code, smiles in zip(products_codes, products_smiles) if code not in excluded_codes]

    # If there are no valid substrates or products left after removal, skip this reaction
    if not substrates_smiles or not products_smiles:
        return
    
    # Compare substrates with products and find the highest similarity
    max_similarity = -1
    best_substrate = None
    best_product = None

    for substrate in substrates_smiles:
        for product in products_smiles:
            if substrate and product:  # Ensure valid SMILES strings
                # Normalize the comparison by using a sorted pair of codes
                substrate_code = [code for code, smile in smiles_dict.items() if smile == substrate][0]
                product_code = [code for code, smile in smiles_dict.items() if smile == product][0]

                # Create a unique comparison identifier (use sorted order to avoid direction bias)
                comparison_pair = tuple(sorted([substrate_code, product_code]))

                # Skip this comparison if it has already been processed
                if comparison_pair in processed_comparisons:
                    continue
                
                # Mark this comparison as processed
                processed_comparisons.add(comparison_pair)

                # Calculate the similarity score
                substrate_embedding = smiles_to_embedding(substrate)
                product_embedding = smiles_to_embedding(product)
                similarity = calculate_similarity(substrate_embedding, product_embedding)
                
                if similarity > max_similarity:
                    max_similarity = similarity
                    best_substrate = substrate
                    best_product = product
    
    # Only add edges with similarity >= 0.80 and prevent self-loops
    if (best_substrate and best_product and max_similarity >= 0.80):
        best_substrate_code = [code for code, smile in smiles_dict.items() if smile == best_substrate][0]
        best_product_code = [code for code, smile in smiles_dict.items() if smile == best_product][0]
        
        # Prevent self-loops by checking that the substrate code is not the same as the product code
        if best_substrate_code != best_product_code:
            # Check if the edge already exists (either direction)
            if not G.has_edge(best_substrate_code, best_product_code):
                print(f'{best_substrate_code} with {best_product_code} with similarity of {max_similarity}')
                G.add_edge(best_substrate_code, best_product_code, weight=max_similarity)

# Process all reactions
df.apply(process_reaction, axis=1)

# Function to find a path of length 8 using DFS
def find_path_dfs(graph, start_node, path_length=8):
    """Find a path of the given length using DFS."""
    def dfs(node, path):
        if len(path) == path_length:
            return path
        for neighbor in graph.neighbors(node):
            if neighbor not in path:
                new_path = path + [neighbor]
                result = dfs(neighbor, new_path)
                if result:
                    return result
        return None
    
    return dfs(start_node, [start_node])

# Example of finding a path of length 8 from a starting node
start_node = list(G.nodes)[0]  # Choose any starting node
path = find_path_dfs(G, start_node, path_length=8)

print("Path of length 8:", path)