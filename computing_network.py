import torch
import pandas as pd
from transformers import RobertaTokenizer, RobertaModel
import networkx as nx
import numpy as np
from sklearn.metrics.pairwise import cosine_similarity

def load_model():
    model_name = "seyonec/ChemBERTa-zinc-base-v1"
    tokenizer = RobertaTokenizer.from_pretrained(model_name)
    model = RobertaModel.from_pretrained(model_name)
    return model_name, tokenizer, model

model_name, tokenizer, model = load_model()

def smiles_to_embedding(smiles):
    inputs = tokenizer(smiles, return_tensors="pt", padding=True, truncation=True, max_length=128)
    with torch.no_grad():
        outputs = model(**inputs)
    embeddings = outputs.last_hidden_state.mean(dim=1)
    return embeddings

def calculate_similarity(embedding1, embedding2):
    return cosine_similarity(embedding1.numpy(), embedding2.numpy())[0][0]

def open_smiles():
    smiles_df = pd.read_csv('smiles.csv')
    smiles_dict = dict(zip(smiles_df['code'], smiles_df['smiles']))
    return smiles_df, smiles_dict

smiles_df, smiles_dict = open_smiles()

def open_excluded_codes():
    with open('excluded_compounds.txt', 'r') as file:
        excluded_codes = set(line.strip() for line in file.readlines())
    return excluded_codes

excluded_codes = open_excluded_codes()

def open_reactions_df():
    return pd.read_csv('gibbs_reactions.csv')

df = open_reactions_df()

def create_graph():
    G = nx.Graph()
    processed_comparisons = set()  # Track processed comparisons
    processed_rows_count = 0
    total_rows = len(df)  # Total rows for progress tracking

    def process_reaction(row):
        nonlocal processed_rows_count  # Use nonlocal to avoid global
        
        substrates_codes = eval(row['substrates_codes'])
        products_codes = eval(row['products_codes'])
        
        # Retrieve SMILES for substrates and products, filtering excluded codes
        substrates_smiles = [smiles_dict.get(code) for code in substrates_codes if code not in excluded_codes]
        products_smiles = [smiles_dict.get(code) for code in products_codes if code not in excluded_codes]

        if not substrates_smiles or not products_smiles:
            return  # Skip if no valid SMILES left

        max_similarity, best_substrate, best_product = -1, None, None

        for substrate in substrates_smiles:
            for product in products_smiles:
                if substrate and product:
                    substrate_code = next(code for code, smile in smiles_dict.items() if smile == substrate)
                    product_code = next(code for code, smile in smiles_dict.items() if smile == product)
                    comparison_pair = tuple(sorted([substrate_code, product_code]))

                    if comparison_pair in processed_comparisons:
                        continue  # Skip processed pairs

                    processed_comparisons.add(comparison_pair)  # Mark as processed
                    substrate_embedding = smiles_to_embedding(substrate)
                    product_embedding = smiles_to_embedding(product)
                    similarity = calculate_similarity(substrate_embedding, product_embedding)
                    
                    if similarity > max_similarity:
                        max_similarity = similarity
                        best_substrate, best_product = substrate, product

        if best_substrate and best_product and max_similarity >= 0.80:
            best_substrate_code = next(code for code, smile in smiles_dict.items() if smile == best_substrate)
            best_product_code = next(code for code, smile in smiles_dict.items() if smile == best_product)

            if best_substrate_code != best_product_code:
                if not G.has_edge(best_substrate_code, best_product_code):
                    G.add_edge(best_substrate_code, best_product_code, weight=max_similarity)

        processed_rows_count += 1
        print(f"\rProcessed {processed_rows_count}/{total_rows} reactions", end='', flush=True)

    # Process all reactions and return the graph
    df.apply(process_reaction, axis=1)
    print("\nGraph creation complete.")
    return G

G = create_graph()

def find_path_dfs(graph, start_node, end_node, path_length=8):
    def dfs(node, path):
        if len(path) == path_length:
            return path if path[-1] == end_node else None
        for neighbor in graph.neighbors(node):
            if neighbor not in path:
                new_path = path + [neighbor]
                result = dfs(neighbor, new_path)
                if result:
                    return result
        return None
    return dfs(start_node, [start_node])

start_node = "C00011"
end_node = "C00024"
path = find_path_dfs(G, start_node, end_node, path_length=7)
print("Path:", path)

kegg_df = pd.read_csv('kegg_compounds.csv')
code_to_name = kegg_df.groupby('code')['name'].apply(lambda names: ';'.join(names)).to_dict()
compound_names = [code_to_name.get(code, "Unknown") for code in path]

for name in compound_names:
    print(name, "->")

