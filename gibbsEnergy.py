from equilibrator_api import ComponentContribution, Reaction, Q_
import os
import pandas as pd
import ast

cc = ComponentContribution()

cc.p_h = Q_(7.4)
cc.p_mg = Q_(3.0)
cc.ionic_strength = Q_("0.25M")
cc.temperature = Q_("298.15K")

if os.path.exists('converted_reactions.csv'):
    df_reactions = pd.read_csv('converted_reactions.csv')
else:
    quit()

new_dict = {'id' : [],
            'enzyme' : [],
            'equation' : [],
            'substrates_names' : [],
            'products_names' : [],
            'substrates_codes' : [],
            'products_codes' : [],
            'gibbs_energy' : []}

for idx, row in df_reactions.iterrows(): 
    id = row['id']
    enzyme = row['enzyme']
    equation = row['equation']
    print(f"Processing reaction {id}: {equation}")
    
    substrates = ast.literal_eval(row['substrates_names'])
    products = ast.literal_eval(row['products_names'])
    sub_codes = ast.literal_eval(row['substrates_codes'])
    prod_codes = ast.literal_eval(row['products_codes'])

    # Retrieve substrates
    sub_comp = []
    for cid in sub_codes:
        compound = cc.get_compound(f"kegg:{cid}")
        if compound is None:
            print(f"Warning: Substrate with KEGG ID {cid} not found.")
        else:
            sub_comp.append(compound)
    
    # Retrieve products
    prod_comp = []
    for cid in prod_codes:
        compound = cc.get_compound(f"kegg:{cid}")
        if compound is None:
            print(f"Warning: Product with KEGG ID {cid} not found.")
        else:
            prod_comp.append(compound)

    # Skip if any compound is missing
    if any(comp is None for comp in sub_comp) or any(comp is None for comp in prod_comp):
        print(f"Skipping reaction {id} due to missing compounds.")
        
        standard_dg_prime = 0.0
        
        new_dict['id'].append(id)
        new_dict['enzyme'].append(enzyme)
        new_dict['equation'].append(equation)
        new_dict['substrates_names'].append(substrates)
        new_dict['products_names'].append(products)
        new_dict['substrates_codes'].append(sub_codes)
        new_dict['products_codes'].append(prod_codes)
        new_dict['gibbs_energy'].append(standard_dg_prime)
        continue

    # Construct the reaction
    compound_dict = {}
    for value in sub_comp:
        compound_dict[value] = -1
    for value in prod_comp:
        compound_dict[value] = 1

    try:
        processed_reaction = Reaction(compound_dict)
        print(processed_reaction)
        standard_dg_prime = cc.standard_dg_prime(processed_reaction)
        #standard_dg_prime = str(standard_dg_prime)
        #standard_dg_prime = standard_dg_prime.split("+/-")[0].replace("(", "").strip()
        standard_dg_prime_value = standard_dg_prime.magnitude  # Extract the numerical value
        print(standard_dg_prime_value)
        print(f"Standard ΔG'° for reaction {id}: {standard_dg_prime}")

        new_dict['id'].append(id)
        new_dict['enzyme'].append(enzyme)
        new_dict['equation'].append(equation)
        new_dict['substrates_names'].append(substrates)
        new_dict['products_names'].append(products)
        new_dict['substrates_codes'].append(sub_codes)
        new_dict['products_codes'].append(prod_codes)
        new_dict['gibbs_energy'].append(standard_dg_prime)

    except Exception as e:
        print(f"Error processing reaction {id}: {e}")
        
        standard_dg_prime = 0.0

        new_dict['id'].append(id)
        new_dict['enzyme'].append(enzyme)
        new_dict['equation'].append(equation)
        new_dict['substrates_names'].append(substrates)
        new_dict['products_names'].append(products)
        new_dict['substrates_codes'].append(sub_codes)
        new_dict['products_codes'].append(prod_codes)
        new_dict['gibbs_energy'].append(standard_dg_prime)

# Save the new DataFrame to CSV
df = pd.DataFrame(new_dict)
df.to_csv('gibbs_reactions.csv', index=False)