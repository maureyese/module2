import os
import pandas as pd
import ast

if os.path.exists('kegg_compounds.csv'):
    df = pd.read_csv('kegg_compounds.csv')
else:
    quit()

if os.path.exists('reactions.csv'):
    df_reactions = pd.read_csv('reactions.csv')
else:
    quit()

new_dict = {'id' : [],
            'enzyme' : [],
            'equation' : [],
            'substrates_names' : [],
            'products_names' : [],
            'substrates_codes' : [],
            'products_codes' : []}

print(df_reactions['substrates'].head(), type(df_reactions['substrates'][0]))
print(df_reactions['products'].head(), type(df_reactions['products'][0]))

# Iterate over reaction rows
for idx, row in df_reactions.iterrows():
    id = row['id']
    enzyme = row['enzyme']
    equation = row['equation']
    
    # Convert strings to lists
    substrates = ast.literal_eval(row['substrates'])
    products = ast.literal_eval(row['products'])

    try:
        sub_codes = []
        for substrate in substrates:
            if substrate in df['name'].values:
                code_value = df.loc[df['name'] == substrate, 'code'].values
                sub_codes.append(code_value[0])
            else:
                raise Exception("No compound found.")
    except Exception as e:
        print(f'Skipping {equation}...')
        continue

    try:
        prod_codes = []
        for product in products:
            if product in df['name'].values:
                code_value = df.loc[df['name'] == product, 'code'].values
                prod_codes.append(code_value[0])
            else:
                raise Exception("No compound found.")
    except Exception as e:
        print(f'Skipping {equation}...')
        continue

    # Debug prints to trace computation
    print(f'Subcodes: {substrates} -> {sub_codes}')
    print(f'Prodcodes: {products} -> {prod_codes}')

    # Append data to the new dictionary
    new_dict['id'].append(id)
    new_dict['enzyme'].append(enzyme)
    new_dict['equation'].append(equation)
    new_dict['substrates_names'].append(substrates)
    new_dict['products_names'].append(products)
    new_dict['substrates_codes'].append(sub_codes)
    new_dict['products_codes'].append(prod_codes)

# Save the new DataFrame to CSV
df = pd.DataFrame(new_dict)
df.to_csv('converted_reactions.csv', index=False)