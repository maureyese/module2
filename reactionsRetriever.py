from Bio.KEGG import REST

result = REST.kegg_list("reaction").read()

import re
import pandas as pd

# Function to remove the stoichiometric value
def remove_stoichiometric_value(compound):
    # Pattern to match either a number followed by a space, 'n', or '(n+#)' at the beginning
    compound = re.sub(r"^(?:\d+\s+|n\s*|\(n\+\d+\)\s*)", "", compound)
    return compound.strip()

reactions = result.split('\n')

id_list = []
enzyme_list = []
equation_list = []
substrates_list = []
products_list = []

count = 0

for i in range(len(reactions)):

    print(f"{count} out of {len(reactions)}")

    reaction = reactions[i]

    if reaction == '':
        continue

    id, full_equation = reaction.split('\t')

    # Split by ';'
    parts = full_equation.split(';')

    # Store all elements except the last one in 'enzyme', and the last element in 'equation'
    enzyme = ';'.join(parts[:-1]).strip()
    equation = parts[-1].strip()

    substrates, products = equation.split(' <=> ')

    #print(f"Original substrates: {substrates}")
    #print(f"Original products: {products}")

    if ' + ' in substrates.strip():
        substrates = substrates.strip().split(' + ')
        substrates = [compound.strip() for compound in substrates]
        substrates = [remove_stoichiometric_value(compound.strip()) for compound in substrates]
        substrates = [compound for compound in substrates if compound != 'e-']
    else:
        substrates = [remove_stoichiometric_value(substrates.strip())]

    if ' + ' in products.strip():
        products = products.strip().split(' + ')
        products = [compound.strip() for compound in products]
        products = [remove_stoichiometric_value(compound.strip()) for compound in products]
        products = [compound for compound in products if compound != 'e-']
    else:
        products = [remove_stoichiometric_value(products.strip())]

    #print(f"Filtered substrates: {substrates}")
    #print(f"Filtered products: {products}")

    id_list.append(id)
    enzyme_list.append(enzyme)
    equation_list.append(equation)
    substrates_list.append(substrates)
    products_list.append(products)

    count += 1

df = pd.DataFrame({
    'id': id_list,
    'enzyme': enzyme_list,
    'equation': equation_list,
    'substrates': substrates_list,
    'products': products_list
})

df.to_csv('reactions.csv', index=False)