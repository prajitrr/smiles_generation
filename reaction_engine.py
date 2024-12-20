from rdkit import Chem
from rdkit.Chem import AllChem
from stqdm import stqdm
import streamlit as st

EXCLUSION_LIST = ["N", "O", "o", "OCCO"]
BATCH_AMINE_SECONDARY_NC = "[N:2]([#6:5])[#6:3].[C:4](=O)[O;D1:1]>>[O:2].[C:4](=O)[N:1]([#6:5])[#6:3]"
BATCH_AMINE_SECONDARY_CN = "[C:4](=O)[O;D1:1].[N:2]([#6:5])[#6:3]>>[O:2].[C:4](=O)[N:1]([#6:5])[#6:3]"




def run_double_reaction(reaction_SMARTS_string, 
                        number_of_reactant_1, 
                        smiles_list, 
                        names_list, 
                        sample_ID_list, 
                        reaction_name):
    product_list = []
    product_name_list = []
    product_ID_list = []
    for i in stqdm(range(number_of_reactant_1)):
            r1_smiles = smiles_list[i]
            first_name = names_list[i]
            first_ID = sample_ID_list[i]
            for num, j in enumerate(smiles_list[number_of_reactant_1:]):
                second_name = names_list[number_of_reactant_1 + num]
                second_ID = sample_ID_list[number_of_reactant_1 + num]
                r2_smiles = Chem.CanonSmiles(j)
                try:
                    r1_mol = Chem.MolFromSmiles(r1_smiles)
                except Exception as e:
                    print(f"Exception: \"{e}\" at row {i + 1}")
                    #print(r1_smiles)
                    pass
                try:
                    r2_mol = Chem.MolFromSmiles(r2_smiles)
                except Exception as e:
                    print(f"Exception: \"{e}\" at row {number_of_reactant_1 + num + 2}")

                reaction = AllChem.ReactionFromSmarts(reaction_SMARTS_string)
                

                product = reaction.RunReactants((r1_mol, r2_mol))
                if product:
                    for g in product:
                        if Chem.MolToSmiles(g[1]) not in EXCLUSION_LIST:
                            try:
                                product_list.append(Chem.CanonSmiles(Chem.MolToSmiles(g[1])))
                                product_name_list.append(f"{first_name}_{second_name}")
                                product_ID_list.append(f"{first_ID}__SEPARATOR__{second_ID}")
                            except Exception as e:
                                print(e)
                elif reaction_name == "Amidation, Amine First" or reaction_name == "Amidation, Amine Second":
                    if reaction_name == "Amidation, Amine First":
                        secondary_reaction = AllChem.ReactionFromSmarts(BATCH_AMINE_SECONDARY_NC)
                    else:
                        secondary_reaction = AllChem.ReactionFromSmarts(BATCH_AMINE_SECONDARY_CN)
                    product = secondary_reaction.RunReactants((r1_mol, r2_mol))
                    if product:
                        for g in product:
                            if Chem.MolToSmiles(g[1]) not in EXCLUSION_LIST:
                                try:
                                    product_list.append(Chem.CanonSmiles(Chem.MolToSmiles(g[1])))
                                    product_name_list.append(f"{first_name}_{second_name}")
                                    product_ID_list.append(f"{first_ID}__SEPARATOR__{second_ID}")
                                except:
                                    pass
                    else:
                        print("***********************")
                        print(r1_smiles)
                        print(r2_smiles)
                        print("***********************")
                else:
                    print("=======================")
                    print(r1_smiles)
                    print(r2_smiles)
                    print("=======================")

    print(len(product_name_list))
    combined_list = zip(product_ID_list, product_name_list, product_list)
    combined_list = list(set(combined_list))
    product_ID_list_unique = [i[0] for i in combined_list]
    product_name_list_unique = [i[1] for i in combined_list]
    product_list_unique = [i[2] for i in combined_list]
    reaction_ran = True

    st.write(f"Reactions finished running. A total number of " 
             + f"{len(product_name_list)}"
             + f" products (may be non-unique) were obtained from approximately {(len(smiles_list) - number_of_reactant_1) * number_of_reactant_1} reactions ran."
             + f" A total number of {len(product_name_list_unique)} unique products were obtained.")
    

    return product_ID_list_unique, product_name_list_unique, product_list_unique, reaction_ran

def run_single_reaction(reaction_SMARTS_string, 
                        smiles_list, 
                        names_list, 
                        sample_ID_list, 
                        reaction_name):
    product_list = []
    product_name_list = []
    sample_ID_list_out = []
    reaction = AllChem.ReactionFromSmarts(reaction_SMARTS_string)
    if "Methylation" or "methylation" in reaction_name:
        secondary_reactant = Chem.MolFromSmiles("C")
    for i in stqdm(range(len(smiles_list))):
        r1_smiles = smiles_list[i]
        try:
            r1_mol = Chem.MolFromSmiles(r1_smiles)
        except:
            st.write(f"Unable to generate molecule from SMILES at row {i+1}. Skipping row.")
            continue
        if r1_mol:
            product = reaction.RunReactants((r1_mol, secondary_reactant))
        if product:
            for mol in product:
                product_list.append(Chem.CanonSmiles(Chem.MolToSmiles(mol[0])))
                if "Methylation" or "methylation" in reaction_name:
                    product_name_list.append("methyl_"+names_list[i])
                else:
                    product_name_list.append(names_list[i])
                sample_ID_list_out.append(sample_ID_list[i])
    reaction_ran = True
    combined_list = zip(sample_ID_list_out, product_name_list, product_list)
    combined_list = list(set(combined_list))
    sample_ID_list_unique = [i[0] for i in combined_list]
    product_name_list_unique = [i[1] for i in combined_list]
    product_list_unique = [i[2] for i in combined_list]

    st.write(f"Reactions finished running. A total number of " 
             + f"{len(product_name_list)}"
             + f" products (may be non-unique) were obtained from approximately {len(smiles_list)} reactions ran."
             + f" A total number of {len(product_name_list_unique)} unique products were obtained.")
    

    return sample_ID_list_unique, product_name_list_unique, product_list_unique, reaction_ran