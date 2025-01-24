from rdkit import Chem
from rdkit.Chem import AllChem
from stqdm import stqdm
import streamlit as st
import numpy as np

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
                if reaction_name == "Amidation, Amine First" or reaction_name == "Amidation, Amine Second":
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
                        pass
                else:
                    pass

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
    else:
        raise Exception("Reaction not recognized.")
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

def run_multi_double_reaction(reaction_SMARTS_string, 
                              number_of_reactant_1, 
                              smiles_list, 
                              names_list, 
                              sample_ID_list, 
                              reaction_name, 
                              multireact_reactant_1, 
                              multireact_reactant_2):
    if not multireact_reactant_1 and not multireact_reactant_2:
        return run_single_reaction(reaction_SMARTS_string, 
                                   smiles_list, 
                                   names_list, 
                                   sample_ID_list, 
                                   reaction_name)
    

    product_list = []
    product_name_list = []
    product_ID_list = []
    reaction_count = 0

    r1_smiles_list = smiles_list[:number_of_reactant_1]
    r2_smiles_list = smiles_list[number_of_reactant_1:]

    r1_names_list = names_list[:number_of_reactant_1]
    r2_names_list = names_list[number_of_reactant_1:]

    r1_ID_list = sample_ID_list[:number_of_reactant_1]
    r2_ID_list = sample_ID_list[number_of_reactant_1:]

    reaction = AllChem.ReactionFromSmarts(reaction_SMARTS_string)

    # Get first reactant from SMARTS String
    reactant_1 = reaction.GetReactants()[0].GetSmarts()

    # Get second reactant from SMARTS String
    reactant_2 = reaction.GetReactants()[1].GetSmarts()
    
    if multireact_reactant_1 or multireact_reactant_2:

        if multireact_reactant_1:
            r1_smarts_match_counts = np.zeros(len(r1_smiles_list))

        if multireact_reactant_2:
            r2_smarts_match_counts = np.zeros(len(r2_smiles_list))

        # Loop over R1 and R2 set of reactants, count number of matches, perform 
        # first set of reactions
        for i in stqdm(range(r1_smiles_list)):
            r1_mol = Chem.MolFromSmiles(r1_smiles_list[i])
            
            
            if multireact_reactant_1:
                num_matches_r1 = len(r1_mol.GetSubstructMatches(reactant_1))
                r1_smarts_match_counts[i] = num_matches_r1

            first_name = r1_names_list[i]
            first_ID = r1_ID_list[i]

            for j in range(len(r2_smiles_list)):
                r2_mol  = Chem.MolFromSmiles(r2_smiles_list[j])


                if multireact_reactant_2 and j == 0:
                    num_matches_r2 = len(r2_mol.GetSubstructMatches(reactant_2))    
                    r2_smarts_match_counts[j] = num_matches_r2
                
                second_name = r2_names_list[j]
                second_ID = r2_ID_list[j]
                                
                product = reaction.RunReactants((r1_mol, r2_mol))
                reaction_count += 1
                
                if product:
                    for g in product:
                        if Chem.MolToSmiles(g[1]) not in EXCLUSION_LIST:
                            try:
                                product_list.append(Chem.CanonSmiles(Chem.MolToSmiles(g[1])))
                                product_name_list.append(f"{first_name}_{second_name}")
                                product_ID_list.append(f"{first_ID}__SEPARATOR__{second_ID}")
                            except Exception as e:
                                print(e)
                if reaction_name == "Amidation, Amine First" or reaction_name == "Amidation, Amine Second":
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



    
        
