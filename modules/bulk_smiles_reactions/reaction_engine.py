from rdkit import Chem
from rdkit.Chem import AllChem
from stqdm import stqdm
import streamlit as st
import numpy as np
import scipy as sp
import pandas as pd
import time

EXCLUSION_LIST = ["N", "O", "o", "OCCO"]
BATCH_AMINE_SECONDARY_NC = "[N:2]([#6:5])[#6:3].[C:4](=O)[O;D1:1]>>[O:2].[C:4](=O)[N:1]([#6:5])[#6:3]"
BATCH_AMINE_SECONDARY_CN = "[C:4](=O)[O;D1:1].[N:2]([#6:5])[#6:3]>>[O:2].[C:4](=O)[N:1]([#6:5])[#6:3]"




def run_double_reaction(reaction_SMARTS_string, 
                        number_of_reactant_1, 
                        smiles_list, 
                        names_list, 
                        sample_ID_list, 
                        reaction_name):
    
    start_time = time.time()

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

    end_time = time.time()
    total_time = end_time - start_time

    if total_time > 60:
        st.write(f"Total time taken: {total_time / 60} minutes.")
    elif total_time < 10e-3:
        st.write(f"Total time taken: {total_time * 10**3} milliseconds.")
    else:
        st.write(f"Total time taken: {end_time - start_time} seconds.")
    
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
    print(multireact_reactant_1)
    print(multireact_reactant_2)
    start_time = time.time()
    if not multireact_reactant_1 and not multireact_reactant_2:
        return run_double_reaction(reaction_SMARTS_string, 
                                   number_of_reactant_1,
                                   smiles_list, 
                                   names_list, 
                                   sample_ID_list, 
                                   reaction_name)
    
    if reaction_name == "Amidation, Amine First":
        secondary_reaction = AllChem.ReactionFromSmarts(BATCH_AMINE_SECONDARY_NC)
    elif reaction_name == "Amidation, Amine Second":
        secondary_reaction = AllChem.ReactionFromSmarts(BATCH_AMINE_SECONDARY_CN)
    else:
        secondary_reaction = None

    product_list = []
    product_name_list = []
    product_ID_list = []
    reaction_count = 0
    total_products = 0

    r1_smiles_list = list(smiles_list[:number_of_reactant_1])
    r2_smiles_list = list(smiles_list[number_of_reactant_1:])
    print(len(r2_smiles_list))
    print(type(r2_smiles_list))

    r1_names_list = list(names_list[:number_of_reactant_1])
    r2_names_list = list(names_list[number_of_reactant_1:])

    r1_ID_list = list(sample_ID_list[:number_of_reactant_1])
    r2_ID_list = list(sample_ID_list[number_of_reactant_1:])

    reaction = AllChem.ReactionFromSmarts(reaction_SMARTS_string)

    # Get first reactant from SMARTS String
    reactant_1 = reaction.GetReactants()[0]

    # Get second reactant from SMARTS String
    reactant_2 = reaction.GetReactants()[1]
    
    if multireact_reactant_1:
        product_matches_r1 = []
        r1_smarts_match_counts = np.zeros(len(r1_smiles_list))

    if multireact_reactant_2:
        product_matches_r2 = []
        r2_smarts_match_counts = np.zeros(len(r2_smiles_list))

    # Loop over R1 and R2 set of reactants, count number of matches, perform 
    # first set of reactions
    for i in stqdm(range(len(r1_smiles_list))):
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
            

            if product:
                reaction_count += 1
                for p in product:
                    product_smiles = Chem.CanonSmiles(Chem.MolToSmiles(p[1]))
                    if product_smiles not in EXCLUSION_LIST:
                        product_list.append(product_smiles)
                        product_name_list.append(f"{first_name}_{second_name}")
                        product_ID_list.append(f"{first_ID}__SEPARATOR__{second_ID}")
                        if multireact_reactant_1:
                            product_matches_r1.append(r1_smarts_match_counts[i])
                        if multireact_reactant_2:
                            product_matches_r2.append(r2_smarts_match_counts[j])



            if secondary_reaction:
                secondary_product = secondary_reaction.RunReactants((r1_mol, r2_mol))
                if secondary_product:
                    reaction_count += 1
                    for s_p in secondary_product:
                        secondary_product_smiles = Chem.CanonSmiles(Chem.MolToSmiles(s_p[1]))
                        if secondary_product_smiles not in EXCLUSION_LIST:
                            product_list.append(secondary_product_smiles)
                            product_name_list.append(f"{first_name}_{second_name}")
                            product_ID_list.append(f"{first_ID}__SEPARATOR__{second_ID}")
                            if multireact_reactant_1:
                                product_matches_r1.append(r1_smarts_match_counts[i])
                            if multireact_reactant_2:
                                product_matches_r2.append(r2_smarts_match_counts[j])
    
    # Calculate upper and lower bounds for the total number of products
    upper_bound = 0
    lower_bound = 0
    r1_upper_bound = 0
    r2_upper_bound = 0
    r1_lower_bound = 0
    r2_lower_bound = 0
    if multireact_reactant_1:
        print(r1_smarts_match_counts)
        print(np.power(len(r2_smiles_list), r1_smarts_match_counts))
        r1_upper_bound = np.sum(np.power(len(r2_smiles_list), r1_smarts_match_counts))
        r1_lower_bound = np.sum(sp.special.comb(r1_smarts_match_counts.astype(np.int64), len(r2_smiles_list), exact=False, repetition=True))

    if multireact_reactant_2:
        r2_upper_bound = np.sum(np.power(len(r1_smiles_list), r2_smarts_match_counts))
        r2_lower_bound = np.sum(sp.special.comb(r2_smarts_match_counts.astype(np.int64), len(r1_smiles_list), exact=False, repetition=True))

    
    upper_bound = int(r1_upper_bound + r2_upper_bound)
    lower_bound = int(r1_lower_bound + r2_lower_bound)
    st.write(f"Upper bound on number of products: {upper_bound}. Lower bound on number of products: {lower_bound}.")
        
        

    # A list of products and matches is obtained

    total_products += len(product_list)
    product_dataframe = pd.DataFrame({"Product": product_list, 
                                      "Sample_ID": product_ID_list, 
                                      "Name": product_name_list})
    if multireact_reactant_1:
        product_dataframe["R1_matches"] = product_matches_r1
    if multireact_reactant_2:
        product_dataframe["R2_matches"] = product_matches_r2
    
    product_dataframe = product_dataframe.drop_duplicates(subset=["Product"])

    if multireact_reactant_1:
        product_dataframe_multi_r1 = product_dataframe[product_dataframe["R1_matches"] > 1]
        product_dataframe_multi_r1.dropna(inplace=True)
        product_dataframe_multi_r1_next = pd.DataFrame().reindex_like(product_dataframe_multi_r1)
        max_r1_matches = product_dataframe_multi_r1["R1_matches"].max()

        for idx in range(2, int(max_r1_matches + 1)):
            product_dataframe_multi_r1.dropna(inplace=True)
            product_dataframe_multi_r1_next.dropna(inplace=True)
    
            for i in product_dataframe_multi_r1.iterrows():
                # print(i[1])
                if (i[1]["Product"] is None or i[1]["Product"] == "" 
                    or i[1]["Sample_ID"] is None or i[1]["Sample_ID"] == "" 
                    or i[1]["Name"] is None or i[1]["Name"] == ""  
                    or i[1]["R1_matches"] is None or i[1]["R1_matches"] == ""):
                    continue
                r1_mol = Chem.MolFromSmiles(i[1]["Product"])
                first_name = i[1]["Name"]
                first_ID = i[1]["Sample_ID"]
                for j in range(len(r2_smiles_list)):
                    r2_mol  = Chem.MolFromSmiles(r2_smiles_list[j])
                    second_name = r2_names_list[j]
                    second_ID = r2_ID_list[j]

                    product = reaction.RunReactants((r1_mol, r2_mol))
                    if product:
                        reaction_count += 1
                        for p in product:
                            product_smiles = Chem.CanonSmiles(Chem.MolToSmiles(p[1]))
                            if product_smiles not in EXCLUSION_LIST:
                                product_list.append(product_smiles)
                                product_name_list.append(f"{first_name}_{second_name}")
                                product_ID_list.append(f"{first_ID}__SEPARATOR__{second_ID}")
                                if i[1]["R1_matches"] != idx:
                                    product_dataframe_multi_r1_next = pd.concat([product_dataframe_multi_r1_next, pd.DataFrame([{"Product": product_smiles, 
                                                                                                              "Sample_ID": f"{first_ID}__SEPARATOR__{second_ID}", 
                                                                                                              "Name": f"{first_name}_{second_name}",
                                                                                                              "R1_matches": i[1]["R1_matches"]}])], ignore_index=True)
                    if secondary_reaction:
                        secondary_product = secondary_reaction.RunReactants((r1_mol, r2_mol))
                        if secondary_product:
                            reaction_count += 1
                            for s_p in secondary_product:
                                secondary_product_smiles = Chem.CanonSmiles(Chem.MolToSmiles(s_p[1]))
                                if secondary_product_smiles not in EXCLUSION_LIST:
                                    product_list.append(secondary_product_smiles)
                                    product_name_list.append(f"{first_name}_{second_name}")
                                    product_ID_list.append(f"{first_ID}__SEPARATOR__{second_ID}")
                                    if i[1]["R1_matches"] != idx:
                                        product_dataframe_multi_r1_next = pd.concat([product_dataframe_multi_r1_next, pd.DataFrame([{"Product": product_smiles, 
                                                                                                                   "Sample_ID": f"{first_ID}__SEPARATOR__{second_ID}", 
                                                                                                                   "Name": f"{first_name}_{second_name}",
                                                                                                                   "R1_matches": i[1]["R1_matches"]}])], ignore_index=True)


            product_dataframe_multi_r1 = product_dataframe_multi_r1_next
            product_dataframe_multi_r1_next = pd.DataFrame().reindex_like(product_dataframe_multi_r1)
            
    
    if multireact_reactant_2:
        # Implement later
        pass
        # product_dataframe_R2_counts = product_dataframe.groupby("R2_matches")
    
    combined_list = zip(product_ID_list, product_name_list, product_list)
    combined_list = list(set(combined_list))
    product_ID_list_unique = [i[0] for i in combined_list]
    product_name_list_unique = [i[1] for i in combined_list]
    product_list_unique = [i[2] for i in combined_list]
    reaction_ran = True

    st.write(f"Reactions finished running. A total number of " 
             + f"{len(product_name_list)}"
             + f" products (may be non-unique) were obtained from {reaction_count} reactions ran."
             + f" A total number of {len(product_name_list_unique)} unique products were obtained.")
    
    end_time = time.time()
    total_time = end_time - start_time
    
    if total_time > 60:
        st.write(f"Total time taken: {total_time / 60} minutes.")
    elif total_time < 10e-3:
        st.write(f"Total time taken: {total_time * 10**3} milliseconds.")
    else:
        st.write(f"Total time taken: {end_time - start_time} seconds.")
    return product_ID_list_unique, product_name_list_unique, product_list_unique, reaction_ran
