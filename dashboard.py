from rdkit import Chem
from rdkit.Chem import AllChem
import pandas as pd
from tqdm import tqdm
import streamlit as st
import urllib.request as urlreq
import urllib.parse
from PIL import *
from stqdm import stqdm

from helpers import *

batch_amine_secondary_NC = "[N:2]([#6:5])[#6:3].[C:4](=O)[O;D1:1]>>[O:2].[C:4](=O)[N:1]([#6:5])[#6:3]"

batch_amine_secondary_CN = "[C:4](=O)[O;D1:1].[N:2]([#6:5])[#6:3]>>[O:2].[C:4](=O)[N:1]([#6:5])[#6:3]"

st.set_page_config(
    page_title="autoSMILES",
    page_icon="./site_files/bulk_smiles_icon.png",
    layout="wide",
)

site_icon = Image.open("./site_files/bulk_smiles_icon.png")
st.markdown(
    """
    <div style="display: flex; align-items: center;">
        <img src="data:image/png;base64,{st.image(site_icon, use_column_width=False)}" alt="Logo" style="margin-right: 10px;">
        <h1 style="margin: 0;">autoSMILES</h1>
    </div>
    """,
    unsafe_allow_html=True
)

SMARTS_RETRIEVAL_URL = "https://smarts.plus/smartsview/download_rest?smarts=INSERT_REACTION_SMARTS;filetype=png;vmode=0;textdesc=0;depsymbols=0;smartsheading=0"

input_reactants = st.file_uploader("Upload input reactant file", 
                                   type=["csv", 'xlsx', "xls", "tsv", "pkl"], 
                                   key='reactants_file')
sample_ID_column_name = st.text_input("Enter Sample ID Column Name", value="unique_sample_id")
compound_name_column_name = st.text_input("Enter Compound Name Column Name", value="compound_name")
SMILES_column_name = st.text_input("Enter SMILES Column Name", value="SMILES")

smiles_list_length = None
if input_reactants is not None:
    input_reactants_name = input_reactants.name
    if input_reactants_name.endswith(".csv"):
        compound_list = pd.read_csv(input_reactants)
    elif input_reactants_name.endswith(".xlsx") or input_reactants.endswith(".xls"):
        compound_list = pd.read_excel(input_reactants)
    elif input_reactants_name.endswith(".tsv"):
        compound_list = pd.read_csv(input_reactants)
    elif input_reactants_name.endswith(".pkl"):
        compound_list = pd.read_pickle(input_reactants)
    else:
        input_reactants_name = pd.read_csv(input_reactants)
    st.write(compound_list)

    smiles_list = compound_list[SMILES_column_name]
    names_list = compound_list[compound_name_column_name]
    sample_ID_list = compound_list[sample_ID_column_name]
    smiles_list_length = len(smiles_list)

number_of_reactant_1 = st.number_input("Enter number of reactant 1", min_value=1, max_value=smiles_list_length, value=1)

reaction_name = st.selectbox("Select reaction to perform",
                                options=["Amidation, Amine First",
                                "Esterification, Alcohol First",
                                "Amidation, Amine Second",
                                "Esterification, Alcohol Second",
                                "Hydroxyl Methylation",
                                "Amine Methylation",
                                "Carboxyl Methylation", 
                                "Custom Reaction"]
                                )

reaction_SMARTS = ["[N;D1:2][#6:3].[C:4](=O)[O;D1:1]>>[O:2].[C:4](=O)[N:1][#6:3]",
                   "[O;D1:2][#6;!$(C=O):3].[C:4](=O)[O;D1:1]>>[O:2].[C:4](=O)[O:1][#6:3]",
                   "[C:4](=O)[O;D1:1].[N;D1:2][#6:3]>>[O:2].[C:4](=O)[N:1][#6:3]",
                   "[C:4](=O)[O;D1:1].[O;D1:2][#6;!$(C=O):3]>>[O:2].[C:4](=O)[O:1][#6:3]",
                   "[O;D1:2][#6;!$(C=O):3]>>C[O;D2:2][#6;!$(C=O):3]",
                   "[N;D1:2][#6:3]>>C[N;D2:2][#6:3]",
                   "[C:4](=O)[O;D1:1]>>[C:4](=O)[O;D2:1]C",
                   "custom_reaction"
        ]

if reaction_name == "Amidation, Amine First":
    reaction = reaction_SMARTS[0]
elif reaction_name == "Esterification, Alcohol First":
    reaction = reaction_SMARTS[1]
elif reaction_name == "Amidation, Amine Second":
    reaction = reaction_SMARTS[2]
elif reaction_name == "Esterification, Alcohol Second":
    reaction = reaction_SMARTS[3]
elif reaction_name == "Hydroxyl Methylation":
    reaction = reaction_SMARTS[4]
elif reaction_name == "Amine Methylation":
    reaction = reaction_SMARTS[5]
elif reaction_name == "Carboxyl Methylation":
    reaction = reaction_SMARTS[6]
elif reaction_name == "Custom Reaction":
    reaction = reaction_SMARTS[7]

if reaction == "custom_reaction":
    custom_reaction = st.text_input("Enter custom reaction SMARTS")
    custom_reaction = custom_reaction.replace("\n", "")
    custom_reaction = custom_reaction.replace(" ", "")
    if custom_reaction != "":
        valid_reaction = False
        try:
            AllChem.ReactionFromSmarts(custom_reaction)
            valid_reaction = True
        except:
            valid_reaction = False

        if valid_reaction:
            st.write("Custom reaction SMARTS is valid.")
            reaction_smarts_url_safe = urllib.parse.quote(custom_reaction)
            smarts_url = SMARTS_RETRIEVAL_URL.replace("INSERT_REACTION_SMARTS", reaction_smarts_url_safe)
            urlreq.urlretrieve(smarts_url, "./smarts_reaction_images/" + f"{custom_reaction}.png")
            reaction_image = Image.open("./smarts_reaction_images/" + f"{custom_reaction}.png")
            reaction_image = trim(reaction_image)
            reaction_image.save("./smarts_reaction_images/" + f"{custom_reaction}.png")
            st.image("./smarts_reaction_images/" + f"{custom_reaction}.png")
        
        else:
            st.write("Custom reaction SMARTS is invalid. Please try again.")

else:
    st.image("./smarts_reaction_images/" + f"{reaction}.png")

#Output Column Filler Values
if "num_rows" not in st.session_state:
    st.session_state.num_rows = 2 

def add_row():
    st.session_state.num_rows += 1

def remove_row():
    if st.session_state.num_rows > 0:
        st.session_state.num_rows -= 1

st.write(f"Enter Any Output File Filler Values:")

for i in range(st.session_state.num_rows):
    col1, col2 = st.columns(2)
    with col1:
        st.text_input(f"Filler Parameter {i+1} Name:", key=f"text_input_{i}_1")
    with col2:
        st.text_input(f"Filler Parameter {i+1} Value:", key=f"text_input_{i}_2")

col_add, col_remove = st.columns(2)
with col_add:
    if st.button("Add Filler Parameter"):
        add_row()
with col_remove:
    if st.button("Remove Filler Parameter"):
        remove_row()

filler_column_names = []
filler_column_values = []
for i in range(st.session_state.num_rows):
    filler_column_names.append(st.session_state[f'text_input_{i}_1'])
    filler_column_values.append(st.session_state[f'text_input_{i}_2'])

filler_column_names = [x for x in filler_column_names if x != ""]
filler_column_values = [x for x in filler_column_values if x != ""]

filler_column_names = [x.strip() for x in filler_column_names]
filler_column_values = [x.strip() for x in filler_column_values]

exclusion_list = ["N", "O", "o", "OCCO"]

product_list = []
product_name_list = []
product_ID_list = []

sample_ID_determinant = st.selectbox("Select which reactant to use as sample ID", options=["first reactant", "second reactant"])
reaction_started = st.button("Start Reaction")
reaction_ran = False

if (reaction_started):
    reaction_SMARTS_string = reaction
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
                    if Chem.MolToSmiles(g[1]) not in exclusion_list:
                        try:
                            product_list.append(Chem.CanonSmiles(Chem.MolToSmiles(g[1])))
                            product_name_list.append(f"{first_name}_{second_name}")
                            product_ID_list.append(f"{first_ID}__SEPARATOR__{second_ID}")
                        except Exception as e:
                            print(e)
            elif reaction_name == "Amidation, Amine First" or reaction_name == "Amidation, Amine Second":
                if reaction_name == "Amidation, Amine First":
                    secondary_reaction = AllChem.ReactionFromSmarts(batch_amine_secondary_NC)
                else:
                    secondary_reaction = AllChem.ReactionFromSmarts(batch_amine_secondary_CN)
                product = secondary_reaction.RunReactants((r1_mol, r2_mol))
                if product:
                    for g in product:
                        if Chem.MolToSmiles(g[1]) not in exclusion_list:
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
    st.write(f"Reaction finished running. A total number of {len(product_name_list)} reactants was obtained from approximately {(smiles_list_length - number_of_reactant_1) * number_of_reactant_1} reactions ran.")
    combined_list = zip(product_ID_list, product_name_list, product_list)
    combined_list = list(set(combined_list))
    product_ID_list_unique = [i[0] for i in combined_list]
    product_name_list_unique = [i[1] for i in combined_list]
    product_list_unique = [i[2] for i in combined_list]
    reaction_ran = True

if reaction_ran:
    output_file_name = st.text_input("Enter output file name", value=input_reactants_name.split(".")[0] + "_output." + input_reactants_name.split(".")[1])

    output_file_headers = filler_column_names.copy()
    output_file_headers.insert(0, sample_ID_column_name)
    output_file_headers.append(SMILES_column_name)
    output_file_headers.append(compound_name_column_name)

    output_frame = pd.DataFrame(columns=output_file_headers)

    product_ID_list_unique = [i.split("__SEPARATOR__")[0] if sample_ID_determinant == "first reactant" else i.split("__SEPARATOR__")[1] for i in product_ID_list_unique]
    output_frame[sample_ID_column_name] = product_ID_list_unique
    output_frame[compound_name_column_name] = product_name_list_unique
    output_frame[SMILES_column_name] = product_list_unique
    for i in range(len(filler_column_names)):
        output_frame[filler_column_names[i]] = filler_column_values[i]

    st.write(output_frame)

    reaction_ran = False


