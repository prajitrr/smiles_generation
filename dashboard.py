import rdkit
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
from reaction_engine import *

MASS_COLUMN_NAME = "Mass"

batch_amine_secondary_NC = "[N:2]([#6:5])[#6:3].[C:4](=O)[O;D1:1]>>[O:2].[C:4](=O)[N:1]([#6:5])[#6:3]"

batch_amine_secondary_CN = "[C:4](=O)[O;D1:1].[N:2]([#6:5])[#6:3]>>[O:2].[C:4](=O)[N:1]([#6:5])[#6:3]"

st.set_page_config(
    page_title="autoSMILES",
    page_icon="./site_files/bulk_smiles_icon.png",
    layout="wide",
)

# site_icon = Image.open("./site_files/bulk_smiles_icon.png")
# st.markdown(
#     """
#     <div style="display: flex; align-items: center;">
#         <img src="data:image/png;base64,{st.image(./site_files/bulk_smiles_icon.png, use_column_width=False)}" alt="Logo" style="margin-right: 10px;">
#         <h1 style="margin: 0;">autoSMILES</h1>
#     </div>
#     """,
#     unsafe_allow_html=True
# )
logo_col, title_col = st.columns([0.3, 4])

# Load and display logo in first column
with logo_col:
    logo = Image.open('./site_files/bulk_smiles_icon.png')  # Replace with your logo path
    st.image(logo, width=100)  # Adjust width as needed

# Display title in second column
with title_col:
    st.title('autoSMILES')


SMARTS_RETRIEVAL_URL = "https://smarts.plus/smartsview/download_rest?smarts=INSERT_REACTION_SMARTS;filetype=png;vmode=0;textdesc=0;depsymbols=0;smartsheading=0"

st.write("Automatic SMILES generation for bulk reactions. Please contact prajkumar@ucsd.edu (Prajit Rajkumar) for any questions or bugs.")

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

    valid_smiles = compound_list[SMILES_column_name].apply(check_valid_smiles)
    if not valid_smiles.all():
        compound_list_indices = compound_list.index[~valid_smiles].tolist()
        compound_list_indices = [str(i) for i in compound_list_indices]
        st.write(f"Invalid SMILES detected at the following position(s): {', '.join(compound_list_indices)}.")
        st.write("Please fix the invalid SMILES string(s) and upload the file again. Note that the position number is based on the numbering in the table above.")
        st.stop()
    else:
        st.write("All SMILES are valid.")

number_of_reactant_1 = st.number_input("Enter number of reactant 1", min_value=1, max_value=smiles_list_length, value=1)

reactions_database = pd.read_csv("./smarts_reaction_database.csv")

reaction_name = st.selectbox("Select reaction to perform",
                                options=reactions_database["reaction_name"].tolist()
                                )

reaction = reactions_database.loc[reactions_database["reaction_name"] == reaction_name, "reaction_smarts"].values[0]
reactant_number = reactions_database.loc[reactions_database["reaction_name"] == reaction_name, "num_reactants"].values[0]

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
    try:
        st.image("./smarts_reaction_images/" + f"{reaction}.png")
    except:
        st.write("Bug accessing reaction image. Please contact prajkumar@ucsd.edu to report bug.")

#Output Column Filler Values
if "num_rows" not in st.session_state:
    st.session_state.num_rows = 2 

st.write(f"Enter Any Output File Filler Values:")

for i in range(st.session_state.num_rows):
    col1, col2 = st.columns(2)
    with col1:
        try:
            st.text_input(f"Filler Parameter {i+1} Name:", key=f"text_input_{i}_1")
        except:
            pass
    with col2:
        try:
            st.text_input(f"Filler Parameter {i+1} Value:", key=f"text_input_{i}_2")
        except:
            pass

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

product_list = []
product_name_list = []
product_ID_list = []

sample_ID_determinant = st.selectbox("Select which reactant to use as sample ID", options=["first reactant", "second reactant"])
reaction_started = st.button("Start Reactions")
reaction_ran = False

if (reaction_started):
    if reactant_number == 2:
        product_ID_list_unique, \
        product_name_list_unique, \
        product_list_unique, \
        reaction_ran = run_double_reaction(reaction, 
                                           number_of_reactant_1, 
                                           smiles_list, 
                                           names_list, 
                                           sample_ID_list, 
                                           reaction_name)
    elif reactant_number == 1:
        product_ID_list_unique, \
        product_name_list_unique, \
        product_list_unique, \
        reaction_ran = run_single_reaction(reaction, 
                                           smiles_list, 
                                           names_list, 
                                           sample_ID_list, 
                                           reaction_name)
if input_reactants is not None:
    output_file_name = st.text_input("Enter output file name", value=input_reactants_name.split(".")[0] + "_output." + input_reactants_name.split(".")[1])

if reaction_ran:

    output_file_headers = filler_column_names.copy()
    output_file_headers.insert(0, sample_ID_column_name)
    output_file_headers.append(compound_name_column_name)
    output_file_headers.append(SMILES_column_name)
    output_file_headers.append(MASS_COLUMN_NAME)
    

    output_frame = pd.DataFrame(columns=output_file_headers)

    product_ID_list_unique = [i.split("__SEPARATOR__")[0] if sample_ID_determinant == "first reactant" else i.split("__SEPARATOR__")[1] for i in product_ID_list_unique]
    output_frame[sample_ID_column_name] = product_ID_list_unique
    output_frame[compound_name_column_name] = product_name_list_unique
    output_frame[SMILES_column_name] = product_list_unique
    output_frame[MASS_COLUMN_NAME] = [rdkit.Chem.rdMolDescriptors.CalcExactMolWt(Chem.MolFromSmiles(i)) for i in product_list_unique]
    
    for i in range(len(filler_column_names)):
        output_frame[filler_column_names[i]] = filler_column_values[i]

    st.write(output_frame)
    csv_output = convert_df(output_frame)
    st.download_button("Download Output File", csv_output, file_name=output_file_name)

    reaction_ran = False


