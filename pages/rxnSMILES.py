import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem

import pandas as pd

from tqdm import tqdm
import requests
import json
import streamlit as st
import urllib.request as urlreq
import urllib.parse
from PIL import *
from stqdm import stqdm

from pathlib import Path
import sys

# from smiles_reaction_modules.helpers import *
# from smiles_reaction_modules.reaction_engine import *


sys.path.append(str(Path(__file__).resolve().parent.parent))

from modules.bulk_smiles_reactions.helpers import *
from modules.bulk_smiles_reactions.reaction_engine import *

SITE_ICON_PATH = "./site_files/bulk_smiles_icon.png"

REACTION_IMAGES_PATH = "./smarts_reaction_images/"

REACTION_DATABASE_PATH = "./databases/smarts_reaction_database.csv"

MASS_COLUMN_NAME = "Exact Mass with MASS_PRECISION Decimal Places"

batch_amine_secondary_NC = "[N:2]([#6:5])[#6:3].[C:4](=O)[O;D1:1]>>[O:2].[C:4](=O)[N:1]([#6:5])[#6:3]"

batch_amine_secondary_CN = "[C:4](=O)[O;D1:1].[N:2]([#6:5])[#6:3]>>[O:2].[C:4](=O)[N:1]([#6:5])[#6:3]"


st.sidebar.header("Bulk Reactor")


logo_col, title_col = st.columns([0.5, 4])

# Load and display logo in first column
with logo_col:
    logo = Image.open(SITE_ICON_PATH)
    st.image(logo, width=100)  # Adjust width as needed

# Display title in second column
with title_col:
    st.title('rxnSMILES')


SMARTS_RETRIEVAL_URL = "https://smarts.plus/smartsview/download_rest?smarts=INSERT_REACTION_SMARTS;filetype=png;vmode=0;textdesc=0;depsymbols=0;smartsheading=0"

SMARTS_RETRIEVAL_URL = "https://api.smarts.plus/smartsView/"
headers = {
    'Content-Type': 'application/json',
    # 'X-API-Key': 'unused'
}

data = {
    "query": {
        "smarts": "",
        "parameters": {
            "file_format": "png",
            "visualization_mode": 0,
            "legend_mode": 0,
            "visualization_of_default_bonds": 0,
            "labels_for_atoms": False,
            "smartstrim_active": False,
            "smarts_string_into_picture": True
        }
    }
}


st.write("Automatic SMILES generation for bulk reactions.")

input_reactants = st.file_uploader("Upload input reactant file", 
                                   type=["csv", 'xlsx', "xls", "tsv", "pkl"], 
                                   key='reactants_file')

smiles_list_length = None
if input_reactants is not None:
    input_reactants_name = input_reactants.name
    if input_reactants_name.endswith(".csv"):
        compound_list = pd.read_csv(input_reactants)
    elif input_reactants_name.endswith(".xlsx") or input_reactants_name.endswith(".xls"):
        compound_list = pd.read_excel(input_reactants)
    elif input_reactants_name.endswith(".tsv"):
        compound_list = pd.read_csv(input_reactants, sep="\t")
    elif input_reactants_name.endswith(".pkl"):
        compound_list = pd.read_pickle(input_reactants)
    else:
        input_reactants_name = pd.read_csv(input_reactants)
    st.write(compound_list)

    # Update to include batch file with parameters and use selectbox if batch not provided
    column_names = compound_list.columns.tolist()
    sample_ID_column_name = st.text_input("Enter Sample ID Column Name", value="unique_sample_id")
    if sample_ID_column_name not in column_names:
        st.write(f"Sample ID column name not found in input file. Please re-enter the correct column name.")
        st.stop()
    compound_name_column_name = st.text_input("Enter Compound Name Column Name", value="compound_name")
    if compound_name_column_name not in column_names:
        st.write(f"Compound name column name not found in input file. Please re-enter the correct column name.")
        st.stop()
    SMILES_column_name = st.text_input("Enter SMILES Column Name", value="SMILES")
    if SMILES_column_name not in column_names:
        st.write(f"SMILES column name not found in input file. Please re-enter the correct column name.")
        st.stop()


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
sample_ID_determinant = st.selectbox("Select which reactant to use as sample ID. Multireactions will use the singular reactant type as the sample ID.", options=["first reactant", "second reactant"])
mass_precision = st.number_input("Enter the number of decimal places for the mass precision. Enter a high value for maximal precision.", min_value=0, value=4)
default_mass_column_name = MASS_COLUMN_NAME.replace("MASS_PRECISION", str(mass_precision))
final_mass_column_name = st.text_input("Enter Mass Column Name", value=default_mass_column_name)
st.write("Select if either of the two sets of reactants should be allowed to undergo multiple reactions.")
# reactant_1_multireact = st.checkbox("Multireact Reactant 1", value=False)
# reactant_2_multireact = st.checkbox("Multireact Reactant 2", value=False)
reactant_1_multireact = False
reactant_2_multireact = False

reactions_database = pd.read_csv(REACTION_DATABASE_PATH)

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
            #reaction_smarts_url_safe = urllib.parse.quote(custom_reaction)
            #smarts_url = SMARTS_RETRIEVAL_URL.replace("INSERT_REACTION_SMARTS", reaction_smarts_url_safe)
            #urlreq.urlretrieve(smarts_url, REACTION_IMAGES_PATH + f"{custom_reaction}.png")
            data["query"]["smarts"] = custom_reaction
            response = requests.post(url, headers=headers, data=json.dumps(data))
            job_id = response.json().get("job_id")
            get_url = f"https://api.smarts.plus/smartsView/?job_id={job_id}"
            image_response = requests.get(get_url)

            if image_response.status_code == 200:
                with open(REACTION_IMAGES_PATH + f"{custom_reaction}.png", "wb") as f:
                    f.write(image_response.content)
            reaction_image = Image.open(REACTION_IMAGES_PATH + f"{custom_reaction}.png")
            reaction_image = trim(reaction_image)
            reaction_image.save(REACTION_IMAGES_PATH + f"{custom_reaction}.png")
            st.image(REACTION_IMAGES_PATH + f"{custom_reaction}.png")
        
        else:
            st.write("Custom reaction SMARTS is invalid. Please try again.")

else:
    try:
        st.image(REACTION_IMAGES_PATH + f"{reaction}.png")
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

reaction_started = st.button("Start Reactions")
reaction_ran = False

if (reaction_started):
    if reactant_number == 2:
        product_ID_list_unique, \
        product_name_list_unique, \
        product_list_unique, \
        reaction_ran = run_multi_double_reaction(reaction, 
                                                 number_of_reactant_1, 
                                                 smiles_list, 
                                                 names_list, 
                                                 sample_ID_list, 
                                                 reaction_name, 
                                                 reactant_1_multireact,
                                                 reactant_2_multireact)
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
    output_file_headers.append(final_mass_column_name)
    

    output_frame = pd.DataFrame(columns=output_file_headers)

    product_ID_list_unique = [i.split("__SEPARATOR__")[0] if sample_ID_determinant == "first reactant" else i.split("__SEPARATOR__")[1] for i in product_ID_list_unique]
    output_frame[sample_ID_column_name] = product_ID_list_unique
    output_frame[compound_name_column_name] = product_name_list_unique
    output_frame[SMILES_column_name] = product_list_unique
    output_frame[final_mass_column_name] = [round(rdkit.Chem.rdMolDescriptors.CalcExactMolWt(Chem.MolFromSmiles(i)), ndigits=mass_precision) for i in product_list_unique]
    
    for i in range(len(filler_column_names)):
        output_frame[filler_column_names[i]] = filler_column_values[i]

    st.write(output_frame)
    csv_output = convert_df(output_frame)
    st.download_button("Download Output File", csv_output, file_name=output_file_name)

    reaction_ran = False


