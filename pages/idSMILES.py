from tqdm import tqdm
from stqdm import stqdm

import streamlit as st
from io import BytesIO
from PIL import Image
import pandas as pd
from rdkit import Chem
import rdkit
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import Descriptors
import xlsxwriter
import os
import numpy as np
from modules.bulk_smiles_reactions.helpers import *

from pathlib import Path
import sys

SITE_ICON_PATH = "./site_files/bulk_smiles_icon.png"
TEMP_SVG_PATH = "./temp_storage/svgs/temporary_mol.svg"
TEMP_PNG_PATH = "./temp_storage/pngs/temporary_mol_NUMBER.png"
TEMP_PNG_FOLDER = "./temp_storage/pngs/"
TEMP_EXCEL_PATH = "./temp_storage/excel/temporary_excel.xlsx"
DATABASE_PATH = "./databases/SMILES_names_master_database.csv"

sys.path.append(str(Path(__file__).resolve().parent.parent))

logo_col, title_col = st.columns([0.5, 4])

# Load and display logo in first column
with logo_col:
    logo = Image.open(SITE_ICON_PATH)
    st.image(logo, width=100)  # Adjust width as needed

# Display title in second column
with title_col:
    st.title('idSMILES')
st.sidebar.header("Compound identifiers & properties")

st.write("Convert between compound names and SMILES and other identifiers as well as retrieve useful molecular properties. Emphasis on standard names of small molecules as well as naming conventions of the Dorrestein Lab.")
#st.write("If you would like to search for identifiers/properties for just individual compounds, please navigate to idSMILES+ on the left sidebar.")
#st.write("Users can also add their own custom names and SMILES strings to the database for future use under the site tab idSMILES+ (located on the left sidebar). Please do not abuse this feature 🥺.")
st.write("Contact Prajit Rajkumar (prajkumar@ucsd.edu) if you have a large number of names you would like to add.")

smiles_database = pd.read_csv(DATABASE_PATH)

input_names = st.file_uploader("Upload input names file", 
                                   type=["csv", 'xlsx', "xls", "tsv", "pkl"], 
                                   key='names_file')
identifier = st.selectbox("Select Compound Identifier", ["Name", "InChI", "SMILES"])
identifier_column_name = st.text_input("Enter Identifier Column Name", value="Names")
st.write('Select which values you would like to output:')
option_1 = st.checkbox('Name')
option_2 = st.checkbox('SMILES')
option_3 = st.checkbox('InChI')
option_4 = st.checkbox('InChIKey')
option_5 = st.checkbox('Molecular Formula')
option_6 = st.checkbox('Molecular Weight')

retrieval_started = False
if input_names is not None:
    input_names_name = input_names.name
    if input_names_name.endswith(".csv"):
        compound_list = pd.read_csv(input_names)
    elif input_names_name.endswith(".xlsx") or input_names_name.endswith(".xls"):
        compound_list = pd.read_excel(input_names)
    elif input_names_name.endswith(".tsv"):
        compound_list = pd.read_csv(input_names, sep="\t")
    elif input_names_name.endswith(".pkl"):
        compound_list = pd.read_pickle(input_names)
    else:
        compound_list = pd.read_csv(input_names)
    st.write(compound_list)

    if identifier_column_name not in compound_list.columns:
        st.error("Names column not found in input file. Please check the column name and try again.")
        st.stop()

    if identifier == "Name":
        smiles_dict = pd.Series(smiles_database.SMILES.values,index=smiles_database.compound_name).to_dict()
    else:
        if identifier == "SMILES":
            valid_smiles = compound_list[identifier_column_name].apply(check_valid_smiles)
            if not valid_smiles.all():
                compound_list_indices = compound_list.index[~valid_smiles].tolist()
                compound_list_indices = [str(i) for i in compound_list_indices]
                st.write(f"Invalid SMILES detected at the following position(s): {', '.join(compound_list_indices)}.")
                st.write("Please fix the invalid SMILES string(s) and upload the file again. Note that the position number is based on the numbering in the table above.")
                st.stop()

        elif identifier == "InChI":
            valid_inchi = compound_list[identifier_column_name].apply(check_valid_inchi)
            if not valid_inchi.all():
                compound_list_indices = compound_list.index[~valid_inchi].tolist()
                compound_list_indices = [str(i) for i in compound_list_indices]
                st.write(f"Invalid InChI detected at the following position(s): {', '.join(compound_list_indices)}.")
                st.write("Please fix the invalid InChI string(s) and upload the file again. Note that the position number is based on the numbering in the table above.")
        
        else:
            st.error("This is a bug and should not be possible. Please contact the developer. 💀")
            st.stop()



        smiles_dict = pd.Series(smiles_database.compound_name.values,index=smiles_database.SMILES).to_dict()

    
    identifier_list = compound_list[identifier_column_name]
    identifier_list_length = len(identifier_list)
    

    

    retrieval_started = st.button("Start Retrieval")
    retrieved = False
    




if retrieval_started:
    st.write("Retrieving identifiers and properties from database...")

    if option_1 and identifier != "Name":
        names_list = []
    if option_2 and identifier != "SMILES":
        smiles_list = []
    if option_3 and identifier != "InChI":
        inchi_list = []
    if option_4:
        inchikey_list = []
    if option_5:
        formula_list = []
    if option_6:
        weight_list = []

    for item in identifier_list:
        if identifier == "InChI":
            molecule = Chem.MolFromInchi(item)

        elif identifier == "SMILES":
            molecule = Chem.MolFromSmiles(item)
        
        
        elif identifier == "Name":
            if item in smiles_dict:
                molecule = Chem.MolFromSmiles(smiles_dict[item])
            else:
                molecule = None
        if molecule is None:
            if option_1 and identifier != "Name":
                names_list.append(np.nan)
            if option_2 and identifier != "SMILES":
                smiles_list.append(np.nan)
            if option_3 and identifier != "InChI":
                inchi_list.append(np.nan)
            if option_4:
                inchikey_list.append(np.nan)
            if option_5:
                formula_list.append(np.nan)
            if option_6:
                weight_list.append(np.nan)
        
        else:
            if option_1 and identifier != "Name":
                if item in smiles_dict:
                    names_list.append(smiles_dict[identifier])
                else:
                    names_list.append(np.nan)
            if option_2 and identifier != "SMILES":
                smiles_list.append(Chem.MolToSmiles(molecule))
            if option_3 and identifier != "InChI":
                inchi_list.append(Chem.MolToInchi(molecule))
            if option_4:
                inchikey_list.append(Chem.MolToInchiKey(molecule))
            if option_5:
                formula_list.append(Chem.rdMolDescriptors.CalcMolFormula(molecule))
            if option_6:
                weight_list.append(rdkit.Chem.rdMolDescriptors.CalcExactMolWt(molecule))

    if option_1 and identifier != "Name":
        compound_list["Name"] = names_list
    if option_2 and identifier != "SMILES":
        compound_list["SMILES"] = smiles_list
    if option_3 and identifier != "InChI":
        compound_list["InChI"] = inchi_list
    if option_4:
        compound_list["InChIKey"] = inchikey_list
    if option_5:
        compound_list["Molecular Formula"] = formula_list
    if option_6:
        compound_list["Molecular Weight"] = weight_list
    
    st.write(compound_list)
    retrieved = True

