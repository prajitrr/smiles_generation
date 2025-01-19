
from tqdm import tqdm
from stqdm import stqdm

import streamlit as st
from io import BytesIO
from PIL import Image
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
import xlsxwriter
import os

from pathlib import Path
import sys

from modules.bulk_smiles_reactions.helpers import *


SITE_ICON_PATH = "./site_files/bulk_smiles_icon.png"
TEMP_SVG_PATH = "./temp_storage/svgs/temporary_mol.svg"
TEMP_PNG_PATH = "./temp_storage/pngs/temporary_mol_NUMBER.png"
TEMP_PNG_FOLDER = "./temp_storage/pngs/"
TEMP_EXCEL_PATH = "./temp_storage/excel/temporary_excel.xlsx"

sys.path.append(str(Path(__file__).resolve().parent.parent))

logo_col, title_col = st.columns([0.5, 4])

# Load and display logo in first column
with logo_col:
    logo = Image.open(SITE_ICON_PATH)
    st.image(logo, width=100)  # Adjust width as needed

# Display title in second column
with title_col:
    st.title('imgSMILES')

st.sidebar.header("Structure Generation")



input_reactants = st.file_uploader("Upload input reactant file", 
                                   type=["csv", 'xlsx', "xls", "tsv", "pkl"], 
                                   key='reactants_file')
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
    smiles_list_length = len(smiles_list)

    valid_smiles = compound_list[SMILES_column_name].apply(check_valid_smiles)
    if not valid_smiles.all():
        compound_list_indices = compound_list.index[~valid_smiles].tolist()
        compound_list_indices = [str(i) for i in compound_list_indices]
        st.write(f"Invalid SMILES detected at the following position(s): {', '.join(compound_list_indices)}.")
        st.write("Please fix the invalid SMILES string(s) and upload the file again. Note that the position number is based on the numbering in the table above.")
        st.stop()
    else:
        num_columns = compound_list.shape[1]
        next_column_letter = chr(ord('A') + num_columns)
        st.write("All SMILES are valid.")



generation_started = st.button("Start Image Generation")
generated = False

if generation_started:

    output = BytesIO()

    with pd.ExcelWriter(output, engine='xlsxwriter') as writer:
        # # Write the DataFrame to Excel
        # compound_list.to_excel(writer, sheet_name='Sheet1', index=False)

        # # Get the xlsxwriter workbook and worksheet objects
        # workbook = writer.book
        # worksheet = writer.sheets['Sheet1']

        # writer.sheets['Sheet1'].write(0, num_columns + 1, "Image")


        output = BytesIO()
        writer = pd.ExcelWriter(output, engine='xlsxwriter')
        compound_list.to_excel(writer, index=False, sheet_name='Sheet1')
        workbook = writer.book
        worksheet = writer.sheets['Sheet1']
        worksheet.write(0, num_columns, "Image")


        # def to_excel(df):
        #     output = BytesIO()
        #     writer = pd.ExcelWriter(output, engine='xlsxwriter')
        #     df.to_excel(writer, index=False, sheet_name='Sheet1')
        #     workbook = writer.book
        #     worksheet = writer.sheets['Sheet1']
        #     format1 = workbook.add_format({'num_format': '0.00'}) 
        #     worksheet.set_column('A:A', None, format1)  
        #     writer.save()
        #     processed_data = output.getvalue()
        #     return processed_data
        # df_xlsx = to_excel(df)

        worksheet.set_column(num_columns, num_columns, 24)
        worksheet.set_default_row(113)

        for i, smiles in enumerate(stqdm(smiles_list)):
            current_path = TEMP_PNG_PATH.replace("NUMBER", str(i))
            mol = Chem.MolFromSmiles(smiles)

            img = Draw.MolToImage(mol, size=(150, 150))
            img.save(current_path, "PNG")



            # # Code for SVG Generation
            # drawer = rdMolDraw2D.MolDraw2DSVG(300,300)
            # drawer.DrawMolecule(mol)
            # drawer.FinishDrawing()
            # drawing = drawer.GetDrawingText()
            # with open(TEMP_SVG_PATH, "w") as f:
            #     f.write(drawing)
            worksheet.insert_image(next_column_letter + str(i + 2), current_path)

        writer.close()
        excel_data = output.getvalue()


    generated = True


    # Insert an image


    # excel_output = xlsxwriter.Workbook(TEMP_EXCEL_PATH)
    # excel_output_sheet = excel_output.add_worksheet()

if generated:
    for file in os.listdir(TEMP_PNG_FOLDER):
        if file.endswith('.png'):
            os.remove(os.path.join(TEMP_PNG_FOLDER, file))
    
    output_file_name = input_reactants_name.split(".")[0] + "_IMAGES.xlsx"
    st.download_button(label='Download Output File',
                                        data=excel_data,
                                        file_name= output_file_name)