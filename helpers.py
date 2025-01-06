from PIL import Image, ImageChops
import streamlit as st
from rdkit import Chem

def trim(im):
    bg = Image.new(im.mode, im.size, im.getpixel((0,0)))
    diff = ImageChops.difference(im, bg)
    diff = ImageChops.add(diff, diff, 2.0, -100)
    bbox = diff.getbbox()
    if bbox:
        return im.crop(bbox)

def add_row():
    st.session_state.num_rows += 1
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

def remove_row():
    if st.session_state.num_rows > 0:
        st.session_state.num_rows -= 1
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

def check_valid_smiles(smiles_string):
    try:
        mol = Chem.MolFromSmiles(smiles_string)
    except:
        return False

    if mol is None:
        return False
    else:
        return True
