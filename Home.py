from tqdm import tqdm
import streamlit as st
import urllib.request as urlreq
import urllib.parse
from PIL import Image
from stqdm import stqdm

SITE_ICON_PATH = "./site_files/bulk_smiles_icon.png"


#st.sidebar.success("Select a workflow above.")

logo_col, title_col = st.columns([0.5, 4])

# Load and display logo in first column
with logo_col:
    logo = Image.open(SITE_ICON_PATH)
    st.image(logo, width=100)  # Adjust width as needed

# Display title in second column
with title_col:
    st.title('autoSMILES')

st.markdown(
    """
    autoSMILES is a comprehensive dashboard for bulk computational chemistry. Contact Prajit Rajkumar (prajkumar@ucsd.edu)
    for further questions or bug fixes, or if you have suggestions for a new workflow to implement. Made as part of the 
    [Dorrestein Lab](https://dorresteinlab.ucsd.edu/) at the UCSD Skaggs School of Pharmacy and Pharmaceutical Sciences.
"""
)