import streamlit as st
from PIL import Image
from util.Functions.gui import create_st_button
from util.pages.common import set_page_container_style
from streamlit_lottie import st_lottie
import requests
import json

st.set_page_config(layout="wide")

def home_page():

    st.markdown(
        """
        <style>
        [data-testid="stSidebar"][aria-expanded="true"]{
            background-color: #48a2cf;
        }
        </style>
        """,
        unsafe_allow_html=True,
    )
    def set_theme(theme):
        if theme == "light":
            st.set_theme("default")
        elif theme == "dark":
            st.set_theme("dark")

    # Theme button
    theme_button = st.button("Toggle Theme")

    # Check if the theme button is clicked
    if theme_button:
        # Check the current theme
        current_theme = st.config.get_option("theme.base")
        # Toggle the theme
        new_theme = "dark" if current_theme == "light" else "light"
        set_theme(new_theme)

    def load_lottiefile(url: str):
        r = requests.get(url)
        if (r.status_code != 200):
            return None
        return r.json()

    lottie_hello = load_lottiefile("https://assets5.lottiefiles.com/packages/lf20_ldf6z8do.json")

    left_col, right_col = st.columns(2)
    with left_col:
        st_lottie(lottie_hello, speed=1, reverse=False, loop=True, quality="high", width=600, height=400, key="initial")
    with right_col:
        right_col.markdown("<h1 style='text-align: center;'>SAMOSA</h1>", unsafe_allow_html=True)
        right_col.markdown("<h3 style='text-align: center;'>A tool to calculate Mutation Prediction</h3>",
                               unsafe_allow_html=True)
        right_col.markdown("<p style='text-align: center;'>Developed By Satyam Sangeet</p>", unsafe_allow_html=True)

    paper_link_dict = {
        "Mutation Prediction Paper": "https://pubs.acs.org/doi/10.1021/acs.jpcb.2c04574"
    }

    st.sidebar.markdown("## Related Paper")
    for link_text, link_url in paper_link_dict.items():
        create_st_button(link_text, link_url, st_col=st.sidebar)

    database_link_dict = {
        "GitHub Page": "https://github.com/psychedelic2007",
        "GISAID": "https://gisaid.org/",
        "NCBI Virus": "https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/sars-cov-2",
        "Influenza Databse": "https://www.ncbi.nlm.nih.gov/genomes/FLU/Database/nph-select.cgi?go=database"
    }

    st.sidebar.markdown("## Database Related Links")
    link_1_col, link_2_col = st.sidebar.columns(2)
    i=0
    link_col_dict = {0: link_1_col, 1: link_2_col}
    for link_text, link_url in database_link_dict.items():
        st_col = link_col_dict[i]
        i +=1
        if(i==len(link_col_dict.keys())):
            i=0
        create_st_button(link_text, link_url, st_col=st_col)

    software_link_dict = {
        "BioPython": "https://biopython.org",
        "Pandas": "https://pandas.pydata.org",
        "NumPy": "https://numpy.org",
        "SciPy": "https://scipy.org",
        "Sklearn": "https://scikit-learn.org/stable/",
        "Matplotlib": "https://matplotlib.org",
        "Seaborn": "https://seaborn.pydata.org",
        "Streamlit": "https://streamlit.io",
        "MAFFT": "https://www.ebi.ac.uk/Tools/msa/mafft/",
        "NetworkX": "https://networkx.org/",
        "Plotly": "https://plotly.com/"
    }

    st.sidebar.markdown("## Software Related Links")
    link_1_col, link_2_col, link_3_col = st.sidebar.columns(3)
    i=0
    link_col_dict = {0: link_1_col, 1: link_2_col, 2: link_3_col}
    for link_text, link_url in software_link_dict.items():
        st_col = link_col_dict[i]
        i +=1
        if(i==len(link_col_dict.keys())):
            i=0
        create_st_button(link_text, link_url, st_col=st_col)
    st.markdown("---")

    #st.write("# Welcome to BioHub! 👋")

    st.markdown(
        """
		### Overview
		
		*SAMOSA* is a website that offers all the necessary calculation required for a protein/genomic sequences

		**👈 Select an analysis from the sidebar** to start your journey of Biological Analysis!
			""")
    # About the Authors
    st.markdown("""
		### Authors
		Please feel free to contact me if you have any issue, comments or questions!
		### Satyam Sangeet
		- Email : satyam85cool@gmail.com
		- [GitHub](https://github.com/psychedelic2007)
		- [WebPage](https://sites.google.com/view/satyamsangeet) 
			""")

    st.write("""***""")

    # Development Field
    st.markdown("""
		## Developed and Maintained by Satyam Sangeet & Dr. Susmita Roy
		[Dr. Roy lab](https://www.drsusmitaroy.com/)
		Copyright (c) 2022 Satyam Sangeet
			""")
