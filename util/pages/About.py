import json
import streamlit as st
from streamlit_lottie import st_lottie
import requests

def about():
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

    def load_lottiefile(url: str):
        r = requests.get(url)
        if (r.status_code != 200):
            return None
        return r.json()

    lottie_hello = load_lottiefile("https://assets4.lottiefiles.com/packages/lf20_qhOsz8P0pl.json")
    st_lottie(lottie_hello,speed=1, reverse=False, loop=True, quality="high")

    st.title("SAMOSA")
    st.write("Hello")
