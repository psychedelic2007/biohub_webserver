import streamlit as st

def theory():

    def local_css(file_name):
        with open(file_name) as f:
            st.markdown(f"<style>{f.read()}</style>", unsafe_allow_html=True)

    local_css("util/pages/style.css")

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
    
    st.title("Theory Corner")
    st.markdown("Welcome to the Theory Corner, where all your burning questions about our web servers are answered! From the basics of server architecture to the "
                "mathematical intricacies of the services that our webserver provides, this page is your one-stop-shop for all things theoretical in the world of "
                "viral bioinformatics. Whether you're a seasoned bioinformatician or just starting out, our webserver and tutorials will help you provide "
                "a platform that will allow you to seamlessly use the mathematical and machine learning techniques for biological exploration. So grab a cup of coffee, "
                "settle in, and get ready to level up your knowledge with the Theory Corner!")
    
    st.write("")
    st.write("***")
