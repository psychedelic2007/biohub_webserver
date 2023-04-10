import json
import streamlit as st
from streamlit_lottie import st_lottie
import requests
import base64
from PIL import Image

def tutorial():

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

    def load_lottiefile(url: str):
        r = requests.get(url)
        if (r.status_code != 200):
            return None
        return r.json()

    lottie_hello = load_lottiefile("https://assets2.lottiefiles.com/private_files/lf30_hbeCPj.json")

    col1,col2 = st.columns([1,2])

    with col1:
        st_lottie(lottie_hello, speed=1, reverse=False, loop=True, quality="high")

    with col2:
        st.title("SAMOSA Tutorial")

    st.write("""***""")

    st.header("Overview")
    st.markdown("Welcome to our streamlit website, dedicated to streamlining the process of viral surveillance and prediction of future mutations."
                "Our website is designed to help biologists predict the epidemiological existence of a species and the growth trajectory that the "
                "organism can take. Our statistical-mechanics guided machine learning model, along with the Mutational Response Function, allows us "
                "to predict the phase transition of any virus from its ancestor to offspring.")

    st.markdown("At our website, users can upload their sequence data in fasta format, and our platform will handle the preprocessing without any "
                "hassle or manual labor. Our platform also provides feature calculation, enabling users to build their own multilayer perceptron "
                "architecture and train the model to predict future mutational residues. This allows users to create custom models that can predict "
                "the future mutations of a virus with great accuracy.")

    st.markdown("The Mutational Response Function is a parameter that we have extensively developed to predict the phase transition of any virus "
                "from its ancestor to offspring. By using this function, users can calculate the probability of a particular mutation occurring, "
                "and this allows them to predict the future trajectory of the virus.")

    st.markdown("Our website is designed to provide a streamlined flow to users, especially biologists, allowing them to predict the future of a "
                "virus with greater accuracy and ease. We are committed to providing the best possible tools to help our users make informed "
                "decisions about viral surveillance and management.")

    st.write("""***""")

    st.header("Tutorial")

    with st.expander("Step 1: Downloading the Dataset"):
        st.markdown("In this step, you will need to download a dataset of protein/genomic sequences either from the "
    "National Center for Biotechnology Information, Virus Database [NCBI Virus](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/sars-cov-2) "
                "or the Global Initiative on Sharing Avian Influenza Data [GISAID](https://gisaid.org/). NCBI and GISAID are two "
    "major databases for accessing genetic sequence data, and you can choose either one based on the specific research question "
    "you are trying to answer.")
        st.markdown("To download the dataset, you will need to navigate to the NCBI Virus or GISAID website and search for the dataset "
        "you want to download. You can filter the search results based on the specific criteria you are interested in, such as "
        "the type of virus or the geographic region. Once you have found the dataset you want to download, you can download it "
                    "in a compressed file format, such as .fasta or .fastq")
        st.markdown("In our previous [work](https://pubs.acs.org/doi/10.1021/acs.jpcb.2c04574) we downloaded the protein data "
                    "from NCBI Virus, in a month-wise manner, by selecting certain filters. Follow the following steps to "
                    "download the data in a month-wise manner ")

        st.markdown("**I) To download the protein data use NCBI Virus Database**")

        st.markdown("**_A_** : Go to [NCBI Virus](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/sars-cov-2) dashboard")
        file_ = open("data/step1_a.gif","rb")
        contents = file_.read()
        data_url = base64.b64encode(contents).decode("utf-8")
        file_.close()
        st.markdown(f'<img src="data:image/gif;base64,{data_url}" alt="step1_a gif">', unsafe_allow_html=True,)

        st.markdown("")
        st.markdown("")

        st.markdown("**_B_** : Click on the Country for which you want to download the data. In this case we select India "
                    "and then scroll up and select **_All Proteins_**")
        file_ = open("data/step1_b.gif", "rb")
        contents = file_.read()
        data_url = base64.b64encode(contents).decode("utf-8")
        file_.close()
        st.markdown(f'<img src="data:image/gif;base64,{data_url}" alt="step1_b gif">', unsafe_allow_html=True,)

        st.markdown("")
        st.markdown("")

        st.markdown("**_C_** : Here you will find the details of all the viral proteins in that country. From the left panel "
                    "select filters to narrow down the search to a specific variant. Click on **_Nucleotide Completeness_** " 
                    "and select **_Complete_**. This step makes sure that you are using a complete sequence and not a partial or "
                    "incomplete sequence")
        file_ = open("data/step1_c.gif", "rb")
        contents = file_.read()
        data_url = base64.b64encode(contents).decode("utf-8")
        file_.close()
        st.markdown(f'<img src="data:image/gif;base64,{data_url}" alt="step1_c gif">', unsafe_allow_html=True,)

        st.markdown("")
        st.markdown("")

        st.markdown("**_D_** : Once you have selected the complete sequence, click on **_Pango Lineage_** to select the "
                    "variant of interest. In this case we are selecting the **_Delta Variant (B.1.617.2)_**. Make sure "
                    "that you enter the variant ID and not the generic name of the variant.")
        file_ = open("data/step1_d.gif", "rb")
        contents = file_.read()
        data_url = base64.b64encode(contents).decode("utf-8")
        file_.close()
        st.markdown(f'<img src="data:image/gif;base64,{data_url}" alt="step1_d gif">', unsafe_allow_html=True, )

        st.markdown("")
        st.markdown("")

        st.markdown("**_E_** : Once you have selected the variant of your choice, then click on **_Proteins_**. Here you "
                    "can select which viral protein you would like to investigate. In this case we have selected the "
                    "**_Surface Glycoprotein_** of **_Delta Variant_**")
        file_ = open("data/step1_e.gif", "rb")
        contents = file_.read()
        data_url = base64.b64encode(contents).decode("utf-8")
        file_.close()
        st.markdown(f'<img src="data:image/gif;base64,{data_url}" alt="step1_e gif">', unsafe_allow_html=True, )

        st.markdown("")
        st.markdown("")

        st.markdown("**_F_** : Click on **_Collection Date_** to display the sequences in a month wise manner")
        file_ = open("data/step1_f.gif", "rb")
        contents = file_.read()
        data_url = base64.b64encode(contents).decode("utf-8")
        file_.close()
        st.markdown(f'<img src="data:image/gif;base64,{data_url}" alt="step1_f gif">', unsafe_allow_html=True,)

        st.markdown("")
        st.markdown("")

        st.markdown("**_G_** : Finally to download the sequences click on **_Download_** and proceed with the default options")
        file_ = open("data/step1_g.gif", "rb")
        contents = file_.read()
        data_url = base64.b64encode(contents).decode("utf-8")
        file_.close()
        st.markdown(f'<img src="data:image/gif;base64,{data_url}" alt="step1_g gif">', unsafe_allow_html=True, )

        st.markdown("")
        st.markdown("")
        st.write("***")

        st.markdown("**II) To download the genomic data use GISAID Database**")
        st.markdown("You need to make an account in GISAID first. Go to [GISAID](https://gisaid.org/) and in the top right corner "
                    "click on **_Log In_**. Create an account and once you log in you will see the home page like below")
        image = Image.open("data/gisiad.png")
        st.image(image)

        st.markdown("**_A_** : Make sure you are on the **_EpiCov_** tab. Click on the Viral structure of SARS-CoV-2 and "
                    "it will take you to the working dashboard")
        






    with st.expander("Step 2: Dataset Preprocessing"):
        st.write("")

    with st.expander("Step 3: Entropy Calculation"):
        st.write("")

    with st.expander("Step 4: Calculating Mutational Response Function"):
        st.write("")

    with st.expander("Step 5: Feature Calculation for Model Training"):
        st.write("")

    with st.expander("Step 6: Model Training & Testing"):
        st.write("")

    with st.expander("Step 7: Model Prediction"):
        st.write("")




