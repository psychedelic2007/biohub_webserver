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
        file_ = open("data/gisaid_step1.gif", "rb")
        contents = file_.read()
        data_url = base64.b64encode(contents).decode("utf-8")
        file_.close()
        st.markdown(f'<img src="data:image/gif;base64,{data_url}" alt="gisaid_step1 gif">', unsafe_allow_html=True, )

        st.markdown("")
        st.markdown("")

        st.markdown("**_B_** : Once you are on the working dashboard navigate to **_Location_** tab. Here you can search "
                    "the genomic sequence based on the countries. But there is a format that you need to use to search. "
                    "First type the name of the **_Continent_** followed by a space and then a **_Forward Slash_**, followed "
                    "by the name of the country. For example if we need to search sequences from India then we will type "
                    "**_Asia / India_**")
        file_ = open("data/gisaid_step2.gif", "rb")
        contents = file_.read()
        data_url = base64.b64encode(contents).decode("utf-8")
        file_.close()
        st.markdown(f'<img src="data:image/gif;base64,{data_url}" alt="gisaid_step2 gif">', unsafe_allow_html=True, )

        st.markdown("")
        st.markdown("")

        st.markdown("**_C_** : After searching for the country of choice, navigate to the right side to select three important "
                    "features: **_Complete_**, **_High Coverage_** and **_Collection Date Complete_**. Selecting **_Complete_** "
                    "assures that the sequences that you will be getting are not incomplete or broken. Selecting **_High Coverage_** "
                    "affirms that the genomic sequences that will be displayed has less amount of unrecognised base pairs i.e. "
                    "low value of **_N's (<1% N and <0.05% mutations)_**. And finally selecting **_Collection Date Complete_** "
                    "allows the server to give you entries with complete collection date.")
        file_ = open("data/gisaid_step3.gif", "rb")
        contents = file_.read()
        data_url = base64.b64encode(contents).decode("utf-8")
        file_.close()
        st.markdown(f'<img src="data:image/gif;base64,{data_url}" alt="gisaid_step3 gif">', unsafe_allow_html=True, )

        st.markdown("")
        st.markdown("")

        st.markdown("**_D_** : Once the filters are selected then you can navigate to **_Lineage_** tab to select the "
                    "variant of your choice. For example if we want to download the genomic sequence of Delta variant then "
                    "in the **_Lineage_** tab we will type the ID of the Delta variant (i.e. B.1.617.2).")
        file_ = open("data/gisaid_step4.gif", "rb")
        contents = file_.read()
        data_url = base64.b64encode(contents).decode("utf-8")
        file_.close()
        st.markdown(f'<img src="data:image/gif;base64,{data_url}" alt="gisaid_step4 gif">', unsafe_allow_html=True, )

        st.markdown("")
        st.markdown("")

        st.markdown("**_E_** : Once you have selected the variant of interest then you can navigate to the **_Collection Date_** "
                    "and collect the sequences in a month wise manner")
        file_ = open("data/gisaid_step5.gif", "rb")
        contents = file_.read()
        data_url = base64.b64encode(contents).decode("utf-8")
        file_.close()
        st.markdown(f'<img src="data:image/gif;base64,{data_url}" alt="gisaid_step5 gif">', unsafe_allow_html=True, )

        st.markdown("")
        st.markdown("")

        st.markdown("**_F_** : After performing the above steps you will get the refined sequences in a month-wise manner. "
                    "To download the sequences you can click on the **_Square Box_** on the left side of **_Virus Name_** tab. "
                    "This will select all the sequences in the entry and then you can navigate to the bottom right and click on "
                    "the **_Download_** button to download the dataset. Once you click on **_Download_** a new window will pop up "
                    "asking you to select the format of the dataset. Click on **_Nucleotide Sequences (FASTA)_** to download the "
                    "dataset in FASTA format and then click on **_Download_**. Then you will be asked to accept the terms and "
                    "conditions of GISAID database. Click on **_I agree to the terms and conditions_** and then click **_Download_**.")
        file_ = open("data/gisaid_step6.gif", "rb")
        contents = file_.read()
        data_url = base64.b64encode(contents).decode("utf-8")
        file_.close()
        st.markdown(f'<img src="data:image/gif;base64,{data_url}" alt="gisaid_step6 gif">', unsafe_allow_html=True, )

        st.markdown("")
        st.markdown("")
        st.write("***")

    with st.expander("Step 2: Dataset Preprocessing"):
        st.markdown("**_A_** : After downloading the dataset it is important to preprocess it as the raw data consists of many junk and broken data."
                   "To preprocess the data click on **_Preprocessing_**tab from the main menu and then upload the sequence file that you have downloaded "
                   "in previous tep in **_FASTA_** format. Click on **_Submit_** after which you will get your preprocessed data. To know more about what exactly "
                   "the server does in preprocessing please visit the **_Theory Corner_**.")
        file_ = open("data/preprocessing_0.gif", "rb")
        contents = file_.read()
        data_url = base64.b64encode(contents).decode("utf-8")
        file_.close()
        st.markdown(f'<img src="data:image/gif;base64,{data_url}" alt="preprocessing_0 gif">', unsafe_allow_html=True, )

        st.markdown("")
        st.markdown("")
        
        st.markdown("**_B_**: If you are unsure about which data file or what type of datafile you need to upload, you click on **_Load Example_** and then "
                   "click on **_Submit_**. Once the preprocessing is done, you can download the file and then view it to understand the file type and its "
                   "architecture.")
        file_ = open("data/preprocessing_1.gif", "rb")
        contents = file_.read()
        data_url = base64.b64encode(contents).decode("utf-8")
        file_.close()
        st.markdown(f'<img src="data:image/gif;base64,{data_url}" alt="preprocessing_1 gif">', unsafe_allow_html=True, )
        st.write("***")

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




