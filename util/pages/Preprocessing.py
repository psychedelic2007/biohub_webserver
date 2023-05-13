import streamlit as st
import zipfile
from Bio import SeqIO
from io import StringIO
from io import BytesIO
import random

def preprocessing():
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

    with open("quotes.txt", "r", encoding="utf-8") as f:
        file_text = f.read()
        quotes = file_text.split("\n\n")

    st.title("Welcome to the Preprocessing Page!")
    st.subheader("Upload your FASTA file")
    uploaded_file = st.file_uploader("Upload", type=["fasta"])
    col1, col2 = st.columns([1,10])
    with col1:
        show_submit = st.button("Submit")
    with col2:
        show_example = st.button("Load Example")

    def load_example_file():
        with open("example.fasta","rb") as f:
            return BytesIO(f.read())

    def preprocess_sequences(records):
        #remove sequences with "X"
        print("Removing the deuplicate sequences...")
        records = [r for r in records if "X" not in r.seq]

        #remove duplicate sequences
        sequences = []
        unique_records = []
        for record in records:
            sequence = str(record.seq)
            if sequence not in sequences:
                unique_records.append(record)
                sequences.append(sequence)
        st.write("Total number of sequences after removing duplicate sequences is:: ", len(unique_records))
        return unique_records
    
    def filter_sequences(sequences, target_length):
        print("Removing unequal length sequences...")
        filtered_sequences = {}
        for accession_id, sequence in sequences.items():
            if(len(sequence) == target_length):
                filtered_sequences[accession_id] = sequence
        st.write("Total number of sequences after making sequences of equal length is:: ", len(filtered_sequences))
        return filtered_sequences

    # Initialize the current index
    if "current_index" not in st.session_state:
        st.session_state.current_index = 0

    st.write("""***""")

    if show_submit:
        with st.spinner("SAMOSA is cleaning your dataset. Please wait...."):
            with zipfile.ZipFile("preprocessed_data.zip", "w") as output_zip:
                if(uploaded_file is not None):
                    stringio = StringIO(uploaded_file.getvalue().decode("utf-8"))
                    records = list(SeqIO.parse(stringio, "fasta"))
                    st.write("Total number of sequences before preprocessing is:: ",len(records))

                    #preprocess the sequences
                    processed_records = preprocess_sequences(records)
                    sequence_filtered = filter_sequences(processed_records)

                    #download the preprocessed file
                    if(len(sequence_filtered)>0):
                        with StringIO() as output:
                            SeqIO.write(sequence_filtered, output, "fasta")
                            processed_file = output.getvalue().encode()

                        st.download_button(label="Download Preprocessed File", data=processed_file, file_name="preprocessed.fasta", mime="application/octet-stream")
        quote = random.choice(quotes)
        st.write(quote)

    elif show_example:
        with st.spinner("SAMOSA is cleaning the example dataset. Please Wait....."):
            with zipfile.ZipFile("example_preprocessed_data.zip", "w") as output_zip:
                if(uploaded_file is None):
                    file_contents = load_example_file()
                    stringio = StringIO(file_contents.getvalue().decode("utf-8"))
                    records = list(SeqIO.parse(stringio, "fasta"))
                    st.write("Total number of sequences before preprocessing is:: ",len(records))

                    processed_records = preprocess_sequences(records)

                    if(len(processed_records)>0):
                        with StringIO() as output:
                            SeqIO.write(processed_records, output, "fasta")
                            processed_file = output.getvalue().encode()

                        st.download_button(label="Download Preprocessed Example File", data=processed_file, file_name="preprocessed_example_file.fasta", mime="application/octet-stream")

        st.write("")
        st.write("""***""")
        quote = random.choice(quotes)
        st.write(quote)
        st.write("""***""")


