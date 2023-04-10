import streamlit as st

def process_fasta_file(file):
    sequences = []
    current_sequence = ""
    with open(file) as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith(">"):
                if current_sequence:
                    sequences.append(current_sequence)
                    current_sequence = ""
            else:
                current_sequence += line.strip()
        sequences.append(current_sequence)
    return sequences

def perform_calculation(sequence):
    # Perform your calculation on the sequence here
    result = len(sequence)
    return result

def main():
    st.title("Fasta Sequence Processing")
    uploaded_file = st.file_uploader("Upload your fasta file", type=["fasta"])
    if uploaded_file is not None:
        sequences = process_fasta_file(uploaded_file)
        results = [perform_calculation(sequence) for sequence in sequences]
        st.write("Calculation Results:", results)

if __name__ == "__main__":
    main()
