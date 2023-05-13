from Bio import SeqIO
import csv

def compare_sequences(sequence1, sequence2):
    comparison = []
    for i in range(len(sequence1)):
        if sequence1[i] == sequence2[i]:
            comparison.append("0")
        else:
            comparison.append("1")
    return comparison

# Read FASTA file
def read_fasta_file(file_path):
    sequences = {}
    for record in SeqIO.parse(file_path, "fasta"):
        sequences[record.id] = str(record.seq)
    return sequences

# Write comparison results to CSV file
def write_csv_file(file_name, positions, comparison):
    with open(file_name, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(["Residue Position", "Target"])
        writer.writerows(zip(positions, comparison))

# Main program
fasta_file_path = "path/to/preprocessed/fasta/file"
sequences = read_fasta_file(fasta_file_path)

sequence_id1 = input("Enter the accession id of the first sequence: ")
sequence_id2 = input("Enter the accession id of the second sequence: ")

sequence1 = sequences.get(sequence_id1)
sequence2 = sequences.get(sequence_id2)

if sequence1 is not None and sequence2 is not None:
    comparison = compare_sequences(sequence1, sequence2)
    csv_file_name = f"{sequence_id1}.csv"
    write_csv_file(csv_file_name, range(1, len(sequence1) + 1), comparison)
    print(f"Comparison results saved to {csv_file_name}")
else:
    print("Invalid sequence accession ids.")
