from Bio import SeqIO
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument("fasta_file", help="Path to the input FASTA file")
parser.add_argument("output", help="Path to the input FASTA file")
parser.add_argument("list_id", help="Path to the input FASTA file")
args = parser.parse_args()

fasta_file = args.fasta_file
output = args.output
list_id = args.list_id

def remove_fasta_headers(input_fasta, output_fasta, headers_to_remove):
    headers_to_remove_set = set(headers_to_remove)
    with open(output_fasta, "w") as out_f:
        SeqIO.write(
            (record for record in SeqIO.parse(input_fasta, "fasta") if record.id not in headers_to_remove_set),
            out_f,
            "fasta")

if os.path.exists(list_id):
    with open(list_id) as f:
        repeat_list = [line.strip() for line in f if line.strip()] 
    remove_fasta_headers(fasta_file, output, repeat_list)
else:
    remove_fasta_headers(fasta_file, output, list_id)