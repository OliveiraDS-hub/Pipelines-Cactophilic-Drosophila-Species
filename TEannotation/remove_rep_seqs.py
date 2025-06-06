from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("fasta_file", help="Path to the input FASTA file")
args = parser.parse_args()

fasta_file = args.fasta_file

def filter_fasta(input_fasta, output_fasta):
    keywords_to_exclude = {"Satellite", "rRNA", "tRNA"}
    
    def keep_record(record):
        return not any(keyword in record.description for keyword in keywords_to_exclude)

    with open(output_fasta, "w") as out_f:
        SeqIO.write(
            (record for record in SeqIO.parse(input_fasta, "fasta") if keep_record(record)),
            out_f,
            "fasta"
        )


filter_fasta(fasta_file, "polished_TEs_s1.fa")