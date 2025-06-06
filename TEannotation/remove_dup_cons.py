import os.path
import argparse
import pandas as pd
import subprocess
from io import StringIO
from Bio import SeqIO
import re
from collections import defaultdict
from collections import Counter
import sys
from Bio import SeqIO

### Take sequences if the same ID family as a dictionary
# Define input FASTA file

# Parse positional argument
parser = argparse.ArgumentParser()
parser.add_argument("fasta_file", help="Path to the input FASTA file")
parser.add_argument("output_file", help="Path to the input FASTA file")
args = parser.parse_args()

fasta_file = args.fasta_file
output_file = args.output_file

blast_result = subprocess.run("blastn -subject {} "
        "-query {} "
        "-outfmt 6 "
        "-qcov_hsp_perc 80 "
        "-perc_identity 80".format(fasta_file, fasta_file),
        shell=True,
        capture_output=True,
        text=True)

seq_id = None
# Check if the command was successful
if blast_result.returncode != 0:
    print("Error running blastn:", blast_result.stderr)
else:
    # Read the BLAST output into a Pandas DataFrame
    blast_columns = ["query", "subject", "pident", "length", "mismatch", "gapopen",
        "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
    blast_out_df = pd.read_csv(StringIO(blast_result.stdout), sep="\t", names=blast_columns)

    ## Remove self matches´
    blast_out_df = blast_out_df[blast_out_df["query"] != blast_out_df["subject"]]

    family_IDs = blast_out_df.iloc[:, 0].drop_duplicates().tolist()

    ## Create .fai from consensus
    subprocess.run('samtools faidx {}'.format(fasta_file), shell = True)

    ## Read the .fai file to get full sequence lengths
    fai_df = pd.read_csv(f"{fasta_file}.fai", sep="\t", header=None, usecols=[0,1],
                            names=["sequence", "seq_length"])

    ## Merge the BLAST results with the .fai file to get full sequence lengths
    blast_out_df_w_length = pd.merge(blast_out_df, fai_df, left_on="query", right_on="sequence")

if len(family_IDs) > 0:
    removal_list = []
    longest_seq_list = []
    for family in family_IDs:
        if not family in removal_list:
            family_matches = blast_out_df_w_length[blast_out_df_w_length["query"] == family]
            unique_subjects = family_matches["subject"].drop_duplicates().tolist()

            unique_subjects.append(family)

            ### Filter only rows where 'subject' matches any of the unique_subjects
            subset_matches = family_matches[family_matches["subject"].isin(unique_subjects)]

            ### Find the entry with the highest seq_length
            longest_subject = subset_matches.loc[subset_matches['seq_length'].idxmax(), 'subject']
            longest_seq_list.append(longest_subject)

            removal_list.extend(unique_subjects)

unique_longest = list(set(longest_seq_list))
# print(f"Longest sequences:\n{len(unique_longest)}")

unique_removal = list(set(removal_list))
# print(f"Sequences to remove:\n{len(unique_removal)}")


list2remove = [item for item in unique_removal if item not in unique_longest]
print(f"Total duplicated consensus to be removed: {len(list2remove)}")



def remove_fasta_headers(input_fasta, output_fasta, headers_to_remove):
    headers_to_remove_set = set(headers_to_remove)  # for fast lookup

    with open(output_fasta, "w") as out_f:
        SeqIO.write(
            (record for record in SeqIO.parse(input_fasta, "fasta") if record.id not in headers_to_remove_set),
            out_f,
            "fasta")
        
remove_fasta_headers(fasta_file, output_file, list2remove)