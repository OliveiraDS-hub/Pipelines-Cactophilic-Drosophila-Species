import subprocess
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("fasta_file", help="Path to the input FASTA file")
args = parser.parse_args()

fasta_f = args.fasta_file

result = subprocess.run(
    'trf {} 2 5 6 75 20 50 500 -m -h'.format(fasta_f),
    shell=True,
    stdout=subprocess.PIPE,
    stderr=subprocess.PIPE,
    text=True)