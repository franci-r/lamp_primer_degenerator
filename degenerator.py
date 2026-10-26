from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio import Align
from funs import *
import pandas as pd
from Bio import motifs
import argparse

# protect from bad alignment 

parser = argparse.ArgumentParser()
parser.add_argument("--input", type=str, required=True)
parser.add_argument("--primers", type=str, required=True)
parser.add_argument("--threshold", type=float, default=15)
parser.add_argument("--output", type=str, default="out")
args = parser.parse_args()

# Read an alignment file (e.g., FASTA, Clustal, PHYLIP)
alignment = AlignIO.read(args.input, "clustal")
threshold = args.threshold

# process all primers
primers = list(SeqIO.parse(args.primers, "fasta"))

all_output = {}

for primer in primers: 
    if primer.id.upper()  in ["F3", "F2", "B1C", "LB"]:
        primer_direction = "F"
        out_df = primer_analysis(primer.seq, primer_direction, threshold , alignment, primer.id)
        all_output[primer.id] = out_df
    elif primer.id.upper()  in ["B3", "F1C", "B2", "LF"]:
        primer_direction = "R"
        out_df = primer_analysis(primer.seq, primer_direction, threshold , alignment, primer.id)
        all_output[primer.id] = out_df
    else: 
        print("Unrecognised primer name")
    

with pd.ExcelWriter(f"{args.output}.xlsx", engine="openpyxl") as writer:
    for sheet_name, df in all_output.items():
        df.to_excel(writer, sheet_name=sheet_name)


