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
parser.add_argument("--primer", type=str, required=True)
parser.add_argument("--primer_dir", type=str, required=True)
parser.add_argument("--threshold", type=float, default=15)
args = parser.parse_args()

# Read an alignment file (e.g., FASTA, Clustal, PHYLIP)
alignment = AlignIO.read(args.input, "clustal")
threshold = args.threshold
primer = args.primer 
primer_direction = args.primer_dir 

out_df = primer_analysis(primer, primer_direction, threshold , alignment)
# color cells 
out_df = out_df.style.map(
    lambda x: "background-color: lightcoral"
    if not isinstance(x, str) and x > threshold
    else ""
)
out_df.to_excel("out.xlsx")

