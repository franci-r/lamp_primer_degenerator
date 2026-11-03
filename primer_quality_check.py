from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio import Align
from funs import *
import pandas as pd
from Bio import motifs
import argparse
import subprocess
from Bio.Data.IUPACData import ambiguous_dna_values
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator


parser = argparse.ArgumentParser()
parser.add_argument("--input", type=str, required=True)
parser.add_argument("--primers", type=str, required=True)
parser.add_argument("--output", type=str, default="out")
args = parser.parse_args()

primers = list(SeqIO.parse(args.primers, "fasta"))

all_output = []

for primer in primers: 
    sequences = list(SeqIO.parse(args.input, "fasta"))
    sequences_rev_comp = [record.reverse_complement(id=True, description=True) for record in sequences]
    if primer.id.upper()  in ["F3", "F2", "B1C", "LB"]:
        out = primer_mismatch_counter(sequences, primer)
        all_output.append(out)
    elif primer.id.upper()  in ["B3", "F1C", "B2", "LF"]:
        out = primer_mismatch_counter(sequences_rev_comp, primer)
        all_output.append(out)
    else: 
        print("Unrecognised primer name")
    

df = pd.concat(all_output, ignore_index=True)
# Faceted boxplots by primer

sns.histplot(
    data=df,
    x="mismatches",
    hue="primer",
    multiple="dodge",
    discrete=True,
    shrink=0.8
)

plt.xlabel("Number of mismatches")
plt.ylabel("Number of sequences")
plt.savefig(f"{args.output}_mismatch.png", dpi=300, bbox_inches="tight")
plt.close()

# Plot 2: Percentages
sns.histplot(
    data=df,
    x="mismatches",
    hue="primer",
    multiple="dodge",
    discrete=True,
    shrink=0.8,
    stat="percent",
    common_norm=False
)

plt.xlabel("Number of mismatches")
plt.ylabel("Percentage of sequences (%)")

# Set y-axis from 0 to 100 with ticks every 10
plt.ylim(0, 100)
plt.yticks(range(0, 101, 10))
plt.gca().xaxis.set_major_locator(MultipleLocator(1))


plt.savefig(f"{args.output}_mismatch_percent.png", dpi=300, bbox_inches="tight")
plt.close()



