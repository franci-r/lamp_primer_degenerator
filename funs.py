# set functions 
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio import Align
from funs import *
import pandas as pd
from Bio import motifs

def primer_finder(sequence, primer):
    name = sequence.id
    sequence = str(sequence.seq).upper().replace("-", "")
    # align
    ali = primer_aligner(sequence, primer)
    # skip for bad alingment 
    if ali is None:
        #print(f"Skipping {name}: primer coverage below 70%")
        return None
    # get location
    #print(name)
    start = ali.aligned[0][0][0]
    end = ali.aligned[0][0][1]  
    # Include ignored terminal nucleotides
    start = start - ali.aligned[1][0][0]
    end = end + (len(primer) - ali.aligned[1][0][1]  )
    # select and output
    return sequence[start:end]

# align primer to sequence
def primer_aligner(sequence, primer):
    # set up aligner
    aligner = Align.PairwiseAligner()
    aligner.mode = "local"
    best = aligner.align(sequence, primer)[0]
    # remove bad alignments 
    counts = best.counts()
    error_rate = (counts.mismatches + counts.gaps) / len(primer)
    # Reject alignments with more than 30% errors
    # reject not full sequences
    aligned_length = best.aligned[0][-1][1] - best.aligned[0][0][0]
    if error_rate <= 0.30 and counts.gaps < 1 and aligned_length == len(primer):
        return best
    else:
        return None


# output df 
def create_output_table(motif_table, seq_number, threshold = 10):
    df = pd.DataFrame(dict(motif_table.counts)).T
    df_perc = df / seq_number * 100
    # create consensus 
    from Bio.Data.IUPACData import ambiguous_dna_values
    degenrate_dic = {frozenset(v): k for k, v in ambiguous_dna_values.items()}
    # reanme cols as new consensus 
    new_names = []
    for col in df_perc.columns:
        pos = "".join(df_perc.loc[df_perc[col] > threshold, col].index)
        new_names.append(degenrate_dic[frozenset(pos)])
    # add nt number
    df_perc.columns = range(1, len(df_perc.columns) + 1)
    # add consensus 
    df_perc.loc["Consensus"] = new_names
    # color 
    df_perc = df_perc.style.map(
        lambda x: "background-color: lightcoral"
        if not isinstance(x, str) and x > threshold
        else ""
    )
    return(df_perc)

# wrapper function
def primer_analysis(primer, primer_direction, threshold, alignment, primer_name):
    primer = primer.upper()
    extracted = []
    for record in alignment:
        name = record.id
        if primer_direction == 'F': 
            fragment = primer_finder(record, primer)
        elif primer_direction == 'R':
            fragment = primer_finder(record.reverse_complement(id=True, description=True), primer)
        else:
            raise ValueError("primer_direction must be 'F' or 'R'")
        # Skip rejected alignments
        if fragment is None:
            continue
        # output
        extracted.append(
                SeqRecord(
                    Seq(fragment),
                    id=name,
                    description=""
                )
            )
    # count % bases at each nt
    print("Number of extracted sequences:", len(extracted))
    motif = motifs.create(list(extracted))
    seq_number = len(extracted)
    # output number of successful finds
    print(f"{primer_name} successfully aligned {seq_number}/{len(alignment)} sequences")
    # output table
    return(create_output_table(motif, seq_number, threshold))

