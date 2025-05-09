import random
import pandas as pd

#Use table of toxic_seed_motifs
toxic_motifs = pd.read_csv("path_to_bad_seeds.csv")

#generates initial library of gRNAs
def generate_random_sequences(num, length=20):

    nucleotides = ['A', 'T', 'C', 'G']

    sequences = set()

    while len(sequences) < num:
        seq = ''.join(random.choices(nucleotides, k=length))
        sequences.add(seq)

    return list(sequences)


#checks if sequence contains a toxic motif
def toxic_motif_check(seq, toxic_motifs):
    for motif in toxic_motifs['seed']:
        if motif in seq:
            return True
    return False

#Determines GC content
def gc_content(seq):
    gcc = seq.count('G') + seq.count('C')
    return (gcc / len(seq)) * 100

# def disallowed_dinucleotides(seq):
#     return "TT" in seq or "AA" in seq


#Running of code
candidate_sequences = generate_random_sequences(100000)

filtered_sequences = []

for sequence in candidate_sequences:
    if toxic_motif_check(sequence, toxic_motifs):
        continue
    # if disallowed_dinucleotides(sequence):
    #     continue
    if 40 <= gc_content(sequence) <= 60:
        filtered_sequences.append(sequence)

print(len(filtered_sequences))

with open("initial_guides_100k.fasta", "w") as output:
    for i, seq in enumerate(filtered_sequences, start=1):
        output.write(f">gRNA_{i}\n{seq}\n")


