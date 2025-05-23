import pandas as pd

def parse_fasta(filename):

    with open(filename, "r") as f:
        header = None
        sequence_lines = []
        for line in f:
            line = line.strip()
            if not line:
                continue  
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(sequence_lines)
                header = line[1:].split()[0]  # Use the first word in the header as the ID
                sequence_lines = []
            else:
                sequence_lines.append(line)
        if header is not None:
            yield header, "".join(sequence_lines)

def write_fasta(records, output_filename):
    with open(output_filename, "w") as f:
        for header, sequence in records:
            f.write(f">{header}\n")
            f.write(f"{sequence}\n")

blast_cols = ["gRNA_name", "seq_id", "percent_ident", "length", "mismatch", "gap_open", "gRNA_start", "gRNA_end",
         "seq_start", "seq_end", "e_value", "bit_score"]

#Read BLAST results for genome and circuit parts
genome_df = pd.read_csv("data/genome_search_results", sep="\t", names=blast_cols)
circuit_df = pd.read_csv("data/circuit_search_results", sep="\t", names=blast_cols)

#computes alignment end position
genome_df["alignment_end"] = genome_df["gRNA_start"] + genome_df["length"] - 1
circuit_df["alignment_end"] = circuit_df["gRNA_start"] + circuit_df["length"] - 1

#Filter 1: Seed match (first 11 nucleotides match, with up to 2 nucleotide mismatch)
seed_hits_genome = genome_df[
    (genome_df["gRNA_start"] <= 11) & (genome_df["length"] > 9) & (genome_df["alignment_end"] <= 11) & (genome_df["mismatch"] <= 2)
]
# print(len(seed_hits_genome))
seed_hits_circuit = circuit_df[
    (circuit_df["gRNA_start"] <= 11) & (circuit_df["length"] > 9) & (circuit_df["alignment_end"] <= 11) & (circuit_df["mismatch"] <= 2)
]

#Filter 2: Signficant off-target matches - exclude matches that have any hit with alignment length between 17 and 20
# nucleotides where the mismatches are less then or equal to alignment length - 17.
sig_hits_genome = genome_df[
    (genome_df["length"] >= 17) & (genome_df["mismatch"] <= (genome_df["length"] - 17))]
sig_hits_circuit = circuit_df[
    (circuit_df["length"] >= 17) & (circuit_df["mismatch"] <= (circuit_df["length"] - 17))]

# Combined the guide names flagged by both filters.
seed_guides = set(seed_hits_genome["gRNA_name"].unique()).union(seed_hits_circuit["gRNA_name"].unique())
sig_guides = set(sig_hits_genome["gRNA_name"].unique()).union(sig_hits_circuit["gRNA_name"].unique())
overlapping_guides = seed_guides.union(sig_guides)

print("no. of guides excluded: ", len(overlapping_guides))

filtered_records = [
    (header, sequence)
    for header, sequence in parse_fasta("data/initial_guides.fasta")
    if header not in overlapping_guides
]
print("no. of guides after filtering:", len(filtered_records))

write_fasta(filtered_records, "data/filtered_guides_230325.fasta")
#
