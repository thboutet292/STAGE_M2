from Bio import SeqIO
import pandas as pd
import sys
from Bio.Seq import Seq

# Vérification des arguments
if len(sys.argv) != 3:
    print("Usage: script_extract_orthologues.py <nom_microsporidie> <i_value>")
    sys.exit(1)

SPECIES = sys.argv[1]
I_VALUE = sys.argv[2]

# Fichiers d'entrée
df_path = f"orthofinder/orthogroup_matches/{SPECIES}_I_{I_VALUE}/resultats_hmmsearch.tsv"
fa_sens = f"genome_complet/genome_prot/{SPECIES}_I_{I_VALUE}_frag.fna"
fa_antisens = f"genome_complet/genome_prot/{SPECIES}_I_{I_VALUE}_frag_rev.fna"
output_fasta = f"data_training/PROT/PROT_{SPECIES}_I_{I_VALUE}"

df = pd.read_csv(df_path, sep="\t")
genome_sens = SeqIO.to_dict(SeqIO.parse(fa_sens, "fasta"))
genome_antisens = SeqIO.to_dict(SeqIO.parse(fa_antisens, "fasta"))

compt_gene = 0
unique_sequences = {}

for index, row in df.iterrows():
    compt_gene += 1
    contig = row["Gène"].strip()
    start, end = int(row["Début"]), int(row["Fin"])
    strand = row["Strand"].strip()

    # Récupération de la séquence protéique (AA)
    if strand == "sens" and contig in genome_sens:
        prot_full = str(genome_sens[contig].seq)
    elif strand == "antisens" and contig in genome_antisens:
        prot_full = str(genome_antisens[contig].seq)
    else:
        continue  # Ignorer si le contig est introuvable

    # Prendre la séquence entre start et end, puis prolonger jusqu'au premier X
    seq_str = prot_full[start-1:end]
    if "X" not in seq_str:
        rest = prot_full[end:]
        for aa in rest:
            seq_str += aa
            if aa == "X":
                break

    # Arrêter à la première occurrence de X (au cas où il y en ait un après le prolongement)
    if "X" in seq_str:
        seq_str = seq_str.split("X", 1)[0]

    # Vérifier que la séquence commence par 'M' et est assez longue
    if seq_str and seq_str[0] == 'M' and len(seq_str) > 50:
        if seq_str not in unique_sequences:
            unique_sequences[seq_str] = f"gene{compt_gene}"
        else:
            if len(seq_str) > len(unique_sequences[seq_str]):
                unique_sequences[seq_str] = f"gene{compt_gene}"

extracted_seqs = [SeqIO.SeqRecord(Seq(seq), id=gene_id, description="") for seq, gene_id in unique_sequences.items()]
if extracted_seqs:
    SeqIO.write(extracted_seqs, output_fasta, "fasta")

