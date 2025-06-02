#!/usr/bin/env python3

import os
import sys
from collections import defaultdict
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO import write

# Vérifier les arguments
if len(sys.argv) != 3:
	print("Usage: script_gene_core.py <nom_microsporidie> <i_value>")
	sys.exit(1)

# Paramètres d'entrées
SPECIES = sys.argv[1]
I_VALUE = sys.argv[2]

# Définition des fichiers d'entrée
GENOME_FILE = f"genome_complet/genome/{SPECIES}.fna"
GENOME_CLEAR_FILE = f"genome_complet/genome_clear/{SPECIES}_clear.fna"
GLIMMER_GFF = f"glimmer/{SPECIES}/{SPECIES}_I_{I_VALUE}/glimmer_{SPECIES}_I_{I_VALUE}.gff"
PRODIGAL_GFF = f"prodigal/{SPECIES}/{SPECIES}_I_{I_VALUE}/result_{SPECIES}_I_{I_VALUE}"
AUGUSTUS_GFF = f"augustus/{SPECIES}/{SPECIES}_I_{I_VALUE}/result_{SPECIES}_I_{I_VALUE}_augustus.gff"
FUNANNOTATE_GFF = f"funannotate/{SPECIES}/{SPECIES}_I_{I_VALUE}/predict_results/{SPECIES}.gff3"

# Vérifier si Augustus et Funannotate sont disponibles
USE_AUGUSTUS = os.path.exists(AUGUSTUS_GFF)
USE_FUNANNOTATE = os.path.exists(FUNANNOTATE_GFF)

# Définition des fichiers de sortie
OUTPUT_CDS = {
	"glimmer": f"glimmer/{SPECIES}/{SPECIES}_I_{I_VALUE}/{SPECIES}_{I_VALUE}_CDS_glimmer.fasta",
	"prodigal": f"prodigal/{SPECIES}/{SPECIES}_I_{I_VALUE}/{SPECIES}_{I_VALUE}_CDS_prodigal.fasta",
	"augustus": f"augustus/{SPECIES}/{SPECIES}_I_{I_VALUE}/{SPECIES}_{I_VALUE}_CDS_augustus.fasta" if USE_AUGUSTUS else None,
	"funannotate": f"funannotate/{SPECIES}/{SPECIES}_I_{I_VALUE}/{SPECIES}_{I_VALUE}_CDS_funannotate.fasta" if USE_FUNANNOTATE else None
}
OUTPUT_PROT = {
	"glimmer": f"glimmer/{SPECIES}/{SPECIES}_I_{I_VALUE}/{SPECIES}_{I_VALUE}_PROT_glimmer.fasta",
	"prodigal": f"prodigal/{SPECIES}/{SPECIES}_I_{I_VALUE}/{SPECIES}_{I_VALUE}_PROT_prodigal.fasta",
	"augustus": f"augustus/{SPECIES}/{SPECIES}_I_{I_VALUE}/{SPECIES}_{I_VALUE}_PROT_augustus.fasta" if USE_AUGUSTUS else None,
	"funannotate": f"funannotate/{SPECIES}/{SPECIES}_I_{I_VALUE}/{SPECIES}_{I_VALUE}_PROT_funannotate.fasta" if USE_FUNANNOTATE else None
}

# Fonction pour lire un fichier FASTA en dictionnaire
def load_genome(fasta_file):
	genome = {}
	chrom_name = ""
	if os.path.exists(fasta_file):
		with open(fasta_file, 'r') as f:
			for line in f:
				line = line.strip()
				if line.startswith(">"):
					chrom_name = line[1:]
					genome[chrom_name] = ""
				else:
					genome[chrom_name] += line.upper()
	return genome if genome else None

# Charger le génome (essayer le genome normal, sinon le genome clear)
genome_dict = load_genome(GENOME_FILE) or load_genome(GENOME_CLEAR_FILE)

# Charger le génome spécifique pour Funannotate
genome_funannotate = load_genome(GENOME_CLEAR_FILE)

# Fonction de traduction d'ADN en protéine
def dna_to_protein(dna_seq):
	return str(Seq(dna_seq).translate(to_stop=True))

# Extraction spécifique pour Glimmer
def glimmer_extraction(gff_file, genome_dict, output_cds, output_prot):
	if not os.path.exists(gff_file):
		return
	
	cds_records, prot_records = [], []
	with open(gff_file, 'r') as f:
		for line in f:
			if not line.startswith("#"):
				cols = line.strip().split()
				if len(cols) >= 7:
					chrom, start, stop, strand = cols[0], int(cols[3]), int(cols[4]), cols[6]
					if chrom in genome_dict:
						gene_seq = genome_dict[chrom][start-1:stop]
						if strand == "-":
							gene_seq = str(Seq(gene_seq).reverse_complement())
						protein_seq = dna_to_protein(gene_seq)
						cds_records.append(SeqRecord(Seq(gene_seq), id=f"{chrom}:{start}-{stop}({strand})", description=""))
						prot_records.append(SeqRecord(Seq(protein_seq), id=f"{chrom}:{start}-{stop}({strand})", description=""))
	
	if output_cds:
		write(cds_records, output_cds, "fasta")
	if output_prot:
		write(prot_records, output_prot, "fasta")

# Extraction des séquences
def extract_sequences(gff_file, genome_dict, output_cds, output_prot):
	if not gff_file or not os.path.exists(gff_file):
		return
	gene_data = defaultdict(list)
	with open(gff_file, 'r') as f:
		for line in f:
			if not line.startswith("#"):
				cols = line.strip().split()
				if cols[2] == "CDS":
					gene_data[cols[0]].append([int(cols[3]), int(cols[4]), cols[6]])
	cds_records, prot_records = [], []
	for chromosome, genes in gene_data.items():
		if chromosome in genome_dict:
			for start, stop, strand in genes:
				gene_seq = genome_dict[chromosome][start-1:stop]
				if strand == "-":
					gene_seq = str(Seq(gene_seq).reverse_complement())
				protein_seq = dna_to_protein(gene_seq)
				cds_records.append(SeqRecord(Seq(gene_seq), id=f"{chromosome}:{start}-{stop}({strand})", description=""))
				prot_records.append(SeqRecord(Seq(protein_seq), id=f"{chromosome}:{start}-{stop}({strand})", description=""))
	if output_cds:
		write(cds_records, output_cds, "fasta")
	if output_prot:
		write(prot_records, output_prot, "fasta")

glimmer_extraction(GLIMMER_GFF, genome_dict, OUTPUT_CDS["glimmer"], OUTPUT_PROT["glimmer"])

gff_files = {"prodigal": PRODIGAL_GFF, "augustus": AUGUSTUS_GFF, "funannotate": FUNANNOTATE_GFF}
for tool, gff_file in gff_files.items():
	if tool == "funannotate":
		extract_sequences(gff_file, genome_dict, OUTPUT_CDS[tool], OUTPUT_PROT[tool])
	elif tool != "glimmer":
		extract_sequences(gff_file, genome_dict, OUTPUT_CDS[tool], OUTPUT_PROT[tool])
		
		
		
		

