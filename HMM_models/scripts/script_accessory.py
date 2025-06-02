import os
import shutil
import re
from collections import defaultdict
import sys
import glob

# Vérification des arguments
def usage():
	print("Usage: python script_accessory.py <nom_microsporidie> <valeur_I>")
	sys.exit(1)

if len(sys.argv) != 3:
	usage()

SPECIES = sys.argv[1]
I_VALUE = sys.argv[2]

# Définir le chemin de base
base_dir = f"orthofinder/results_{SPECIES}_I_{I_VALUE}"

# Trouver le dossier Results_*
results_dir = glob.glob(os.path.join(base_dir, "Results_*"))
if not results_dir:
	raise FileNotFoundError(f"Aucun dossier Results_* trouvé dans {base_dir}")
results_dir = results_dir[0]  # Prendre le premier résultat

# Définir les chemins
ORTHOGROUP_DIR = os.path.join(results_dir, "Orthogroup_Sequences")
CORE_SOURCE_DIR = os.path.join(results_dir, "Single_Copy_Orthologue_Sequences")
PANGENOME_DIR = os.path.join(base_dir, "pangenome")
ACCESSORY_DIR = os.path.join(base_dir, "accessory_genome")

# Vérifier si la microsporidie est dans la base de données
DATABASE_DIR = "orthofinder/database/Proteomes"
PROTEOME_DIR = f"orthofinder/database/Proteomes_no_{SPECIES}" if os.path.exists(f"{DATABASE_DIR}/Proteomes_{SPECIES}.fasta") else DATABASE_DIR

# Création des dossiers de sortie
os.makedirs(ACCESSORY_DIR, exist_ok=True)
os.makedirs(PANGENOME_DIR, exist_ok=True)

# Fonction pour normaliser les noms des gènes
def normalize_gene_name(gene):
	return re.sub(r'[:()]+', '_', gene)  # Remplace ":" et "()" par "_"

# Étape 1 : Récupérer l'association gène -> protéome
gene_to_proteome = {}

for proteome_file in os.listdir(PROTEOME_DIR):
	if proteome_file.endswith(".fasta"):
		proteome_name = os.path.splitext(proteome_file)[0]
		with open(os.path.join(PROTEOME_DIR, proteome_file), "r") as f:
			for line in f:
				if line.startswith(">"):
					gene_id = line.strip().split()[0][1:]
					gene_id = normalize_gene_name(gene_id)
					gene_to_proteome[gene_id] = proteome_name

print(f"Nombre total de gènes indexés : {len(gene_to_proteome)}")

# Étape 2 : Sélectionner l'accessory genome
accessory_count = 0
ignored_count = 0

for file in os.listdir(ORTHOGROUP_DIR):
	if file.endswith(".fa"):
		file_path = os.path.join(ORTHOGROUP_DIR, file)
		proteome_counts = defaultdict(int)
		total_sequences = 0
		valid_genes = True

		with open(file_path, "r") as f:
			for line in f:
				if line.startswith(">"):
					total_sequences += 1
					gene_id = line.strip().split()[0][1:]
					gene_id = normalize_gene_name(gene_id)
					
					if gene_id in gene_to_proteome:
						proteome_name = gene_to_proteome[gene_id]
						proteome_counts[proteome_name] += 1
					else:
						valid_genes = False
						break

		if not valid_genes or not (2 <= total_sequences <= len(os.listdir(PROTEOME_DIR))):
			ignored_count += 1
			continue

		if 2 <= len(proteome_counts) < len(os.listdir(PROTEOME_DIR)):
			shutil.copy(file_path, os.path.join(ACCESSORY_DIR, file))
			accessory_count += 1
			


# Étape 3 : Fusionner accessory genome et core genome dans pangenome
for file in os.listdir(ACCESSORY_DIR):
	shutil.copy(os.path.join(ACCESSORY_DIR, file), os.path.join(PANGENOME_DIR, file))

for file in os.listdir(CORE_SOURCE_DIR):
	if file.endswith(".fa"):
		shutil.copy(os.path.join(CORE_SOURCE_DIR, file), os.path.join(PANGENOME_DIR, file))

shutil.rmtree(ACCESSORY_DIR)

# Résumé
print(f"\nTri terminé :")
print(f"\t- {len(os.listdir(CORE_SOURCE_DIR))} fichiers core genome ajoutés depuis {CORE_SOURCE_DIR}.")
print(f"\t- Total dans pangenome : {len(os.listdir(PANGENOME_DIR))} fichiers.")
print(f"\t- {ignored_count} fichiers ignorés.")

