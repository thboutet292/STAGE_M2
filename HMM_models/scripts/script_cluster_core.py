#!/usr/bin/env python3

from collections import defaultdict
import re
import sys
import os
import csv
import matplotlib.pyplot as plt
import numpy as np
from venn import venn
import pandas as pd

# Fonction pour récupérer les séquences des gènes prédits
def gene_predict(fichier):
		dico = defaultdict(str)
		with open(fichier, 'r') as f:
				for line in f:
						line = line.rstrip()
						if line.startswith('>'):
								contig = line[1:]
						else:
								dico[contig] += line
		return dico

# Fonction pour parser les fichiers de cluster de CD-HIT-2D
def parse_cdhit(fichier):
		dico = defaultdict(list)
		with open(fichier, 'r') as f:
				for line in f:
						line = line.rstrip()
						if line.startswith('>'):
								cluster = line[1:]
						else:
								dico[cluster].append(line.split()[2][1:-3])  # Enlever le chevron et les points à la fin
								dico[cluster].append(line.split()[1][0:-1])  # Longueur en aa
		return dico

# Fusion des dictionnaires pour obtenir un seul contenant les gènes clusterisés et non-clusterisés
def combined_dico(clusters, supclusters):
		for cluster, liste in clusters.items():
				if len(liste) == 2:
						for supcluster, gene_predits in supclusters.items():
								if len(gene_predits) > 2 and liste[0] == gene_predits[2]:
										clusters[cluster].append(gene_predits[0])
										clusters[cluster].append(gene_predits[1])
		return clusters

# Compter le nombre de gènes prédits et clusterisés
def compt_gene(dico):
		compt = 0
		for _, liste in dico.items():
				if len(liste) > 2:
						compt += 1
		return compt

# Trouver les gènes correctement prédits (même taille)
def same_gene(clusters):
		compt = 0
		for _, liste in clusters.items():
				if len(liste) == 4:
						if liste[1] == liste[3]:
								compt += 1
		return compt

# Vérification des arguments
if len(sys.argv) != 3:
		print("Usage: script_cluster_core.py <nom_microsporidie> <i_value>")
		sys.exit(1)

SPECIES = sys.argv[1]
I_VALUE = sys.argv[2]
RESULT_DIR = f"result/{SPECIES}/{SPECIES}_{I_VALUE}"

# Créer les répertoires nécessaires
def create_directory(directory):
		os.makedirs(directory, exist_ok=True)
		
PROTEOME_FILE = f"Proteomes/Proteomes_{SPECIES}.txt"
dico_microannot = gene_predict(PROTEOME_FILE)
VENN_DIR = f"{RESULT_DIR}/Venn"
create_directory(VENN_DIR)
ERRORS_DIR = f"{RESULT_DIR}/{SPECIES}_{I_VALUE}_errors.csv"
os.makedirs(os.path.dirname(ERRORS_DIR), exist_ok=True)
BAR_DIR = f"{RESULT_DIR}/Bar_tools"
os.makedirs(os.path.dirname(BAR_DIR), exist_ok=True)

#AJOUTER LES OUTILS EN FONCTIONS (SI AUGUSTUS GALERE FAIRE QUE GLIMMER ET PRODIGAL)
TOOLS =  ["glimmer", "prodigal"]

clusters = {}
supclusters = {}

for tool in TOOLS:
		clusters[tool] = parse_cdhit(f"{tool}/{SPECIES}/{SPECIES}_I_{I_VALUE}/cluster/cluster_{SPECIES}_{tool}.clstr")
		supclusters[tool] = parse_cdhit(f"{tool}/{SPECIES}/{SPECIES}_I_{I_VALUE}/cluster/cluster_sup_{SPECIES}_{tool}.clstr")

complete_clusters = {tool: combined_dico(clusters[tool], supclusters[tool]) for tool in TOOLS}


# Extraction des gènes prédits et non prédits
def predicted_genes(clusters):
		liste_genes_predicted = []
		liste_genes_unpredicted = []
		for cluster, liste in clusters.items():		
				if len(liste) > 2:  # Si mon cluster contient un gène prédit
						liste_genes_predicted.append(liste[0])
				else:  # Sinon, il est non prédit
						liste_genes_unpredicted.append(liste[0])
		return liste_genes_predicted, liste_genes_unpredicted

tool_predicted = {tool: predicted_genes(complete_clusters[tool])[0] for tool in TOOLS}

# Fonction pour compter le nombre de séquences dans un fichier FASTA
def count_fasta_sequences(fasta_file):
		count = 0
		with open(fasta_file, 'r') as f:
				for line in f:
						if line.startswith('>'):
								count += 1
		return count

# Compter le nombre de séquences dans chaque fichier FASTA et les afficher
fasta_counts = {}
for tool in TOOLS:
		fasta_file = f"{tool}/{SPECIES}/{SPECIES}_I_{I_VALUE}/{SPECIES}_{I_VALUE}_PROT_{tool}.fasta"
		if os.path.exists(fasta_file):
				fasta_counts[tool] = count_fasta_sequences(fasta_file)
		else:
				fasta_counts[tool] = 0

def venn_tools(tool_predicted, microsporidie):
		sets = {tool.capitalize(): set(genes) for tool, genes in tool_predicted.items()}
		plt.figure(figsize=(10, 10))
		venn(sets, fmt="{size}")
		plt.title(f"Venn diagram for {microsporidie} ({I_VALUE})", fontweight='bold')
		plt.savefig(f"{RESULT_DIR}/Venn/all_tools_{microsporidie}.png")
		plt.close()

venn_tools(tool_predicted, SPECIES)

# Fonction pour compter le nombre de séquences dans un fichier FASTA
def count_fasta_sequences(fasta_file):
		count = 0
		with open(fasta_file, 'r') as f:
				for line in f:
						if line.startswith('>'):
								count += 1
		return count

def calculate_and_save_metrics(tools, proteome_file, species, result_dir, database_ref):
		proteome_count = count_fasta_sequences(proteome_file)
		metrics_data = []

		print("Nombre de gènes prédits retrouvés dans le protéomes")
		for tool in tools:
				fasta_file = f"{tool}/{species}/{species}_I_{I_VALUE}/{species}_{I_VALUE}_PROT_{tool}.fasta"
				total_predicted = count_fasta_sequences(fasta_file) if os.path.exists(fasta_file) else 0

				tp = len(tool_predicted.get(tool, []))  # True Positives
				fn = proteome_count - tp  # False Negatives
				fp = total_predicted - tp  # False Positives

				sensitivity = (tp / proteome_count) * 100 if proteome_count > 0 else 0
				specificity = (tp / (tp + fp)) * 100 if (tp + fp) > 0 else 0

				print(f"{tool}\t {tp}")
				print(f"Sensibilité de {tool}\t {sensitivity:.2f}%")
				print(f"Spécificité de {tool}\t {specificity:.2f}%")

				metrics_data.append([
						tool.capitalize(),
						total_predicted,
						tp,
						round(sensitivity, 2),
						round(specificity, 2),
						fp,
						fn,
						database_ref,
						proteome_count  # Ajout du nombre de gènes dans la base de données
				])

		print(f"{species}\t {proteome_count}")

		# Sauvegarde des métriques dans un CSV en incluant la référence de la base de données et le nombre de gènes
		metrics_df = pd.DataFrame(
				metrics_data,
				columns=["Tool", "Total Predicted", "True Positives", "Sensitivity (%)", "Specificity (%)", "False Positives", "False Negatives", "Database Reference", "Database Gene Count"]
		)

		os.makedirs(result_dir, exist_ok=True)
		metrics_path = f"{result_dir}/metrics_summary_{SPECIES}.csv"
		metrics_df.to_csv(metrics_path, index=False)
		print(f"\nMétriques sauvegardées dans : {metrics_path}")
		print(metrics_df.to_string(index=False))

# Appel de la fonction avec les paramètres du script et la base de données en référence
database_reference = os.path.basename(PROTEOME_FILE)  # Utilise le nom du fichier comme référence de la base de données
calculate_and_save_metrics(TOOLS, PROTEOME_FILE, SPECIES, RESULT_DIR, database_reference)

def extract_cluster_errors(cluster_dict, output_file):
	errors = []
	all_clusters = []
	single_gene_clusters = 0

	if os.path.isdir(output_file):
		print(f"Erreur : {output_file} est un dossier, pas un fichier.")
		return

	if not cluster_dict:
		print("Aucun cluster trouvé dans complete_clusters.")
		return

	for tool, clusters in sorted(cluster_dict.items()):  # S'assurer que les clusters sont traités dans l'ordre
		for cluster_id, values in sorted(clusters.items(), key=lambda x: int(x[0].split()[-1])):
			gene_names = values[::2]  # Récupérer les noms de contigs (indices pairs)
			aa_lengths = [int(v.replace("aa", "")) for v in values[1::2]]  # Récupérer les tailles en aa (indices impairs)

			# Vérifier si c'est un cluster avec un seul gène
			if len(values) <= 2:
				single_gene_clusters += 1
				all_clusters.append((tool, cluster_id, "; ".join(gene_names), "Gène non prédit"))
			else:
				# Vérifier si au moins deux gènes ont les mêmes coordonnées
				unique_coords = set()
				matching_genes = False
				for gene in gene_names:
					match = re.match(r'([^:]+):c?\(?([\d]+)-([\d]+)\)?(?:\([-+]\))?', gene)
					if match:
						contig, start, end = match.groups()
						start, end = int(start), int(end)
						if (start, end) in unique_coords:
							matching_genes = True
							break
						unique_coords.add((start, end))

				if matching_genes:
					aa_status = "OK"
				else:
					# Vérifier si toutes les longueurs AA sont identiques
					if len(set(aa_lengths)) == 1:
						aa_status = "OK"
					else:
						aa_status = "Différent"

					# Identifier la différence si les AA sont différents
					if aa_status == "Différent" and len(gene_names) > 1:
						gene_diffs = []
						for i in range(len(gene_names) - 1):
							match1 = re.match(r'([^:]+):c?\(?([\d]+)-([\d]+)\)?(?:\([-+]\))?', gene_names[i])
							match2 = re.match(r'([^:]+):c?\(?([\d]+)-([\d]+)\)?(?:\([-+]\))?', gene_names[i + 1])
							if match1 and match2:
								contig1, start1, end1 = match1.groups()
								contig2, start2, end2 = match2.groups()
								start1, end1, start2, end2 = map(int, [start1, end1, start2, end2])

								if "(-)" in gene_names[i] or "(-)" in gene_names[i + 1]:
									if end1 != end2:
										gene_diffs.append(f"Différence au début: {end1} vs {end2}")
								else:
									if start1 != start2:
										gene_diffs.append(f"Différence au début: {start1} vs {start2}")
									if end1 != end2:
										gene_diffs.append(f"Différence à la fin: {end1} vs {end2}")

							aa_status += " (" + ", ".join(gene_diffs) + ")"

				all_clusters.append((tool, cluster_id, "; ".join(gene_names), aa_status))

	print(f"Clusters avec un seul gène: {single_gene_clusters}")
	print(f"Total clusters traités: {len(all_clusters)}")

	if not all_clusters:
		print("Aucune donnée à enregistrer.")
		return

	os.makedirs(os.path.dirname(output_file), exist_ok=True)

	with open(output_file, 'w', newline='') as csvfile:
		csv_writer = csv.writer(csvfile, delimiter='\t')
		csv_writer.writerow(["Outil", "Cluster", "Gènes", "État AA"])
		for entry in all_clusters:
			csv_writer.writerow([entry[0], entry[1], entry[2], entry[3]])

	print(f"{len(all_clusters)} clusters sauvegardés dans {output_file}")
	return all_clusters

extract_cluster_errors(complete_clusters, ERRORS_DIR)



output_bar = os.path.join(BAR_DIR, "comparaison_outils2.png")

# Charger le fichier en précisant le séparateur
df = pd.read_csv(ERRORS_DIR, sep="\t")

# Nettoyer la colonne "État AA" pour uniformiser les catégories
df["État AA"] = df["État AA"].str.replace(r"Différent.*", "Misspredicted", regex=True)
df["État AA"] = df["État AA"].replace({"OK": "Predicted", "Gène non prédit": "Unpredicted"})

# Filtrer pour Glimmer et Prodigal
df_glimmer = df[df['Outil'] == 'glimmer']
df_prodigal = df[df['Outil'] == 'prodigal']

# Catégories d'intérêt et couleurs associées
categories = ['Predicted', 'Misspredicted', 'Unpredicted']
colors = ['#388E3C', '#FBC02D', '#D32F2F']

def get_proportions(df):
    counts = df['État AA'].value_counts(normalize=True) * 100
    return [counts.get(cat, 0) for cat in categories]

# Obtenir les proportions
glimmer_props = get_proportions(df_glimmer)
prodigal_props = get_proportions(df_prodigal)

# Création du graphique empilé
fig, ax = plt.subplots(figsize=(8, 6))
x = np.arange(2)  # Une barre pour Glimmer et une pour Prodigal
bottom = np.zeros(2)

for i, (cat, color) in enumerate(zip(categories, colors)):
    values = [glimmer_props[i], prodigal_props[i]]
    bars = ax.bar(x, values, label=cat, bottom=bottom, color=color)
    bottom += values  # Mise à jour du bas de la barre
    
    # Ajouter les annotations
    for bar in bars:
        height = bar.get_height()
        if height > 0:
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_y() + height/2, f'{height:.1f}%', ha='center', va='center', fontsize=10, color='black', fontweight='bold')

ax.set_xticks(x)
ax.set_xticklabels(['Glimmer', 'Prodigal'])
ax.set_ylabel("Proportion (%)")
ax.set_title(f"Comparaison des outils de prédictions pour {SPECIES} ({I_VALUE})")
ax.legend()
ax.set_ylim(0, 100)
ax.grid(axis="y", linestyle="--", alpha=0.7)

# Sauvegarder le graphique
plt.savefig(output_bar, dpi=200, bbox_inches='tight')




