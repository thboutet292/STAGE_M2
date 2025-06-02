import os
import sys

def parse_glimmer_result(input_file, output_file):
	""" Convertit un fichier Glimmer en GFF3 """
	gene_id = 0

	if not os.path.exists(input_file):
		print(f"Erreur : Fichier introuvable {input_file}")
		return

	with open(input_file, "r") as f_in, open(output_file, "w") as f_out:
		# Ajouter l'en-tête GFF3
		f_out.write("##gff-version 3\n")

		for line in f_in:
			line = line.strip()
			if line.startswith(">"):
				chromosome = line[1:]
			else:
				fields = line.split()
				gene_id += 1
				start, end = map(int, (fields[1], fields[2]))

				strand = "-" if start > end else "+"
				start, end = min(start, end), max(start, end)

				f_out.write(f"{chromosome}\tGlimmer\tgene\t{start}\t{end}\t.\t{strand}\t.\tID=gene{gene_id};Name=gene{gene_id}\n")

	print(f"Fichier GFF3 généré : {output_file}")


# Vérifier les arguments
if len(sys.argv) != 3:
	print("Usage: python script_gff.py <nom_microsporidie> <i_value>")
	sys.exit(1)

# Récupérer le nom de la microsporidie et la valeur de I depuis l'argument
species = sys.argv[1]
I_VALUE = sys.argv[2]

# Définition des chemins d'entrée et de sortie
input_path = f"glimmer/{species}/{species}_I_{I_VALUE}/result_{species}_I_{I_VALUE}.predict"
output_path = f"glimmer/{species}/{species}_I_{I_VALUE}/glimmer_{species}_I_{I_VALUE}.gff"

# Exécution de la conversion
print(f"Traitement de {species}...")
parse_glimmer_result(input_path, output_path)

print("Fin du traitement.")

