#!/usr/bin/env python3
import sys
import os
import re
from Bio import SeqIO
from Bio.Seq import Seq

def get_short_id(s):
	"""
	Extrait la partie courte de l'identifiant d'un record.
	Exemples :
	  - Pour un header FASTA comme ">CM017447.1 Antonospora locustae strain CLX chromosome 1, whole genome shotgun sequence",
	    la fonction retournera "CM017447.1".
	  - Pour un record EMBL dont l'ID est "CM017455.1ANTONOSPOR", elle retournera "CM017455.1".
	"""
	m = re.search(r'(CM\d+\.\d+)', s)
	if m:
		return m.group(1)
	else:
		return s.split()[0]

def contains_polyA_signal(window):
	"""
	Parcourt la fenêtre (chaîne de nucléotides) à la recherche d'une sous-séquence qui 
	correspond à l'un des motifs autorisés ("AATAAA" ou "AATTAAA") avec au plus 1 mutation.
	La comparaison se fait uniquement sur des candidats dont la longueur est égale à celle du motif.
	Si un tel candidat est trouvé, la fonction retourne un tuple contenant :
		- le candidat (la sous-séquence trouvée),
		- son indice relatif dans la fenêtre,
		- et le motif autorisé correspondant.
	Sinon, elle retourne (None, None, None).
	"""
	allowed_signals = ["AATAAA", "AATTAAA"]							# liste des motifs polyA autorisés
	for motif in allowed_signals:             						# boucle sur chaque motif autorisé
		motif_len = len(motif)                						# calcule la longueur du motif
		for i in range(len(window) - motif_len + 1):  					# parcourt la fenêtre pour extraire une sous-séquence
			candidate = window[i:i+motif_len] 					# extrait la sous-séquence de la longueur du motif
			mismatches = sum(1 for a, b in zip(candidate, motif) if a != b)  	# compte le nombre de différences
			if mismatches <= 1:           						# tolérance d'une mutation
				return candidate, i, motif  					# retourne le candidat, sa position relative et le motif utilisé
	return None, None, None                   						# sinon, aucun motif n'a été trouvé

def find_orfs(seq, strand="+", min_length=80):
	"""
	Recherche des ORFs dans une séquence donnée pour un cadre de lecture.
	Un ORF commence par ATG et se termine par l'un des codons stop (TAA, TAG, TGA).
	Seuls les ORFs d'au moins `min_length` nt sont retenus, et il faut que dans les 20
	nucléotides en amont de l'ATG, il y ait soit le signal "GGG" ou "CCC" (signal fort),
	soit une région riche en AT (>= 80% sur 15 nt, signal faible).
	Pour les petits ORFs (< 244 nt), une recherche d'un signal polyA est effectuée dans une
	fenêtre allant de -5 à +60 nt autour du codon stop. Si le signal est trouvé, le motif et
	ses coordonnées absolues dans la séquence sont enregistrés.
	Retourne une liste de dictionnaires contenant les informations de chaque ORF.
	"""
	start_codon = "ATG"
	stop_codons = {"TAA", "TAG", "TGA"}
	orfs = [] 							
	
	# Parcourt chacun des 3 cadres de lecture pour le brin courant
	for frame in range(3):
		i = frame
		while i <= len(seq) - 3:
			codon = seq[i:i+3]
			if codon == start_codon:
				# Recherche du premier codon stop dans le même cadre
				for j in range(i+3, len(seq)-2, 3):
					curr_codon = seq[j:j+3]
					if curr_codon in stop_codons:
						orf_length = j + 3 - i
						if orf_length >= min_length:
							if i >= 20:
								upstream = seq[i-20:i]
								# Vérifie la présence d'un signal fort ("GGG" ou "CCC")
								if ("GGG" in upstream) or ("CCC" in upstream):
									valid = True
									signal_type = "strong"
								else:
									at_count = upstream.count("A") + upstream.count("T")
									if at_count / 15 >= 0.8:
										valid = True
										signal_type = "weak"
									else:
										valid = False
								if valid:
									orf_data = {
										"strand": strand,
										"frame": frame,
										"start": i,
										"end": j+3,
										"sequence": seq[i:j+3],
										"upstream": upstream,
										"signal_type": signal_type
									}
									# Pour les petits ORFs (< 244 nt), recherche du signal polyA
									if orf_length < 244:
										polyA_window_start = max(j - 5, 0)
										polyA_window_end = min(j + 60, len(seq))
										polyA_window = seq[polyA_window_start:polyA_window_end]
										polyA_signal, relative_index, motif_used = contains_polyA_signal(polyA_window)
										if polyA_signal is None:
											break
										else:
											polyA_abs_start = polyA_window_start + relative_index
											polyA_abs_end = polyA_abs_start + len(polyA_signal)
											orf_data["polyA_signal"] = polyA_signal
											orf_data["polyA_coords"] = (polyA_abs_start, polyA_abs_end)
									orfs.append(orf_data)
						break
			i += 3
	return orfs

def scan_genome(input_fasta):
	"""
	Parcourt chaque record du fichier FASTA et recherche les ORFs sur les deux brins.
	Pour le brin direct, la séquence est utilisée telle quelle.
	Pour le brin inverse, on calcule le reverse complement et on ajuste les coordonnées
	pour retrouver leur position d'origine.
	Retourne une liste de dictionnaires, chacun contenant l'ID du record (identifiant court)
	et ses ORFs détectées.
	"""
	results = []
	for record in SeqIO.parse(input_fasta, "fasta"):
		short_id = get_short_id(record.id)
		seq = str(record.seq).upper()
		record_results = {"id": short_id, "orfs": []}
		
		# ORFs sur le brin direct
		orfs_forward = find_orfs(seq, strand="+")
		record_results["orfs"].extend([
			{
				"record_id": short_id,
				"strand": orf["strand"],
				"frame": orf["frame"],
				"start": orf["start"],
				"end": orf["end"],
				"sequence": orf["sequence"],
				"upstream": orf["upstream"],
				"signal_type": orf["signal_type"],
				"polyA_signal": orf.get("polyA_signal"),
				"polyA_coords": orf.get("polyA_coords")
			} for orf in orfs_forward
		])
		
		# ORFs sur le brin inverse
		rev_seq = str(Seq(seq).reverse_complement())
		orfs_reverse = find_orfs(rev_seq, strand="-")
		for orf in orfs_reverse:
			L = len(seq)
			orig_start = L - orf["end"]
			orig_end = L - orf["start"]
			record_results["orfs"].append({
				"record_id": short_id,
				"strand": orf["strand"],
				"frame": orf["frame"],
				"start": orig_start,
				"end": orig_end,
				"sequence": orf["sequence"],
				"upstream": orf["upstream"],
				"signal_type": orf["signal_type"],
				"polyA_signal": orf.get("polyA_signal"),
				"polyA_coords": orf.get("polyA_coords")
			})
		results.append(record_results)
	return results

def parse_embl_cds_from_directory(embl_dir):
	"""
	Parcourt le répertoire donné et parse chaque fichier ayant l'extension ".embl".
	Retourne un dictionnaire dont les clés sont les IDs courts des records et les valeurs sont 
	des listes de tuples (start, end, strand) pour chaque CDS annotée.
	"""
	annotated_cdss = {}
	for filename in os.listdir(embl_dir):
		if filename.endswith(".embl"):
			filepath = os.path.join(embl_dir, filename)
			for record in SeqIO.parse(filepath, "embl"):
				short_id = get_short_id(record.id)
				cds_list = []
				for feature in record.features:
					if feature.type == "CDS":
						cds_start = int(feature.location.start)
						cds_end = int(feature.location.end)
						cds_strand = feature.location.strand
						cds_list.append((cds_start, cds_end, cds_strand))
				if short_id in annotated_cdss:
					annotated_cdss[short_id].extend(cds_list)
				else:
					annotated_cdss[short_id] = cds_list
	return annotated_cdss

def filter_nested_orfs(orfs):
	"""
	Filtre les ORFs imbriquées (sur le même brin et même cadre).
	Pour deux ORFs imbriquées, on applique la règle suivante :
	  - Si l'ORF en aval est complètement contenue dans l'ORF en amont et que la différence 
	    entre leurs positions de départ est d'au moins 30 nt (plus de 10 acides aminés), on garde l'ORF en amont,
	    même si son signal est AT-rich (weak).
	  - Sinon, on conserve par défaut le plus grand, sauf si le plus grand a un signal weak et
	    le plus petit un signal strong (dans ce cas, on garde le plus petit).
	Retourne la liste des ORFs filtrées.
	"""
	filtered = orfs[:]
	changed = True
	while changed:
		changed = False
		remove_indices = set()
		for i in range(len(filtered)):
			for j in range(i+1, len(filtered)):
				if filtered[i]["strand"] == filtered[j]["strand"] and filtered[i]["frame"] == filtered[j]["frame"]:
					if filtered[i]["start"] >= filtered[j]["start"] and filtered[i]["end"] <= filtered[j]["end"]:
						delta = filtered[i]["start"] - filtered[j]["start"]
						if delta >= 30:
							remove_indices.add(i)
						else:
							if filtered[j]["signal_type"] == "weak" and filtered[i]["signal_type"] == "strong":
								remove_indices.add(j)
							else:
								remove_indices.add(i)
					elif filtered[j]["start"] >= filtered[i]["start"] and filtered[j]["end"] <= filtered[i]["end"]:
						delta = filtered[j]["start"] - filtered[i]["start"]
						if delta >= 30:
							remove_indices.add(j)
						else:
							if filtered[i]["signal_type"] == "weak" and filtered[j]["signal_type"] == "strong":
								remove_indices.add(i)
							else:
								remove_indices.add(j)
		if remove_indices:
			filtered = [filtered[k] for k in range(len(filtered)) if k not in remove_indices]
			changed = True
	return filtered

def filter_overlapping_small_orfs(orfs):
	"""
	Filtre les petits ORFs (< 244 nt) qui chevauchent (peu importe le brin ou le cadre)
	avec un grand ORF (>= 244 nt). Les petits ORFs chevauchants sont éliminés.
	Retourne la liste des petits ORFs non chevauchants.
	"""
	large_orfs = [orf for orf in orfs if len(orf["sequence"]) >= 244]
	filtered = []
	for orf in orfs:
		if len(orf["sequence"]) < 244:
			overlap = False
			for large in large_orfs:
				if orf["start"] < large["end"] and orf["end"] > large["start"]:
					overlap = True
					break
			if not overlap:
				filtered.append(orf)
	return filtered

def write_output(orfs_results, cds_filename, prot_filename):
	"""
	Enregistre les séquences CDS et leur traduction en protéines dans deux fichiers FASTA.
	Pour chaque ORF retenu, l'en-tête inclut le nom du contig (identifiant court) et ses coordonnées.
	Seules les ORFs de petite taille (CDS < 244 nt et protéine <= 81 acides aminés) sont écrites.
	"""
	cds_records = []
	prot_records = []
	for record in orfs_results:
		rec_id = record["id"]
		filtered_orfs = filter_nested_orfs(record["orfs"])
		filtered_orfs = filter_overlapping_small_orfs(filtered_orfs)
		for orf in filtered_orfs:
			if orf["strand"] == "+":
				header = f">{rec_id}:{orf['start']+1}-{orf['end']}"
			else:
				header = f">{rec_id}:c({orf['start']+1}-{orf['end']})"
			nt_seq = orf["sequence"]
			prot_seq = str(Seq(nt_seq).translate(to_stop=False))
			if prot_seq.endswith("*"):
				prot_seq = prot_seq[:-1]
			if len(nt_seq) < 244 and len(prot_seq) <= 81:
				cds_records.append(f"{header}\n{nt_seq}")
				prot_records.append(f"{header}\n{prot_seq}")
	with open(cds_filename, "w") as f_cds:
		f_cds.write("\n".join(cds_records) + "\n")
	with open(prot_filename, "w") as f_prot:
		f_prot.write("\n".join(prot_records) + "\n")

def main():
	"""
	Usage: python script_petits_genes.py <input_fasta> <embl_directory> <output_CDS> <output_PROT>
	
	Le script procède ainsi :
	  - Parse le répertoire contenant les fichiers EMBL (MicroAnnot) pour récupérer les CDS annotées.
	  - Scanne le fichier FASTA du génome pour détecter les ORFs candidates sur les 6 cadres.
	  - Filtre les sCDS qui chevauchent (même partiellement) une CDS annotée.
	    Un message de filtrage est affiché uniquement pour les petits gènes (longueur < 244 nt).
	  - Écrit les résultats dans deux fichiers FASTA (séquences CDS et traduction protéique).
	"""
	if len(sys.argv) != 5:
		sys.exit("Usage: python script_petits_genes.py <input_fasta> <embl_directory> <output_CDS> <output_PROT>")
	input_fasta = sys.argv[1]
	embl_directory = sys.argv[2]
	output_cds = sys.argv[3]
	output_prot = sys.argv[4]
	
	# Récupère les CDS annotées depuis le répertoire contenant les fichiers EMBL
	annotated_cdss = parse_embl_cds_from_directory(embl_directory)
	
	# Détecte les ORFs candidates dans le génome (sur 6 cadres)
	orfs_results = scan_genome(input_fasta)
	
	# Pour chaque record, élimine les ORFs qui chevauchent une CDS annotée (même partiellement)
	for record in orfs_results:
		rec_id = record["id"]
		if rec_id in annotated_cdss:
			cds_intervals = annotated_cdss[rec_id]
			new_orfs = []
			for orf in record["orfs"]:
				overlap = False
				for cds in cds_intervals:
					cds_start, cds_end, _ = cds
					# Si la région de l'ORF chevauche la CDS annotée
					if orf["start"] < cds_end and orf["end"] > cds_start:
						# Affiche un message de filtrage uniquement pour les petits gènes (< 244 nt)
						"""
						if (orf["end"] - orf["start"]) < 244:
							print(f"Filtrage : sCDS {rec_id} {orf['start']+1}-{orf['end']} filtré (chevauche CDS annotée {cds_start+1}-{cds_end}).")
						"""
						overlap = True
						break
				if not overlap:
					new_orfs.append(orf)
			record["orfs"] = new_orfs
	write_output(orfs_results, output_cds, output_prot)

if __name__ == "__main__":
	main()

