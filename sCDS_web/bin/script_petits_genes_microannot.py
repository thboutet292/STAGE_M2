#!/usr/bin/env python3
import sys	# Pour la gestion des arguments de la ligne de commande
import os	# Pour la gestion des fichiers et répertoires
import re	# Pour les expressions régulières
from Bio import SeqIO	# Pour lire les fichiers FASTA/EMBL via Biopython
from Bio.Seq import Seq	# Pour manipuler des séquences biologiques avec Biopython

def get_short_id(header):
	"""
	Renvoie l'identifiant du contig en coupant tout ce qui suit un point (.)
	et en ne prenant que le premier "mot" (avant espace ou point-virgule).
	Exemple : 'GL877538.1', 'GL877538', 'GL877538.2 description' --> 'GL877538'
	"""
	id_simple = header.strip().split()[0].replace(';', '')
	return id_simple.split('.')[0]





def contains_polyA_signal(window):
	"""
	Parcourt la fenêtre pour trouver AATAAA ou AATTAAA avec ≤1 mutation.
	Retourne (signal, index_relatif, motif) ou (None, None, None).
	"""
	allowed = ["AATAAA", "AATTAAA"]						# motifs polyA autorisés
	for motif in allowed:							# pour chaque motif
		L = len(motif)							# on garde la longueur du motif
		for i in range(len(window) - L + 1):				# pour chaque position possible
			cand = window[i:i+L]					# prend la séquence locale
			mismatches = sum(1 for a,b in zip(cand,motif) if a!=b)	# compte les différences entre la séquence gardé et les motifs polyA
			if mismatches <= 1:					# accepte le motif avec jusqu'à une erreur d'autorisé (donc soit motif polyA identique, soit dégénéré d'un nt) 
				return cand, i, motif				# retourne la séquence trouvée, son index, et le motif d'origine
	return None, None, None							# si rien n'est trouvé



def find_orfs(seq, strand="+", min_length=80):
	"""
	Recherche d’ORFs :
	 - ATG jusqu'à un codon stop (TAA,TAG,TGA) 
	 - upstream ≥20nt avec signal fort (GGG/CCC) ou AT-rich ≥80%
	 - pour ORF <244nt, recherche polyA dans [-5,+60] autour du stop
	"""
	starts = "ATG"								# codon start
	stops = {"TAA","TAG","TGA"}						# codons stop
	orfs = []								# liste de sortie
	
	for frame in range(3):							# pour chaque cadre de lecture possible
		i = frame
		while i <= len(seq)-3:						# tant qu'on peut prendre un triplet
			if seq[i:i+3] == starts and i>=20:			# cherche un codon start avec >=20 nt en amont
				up = seq[i-20:i]				# prend la séquence en amont
				if ("GGG" in up) or ("CCC" in up):		# si signal fort GGG ou CCC
					valid, stype = True, "strong"		# marque comme valide (signal fort)
				else:
					atc = up.count("A")+up.count("T")	# compte le nombre d'A/T dans les 20 nt en amonts du codon start
					if atc/20 >= 0.8:			# si AT-rich (≥80%)
						valid, stype = True, "weak"	# marque comme valide (signal faible)
					else:
						valid = False			# sinon, non valide car pas de signal trouvé
						
				if valid:							# si l'orf est valide
					for j in range(i+3, len(seq)-2, 3):			# parcours par codons à partir du start
						if seq[j:j+3] in stops:				# si codon stop trouvé
							length = j+3 - i			# resort la longueur de l'ORF
							if length >= min_length:		# vérifie que la séquence fait au moins la taille minimum définit dans la commande (ici 80 nt)
							
								orf = {					# garde des infos sur l'orf 
									"strand": strand,		# le brin ('+' ou '-')
									"frame": frame,			# cadre de lecture
									"start": i,			# position start 
									"end": j+3,			# position end (codon stop inclus)
									"sequence": seq[i:j+3],		# séquence
									"upstream": up,			# séquence région upstream
									"signal_type": stype		# type de signal détecté
								}
								if length < 244:			# si c'est une petite ORF
									w0 = max(j-5, 0)		# début fenêtre polyA
									w1 = min(j+60, len(seq))	# fin fenêtre polyA
									win = seq[w0:w1]		# séquence de la fenêtre
									sig, ridx, mot = contains_polyA_signal(win)	# cherche motifs polyA (via la fonction 'contains_polyA_signal', définit avant)
									if sig is None:			# si aucun résultat
										break			# ignore cette ORF
									orf["polyA_signal"] = sig	#sinon garde le signal
									orf["polyA_coords"] = (w0+ridx, w0+ridx+len(sig))
								orfs.append(orf)			# ajoute l'ORF à la liste
							break  						# sort de la boucle, ne considère qu'un stop par start
			
			i += 3	# avance de 3 nt (un codon)
	return orfs		# retourne la liste d'ORFs correspondant aux critères trouvées


def scan_genome(input_fasta):
	"""
	Lit le FASTA, détecte ORFs sur les deux brins.
	Retourne liste de {"id": accession, "orfs": [...]}.
	"""
	results = []					# liste des résultats pour tous les contigs
	for rec in SeqIO.parse(input_fasta, "fasta"):	# boucle sur tous les contigs/sequences du fasta
		rid = get_short_id(rec.id)		# obtient l'identifiant court (via la fonction 'get_short_id' définit précedemment)
		seq = str(rec.seq).upper()		# met la séquence en majuscules
		entry = {"id": rid, "orfs": []}		# dico pour ce contig
		
		# brin +
		for o in find_orfs(seq, "+"):				# applique la fonction 'find_orfs' pour le brin sens 
			entry["orfs"].append({**o, "record_id": rid})	# ajoute les ORFs trouvées sur +
		
		# brin -
		rc = str(Seq(seq).reverse_complement())			# Calcule le brin -
		L = len(seq)	
		for o in find_orfs(rc, "-"):				# applique la fonction 'find_orfs' pour le brin sens
			s0 = L - o["end"]				# Transforme coord brin - en coord brin +
			e0 = L - o["start"]
			o2 = o.copy()
			o2.update({"start": s0, "end": e0, "record_id": rid})
			entry["orfs"].append(o2)
		results.append(entry)					# ajoute le résultat du contig
	return results	


def parse_embl_cds_from_directory(embl_dir):
	"""
	Parcourt tous les fichiers .embl du répertoire,
	extrait pour chaque CDS le tuple (start,end,strand),
	et indexe par accession extraite avec get_short_id().
	"""
	annotated = {}							# dico des CDS du fichier EMBL
	for fn in os.listdir(embl_dir):					# parcours les fichiers du répertoire
		if not fn.lower().endswith(".embl"):			# cherche que les .embl
			continue
		path = os.path.join(embl_dir, fn)			# chemin complet du fichier
		for rec in SeqIO.parse(path, "embl"):			# pour chaque entrée dans le fichier EMBL
			acc = get_short_id(rec.id)           		# utilisation de la même fonction 'get_short_id' pour prendre le même nom entre le fichier du génome et les .embl
			cdss = []					# liste de CDS trouvées pour ce contig
			for feat in rec.features:			# parcours les features
				if feat.type == "CDS":			# si la feature est de type CDS
					s = int(feat.location.start)	# récupère position start
					e = int(feat.location.end)	# position stop
					st = feat.location.strand	# sens du brin
					cdss.append((s, e, st))		# ajoute ce CDS à la liste
			if cdss:
				annotated.setdefault(acc, []).extend(cdss)	# ajoute les cds au dico
	return annotated


def filter_nested_orfs(orfs):
	"""
	Filtre ORFs imbriquées 
	"""
	filt = orfs[:]	# Copie la liste d’ORFs initiale pour ne pas modifier l’originale
	changed = True	# Booléen pour indiquer si une modification a été faite lors de la boucle
	while changed:	# Boucle tant qu’au moins un ORF est retiré à chaque passage
		changed = False	# Réinitialise le booléen à False au début de chaque boucle
		to_remove = set()	# Initialise un ensemble pour stocker les indices à supprimer
		for i in range(len(filt)):	# Boucle sur chaque ORF par indice i
			for j in range(i+1, len(filt)):	# Boucle sur chaque ORF suivant i, par indice j
				a, b = filt[i], filt[j]	# Récupère les deux ORFs à comparer
				if a["strand"]==b["strand"] and a["frame"]==b["frame"]:	# Compare uniquement les ORFs sur le même brin et frame
					# a dans b ?
					if a["start"]>=b["start"] and a["end"]<=b["end"]:	# Si a est inclus dans b
						d = a["start"] - b["start"]	# Calcule la distance entre le début de a et b
						if d>=30:	# Si cette distance est supérieure ou égale à 30 nt, on retire a
							to_remove.add(i)
						else:
							if b["signal_type"]=="weak" and a["signal_type"]=="strong":	# Si b est un signal faible et a un signal fort, retire b
								to_remove.add(j)
							else:	# Sinon retire a
								to_remove.add(i)
					# b dans a ?
					elif b["start"]>=a["start"] and b["end"]<=a["end"]:	# Si b est inclus dans a
						d = b["start"] - a["start"]	# Calcule la distance entre le début de b et a
						if d>=30:	# Si cette distance est supérieure ou égale à 30 nt, on retire b
							to_remove.add(j)
						else:
							if a["signal_type"]=="weak" and b["signal_type"]=="strong":	# Si a est faible et b fort, retire a
								to_remove.add(i)
							else:	# Sinon retire b
								to_remove.add(j)
		if to_remove:	# Si des indices ont été identifiés à retirer
			filt = [x for idx,x in enumerate(filt) if idx not in to_remove]	# Supprime les ORFs correspondants
			changed = True	# Indique qu’on a modifié la liste, donc on relance la boucle
	return filt	# Retourne la liste filtrée

# ------------------------------------------------------------
# filter_overlapping_small_orfs : supprime petits ORFs chevauchant grands
# ------------------------------------------------------------
def filter_overlapping_small_orfs(orfs):
	"""
	Supprime tout ORF <244nt chevauchant un ORF ≥244nt
	(tout brin ou frame).
	"""
	large = [o for o in orfs if len(o["sequence"])>=244]	# Liste des grands ORFs (≥244 nt)
	smalls = []	# Initialise la liste des petits ORFs à conserver
	for o in orfs:	# Parcourt tous les ORFs
		if len(o["sequence"])<244:	# Si l’ORF est petit (<244 nt)
			overlap = False	# Initialise un indicateur de chevauchement
			for L in large:	# Parcourt chaque grand ORF
				if not (o["end"]<=L["start"] or o["start"]>=L["end"]):	# Si il y a chevauchement
					overlap = True	# Met l’indicateur à True
					break	# Sort de la boucle dès qu’un chevauchement est trouvé
			if not overlap:	# Si pas de chevauchement
				smalls.append(o)	# Ajoute le petit ORF à la liste à garder
	return smalls	# Retourne la liste des petits ORFs qui ne chevauchent pas de grands


def write_output(orfs_results, cds_file, prot_file):
	"""
	Génère deux FASTA :
	  - cds_file : séquences nucléotides (<244nt)
	  - prot_file: protéines (≤81 AA)
	"""
	cds_out = []									# liste pour séquence en nt 
	prot_out = []									# liste pour séquence en aa
	
	for rec in orfs_results:							# voucle sur chaque orf
		cid = rec["id"]								# récupère l’identifiant court du contig
		fo = filter_nested_orfs(rec["orfs"])					# applique le filtre pour supprimer les ORFs imbriqués (via 'filter_nested_orfs')
		fo = filter_overlapping_small_orfs(fo)					# applique le filtre pour supprimer les petits ORFs qui chevauchent les grands (via 'filter_overlapping_orfs')
		for o in fo:								# boucle sur chaque ORF retenu
			nt = o["sequence"]						# récupère la séquence nucléotidique de l’ORF
			prot = str(Seq(nt).translate(to_stop=False)).rstrip("*")	# traduit la séquence en protéine (AA), retire l’étoile finale éventuelle (codon stop)
			
			if len(nt)<244 and len(prot)<=81:				# ne conserve que les CDS de moins de 244 nt et protéines ≤81 AA
				if o["strand"]=="+":					# si ORF sur le brin +
					hdr = f">{cid}:{o['start']+1}-{o['end']} | {o['signal_type']}"	# entête fasta (1-based)
				else:	# Si ORF sur le brin -
					hdr = f">{cid}:c({o['start']+1}-{o['end']}) | {o['signal_type']}"	# entête fasta pour le brin -
					
				cds_out.append(f"{hdr}\n{nt}")				# ajoute la séquence nucléotidique à la liste
				prot_out.append(f"{hdr}\n{prot}")			# ajoute la séquence protéique à la liste
	with open(cds_file, "w") as f:							# ouvre le fichier CDS en écriture
		f.write("\n".join(cds_out)+"\n")					# écrit toutes les séquences nucléotidiques dans le fichier
	with open(prot_file, "w") as f:							# ouvre le fichier protéines en écriture
		f.write("\n".join(prot_out)+"\n")					# écrit toutes les séquences protéiques dans le fichier


def main():
	if len(sys.argv)!=5:									# vérifie que 4 arguments sont bien passés au script
		sys.exit("Usage: python script.py <fastafile> <embl_dir> <out_CDS> <out_PROT>")	# affiche l’aide et quitte si pas le bon nombre d’arguments

	fasta, embl_dir, out_cds, out_prot = sys.argv[1:]					# récupère les 4 arguments : fasta, dossier EMBL, fichier CDS, fichier protéines

	annotated = parse_embl_cds_from_directory(embl_dir)					# appelle la fonction 'parse_embl_cds_from_directory' pour lire les CDS annotées dans tous les fichiers EMBL du dossier
	
	orfs_res = scan_genome(fasta)								# appelle la fonction 'scan_genome' pour détecter les ORFs sur le génome FASTA
	
	for rec in orfs_res:									# parcourt chaque résultat ORF par contig
		cid = rec["id"]									# récupère l’identifiant du contig
		cdss = annotated.get(cid, [])							# récupère les CDS EMBL pour ce contig
		if not cdss:									# si aucune CDS EMBL n’est trouvée
			print(f"[DEBUG] Pas de CDS EMBL pour {cid}")				# affiche un message debug
			continue								# passe au contig suivant
		#print(f"[DEBUG] Contig {cid} — {len(cdss)} CDS annotées : {cdss}")		# affiche un message debug sur le nombre de CDS EMBL
		kept = []									# initialise une liste pour les ORFs à conserver
		for o in rec["orfs"]:								# parcourt chaque ORF détecté
			#print(f"[DEBUG]   Test ORF {o['start']+1}-{o['end']} (len={o['end']-o['start']})")	# affiche la position de l’ORF
			over = False								# booléen pour savoir si l’ORF chevauche une CDS EMBL
			for (s,e,_) in cdss:							# parcourt chaque CDS EMBL du contig
				#print(f"[DEBUG]     VS CDS {s+1}-{e}")				# affiche la position de la CDS EMBL comparée
				if not (o["end"]<=s or o["start"]>=e):				# si il y a chevauchement entre ORF et CDS EMBL
					over = True						# met le booléen à True
					if (o["end"]-o["start"])<244:				# si l’ORF est petit
						print(f"Filtrage : petit CDS {cid} {o['start']+1}-{o['end']} chevauche CDS {s+1}-{e}")	# affiche un message de filtrage
					break							# sort de la boucle sur les CDS EMBL dès qu’un chevauchement est trouvé
			if not over:								# si pas de chevauchement, on conserve l’ORF
				kept.append(o)
		rec["orfs"] = kept								# on met à jour la liste d’ORFs du contig avec ceux conservés

	
	write_output(orfs_res, out_cds, out_prot)	# appelle la fonction 'write_output' pour écrire les fichiers FASTA finaux

if __name__=="__main__":	
	main()


