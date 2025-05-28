#!/usr/bin/env python3

import sys			 
import os			
import csv			
import shutil		
import re			
from Bio import SeqIO			
from Bio.SeqRecord import SeqRecord	

# Seuils de qualité et paramètres globaux
EVALUE_BLAST   = 1e-5			# E-value maximal pour BLASTP
IDENTITY_BLAST = 50.0			# Identité (%) minimale pour BLASTP
COVERAGE_BLAST = 0.7			# Couverture minimale pour BLASTP

EVALUE_DB      = 1e-5			# E-value maximal pour tBLASTn
IDENTITY_DB    = 80.0			# Identité (%) minimale pour tBLASTn
COVERAGE_DB    = 0.0			# Couverture minimale pour tBLASTn

IDENTITY_DUP   = 100.0			# Identité (%) pour duplication exacte
COVERAGE_DUP   = 1.0			# Couverture pour duplication exacte

TRANSLATION_TABLE = 1			# Table de traduction standard Biopython

DEFAULT_BDD1 = "data//database_ortholog/sCDS_ortholog.fasta"
DEFAULT_BDD2 = "data/database_not_ortholog/sCDS_not_ortholog.fasta"
BDD1_PATH    = os.environ.get("BDD1_PATH", DEFAULT_BDD1)	# Chemin BDD1 (orthologues)
BDD2_PATH    = os.environ.get("BDD2_PATH", DEFAULT_BDD2)	# Chemin BDD2 (non-orthologues)

proteins_seen = set()		# Pour éviter d'enregistrer deux fois la même séquence protéique

def normalize_id(s):
	"""Garde tout avant le premier '.', '_' ou ':' pour extraire le nom du contig."""
	return re.split(r"[._:]", s, 1)[0]

def backup(fp):
	"""Crée une sauvegarde fp_backup.fasta si le fichier existe."""
	if os.path.isfile(fp):
		bak = fp.replace('.fasta', '_backup.fasta')
		shutil.copy(fp, bak)
		print(f"[INFO] Backup created: {bak}")
	else:
		print(f"[WARN] {fp} not found, skipped backup.")

def load_fasta(fp):
	"""Charge un fichier FASTA et retourne un dict id -> SeqRecord."""
	return { rec.id: rec for rec in SeqIO.parse(fp, 'fasta') }

def parse_blast(tab):
	"""Parse BLASTP outfmt6 et retourne un dict qid -> liste de hits."""
	hits = {}
	with open(tab) as f:
		reader = csv.reader(f, delimiter='\t')
		for row in reader:
			if len(row) < 12: continue	# Vérifie le bon format BLAST
			qid, sid = row[0], row[1]	# Query et Subject
			pid = float(row[2])			# % identité
			aln = int(row[3])			# longueur alignement
			ev = float(row[10])			# e-value
			hits.setdefault(qid, []).append({
				'sid': sid,
				'sid_norm': normalize_id(sid),
				'pid': pid,
				'alen': aln,
				'ev': ev
			})
	return hits

def parse_tblastn(tab):
	"""Parse tBLASTn outfmt6 -> dict qid -> liste de hits avec contig_range."""
	hits = {}
	with open(tab) as f:
		reader = csv.reader(f, delimiter='\t')
		for row in reader:
			if len(row) < 12: continue
			qid = row[0]				# Query
			full_s = row[1]				# Subject
			pid = float(row[2])			# % identité
			aln = int(row[3])			# longueur alignement
			ev = float(row[10])			# e-value
			ss, se = int(row[8]), int(row[9])	# Start/end subject
			if ss <= se:
				cr = f"{normalize_id(full_s)}:{ss}-{se}"
			else:
				cr = f"{normalize_id(full_s)}:c({se}-{ss})"
			hits.setdefault(qid, []).append({
				'contig_range': cr,
				'contig': normalize_id(full_s),
				'pid': pid,
				'alen': aln,
				'ev': ev
			})
	return hits

def cover(aln, length):
	"""Retourne la couverture de l'alignement (aln/length) ou 0 si length == 0."""
	return aln/length if length else 0

def load_embl_cds(embl_dir):
	"""Parse tous les .embl pour extraire les CDS annotées par contig."""
	annot = {}
	for fn in os.listdir(embl_dir):
		if not fn.endswith('.embl'): continue
		path = os.path.join(embl_dir, fn)
		for rec in SeqIO.parse(path, 'embl'):
			cont = normalize_id(rec.id)
			for feat in rec.features:
				if feat.type=='CDS' and feat.location:
					a = int(feat.location.start)+1	# Début (1-based)
					b = int(feat.location.end)		# Fin (1-based, inclusif)
					annot.setdefault(cont, []).append((a,b))
	return annot

def inside_annotated(contig, s, e, annot):
	"""Vérifie si [s, e] est inclus dans une CDS annotée pour ce contig."""
	for a, b in annot.get(contig, []):
		if s >= a and e <= b:
			return True, (a, b)
	return False, None

def find_longest_M_window_around(seq, start, window_nt=15):
	"""
	Cherche tous les codons ATG dans une fenêtre définit autour de start.
	Pour chaque ATG, traduit jusqu'au premier stop, et retourne la position (offset)
	qui donne la protéine la plus longue. Renvoie (best_offset, best_prot) ou (None, None) si rien trouvé.
	"""
	start0 = start - 1	# start (1-based) vers index 0-based
	best_offset = None
	best_prot = ""
	seq_len = len(seq)
	# Fenêtre centrée sur start0 : de start0-window_nt à start0+window_nt
	for offset in range(-window_nt, window_nt+1, 3):
		pos = start0 + offset
		if pos < 0 or pos+3 > seq_len:
			continue
		codon = str(seq[pos:pos+3])
		if codon == "ATG":
			prot_seq = seq[pos:].translate(table=TRANSLATION_TABLE)
			prot_str = str(prot_seq)
			if "*" in prot_str:
				prot_str = prot_str[:prot_str.index("*")]
			if len(prot_str) > len(best_prot):
				best_prot = prot_str
				best_offset = offset
	if best_offset is not None:
		return best_offset, best_prot
	return None, None

def main():
	# Vérification du nombre d'arguments
	if len(sys.argv)!=8:
		sys.exit("Usage: update_bdd.py blast1 blast2 tblastn prot_fa genome_fa out_fa embl_dir")
	# Récupération des chemins des fichiers d'entrée/sortie
	b1_tab, b2_tab, tbl_tab, prot_fa, genome_fa, out_fa, embl_dir = sys.argv[1:]

	# Sauvegarde des bases avant modification
	backup(BDD1_PATH)
	backup(BDD2_PATH)

	# Chargement des données et annotations
	pred   = load_fasta(prot_fa)
	bdd1   = load_fasta(BDD1_PATH)
	bdd2   = load_fasta(BDD2_PATH)
	genome = { normalize_id(r.id): r for r in SeqIO.parse(genome_fa,'fasta') }
	contigs = set(genome.keys())
	annot   = load_embl_cds(embl_dir)

	# Si le génome est déjà dans la BDD1, export uniquement
	if contigs & { normalize_id(k) for k in bdd1 }:
		print("[INFO] Genome already in BDD1   ->  exporting existing entries")
		with open(out_fa,'w') as out:
			for rid, rec in bdd1.items():
				if normalize_id(rid) in contigs:
					out.write(f">{rec.id} | BDD1\n{rec.seq}\n")
		return

	# Parsing des résultats BLAST
	blast1  = parse_blast(b1_tab)
	blast2  = parse_blast(b2_tab)
	tblastn = parse_tblastn(tbl_tab)

	new1, new_q2, new_s2, new_t = [], [], [], []
	source_map, extracted = {}, {}
	added2 = 0

	# Validation des prédits avec BLASTP
	for qid, rec in pred.items():
		L = len(rec.seq)
		for h in blast1.get(qid, []):
			cov = cover(h['alen'],L)
			if h['pid']==IDENTITY_DUP and cov==COVERAGE_DUP: continue
			if h['ev']<EVALUE_BLAST and h['pid']>=IDENTITY_BLAST and cov>=COVERAGE_BLAST:
				if qid not in bdd1 and len(str(rec.seq)) > 6 and len(str(rec.seq)) < 81:	
					bdd1[qid]=rec
					new1.append(qid)
					source_map[qid]="BDD1"
				break
		for h in blast2.get(qid, []):
			if h['sid_norm'] in contigs:
				print(f"[DEBUG] Skip BLAST2 {qid}  ->  {h['sid']} (same genome)")
				continue
			cov = cover(h['alen'],L)
			if h['pid']==IDENTITY_DUP and cov==COVERAGE_DUP: continue
			if h['ev']<EVALUE_BLAST and h['pid']>=IDENTITY_BLAST and cov>=COVERAGE_BLAST:
				sid = h['sid']
				if sid in bdd2 and sid not in bdd1 and len(str(bdd2[sid].seq)) > 6 and len(str(bdd2[sid].seq)) < 81 :	
					bdd1[sid]=bdd2.pop(sid)
					new_s2.append(sid)
					source_map[sid]="BDD2"
				if qid not in bdd1 and len(str(rec.seq)) > 6 and len(str(rec.seq)) < 81:	
					bdd1[qid]=rec
					new_q2.append(qid)
					source_map[qid]="BDD2"
				break

	# Extraction de nouveaux sCDS depuis tBLASTn
	for qid, hits in tblastn.items():
		if normalize_id(qid) in contigs:
			print(f"[DEBUG] Skip tBLASTn {qid} (query from same genome)")
			continue
		for h in hits:
			seq = genome[h['contig']].seq
			cov = cover(h['alen'], len(seq))
			if h['pid'] == IDENTITY_DUP and cov == COVERAGE_DUP:
				continue
			if h['ev'] >= EVALUE_DB or h['pid'] < IDENTITY_DB or cov < COVERAGE_DB:
				continue
			range_str = h['contig_range'].split(':', 1)[1]
			if range_str.startswith("c(") and range_str.endswith(")"):
				start, end = map(int, range_str[2:-1].split('-'))
				is_reverse = True
			else:
				start, end = map(int, range_str.split('-'))
				is_reverse = False
			is_in, ab = inside_annotated(h['contig'], start, end, annot)
			if is_in:
				print(f"[DEBUG] Skip tBLASTn {h['contig_range']} ({start}-{end}) : inside annotated CDS {ab[0]}-{ab[1]} on {h['contig']}")
				continue

			# Extraction sens +
			if not is_reverse:
				seq_cds = seq[start-1:]
				found_stop = False
				for i in range(0, len(seq_cds)-2, 3):
					codon = str(seq_cds[i:i+3])
					if codon in ("TAA", "TAG", "TGA"):
						found_stop = True
						end_cds = start - 1 + i + 3
						prot_seq = seq[start-1:end_cds].translate(table=TRANSLATION_TABLE)
						prot_str = str(prot_seq)
						# Si pas de M au début
						if not prot_str or prot_str[0] != "M":
							offset, best_prot = find_longest_M_window_around(seq, start, window_nt=60)
							if offset is not None and best_prot and best_prot[0]=="M":
								new_start = start + offset
								stop_found = False
								seq_cds2 = seq[new_start-1:]
								for j in range(0, len(seq_cds2)-2, 3):
									codon2 = str(seq_cds2[j:j+3])
									if codon2 in ("TAA", "TAG", "TGA"):
										stop_found = True
										end_cds2 = new_start - 1 + j + 3
										header = f"{h['contig']}:{new_start}-{end_cds2} | Mwindow"
										# FILTRE TAILLE ICI
										if best_prot and best_prot not in proteins_seen and len(best_prot) > 6 and len(best_prot) < 81:
											rec2 = SeqRecord(best_prot, id=header, description="")
											extracted[header] = rec2
											new_t.append(header)
											source_map[header] = "tBLASTn"
											proteins_seen.add(best_prot)
										else:
											print(f"[INFO] Skip duplicate protein (window) or length: {header}")
										break
								if not stop_found:
									print(f"[INFO] Skip {h['contig']}:{new_start}-... : no stop codon after best M")
							else:
								print(f"[INFO] Skip {h['contig']}:{start}-{end_cds}: no M found in window")
							break
						if '*' not in prot_str:
							print(f"[INFO] Skip {h['contig']}:{start}-{end_cds}: no stop in translated sequence")
							break
						stop_index = prot_str.index('*')
						prot_str = prot_str[:stop_index]
						header = h['contig_range']
						# FILTRE TAILLE ICI
						if prot_str and prot_str not in proteins_seen and len(prot_str) > 6 and len(prot_str) < 81:
							rec2 = SeqRecord(prot_str, id=header, description=f"tBLASTn|from:{qid}")
							extracted[header] = rec2
							new_t.append(header)
							source_map[header] = "tBLASTn"
							proteins_seen.add(prot_str)
						else:
							print(f"[INFO] Skip duplicate protein: {header}")
						break
				if not found_stop:
					print(f"[INFO] Skip {h['contig']}:{start}-... : no stop codon found downstream")
			# Extraction sens - (reverse complement)
			else:
				seq_rc = seq.reverse_complement()
				genome_len = len(seq)
				rc_start = genome_len - end + 1
				seq_cds = seq_rc[rc_start-1:]
				found_stop = False
				for i in range(0, len(seq_cds)-2, 3):
					codon = str(seq_cds[i:i+3])
					if codon in ("TAA", "TAG", "TGA"):
						found_stop = True
						end_cds = rc_start - 1 + i + 3
						prot_seq = seq_rc[rc_start-1:end_cds].translate(table=TRANSLATION_TABLE)
						prot_str = str(prot_seq)
						if not prot_str or prot_str[0] != "M":
							offset, best_prot = find_longest_M_window_around(seq_rc, rc_start, window_nt=60)
							if offset is not None and best_prot and best_prot[0]=="M":
								new_rc_start = rc_start + offset
								if new_rc_start < 1 or new_rc_start > genome_len:
									break
								stop_found = False
								seq_cds2 = seq_rc[new_rc_start-1:]
								for j in range(0, len(seq_cds2)-2, 3):
									codon2 = str(seq_cds2[j:j+3])
									if codon2 in ("TAA", "TAG", "TGA"):
										stop_found = True
										end_cds2 = new_rc_start - 1 + j + 3
										new_genome_start = genome_len - (new_rc_start - 1)
										new_genome_end   = genome_len - (end_cds2 - 1) + 1
										header = f"{h['contig']}:c({min(new_genome_end,new_genome_start)}-{max(new_genome_end,new_genome_start)}) | Mwindow"
										# FILTRE TAILLE ICI
										if best_prot and best_prot not in proteins_seen and len(best_prot) > 6 and len(best_prot) < 81:
											rec2 = SeqRecord(best_prot, id=header, description="")
											extracted[header] = rec2
											new_t.append(header)
											source_map[header] = "tBLASTn"
											proteins_seen.add(best_prot)
										else:
											print(f"[INFO] Skip duplicate protein (window) or length: {header}")
										break
								if not stop_found:
									print(f"[INFO] Skip {h['contig']}:c({rc_start}-...): no stop codon after best M")
							else:
								print(f"[INFO] Skip {h['contig']}:c({start}-{end_cds}): no M found in window")
							break
						if '*' not in prot_str:
							print(f"[INFO] Skip {h['contig']}:c({start}-{end_cds}): no stop in translated sequence")
							break
						stop_index = prot_str.index('*')
						prot_str = prot_str[:stop_index]
						header = h['contig_range']
						# FILTRE TAILLE ICI
						if prot_str and prot_str not in proteins_seen and len(prot_str) > 6 and len(prot_str) < 81:
							rec2 = SeqRecord(prot_str, id=header, description=f"tBLASTn|from:{qid}")
							extracted[header] = rec2
							new_t.append(header)
							source_map[header] = "tBLASTn"
							proteins_seen.add(prot_str)
						else:
							print(f"[INFO] Skip duplicate protein: {header}")
						break
				if not found_stop:
					print(f"[INFO] Skip {h['contig']}:c({start}-{end})-... : no stop codon found downstream")
			break  # un seul hit par query

	# Ajout des prédits non trouvés dans BDD2, AVEC FILTRE TAILLE
	for qid, rec in pred.items():
		if qid not in bdd1 and qid not in bdd2 and len(str(rec.seq)) > 6 and len(str(rec.seq)) < 81:
			bdd2[qid] = rec
			added2 += 1

	# Écriture du fichier FASTA de sortie (sécurité, refiltre la taille)
	with open(out_fa, 'w') as out:
		seqs_written = set()
		for q in [x for x in (new1+new_q2+new_s2) if x in pred]:
			seqstr = str(pred[q].seq)
			if seqstr not in seqs_written :
				out.write(f">{q} | {source_map[q]}\n{seqstr}\n")
				seqs_written.add(seqstr)
		for cr in new_t:
			rec2 = extracted[cr]
			seqstr = str(rec2.seq)
			if seqstr not in seqs_written :
				out.write(f">{cr} | tBLASTn\n{seqstr}\n")
				seqs_written.add(seqstr)

	# Sauvegarde des BDD mises à jour
	SeqIO.write(bdd1.values(), open(BDD1_PATH,'w'), 'fasta')
	SeqIO.write(bdd2.values(), open(BDD2_PATH,'w'), 'fasta')

	# Affichage résumé
	total1 = len(new1)+len(new_q2)+len(new_s2)+len(new_t)
	print(f"\n[SUMMARY] Added to BDD1\t: {total1}")
	print(f"  - BLASTP  ->  BDD1 (orig)     : {len(new1)}")
	print(f"  - BLASTP  ->  BDD1 (via BDD2) : {len(new_q2)+len(new_s2)}")
	print(f"  - tBLASTn  ->  BDD1           : {len(new_t)}")
	print(f"[SUMMARY] Added to BDD2\t\t: {added2}")
	print("[Done]")

if __name__=='__main__':
	main()

