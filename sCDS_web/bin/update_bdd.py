#!/usr/bin/env python3										
import sys, os, csv, shutil, re								
from Bio import SeqIO									
from Bio.SeqRecord import SeqRecord							

# DEFINITION DES SEUILS
EVALUE_BLAST   = 1e-5										# seuil E-value pour BLASTP
IDENTITY_BLAST = 30.0										# seuil identité % pour BLASTP
COVERAGE_BLAST = 0.5										# seuil couverture alignement pour BLASTP

EVALUE_DB      = 1e-10										# seuil E-value pour tBLASTn
IDENTITY_DB    = 80.0										# seuil identité % pour tBLASTn
COVERAGE_DB    = 0.0										# seuil couverture pour tBLASTn

IDENTITY_DUP   = 100.0										# identité % pour considérer duplication exacte
COVERAGE_DUP   = 1.0										# couverture pour duplication exacte

TRANSLATION_TABLE = 1										# table de traduction standard

DEFAULT_BDD1 = "data//database_ortholog/sCDS_ortholog.fasta"					# chemin par défaut BDD1
DEFAULT_BDD2 = "data/database_not_ortholog/sCDS_not_ortholog.fasta"				# chemin par défaut BDD2
BDD1_PATH    = os.environ.get("BDD1_PATH", DEFAULT_BDD1)					# chemin BDD1 paramétrable
BDD2_PATH    = os.environ.get("BDD2_PATH", DEFAULT_BDD2)					# chemin BDD2 paramétrable

# DEFINITION DES FONCTIONS 

def normalize_id(s):										# normalise un identifiant de contig
	"""Garde tout avant premier '.', '_' ou ':' pour extraire le contig."""
	return re.split(r"[._:]", s, 1)[0]							# découpe et renvoie la première partie

def backup(fp):											# crée un backup d'un FASTA existant
	"""Crée un backup du FASTA s’il existe."""
	if os.path.isfile(fp):									# si le fichier existe
		bak = fp.replace('.fasta', '_backup.fasta')					# nom du backup
		shutil.copy(fp, bak)								# copie
		print(f"[INFO] Backup created: {bak}")						
	else:
		print(f"[WARN] {fp} not found, skipped backup.")				# avertissement si absent

def load_fasta(fp):										# charge un FASTA dans un dict
	"""Charge un FASTA, clef = rec.id."""
	return { rec.id: rec for rec in SeqIO.parse(fp, 'fasta') }				# comprehension dict id -> SeqRecord

def parse_blast(tab):										# parse un BLASTP outfmt6
	"""Parse BLAST tabulé (outfmt 6)."""
	hits = {}										# dict pour stocker hits
	with open(tab) as f:									# ouvre le fichier
		reader = csv.reader(f, delimiter='\t')						# CSV tabulé
		for row in reader:								# pour chaque ligne
			if len(row) < 12: continue						# on s’assure du nombre de colonnes
			qid, sid = row[0], row[1]						# requête et subject
			pid = float(row[2]); aln = int(row[3]); ev = float(row[10])		# identity, alignment length, E-value
			hits.setdefault(qid, []).append({					# ajoute dans la liste pour qid
				'sid': sid,
				'sid_norm': normalize_id(sid),
				'pid': pid, 'alen': aln, 'ev': ev
			})
	return hits										# renvoie dict qid -> list(hits)

def parse_tblastn(tab):										# parse un tBLASTn outfmt6
	"""Parse TBLASTN tabulé, génère contig_range et contig."""
	hits = {}										# dict pour stocker hits
	with open(tab) as f:									# ouvre le fichier
		reader = csv.reader(f, delimiter='\t')						# CSV tabulé
		for row in reader:								# pour chaque ligne
			if len(row) < 12: continue						# vérification colonnes
			qid = row[0]								# requête
			full_s = row[1]								# sujet complet (avec positions)
			pid = float(row[2]); aln = int(row[3]); ev = float(row[10])		# identity, alen, E-value
			ss, se = int(row[8]), int(row[9])					# start and end positions subject
			start, end = min(ss, se), max(ss, se)					# ordonne start < end
			cont = normalize_id(full_s)						# identifiant normalisé du contig
			cr = f"{cont}:{start}-{end}"						# assemble contig_range
			hits.setdefault(qid, []).append({					# ajoute hit
				'contig_range': cr,
				'contig': cont,
				'pid': pid,
				'alen': aln,
				'ev': ev
			})
	return hits										# renvoie dict qid -> list(hits)

def cover(aln, length):										# calcule couverture
	"""Retourne la couverture."""
	return aln / length if length else 0							# aln/length ou 0 si length=0

def load_embl_cds(embl_dir):									# parse tous les .embl pour leurs CDS
	"""Charge toutes les CDS annotées (start,end) par contig normalisé."""
	annot = {}										# dict contig -> list((start,end))
	for fn in os.listdir(embl_dir):								# pour chaque fichier dans le répertoire
		if not fn.endswith('.embl'): continue						# ignore les autres
		path = os.path.join(embl_dir, fn)						# chemin complet
		for rec in SeqIO.parse(path, 'embl'):						# parse EMBL
			cont = normalize_id(rec.id)						# ID de contig normalisé
			for feat in rec.features:						# pour chaque feature
				if feat.type == 'CDS' and feat.location:			# si c’est une CDS
					a = int(feat.location.start) + 1			# coord start +1
					b = int(feat.location.end)				# coord end
					annot.setdefault(cont, []).append((a, b))		# ajoute intervalle
	return annot										# renvoie dict contig -> list(intervals)

def inside_annotated(contig, s, e, annot):							# teste inclusion dans CDS annotée
	"""Vrai si [s,e] est inclus dans une CDS annotée du contig."""
	for a, b in annot.get(contig, []):							# pour chaque intervalle annoté
		if s >= a and e <= b:								# si entièrement contenu
			return True								# retourne True
	return False										# sinon False

# FONCTION PRINCIPALE DU SCRIPT 

def main():												# fonction principale
	if len(sys.argv) != 8:										# si nombre d’arguments incorrect
		sys.exit("Usage: update_bdd.py blast1 blast2 tblastn prot_fa genome_fa out_fa embl_dir")
	b1_tab, b2_tab, tbl_tab, prot_fa, genome_fa, out_fa, embl_dir = sys.argv[1:]			# assignation

	# 0) backups
	backup(BDD1_PATH)										# backup BDD1
	backup(BDD2_PATH)										# backup BDD2

	# 1) charger FASTA et annotations
	pred   = load_fasta(prot_fa)									# sCDS prédits
	bdd1   = load_fasta(BDD1_PATH)									# base orthologue
	bdd2   = load_fasta(BDD2_PATH)									# base non-orthologue
	genome = { normalize_id(r.id): r for r in SeqIO.parse(genome_fa, 'fasta') }			# dict contig -> SeqRecord
	contigs = set(genome.keys())									# ensemble des contigs du génome
	annot   = load_embl_cds(embl_dir)								# annotations EMBL

	# 1b) si ce génome est déjà dans BDD1  -> exporter existants
	if contigs & { normalize_id(k) for k in bdd1 }:		
		print("[INFO] Genome already in BDD1  ->  exporting existing entries")
		with open(out_fa, 'w') as out:					
			for rid, rec in bdd1.items():				
				if normalize_id(rid) in contigs:		
					out.write(f">{rec.id} | BDD1\n{rec.seq}\n")			# écrit dans out_fa
		return												

	# 2) parser BLAST
	blast1  = parse_blast(b1_tab)									# BLASTP vs BDD1
	blast2  = parse_blast(b2_tab)									# BLASTP vs BDD2
	tblastn = parse_tblastn(tbl_tab)								# tBLASTn vs génome

	new1, new_q2, new_s2, new_t = [], [], [], []							# listes de nouvelles séquences
	source_map = {}											# pour traçabilité
	extracted  = {}											# pour tBLASTn extraits
	added2     = 0											# compteur BDD2

	# 3) BLASTP vs BDD1 & BDD2
	for qid, rec in pred.items():									# pour chaque sCDS prédit
		L = len(rec.seq)									# longueur protéine
		# -> BDD1
		for h in blast1.get(qid, []):								# hits BLAST1
			cov = cover(h['alen'], L)						
			if h['pid'] == IDENTITY_DUP and cov == COVERAGE_DUP: continue			# duplication exacte
			if h['ev'] < EVALUE_BLAST and h['pid'] >= IDENTITY_BLAST and cov >= COVERAGE_BLAST:
				if qid not in bdd1:						
					bdd1[qid] = rec						
					new1.append(qid)						
					source_map[qid] = "BDD1"				
				break										# stop après premier match
		# -> BDD2
		#  ->  Boucle sur les hits BLAST2 pour qid (comparaison contre BDD2)
		for h in blast2.get(qid, []):												
			# Si le contig normalisé du hit est dans le même génome, on l’ignore
			if h['contig'] in contigs:											
				print(f"[DEBUG] Skip BLAST2 {qid} -> {h['sid']} (same genome)")		
				continue														
			# Calcul de la couverture de l’alignement
			cov = cover(h['alen'], L)											
			# Si identité et couverture correspondent à une duplication exacte, on skip
			if h['pid'] == IDENTITY_DUP and cov == COVERAGE_DUP:					
				continue														
			# Si le hit passe les seuils E-value, identité et couverture
			if h['ev'] < EVALUE_BLAST and h['pid'] >= IDENTITY_BLAST and cov >= COVERAGE_BLAST:
				# On récupère l’ID du sujet
				sid = h['sid']													
				# Si ce sujet existe dans BDD2 mais pas encore dans BDD1, on le migre
				if sid in bdd2 and sid not in bdd1:								
					bdd1[sid] = bdd2.pop(sid)									
					new_s2.append(sid)											
					source_map[sid] = "BDD2"									
				# On ajoute la requête qid dans BDD1 si elle n’y est pas déjà
				if qid not in bdd1:												
					bdd1[qid] = rec											
					new_q2.append(qid)											
					source_map[qid] = "BDD2"									
				# On sort de la boucle BLAST2 après avoir trouvé un hit valide
				break															

	# 4) TBLASTn  ->  Parcours des résultats tBLASTn pour chaque requête
	for qid, hits in tblastn.items():											
		for h in hits:															
			# Ignore les hits issus du même génome
			if h['contig'] in contigs:											
				print(f"[DEBUG] Skip tBLASTn {h['contig_range']} (same genome)")	
				continue														
			# On récupère la séquence du contig
			seq = genome[h['contig']].seq											
			# Calcul de la couverture de l’alignement sur le génome
			cov = cover(h['alen'], len(seq))										
			# Skip si duplication exacte
			if h['pid'] == IDENTITY_DUP and cov == COVERAGE_DUP:						
				continue														
			# Skip si E-value, identité ou couverture hors seuil
			if h['ev'] >= EVALUE_DB or h['pid'] < IDENTITY_DB or cov < COVERAGE_DB:	
				continue														
			# Extraction des positions start/end depuis contig_range
			start, end = map(int, h['contig_range'].split(':',1)[1].split('-',1))
			# Skip si l’intervalle est déjà annoté comme CDS
			if inside_annotated(h['contig'], start, end, annot):					
				print(f"[INFO] Skip {h['contig_range']}: inside annotated CDS")	
				continue														
			# Vérifie la présence d’un codon stop après la région
			stop = seq[end:end+3]													
			if len(stop) != 3 or stop.translate(table=TRANSLATION_TABLE) == '':		
				print(f"[ERROR] No stop at {h['contig_range']}, skip")		
				continue														
			# Traduit la région pour obtenir la protéine
			prot = seq[start-1:end].translate(table=TRANSLATION_TABLE)			
			# Skip si pas de M au début
			if not prot or prot[0] != 'M':										
				print(f"[INFO] Skip {h['contig_range']}: no start-M")		
				continue														
			# Si on n’a pas encore extrait ce contig_range, on l’ajoute
			if h['contig_range'] not in extracted:								
				rec2 = SeqRecord(prot, id=h['contig_range'], description=f"tBLASTn|from:{qid}")	
				extracted[h['contig_range']] = rec2							
				new_t.append(h['contig_range'])								
				source_map[h['contig_range']] = "tBLASTn"						
			# On sort de la boucle des hits pour cette requête
			break															

	# 5) Ajout des prédits restants dans BDD2
	for qid, rec in pred.items():											
		if qid not in bdd1 and qid not in bdd2:								
			bdd2[qid] = rec													
			added2 += 1													

	# 6) Écriture du FASTA de sortie (out_fa)
	seen = set()															
	with open(out_fa, 'w') as out:											
		# a) Hits validés BLASTP (BDD1 & BDD2)
		for q in new1 + new_q2 + new_s2:										
			if q in seen: continue											
			seen.add(q)														
			out.write(f">{q} | {source_map[q]}\n{pred[q].seq}\n")			
		# b) Séquences extraites par tBLASTn
		for cr, rec2 in extracted.items():									
			if cr in seen: continue											
			seen.add(cr)														
			out.write(f">{cr} | tBLASTn\n{rec2.seq}\n")						

	# 7) sauver BDD1 & BDD2
	SeqIO.write(bdd1.values(), open(BDD1_PATH, 'w'), 'fasta')	# mise à jour BDD1
	SeqIO.write(bdd2.values(), open(BDD2_PATH, 'w'), 'fasta')	# mise à jour BDD2

	# 8) résumé
	total1 = len(new1) + len(new_q2) + len(new_s2) + len(new_t)	# total ajoutés BDD1
	print(f"\n[SUMMARY] Added to BDD1: {total1}")					
	print(f"  - BLASTP  ->  BDD1 (orig)     : {len(new1)}")				
	print(f"  - BLASTP  ->  BDD1 (via BDD2) : {len(new_q2)+len(new_s2)}")
	print(f"  - tBLASTn  ->  BDD1          : {len(new_t)}")
	print(f"[SUMMARY] Added to BDD2: {added2}")						
	print("[Done]")											

if __name__=='__main__':											
	main()													# exécution du main si script appelé directement

