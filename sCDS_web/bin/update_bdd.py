#!/usr/bin/env python3										
import sys, os, csv, shutil, re							
from Bio import SeqIO								
from Bio.SeqRecord import SeqRecord						

# DEFINITION DES SEUILS
EVALUE_BLAST   = 1e-5								# Seuil E-value maximal pour BLASTP
IDENTITY_BLAST = 30.0								# Seuil identité (%) minimal pour BLASTP
COVERAGE_BLAST = 0.5								# Seuil couverture minimal pour BLASTP

EVALUE_DB      = 1e-10								# Seuil E-value maximal pour tBLASTn
IDENTITY_DB    = 80.0								# Seuil identité (%) minimal pour tBLASTn
COVERAGE_DB    = 0.0								# Seuil couverture minimal pour tBLASTn

IDENTITY_DUP   = 100.0								# Identité (%) pour détecter duplications exactes
COVERAGE_DUP   = 1.0								# Couverture (%) pour duplications exactes

TRANSLATION_TABLE = 1								# Table de traduction standard pour Biopython

DEFAULT_BDD1 = "data//database_ortholog/sCDS_ortholog.fasta"			# Chemin par défaut vers la base orthologue
DEFAULT_BDD2 = "data/database_not_ortholog/sCDS_not_ortholog.fasta"		# Chemin par défaut vers la base non-orthologue
BDD1_PATH    = os.environ.get("BDD1_PATH", DEFAULT_BDD1)			# Permet de surcharger le chemin de BDD1 via variable d’environnement
BDD2_PATH    = os.environ.get("BDD2_PATH", DEFAULT_BDD2)			# Permet de surcharger le chemin de BDD2 via variable d’environnement

# DEFINITION DES FONCTIONS

def normalize_id(s):										# Fonction pour normaliser un identifiant de contig
	"""Garde tout avant le premier '.', '_' ou ':' afin d’extraire le nom du contig."""
	return re.split(r"[._:]", s, 1)[0]							# Découpe la chaîne et retourne la première partie

def backup(fp):											# Fonction pour créer une sauvegarde d’un FASTA
	"""Copie fp vers fp_backup.fasta s’il existe."""
	if os.path.isfile(fp):									# Vérifie si le fichier existe
		bak = fp.replace('.fasta', '_backup.fasta')					# Construit le nom du fichier de sauvegarde
		shutil.copy(fp, bak)								# Copie le fichier original vers la sauvegarde
		print(f"[INFO] Backup created: {bak}")						# Informe l’utilisateur
	else:
		print(f"[WARN] {fp} not found, skipped backup.")				# Avertissement si le fichier n’existe pas

def load_fasta(fp):										# Fonction pour charger un fichier FASTA
	"""Charge un FASTA, retourne un dict id -> SeqRecord."""
	return { rec.id: rec for rec in SeqIO.parse(fp, 'fasta') }				# Parcours le FASTA et crée un dict

def parse_blast(tab):										# Fonction pour parser la sortie BLASTP (outfmt6)
	"""Parse BLASTP outfmt6  ->  dict qid -> liste de hits."""
	hits = {}										# Initialise le dict des hits
	with open(tab) as f:									# Ouvre le fichier BLASTP
		reader = csv.reader(f, delimiter='\t')						# Lit en tant que CSV tabulé
		for row in reader:								# Pour chaque ligne du fichier
			if len(row) < 12: continue						# Ignore si moins de 12 colonnes
			qid, sid = row[0], row[1]						# Récupère qid et sid
			pid = float(row[2]); aln = int(row[3]); ev = float(row[10])		# Récupère %identity, longueur alignement, E-value
			hits.setdefault(qid, []).append({					# Ajoute un hit à la liste pour qid
				'sid': sid,							# Champ subject id complet
				'sid_norm': normalize_id(sid),					# Champ subject id normalisé
				'pid': pid,							# % identity
				'alen': aln,							# alignment length
				'ev': ev							# E-value
			})
	return hits										

def parse_tblastn(tab):										# Fonction pour parser la sortie tBLASTn (outfmt6)
	"""Parse tBLASTn outfmt6  ->  dict qid -> liste de hits avec contig_range."""
	hits = {}										# Initialise le dict des hits
	with open(tab) as f:									# Ouvre le fichier tBLASTn
		reader = csv.reader(f, delimiter='\t')						# Lit en tant que CSV tabulé
		for row in reader:								# Pour chaque ligne
			if len(row) < 12: continue						# Ignore si moins de 12 colonnes
			qid = row[0]								# Récupère qid
			full_s = row[1]								# Récupère le subject complet
			pid = float(row[2]); aln = int(row[3]); ev = float(row[10])		# Récupère %identity, alen, E-value
			ss, se = int(row[8]), int(row[9])					# Récupère start et end sur le sujet
			start, end = min(ss,se), max(ss,se)					# Garantit start < end
			cont = normalize_id(full_s)						# Normalise l’identifiant du contig
			cr = f"{cont}:{start}-{end}"						# Crée la chaîne contig_range
			hits.setdefault(qid, []).append({					# Ajoute un hit
				'contig_range': cr,							# Intervalle contig_range
				'contig': cont,							# Nom du contig
				'pid': pid,							# % identity
				'alen': aln,							# alignment length
				'ev': ev							# E-value
			})
	return hits										# Retourne le dict des hits

def cover(aln, length):										# Fonction pour calculer la couverture d’un alignement
	"""Retourne aln/length, ou 0 si length == 0."""
	return aln/length if length else 0							# Renvoie la fraction ou 0

def load_embl_cds(embl_dir):									# Fonction pour parser les CDS annotées dans EMBL
	"""Parse tous les .embl pour extraire CDS annotées par contig."""
	annot = {}										# Initialise le dict des annotations
	for fn in os.listdir(embl_dir):								# Parcourt tous les fichiers du répertoire
		if not fn.endswith('.embl'): continue						# Ignore les fichiers non .embl
		path = os.path.join(embl_dir, fn)						# Construit le chemin complet
		for rec in SeqIO.parse(path, 'embl'):						# Parse chaque record EMBL
			cont = normalize_id(rec.id)						# Normalise l’ID du contig
			for feat in rec.features:						# Parcourt les features
				if feat.type=='CDS' and feat.location:				# Si la feature est de type CDS
					a = int(feat.location.start)+1				# Récupère la position start (+1 pour conversion 0 -> 1)
					b = int(feat.location.end)				# Récupère la position end
					annot.setdefault(cont, []).append((a,b))		# Ajoute l’intervalle annoté
	return annot										# Retourne dict contig -> [(start,end),...]

def inside_annotated(contig, s, e, annot):							# Fonction pour tester si un intervalle chevauche une annotation
	"""Retourne True si [s,e] est contenu dans une CDS annotée du contig."""
	for a,b in annot.get(contig, []):							# Parcourt les intervalles annotés pour ce contig
		if s>=a and e<=b:								# Vérifie si l’intervalle testé est inclus
			return True								# Retourne True si inclusion
	return False										# Sinon retourne False

# ——— Main ——————————————————————————————————————————————

def main():											# Fonction principale du script
	# Vérification des arguments
	if len(sys.argv)!=8:									# Si le nombre d’arguments est incorrect
		sys.exit("Usage: update_bdd.py blast1 blast2 tblastn prot_fa genome_fa out_fa embl_dir")
	b1_tab, b2_tab, tbl_tab, prot_fa, genome_fa, out_fa, embl_dir = sys.argv[1:]		# Assigne les chemins reçus

	# 0) Création de sauvegardes
	backup(BDD1_PATH)									# Backup de BDD1 avant modification
	backup(BDD2_PATH)									# Backup de BDD2 avant modification

	# 1) Chargement des données et annotations
	pred   = load_fasta(prot_fa)								# Charge les prédits sCDS
	bdd1   = load_fasta(BDD1_PATH)								# Charge la base orthologue
	bdd2   = load_fasta(BDD2_PATH)								# Charge la base non-orthologue
	genome = { normalize_id(r.id): r for r in SeqIO.parse(genome_fa,'fasta') }		# Charge le génome en dict contig -> SeqRecord
	contigs = set(genome.keys())								# Ensemble des noms de contigs du génome
	annot   = load_embl_cds(embl_dir)							# Charge les annotations EMBL

	# 1b) Si le génome est déjà présent dans BDD1, on exporte simplement ces entrées
	if contigs & { normalize_id(k) for k in bdd1 }:						# Intersection entre contigs et BDD1
		print("[INFO] Genome already in BDD1   ->  exporting existing entries")
		with open(out_fa,'w') as out:							# Ouvre le fichier de sortie
			for rid, rec in bdd1.items():						# Parcourt les enregistrements de BDD1
				if normalize_id(rid) in contigs:				# Si le contig correspond
					out.write(f">{rec.id} | BDD1\n{rec.seq}\n")		# Écrit en FASTA
		return										# Termine l’exécution

	# 2) Parsing des résultats BLAST
	blast1  = parse_blast(b1_tab)								# Résultats BLASTP vs BDD1
	blast2  = parse_blast(b2_tab)								# Résultats BLASTP vs BDD2
	tblastn = parse_tblastn(tbl_tab)							# Résultats tBLASTn vs génome

	# 3) Initialisation des listes de nouvelles séquences
	new1, new_q2, new_s2, new_t = [], [], [], []						# Listes pour stocker les IDs ajoutés
	source_map, extracted = {}, {}								# Dictionnaires pour suivi et extraits
	added2 = 0										# Compteur pour BDD2

	# 4) Boucle BLASTP pour valider les prédits
	for qid, rec in pred.items():								# Pour chaque protéine prédite
		L = len(rec.seq)								# Longueur de la protéine
		#  ->  recherche dans BDD1
		for h in blast1.get(qid, []):									# Parcourt les hits BLAST1
			cov = cover(h['alen'],L)								# Calcule la couverture
			if h['pid']==IDENTITY_DUP and cov==COVERAGE_DUP: continue				# Ignore duplications exactes
			if h['ev']<EVALUE_BLAST and h['pid']>=IDENTITY_BLAST and cov>=COVERAGE_BLAST:		# Si hit valide
				if qid not in bdd1:								# Si non déjà dans BDD1
					bdd1[qid]=rec								# Ajoute dans BDD1
					new1.append(qid)							# Enregistre l’ID
					source_map[qid]="BDD1"							# Mémorise la source
				break
													
		#  ->  recherche dans BDD2
		for h in blast2.get(qid, []):									# Parcourt les hits BLAST2
			# Ignore si provenant du même génome
			if h['sid_norm'] in contigs:								# Compare le contig normalisé
				print(f"[DEBUG] Skip BLAST2 {qid}  ->  {h['sid']} (same genome)")
				continue									# Continue au hit suivant
			cov = cover(h['alen'],L)								# Calcule la couverture
			if h['pid']==IDENTITY_DUP and cov==COVERAGE_DUP: continue				# Ignore duplications exactes
			if h['ev']<EVALUE_BLAST and h['pid']>=IDENTITY_BLAST and cov>=COVERAGE_BLAST:		# Si hit valide
				sid = h['sid']									# Récupère l’ID du sujet
				if sid in bdd2 and sid not in bdd1:						# Si sujet dans BDD2 mais pas BDD1
					bdd1[sid]=bdd2.pop(sid)							# Déplace de BDD2 vers BDD1
					new_s2.append(sid)							# Enregistre l’ID
					source_map[sid]="BDD2"							# Mémorise la source
				if qid not in bdd1:								# Si qid pas encore dans BDD1
					bdd1[qid]=rec								# Ajoute qid
					new_q2.append(qid)							# Enregistre l’ID
					source_map[qid]="BDD2"							# Mémorise la source
				break										# Sort de la boucle BLAST2

	# 5) Boucle tBLASTn pour extraire de nouvelles protéines
	for qid, hits in tblastn.items():									# Pour chaque requête tBLASTn
		for h in hits:											# Parcourt chaque hit
			if h['contig'] in contigs:								# Si hit sur même génome
				print(f"[DEBUG] Skip tBLASTn {h['contig_range']} (same genome)")
				continue									# Continue au hit suivant
			seq = genome[h['contig']].seq								# Séquence du contig
			cov = cover(h['alen'],len(seq))								# Calcule la couverture
			if h['pid']==IDENTITY_DUP and cov==COVERAGE_DUP: continue				# Ignore duplications exactes
			if h['ev']>=EVALUE_DB or h['pid']<IDENTITY_DB or cov<COVERAGE_DB: continue		# Filtre par seuils
			start,end = map(int, h['contig_range'].split(':',1)[1].split('-',1))			# Extrait les positions
			if inside_annotated(h['contig'],start,end,annot):					# Si chevauche annotation
				print(f"[INFO] Skip {h['contig_range']}: inside annotated CDS")
				continue							
				
			stop = seq[end:end+3]									# Codon stop après intervalle
			if len(stop)!=3 or stop.translate(table=TRANSLATION_TABLE)=='':				# Vérifie codon stop valide
				print(f"[ERROR] No stop at {h['contig_range']}, skip")
				continue							
				
			prot = seq[start-1:end].translate(table=TRANSLATION_TABLE)				# Traduit la région
			if not prot or prot[0]!='M':								# Vérifie la présence d’un M initial
				print(f"[INFO] Skip {h['contig_range']}: no start-M")	
				continue							
				
			if h['contig_range'] not in extracted:							# Si pas déjà extrait
				rec2 = SeqRecord(prot, id=h['contig_range'], description=f"tBLASTn|from:{qid}")
				extracted[h['contig_range']] = rec2						# Stocke le SeqRecord
				new_t.append(h['contig_range'])							# Enregistre l’intervalle
				source_map[h['contig_range']] = "tBLASTn"					# Mémorise la source
			break									

	# 6) Ajout des prédits non trouvés  ->  BDD2
	for qid, rec in pred.items():										# Pour chaque prédiction
		if qid not in bdd1 and qid not in bdd2:								# Si non dans aucune base
			bdd2[qid] = rec										# Ajoute dans BDD2
			added2 += 1										# Incrémente compteur

	# 7) Écriture du fichier FASTA de sortie
	with open(out_fa,'w') as out:										# Ouvre le fichier de sortie
		# a) Proteines validées par BLASTP (uniquement issues de pred)
		for q in [x for x in (new1+new_q2+new_s2) if x in pred]:					# Filtre pour garder q dans pred
			out.write(f">{q} | {source_map[q]}\n{pred[q].seq}\n")					# Écrit le FASTA
		# b) Proteines extraites par tBLASTn
		for cr in new_t:										# Pour chaque contig_range extrait
			rec2 = extracted[cr]					# Récupère le SeqRecord
			out.write(f">{cr} | tBLASTn\n{rec2.seq}\n")		# Écrit le FASTA

	# 8) Sauvegarde des bases mises à jour
	SeqIO.write(bdd1.values(), open(BDD1_PATH,'w'), 'fasta')		# Écrit BDD1 mise à jour
	SeqIO.write(bdd2.values(), open(BDD2_PATH,'w'), 'fasta')		# Écrit BDD2 mise à jour

	# 9) Affichage du résumé
	total1 = len(new1)+len(new_q2)+len(new_s2)+len(new_t)			# Nombre total ajouté à BDD1
	print(f"\n[SUMMARY] Added to BDD1	: {total1}")
	print(f"  - BLASTP  ->  BDD1 (orig)     : {len(new1)}")			# Détails BLASTP originaire
	print(f"  - BLASTP  ->  BDD1 (via BDD2) : {len(new_q2)+len(new_s2)}")	# Détails BLASTP via BDD2
	print(f"  - tBLASTn  ->  BDD1           : {len(new_t)}")		# Détails tBLASTn
	print(f"[SUMMARY] Added to BDD2		: {added2}")			# Nombre ajouté à BDD2
	print("[Done]")								

if __name__=='__main__':									
	main()									

