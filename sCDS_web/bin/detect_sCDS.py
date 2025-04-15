#!/usr/bin/env python3
import sys
import os
import re
import tempfile
from flask import Flask, render_template, request, redirect, url_for, send_file
from werkzeug.utils import secure_filename
from Bio import SeqIO
from Bio.Seq import Seq

# Détermine le chemin absolu du dossier templates (qui est dans le répertoire parent de bin)
template_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "templates")

# Crée l'application Flask en spécifiant le dossier des templates
app = Flask(__name__, template_folder=template_dir)

# Utilisation du dossier temporaire pour stocker les uploads
app.config["UPLOAD_FOLDER"] = tempfile.gettempdir()

###############################################################################
# Fonctions de traitement (adaptées de votre script)

def get_short_id(s):
	"""
	Extrait la partie courte de l'identifiant d'un record.
	Exemples :
	  - Pour un header FASTA comme ">CM017447.1 Antonospora locustae ...", retourne "CM017447.1"
	  - Pour un record EMBL comme "CM017455.1ANTONOSPOR", retourne "CM017455.1"
	"""
	m = re.search(r'(CM\d+\.\d+)', s)
	if m:
		return m.group(1)
	else:
		return s.split()[0]

def contains_polyA_signal(window):
	allowed_signals = ["AATAAA", "AATTAAA"]
	for motif in allowed_signals:
		motif_len = len(motif)
		for i in range(len(window) - motif_len + 1):
			candidate = window[i:i+motif_len]
			mismatches = sum(1 for a, b in zip(candidate, motif) if a != b)
			if mismatches <= 1:
				return candidate, i, motif
	return None, None, None

def find_orfs(seq, strand="+", min_length=80):
	start_codon = "ATG"
	stop_codons = {"TAA", "TAG", "TGA"}
	orfs = []
	for frame in range(3):
		i = frame
		while i <= len(seq) - 3:
			codon = seq[i:i+3]
			if codon == start_codon:
				for j in range(i+3, len(seq)-2, 3):
					curr_codon = seq[j:j+3]
					if curr_codon in stop_codons:
						orf_length = j + 3 - i
						if orf_length >= min_length:
							if i >= 20:
								upstream = seq[i-20:i]
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
	results = []
	for record in SeqIO.parse(input_fasta, "fasta"):
		short_id = get_short_id(record.id)
		seq = str(record.seq).upper()
		record_results = {"id": short_id, "orfs": []}
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

###############################################################################
# Fonction de traitement des fichiers uploadés
def process_inputs(fasta_path, embl_dir, output_cds_path, output_prot_path):
	annotated_cdss = parse_embl_cds_from_directory(embl_dir)
	orfs_results = scan_genome(fasta_path)
	for record in orfs_results:
		rec_id = record["id"]
		if rec_id in annotated_cdss:
			cds_intervals = annotated_cdss[rec_id]
			new_orfs = []
			for orf in record["orfs"]:
				overlap = False
				for cds in cds_intervals:
					cds_start, cds_end, _ = cds
					if orf["start"] < cds_end and orf["end"] > cds_start:
						# N'affiche le message que pour les petits gènes (< 244 nt)
						if (orf["end"] - orf["start"]) < 244:
							print(f"Filtrage : sCDS {rec_id} {orf['start']+1}-{orf['end']} filtré (chevauche CDS annotée {cds_start+1}-{cds_end}).")
						overlap = True
						break
				if not overlap:
					new_orfs.append(orf)
			record["orfs"] = new_orfs
	write_output(orfs_results, output_cds_path, output_prot_path)

###############################################################################
# Routes de l'application Web

@app.route("/", methods=["GET", "POST"])
def index():
	if request.method == "POST":
		# Récupère le fichier FASTA
		fasta_file = request.files.get("fasta")
		# Récupère les fichiers EMBL (sélection multiple)
		embl_files = request.files.getlist("embl_files")
		# Récupère les noms des fichiers de sortie
		output_cds = request.form.get("output_cds")
		output_prot = request.form.get("output_prot")
		
		# Sauvegarde le fichier FASTA
		fasta_filename = secure_filename(fasta_file.filename)
		fasta_path = os.path.join(app.config["UPLOAD_FOLDER"], fasta_filename)
		fasta_file.save(fasta_path)
		
		# Crée un dossier temporaire pour stocker les fichiers EMBL
		embl_dir = os.path.join(app.config["UPLOAD_FOLDER"], "embl_upload")
		if not os.path.exists(embl_dir):
			os.makedirs(embl_dir)
		for file in embl_files:
			filename = secure_filename(file.filename)
			file.save(os.path.join(embl_dir, filename))
		
		# Définit les chemins pour les fichiers de sortie
		output_cds_path = os.path.join(app.config["UPLOAD_FOLDER"], secure_filename(output_cds))
		output_prot_path = os.path.join(app.config["UPLOAD_FOLDER"], secure_filename(output_prot))
		
		# Lance le traitement
		process_inputs(fasta_path, embl_dir, output_cds_path, output_prot_path)
		
		# Redirige vers la page de téléchargement
		return redirect(url_for("download", cds_file=output_cds, prot_file=output_prot))
	return render_template("index.html")

@app.route("/download")
def download():
	cds_file = request.args.get("cds_file")
	prot_file = request.args.get("prot_file")
	return f"""
		<h1>Résultats générés</h1>
		<ul>
			<li><a href="/files/{cds_file}">Télécharger le fichier CDS</a></li>
			<li><a href="/files/{prot_file}">Télécharger le fichier protéines</a></li>
		</ul>
	"""

@app.route("/files/<filename>")
def files(filename):
	filepath = os.path.join(app.config["UPLOAD_FOLDER"], filename)
	return send_file(filepath, as_attachment=True)

if __name__ == "__main__":
	app.run(debug=True, port=5001)
