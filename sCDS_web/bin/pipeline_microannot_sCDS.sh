#!/usr/bin/env bash
#
# pipeline_microannot_sCDS.sh : complete sCDS detection  ->  validation  ->  update
#
# Usage:
#   bash pipeline_microannot_sCDS.sh <genome_fasta> <embl_directory> <output_fasta> [--use-updated]
#
#   --use-updated : reuse the last‐updated DBs in data/ instead of copying fresh ones

export PATH="/home/thomas/anaconda3/envs/solo/bin:$PATH"

set -euo pipefail

if [ "$#" -lt 3 ] || [ "$#" -gt 4 ]; then
    echo "Usage: $0 <genome_fasta> <embl_directory> <output_fasta> [--use-updated]"
    exit 1
fi

genome_fasta="$1"
embl_directory="$2"
output_fasta="$3"
use_updated=false
if [ "${4:-}" = "--use-updated" ]; then
    use_updated=true
fi

echo "Genome FASTA:     $genome_fasta"
echo "EMBL directory:   $embl_directory"
echo "Output FASTA:     $output_fasta"
echo "Reuse updated DBs? $use_updated"
echo

# 2) Prepare working DBs
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
DATA_DIR="$SCRIPT_DIR/../data"
WORKDIR="$(pwd)"
WORK_BDD1="$WORKDIR/data/work_sCDS_ortholog.fasta"
WORK_BDD2="$WORKDIR/data/work_sCDS_not_ortholog.fasta"
mkdir -p "$WORKDIR/data"

if [ "$use_updated" = false ]; then
    echo "Step 0: Initializing fresh DBs…"
    cp "$DATA_DIR/database_ortholog/sCDS_ortholog.fasta"   "$WORK_BDD1"
    cp "$DATA_DIR/database_not_ortholog/sCDS_not_ortholog.fasta" "$WORK_BDD2"
else
    echo "Step 0: Reusing existing work DBs."
fi

# 3) Detect initial sCDS
echo "Step 1: sCDS detection…"
python3 "$SCRIPT_DIR/script_petits_genes_microannot.py" \
    "$genome_fasta" "$embl_directory" \
    "temp_CDS.fasta" "temp_prot.fasta"

# 4) Build BLAST DBs
echo "Step 2: makeblastdb BDD1 & BDD2…"
makeblastdb -in "$WORK_BDD1" -dbtype prot -parse_seqids -out "$WORKDIR/data/BDD1_db"
makeblastdb -in "$WORK_BDD2" -dbtype prot -parse_seqids -out "$WORKDIR/data/BDD2_db"

# 5) BLASTP against both
echo "Step 3: BLASTP vs BDD1…"
blastp -query "temp_prot.fasta" -db "$WORKDIR/data/BDD1_db" \
       -outfmt 6 -evalue 1e-5 -num_threads 6 -out "result/blast_BDD1.txt"
echo "Step 4: BLASTP vs BDD2…"
blastp -query "temp_prot.fasta" -db "$WORKDIR/data/BDD2_db" \
       -outfmt 6 -evalue 1e-5 -num_threads 6 -out "result/blast_BDD2.txt"

# 6) makeblastdb genome
echo "Step 5: makeblastdb genome…"
makeblastdb -in "$genome_fasta" -dbtype nucl -parse_seqids -out "$WORKDIR/data/genome_db"

# 7) TBLASTN merged DB vs genome
echo "Step 6: TBLASTN merged DB…"
cat "$WORK_BDD1" "$WORK_BDD2" > "$WORKDIR/data/merged_BDD.fasta"
tblastn -query "$WORKDIR/data/merged_BDD.fasta" -db "$WORKDIR/data/genome_db" \
       -outfmt 6 -evalue 1e-5 -num_threads 6 -out "result/tblastn_db.txt"

# 8) Run update_bdd.py  ->  writes directly to $output_fasta
echo "Step 7: Updating DB (and producing $output_fasta)…"
export BDD1_PATH="$WORK_BDD1" BDD2_PATH="$WORK_BDD2"
python3 "$SCRIPT_DIR/update_bdd.py" \
    "result/blast_BDD1.txt" \
    "result/blast_BDD2.txt" \
    "result/tblastn_db.txt" \
    "temp_prot.fasta" \
    "$genome_fasta" \
    "$output_fasta" \
    "$embl_directory"

echo "[PIPELINE] Done : $output_fasta"


# 9) CLEANUP temporary & BLAST‐DB files
rm -f temp_CDS.fasta temp_prot.fasta
rm -f "$WORKDIR/data/merged_BDD.fasta"
rm -f "$WORKDIR/data/BDD1_db".*
rm -f "$WORKDIR/data/BDD2_db".*
rm -f "$WORKDIR/data/genome_db".*


