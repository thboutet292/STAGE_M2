#!/bin/bash

# Vérifier les arguments
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <nom_microsporidie> <genome_fasta> <valeur_I>"
    echo "Exemple: $0 E_cuniculi chemin/vers/E_cuniculi.fa 1.0"
    exit 1
fi

# Paramètres d'entrée
SPECIES="$1"
GENOME_FILE="$2"
I_VALUE="$3"

# Définition des chemins
GENOME_DIR="${GENOME_FILE}"
PROT_DIR="genome_complet/genome_prot"
FRAG_DIR="${PROT_DIR}/${SPECIES}_I_${I_VALUE}_frag"
HMM_RESULTS_DIR="orthofinder/hmm_results/${SPECIES}_I_${I_VALUE}"
ORTHOGROUP_MATCH_DIR="orthofinder/orthogroup_matches/${SPECIES}_I_${I_VALUE}"
TRAINING_FILE="data_training/CDS/CDS_${SPECIES}_I_${I_VALUE}"
GLIMMER_ICM_FILE="glimmer/icm_file/icm_file_${SPECIES}_I_${I_VALUE}"
PRODIGAL_TRAINING_FILE="prodigal/training_file/data_training_${SPECIES}_I_${I_VALUE}"
PROTEOME_REF="Proteomes/Proteomes_${SPECIES}.txt"

# Définition des répertoires de sortie
# Fichiers des outils 
GLIMMER_DIR="glimmer/${SPECIES}/${SPECIES}_I_${I_VALUE}"
PRODIGAL_DIR="prodigal/${SPECIES}/${SPECIES}_I_${I_VALUE}"
AUGUSTUS_DIR="augustus/${SPECIES}/${SPECIES}_I_${I_VALUE}"
FUNANNOTATE_DIR="funannotate/${SPECIES}/${SPECIES}_I_${I_VALUE}"

#Fichiers contenant les résultats des cd-hit-2d
CLUSTER_GLM_DIR="$GLIMMER_DIR/cluster"
CLUSTER_PROD_DIR="$PRODIGAL_DIR/cluster"
CLUSTER_AUG_DIR="$AUGUSTUS_DIR/cluster"
CLUSTER_FUN_DIR="$FUNANNOTATE_DIR/cluster"

#Fichiers pour stocker les résultats finaux
RESULT_DIR="result/${SPECIES}"

mkdir -p "$CLUSTER_GLM_DIR" "$CLUSTER_PROD_DIR" "$CLUSTER_AUG_DIR" "$CLUSTER_FUN_DIR"
mkdir -p "$RESULT_DIR"

# Création des répertoires nécessaires
mkdir -p "$PROT_DIR" "$FRAG_DIR" "$HMM_RESULTS_DIR" "$ORTHOGROUP_MATCH_DIR"
mkdir -p "$GLIMMER_DIR" "$PRODIGAL_DIR" "$AUGUSTUS_DIR" "$FUNANNOTATE_DIR"

# Vérifier si la microsporidie est dans la base de données
DATABASE_DIR="orthofinder/database/Proteomes"
if [ -f "$DATABASE_DIR/Proteomes_${SPECIES}.fasta" ]; then
    echo "$SPECIES est déjà dans la base de données, exclusion en cours..."
    DATABASE_NO_SPECIES="orthofinder/database/Proteomes_no_${SPECIES}"
    mkdir -p "$DATABASE_NO_SPECIES"
    for proteome in "$DATABASE_DIR"/*.fasta; do
        if [[ "$proteome" != *"${SPECIES}.fasta" ]]; then
            cp "$proteome" "$DATABASE_NO_SPECIES/"
        fi
    done
    DATABASE_PATH="$DATABASE_NO_SPECIES"
else
    echo "$SPECIES n'est pas dans la base de données, utilisation de la base complète."
    DATABASE_PATH="$DATABASE_DIR"
fi

# Lancer OrthoFinder
OUTPUT_DIR="orthofinder/results_${SPECIES}_I_${I_VALUE}"
if [ ! -d "$OUTPUT_DIR" ]; then
    echo "Lancement d'OrthoFinder avec -I $I_VALUE..."
    orthofinder -f "$DATABASE_PATH" -t 4 -a 4 -I "$I_VALUE" -o "$OUTPUT_DIR"
else
    echo "OrthoFinder a déjà été exécuté."
fi


# Analyse du core ou pangenome
python scripts/script_accessory.py "$SPECIES" "$I_VALUE"


# Alignement des orthogroupes
ALIGNED_ORTHOGROUP_DIR="orthofinder/aligned_orthogroups/${SPECIES}_I_$I_VALUE"
mkdir -p "$ALIGNED_ORTHOGROUP_DIR"
if [ ! -f "$ALIGNED_ORTHOGROUP_DIR/alignments_done.flag" ]; then
    for file in "$OUTPUT_DIR/pangenome"/*.fa; do
        orthogroup_name=$(basename "$file" .fa)
        mafft --auto "$file" > "$ALIGNED_ORTHOGROUP_DIR/${orthogroup_name}_aligned.fa"
    done
    touch "$ALIGNED_ORTHOGROUP_DIR/alignments_done.flag"
fi



# Génération des modèles HMM
HMM_MODELS_DIR="orthofinder/hmm_models/${SPECIES}_I_$I_VALUE"
mkdir -p "$HMM_MODELS_DIR"
if [ ! -f "$HMM_MODELS_DIR/hmmbuild_done.flag" ]; then
    for file in "$ALIGNED_ORTHOGROUP_DIR"/*.fa; do
        orthogroup_name=$(basename "$file" _aligned.fa)
        hmmbuild "$HMM_MODELS_DIR/${orthogroup_name}.hmm" "$file"
    done
    touch "$HMM_MODELS_DIR/hmmbuild_done.flag"
fi

# Génération du brin reverse
REV_GENOME_DIR="genome_complet/rev"
mkdir -p "$REV_GENOME_DIR"

seqkit seq -r -p -t dna "$GENOME_DIR" > "$REV_GENOME_DIR/${SPECIES}_rev.fna"


# Traduction en 6 cadres et fragmentation du génome
transeq -sequence "$GENOME_DIR" -outseq "$PROT_DIR/${SPECIES}_forward.fna" -frame 6 -clean -trim
transeq -sequence "$REV_GENOME_DIR/${SPECIES}_rev.fna" -outseq "$PROT_DIR/${SPECIES}_reverse.fna" -frame 6 -clean -trim



splitter -sequence "$PROT_DIR/${SPECIES}_forward.fna" -size 100000 -overlap 3000 -outseq "${FRAG_DIR}.fna"
splitter -sequence "$PROT_DIR/${SPECIES}_reverse.fna" -size 100000 -overlap 3000 -outseq "${FRAG_DIR}_rev.fna"



# Exécution de HMMER sur les fragments des deux brins
mkdir -p "$HMM_RESULTS_DIR/forward" "$HMM_RESULTS_DIR/reverse"
for hmm in "$HMM_MODELS_DIR"/*.hmm; do
    orthogroup_name=$(basename "$hmm" .hmm)
    hmmsearch --domtblout "$HMM_RESULTS_DIR/forward/${orthogroup_name}_hits.txt" "$hmm" "$FRAG_DIR.fna"
    hmmsearch --domtblout "$HMM_RESULTS_DIR/reverse/${orthogroup_name}_hits.txt" "$hmm" "${FRAG_DIR}_rev.fna"
done


mkdir -p "$ORTHOGROUP_MATCH_DIR"
echo -e "Gène\tDébut\tFin\tStrand" > "$ORTHOGROUP_MATCH_DIR/resultats_hmmsearch.tsv"

for file in $HMM_RESULTS_DIR/forward/*_hits.txt; do
    [ -e "$file" ] || continue  # Vérifie si le fichier existe
    awk 'BEGIN {OFS="\t"} $1 !~ /^#/ && $5 < 1e-5 && $4 > 50 {print $1, $18, $19, "sens"}' "$file" >> "$ORTHOGROUP_MATCH_DIR/resultats_hmmsearch.tsv"
done

for file in $HMM_RESULTS_DIR/reverse/*_hits.txt; do
    [ -e "$file" ] || continue  # Vérifie si le fichier existe
    awk 'BEGIN {OFS="\t"} $1 !~ /^#/ && $5 < 1e-5 && $4 > 50 {print $1, $18, $19 , "antisens"}' "$file" >> "$ORTHOGROUP_MATCH_DIR/resultats_hmmsearch.tsv"
done



python scripts/script_extract_orthologues.py "$SPECIES" "$I_VALUE"

bash scripts/script_trad_CDS.sh "data_training/PROT/PROT_${SPECIES}_I_${I_VALUE}" "$GENOME_DIR" "$SPECIES" "$I_VALUE"


# Faire le fichier d'entrainement de GLIMMER
build-icm -r "$GLIMMER_ICM_FILE" < "$TRAINING_FILE"

# Run Glimmer
glimmer3 -g 240 --start_codons atg "$GENOME_DIR" "$GLIMMER_ICM_FILE" "$GLIMMER_DIR/result_${SPECIES}_I_$I_VALUE"
python glimmer/script_gff.py "$SPECIES" "$I_VALUE"
    
# Fichier d'entrainement de PRODIGAL 
prodigal -i "$TRAINING_FILE" -t "$PRODIGAL_TRAINING_FILE" -p single

# Run PRODIGAL
prodigal -i "$GENOME_DIR" -t "$PRODIGAL_TRAINING_FILE" -f gff > "$PRODIGAL_DIR/result_${SPECIES}_I_$I_VALUE"


# Extraction des CDS et protéines
python scripts/script_gene_core.py "$SPECIES" "$I_VALUE"

# Exécution de CD-HIT-2D pour tous les outils
for tool in glimmer prodigal; do
    CLUSTER_DIR="$(eval echo \$${tool^^}_DIR)/cluster"
    PROT_FILE="$(eval echo \$${tool^^}_DIR)/${SPECIES}_${I_VALUE}_PROT_${tool}.fasta"
    
    cd-hit-2d -i "$PROTEOME_REF" -i2 "$PROT_FILE" -d 0 -o "$CLUSTER_DIR/cluster_${SPECIES}_${tool}" -c 1 -A 1
    cd-hit-2d -i "$PROT_FILE" -i2 "$PROTEOME_REF" -d 0 -o "$CLUSTER_DIR/cluster_sup_${SPECIES}_${tool}" -c 1 -A 1
done

# Analyse des clusters
python scripts/script_cluster_core.py "$SPECIES" "$I_VALUE"

echo "Pipeline terminé pour $SPECIES $I_VALUE !"



