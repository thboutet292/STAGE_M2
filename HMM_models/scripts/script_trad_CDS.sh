#!/bin/bash

# Vérifier que quatre arguments sont passés
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <proteins.fasta> <genome.fna> <name> <i_value>"
    exit 1
fi

PROTEINS=$1
GENOME=$2
NAME=$3
I_VALUE=$4

DB_NAME="blast_genome_db"
BLAST_OUTPUT="blast_results.txt"
BED_FILE="cds_regions.bed"
CDS_FORWARD="cds_forward.fasta"
CDS_REVERSE="cds_reverse.fasta"
FILTERED_CDS="data_training/CDS/CDS_${NAME}_I_${I_VALUE}"

# Étape 1: Créer la base de données BLAST
makeblastdb -in $GENOME -dbtype nucl -out $DB_NAME

# Étape 2: Exécuter TBLASTN pour aligner les protéines sur l'ADN
tblastn -query $PROTEINS -db $DB_NAME -out $BLAST_OUTPUT \
        -outfmt "6 qseqid sseqid pident qstart qend sstart send evalue bitscore sseq" \
        -evalue 1e-5 -max_target_seqs 1

# Vérifier si BLAST a bien généré un fichier de résultats
if [ ! -s $BLAST_OUTPUT ]; then
    echo "Erreur: BLAST n'a retourné aucun alignement."
    exit 1
fi

# Étape 3: Filtrer les alignements à >90% d'identité et extraire les régions en format BED
awk '$3 == 100 {if ($6 < $7) print $2 "\t" $6-1 "\t" $7 "\t" $1; else print $2 "\t" $7-1 "\t" $6 "\t" $1}' $BLAST_OUTPUT > $BED_FILE

# Vérifier si des hits ont été trouvés
if [ ! -s $BED_FILE ]; then
    echo "Aucun alignement avec >90% d'identité trouvé."
    exit 1
fi

# Étape 4: Extraire les séquences nucléotidiques correspondantes avec bedtools (brin direct et reverse)
bedtools getfasta -fi $GENOME -bed $BED_FILE -fo $CDS_FORWARD
bedtools getfasta -s -fi $GENOME -bed $BED_FILE -fo $CDS_REVERSE

# Vérifier si des séquences ont été extraites
if [ ! -s $CDS_FORWARD ] && [ ! -s $CDS_REVERSE ]; then
    echo "Erreur: Aucune CDS extraite."
    exit 1
fi

# Étape 5: Filtrer les CDS qui commencent par ATG dans les deux fichiers
FILTERED_FORWARD="filtered_forward.fasta"
FILTERED_REVERSE="filtered_reverse.fasta"

awk 'BEGIN {RS=">"; ORS=""} NR>1 {header=$1; seq=toupper($2); if (substr(seq,1,3) == "ATG") print ">" header "\n" seq "\n"}' $CDS_FORWARD > $FILTERED_FORWARD

awk 'BEGIN {RS=">"; ORS=""} NR>1 {header=$1; seq=toupper($2); if (substr(seq,1,3) == "ATG") print ">" header "\n" seq "\n"}' $CDS_REVERSE > $FILTERED_REVERSE

# Fusionner les CDS des deux brins
cat $FILTERED_FORWARD $FILTERED_REVERSE > $FILTERED_CDS

# Vérifier si des séquences filtrées existent
if [ ! -s $FILTERED_CDS ]; then
    echo "Erreur: Aucune CDS valide (commençant par ATG) n'a été trouvée."
    exit 1
fi

# Supprimer les fichiers temporaires
rm -f $BLAST_OUTPUT $BED_FILE $DB_NAME.* $CDS_FORWARD $CDS_REVERSE $FILTERED_FORWARD $FILTERED_REVERSE

# Afficher le nombre de séquences extraites
echo "Nombre de CDS extraites : $(grep -c ">" $FILTERED_CDS)"
echo "Les CDS ont été enregistrées dans $FILTERED_CDS"

