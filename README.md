# STAGE M2 BIO-INFORMATIQUE, Thomas BOUTET

Ce stage avait deux objectifs principaux à savoir la création de modèles HMM à partir de gènes orthologues et la création d'un site web qui permet de détecter les petites régions codantes (sCDS) chez les Microsporidies.

## Table des matières

- [Modèles HMM](#modeles-hmm)  
- [Interface Web sCDS_web](#interface-web-scds_web)

## Modeles HMM

*EN TRAVAUX*

---

## Interface Web sCDS_web

L’interface **sCDS_web** fournit un outil en ligne pour détecter des petites séquences codantes (sCDS) dans des génomes de microsporidies. Elle a pour but d'être une suite à l'outil MicroAnnot (https://www.microannot.org/) en reprenant les fichiers de sortie de celui-ci (EMBL), ainsi que le fichier du génome de la microsporidie à annoter au format FASTA.

### Arborescence

```plaintext
sCDS_web/
├── bin/
│   ├── detect_sCDS.py
│   ├── pipeline_microannot_sCDS.sh
│   ├── script_petits_genes_microannot.py
│   └── update_bdd.py
├── data/
│   ├── database_ortholog/
│   │   └── sCDS_ortholog.fasta
│   ├── database_not_ortholog/
│   │   └── sCDS_not_ortholog.fasta
│   ├── work_sCDS_ortholog.fasta
│   ├── work_sCDS_not_ortholog.fasta
│   ├── work_sCDS_ortholog_backup.fasta
│   └── work_sCDS_not_ortholog_backup.fasta
├── result/
├── static/
│   └── images/
│       └── background.png
└── templates/
    ├── index.html
    └── databases.html
```

#### 1. Installer Miniconda (si nécessaire)

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh
bash ~/miniconda.sh -b -p $HOME/miniconda
export PATH="$HOME/miniconda/bin:$PATH"
```

#### 2. Fermer et rouvrir le terminal ou recharger le profil

```bash
source ~/.bashrc
```

#### 3. Créer et activer l’environnement Conda dédié

```bash
conda create -n scds_web python=3.8 flask biopython blast -c conda-forge -c bioconda -y
conda activate scds_web
```

#### 4. Cloner le dépôt et se placer dans `sCDS_web`

```bash
git clone https://github.com/thboutet292/STAGE_M2.git
cd STAGE_M2/sCDS_web
```

#### 5. Lancer l’application Flask

```bash
python bin/detect_sCDS.py
```
> [!NOTE]  
> Ouvrez ensuite la page web via le lien indiqué dans le terminal, il devrait ressembler à 'http://localhost:5000'


> [!IMPORTANT]  
> Si l'utilisateur préfère utiliser directement l'outil depuis le terminal cela est possible via :

```bash
bash bin/pipeline_microannot_sCDS.sh <genome_fasta> <embl_directory> <output_fasta> [--use-updated]
```
+ <genome_fasta> = chemin vers le fichier FASTA du génome
+ <embl_directory> = chemin vers le répertoire contenant les fichiers EMBL de MicroAnnot correspondant au génome
+ <output_fasta> = nom du fichier de sortie
+ --use-updated = option pour utiliser les bases de données mises à jour, par default désactivé
 



