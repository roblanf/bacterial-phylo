#!/usr/bin/env bash

# A script to get all the raw data for these analyses
# 4 files per release: full alignment, taxonomy file, tree from fasttree, and short (5K) alignment. The short alignment is what is used to make the tree.
# 4 releases: 202, 207, 214, 220 (in order of size)

# full alignments (120 loci)
wget https://data.ace.uq.edu.au/public/gtdb/temporary_folder/full_msas/gtdb_r202_bac120.faa.gz
wget https://data.ace.uq.edu.au/public/gtdb/temporary_folder/full_msas/gtdb_r207_bac120.faa.gz
wget https://data.ace.uq.edu.au/public/gtdb/temporary_folder/full_msas/gtdb_r214_bac120.faa.gz
wget https://data.ace.uq.edu.au/public/gtdb/temporary_folder/full_msas/gtdb_r220_bac120.faa.gz

# taxonomy files
wget https://data.gtdb.ecogenomic.org/releases/release202/202.0/bac120_taxonomy_r202.tsv.gz
wget https://data.gtdb.ecogenomic.org/releases/release207/207.0/bac120_taxonomy_r207.tsv.gz
wget https://data.gtdb.ecogenomic.org/releases/release214/214.1/bac120_taxonomy_r214.tsv.gz
wget https://data.gtdb.ecogenomic.org/releases/release220/220.0/bac120_taxonomy_r220.tsv.gz

# trees from fasttree
wget https://data.gtdb.ecogenomic.org/releases/release202/202.0/bac120_r202.tree.tar.gz
wget https://data.gtdb.ecogenomic.org/releases/release207/207.0/bac120_r207.tree.tar.gz
wget https://data.gtdb.ecogenomic.org/releases/release214/214.1/bac120_r214.tree.tar.gz
wget https://data.gtdb.ecogenomic.org/releases/release220/220.0/bac120_r220.tree.gz

# short (5K) alignments
wget https://data.gtdb.ecogenomic.org/releases/release202/202.0/genomic_files_reps/bac120_msa_reps_r202.tar.gz
wget https://data.gtdb.ecogenomic.org/releases/release207/207.0/genomic_files_reps/bac120_msa_reps_r207.tar.gz
wget https://data.gtdb.ecogenomic.org/releases/release214/214.1/genomic_files_reps/bac120_msa_reps_r214.tar.gz
wget https://data.gtdb.ecogenomic.org/releases/release220/220.0/genomic_files_reps/bac120_msa_reps_r220.faa.gz


# now unzip the lot
gunzip *.gz

# and untar the three remaining trees
tar -xvf bac120_r202.tree.tar
tar -xvf bac120_r207.tree.tar
tar -xvf bac120_r214.tree.tar

# I don't like the formatting of the tsv files, so let's fix that next
awk 'BEGIN { FS=OFS="\t" } { gsub(";", "\t", $2) } 1' bac120_taxonomy_r202.tsv > bac120_taxonomy_r202_clean.tsv
awk 'BEGIN { FS=OFS="\t" } { gsub(";", "\t", $2) } 1' bac120_taxonomy_r207.tsv > bac120_taxonomy_r207_clean.tsv
awk 'BEGIN { FS=OFS="\t" } { gsub(";", "\t", $2) } 1' bac120_taxonomy_r214.tsv > bac120_taxonomy_r214_clean.tsv
awk 'BEGIN { FS=OFS="\t" } { gsub(";", "\t", $2) } 1' bac120_taxonomy_r220.tsv > bac120_taxonomy_r220_clean.tsv
