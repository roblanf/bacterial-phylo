# input files
phylum="p__Patescibacteria"
taxonomy="../../raw_data/bac120_taxonomy_r202_clean.tsv"
alignment="../../raw_data/bac120_msa_reps_r202.faa"
tree="../../raw_data/bac120_r202.tree"

# Get the taxa we want - just for the phylum of interest
echo "Extracting taxa for $phylum"
awk -F"\t" -v phylum="$phylum" '$3 == phylum { print $1 }' "$taxonomy" > taxa.txt
N=$(wc -l < taxa.txt)
echo "Found $N taxa in $phylum"

# Get the alignment and starting subtree
echo "Subsetting alignment and tree"
faSomeRecords $alignment taxa.txt alignment.faa

# get the subtree
nw_labels -I $tree | sort > all_tips_sorted.txt # gets all tip labels
sort taxa.txt > taxa_sorted.txt
comm -23 all_tips_sorted.txt taxa_sorted.txt > tips_to_prune.txt # gets the tips not in my list
nw_prune $tree $(cat tips_to_prune.txt) > subtree.nwk
nw_topology -I subtree.nwk > subtree_topology.nwk # discard all unnecessary labels and branch lengths
sed -E 's/^\((.*)\);$/\1;/' subtree_topology.nwk > subtree_topology_clean.nwk


# now we get a subtree and sub-alignment of 100 random sequences
# these are alignment_100.faa, and subtree_100_topology_clean.nwk
shuf -n 100 taxa.txt > taxa_100.txt
faSomeRecords $alignment taxa_100.txt alignment_100.faa
sort taxa_100.txt > taxa_100_sorted.txt
comm -23 all_tips_sorted.txt taxa_100_sorted.txt > tips_to_prune_100.txt # gets the tips not in my list
nw_prune $tree $(cat tips_to_prune_100.txt) > subtree_100.nwk
nw_topology -I subtree_100.nwk > subtree_100_topology.nwk # discard all unnecessary labels and branch lengths
sed -E 's/^\((.*)\);$/\1;/' subtree_100_topology.nwk > subtree_100_topology_clean.nwk



# IQ-TREE analyses

#### C20 on the full dataset
iqtree -s alignment.faa -t subtree_topology_clean.nwk -m C20 -nt 60 -safe -pre C20fixed

#### Estimate PMSF on the full dataset: bacteria
# 1. Basic analysis to get the tree
iqtree -s alignment.faa -t subtree_topology_clean.nwk -m Q.pfam+G -nt 60 -safe -pre basic
# get the site frequencies
iqtree -s alignment.faa -ft basic.treefile -m Q.pfam+C60+G -nt 60 -safe -pre qpfam_c60_g_sf -n 0
# apply them with +G
iqtree -s alignment.faa -fs qpfam_c60_g_sf.sitefreq -t basic.treefile -m Q.pfam+C60+G -nt 60 -safe -pre qpfam_c60_g_tree_noboot
# apply them with +G and UFBOOT
iqtree -s alignment.faa -fs qpfam_c60_g_sf.sitefreq -t basic.treefile -m Q.pfam+C60+G -nt 60 -safe -pre qpfam_c60_g_tree_ufboot -bb 1000


#### Estimate PMSF on a 100 taxon dataset: bac3
# 1. get the tree with C60
iqtree -s alignment_100.faa -t subtree_100_topology_clean.nwk -m Q.pfam+C60+R8 -nt 60 -safe -pre 100_tree
# turns out that the R8 distributino looks very gamma-like, so I will switch to gamma
# get the site frequencies
iqtree -s alignment_100.faa -ft 100_tree.treefile -m Q.pfam+C60+G -nt 60 -safe -pre qpfam_c60_g_sf_100 -n 0
# apply them to the FULL alignment (using the same starting tree as before)
iqtree -s alignment.faa -fs qpfam_c60_g_sf_100.sitefreq -t basic.treefile -m Q.pfam+C60+G -nt 60 -safe -pre qpfam_c60_g_tree_noboot_100
# apply them to the FULL alignment with UFBOOT (using the same starting tree as before)
iqtree -s alignment.faa -fs qpfam_c60_g_sf_100.sitefreq -t basic.treefile -m Q.pfam+C60+G -nt 60 -safe -pre qpfam_c60_g_tree_ufboot_100 -bb 1000


#### GTRpmix on 100 taxon dataset, attempt with mwopt: bac2
# 1. Get the tree with C60 - take it from 100_tree above
# 2. Estimate GTR+C60 with mwopt
iqtree -s alignment_100.faa -m GTR20+C60+G4 --link-exchange --init-exchange q.pfam -te 100_tree.treefile -me 0.99 -nt 60 -safe -mwopt -pre GTR_c60_g_mwopt
# 3. get the site frequencies
iqtree -s alignment_100.faa -ft 100_tree.treefile -m GTRPMIX+C60+G4 -mdef GTR_c60_g_mwopt.GTRPMIX.nex -nt 60 -safe -pre GTR_c60_g_mwopt_sf_100 -n 0
# 4. apply them to the FULL alignment (using the same starting tree as before)
iqtree -s alignment.faa -fs GTR_c60_g_mwopt_sf_100.sitefreq -t basic.treefile -m GTRPMIX+C60+G4 -mdef GTR_c60_g_mwopt.GTRPMIX.nex -nt 60 -safe -pre GTR_c60_g_mwopt_tree_noboot_100



#### GTRpmix on 100 taxon dataset, site freqs on the whole thing? bacteria
# 1. Get the tree with C60 - take it from 100_tree above
# 2. Estimate GTR+C60 with mwopt: as above
# 3. get the site frequencies
iqtree -s alignment.faa -ft basic.treefile -m GTRPMIX+C60+G4 -mdef GTR_c60_g_mwopt.GTRPMIX.nex -nt 60 -safe -pre GTR_c60_g_mwopt_sf -n 0

# 4. apply them to the FULL alignment (using the same starting tree as before)
iqtree -s alignment.faa -fs GTR_c60_g_mwopt_sf.sitefreq -t basic.treefile -m GTRPMIX+C60+G4 -mdef GTR_c60_g_mwopt.GTRPMIX.nex -nt 60 -safe -pre GTR_c60_g_mwopt_tree_noboot_100

















#### Estimate a GTR20+F20 PMSF on the 100 taxon dataset: bac4
# 1. Estimate GTR20+F20
iqtree -s alignment_100.faa -te subtree_100_topology_clean.nwk -m GTR20+F20+R8 -nt 60 -safe -pre 100_GTR20F20

# 2. get the exchangeability matrix
iqtree2 -s alignment_100.faa -m GTR20+C60+R8 --link-exchange -te 100_tree.tree -me 0.99 -safe --init-exchange q.pfam -pre 100_pmix




#### C20 model no PMSF: bac2
iqtree -s alignment.faa -t subtree_topology_clean.nwk -m C20 -nt 60 -safe -pre C20fixed


