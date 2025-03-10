# input files
phylum="p__Patescibacteria"
taxonomy="../../raw_data/bac120_taxonomy_r202_clean.tsv"
alignment="../../raw_data/bac120_msa_reps_r202.faa"
tree="../../raw_data/bac120_r202.tree"
log="log_r202.txt"
subsample=128
threads=60

echo "analysing phylum: $phylum" > $log
echo "taxonomy file: $taxonomy" >> $log
echo "alignment file: $alignment" >> $log
echo "tree file: $tree" >> $log


echo "" >> $log
echo "" >> $log
echo "1. Extracting taxa for $phylum from $taxonomy" >> $log

awk -F"\t" -v phylum="$phylum" '$3 == phylum { print $1 }' "$taxonomy" > taxa.txt
N=$(wc -l < taxa.txt)
echo "Found $N taxa in $phylum in $taxonomy" >> $log


echo "" >> $log
echo "" >> $log
echo " 2. Subsetting alignment and tree" >> $log

faSomeRecords $alignment taxa.txt alignment.faa
N_taxa_aln=$(grep ">" alignment.faa | wc -l)
echo "Found $N_taxa_aln taxa from $phylum in $alignment" >> $log

# get the subtree of just the phylum
nw_labels -I $tree | sort > all_tips_sorted.txt # gets all tip labels
sort taxa.txt > taxa_sorted.txt
comm -23 all_tips_sorted.txt taxa_sorted.txt > tips_to_prune.txt # gets the tips not in my list
nw_prune $tree $(cat tips_to_prune.txt) > subtree.nwk
nw_topology -I subtree.nwk > subtree_topology.nwk # discard all unnecessary labels and branch lengths
sed -E 's/^\((.*)\);$/\1;/' subtree_topology.nwk > subtree_topology_clean.nwk
nw_labels -I subtree_topology_clean.nwk | sort > subtree_tips.txt
N_taxa_tree=$(wc -l subtree_tips.txt | cut -d ' ' -f1)
echo "Found $N_taxa_tree taxa from $phylum in $tree" >> $log

# make nice files
cp subtree_topology_clean.nwk tree.nex
cp subtree_tips.txt species.txt




# now we get a subtree and sub-alignment of 128 random sequences

echo "" >> $log
echo "" >> $log
echo " 4. Subsampling $subsample taxa from alignment and tree for model estimation" >> $log

shuf -n $subsample species.txt > sub_taxa.txt
faSomeRecords $alignment sub_taxa.txt sub_alignment.faa
sort sub_taxa.txt > sub_taxa_sorted.txt
comm -23 all_tips_sorted.txt sub_taxa_sorted.txt > tips_to_prune.txt # gets the tips not in my list
nw_prune $tree $(cat tips_to_prune.txt) > sub_tree.nwk
nw_topology -I sub_tree.nwk > sub_tree_topology.nwk # discard all unnecessary labels and branch lengths
sed -E 's/^\((.*)\);$/\1;/' sub_tree_topology.nwk > sub_tree_topology_clean.nwk
cp sub_tree_topology_clean.nwk sub_tree.nex
nw_labels -I sub_tree.nex | sort > sub_tree_tips.txt

N_taxa_sub_aln=$(grep ">" sub_alignment.faa | wc -l)
N_taxa_sub_tree=$(wc -l sub_tree_tips.txt | cut -d ' ' -f1)
echo "sub_alignmnet.faa has $N_taxa_sub_aln taxa" >> $log
echo "sub_tree.nex has $N_taxa_sub_tree taxa" >> $log


echo "" >> $log
echo "" >> $log
echo " 5. Files for analysis" >> $log
echo "Complete alignment of $N_taxa_aln taxa: alignment.faa" >> $log
echo "Complete tree of $N_taxa_tree taxa: tree.nex" >> $log
echo "Sub alignment of $N_taxa_sub_aln taxa: sub_alignment.faa" >> $log
echo "Sub tree of $N_taxa_sub_tree taxa: sub_tree.nex" >> $log


echo "" >> $log
echo "" >> $log
echo " 5. Running IQ-TREE" >> $log

#### GTRpmix on 250 taxon dataset, site freqs on the whole thing: bacteria
# 1. Get the tree with C60
/usr/bin/time -v -o 01_Qpfam_C60G_sub_tree.txt iqtree -s sub_alignment.faa -t sub_tree.nex -m Q.pfam+C60+G -nt $threads -safe -pre 01_Qpfam_C60G_sub_tree

# 2. Estimate GTR+C60 with mwopt on a 250 subset
/usr/bin/time -v -o 02_GTR_c60_g_mwopt.txt iqtree -s sub_alignment.faa -m GTR20+C60+G4 --link-exchange --init-exchange q.pfam -te 01_Qpfam_C60G.treefile -me 0.99 -nt $threads -safe -mwopt -pre 02_GTR_c60_g_mwopt

# 3. get the site frequencies using the fasttree tree
/usr/bin/time -v -o 03_sitefreqs_iteration1.txt iqtree -s alignment.faa -ft tree.nex -m GTRPMIX+C60+G4 -mdef 02_GTR_c60_g_mwopt.GTRPMIX.nex -nt $threads -safe -pre 03_sitefreqs_iteration1 -n 0

# 4. update the tree with this PMSF model
/usr/bin/time -v -o 04_PMSF_tree_iteration1.txt iqtree -s alignment.faa -fs 03_sitefreqs_iteration1.sitefreq -t tree.nex -m GTRPMIX+C60+G4 -mdef 02_GTR_c60_g_mwopt.GTRPMIX.nex -nt $threads -safe -pre 04_PMSF_tree_iteration1

# 5. update the site freqs with the new tree
/usr/bin/time -v -o 05_sitefreqs_iteration2.txt iqtree -s alignment.faa -ft 04_PMSF_tree_iteration1.treefile -m GTRPMIX+C60+G4 -mdef 02_GTR_c60_g_mwopt.GTRPMIX.nex -nt $threads -safe -pre 05_sitefreqs_iteration2 -n 0

# 6. Full analysis with UFBOOT
/usr/bin/time -v -o 06_PMSF_tree_ufboot.txt iqtree -s alignment.faa -fs 05_sitefreqs_iteration2.sitefreq -t 05_sitefreqs_iteration2.treefile -m GTRPMIX+C60+G4 -mdef 02_GTR_c60_g_mwopt.GTRPMIX.nex -nt $threads -safe -pre 06_PMSF_tree_ufboot -bb 1000


