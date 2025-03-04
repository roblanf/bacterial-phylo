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



# IQ-TREE analyses


#### Estimate PMSF on the full dataset
# 1. Basic analysis
iqtree -s alignment.faa -t subtree_topology_clean.nwk -m Q.pfam+G -nt 60 -safe -pre basic


#### Estimate PMSF on a 100 taxon dataset



#### C20 model no PMSF

iqtree -s alignment.faa -t subtree_topology_clean.nwk -m C20 -nt 60 -safe -pre C20fixed


# C20 estimate weights
iqtree -s alignment.faa -t subtree_topology_clean.nwk -m Q.pfam+C20+R8 -mwopt -nt 120 -safe -pre C20estimated

# CAT-PMSF
# first optimise the branch length and C20 weights on a fixed starting tree
iqtree -s alignment.faa -te subtree_topology_clean.nwk -m Q.pfam+C20+R8 -mwopt -nt 60 -safe -pre C20_R8_fixed

# optionally we could infer the GTR matrix here, but that will likely take too long...

# then fix the model and get the PMSF site profiles
iqtree -s alignment.faa -ft C20_R8_fixed.tree -m Q.pfam+C20{}+R8{} -mwopt -nt 60 -safe -pre CATPMSF_C20