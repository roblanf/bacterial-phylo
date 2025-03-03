# input files
phylum="p__Patescibacteria"
taxonomy="../../raw_data/bac120_taxonomy_r202_clean.tsv"
alignment="../../raw_data/gtdb_r202_bac120.faa"
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
comm -23 all_tips_sorted.txt tax_sorted.txt > tips_to_prune.txt # gets the tips not in my list
nw_prune $tree $(cat tips_to_prune.txt) > subtree.nwk
nw_topology -I subtree.nwk > subtree_topology.nwk # discard all unnecessary labels and branch lengths
sed -E 's/^\((.*)\);$/\1;/' subtree_topology.nwk > subtree_topology_clean.nwk



# IQ-TREE analyses

# Basic
iqtree -s alignment.faa -t subtree_topology_clean.nwk -m Q.pfam+R8 -nt 120 -safe -pre basic

# C20
iqtree -s alignment.faa -t subtree_topology_clean.nwk -m Q.pfam+C20+R8 -nt 120 -safe -pre C20fixed

# C20 estimate weights
iqtree -s alignment.faa -t subtree_topology_clean.nwk -m Q.pfam+C20+R8 -mwopt -nt 120 -safe -pre C20estimated

