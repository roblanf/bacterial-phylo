# bacterial-phylo

This repo details some attempts to improve GTDB's trees.

# The Data:`/raw_data`

The data comprise:

4 releases of GTDB:
* 202
* 207
* 214
* 220

4 raw data files per release (e.g. from 202)
* full alignment: `gtdb_r202_bac120.faa.gz`
* taxonomy spreadsheet: `bac120_taxonomy_r202.tsv.gz`
* fasttree file: `bac120_r202.tree.tar.gz`
* 5K alignment: `bac120_msa_reps_r202.tar.gz`

Note that the tree file is estimated from the 5K alignment by GTDB with FastTree2.1

To get the raw data and uncompress it:

```{bash}
cd raw_data
bash download_data.sh
```


# The analysis

Each phylum has its own folder, e.g. `/patescibacteria`, with a readme. 