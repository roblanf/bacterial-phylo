# Patescibacteria

Details of the analysis are in the analysis folder. Here's a summary

## Results

AICs are approx because they don't account for parameters estimated in earlier steps. These can be ~300 parameters (rate matrix, plus ~60 profile weights). 

| ID  | LnL   | AIC(approx) | Tree Length | % internal | Total Time (tree) | Total Time (UFBoot) | Summary  |
|-----|-------|-----------|----------|-------|-----|-----|-----|
| 1 | -5457790 | 10921981 | 487 | 27.5678% | 24h:1m:14s | NA |  `-m C20`, that's it |
| 2 | -5327202 | 10660842 | 713 | 27.5436% | 9h:21m:28s | 11h:30m:51s |  PMSF with Q.pfam+C60+G on the full dataset |
| 3 | -5381055 | 10768548 | 668 | 28.0951% | 7h:54m:1s | 11h:48m:28s | PMSF as for 2 but site freqs all done from 100 random taxa |
| 4 | -5306068 | 10618575 | 525 | 27.8594% | 11h:13m:38s |  |  |
| 5 |  |  |  |  |  |  |  |
| 3 |  |  |  |  |  |  |  |
| 3 |  |  |  |  |  |  |  |
| 3 |  |  |  |  |  |  |  |
| 3 |  |  |  |  |  |  |  |
| 3 |  |  |  |  |  |  |  |


## The analyses

### 1. C20

Just C20 on the full dataset. Nothing complex. But a good starting point.
Starting tree is the FastTree tree.

``` {bash}
#### C20 on the full dataset
iqtree -s alignment.faa -t subtree_topology_clean.nwk -m C20 -nt 60 -safe -pre C20fixed
```

Time: 24h:1m:14s

This takes a while! Because C20 means a lot of likelihoods.

### 2. PMSF with Q.pfam+C60+G

Standard PMSF with the most complex model we can realistically do on a huge dataset. I.e. no optimising weights.

```{bash}
#### Estimate PMSF on the full dataset: bacteria
# 1. Basic analysis to get the tree
iqtree -s alignment.faa -t subtree_topology_clean.nwk -m Q.pfam+G -nt 60 -safe -pre basic
# get the site frequencies
iqtree -s alignment.faa -ft basic.treefile -m Q.pfam+C60+G -nt 60 -safe -pre qpfam_c60_g_sf -n 0
# apply them with +G
iqtree -s alignment.faa -fs qpfam_c60_g_sf.sitefreq -t basic.treefile -m Q.pfam+C60+G -nt 60 -safe -pre qpfam_c60_g_tree_noboot
# apply them with +G and UFBOOT
iqtree -s alignment.faa -fs qpfam_c60_g_sf.sitefreq -t basic.treefile -m Q.pfam+C60+G -nt 60 -safe -pre qpfam_c60_g_tree_ufboot -bb 1000

```

Timing:

1. 4h:23m:8s
2. 0h:2m:1s
3. 4h:56m:19s
4. 7h:5m:42s

Observations: it's the tree search which takes all the time. And C60 PMSF is almost just as quick as Q.pfam+G, as expected. Also, UFBOOT roughly doubles the time. 



### 3. PMSF with parameters from a 100 taxon subset

The point here is a fairly direct comparison with the above. I estimate the site frequencies themselves on a subset of 100 taxa, then apply them to the full dataset. 

```{bash}
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
```

Note that the timing for step 1 is long, and off, because I tried an R8 model. 

1. 2h:42m:58s
2. 0h:0m:3s
3. 5h:11m:0s
4. 9h:5m:27s

Observations: not sure why this PMSF analysis took so much longer. Worse model fit? As expected, the lnL is lower. But it's not a LOT lower. So that's encouraging! This suggests we could try doing site profiles from super fancy models on a subset. Let's see how that goes.


### 4. PMSF with parameters from a 100 taxon subset, using `-mwopt` and `GTRPMIX`

The idea is that maybe a super fancy model on a small subset gives better site frequencies than a simpler model on the full data. Let's see.


```{bash}
#### GTRpmix on 100 taxon dataset, attempt with mwopt: bac2
# 1. Get the tree with C60 - take it from 100_tree above
# 2. Estimate GTR+C60 with mwopt
iqtree -s alignment_100.faa -m GTR20+C60+G4 --link-exchange --init-exchange q.pfam -te 100_tree.treefile -me 0.99 -nt 60 -safe -mwopt -pre GTR_c60_g_mwopt
# 3. get the site frequencies
iqtree -s alignment_100.faa -ft 100_tree.treefile -m GTRPMIX+C60+G4 -mdef GTR_c60_g_mwopt.GTRPMIX.nex -nt 60 -safe -pre GTR_c60_g_mwopt_sf_100 -n 0
# 4. apply them to the FULL alignment (using the same starting tree as before), then try ufboot
iqtree -s alignment.faa -fs GTR_c60_g_mwopt_sf_100.sitefreq -t basic.treefile -m GTRPMIX+C60+G4 -mdef GTR_c60_g_mwopt.GTRPMIX.nex -nt 60 -safe -pre GTR_c60_g_mwopt_tree_noboot_100
iqtree -s alignment.faa -fs GTR_c60_g_mwopt_sf_100.sitefreq -t basic.treefile -m GTRPMIX+C60+G4 -mdef GTR_c60_g_mwopt.GTRPMIX.nex -nt 60 -safe -pre GTR_c60_g_mwopt_tree_ufboot_100 -bb 1000

```

Timing:

1. 2h:42m:58s
2. 2h:53m:12s
3. 0h:0m:2s
4. 5h:37m:26s
5. 

Observations: best likelihood yet! So this works. Also, the tree lenght is a lot shorter, so I'm not sure what to think of that. Less LBA perhaps? 

