# Patescibacteria: analysis

Based on the observations from `patescibacteria_modelfit`, now we apply a PMSF model to each dataset.

We do this as follows:

1. Take a random subset of 128 taxa, estimate a tree with `Q.pfam+C60+G`, using the FastTree tree as the starting tree
2. Fix that tree, and estimate the parameters of a `GTR20+C60+G -mwopt` model on the 128 taxon subset
3. Use that model and the FastTree tree to estimate PMSF site profiles on the complete dataset
4. Update the complete tree with GTR20+C60+G PMSF model 
5. Re-estimate the site profiles with the new tree
6. Update the complete tree with GTR20+C60+G PMSF from step 5, with UFBOOT

## The analysis

All are in four scripts, which are identical except for the datasets they point to:

```{bash}

bash r202_dataset_analysis.sh
```



