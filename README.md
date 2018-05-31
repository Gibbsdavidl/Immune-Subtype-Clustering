# Immune-Subtype-Clustering


This repo contains the code necessary to reproduce the clusters found in "The Immune Landscape of Cancer".

The signature_mclust_ensemble.R file is a script that reads the signature scores in five_signature_mclust_ensemble_results.tsv.gz (same as found in the manuscript supplement) and generates the ensemble of mclust models, and finally compares the clustering to the cluster labels reported in the manuscript.

The ImmuneSigs68_function.R file contains the code necessary to produce the signature scores.


