# Immune-Subtype-Clustering


This repo contains the code necessary to reproduce the signature scores and clusters found in "The Immune Landscape of Cancer".

NOTE! This method is very sensitive to the software pipeline used in quantifying genes. It does not work well with FPKM, RPKM, TPM, etc.

***A new, more robust method is found at; https://github.com/Gibbsdavidl/ImmuneSubtypeClassifier***

In this repo:
The signature_mclust_ensemble.R file is a script that reads the signature scores in five_signature_mclust_ensemble_results.tsv.gz (same as found in the manuscript supplement) and generates the ensemble of mclust models, and finally compares the clustering to the cluster labels reported in the manuscript.

The ImmuneSigs68_function.R file contains the code necessary to produce the signature scores.


