

# generating the signature scores. #

# running sigs
library(stringr)
load('ebpp.rda')  # the pancan batch corrected expression data with symbols in the rows.
load('data/comparative_immuneSigs_geneLists4.rda')

data2 <- ebpp[,1:100]  # smaller test set

source('ImmuneSigs68_function.R')
scrs <- ImmuneSigs_function(data2, sigs1_2_eg2,sigs12_weighted_means,sigs12_module_weights,sigs1_2_names2,sigs1_2_type2)

library(superheat)
superheat(scrs, left.label.text.size=2)

