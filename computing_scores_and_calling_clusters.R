
# David L Gibbs
# dgibbs@systemsbiology.org

# this script will form the backend of the shiny app.


# this function computes scores given some expression data.
newScores <- function(newdata, logflag) {
	
  source('src/ImmuneSigs68_function.R')
  load('data/comparative_immuneSigs_geneLists4.rda')
  zscore.cols2<-function(x){
    return((apply(x, 2, function(x) (x - median(na.omit(x)))/sd(na.omit(x)))))
  }

  # 1 Recomputed non-Z-scored scores from the EBPP matrix.
  ebppScores <- load("data/ebpp_scores.rda")

  # 2 we get some new data in.. require:
  #    it is FPKM  
  #    and   log2 transformed 
  #    and   gene symbols as row names
  #    and   median centered
  dat <- newdata  # needs row names as symbols

  # just in case we need a log transform
  if (logflag) {
    datlog2 <- log2(dat+1)
  } else {
    datlog2 <- dat
  }

  ### median scaled for each gene
  datmeds <- apply(datlog2, 1, median, na.rm=T)  
  datscaled <- sweep(datlog2,1,datmeds,'-')
  datScores <- ImmuneSigs_function(datscaled, sigs1_2_eg2,sigs12_weighted_means,sigs12_module_weights,sigs1_2_names2,sigs1_2_type2)

  # then batch correction between scores...
  library(sva)
  df <- as.matrix(cbind(datScores[rownames(ebppScores),], ebppScores))
  batch <- c(rep(1,ncol(datScores)), rep(2,ncol(ebppScores)))
  modcombat = model.matrix(~1, data=as.data.frame(t(df)))
  combat_edata = ComBat(dat=df, batch=batch, mod=modcombat, par.prior=TRUE,prior.plots=FALSE)

  # and we subset the 5 scores used in clustering
  idx <- c("LIexpression_score", "CSF1_response", "TGFB_score_21050467", "Module3_IFN_score", "CHANG_CORE_SERUM_RESPONSE_UP")
  scores <- t(combat_edata[idx,])
  zscores <- zscore.cols2(scores)

  # load the clustering model trained on all pancan data.
  load("data/wolf_set_model1.rda")
  calls <- consensusEnsemble(mods1, zscores)
  maxcalls <- apply(calls$.Data, 1, function(a) which(a == max(a))[1])

  # and get the reported scores from the manuscript
  wolf <- read.table("data/five_signature_mclust_ensemble_results.tsv", sep='\t', header=T, stringsAsFactors = F)
  wolfscrs <- wolf[,c(5:9)]
  wolfNames <- str_replace_all(wolf$AliquotBarcode, '\\.', '-')

  # Then we make sure the pancan cluster labels have not changed *much*  *how much is OK?*
  idx <- match(table=rownames(scores), x = wolfNames)
  all(rownames(scores)[idx] == wolfNames)
  table(New=maxcalls[idx], Wolf=wolf$ClusterModel1)
}