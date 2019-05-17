
library(stringr)
library(xgboost)

# used to subset the genes
testfun <- function(a, y) {
  ranka <- rank(a)
  testres <- (sum(ranka[y == 0]) / sum(y == 0)) - (sum(ranka[y == 1]) / sum(y == 1))
  return(testres)
}


binGenes <- function(a) {
  brks <- quantile(as.numeric(a), probs=c(0, 0.25, 0.5, 0.75, 1.0), na.rm = T) ##### within one sample
  Xbin <- sapply(a, function(x) .bincode(x = x, breaks = brks))
  Xbin <- as.numeric(Xbin)
  Xbin
}


load('geneFeatureRanking.rda')  # loads preResList
idxList <- lapply(preResList, function(preRes)  which( (preRes < quantile(preRes, 0.05)) | (preRes > quantile(preRes, 0.95)) ))
names(idxList) <- names(preResList)
#allgenes <- unique(unlist(lapply(idxList, function(a) names(a))))  # used to subset full EBpp matrix, 0.10 adn 0.9


#load('/home/davidgibbs/Work/Cytokine-networks/Data/ebpp.rda')
#ebpp <- read.table('ebppSubset.tsv.bz2', header = T, sep = '\t', stringsAsFactors = F)
#colnames(ebpp) <- str_replace_all(colnames(ebpp), pattern = '\\.', replacement = '-')

# get cluster labels
# then make the rownames the aliquot barcodes to match up with TCGA calls
reportedScores <- read.table('five_signature_mclust_ensemble_results.tsv.gz', sep='\t', header=T, stringsAsFactors = F)
rownames(reportedScores) <- str_replace_all(reportedScores$AliquotBarcode, pattern = '\\.', replacement = '-')


#X <- ebpp[allgenes, rownames(reportedScores)] # get samples into correct order matching labels.
#Xmat <- as.matrix(X)  # this is ebpp_mini_subset.rda
#Xbin <- apply(Xmat, 2, function(x) binGenes(x))  # dim is samples in columns
load('ebpp_mini_subset_binned.rda')

all(colnames(Xbin) == rownames(reportedScores))
#[1] TRUE

Y <- reportedScores[,"ClusterModel1"]

preResList <- list()
modelList <- list()

for (cluster in 1:6) {
  print(cluster)

  Ybin <- ifelse(Y == cluster, yes = 1, no=0)

  print('building data sets')
  gidx <- names(idxList[[cluster]])
  Xsub <- t(Xbin[gidx,])

  #nrounds <- 15
  #param <- list(max_depth=2, eta=0.5, silent=1, nthread=4, objective='binary:logistic')

  #cat('running cross validation\n')
  # do cross validation, this will print result out as
  # [iteration]  metric_name:mean_value+std_value
  # std_value is standard deviation of the metric
  #xgb.cv(param, dtrain, nrounds, nfold=5, metrics={'error'})

  print('fitting model')
  bst <- xgboost(data = Xsub, label = Ybin, max_depth = 2, eta = 0.5, nrounds = 33, nthread = 5, objective = "binary:logistic")
  pred <- predict(bst, Xsub)
  err <- mean(as.numeric(pred > 0.5) != Ybin)
  print(table((pred > 0.5), Ybin))

  modelList[[cluster]] <- bst
  preResList[[cluster]] <- preRes

  print('making plot')
  library(plotROC )
  df <- data.frame(predictions=pred, labels=Ybin)
  rocplot <- ggplot(df, aes(m = predictions, d = labels))+ geom_roc(n.cuts=20,labels=FALSE)
  rocplot <- rocplot + style_roc(theme = theme_grey) + geom_rocci(fill="pink") + ggtitle(paste0('C', cluster, ', ROC'))
  ggsave(paste0('cluster_',cluster,'.png'), plot = rocplot)

}
#save(preResList, file='preResList.rda')
save(modelList, file='modelList.rda')



predOne <- function(ex, modelList) {
  # ex should already be binned.
  # model wants data to be samples x genes
  res0 <- lapply(modelList, function(ml) {
    gs  <- ml$feature_names
    ex2 <- matrix(data=ex[gs], ncol=length(gs), nrow=1)
    colnames(ex2) <- gs
    predict(ml, ex2)
  })
  unlist(res0)
}

trainingPred <- t(apply(Xbin, 2, function(a) predOne(a, modelList) ))
colnames(trainingPred) <- c('C1', 'C2', 'C3', 'C4', 'C5', 'C6')
bestPred <- apply(trainingPred, 1, function(a) which(a == max(a)) )
trainingRes <- as.data.frame(cbind(Subtype=Y, bestPred, trainingPred))

t1 <- (table(Subtype=trainingRes$Subtype, Prediction=trainingRes$bestPred))

diag(t1) / rowSums(t1)

