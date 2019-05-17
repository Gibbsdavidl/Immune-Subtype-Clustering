
library(stringr)
library(xgboost)

testfun <- function(a, y) {
  ranka <- rank(a)
  testres <- (sum(ranka[y == 0]) / sum(y == 0)) - (sum(ranka[y == 1]) / sum(y == 1))
  return(testres)
}

load('/home/davidgibbs/Work/Cytokine-networks/Data/ebpp.rda')
#ebpp <- read.table('ebppSubset.tsv.bz2', header = T, sep = '\t', stringsAsFactors = F)
#colnames(ebpp) <- str_replace_all(colnames(ebpp), pattern = '\\.', replacement = '-')

# get cluster labels
# then make the rownames the aliquot barcodes to match up with TCGA calls
reportedScores <- read.table('five_signature_mclust_ensemble_results.tsv.gz', sep='\t', header=T, stringsAsFactors = F)
rownames(reportedScores) <- str_replace_all(reportedScores$AliquotBarcode, pattern = '\\.', replacement = '-')


X <- ebpp[, rownames(reportedScores)]
Xmat <- as.matrix(X)
Y <- reportedScores[,"ClusterModel1"]

preResList <- list()
modelList <- list()

for (cluster in 1:6) {
  print(cluster)

  Ybin <- ifelse(Y == cluster, yes = 1, no=0)

  print('gen gene rank scores')
  preRes <- apply(Xmat, 1, FUN=function(a) testfun(a,Ybin))

  print('building data sets')
  idx <- which( (preRes < quantile(preRes, 0.10)) | (preRes > quantile(preRes, 0.90)) )
  Xsub <- Xmat[idx,]
  Xsub[is.na(Xsub)] <- 0
  Xscl <- t(scale(Xsub)) # scale each sample
  brks <- quantile(as.numeric(Xscl), probs=c(0, 0.25, 0.5, 0.75, 1.0))
  Xbin <- apply(Xscl, 2, function(x) .bincode(x = x, breaks = brks))
  Xbin <- apply(Xbin, 2, as.numeric)

  #nrounds <- 15
  #param <- list(max_depth=2, eta=0.5, silent=1, nthread=4, objective='binary:logistic')

  #cat('running cross validation\n')
  # do cross validation, this will print result out as
  # [iteration]  metric_name:mean_value+std_value
  # std_value is standard deviation of the metric
  #xgb.cv(param, dtrain, nrounds, nfold=5, metrics={'error'})

  print('fitting model')
  bst <- xgboost(data = Xbin, label = Ybin, max_depth = 2, eta = 0.5, nrounds = 33, nthread = 5, objective = "binary:logistic")
  pred <- predict(bst, Xbin)
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
save(preResList, file='preResList.rda')
save(modelList, file='modelList.rda')

