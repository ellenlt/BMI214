require(class)
require(ROCR)
require(e1071)
require(FSelector)
require(foreign)
require(arules)

# ------------------------------------------------
# ------------------------------------------------
# LEUKEMIA DATA
# ------------------------------------------------
# ------------------------------------------------
# 72 leukemia patients (rows) and expression values for 150 genes (columns)

leukemia <-read.table("http://bmi214.stanford.edu/files/a2/leukemia.csv",header=TRUE,sep=",")
summary(leukemia)

# ------------------------------------------------
# K NEAREST NEIGHBORS
# ------------------------------------------------
## Train the classifier, with leave-one-out cross-validation and k=5:
knn.pred <- knn.cv(k=5,
                   cl=leukemia[,c("leukemia_type")],
                   train=leukemia[,names(leukemia)!="leukemia_type"],
                   prob=TRUE)
## Generate the confusion matrix for the classifier:
knn.actual <- leukemia[,c("leukemia_type")]
confusion.table <- table(knn.pred, knn.actual)
confusion.table

## Visualize ROC curve - TP (sensitivity) vs FP (1-specificity)
knn.prob <- attr(knn.pred,"prob")
knn.prob[knn.pred == "AML"] = 1 - knn.prob[knn.pred == "AML"]
knn.rocr.pred <- prediction(knn.prob,
                            leukemia[,c("leukemia_type")],
                            label.ordering=c("AML","ALL"))
knn.perf <- performance(knn.rocr.pred,"tpr","fpr")
plot(knn.perf)

# ------------------------------------------------
# CROSS VALIDATION
# ------------------------------------------------
## This is function that we will use to run cross validation
## for several methods:
cross.validation <- function(data, response, FUN, num.folds=5) {
  
  fold.size = floor(nrow(data)/num.folds)
  
  ## Loop for each fold:
  for (i in 1:num.folds) {
    
    ## Determine the indexes for test and train partitions:
    start.index <- (i-1)*fold.size + 1
    if(i == num.folds) end.index <- nrow(data) else end.index <-
      start.index +
      fold.size - 1
    excl.index <- c(start.index:end.index)
    test <- data[excl.index,]
    train <- data[-excl.index,]
    
    
    ## Test:
    if (i==1) {
      weights <- FUN(as.formula(paste(response,'~.',sep='')),test)
    } else {
      weights[,1] <- weights[,1] +
        FUN(as.formula(paste(response,'~.',sep='')),test)[,1]
    }
  }
  weights[,1] <- weights[,1]/num.folds
  return(weights[order(weights[,1], decreasing=TRUE),,drop=FALSE])
  
}

## zeroR returns a data frame, where the rows represent:
## TP, TN, FP, FN
classOne<-"1"  #ALL
classTwo<-"2"  #AML
zeroR <- function(formula, data) {
  class.summary <- summary(data[,c("leukemia_type")])
  classes <- labels(class.summary)
  if ( class.summary[as.numeric(classOne)] > class.summary[as.numeric(classTwo)] ) {
    other.class <- classes[as.numeric(classTwo)]
    predictions <- rep(classes[as.numeric(classOne)],times=nrow(data))
  } else {
    other.class <- classes[as.numeric(classOne)]
    predictions <- rep(classes[as.numeric(classTwo)],times=nrow(data))
  }
  confusion.table <- table(predictions, data[,c("leukemia_type")])
  results <- matrix(nrow=4,data=0)
  results[1,] <- confusion.table[c(predictions[as.numeric(classOne)]),c(predictions[as.numeric(classOne)])]
  results[3,] <- confusion.table[c(predictions[as.numeric(classOne)]),other.class]
  results.frame <- data.frame(results,row.names=c("TP","TN","FP","FN"))
  return(results.frame)
}

cross.validation(leukemia,"",zeroR)

# ------------------------------------------------
# NAIVE BAYES
# ------------------------------------------------
# Fitting Naive Bayes to predict leukemia type
nb.fit <- naiveBayes(leukemia[,names(leukemia)!="leukemia_type"],
                     leukemia[,c("leukemia_type")])
table(predict(nb.fit,leukemia[,names(leukemia)!="leukemia_type"]),
      leukemia[,c("leukemia_type")])

# ------------------------------------------------
# FEATURE SELECTION
# ------------------------------------------------
# Choose subset of five genes (arbitrarily the first five)
# and using KNN to assess predictive power
leukemia.subset <- leukemia[,c(1:5,ncol(leukemia))]
feature.col.index <- names(leukemia.subset) != "leukemia_type"
leukemia.subset.knn.pred <- knn.cv(k=5,
                                   cl=leukemia.subset[,c("leukemia_type")],
                                   train=leukemia.subset[,feature.col.index])
leukemia.subset.confusion.table <- table(leukemia.subset.knn.pred,
                                         leukemia.subset[,c("leukemia_type")])
leukemia.subset.confusion.table

## Feature selection by information gain
weights <- information.gain("leukemia_type ~.",leukemia)
weights[with(weights, order(-attr_importance)),]
rownames(weights)[with(weights, order(-attr_importance))]

# ------------------------------------------------
# ------------------------------------------------
# YEAST EXPRESSION DATA
# ------------------------------------------------
# ------------------------------------------------
# Expression profiles of 2467 yeast genes across 79 experiments
# Rows are genes, columns are experiments
yeast <-
  read.table("http://bmi214.stanford.edu/files/a2/yeast.dat.csv",
             header=TRUE,sep=",")
summary(yeast)

# ------------------------------------------------
# K-MEANS CLUSTERING
# ------------------------------------------------
yeast.clust <- kmeans(yeast[, - ncol(yeast)], 2)
yeast.clust
sum(yeast.clust$cluster==1)
sum(yeast.clust$cluster==2)

## Multiple random initializations
yeast.clust.15 <- kmeans(yeast[, - ncol(yeast)], 2,nstart=100)
sum(yeast.clust.15$cluster==1)
sum(yeast.clust.15$cluster==2)

## Compare actual data to computed clusters
yeast.cluster.vs.label <- table(yeast.clust$cluster, yeast[,c("ribosomal")])
yeast.cluster.vs.label

# ------------------------------------------------
# K NEAREST NEIGHBORS
# ------------------------------------------------
yeast.knn <- knn.cv(k=5,cl=yeast[,c("ribosomal")],train=yeast[,-ncol(yeast)])
yeast.knn.table <- table(yeast.knn,yeast[,c("ribosomal")])
yeast.knn.table


# ------------------------------------------------
# ------------------------------------------------
# GENOTENUREITIS DATA
# ------------------------------------------------
# ------------------------------------------------
# Contains genotype-phenotype data for 557 human subjects.
# Each row is subject, first column is ID, next 188 columns are SNPs
# 11: both copies are allele 1; 12: one copy allele 1, one copy allele 2
# 22: both copies are allele 2

genotenureitis <-read.arff("http://bmi214.stanford.edu/files/a2/genotenureitus1.arff")
summary(genotenureitis)
# Filter out irep attribute because it is not useful
genotenureitis.filtered <-genotenureitis[, colnames(genotenureitis) != 'irep']

# ------------------------------------------------
# ASSOCIATION RULE MINING
# ------------------------------------------------
rules <- apriori(genotenureitis.filtered,
                 parameter=list(supp=0.5,
                                conf=0.9,
                                target="rules",
                                maxlen=4))
summary(rules)

my.subset <- subset(rules,subset=confidence == 1 & support > 0.8)
summary(my.subset)
inspect(my.subset)

# Inspect gotgrants phenotype to see which SNPs best predict this phenotype
cross.validation(genotenureitis.filtered,'gotgrants',gain.ratio)
# chi squared
cross.validation(genotenureitis.filtered,'gotgrants',chi.squared)
# information gain
cross.validation(genotenureitis.filtered,'gotgrants',information.gain)
# symmetrical uncertainty
cross.validation(genotenureitis.filtered,'gotgrants',symmetrical.uncertainty)

# Inspect pctdrivel phenotype
# information gain ratio
cross.validation(genotenureitis.filtered,'pctdrivel',gain.ratio)
# chi squared
cross.validation(genotenureitis.filtered,'pctdrivel',chi.squared)
# information gain
cross.validation(genotenureitis.filtered,'pctdrivel',information.gain)
# symmetrical uncertainty
cross.validation(genotenureitis.filtered,'pctdrivel',symmetrical.uncertainty)
