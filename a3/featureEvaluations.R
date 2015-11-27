#Load the required libraries
require(e1071)
require(FSelector)

#Load and examine the input data
#Replace PATH with your path
featureData<-read.csv("~/BMI214/a3/features.csv")
print(head(featureData))

#Select testing and training data sets using an 80% training, 20% testing split
#NOTE: you must run the following three lines of code in EXACTLY this order to guarantee the right answer.
#				Otherwise, your random samples will differ from the our solutions
set.seed(1028)
trainingRows <- sample(1:nrow(featureData),nrow(featureData)*.8)
testingRows <- which(!c(1:nrow(featureData)) %in% trainingRows)

#Fit a Naive Bayes model using the training data
nb.fit <- naiveBayes(SITE~., data=featureData[trainingRows,])

# Generate the confusion matrix for the Naive Bayes data
#		The rows of the output are the predictions, the columns are the actual data
print(table(predict(nb.fit,featureData[testingRows,names(featureData)!="SITE"]), featureData[testingRows,"SITE"]))

#Examine the information gain
info.gain <- information.gain("SITE ~ .", featureData)

#Rank the attributes by information gain
info.gain.ordered <- data.frame("Attribute importance"=info.gain[order(-info.gain$attr_importance),], 
																"Name"=rownames(info.gain)[order(-info.gain$attr_importance)])

#Examine the attributes
print(info.gain.ordered)
