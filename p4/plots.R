setwd("/Users/ellen/BMI214/p4/")

all_tanimoto <- read.csv("~/BMI214/p4/output.csv", header=F)
shared_tanimoto <- all_tanimoto[all_tanimoto$V4==1,]
notshared_tanimoto <- all_tanimoto[all_tanimoto$V4==0,]

hist(all_tanimoto$V3,main="05708608 All",xlab="Tanimoto Score")
hist(shared_tanimoto$V3,main="05708608 Shared",xlab="Tanimoto Score")
hist(notshared_tanimoto$V3,main="05708608 Not Shared",xlab="Tanimoto Score")
