import sys
from random import shuffle
import math
import itertools

# Load correctly formatted data. Tab-delimited, 1 column per sample, 1 row per feature
posFile = open(sys.argv[1]).readlines()	# Positive samples
negFile = open(sys.argv[2]).readlines()	# Negative samples
k = int(sys.argv[3])	# Number of neighbors to consider
p = float(sys.argv[4])	# Minimum fraction of neighbors needed to classify a sample as positive
n = int(sys.argv[5])	# Number of folds for n-fold validation

# run
# Runs the k-nearest neighbors algorithm
# Input: none
# Output: prints out parameters, accuracy, sensitivity, and specificity
#		  to console and a file called knn.out 
def run():
	pos = processData(posFile, 1)
	neg = processData(negFile, 0)
	metrics = nfold(pos, neg)
	evaluatePerformance(metrics)

# processData
# Input: data with 1 line per feature and tab-delimited samples
# Output: matrix with 1 row per sample and 1 column per feature.
# 		  The class (pos, neg) is denoted by a binary number in the first column
def processData(data, label):
	m = []
	for f in range(len(data)):
		feature = data[f].split()
		m.append(feature)
	# Swap rows and columns
	result = [list(i) for i in zip(*m)]
	for sample in result:
		sample[0] = label
	return result

# nfold
# Implements n-fold cross-validation
# Inputs: two matrices containing m rows of samples and n columns of features
# 		  the first matrix contains the positive-labeled data, the second contains negative-labeled data
# Output: a dictionary containing the total number of TP, FP, TN, and FN found across the n-folds cross-validation
def nfold(pos, neg):
	# Randomly divide and combine data into n groups
	# with equal proportions of pos and neg in each group
	posSets = divideData(n, pos)
	negSets = divideData(n, neg)
	dataSets = combineData(posSets, negSets)
	# Iterate over all sets. At each iteration treat current set s
	# as unlabeled data and other n-1 sets as labeled data
	TP=0; FP=0; TN=0; FN=0;
	for s in dataSets:
		# Merge all the other "labeled" sets
		#labeled = dataSets.copy()
		labeled = dataSets[:]
		labeled.remove(s)
		labeled = list(itertools.chain(*labeled))
		predictions = knn(s, labeled)
		answer = [i[0] for i in s]
		TP += sum(i==1 and j==1 for i,j in zip(predictions, answer))
		FP += sum(i==1 and j==0 for i,j in zip(predictions, answer))
		TN += sum(i==0 and j==0 for i,j in zip(predictions, answer))
		FN += sum(i==0 and j==1 for i,j in zip(predictions, answer))
	return {'TP':TP, 'FP':FP, 'TN':TN, 'FN':FN}

# evaluatePerformance
# Calculates accuracy, sensitivity, and specificity
# Input: a dictionary containing the total number of TP,FP,TN,FN found during n-fold cross-validation
# Output: prints out parameters, accuracy, sensitivity, and specificity
#		  to console and a file called knn.out 
def evaluatePerformance(metrics):
	output_file = open("knn.out", 'w+')
	output_file.write("k: %s\np: %.2f\nn: %s\n" % (k, p, n))
	print("k: %s\np: %.2f\nn: %s" % (k, p, n))
	
	sensitivity = float(metrics['TP']) / float(metrics['TP']+metrics['FN'])
	specificity = float(metrics['TN']) / float(metrics['TN']+metrics['FP'])
	total = sum(metrics.values())
	accuracy = float(metrics['TP']+metrics['TN']) / float(total)
	
	output_file.write("accuracy: %.2f\nsensitivity: %.2f\nspecificity: %.2f" % (accuracy, sensitivity, specificity))
	output_file.close()
	print("accuracy: %.2f\nsensitivity: %.2f\nspecificity: %.2f" % (accuracy, sensitivity, specificity))


# knn
# Implements K-nearest neighbors based on Euclidean distance
# Inputs: two matrices of labeled data and unlabeled data,
# 		  with samples as rows and features as columns
# Outputs: vector of labels where the ith entry is the label
#		   for the ith unlabeled sample
def knn(unlabeled, labeled):
	result = []
	for u in unlabeled:
		nearest = []
		
		# Compute Euclidean distances to all labeled data
		for l in labeled:
			dist = euclidean(u[1:],l[1:])
			label = l[0]
			nearest.append( [dist, label] )
			
		# Find top k Euclidean distances
		nearest = sorted(nearest, key = lambda pair: pair[0])
		nearest = nearest[0:k]
		# Label using majority vote
		numPos = sum(x[1] for x in nearest)
		if float(numPos)/float(k) >= p: result.append(1)
		else: result.append(0)
	return result

# euclidean
# Computes euclidean distance between two vectors
# Input: two vectors containing strings of numbers
# Output: the euclidean distance (a single number) between the two vectors
def euclidean(v1, v2):
	result = 0
	for i in range(len(v1)):
		result += (int(v1[i])-int(v2[i]))**2
	result = math.sqrt(result)
	return result

# divideData
# Distributes data randomly into n sets as evenly as possible
# Inputs:
#	1) Specified number of sets to divide into (n)
#	2) The data that is to be divided, where rows are samples
# Output: a list containing n datasets, where each element of the list is
#		  a dataset (matrix)
def divideData(n, data):
	result = []
	shuffle(data)
	setSize = math.trunc(len(data)/n)
	remainder = len(data)%n
	for i in range(0, len(data)-remainder, setSize):
		result.append(data[i:i+setSize])
	# Distribute remaining samples evenly among groups
	for i in range(remainder):
		result[-(1+i)].append(data[-(1+i)])
	return result

# combineData
# Combines subsets of data from two sources s1, s2 and returns
# resulting hybrid subsets in a list
# Input:
# 		1) s1, a list of datasets. Each element in this list contains a set of data (matrices)
#		2) s2, another list of datasets that is the same length as s1. However, the size of the sets of data
#			need not be the same.
# Output: a single list where the ith element is a dataset (matrix) that includes both the ith dataset
#			from s1 and the ith dataset from s2
def combineData(s1, s2):
	result = []
	for i in range(len(s1)):
		result.append(s1[i]+s2[i])
	return result

run()