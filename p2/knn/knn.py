import sys
from random import shuffle
import math
import itertools
"""
# Correctly formatted data. Tab-delimited, 1 column per sample, 1 row per feature
posData = open(sys.argv[1]).readlines()	# Positive samples
negData = open(sys.argv[2]).readlines()	# Negative samples
k = sys.argv[3]	# Number of neighbors to consider
p = sys.argv[4]	# Minimum fraction of neighbors needed to classify a sample as positive
n = sys.argv[5]	# Number of folds for n-fold validation
"""
	
# Correctly formatted data. Tab-delimited, 1 column per sample, 1 row per feature
posFile = open("ALL.dat").readlines()	# Positive samples
negFile = open("AML.dat").readlines()	# Negative samples
k = 4	# Number of neighbors to consider
p = 0.51	# Minimum fraction of neighbors needed to classify a sample as positive
n = 8	# Number of folds for n-fold validation
	
def run():
	print("k: %s\np: %.2f\nn: %s" % (k, p, n))
	pos = processData(posFile, 1)
	neg = processData(negFile, 0)
	metrics = nfold(pos, neg)
	evaluatePerformance(metrics)
	
# Transforms data with 1 line per feature and tab-delimited samples
# into a matrix with 1 row per sample and 1 column per feature.
# The class (pos, neg) is denoted by a binary number in the first column
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

# Implements n-fold cross-validation
# Returns the total number of TP, FP, TN, and FN in a dictionary
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
		labeled = dataSets.copy()
		labeled.remove(s)
		labeled = list(itertools.chain(*labeled))
		predictions = knn(s, labeled)
		print(predictions)
		print([i[0] for i in s])
		print("")
		answer = [i[0] for i in s]
		TP += sum(i==1 and j==1 for i,j in zip(predictions, answer))
		FP += sum(i==1 and j==0 for i,j in zip(predictions, answer))
		TN += sum(i==0 and j==0 for i,j in zip(predictions, answer))
		FN += sum(i==0 and j==1 for i,j in zip(predictions, answer))
	return {'TP':TP, 'FP':FP, 'TN':TN, 'FN':FN}

# Calculates accuracy, sensitivity, and specificity from a 
# dictionary of TP/FP/TN/FN and prints out results
def evaluatePerformance(metrics):
	sensitivity = metrics['TP'] / (metrics['TP']+metrics['FN'])
	specificity = metrics['TN'] / (metrics['TN']+metrics['FP'])
	total = sum(metrics.values())
	accuracy = (metrics['TP']+metrics['TN']) / total
	print("TP: %s, FP: %s, TN: %s, FN: %s, total: %s" % (metrics['TP'], metrics['FP'], metrics['TN'], metrics['FN'], total))
	print("accuracy: %s\nsensitivity: %s\nspecificity: %s" % (accuracy, sensitivity, specificity))

# Implements K-nearest neighbors based on Euclidean distance
# and returns labels in a vector 
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
		if numPos/k >= p: result.append(1)
		else: result.append(0)
	return result

# Computes euclidean distance between two vectors
def euclidean(v1, v2):
	result = 0
	for i in range(len(v1)):
		result += (int(v1[i])-int(v2[i]))**2
	result = math.sqrt(result)
	return result

# Distributes data randomly into n sets as evenly as possible and 
# returns datasets in a list.  
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

# Combines subsets of data from two sources and returns
# resulting hybrid subsets in a list
def combineData(s1, s2):
	result = []
	for i in range(len(s1)):
		result.append(s1[i]+s2[i])
	return result

run()