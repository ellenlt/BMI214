import sys
import random
import math

# Correctly formatted data. Tab-delimited, 1 row per gene, 1 column per experimental condition
k = int(sys.argv[1])    # Number of centroids
dataFile = open(sys.argv[2]).readlines()    # Expression data
maxIters = int(sys.argv[3])    # Max number of iterations
# Optional parameter specifying a file with predetermined centroids
if len(sys.argv) > 4:
    centroidsFile = open(sys.argv[4]).readlines()
    if k > len(centroidsFile): k = len(centroidsFile)   # Choose first k of the provided centroids 

"""
# Correctly formatted data. Tab-delimited, 1 row per sample, 1 column per condition
k = 3    # Number of centroids
dataFile = open("test.dat").readlines()    # Data
maxIters = 100    # Max number of iterations
#centroidsFile = open("yeast_centroids.txt").readlines()    # Optional parameter specifying file with predetermined centroids
"""

def run():
    data = processData(dataFile)
    # Initialize k centroids in matrix
    if 'centroidsFile' in globals(): centroids = processData(centroidsFile[:k])
    else: centroids = generateCentroids(data)
    # Assign genes to clusters
    results = kmeans(data, centroids)
    printResults(results)

def printResults(results):
    output_file = open("kmeans.out", 'w+')
    print("iterations: %d" % results['iterations'])
    assignments = results['assignments']
    for a in range(len(assignments)):
        output_file.write("%s\t%s\n" % (a+1, assignments[a]))
    output_file.close()

# kmeans
# Implements k-means clustering based on Euclidean distance
# Input: 1) data with samples as rows and conditions/features as columns,
#        2) list of initial centroids
# Output: a dictionary containing:
#     1) A list, where the first element is the number of iterations i
#         (after the ith iteration, clusters don't change), and the
#         other elements i are cluster assignments for sample i.
#     2) The number of iterations it took to reach convergence.
def kmeans(data, centroids):
    clusterAssignments = []
    iterations = 0
    # Reassign genes to clusters and recalculate centroids until
    # you reach either convergence, or max number of iterations
    while True:
        newClusterAssignments = assignClusters(data, centroids)
        if newClusterAssignments == clusterAssignments or iterations == maxIters:
            clusterAssignments = newClusterAssignments
            break
        clusterAssignments = newClusterAssignments
        centroids = updateCentroids(clusterAssignments, data)
        iterations += 1
    return {'iterations': iterations, 'assignments': clusterAssignments}

# assignClusters
# Assigns each sample to a cluster based on closest Euclidean distance
# Input: 1) data with samples as rows and conditions/features as columns,
#        2) list of initial centroids
# Output: vector with numeric cluster assignments for each sample, in the
#         same order as the samples appear in the data
def assignClusters(data, centroids):
    result = []
    for sample in data:
        # Find closest centroid to one particular sample
        minDist = 10000000000000000000
        for centroid in centroids:
            currDist = euclidean(centroid, sample)
            # Update minimum distance and closest cluster
            if currDist < minDist:
                minDist = currDist
                cluster = centroids.index(centroid)+1   # Number first centroid as 1
        # Store in order
        result.append(cluster)
    return result

# updateCentroids
# Updates centroid locations to be the average of a cluster.
# Input: 1) vector containing assigned cluster numbers for each sample, in same order
#            as the appearance of the samples in the data
#        2) data with samples as rows and conditions/features as columns
# Output: vector of updated centroids
def updateCentroids(clusterAssignments, data):
    result = []
    # For each cluster, update centroid
    for c in range(1,k+1):
        # Aggregates indices of all samples in cluster
        indicesOfClusterMembers = [i for i, assignment in enumerate(clusterAssignments) if assignment==c]
        # Cluster containing samples at specified indices
        cluster = [data[i] for i in indicesOfClusterMembers]
        # Compute average of samples in cluster
        clusterCenter = [float(sum(col))/len(col) for col in zip(*cluster)]
        result.append(clusterCenter)
    return result

# generateCentroids
# Input: data with samples as rows and conditions/features as columns
# Output: list of k randomly generated centroids within the range of the data
def generateCentroids(data):
    centroids = [[] for c in range(k)]
    # For each experimental condition... 
    for e in range(len(data[0])):
        # ...find range of values
        maxVal = int(max([i[0][e] for i in zip(data)]))
        minVal = int(min([i[0][e] for i in zip(data)]))
        # ...and generate k random values (1 per centroid)
        for c in centroids:
            c.append(random.randrange(minVal, maxVal+1))
    return centroids

# processData
# Input: data with 1 row per sample and tab-delimited features
# Output: matrix with 1 row per sample and 1 column per feature.
def processData(data):
    result = []
    for i in range(len(data)):
        sample = data[i].split()
        sample = [float(x) for x in sample]
        result.append(sample)
    return result

# Computes Euclidean distance between two vectors
def euclidean(v1, v2):
    result = 0
    for i in range(len(v1)):
        result += (v1[i]-v2[i])**2
    result = math.sqrt(result)
    return result

run()
    
    