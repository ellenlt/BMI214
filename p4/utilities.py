import csv
import numpy as np
import random
import time
start_time = time.time()

# Reads in a csv file into a Numpy array 
# Input: filename of csv file
# Output: Numpy array of strings. Omits header line of csv file
def readCsvData(filename):
    targets = np.loadtxt(filename, dtype='string', delimiter=',', skiprows=1)
    return targets

# Reads in a csv file of drug information into a Python 2D list
# Input: filename of csv file containing drug information. One row per drug.
#        Columns are the DrugBank Database ID, generic name, and 
#        2D fingerprint (space-delimited features)
# Output: drugIDs: Numpy array of strings containing the DrugBank Database ID
#         fingerprints: Python integers list storing the fingerprints
def readDrugData(filename):
    drugIDs = []
    fingerprints = []
    with open(filename, 'rb') as drugsCsv:
        reader = csv.reader(drugsCsv, delimiter=',')
        next(reader, None)
        for row in reader:
            #fingerprints.append([int(i) for i in row[2].split()])
            fingerprints.append(set(row[2].split()))
            drugIDs.append(row[0])
    return [np.array(drugIDs), fingerprints]

# Computes the Tanimoto score/coefficient for two molecules
# Input: two lists of integers that represent the fingerprints of the two molecules
# Output: the Tanimoto coefficient as a float
def tanimotoScore(fingerprint1, fingerprint2):
    #intersection = list(set(fingerprint1) & set(fingerprint2))
    #union = list(set(fingerprint1) | set(fingerprint2))
    score = float(len(fingerprint1 & fingerprint2))/float(len(fingerprint1 | fingerprint2))
    return score


# Calculates Tanimoto scores for all drugs pairs and prints to a csv file
# Inputs: 
#        drugs: a 2D Python list containing drug information. One row per drug.
#               Columns are the DrugBank Database ID, generic name, and 
#               2D fingerprint (space-delimited features)
#        targets: a Numpy array of strings containing drug target information, where
#                the columns are the drug database ID, target UniProt accession number and target ID/name.
#        outputfilename: optional parameter; filename of the output csv file
# Output: 
#        If outputfilename is specified,
#        writes to a csv file the Tanimoto similarity scores for all drug pairs where the
#        first two columns are the DrugBank Database IDs, the third is the Tanimoto score,
#        and the fourth column indicates whether the drugs share a target or not (1 or 0)
#        Returns the first three columns as a Numpy array of strings
def computeAllTanimotoScores(drugIDs, fingerprints, targets, outputfilename=False):
    tanimotoScores = []
    #if outputfilename != False: outputFile = open(outputfilename, 'w+')
    #possibleDrugPairs = []
    pasttime=time.time()
    for drugIdx1 in range(len(drugIDs)):
        #print("--- %s seconds --- New drug and all its combos" % (time.time() - pasttime))
        pasttime=time.time()
        for drugIdx2 in range(drugIdx1+1,len(drugIDs)):
            #print("--- %s seconds --- Comparing new pair of drugs" % (time.time() - start_time))
            id1 = drugIDs[drugIdx1]
            id2 = drugIDs[drugIdx2]
            score = tanimotoScore(fingerprints[drugIdx1], fingerprints[drugIdx2])
            #print("--- %s seconds --- Finished computing tanimoto score" % (time.time() - start_time))
            #if outputfilename != False: outputFile.write("%s,%s,%.6f,%d\n" % (id1, id2, score, sharedTarget(id1, id2, targets)))  
            #tanimotoScores.append([id1, id2, score])
            tanimotoScores.append([id1, id2, format(score,'.6f'), sharedTarget(id1, id2, targets)]) 
    tanimotoScores = np.array(tanimotoScores)
    if outputfilename != False: np.savetxt(outputfilename, tanimotoScores, fmt='%s', delimiter=',')
    #if outputfilename != False: outputFile.close()
    return np.array(tanimotoScores)


# Determines whether two drugs have a shared target
# Input: id1 and id2: the two DrugBank Database IDs as strings, as well as
#        targets: a Numpy array of strings containing drug target information, where
#                the columns are the drug database ID, target UniProt accession number and target ID/name.
# Output: 0 if no shared target, 1 if shared target
def sharedTarget(id1, id2, targets):
    targets1 = targets[targets[:,0] == id1][:,1]
    targets2 = targets[targets[:,0] == id2][:,1]
    sharedTargets = np.intersect1d(targets1, targets2, assume_unique=True)
    if len(sharedTargets) > 0: return 1
    else: return 0


# Calculates bootstrap p-value for two proteins
# Input: 
#        n: number of iterations
#        r: seed for random number generator
#        drugs: Python 2D list containing drug information. One row per drug.
#            Columns are the DrugBank Database ID, generic name, and 2D fingerprint (space-delimited features)
#        targets: Numpy array of strings containing drug target information, where
#                the columns are the drug database ID, target UniProt accession number and target ID/name.
#        proteinA and proteinB: UniProt accession numbers as strings
#        tanimotoScores: optional parameter; Numpy array of strings containing precomputed Tanimoto scores 
#                        for drug pairs where the first two columns are the DrugBank Database IDs
#                        and the third is the Tanimoto score
# Output: bootstrap p-value as a float
def calculateBootstrapPValue(n, r, drugIDs, fingerprints, targets, proteinA, proteinB, tanimotoScores=None):
    ligandSetA = findLigands(proteinA, targets)
    ligandSetB = findLigands(proteinB, targets)
    realSummaryScore = computeTanimotoSummaryScore(ligandSetA, ligandSetB, drugIDs, fingerprints, tanimotoScores)
    
    count = 0   # The number of times a random summary Tanimoto score exceeds the observed
    random.seed(r)
    #drugsList = [row[0] for row in drugs]   # List of all drug IDs from which to sample
    
    for i in range(n):
        randomSetA = random.sample(drugIDs, len(ligandSetA))
        randomSetB = random.sample(drugIDs, len(ligandSetB))
        randomSummaryScore = computeTanimotoSummaryScore(randomSetA, randomSetB, drugIDs, fingerprints, tanimotoScores)
        if randomSummaryScore >= realSummaryScore: count = count + 1
    bootstrapPVal = float(count)/float(n)
    return bootstrapPVal

# Computes the Tanimoto summary score for two proteins to assess their similarity
# Inputs:
#        ligandSetA and ligandSetB: two lists of the DrugBank Database IDs (strings) of all drugs that
#                                    bind to proteins A and B, respectively
#        drugs: Python 2D list containing the same information as above. First two columns are
#                 strings while third column stores fingerprint as a vector of integers
#        tanimotoScores: optional parameter; Numpy array of strings containing precomputed Tanimoto scores 
#                        for drug pairs where the first two columns are the DrugBank Database IDs
#                        and the third is the Tanimoto score
# Output: Tanimoto summary score (float)
def computeTanimotoSummaryScore(ligandSetA, ligandSetB, drugIDs, fingerprints, tanimotoScores=None):
    summaryScore = 0.0    # Initialize summary score to be 0
    # If no precomputed scores available, process drugs data for easy computation later on
    #if tanimotoScores is None:
    #    drugsIDs = np.array([row[0] for row in drugs])
    #    fingerprints = [row[2] for row in drugs]
    # Compare all ligands in set A to those in set B
    for ligandA in ligandSetA:
        for ligandB in ligandSetB:
            # To save computation time, hard code score if ligands are identical
            if ligandA == ligandB: score = 1.0
            # If precomputed score is available, use it
            elif tanimotoScores is not None:
                print 1
                score = float(tanimotoScores[((tanimotoScores[:,0] == ligandA) & (tanimotoScores[:,1] == ligandB)) | ((tanimotoScores[:,0] == ligandB) & (tanimotoScores[:,1] == ligandA))][:,2][0])
            else:
                ligandAIdx = np.where(drugIDs==ligandA)[0]
                ligandBIdx = np.where(drugIDs==ligandB)[0]
                score = tanimotoScore(fingerprints[ligandAIdx], fingerprints[ligandBIdx])
            if score > 0.5: summaryScore = summaryScore + score
    return summaryScore

# Finds all drugs that bind to a given protein (the protein's ligand set)
# Input: proteinID: UniProt accession number as a string
#        targets: Numpy array of strings containing drug target information, where
#                the columns are the drug database ID, target UniProt accession number and target ID/name.
# Output: a list of DrugBank Database IDs (strings) of all drugs that bind to the given protein
def findLigands(proteinID,targets):
    ligands = targets[targets[:,1]==proteinID][:,0]
    return ligands