# This file contains utility functions for tanimoto.py, pvalue.py, and networkgen.py

import csv
import numpy as np
import random

# Reads / stores a csv file containing drug target information
# Inputs: 
#       targets.csv: file containing a drug target information. One row per drug.
#           Columns are the drug's database ID, 
#           and the drug target's UniProt accession number and ID/name
#        drugIDs: array of strings containing the DrugBank Database ID
# Outputs:
#        targetSets: dictionary mapping drug database IDs (strings) -> sets of
#                    drug targets/proteins (the accession number as a string)
#        ligandSets: dictionary mapping protein/drug target accession numbers (strings) ->
#                    ligand sets or drugs that they bind to (the drug database ID as a string)
def readTargetData(filename, drugIDs):
    data = np.loadtxt(filename, dtype='string', delimiter=',', skiprows=1)
    targetSets = {}
    for drugID in drugIDs:
        targetSets[drugID]=set(data[data[:,0] == drugID][:,1])
    ligandSets = {}
    ligandAccessionNums = list(np.unique(data[:,1]))
    for ligandAccessionNum in ligandAccessionNums:
        ligandSets[ligandAccessionNum] = set(data[data[:,1] == ligandAccessionNum][:,0])
    return [targetSets, ligandSets]
    
# Reads / stores a csv file of drug information
# Input:
#        filename: filename of csv file containing drug information. One row per drug.
#                Columns are the DrugBank Database ID, generic name, and 
#                2D fingerprint (space-delimited features)
# Output:
#        fingerprints: dictionary that maps each drug's database ID -> its fingerprint (set of strings)
def readDrugData(filename):
    fingerprints = {}
    with open(filename, 'rb') as drugsCsv:
        reader = csv.reader(drugsCsv, delimiter=',')
        next(reader, None)
        for row in reader:
            fingerprints[row[0]]=set(row[2].split())
    return fingerprints

# Computes the Tanimoto score/coefficient for two molecules
# Input: two sets of strings that represent the fingerprints of the two molecules
# Output: the Tanimoto coefficient as a float
def tanimotoScore(fingerprint1, fingerprint2):
    intersectionSize = len(fingerprint1 & fingerprint2)
    unionSize = len(fingerprint1) + len(fingerprint2) - intersectionSize
    score = float(intersectionSize)/float(unionSize)
    return score


# Calculates Tanimoto scores for all drugs pairs and prints to a csv file
# Inputs: 
#        fingerprints: dictionary that maps each drug's database ID -> its fingerprint (set of strings)
#        targetSets: dictionary mapping drug database IDs (strings) -> sets of
#                    drug targets/proteins (the accession number as a string)
#        outputfilename: optional parameter; filename of the output csv file
# Output: 
#        If outputfilename is specified, writes the Tanimoto similarity scores for all drug pairs
#            to a csv file. The first two columns are the DrugBank Database IDs, 
#            the third is the Tanimoto score (6-point precision), and the fourth column indicates 
#            whether the drugs share a target or not (1 or 0). Sorted alphabetically by col1 then col2.
#            Returns the above information as a Numpy array of strings
#        If no outputfilename is provided, returns the above information as a Numpy array of strings
#            but omitting the 4th column and with unspecified float precision
def computeAllTanimotoScores(fingerprints, targetSets, outputfilename=False):
    tanimotoScores = []
    drugIDs = fingerprints.keys()
    drugIDs.sort()
    for drugIdx1 in range(len(drugIDs)):
        for drugIdx2 in range(drugIdx1+1,len(drugIDs)):
            drugID1 = drugIDs[drugIdx1]
            drugID2 = drugIDs[drugIdx2]
            score = tanimotoScore(fingerprints[drugID1], fingerprints[drugID2])
            if outputfilename is False:
                tanimotoScores.append([drugID1, drugID2, score])
            else:
                target = sharedTarget(drugID1, drugID2, targetSets)
                tanimotoScores.append([drugID1, drugID2, format(score,'.6f'), target]) 
    tanimotoScores = np.array(tanimotoScores)
    if outputfilename != False: np.savetxt(outputfilename, tanimotoScores, fmt='%s', delimiter=',')
    return tanimotoScores


# Determines whether two drugs have a shared target
# Input: 
#        id1 and id2: the two DrugBank Database IDs as strings
#        targetSets: dictionary mapping drug database IDs (strings) -> sets of
#                    drug targets/proteins (the accession number as a string)
# Output: 0 if no shared target, 1 if shared target
def sharedTarget(id1, id2, targetSets):
    if len(targetSets[id1] & targetSets[id2])>0: return 1
    else: return 0


# Calculates bootstrap p-value for two proteins
# Input: 
#        n: number of iterations
#        r: seed for random number generator
#        fingerprints: dictionary that maps each drug's database ID -> its fingerprint (set of strings)
#        ligandSets: dictionary mapping protein/drug target accession numbers (strings) ->
#                    ligand sets or drugs that they bind to (the drug database ID as a string)
#        proteinA and proteinB: UniProt accession numbers as strings
# Output: bootstrap p-value as a float
def calculateBootstrapPValue(n, r, fingerprints, ligandSets, proteinA, proteinB):
    ligandSetA = ligandSets[proteinA]
    ligandSetB = ligandSets[proteinB]
    realSummaryScore = computeTanimotoSummaryScore(ligandSetA, ligandSetB, fingerprints)
    count = 0   # The number of times a random summary Tanimoto score exceeds the observed
    random.seed(r)
    drugIDs = fingerprints.keys()
    drugIDs.sort()
    
    for i in range(n):
        randomSetA = random.sample(drugIDs, len(ligandSetA))
        randomSetB = random.sample(drugIDs, len(ligandSetB))
        randomSummaryScore = computeTanimotoSummaryScore(randomSetA, randomSetB, fingerprints)
        if randomSummaryScore >= realSummaryScore: count = count + 1
    bootstrapPVal = float(count)/float(n)
    return bootstrapPVal

# Computes the Tanimoto summary score for two proteins by comparing the ligand sets of the
# two proteins and summing the Tanimoto scores for all ligands that exceed 0.5
# Inputs:
#        ligandSetA and ligandSetB: the ligand sets of proteins A and B, respectively.
#                                   aka, the DrugBank Database IDs (strings) of
#                                   all drugs that bind each protein
#        fingerprints: dictionary that maps each drug's database ID -> its fingerprint (set of strings)
# Output: Tanimoto summary score (float)
def computeTanimotoSummaryScore(ligandSetA, ligandSetB, fingerprints):
    summaryScore = 0.0    # Initialize summary score to be 0
    # Compare all ligands in set A to those in set B
    for ligandA in ligandSetA:
        for ligandB in ligandSetB:
            # To save computation time, hard code score if ligands are identical
            if ligandA == ligandB: score = 1.0
            else: score = tanimotoScore(fingerprints[ligandA], fingerprints[ligandB])
            if score > 0.5: summaryScore = summaryScore + score
    return summaryScore

