"""
This program generates the pairwise Tanimoto similarity scores for all pairs of drugs
Example function call:
    tanimoto.py drugs.csv targets.csv outputfile.csv
Inputs:
    drugs.csv: file containing drug information. One row per drug.
        Columns are the DrugBank Database ID, generic name, and 
        2D fingerprint (space-delimited features)
    targets.csv: file containing a drug target information. One row per drug.
        Columns are the drug's database ID, 
        and the drug target's UniProt accession number and ID/name
    outputfile.csv: filename of the output
Output:
    csv file containing the Tanimoto similarity scores for all drug pairs
    First two columns are the DrugBank Database IDs, the third is the Tanimoto score (6-point float),
    are the fourth is a binary number indicating whether the drugs share a target (1) or not (0)
    The list is sorted by col1 and then by col 2
"""
import sys
import utilities as util
#import time
#start_time = time.time()

'''
Inputs:
    drugsfile: filename of csv containing drug information
    targetsfile: filename of csv containing a drug target information
    outputfile: filename of the output csv file
 Outputs:
    A csv file containing the Tanimoto similarity scores for all drug pairs where the
    first two columns are the DrugBank Database IDs, the third is the Tanimoto score (6-point float),
    and the fourth is a binary number indicating whether the drugs share a target (1) or not (0)
    The list is sorted by col1 and then by col 2
'''
def tanimotoScoresForAllDrugPairs(drugsfile, targetsfile, outputfile):
    fingerprints = util.readDrugData(drugsfile)
    [targetSets, ligandSets] = util.readTargetData(targetsfile, fingerprints.keys())
    util.computeAllTanimotoScores(fingerprints, targetSets, outputfile)

tanimotoScoresForAllDrugPairs(sys.argv[1],sys.argv[2],sys.argv[3])
#print("--- %s seconds --- Complete" % (time.time() - start_time))