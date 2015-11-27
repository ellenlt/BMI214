'''
This program generates the bootstrap p-value for the comparison of two proteins.
Example function call:
    pvalue.py -n <INT> -r <INT> <drugs.csv> <targets.csv> <proteinA> <proteinB>
Inputs:
    n: number of bootstrap iterations (default = 100)
    r: seed for pseudo-random number generator (default = 214)
    drugs.csv: file containing drug information. One row per drug.
        Columns are the DrugBank Database ID, generic name, and 
        2D fingerprint (space-delimited features)
    targets.csv: file containing a drug target information. One row per drug.
        Columns are the drug's database ID, 
        and the drug target's UniProt accession number and ID/name
    proteinA and proteinB: UniProt accession numbers of the two proteins that are being compared
Output:
    Bootstrap p-value as a float
'''

import utilities as util
import argparse
#import time
#start_time = time.time()

# Prints the bootstrap p value for two proteins
# Input: command line arguments (see header comments)
# Output: bootstrap p value as a float
def printBootstrapPValue():
    params = initParams()
    print util.calculateBootstrapPValue(params['n'], params['r'],
                                   params['fingerprints'], params['ligandSets'],
                                   params['proteinA'], params['proteinB'])

# Initializes all parameters needed
# Input: none; uses arguments supplied to the command line (see header comments)
# Output: a dictionary of parameters containing:
#        n: number of iterations
#        r: seed for random number generator
#        fingerprints: dictionary that maps DrugBank Database IDs -> drug fingerprint (set of strings)
#        ligandSets: dictionary mapping protein/drug target accession numbers (strings) ->
#                    ligand sets or drugs that they bind to (the drug database ID as a string)
#        proteinA and proteinB: UniProt accession numbers as strings
def initParams():
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', metavar='n', type=int, default='100')
    parser.add_argument('-r', metavar='r', type=int, default='214')
    parser.add_argument('drugs', type=str)
    parser.add_argument('targets', type=str)
    parser.add_argument('proteinA', type=str)
    parser.add_argument('proteinB', type=str)
    
    args = parser.parse_args()
    fingerprints = util.readDrugData(args.drugs)
    [targetSets, ligandSets] = util.readTargetData(args.targets, fingerprints.keys())
    # Return a dictionary with all parameters
    return {'n':args.n,'r':args.r, 
            'fingerprints':fingerprints, 'ligandSets':ligandSets, 
            'proteinA':args.proteinA, 'proteinB':args.proteinB}

printBootstrapPValue()
#print("--- %s seconds ---" % (time.time() - start_time))