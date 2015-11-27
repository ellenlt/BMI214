'''
This program generates file necessary to visualize a protein network
where nodes are proteins and edges are drawn between protein nodes with
a statistically significant bootstrap p-value (<=0.05)
Example function call:
    networkgen.py drugs.csv targets.csv protein_nodes.csv
Inputs:
    drugs.csv: file containing drug information. One row per drug.
        Columns are the DrugBank Database ID, generic name, and 
        2D fingerprint (space-delimited features)
    targets.csv: file containing a drug target information. One row per drug.
        Columns are the drug's database ID, 
        and the drug target's UniProt accession number and ID/name
    protein_nodes.csv: file containing protein nodes to be included in the network
        Columns are the UniProt accession numbers, ID/name, and indications
Outputs:
    network.sif: links pairs of UniProt accession numbers with edges.
                    Each row has the format "<accession#1> edge <accession#2>"
                    accession#1 is alphabetically less than accession#2, 
                    and the first column is a sorted list.
    name.nodeAttr: labels the protein nodes with their UniProt name/ID. First row reads "name" and
                    subsequent rows have the format: <accession#> = <name/ID>
    indication.nodeAttr: labels the protein nodes with their associated indications.
                        First row reads "indication" and subsequent rows have the format:
                        <accession#> = <indication>
'''
import sys
import utilities as util
import time
start_time = time.time()

def run(drugsfile, targetsfile, proteinnodesfile):
    drugs = util.readDrugData(drugsfile)
    targets = util.readCsvData(targetsfile)
    proteinNodes = util.readCsvData(proteinnodesfile)
    # Sort protein names alphabetically
    proteinNodes = proteinNodes[proteinNodes[:,0].argsort()]
    
    # Precompute all Tanimoto scores to save computation time
    tanimotoScores = util.computeAllTanimotoScores(drugs, targets)
    
    # Initialize output files and write necessary headers
    sifFile = open('network.sif', 'w+')
    nameFile = open('name.nodeAttr', 'w+'); nameFile.write("name\n")
    indicationFile = open('indication.nodeAttr', 'w+'); indicationFile.write("indication\n")
    
    for i in range(len(proteinNodes)-1):
        for j in range(i+1, len(proteinNodes)):
            proteinA = proteinNodes[i,0]
            proteinB = proteinNodes[j,0]
            p = util.calculateBootstrapPValue(100, 214, drugs, targets, proteinA, proteinB, tanimotoScores)
            if p<= 0.05:
                sifFile.write("%s edge %s\n" % (proteinA, proteinB))
                nameFile.write("%s = %s\n" % (proteinA, proteinNodes[i,1]))
                nameFile.write("%s = %s\n" % (proteinB, proteinNodes[j,1]))
                indicationFile.write("%s = %s\n" % (proteinA, proteinNodes[i,2]))
                indicationFile.write("%s = %s\n" % (proteinB, proteinNodes[j,2]))
    
    sifFile.close()
    nameFile.close()    
    indicationFile.close()
            
    
run(sys.argv[1],sys.argv[2],sys.argv[3])
print("--- %s seconds ---" % (time.time() - start_time))