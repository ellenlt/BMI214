# This program generates file necessary to visualize a protein network
# where nodes are proteins and edges are drawn between protein nodes with
# a statistically significant bootstrap p-value (<=0.05)
# Example function call:
#     networkgen.py drugs.csv targets.csv protein_nodes.csv
# Inputs:
#     drugs.csv: file containing drug information. One row per drug.
#         Columns are the DrugBank Database ID, generic name, and 
#         2D fingerprint (space-delimited features)
#     targets.csv: file containing a drug target information. One row per drug.
#         Columns are the drug's database ID, 
#         and the drug target's UniProt accession number and ID/name
#     protein_nodes.csv: file containing protein nodes to be included in the network
#         Columns are the UniProt accession numbers, ID/name, and indications
# Outputs:
#     network.sif: links pairs of UniProt accession numbers with edges.
#                     Each row has the format "<accession#1> edge <accession#2>"
#                     accession#1 is alphabetically less than accession#2, 
#                     and the first column is a sorted list.
#     name.nodeAttr: labels the protein nodes with their UniProt name/ID. First row reads "name" and
#                     subsequent rows have the format: <accession#> = <name/ID>
#     indication.nodeAttr: labels the protein nodes with their associated indications.
#                         First row reads "indication" and subsequent rows have the format:
#                         <accession#> = <indication>

import sys
import chemoUtils as util
import numpy as np

def generateNetworkVisualizationFiles(drugsfile, targetsfile, proteinnodesfile):
    fingerprints = util.readDrugData(drugsfile)
    [targetSets, ligandSets] = util.readTargetData(targetsfile, fingerprints.keys())
    # Unsorted protein nodes data
    proteinNodes = np.loadtxt(proteinnodesfile, dtype='string', delimiter=',', skiprows=1)
    # Protein nodes data, sorted alphabetically by protein name
    sortedProteinNodes = proteinNodes[proteinNodes[:,0].argsort()]
    # Initialize data structures
    networkData=[]; nameData=[]; indicationData=[]; networkNodes=set()
    
    # Iterate over each node pair to determine if their bootstrap p-value merits an edge connection
    for i in range(len(sortedProteinNodes)):
        for j in range(i+1, len(sortedProteinNodes)):
            proteinA = sortedProteinNodes[i,0]
            proteinB = sortedProteinNodes[j,0]
            bootstrapPVal = util.calculateBootstrapPValue(100, 214, fingerprints, ligandSets, proteinA, proteinB)
            
            if bootstrapPVal <= 0.05:
                networkData.append([proteinA, proteinB])    # Create edge between connected nodes
                networkNodes.update([proteinA, proteinB])   # Update set of nodes in network
    
    # Annotate each node with its name and indication
    for protein in proteinNodes:
        if protein[0] in networkNodes:
            nameData.append([protein[0],protein[1]])
            indicationData.append([protein[0],protein[2]])
    
    # Output files        
    np.savetxt('network.sif', np.array(networkData), fmt='%s edge %s')
    np.savetxt('name.nodeAttr', np.array(nameData), fmt='%s = %s', header="name", comments='')
    np.savetxt('indication.nodeAttr', np.array(indicationData), fmt='%s = %s', header="indication", comments='')
            
    
generateNetworkVisualizationFiles(sys.argv[1],sys.argv[2],sys.argv[3])
