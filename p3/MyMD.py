import numpy as np
import scipy.spatial.distance as dist
import argparse

# Input: none
# Output: none
# Initializes all parameters
def initParams():
    parser = argparse.ArgumentParser(description='Process inputs to program')
    parser.add_argument('--iF', metavar='iF', type=str, help='Input filename, including directory and extension')
    parser.add_argument('--kB', metavar='kB', type=float, help='Spring constant for bonds', default='4000.0')
    parser.add_argument('--kN', metavar='kN', type=float, help='Spring constant for non-bonds', default='400.0')
    parser.add_argument('--nbCutoff', metavar='nbCutoff', type=float, help='Distance within which non-bonded atoms have interactions', default='0.50')
    parser.add_argument('--m', metavar='m', type=float, help='Atom mass for all atoms', default='12.0')
    parser.add_argument('--dt', metavar='dt', type=float, help='Length of timestep', default='0.001')
    parser.add_argument('--n', metavar='n', type=int, help='Number of timesteps to iterate', default='1000')
    parser.add_argument('--out', metavar='out', type=str, help='Prefix of the filename to output')
    
    args = parser.parse_args()
    
    global iF; global kB; global kN; global nbCutoff; global m; global dt; global n; global out
    # Input filename, including directory and file extension. Input file should be an .rvc file containing
    # a header row as the first row. The rest of the rows represent atoms.
    # The columns represent the atom's ID (1-indexed, in sequential order),
    # x, y, z coordinates, x, y, z velocities, and ID's of all connected atoms
    iF = args.iF
    kB = args.kB    # Spring constant for bonds
    kN = args.kN      # Spring constant for nonbonds
    nbCutoff = args.nbCutoff # Distance within which non-bonded atoms have interactions
    m = args.m        # Atom mass for all atoms
    dt = args.dt      # Length of timestep
    n = args.n        # Number of timesteps to iterate
    if args.out is None: out = iF.split('.rvc')[0]  # Prefix of filename to output
    else: out = args.out

def generateSimulation():
    initParams()
    molecule = initMolecule()
    # numpy vector for force and energy
    # Iterate through n timesteps
        # Initialize force and energy vectors to 0
        # Update velocities (t+dt/2) on each atom
        # Update positions
        # For bonded pairs and nonbonded pairs:
            # Calculate and update energy
            # Calculate forces and update in each dimension
        # Update velocities (t+dt) and calculate kinetic energy
    # Every 1 timesteps, calculate E_tot = KE + PEbonds + PEnonbonds as sanity check
    # Every 10 timesteps, write appropriate output to rvc, crd, files needed to generate euc 
    # Output: .erg and .rvc files
    

def initMolecule():
    global iF; global nbCutoff; global positions; global velocities; global forces
    data = open(iF).readlines()
    numAtoms = len(data)-1
    
    # Three Numpy matrices containing the spatial positions, velocities, and total forces on each atom.
    # Each row represents an atom; the 3 columns contain the x, y, and z coordinates/velocities/forces
    positions = np.empty([numAtoms, 3])
    velocities = np.empty([numAtoms, 3])
    forces = np.zeros([numAtoms, 3])   # Set initial forces to be zero
    # Two 2D lists (ultimately converted to Numpy matrices) containing information for bonded and nonbonded atom pairs.
    # Each row represents a bond (or nonbond). The 6 columns contain the ID of atom 1, the ID of atom 2 (the lower ID is always first),
    # distance between the two atoms, reference/initial distance between the two atoms,
    # the potential energy of the bond (or nonbond), and the force between the bonded (or nonbonded) atoms
    bonds = []
    nonbonds = []
    
    # Boolean Numpy array keeps track of atoms whose positions and velocities have already been stored
    atomHasAlreadyBeenStored = np.zeros(numAtoms, dtype=bool)
    
    for currAtomID in range(1,numAtoms+1):
        if not atomHasAlreadyBeenStored[currAtomID-1]:
            connectivities = storePositionAndVelocity(currAtomID, data)
            atomHasAlreadyBeenStored[currAtomID-1] = True
        for connectedAtomID in connectivities:
            # Only need to store each bond/nonbonded pair of atoms once
            if currAtomID < connectedAtomID:
                # Classify them as bonded/nonbonded and store their info in the appropriate data-structure
                storePositionAndVelocity(connectedAtomID, data)
                atomHasAlreadyBeenStored[connectedAtomID-1] = True
                distance = dist.euclidean(positions[currAtomID], positions[connectedAtomID])
                if distance > nbCutoff
                #if(distance)
                
    # numpy array for positions
    # numpy array for velocities
    # store reference bond distance for each pair of atams
    # numpy vector with indices of bonded and nonbonded pairs
    # you know initial positions and velocities
    # find and store nonbonded pairs. future iterations: only bonded atoms have bonded interactions.
    # nonbonded either don't interact or have nonbonded interactions
    # 
    
# Input: The ID number of an atom, and a 1D array of strings where each row is a row of the input file
# Output: Stores the atom's x, y, and z position and velocities as floats in global Numpy matrices ("positions"
#        and "velocities"). Returns an array of the ID's of connected atoms as integers.
def storePositionAndVelocity(atomID, data):
    global positions; global velocities
    atomInfo = data[atomID].split()
    positions[atomID-1] = map(float, atomInfo[1:4])
    velocities[atomID-1] = map(float, atomInfo[4:7])
    return map(int, atomInfo[7:])

generateSimulation()