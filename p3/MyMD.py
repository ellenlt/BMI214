import numpy as np
import scipy.spatial.distance as dist
import argparse

# initParams gets and processes the desired parameters from the console, and sets defaults as needed
# Input: none
# Output: a dictionary containing:
#        - iF: Input filename, including directory and extension (string)
#        - kB: Spring constant for bonds (float)
#        - kN: Spring constant for non-bonds (float)
#        - nbCutoff: Distance within which non-bonded atoms have interactions (float)
#        - m: Atom mass for all atoms (float)
#        - dt: Length of timestep (float)
#        - n: Number of timesteps to iterate (int)
#        - out: Prefix of the filename to output (string)
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
    
    # Input filename, including directory and file extension. Input file should be an .rvc file containing
    # a header row as the first row. The rest of the rows represent atoms.
    # The columns represent the atom's ID (1-indexed, in sequential order),
    # x, y, z coordinates, x, y, z velocities, and ID's of all connected atoms
    iF = args.iF
    kB = args.kB                # Spring constant for bonds
    kN = args.kN                # Spring constant for nonbonds
    nbCutoff = args.nbCutoff    # Distance within which non-bonded atoms have interactions
    m = args.m                  # Atom mass for all atoms
    dt = args.dt                # Length of timestep
    n = args.n                  # Number of timesteps to iterate
    if args.out is None: out = iF.split('.rvc')[0]  # Prefix of filename to output
    else: out = args.out
    
    # Return a dictionary with all parameters
    return {'iF':iF,'kB':kB, 'kN':kN, 'nbCutoff':nbCutoff, 'm':m, 'dt':dt, 'n':n, 'out':out}


def generateSimulation():
    params = initParams()
    initialInfo = initMolecule(params)
    stepThruSimulation(params, initialInfo)
    
# Inputs:
#    1) a dictionary containing:
#        - kB: Spring constant for bonds (float)
#        - kN: Spring constant for non-bonds (float)
#        - m: Atom mass for all atoms (float)
#        - dt: Length of timestep (float)
#        - n: Number of timesteps to iterate (int)
#        - out: Prefix of the filename to output (string)
#    2) dictionary containing the following data structures:
#         - positions: nx3 Numpy ndarray containing the x, y, and z coordinates of n atoms
#         - velocities: nx3 Numpy ndarray containing the x, y, and z velocities of n atoms
#         - bonds: nx9 Numpy ndarray where each row represents a bond and the columns are:
#                    0) ID of atom 1 (the lower ID # is first)
#                    1) ID of atom 2
#                    2) bond length
#                    3) reference/initial bond length
#                    4) bond potential energy
#                    5) magnitude of force between bonded atoms
#                    6) x-component of force exerted on atom 1 by atom 2
#                    7) y-component of force exerted on atom 1 by atom 2
#                    8) z-component of force exerted on atom 1 by atom 2
#         - nonbonds: identical to "bonds" above, except values are for non-bonded atom pairs
#         - header: first line of the data input file
def stepThruSimulation(params, molecule):
    # Unpack dictionary inputs into local variables
    positions=molecule['positions']; velocities=molecule['velocities']; bonds=molecule['bonds'];nonbonds=molecule['nonbonds']
    dt=params['dt']; n=params['n']; kB=params['kB']; kN=params['kN']; m=params['m'];out=params['out']
    
    numAtoms = len(velocities)
    # nx3 Numpy ndarrays containing the x, y, and z components of force and acceleration for all n atoms
    forces = np.zeros([numAtoms, 3]); accelerations = np.zeros([numAtoms, 3])
    
    # Create output files
    ergFile = open(out + ".erg", 'w+')
    ergFile.write("# step\tE_k\tE_b\tE_nB\tE_tot")
    rvcFile = open(out + ".rvc", 'w+')
    rvcFile.write(molecule['header'])
    
    for i in range(n):
        # Update velocities at t+0.5dt for each atom
        velocities = velocities + accelerations*0.5*dt  
        # Update positions
        positions = positions + velocities*dt + accelerations*dt**2 
        # Update potential energies for each bond and nonbond
        bonds[:,4] = 0.5*kB*(bonds[:,2]-bonds[:,3])**2
        nonbonds[:,4] = 0.5*kB*(nonbonds[:,2]-nonbonds[:,3])**2
        # Calculate the magnitude of the force for each bond and nonbond
        bonds[:,5] = kB*(bonds[:,2]-bonds[:,3])
        nonbonds[:,5] = kB*(nonbonds[:,2]-nonbonds[:,3])
        # Translate bond forces into x, y, and z components (force exerted on the first atom in the bond by the second atom)
        pos1 = positions[map(int, bonds[:,0]-1)]       # mx3 ndarray of x, y, z coordinates of first atom in all m bonds
        pos2 = positions[map(int, bonds[:,1]-1)]       # mx3 ndarray of x, y, z coordinates of second atom in all m bonds
        bonds[:,6] = bonds[:,5] * (pos1[:,0] - pos2[:,0]) / bonds[:,2]  # Update x-force
        bonds[:,7] = bonds[:,5] * (pos1[:,1] - pos2[:,1]) / bonds[:,2]  # Update y-force
        bonds[:,8] = bonds[:,5] * (pos1[:,2] - pos2[:,2]) / bonds[:,2]  # Update z-force
        # Similarly, translate nonbond forces into x, y, and z components
        pos1 = positions[map(int, nonbonds[:,0]-1)]
        pos2 = positions[map(int, nonbonds[:,1]-1)]
        nonbonds[:,6] = nonbonds[:,5] * (pos1[:,0] - pos2[:,0]) / nonbonds[:,2]
        nonbonds[:,7] = nonbonds[:,5] * (pos1[:,1] - pos2[:,1]) / nonbonds[:,2]
        nonbonds[:,8] = nonbonds[:,5] * (pos1[:,2] - pos2[:,2]) / nonbonds[:,2]
        
        # Calculate the total force on each atom in the molecule
        for atomID in range(1,numAtoms+1):
            forces[atomID-1] = totalForceOnAtom(atomID, bonds, nonbonds)
        
        accelerations = forces/m
        velocities = velocities + 0.5*accelerations*dt
        totalKineticEnergy = (0.5*m*velocities**2).sum()
        totalEnergyOfBonds = sum(bonds[:,4])
        totalEnergyOfNonbonds = sum(nonbonds[:,4])
        totalEnergy = totalKineticEnergy + totalEnergyOfBonds + totalEnergyOfNonbonds
        ergFile.write("%d\t%.1f\t%.1f\t%.1f\t%.1f" % (i, totalKineticEnergy, totalEnergyOfBonds, totalEnergyOfNonbonds, totalEnergy))
        
 #       if i%10==0:
 #           rvcFile.write("%d\t%.1f\t%.1f\t%.1f\t%.1f" % (i, totalKineticEnergy, totalEnergyOfBonds, totalEnergyOfNonbonds, totalEnergy))
    # Every 1 timesteps, calculate E_tot = KE + PEbonds + PEnonbonds as sanity check
    # Every 10 timesteps, write appropriate output to rvc, crd, files needed to generate euc 
    # Output: .erg and .rvc files

# Calculates x, y, and z components of the total force on a given atom 
# by summing over the individual forces from each bond/nonbond that the atom is involved with
# Inputs:
#        - atomID: ID number of atom (int)
#        - bonds: nx9 Numpy ndarray where each row represents a bond and the columns are:
#                    0) ID of atom 1 (the lower ID # is first)
#                    1) ID of atom 2
#                    2) bond length
#                    3) reference/initial bond length
#                    4) bond potential energy
#                    5) magnitude of force between bonded atoms
#                    6) x-component of force exerted on atom 1 by atom 2
#                    7) y-component of force exerted on atom 1 by atom 2
#                    8) z-component of force exerted on atom 1 by atom 2
#         - nonbonds: identical to "bonds" above, except values are for non-bonded atom pairs
# Output: 1x3 ndarray containing the x, y, and z components of the total force on the given atom
def totalForceOnAtom(atomID, bonds, nonbonds):
    # Find bonds that contain the atom. Take the negative sign of the force for bonds
    # where the atom is listed as second, to account for directionality of force
    forceForBondsWhereAtomIs1st = bonds[bonds[:,0]==float(atomID)][:,6:9]
    forceForBondsWhereAtomIs2nd = -bonds[bonds[:,1]==float(atomID)][:,6:9]
    # Vector containing three elements: the total x, y, and z forces
    # Set forces to zero if atom has no bonds
    allBondForces = np.array(forceForBondsWhereAtomIs1st.tolist() + forceForBondsWhereAtomIs2nd.tolist())
    if len(allBondForces) > 0: totalForcesDueToBonds = allBondForces.sum(axis=0)
    else: totalForcesDueToBonds = np.array([0,0,0])
            
    # Find bonds that contain the atom. Take the negative sign of the force for bonds
    # where the atom is listed as second, to account for directionality
    forceForNonbondsWhereAtomIs1st = nonbonds[nonbonds[:,0]==float(atomID)][:,6:9]
    forceForNonbondsWhereAtomIs2nd = -nonbonds[nonbonds[:,1]==float(atomID)][:,6:9]
    # Vector containing three elements: the total x, y, and z forces
    # Set forces to zero if atom has no bonds
    allNonbondForces = np.array(forceForNonbondsWhereAtomIs1st.tolist() + forceForNonbondsWhereAtomIs2nd.tolist())
    if len(allNonbondForces) > 0: totalForcesDueToNonbonds = allNonbondForces.sum(axis=0)
    else: totalForcesDueToNonbonds = np.array([0,0,0])
    
    return (totalForcesDueToBonds + totalForcesDueToNonbonds)

# initMolecule: stores initial information about the atoms in a molecule in the appropriate data structures
# Input: dictionary containing the following parameters:
#        - iF: Name of input datafile, including directory and extension (string)
#        - nbCutoff: Distance within which non-bonded atoms have interactions (float)
# Output: dictionary containing the following data structures:
#         - positions: nx3 Numpy ndarray containing the x, y, and z coordinates of n atoms
#         - velocities: nx3 Numpy ndarray containing the x, y, and z velocities of n atoms
#         - bonds: nx9 Numpy ndarray where each row represents a bond and the columns are:
#                    0) ID of atom 1 (the lower ID # is first)
#                    1) ID of atom 2
#                    2) bond length
#                    3) reference/initial bond length
#                    4) bond potential energy
#                    5) magnitude of force between bonded atoms
#                    6) x-component of force exerted on atom 1 by atom 2
#                    7) y-component of force exerted on atom 1 by atom 2
#                    8) z-component of force exerted on atom 1 by atom 2
#         - nonbonds: identical to "bonds" above, except values are for non-bonded atom pairs
#         - header: first line of the data input file
def initMolecule(params):
    data = open(params['iF']).readlines()
    header = data[0]
    numAtoms = len(data)-1
    # x, y, and z coordinates and velocities
    positions = np.empty([numAtoms, 3]); velocities = np.empty([numAtoms, 3])
    # Information about bonded and nonbonded atom pairs
    bonds = []; nonbonds = []
    # Array where each row represents an atom (in order of ascending ID)
    # and contains the ID #s of all atoms it is bonded to. Zero-indexed
    connectivities = []
    
    # For each atom (starting from atom ID #1), store its position, velocity, and atoms it is connected to
    for atomID in range(1,numAtoms+1):
        atomInfo = data[atomID].split()
        positions[atomID-1] = map(float, atomInfo[1:4])
        velocities[atomID-1] = map(float, atomInfo[4:7])
        connectivities.append(map(int, atomInfo[7:]))
    
    # For each atom, iterate over all the other atoms (to avoid double-counting, only look at atoms following the
    # atom of interest) to identify bonded/nonbonded pairs and store them in the appropriate data structures
    for currAtomID in range(1,numAtoms+1):
        for otherAtomID in range(currAtomID+1,numAtoms+1):
            distance = dist.euclidean(positions[currAtomID-1], positions[otherAtomID-1])
            if otherAtomID in connectivities[currAtomID-1]:
                bonds.append([currAtomID, otherAtomID, distance, distance, 0, 0, 0, 0, 0])
            else:
                if distance < params['nbCutoff']:
                    nonbonds.append([currAtomID, otherAtomID, distance, distance, 0, 0, 0, 0, 0])
    
    # Convert to Numpy ndarrays
    bonds = np.array(bonds); nonbonds = np.array(nonbonds)
    return {'positions':positions, 'velocities':velocities, 'bonds':bonds, 'nonbonds':nonbonds, 'header':header}
        
generateSimulation()