import sys
import numpy as np
import argparse
#import time
#start_time = time.time()

# generateSimulation is the wrapper function that runs through a molecular dynamics
# simulation and writes the resulting .rvc and .erg files
# Input: parameters from console (see initParams function)
# Output: prints resulting .rvc and .erg files
def generateSimulation():
    params = initParams()
    initialInfo = initMolecule(params)
    stepThruSimulation(params, initialInfo)

# initParams gets and processes the desired parameters from the console, and sets defaults as needed
# Input: none
# Output: a dictionary containing:
#        - iF: Input filename, including directory and extension (string).
#                Input file should be an .rvc file containing
#                 a header row as the first row. The rest of the rows represent atoms.
#                 The columns represent the atom's ID (1-indexed, in sequential order),
#                 x, y, z coordinates, x, y, z velocities, and ID's of all connected atoms
#        - kB: Spring constant for bonds (float)
#        - kN: Spring constant for non-bonds (float)
#        - nbCutoff: Distance within which non-bonded atoms have interactions (float)
#        - m: Atom mass for all atoms (float)
#        - dt: Length of timestep (float)
#        - n: Number of timesteps to iterate (int)
#        - out: Prefix of the filename to output (string)
def initParams():
    parser = argparse.ArgumentParser(description='Process inputs to program')
    parser.add_argument('--iF', metavar='iF', type=str)
    parser.add_argument('--kB', metavar='kB', type=float, default='40000.0')
    parser.add_argument('--kN', metavar='kN', type=float, default='400.0')
    parser.add_argument('--nbCutoff', metavar='nbCutoff', type=float, default='0.50')
    parser.add_argument('--m', metavar='m', type=float, default='12.0')
    parser.add_argument('--dt', metavar='dt', type=float, default='0.001')
    parser.add_argument('--n', metavar='n', type=int, default='1000')
    parser.add_argument('--out', metavar='out', type=str)
    
    args = parser.parse_args()
    # Assign arguments to local variables
    iF = args.iF; kB = args.kB; kN = args.kN; nbCutoff = args.nbCutoff; m = args.m; dt = args.dt; n = args.n                  # Number of timesteps to iterate
    # If --out is unspecified, set default to be the prefix of the input filename
    if args.out is None: out = iF.split('.rvc')[0]
    else: out = args.out
    
    # Return a dictionary with all parameters
    return {'iF':iF,'kB':kB, 'kN':kN, 'nbCutoff':nbCutoff, 'm':m, 'dt':dt, 'n':n, 'out':out}

# stepThruSimulation steps through specified number of iterations and updates the positions, 
# velocities, and energies of the atoms to eventually generate the .rvc and .erg files
# Inputs:
#    2) a dictionary containing:
#        - iF: Input filename, including directory and extension (string)
#        - kB: Spring constant for bonds (float)
#        - kN: Spring constant for non-bonds (float)
#        - nbCutoff: Distance within which non-bonded atoms have interactions (float)
#        - m: Atom mass for all atoms (float)
#        - dt: Length of timestep (float)
#        - n: Number of timesteps to iterate (int)
#        - out: Prefix of the filename to output (string)
#    2) dictionary containing the following data structures:
#         - positions: nx3 Numpy ndarray of floats containing the x, y, and z coordinates of n atoms
#         - velocities: nx3 Numpy ndarray of floats containing the x, y, and z velocities of n atoms
#         - bonds: nx9 Numpy ndarray of floats where each row represents a bond and the columns are:
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
#         - connectivities: nx1 tab-delimited array of strings where each row contains the ID's of all atoms bonded to that atom
#         - temp: temperature at which the initial rvc's were measured (string: "T=...")
#         - name: name of protein/molecule that is being simulated (string: "proteinName:")
def stepThruSimulation(params, molecule):
    # Unpack dictionary inputs into local variables (at time = t)
    positions=molecule['positions']; velocities=molecule['velocities']; bonds=molecule['bonds'];nonbonds=molecule['nonbonds'];
    connectivities=molecule['connectivities']; temp=molecule['temp']; proteinName=molecule['name']
    dt=params['dt']; n=params['n']; kB=params['kB']; kN=params['kN']; m=params['m'];out=params['out'];nbCutoff=params['nbCutoff']
    
    numAtoms = len(velocities)
    initialEnergy = 0
    
#    print 0
#    print "17-96        298-385"
#    print "%.4f\t%.4f" % (np.linalg.norm(positions[17-1] - positions[96-1]), np.linalg.norm(positions[298-1] - positions[385-1]))
    
    # nx3 Numpy ndarrays containing the x, y, z components of force and acceleration for all n atoms
    # Initialize to zero. (time = t)
    forces = np.zeros([numAtoms, 3]); accelerations = np.zeros([numAtoms, 3])
        
    # Create output files
    ergFile = open(out + "_out.erg", 'w+'); ergFile.write("# step\tE_k\tE_b\tE_nB\tE_tot\n")
    rvcFile = open(out + "_out.rvc", 'w+')
    rvcFile.write("# %s\tkB=%.1f\tkN=%.1f\tnbCutoff=%.2f\tdt=%.4f\tmass=%.1f\t%s\n" % (proteinName,kB,kN,nbCutoff,dt,m,temp))
    inputData = open(params['iF']).readlines()    
    rvcFile.writelines(inputData[1:])
     
    for i in range(1,n+1):      # Iterate n times
        # Update velocities at time = t+0.5dt*i for each atom
        # v_(t+0.5dt) = v_t + 1/2 * a_t * dt
        velocities = velocities + accelerations*0.5*dt  
        # Update positions at time = t+dt*i for each atom
        # r_(t+dt) = r_t + v_(t+0.5dt) * dt
        positions = positions + velocities*dt
        # Updates current distance between each atom at time = t+dt*i
        bonds[:,2] = [np.linalg.norm(a-b) for (a,b) in zip(positions[bonds[:,0].astype(int)-1], positions[bonds[:,1].astype(int)-1])]
        nonbonds[:,2] = [np.linalg.norm(a-b) for (a,b) in zip(positions[nonbonds[:,0].astype(int)-1], positions[nonbonds[:,1].astype(int)-1])]
        # Updates potential energies for bonds and nonbonds (time = t+dt*i)
        # PE = 1/2*spring constant*(current distance - initial distance)^2
        bonds[:,4] = 0.5*kB*(bonds[:,2]-bonds[:,3])**2; nonbonds[:,4] = 0.5*kN*(nonbonds[:,2]-nonbonds[:,3])**2
        # Calculate the magnitude of the force for each bond and nonbond (time = t+dt*i)
        bonds[:,5] = kB*(bonds[:,2]-bonds[:,3]); nonbonds[:,5] = kN*(nonbonds[:,2]-nonbonds[:,3])
        # Translate bond forces F into F_x, F_y, and F_z components
        # F_x on atom A by B = magnitude of F * (distance between A and B along x axis)/(overall distance between A and B)
        coord1 = positions[bonds[:,0].astype(int)-1]       # mx3 ndarray of x, y, z coordinates of first atom in all m bonds
        coord2 = positions[bonds[:,1].astype(int)-1]       # mx3 ndarray of x, y, z coordinates of second atom in all m bonds
        bonds[:,6] = bonds[:,5] * (coord2[:,0] - coord1[:,0]) / bonds[:,2]  # Update x-force
        bonds[:,7] = bonds[:,5] * (coord2[:,1] - coord1[:,1]) / bonds[:,2]  # Update y-force
        bonds[:,8] = bonds[:,5] * (coord2[:,2] - coord1[:,2]) / bonds[:,2]  # Update z-force
        # Do the same for nonbond forces
        coord1 = positions[nonbonds[:,0].astype(int)-1]; coord2 = positions[nonbonds[:,1].astype(int)-1]
        nonbonds[:,6] = nonbonds[:,5] * (coord2[:,0] - coord1[:,0]) / nonbonds[:,2]
        nonbonds[:,7] = nonbonds[:,5] * (coord2[:,1] - coord1[:,1]) / nonbonds[:,2]
        nonbonds[:,8] = nonbonds[:,5] * (coord2[:,2] - coord1[:,2]) / nonbonds[:,2]
        
        # Calculate the total force on each atom in the molecule (time = t+dt*i)
        for atomID in range(1,numAtoms+1): forces[atomID-1] = totalForceOnAtom(atomID, bonds, nonbonds)
        # Update accelerations and velocities on each atom (time = t+dt*i)
        accelerations = forces/m; velocities = velocities + 0.5*accelerations*dt
#        if i%10==0:
#            print i
#            print "17-96        298-385"
#            print "%.4f\t%.4f" % (np.linalg.norm(positions[17-1] - positions[96-1]), np.linalg.norm(positions[298-1] - positions[385-1]))
      
        # Every 1 timestep, calculate energies and write output to .erg file
        totalKineticEnergy = (0.5*m*velocities**2).sum()
        totalEnergyOfBonds = sum(bonds[:,4])
        totalEnergyOfNonbonds = sum(nonbonds[:,4])
        totalEnergy = totalKineticEnergy + totalEnergyOfBonds + totalEnergyOfNonbonds
        # Terminate code if energies overflow
        if i==1: initialEnergy=totalEnergy
        if totalEnergy/initialEnergy > 10**1 or totalEnergy/initialEnergy < 10**-1: sys.exit()
        # Otherwise, write output to .erg file
        ergFile.write("%d\t%.1f\t%.1f\t%.1f\t%.1f\n" % (i, totalKineticEnergy, totalEnergyOfBonds, totalEnergyOfNonbonds, totalEnergy))
              
        # Every 10 timesteps, write appropriate output to .rvc files 
        if i%10==0:
            rawOutput = np.c_[positions, velocities]
            output = np.array(["%.4f" % w for w in rawOutput.reshape(rawOutput.size)])
            output = output.reshape(rawOutput.shape)
            output = np.c_[range(1,numAtoms+1), output, connectivities]
            np.savetxt(rvcFile, output, fmt='%s', delimiter='\t', newline='\n', header='#At time step %d,energy = %.3fkJ' % (i, totalEnergy), comments='')
            
    ergFile.close()
    rvcFile.close()

   
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
    # Find all bond & nonbond interactions the atom is involved in. To account for directionality of force,
    # take the negative of the force for interactions where the atom is listed as second.
    forceForBondsWhereAtomIs1st = bonds[bonds[:,0].astype(int)==atomID][:,6:9]
    forceForBondsWhereAtomIs2nd = -bonds[bonds[:,1].astype(int)==atomID][:,6:9]
    forceForNonbondsWhereAtomIs1st = nonbonds[nonbonds[:,0].astype(int)==atomID][:,6:9]
    forceForNonbondsWhereAtomIs2nd = -nonbonds[nonbonds[:,1].astype(int)==atomID][:,6:9]
    # Sum up x, y, and z components of forces for all bond and nonbond interactions
    totalForces = np.array(forceForBondsWhereAtomIs1st.tolist() + forceForBondsWhereAtomIs2nd.tolist() + forceForNonbondsWhereAtomIs1st.tolist() + forceForNonbondsWhereAtomIs2nd.tolist())
    if len(totalForces) > 0: totalForces = totalForces.sum(axis=0)
    else: totalForces = np.array([0,0,0])
    # totalForces is a 3-element matrix containing the total x, y, and z forces for the given atom
    # If atom is not involved in any interactions with other atoms, set totalForces=0
    return totalForces

# initMolecule: stores initial information about the atoms in a molecule in the appropriate data structures
# Input: dictionary containing the following parameters:
#        - iF: Name of input datafile, including directory and extension (string)
#        - nbCutoff: Distance within which non-bonded atoms have interactions (float)
# Output: dictionary containing the following data structures:
#         - positions: nx3 Numpy ndarray of floats containing the x, y, and z coordinates of n atoms
#         - velocities: nx3 Numpy ndarray of floats containing the x, y, and z velocities of n atoms
#         - bonds: nx9 Numpy ndarray of floats where each row represents a bond and the columns are:
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
#         - connectivities: nx1 tab-delimited array of strings where each row contains the ID's of all atoms bonded to that atom
#         - temp: temperature at which the initial rvc's were measured (string: "T=...")
#         - name: name of protein/molecule that is being simulated (string: "proteinName:")
def initMolecule(params):
    data = open(params['iF']).readlines()
    # Name of protein and temperature of simulation
    temp = (data[0].split())[-1];    proteinName = (data[0].split())[1]
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
            distance = np.linalg.norm(positions[currAtomID-1] - positions[otherAtomID-1])
            #print "%d--%d: %f" %(currAtomID, otherAtomID, distance)
            if otherAtomID in connectivities[currAtomID-1]:
                # Store atom 1's ID, atom 2's ID, current distance, initial distance, PE, force, f_x, f_y, f_z
                bonds.append([currAtomID, otherAtomID, distance, distance, 0, 0, 0, 0, 0])
            else:
                if distance < params['nbCutoff']:
                    nonbonds.append([currAtomID, otherAtomID, distance, distance, 0, 0, 0, 0, 0])
    
    # Convert to Numpy ndarrays
    bonds = np.array(bonds); nonbonds = np.array(nonbonds);
    connectivities = ['\t'.join(str(connectedAtom) for connectedAtom in atom) for atom in connectivities]
    return {'positions':positions, 'velocities':velocities, 'bonds':bonds, 'nonbonds':nonbonds, 'connectivities':connectivities, 'temp':temp, 'name':proteinName}
        
generateSimulation()
#print("--- %s seconds ---" % (time.time() - start_time))