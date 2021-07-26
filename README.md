# MALA Data Acquistion (MALADA)

## Scope

MALADA provides helpful functions to automate the data acquisition process for the Materials Learning Algorithms (MALA) package. It helps the user create and run the typical MALA data acquisition pipelin, consisting of creation of suitable supercells from crystal structures, determination of suitable calculation parameters, molecular dynamics (MD) simulations, density functional theory simulations (DFT) and post-processing of the latter to acquire the local density of states (LDOS). At the end of a run, MALA can be used to preprocess the data and build surrogate models. 

## Usage

### Providers

MALADA is a python package and works by writing run scripts. 
It is based on "providers" that provide certain files for the next step. These can be directly linked in the python script or via files. The following providers exist, reflecting the steps of a data acquisition pipeline: 

1. CrystalStructureProvider
    - input: Element name, element ID in MaterialsProject
    - output: .cif file containing the crystal structure

2. SupercellProvider
    - input: .cif file containing the crystal structure, number of atoms
    - output: .vasp file containing the positions of the atoms
    
3. DFTConvergenceProvider
    - input: temperature at which to perform the simulation, .vasp file containing the positions of the atoms
    - output: .xml file containing converged DFT parameters (k-grid and cutoff energy for basis set)
    
4. MDPerformanceProvider
    - input: .pkl file containing converged DFT parameters (k-grid and cutoff energy for basis set), temperature at which to perform the simulation, .vasp file containing the positions of the atoms
    - output: .xml file containing optimal parameters for MD simulation
    
5. MDProvider
    - input: .pkl file containing converged DFT parameters (k-grid and cutoff energy for basis set), temperature at which to perform the simulation, .vasp file containing the positions of the atoms, .xml file containing optimal parameters for MD simulation
    - output: .traj file containing trajectory of the MD run, .npy files containing temperature

6. SnapshotProvider
    - input: .traj file containing trajectory of the MD run
    - output: .npy file containing potential snapshots for DFT calculation
    
7. LDOSConvergenceProvider
    - input: .npy file containing potential snapshots for DFT calculation, .xml file containing converged DFT parameters (k-grid and cutoff energy for basis set)
    - output: .xml file containing updated DFT parameters for post processing
    
8. DFTProvider:
    - input: .xml file cotnaining DFT parameters (for both energy and LDOS), .npy file containing potential snapshots for DFT calculation
    - output: DFT calculation output files and LDOS cube files


### Global parameters

Global parameters such as number of atoms, temperature and runner information (which run system, e.g. slurm, how many cores, etc.)
