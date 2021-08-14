# List of providers

The following is a lit of providers. 
"Input" in this context always refers to the variables that have to be given to the `provide()` function. "Output" values are attributes of the classes,
that can be accessed after running. "Optional inputs" are inputs given upon creation, that allow to build complex workflows, by giving a provider half or fully processed data. Generally speaking, the following interactions with a provider are possible:

```python 
some_provider = SomeProvider(optional_param1=some_value)
some_provider.provide()
print(some_provider.some_output_quantity)
```

[Global parameters](global_parameters.md) will not be counted as inputs in this context. Every provider has access to them - they define how a pipeline is built.

## CrystalStructureProvider

Provides a crystal structure.

- input: None
- output: 
  - .cif file containing the crystal structure

## SupercellProvider

Provides a supercell in VASP format.

- input: .cif file containing the crystal structure
- output: .vasp file containing the positions of the atoms
    
## DFTConvergenceProvider

Calculates optimal DFT parameters (cutoff energy, k-grid).

- input: 
  - vasp file containing the positions of the atoms
- output: 
  - .xml file containing converged DFT parameters (k-grid and cutoff energy for basis set)
    
## MDPerformanceProvider

Determines optimal MD performance for large scale calculations.

- input: 
  - .xml file containing converged DFT parameters (k-grid and cutoff energy for basis set), .vasp file containing the positions of the atoms
- output: .xml file containing optimal parameters for MD simulation
    
## MDProvider

Calculates DFT-MD trajectory.

- input: .xml file containing converged DFT parameters (k-grid and cutoff energy for basis set), .vasp file containing the positions of the atoms, .xml file containing optimal parameters for MD simulation
- output: .traj file containing trajectory of the MD run, .npy file containing temperatures

## SnapshotProvider

Parses MD trajectory for suitable MD snapshots. 

- input: .traj file containing trajectory of the MD run, .npy file containing temperatures
- output: .traj file containing potential snapshots for DFT calculation
    
## LDOSConvergenceProvider

Determines optimal parameters for LDOS calculation.

 - input: .traj file containing potential snapshots for DFT calculation, .xml file containing converged DFT parameters (k-grid and cutoff energy for basis set)
 - output: .xml file containing updated LDOS parameters for post processing
    
## DFTProvider:

 - input: xml file containing converged DFT parameters (k-grid and cutoff energy for basis set), .xml file containing updated LDOS parameters for post processing, .traj file containing potential snapshots for DFT calculation
 - output: DFT calculation output files and LDOS cube files





