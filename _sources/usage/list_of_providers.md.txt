# List of providers

The following is a lit of providers. 
"Input" in this context always refers to the variables that have to be given to the `provide()` function. Each provider requires a path to operate in as input, which will be omitted further down. "Output" values are attributes of the classes,
that can be accessed after running. "Optional inputs" are inputs given upon creation, that allow to build complex workflows, by giving a provider half or fully processed data. Generally speaking, the following interactions with a provider are possible:

```python 
some_provider = SomeProvider(optional_param1=some_value)
some_provider.provide(provier_path)
print(some_provider.some_output_quantity)
```

[Global parameters](global_parameters.md) will not be counted as inputs in this context. Every provider has access to them - they define how a pipeline is built.

## CrystalStructureProvider

Provides a crystal structure.

- input: None
- output: 
  - cif_file: .cif file containing the crystal structure
- optional inputs:
  - external_cif_file: path to a .cif file previously downloaded that is to be used by MALADA. If given, no calculation/download will be performed.

## SupercellProvider

Provides a supercell in VASP format.

- input: 
  - cif_file: .cif file containing the crystal structure
- output: 
  - supercell_file: .vasp file containing the positions of the atoms
- optional inputs:
  - external_supercell_file: path to a.vasp file containing the positions of the atoms to be used by MALADA. If given, no calculation will be performed. 

    
## DFTConvergenceProvider

Calculates optimal DFT parameters (cutoff energy, k-grid).

- input: 
  - supercell_file: .vasp file containing the positions of the atoms
- output: 
  - convergence_results_file: .xml file containing converged DFT parameters (k-grid and cutoff energy for basis set)
- optional inputs:
  - external_convergence_results: path to a .xml file containing converged DFT parameters, if given, no DFT calculations will be performed
  - external_convergence_folder: path to a folder containing already calculated DFT results from which convergence can be determined; if given, no (or less) DFT calculations will be performed
  - predefined_kgrid: If an optimal k-grid is known, it can be provided and then no k-grid convergence calculations will be performed
  - predefined_cutoff: If an optimal cutoff energy is known, it can be provided and then no cutoff energy convergence calculations will be performed
    
## MDPerformanceProvider

Determines optimal MD performance for large scale calculations.

- input: 
  - dft_convergence_file: .xml file containing converged DFT parameters (k-grid and cutoff energy for basis set)
  - supercell_file: .vasp file containing the positions of the atoms
- output: 
  - md_performance_xml: .xml file containing optimal parameters for MD simulation
- optional inputs:
  - external_performance_file: path to a .xml file containing optimal parameters for MD simulation, if given, no calculations will be performed

## MDProvider

Calculates DFT-MD trajectory.

- input: 
  - dft_convergence_file: .xml file containing converged DFT parameters (k-grid and cutoff energy for basis set)
  - supercell_file: .vasp file containing the positions of the atoms
  - md_performance_xml: .xml file containing optimal parameters for MD simulation
- output: 
  - trajectory_file: .traj file containing trajectory of the MD run 
  - temperature_file: .npy file containing temperatures
- optional inputs:
  - external_trajectory: path to a .traj file containing trajectory of the MD run, if given in conjucntion with external_temperatures, no calculations will be performed
  - external_temperatures: path to a .npy file containing temperatures, if given in conjucntion with external_trajectory, no calculations will be performed
  - external_run_folder: path to a directory in which a MD simulation was performed, if given, only analysis will be performed on this MD run
  
## SnapshotProvider

Parses MD trajectory for suitable MD snapshots. 

- input: 
  - trajectory_file: .traj file containing trajectory of the MD run 
  - temperature_file: .npy file containing temperatures
- output:
  - snapshot_file: .traj file containing potential snapshots for DFT calculation
- optional inputs:
  - external_snapshots: path to a .traj file containing potential snapshots for DFT calculation, if given, no analysis will be performed on MD data 
  
## LDOSConvergenceProvider

Determines optimal parameters for LDOS calculation.

- input: 
  - snapshot_file: .traj file containing potential snapshots for DFT calculation
  - dft_convergence_file: .xml file containing converged DFT parameters (k-grid and cutoff energy for basis set)
- output: 
  - ldos_configuration_file: .xml file containing updated LDOS parameters for post processing
- optional inputs:
  - external_ldos_configuration: path to a .xml file containing updated LDOS parameters for post processing, if given, no calculation will be performed
      
## DFTProvider:

- input: 
  - dft_convergence_file: .xml file containing converged DFT parameters (k-grid and cutoff energy for basis set)
  - ldos_configuration_file: .xml file containing updated LDOS parameters for post processing
  - snapshot_file: .traj file containing potential snapshots for DFT calculation
- output: 
  - calculation_folders: directory containing DFT calculation output files and LDOS cube files
- optional inputs:
  - external_calculation_folders: directory containing DFT calculation output files and LDOS cube files, if given, no calculations will be provided, only a consistency check
