# Global parameters

Each object in the MALADA workflow needs a `Parameters` object provided at time of instantiation. It holds information such as number of atoms, temperature and runner information (which run system, e.g. slurm, how many cores, etc.). Some of these are optional, and default values can be used. Others have to be explicitly set for the workflow not to fail. For a detailed description of each parameter, see the API documentation. 

## Required parameters

- temperature 
- number_of_atoms 
- element

## Optional parameters

- crystal_structure
- base_folder
- dft_calculator
- md_calculator
- pseudopotential
- run_system
- dft_slurm
- md_slurm
- dft_conv_accuracy_meVperatom
- maximum_cutoff_try
- maximum_kpoint_try
- dft_scf_accuracy_per_atom_Ry
- md_at_gamma_point
- maximum_number_of_timesteps 
- time_step_fs
- md_thermostat_controller 
- snapshot_parsing_beginning 
- snapshot_parsing_temperature_tolerance_percent
- snapshot_parsing_criterion 
- number_of_snapshots 
- distance_metric_snapshots
- distance_metric_snapshots_cutoff
