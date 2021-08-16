# Current limitations 

MALADA is currently in its very early development stages. Please feel free to report any bug or oddity you encounter. The following is a list of known limitations users currently have to work around, please excuse the inconveniences. 

- Run system: The run system is currently rudimentary, and uses system calls (bash), or the automatic creation of slurm batch scripts. This of course comes with some problems and will fail on some machines. In the future, a more stable run system will be provided drawing on ASE and python-slurm interfaces.
- Automatic .cif download: This feature is currently not implemented and users will have to provide their own crystal structures.
- MD performance detection: This feature is currently not implemented, and users will have to provide their own performance parameter
- LDOS convergence: This feature is currently not implemented, and users will have to provide their own parameters
- Species: MALADA can currently only handle one chemical species at a time 
- Temperature: MALADA has not been tested with higher temperatures, there could be some automated parameter estimations that fail
- Mass densities: MALADA currently has no feature to automatically work with higher mass densities, and will use the one from the .cif file
- MD trajectory parsing: This is currently limited to randomly sampling snapshots; other metrics will follow soon
- Also see https://github.com/mala-project/malada/issues
