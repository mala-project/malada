# General Usage

The idea of MALADA is to automate as many steps of the data generation for MALA workflow as posible. It is not meant to be a universal tool for DFT-MD workflows, and is specifically tailored to this application. In order to create training data for MALA networks, the following steps need to be undertaken:

1. Crystal structure determination (What element in what structure do you want to model?)
2. Supercell construction (How many atoms will be in your supercell?)
3. DFT convergence (Which cutoff energy and k-grid should be used for DFT calculations?)
4. DFT-MD performance test (How do you need to parallelize a DFT-MD run so that you can get results in reasonable time?)
5. DFT-MD simulations (Calculate a trajectory with your given parameters)
6. Snapshot selection (Filter this trajectory)
7. LDOS parameter determination (How should the LDOS be discretized? On which grid should it be discretized?)
8. DFT simulations, LDOS calculation (Perform a DFT calculation and calculate the LDOS from the final wave functions)

To automate this, MALADA provides a [pipeline](pipeline.md), that can be controlled using the [global parameters](global_parameters.md).
