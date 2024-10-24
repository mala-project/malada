import malada

"""
This example shows how MALADA can be used to run a simple pipeline for 
creating atomic position / LDOS data for machine learning workflows. 
The example system is a very limited system of Be atoms.
The example sets parameter in such a way that calculations will be performed
fast - not necessarily in the most physical sense. 
"""
###########################################
# First, we have to define the parameters.
###########################################

# These describe the general quantities of the system.
params = malada.Parameters()
params.number_of_atoms = 2
params.temperature = 298
params.element = "Be"
params.base_folder = "ex01"
params.pseudopotential["valence_electrons"] = 2
params.pseudopotential["name"] = "Be.pbe-n-rrkjus_psl.1.0.0.UPF"

# This path has to be altered updated with your path.
params.pseudopotential["path"] = (
    "/home/fiedlerl/tools/pslibrary/" "pbe/PSEUDOPOTENTIALS/"
)

# These are technical parameters for DFT and MD.
params.maximum_kpoint_try = 3
params.snapshot_parsing_beginning = 5
params.number_of_snapshots = 2

# All of these parameters are chosen too small/large to accomodate for
# fast execution.
# Generally we would want more timesteps, a smaller accuracy threshold,
# a larger distance cutoff and a smaller temperature tolerance.
params.maximum_number_of_timesteps = 20
params.dft_scf_accuracy_per_atom_Ry = 1e-1
params.dft_conv_accuracy_meVperatom = 10
params.distance_metric_snapshots_cutoff = 0.0001
params.snapshot_parsing_temperature_tolerance_percent = 10

# Currently, both the crystal structure and LDOS configuration have to be
# provided by the user.
crystal_structure = malada.CrystalStructureProvider(
    params, external_cif_file="./ex01_inputs/Be_bcc.cif"
)

ldosconvergence = malada.LDOSConvergenceProvider(
    params, external_ldos_configuration="./ex01_inputs/ldosconf.xml"
)

# After defining parameters and important files, the pipeline can be
# instantiated and run.
pipeline = malada.DataPipeline(
    params,
    crystal_structure_provider=crystal_structure,
    ldos_configuration_provider=ldosconvergence,
)
pipeline.run()
