"""Collection of all parameters for the pipeline."""

import os
from .slurmparams import SlurmParameters


class Parameters:
    """
    Holds parameters needed for constructing a data acqusition pipeline.

    Attributes
    ----------
    temperature : float
        Temperature in K.

    number_of_atoms : int
        Number of atoms in supercell.

    crystal_structure : string
        Crystal structure, this is used for supercell creation because
        it determines the transformations.

    element : string
        Name of the element for which to generate data.

    base_folder : string
        Base folder to operate in.

    dft_calculator : string
        Name of the DFT calculator, currently only QuantumESPRESSO and VASP
        are supported.

    md_calculator : string
        Name of the DFT-MD calculator, currently only QuantumESPRESSO and VASP
        are supported.

    pseudopotential : dict
        Dictionary for the pseudopotential. Needs to contain "path" for
        path to pseudopotential, "name" for name, and "valence_electrons"
        for number of valence electrons.

    run_system : string
        Run system used during pipeline, currently only "bash" and
        "slurm_creator" are supported.

    dft_slurm : malada.SlurmParameters
        Slurm parameters used for DFT calculations

    md_slurm : malada.SlurmParameters
        Slurm parameters used for DFT-MD calculations. Please note that
        these should not be set here, as they will be overwritten in
        the pipeline.

    dft_conv_accuracy_meVperatom : float
        Accuracy of the DFT convergence calculations in meV/atom.

    maximum_cutoff_try : int
        Maxmimum number of cutoff convergence iterations. If the optimal
        cutoff is not found with initial values, calculation
        of it it is retried maximum_cutoff_try number of times.

    maximum_kpoint_try : int
        Maxmimum number of k-grid convergence iterations. If the optimal
        k-grid is not found with initial values, calculation
        of it it is retried maximum_kpoint_try number of times.

    dft_scf_accuracy_per_atom_Ry : float
        Accuracy of all DFT calculation in terms of total energy/per atom.

    md_at_gamma_point : bool
        If True, MD simulations will be run at the gamma point.

    maximum_number_of_timesteps : int
        Maximum number of timesteps for MD calculation.

    time_step_fs : float
        Timestep for MD calculation in fs.

    md_thermostat_controller : float
        Thermostat controller for MD calculation. This is "NRAISE" for QE
        and "SMASS" for VASP.

    snapshot_parsing_beginning : int
        The snapshot after which snapshot parsing is started. If < 0 automatic
        detection of this snapshot will be performed (currently not supported)

    snapshot_parsing_temperature_tolerance_percent : float
        Maximum deviation of temperature between snapshot and desired
        temperature for snapshot to be considered for DFT calculation
        (in percent)

    snapshot_parsing_criterion : string
        Criterion with which to parse the snapshots, currently only
        "random" is supported.

    number_of_snapshots : int
        Number of snapshots the pipeline shoudl generate.

    distance_metric_snapshots : string
        Distance metric to determine how alike two atomic snapshots are.
        Currently only "realspace", i.e. the realspace distance between atoms
        is calculated.

    distance_metric_snapshots_cutoff : float
        Minimum distance in terms of distance_metric_snapshots after which
        after which two snapshots are considered distinct by the algorithm.

    number_of_bands_factor : float
        Determines how many more many more bands then number of electrons will
        be used. E.g. if =0.05, the overall number of bands will be the number
        of electrons times 1.05. Has to be scaled with temperature. Default
        is 0.05, which should be ok up to ~2500K.

    dft_assume_two_dimensional : bool
        If True, two dimensional DFT calculations will be performed.
        To that end, the relevant QE parameters will be selected.
    """

    def __init__(self):
        # General information about the system.
        self.temperature = 0
        self.stretch_factor = None
        self.mass_density = None
        self.WS_radius = None
        self.number_of_atoms = 0
        self.crystal_structure = "bcc"
        self.element = None

        # Information about running the pipeline.
        self.base_folder = "./"
        self.dft_calculator = "qe"
        self.md_calculator = "qe"
        # TODO: Get number of electrons directly from file.
        self.pseudopotential = {
            "path": None,
            "valence_electrons": 0,
            "name": None,
        }
        self.run_system = "bash"
        self.mp_api_file = os.path.expanduser("~") + "/malada/.mp_api/.api_key"
        self.dft_slurm = SlurmParameters()
        self.md_slurm = SlurmParameters()

        # Information about DFT and MD calculations.
        self.dft_conv_accuracy_meVperatom = 1
        self.maximum_cutoff_try = 2
        self.maximum_kpoint_try = 2
        self.dft_scf_accuracy_per_atom_Ry = 1e-6
        self.md_at_gamma_point = True
        self.maximum_number_of_timesteps = 10000
        self.time_step_fs = 1
        self.md_thermostat_controller = 0.001
        self.number_of_bands_factor = 0.05
        self.dft_calculate_stress = True
        self.dft_calculate_force = True
        self.dft_use_inversion_symmetry = False
        self.dft_mixing_beta = 0.1
        self.dft_assume_two_dimensional = False
        self.twodimensional_cutting_tolerance = 1e-3

        # Information about MD parsing.
        self.snapshot_parsing_beginning = -1
        self.snapshot_parsing_temperature_tolerance_percent = 1
        self.number_of_snapshots = 10
        # Technical parameters; for certain types of crystals, without
        # such a tolerance, the RDF calculation might fail. Also we need
        # to specify the number of bins. Usually, we do not have to temper
        # with these values, the defaults should give pretty steady results.
        self.distance_metric_snapshots_rdf_tolerance = 0.0001
        self.distance_metric_snapshots_rdf_bins = 500
        # The distance metric is denoised prior to analysis using a certain
        # width. This should be adjusted if there is reason to believe
        # the trajectory will be noise for some reason.
        self.distance_metrics_denoising_width = 100
        # Number of time steps that have to consecutively below the average
        # of the distance metric curve, before we consider the trajectory
        # to be equilibrated.
        # Usually does not have to be changed.
        self.distance_metrics_below_average_counter = 50
        # The analysis of the trajectory builds on the assumption that at some
        # point of the trajectory, the system is equilibrated.
        # For this, we need to provide the fraction of the trajectory (counted
        # from the end). Usually, 10% is a fine assumption. This value usually
        # does not need to be changed.
        self.distance_metrics_estimated_equilibrium = 0.1
        self.distance_metric_snapshots_cutoff = -0.1
