class Parameters:
    """Holds parameters needed for constructing a data acqusition pipeline."""
    def __init__(self):
        self.temperature = 0
        self.run_system = "bash"
        self.number_of_atoms = 0
        self.crystal_structure = "bcc"
        self.element = None
        self.base_folder = "./"
        self.dft_calculator = "qe"
        # TODO: Get number of electrons directly from file.
        self.pseudopotential = {"path": None, "valence_electrons": 0,
                                "name": None}
        self.dft_accuracy_meVperatom = 1
        self.maximum_cutoff_try = 2
        self.maximum_kpoint_try = 2
        self.dft_scf_accuracy_per_atom = 1e-6
