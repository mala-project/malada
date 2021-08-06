from malada import CrystalStructureProvider, SuperCellProvider, \
                   DFTConvergenceProvider, MDPerformanceProvider, MDProvider, \
                   SnapshotsProvider
import os


class DataPipeline:
    """Uses providers to run an entire data pipeline, ending with the LDOS."""

    def __init__(self, parameters,
                 crystal_structure_provider: CrystalStructureProvider = None,
                 supercell_provider: SuperCellProvider = None,
                 dft_convergence_provider: DFTConvergenceProvider = None,
                 md_performance_provider: MDPerformanceProvider = None,
                 md_provider: MDProvider = None,
                 snapshots_provider: SnapshotsProvider = None,
                 ):
        self.parameters = parameters

        # Create the providers object that were not otherwise specified.
        if crystal_structure_provider is None:
            self.crystal_structure_provider = \
                CrystalStructureProvider(self.parameters)
        else:
            self.crystal_structure_provider = crystal_structure_provider
        if supercell_provider is None:
            self.supercell_provider = \
                SuperCellProvider(self.parameters)
        else:
            self.supercell_provider = supercell_provider
        if dft_convergence_provider is None:
            self.dft_convergence_provider = \
                DFTConvergenceProvider(self.parameters)
        else:
            self.dft_convergence_provider = dft_convergence_provider
        if md_performance_provider is None:
            self.md_performance_provider = \
                MDPerformanceProvider(self.parameters)
        else:
            self.md_performance_provider = md_performance_provider
        if md_provider is None:
            self.md_provider = \
                MDProvider(self.parameters)
        else:
            self.md_provider = md_provider
        if snapshots_provider is None:
            self.snapshots_provider = SnapshotsProvider(self.parameters)
        else:
            self.snapshots_provider = snapshots_provider

    def run(self):
        # Step one: Get the crystal structure.
        print("Getting the crystal structure...")
        path00 = os.path.join(self.parameters.base_folder,
                              "00_crystal_structure")
        if not os.path.exists(path00):
            os.makedirs(path00)
        self.crystal_structure_provider.provide(path00)

        print("Getting the crystal structure: Done.")

        # Step two: Build the supercell.
        print("Building supercell...")
        path01 = os.path.join(self.parameters.base_folder,
                              "01_crystal_structure")
        if not os.path.exists(path01):
            os.makedirs(path01)
        self.supercell_provider.provide(path01, self.crystal_structure_provider.cif_file)
        print("Building supercell: Done.")

        # Step three: Converge DFT parameters.
        print("Converging DFT parameters...")
        path02 = os.path.join(self.parameters.base_folder, "02_dft_convergence")
        if not os.path.exists(path02):
            os.makedirs(path02)
        self.dft_convergence_provider.provide(path02,
                                              self.supercell_provider.supercell_file)
        print("Converging DFT parameters: Done.")

        # Step four: Get optimal MD run parameters.
        print("Testing MD performance...")
        path03 = os.path.join(self.parameters.base_folder, "03_md_performance")
        if not os.path.exists(path03):
            os.makedirs(path03)
        self.md_performance_provider.provide(path03,
                                             self.dft_convergence_provider.
                                             convergence_results_file)
        print("Testing MD performance: Done.")

        # Step five: Get/Calculate a MD trajectory.
        print("Getting MD trajectory...")
        path04 = os.path.join(self.parameters.base_folder, "04_md")
        if not os.path.exists(path04):
            os.makedirs(path04)
        self.md_provider.provide(path04, self.supercell_provider.supercell_file,
                                 self.dft_convergence_provider.
                                 convergence_results_file, self.
                                 md_performance_provider.md_performance_xml)
        print("Getting MD trajectory: Done.")

        # Step six: Parsing MD trajectory for snapshots.
        print("Parsing snapshots from MD trajectory...")
        path05 = os.path.join(self.parameters.base_folder, "05_snapshots")
        if not os.path.exists(path05):
            os.makedirs(path05)
        self.snapshots_provider.provide(path05,
                                        self.md_provider.trajectory_file,
                                        self.md_provider.temperature_file)
        print("Parsing snapshots from MD trajectory: Done.")

