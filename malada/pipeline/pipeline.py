from malada import CrystalStructureProvider, SuperCellProvider
import os

class DataPipeline:
    """Uses providers to run an entire data pipeline, ending with the LDOS."""

    def __init__(self, parameters,
                 crystal_structure_provider:CrystalStructureProvider=None,
                 supercell_provider:SuperCellProvider=None,
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

    def run(self):
        # Step one: Get the crystal structure.
        print("Getting the crystal structure...")
        path00 = os.path.join(self.parameters.base_folder, "00_crystal_structure")
        if not os.path.exists(path00):
            os.makedirs(path00)
        self.crystal_structure_provider.provide(path00)

        print("Getting the crystal structure: Done.")

        # Step two: Build the supercell.
        print("Building supercell...")
        path01 = os.path.join(self.parameters.base_folder, "01_crystal_structure")
        if not os.path.exists(path01):
            os.makedirs(path01)
        self.supercell_provider.provide(path01, self.crystal_structure_provider.cif_file)
        print("Building supercell: Done.")

