from malada import CrystalStructureProvider


class DataPipeline:
    """Uses providers to run an entire data pipeline, ending with the LDOS."""

    def __init__(self, parameters, crystal_structure_provider=None):
        self.parameters = parameters

        # Create the providers object that were not otherwise specified.
        if crystal_structure_provider is None:
            self.crystal_structure_provider = \
                CrystalStructureProvider(self.parameters)
        else:
            self.crystal_structure_provider = crystal_structure_provider

    def run(self):
        # Step one: Get the crystal structure.
        print("Getting the crystal structure...")
        self.crystal_structure_provider.provide()
        print("Getting the crystal structure: Done.")

        # Step two: Build the supercell.
        print("Building supercell...")
        self.crystal_structure_provider.provide()
        print("Building supercell: Done.")

