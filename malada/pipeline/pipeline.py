from malada import CrystalStructureProvider

class DataPipeline:
    """Uses providers to run an entire data pipeline, ending with the LDOS."""

    def __init__(self, parameters, crystal_structure_provider=None):
        self.parameters = parameters

        # Create the providers object that were not otherwise specified.
        self.crystal_structure_provider = \
            CrystalStructureProvider(self.parameters)


    def run(self):
        # Step one: Get the crystal structure.
        print("Getting the crystal structure...")

        self.crystal_structure_provider.provide()
        print(self.crystal_structure_provider.cif_file)

        print("Getting the crystal structure: Done.")
