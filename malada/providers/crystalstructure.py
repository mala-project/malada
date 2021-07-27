from .provider import Provider

class CrystalStructureProvider(Provider):
    """Provides crystal structure in the form of a .cif file."""
    def __init__(self, parameters, cif_file=None):
        super(CrystalStructureProvider, self).__init__(parameters)
        self.cif_file = cif_file


    def provide(self):
        if self.cif_file is None:
            raise Exception("Currently there is no way to provide a cif file"
                            "on the fly.")
        else:
            print("Getting <<crystal_structure>>.cif file from disc.")
