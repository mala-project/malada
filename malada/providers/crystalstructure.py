from .provider import Provider
import os
from shutil import copyfile

class CrystalStructureProvider(Provider):
    """Provides crystal structure in the form of a .cif file."""
    def __init__(self, parameters, external_cif_file=None):
        super(CrystalStructureProvider, self).__init__(parameters)
        self.external_cif_file = external_cif_file
        self.cif_file = None

    def provide(self, provider_path):
        file_name = self.parameters.element + \
                    "_" + self.parameters.crystal_structure + ".cif"
        self.cif_file = os.path.join(provider_path,file_name)
        if self.external_cif_file is None:
            raise Exception("Currently there is no way to provide a cif file"
                            "on the fly.")
        else:

            copyfile(self.external_cif_file, self.cif_file)
            print("Getting <<crystal_structure>>.cif file from disc.")
