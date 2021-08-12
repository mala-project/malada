"""Provider for crystal structures."""
from .provider import Provider
import os
from shutil import copyfile


class CrystalStructureProvider(Provider):
    """
    Provides crystal structure in the form of a .cif file.

    Currently limited to copying user input. In the future this will enable
    automatic download.

    Parameters
    ----------
    parameters : malada.utils.parametes.Parameters
        Parameters used to create this object.

    external_cif_file : string
        Path to cif file provided by user. In the current state of the
        code, the pipeline will fail if this is None.
    """

    def __init__(self, parameters, external_cif_file=None):
        super(CrystalStructureProvider, self).__init__(parameters)
        self.external_cif_file = external_cif_file
        self.cif_file = None

    def provide(self, provider_path):
        """
        Provide a crystal structure in the form of a cif file.

        Parameters
        ----------
        provider_path : string
            Path in which to operate in.
        """
        file_name = self.parameters.element + \
                    "_" + self.parameters.crystal_structure + ".cif"
        self.cif_file = os.path.join(provider_path,file_name)
        if self.external_cif_file is None:
            raise Exception("Currently there is no way to provide a cif file"
                            "on the fly.")
        else:

            copyfile(self.external_cif_file, self.cif_file)
            print("Getting <<crystal_structure>>.cif file from disc.")
