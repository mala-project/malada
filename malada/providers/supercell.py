"""Provider for creation of supercell from crystal structure."""
from malada.utils import structure_to_transformation
from .provider import Provider
import os
import ase.io
import ase.build
from shutil import copyfile


class SuperCellProvider(Provider):
    """
    Builds a supercell file (vasp format).

    Parameters
    ----------
    parameters : malada.utils.parametes.Parameters
        Parameters used to create this object.
    """

    def __init__(self, parameters, external_supercell_file=None):
        super(SuperCellProvider, self).__init__(parameters)
        self.external_supercell_file = external_supercell_file
        self.supercell_file = None

    def provide(self, provider_path, cif_file):
        """
        Provide supercell file in VASP format.

        Parameters
        ----------
        provider_path : string
            Path in which to operate in.

        cif_file : string
            Path to cif file used for supercell creation.
        """
        file_name = self.parameters.element + "_" + str(self.parameters.number_of_atoms) \
                    + "_" + self.parameters.crystal_structure + ".vasp"
        self.supercell_file = os.path.join(provider_path, file_name)
        if self.external_supercell_file is None:
            primitive_cell = ase.io.read(cif_file, format="cif")
            transformation_matrix = structure_to_transformation\
                                    [self.parameters.crystal_structure]\
                                    [self.parameters.number_of_atoms]

            super_cell = ase.build.make_supercell(primitive_cell,
                                                  transformation_matrix)
            ase.io.write(self.supercell_file,
                         super_cell, format="vasp", long_format=True)
        else:
            copyfile(self.external_supercell_file, self.supercell_file)
            print("Getting <<supercell>>.vasp file from disc.")
