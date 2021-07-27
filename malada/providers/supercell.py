from .provider import Provider
import ase.io

class SuperCellProvider(Provider):
    """Provides crystal structure in the form of a .cif file."""
    def __init__(self, parameters, supercell_file=None):
        super(SuperCellProvider, self).__init__(parameters)
        self.supercell_file = supercell_file

    def provide(self, cif_file):
        if self.supercell_file is None:
            primitive_cell = ase.io.read(cif_file, format="cif")
            if self.parameters.number_of_atoms == 256:
                transformation_matrix = [[8, 0, 0],
                                         [0, 4, 0],
                                         [0, 0, 4]]

            if number_of_atoms == 128:
                transformation_matrix = [[4, 0, 0],
                                         [0, 4, 0],
                                         [0, 0, 4]]
            print(element, structure, number_of_atoms)
            print(element + "_" + structure + "_" + str(
                number_of_atoms) + "_atoms.vasp")
            super_cell = ase.build.make_supercell(primitive_cell,
                                                  transformation_matrix)
            ase.io.write(element + "_" + structure + "_" + str(
                number_of_atoms) + "_atoms.vasp",
                         super_cell, format="vasp", long_format=True)
        else:
            print("Getting <<supercell>>.vasp file from disc.")
