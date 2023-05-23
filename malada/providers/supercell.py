"""Provider for creation of supercell from crystal structure."""
from malada.utils import structure_to_transformation
from .provider import Provider
import os
import ase.io
import ase.build
from ase.units import m, kg, Bohr
from shutil import copyfile
import numpy as np


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
            super_cell = self.get_compressed_cell(super_cell,
                                                  stretch_factor = self.parameters.stretch_factor,
                                                  density = self.parameters.mass_density,
                                                  radius = self.parameters.WS_radius)
            ase.io.write(self.supercell_file,
                         super_cell, format="vasp", long_format=True)
        else:
            copyfile(self.external_supercell_file, self.supercell_file)
            print("Getting <<supercell>>.vasp file from disc.")

    @staticmethod
    def get_compressed_cell(
        supercell,
        stretch_factor=None,
        density=None,
        radius=None,
        units_density="g/(cm^3)",
    ):

        if sum([stretch_factor is None, density is None, radius is None]) == 3:
            raise ValueError(
                "At least one of stretch_factor, density and radius must be speficied"
            )
        elif sum([stretch_factor is None, density is None, radius is None]) < 2:
            print("Warning: More than one of stretch factor, density, "
                  "and radius is specified.\nRadius takes first priority, "
                  "then density, finally stretch factor.")
        if density is not None:
            density_ambient = SuperCellProvider.get_mass_density(
                supercell, unit=units_density
            )
            stretch_factor = (density_ambient / density) ** (1.0 / 3.0)
        if radius is not None:
            radius_ambient = SuperCellProvider.get_wigner_seitz_radius(supercell)
            stretch_factor = radius / radius_ambient

        supercell.set_cell(supercell.get_cell() * stretch_factor, scale_atoms=True)

        return supercell

    @staticmethod
    def get_number_density(supercell, unit="Angstrom^3"):
        nr_atoms = len(supercell)
        volume_atoms = supercell.get_volume()
        nr_electrons = 0
        for i in range(0, nr_atoms):
            nr_electrons += supercell[i].number
        number_density = nr_electrons/volume_atoms
        angstrom3_in_m3 = 1 / (m * m * m)
        angstrom3_in_cm3 = angstrom3_in_m3 * 1000000
        if unit == "Angstrom^3":
            return number_density
        elif unit == "cm^3":
            return number_density/angstrom3_in_cm3
        else:
            raise Exception("Unit not implemented")

    @staticmethod
    def get_wigner_seitz_radius(supercell):
        number_density = SuperCellProvider.\
            get_number_density(supercell, unit="Angstrom^3")
        rs_angstrom = (3 / (4 * np.pi * number_density)) ** (1 / 3)
        return rs_angstrom / Bohr

    @staticmethod
    def get_mass_density(supercell, unit="g/(cm^3)"):
        nr_atoms = len(supercell)
        mass_atom = supercell[0].mass
        mass_atoms = nr_atoms * mass_atom
        volume_atoms = supercell.get_volume()
        mass_density = mass_atoms / volume_atoms
        u_angstrom3_in_kg_m3 = kg / (m * m * m)
        u_angstrom3_in_g_cm3 = u_angstrom3_in_kg_m3 * 1000
        if unit == "kg_m^3":
            return mass_density/u_angstrom3_in_kg_m3
        elif unit == "g/(cm^3)":
            return mass_density/u_angstrom3_in_g_cm3
        elif unit == "u_Angstrom^3":
            return mass_density
        else:
            raise Exception("Unit not implemented")


