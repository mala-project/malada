"""Provider for creation of supercell from crystal structure."""

from .provider import Provider
import os
import ase.io
import ase.build
from ase.units import m, kg, Bohr
from shutil import copyfile
import numpy as np

try:
    from mp_api.client import MPRester
except:
    pass


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
        file_name = (
            self.parameters.element
            + "_"
            + str(self.parameters.number_of_atoms)
            + "_"
            + self.parameters.crystal_structure
            + ".vasp"
        )
        self.supercell_file = os.path.join(provider_path, file_name)
        if self.external_supercell_file is None:
            try:
                primitive_cell = ase.io.read(cif_file, format="cif")
            except FileNotFoundError:
                self.generate_cif(cif_file)
                primitive_cell = ase.io.read(cif_file, format="cif")

            transformation_matrix = self.get_transformation_matrix(cif_file)

            super_cell = ase.build.make_supercell(
                primitive_cell, transformation_matrix
            )
            super_cell = self.get_compressed_cell(
                super_cell,
                stretch_factor=self.parameters.stretch_factor,
                density=self.parameters.mass_density,
                radius=self.parameters.WS_radius,
            )
            ase.io.write(
                self.supercell_file,
                super_cell,
                format="vasp",
                long_format=True,
            )
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
        elif (
            sum([stretch_factor is None, density is None, radius is None]) < 2
        ):
            print(
                "Warning: More than one of stretch factor, density, "
                "and radius is specified.\nRadius takes first priority, "
                "then density, finally stretch factor."
            )
        if density is not None:
            density_ambient = SuperCellProvider.get_mass_density(
                supercell, unit=units_density
            )
            stretch_factor = (density_ambient / density) ** (1.0 / 3.0)
        if radius is not None:
            radius_ambient = SuperCellProvider.get_wigner_seitz_radius(
                supercell
            )
            stretch_factor = radius / radius_ambient

        supercell.set_cell(
            supercell.get_cell() * stretch_factor, scale_atoms=True
        )

        return supercell

    @staticmethod
    def get_number_density(supercell, unit="Angstrom^3"):
        nr_atoms = len(supercell)
        volume_atoms = supercell.get_volume()
        nr_electrons = 0
        for i in range(0, nr_atoms):
            nr_electrons += supercell[i].number
        number_density = nr_electrons / volume_atoms
        angstrom3_in_m3 = 1 / (m * m * m)
        angstrom3_in_cm3 = angstrom3_in_m3 * 1000000
        if unit == "Angstrom^3":
            return number_density
        elif unit == "cm^3":
            return number_density / angstrom3_in_cm3
        else:
            raise Exception("Unit not implemented")

    @staticmethod
    def get_wigner_seitz_radius(supercell):
        number_density = SuperCellProvider.get_number_density(
            supercell, unit="Angstrom^3"
        )
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
            return mass_density / u_angstrom3_in_kg_m3
        elif unit == "g/(cm^3)":
            return mass_density / u_angstrom3_in_g_cm3
        elif unit == "u_Angstrom^3":
            return mass_density
        else:
            raise Exception("Unit not implemented")

    def generate_cif(self, cif_file):
        # read in the api key
        try:
            with open(self.parameters.mp_api_file, "r") as f:
                api_key = f.readlines()[0].strip()
        except FileNotFoundError:
            raise Exception(
                "File containing materials project API key not found.\n"
                "If you want to generate structures automatically, please "
                "register an account on the Materials Project website.\n"
                "Then save the API key in a file and set the file location "
                "using parameters.mp_api_file."
            )
        # initialize materials project api
        mpr = MPRester(api_key)
        # dictionary relating crystal structures to their space groups
        # TODO should this be moved elsewhere?
        lattice_space_groups = {"fcc": 225, "bcc": 229, "hcp": 194}
        # search materials project database for the desired element
        structures = mpr.summary.search(
            chemsys=self.parameters.element,
            fields=[
                "symmetry",
                "material_id",
                "theoretical",
                "energy_above_hull",
            ],
        )
        # narrow down by correct structure
        structures = list(
            filter(
                lambda x: x.symmetry.number
                == lattice_space_groups[self.parameters.crystal_structure],
                structures,
            )
        )
        if len(structures) == 0:
            raise Exception(
                "No structure found for element "
                + self.parameters.element
                + " with crystal structure"
                + self.parameters.crystal_structure
                + "\n."
                + "Please reconsider parameters or provide your own input files."
            )
        # narrow down by structures which are experimentally verified
        structures_expt = list(
            filter(lambda x: x.theoretical == False, structures)
        )
        # choose minimal energy structure
        if len(structures_expt) == 0:
            best_structure = min(structures, key=lambda x: x.energy_above_hull)
        elif len(structures_expt) == 1:
            best_structure = structures_expt[0]
        elif len(structures_expt) > 1:
            best_structure = min(
                structures_expt, key=lambda x: x.energy_above_hull
            )
        # get the ID of the best structure
        best_structure_id = best_structure.material_id
        structure = mpr.get_structure_by_material_id(
            best_structure_id, conventional_unit_cell=True
        )
        # write to the cif file
        structure.to(filename=cif_file)

    def get_transformation_matrix(self, cif_file):
        # extract the number of atoms in the unit cell
        nsites = self.get_nsites_from_cif(cif_file)
        if nsites is None:
            raise Exception(
                "cif file is not formatted correctly.\n"
                "It must contain _cell_formula_units_Z."
            )
        atoms_over_nsites = self.parameters.number_of_atoms / nsites
        if np.log2(atoms_over_nsites) % 1 != 0:
            raise Exception(
                "Number of atoms with this crystal structure "
                "is not supported. The ratio of these two values "
                "must be a power of 2."
            )
        # matrix which creates supercell
        transform_ratios = [1, 1, 1]
        i = 0
        while atoms_over_nsites > 1:
            transform_ratios[i] *= 2
            atoms_over_nsites /= 2
            i += 1
            if i > 2:
                i = 0

        transformation_matrix = [
            [transform_ratios[0], 0, 0],
            [0, transform_ratios[1], 0],
            [0, 0, transform_ratios[2]],
        ]

        return transformation_matrix

    @staticmethod
    def get_nsites_from_cif(cif_file):
        with open(cif_file, "r") as file:
            for line in file:
                if line.startswith("_cell_formula_units_Z"):
                    return int(line.split()[-1])
        # returns None if the line isn't found
        return None
