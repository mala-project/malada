import os
import numpy as np

from shutil import copyfile
import ase
import ase.io

from malada.utils.convergence_guesses import *
from malada.utils.custom_converter import *
from malada.utils.vasp_utils import VaspUtils
from .provider import Provider


class DFTConvergenceProvider(Provider):
    """For a given supercell and calculator, determine convergence."""

    def __init__(self, parameters, external_convergence_results = None,
                 external_convergence_folder = None):
        super(DFTConvergenceProvider, self).__init__(parameters)
        self.parameters = parameters
        self.external_convergence_results = external_convergence_results
        self.external_convergence_folder = external_convergence_folder
        self.convergence_results_file = None
        self.converged_cutoff = None
        self.converged_kgrid = None

    def provide(self, provider_path, supercell_file):
        file_name = self.parameters.element + \
                    "_" + str(self.parameters.number_of_atoms) + \
                    "_" + self.parameters.crystal_structure +\
                    "_" + str(self.parameters.temperature) +\
                    "_" + self.parameters.dft_calculator+".conv.xml"
        self.supercell_file = os.path.join(provider_path, file_name)

        # Check if there exist results or if we have to work from scratch.
        if self.external_convergence_results is None:
            if self.external_convergence_folder is None:
                # In this case we have to run convergence calculations.
                # First: cutoffs.
                self.__create_dft_convergence_inputs("cutoff", supercell_file,
                                                     provider_path)


            else:
                pass # TODO: Write analysis here.
        else:
            copyfile(self.external_convergence_results,
                     self.convergence_results_file)
            print("Getting <<convergence_results>>.pkl file from disc.")

    def __check_input_correctness(self, atoms):
        # Check if what we read made sense (assuming there was specified
        # user input.)
        if self.parameters.element is not None \
                and self.parameters.number_of_atoms is not None:
            if not np.all(
                    atoms.get_atomic_numbers() ==
                    atoms.get_atomic_numbers()):
                raise Exception(
                    "Mismatch between user input and provided file.")
            if atoms.get_chemical_symbols()[0] != \
                    self.parameters.element:
                raise Exception(
                    "Mismatch between user input and provided file.")
            if len(atoms) != self.parameters.number_of_atoms:
                raise Exception(
                    "Mismatch between user input and provided file.")


    def __create_dft_convergence_inputs(self, parameters_to_converge, posfile,
                                        working_directory):
        atoms_Angstrom = ase.io.read(posfile, format="vasp")
        self.__check_input_correctness(atoms_Angstrom)

        # This is needed for QE.
        scaling = int(np.round(
            atoms_Angstrom.cell.cellpar()[0] / atoms_Angstrom.cell.cellpar()[
                1]))

        # Determine parameters to converge.
        converge_list = []
        if parameters_to_converge == "kpoints":
            if self.parameters.dft_calculator == "qe":
                converge_list = kpoints_guesses[self.parameters.element]
            elif self.parameters.dft_calculator == "vasp":
                converge_list = kpoints_guesses[self.parameters.element]
            cutoff = self.converged_cutoff
        else:
            if self.parameters.dft_calculator == "qe":
                converge_list = cutoff_guesses_qe[self.parameters.element]
            elif self.parameters.dft_calculator == "vasp":
                converge_list = cutoff_guesses_vasp[self.parameters.element]
            kpoints = (1, 1, 1)

        # Create files for submission.
        for entry in converge_list:
            if parameters_to_converge == "kpoints":
                kpoints = entry
            else:
                cutoff = entry
            this_folder = self.__make_convergence_folder(kpoints, cutoff,
                                                         working_directory)

            qe_pseudopotentials = {self.parameters.element :
                                   self.parameters.pseudopotential["path"]}
            nbands = int(self.parameters.number_of_atoms *
                         self.parameters.pseudopotential["valence_electrons"]
                         * 1.05)

            # Get the run parameters.
            # TODO: Find some metric for the mixing!
            mixing = 0.1
            qe_input_data = {
                "occupations": 'smearing',
                "calculation": 'scf',
                "restart_mode": 'from_scratch',
                "prefix": self.parameters.element,
                "pseudo_dir": self.parameters.pseudopotential["path"],
                "outdir": "temp",
                "ibrav": 0,
                "smearing": 'fermi-dirac',
                "degauss": round(kelvin_to_rydberg(
                                 self.parameters.temperature), 7),
                "ecutrho": cutoff * 4,
                "ecutwfc": cutoff,
                "tstress": True,
                "tprnfor": True,
                "nbnd": nbands,
                "mixing_mode": "plain",
                "mixing_beta": mixing,
                "conv_thr": 1e-6 * self.parameters.number_of_atoms,
                "nosym": True,
                "noinv": True
            }

            vasp_input_data = {
                "ISTART": 0,
                "ENCUT": str(cutoff) + " eV",
                "EDIFF": 1e-6 * self.parameters.number_of_atoms *
                                ase.units.Rydberg,
                "ISMEAR": -1,
                "SIGMA": round(kelvin_to_eV(
                               self.parameters.temperature) *
                               ase.units.Rydberg, 7),
                "ISYM": 0,
                "NBANDS": nbands,
                "LREAL": "A",
            }

            # Write the convergence files.
            if self.parameters.dft_calculator == "qe":
                ase.io.write(this_folder + self.parameters.element
                             + ".pw.scf.in",
                             atoms_Angstrom,
                             "espresso-in", input_data=qe_input_data,
                             pseudopotentials= \
                             qe_pseudopotentials, kpts=kpoints)

            elif self.parameters.dft_calculator == "vasp":
                ase.io.write(this_folder + "POSCAR", atoms_Angstrom,
                             "vasp")
                VaspUtils.write_to_incar(this_folder, "INCAR", vasp_input_data)
                VaspUtils.write_to_kpoints(this_folder, "KPOINTS", kpoints)

    def __make_convergence_folder(self, kgrid, cutoff, base_folder):
            if self.parameters.dft_calculator == "qe":
                folder_name = str(cutoff) + "Ry_k" + str(kgrid[0]) + str(
                    kgrid[1]) + str(kgrid[2]) + "/"
            elif self.parameters.dft_calculator == "vasp":
                folder_name = str(cutoff) + "eV_k" + str(kgrid[0]) + str(
                    kgrid[1]) + str(kgrid[2]) + "/"
            this_folder = os.path.join(base_folder, folder_name)
            if not os.path.exists(this_folder):
                os.makedirs(this_folder)
            return this_folder

