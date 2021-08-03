import os
import numpy as np
import glob

from shutil import copyfile
import ase
import ase.io
from ase.units import Rydberg

import malada
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

        # Instantiate a runner.
        dft_runner = malada.BashRunner()

        # Check if there exist results or if we have to work from scratch.
        if self.external_convergence_results is None:
            if self.external_convergence_folder is None:
                # In this case we have to run convergence calculations.
                # First: cutoffs.
                # Create, run and analyze.
                cutoff_try = 0
                kpoint_try = 0
                cutoff_folders = self.\
                    __create_dft_convergence_inputs("cutoff", supercell_file,
                                                    provider_path, cutoff_try)
                # for cutoff_folder in cutoff_folders:
                #     print("Running DFT in", cutoff_folder)
                #     dft_runner.run_folder(cutoff_folder,
                #                           self.parameters.dft_calculator)
                self.converged_cutoff = self.__analyze_convergence_runs(provider_path, "cutoff",
                                                                        fixed_kpoints=(1, 1, 1))

                while (self.converged_kgrid is None and kpoint_try
                       < self.parameters.maximum_kpoint_try):
                    kpoints_folders = self. \
                        __create_dft_convergence_inputs("kpoints", supercell_file,
                                                        provider_path, kpoint_try)
                    for kpoint_folder in kpoints_folders:
                        print("Running DFT in", kpoint_folder)
                        dft_runner.run_folder(kpoint_folder,
                                              self.parameters.dft_calculator)
                    self.converged_kgrid = self.__analyze_convergence_runs(provider_path,
                                                                           "kpoints",
                                                                            fixed_cutoff=self.converged_cutoff)
                    if self.converged_kgrid is None:
                        print("Could not find an aedaquate k-grid, trying again"
                              "with larger k-grids. This will be try nr. "
                              +str(kpoint_try+1))
                    kpoint_try += 1
                print(self.converged_cutoff, self.converged_kgrid)
            else:
                # In this case we just have to run the analysis on the folder.
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
                                        working_directory, try_number):
        atoms_Angstrom = ase.io.read(posfile, format="vasp")
        self.__check_input_correctness(atoms_Angstrom)

        # This is needed for QE.
        scaling = int(np.round(
            atoms_Angstrom.cell.cellpar()[0] / atoms_Angstrom.cell.cellpar()[
                1]))

        # Determine parameters to converge.
        converge_list = []
        if parameters_to_converge == "kpoints":
            if try_number == 0:
                for grids in kpoints_guesses[self.parameters.element]:
                    converge_list.append(self.__k_edge_length_to_grid(grids))
            else:
                spacing = kpoints_guesses[self.parameters.element][-1]-\
                          kpoints_guesses[self.parameters.element][-2]
                converge_list.append(self.__k_edge_length_to_grid(kpoints_guesses[self.parameters.element][-1]+spacing*try_number))
                converge_list.append(self.__k_edge_length_to_grid(kpoints_guesses[self.parameters.element][-1]+spacing*(try_number+1)))
            cutoff = self.converged_cutoff
        else:
            if self.parameters.dft_calculator == "qe":
                converge_list = cutoff_guesses_qe[self.parameters.element]
            elif self.parameters.dft_calculator == "vasp":
                converge_list = cutoff_guesses_vasp[self.parameters.element]
            kpoints = (1, 1, 1)

        # Create files for submission.
        converge_folder_list = []
        for entry in converge_list:
            if parameters_to_converge == "kpoints":
                kpoints = entry
            else:
                cutoff = entry
            this_folder = self.__make_convergence_folder(kpoints, cutoff,
                                                         working_directory)
            converge_folder_list.append(this_folder)
            qe_pseudopotentials = {self.parameters.element :
                                   self.parameters.pseudopotential["name"]}
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
                "conv_thr": self.parameters.dft_scf_accuracy_per_atom * self.parameters.number_of_atoms,
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

        return converge_folder_list

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

    def __analyze_convergence_runs(self, base_folder, parameter_to_converge,
                                  fixed_cutoff=None, fixed_kpoints=None):

        # First, parse the results out of the output files.
        result_list = []
        if parameter_to_converge == "cutoff":
            converge_list = glob.glob(os.path.join(base_folder,
                                      "*Ry_k"+str(fixed_kpoints[0]) +
                                      str(fixed_kpoints[1]) +
                                      str(fixed_kpoints[2])))
        elif parameter_to_converge == "kpoints":
            converge_list = glob.glob(os.path.join(base_folder,
                                                   str(fixed_cutoff)+"Ry_k*"))
        else:
            raise Exception("Unknown convergence parameter.")

        for entry in converge_list:
            if self.parameters.dft_calculator == "qe":
                if parameter_to_converge == "cutoff":
                   argument = int(os.path.basename(entry).split("Ry")[0])
                elif parameter_to_converge == "kpoints":
                    argument = (int((os.path.basename(entry).split("k")[1])[0]),
                                int((os.path.basename(entry).split("k")[1])[1]),
                                int((os.path.basename(entry).split("k")[1])[2]))
                convergence_file_lists = glob.glob(os.path.join(entry,
                                                   "*.out"))
                if len(convergence_file_lists) != 1:
                    print(convergence_file_lists)
                    raise Exception("Run folder with ambigous content.")
                energy = self.__get_qe_energy(convergence_file_lists[0])

            elif self.parameters.dft_calculator == "vasp":
                argument = int(os.path.basename(entry).split("eV")[0])
                energy = self.__get_vasp_energy(os.path.join(entry,
                                                               "OUTCAR"))
            else:
                raise Exception("Unknown calculator chosen.")
            result_list.append([argument, energy])
        result_list = sorted(result_list, key=lambda d: [d[0]])

        # Next, analyze the convergence.
        for i in range(1, len(result_list)):
            if np.abs(result_list[i][1] - result_list[i-1][1]) < self.\
                    parameters.dft_accuracy_meVperatom:
                if best_param is not None:
                    break
                else:
                    best_param = result_list[i - 1][0]
            else:
                best_param = None
        return best_param

    def __get_qe_energy(self, file):
        energy = np.nan
        convergence_achieved = True
        with open(file, "r") as f:
            ll = f.readlines()
            for l in ll:
                if "total energy" in l and "is F=E-TS" not in l:
                    energy = float((l.split('=')[1]).split('Ry')[0])
                if "convergence NOT achieved" in l or "oom-kill" in l or\
                        "CANCELLED" in l or "BAD TERMINATION" in l:
                    convergence_achieved = False
        if convergence_achieved is False:
            raise Exception("Convergence was not achieved at", file)

        return float((energy * Rydberg * 1000)/self.parameters.number_of_atoms)

    def __get_vasp_energy(self, file):
        energy = np.nan
        convergence_achieved = True
        found_end = False
        with open(file, "r") as f:
            ll = f.readlines()
            for l in ll:
                if "FREE ENERGIE" in l:
                    found_end = True
                if "TOTEN" in l and found_end:
                    energy = float((l.split('=')[1]).split('eV')[0])
        if convergence_achieved is False:
            raise Exception("Convergence was not achieved at", file)

        return float((energy * 1000)/self.parameters.number_of_atoms)

    def __k_edge_length_to_grid(self, edge_length):
        # TODO: Reflect geometry here.
        return (edge_length, edge_length, edge_length)
