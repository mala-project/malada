"""Provider for optimized DFT calculation parameters."""
import os
import numpy as np
import glob

from shutil import copyfile
import ase
import ase.io
from ase.units import Rydberg
from xml.etree.ElementTree import Element, SubElement, tostring
from xml.dom import minidom

import malada
from malada.utils.convergence_guesses import *
from malada.utils.custom_converter import *
from malada.utils.vasp_utils import VaspUtils
from .provider import Provider


class DFTConvergenceProvider(Provider):
    """
    For a given supercell and calculator, determine convergence.

    Parameters
    ----------
    parameters : malada.utils.parametes.Parameters
        Parameters used to create this object.

    external_convergence_results : string
        Path to xml file containing previously calculated convergence
        results. If not None, no DFT caclulations will be done.

    external_convergence_folder : string
        Path to a folder containing already calculated DFT convergence
        results. If not none, no DFT calculations will be done.

    predefined_kgrid : tuple
        Tuple in the form (kx,ky,kz). If not None, this k grid will be
        used and no attempt will be made to find a more optimal one.

    predefined_cutoff : float
        Kinetic energy cutoff. If not None, this cutoff will be used and
        no attempt will be made to find a more optimal one.
    """

    def __init__(self, parameters, external_convergence_results=None,
                 external_convergence_folder=None,
                 predefined_kgrid=None, predefined_cutoff=None):
        super(DFTConvergenceProvider, self).__init__(parameters)
        self.parameters = parameters
        self.external_convergence_results = external_convergence_results
        self.external_convergence_folder = external_convergence_folder
        self.convergence_results_file = None
        self.converged_cutoff = predefined_cutoff
        self.converged_kgrid = predefined_kgrid

    def provide(self, provider_path, supercell_file, create_submit_script=True):
        """
        Provide DFT parameters converged to within user specification.

        The cutoff energy (=basis set size) and k-grid will be optimized.

        Parameters
        ----------
        provider_path : string
            Path in which to operate in.

        supercell_file
        """
        file_name = self.parameters.element + \
                    str(self.parameters.number_of_atoms) + \
                    "_" + self.parameters.crystal_structure +\
                    "_" + str(self.parameters.temperature) +\
                    "_" + self.parameters.dft_calculator+".conv.xml"
        self.convergence_results_file = os.path.join(provider_path, file_name)

        # Instantiate a runner.
        dft_runner = malada.RunnerInterface(self.parameters)

        # Check if there exist results or if we have to work from scratch.
        if self.external_convergence_results is None:
            if self.external_convergence_folder is None:
                # In this case we have to run convergence calculations.
                # First: cutoffs.
                # Create, run and analyze.
                cutoff_try = 0
                kpoint_try = 0
                while (self.converged_cutoff is None and cutoff_try
                       < self.parameters.maximum_cutoff_try):

                    # First, we create the inputs.
                    cutoff_folders = self.\
                        __create_dft_convergence_inputs("cutoff", supercell_file,
                                                        provider_path, cutoff_try)

                    # Then we run.
                    for cutoff_folder in cutoff_folders:
                        print("Running DFT in", cutoff_folder)
                        dft_runner.run_folder(cutoff_folder,
                                              "dft")

                    # Afterwards, we check if running was succesful.
                    if self.__check_run_success(provider_path, "cutoff", fixed_kpoints=(1,1,1)):
                        self.converged_cutoff = self.__analyze_convergence_runs(provider_path, "cutoff",
                                                                                supercell_file,
                                                                                fixed_kpoints=(1, 1, 1))
                    else:
                        if self.parameters.run_system == "slurm_creator":
                            all_submit_file = open(
                                os.path.join(provider_path, "submit_all.sh"), mode='w')
                            all_submit_file.write("#!/bin/bash\n\n")
                            for cutoff_folder in cutoff_folders:
                                all_submit_file.write("cd "+cutoff_folder.split("/")[-2]+"\n")
                                all_submit_file.write("sbatch submit.slurm\n")
                                all_submit_file.write("cd ..\n")
                            all_submit_file.close()
                            print("Run scripts created, please run via slurm.\n"
                                  "Quitting now.")
                            quit()
                        else:
                            raise Exception("DFT calculations failed.")
                    if self.converged_cutoff is None:
                        print("Could not find an aedaquate cutoff energy, "
                              "trying again with larger cutoff energies. "
                              "This will be try nr. "
                              +str(cutoff_try+2))
                    cutoff_try += 1

                # Second: cutoffs.
                # Create, run and analyze.
                while (self.converged_kgrid is None and kpoint_try
                       < self.parameters.maximum_kpoint_try):
                    kpoints_folders = self. \
                        __create_dft_convergence_inputs("kpoints", supercell_file,
                                                        provider_path, kpoint_try)
                    for kpoint_folder in kpoints_folders:
                        print("Running DFT in", kpoint_folder)
                        dft_runner.run_folder(kpoint_folder,
                                              "dft")

                    # Afterwards, we check if running was succesful.
                    if self.__check_run_success(provider_path, "kpoints", fixed_cutoff=self.converged_cutoff):
                        self.converged_kgrid = self.__analyze_convergence_runs(
                            provider_path,
                            "kpoints", supercell_file,
                            fixed_cutoff=self.converged_cutoff)
                    else:
                        if self.parameters.run_system == "slurm_creator":
                            all_submit_file = open(
                                os.path.join(provider_path, "submit_all.sh"), mode='w')
                            all_submit_file.write("#!/bin/bash\n\n")
                            for kpoint_folder in kpoints_folders:
                                all_submit_file.write("cd "+kpoint_folder.split("/")[-2]+"\n")
                                all_submit_file.write("sbatch submit.slurm\n")
                                all_submit_file.write("cd ..\n")
                            all_submit_file.close()

                            print("Run scripts created, please run via slurm.\n"
                                  "Quitting now.")
                            quit()
                        else:
                            raise Exception("DFT calculations failed.")

                    if self.converged_kgrid is None:
                        print("Could not find an aedaquate k-grid, trying again"
                              "with larger k-grids. This will be try nr. "
                              +str(kpoint_try+2))
                    kpoint_try += 1
            else:
                # TODO: Add a consistency check here. It could be the provided
                # folder does not match!
                print("Reading precalculated convergence results.")
                # Analyze the convergence analsyis.
                # First: cutoffs.
                if self.converged_cutoff is None:
                    self.converged_cutoff = self.__analyze_convergence_runs(self.external_convergence_folder, "cutoff",
                                                                            supercell_file, fixed_kpoints=(1, 1, 1))
                if self.converged_cutoff is None:
                    raise Exception("Provided convergence data not sufficient,"
                                    "please perform additional calculations"
                                    " with larger cutoff energies.")

                # Second: kgrid.
                if self.converged_kgrid is None:
                    self.converged_kgrid = self.__analyze_convergence_runs(self.external_convergence_folder,
                                                                           "kpoints",supercell_file,
                                                                            fixed_cutoff=self.converged_cutoff)
                if self.converged_kgrid is None:
                    raise Exception("Provided convergence data not sufficient,"
                                    "please perform additional calculations"
                                    " with larger k-grids.")
            # Print the output.
            unit = "eV" if self.parameters.dft_calculator == "vasp" else "Ry"
            print("Converged energy cutoff: ", self.converged_cutoff, unit)
            print("Converged k-grid: ", self.converged_kgrid)

            # Write the output to xml.
            self.__write_to_xml()
        else:
            copyfile(self.external_convergence_results,
                     self.convergence_results_file)
            print("Getting <<convergence_results>>.xml file from disc.")


    def __write_to_xml(self):
        top = Element('dftparameters')
        calculationparameters = SubElement(top, "calculationparameters")
        cpnode = SubElement(calculationparameters, "element",
                            {"type": "string"})
        cpnode.text = self.parameters.element
        cpnode = SubElement(calculationparameters, "number_of_atoms",
                          {"type": "int"})
        cpnode.text = str(self.parameters.number_of_atoms)
        cpnode = SubElement(calculationparameters, "crystal_structure",
                                  {"type": "string"})
        cpnode.text = self.parameters.crystal_structure
        cpnode = SubElement(calculationparameters, "temperature",
                                  {"type": "float"})
        cpnode.text = str(self.parameters.temperature)
        cpnode = SubElement(calculationparameters, "dft_calculator",
                                  {"type": "string"})
        cpnode.text = self.parameters.dft_calculator
        cutoff = SubElement(top, "cutoff_energy", {'type': "int"})
        cutoff.text = str(self.converged_cutoff)
        kpoints = SubElement(top, "kpoints")
        kx = SubElement(kpoints, "kx", {'type': "int"})
        kx.text = str(self.converged_kgrid[0])
        ky = SubElement(kpoints, "ky", {'type': "int"})
        ky.text = str(self.converged_kgrid[1])
        kz = SubElement(kpoints, "kz", {'type': "int"})
        kz.text = str(self.converged_kgrid[2])
        rough_string = tostring(top, 'utf-8')
        reparsed = minidom.parseString(rough_string)
        with open(self.convergence_results_file, "w") as f:
            f.write(reparsed.toprettyxml(indent="  "))


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
                    converge_list.append(self.__k_edge_length_to_grid(grids, atoms_Angstrom))
            else:
                spacing = kpoints_guesses[self.parameters.element][-1]-\
                          kpoints_guesses[self.parameters.element][-2]
                converge_list.append(self.__k_edge_length_to_grid(kpoints_guesses[self.parameters.element][-1]+spacing*(try_number-1)*2+1*spacing, atoms_Angstrom))
                converge_list.append(self.__k_edge_length_to_grid(kpoints_guesses[self.parameters.element][-1]+spacing*(try_number-1)*2+2*spacing, atoms_Angstrom))
            cutoff = self.converged_cutoff
        else:
            if self.parameters.dft_calculator == "qe":
                if try_number == 0:
                    converge_list = cutoff_guesses_qe[self.parameters.element]
                else:
                    spacing = cutoff_guesses_qe[self.parameters.element][-1] - \
                              cutoff_guesses_qe[self.parameters.element][-2]
                    converge_list.append(
                        cutoff_guesses_qe[self.parameters.element][
                            -1] + spacing * (try_number - 1) * 2 + 1*spacing)
                    converge_list.append(
                        cutoff_guesses_qe[self.parameters.element][
                            -1] + spacing * (try_number - 1) * 2 + 2*spacing)


            elif self.parameters.dft_calculator == "vasp":
                if try_number == 0:
                    converge_list = cutoff_guesses_vasp[self.parameters.element]
                else:
                    spacing = cutoff_guesses_vasp[self.parameters.element][-1] - \
                              cutoff_guesses_vasp[self.parameters.element][-2]
                    converge_list.append(
                        cutoff_guesses_vasp[self.parameters.element][
                            -1] + spacing * (try_number - 1) * 2 + 1*spacing)
                    converge_list.append(
                        cutoff_guesses_vasp[self.parameters.element][
                            -1] + spacing * (try_number - 1) * 2 + 2*spacing)
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
            if self.parameters.dft_calculator == "qe":
                qe_pseudopotentials = {self.parameters.element :
                                       self.parameters.pseudopotential["name"]}

            # TODO: Find a better way to calculate this
            # With increasing temperature, this breaks of fairly quickly.
            nbands = self._get_number_of_bands()

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
                "conv_thr": self.parameters.dft_scf_accuracy_per_atom_Ry * self.parameters.number_of_atoms,
                "nosym": True,
                "noinv": True,
                "verbosity": "high"
            }

            vasp_input_data = {
                "ISTART": 0,
                "ENCUT": str(cutoff) + " eV",
                "EDIFF": self.parameters.dft_scf_accuracy_per_atom_Ry*
                         self.parameters.number_of_atoms *
                                ase.units.Rydberg,
                "ISMEAR": -1,
                "SIGMA": round(kelvin_to_eV(
                               self.parameters.temperature), 7),
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
                VaspUtils.write_to_incar(this_folder, vasp_input_data)
                VaspUtils.write_to_kpoints(this_folder, kpoints)
                VaspUtils.write_to_potcar_copy(this_folder,
                                               os.path.join(self.parameters.pseudopotential["path"],
                                                            self.parameters.pseudopotential["name"]))


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

    def __check_run_success(self, base_folder, parameter_to_converge,
                            fixed_cutoff=None, fixed_kpoints=None):
        # First, parse the results out of the output files.
        if parameter_to_converge == "cutoff":
            if self.parameters.dft_calculator == "qe":
                converge_list = glob.glob(os.path.join(base_folder,
                                          "*Ry_k"+str(fixed_kpoints[0]) +
                                          str(fixed_kpoints[1]) +
                                          str(fixed_kpoints[2])))
            elif self.parameters.dft_calculator == "vasp":
                converge_list = glob.glob(os.path.join(base_folder,
                                          "*eV_k"+str(fixed_kpoints[0]) +
                                          str(fixed_kpoints[1]) +
                                          str(fixed_kpoints[2])))
        elif parameter_to_converge == "kpoints":
            if self.parameters.dft_calculator == "qe":
                converge_list = glob.glob(os.path.join(base_folder,
                                                       str(fixed_cutoff)+"Ry_k*"))
            elif self.parameters.dft_calculator == "vasp":
                converge_list = glob.glob(os.path.join(base_folder,
                                                       str(fixed_cutoff)+"eV_k*"))
        else:
            raise Exception("Unknown convergence parameter.")

        for entry in converge_list:
            if self.parameters.dft_calculator == "qe":
                convergence_file_lists = glob.glob(os.path.join(entry,
                                                   "*.out"))
                if len(convergence_file_lists) != 1:
                    return False
            elif self.parameters.dft_calculator == "vasp":
                if not os.path.isfile(os.path.join(entry, "OUTCAR")):
                    return False
            else:
                raise Exception("Unknown calculator chosen.")
        return True

    def __analyze_convergence_runs(self, base_folder, parameter_to_converge,
                                   supercell_file,
                                   fixed_cutoff=None, fixed_kpoints=None):
        atoms_Angstrom = ase.io.read(supercell_file, format="vasp")

        # First, parse the results out of the output files.
        result_list = []
        if parameter_to_converge == "cutoff":
            if self.parameters.dft_calculator == "qe":
                converge_list = glob.glob(os.path.join(base_folder,
                                          "*Ry_k"+str(fixed_kpoints[0]) +
                                          str(fixed_kpoints[1]) +
                                          str(fixed_kpoints[2])))
            elif self.parameters.dft_calculator == "vasp":
                converge_list = glob.glob(os.path.join(base_folder,
                                          "*eV_k"+str(fixed_kpoints[0]) +
                                          str(fixed_kpoints[1]) +
                                          str(fixed_kpoints[2])))
        elif parameter_to_converge == "kpoints":
            if self.parameters.dft_calculator == "qe":
                converge_list = glob.glob(os.path.join(base_folder,
                                                       str(fixed_cutoff)+"Ry_k*"))
            elif self.parameters.dft_calculator == "vasp":
                converge_list = glob.glob(os.path.join(base_folder,
                                                       str(fixed_cutoff)+"eV_k*"))
        else:
            raise Exception("Unknown convergence parameter.")

        for entry in converge_list:
            energy = None
            if parameter_to_converge == "cutoff":
                if self.parameters.dft_calculator == "qe":
                    argument = int(os.path.basename(entry).split("Ry")[0])
                elif self.parameters.dft_calculator == "vasp":
                    argument = int(os.path.basename(entry).split("eV")[0])
                else:
                    raise Exception("Unknown calculator chosen.")

            elif parameter_to_converge == "kpoints":
                k_argument = os.path.basename(entry).split("k")[1]
                if len(k_argument) == 3:
                    argument = (
                    int((os.path.basename(entry).split("k")[1])[0]),
                    int((os.path.basename(entry).split("k")[1])[1]),
                    int((os.path.basename(entry).split("k")[1])[2]))
                else:
                    # All three dimensions are double digits.
                    if len(k_argument) == 6:
                        argument = (
                            int((os.path.basename(entry).split("k")[1])[0:2]),
                            int((os.path.basename(entry).split("k")[1])[2:4]),
                            int((os.path.basename(entry).split("k")[1])[4:6]))
                    else:
                        # One of the k-dimensions is a double digit.
                        # We will attempt to recover by making an educated
                        # guess which one it is.
                        scaling_x = int(np.ceil(
                            atoms_Angstrom.cell.cellpar()[0] /
                            atoms_Angstrom.cell.cellpar()[
                                2]))
                        scaling_y = int(np.ceil(
                            atoms_Angstrom.cell.cellpar()[1] /
                            atoms_Angstrom.cell.cellpar()[
                                2]))
                        argument_x = int(
                            (os.path.basename(entry).split("k")[1])[0])
                        argument_y = int(
                            (os.path.basename(entry).split("k")[1])[1])
                        argument_z = int(
                            (os.path.basename(entry).split("k")[1])[2])

                        if scaling_x != 1:
                            argument_x = int(
                                (os.path.basename(entry).split("k")[1])[0:2])
                            argument_y = int(
                                (os.path.basename(entry).split("k")[1])[2:3])

                        if scaling_y != 1:
                            if scaling_x == 1:
                                argument_y = int(
                                    (os.path.basename(entry).split("k")[1])[1:3])
                            else:
                                argument_y = int(
                                    (os.path.basename(entry).split("k")[1])[2:4])
                        argument = (argument_x, argument_y, argument_z)

            else:
                raise Exception("Unknown parameter chosen.")

            if self.parameters.dft_calculator == "qe":
                convergence_file_lists = glob.glob(os.path.join(entry,
                                                   "*.out"))
                if len(convergence_file_lists) > 1:
                    print(convergence_file_lists)
                    raise Exception("Run folder with ambigous content.")
                if len(convergence_file_lists) > 0:
                    energy = self.__get_qe_energy(convergence_file_lists[0])

            elif self.parameters.dft_calculator == "vasp":
                energy = self.__get_vasp_energy(os.path.join(entry,
                                                                   "OUTCAR"))

            else:
                raise Exception("Unknown calculator chosen.")
            if energy is not None:
                result_list.append([argument, energy])
        result_list = sorted(result_list, key=lambda d: [d[0]])

        # Next, analyze the convergence.
        best_param = None
        print("Convergence results for: ", parameter_to_converge)
        for i in range(1, len(result_list)):
            print(result_list[i-1][0], result_list[i][0],
                  np.abs(result_list[i][1] - result_list[i-1][1]))
            if np.abs(result_list[i][1] - result_list[i-1][1]) < self.\
                    parameters.dft_conv_accuracy_meVperatom:
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
                if "total energy" in l and "is F=E-TS" not in l and \
                    "is the sum of" not in l:
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

    def __k_edge_length_to_grid(self, edge_length, atoms_Angstrom):
        # TODO: Reflect geometry here.
        scaling_x = int(np.ceil(
            atoms_Angstrom.cell.cellpar()[0] / atoms_Angstrom.cell.cellpar()[
                2]))
        scaling_y= int(np.ceil(
            atoms_Angstrom.cell.cellpar()[1] / atoms_Angstrom.cell.cellpar()[
                2]))
        return (edge_length*scaling_x, edge_length*scaling_y, edge_length)
