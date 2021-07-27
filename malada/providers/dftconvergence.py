import os
from malada.utils.convergence_guesses import *

class DFTConvergenceProvider:
    """For a given supercell and calculator, determine convergence."""

    def __init__(self, parameters, external_convergence_results,
                 external_convergence_folder):
        super(DFTConvergenceProvider, self).__init__(parameters)
        self.parameters = parameters
        self.external_convergence_results = external_convergence_results
        self.external_convergence_folder = external_convergence_folder

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
            copyfile(self.external_cif_file, self.supercell_file)
            print("Getting <<supercell>>.vasp file from disc.")

    def __create_dft_convergence_inputs(self, parameters_to_converge, posfile,
                                        working_directory):
        atoms_Angstrom = ase.io.read(posfile, format="vasp")

        # Check if what we read made sense (assuming there was specified user input.)
        if element is not None and number_of_atoms is not None:
            if not np.all(
                    atoms_Angstrom.get_atomic_numbers() ==
                    atoms_Angstrom.get_atomic_numbers()):
                raise Exception(
                    "Mismatch between user input and provided file.")
            if atoms_Angstrom.get_chemical_symbols()[0] != \
                    self.parameters.element:
                raise Exception(
                    "Mismatch between user input and provided file.")
            if len(atoms_Angstrom) != self.parameters.number_of_atoms:
                raise Exception(
                    "Mismatch between user input and provided file.")
        # This is needed for QE.
        scaling = int(np.round(
            atoms_Angstrom.cell.cellpar()[0] / atoms_Angstrom.cell.cellpar()[
                1]))

        # Determine parameters to converge.
        if parameter_to_converge == "kpoints":
            if self.parameters.element == "Fe":
                basic_list = []
                converge_list = []
                for i in basic_list:
                    converge_list.append((i * scaling, i, i))
            if self.parameters.element == "Be":
                basic_list = []
                converge_list = []
                for i in basic_list:
                    converge_list.append((i * scaling, i, i))
        else:
            if self.parameters.element == "Fe":
                if calculator == "qe":
                    converge_list = cutoff_guesses_qe["Fe"]
                elif calculator == "vasp":
                    converge_list = cutoff_guesses_vasp["Fe"]
                kpoints = (1, 1, 1)
            if self.parameters.element == "Be":
                converge_list = [40, 50, 60, 70]
                kpoints = (1, 1, 1)
        # Create files for submission.
        for entry in converge_list:
            if parameter_to_converge == "kpoints":
                kpoints = entry
            else:
                cutoff = entry
                kpoints = (1, 1, 1)
            this_folder = self.__make_convergence_folder(kpoints, cutoff, base_folder,
                                                  calculator)

            qe_pseudopotentials = {self.parameters.element :
                                   self.parameters.pseudopotential["path"]}
            nbands = int(self.parameters.number_of_atoms*
                         self.parameters.pseudopotential["valence_electrons"]
                         *1.05)

            mixing = 0.1
            qe_input_data = {
                "occupations": 'smearing',
                "calculation": 'scf',
                "restart_mode": 'from_scratch',
                "prefix": element,
                "pseudo_dir": pspath,
                "outdir": base_data_folder,
                "ibrav": 0,
                "smearing": 'fermi-dirac',
                "degauss": round(get_smearing_from_temperature(temperature),
                                 7),
                "ecutrho": cutoff * 4,
                "ecutwfc": cutoff,
                "tstress": True,
                "tprnfor": True,
                "nbnd": nbands,
                "mixing_mode": "plain",
                "mixing_beta": mixing,
                "conv_thr": 1e-6 * number_of_atoms,
                "nosym": True,
                "noinv": True
            }

            vasp_input_data = {
                "ISTART": 0,
                "ENCUT": str(cutoff) + " eV",
                "EDIFF": 1e-6 * number_of_atoms * ase.units.Rydberg,
                "ISMEAR": -1,
                "SIGMA": round(get_smearing_from_temperature(
                    temperature) * ase.units.Rydberg, 7),
                "ISYM": 0,
                "NBANDS": nbands,
                "LREAL": "A",
            }

            if calculator == "qe":
                ase.io.write(this_folder + element + ".pw.scf.in",
                             atoms_Angstrom,
                             "espresso-in", input_data=qe_input_data,
                             pseudopotentials= \
                                 qe_pseudopotentials, kpts=kpoints)
            elif calculator == "vasp":
                ase.io.write(this_folder + "POSCAR", atoms_Angstrom,
                             "vasp")
                write_to_incar(this_folder, "INCAR", vasp_input_data)
                write_to_kpoints(this_folder, "KPOINTS", kpoints)

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

