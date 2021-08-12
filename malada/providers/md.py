"""Provider for DFT-MD calculations."""

from .provider import Provider
import os
from shutil import copyfile
import xml.etree.ElementTree as ET
import ase
from ase.units import Rydberg
import ase.io
from ..utils import kelvin_to_rydberg, second_to_rydberg_time, kelvin_to_eV
from malada.utils.vasp_utils import VaspUtils
import malada


class MDProvider(Provider):
    """
    Performs a DFT-MD calculation and provides an ASE trjactory.

    This can be done with a number of runners.

    Parameters
    ----------
    parameters : malada.utils.parametes.Parameters
        Parameters used to create this object.

    external_trajectory : string
        Path to a file containing an ASE trajectory.
        If a MD trajectory has already been calculated, it can be provided
        here. In this case no MD calculation or processing is done whatsoever.
        Will be ignored if external_temperatures is None.

    external_temperatures : string
        Path to a file containing a numpy array with temperature values.
        If an MD trajectory has already been calculated, the corresponding
        temperatures can be provided here. In this case no MD calculation or
        processing is done whatsoever.
        Will be ignored if external_trajectory is None.

    external_run_folder : string
        Path to a folder containing a finished MD calculation.
        If provided, this MD calculation will be preprocessed by this
        file.
    """

    def __init__(self, parameters, external_trajectory=None,
                 external_temperatures=None, external_run_folder=None):
        super(MDProvider, self).__init__(parameters)
        self.external_trajectory = external_trajectory
        self.external_temperatures = external_temperatures
        self.external_run_folder = external_run_folder
        self.trajectory_file = None
        self.temperature_file = None

    def provide(self, provider_path, supercell_file, dft_convergence_file,
                md_performance_file):
        """
        Provide a MD trajectory and temperatures from previous steps.

        Depending on arguments at creation, no actual MD calculation
        may be necessary.

        Parameters
        ----------
        provider_path : string
            Path in which to operate in.

        supercell_file : string
            Path to file containing structue of the supercell in VASP format.

        dft_convergence_file : string
            Path to xml file containing the DFT convergence parameter.

        md_performance_file : string
            Path to xml file containing the optimal run parameters for
            MD.
        """
        # Creating file for output.
        file_name = self.parameters.element + \
                    str(self.parameters.number_of_atoms) + \
                    "_" + self.parameters.crystal_structure +\
                    "_" + str(self.parameters.temperature)
        self.trajectory_file = os.path.join(provider_path, file_name+".traj")
        self.temperature_file = os.path.join(provider_path, file_name +
                                             ".temp.npy")

        # MD performance has to be loaded if SLURM is used.
        # TODO: Fix this, the MD parameters should also contain the thermostat
        # controller
        if self.parameters.run_system == "slurm_creator":
            self.parameters.md_slurm = malada.SlurmParameters.\
                                       from_xml(md_performance_file)
        if self.external_trajectory is None or self.external_temperatures is None:
            if self.external_run_folder is None:
                # Here, we have to perform the MD first.
                # First, create MD inputs.
                self.__create_md_run(dft_convergence_file,
                                     supercell_file, provider_path)

                # Run the MD calculation.
                mdrunner = malada.RunnerInterface(
                    self.parameters)
                mdrunner.run_folder(provider_path,"md",
                                    qe_input_type="*.pw.md.in")
                folder_to_parse = provider_path

                # If SLURM is used as run tool, we stop here while the user
                # runs.
                if self.parameters.run_system == "slurm_creator":
                    print("Created run scripts. Please run via slurm.")
                    print("Quitting...")
                    quit()
            else:
                folder_to_parse = self.external_run_folder

            # Now we parse the results.
            self._qe_out_to_trajectory(folder_to_parse, self.trajectory_file)
            self._qe_out_to_temperature(folder_to_parse, self.temperature_file)

        else:
            copyfile(self.external_trajectory, self.trajectory_file)
            copyfile(self.external_temperatures, self.temperature_file)
            print("Getting <<trajectory>>.traj and <<temperatures>>.temp.npy"
                  " files from disc.")

    def __create_md_run(self, dft_convergence_file,
                        posfile, base_path):
        cutoff, kgrid = self._read_convergence(dft_convergence_file)
        if self.parameters.md_at_gamma_point:
            kgrid = (1, 1, 1)
            if self.parameters.dft_calculator == "qe":
                kgrid = None

        qe_pseudopotentials = {self.parameters.element:
                                   self.parameters.pseudopotential["name"]}
        nbands = self._get_number_of_bands()

        qe_input_data = {
            "occupations": 'smearing',
            "calculation": 'md',
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
            "mixing_beta": 0.1,
            "conv_thr": self.parameters.dft_scf_accuracy_per_atom_Ry *
                        self.parameters.number_of_atoms,
            # "verbosity" : "high", # This is maybe a bit high
            "nosym": True,
            # MD variables - these are not final in any way.
            # Choice of thermostat. We want a constant temperature throughout.
            # Looking at this brief note here:
            # https://www2.mpip-mainz.mpg.de/~andrienk/journal_club/thermostats.pdf
            # it is apparent that we would need Nose-Hover, but the next
            # best thing QE gives us is berendsen.
            "ion_temperature": 'berendsen',
            "tempw": self.parameters.temperature,
            # Time step: In Ryberg atomic units.
            # For now: 1 fs
            "dt": second_to_rydberg_time(1e-15*self.parameters.time_step_fs),
            # Number of steps.
            # For now: 500.
            "nstep": self.parameters.maximum_number_of_timesteps,
            # Control how the ionic dynamics are calculated.
            # I think we want verlet, which is also the default.
            "ion_dynamics": "verlet",
            # I think this helps with performance.
            "wfc_extrapolation": 'second order',
            # This controls the velocity rescaling performed by the Berendsen
            # thermostat. It corresponds to tau/dt in e.g. eq. 1.12
            # of https://link.springer.com/content/pdf/10.1007%2F978-3-540-74686-7_1.pdf
            # By default it is 1, but I think we have to increase it a bit.
            "nraise": 1#self.parameters.md_thermostat_controller,
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
            # MD variables - these are not final in any way.
            # Interpolation between time steps.
            "IWAVPR": 2,
            "IWAVPRE": 12,
            # Precision. For our purposes, normal should suffice.
            "PREC": "NORMAL",
            # Mixing algorithm, 48 is RMM-DIIS with preconditioning
            "IALGO": 48,
            # Minimum number of electronic steps, the default is 2 and should
            # be sligthly higher for MD runs.
            "NELMIN": 6,
            # These are some missing parameters that Kushal and Mani used
            # and I am copying them because I assume they are beneficial for
            # performance.
            "BMIX": 2,
            "MAXMIX": 50,
            # Selects MD (and Nose Hoover thermostat for NVT ensemble)
            "IBRION": 0,
            "MDALGO": 2,
            # We want to do 10000 time steps
            "NSW": 10000,
            # This enables the Nose Hoover thermostat, it is controlled by this
            # parameter
            "SMASS": self.parameters.md_thermostat_controller,

            # Time step in fs, we want 1 fs.
            "POTIM": self.parameters.time_step_fs,

            # TEBEG: Initial temperature.
            "TEBEG": self.parameters.temperature,

            # This surpresses large disk operations at the end of the run
            "LWAVE": ".FALSE.",
            "LCHARG": ".FALSE.",
            "LPLANE": ".FALSE.",

        }

        if self.parameters.dft_calculator == "qe":
            atoms_Angstrom = ase.io.read(posfile, format="vasp")
            ase.io.write(os.path.join(base_path, self.parameters.element
                                      + ".pw.md.in"),
                         atoms_Angstrom,
                         "espresso-in", input_data=qe_input_data,
                         pseudopotentials= \
                             qe_pseudopotentials, kpts=kgrid)
        elif self.parameters.dft_calculator == "vasp":
            atoms_Angstrom = ase.io.read(posfile, format="vasp")
            ase.io.write(os.path.join(base_path, "POSCAR"), atoms_Angstrom,
                         "vasp")
            ase.io.write(os.path.join(base_path, "POSCAR_original"), atoms_Angstrom,
                         "vasp")
            VaspUtils.write_to_incar(base_path, vasp_input_data)
            VaspUtils.write_to_kpoints(base_path, kgrid)
            VaspUtils.write_to_potcar_copy(base_path,
                                           os.path.join(
                                               self.parameters.pseudopotential[
                                                   "path"],
                                               self.parameters.pseudopotential[
                                                   "name"]))


