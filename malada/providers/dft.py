"""Provider for DFT calculations to get energies and LDOS."""

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


class DFTProvider(Provider):
    """
    Performs a DFT calculation and provides an DFT output plus LDOS.

    Parameters
    ----------
    parameters : malada.utils.parametes.Parameters
        Parameters used to create this object.

    external_calculation_folders:
        Path to folders containing already (half) finished calculations.
        If not None, MALADA will try to assess which calculations are still
        missing and then perform those.
    """

    def __init__(self, parameters, external_calculation_folders=None):
        super(DFTProvider, self).__init__(parameters)
        self.external_calculation_folders = external_calculation_folders
        self.calculation_folders = None

    def provide(
        self,
        provider_path,
        dft_convergence_file,
        ldos_convergence_file,
        possible_snapshots_file,
        do_postprocessing=True,
        numbering_starts_at=0,
        parsing_starts_at=0,
        ignore_atom_number=False,
    ):
        """
        Provide a set of DFT calculations on predefined snapshots.

        This includes DFT energies and snapshots.

        Parameters
        ----------
        provider_path : string
            Path in which to operate in.

        dft_convergence_file : string
            Path to xml file containing the DFT convergence parameter.

        ldos_convergence_file : string
            Path to xml file containing the LDOS convergence parameter.
            This means a different set of DFT parameters needed for
            LDOS calculation.

        possible_snapshots_file : string
            Path to a file containing an ASE trajectory containing atomic
            snapshots for DFT/LDOS calculation.

        numbering_starts_at : int
            Number from which to start the numbering of the snapshots.
            Zero by default. Useful for incremental data generation.
            THIS WILL ONLY AFFECT THE NUMBERING! NOT THE SNAPSHOTS USED.
            E.g. if this value is set to 10, the 1st snapshot of the
            selected snapshots will be set as number 10. If you want to add
            more snapshots to an existing data set, use parsing_starts_at.

        parsing_starts_at : int
            Overwrites numbering_starts_at. Number from which to start
            processing snapshots.

        ignore_atom_number : bool
            If True, a convergence file for a different number of atoms
            can be loaded.
        """
        if parsing_starts_at > 0:
            numbering_starts_at = parsing_starts_at
        self.calculation_folders = provider_path
        if self.external_calculation_folders is None:
            # Here we have to perform the actucal calculation.
            if self.parameters.dft_calculator != "qe":
                raise Exception(
                    "Currently only QE is supported for " "DFT calculations."
                )
            dft_runner = malada.RunnerInterface(self.parameters)
            all_valid_snapshots = ase.io.Trajectory(possible_snapshots_file)
            for i in range(0, self.parameters.number_of_snapshots):
                snapshot_number = i + numbering_starts_at
                snapshot_path = os.path.join(
                    provider_path, "snapshot" + str(snapshot_number)
                )
                self.__create_dft_run(
                    dft_convergence_file,
                    ldos_convergence_file,
                    all_valid_snapshots[
                        snapshot_number
                        - numbering_starts_at
                        + parsing_starts_at
                    ],
                    snapshot_path,
                    "snapshot" + str(snapshot_number),
                    do_postprocessing,
                    ignore_atom_number,
                )
            for i in range(0, self.parameters.number_of_snapshots):
                # Run the individul files.
                snapshot_number = i + numbering_starts_at
                snapshot_path = os.path.join(
                    provider_path, "snapshot" + str(snapshot_number)
                )
                print("Running DFT in", snapshot_path)
                if do_postprocessing:
                    dft_runner.run_folder(snapshot_path, "dft+pp")
                else:
                    dft_runner.run_folder(snapshot_path, "dft")
        else:
            # Here, we have to do a consistency check.
            pass

    def __create_dft_run(
        self,
        dft_convergence_file,
        ldos_convergence_file,
        atoms,
        snapshot_path,
        snapshot_name,
        do_postprocessing,
        ignore_atom_number,
    ):
        # Get cluster info
        # TODO: Use DFT kgrid for scf and LDOS kgrid for NSCF calculations.
        cutoff, kgrid = self._read_convergence(
            dft_convergence_file, ignore_atom_number
        )
        ldos_params, kgrid = self.__read_ldos_convergence(
            ldos_convergence_file, ignore_atom_number
        )
        # Create folder
        if not os.path.exists(snapshot_path):
            os.makedirs(snapshot_path)

        qe_pseudopotentials = {
            self.parameters.element: self.parameters.pseudopotential["name"]
        }
        nbands = self._get_number_of_bands()
        outdir = "temp"
        if (
            self.parameters.dft_scf_accuracy_per_atom_Ry >= 1e-6
            and self.parameters.dft_assume_two_dimensional
        ):
            print(
                "Large DFT accuracy threshold detected. When running "
                "two-dimensional calculations, which include areas of low "
                "electronic density, smaller accuracy thresholds are "
                "recommended. Consider setting dft_scf_accuracy_per_atom_Ry "
                "to, e.g., 1-e9."
            )
        qe_input_data = {
            "occupations": "smearing",
            "calculation": "scf",
            "restart_mode": "from_scratch",
            "verbosity": "high",
            "prefix": self.parameters.element,
            "pseudo_dir": self.parameters.pseudopotential["path"],
            "outdir": outdir,
            "ibrav": 0,
            "smearing": "fermi-dirac",
            "degauss": round(
                kelvin_to_rydberg(self.parameters.temperature), 7
            ),
            "ecutrho": cutoff * 4,
            "ecutwfc": cutoff,
            "tstress": True,
            "tprnfor": True,
            "nbnd": nbands,
            "mixing_mode": "plain",
            "mixing_beta": self.parameters.dft_mixing_beta,
            "conv_thr": self.parameters.dft_scf_accuracy_per_atom_Ry
            * self.parameters.number_of_atoms,
            # "verbosity" : "high", # This is maybe a bit high
            "nosym": True,
            "noinv": True,
        }
        if self.parameters.dft_use_inversion_symmetry is True:
            qe_input_data["noinv"] = False
        if self.parameters.dft_calculate_stress is False:
            qe_input_data["tstress"] = False
        if self.parameters.dft_calculate_force is False:
            qe_input_data["tprnfor"] = False
        if self.parameters.dft_assume_two_dimensional:
            qe_input_data["assume_isolated"] = "2D"

        id_string = self.parameters.element + "_" + snapshot_name
        ase.io.write(
            os.path.join(snapshot_path, id_string + ".pw.scf.in"),
            atoms,
            "espresso-in",
            input_data=qe_input_data,
            pseudopotentials=qe_pseudopotentials,
            kpts=kgrid,
        )

        emin = ldos_params["ldos_offset_eV"]
        deltae = ldos_params["ldos_spacing_eV"]
        emax = (
            ldos_params["ldos_length"] * ldos_params["ldos_spacing_eV"]
            + ldos_params["ldos_offset_eV"]
        )
        smearing_factor = ldos_params["smearing_factor"]

        if do_postprocessing:
            # DOS file
            dos_file = open(
                os.path.join(snapshot_path, id_string + ".dos.in"), mode="w"
            )
            dos_file.write("&dos\n")
            dos_file.write(" outdir='" + outdir + "',\n")
            dos_file.write(" prefix='" + self.parameters.element + "',\n")
            dos_file.write(" Emin=" + str(emin) + ",\n")
            dos_file.write(" Emax=" + str(emax) + ",\n")
            dos_file.write(" DeltaE=" + str(deltae) + ",\n")
            dos_file.write(
                " degauss=" + str(deltae * smearing_factor / Rydberg) + ",\n"
            )
            dos_file.write(" fildos='" + id_string + ".dos'\n")
            dos_file.write("/\n")

            # LDOS file
            ldos_file = open(
                os.path.join(snapshot_path, id_string + ".pp.ldos.in"),
                mode="w",
            )
            ldos_file.write("&inputpp\n")
            ldos_file.write(" outdir='" + outdir + "',\n")
            ldos_file.write(" prefix='" + self.parameters.element + "',\n")
            ldos_file.write(" plot_num=3,\n")
            ldos_file.write(" emin=" + str(emin) + ",\n")
            ldos_file.write(" emax=" + str(emax) + ",\n")
            ldos_file.write(" delta_e=" + str(deltae) + ",\n")
            ldos_file.write(
                " degauss_ldos=" + str(deltae * smearing_factor) + ",\n"
            )
            ldos_file.write(" use_gauss_ldos=.true.\n")
            ldos_file.write("/\n")
            ldos_file.write("&plot\n")
            ldos_file.write(" iflag=3,\n")
            ldos_file.write(" output_format=6,\n")
            ldos_file.write(" fileout='" + id_string + "_ldos.cube',\n")
            ldos_file.write("/\n")

            # Density file
            dens_file = open(
                os.path.join(snapshot_path, id_string + ".pp.dens.in"),
                mode="w",
            )
            dens_file.write("&inputpp\n")
            dens_file.write(" outdir='" + outdir + "',\n")
            dens_file.write(" prefix='" + self.parameters.element + "',\n")
            dens_file.write(" plot_num=0,\n")
            dens_file.write("/\n")
            dens_file.write("&plot\n")
            dens_file.write(" iflag=3,\n")
            dens_file.write(" output_format=6,\n")
            dens_file.write(" fileout='" + id_string + "_dens.cube',\n")
            dens_file.write("/\n")

            dos_file.close()
            ldos_file.close()
            dens_file.close()

    def __read_ldos_convergence(self, filename, ignore_atom_number):

        # Parse the XML file and first check for consistency.
        filecontents = ET.parse(filename).getroot()
        dftparams = filecontents.find("calculationparameters")
        number_of_atoms_check = (
            (
                int(dftparams.find("number_of_atoms").text)
                != self.parameters.number_of_atoms
            )
            if (ignore_atom_number is False)
            else False
        )
        if (
            dftparams.find("element").text != self.parameters.element
            or dftparams.find("crystal_structure").text
            != self.parameters.crystal_structure
            or dftparams.find("dft_calculator").text
            != self.parameters.dft_calculator
            or float(dftparams.find("temperature").text)
            != self.parameters.temperature
            or number_of_atoms_check
        ):
            raise Exception("Incompatible convergence parameters provided.")
        ldos_params = {
            "ldos_offset_eV": float(filecontents.find("ldos_offset_eV").text),
            "ldos_spacing_eV": float(
                filecontents.find("ldos_spacing_eV").text
            ),
            "ldos_length": int(filecontents.find("ldos_length").text),
            "smearing_factor": float(
                filecontents.find("smearing_factor").text
            ),
        }
        kpoints = filecontents.find("kpoints")
        kgrid = (
            int(kpoints.find("kx").text),
            int(kpoints.find("ky").text),
            int(kpoints.find("kz").text),
        )

        return ldos_params, kgrid
