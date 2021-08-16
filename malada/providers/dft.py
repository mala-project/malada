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

    def provide(self, provider_path, dft_convergence_file,
                ldos_convergence_file, possible_snapshots_file):
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
        """
        self.calculation_folders = provider_path
        if self.external_calculation_folders is None:
            # Here we have to perform the actucal calculation.
            if self.parameters.dft_calculator != "qe":
                raise Exception("Currently only QE is supported for "
                                "DFT calculations.")
            dft_runner = malada.RunnerInterface(self.parameters)
            all_valid_snapshots = ase.io.Trajectory(possible_snapshots_file)
            for i in range(0, self.parameters.number_of_snapshots):
                snapshot_path = os.path.join(provider_path,"snapshot"+str(i))
                self.__create_dft_run(dft_convergence_file,
                                      ldos_convergence_file,
                                      all_valid_snapshots[i],
                                      snapshot_path,
                                      "snapshot"+str(i))
            for i in range(0, self.parameters.number_of_snapshots):
                # Run the individul files.
                snapshot_path = os.path.join(provider_path,"snapshot"+str(i))
                print("Running DFT in", snapshot_path)
                dft_runner.run_folder(snapshot_path, "dft+pp")
        else:
            # Here, we have to do a consistency check.
            pass

    def __create_dft_run(self, dft_convergence_file, ldos_convergence_file,
                        atoms, snapshot_path, snapshot_name):
        # Get cluster info
        # TODO: Use DFT kgrid for scf and LDOS kgrid for NSCF calculations.
        cutoff, kgrid = self._read_convergence(dft_convergence_file)
        ldos_params, kgrid = self.__read_ldos_convergence(ldos_convergence_file)
        # Create folder
        if not os.path.exists(snapshot_path):
            os.makedirs(snapshot_path)

        qe_pseudopotentials = {self.parameters.element:
                                   self.parameters.pseudopotential["name"]}
        nbands = self._get_number_of_bands()
        outdir = "temp"
        qe_input_data = {
            "occupations": 'smearing',
            "calculation": 'scf',
            "restart_mode": 'from_scratch',
            "verbosity": 'high',
            "prefix": self.parameters.element,
            "pseudo_dir": self.parameters.pseudopotential["path"],
            "outdir": outdir,
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
            "conv_thr": self.parameters.dft_scf_accuracy_per_atom_Ry * self.parameters.number_of_atoms,
            # "verbosity" : "high", # This is maybe a bit high
            "nosym": True,
            "noinv": True

        }
        id_string = self.parameters.element+"_"+snapshot_name
        ase.io.write(os.path.join(snapshot_path,id_string+".pw.scf.in"),
            atoms,
            "espresso-in", input_data=qe_input_data, pseudopotentials= \
                qe_pseudopotentials, kpts=kgrid)

        emin = ldos_params["ldos_offset_eV"]
        deltae = ldos_params["ldos_spacing_eV"]
        emax = ldos_params["ldos_length"]*ldos_params["ldos_spacing_eV"]+ldos_params["ldos_offset_eV"]
        smearing_factor = ldos_params["smearing_factor"]

        # DOS file
        dos_file = open(
            os.path.join(snapshot_path,id_string+".dos.in"),
            mode='w')
        dos_file.write("&dos\n")
        dos_file.write(" outdir='" + outdir+"',\n")
        dos_file.write(" prefix='" + self.parameters.element + "',\n")
        dos_file.write(" Emin=" + str(emin) + ",\n")
        dos_file.write(" Emax=" + str(emax) + ",\n")
        dos_file.write(" DeltaE=" + str(deltae) + ",\n")
        dos_file.write(
            " degauss=" + str(deltae * smearing_factor / Rydberg) + ",\n")
        dos_file.write(" fildos='" + id_string + ".dos'\n")
        dos_file.write("/\n")

        # LDOS file
        ldos_file = open(os.path.join(snapshot_path,id_string+".pp.ldos.in"),
            mode='w')
        ldos_file.write("&inputpp\n")
        ldos_file.write(" outdir='" + outdir+"',\n")
        ldos_file.write(" prefix='" + self.parameters.element + "',\n")
        ldos_file.write(" plot_num=3,\n")
        ldos_file.write(" emin=" + str(emin) + ",\n")
        ldos_file.write(" emax=" + str(emax) + ",\n")
        ldos_file.write(" delta_e=" + str(deltae) + ",\n")
        ldos_file.write(
            " degauss_ldos=" + str(deltae * smearing_factor) + ",\n")
        ldos_file.write("/\n")
        ldos_file.write("&plot\n")
        ldos_file.write(" iflag=3,\n")
        ldos_file.write(" output_format=6,\n")
        ldos_file.write(
            " fileout='" + id_string + "_ldos.cube',\n")
        ldos_file.write("/\n")

        # Density file
        dens_file = open(os.path.join(snapshot_path,id_string+".pp.dens.in"),
            mode='w')
        dens_file.write("&inputpp\n")
        dens_file.write(" outdir='" + outdir+"',\n")
        dens_file.write(" prefix='" + self.parameters.element + "',\n")
        dens_file.write(" plot_num=0,\n")
        dens_file.write("/\n")
        dens_file.write("&plot\n")
        dens_file.write(" iflag=3,\n")
        dens_file.write(" output_format=6,\n")
        dens_file.write(
            " fileout='" + id_string + "_dens.cube',\n")
        dens_file.write("/\n")

        dos_file.close()
        ldos_file.close()
        dens_file.close()

    def __read_ldos_convergence(self, filename):

        # Parse the XML file and first check for consistency.
        filecontents = ET.parse(filename).getroot()
        dftparams = filecontents.find("calculationparameters")
        if dftparams.find("element").text != self.parameters.element or \
           dftparams.find("crystal_structure").text != self.parameters.crystal_structure or \
           dftparams.find("dft_calculator").text != self.parameters.dft_calculator or \
           float(dftparams.find("temperature").text) != self.parameters.temperature or \
           int(dftparams.find("number_of_atoms").text) != self.parameters.number_of_atoms:
            raise Exception("Incompatible convergence parameters provided.")
        ldos_params = {"ldos_offset_eV": float(filecontents.find("ldos_offset_eV").text),
                       "ldos_spacing_eV": float(filecontents.find("ldos_spacing_eV").text),
                       "ldos_length": int(filecontents.find("ldos_length").text),
                       "smearing_factor": float(filecontents.find("smearing_factor").text),}
        kpoints = filecontents.find("kpoints")
        kgrid = (int(kpoints.find("kx").text),int(kpoints.find("kx").text),
                 int(kpoints.find("kx").text))

        return ldos_params, kgrid
