"""Base class for all pipeline providers."""

from abc import ABC, abstractmethod
import xml.etree.ElementTree as ET
import glob
import os
import ase
import ase.io
import numpy as np


class Provider:
    """
    Abstract base class for defining providers subclasses.

    Apart from the constructor, each provider should have a provide() method.
    """

    def __init__(self, parameters):
        self.parameters = parameters

    @abstractmethod
    def provide(self, provider_path):
        """
        Use output from previous step to provide input for the next.

        Parameters
        ----------
        provider_path : string
            Path in which to operate in.
        """
        pass

    def _read_convergence(self, filename, ignore_atom_number=False):

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
        element_list = (
            [self.parameters.element]
            if isinstance(self.parameters.element, str)
            else self.parameters.element
        )
        element_string = ""
        for element in element_list:
            element_string += element

        element_string_saved = ""
        if isinstance(self.parameters.element, str):
            element_string_saved = dftparams.find("element").text
        else:
            for element_type in range(0, len(self.parameters.element)):
                element_string_saved += dftparams.find("element"+str(element_type)).text

        if (
                element_string != element_string_saved
            or dftparams.find("crystal_structure").text
            != self.parameters.crystal_structure
            or dftparams.find("dft_calculator").text
            != self.parameters.dft_calculator
            or float(dftparams.find("temperature").text)
            != self.parameters.temperature
            or number_of_atoms_check
        ):
            raise Exception("Incompatible convergence parameters provided.")

        cutoff_energy = int(filecontents.find("cutoff_energy").text)
        kpoints = filecontents.find("kpoints")
        kgrid = (
            int(kpoints.find("kx").text),
            int(kpoints.find("kx").text),
            int(kpoints.find("kx").text),
        )

        return cutoff_energy, kgrid

    def _qe_out_to_trajectory(self, out_folder, file_name):
        file_list = glob.glob(os.path.join(out_folder, "*.out"))
        ordered_file_list = sorted(file_list)

        # I know this breaks down if one of the out files is for any reason incorrect.
        i_actual = 0
        for posfile in ordered_file_list:
            current_atoms = ase.io.read(
                posfile, index=":", format="espresso-out"
            )
            for i in range(0, len(current_atoms)):
                if i_actual < self.parameters.maximum_number_of_timesteps:
                    atoms_to_write = self.enforce_pbc(current_atoms[i])
                    if i_actual == 0:
                        traj_writer = ase.io.trajectory.TrajectoryWriter(
                            file_name, mode="w"
                        )
                        traj_writer.write(atoms=atoms_to_write)
                        i_actual += 1
                    else:
                        traj_writer = ase.io.trajectory.TrajectoryWriter(
                            file_name, mode="a"
                        )
                        if i > 0:
                            traj_writer.write(atoms=atoms_to_write)
                            i_actual += 1
                else:
                    break

    def _qe_out_to_temperature(self, out_folder, file_name):
        file_list = glob.glob(os.path.join(out_folder, "*.out"))
        ordered_file_list = sorted(file_list)
        temps = []

        # I know this breaks down if one of the out files is for any reason incorrect.
        i = 0
        for file_to_open in ordered_file_list:
            posfile = open(file_to_open)
            for line in posfile.readlines():
                if i < self.parameters.maximum_number_of_timesteps:
                    if (
                        "temperature" in line
                        and "=" in line
                        and "Starting" not in line
                    ):
                        temp = float((line.split("=")[1]).split("K")[0])
                        temps.append(temp)
                        i += 1
                else:
                    break
        np.save(file_name, temps)

    def _qe_out_to_timing(self, out_folder, file_name):
        file_list = glob.glob(os.path.join(out_folder, "*.out"))
        ordered_file_list = sorted(file_list)

        # I know this breaks down if one of the out files is for any reason incorrect.
        i = 0
        times = []
        for file_to_open in ordered_file_list:
            posfile = open(file_to_open)
            time_step_found = False
            last_time = 0
            for line in posfile.readlines():
                if i < self.parameters.maximum_number_of_timesteps:
                    if "Ekin + Etot (const)" in line:
                        time_step_found = True
                    if time_step_found and "cpu time" in line:
                        current_time = float(line.split()[-2])
                        times.append(current_time - last_time)
                        last_time = current_time
                        time_step_found = False
                        i += 1
                else:
                    break
        np.save(file_name, times)

    def _vasp_out_to_temperature(self, out_folder, file_name):
        file_list = glob.glob(os.path.join(out_folder, "slurm-*/"))
        ordered_file_list = sorted(file_list)
        temps = []
        # I know this breaks down if one of the out files is for any reason incorrect.
        i = 0
        for file_to_open in ordered_file_list:
            posfile = open(os.path.join(file_to_open, "REPORT"))
            for line in posfile.readlines():
                if i < self.parameters.maximum_number_of_timesteps:
                    if "tmprt" in line:
                        temp = float(line.split()[2])
                        temps.append(temp)
                        i += 1
                else:
                    break
        np.save(file_name, temps)

    def _vasp_out_to_timing(self, out_folder, file_name):
        file_list = glob.glob(os.path.join(out_folder, "slurm-*/"))
        ordered_file_list = sorted(file_list)
        temps = []
        # I know this breaks down if one of the out files is for any reason incorrect.
        i = 0
        times = []
        for file_to_open in ordered_file_list:
            posfile = open(os.path.join(file_to_open, "OUTCAR"))
            for line in posfile.readlines():
                if i < self.parameters.maximum_number_of_timesteps:
                    if "LOOP+" in line:
                        current_time = float(line.split()[-1])
                        times.append(current_time)
                        i += 1
                else:
                    break
        np.save(file_name, times)

    def _vasp_out_to_trajectory(self, out_folder, file_name):
        file_list = glob.glob(os.path.join(out_folder, "slurm-*/"))
        ordered_file_list = sorted(file_list)
        # I know this breaks down if one of the out files is for any reason incorrect.
        i_actual = 0
        for file_to_open in ordered_file_list:
            current_atoms = ase.io.read(
                os.path.join(file_to_open, "OUTCAR"),
                index=":",
                format="vasp-out",
            )
            for i in range(0, len(current_atoms)):
                if i_actual < self.parameters.maximum_number_of_timesteps:
                    atoms_to_write = self.enforce_pbc(current_atoms[i])
                    if i_actual == 0:
                        traj_writer = ase.io.trajectory.TrajectoryWriter(
                            file_name, mode="w"
                        )
                        traj_writer.write(atoms=atoms_to_write)
                        i_actual += 1
                    else:
                        traj_writer = ase.io.trajectory.TrajectoryWriter(
                            file_name, mode="a"
                        )
                        if i > 0:
                            traj_writer.write(atoms=atoms_to_write)
                            i_actual += 1
                else:
                    break

    def _get_number_of_bands(self):
        number_of_bands = int(
            self.parameters.number_of_atoms
            * self.parameters.pseudopotential["valence_electrons"]
            * (1.0 + self.parameters.number_of_bands_factor)
        )
        return number_of_bands

    @staticmethod
    def enforce_pbc(atoms):
        """
        Explictly enforeces the PBC on an ASE atoms object.

        QE (and potentially other codes?) do that internally. Meaning that the
        raw positions of atoms (in Angstrom) can lie outside of the unit cell.
        When setting up the DFT calculation, these atoms get shifted into
        the unit cell. Since we directly use these raw positions for the
        descriptor calculation, we need to enforce that in the ASE atoms
        objects, the atoms are explicitly in the unit cell.

        Parameters
        ----------
        atoms : ase.atoms
            The ASE atoms object for which the PBC need to be enforced.

        Returns
        -------
        new_atoms : ase.atoms
            The ASE atoms object for which the PBC have been enforced.
        """
        new_atoms = atoms.copy()
        new_atoms.pbc = True
        new_atoms.set_scaled_positions(new_atoms.get_scaled_positions())

        # This might be unecessary, but I think it is nice to have some sort of
        # metric here.
        rescaled_atoms = 0
        for i in range(0, len(atoms)):
            if False in (
                np.isclose(
                    new_atoms[i].position, atoms[i].position, atol=0.001
                )
            ):
                rescaled_atoms += 1
        return new_atoms
