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

    Each provider should have at least three methods, the constructor
    included:
        - constructor
        - from_file() from providing its results from a user-defined file
        - provide() function to perform its calculation and provide the
        results.
    """

    def __init__(self, parameters):
        self.parameters = parameters

    @abstractmethod
    def provide(self, provider_path):
        pass

    def _read_convergence(self, filename):

        # Parse the XML file and first check for consistency.
        filecontents = ET.parse(filename).getroot()
        dftparams = filecontents.find("calculationparameters")
        if dftparams.find("element").text != self.parameters.element or \
           dftparams.find("crystal_structure").text != self.parameters.crystal_structure or \
           dftparams.find("dft_calculator").text != self.parameters.dft_calculator or \
           float(dftparams.find("temperature").text) != self.parameters.temperature or \
           int(dftparams.find("number_of_atoms").text) != self.parameters.number_of_atoms:
            raise Exception("Incompatible convergence parameters provided.")

        cutoff_energy = int(filecontents.find("cutoff_energy").text)
        kpoints = filecontents.find("kpoints")
        kgrid = (int(kpoints.find("kx").text),int(kpoints.find("kx").text),
                 int(kpoints.find("kx").text))

        return cutoff_energy, kgrid

    def _qe_out_to_trajectory(self, out_folder, file_name):
        file_list = glob.glob(os.path.join(out_folder, "*.out"))
        ordered_file_list = sorted(file_list)

        # I know this breaks down if one of the out files is for any reason incorrect.
        i_actual = 0
        for posfile in ordered_file_list:
            current_atoms = ase.io.read(posfile,index = ':', format="espresso-out")
            for i in range(0, len(current_atoms)):
                if i_actual < self.parameters.maximum_number_of_timesteps:
                    if i_actual == 0:
                        traj_writer = ase.io.trajectory.TrajectoryWriter(file_name, mode='w')
                        traj_writer.write(atoms=current_atoms[i])
                        i_actual += 1
                    else:
                        traj_writer = ase.io.trajectory.TrajectoryWriter(file_name, mode='a')
                        if i > 0:
                            traj_writer.write(atoms=current_atoms[i])
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
                    if "temperature" in line and "=" in line and "Starting" not in line:
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
                        times.append(current_time-last_time)
                        last_time = current_time
                        time_step_found = False
                        i += 1
                else:
                    break
        np.save(file_name)
