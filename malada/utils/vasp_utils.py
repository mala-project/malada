"""Utilities for setting up VASP calculations."""
import os


class VaspUtils:
    """Collection of functions to write to VASP files other then POSCAR."""

    @staticmethod
    def write_to_incar(folder, vasp_array):
        """
        Write a dictionary of VASP input quantities to an INCAR file.

        Currently, there is NO sanity checking done here. Meaning that
        wrong INCAR files can be reated very easily.

        Parameters
        ----------
        folder : string
            Path to where to write the INCAR file.

        vasp_array : dict
            Dictionary of VASP quantities.
        """
        file_handle = open(os.path.join(folder, "INCAR"), "w")
        for entry in vasp_array:
            file_handle.write(entry + "=" + str(vasp_array[entry]) + "\n")
        file_handle.close()

    @staticmethod
    def write_to_kpoints(folder, kgrid):
        """
        Write a k-grid to VASP KPOINTS file.

        Parameters
        ----------
        folder : string
            Path to where to write the INCAR file.

        kgrid : tuple
            k-grid to write in form of (kx, ky, kz).
        """
        file_handle = open(os.path.join(folder, "KPOINTS"), "w")
        file_handle.write("Automatic mesh\n")
        file_handle.write("0\n")
        file_handle.write("Gamma\n")
        file_handle.write(
            str(kgrid[0]) + " " + str(kgrid[1]) + " " + str(kgrid[2]) + "\n")
        file_handle.write("0. 0. 0.\n")
        file_handle.close()

    # TODO: Find a more elegant way to do this...
    @staticmethod
    def write_to_potcar_copy(folder, pspstring):
        """
        Create a file that will copy the POTCAR file to the run folder.

        VASP always looks for the POTCAR (PSP) in the current folder.
        Therefore, before running, we have to copy it.

        Parameters
        ----------
        folder : string
            Path tto where to create this script.

        pspstring : string
            Path to POTCAR file.
        """
        file_handle = open(os.path.join(folder, "potcar_copy.sh"), "w")
        file_handle.write("#!/bin/bash\n")
        file_handle.write("cp "+ pspstring+" POTCAR\n")
        file_handle.close()
