import os

class VaspUtils:
    """Collection of functions to write to VASP files other then POSCAR."""

    @staticmethod
    def write_to_incar(folder, vasp_array):
        file_handle = open(os.path.join(folder, "INCAR"), "w")#
        for entry in vasp_array:
            file_handle.write(entry + "=" + str(vasp_array[entry]) + "\n")
        file_handle.close()

    @staticmethod
    def write_to_kpoints(folder, kgrid):
        file_handle = open(os.path.join(folder, "KPOINTS"), "w")#
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
        file_handle = open(os.path.join(folder, "potcar_copy.sh"), "w")#
        file_handle.write("#!/bin/bash\n")
        file_handle.write("cp "+ pspstring+" POTCAR\n")
        file_handle.close()
