class VaspUtils:
    """Collection of functions to write to VASP files other then POSCAR."""

    @staticmethod
    def write_to_incar(folder, filename, vasp_array):
        file_handle = open(folder + filename, "w")
        for entry in vasp_array:
            file_handle.write(entry + "=" + str(vasp_array[entry]) + "\n")

    @staticmethod
    def write_to_kpoints(folder, filename, kgrid):
        file_handle = open(folder + filename, "w")
        file_handle.write("Automatic mesh\n")
        file_handle.write("0\n")
        file_handle.write("Gamma\n")
        file_handle.write(
            str(kgrid[0]) + " " + str(kgrid[1]) + " " + str(kgrid[2]) + "\n")
        file_handle.write("0. 0. 0.\n")
