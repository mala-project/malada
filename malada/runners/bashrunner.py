"""Bash based runner. This will be replaced by an ASE based system."""

import glob
import os
import ase.io
import subprocess

from .runner import Runner


class BashRunner(Runner):
    """
    Runs a simulation.

    This is a really ugly, system calls based runner, that will soon be changed
    out for a better, ASE based runner.
    """

    def __init__(self, parameters):
        super(BashRunner, self).__init__(parameters)

    def run_folder(self, folder, calculation_type):
        """
        Run a folder on the command line.

        Parameters
        ----------
        folder : string
            Folder to run in.

        calculation_type : string
            Type of calculation, currently supported are "dft" and "md".
        """
        if calculation_type == "dft" or calculation_type == "dft+pp":
            calculator_type = self.parameters.dft_calculator
        elif calculation_type == "md":
            calculator_type = self.parameters.md_calculator
        else:
            raise Exception("Unknown calculation type encountered.")
        job_name = os.path.basename(os.path.normpath(folder))
        if calculator_type == "qe":
            if calculation_type == "dft":
                filelist = glob.glob(os.path.join(folder, "*.pw.scf.in"))
                if len(filelist) != 1:
                    print(filelist, folder)
                    raise Exception("Run folder with ambigous content.")
                filename = os.path.basename(filelist[0])
                run_process = subprocess.Popen(
                    "pw.x -in " + filename + " > " + job_name + ".out",
                    cwd=folder,
                    shell=True,
                )
                run_process.wait()
            if calculation_type == "md":
                filelist = glob.glob(os.path.join(folder, "*.pw.md.in"))
                if len(filelist) != 1:
                    print(filelist, folder)
                    raise Exception("Run folder with ambigous content.")
                filename = os.path.basename(filelist[0])
                run_process = subprocess.Popen(
                    "pw.x -in " + filename + " > " + job_name + ".out",
                    cwd=folder,
                    shell=True,
                )
                run_process.wait()
            elif calculation_type == "dft+pp":
                scf_file = glob.glob(os.path.join(folder, "*.pw.scf.in"))
                ldos_file = glob.glob(os.path.join(folder, "*.pp.ldos.in"))
                dens_file = glob.glob(os.path.join(folder, "*.pp.dens.in"))
                dos_file = glob.glob(os.path.join(folder, "*.dos.in"))
                if (
                    len(scf_file) != 1
                    or len(ldos_file) != 1
                    or len(dens_file) != 1
                    or len(dos_file) != 1
                ):
                    print(scf_file, ldos_file, dens_file, dos_file, folder)
                    raise Exception("Run folder with ambigous content.")
                run_process = subprocess.Popen(
                    "pw.x -in "
                    + os.path.basename(scf_file[0])
                    + " > "
                    + job_name
                    + ".out;"
                    + "pp.x -in "
                    + os.path.basename(dens_file[0])
                    + ";"
                    + "dos.x -in "
                    + os.path.basename(dos_file[0])
                    + ";"
                    + "pp.x -in "
                    + os.path.basename(ldos_file[0])
                    + ";",
                    cwd=folder,
                    shell=True,
                )
                run_process.wait()

        elif calculator_type == "vasp":
            raise Exception("VASP currently not implemented.")
        else:
            raise Exception("Calculator type unknown.")
