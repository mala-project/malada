"""Runner that only creates SLURM input files."""

from .runner import Runner
from malada import SlurmParameters
import os
import glob


class SlurmCreatorRunner(Runner):
    """
    Create slurm input scripts without running them.

    Parameters
    ----------
    parameters : malada.utils.parametes.Parameters
        Parameters used to create this object.
    """

    def __init__(self, parameters):
        super(SlurmCreatorRunner, self).__init__(parameters)

    def run_folder(self, folder, calculation_type):
        """
        Run a folder (=create submit.slurm for it).

        Parameters
        ----------
        folder : string
            Folder to run in.

        calculation_type : string
            Type of calculation, currently supported are "dft" and "md".
        """
        # Write a submit.slurm file.
        slurm_params: SlurmParameters
        if calculation_type == "dft" or calculation_type == "dft+pp":
            calculator_type = self.parameters.dft_calculator
            slurm_params = self.parameters.dft_slurm
        elif calculation_type == "md":
            calculator_type = self.parameters.md_calculator
            slurm_params = self.parameters.md_slurm
        else:
            raise Exception("Unknown calculation type encountered.")

        job_name = os.path.basename(os.path.normpath(folder))
        submit_file = open(os.path.join(folder, "submit.slurm"), mode="w")
        submit_file.write("#!/bin/bash\n")
        submit_file.write("#SBATCH --nodes=" + str(slurm_params.nodes) + "\n")
        submit_file.write(
            "#SBATCH --ntasks-per-node="
            + str(slurm_params.tasks_per_node)
            + "\n"
        )
        submit_file.write("#SBATCH --job-name=" + job_name + "\n")
        if calculation_type != "md":
            submit_file.write("#SBATCH --output=" + job_name + ".out\n")
        submit_file.write(
            "#SBATCH --time=" + str(slurm_params.execution_time) + ":00:00\n"
        )
        submit_file.write(slurm_params.partition_string)
        submit_file.write("\n")
        submit_file.write(slurm_params.module_loading_string)
        submit_file.write("\n")
        if calculator_type == "qe":
            # TODO: Fix this.
            if calculation_type == "dft":
                # Get the filename.
                filelist = glob.glob(os.path.join(folder, "*.pw.scf.in"))
                if len(filelist) != 1:
                    print(filelist, folder)
                    raise Exception("Run folder with ambigous content.")
                filename = os.path.basename(filelist[0])
                submit_file.write(
                    slurm_params.mpi_runner
                    + " "
                    + slurm_params.get_mpirunner_process_params()
                    + " "
                    + str(slurm_params.nodes * slurm_params.tasks_per_node)
                    + " "
                    + slurm_params.scf_executable
                    + " -in "
                    + filename
                    + "\n"
                )
            elif calculation_type == "dft+pp":
                # Get the filenames.
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
                submit_file.write(
                    slurm_params.mpi_runner
                    + " "
                    + slurm_params.get_mpirunner_process_params()
                    + " "
                    + str(slurm_params.nodes * slurm_params.tasks_per_node)
                    + " "
                    + slurm_params.scf_executable
                    + " -in "
                    + os.path.basename(scf_file[0])
                    + " \n"
                )
                submit_file.write(
                    slurm_params.mpi_runner
                    + " "
                    + slurm_params.get_mpirunner_process_params()
                    + " "
                    + str(slurm_params.nodes * slurm_params.tasks_per_node)
                    + " "
                    + slurm_params.pp_executable
                    + " -in "
                    + os.path.basename(dens_file[0])
                    + " \n"
                )
                submit_file.write(
                    slurm_params.mpi_runner
                    + " "
                    + slurm_params.get_mpirunner_process_params()
                    + " "
                    + str(slurm_params.nodes * slurm_params.tasks_per_node)
                    + " "
                    + slurm_params.dos_executable
                    + " -in "
                    + os.path.basename(dos_file[0])
                    + " \n"
                )
                submit_file.write(
                    slurm_params.mpi_runner
                    + " "
                    + slurm_params.get_mpirunner_process_params()
                    + " "
                    + str(slurm_params.nodes * slurm_params.tasks_per_node)
                    + " "
                    + slurm_params.pp_executable
                    + " -in "
                    + os.path.basename(ldos_file[0])
                    + " \n"
                )
            elif calculation_type == "md":
                md_file = glob.glob(os.path.join(folder, "*.pw.md.in"))
                if len(md_file) != 1:
                    print(md_file, folder)
                    raise Exception("Run folder with ambigous content.")
                submit_file.write(
                    slurm_params.mpi_runner
                    + " "
                    + slurm_params.get_mpirunner_process_params()
                    + " "
                    + str(slurm_params.nodes * slurm_params.tasks_per_node)
                    + " "
                    + slurm_params.scf_executable
                    + " -in "
                    + os.path.basename(md_file[0])
                    + " \n"
                )
        elif calculator_type == "vasp":
            if calculation_type == "dft":
                submit_file.write("bash potcar_copy.sh\n")
                submit_file.write(
                    slurm_params.mpi_runner
                    + " "
                    + slurm_params.get_mpirunner_process_params()
                    + " "
                    + str(slurm_params.nodes * slurm_params.tasks_per_node)
                    + " "
                    + slurm_params.scf_executable
                    + " \n"
                )
            elif calculation_type == "md":
                submit_file.write("mkdir slurm-$SLURM_JOB_ID\n")
                submit_file.write("cp INCAR slurm-$SLURM_JOB_ID\n")
                submit_file.write("cp POSCAR slurm-$SLURM_JOB_ID\n")
                submit_file.write("cp KPOINTS slurm-$SLURM_JOB_ID\n")
                submit_file.write("cp potcar_copy.sh slurm-$SLURM_JOB_ID\n")
                submit_file.write("cd slurm-$SLURM_JOB_ID\n")
                submit_file.write("bash potcar_copy.sh\n")
                submit_file.write(
                    slurm_params.mpi_runner
                    + " "
                    + slurm_params.get_mpirunner_process_params()
                    + " "
                    + str(slurm_params.nodes * slurm_params.tasks_per_node)
                    + " "
                    + slurm_params.scf_executable
                    + " \n"
                )
        submit_file.write(slurm_params.cleanup_string)
        submit_file.write("\n")

        submit_file.close()
