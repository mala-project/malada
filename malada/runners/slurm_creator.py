"""Runner that only creates SLURM input files."""
from .runner import Runner
from malada import SlurmParameters
import os


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

    def run_folder(self, folder, calculation_type,
                   qe_input_type="*.pw.scf.in"):
        """
        Run a folder (=create submit.slurm for it).

        Parameters
        ----------
        folder : string
            Folder to run in.

        calculation_type : string
            Type of calculation, currently supported are "dft" and "md".

        qe_input_type : string
            Details the type of DFT calculation to be performed by QE.

        """
        # Write a submit.slurm file.
        slurm_params : SlurmParameters
        if calculation_type == "dft":
            calculator_type = self.parameters.dft_calculator
            slurm_params = self.parameters.dft_slurm
        if calculation_type == "md":
            calculator_type = self.parameters.md_calculator
            slurm_params = self.parameters.md_slurm

        job_name = os.path.basename(os.path.normpath(folder))
        submit_file = open(folder+"submit.slurm", mode='w')
        submit_file.write("#!/bin/bash\n")
        submit_file.write("#SBATCH --nodes="+str(slurm_params.nodes)+"\n")
        submit_file.write("#SBATCH --ntasks-per-node="+str(slurm_params.tasks_per_node)+"\n")
        submit_file.write("#SBATCH --job-name="+job_name+"\n")
        if calculation_type != "md":
            submit_file.write("#SBATCH --output="+job_name+".out\n")
        submit_file.write("#SBATCH --time="+str(slurm_params.execution_time)+":00:00\n")
        submit_file.write(slurm_params.partition_string)
        submit_file.write("\n")
        submit_file.write(slurm_params.module_loading_string)
        submit_file.write("\n")
        if calculator_type == "qe":
            # TODO: Fix this.
            submit_file.write(slurm_params.mpi_runner+" -np "+
                              str(slurm_params.nodes*slurm_params.tasks_per_node)+
                              " pw.x -in "+self.parameters.element+".pw.scf.in \n")
        elif calculator_type == "vasp":
            if calculation_type == "dft":
                submit_file.write("bash potcar_copy.sh\n")
                submit_file.write(slurm_params.mpi_runner+" -np "+
                                  str(slurm_params.nodes*slurm_params.tasks_per_node)+" "+
                                  slurm_params.scf_executable+" \n")
            elif calculation_type == "md":
                submit_file.write("mkdir slurm-$SLURM_JOB_ID\n")
                submit_file.write("cp INCAR slurm-$SLURM_JOB_ID\n")
                submit_file.write("cp POSCAR slurm-$SLURM_JOB_ID\n")
                submit_file.write("cp KPOINTS slurm-$SLURM_JOB_ID\n")
                submit_file.write("cp potcar_copy.sh slurm-$SLURM_JOB_ID\n")
                submit_file.write("cd slurm-$SLURM_JOB_ID\n")
                submit_file.write("bash potcar_copy.sh\n")
                submit_file.write(slurm_params.mpi_runner+" -np "+
                                  str(slurm_params.nodes*slurm_params.tasks_per_node)+" "+
                                  slurm_params.scf_executable+" \n")

        submit_file.close()
