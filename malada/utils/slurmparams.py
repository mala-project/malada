
class SlurmParameters:
    """
    Class holding slurm parameters.
    """

    def __init__(self):
        self.scf_executable = "pw.x"
        self.module_loading_string = ""
        self.execution_time = ""
        self.partition_string = ""
        self.mpi_runner = "mpirun"

        # TODO: Maybe some sanity checks here?
        self.tasks_per_node = 0
        self.nodes = 0
