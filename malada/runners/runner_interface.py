"""Interface to automate creation of Runners."""
from .bashrunner import BashRunner
from .slurm_creator import SlurmCreatorRunner


def RunnerInterface(parameters):
    """
    Get the correct runner for the parameters provided.

    Parameters
    ----------
    parameters : malada.utils.parametes.Parameters
        Parameters used to create this object.

    Returns
    -------
    runner : malada.Runner

    """
    runner = None
    if parameters.run_system == "bash":
        runner = BashRunner(parameters)
    if parameters.run_system == "slurm_creator":
        runner = SlurmCreatorRunner(parameters)

    if runner is not None:
        return runner
    else:
        raise Exception("Unknown runner type.")

