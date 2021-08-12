"""Base class for all runners."""
from abc import ABC, abstractmethod


class Runner:
    """
    Abstract base class to run a simulation.

    Parameters
    ----------
    parameters : malada.utils.parametes.Parameters
        Parameters used to create this object.
    """

    def __init__(self, parameters):
        self.parameters = parameters
        pass

    @abstractmethod
    def run_folder(self, folder, calculation_type, qe_input_type=".pw.scf.in"):
        """
        Run a simulation in a folder.

        Parameters
        ----------
        folder : string
            Folder to run in.

        calculation_type : string
            Type of calculation, currently supported are "dft" and "md".

        qe_input_type : string
            Details the type of DFT calculation to be performed by QE.
        """
        pass
