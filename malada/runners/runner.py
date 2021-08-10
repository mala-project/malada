from abc import ABC, abstractmethod


class Runner:
    """
    Runs a simulation.
    """

    def __init__(self, parameters):
        self.parameters = parameters
        pass

    @abstractmethod
    def run_folder(self, folder, calculator_type, qe_input_type=".pw.scf.in"):
        pass
