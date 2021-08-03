from abc import ABC, abstractmethod


class Runner:
    """
    Runs a simulation.
    """

    # TODO: Add some parameters to characterize the running.
    def __init__(self):
        pass

    @abstractmethod
    def run_folder(self, folder, calculator_type, qe_input_type=".pw.scf.in"):
        pass
