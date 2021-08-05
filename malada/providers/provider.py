from abc import ABC, abstractmethod
import xml.etree.ElementTree as ET

class Provider:
    """
    Abstract base class for defining providers subclasses.

    Each provider should have at least three methods, the constructor
    included:
        - constructor
        - from_file() from providing its results from a user-defined file
        - provide() function to perform its calculation and provide the
        results.
    """

    def __init__(self, parameters):
        self.parameters = parameters

    @abstractmethod
    def provide(self, provider_path):
        pass

    def _read_convergence(self, filename):

        # Parse the XML file and first check for consistency.
        filecontents = ET.parse(filename).getroot()
        dftparams = filecontents.find("calculationparameters")
        if dftparams.find("element").text != self.parameters.element or \
           dftparams.find("crystal_structure").text != self.parameters.crystal_structure or \
           dftparams.find("dft_calculator").text != self.parameters.dft_calculator or \
           float(dftparams.find("temperature").text) != self.parameters.temperature or \
           int(dftparams.find("number_of_atoms").text) != self.parameters.number_of_atoms:
            raise Exception("Incompatible convergence parameters provided.")

        cutoff_energy = int(filecontents.find("cutoff_energy").text)
        kpoints = filecontents.find("kpoints")
        kgrid = (int(kpoints.find("kx").text),int(kpoints.find("kx").text),
                 int(kpoints.find("kx").text))

        return cutoff_energy, kgrid
