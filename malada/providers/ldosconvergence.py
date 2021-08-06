import os
import numpy as np
import glob

from shutil import copyfile
import ase
import ase.io
from ase.units import Rydberg
from xml.etree.ElementTree import Element, SubElement, tostring
from xml.dom import minidom

import malada
from malada.utils.convergence_guesses import *
from malada.utils.custom_converter import *
from malada.utils.vasp_utils import VaspUtils
from .provider import Provider


class LDOSConvergenceProvider(Provider):
    """Determine number of k points needed for a smooth LDOS."""

    def __init__(self, parameters, external_ldos_configuration=None):
        super(LDOSConvergenceProvider, self).__init__(parameters)
        self.parameters = parameters
        self.external_ldos_configuration = external_ldos_configuration
        self.ldos_configuration_file = None

    def provide(self, provider_path, snapshot_file, dft_convergence_file):
        file_name = self.parameters.element + \
                    str(self.parameters.number_of_atoms) + \
                    "_" + self.parameters.crystal_structure +\
                    "_" + str(self.parameters.temperature) +\
                    "_" + self.parameters.dft_calculator+".ldosconf.xml"
        self.ldos_configuration_file = os.path.join(provider_path, file_name)
        # Check if there exist results or if we have to work from scratch.
        if self.external_ldos_configuration is None:
            raise Exception("Currently there is no way to determine LDOS"
                            "configuration automatically.")
        else:
            copyfile(self.external_ldos_configuration,
                     self.ldos_configuration_file)
            print("Getting <<ldos_configuration>>.xml file from disc.")


