"""Provider for optimal LDOS calculation parameters."""
import os
from shutil import copyfile
from .provider import Provider


class LDOSConvergenceProvider(Provider):
    """
    Determine number of k points and energy levels needed for a smooth LDOS.

    Parameters
    ----------
    parameters : malada.utils.parametes.Parameters
        Parameters used to create this object.

    external_ldos_configuration : string
        Path xml file containing k grid for LDOS creation and energy levels.
        If not None, no DFT calculations will be performed.

    """

    def __init__(self, parameters, external_ldos_configuration=None):
        super(LDOSConvergenceProvider, self).__init__(parameters)
        self.parameters = parameters
        self.external_ldos_configuration = external_ldos_configuration
        self.ldos_configuration_file = None

    def provide(self, provider_path, snapshot_file, dft_convergence_file):
        """
        Provide correct number of k points and energy levels to calculate LDOS.

        The results of LDOS based (ML) workflows is heaviliy dependent on
        a correct choice of energy levels and the k-grid. If those are
        chosen wrongly, any neural network training will be hard and yet
        give insufficient accuracies.

        Parameters
        ----------
        provider_path : string
            Path in which to operate in.

        dft_convergence_file : string
            Path to xml file containing the DFT convergence parameter.

        snapshot_file : string
            Path to a file containing an ASE trajectory containing atomic
            snapshots for DFT/LDOS calculation.
        """
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


