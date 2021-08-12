"""Provider for optimal MD performance parameters."""
import os
from shutil import copyfile
from xml.etree.ElementTree import Element, SubElement, tostring
from xml.dom import minidom
from .provider import Provider


class MDPerformanceProvider(Provider):
    """
    Determine parallelization parameters for optimal MD performance.

    Parameters
    ----------
    parameters : malada.utils.parametes.Parameters
        Parameters used to create this object.

    external_performance_file : string
        Path to file containing optimal DFT-MD performance. If not None,
        no DFT-MD calculations will be done.
    """

    def __init__(self, parameters, external_performance_file=None):
        super(MDPerformanceProvider, self).__init__(parameters)
        self.md_performance_xml = None
        self.external_performance_file = external_performance_file

    def provide(self, provider_path, dft_convergence_file):
        """
        Provide parallelization parameters for optimal DFT-MD performance.

        This is not necessary when running in serial.

        Parameters
        ----------
        provider_path : string
            Path in which to operate in.

        dft_convergence_file : string
            Path to xml file containing the DFT convergence parameter.
        """
        file_name = self.parameters.element + \
                    str(self.parameters.number_of_atoms) + \
                    "_" + self.parameters.crystal_structure +\
                    "_" + str(self.parameters.temperature) +\
                    "_" + self.parameters.dft_calculator+".mdperformance.xml"
        self.md_performance_xml = os.path.join(provider_path, file_name)
        if self.external_performance_file is None:
            if self.parameters.run_system == "bash":
                # In this case, we can simply write an empty xml file
                print("Using bash based run system, no MD performance "
                      "optimization possible or necessary.")
                top = Element('mdperformanceparameters')
                dummy = SubElement(top, "nraise", {'type': "float"})
                dummy.text = "0.001"
                rough_string = tostring(top, 'utf-8')
                reparsed = minidom.parseString(rough_string)
                with open(self.md_performance_xml, "w") as f:
                    f.write(reparsed.toprettyxml(indent="  "))
            else:
                raise Exception("Currently there is no way to evaluate MD "
                                "performance on the fly.")
        else:
            copyfile(self.external_performance_file, self.md_performance_xml)
            print("Getting <<md_performance>>.xml file from disc.")
