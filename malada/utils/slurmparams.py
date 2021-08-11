from xml.etree.ElementTree import Element, SubElement, tostring, parse
from xml.dom import minidom


class SlurmParameters:
    """
    Class holding slurm parameters.
    """

    def __init__(self):
        self.scf_executable = "pw.x"
        self.module_loading_string = ""
        self.execution_time = 0
        self.partition_string = ""
        self.mpi_runner = "mpirun"

        # TODO: Maybe some sanity checks here?
        self.tasks_per_node = 0
        self.nodes = 0

    def save(self, filename):
        top = Element('slurmparameters')
        node = SubElement(top, "scf_executable",
                          {"type": "string"})
        node.text = self.scf_executable
        node = SubElement(top, "module_loading_string",
                          {"type": "string"})
        node.text = self.module_loading_string
        node = SubElement(top, "execution_time",
                          {"type": "int"})
        node.text = str(self.execution_time)
        node = SubElement(top, "partition_string",
                          {"type": "string"})
        node.text = self.partition_string
        node = SubElement(top, "mpi_runner",
                          {"type": "string"})
        node.text = self.mpi_runner
        node = SubElement(top, "tasks_per_node",
                          {"type": "int"})
        node.text = str(self.tasks_per_node)
        node = SubElement(top, "nodes",
                          {"type": "int"})
        node.text = str(self.nodes)
        rough_string = tostring(top, 'utf-8')
        reparsed = minidom.parseString(rough_string)
        with open(filename, "w") as f:
            f.write(reparsed.toprettyxml(indent="  "))

    @classmethod
    def from_xml(cls, filename):
        new_object = SlurmParameters()
        filecontents = parse(filename).getroot()
        new_object.scf_executable = filecontents.find("scf_executable").text
        new_object.module_loading_string = filecontents.find("module_loading_string").text
        new_object.execution_time = int(filecontents.find("execution_time").text)
        new_object.partition_string = filecontents.find("partition_string").text
        new_object.tasks_per_node = int(filecontents.find("tasks_per_node").text)
        new_object.nodes = int(filecontents.find("nodes").text)
        return new_object


