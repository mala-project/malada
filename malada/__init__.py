"""
MALA Data Acquistion.

Can be used to automate a data acquistion pipeline for MALA.
"""


from .version import __version__
from .utils import Parameters
from .providers import CrystalStructureProvider, SuperCellProvider,\
                       DFTConvergenceProvider, MDPerformanceProvider,\
                       MDProvider, SnapshotsProvider, LDOSConvergenceProvider
from .pipeline import DataPipeline
from .runners import Runner, BashRunner

