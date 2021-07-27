from abc import ABC, abstractmethod

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
