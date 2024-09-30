# Written in 2018 by Dr. Robert Pantazes. This file defines a Python error for
# issues that arise in Protein - structure objects

class ProteinError (Exception):
    """An error for problems in protein - structure objects."""

    def __init__ (self, error = ''):
        """The initialization of the ProteinError class."""
        if isinstance(error, str):
            self.error = error
        else:
            self.error = ''

    def __str__ (self):
        if not isinstance(self.error, str) or len(self.error) == 0:
            return ''
        else:
            return "\n" + self.error
