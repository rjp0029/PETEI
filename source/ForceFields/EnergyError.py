# Written in 2018 by Dr. Robert Pantazes. This file defines a Python error for
# issues that arise in force fields or energy calculations.

class EnergyError (Exception):
    """An error for problems in force fields or energy calculations."""

    def __init__ (self, error = ''):
        """The initialization of the Energy Error class."""
        if isinstance(error, str):
            self.error = error
        else:
            self.error = ''

    def __str__ (self):
        """A String representation of an Energy Error"""
        if not isinstance(self.error, str) or len(self.error) == 0:
            return ''
        else:
            return "\n" + self.error
