# Written in 2018 by Dr. Robert Pantazes. Copyright owned by Auburn
# University.
# This file contains the implementation of a Python Error class for use in
# designing binding proteins.

class BindingProteinError (Exception):
    """An error for problems with designing binding proteins"""

    def __init__ (self, error = ''):
        """The initialization of the BindingProteinError class"""
        if isinstance(error, str):
            self.error = error
        else:
            self.error = ''

    def __str__ (self):
        """A string representation of a Binding Protein Error"""
        if not isinstance(self.error, str) or len(self.error) == 0:
            return ''
        else:
            return "\n" + self.error
