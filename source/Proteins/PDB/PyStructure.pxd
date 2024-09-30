# Written in 2018 by Dr. Robert Pantazes. Copyright owned by Auburn
# University.
# This file contains the C-level declaration of a Python extension-type
# wrapper for the C++ structure class. That class is only used in the context
# of PDB files

# Cython compiler directives
# cython: c_string_type=str
# cython: c_string_encoding=ascii
# cython: language_level=3
# distutils: language=c++

# Include the Structure class
from Proteins.Proteins cimport Structure
# And a boolean value
from libcpp cimport bool

# C-declare the extension type
cdef class PyStructure:

    # A pointer to the object
    cdef Structure * _thisptr

    # A flag to say whether or not deletion should occur
    cdef bool delete_flag

    # An error checking function
    cdef void _check_for_null (self)
