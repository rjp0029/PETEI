# Written in 2018 by Dr. Robert Pantazes. Copyright owned by Auburn
# University.
# This file contains the C-level declaration of a Python extension-type
# wrapper for the C++ PDB class, which is intended for processing PDB
# formatted files

# Cython compiler directives
# cython: c_string_type=str
# cython: c_string_encoding=ascii
# cython: language_level=3
# distutils: language=c++

# Include all of the necessary classes
from Proteins.Proteins cimport PDB

# Declare the class
cdef class PyPDB:

    # A pointer to the object
    cdef PDB * _thisptr
