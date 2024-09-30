# Written in 2018 by Dr. Robert Pantazes. Copyright owned by Auburn
# University.
# This file contains the declaration of the Python Atom class for packaging a
# C++ Atom class intended for working with PDB files

# Include some Cython compiler directives
# cython: c_string_type=str
# cython: c_string_encoding=ascii
# cython: language_level=3
# distutils: language=c++

# Include the C++ Atom class
from Proteins.Proteins cimport Atom
# And the boolean class
from libcpp cimport bool

# Declare the class
cdef class PyAtom:

    # A pointer to the C++ atom
    cdef Atom * _thisptr
    # A boolean flag to make sure that PyAtoms returned by PyResidues do NOT
    # delete the information in the Residue
    cdef bool delete_flag
    # An error-checking function
    cdef void _check_for_null_pointer (self)

    # C++ and Python declared functions
    cpdef float calculate_distance (self, PyAtom other, bool squared=*)
    cpdef size_t calculate_vdw_overlap (self, PyAtom other)
    cpdef float calculate_energy (self, PyAtom other)
    cpdef float elec_charge (self)
