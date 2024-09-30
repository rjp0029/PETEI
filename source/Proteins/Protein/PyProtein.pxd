# Written in 2018 by Dr. Robert Pantazes. Copyright owned by Auburn
# University.
# This file contains the c-level declaration of a Python extension-type
# wrapper for a C++ Protein class.

# Cython compiler directives
# cython: c_string_type=str
# cython: c_string_encoding=ascii
# cython: language_level=3
# distutils: language=c++

# Include the C++ Protein class
from Proteins.Proteins cimport Protein
# Include the Py-wrappers for the Atom and Residue classes
from Proteins.Atom.PyAtom cimport PyAtom
from Proteins.Residue.PyResidue cimport PyResidue
# Include the C++ true / false boolean value
from libcpp cimport bool

# C-declare the class
cdef class PyProtein:

    # A pointer to a C++ Protein
    cdef Protein * _thisptr

    # A flag to indicate when a PyProtein should not delete its own
    # information because it came from another source instead of generating
    # itself
    cdef bool delete_flag

    # An error checking function
    cdef void _check_for_null (self)

    # C++ and Python declared functions
    cpdef size_t calculate_atom_overlap (self, PyAtom atom)
    cpdef size_t calculate_residue_overlap (self, PyResidue residue, 
                                            bool rotamer =*)
    cpdef size_t calculate_protein_overlap (self, PyProtein other)
    cpdef float calculate_atom_energy (self, PyAtom atom)
    cpdef float calculate_residue_energy (self, PyResidue residue, 
                                          bool rotamer =*)
    cpdef float calculate_protein_energy (self, PyProtein other)

    cpdef PyProtein duplicate (self)
