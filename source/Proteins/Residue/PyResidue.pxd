# Written in 2018 by Dr. Robert Pantazes. Copyright owned by Auburn
# University.
# This file contains the C-level declarations of the PyResidue class for
# compilation in Cython. A Residue is a container of Atoms in a Protein

# Cython compiler directives
# cython: c_string_type=str
# cython: c_string_encoding=ascii
# cython: language_level=3
# distutils: language=c++

# Include the C++ Residue and Atom classes
from Proteins.Proteins cimport Residue
from Proteins.Proteins cimport Atom
from Proteins.Atom.PyAtom cimport PyAtom
from libcpp cimport bool

# Declare the class
cdef class PyResidue:

    # The class contains a pointer to a C++ Residue
    cdef Residue * _thisptr
    # A flag to make sure that PyResidues generated by PyProteins do not
    # delete the Protein's information
    cdef bool delete_flag

    # An error checking function
    cdef void _check_for_null (self)
    # Functions for accessing Atoms in the Residue
    cdef Atom * pointer_by_index (self, int)
    cdef Atom * pointer_by_name (self, str)
    
    # C++ and Python declared functions
    cpdef size_t calculate_atom_overlap (self, PyAtom atom, bool rotamer =*)
    cpdef size_t calculate_residue_overlap (self, PyResidue other, 
                 bool selfRotamer =*, bool otherRotamer =*)
    cpdef float calculate_atom_energy (self, PyAtom atom, bool rotamer =*)
    cpdef float calculate_residue_energy (self, PyResidue other, 
                bool selfRotamer =*, bool otherRotamer =*)
