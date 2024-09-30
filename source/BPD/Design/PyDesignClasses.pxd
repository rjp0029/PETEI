# Written in 2018 by Dr. Robert Pantazes. Copyright owned by Auburn
# University.
# This file contains Cython declarations of C++ classes that are useful for
# designing binding proteins using a database of protein pieces

# Cython compiler directives
# cython: c_string_type=str
# cython: c_string_encoding=ascii
# cython: language_level=3
# distutils: language=c++

# Include the Atom C++ class
from Proteins.Proteins cimport Atom
# Include the PyProtein and PyResidue classes
from Proteins.Protein.PyProtein cimport PyProtein
from Proteins.Residue.PyResidue cimport PyResidue
# Include the C++ design classes
from BPD.Design.DesignClasses cimport Position
from BPD.Design.DesignClasses cimport Interaction
from BPD.Design.DesignClasses cimport Solution
# Include the vector, string and bool classes
from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp cimport bool

# Declare the LoopStructure class
cdef class LoopStructure:

    # The attributes contained in the class
    # The Protein structure
    cdef public PyProtein protein
    # Whether or not the structure is usable in a solution
    cdef public bool usable
    # The name of the loop the structure belongs to
    cdef readonly str loop
    # The number of the structure (starts at 1)
    cdef readonly size_t number
    # Boolean values storing compatibility information with every other loop
    # structure
    cdef vector[vector[bool]] compatibility
    # A dictionary for converting loop names into indexes
    cdef dict loopIndexes
    # Whether or not the structure is compatible with another
    cpdef bool compatible (self, str, size_t)

# Declare the Loop class
cdef class Loop:

    # A list of structures in the loop
    cdef public list structures
    # The name of the loop
    cdef readonly str name
    # The number of structures in the loop
    cdef readonly size_t size

# Declare the Interaction Type class
cdef class InteractionType:

    # A vector of the interactions
    cdef vector[Interaction] interactions
    # The scoring threshold
    cdef float threshold
    # The number of the interaction
    cdef readonly size_t number
    # The functions that find the positions from the interactions use these
    # three variables
    cdef string res
    cdef vector[Atom *] atoms
    cdef char ProtChain

    # C-defined methods of the Interaction Type class
    cdef void get_atoms (self, PyResidue, list, bool)
    cdef void atom_swap (self, size_t, size_t)
    cdef void evaluate (self, vector[Position]&)
    cdef void F1 (self, PyResidue, vector[Position]&)
    cdef void F2 (self, PyResidue, vector[Position]&)
    cdef void F3 (self, PyResidue, vector[Position]&)
    cdef void F4 (self, PyResidue, vector[Position]&)
    cdef void F5 (self, PyResidue, vector[Position]&)
    cdef void F6 (self, PyResidue, vector[Position]&)
    cdef void F7 (self, PyResidue, vector[Position]&)
    cdef void F8 (self, PyResidue, vector[Position]&)
    cdef void F9 (self, PyResidue, vector[Position]&)
    cdef void F10 (self, PyResidue, vector[Position]&)
    cdef void F11 (self, PyResidue, vector[Position]&)
    cdef void F12 (self, PyResidue, vector[Position]&)
    cdef void F13 (self, PyResidue, vector[Position]&)
    cdef void F14 (self, PyResidue, vector[Position]&)
    cdef void F15 (self, PyResidue, vector[Position]&)
    cdef void F16 (self, PyResidue, vector[Position]&)
    cdef void F17 (self, PyResidue, vector[Position]&)
    cdef void F18 (self, PyResidue, vector[Position]&)
    cdef void F19 (self, PyResidue, vector[Position]&)
    cdef void F20 (self, PyResidue, vector[Position]&)
    cdef void F21 (self, PyResidue, vector[Position]&)
    cdef void F22 (self, PyResidue, vector[Position]&)
    cdef void F23 (self, PyResidue, vector[Position]&)
    cdef void find_positions (self, vector[Position]&, list)

# Declare the PySolution class
cdef class PySolution:

    # A pointer to the C++ solution object
    cdef Solution * _thisptr

    # The C-level class methods
    cdef void _check_for_null (self) except +
    cdef void allocate (self, vector[Position]&, vector[size_t]&, float,
                        float, float)
    cdef size_t size (self) except +
    cdef Position * position (self, size_t) except +
    cpdef float X (self)
    cpdef float Y (self)
    cpdef float Z (self)
    cpdef float XA (self)
    cpdef float YA (self)
    cpdef float ZA (self)
    cdef dict list_loops (self)
    cpdef float similarity (self, PySolution)

# Declare the solution storage class
cdef class SolutionStorage:

    # The information stored in the class
    cdef public dict solutions
    cdef public list order

    # The C-level methods of the class
    cdef void clean_up (self)
    cdef void store_solutions (self, vector[Position]&, 
              vector[vector[vector[size_t]]]&, float, float, float)


