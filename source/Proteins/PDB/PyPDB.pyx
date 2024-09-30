# Written in 2018 by Dr. Robert Pantazes. Copyright owned by Auburn
# University.
# This file contains the implementation of the PyPDB class, a Python
# extension-type wrapper for the C++ PDB class

# Cython compiler directives
# cython: c_string_type=str
# cython: c_string_encoding=ascii
# cython: language_level=3
# distutils: language=c++

# Include necessary modules and classes
import os
import sys

from Proteins.Proteins cimport PDB
from Proteins.Proteins cimport Protein
from Proteins.Proteins cimport Structure
from Proteins.Protein.PyProtein cimport PyProtein
from Proteins.PDB.PyStructure cimport PyStructure
from libcpp.string cimport string
from libcpp cimport bool

# Implement the class
cdef class PyPDB:

    def __cinit__ (self, str input1, str input2 = ""):
        """C-level initialization of a PyPDB object"""
        cdef string string1 = input1
        cdef string string2 = input2
        self._thisptr = new PDB (string1, string2)

    def __dealloc__ (self):
        """Deallocate class memory"""
        del self._thisptr

    def __init__ (self, str input1, str input2 = ""):
        """Initialization of a PyPDB object"""
        # Don't do anything here, everything is done in cinit
        pass

    def name (self):
        """The name of the PDB file"""
        return self._thisptr.name()

    def folder (self):
        """The location of the PDB file"""
        return self._thisptr.folder()

    def lines (self):
        """The number of lines loaded from the PDB file"""
        return self._thisptr.lines()

    def line (self, int i):
        """A specific line loaded from the file"""
        return self._thisptr.line(i)

    def kind (self):
        """The kind of experiment that generated the structure"""
        return self._thisptr.type()

    def resolution (self):
        """The experimental resolution of the PDB file"""
        return self._thisptr.resolution()

    def proteins (self):
        """The number of proteins in the file"""
        return self._thisptr.proteins()

    def protein (self, int i):
        """A specific protein from the file"""
        cdef Protein * ptr = self._thisptr.protein(i)
        cdef PyProtein prot = PyProtein ()
        prot._thisptr = ptr
        prot.delete_flag = False
        return prot

    def structures (self):
        """The number of protein structures in the file"""
        return self._thisptr.structures()

    def structure (self, int i):
        """Access to a specific protein structure group"""
        cdef Structure * ptr = self._thisptr.structure(i)
        cdef PyStructure struc = PyStructure ()
        struc._thisptr = ptr
        struc.delete_flag = False
        return struc

    def has_problem (self):
        """Whether or not the PDB file is appropriate for use"""
        return self._thisptr.obsolete()

    def __str__ (self):
        """A String representation of the file's contents"""
        return self._thisptr.str()
