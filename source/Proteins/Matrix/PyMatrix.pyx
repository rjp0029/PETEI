# Written in 2018 by Dr. Robert Pantazes. Copyright owned by Auburn
# University.
# This file contains the implementation of the Python-level functions of the
# PyMatrix class. This class is a wrapper of the C++ Matrix class. 

# Cython compiler directives
# cython: c_string_type=str
# cython: c_string_encoding=ascii
# cython: language_level=3
# distutils: language=c++

# Include the "standard" Python modules that I always include
import os
import sys
# Include the "Numbers for Python" module that is needed for calculating
# rotation matrices
import numpy
# Include the Protein Error class
from Proteins.ProteinError import ProteinError

# Include the C++ class definition
from Proteins.Proteins cimport Matrix
# Include the C++ protein-structure definitions
from Proteins.Proteins cimport Atom
from Proteins.Atom.PyAtom cimport PyAtom
from Proteins.Protein.PyProtein cimport PyProtein
# Include the coordinate typedefinition
from Proteins.Proteins cimport coor
# Include the C++ boolean variable
from libcpp cimport bool

# C-declare the extension type
cdef class PyMatrix:
    """An extension type for a C++ Matrix data type."""

    def __cinit__ (self):
        """The c-level initialization of the PyMatrix class."""
        # Dynamically allocate a new Matrix for this object 
        self._thisptr = new Matrix ()

    def __dealloc__ (self):
        """Deallocate the Matrix's memory"""
        if self._thisptr != NULL:
            del self._thisptr

    def __init__ (self):
        """Initialize a PyMatrix"""
        # Because of some of the rules associated with cython, it isn't
        # possible to set up this function to do all of the different
        # initializations I might want. They will be separated into other
        # functions explicitly. So this function does nothing
        pass

    def dimension_allocation (self, size_t r, size_t c):
        """Allocate the Matrix's dimensions explicitly"""
        self._thisptr.allocate(r, c)

    def PyAtom_allocation (self, PyAtom atom):
        """Allocate the Matrix's dimensions from an Atom."""
        self._thisptr.allocate(atom._thisptr)

    def PyAtoms_allocation_for_centering (self, atomList):
        """Create a PyMatrix that will center the provided Atoms"""
        # Allocate a matrix with 1 row and ATOMSIZE (3) coordinates
        self._thisptr.allocate(1, len(atomList[0]))
        # Loop through the coordinates, calculating averages as it goes
        cdef float average
        for i in range(len(atomList[0])):
            average = 0.0
            for atom in atomList:
                average += atom[i]
            average /= float(len(atomList))
            self._thisptr.set(0, i, average)

    def rows (self):
        """The number of rows in the matrix"""
        return self._thisptr.rows()

    def columns (self):
        """The number of columns in the matrix"""
        return self._thisptr.columns()

    def __call__ (self, size_t r, size_t c):
        """Access to a specific entry in the matrix"""
        return self._thisptr[0](r, c)

    def __str__ (self):
        """A String representation of the matrix"""
        return self._thisptr.str()

    cdef void SVD_for_RMSD (self):
        """Carry out a singular value decomposition to minimize RMSD"""
        # Get the number of rows (i.e. Atom coordinates). This function will
        # only ever be called from within the class, so it is assumed the
        # matrix is known to be square
        cdef size_t r = self._thisptr.rows()
        # Make a r x r numpy array
        mat = numpy.zeros(shape=(r, r))
        # Store the values in the matrix
        cdef size_t i, j
        for i in range(r):
            for j in range(r):
                mat[i][j] = self._thisptr[0](i, j)
        # Do the linear value decomposition
        v, s, w = numpy.linalg.svd(mat)
        if numpy.linalg.det(v) * numpy.linalg.det(w) < 0:
            for i in range(r):
                v[i][r-1] = -v[i][r-1]
        # Calculate the rotation matrix
        mat = numpy.dot(v, w)
        for i in range(r):
            for j in range(r):
                self._thisptr.set(i, j, mat[i][j])

    def minimized_RMSD (self, PyProtein p1, PyProtein p2, 
                        bool backboneOnly = False):
        """Generate a matrix to minimize the RMSD between the 2 proteins"""
        # The matrix must be applied to the second Protein
        # Create a Matrix from the first Protein, transposed
        cdef Matrix first = Matrix(p1._thisptr, backboneOnly).transpose()
        # Create a Matrix from the second Protein
        cdef Matrix second = Matrix (p2._thisptr, backboneOnly)
        # Calcualte the product of the two matrices
        cdef Matrix prod = first.product(&second)
        # Allocate that in the pointer
        self._thisptr.allocate(&prod)
        # Do the singular value decomposition 
        self.SVD_for_RMSD()

    def minimized_RMSD_atoms (self, atomList1, atomList2):
        """Calculate a rotation matrix to minimize the RMSD between the atoms"""
        # Make sure that the atom lists are the same non-zero length
        if len(atomList1) == 0 or len(atomList1) != len(atomList2):
            error = "Invalid inputs to the minimized_RMSD_atoms function.\n"
            raise ProteinError(error)
        # Create two matrices using the atoms
        cdef size_t i, j, k, l
        cdef coor value
        k = len(atomList1)
        l = len(atomList1[0])
        cdef Matrix first = Matrix(l, k)
        for i in range(k):
            for j in range(l):
                value = atomList1[i][j]
                first.set(j, i, value)
        cdef Matrix second = Matrix(k, l)
        for i in range(k):
            for j in range(l):
                value = atomList2[i][j]
                second.set(i, j, value)
        # Calculate the product of the two matrices
        cdef Matrix prod = first.product(&second)
        # Allocate that into this object
        self._thisptr.allocate(&prod)
        # Do the SVD
        self.SVD_for_RMSD()

    def specified_rotation (self, coor angle, list values):
        """Generate a rotation matrix using Rodriguez's Rotation Formula"""
        # If the list is not the proper length, raise an error
        if len(values) != 3:
            text = "The specified_rotation function requires a list of 3 "
            text += "values, not " + str(len(list)) + ".\n"
            raise ProteinError (text)
        # c-define the necessary array
        cdef coor nums [3]
        # Go through the values and store them
        cdef int i, j
        j = len(values)
        for i in range(j):
            # Make sure the data type is acceptable
            if not isinstance(values[i], float):
                text = "The specified rotation function only works with "
                text += "float variables, not " + type(values[i]) + ".\n";
                raise ProteinError (text)
            # Store the value in the array
            nums[i] = values[i]
        # Allocate the matrix
        self._thisptr.allocate(angle, nums)
