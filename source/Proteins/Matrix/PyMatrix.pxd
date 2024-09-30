# Written in 2018 by Dr. Robert Pantazes. Copyright owned by Auburn
# University.
# This file contains the declaration of the Python Matrix class for linear
# algebra calculations involving C++ Protein structure data-types

# Cython compiler directives
# cython: c_string_type=str
# cython: c_string_encoding=ascii
# cython: language_level=3
# distutils: language=c++

# Include the C++ Matrix class
from Proteins.Proteins cimport Matrix

# Declare the class
cdef class PyMatrix:

    # A Pointer to the C++ matrix object
    cdef Matrix * _thisptr

    # the declaration of a function that does an SVD decomposition of a Matrix
    cdef void SVD_for_RMSD (self)
