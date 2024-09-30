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
# Include the vector, string and bool classes
from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp cimport bool

# Include contents from the INTERaction namespace in the DesignClasses header
# file
cdef extern from "DesignClasses.h" namespace "INTER":

    # Declare the Spot class
    cdef cppclass Spot:
        Spot ()
        Spot (const float, const float, const float)
        Spot (const Spot&)
        size_t size () const
        float operator[] (const size_t) except +
        Spot calculate_movement (const Atom *) const
        float calculate_distance (const Atom *, const Spot *) const
        string str () const

    # Declare the Plane class
    cdef cppclass Plane:
        Plane ()
        Plane (const float, const float, const float, const float)
        Plane (const Plane&)
        size_t size () const
        float operator[] (const size_t) except +
        float calculate_distance (const Atom *, const Spot *) const

    cdef cppclass Interaction:
        Interaction ()
        Interaction (const size_t, const string&)
        Interaction (const Interaction&)
        string loop () const
        size_t structure () const
        size_t position () const
        size_t rotamer () const
        size_t kind () const
        Spot movement (const vector[Atom *]&) const
        bool evaluate (const vector[Atom *]&, const float) const

    cdef cppclass Position:
        Position ()
        Position (Interaction *, const vector[Atom *]&, const string&, 
                  const char)
        Position (const Position&)
        float score () const
        string residue () const
        char protein () const
        string loop () const
        size_t structure () const
        size_t position () const
        size_t rotamer () const
        size_t interaction () const
        Spot * movement ()
        float operator[] (const size_t) const
        bool distance_compatible (const Position *) const
        bool compatible (const Position *, const float) const
        string str () const

    cdef cppclass Solution:
        Solution ()
        Solution (vector[Position]&, vector[size_t]&, const float, 
                  const float, const float)
        Solution (const Solution&)
        size_t size () const
        float X () const
        float Y () const
        float Z () const
        float XA () const
        float YA () const
        float ZA () const
        Position * operator[] (const size_t)

