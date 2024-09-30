# Written in 2018 by Dr. Robert Pantazes. Copyright owned by Auburn
# University.
# This file contains the implementation of the Python-level functions of the
# PyAtom class. This class is a wrapper for a C++ Atom class for working with
# PDB-formatted data

# Cython compiler directives
# cython: c_string_type=str
# cython: c_string_encoding=ascii
# cython: language_level=3
# distutils: language=c++

# Include the "standard" Python modules
import os
import sys

# Include the C++ Atom class
from Proteins.Proteins cimport Atom
# Include the C++ string class
from libcpp.string cimport string
# Include the Python Protein error class
from Proteins.ProteinError import ProteinError
# Include the PyResidue and PyProtein classes
from Proteins.Residue.PyResidue cimport PyResidue
from Proteins.Protein.PyProtein cimport PyProtein
# And the Residue class
from Proteins.Proteins cimport Residue
# And the Matrix class
from Proteins.Matrix.PyMatrix cimport PyMatrix
from Proteins.Proteins cimport Matrix
# Include true / false values
from libcpp cimport bool

# C-declare the extension type
cdef class PyAtom:
    """An extension type for a C++ class of a PDB Atom."""

    def __cinit__ (self, text = None):
        """The c-level initialization of the PyAtom class."""
        # Set the self-pointer (declared in the pxd file) to Null
        self._thisptr = NULL
        # Set the delete flag to True
        self.delete_flag = True
        # C-declare a string, since it can't be done in the if statement
        cdef string data
        # If the text is a Python String, create a new Atom using it
        if isinstance (text, str):
            # Create a new Atom using it
            data = text
            self._thisptr = new Atom (data)

    def __dealloc__ (self):
        """Deallocate the C-level pointer information"""
        if self._thisptr != NULL and self.delete_flag:
            del self._thisptr

    def __init__ (self, text = None):
        """The standard initialization function of the PyAtom class."""
        # Because this class is primarily a wrapper for the C++ class and the
        # cinit function does the initial memory allocation, this function
        # does not actually do anything
        pass

    cdef void _check_for_null_pointer (self):
        """If the pointer is null, raise a Protein Error."""
        if self._thisptr == NULL:
            error = "No PyAtom operation may be performed on an unallocated "
            error += "PyAtom.\n"
            raise ProteinError (error)

    def __str__ (self):
        """Create a string represenation of the Atom"""
        # Make sure the pointer isn't NULL
        self._check_for_null_pointer()
        # Since it isn't, get the string of information
        cdef str text = self._thisptr.str()
        return text

    def __len__ (self):
        """The number of coordinates in an Atom"""
        self._check_for_null_pointer()
        return self._thisptr.size()

    def kind (self):
        """ATOM or HETATM"""
        self._check_for_null_pointer ()
        cdef str answer = self._thisptr.type()
        return answer

    def name (self):
        """The name of the Atom"""
        self._check_for_null_pointer ()
        cdef str answer = self._thisptr.name()
        return answer

    def residue (self):
        """The name of the Atom's residue"""
        self._check_for_null_pointer ()
        cdef str answer = self._thisptr.residue()
        return answer

    def element (self):
        """The type of element that the Atom is."""
        self._check_for_null_pointer ()
        cdef str answer = self._thisptr.element()
        return answer

    def charge (self):
        """The charge of the Atom."""
        self._check_for_null_pointer ()
        cdef str answer = self._thisptr.charge()
        return answer

    def number (self):
        """The Atom's number."""
        self._check_for_null_pointer()
        return self._thisptr.number()

    def residue_number (self):
        """The Atom's Residue's number."""
        self._check_for_null_pointer()
        return self._thisptr.residue_number ()

    def occupancy (self):
        """The Atom's occupancy"""
        self._check_for_null_pointer()
        return self._thisptr.occupancy()

    def temperature (self):
        """The Atom's temperature."""
        self._check_for_null_pointer()
        return self._thisptr.temperature()

    def alternative_location (self):
        """The Atom's Residue's alternative location information."""
        self._check_for_null_pointer()
        cdef char L = self._thisptr.alternative_location()
        cdef string data
        data.push_back(L)
        cdef str answer = data
        return answer

    def insertion_code (self):
        """The Atom's Residue's insertion code."""
        self._check_for_null_pointer()
        cdef char L = self._thisptr.insertion_code()
        cdef string data
        data.push_back(L)
        cdef str answer = data
        return answer

    def protein (self):
        """The name of the Atom's protein."""
        self._check_for_null_pointer()
        cdef char L = self._thisptr.protein()
        cdef string data
        data.push_back(L)
        cdef str answer = data
        return answer

    cpdef float elec_charge (self):
        """Get the Atom's electric charge"""
        self._check_for_null_pointer()
        return self._thisptr.get_elec_charge()

    def __getitem__ (self, value):
        """Access a coordinate in the Atom."""
        # Confirm that the Atom is valid
        self._check_for_null_pointer()
        # Cdef an integer for use in the Atom's access operator
        cdef size_t i
        # If the value is an integer
        if isinstance(value, int):
            i = value
        # Or an appropriate character
        elif value in ['x', 'X']:
            i = 0
        elif value in ['y', 'Y']:
            i = 1
        elif value in ['z', 'Z']:
            i = 2
        # Otherwise throw an error
        else:
            text = str(value) + " is not a valid Atom coordinate specifier.\n"
            raise ProteinError (text)
        # Return the appropriate coordinate from the Atom
        return self._thisptr[0][i]

    cpdef float calculate_distance (self, PyAtom other, bool squared = False):
        """Calculate the distance between two Atoms."""
        # Check this Atom
        self._check_for_null_pointer()
        # Check the other Atom's status
        other._check_for_null_pointer()
        # Calculate the distance
        return self._thisptr.calculate_distance(other._thisptr, squared)

    cpdef size_t calculate_vdw_overlap (self, PyAtom other):
        """Calculate the VDW overlaps between 2 Atoms (0 or 1)"""
        # Validate the pointers
        self._check_for_null_pointer()
        other._check_for_null_pointer()
        # Do the calculation
        return self._thisptr.calculate_vdw_overlap (other._thisptr)

    cpdef float calculate_energy (self, PyAtom other):
        """Calculate the energy between two Atoms."""
        # Validate the pointers
        self._check_for_null_pointer()
        other._check_for_null_pointer ()
        # Calculate the energy
        return self._thisptr.calculate_energy (other._thisptr)

    def parameterize (self, parameters, bool VDW, bool ELEC, bool LK, 
                      bool GB=False, bool VDW_H = True):
        """Assign non-bonded energy parameters to the Atom"""
        # Make sure the atom is valid
        self._check_for_null_pointer()
        # Use this boolean value to indicate what should be done for certain
        # tasks
        cdef bool doThis = True
        # And use this float value for the parameters
        cdef float value
        # Whether or not electrostatics calculations should be done
        doThis = ELEC
        self._thisptr.set_calculate_elec(doThis)
        # If they should be calculated, set the charge
        if ELEC:
            try:
                value = parameters.elec_charge
                self._thisptr.set_elec_charge(value)
            except AttributeError:
                error = "Missing electrostatics information for this Atom:\n"
                error += self.__str__()
                raise ProteinError (error)
        # If VDW should be calculated
        doThis = VDW
        # Get the Atom's name
        cdef str label = self._thisptr.name()
        # If the atom is a Hydrogen and VDW calculations should be done, only
        # include VDW parameters if Hydrogens should be used in VDW
        if label[0] == "H" and VDW and not VDW_H:
            doThis = VDW_H
        self._thisptr.set_calculate_vdw(doThis)
        # Set the VDW information
        if doThis:
            try:
                value = parameters.vdwE
                self._thisptr.set_vdw_epsilon(value)
                value = parameters.vdwR
                self._thisptr.set_vdw_radius(value)
                value = parameters.vdw14E
                self._thisptr.set_vdw_1_4_epsilon(value)
                value = parameters.vdw14R
                self._thisptr.set_vdw_1_4_radius(value)
                # Currently, the Lennard-Jones softening term is hard-coded to
                # a value of 1 (i.e. do not soften the radius)
                value = 1.0
                self._thisptr.set_vdw_softening(value)
                # Softening values are not currently set
            except AttributeError:
                error = "Missing VDW parameters for this Atom:\n"
                error += self.__str__()
                raise ProteinError (error)
        # Set the LK values
        doThis = LK
        self._thisptr.set_calculate_lk(doThis)
        if LK:
            try:
                value = parameters.lkR
                self._thisptr.set_lk_radius(value)
                value = parameters.lkL
                self._thisptr.set_lk_lambda(value)
                value = parameters.lkV
                self._thisptr.set_lk_volume(value)
                value = parameters.lkG
                self._thisptr.set_lk_gibbs(value)
            except AttributeError:
                error = "Missing LK Solvation data for this Atom:\n"
                error += self.__str__()
                raise ProteinError (error)
        # Generalize Borne Implicit solvation is not currently supported
        doThis = False
        self._thisptr.set_calculate_gb(doThis)

    def is_backbone_atom (self):
        """Determine whether or not the Atom is a backbone atom"""
        self._check_for_null_pointer()
        return self._thisptr.is_backbone_atom()
    
    def disable_VDW (self):
        """Turn off VDW calculations for the Atom"""
        # Make sure the atom exists
        self._check_for_null_pointer()
        # turn off the calculate_VDW value
        cdef bool off = False
        self._thisptr.set_calculate_vdw(off)

    def move (self, PyMatrix mat, bool subtract = True):
        """Move the Atom"""
        # Validate the pointer
        self._check_for_null_pointer()
        # Convert the boolean value into the proper character
        cdef str temp = "-"
        if not subtract:
            temp = "+"
        cdef string use = temp
        cdef char how = use[0]
        self._thisptr.move(mat._thisptr, how)

    def rotate (self, PyMatrix mat):
        """Rotate the Atom"""
        self._check_for_null_pointer()
        self._thisptr.rotate(mat._thisptr)
