# Written in 2018 by Dr. Robert Pantazes. Copyright owned by Auburn
# University.
# This file contains the implementation of a Python extension-type wrapper of
# a C++ Protein class. This class is intended to be a primary worker in
# computational protein engineering tasks involving PDB formatted data

# Cython compiler directives
# cython: c_string_type=str
# cython: c_string_encoding=ascii
# cython: language_level=3
# distutils: language=c++

# Include the standard Python modules I always include in Python codes
import os
import sys
# Include the Protein Error class
from Proteins.ProteinError import ProteinError

# C-import C++ classes
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp cimport bool

# C-import the C++ classes and the extension types for the other PDB data
# levels
from Proteins.Proteins cimport Matrix
from Proteins.Proteins cimport Atom
from Proteins.Proteins cimport Residue
from Proteins.Proteins cimport Protein
from Proteins.Matrix.PyMatrix cimport PyMatrix
from Proteins.Atom.PyAtom cimport PyAtom
from Proteins.Residue.PyResidue cimport PyResidue
from Proteins.Proteins cimport coor

# Define the PyProtein class
cdef class PyProtein:

    def __cinit__ (self, fileName = None, path = None):
        """C-initialize a PyProtein object."""
        # Set the pointer to null
        self._thisptr = NULL
        # Set the delete flag to false
        self.delete_flag = False
        # C-declare two C++ strings
        cdef string part1, part2
        # If both inputs are None
        if fileName == None and path == None:
            # Be done
            return
        # If the path does not exist but the file name is a string
        elif isinstance(fileName, str) and path == None:
            part1 = fileName
            self._thisptr = new Protein (part1, part2)
            self.delete_flag = True
        # If both are strings
        elif isinstance(fileName, str) and isinstance(path, str):
            part1 = fileName
            part2 = path
            self._thisptr = new Protein (part1, part2)
            self.delete_flag = True
        # If the fileName is instead a True value, create an empty Protein
        elif fileName == True:
            self._thisptr = new Protein ()
            self.delete_flag = True
        # Otherwise raise an error
        else:
            text = "A PyProtein can only be initialized directly from a file.\n"
            text += "Inputs of " + type(fileName) + " and "
            text += type(path) + " are not acceptable.\n"
            raise ProteinError (text)

    def __dealloc__ (self):
        """Deallocate dynamically allocated memory"""
        if self._thisptr != NULL and self.delete_flag:
            del self._thisptr

    def __init__ (self, fileName = None, path = None):
        """The PyProtein initialization function."""
        # Initialization is done in cinit. This function does nothing
        pass

    cdef void _check_for_null (self):
        """Check to see if the Protein is a null pointer"""
        if self._thisptr == NULL:
            text = "Operations may not be performed on NULL PyProteins.\n"
            raise ProteinError (text)

    def format (self, bool internal = False):
        """Create a formatted string of text of the Protein's Residues"""
        self._check_for_null()
        return self._thisptr.str(internal)

    def __str__ (self):
        """Create a formatted string of text of the Protein's contents"""
        return self.format(False)

    def output (self, str fileName, bool internal = False):
        """Output the Protein to a file"""
        self._check_for_null()
        # Open the file
        f = open(fileName, "w")
        # Write the formatted text
        f.write(self._thisptr.str(internal))
        f.write("END\n")
        # Close the file
        f.close()

    def load (self, str fileName, str path = ""):
        """Load the protein from a file"""
        self._check_for_null()
        self._thisptr.load(fileName, path)

    def number_of_atoms (self):
        """Get the number of Atoms in the Protein"""
        self._check_for_null()
        return self._thisptr.number_of_atoms()

    def last_atom_number (self):
        """Get the number of the last Atom in the Protein."""
        self._check_for_null ()
        return self._thisptr.last_atom_number()

    def last_residue_number (self):
        """Get the number of the last residue in the protein."""
        self._check_for_null ()
        return self._thisptr.last_residue_number()

    def __len__ (self):
        """The number of residues in the protein"""
        self._check_for_null()
        return self._thisptr.size()

    def name (self):
        """The name of the protein"""
        self._check_for_null()
        cdef char L = self._thisptr.name()
        cdef string temp
        temp.push_back(L)
        cdef str output = temp
        return output

    def renumber_residues (self, long n):
        """Renumber the residues in the protein."""
        self._check_for_null()
        return self._thisptr.renumber_residues(n)

    def renumber_atoms (self, long n):
        """Renumber the atoms in the protein."""
        self._check_for_null()
        return self._thisptr.renumber_atoms(n)

    def __getitem__ (self, long n):
        """Access a Residue in the Protein"""
        self._check_for_null()
        # Create an empty PyResidue
        cdef PyResidue res = PyResidue ()
        # Get a pointer to a residue from the protein
        cdef char L = 32
        cdef bool how = True
        cdef Residue * ptr = self._thisptr[0](n, L, how)
        res._thisptr = ptr
        # Set the delete flag of the residue to false
        res.delete_flag = False
        return res

    def move (self, PyMatrix matrix, bool subtract = True):
        """Move the Protein using the information in a matrix"""
        self._check_for_null()
        self._thisptr.move(matrix._thisptr, subtract)

    def rotate (self, PyMatrix matrix):
        """Rotate the Protein using the information in the matrix."""
        self._check_for_null()
        self._thisptr.rotate(matrix._thisptr)

    def center (self, bool backboneOnly = False):
        """Move a protein so that its center of mass is at the origin"""
        self._check_for_null()
        # Make an empty PyMatrix
        cdef PyMatrix output = PyMatrix()
        # Center the protein and get the matrix result
        cdef Matrix result = self._thisptr.center(backboneOnly)
        # use that to fill in the values in the matrix that will be returned
        output._thisptr.allocate(&result)
        return output

    def purify (self):
        """Return an amino-acid only version of this protein."""
        self._check_for_null()
        # allocate a vector of Residues
        cdef vector[Residue] residues
        # A pointer to residues
        cdef Residue * res
        # Looping variables
        cdef size_t i, j
        j = self._thisptr.size()
        # Loop through the protein's residues
        for i in range(j):
            # Get the Residue
            res = self._thisptr[0](i, 32, True)
            # If the residue is an amino acid with atoms
            if res.is_amino_acid() and res.size() > 0:
                # Store the Residue in the vector
                residues.push_back(res[0])
        # Create a new PyProtein that allocates a blank protein
        cdef PyProtein prot = PyProtein (True)
        # If there are residues
        if residues.size() > 0:
            # Load the residues into that protein
            prot._thisptr.load(residues)
        return prot

    def set_name (self, str text):
        """Set the name of the protein"""
        self._check_for_null()
        cdef string temp = text
        cdef char L
        if len(text) > 0:
            L = temp[0]
        self._thisptr.set_name(L)

    def charmm_his_prep (self):
        """Change residue names from HIS to HSD for CHARMM."""
        self._check_for_null ()
        self._thisptr.charmm_his_prep ()

    def charmm_his_fix (self):
        """Change residue names from HSD to HIS after CHARMM."""
        self._check_for_null ()
        self._thisptr.charmm_his_fix ()

    def get_sequence (self):
        """Get a string of the Protein's sequence."""
        self._check_for_null()
        cdef str text
        text = self._thisptr.get_sequence()
        return text

    def calculate_dihedrals (self, bool omega = False):
        """Calculate the Protein's dihedral angles."""
        self._check_for_null()
        self._thisptr.calculate_dihedrals(omega)

    def calculate_RMSD (self, PyProtein other, bool backboneOnly = False):
        """Calculate the RMSD between two PyProteins"""
        # Make sure both proteins have valid pointers
        self._check_for_null()
        other._check_for_null()
        # Calculate the rotation matrix between the proteins
        cdef PyMatrix matrix = PyMatrix()
        matrix.minimized_RMSD (self, other, backboneOnly)
        # Rotate the other protein
        other._thisptr.rotate(matrix._thisptr)
        # Calculate the RMSD value
        return self._thisptr.calculate_RMSD(other._thisptr, backboneOnly)

    def parameterize (self, parameters, bool VDW, bool ELEC, bool LK, 
                      bool GB = False, bool VDW_H = True):
        """Assign parameters to the Protein's Atoms for energy calculations"""
        # Check that the protein is not empty
        self._check_for_null()
        # Use this PyResidue for looping through the Protein's residues
        cdef PyResidue residue = PyResidue ()
        # Make sure it does not delete information when done
        residue.delete_flag = False
        # Use these values for indicating if this is the first or last residue
        # in the protein
        cdef bool NTER, CTER
        # Variables for accessing residues
        cdef char L = 32
        cdef bool how = True
        # Use these variables for looping through the residues
        cdef long i, j
        j = self._thisptr.size()
        # Loop through the residues
        for i in range(j):
            # If this is the first residue, NTER is True
            if i == 0:
                NTER = True
            else:
                NTER = False
            # If this is the last residue, CTER is True. NTER and CTER are not
            # exclusive of one another
            if i == j-1:
                CTER = True
            else:
                CTER = False
            # Put the pointer to the current residue in the PyResidue
            residue._thisptr = self._thisptr[0](i, L, how)
            # Parameterize the residue
            residue.parameterize(parameters, NTER, CTER, VDW, ELEC, LK, GB,
                                 VDW_H)

    def disable_VDW (self):
        """Turn off VDW energy calculations for the Protein's atoms"""
        # Make sure the protein is valid
        self._check_for_null()
        # Looping variables
        cdef size_t i, j, k, l
        # A residue pointer
        cdef Residue * residue
        # Residue access variables
        cdef bool how = True
        cdef char L = 32
        # An atom pointer
        cdef Atom * atom
        # The boolean off value
        cdef bool off = False
        # Loop through the residues
        k = self._thisptr.size()
        for i in range(k):
            residue = self._thisptr[0](i, L, how)
            # Loop through the Residue's atoms
            l = residue.size()
            for j in range(l):
                atom = residue[0][j]
                atom.set_calculate_vdw(off)

    cpdef size_t calculate_atom_overlap (self, PyAtom atom):
        """Count the VDW overlaps with an Atom"""
        # Check the pointers
        self._check_for_null()
        atom._check_for_null_pointer()
        # Use the C++ function
        return self._thisptr.calculate_vdw_overlap(atom._thisptr)

    cpdef size_t calculate_residue_overlap (self, PyResidue residue, 
                                            bool rotamer = False):
        """Count the VDW overlaps with a Residue"""
        # Validate the pointers
        self._check_for_null()
        residue._check_for_null()
        # Use the C++ function
        return self._thisptr.calculate_vdw_overlap(residue._thisptr, rotamer)

    cpdef size_t calculate_protein_overlap (self, PyProtein other):
        """Count the VDW clashes between two Proteins"""
        # Validate the pointers
        self._check_for_null()
        other._check_for_null()
        # Use the C++ function
        return self._thisptr.calculate_vdw_overlap (other._thisptr)

    cpdef float calculate_atom_energy (self, PyAtom atom):
        """Calculate the energy of the Protein with an Atom"""
        # Check the pointers
        self._check_for_null()
        atom._check_for_null_pointer()
        # Use the C++ function
        return self._thisptr.calculate_energy(atom._thisptr)

    cpdef float calculate_residue_energy (self, PyResidue residue, 
                                          bool rotamer = False):
        """Calculate the energy of the Protein with a Residue"""
        # Check the pointers
        self._check_for_null()
        residue._check_for_null()
        # Use the C++ function
        return self._thisptr.calculate_energy(residue._thisptr, rotamer)

    cpdef float calculate_protein_energy (self, PyProtein other):
        """Calculate the energy of the Protein with another Protein"""
        # Validate the pointers
        self._check_for_null()
        other._check_for_null()
        # Use the C++ function
        return self._thisptr.calculate_energy(other._thisptr)

    cpdef PyProtein duplicate (self):
        """Create a duplicated copy of the Protein"""
        # Validate the pointer
        self._check_for_null()
        # Create an empty PyProtein
        cdef PyProtein output = PyProtein()
        # Set it's pointer to a new Protein created by duplicating this one
        output._thisptr = new Protein (self._thisptr.duplicate())
        # Set the delete flag of the output to True
        output.delete_flag = True
        return output
