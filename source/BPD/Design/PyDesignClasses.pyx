# Written in 2018 by Dr. Robert Pantazes. Copyright owned by Auburn
# University.
# This file contains the Python level implementation of classes that are
# useful for searching for binding proteins designed to bind a user-specified
# epitope

# Cython compiler directives
# cython: c_string_type=str
# cython: c_string_encoding=ascii
# cython: language_level=3
# distutils: language=c++

# Include pure python modules
import os
import sys
# Include Binding Protein Design contents
from BPD.BindingProteinError import BindingProteinError
from BPD.General import parameterize_protein
from BPD.Design.DesignClasses cimport Position
from BPD.Design.DesignClasses cimport Interaction
# Include the Protein classes
from Proteins.Protein.PyProtein cimport PyProtein
from Proteins.Residue.PyResidue cimport PyResidue
from Proteins.Proteins cimport Residue
from Proteins.Proteins cimport Atom
# Include C++ classes
from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp cimport bool

cdef class LoopStructure:
    """The Structure class is a container for a binding loop structure"""

    def __init__ (self, DBFolder, LoopName, StructNum, instructions, LoopInfo, 
                  parameters):
        """Initialize a Structure object"""
        # Store the parameters about the Structure
        self.loop = LoopName
        self.number = StructNum
        # Indicate that the structure is usable
        self.usable = True
        # The folder the structure is located in
        folder = DBFolder + self.loop + "/Structures/"
        # The name of the structure file
        fileName = instructions["database"] + "_" + self.loop + "_" \
                 + str(StructNum) + ".pdb"
        # Load the Protein's structure
        self.protein = PyProtein (fileName, folder)
        # Parameterize the protein
        parameterize_protein (self.protein, parameters)
        # Calculate it's dihedral angles
        self.protein.calculate_dihedrals()
        # Create the dictionary of loop indexes
        self.loopIndexes = {}
        for value in LoopInfo:
            self.loopIndexes[value] = LoopInfo[value][0]
        # Properly size the compatibility vector
        cdef size_t I1, I2, J1, J2
        J1 = len(self.loopIndexes)
        self.compatibility.resize(J1)
        for value in self.loopIndexes:
            # If this is the same loop, skip it
            if value == self.loop:
                continue
            # Get the index of that loop
            I1 = self.loopIndexes[value]
            # Get the number of structures in that loop
            J2 = LoopInfo[value][1]
            # Resize the vector
            self.compatibility[I1].resize(J2)
            # Set each value to False
            for I2 in range(J2):
                self.compatibility[I1][I2] = False
        # The compatibility information is stored in this file
        fileName = DBFolder + self.loop \
                 + "/Compatible/" + instructions["database"] + "_" \
                 + self.loop + "_" + str(self.number) + ".txt"
        # Open the file
        f = open(fileName, "r")
        # Define certain variables
        cdef str line, LN
        cdef list items
        # Go through the lines of the file
        for line in f:
            # Split it into items
            items = line.split()
            # If there are exactly 2 (and there should be)
            if len(items) == 2:
                I1 = self.loopIndexes[items[0]]
                I2 = int(items[1]) - 1
                # Indicate that this Structure is compatible with that
                # structure
                self.compatibility[I1][I2] = True

    cpdef bool compatible (self, str L, size_t N):
        """Determine whether or not this structre is compatible with another"""
        # Get the loop index of the specified loop
        cdef size_t I = self.loopIndexes[L]
        # Return the appropriate value from the compatibility vector
        return self.compatibility[I][N-1]

cdef class Loop:
    """A container of information about a Binding Loop"""

    def __init__ (self, label, instructions, DBFolder, LoopInfo, parameters):
        """Initialize a Loop object"""
        # Set the name of the loop
        self.name = label
        # Get the number of structures that this loop has
        self.size = LoopInfo[self.name][1]
        # Set the structures up as an empty list
        self.structures = []
        # Store each structure
        cdef size_t I
        for I in range(1, self.size+1):
            self.structures.append(LoopStructure(DBFolder, self.name, I, 
                                   instructions, LoopInfo, parameters))

    def __getitem__ (self, I):
        """Access to a particular structure"""
        return self.structures[I]

    def __len__ (self):
        """The number of structures in the Loop"""
        return self.size

cdef class InteractionType:
    """The collection of all Interactions of a certain type"""

    def __init__ (self, instructions, N, DBFolder):
        """Initialize an Interaction Type object"""
        # Set the number of the interaction
        self.number = N
        # The distance threshold
        self.threshold = instructions["threshold"]
        # The file containing the interactions of this type
        fileName = DBFolder + "Interactions/" \
                 + instructions["database"] + "_Interactions_" \
                 + str(self.number) + ".txt"
        # Open the file
        f = open(fileName, "r")
        # Interactions are loaded from a string
        cdef string text
        cdef str line
        # Go through the contents of the file
        for line in f:
            # Store it as a C++ string
            text = line
            # Store an interaction
            self.interactions.push_back(Interaction(self.number, text))

    cdef void get_atoms (self, PyResidue residue, list atomNames, 
                         bool errorFlag):
        """Get the specified Atoms from the Residue"""
        # Clear the atoms vector
        self.atoms.clear()
        # Looping variables
        cdef size_t I1, I2, J1, J2
        J1 = len(atomNames)
        J2 = residue._thisptr.size()
        # A flag to indicate whether or not something was found
        cdef bool found
        # A string of the atom name
        cdef string name
        # An Atom pointer
        cdef Atom * atom
        # Go through the atom names
        for I1 in range(J1):
            # Get the name of the atom
            name = atomNames[I1]
            # Indicate that this atom has not yet been found
            found = False
            # Go through the Atoms in the Residue
            for I2 in range(J2):
                atom = residue._thisptr[0][I2]
                # If it is the proper atom
                if atom.name() == name:
                    # Indicate that the atom was found
                    found = True
                    # Store the atom
                    self.atoms.push_back(atom)
                    # Stop the search for this atom
                    break
            # If the atom was not identified
            if not found:
                # If it is OK for some of the atoms to be missing, just end
                # the function
                if not errorFlag:
                    return
                # Otherwise, raise an error
                else:
                    error = "Failure to find each of these atoms in this "
                    error += "Residue: " + str(atomNames) + "\n" +str(residue)
                    raise BindingProteinError (error)

    cdef void evaluate (self, vector[Position]& positions):
        """Evaluate a set of Atoms and store the identified Positions"""
        # Looping variables
        cdef size_t I, J
        J = self.interactions.size()
        # Go through the interactions
        for I in range(J):
            # If the interaction is acceptable
            if self.interactions[I].evaluate(self.atoms, self.threshold):
                # Store a Position
                positions.push_back(Position(&(self.interactions[I]),
                                             self.atoms, self.res,
                                             self.ProtChain))

    cdef void F1 (self, PyResidue residue, vector[Position]& positions):
        """Calculations for interaction 1"""
        # Set the residue string and protein chain
        self.res = residue.get_label()
        self.ProtChain = residue._thisptr.protein()
        # Try to get the backbone O=C atoms. It is OK if that fails - terminal
        # residues won't have an O
        self.get_atoms (residue, ['O','C'], False)
        # If there are exactly 2 of them
        if self.atoms.size() == 2:
            self.evaluate(positions)
        # The lists of double-bonded oxygens in amino acid side chains
        cdef list data
        cdef list group
        # The name of the residue
        cdef str name = residue.name()
        if name == "ASN":
            data = [['OD1', 'CG']]
        elif name == "ASP":
            data = [['OD1', 'CG'], ['OD2', 'CG']]
        elif name == "GLN":
            data = [['OE1', 'CD']]
        elif name == "GLU":
            data = [['OE1', 'CD'], ['OE2', 'CD']]
        else:
            return
        for group in data:
            self.get_atoms(residue, group, True)
            self.evaluate(positions)

    cdef void F2 (self, PyResidue residue, vector[Position]& positions):
        """Calculations for interaction 2"""
        # Set the residue string
        self.res = residue.get_label ()
        self.ProtChain = residue._thisptr.protein()
        # the lists of atoms
        cdef list data, group
        # The name of the residue
        cdef str name = residue.name()
        if name == "SER":
            data = [['OG', 'CB']]
        elif name == "THR":
            data = [['OG1', 'CB']]
        elif name == "TYR":
            data = [['OH', 'CZ']]
        else:
            return
        for group in data:
            self.get_atoms (residue, group, True)
            self.evaluate(positions)

    cdef void F3 (self, PyResidue residue, vector[Position]& positions):
        """Calculations for interaction 3"""
        return

    cdef void F4 (self, PyResidue residue, vector[Position]& positions):
        """Calculations for interaction 4"""
        return

    cdef void F5 (self, PyResidue residue, vector[Position]& positions):
        """Calculations for interaction 5"""
        self.F1 (residue, positions)

    cdef void F6 (self, PyResidue residue, vector[Position]& positions):
        """Calculations for interaction 6"""
        self.F2 (residue, positions)

    cdef void F7 (self, PyResidue residue, vector[Position]& positions):
        """Calculations for interaction 7"""
        # This interaction is an acceptance of H-N hydrogen bonds in the
        # epitope. Every residue except proline has at least one of those
        # (unless it is a terminal amino acid)
        self.res = residue.get_label()
        self.ProtChain = residue._thisptr.protein()
        cdef str name = residue.name()
        if name != "PRO":
            self.get_atoms (residue, ['HN', 'N'], False)
            if self.atoms.size() == 2:
                self.evaluate(positions)
        # The lists of H-N atoms in amino acid side chains
        cdef list data, group
        if name == "ARG":
            data = [['HE', 'NE'], ['HH11', 'NH1'], ['HH12', 'NH1'],
                    ['HH22', 'NH2'], ["HH21", "NH2"]]
        elif name == "ASN":
            data = [['HD21', 'ND2'], ['HD22', 'ND2']]
        elif name == "GLN":
            data = [['HE21', 'NE2'], ['HE22', 'NE2']]
        elif name == "HIS" or name == "HSD":
            data = [['HD1', 'ND1']]
        elif name == "LYS":
            data = [['HZ1', 'NZ'], ['HZ2', 'NZ'], ['HZ3', 'NZ']]
        elif name == "TRP":
            data = [['HE1', 'NE1']]
        else:
            return
        for group in data:
            self.get_atoms (residue, group, True)
            self.evaluate(positions)

    cdef void F8 (self, PyResidue residue, vector[Position]& positions):
        """Calculations for interaction 8"""
        self.res = residue.get_label ()
        self.ProtChain = residue._thisptr.protein()
        cdef str name = residue.name()
        cdef list data, group
        if name == "SER":
            data = [['HG1', 'OG']]
        elif name == "THR":
            data = [['HG1', 'OG1']]
        elif name == 'TYR':
            data = [['HH', 'OH']]
        else:
            return
        for group in data:
            self.get_atoms (residue, group, True)
            self.evaluate(positions)

    cdef void atom_swap (self, size_t i1, size_t i2):
        """Swap the pointers of two atoms in the atoms vector"""
        cdef Atom * temp = self.atoms[i1]
        self.atoms[i1] = self.atoms[i2]
        self.atoms[i2] = temp

    cdef void F9 (self, PyResidue residue, vector[Position]& positions):
        """Calculations for interaction 9"""
        # Get the necessary information
        self.res = residue.get_label ()
        self.ProtChain = residue._thisptr.protein()
        cdef str name = residue.name()
        cdef list data, group
        # Interaction 9 only includes acidic amino acids
        if name == "ASP":
            data = [['CG', 'OD1', 'OD2']]
        elif name == "GLU":
            data = [['CD', 'OE1', 'OE2']]
        else:
            return
        # go through the atom sets
        for group in data:
            # Get the atoms
            self.get_atoms (residue, group, True)
            # Evaluate them
            self.evaluate(positions)
            # In this case, the positions of all 3 atoms have to be exactly
            # right, but the oxygens are interchangable. So swap them and
            # check again
            self.atom_swap(1, 2)
            self.evaluate(positions)

    cdef void F10 (self, PyResidue residue, vector[Position]& positions):
        """Calculations for interaction 10"""
        # This is similar to F9, but without the swapped atoms
        self.res = residue.get_label ()
        self.ProtChain = residue._thisptr.protein()
        cdef str name = residue.name()
        cdef list data, group
        if name == "ASP":
            data = [['CG', 'OD1', 'OD2']]
        elif name == "GLU":
            data = [['CD', 'OE1', 'OE2']]
        else:
            return
        for group in data:
            self.get_atoms (residue, group, True)
            self.evaluate(positions)

    cdef void F11 (self, PyResidue residue, vector[Position]& positions):
        """Calculations for interaction 11"""
        # If this is not an ARG residue, skip it
        if residue.name() != "ARG":
            return
        # Set variables
        self.res = residue.get_label()
        self.ProtChain = residue._thisptr.protein()
        cdef list data, group
        data = [['CZ', 'NH1', 'NH2'], ['CZ', 'NE', 'NH1'], ['CZ', 'NE','NH2']]
        for group in data:
            # Get the atoms and evaluate them
            self.get_atoms (residue, group, True)
            self.evaluate(positions)
            # Swap the two nitrogens and repeat the evaluation
            self.atom_swap (1, 2)
            self.evaluate(positions)

    cdef void F12 (self, PyResidue residue, vector[Position]& positions):
        """Calculations for Interaction 12"""
        # This is only for Lysine ersidues
        if residue.name() != "LYS":
            return
        self.res = residue.get_label()
        self.ProtChain = residue._thisptr.protein()
        self.get_atoms (residue, ['NZ', 'CE'], True)
        self.evaluate(positions)

    cdef void F13 (self, PyResidue residue, vector[Position]& positions):
        """Calculations for Interaction 13"""
        # Get the preliminary variables
        self.res = residue.get_label ()
        self.ProtChain = residue._thisptr.protein()
        cdef str name = residue.name()
        cdef list data
        # This is only for aromatic residues
        if name == "PHE":
            data = ['CZ', 'CE1', 'CD1', 'CG', 'CD2', 'CE2']
        elif name == "TRP":
            data = ['CH2', 'CZ3', 'CE3', 'CD2', 'CE2', 'CZ2']
        elif name == "TYR":
            data = ['CZ', 'CE1', 'CD1', 'CG', 'CD2', 'CE2']
        else:
            return
        # Get the atoms and evaluate them
        self.get_atoms (residue, data, True)
        self.evaluate(positions)

    cdef void F14 (self, PyResidue residue, vector[Position]& positions):
        """Calculations for Interaction 14"""
        self.F13 (residue, positions)

    cdef void F15 (self, PyResidue residue, vector[Position]& positions):
        """Calculations for Interaction 15"""
        # This is only for ARG residues
        if residue.name () != "ARG":
            return
        self.res = residue.get_label()
        self.ProtChain = residue._thisptr.protein()
        self.get_atoms (residue, ['CZ', 'NE', 'NH1', 'NH2'], True)
        self.evaluate(positions)

    cdef void F16 (self, PyResidue residue, vector[Position]& positions):
        """Calculations for interaction 16"""
        # This is only for Lysine residues"""
        if residue.name() != "LYS":
            return
        self.res = residue.get_label()
        self.ProtChain = residue._thisptr.protein()
        self.get_atoms (residue, ['NZ', 'CE'], True)
        self.evaluate(positions)

    cdef void F17 (self, PyResidue residue, vector[Position]& positions):
        """Calculations for interaction 17"""
        # Set the preliminary variables
        self.res = residue.get_label()
        self.ProtChain = residue._thisptr.protein()
        cdef str name = residue.name()
        cdef list data
        # This is only for aromatic residues, but it needs to do some
        # additional work compared to the F13 function
        if name == "PHE":
            data = ['CZ', 'CE1', 'CD1', 'CG', 'CD2', 'CE2']
        elif name == "TRP":
            data = ['CH2', 'CZ3', 'CE3', 'CD2', 'CE2', 'CZ2']
        elif name == "TYR":
            data = ['CZ', 'CE1', 'CD1', 'CG', 'CD2', 'CE2']
        else:
            return
        # Get the atoms and evaluate them
        self.get_atoms (residue, data, True)
        self.evaluate(positions)
        # Each of the 6 atoms needs to be checked as "tip" atom
        cdef size_t I, J
        J = 0
        for I in range(1, 6):
            self.atom_swap(J, I)
            self.evaluate(positions)

    cdef void F18 (self, PyResidue residue, vector[Position]& positions):
        """Calculations for Interaction 18"""
        self.F13 (residue, positions)

    cdef void F19 (self, PyResidue residue, vector[Position]& positions):
        """Calculations for Interaction 19"""
        self.F13 (residue, positions)

    cdef void F20 (self, PyResidue residue, vector[Position]& positions):
        """Calculations for Interaction 20"""
        self.F13 (residue, positions)

    cdef void F21 (self, PyResidue residue, vector[Position]& positions):
        """Calculations for Interaction 21"""
        self.F13 (residue, positions)

    cdef void F22 (self, PyResidue residue, vector[Position]& positions):
        """Calculations for Interaction 22"""
        # This is only for histidine residues
        if residue.name() not in ['HIS', 'HSD']:
            return
        # Do the calculations
        self.res = residue.get_label ()
        self.ProtChain = residue._thisptr.protein()
        self.get_atoms (residue, ['CG', 'CD2', 'NE2', 'CE1', 'ND1'], True)
        self.evaluate(positions)

    cdef void F23 (self, PyResidue residue, vector[Position]& positions):
        """Calculations for Interaction 23"""
        # This is only for histidine residues
        if residue.name() not in ['HIS', 'HSD']:
            return
        # Do the calculations for NE2 to point at the aromatic ring
        self.res = residue.get_label ()
        self.ProtChain = residue._thisptr.protein()
        self.get_atoms (residue, ['NE2', 'CE1', 'ND1', 'CG', 'CD2'], True)
        self.evaluate(positions)
        # Do the calculations for CE1 to point at the aromatic ring
        self.atom_swap (0, 1)
        self.evaluate(positions)

    cdef void find_positions (self, vector[Position]& positions, list epitope):
        """Find the positions of this interaction type"""
        # Go through the epitope residues
        cdef size_t I, J
        J = len(epitope)
        cdef PyResidue residue
        cdef str name
        for I in range(J):
            residue = epitope[I]
            # Depending on the interaction's number, call the proper function
            if self.number == 1:
                self.F1(residue, positions)
            elif self.number == 2:
                self.F2 (residue, positions)
            elif self.number == 3:
                self.F3 (residue, positions)
            elif self.number == 4:
                self.F4 (residue, positions)
            elif self.number == 5:
                self.F5 (residue, positions)
            elif self.number == 6:
                self.F6 (residue, positions)
            elif self.number == 7:
                self.F7 (residue, positions)
            elif self.number == 8:
                self.F8 (residue, positions)
            elif self.number == 9:
                self.F9 (residue, positions)
            elif self.number == 10:
                self.F10 (residue, positions)
            elif self.number == 11:
                self.F11 (residue, positions)
            elif self.number == 12:
                self.F12 (residue, positions)
            elif self.number == 13:
                self.F13 (residue, positions)
            elif self.number == 14:
                self.F14 (residue, positions)
            elif self.number == 15:
                self.F15 (residue, positions)
            elif self.number == 16:
                self.F16 (residue, positions)
            elif self.number == 17:
                self.F17 (residue, positions)
            elif self.number == 18:
                self.F18 (residue, positions)
            elif self.number == 19:
                self.F19 (residue, positions)
            elif self.number == 20:
                self.F20 (residue, positions)
            elif self.number == 21:
                self.F21 (residue, positions)
            elif self.number == 22:
                self.F22 (residue, positions)
            elif self.number == 23:
                self.F23 (residue, positions)
            else:
                error = "Unrecognized interaction type: "+str(self.number)+"\n"
                raise BindingProteinError (error)

    def __len__ (self):
        """The number of interactions of this type"""
        return self.interactions.size()

cdef class PySolution:
    """An extension type wrapper of the C++ Solution class"""

    def __cinit__ (self):
        """C-level initialization of the PySolution object"""
        # Set the pointer to NULL
        self._thisptr = NULL

    def __dealloc__ (self):
        """Delete the Solution"""
        if self._thisptr != NULL:
            del self._thisptr
            self._thisptr = NULL

    def __init__ (self):
        """Initialize a PySolution object"""
        pass

    cdef void _check_for_null (self) except +:
        """Check to see if the Solution is valid"""
        if self._thisptr == NULL:
            error = "Operations cannot be performed on an unallocated "
            error += "PySolution object\n"
            raise BindingProteinError (error)

    cdef void allocate (self, vector[Position]& positions, 
                        vector[size_t]& spots, float n1, float n2, float n3):
        """Allocate information in the Solution object"""
        # If there is already a solution stored, delete it
        if self._thisptr != NULL:
            del self._thisptr
            self._thisptr = NULL
        # Create a new solution using the provided information
        self._thisptr = new Solution (positions, spots, n1, n2, n3)

    cdef size_t size (self) except +:
        """The number of positions in the solution"""
        self._check_for_null()
        return self._thisptr.size()

    def __len__ (self):
        """The number of positions in the solution"""
        self._check_for_null()
        return self._thisptr.size()

    def __str__ (self):
        """The summary of the Solution's positions"""
        self._check_for_null()
        cdef size_t I, J
        J = self._thisptr.size()
        cdef str output = ""
        for I in range(J):
            output += self._thisptr[0][I].str()
        return output

    cdef Position * position (self, size_t i) except +:
        """Get the pointer to a Position in the solution"""
        self._check_for_null ()
        return self._thisptr[0][i]

    cpdef float X (self):
        """Return the X coordinate of the Solution's movement"""
        self._check_for_null ()
        return self._thisptr.X()

    cpdef float Y (self):
        """Return the Y coordinate of the Solution's movement"""
        self._check_for_null ()
        return self._thisptr.Y()

    cpdef float Z (self):
        """Return the Z coordinate of the Solution's movement"""
        self._check_for_null ()
        return self._thisptr.Z ()

    cpdef float XA (self):
        """Return the angle of rotation around the X axis"""
        self._check_for_null ()
        return self._thisptr.XA()

    cpdef float YA (self):
        """Return the angle of rotation around the Y axis"""
        self._check_for_null ()
        return self._thisptr.YA()

    cpdef float ZA (self):
        """Return the angle of rotation around the Z axis"""
        self._check_for_null ()
        return self._thisptr.ZA ()

    cdef dict list_loops (self):
        """List the loop structures the solution uses"""
        # Store the values here
        cdef dict answers = {}
        # Loop and structure variables
        cdef str loop
        # Go through the Positions
        cdef size_t I, J
        J = self._thisptr.size()
        cdef Position * ptr
        for I in range(J):
            ptr = self._thisptr[0][I]
            # Get the loop and structure information for the position
            loop = ptr.loop()
            # If this is a new loop in the answers
            if loop not in answers:
                answers[loop] = ptr.structure()
        return answers

    cpdef float similarity (self, PySolution other):
        """Calculate the similarity between 2 PySolutions"""
        # Get the binding loop information for each solution
        cdef dict selfLoops = self.list_loops()
        cdef dict otherLoops = other.list_loops()
        # Calculate the similarity between this solution and the other
        # solution
        cdef float N1 = 0.0
        # how much each value is the same
        cdef float increment = 1.0 / float(len(selfLoops))
        # Go through the loops in the self dictionary
        for L in selfLoops:
            # If the loop is in the other dictionary and the structure number
            # is the same
            if L in otherLoops and selfLoops[L] == otherLoops[L]:
                N1 += increment
        # Calculate the similarity between the other solution and this
        # solution
        cdef float N2 = 0.0
        increment = 1.0 / float(len(otherLoops))
        for L in otherLoops:
            if L in selfLoops and otherLoops[L] == selfLoops[L]:
                N2 += increment
        # Return which ever value is larger
        if N1 > N2:
            return N1
        return N2

cdef class SolutionStorage:
    """A class for storing Solution objects"""

    def __init__ (self):
        """Initialize a Solution Storage object"""
        # Set up the class attributes
        self.solutions = {}
        self.order = []

    def store_order (self, size_t N):
        """Modify the storage to include a new position count"""
        # If the value is already in the order, be done
        if N in self.order:
            return
        # Create an empty list of this length in the solution
        self.solutions[N] = []
        # Put it in the right spot in the order
        if len(self.order) == 0:
            self.order.append(N)
        else:
            found = False
            for i in range(len(self.order)):
                if N > self.order[i]:
                    self.order.insert(i, N)
                    found = True
                    break
            if not found:
                self.order.append(N)

    cdef void clean_up (self):
        """Possibly remove identified solutions from the container"""
        # The total number of solutions identified so far
        count = 0
        # Go through the values in the order
        for N in self.order:
            # If the solution dictionary is not None for this value
            if self.solutions[N] != None:
                # If the count is already too large, delete the list and set
                # it to None
                if count >= 100000:
                    self.solutions[N] = None
                # Otherwise, add this list's quantity to the count
                else:
                    count += len(self.solutions[N])

    cdef void store_solutions (self, vector[Position]& positions,
              vector[vector[vector[size_t]]]& findings, 
              float n1, float n2, float n3):
        """Store solutions in the Solution Storage object"""
        # Looping variables
        cdef size_t I1, I2, J1, J2
        # Capacity / size variables
        cdef size_t N
        # A solution object
        cdef PySolution answer
        # Go through the findings
        J1 = findings.size()
        for I1 in range(J1):
            # Get the number of entries in this count number
            J2 = findings[I1].size()
            # If there aren't any, skip them
            if J2 == 0:
                continue
            # Get the number of positions in these solutions
            N = findings[I1][0].size()
            # If N is not in self.order
            if N not in self.order:
                self.store_order(N)
            # If the entry for N is None, skip this size
            if self.solutions[N] == None:
                continue
            # Go through the solutions of this many positions
            for I2 in range(J2):
                # Create a solution
                answer = PySolution ()
                # Allocate information into the answer
                answer.allocate(positions, findings[I1][I2], n1, n2, n3)
                # Store the solution
                self.solutions[N].append(answer)
        # Clean up information in the solution if needed
        self.clean_up()

    def __len__ (self):
        """The number of solutions identified"""
        count = 0
        for N in self.order:
            if self.solutions[N] != None:
                count += len(self.solutions[N])
        return count

    def __str__ (self):
        """A string summary of the identified solutions"""
        message = ""
        for N in self.order:
            if self.solutions[N] != None:
                message += "Solutions of " + str(N) + " Positions: " \
                         + str(len(self.solutions[N])) + "\n"
        message += str(self.__len__()) + " total Solutions found\n"
        return message

