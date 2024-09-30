# Written in 2018 by Dr. Robert Pantazes. Copyright owned by Auburn
# University.
# This file contains the implementation of the PyResidue class. This class is
# a container of Atoms in a Protein.

# Cython compiler directives
# cython: c_string_type=str
# cython: c_string_encoding=ascii
# cython: language_leve=3
# distutils: language=c++

# Include Python modules
import os
import sys
import math

# Include the error class 
from Proteins.ProteinError import ProteinError

# Include C++ classes
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp cimport bool

# Include the C++ and Cython classes
from Proteins.Proteins cimport Atom
from Proteins.Proteins cimport Matrix
from Proteins.Proteins cimport Residue
from Proteins.Proteins cimport Protein
from Proteins.Atom.PyAtom cimport PyAtom
from Proteins.Matrix.PyMatrix cimport PyMatrix
from Proteins.Protein.PyProtein cimport PyProtein

# Declare the extension type
cdef class PyResidue:
    """An extension type for the C++ Residue class."""

    def __cinit__ (self, inputs = None):
        """Do c-level initialization of the class"""
        # Allocate a null pointer
        self._thisptr = NULL
        # Set the delete flag to true
        self.delete_flag = True
        # If there are no inputs, be done
        if inputs == None:
            return
        # If the inputs are not a list, raise an error
        elif not isinstance (inputs, list):
            text = "A PyResidue must be initialized from a list.\n"
            text += type(inputs) + " is not acceptable.\n"
            raise ProteinError (text)
        # If the list is empty, raise an error
        elif len(inputs) == 0:
            text = "A PyResidue cannot be initialized from an empty list.\n"
            raise ProteinError (text)
        # If the first entry in the list is not a string, raise an error
        elif not isinstance(inputs[0], str):
            text = "A PyResidue must initialized from a list of strings.\n"
            text += type(inputs[0]) + " is not acceptable.\n"
            raise ProteinError (text)
        # C-define a vector of Atoms and a vector of Atom pointers
        cdef vector[Atom] atoms
        cdef vector[Atom *] ptrs
        # C-define looping variables
        cdef int i, N
        # C-define a string
        cdef string line
        # Get the number of lines
        N = len(inputs)
        # Reserve enough space in the vectors
        atoms.reserve(N)
        ptrs.reserve(N)
        # Loop through the lines
        for i in range(N):
            # Store the contents of the list as a string
            line = inputs[i]
            # Create and store an Atom
            atoms.push_back(Atom(line))
        # Loop through the Atoms and store them as pointers
        for i in range(N):
            ptrs.push_back(&(atoms[i]))
        # Create a new Residue using this list of Atom pointers
        self._thisptr = new Residue (ptrs)

    def __dealloc__ (self):
        """Deallocate c-level pointer information"""
        if self._thisptr != NULL and self.delete_flag:
            del self._thisptr

    def __init__ (self, inputs = None):
        """Initialize a PyResidue object"""
        # Initialization is done in the c-initialization function, so this
        # function does nothing
        pass

    cdef void _check_for_null (self):
        """Check to see if the pointer is NULL and raise an error if so."""
        if self._thisptr == NULL:
            text = "Operations cannot be performed on a NULL PyResidue.\n"
            raise ProteinError (text)
    
    def __len__ (self):
        """The number of Atoms in the Residue"""
        self._check_for_null ()
        return self._thisptr.size()

    def name (self):
        """The name of the Residue"""
        self._check_for_null()
        return self._thisptr.name()

    def number (self):
        """The Residue's numbering information"""
        self._check_for_null ()
        return self._thisptr.number ()

    def insertion_code (self):
        """The Residue's insertion code"""
        self._check_for_null ()
        cdef char L = self._thisptr.insertion_code()
        cdef string temp
        temp.push_back(L)
        cdef string answer = temp
        return answer

    def get_label (self):
        """A string of the Residue's number and insertion code"""
        self._check_for_null()
        cdef char L = self._thisptr.insertion_code()
        cdef string code
        code.push_back(L)
        cdef str answer = str(self._thisptr.number())
        answer += code
        return answer.strip()

    def protein (self):
        """The protein the Residue belongs to"""
        self._check_for_null ()
        cdef char L = self._thisptr.protein()
        cdef string temp
        temp.push_back(L)
        cdef string answer = temp
        return answer

    cdef Atom * pointer_by_index (self, int i):
        """Access to a Residue's Atom by index"""
        cdef size_t j
        if (i >= 0):
            j = i
        elif (self._thisptr.size() + i >= 0):
            j = self._thisptr.size() + i
        else:
            j = -i
        return self._thisptr[0][j]

    cdef Atom * pointer_by_name (self, str n):
        """Access to a Residue's Atom by name"""
        cdef string name1 = n
        return self._thisptr[0][name1]
    
    def __getitem__ (self, value):
        """Access to the Atoms of the Residue"""
        self._check_for_null ()
        # The PyAtom that will contain the pointer
        cdef PyAtom atom = PyAtom ()
        # Set its pointer value
        if isinstance(value, int):
            atom._thisptr = self.pointer_by_index(value)
        elif isinstance(value, str):
            atom._thisptr = self.pointer_by_name(value)
        else:
            text = type(value) + " is not a valid data type for accessing "
            text += "Atoms in a PyResidue.\n"
            raise ProteinError (text)
        # Set the delete flag of the atom to false
        atom.delete_flag = False
        return atom

    def set_number (self, long n, str Insert):
        """Set the Residue's numbering information"""
        self._check_for_null()
        # Get a character from the Insert information
        cdef char L = 32
        cdef str temp = Insert
        if len(Insert) > 0:
            L = temp[0]
        cdef bool internal = False
        # Set the Residue's numbering information
        self._thisptr.set_number(n, L, internal)

    def renumber_atoms (self, long n):
        """Renumber the Atoms in the Residue"""
        self._check_for_null ()
        return self._thisptr.renumber_atoms(n)

    def last_atom_number (self):
        """The number of the last atom in the Residue"""
        self._check_for_null ()
        return self._thisptr.last_atom_number()

    def format (self, bool internal = False):
        """Create a string of formatted text of the Residue's Atoms"""
        self._check_for_null()
        cdef bool how = internal
        return self._thisptr.str(how)

    def __str__ (self):
        """Create a string of formatted text of the Residue's Atoms"""
        return self.format(False)

    def move (self, PyMatrix mat, bool subtract = True):
        """Move the Residue using the matrix"""
        self._check_for_null ()
        self._thisptr.move(mat._thisptr, subtract)

    def rotate (self, PyMatrix mat):
        """Rotate the Residue"""
        self._check_for_null ()
        self._thisptr.rotate(mat._thisptr)

    def center (self):
        """Move the Residue so that its center of mass is at the origin"""
        self._check_for_null ()
        cdef PyMatrix mat = PyMatrix ()
        cdef Matrix result = self._thisptr.center()
        mat._thisptr.allocate(&result)
        return mat

    def position (self, PyMatrix mat, bool forward = True):
        """Move the Residue into a starndard position"""
        self._check_for_null ()
        self._thisptr.position(mat._thisptr[0], forward)

    def missing_atoms (self):
        """Find out how many missing atoms a residue has"""
        self._check_for_null()
        return self._thisptr.missing_atoms()

    def is_present (self):
        """Whether or not the Residue is present in structural data"""
        self._check_for_null()
        return self._thisptr.is_present()

    def is_amino_acid (self):
        """Whether or not the Residue is an amino acid"""
        self._check_for_null()
        return self._thisptr.is_amino_acid()

    def phi (self):
        """The Residue's phi dihedral angle"""
        self._check_for_null()
        return self._thisptr.phi()

    def psi (self):
        """The Residue's psi dihedral angle"""
        self._check_for_null()
        return self._thisptr.psi()

    def omega (self):
        """The Residue's omega dihedral angle"""
        self._check_for_null()
        return self._thisptr.omega()

    def charmm_his_fix (self):
        """Correct the names of histidine residues"""
        self._check_for_null()
        self._thisptr.charmm_his_fix()

    def rotamer_library_path (self, start):
        """Generate a path to the proper portion of the rotamer library"""
        # Make sure this rotamer isn't empty
        self._check_for_null()
        # Use a try statement in case either the phi or psi dihedral angle is
        # unavailable
        #try:
        #    # Get the phi angle rounded to the nearest 10
        #    f = format(int(round(self._thisptr.phi(), -1)), "+") 
        #    s = format(int(round(self._thisptr.psi(), -1)), "+")
        #    # Make the path
        #    path = start + f + s + "/"
        #    return path
        # If there is an error, reference the folder for the ends of a protein
        #except (RuntimeError, ValueError):
        #    path = start + "independ/"
        #    return path
        return start

    def collect_rotamers (self, rotLibPath, which = None, specific = None):
        """Collect and position the rotamers for this residue"""
        # Get the path to the proper portion of the rotamer library. That
        # function will check to make sure that this is not a null pointer
        path = self.rotamer_library_path (rotLibPath)
        # Make a list of the amino acids that should be used
        if which == None:
            tag = self.name()
            if tag == "HIS":
                tag = "HSD"
            which = [tag]
        # Generate the matrix for positioning rotamers into this protein
        here = PyMatrix()
        self.position(here, True)
        self.position(here, False)
        # Store the residues here
        results = []
        # Loop through the amino acids
        for amino in which:
            # If this is a proline, there is no Proline in the rotamer
            # library. In that case, use the one that is hard coded here. You
            # can refer to the IPRO Suite paper for details of how it was
            # identified. The Ramachandran values come from data in S.C.
            # Lovell, et al., "Structure validation by Ca geometry: phi, psi
            # and Cb deviation," Proteins, 50:437-450 (2003). The residue
            # itself is residue 14 in chain B of 1CFQ.pdb
            if amino == "PRO":
                # Whether or not a rotamer or proline can be collected
                use = False
                # Try to make a decision based on dihedral angles
                try:
                    f = int(round(self._thisptr.phi(), -1))
                    s = int(round(self._thisptr.psi(), -1))
                    if f == -40 and ((-60 <= s <= -30) or (110 <= s <= 120)):
                        use = True
                    elif f == -50 and ((-60 <= s <= -30) or (110 <= s <= 150)):
                        use = True
                    elif f == -60 and ((-60 <= s <= 0) or (100 <= s <= 170)):
                        use = True
                    elif f == -70 and -40 <= s <= 170:
                        use = True
                    elif f == -80 and ((-30 <= s <= 180) or (-180 <= s <= -170)):
                        use = True
                    elif f == -90 and (s == -180 or (-20 <= s <= 30) or \
                            (120 <= s <= 180)):
                        use = True
                    elif f == -100 and ((-10 <= s <= 20) or (120 <= s <= 170)):
                        use = True
                # If that fails because the dihedral angles aren't available,
                # permit the rotamer
                except RuntimeError:
                    use = True
                # Also permit the rotamer if the residue is currently an amino
                # acid
                if self.name() == "PRO":
                    use = True
                if len(which) == 1:
                    use = True
                # If the rotamer should not be used
                if not use:
                    continue
                # Since the rotamer should be used
                text = """
ATOM      1  N   PRO P   1       0.000   0.000   1.451  1.00 52.87      MLP
ATOM      2  CD  PRO P   1       1.380  -0.040   1.994  1.00 52.80      MLP
ATOM      3  HD1 PRO P   1       0.768   0.827   2.331  1.00  0.00      MLP
ATOM      4  HD2 PRO P   1       2.021  -0.413   2.827  1.00  0.00      MLP
ATOM      5  CA  PRO P   1       0.000   0.000   0.000  1.00 51.90      MLP
ATOM      6  HA  PRO P   1      -0.455   0.913  -0.367  1.00  0.00      MLP
ATOM      7  CB  PRO P   1       1.468   0.000  -0.396  1.00 52.39      MLP
ATOM      8  HB1 PRO P   1       1.843   1.048  -0.400  1.00  0.00      MLP
ATOM      9  HB2 PRO P   1       1.643  -0.449  -1.394  1.00  0.00      MLP
ATOM     10  CG  PRO P   1       2.236   0.357   0.825  1.00 52.57      MLP
ATOM     11  HG1 PRO P   1       3.233   0.784   0.599  1.00  0.00      MLP
ATOM     12  HG2 PRO P   1       2.357  -0.550   1.460  1.00  0.00      MLP
ATOM     13  C   PRO P   1      -0.731  -1.192  -0.590  1.00 50.77      MLP
ATOM     14  O   PRO P   1      -0.717  -2.296  -0.043  1.00 50.66      MLP"""
                # Split the text into lines
                lines = text.split("\n")[1:]
                # Store a Residue and put it into the protein's backbone. It
                # does not need to have the original positioning done because
                # it is already in that state
                results.append(PyResidue(lines))
                results[-1].position(here, False)
            # If it is any other amino acid
            else:
                # If a specific rotamer is needed
                if specific != None:
                    I1 = specific
                    I2 = I1+1
                else:
                    I1 = 1
                    I2 = 1001
                # Loop through rotamer numbers until a file is not found
                for INDEX in range(I1, I2):
                    # The name of the file
                    fileName = path + amino.lower() + str(INDEX) + ".pdb"
                    # Open it
                    try:
                        f = open(fileName, "r")
                    # If that fails, break the search
                    except IOError:
                        break
                    # Store the Atom lines here
                    atoms = []
                    for line in f:
                        if line.startswith("ATOM"):
                            atoms.append(line.strip())
                    # Close the file
                    f.close()
                    # Store a PyResidue
                    results.append(PyResidue(atoms))
                    # And properly position it
                    results[-1].position(PyMatrix(), True)
                    results[-1].position(here, False)
        # Return the residues list
        return results

    def closest_rotamer (self, rotLibPath):
        """Identify the most similar rotamer to this residue"""
        # Make sure this is a valid residue
        self._check_for_null()
        # Collect the rotamers for this residue
        rotamers = self.collect_rotamers(rotLibPath)
        # Determine which is most similar to this residue
        cdef size_t bestIndex = len(rotamers) + 1
        cdef float bestRMSD = 100000000.0
        cdef float RMSD
        # A PyResidue
        cdef PyResidue rotamer
        # Two Atom pointers
        cdef Atom * atom1
        cdef Atom * atom2
        # A true boolean value
        cdef bool how = True
        # Loop through the rotamers
        cdef size_t i, j, k, l, m, n
        cdef float count
        l = bestIndex - 1
        m = self._thisptr.size()
        for i in range(l):
            # Get the PyResidue
            rotamer = rotamers[i]
            # The number of Atoms in the rotamer
            n = rotamer._thisptr.size()
            # Set the RMSD for this rotamer to 0
            RMSD = 0
            # The number of atoms used starts at 0
            count = 0
            # Loop through the Residue's Atoms
            for j in range(m):
                # Get a pointer to the atom
                atom1 = self._thisptr[0][j]
                # find the corresponding atom in the other residue
                for k in range(n):
                    atom2 = rotamer._thisptr[0][k]
                    # If they have the same name
                    if atom1.name() == atom2.name():
                        # Increment the counter
                        count += 1
                        # Calculate the square of the distance between the
                        # atoms
                        RMSD += atom1.calculate_distance(atom2, how)
                        # Stop searching the atoms in the rotamer
                        break
            # Calculate the final RMSD
            if count > 0:
                RMSD = math.sqrt(RMSD / count)
                # If this is the best one found so far, store that information
                if RMSD < bestRMSD:
                    bestRMSD = RMSD
                    bestIndex = i
        # If no best rotamer was identified
        if bestIndex > len(rotamers):
            error = "Unexpected failure in the closest rotamer function.\n"
            raise ProteinError (error)
        # Make a vector of the best rotamer's atoms
        cdef vector[Atom *] use
        rotamer = rotamers[bestIndex]
        use.reserve(rotamer._thisptr.size())
        # Loop through said atoms
        for i in range(rotamer._thisptr.size()):
            # Store the rotamer's atom's pointer
            use.push_back(rotamer._thisptr[0][i])
        # Use those atom's in this residue's load function
        cdef bool sidechainOnly = True
        cdef bool complete = False
        self._thisptr.load(use, sidechainOnly, complete)

    def specific_rotamer (self, rotLibPath, I):
        """Patch a specific rotamer into this one"""
        # Make sure the pointer is valid
        self._check_for_null()
        # Collect the specific rotamer that is needed
        rotamers = self.collect_rotamers(rotLibPath, None, I)
        # Make a vector of Atom pointers
        cdef vector[Atom *] use
        # Get the specific rotamer
        cdef PyResidue rotamer = rotamers[0]
        use.reserve(rotamer._thisptr.size())
        # Store the atom pointers
        cdef size_t i1, j1
        j1 = rotamer._thisptr.size()
        for i1 in range(j1):
            use.push_back(rotamer._thisptr[0][i1])
        # Use those atoms in the residue
        cdef bool sidechainOnly = True
        cdef bool complete = False
        self._thisptr.load(use, sidechainOnly, complete)

    def load (self, PyResidue other):
        """Patch the Atoms from another PyResidue into this residue"""
        # make sure both residues are usable
        self._check_for_null()
        other._check_for_null()
        # Make a vector of atom pointers from the other residue
        cdef vector[Atom *] use
        cdef size_t i, j
        j = other._thisptr.size()
        use.reserve(j)
        for i in range(j):
            use.push_back(other._thisptr[0][i])
        # Only the side chain is being patched in and it is not a complete
        # loading of all residue information
        cdef bool sidechainOnly = True
        cdef bool complete = False
        # Load the information
        self._thisptr.load(use, sidechainOnly, complete)

    def parameterize (self, parameters, bool NTER, bool CTER, bool VDW, 
                      bool ELEC, bool LK, bool GB = False, bool VDW_H = True):
        """Assign energy parameters to the Residue's Atoms."""
        # Make sure that the Residue isn't empty
        self._check_for_null()
        # The Residue's Atoms will be worked with through this PyAtom
        cdef PyAtom atom = PyAtom ()
        # Set it's delete flag to 0 so the Atom isn't deleted when the
        # function is finished
        atom.delete_flag = False
        # A string for the Atom's name
        cdef str atomName
        # A string of the Residue's name
        cdef str resName = self._thisptr.name()
        # Loop through the Residue's Atoms
        cdef size_t i, j
        j = self._thisptr.size()
        for i in range(j):
            # Store this C++ Atom in the PyAtom so the parameterize function
            # can be used
            atom._thisptr = self._thisptr[0][i]
            # Get the Atom's name
            atomName = atom._thisptr.name()
            # Determine the proper residue label to get information from the
            # parameters dictionary
            if NTER and resName == "GLY" and atomName in ['N', 'HT1', 'HT2', \
                    'HT3', 'CA', 'HA1', 'HA2']:
                res = "GLYP"
            elif NTER and resName == "PRO" and atomName in ['N', 'HN1', 'HN2',\
                    'CA', 'HA', 'CD', 'HD1', 'HD2']:
                res = "PROP"
            elif NTER and (self._thisptr.is_amino_acid() or resName == "HSD") \
                    and atomName in ['N', 'HT1', 'HT2', 'HT3', 'CA', 'HA']:
                res = "NTER"
            elif CTER and (self._thisptr.is_amino_acid() or resName == "HSD") \
                    and atomName in ['C', 'OT1', 'OT2']:
                res = "CTER"
            elif resName == "HIS":
                res = "HSD"
            else:
                res = resName
            # Parameterize the Atom
            try:
                atom.parameterize(parameters[res][atomName], VDW, ELEC, LK,
                                  GB, VDW_H)
            # If the information isn't in the parameters dictionary, raise an
            # error
            except KeyError:
                error = "Missing parameter information for Atom " + atomName
                error += " in Residue " + res + ".\n"
                raise ProteinError (error)

    cpdef size_t calculate_atom_overlap (self, PyAtom atom, 
                                         bool rotamer = False):
        """Count the VDW overlaps between the Residue and an atom"""
        # Validate the pointers
        self._check_for_null()
        atom._check_for_null_pointer()
        # Use the C++ function
        return self._thisptr.calculate_vdw_overlap (atom._thisptr, rotamer)

    cpdef size_t calculate_residue_overlap (self, PyResidue other, 
                 bool selfRotamer = False, bool otherRotamer = False):
        """Count the VDW overlaps between 2 Residues"""
        # Validate the pointers
        self._check_for_null()
        other._check_for_null()
        # Use the C++ function
        return self._thisptr.calculate_vdw_overlap (other._thisptr,
                             selfRotamer, otherRotamer)

    cpdef float calculate_atom_energy (self, PyAtom atom, 
                                       bool rotamer = False):
        """Calculate the energy between the residue and an atom"""
        # Validate the pointers
        self._check_for_null()
        atom._check_for_null_pointer()
        # Use the C++ function
        return self._thisptr.calculate_energy(atom._thisptr, rotamer)

    cpdef float calculate_residue_energy (self, PyResidue other,
                bool selfRotamer = False, bool otherRotamer = False):
        """Calculate the energy between two residues"""
        # Validate the pionters
        self._check_for_null()
        other._check_for_null()
        # Use the C++ function
        return self._thisptr.calculate_energy(other._thisptr, selfRotamer,
                                              otherRotamer)

    def disable_VDW (self):
        """Turn off VDW calculations for the Residue's atoms"""
        # Make sure the Residue is valid
        self._check_for_null()
        # The boolean off value
        cdef bool off = False
        # Looping variables
        cdef size_t i, j
        # An Atom pointer
        cdef Atom * atom
        # Loop through the Residue's atoms
        j = self._thisptr.size()
        for i in range(j):
            atom = self._thisptr[0][i]
            atom.set_calculate_vdw(off)

    def duplicate (self):
        """Create another PyResidue that matches this one"""
        # Make sure the residue is valid
        self._check_for_null()
        # Create an empty PyResidue
        cdef PyResidue other = PyResidue ()
        # Allocate it's Pointer using the duplicate method of this Residue
        other._thisptr = new Residue (self._thisptr.duplicate())
        # Return the residue
        return other
