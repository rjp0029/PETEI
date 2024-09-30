# Written in 2018 by Dr. Robert Pantazes. Copyright owned by Auburn
# University.
# This file contains the classes and functions necessary to search a
# collection of Protein Data Bank files for pieces that are structurally
# compatible with the framework of the database of parts being generated.

# Cython compiler directives
# cython: c_string_type=str
# cython: c_string_encoding=ascii
# cython: language_level=3
# distutils: language=c++

# Include pure Python modules
import os
import sys
# Include Binding Protein Design functions
from BPD.BindingProteinError import BindingProteinError
from BPD.General import parameterize_protein
from BPD.General import parameterize_residue
from BPD.Instructions import time_stamp
from BPD.Database.Preliminary import database_folder_name
from BPD.Database.Preliminary import summary_file_name
# Include the function that loads energy calculation parameters
from ForceFields.Parameters import load_parameters
# Include the Py-wrappers of the C++ PDB structure classes
from Proteins.PDB.PyPDB cimport PyPDB
from Proteins.Protein.PyProtein cimport PyProtein
from Proteins.Residue.PyResidue cimport PyResidue
from Proteins.Matrix.PyMatrix cimport PyMatrix
from Proteins.Atom.PyAtom cimport PyAtom
# Also include the C++ classes
from Proteins.Proteins cimport Protein
from Proteins.Proteins cimport Residue
from Proteins.Proteins cimport Matrix
from Proteins.Proteins cimport Atom
# Include some useful C++ classes
from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector

# Create a C++ structure for a distance bound between atoms in a binding loop.
# It stores the average expected distances, an Upper Bound and a Lower Bound
cdef struct Bound:
    float average
    float UB
    float LB

cdef float bound_compatibility (Bound& bound, Atom * atom1, Atom * atom2):
    """Calculate whether or not a pair of atoms is compatible with a Bound"""
    # A boolean value is needed for the distance calculation
    cdef bool squared = False
    # The distance between the atoms
    cdef float D = atom1.calculate_distance(atom2, squared)
    # If the distance is not acceptable, return a negative value to indicate
    # this
    if D < bound.LB or D > bound.UB:
        return -1.0
    # Otherwise, return the square of the average minus the distance
    return (bound.average - D)**2

cdef vector[vector[string]] get_atom_names (dict instructions):
    """Get the lists of atom names that need distance calculations"""
    # Store the names here
    cdef vector[vector[string]] atomNames
    # Set it up to contain a proper number of values
    cdef size_t N = len(instructions["atoms"])
    atomNames.resize(2)
    atomNames[0].reserve(N)
    atomNames[1].reserve(N)
    # Have strings for converting the Python values to c++ strings
    cdef string name1, name2
    # Go through the information in the instructions
    for pair in instructions["atoms"]:
        name1 = pair[0]
        name2 = pair[1]
        atomNames[0].push_back(name1)
        atomNames[1].push_back(name2)
    # Return the atom names
    return atomNames

cdef bool get_atom_pointers(Residue * res, vector[string]& labels, 
                            vector[Atom *]& atoms):
    """Retrieve the pointers to the relevant Atoms from the Residue"""
    # Set the vector up to hold the proper number of atoms
    atoms.clear()
    atoms.reserve(labels.size())
    # Have looping variables
    cdef size_t I1, I2, J1, J2
    J1 = labels.size()
    J2 = res.size()
    # Have variables for storing pointers and strings
    cdef Atom * atom
    cdef string label, name
    # Special care needs to be taken for handling glycines
    cdef str temp = "GLY"
    cdef string RESNAME = temp
    temp = "CB"
    cdef string ATOMNAME = temp
    temp = "HA1"
    cdef string REPLACE = temp
    # Have a boolean value to indicate whether or not the function is
    # successful
    cdef bool found
    # Go through the labels
    for I1 in range(J1):
        # Store the label
        label = labels[I1]
        # Possibly update it if the residue is Glycine
        if res.name() == RESNAME and label == ATOMNAME:
            label = REPLACE
        # Set the flag to False
        found = False
        # Go through the residue atoms
        for I2 in range(J2):
            # Get the atom and it's name
            atom = res[0][I2]
            name = atom.name()
            # If the name is right, store the atom pointer
            if name == label:
                atoms.push_back(atom)
                found = True
                break
        # If the atom wasn't found, the function can auto-return
        if not found:
            break
    return found

cdef list make_PyAtoms (vector[Atom *]& atoms):
    """Make a list of PyAtoms from a vector of Atom Pointers"""
    # Store the PyAtoms here
    cdef list results = []
    cdef PyAtom atom
    # Looping variables
    cdef size_t I1, J1
    J1 = atoms.size()
    # Go through the Atoms
    for I1 in range(J1):
        atom = PyAtom ()
        atom._thisptr = atoms[I1]
        atom.delete_flag = False
        results.append(atom)
    # Return the list of PyAtoms
    return results

cdef class Loop:
    """A container of information about a binding loop"""

    # The attributes of the class
    cdef vector[Bound] bounds
    cdef list atoms
    cdef PyMatrix movement
    cdef vector[string] names
    cdef vector[float] scores
    cdef int limit
    cdef readonly str folder, name
    cdef public size_t count

    def __init__ (self, label, instructions, data, model):
        """Initialize a Loop object"""
        # Initialize the lists and matrix
        self.atoms = []
        self.movement = PyMatrix ()
        # Store some information directly
        self.name = label
        self.limit = instructions["structures"]
        # Store appropriate space in the vectors
        self.names.reserve(self.limit+1)
        self.scores.reserve(self.limit+1)
        # Create the folder for the loop
        self._make_folder(instructions)
        # Get the distance bound information
        self._get_bounds (instructions, data)
        # Get the atoms that are needed from the model for this binding loop
        self._get_atoms (model, instructions)
        # Move those atoms so their center of mass is at the origin
        self.movement.PyAtoms_allocation_for_centering (self.atoms)
        for I, atom in enumerate(self.atoms):
            # Atoms can be duplicated in the lists, but we don't want to move
            # them twice (took me a long time to track down this error). The
            # following steps makes sure that is done correctly
            unique = True
            if I > 0:
                for J in range(I):
                    atom2 = self.atoms[J]
                    if atom2.name() == atom.name() and \
                       atom2.residue() == atom.residue() and \
                       atom2.residue_number() == atom.residue_number():
                           unique = False
                           break
            # If the atom is unique, move it
            if unique:
                atom.move(self.movement)

    def _make_folder (self, instructions):
        """Make the folder for the binding loop"""
        self.folder = database_folder_name(instructions) + self.name + "/"
        try:
            os.mkdir(self.folder)
        except OSError:
            error = "Failure to create the " + self.name + " binding loop "
            error += "folder. It likely already exists.\n"
            raise BindingProteinError (error)
        # The folder should actually lead to the PDB folder for the loop
        self.folder += "PDB/"
        os.mkdir(self.folder)

    def _get_bounds (self, instructions, data):
        """Gather the Bound information for the binding loop"""
        # Get the number of bounds that are expected
        cdef size_t N = len(instructions["atoms"])
        # reserve that much memory
        self.bounds.reserve(N)
        # The distance variables
        cdef float n1, n2
        # A flag to indicate if things weren't found
        flag = False
        # A Bound structure
        cdef Bound bound
        # Go through the lines in the data
        for i, line in enumerate(data):
            # If it is the proper line
            if line.startswith("Binding Loop " + self.name):
                # The flag is True
                flag = True
                # Go through the distance information
                for j in range(1, N + 1):
                    # Split the line on whitespace
                    pieces = data[i+j].split()
                    # Get the average and standard deviation
                    n1 = float(pieces[-3])
                    n2 = float(pieces[-1])
                    # n2 should be a standard deviation, but if only a single
                    # structure was used it will have a value of 0. Use 0.5 as
                    # the default in those cases
                    if n2 < 1e-6:
                        n2 = 0.225
                    # Create the bound
                    bound = Bound (n1, n1+n2, n1-n2)
                    # Store the bound
                    self.bounds.push_back(bound)
                # Stop the search once the values are found
                break
        # Raise an error if there was a problem
        if not flag:
            error = "Failure to identify the distance parameters for binding "
            error += "loop " + self.name + "\n"
            raise BindingProteinError (error)

    def _get_atoms (self, model, instructions):
        """Get the atoms needed for positioning the binding loop"""
        # Find the information about this binding loop in the instructions
        use = None
        for loop in instructions["loops"]:
            if loop[0] == self.name:
                use = loop
                break
        # If it wasn't found, raise an error
        if use == None:
            error = "Failure to identify the " + self.name + " binding loop "
            error += "instructions.\n"
            raise BindingProteinError (error)
        # Store the relevant residues here
        cdef PyResidue startRes
        cdef PyResidue endRes
        # Go through the model proteins
        for i in range(model.proteins()):
            protein = model.protein(i)
            if protein.name() != use[1]:
                continue
            for j in range(len(protein)):
                res = protein[j]
                if res.get_label() == use[2]:
                    startRes = res
                elif res.get_label() == use[3]:
                    endRes = res
                    break
            break
        # Get the lists of needed atom names
        cdef vector[vector[string]] atomNames = get_atom_names(instructions)
        # Have a vector of the needed atoms
        cdef vector[Atom *] atoms
        # Get them for the first residue
        cdef bool found = get_atom_pointers(startRes._thisptr, atomNames[0],
                                            atoms)
        # If they weren't found raise an error
        if not found:
            error = "Failure to find all model atoms in binding loop "
            error += self.name + "\n"
            raise BindingProteinError (error)
        # Store the PyAtoms
        self.atoms.extend(make_PyAtoms(atoms))
        # Do that again for the ending residue
        found = get_atom_pointers(endRes._thisptr, atomNames[1], atoms)
        if not found:
            error = "Failure to find all model atoms in binding loop "
            error += self.name + "\n"
            raise BindingProteinError (error)
        self.atoms.extend(make_PyAtoms(atoms))
        # Make sure there are the proper number of atoms
        if len(self.atoms) != 2*len(instructions["atoms"]):
            error = "Failure to find all model atoms in binding loop "
            error += self.name + "\n"
            raise BindingProteinError (error)

    cdef float compatible (self, vector[Atom *]& atoms1, 
                                 vector[Atom *]& atoms2):
        """Determine if a set of atoms are compatible with the binding loop"""
        # Looping variables
        cdef size_t I1, J1
        J1 = self.bounds.size()
        # An overall performance
        cdef float overall = 0.0
        # And a score for a specific atom pair
        cdef float score 
        # Go through the atom pairs
        for I1 in range(J1):
            score = bound_compatibility(self.bounds[I1], atoms1[I1],
                                        atoms2[I1])
            # If the score is negative, the atoms are not compatible with the
            # binding loop and a negative value should indicate this
            if score < 0:
                return score
            # Otherwise, add the score to the overall comparison
            overall += score
        # Return the overall comparison value
        return overall

    cdef int find_index (self, float score):
        """Find the index for where this result goes in the Loop's results"""
        # Looping variables
        cdef int i, j
        j = self.scores.size()
        # If there are currently values
        if j > 0:
            # If there are already too many scores and this score is bad,
            # don't do this search
            if j == self.limit and score > self.scores[j-1]:
                return -1
            # Go through the current values
            for i in range(j):
                if score < self.scores[i]:
                    return i
        return j

    cdef list prepare_residues (self, PyProtein protein, size_t I1, size_t I2,
                                vector[vector[string]]& atomNames, parameters,
                                str rotlib):
        """Position the possible structure where the binding loop goes"""
        # Store the identified residues here
        cdef list residues = []
        # We need a PyResidue object
        cdef PyResidue residue
        # A looping variable
        cdef size_t i
        # The list of needed PyAtoms
        cdef list atoms = []
        # Variables for getting those atoms
        cdef vector[Atom *] atomPtrs
        cdef bool flag
        # Go through the protein's residues
        for i in range(I1, I2+1):
            # Get a duplicated copy of the residue
            residue = protein[i].duplicate()
            # Patch in it's closest rotamer
            residue.closest_rotamer(rotlib)
            # Parameterize it
            parameterize_residue(residue, parameters)
            # Fix the names of HIS residues
            residue.charmm_his_fix()
            # Get the atoms that are needed
            if i == I1:
                flag = get_atom_pointers(residue._thisptr, atomNames[0],
                                         atomPtrs)
                atoms.extend(make_PyAtoms(atomPtrs))
            elif i == I2:
                flag = get_atom_pointers(residue._thisptr, atomNames[1],
                                         atomPtrs)
                atoms.extend(make_PyAtoms(atomPtrs))
            # Store the residue
            residues.append(residue)
        # Move the Residues
        cdef PyMatrix matrix = PyMatrix()
        matrix.PyAtoms_allocation_for_centering(atoms)
        for i in range((I2-I1)+1):
            residues[i].move(matrix)
        # Calculate the rotation matrix
        matrix.minimized_RMSD_atoms (self.atoms, atoms)
        # Rotate the residues and put them in the framework
        for i in range((I2-I1)+1):
            residue = residues[i]
            residue.rotate(matrix)
            residue.move(self.movement, False)
        return residues

    cdef void store_result (self, str name, float score, int index):
        """Store information about a result in the binding loop"""
        # Convert the name to a C++ string
        cdef string use = name
        # If the index is equal to the current length of the vectors
        if index == self.names.size():
            self.names.push_back(use)
            self.scores.push_back(score)
            return
        # Otherwise, insert needs to be used
        self.names.insert(self.names.begin()+index, 1, use)
        self.scores.insert(self.scores.begin()+index, 1, score)
        # If there are too many values, delete the last
        if self.names.size() <= self.limit:
            return
        # Delete the last score
        self.scores.pop_back()
        # Get the last value in the names vector
        name = self.names[self.limit]
        # Delete that file
        os.remove(self.folder + name)
        # Delete the entry from the vector
        self.names.pop_back()

    def __len__ (self):
        """The number of identified binding loop structures"""
        return self.names.size()

def load_model (instructions):
    """Load the model structure for the database"""
    # It is in this file in this folder
    fileName = instructions["name"] + "_Model.pdb"
    folder = database_folder_name(instructions) + "framework/"
    # Load the PyPDB object
    model = PyPDB(fileName, folder)
    return model

def load_frameworks (instructions, summaryLines, parameters):
    """Load the framework structures"""
    # Get the number of structures expected
    N = None
    for line in summaryLines:
        if line.startswith("Framework Pieces Identified:"):
            N = int(line.split()[3])
            break
    if N == None:
        error = "Error finding the number of framework pieces\n"
        raise BindingProteinError (error)
    # Store the frameworks here
    frameworks = []
    # Go through the numbers
    for I in range(1, N+1):
        # The structure is in this file / folder
        fileName = instructions["name"] + "_framework_" + str(I) + ".pdb"
        folder = database_folder_name(instructions) + "framework/"
        # Load the PyProtein
        protein = PyProtein(fileName, folder)
        # Parameterize and store it
        parameterize_protein(protein, parameters)
        frameworks.append(protein)
    return frameworks

def prepare_for_search (instructions):
    """Prepare to do a search of PDB files for binding loop structures"""
    # Get the name of the summary file
    summaryFileName = summary_file_name (instructions)
    # Read in the contents of that file
    f = open(summaryFileName, "r")
    summaryLines = f.readlines()
    f.close()
    # Open the summary file for writing with appending at the end of the file
    summary = open(summaryFileName, "a")
    # Write a message that this is starting
    message = "Searching for binding loop structures started on "+time_stamp()
    summary.write(message)
    summary.flush()
    # Load the energy parameters
    parameters = load_parameters (instructions["topology"],
                                  instructions["parameters"],
                                  instructions["solvation"])
    # Load the model and frameworks
    model = load_model (instructions)
    frameworks = load_frameworks (instructions, summaryLines, parameters)
    # Find the order of the binding loops
    loopOrder = []
    for line in summaryLines:
        if line.startswith("Binding Loop "):
            loopOrder.append(line.split()[2])
    # Create the loop objects
    loops = []
    for L in loopOrder:
        loops.append(Loop(L, instructions, summaryLines, model))
    # Return the generated objects
    return summary, parameters, model, frameworks, loops

cdef bool assess_framework_compatibility (list residues, list frameworks):
    """Determine if a binding loop is framework - compatible"""
    # Looping variables
    cdef size_t I1, I2, J1, J2
    J1 = len(residues)
    J2 = len(frameworks)
    # A flag to indicate whether or not the loop is compatible with the
    # framework
    cdef bool result = False
    # A flag to indicate whether or not a residue is the first or last in the
    # list. In that case, it should be treated as a "rotamer" for calculations
    # so that it's backbone atoms aren't misused
    cdef bool rotamer
    # Python wrapper types
    cdef PyResidue residue
    cdef PyProtein framework
    # Count the number of VDW clashes that the loop has with the framework
    cdef size_t clashes = 0
    # Go through the binding loop residues
    for I1 in range(J1):
        # Get the PyResidue
        residue = residues[I1]
        # Determine if it should be treated as a rotamer or not
        if I1 == 0 or I1 == J1-1:
            rotamer = True
        else:
            rotamer = False
        # Go through the framework pieces
        for I2 in range(J2):
            framework = frameworks[I2]
            # Get the number of VDW clashes
            clashes += framework.calculate_residue_overlap(residue, rotamer)
            # 90% of experimentally determined CDRs treated the same way these
            # loops are have 16 or fewer VDW clashes
            if clashes > 16:
                return result
    # Since the VDW clashes are sufficiently few, next calculate the energies
    # between the binding loop and the framework pieces
    cdef float energy = 0.0
    for I1 in range(J1):
        residue = residues[I1]
        if I1 == 0 or I1 == J1-1:
            rotamer = True
        else:
            rotamer = False
        for I2 in range(J2):
            framework = frameworks[I2]
            energy += framework.calculate_residue_energy(residue, rotamer)
    # 90% of experimentally determined CDRs have energies <= 237.50 using
    # these energy functions
    if energy <= 237.50:
        result = True
    return result

cdef void loop_determination (Loop loop, PyProtein protein, size_t I1, 
                              size_t I2, vector[vector[string]]& atomNames,
                              parameters, str rotlib, size_t index, 
                              float score, list frameworks, str PDBName):
    """Determine whether or not a loop structure is actually acceptable"""
    # Prepare the residues for comparison with the frameworks
    cdef list residues = loop.prepare_residues(protein, I1, I2, atomNames,
                                               parameters, rotlib)
    # Determine if the residues are compatible with the framework or not
    cdef bool accept = assess_framework_compatibility(residues, frameworks)
    # If the result is not acceptable, be done
    if not accept:
        return
    # Since it is acceptable, the binding loop should be output and the
    # information should be stored in the binding loop
    # Increment the count in the loop
    loop.count += 1
    # Create the name for this binding loop file
    label = PDBName.split(".")[0]
    fileName = label + "_" + protein.name() + "_" + str(loop.count) + ".pdb"
    # Store the information in the binding loop
    loop.store_result(fileName, score, index)
    # Create the file
    f = open(loop.folder + fileName, "w")
    # Write a Remark at the top of the file
    remark = "REMARK    Residues " + residues[0].get_label() + " to " \
           + residues[-1].get_label() + " in Protein " + protein.name() \
           + " of PDB file " + label + "\n"
    f.write(remark)
    # Include each residue
    cdef size_t i, j
    j = len(residues)
    atomNum = 1
    for i in range(j):
        atomNum = residues[i].renumber_atoms(atomNum)
        f.write(str(residues[i]))
    # Close the file
    f.close()

cdef void search_PDB_file (instructions, str fileName, str fileFolder, 
        list loops, vector[vector[string]]& atomNames, parameters, 
        list frameworks):
    """Search a PDB file's Proteins for possible binding loop structures"""
    # Load the PDB file
    try:
        pdb = PyPDB (fileName, fileFolder)
    # Skip files that don't load correctly
    except RuntimeError:
        return
    # Extract necessary information from the instructions
    cdef size_t minL = instructions["lengths"][0]
    cdef size_t maxL = instructions["lengths"][1]
    cdef str rotlib = instructions["rotamers"]
    # Looping variables
    cdef size_t I1, I2, I3, I4, J1, J2, J3, J4
    # The number of Proteins in the pdb file
    J1 = pdb.proteins()
    # The number of binding loops
    J4 = len(loops)
    # Vectors of needed atoms
    cdef vector[Atom *] atoms1
    cdef vector[Atom *] atoms2
    # A PyProtein
    cdef PyProtein protein
    # C++ pointer objects
    cdef Protein * protptr
    cdef Residue * resptr1
    cdef Residue * resptr2
    # Necessary variables for accessing residue pointers
    cdef char L = 32
    cdef bool how = True
    # A Loop object
    cdef Loop loop
    # Whether or not an activity was successful
    cdef bool accept
    # A compatibility score for a loop
    cdef float score
    # And an index
    cdef int index
    # Go through the proteins in the PDB file
    for I1 in range(J1):
        # Get the protein
        protein = pdb.protein(I1)
        # The protein's pointer
        protptr = protein._thisptr
        # The number of Residues in the protein. The search will actually only
        # consider the second residue to the second to last (don't want to
        # include terminal amino acids), so 1 is subtracted from that.
        J2 = protptr.size()-1
        # If there are too few residues, skip this protein
        if J2 < minL+1:
            continue
        # Reset the counts for the binding loops
        for I4 in range(J4):
            loop = loops[I4]
            loop.count = 0
        # Go through the positions that can be the start of a binding loop
        for I2 in range(1, J2 - (minL-1)):
            # Get the Residue's pointer
            resptr1 = protptr[0](I2, L, how)
            # If the residue is not present, continue the search
            if not resptr1.is_present():
                continue
            # If the previous residue is not present, continue the search
            if not protptr[0](I2-1, L, how).is_present():
                continue
            # Get the atoms from the residue
            accept = get_atom_pointers(resptr1, atomNames[0], atoms1)
            # If they weren't all found, skip this residue
            if not accept:
                continue
            # Go through the residues that are binding loop compatible with
            # this residue
            J3 = I2 + maxL
            if J3 > J2:
                J3 = J2
            for I3 in range(I2+1, J3):
                # Get this residue's pointer
                resptr2 = protptr[0](I3, L, how)
                # If it isn't present, stop the search
                if not resptr2.is_present():
                    break
                # If the next residue is not present, stop the search
                if not protptr[0](I3+1, L, how).is_present():
                    break
                # If it is not yet enough residues, continue the search
                if (I3-I2) + 1 < minL:
                    continue
                # Get the atoms from the residue
                accept = get_atom_pointers(resptr2, atomNames[1], atoms2)
                # If they weren't all found, continue the search
                if not accept:
                    continue
                # Go through the binding loops
                for I4 in range(J4):
                    loop = loops[I4]
                    # Calculate the score for these atoms with the binding
                    # loop
                    score = loop.compatible(atoms1, atoms2)
                    # If the score is negative, the atoms aren't compatible
                    if score < 0:
                        continue
                    # Get an index for the structure in the binding loop
                    index = loop.find_index(score)
                    # If the index is negative, the structure shouldn't be
                    # saved
                    if index < 0:
                        continue
                    # If the score and index are acceptable, do the final
                    # determination of whether or not the loop is OK and save
                    # it's structure
                    loop_determination(loop, protein, I2, I3, atomNames,
                            parameters, rotlib, index, score, frameworks,
                            fileName)

def PDB_Search (instructions):
    """Search PDB files for binding loop structures"""
    # Do the preliminary preparation
    summary, parameters, model, frameworks, loops = prepare_for_search(instructions)
    # Load the PDB file information
    SourceFile = instructions["PDB"]
    # Try to open that file
    try:
        f = open(SourceFile, "r")
    except IOError:
        error = "Failure to open the PDB File List: " + SourceFile + "\n"
        raise BindingProteinError(error)
    # Store the information here
    fileData = []
    for line in f:
        items = line.split()
        if len(items) == 2:
            fileData.append(items)
    f.close()
    # If there are no files, raise an error
    if len(fileData) == 0:
        error = "The PDB File List did not specify any PDB files\n"
        raise BindingProteinError (error)
    # Get the names of the atom pairs
    cdef vector[vector[string]] atomNames = get_atom_names(instructions)
    # Go through every file
    for pair in fileData:
        # Search the PDB file
        search_PDB_file(instructions, pair[0], pair[1], loops, atomNames,
                parameters, frameworks)
    # Make sure that every loop has at least 1 compatible structure
    error = ""
    for loop in loops:
        if len(loop) == 0:
            error += "No " + loop.name + " binding loop structures were "
            error += "identified.\n"
    if len(error) > 0:
        raise BindingProteinError(error)
    # Write to the summary file that this is finished
    message = ""
    for loop in loops:
        message += loop.name + " Binding Loop: " + str(len(loop))
        message += " structures identified\n"
    message += "Searching for binding loop structures ended on "
    message += time_stamp() + "\n"
    summary.write(message)
    summary.close()
