# Written in 2018 by Dr. Robert Pantazes. Copyright owned by Auburn
# University.
# This file contains the classes and functions that are needed to determine
# whether or not structures of different binding loops are compatible with one
# another or not. The results of these functions and files are binding loop
# structures that are guaranteed to be able to be part of a solution.

# Cython compiler directives
# cython: c_string_type=str
# cython: c_string_encoding=ascii
# cython: language_level=3
# distutils: language=c++

# Include pure Python modules
import os
import sys
# Include Binding Protein Design contents
from BPD.BindingProteinError import BindingProteinError
from BPD.General import parameterize_protein
from BPD.Instructions import time_stamp
from BPD.Database.Preliminary import database_folder_name
from BPD.Database.Preliminary import summary_file_name
# Include the function for loading energy parameters
from ForceFields.Parameters import load_parameters
# Include the PyProtein and PyResidue classes
from Proteins.Protein.PyProtein cimport PyProtein
from Proteins.Residue.PyResidue cimport PyResidue
# Include the C++ boolean class
from libcpp cimport bool

# Define a Structure class here, where a Structure is a specific Protein in a
# specific Binding Loop. 
cdef class Structure:

    # The class attributes are:
    # The protein structure
    cdef readonly PyProtein protein
    # Whether or not the protein's structure is currently loaded
    cdef readonly bool loaded
    # The name of this Structure
    cdef readonly str name
    # Lists of the Structures belonging to other Binding Loops that this
    # structure is compatible with 
    cdef public dict compatible
    # Whether or not this sturcture is still known to be part of possible
    # solutions
    cdef public bool useable

    def __init__ (self, label, otherLoops):
        """Initialize the Structure class object"""
        # Set the name of the structure
        self.name = label
        # Set up the compatible dictionary
        self.compatible = {}
        for L in otherLoops:
            self.compatible[L] = []
        # The protein is not initially loaded
        self.loaded = False
        # Set the protein as a blank PyProtein. This shouldn't be necessary,
        # but do it to be safe
        self.protein = PyProtein()
        # The Structure starts out able to be used in solutions
        self.useable = True

    def load (self, folder, parameters):
        """Load the PyProtein structure"""
        # If the protien is already loaded, raise an error
        if self.loaded:
            error = "PyProtein " + self.name + " is already loaded - it "
            error += "can't be loaded again\n"
            raise BindingProteinError (error)
        # Load the PyProtein
        self.protein = PyProtein (self.name + ".pdb", folder)
        # Parameterize the protein
        parameterize_protein(self.protein, parameters)
        # Indicate that the protein is loaded
        self.loaded = True

    def clean_up (self):
        """Delete the existing Protein's pointer"""
        # Only do this if the Protein is loaded
        if self.loaded:
            # Delete the Protein pointer
            del self.protein._thisptr
            # Set the pointer to null
            self.protein._thisptr = NULL
            # Indicate that the protein is no longer loaded
            self.loaded = False

    cpdef bool compare (self, Structure other):
        """Compare two Structures to see if they are compatible"""
        # The answer is this
        cdef bool answer = False
        # Get the number of VDW clashes between the two structures
        cdef size_t clashes = self.protein.calculate_protein_overlap(other.protein)
        # 90% of naturally occurring CDR pairs have 3 or fewer VDW clashes
        if clashes > 3:
            return answer
        # Get the energy between the two structures
        cdef float energy = self.protein.calculate_protein_energy(other.protein)
        # 90% of naturally occurring CDR pairs have interaction energies <=
        # 55.6
        if energy <= 55.6:
            answer = True
        return answer

    def update_label (self, int I, str L, Structure other):
        """Update the label information in the other Structure"""
        if L in other.compatible:
            try:
                index = other.compatible[L].index(self.name)
                other.compatible[L][index] = I
            except ValueError:
                pass

    def remove_self (self, str L, Structure other):
        """Remove this Structure from being compatible with another"""
        index = other.compatible[L].index(self.name)
        waste = other.compatible[L].pop(index)

    def renumber (self, int n1, int n2, str L):
        """Renumber the protein's residues"""
        # Letters for intermediate positions
        cdef str Letters = ' ABCDEFGHIJKLMNOPQRSTUVWXYZ'
        # Set the Protein's name
        self.protein.set_name(L)
        # Renumber the atoms in the binding loop structure
        atomNum = self.protein.renumber_atoms(1)
        # Information for numbering from the start of the protein
        cdef int FRONT = -1
        cdef int NUM1 = n1
        cdef size_t L1 = 0
        # Information for numbering from the end of the protein
        cdef int END = len(self.protein)
        cdef int NUM2 = n2
        cdef size_t L2 = 0
        # A PyResidue object
        cdef PyResidue res
        # Go through all of the residues, alternating between those at the
        # front and those at the end of the protein
        while FRONT < END - 1:
            # Increment the residue counter from the start of the protein
            FRONT += 1
            # Get the residue
            res = self.protein[FRONT]
            # Determine how to renumber the residue
            if NUM1 < NUM2 - 1:
                res.set_number(NUM1, ' ')
                NUM1 += 1
            else:
                res.set_number(NUM1, Letters[L1])
                L1 += 1
            # If there are more residues to number from the end of the loop
            if END > FRONT + 1:
                END -= 1
                res = self.protein[END]
                if NUM2 > NUM1 + 1:
                    res.set_number(NUM2, ' ')
                    NUM2 -= 1
                else:
                    res.set_number(NUM2, Letters[L2])
                    L2 += 1

    def output (self, int NUM, str LoopFolder, str nameStart, 
                int n1, int n2, str L, list loopOrder, 
                str sourceFolder):
        """Output the Structure's information"""
        # If the protein is not currently loaded
        if not self.loaded:
            self.protein = PyProtein(self.name + ".pdb", sourceFolder)
            self.loaded = True
        # Get the REMARK from the start of that file listing the file, chain
        # and residues of the protein piece
        f = open(sourceFolder + self.name + ".pdb", "r")
        REMARK = f.readline().strip()
        f.close()
        # Renumber the protein's residues
        self.renumber(n1, n2, L)
        # The structure should be written to this file
        fileName = LoopFolder + "Structures/" + nameStart + str(NUM) + ".pdb"
        # Output the Protein's structure
        f = open(fileName, "w")
        f.write(REMARK + "\n")
        f.write(str(self.protein))
        f.close()
        # Remove the protein from memory
        self.clean_up()
        # The compatibility information should be written to this file
        fileName = LoopFolder + "Compatible/" + nameStart + str(NUM) + ".txt"
        f = open(fileName, "w")
        # Go through the binding loop order
        for loop in loopOrder:
            # If loop isn't in the compatible dictionary, skip it
            if loop not in self.compatible:
                continue
            for N in self.compatible[loop]:
                f.write(loop + " " + str(N) + "\n")
        # close the file
        f.close()

cdef class Loop:
    """All the structures of a binding loop"""

    # The name of the loop
    cdef readonly str name
    # The various folders of the binding loop
    cdef readonly str pdbFolder
    cdef readonly str clusterFolder
    cdef readonly str mainFolder
    # The name of the protein the loop goes in 
    cdef readonly str protein
    # The numbers of the first and last residues in the binding loop
    cdef readonly int firstRes
    cdef readonly int lastRes
    # The list of structures that the Binding Loop contains
    cdef public list structures

    def __init__ (self, label, instructions):
        """Initialize the Binding Loop object"""
        # Set the name of the loop
        self.name = label
        # Get the folder for the database
        DBFolder = database_folder_name(instructions)
        # Set the names of the folders
        self.pdbFolder = DBFolder + self.name + "/PDB/"
        self.clusterFolder = DBFolder + self.name + "/Clustering/"
        self.mainFolder = DBFolder + self.name + "/"
        # Try to make the structure folder and raise an error if that fails
        try:
            os.mkdir(self.mainFolder + "Structures/")
        except OSError:
            error = "Failure to make the " + self.name + " Binding Loop "
            error += "Structures folder. It likely already exists.\n"
            raise BindingProteinError (error)
        # Try to make the compatibility folder and raise an error if that
        # fails
        try:
            os.mkdir(self.mainFolder + "Compatible/")
        except OSError:
            error = "Failure to make the " + self.name + " Binding Loop "
            error += "Compatible folder. It likely already exists.\n"
            raise BindingProteinError (error)
        # Find the information about this binding loop in the instructions
        found = False
        otherLoops = []
        for data in instructions["loops"]:
            if data[0] == self.name:
                found = True
                self.protein = data[1]
                try:
                    self.firstRes = int(data[2])
                except ValueError:
                    self.firstRes = int(data[2][:-1]) + 1
                try:
                    self.lastRes = int(data[3])
                except ValueError:
                    self.lastRes = int(data[3][:-1]) - 1
            else:
                otherLoops.append(data[0])
        # If the information wasn't found, raise an error
        if not found:
            error = "Failure to find the " + self.name + " instructions.\n"
            raise BindingProteinError (error)
        # Use the clustering information to get the structures
        self.structures = []
        # Go through the possible loop lengths
        for N in range(instructions["lengths"][0], instructions["lengths"][1]+1):
            # The file is:
            fileName = self.clusterFolder + self.name + "_" + str(N)
            fileName += "_Clusters.txt"
            # Try to open the file
            try:
                f = open(fileName, "r")
            # It is OK if that fails - it just means there were no loops of
            # that length identified
            except IOError:
                continue
            # This flag indicates whether or not a line includes a structure
            # listing
            flag = False
            # This list includes sequences in a cluster that are already
            # included
            sequences = []
            # And this flag indicates whether or not those sequences are used
            # to eliminate members of a cluster
            RemoveDuplicates = False
            # Go through the contents of the file
            for line in f:
                # If the line indicates that the flag should be turned off
                if line.startswith("Cluster "):
                    flag = False
                # Or if it should be turned on
                elif line.startswith("Structures"):
                    flag = True
                    # Empty the sequences list
                    sequences = []
                # If this line might contain a sequence
                elif flag:
                    # Split the line into pieces
                    items = line.split()
                    # If there are at least three, store the information
                    if len(items) > 3:
                        name = items[0]
                        sequence = items[2:]
                        # If duplicated sequences should be excluded
                        if RemoveDuplicates and sequence in sequences:
                            continue
                        # Store the sequence
                        sequences.append(sequence)
                        # Create and store a structure
                        self.structures.append(Structure(name, otherLoops))
            # Close the file
            f.close()
        # If there are no structures, raise an error
        if len(self.structures) == 0:
            error = "No " + self.name + " Binding Loop Structures were "
            error += "created.\n"
            raise BindingProteinError (error)

    cdef void load_structures (self, parameters):
        """Load the Loop's structures"""
        # Looping variables
        cdef size_t I, J
        J = len(self.structures)
        # Go through the structures
        for I in range(J):
            self.structures[I].load(self.pdbFolder, parameters)

    cdef void clean_up (self):
        """Delete the PyProtein portions of the Loop's structures"""
        # Looping variables
        cdef size_t I, J
        J = len(self.structures)
        for I in range(J):
            self.structures[I].clean_up()

    cdef void compare (self, Loop other):
        """Compare two Loops' structures to find compatibilities"""
        # Looping variables
        cdef size_t I1, I2, J1, J2
        J1 = len(self.structures)
        J2 = len(other.structures)
        # Structure objects
        cdef Structure protein1
        cdef Structure protein2
        # A boolean value for when proteins are compatible
        cdef bool accept
        # Go through the structure pairs
        for I1 in range(J1):
            protein1 = self.structures[I1]
            for I2 in range(J2):
                protein2 = other.structures[I2]
                # Determine if they are compatible or not
                accept = protein1.compare(protein2)
                # If they are mutually acceptable
                if accept:
                    protein1.compatible[other.name].append(protein2.name)
                    protein2.compatible[self.name].append(protein1.name)

    def __len__ (self):
        """The number of structures in the Loop"""
        return len(self.structures)

def prepare (instructions):
    """Do the initial preparations / calculations"""
    # Open the summary file
    summary = open(summary_file_name(instructions), "a")
    # Write a message saying what is happening
    message = "Finding compatible binding loop structures started on "
    message += time_stamp()
    summary.write(message)
    summary.flush()
    # Create the Loops
    loops = []
    for data in instructions["loops"]:
        loops.append(Loop(data[0], instructions))
    # Load the energy calculation parameters
    parameters = load_parameters(instructions["topology"],
                                 instructions["parameters"], 
                                 instructions["solvation"])
    # Return the generated variables
    return summary, loops, parameters

def find_compatible_loops (list loops, parameters):
    """Find the compatible binding loop structures"""
    # Looping variables
    cdef size_t I1, I2, J1, J2
    J2 = len(loops)
    J1 = J2 - 1
    # Loop variables
    cdef Loop loop1, loop2
    # Go through the loop pairs
    for I1 in range(J1):
        loop1 = loops[I1]
        # Load the loop's structures
        loop1.load_structures(parameters)
        # Go through the subsequent structures
        for I2 in range(I1+1, J2):
            loop2 = loops[I2]
            # Load the structures
            loop2.load_structures(parameters)
            # Compare the two loops
            loop1.compare(loop2)
            # Delete the PyProteins from the second loop
            loop2.clean_up()
        # Delete the PyProteins from the first loop
        loop1.clean_up()

def purge_unusable_structures (list loops):
    """Remove structures that can't be part of solutions"""
    # Do this iteratively until changes stop occurring
    cdef bool changes = True
    # And indicate whether or not ANY changes occurred
    cdef bool anyChanges = False
    # A flag to indicate if problems are identified for a structure
    cdef bool problems = False
    # A list of protein names
    cdef list names
    # Necessary variables
    cdef Loop loop1, loop2
    cdef Structure protein1, protein2
    cdef I1, I2, I3, I4, J1, J2, J3
    J1 = len(loops)
    # Use a while loop
    while changes:
        changes = False
        # Go through the binding loops
        for I1 in range(J1):
            loop1 = loops[I1]
            # Go through the loop's structures BACKWARDS
            J2 = len(loop1.structures)
            for I2 in range(J2, 0, -1):
                protein1 = loop1.structures[I2-1]
                # If the protein has been deactivated, it should be removed
                # from consideration
                problems = (not protein1.useable)
                # If it doesn't have known problems already, check that it is
                # compatible with at least one structure in every other
                # binding loop
                if not problems:
                    # Go through the protein's compatibility dictionary
                    for L in protein1.compatible:
                        # If there are no compatible structures in a binding loop
                        if len(protein1.compatible[L]) == 0:
                            problems = True
                            break
                # If there are no problems, go to the next protein
                if not problems:
                    continue
                # If there are problems, the overall while loop should be
                # repeated
                changes = True
                anyChanges = True
                # Go through all of the binding loops
                for I3 in range(J1):
                    # Skip the loop if it is the same as the one the structure
                    # is in
                    if I3 == I1:
                        continue
                    # Get the binding loop
                    loop2 = loops[I3]
                    # Get the names of the structures in loop2 that protein1
                    # is compatible with
                    names = protein1.compatible[loop2.name]
                    # If there aren't any, skip this loop
                    if len(names) == 0:
                        continue
                    # Get the number of proteins in the binding loop
                    J3 = len(loop2)
                    for I4 in range(J3):
                        # Get the protein
                        protein2 = loop2.structures[I4]
                        # If this protein is compatible with protein1,
                        # eliminate that information
                        if protein2.name in names:
                            protein1.remove_self(loop1.name, protein2)
                # Pop this protein out of the structures list
                protein1 = loop1.structures.pop(I2-1)
            # If the loop's structures list is now empty, raise an error
            if len(loop1.structures) == 0:
                error = "No " + loop1.name + " structures were identified "
                error += "that were compatible with other loop structures\n"
                raise BindingProteinError (error)
    return anyChanges

cdef bool check_for_solution (list structures, list loopOrder, 
                              size_t loopIndex, list loops):
    # The answer to the function
    cdef bool answer = False
    # The number of loops in the order
    cdef size_t LON = len(loopOrder)
    # If the loop index has gone too far, the function was successful
    if loopIndex >= LON:
        answer = True
        return answer
    # The loop that needs to be considered this recursion
    cdef Loop loop
    cdef size_t I1, I2, J1, J2
    J1 = len(loops)
    for I1 in range(J1):
        loop = loops[I1]
        if loop.name == loopOrder[loopIndex]:
            break
    # Structure objects
    cdef Structure protein1, protein2
    J2 = len(structures)
    # A flag to indicate whether or not a structure is acceptable
    cdef bool allowed = False
    # A list of the structures that are compatible with all the current
    # structures
    cdef list identified = []
    # Go through the structures of the loop
    J1 = len(loop.structures)
    for I1 in range(J1):
        protein1 = loop.structures[I1]
        # If the protein is not useable, skip it
        if not protein1.useable:
            continue
        # Indicate that this protein is preliminarily allowed
        allowed = True
        # Go through the current structures
        for I2 in range(J2):
            protein2 = structures[I2]
            # If the protein is not compatible with this structure, it is not
            # allowed
            if protein1.name not in protein2.compatible[loop.name]:
                allowed = False
                break
        # If the protein is allowed, store it in the identified structures
        if allowed:
            identified.append(protein1)
    # If there are no identified structures, no solution is possible
    J1 = len(identified)
    if J1 == 0:
        return answer
    # Continue the recursion
    cdef list newStructures
    for I1 in range(J1):
        newStructures = []
        newStructures.extend(structures)
        newStructures.append(identified[I1])
        answer = check_for_solution(newStructures, loopOrder, loopIndex+1,
                                    loops)
        # If a successful answer was found, return that
        if answer:
            return answer
    # Return that no successful answer was identified
    return answer

def compatibility_checking (list loops):
    """Check that every structure can be part of a solution"""
    # Looping variables
    cdef size_t I1, I2, J1, J2
    J1 = len(loops)
    # Various variables needed for the solution
    cdef Loop loop1, loop2
    cdef Structure protein
    # A list of loop names
    cdef list loopOrder
    # Whether or not a particular structure is acceptable
    cdef bool acceptable
    # Do a preliminary purging of the binding loop structures
    cdef bool changes = purge_unusable_structures(loops)
    # Do the following steps until changes stop happening
    changes = True
    while changes:
        changes = False
        # Go through each binding loop
        for I1 in range(J1):
            loop1 = loops[I1]
            # Make a list of the names of the other binding loops
            loopOrder = []
            for I2 in range(J1):
                if I2 != I1:
                    loopOrder.append(loops[I2].name)
            # Go through the proteins of the binding loop
            J2 = len(loop1.structures)
            for I2 in range(J2):
                protein = loop1.structures[I2]
                # Determine whether or not the protein is acceptable
                acceptable = check_for_solution ([protein], loopOrder, 0,
                                                 loops)
                # If it is not acceptable, mark it as unusable
                if not acceptable:
                    changes = True
                    protein.useable = False
        # If there were changes identified, purge unuseable structures
        if changes:
            acceptable = purge_unusable_structures(loops)

def update_labels (list loops):
    """Update the naming information in the compatible dictionaries"""
    # Looping variables
    cdef int I1, I2, I3, I4, J1, J2, J3, J4
    # Data types
    cdef Loop loop1, loop2
    cdef Structure protein1, protein2
    # Go through the binding loops
    J1 = len(loops)
    for I1 in range(J1):
        loop1 = loops[I1]
        # Go through the structures of that loop
        J2 = len(loop1.structures)
        for I2 in range(J2):
            protein1 = loop1.structures[I2]
            # Go through all other proteins in all other loops
            for I3 in range(J1):
                if I3 == I1:
                    continue
                loop2 = loops[I3]
                J4 = len(loop2.structures)
                for I4 in range(J4):
                    protein2 = loop2.structures[I4]
                    # Update the labelling information in protein 2
                    protein1.update_label(I2+1, loop1.name, protein2)

def Compatibility (instructions):
    """Identify which binding loop structures can be part of solutions"""
    # Do the preliminary calculations
    summary, loops, parameters = prepare(instructions)
    # Find out which binding loop structures are compatible with one another
    find_compatible_loops (loops, parameters)
    # Iteratively find structures that are able to be part of solutions in
    # designed binding proteins
    compatibility_checking(loops)
    # Update the labelling information for the remaining structures
    update_labels (loops)
    # Output the binding loop structures
    for loop in loops:
        # Make a name for the start of every file
        nameStart = instructions["name"] + "_" + loop.name + "_"
        # Make list of the other loops
        otherLoops = []
        for other in loops:
            if other.name != loop.name:
                otherLoops.append(other.name)
        # Go through the structures of the binding loop
        for i, protein in enumerate(loop.structures):
            protein.output(i+1, loop.mainFolder, nameStart, loop.firstRes,
                           loop.lastRes, loop.protein, otherLoops, 
                           loop.pdbFolder)
    # Update the summary file with what was found
    message = ""
    for loop in loops:
        message += loop.name + " Binding Loop: " + str(len(loop))
        message += " usable structures identified\n"
    message += "Finding compatible binding loop structures ended on "
    message += time_stamp() + "\n"
    summary.write(message)
    summary.flush()
