# Written in 2018 by Dr. Robert Pantazes. Copyright owned by Auburn
# University.
# This file contains General Functions used in creating databases of binding
# protein pieces and designing binding proteins with those pieces. 

# Cython compiler directives
# cython: c_string_type=str
# cython: c_string_encoding=ascii
# cython: language_level=3
# distutils: language=c++

# Include Python modules
import os
import sys
# Include the error class
from BPD.BindingProteinError import BindingProteinError
# C-import classes, too
from Proteins.Protein.PyProtein cimport PyProtein
from Proteins.Residue.PyResidue cimport PyResidue
from Proteins.Atom.PyAtom cimport PyAtom
from Proteins.Proteins cimport Protein
from Proteins.Proteins cimport Residue
from Proteins.Proteins cimport Atom
from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector

def get_residue (PyProtein prot, str label):
    """Get the specified residue from the protein"""
    # Looping variables
    cdef size_t i, j
    # A PyResidue object
    cdef PyResidue res
    # The number of residues in the protein
    j = prot._thisptr.size()
    # Go through the residues in the protein
    for i in range(j):
        # Get the residue
        res = prot[i]
        # If it is the right residue, return it
        if res.get_label() == label:
            return res
    # If the residue was not found, raise an error
    error = "Residue " + label + " was not found in protein "
    error += prot.name() + "\n"
    raise BindingProteinError (error)


def get_atoms (PyResidue res, list atomNames, errorFlag = True):
    """Get the specified atoms from the residue"""
    # Store the atoms here
    cdef list atoms = []
    # Looping variables
    cdef size_t i1, i2, j1, j2
    j2 = res._thisptr.size()
    j1 = len(atomNames)
    # the name of each atom
    cdef string name
    cdef str atomName
    # A flag to indicate whether or not the atom was found
    cdef bool flag
    # Go through the atom names
    for i1 in range(j1):
        # Get the name of the atom
        atomName = atomNames[i1]
        # If this is a GLY and the atom is CB, use HA1
        if res.name() == "GLY" and atomName == "CB":
            atomName = "HA1"
        # Store the name of the atom
        name = atomName
        # Indicate whether or not the atom is found
        flag = False
        # Go through the atoms in the residue
        for i2 in range(j2):
            # If this is the right atom
            if res._thisptr[0][i2].name() == name:
                # Store the PyAtom
                atoms.append(res[i2])
                flag = True
                break
        # If the atom was not found, raise an error
        if not flag:
            # If an error should not be created, just return that the function
            # failed
            if not errorFlag:
                return False
            # Otherwise, make an error
            error = "This Residue does not contain each of these Atoms:"
            for atomName in atomNames:
                error += " " + atomName
            error += "\n"
            error += str(res)
            raise BindingProteinError (error)
    # Return the identified atoms
    return atoms

def parameterize_residue (PyResidue residue, parameters = None):
    """Assign VDW parameters to the Atoms in a Residue"""
    # A PyAtom object
    cdef PyAtom atom
    # Looping variables
    cdef size_t I, J
    J = len(residue)
    # The name of an atom
    cdef str name
    # Go through the Residue's atoms
    for I in range(J):
        # Get the Atom
        atom = residue[I]
        # Get the first letter of the atom name
        name = atom.name()[0]
        # If it is a hydrogen
        if name == "H":
            atom._thisptr.set_vdw_radius(1.20)
        # If it is a Carbon
        elif name == "C":
            atom._thisptr.set_vdw_radius(1.70)
        # Or if it is a nitrogen
        elif name == "N":
            atom._thisptr.set_vdw_radius(1.55)
        # Or if it is an oxygen
        elif name == "O":
            atom._thisptr.set_vdw_radius(1.52)
        # If it is a sulfur
        elif name == "S":
            atom._thisptr.set_vdw_radius(1.80)
        # For anything else, raise an error
        else:
            error = "Failure to parameterize this Atom:\n" + str(atom)
            raise BindingProteinError(error)

def parameterize_protein (PyProtein protein, parameters = None):
    """Assign VDW parameters to a Protein's Residues"""
    # A residue in the protein
    cdef PyResidue residue
    # Looping variables
    cdef size_t I, J
    J = protein._thisptr.size()
    # Go through the residues
    for I in range(J):
        # Get the residue
        residue = protein[I]
        # Parameterize the Residue
        parameterize_residue(residue, parameters)

cdef void allocate_atoms (Residue * residue, vector[Atom *]& backbone,
                          vector[Atom *]& sidechain):
    """Group a Residue's non-hydrogen atoms into backbone or sidechain"""
    # Looping variables
    cdef size_t I, J
    J = residue.size()
    # Clear the two vectors
    backbone.clear()
    sidechain.clear()
    # An atom pointer
    cdef Atom * atom
    # Go through the atoms
    for I in range(J):
        # Get the atom pointer
        atom = residue[0][I]
        # If the atom is a hydrogen, skip it
        if atom.is_hydrogen():
            continue
        # Store the atom in the proper vector
        if atom.is_backbone_atom():
            backbone.push_back(atom)
        else:
            sidechain.push_back(atom)

cdef bool vdw_clash (Atom * atom1, Atom * atom2, bool SC = False):
    """Determine whether or not two Atoms have a VDW clash"""
    # Set a multiplier for the clashes depending on the SC variable
    cdef float multiplier = 0.8
    if SC:
        multiplier = 0.5
    # Calculate the distance between the two atoms
    cdef bool squared = False
    cdef float dis = atom1.calculate_distance(atom2, squared)
    # The answer to this function
    cdef bool answer = False
    # If the atoms are closer than 80% of their VDW radii
    if dis < multiplier * (atom1.get_vdw_radius() + atom2.get_vdw_radius()):
        answer = True
    return answer

cdef void residue_vdw_clashes (list proteins, list antigens, 
                               vector[size_t]& clashes, size_t method):
    """Find how many residues in the proteins have VDW clashes"""
    # Declare the necessary objects
    cdef PyProtein protein, antigen
    cdef PyResidue protRes, antRes
    cdef vector[Atom *] protBB, protSC, antBB, antSC
    cdef vector[size_t] counts
    cdef size_t I1, I2, I3, I4, I5, I6, J1, J2, J3, J4, J5, J6
    # The number of different clash types there can be
    cdef size_t N = 3
    # Make sure the clashes are set to 0
    clashes.clear()
    clashes.resize(N)
    for I1 in range(N):
        clashes[I1] = 0
    # Make sure the counts are set to a size of 3, too
    counts.resize(N)
    # Set the numbers of proteins and antigens
    J1 = len(proteins)
    J3 = len(antigens)
    # Go through the proteins
    for I1 in range(J1):
        # Get the protein
        protein = proteins[I1]
        # The number of residues in the protein
        J2 = protein._thisptr.size()
        # Go through the residues
        for I2 in range(J2):
            # Get the Residue
            protRes = protein[I2]
            # Allocate its atoms into backbone and sidechain groups
            allocate_atoms (protRes._thisptr, protBB, protSC)
            # Set the counts for this residue to 0
            for I3 in range(N):
                counts[I3] = 0
            # Go through the antigen residues
            for I3 in range(J3):
                antigen = antigens[I3]
                J4 = antigen._thisptr.size()
                for I4 in range(J4):
                    antRes = antigen[I4]
                    # Allocate it's atoms
                    allocate_atoms (antRes._thisptr, antBB, antSC)
                    # Always count the number of clashing backbone atoms
                    J5 = protBB.size()
                    J6 = antBB.size()
                    for I5 in range(J5):
                        for I6 in range(J6):
                            # If the atoms have a VDW clash, increment the
                            # first count
                            if vdw_clash (protBB[I5], antBB[I6]):
                                counts[0] += 1
                    # If the method is at least 2, count the number of side
                    # chain to backbone VDW clashes
                    if method > 1:
                        J6 = antSC.size()
                        for I5 in range(J5):
                            for I6 in range(J6):
                                if vdw_clash(protBB[I5], antSC[I6], True):
                                    counts[1] += 1
                        J5 = protSC.size()
                        J6 = antBB.size()
                        for I5 in range(J5):
                            for I6 in range(J6):
                                if vdw_clash(protSC[I5], antBB[I6], True):
                                    counts[1] += 1
                    # If the method is at least 3, count the number of
                    # sidechain - sidechain VDW clashes
                    if method > 2:
                        J6 = antSC.size()
                        for I5 in range(J5):
                            for I6 in range(J6):
                                if vdw_clash(protSC[I5], antSC[I6], True):
                                    counts[2] += 1
            # If this protein residue has clashes, store that information
            for I3 in range(N):
                if counts[I3] > 0:
                    clashes[I3] += 1

cdef void get_positive_atoms (list proteins, vector[Atom *]& atoms):
    """Get the positively charged atoms from a set of proteins"""
    # Looping variables
    cdef size_t I1, I2, I3, J1, J2, J3
    J1 = len(proteins)
    # Py objects
    cdef PyProtein protein
    cdef PyResidue residue
    cdef PyAtom atom
    # Go through the proteins
    for I1 in range(J1):
        protein = proteins[I1]
        # Go through the residues
        J2 = protein._thisptr.size()
        for I2 in range(J2):
            residue = protein[I2]
            # If the residue is ARG, get the NH1 and NH2 atoms
            if residue.name() == "ARG":
                J3 = residue._thisptr.size()
                for I3 in range(J3):
                    atom = residue[I3]
                    if atom.name() in ['NH1', 'NH2']:
                        atoms.push_back(atom._thisptr)
            # If the residue is LYS, get the NZ atom
            elif residue.name() == "LYS":
                J3 = residue._thisptr.size()
                for I3 in range(J3):
                    atom = residue[I3]
                    if atom.name() == "NZ":
                        atoms.push_back(atom._thisptr)
                        break

cdef void get_negative_atoms (list proteins, vector[Atom *]& atoms):
    """Get the negatively charged atoms from a set of proteins"""
    # This is exactly the same as the previous function, except it collects
    # OD1/2 atoms from ASP and OE1/2 atoms from GLU
    cdef size_t I1, I2, I3, J1, J2, J3
    J1 = len(proteins)
    cdef PyProtein protein
    cdef PyResidue residue
    cdef PyAtom atom
    for I1 in range(J1):
        protein = proteins[I1]
        J2 = protein._thisptr.size()
        for I2 in range(J2):
            residue = protein[I2]
            if residue.name() == "ASP":
                J3 = residue._thisptr.size()
                for I3 in range(J3):
                    atom = residue[I3]
                    if atom.name() in ['OD1', 'OD2']:
                        atoms.push_back(atom._thisptr)
            elif residue.name() == "GLU":
                J3 = residue._thisptr.size()
                for I3 in range(J3):
                    atom = residue[I3]
                    if atom.name() in ['OE1', 'OE2']:
                        atoms.push_back(atom._thisptr)

cdef size_t charge_clashes (list proteins, list antigens):
    """Find the number of charge / charge clashes between the proteins"""
    # Get the positively and negatively charged atoms from the proteins
    cdef vector[Atom *] protPos
    get_positive_atoms (proteins, protPos)
    cdef vector[Atom *] protNeg
    get_negative_atoms (proteins, protNeg)
    # Do the same with the antigens
    cdef vector[Atom *] antPos
    get_positive_atoms (antigens, antPos)
    cdef vector[Atom *] antNeg
    get_negative_atoms (antigens, antNeg)
    # The answer from this function
    cdef size_t answer = 0
    # Looping variables
    cdef size_t I1, I2, J1, J2
    # A distance value
    cdef float dis
    # A boolean value for distance calculations
    cdef bool squared = False
    # Get the numbers of positively charged atoms
    J1 = protPos.size()
    J2 = antPos.size()
    # If there are positively charged atoms in both lists, loop through them
    if J1 > 0 and J2 > 0:
        for I1 in range(J1):
            for I2 in range(J2):
                dis = protPos[I1].calculate_distance(antPos[I2], squared)
                if dis < 4.2:
                    answer += 1
    # Do the same for the negatively charged atoms
    J1 = protNeg.size()
    J2 = antNeg.size()
    if J1 > 0 and J2 > 0:
        for I1 in range(J1):
            for I2 in range(J2):
                dis = protNeg[I1].calculate_distance(antNeg[I2], squared)
                if dis < 4.2:
                    answer += 1
    return answer

def loop_framework_compatibility (PyProtein loop, list frameworks):
    """Determine whether or not a loop is compatible with the frameworks"""
    # Find the number of loop residues with VDW clashes with the frameworks
    cdef vector[size_t] clashes
    residue_vdw_clashes ([loop], frameworks, clashes, 1)
    # If a failure condition is identified, return that
    # There should be 2 loop residues that have backbone clashes with the
    # frameworks (the attachment points)
    if clashes[0] > 2:
        return 1
    # The comparison values should be modified if / when the method for this
    # function is changed from 1 to another value
    elif clashes[1] > 0:
        return 2
    elif clashes[2] > 0:
        return 3
    # Find the number of charge / charge clashes
    cdef size_t charges = charge_clashes ([loop], frameworks)
    # This threshold should be set based on an analysis of data that hasn't
    # been performed yet
    if charges > 0:
        return 4
    # If the function was successful, return 0
    return 0

def loop_loop_compatibility (PyProtein loop1, PyProtein loop2):
    """Determine whether or not 2 binding loops are compatible"""
    # Do the same things as the previous function, just with tweaked methods
    # and cutoffs
    cdef vector[size_t] clashes
    residue_vdw_clashes ([loop1], [loop2], clashes, 1)
    if clashes[0] > 0:
        return 1
    elif clashes[1] > 0:
        return 2
    elif clashes[2] > 0:
        return 3
    cdef size_t charges = charge_clashes ([loop1], [loop2])
    if charges > 0:
        return 4
    return 0

def loop_antigen_compatibility (PyProtein loop, list antigens):
    """Determine whether or not a loop is compatible with the antigens"""
    cdef vector[size_t] clashes
    residue_vdw_clashes ([loop], antigens, clashes, 2)
    if clashes[0] > 0:
        return 1
    elif clashes[1] > 0:
        return 2
    elif clashes[2] > 0:
        return 3
    cdef size_t charges = charge_clashes ([loop], antigens)
    if charges > 0:
        return 4
    return 0

def framework_antigen_compatibility (list frameworks, list antigens):
    """Check for framework - antigen compatibility"""
    cdef vector[size_t] clashes
    residue_vdw_clashes (frameworks, antigens, clashes, 1)
    if clashes[0] > 0:
        return 1
    elif clashes[1] > 0:
        return 2
    elif clashes[2] > 0:
        return 3
    cdef size_t charges = charge_clashes (frameworks, antigens)
    if charges > 0:
        return 4
    return 0
