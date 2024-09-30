# Written in 2018 by Dr. Robert Pantazes. Copyright owned by Auburn
# University.
# This file contains the functions necessary for doing all the setup /
# preparation calculations for a Binding Protein Design problem

# Include standard Python modules
import os
import sys
import math
import numpy
# Include Protein structures
from Proteins.PDB.PyPDB import PyPDB
from Proteins.Matrix.PyMatrix import PyMatrix
# Include the error class
from BPD.BindingProteinError import BindingProteinError
# Include the Loop and Interaction Type classes
from BPD.Design.PyDesignClasses import Loop
from BPD.Design.PyDesignClasses import InteractionType
# Include the function that gets interaction parameters
from ForceFields.Parameters import load_parameters
# And the ability to parameterizer proteins
from BPD.General import parameterize_protein
# And the time stamp function
from BPD.Instructions import time_stamp
from Proteins.Protein.PyProtein import PyProtein

def results_folder_name (instructions):
    """Get the name of the folder where the results are stored"""
    return "./results/" + instructions["name"] + "/"

def database_folder_name (instructions):
    """The name of the folder where the database is located"""
    return "./databases/" + instructions["database"] + "/"

def get_database_summary (instructions):
    """Get the summary of what was done when the database was made"""
    # The name of the summary file
    fileName = database_folder_name(instructions) + instructions["database"]
    fileName += "_Summary.txt"
    # Open that file, raising an error if there is an issue
    try:
        f = open(fileName, "r")
    except IOError:
        error = "Failure to load the " + instructions["database"] 
        error += " summary file\n"
        raise BindingProteinError (error)
    # Read in the lines from that file
    lines = f.readlines()
    # close the file and return the liens
    f.close()
    return lines

def get_database_information (instructions):
    """Get information about the database"""
    # Get the lines from the database summary file
    lines = get_database_summary (instructions)
    # Determine the number of framework pieces, the order of the binding
    # loops, the protein they belong to, and the number of structures per
    # binding loop
    loops = []
    FN = None
    for line in lines:
        # The binding loop specifications and order
        if line.startswith("Binding Loop"):
            items = line.split()
            loops.append([items[2], items[5]])
        # The number of framework pieces
        elif line.startswith("Framework Pieces"):
            FN = int(line.split()[3])
        # The number of structures for each binding loop
        elif "usable" in line:
            items = line.split()
            L = items[0]
            N = int(items[3])
            for i in range(len(loops)):
                if loops[i][0] == L:
                    loops[i].append(N)
                    break
    # Determine the order of the framework pieces and binding loops
    Order = []
    # Include the first framework piece
    Order.append([False, 0])
    FI = 0
    # Go through the binding loops
    for i in range(len(loops)):
        # If this is not the first binding loop and it is in a different
        # protein than the previous loop, include the next framework structure
        if i > 0 and loops[i][1] != loops[i-1][1]:
            FI += 1
            Order.append([False, FI])
        # Include the Binding loop 
        Order.append([True, loops[i][0]])
        # Include the next framework structure
        FI += 1
        Order.append([False, FI])
    # Create a dictionary that lists the indexes of each binding loop
    LoopInfo = {}
    for i, loop in enumerate(loops):
        LoopInfo[loop[0]] = [i, loop[2]]
    # Return this information
    return FN, Order, LoopInfo

def load_frameworks (instructions, FN, parameters):
    """Load and parameterize the framework pieces"""
    # Get the database folder
    folder = database_folder_name(instructions) + "framework/"
    # Store the frameworks in this list
    frameworks = []
    # Loop through the structure numbers
    for I in range(1, FN+1):
        # The name of the file
        fileName = instructions["database"] + "_framework_" + str(I) + ".pdb"
        # Load the PyProtein
        framework = PyProtein (fileName, folder)
        # parameterize it
        parameterize_protein(framework, parameters)
        # Store it
        frameworks.append(framework)
    return frameworks

def load_antigens (instructions, parameters):
    """Load the antigen structures"""
    # Currently, the Instructions code only permits a single protein chain to
    # be an antigen. This code will need to be changed in the future if that
    # ever changes, but this function already returns a list of PyProteins.
    # The rest of the program won't have to be modified if multiple antigen
    # proteins are permitted in the future
    antigens = []
    # Load the antigen file
    pdb = PyPDB(instructions["antigen"], "")
    # Get the chains in the epitope
    chains = []
    for res in instructions["epitope"]:
        if res[0] not in chains:
            chains.append(res[0])
    # Get the protein chains
    for L in chains:
        flag = False
        for i in range(pdb.proteins()):
            protein = pdb.protein(i)
            if protein.name() == L:
                flag = True
                antigens.append(protein)
                break
        if not flag:
            error = "Failure to find chain " + L + " in PDB file "
            error += instructions["antigen"] + "\n"
            raise BindingProteinError (error)
    # Parameterize every protein
    for protein in antigens:
        parameterize_protein(protein, parameters)
    return antigens, pdb

def get_epitope (instructions, antigens):
    """Get the epitope residues"""
    # Store them here
    epitope = []
    # Go through the epitope information
    for pair in instructions["epitope"]:
        # Set the labelling information for the epitope residue
        L = pair[0]
        res = pair[1]
        # Whether or not the residue is found
        Found = False
        # Go through the antigens
        for antigen in antigens:
            # If this is not the right protein, skip it
            if antigen.name() != L:
                continue
            # Go through the residues
            for i in range(len(antigen)):
                # Get the residue
                residue = antigen[i]
                # If this is the right residue
                if residue.get_label() == res:
                    Found = True
                    epitope.append(residue)
                    break
            if Found:
                break
        # If the residue was not found, raise an error
        if not Found:
            error = "Failure to find epitope residue " + res + " in protein "
            error += L + "\n"
            raise BindingProteinError (error)
    # Return the epitope residues
    return epitope

def calculate_best_fit_plane (atoms):
    """Calculate a best fit plane for a set of Atoms"""
    # Refer to the same function in the Database methods for an explanation
    tmpA = []
    tmpB = []
    for atom in atoms:
        tmpA.append([atom[0], atom[1], 1.0])
        tmpB.append(atom[2])
    B = numpy.matrix(tmpB).T
    A = numpy.matrix(tmpA)
    fit = (A.T * A).I * A.T * B
    return fit

def position_antigen (instructions, antigens):
    """Position the antigen so the epitope points down"""
    # Store all epitope atoms here
    epitope = []
    # And all non-epitope atoms here
    other = []
    # go through the antigen chains
    for protein in antigens:
        # Go through the residues
        for I in range(len(protein)):
            residue = protein[I]
            # If this residue is part of the epitope
            if [residue.protein(), residue.get_label()] in instructions["epitope"]:
                for J in range(len(residue)):
                    epitope.append(residue[J])
            else:
                for J in range(len(residue)):
                    other.append(residue[J])
    # Move the complex so that the epitope atoms are centered around the
    # origin
    matrix = PyMatrix ()
    matrix.PyAtoms_allocation_for_centering(epitope)
    for protein in antigens:
        protein.move(matrix)
    # If there are more epitope atoms than non-epitope atoms, assume the user
    # properly positioned their antigen to begin with
    if len(other) < len(epitope):
        return
    # Experience has shown that the following steps need to be repeated a few
    # times to get a good fit
    for I in range(4):
        fit = calculate_best_fit_plane(epitope)
        angle = math.fabs(math.acos(1.0 /
                math.sqrt(math.pow(fit[0][0], 2) + 1)))
        if fit[0][0] < 0:
            angle = -angle
        matrix.specified_rotation (angle, [0.0, 1.0, 0.0])
        for protein in antigens:
            protein.rotate(matrix)
        fit = calculate_best_fit_plane(epitope)
        angle = math.fabs(math.acos(1.0 /
                math.sqrt(math.pow(fit[1][0], 2) + 1.0)))
        if fit[1][0] > 0:
            angle = -angle
        matrix.specified_rotation(angle, [1.0, 0.0, 0.0])
        for protein in antigens:
            protein.rotate(matrix)
    # Calculate the average Z coordinate of the epitope atoms
    EZ = 0.0
    for atom in epitope:
        EZ += atom[2]
    EZ /= float(len(epitope))
    # Calculate the average Z coordinate of the other atoms
    OZ = 0.0
    for atom in other:
        OZ += atom[2]
    OZ /= float(len(other))
    # If the average Z coordinate of the epitope atoms is greater than that of
    # the other atoms, flip the system
    if EZ > OZ:
        matrix.specified_rotation (math.radians(180.0), [0.0, 1.0, 0.0])
        for protein in antigens:
            protein.rotate(matrix)
    # For visualization purposes, raise the antigens up by 20 angstroms. 
    matrix.PyAtoms_allocation_for_centering([[0.0, 0.0, 20.0]])
    for protein in antigens:
        protein.move(matrix, False)

def summary_file_name (instructions):
    """Get the name of the summary file"""
    fileName = results_folder_name (instructions) + instructions["name"]
    fileName += "_Summary.txt"
    return fileName

def preparation (instructions, overwrite = False):
    """Prepare to do the Design Calculations"""
    # I don't yet want to allow the distance threshold to be a user-settable
    # value, but I might someday. For now, store the value in the instructions
    # as this function starts
    instructions["threshold"] = 1.0 / 3.0
    # Create a message saying what is happening. It will be written to the
    # summary file once that is created. It starts with the instructions for
    # the calculations
    message = "Name: " + instructions["name"] + "\n"
    message += "User: " + instructions["who"] + "\n\n"
    message += "Database: " + instructions["database"] + "\n"
    message += "Maximum Designs: " + str(instructions["designs"]) + "\n"
    message +="Minimum Interactions: "+str(instructions["interactions"])+"\n\n"
    message += "Antigen File: " + instructions["antigen"] + "\n"
    for pair in instructions["epitope"]:
        message += "Epitope Residue: " + pair[0] + " " + pair[1] + "\n"
    message += "\nRotamer Library: " + instructions["rotamers"] + "\n"
    for data in instructions["topology"]:
        message += "CHARMM Topology File: " + data + "\n"
    for data in instructions["parameters"]:
        message += "CHARMM Parameter File: " + data + "\n"
    for data in instructions["solvation"]:
        message += "LK Solvation File: " + data + "\n"
    message += "\nThe " + instructions["name"] + " Binding Protein Design "
    message += "calculations started on " + time_stamp()
    # Create the folder for the results
    folder = results_folder_name (instructions)
    # Try to make that folder, and raise an error if it already exists
    try:
        os.mkdir(folder)
    except OSError:
        if overwrite:
            command = "rm -rf " + folder
            os.system(command)
            os.mkdir(folder)
        else:
            error = "There is already a " + instructions["name"] + " results "
            error += "folder. Select a different name.\n"
            raise BindingProteinError (error)
    # Create the summary file
    summary = open(summary_file_name (instructions), "w")
    summary.write(message)
    summary.flush()
    # Load the energy calculation parameters
    parameters = load_parameters (instructions["topology"],
                                  instructions["parameters"],
                                  instructions["solvation"])
    # Load the information about the database
    FN, Order, LoopInfo = get_database_information (instructions)
    # Load the framework structures
    frameworks = load_frameworks (instructions, FN, parameters)
    # Load the antigens
    antigens, pdb = load_antigens (instructions, parameters)
    # Move them so the epitope is pointed down at the binding pocket
    position_antigen (instructions, antigens)
    # Get the epitope residues
    epitope = get_epitope (instructions, antigens)
    # Get the name of the database's folder
    DBFolder = database_folder_name(instructions)
    # Create the loops, which loads the structures
    loops = {}
    for label in LoopInfo:
        loops[label] = Loop(label, instructions, DBFolder, LoopInfo, parameters)
    # Load the Interactions
    interactions = []
    for I in range(1, 24):
        if I == 3 or I == 4:
            continue
        interactions.append(InteractionType(instructions, I, DBFolder))
    # Indicate when the database information stopped being loaded into memory
    message = "Loading database information ended on " + time_stamp() + "\n"
    summary.write(message)
    summary.flush()
    # Return the information generated by this function
    return summary, Order, frameworks, antigens, pdb, epitope, loops, \
           interactions, parameters
