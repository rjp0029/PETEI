# Written in 2018 by Dr. Robert Pantazes. Copyright owned by Auburn
# University.
# This file contains functions used in getting user-specified instructions for
# designing binding proteins.

# Include the standard python modules
import os
import sys
# Include the BindingProteinError class
from BPD.BindingProteinError import BindingProteinError
# Include the datetime module
import datetime

def time_stamp ():
    """Create a string listing the current time."""
    # get the current time
    a = datetime.datetime.now()
    # Create a summary string
    text = str(a.month) + "/" + str(a.day) + "/" + str(a.year)
    text += " at " + str(a.hour).rjust(2, '0')
    text += ":" + str(a.minute).rjust(2, '0')
    text += ":" + str(a.second).rjust(2, '0') + "\n"
    return text

def load_instructions (fileName):
    """Load the instructions from the user-input file"""
    # Try to open the file
    try:
        f = open(fileName, "r")
    except IOError:
        error = "Failure to open the instructions file: " + fileName
        raise BindingProteinError (error + "\n")
    # Store the contents of the file here, stripped of whitespace
    contents = []
    for line in f:
        contents.append(line.strip())
    # Close the file
    f.close()
    # If the file was empty, raise an error
    if len(contents) == 0:
        error = fileName + " did not contain any content.\n"
        raise BindingProteinError (error)
    # Return those values
    return contents

def get_instructions (fileContents):
    """Get the inputs that give commands from the user"""
    # Store the commands here
    commands = []
    # Loop through the contents of the file
    for line in fileContents:
        # Skip lines that are empty or start with a comment
        if len(line) == 0 or line.startswith("#"):
            continue
        # Store the line
        commands.append(line)
    return commands

def find_instructions (commands, keyPhrase, singular = True):
    """Find instructions that specify specific information"""
    # Store those commands here
    want = []
    # Loop through the commands
    for line in commands:
        # If the line starts with the key phrase + ":"
        if line.upper().startswith(keyPhrase.upper() + ":"):
            # Get the end of the line, stripped of whitespace
            command = line[len(keyPhrase)+1:].strip()
            # If the command is not empty, store it
            if len(command) > 0:
                want.append(command)
    # If no commands were found, raise an error
    if len(want) == 0:
        error = "No " + keyPhrase + " commands were identified.\n"
        raise BindingProteinError(error)
    # If there is only supposed to be a single command, but multiple ones were
    # identified, raise an error
    if singular and len(want) > 1:
        error = "Multiple " + keyPhrase + " commands were identified.\n"
        raise BindingProteinError (error)
    # If there is only supposed to be a single command, return that
    if singular:
        return want[0]
    # Otherwise, return the list
    return want

# The user and name are general commands that provide information about a run
# of the Binding Protein Design software

def get_user (commands):
    """Get the person requesting the calculations"""
    user = find_instructions(commands, "User")
    return user

def get_name (commands, design):
    """Get the name of the calculations."""
    # Identify the name of the calculations
    name = find_instructions(commands, "Name", True)
    # Make sure that the name is a single string
    pieces = name.split()
    if len(pieces) > 1:
        error = "The name for a process may not contain whitespace.\n"
        error += name + " is not acceptable.\n"
        raise BindingProteinError (error)
    # Return the name of the calculations
    return name

# The first step in database generation is an analysis of MD structures. These
# next commands specify that information

def get_structure_folder (commands):
    """Identify the folder where the MD structure results are located"""
    # Get the value
    folder = find_instructions(commands, "MD Folder", True)
    # Make sure it is a single value
    pieces = folder.split()
    if len(pieces) > 1:
        error = "The folder with MD results cannot contain whitespace.\n"
        error += name + " is not acceptable.\n"
        raise BindingProteinError (error)
    # Make sure it is accessible
    current = os.getcwd()
    try:
        os.chdir(folder)
        os.chdir(current)
    except OSError:
        error = "The specified MD folder was not accessible.\n"
        raise BindingProteinError (error)
    # Make sure the folder ends with the proper character
    if not folder.endswith("/"):
        folder += "/"
    return folder

def get_atom_pairs (commands):
    """Identify the atoms that specify the binding loop ends"""
    # Get the values
    pairs = find_instructions(commands, "Atom Pair", False)
    # The results will be stored here
    results = []
    # Loop through the information
    for pair in pairs:
        # split on whitespace
        items = pair.split()
        # Make sure there are 2
        if len(items) != 2:
            error = "Atom Pairs must contain exactly two labels.\n"
            error += pair + " is not acceptable.\n"
            raise BindingProteinError (error)
        # Make sure each value is OK
        for atom in items:
            if atom not in ['N', 'CA', 'C', 'CB', 'O']:
                error = atom + " is not currently an acceptable atom specifier"
                error += " for the start and end points of binding loops.\n"
                raise BindingProteinError (error)
        # Make sure the information is unique
        if items in results:
            error = "There are duplicated " + pair + " Atom pair designations.\n"
            raise BindingProteinError (error)
        # Store the pair of atoms
        results.append(items)
    return results

def get_binding_loops (commands):
    """Identify the information specifying binding loops"""
    # Get the information
    data = find_instructions(commands, "Binding Loop", False)
    # Store the organized in this list
    results = []
    # Loop through the data
    for info in data:
        # Split it into whitespace
        parts = info.split()
        # There must be exactly 4 entries
        if len(parts) != 4:
            error = "This is not a properly formatted binding loop "
            error += "specification: " + info + "\n"
            raise BindingProteinError (error)
        # Make sure it is unique
        if parts in results:
            error = "There are duplicated Binding Loop designations.\n"
            raise BindingProteinError (error)
        results.append(parts)
    return results

# The second step of database generation analyzes the PDB for structures that
# match the geometric features of the binding loops analyzed in Step 1. The
# next set of commands provides that information

def get_loop_lengths (commands):
    """Identify the minimum and maximum permitted binding loop lengths"""
    # Get the information
    data = find_instructions (commands, "Loop Lengths", True)
    # Split it into pieces
    items = data.split()
    if len(items) != 2:
        error = "Exactly 2 Loop Lengths must be specified.\n"
        error += data + " is not acceptable.\n"
        raise BindingProteinError (error)
    # Store the lengths here
    lengths = []
    for value in items:
        try:
            n = int(value)
            if n <= 0:
                raise ValueError
            lengths.append(n)
        except ValueError:
            error = "Loop lengths must be positive integers.\n"
            error += value + " is not acceptable.\n"
            raise BindingProteinError (error)
    # Make sure the second value is greater than or equal to the first by
    # sorting them
    lengths.sort()
    return lengths

def get_PDB_file_list (commands):
    """Identify the file that lists PDB files for analysis"""
    # The name of the file
    fileName = find_instructions(commands, "PDB File List", True)
    # Make sure it has no spaces in it
    pieces = fileName.split()
    if len(pieces) > 1:
        error = "The PDB File List must be a string without whitespace.\n"
        error += fileName + " is not acceptable.\n"
        raise BindingProteinError(error)
    # Make sure the file can be opened
    try:
        f = open(fileName, "r")
        f.close()
    except IOError:
        error = "This file could not be opened: " + fileName + "\n"
        raise BindingProteinError (error)
    return fileName

def get_max_structures (commands):
    """Identify the maximum number of binding loop structures to retain"""
    data = find_instructions (commands, "Maximum Structures", True)
    try:
        n = int(data)
        if n < 1:
            raise ValueError
    except ValueError:
        error = "The Maximum Structures value must be a positive integer.\n"
        error += data + " is not acceptable.\n"
        raise BindingProteinError (error)
    # Throw out a warning if the number is too large
    if n > 5000:
        message = data + " is an unusually large maximum number of binding "
        message += "loop structures.\n This may slow down database creation "
        message += "and subsequent design calculations.\n"
        print(message)
        sys.stdout.flush()
    return n

# The third step of database generation clusters the structures identified in
# Step 2 to identify groups of similar structures. The maximum RMSD they are
# permitted to have must be known

def get_max_RMSD (commands):
    """Get the maximum RMSD allowed during structure clustering"""
    # Get the value
    value = find_instructions (commands, "Clustering RMSD Max", True)
    # Make sure it is a single value
    pieces = value.split()
    if len(pieces) > 1:
        error = "The Clustering RMSD Max value must contain exactly one "
        error += "number.\n" + value + " is not acceptable.\n"
        raise BindingProteinError(error)
    # Make sure it is a positive number
    try:
        answer = float(value)
        if answer <= 0.0:
            raise ValueError
    except ValueError:
        error = "A Clustering RMSD Max value must be a positive number.\n"
        error += value + " is not acceptable.\n"
        raise BindingProteinError (error)
    return answer

# The fourth step of database generation involves positioning the structures,
# removing structures with redundant sequences, patching in rotamers,
# calculating structure - structure energies, and calculating atom distances
# in pairs of acceptable structures. To do all of that requires knowing where
# the Rotamer Library is located and identifying files that contain parameters
# for energy calculations

def get_topology_files (commands):
    """Identify the files that contain CHARMM topology information"""
    # Store them here
    files = find_instructions(commands, "CHARMM Topology File", False)
    # Check each of them
    for value in files:
        # Make sure it is a single string
        pieces = value.split()
        if len(pieces) > 1:
            error = "A CHARMM Topology File may not contain whitespace.\n"
            error += value + " is not acceptable.\n"
            raise BindingProteinError (error)
        # Make sure the file is readable
        try:
            f = open(value, "r")
            f.close()
        except IOError:
            error = "This file could not be opened: " + value + "\n"
            raise BindingProteinError(error)
    return files

def get_parameter_files (commands):
    """Identify the files that contain CHARMM parameter information"""
    # This is effectively identical to the topology file function
    files = find_instructions(commands, "CHARMM Parameter File", False)
    for value in files:
        pieces = value.split()
        if len(pieces) > 1:
            error = "A CHARMM Parameter File may not contain whitespace.\n"
            error += value + " is not acceptable.\n"
            raise BindingProteinError (error)
        try:
            f = open(value, "r")
            f.close()
        except IOError:
            error = "This file could not be opened: " + value + "\n"
            raise BindingProteinError (error)
    return files

def get_solvation_files (commands):
    """Identify files that contain LK Implicit Solvation information"""
    files = find_instructions(commands, "LK Solvation File", False)
    for value in files:
        pieces = value.split()
        if len(pieces) > 1:
            error = "A LK Solvation File may not contain whitespace.\n"
            error += value + " is not acceptable.\n"
            raise BindingProteinError (error)
        try:
            f = open(value, "r")
            f.close()
        except IOError:
            error = "This file coult not be opened: " + value + "\n"
            raise BindingProteinError (error)
    return files

def get_rotamer_library (commands):
    """Identify the location of the rotamer library"""
    # The folder is here
    folder = find_instructions(commands, "Rotamer Library", True)
    # Make sure it does not contain whitespace
    pieces = folder.split()
    if len(pieces) > 1:
        error = "The Rotamer Library specification may not contain whitespace"
        error += ".\n" + folder + " is not acceptable.\n"
        raise BindingProteinError (error)
    # Make sure it can be accessed
    current = os.getcwd()
    try:
        os.chdir(folder)
        os.chdir(current)
    except OSError:
        error = "Unable to access the rotamer library location.\n"
        raise BindingProteinError (error)
    # Make sure the folder ends with the proper character
    if not folder.endswith("/"):
        folder += "/"
    return folder

def load_database_instructions (fileName):
    """Load the instructions for creating a database for protein design."""
    # Get the contents of the file
    allLines = load_instructions(fileName)
    # Get the command lines
    commands = get_instructions(allLines)
    # Access the appropriate information
    instructions = {}
    # General Commands
    instructions["name"] = get_name (commands, False)
    instructions["who"] = get_user (commands)
    # Commands for Step 1
    instructions["MD"] = get_structure_folder (commands)
    instructions["atoms"] = get_atom_pairs (commands)
    instructions["loops"] = get_binding_loops (commands)
    #instructions["rotamers"] = get_rotamer_library(commands)
    instructions["rotamers"] = "./resources/rotamers/"
    # Commands for Step 2
    instructions["PDB"] = get_PDB_file_list (commands)
    instructions["lengths"] = get_loop_lengths (commands)
    instructions["structures"] = get_max_structures (commands)
    instructions["topology"] = get_topology_files (commands)
    instructions["parameters"] = get_parameter_files(commands)
    instructions["solvation"] = get_solvation_files (commands)
    # Commands for Step 3
    instructions["RMSD"] = get_max_RMSD (commands)
    # Give the necessary information back to the calling function
    return instructions

# The next set of functions are those that are required for using a database
# to design binding proteins
def get_database (commands):
    """Identify the database folder"""
    # Get the label
    name = find_instructions(commands, "Database", True)
    # Make sure it is a single value
    pieces = name.split()
    if len(pieces) > 1:
        error = "The database specification cannot contain whitespace.\n"
        error += name + " is not acceptable.\n"
        raise BindingProteinError (error)
    # Try and access it as a folder
    current = os.getcwd()
    try:
        os.chdir("databases/" + name + "/")
        os.chdir(current)
    except OSError:
        error = "The specified database was not accessible.\n"
        raise BindingProteinError (error)
    return name

def get_maximum_designs (commands):
    """Identify how many designs should be created"""
    # Get the information
    value = find_instructions(commands, "Maximum Designs", True)
    # Try converting it to an integer
    try:
        answer = int(value)
        if answer < 1:
            raise ValueError
    except ValueError:
        error = "The Maximum Designs number must be a positive integer.\n"
        error += value + " is not acceptable.\n"
        raise BindingProteinError(error)
    return answer

def get_minimum_interactions (commands):
    """Identify the minimum number of well-positioned interactions"""
    value = find_instructions (commands, "Minimum Interactions", True)
    try:
        answer = int(value)
        if answer < 3:
            raise ValueError
    except ValueError:
        error = "The Minimum Interactions number must be an integer greater " 
        error += "than 3.\n" + value + " is not acceptable.\n"
        raise BindingProteinError (error)
    return answer

def get_antigen (commands):
    """Identify the file that contains the antigen"""
    fileName = find_instructions (commands, "Antigen File", True)
    pieces = fileName.split()
    if len(pieces) > 1:
        error = "The name of the Antigen File must be a string with no "
        error += "whitespace.\n" + fileName + " is not acceptable.\n"
        raise BindingProteinError(error)
    try:
        f = open(fileName,"r")
        f.close()
    except IOError:
        error = "This file could not be opened: " + fileName + "\n"
        raise BindingProteinError (error)
    return fileName

def get_epitope (commands):
    """Get the specifications of the epitope residues"""
    # Get the values
    results = find_instructions(commands, "Epitope Residue", False)
    # If there are not at least five values, raise an error
    if len(results) < 5:
        error = "Binding Protein Design currently requires an epitope of "
        error += "at least 5 residues.\n"
        raise BindingProteinError (error)
    # Store the epitope information here
    epitope = []
    # Go through the data
    for value in results:
        # Make sure the entry has 2 specifications
        pieces = value.split()
        if len(pieces) != 2:
            error = "An Epitope Residue specification should include "
            error += "exactly 2 pieces of data.\n"
            error += value + " is not acceptable.\n"
            raise BindingProteinError (error)
        # The first should be a single character
        if len(pieces[0]) != 1:
            error = "Protein chains in Epitope Residue specifications "
            error += "must be exactly 1 character in length.\n"
            error += pieces[0] + " is not acceptable.\n"
            raise BindingProteinError (error)
        # Make sure the information is unique
        if pieces in epitope:
            error = "Every Epitope Residue specification must be unique.\n"
            raise BindingProteinError(error)
        # Make sure every residue is in the same protein
        for previous in epitope:
            if previous[0] != pieces[0]:
                error = "Currently, every Epitope Residue must be in the "
                error += "same protein chain. This may be changed in a "
                error += "future update.\n"
                raise BindingProteinError (error)
        epitope.append(pieces)
    return epitope

def load_design_instructions (fileName):
    """Load the instructions for designing Binding Proteins"""
    # Get the contents of the file
    allLines = load_instructions (fileName)
    # Get the commands
    commands = get_instructions(allLines)
    # Store the instructions here
    instructions = {}
    # Get the information
    instructions["name"] = get_name (commands, True)
    instructions["who"] = get_user (commands)
    instructions["database"] = get_database(commands)
    instructions["designs"] = get_maximum_designs(commands)
    instructions["interactions"] = get_minimum_interactions(commands)
    instructions["antigen"] = get_antigen(commands)
    instructions["epitope"] = get_epitope(commands)
    instructions["topology"] = get_topology_files(commands)
    instructions["parameters"] = get_parameter_files(commands)
    instructions["solvation"] = get_solvation_files(commands)
    #instructions["rotamers"] = get_rotamer_library (commands)
    instructions["rotamers"] = "./resources/rotamers/"
    return instructions
