# Written in 2018 by Dr. Robert Pantazes. Copyright owned by Auburn
# University. 
# This file contains functions that set up a folder for storing a database of
# protein pieces that can be assembled into a binding protein structure

# Include Python modules
import os
import sys
# Include the error class
from BPD.BindingProteinError import BindingProteinError
# Include the ability to list the time something occurs
from BPD.Instructions import time_stamp

def database_folder_name (instructions):
    """The name of the folder where the database is being stored"""
    return "databases/" + instructions["name"] + "/"

def create_database_folder (instructions):
    """Make the folder to store the database"""
    # The name of the folder
    folder = database_folder_name(instructions)
    # Try to make that folder
    try:
        os.mkdir(folder)
    # If that failed, raise an error
    except OSError:
        error = "Failure to create the database folder. It likley already " \
              + "exists:\n" + folder + "\n"
        raise BindingProteinError (error)

def make_instructions_file (instructions):
    """Make a file listing the user-specified instructions"""
    # Store the information here
    output = "# General Information\n"
    output += "User: " + instructions["who"] + "\n"
    output += "Name: " + instructions["name"] + "\n\n"
    output += "# Information for assessing the framework protein(s)\n"
    output += "MD Folder: " + instructions["MD"] + "\n"
    for pair in instructions["atoms"]:
        output += "Atom Pair: " + pair[0] + " " + pair[1] + "\n"
    output += "\n"
    for loop in instructions["loops"]:
        output += "Binding Loop:"
        for value in loop:
            output += " " + value
        output += "\n"
    output += "\nRotamer Library: " + instructions["rotamers"] + "\n\n"
    output += "# Information for searching PDB structures for binding loops\n"
    output += "PDB File List: " + instructions["PDB"] + "\n"
    output += "Loop Lengths: " + str(instructions["lengths"][0]) + " "
    output += str(instructions["lengths"][1]) + "\n"
    output += "Maximum Structures: " + str(instructions["structures"]) + "\n"
    for entry in instructions["topology"]:
        output += "CHARMM Topology File: " + entry + "\n"
    for entry in instructions["parameters"]:
        output += "CHARMM Parameter File: " + entry + "\n"
    for entry in instructions["solvation"]:
        output += "LK Solvation File: " + entry + "\n"
    output += "\n# Information for clustering binding loop structures\n"
    output += "Clustering RMSD Max: " + str(instructions["RMSD"]) + "\n"
    # Create the file
    fileName = database_folder_name (instructions) + instructions["name"]
    fileName += "_Instructions.txt"
    f = open(fileName, "w")
    f.write(output)
    f.close()

def summary_file_name (instructions):
    """The name of the summary file"""
    fileName = database_folder_name(instructions) + instructions["name"]
    fileName += "_Summary.txt"
    return fileName

def start_summary_file (instructions):
    """Start the summary file"""
    # Get the name of the file
    fileName = summary_file_name (instructions)
    # Open the file for writing
    f = open(fileName, "w")
    # Write this message
    message = "Creating the " + instructions["name"] + " binding protein "
    message += "database started on " + time_stamp() + "\n"
    f.write(message)
    # Actually close the file - the individual steps are time consuming and it
    # is easier to debug / re-work this program if each step independently
    # deals with the file
    f.close()

def Preliminary (instructions):
    """Do the preliminary preparation for the database creation"""
    create_database_folder (instructions)
    make_instructions_file (instructions)
    start_summary_file (instructions)
