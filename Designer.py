# Written in 2018 by Dr. Robert Pantazes. Copyright owned by Auburn
# University. 
# This file contains the code that creates a database of protein pieces for
# designing binding proteins.

# Include the OS and SYS python modules
import os 
import sys
# Include the path to the Protein ANalysis Toolz 
sys.path.append("./source/")
# Include the necessary functions
from BPD.Instructions import load_design_instructions
from BPD.Instructions import time_stamp
from BPD.BindingProteinError import BindingProteinError
from BPD.Design.Prepare import preparation
from BPD.Design.Searching import find_possible_solutions
from BPD.Design.Assessing import find_good_solutions

# Load the instructions
if len(sys.argv) < 2:
    error = "At least one command line argument must be specified.\n"
    raise BindingProteinError (error)
instructions = load_design_instructions(sys.argv[1])
# Set the overwrite value to False
overwrite = False
# If there is a third value and it specifies that existing values should be
# overwritten
if len(sys.argv) >= 3 and sys.argv[2].lower() == "overwrite":
    overwrite = True

# Prepare for the design calculations
summary, order, frameworks, antigens, pdb, epitope, loops, interactions, \
parameters = preparation(instructions, overwrite)

# Set the numbers of possible and accepted solutions to 0
PS = 0
AS = 0
# The previously identified solutions
previous = []

# Rotations Counter
RC = 360

# Rotate around the Z axis
for zI in range(RC):

    message = "Rotation " + str(zI+1) + " of " + str(RC) + "\n"
    summary.write(message)

    za = 360 * float(zI) / RC
    
    # Find possible solutions at this angle of rotation
    solutions = find_possible_solutions (instructions, loops, interactions,
                                  epitope, antigens, summary, 0.0, 0.0, za)

    # Find any good solutions in that set
    PS, AS = find_good_solutions (instructions, solutions, loops, order,
             frameworks, antigens, summary, parameters, PS, AS, previous)

    # If the number of accepted solutions indicates the program should end
    if AS >= instructions["designs"]:
        break

# Write the message to the summary file saying when this finished
message = "The " + instructions["name"] + " Binding Protein Design "
message += "calculations ended on " + time_stamp()
summary.write(message)
summary.close()
