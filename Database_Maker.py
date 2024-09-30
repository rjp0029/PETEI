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
from BPD.Instructions import load_database_instructions
from BPD.BindingProteinError import BindingProteinError
from BPD.Database.Preliminary import Preliminary
from BPD.Database.Framework import Get_Framework
from BPD.Database.LoopSearch import PDB_Search
from BPD.Database.Cluster import Clustering
from BPD.Database.Compatible import Compatibility
from BPD.Database.Interactions import Find_Interactions

# Load the instructions
if len(sys.argv) != 2:
    error = "Exactly one command line argument must be specified.\n"
    raise BindingProteinError (error)
instructions = load_database_instructions(sys.argv[1])
# Do the preliminary database making
Preliminary(instructions)
# Get the framework
#Get_Framework(instructions)
# Search the PDB files for compatible structures
#PDB_Search(instructions)

f = open("./databases/1TTG/1TTG_Summary.txt", "r")
lines = f.readlines()
f.close()

message = "Identifying and positioning the framework copied from 1TTG\n"
for i in range(3, 28):
    message += lines[i].strip() + "\n"
message += "Identifying and positioning the framework ended\n\n"
message += "Searching for binding loop structures copied from 1TTG\n"
message += "H1 Binding Loop: 1000 structures identified\n"
message += "H2 Binding Loop: 1000 structures identified\n"
message += "H3 Binding Loop: 1000 structures identified\n"
message += "Searching for binding loop structures ended\n\n"

f = open("./databases/1TTGu/1TTGu_Summary.txt", "a")
f.write(message)
f.close()

OUTPATH = "./databases/1TTGu/"
INPATH = "./databases/1TTG/"

framework = OUTPATH + "framework/"
os.mkdir(framework)
for term in ['Model', 'framework_1', 'framework_2', 'framework_3', 'framework_4']:
    command = "cp " + INPATH + "framework/1TTG_" + term + ".pdb " + OUTPATH
    command += "framework/1TTGu_" + term + ".pdb"
    os.system(command)

for LOOP in ['H1', 'H2', 'H3']:
    os.mkdir(OUTPATH + LOOP + "/")
    os.mkdir(OUTPATH + LOOP + "/PDB/")
    command = "cp " + INPATH + LOOP + "/PDB/* " + OUTPATH + LOOP + "/PDB/"
    os.system(command)

# Cluster the identified structures into similar groups
Clustering(instructions)
# Determine which binding loop structures are actually able to be parts of
# solutions
Compatibility (instructions)
# Find the positions of the optimal interactions for the database
Find_Interactions (instructions)
