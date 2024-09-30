# Written in 2018 by Dr. Robert Pantazes. Copyright owned by Auburn
# University.
# This file contains the functions and classes needed to assess structures and
# identify the framework for a database of binding protein pieces.

# Include pure Python modules
import os
import sys
import math
import numpy
# Include Binding Protein Design contents
from BPD.BindingProteinError import BindingProteinError
from BPD.General import get_residue
from BPD.General import get_atoms
from BPD.Database.Preliminary import summary_file_name
from BPD.Database.Preliminary import database_folder_name
from BPD.Instructions import time_stamp
# Include the PDB file class and the Matrix class
from Proteins.PDB.PyPDB import PyPDB
from Proteins.Matrix.PyMatrix import PyMatrix

class Loop (object):
    """Information about a specific binding loop"""
    
    def __init__ (self, data, instructions):
        """Initialize a Loop object"""
        # Store the name of the loop
        self.name = data[0]
        # Store the protein the loop belongs to
        self.protein = data[1]
        # The label of the first residue
        self.start = data[2]
        # And the label of the last residue
        self.end = data[3]
        # The distance information will be stored here
        self.distances = []
        # Set up the distances list to store the proper number of values
        for i in range(len(instructions["atoms"])):
            self.distances.append([])

    def __gt__ (self, other):
        """Determine whether or not this loop comes after another"""
        # Get the number of each starting residue
        try:
            n1 = int(self.start)
        except ValueError:
            n1 = int(self.start[:-1])
        try:
            n2 = int(other.start)
        except ValueError:
            n2 = int(other.start[:-1])
        # Return the value of comparing them
        return n1 > n2

    def evaluate (self):
        """Calculate the average and SD of the atom distances"""
        # Do this for every distance pair
        for I in range(len(self.distances)):
            # First calculate the average
            average = 0.0
            for value in self.distances[I]:
                average += value
            average /= float(len(self.distances[I]))
            # Now calculate the standard deviation
            sd = 0.0
            for value in self.distances[I]:
                sd += math.pow(average - value, 2)
            sd = math.sqrt(sd / float(len(self.distances[I])))
            # Store that information in the binding loop
            self.distances[I] = [average, sd]

def prepare_loops (instructions):
    """Generate ordered loop information"""
    # Store the loops here
    preliminary = []
    # Go through the relevant information and make the loops
    for data in instructions["loops"]:
        preliminary.append(Loop(data, instructions))
    # Store them by the protein chain they belong to
    chainOrder = []
    chains = {}
    for loop in preliminary:
        if loop.protein not in chainOrder:
            chainOrder.append(loop.protein)
            chains[loop.protein] = []
        chains[loop.protein].append(loop)
    # Store the final loops here
    loops = []
    # Go through the protein chains
    for L in chainOrder:
        # sort the loops in this chain
        if len(chains[L]) > 1:
            done = False
            while not done:
                done = True
                for i in range(len(chains[L]) - 1):
                    if chains[L][i] > chains[L][i+1]:
                        done = False
                        temp = chains[L][i]
                        chains[L][i] = chains[L][i+1]
                        chains[L][i+1] = temp
        # Store the loops
        for loop in chains[L]:
            loops.append(loop)
    return loops

def get_atom_sets (instructions):
    """Get the lists of atoms for the starts and ends of binding loops"""
    # Store the atoms that are used at the start of the binding loops here
    startAtoms = []
    # And those at the end here
    endAtoms = []
    # Go through the information
    for pair in instructions["atoms"]:
        startAtoms.append(pair[0])
        endAtoms.append(pair[1])
    return startAtoms, endAtoms

def calculate_loop_distances (instructions, loops, fileNames):
    """Calculate the loop distance information"""
    # Get the sets of atoms that need to be used
    startNames, endNames = get_atom_sets(instructions)
    # Go through every file in the folder
    for fileName in fileNames:
        # Try to load the PDB file
        try:
            pdb = PyPDB (fileName, instructions["MD"])
        # If there was an error, tell the user
        except RuntimeError:
            error = "Error loading PDB file: " + fileName + "\n"
            error += "Most likely, not all files in the MD folder are PDB "
            error += "formatted\n"
            raise BindingProteinError (error)
        # Go through the loops
        for loop in loops:
            # Find the proper protein in the PDB file
            protein = None
            for i in range(pdb.proteins()):
                prot = pdb.protein(i)
                if prot.name() == loop.protein:
                    protein = prot
                    break
            # If the protein was not found, raise an error
            if protein == None:
                error = "Failure to find Protein " + loop.protein + " in PDB "
                error += "file " + fileName + "\n"
                raise BindingProteinError (error)
            # Get the necessary atoms for the distance calculations
            try:
                start = get_residue(protein, loop.start)
                end = get_residue(protein, loop.end)
                startAtoms = get_atoms(start, startNames)
                endAtoms = get_atoms(end, endNames)
            # If that failed, raise an error
            except BindingProteinError as e:
                error = e.error + "This occurred in file " + fileName + "\n"
                raise BindingProteinError (error)
            # Do the distance calculations
            for i in range(len(startAtoms)):
                D = startAtoms[i].calculate_distance(endAtoms[i])
                loop.distances[i].append(D)
    # Calculate the average and standard deviations for the loop distance
    # information
    for loop in loops:
        loop.evaluate()

def summarize_loop_information (instructions, loops, summary):
    """Summarize the calculated loop parameters"""
    # Store the message here
    message = ""
    # Go through the binding loops
    for loop in loops:
        message += "Binding Loop " + loop.name + " in Protein " + loop.protein
        message += "\n"
        # Go through the atom pairs
        for i, pair in enumerate(instructions["atoms"]):
            message += pair[0] + " - " + pair[1] + " distance: "
            message += format(loop.distances[i][0], '.5f') + " +/- "
            message += format(loop.distances[i][1], '.5f') + "\n"
    # Write the message to the summary file
    summary.write(message)
    summary.flush()

def select_model_structure (instructions, loops, fileNames):
    """Select the structure that will serve as the model framework"""
    # Get the sets of atoms that are the start and ends of the binding loops
    startNames, endNames = get_atom_sets (instructions)
    # The best PDB file goes here
    best = None
    # And the best score goes here
    bestScore = None
    # Go through the files. Unlike the calculate loop distances function, this
    # function does not error check the necessary information because any
    # errors will have been detected in the calculate loop distances function
    for fileName in fileNames:
        # Get the PDB file object
        pdb = PyPDB (fileName, instructions["MD"])
        # The deviation or score for this file
        deviation = 0.0
        # Go through the binding loops
        for loop in loops:
            # Find the protein
            for i in range(pdb.proteins()):
                protein = pdb.protein(i)
                if protein.name() == loop.protein:
                    break
            # Get the atoms
            startRes = get_residue(protein, loop.start)
            endRes = get_residue(protein, loop.end)
            startAtoms = get_atoms(startRes, startNames)
            endAtoms = get_atoms(endRes, endNames)
            # Go through each atom pair
            for i in range(len(startAtoms)):
                # Calculate the distance between the atoms
                D = startAtoms[i].calculate_distance(endAtoms[i])
                # Calculate the deviation of the distance with the average
                deviation += math.pow(loop.distances[i][0] - D, 2)
        # Possibly store this result
        if best == None or deviation < bestScore:
            best = pdb
            bestScore = deviation
    # Return the selected model complex
    return best

def calculate_best_fit_plane (atoms):
    """Calculate a best fit plane for a set of Atoms"""
    # This code is based on the discussion and sample Python code for this
    # problem found at:
    # math.stackexchange.com/questions/99299/best-fitting-plane-given-a-set-of-points?utm_medium=organic&utm_source=google_rich_qa&utm_campaign=google_rich_qa
    # Collect the coordinates in two lists
    tmpA = []
    tmpB = []
    for atom in atoms:
        tmpA.append([atom[0], atom[1], 1.0])
        tmpB.append(atom[2])
    # Make numpy matrices
    B = numpy.matrix(tmpB).T
    A = numpy.matrix(tmpA)
    # Calculate the best fit plane
    fit = (A.T * A).I * A.T * B
    return fit

def calculate_average_loop_Z(proteins, loops):
    """Calculate the average Z coordinate of the atoms in the binding loops"""
    # Store the average here
    average = 0.0
    # The number of contributing atoms
    used = 0
    # Go through the binding loops
    for loop in loops:
        # Go through the proteins
        for protein in proteins:
            # If it isn't the right protein, skip it
            if protein.name() != loop.protein:
                continue
            # Use this flag to indicate whether or not to include the atoms of
            # the residue
            include = False
            # go through the protein's residues
            for i in range(len(protein)):
                # Get the residue
                res = protein[i]
                # If this residue is the start of the binding loop
                if res.get_label() == loop.start:
                    include = True
                # If the residue's atoms should be included in the average
                if include:
                    for j in range(len(res)):
                        average += res[j][2]
                        used += 1
                # If this is the last residue in the binding loop
                if res.get_label() == loop.end:
                    include = False
                    break
    # Return the average value
    return average / float(used)

def prepare_model (instructions, loops, model):
    """Prepare the model for use in a database for designing proteins"""
    # Get the proteins that hold the framework
    proteins = []
    chains = []
    for loop in loops:
        if loop.protein not in chains:
            for i in range(model.proteins()):
                protein = model.protein(i)
                if protein.name() == loop.protein:
                    proteins.append(protein)
                    chains.append(loop.protein)
                    break
    # Assign rotamers to these proteins
    for protein in proteins:
        for i in range(len(protein)):
            protein[i].closest_rotamer(instructions["rotamers"])
    # Collect all of the atoms that are used as the attachment points for the
    # binding loops
    atoms = []
    #startNames, endNames = get_atom_sets(instructions)
    for loop in loops:
        for protein in proteins:
            if protein.name() == loop.protein:
                startRes = get_residue(protein, loop.start)
                atoms.extend(get_atoms(startRes, ['N', 'CA', 'C', 'CB']))
                endRes = get_residue(protein, loop.end)
                atoms.extend(get_atoms(endRes, ['N', 'CA', 'C', 'CB']))
    # Create a PyMatrix for centering those atoms around the origin
    matrix = PyMatrix ()
    matrix.PyAtoms_allocation_for_centering(atoms)
    # Move the model
    for protein in proteins:
        protein.move(matrix)
    # Testing demonstrated that this approach to orienting the best fit plane
    # required multiple iterations to get to a reasonable answer. There is
    # therefore probably a more efficient way to do this, but this approach is
    # effective
    for I in range(4):
        fit = calculate_best_fit_plane(atoms)
        # Calculate the angle of rotation around the Y axis so that the X
        # contribution becomes negligible
        angle = math.fabs(math.acos(1.0 / 
                math.sqrt(math.pow(fit[0][0], 2) + 1)))
        if fit[0][0] < 0:
            angle = -angle
        # Create the rotation matrix to do this
        matrix.specified_rotation(angle, [0.0, 1.0, 0.0])
        # Rotate the proteins
        for protein in proteins:
            protein.rotate(matrix)
        # Calculate the new best fit plane
        fit = calculate_best_fit_plane(atoms)
        # Calculate the angle of rotation around the X axis so that the Y
        # contribution becomes negligible
        angle = math.fabs(math.acos( 1.0 / 
                math.sqrt(math.pow(fit[1][0], 2) + 1.0)))
        if fit[1][0] > 0:
            angle = -angle
        matrix.specified_rotation(angle, [1.0, 0.0, 0.0])
        for protein in proteins:
            protein.rotate(matrix)
    # Calculate the average Z coordinate of the loop atoms
    Z = calculate_average_loop_Z (proteins, loops)
    # If that value is negative, rotate the system 180 around the X-axis
    if Z < 0:
        angle = math.radians(180.0)
        matrix.specified_rotation(angle, [1.0, 0.0, 0.0])
        for protein in proteins:
            protein.rotate(matrix)

def output_framework(instructions, loops, model, summary):
    """Output the framework for the database"""
    # First, create the framework folder
    folder = database_folder_name(instructions) + "framework/"
    # Make the folder
    os.mkdir(folder)
    # Collect the proteins that contain the framework
    proteins = []
    order = []
    for loop in loops:
        if loop.protein not in order:
            order.append(loop.protein)
            for i in range(model.proteins()):
                protein = model.protein(i)
                if protein.name() == loop.protein:
                    proteins.append(protein)
                    break
    # Output those proteins to a file
    atomNum = 1
    f = open(folder + instructions["name"] + "_Model.pdb", "w")
    for protein in proteins:
        atomNum = protein.renumber_atoms(atomNum)
        f.write(str(protein))
    f.write("END\n")
    f.close()
    # Get the framework pieces for the complex
    frameworks = []
    framework = []
    LoopFlag = False
    # Go through the proteins
    for protein in proteins:
        # Go through the residues
        for i in range(len(protein)):
            res = protein[i]
            # Determine if this residue is the start of a binding loop
            for loop in loops:
                if loop.protein == protein.name() and \
                   loop.start == res.get_label():
                       LoopFlag = True
                       if len(framework) > 0:
                           frameworks.append(framework)
                           framework = []
                       break
            # If the residue is not part of a loop, include it in the
            # framework
            if not LoopFlag:
                framework.append(res)
            # Determine if the residue is the end of a binding loop
            for loop in loops:
                if loop.protein == protein.name() and \
                   loop.end == res.get_label():
                       LoopFlag = False
        # If there is a framework, store it
        if len(framework) > 0:
            frameworks.append(framework)
            framework = []
    # Output those frameworks
    for I in range(1, len(frameworks)+1):
        # The name of the file
        fileName = folder + instructions["name"] + "_framework_" + str(I)
        fileName += ".pdb"
        # Open that file
        f = open(fileName, "w")
        # Write the framework residues
        for res in frameworks[I-1]:
            f.write(str(res))
        f.close()
    # Tell the summary file how many framework pieces there are
    message = "Framework Pieces Identified: " + str(len(frameworks)) + "\n"
    summary.write(message)
    summary.flush()

def Get_Framework (instructions):
    """Get the Framework for the database being created"""
    # Get the name of the database folder
    folder = database_folder_name (instructions)
    # Get the name of the summary file
    summaryName = summary_file_name (instructions)
    # Open that file for writing, but append at the end of the file
    summary = open(summaryName, "a")
    # Write that this is starting
    message = "Identifying and positioning the framework started on "
    message += time_stamp()
    summary.write(message)
    summary.flush()
    # Prepare the Binding Loops
    loops = prepare_loops(instructions)
    # Get the names of the PDB files
    try:
        fileNames = os.listdir(instructions["MD"])
    except OSError:
        error = "Error accessing the MD folder: " + instructions["MD"] + "\n"
        raise BindingProteinError (error)
    # Calculate the loop parameter information
    calculate_loop_distances(instructions, loops, fileNames)
    # Summarize the binding loop information
    summarize_loop_information(instructions, loops, summary)
    # Find the model complex
    model = select_model_structure (instructions, loops, fileNames)
    # Prepare the model by patching in rotamers and properly positioning it
    prepare_model (instructions, loops, model)
    # Output the model and framework pieces
    output_framework(instructions, loops, model, summary)
    # Tell the summary file that this task is finished
    message = "Identifying and positioning the framework ended on "
    message += time_stamp() + "\n"
    summary.write(message)
    summary.close()
