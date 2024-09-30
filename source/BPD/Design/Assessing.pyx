# Written in 2018 by Dr. Robert Pantazes. Copyright owned by Auburn
# University.
# This file contains the methods necessary to assess possible solutions to
# binding protein design problems to determine whether or not they are
# satisfactory

# Cython compiler directives
# cython: c_string_type=str
# cython: c_string_encoding=ascii
# cython: language_level=3
# distutils: language=c++

# Include pure python modules
import os
import sys
import time
import math
# Include Binding Protein Design contents
from BPD.BindingProteinError import BindingProteinError
from BPD.Instructions import time_stamp
from BPD.General import parameterize_residue
from BPD.General import loop_antigen_compatibility
from BPD.General import framework_antigen_compatibility
from BPD.Design.DesignClasses cimport Interaction
from BPD.Design.DesignClasses cimport Position
from BPD.Design.DesignClasses cimport Solution
from BPD.Design.PyDesignClasses cimport LoopStructure
from BPD.Design.PyDesignClasses cimport InteractionType
from BPD.Design.PyDesignClasses cimport PySolution
from BPD.Design.PyDesignClasses cimport SolutionStorage
from BPD.Design.Prepare import results_folder_name
# Include the Protein classes
from Proteins.Protein.PyProtein cimport PyProtein
from Proteins.Residue.PyResidue cimport PyResidue
from Proteins.Matrix.PyMatrix cimport PyMatrix
from Proteins.Proteins cimport Protein
from Proteins.Proteins cimport Residue
from Proteins.Proteins cimport Atom
# Include C++ classes
from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp cimport bool

cdef list position_antigens (PySolution solution, list antigens):
    """Position the antigens as they need to be for this solution"""
    # Create the matrices for moving the antigens
    cdef PyMatrix movement = PyMatrix ()
    X = solution.X()
    Y = solution.Y()
    Z = solution.Z()
    movement.PyAtoms_allocation_for_centering([[X, Y, Z]])
    cdef PyMatrix r1 = PyMatrix ()
    r1.specified_rotation(solution.XA(), [1.0, 0.0, 0.0])
    cdef PyMatrix r2 = PyMatrix ()
    r2.specified_rotation(solution.YA(), [0.0, 1.0, 0.0])
    cdef PyMatrix r3 = PyMatrix ()
    r3.specified_rotation(solution.ZA(), [0.0, 0.0, 1.0])
    # This is the list of proteins returned by the function
    cdef list answer = []
    # Go through the antigens
    cdef PyProtein protein, antigen
    for antigen in antigens:
        protein = antigen.duplicate()
        protein.rotate(r1)
        protein.rotate(r2)
        protein.rotate(r3)
        protein.move(movement)
        answer.append(protein)
    return answer

cdef bool preliminary_check (PySolution solution, dict loops, list previous):
    """Check to see whether or not the specified structures are usable"""
    # Information about the positions
    cdef Position * ptr
    cdef size_t J = solution.size()
    cdef str L
    cdef size_t SN, I
    # The answer to this function
    cdef bool answer = True
    # Go through the positions
    for I in range(J):
        # The pointer
        ptr = solution.position(I)
        # The specified loop
        L = ptr.loop()
        # The structure number
        SN = ptr.structure()
        # If that structure is not usable
        if not loops[L].structures[SN-1].usable:
            answer = False
            break
    # If the solution isn't possible because of the loop structures, be done
    if not answer:
        return answer
    # Go through the previously identified solutions
    J = len(previous)
    cdef PySolution other
    for I in range(J):
        other = previous[I]
        # Calculate the similarity score between the two solutions
        if solution.similarity(other) > 0.6:
            answer = False
            break
    return answer

cdef dict collect_loops (instructions, PySolution solution, dict loops, 
                         parameters):
    """Collect the loop structures that are listed in a solution"""
    # The value returned by this function
    cdef dict answer = {}
    # Set up that dictionary
    for value in loops:
        answer[value] = None
    # Go through the positions in the solution
    cdef Position * ptr
    cdef size_t J = solution.size()
    # The loop name, structure number, position index and rotamer number of
    # the interaction
    cdef str L
    cdef size_t SN, PI, RN
    # A looping variable
    cdef size_t I
    for I in range(J):
        # Get the Position's pointer
        ptr = solution.position(I)
        # Get the loop specifier
        L = ptr.loop()
        # If that value is None in the dictionary
        if answer[L] == None:
            # Get the structure number
            SN = ptr.structure()
            # Store the appropriate structure from the loops
            answer[L] = [loops[L].structures[SN-1]]
            # Also store a duplicated copy of the Protein
            answer[L].append(answer[L][0].protein.duplicate())
        # Get the index of the position and the rotamer number
        RN = ptr.rotamer()
        # If it indicates a backbone interaction, skip the rest of this loop
        if RN == 0:
            continue
        # Get the position index
        PI = ptr.position()
        # Load the specific rotamer into the position
        answer[L][1][PI].specific_rotamer(instructions["rotamers"], RN)
        # parameterize that residue
        parameterize_residue(answer[L][1][PI], parameters)
    return answer

cdef void recursive_loop_finding (dict useLoops, dict loops, list unknown, 
        dict acceptable, size_t INDEX, PySolution solution, list antigens):
    """Recursively search for missing binding loop structures"""
    # If the index indicates that all values have been checked
    if INDEX >= len(unknown):
        return
    # Get the next binding loop to search for
    L = unknown[INDEX]
    # Get the number of structures in that binding loop
    cdef size_t I, J
    J = len(loops[L].structures)
    # Find the list of structures from that binding loop that can still be
    # part of a solution
    previous = []
    for I in range(1, J+1):
        # If the structure is usable and known to be acceptable, store the
        # index
        if loops[L].structures[I-1].usable and acceptable[L][I] != False:
            previous.append(I)
    # Now recursively find the values that are compatible with the current
    # solution
    # Go through the use loops information
    for loop in useLoops:
        # If that loop hasn't yet been assigned, skip it
        if useLoops[loop] == None:
            continue
        # Set up an empty list
        structures = []
        # Go through the structures in the target loop
        for I in previous:
            # If this structure is compatible with that one, store the index
            if useLoops[loop][0].compatible(L, I):
                structures.append(I)
        # Replace the previous list
        previous = structures
    # If there are no possible structures, end the search
    if len(structures) == 0:
        return
    # Create a string of this loop's name
    cdef string loopName = L
    # A LoopStructure object
    cdef LoopStructure structure
    # Go through each acceptable structure
    for I in structures:
        # Get the structure
        structure = loops[L].structures[I]
        # Determine whether or not the structure is compatible with the
        # antigens. 
        if acceptable[L][I] == None:
            # Do the calculations
            how = loop_antigen_compatibility (structure.protein, antigens)
            # If it is acceptable, store that so the calculations don't have
            # to be repeated in the future
            if how == 0:
                acceptable[L][I] = True
            # If it isn't acceptable, store that and continue the for loop
            else:
                acceptable[L][I] = False
                continue
        # Create a new dictionary of solutions
        newLoops = {}
        # Store information in it from the use loops
        for loop in useLoops:
            newLoops[loop] = useLoops[loop]
        # Update the information for this loop. Include a 3rd entry so it is
        # clear to the output solution function that this was not a chosen
        # portion of the result
        newLoops[L] = [structure, structure.protein, False]
        # Call the function recursively
        recursive_loop_finding (newLoops, loops, unknown, acceptable, INDEX+1,
                                solution, antigens)
        # Determine if a viable solution has been found or not
        viable = True
        for loop in newLoops:
            if newLoops[loop] == None:
                viable = False
                break
        # If a viable solution has been found, update the use loops and end
        # the function
        if viable:
            for loop in useLoops:
                if useLoops[loop] == None:
                    useLoops[loop] = newLoops[loop]
            return

cdef size_t find_missing_binding_loops (dict useLoops, dict loops,
                                        PySolution solution, list antigens):
    """Find binding loop structures that are missing from the solution"""
    # The answer provided by the function
    cdef size_t answer = 0
    # Identify which loops are missing
    unknown = []
    for L in useLoops:
        if useLoops[L] == None:
            unknown.append(L)
    # If there are no unknown loops, be done
    if len(unknown) == 0:
        return answer
    # Create a dictionary storing information about each possible binding loop
    # structure
    acceptable = {}
    for L in unknown:
        acceptable[L] = {}
        for I in range(len(loops[L].structures)):
            acceptable[L][I+1] = None
    # Use the recursive function to try and find answers
    recursive_loop_finding (useLoops, loops, unknown, acceptable, 0, solution,
                            antigens)
    # If any value in the use loops dictionary is still None, no acceptable
    # answer was found
    for L in useLoops:
        if useLoops[L] == None:
            answer = 1
            break
    return answer

def output_solution (instructions, order, frameworks, useLoops, antigens, AS,
                     summary):
    """Output a successfully identified solution"""
    # Write information to the summary
    message = "Successful solution " + str(AS) + " was identified at "
    message += time_stamp() + "\n"
    summary.write(message)
    summary.flush()
    # The file containing the result is
    fileName = results_folder_name(instructions) + instructions["name"] \
             + "_Solution_" + str(AS) + ".pdb"
    # Open the file for writing
    f = open(fileName, "w")
    # Make sure atoms are numbered sequentially as the analysis progresses
    atomNum = 0
    # Go through the order information
    for pair in order:
        # Get the proper protein piece
        if pair[0]:
            protein = useLoops[pair[1]][1]
        else:
            protein = frameworks[pair[1]]
        # Renumber it's atoms
        atomNum = protein.renumber_atoms(atomNum)
        # Write the structure to the file
        f.write(str(protein))
    # Write the antigens to the file
    for antigen in antigens:
        atomNum = antigen.renumber_atoms(atomNum)
        f.write(str(antigen))
    # End the file
    f.write("END\n")
    f.close()
    # Deactivate the structures that had key interactions as part of this
    # solution
    #
    # This has been disabled because there is now the ability to check and
    # make sure that a solution is at least 40% different than previously
    # accepted solutions
    #for L in useLoops:
    #    if len(useLoops[L]) == 2:
    #        useLoops[L][0].usable = False

def assess_solution (instructions, PySolution solution, loops, order,
                     frameworks, antRef, summary, PS, AS, parameters,
                     previous):
    """Assess a possible solution"""
    # Do a preliminary check to see if the structures are still usable
    cdef bool preliminary = preliminary_check (solution, loops, previous)
    # If not, return the provided indexes directly
    if not preliminary:
        return PS, AS
    # Increment the preliminary solution number
    PS += 1
    # Create a message saying what is happening.
    message = "Assessing preliminary solution " + str(PS) + " started on "
    message += time_stamp() + "This solution includes:\n" + str(solution)
    # Because the code is now fairly well established, summary messages are
    # written all at once and ONLY for successful designs
    #summary.write(message)
    # Collect the necessary loop structures
    cdef dict useLoops = collect_loops (instructions, solution, loops,
                                        parameters)
    # Get the re-positioned antigens
    antigens = position_antigens (solution, antRef)
    # A variable indicating how assessment of a framework or binding loop
    # structure went
    cdef size_t performance
    # Check each binding loop structure
    cdef string loopName
    for L in useLoops:
        if useLoops[L] != None:
            loopName = L
            performance = loop_antigen_compatibility (useLoops[L][1],
                                                      antigens)
            # Write information to the summary file depending on how things
            # went
            if performance == 0:
                message += "The " + L + " binding loop was acceptable\n"
                #summary.write(message)
                continue
            # Any other value indicates failure
            #message += "The " + L + " binding loop was not acceptable because of "
            #if performance == 1:
            #    message += "backbone / backbone VDW clashes\n"
            #elif performance == 2:
            #    message += "backbone / sidechain VDW clashes\n"
            #elif performance == 3:
            #    message += "sidechain / sidechain VDW clashes\n"
            #elif performance == 4:
            #    message += "charge / charge clashes\n"
            #summary.write(message + "\n")
            return PS, AS
    # Check the framework structures
    performance = framework_antigen_compatibility (frameworks, antigens)
    if performance == 0:
        message += "The framework was acceptable\n"
        #summary.write(message)
    else:
        #message = "The framework was not acceptable because of "
        #if performance == 1:
        #    message += "backbone / backbone VDW clashes\n"
        #elif performance == 2:
        #    message += "backbone / sidechain VDW clashes\n"
        #elif performance == 3:
        #    message += "sidechain / sidechain VDW clashes\n"
        #elif performance == 4:
        #    message += "charge / charge clashes\n"
        #summary.write(message + "\n")
        return PS, AS
    # Search for missing binding loop structures
    performance = find_missing_binding_loops (useLoops, loops, solution,
                                              antigens)
    # If the missing loops could not be found
    if performance == 1:
        #message = "The solution failed because no acceptable non-interaction "
        #message += "binding loops could be identified\n\n"
        #summary.write(message)
        return PS, AS
    # Confirm that not more than 3 positively charged residues are used
    posChargeCount = 0
    for pair in order:
        if pair[0]:
            protein = useLoops[pair[1]][1]
            for i in range(len(protein)):
                res = protein[i]
                name = res.name()
                if name in ['ARG', 'LYS']:
                    posChargeCount += 1
    if posChargeCount > 3:
        return PS, AS
    # If the function has reached this point, the solution is viable. 
    # Increment the Acceptable Solution counter
    AS += 1
    # Store this solution in the previous list
    previous.append(solution)
    # Write the message to the summary file
    summary.write(message)
    # Output the solution's complex
    output_solution (instructions, order, frameworks, useLoops, antigens, AS,
                     summary)
    return PS, AS

def find_good_solutions (instructions, solutions, loops, order, frameworks, 
                         antigens, summary, parameters, PS, AS, previous):
    """Find good solutions to the binding protein design problem"""
    # Go through the keys in the order
    for N in solutions.order:
        # Go through the solutions of this number
        for solution in solutions.solutions[N]:
            # assess the solution
            PS, AS = assess_solution (instructions, solution, loops, order,
                         frameworks, antigens, summary, PS, AS, parameters,
                         previous)
            # If the number of acceptable solutions surpasses the target goal,
            # end the function
            if AS >= instructions["designs"]:
                return PS, AS
    return PS, AS
