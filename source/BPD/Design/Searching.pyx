# Written in 2018 by Dr. Robert Pantazes. Copyright owned by Auburn
# University.
# This file contains the methods necessary to identify possible binding
# protein design solutions given a position of an antigen

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
from BPD.Design.DesignClasses cimport Interaction
from BPD.Design.DesignClasses cimport Position
from BPD.Design.DesignClasses cimport Solution
from BPD.Design.PyDesignClasses cimport LoopStructure
from BPD.Design.PyDesignClasses cimport InteractionType
from BPD.Design.PyDesignClasses cimport SolutionStorage
# Include the Matrix and Protein classes
from Proteins.Matrix.PyMatrix cimport PyMatrix
from Proteins.Protein.PyProtein cimport PyProtein
# Include C++ classes
from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp cimport bool

cdef void find_positions (list interactions, list epitope, 
          vector[Position]& positions):
    """Find all compatible positions given the current epitope"""
    # Looping variables
    cdef size_t I, J
    J = len(interactions)
    # data types
    cdef InteractionType interaction
    # Go through the interaction types
    for I in range(J):
        interaction = interactions[I]
        interaction.find_positions(positions, epitope)

cdef bool compatible_positions (Position * p1, Position * p2, dict loops): 
    """Determine whether or not a pair of positions are compatible"""
    # The answer from the function
    cdef bool answer
    # Do considerations based on the exact positions / rotamers involved
    # If the two positions are in different loop structures, make sure those
    # structures are compatible. Start by getting python str representations
    # of the loop names
    cdef str L1 = p1.loop()
    cdef str L2 = p2.loop()
    # This check needs a Structure
    cdef LoopStructure structure
    if L1 != L2:
        structure = loops[L1].structures[p1.structure()-1]
        answer = structure.compatible(L2, p2.structure())
        return answer
    # At this point, we know the positions are in the same loop. If they
    # aren't in the same structure, they aren't compatible
    if p1.structure() != p2.structure():
        answer = False
        return answer
    # If they are in the same structure, they are compatible if they are at
    # different positions OR are the same rotamer OR either rotamer is a
    # backbone interaction
    if p1.position() == p2.position() and p1.rotamer() != 0 \
       and p2.rotamer() != 0 and p1.rotamer() != p2.rotamer():
        answer = False
        return answer
    # Special consideration is needed for interactions 9 - 12, which are
    # strong hydrogen bonds between ASP / GLU and ARG / LYS. Don't allow those
    # with instances where a single H-bond is occurring.
    if (p1.interaction() <= 8 and (9 <= p2.interaction() <= 12)) or \
       (p2.interaction() <= 8 and (9 <= p1.interaction() <= 12)):
        # If neither interaction is a backbone interaction
        if p1.rotamer() != 0 and p2.rotamer() != 0:
            # If the interactions are with the same residues, they are not
            # compatible
            if p1.position() == p2.position() and p1.residue() == p2.residue():
                answer = False
                return answer
    # If the function has reached this point, the positions are compatible
    answer = True
    return answer

cdef size_t calculate_index (float coordinate, float threshold, float minimum):
    """Calculate an index for the coordinate based on the minimum value"""
    # First calculate it as an integer, because I know how to do that directly
    cdef int value = int((coordinate-minimum)/threshold)
    # Convert that to a size_t
    cdef size_t answer = value
    return answer

cdef void order_positions (vector[Position]& positions, float threshold,
          vector[vector[vector[vector[size_t]]]]& ordered, 
          vector[float]& minimums):
    """Organize the positions based on their X, Y and Z coordinates"""
    # Find the minimum and maximum values
    cdef vector[float] maximums
    # Initialize the values in the vectors
    cdef size_t I1, I2, I3, J1
    for I1 in range(3):
        maximums.push_back(-1e50)
        minimums.push_back(1e50)
    # Go through the positions
    cdef float value
    J1 = positions.size()
    for I1 in range(J1):
        for I2 in range(3):
            value = positions[I1][I2]
            if value > maximums[I2]:
                maximums[I2] = value
            if value < minimums[I2]:
                minimums[I2] = value
    # Determine the number of X, Y and Z bins needed
    cdef size_t X = calculate_index (maximums[0], threshold, minimums[0])+1
    cdef size_t Y = calculate_index (maximums[1], threshold, minimums[1])+1
    cdef size_t Z = calculate_index (maximums[2], threshold, minimums[2])+1
    # Create the vector of size_t values that is needed
    ordered.resize(X)
    for I1 in range(X):
        ordered[I1].resize(Y)
        for I2 in range(Y):
            ordered[I1][I2].resize(Z)
    # Store each position in the proper bin
    for I1 in range(J1):
        X = calculate_index (positions[I1][0], threshold, minimums[0])
        Y = calculate_index (positions[I1][1], threshold, minimums[1])
        Z = calculate_index (positions[I1][2], threshold, minimums[2])
        ordered[X][Y][Z].push_back(I1)

cdef void first_level_compatibility (instructions, 
        vector[vector[size_t]]& solutions, vector[Position]& positions, 
        dict loops) except +:
    """Do the first round of finding compatible movement combinations"""
    # The minimum number of interactions required in a solution
    cdef size_t minNum = instructions["interactions"]
    # If there aren't enough positions to meet that threshold, end the
    # function
    if positions.size() < minNum:
        return 
    # Get the movement threshold from the instructions
    cdef float threshold = instructions["threshold"]
    # Get the ordered position information
    cdef vector[vector[vector[vector[size_t]]]] ordered 
    cdef vector[float] minimums
    order_positions(positions, threshold, ordered, minimums)
    # The number of X, Y and Z bins
    cdef size_t X = ordered.size()
    cdef size_t Y = ordered[0].size()
    cdef size_t Z = ordered[0][0].size()
    # Looping variables
    cdef size_t I1, I2, I3, I4, I5
    # A vector of compatible positions
    cdef vector[size_t] group
    # Pointers to positions
    cdef Position * p1
    cdef Position * p2
    # Go through the X, Y and Z bins
    for I1 in range(X):
        # Information about the next X bin
        L1 = I1+1
        if L1 == X:
            L1 = I1
        for I2 in range(Y):
            # Information about the next Y bin
            L2 = I2 + 1
            if L2 == Y:
                L2 = I2
            for I3 in range(Z):
                # The number of values in this bin
                P = ordered[I1][I2][I3].size()
                # If the bin is empty, skip it
                if P == 0:
                    continue
                # Information about the next Z bin
                L3 = I3 + 1
                if L3 == Z:
                    L3 = I3
                # Go through the positions in this bin
                for I4 in range(P):
                    # Get the pointer to the position
                    p1 = &(positions[ordered[I1][I2][I3][I4]])
                    # Clear the group vector
                    group.clear()
                    # Store the index of this position
                    group.push_back(ordered[I1][I2][I3][I4])
                    # Go through all combinations of bins that might contain
                    # positions that could match with this one
                    for J1 in range(I1, L1+1):
                        for J2 in range(I2, L2+1):
                            for J3 in range(I3, L3+1):
                                # If this is the same X, Y and Z bin, only
                                # search values that come after I4.
                                if J1 == I1 and J2 == I2 and J3 == I3:
                                    # If I4 is the last value in P, skip this
                                    # step
                                    if I4 == P-1:
                                        continue
                                    P2 = I4+1
                                # Otherwise, consider all Positions in the bin
                                else:
                                    P2 = 0
                                P3 = ordered[J1][J2][J3].size()
                                # Go through the positions in that bin
                                for J4 in range(P2, P3):
                                    # Get the pointer to the position
                                    p2 = &(positions[ordered[J1][J2][J3][J4]])
                                    # If the pointers can tell they aren't
                                    # compatible, continue the search
                                    if not p1.compatible(p2, threshold):
                                        continue
                                    # If they are compatible based on whether
                                    # or not the structures are compatible
                                    if compatible_positions (p1, p2, loops):
                                        # Store the index of the position
                                        group.push_back(ordered[J1][J2][J3][J4])
                    # If the group contains enough values, store it in the
                    # solutions
                    if group.size() >= minNum:
                        solutions.push_back(group)
    # Be sure to clear the ordered vector, just to be on the safe side about
    # memory usage
    ordered.clear()

cdef bool previously_seen (vector[size_t]& known, vector[size_t]& novel):
    """Determine if the novel vector is part of the known vector"""
    # The answer from the function starts out as they are the same
    cdef bool same = True
    # For individual values, indicate whether or not they are found
    cdef bool found
    # Looping variables
    cdef size_t I1, I2, J1, J2
    # stored indexes
    cdef size_t n1, n2
    # The lengths of the vectors
    J1 = novel.size()
    J2 = known.size()
    # Go through the novel vector
    for I1 in range(J1):
        # Get the value
        n1 = novel[I1]
        # Indicate that this value has not yet been found
        found = False
        # Go through the known vector
        for I2 in range(J2):
            # Get the value
            n2 = known[I2]
            # If it is a match, indicate that the value was found
            if n1 == n2:
                found = True
                break
        # If the value was not found, the novel vector is in fact novel
        if not found:
            same = False
            break
    # Return whether or not the vector was previously seen
    return same

cdef void store_solution (vector[size_t]& possible, 
          vector[vector[size_t]]& solutions, size_t minNum):
    """Possibly store the provided solution"""
    # If the length of the solution is not long enough, don't store the
    # solution
    if possible.size() < minNum:
        return
    # Whether or not this possible solution has previously been seen
    cdef bool known = False
    # Looping variables
    cdef size_t I1, J1
    J1 = solutions.size()
    # If there are previous solutions to check against
    if J1 > 0:
        # Loop through them
        for I1 in range(J1):
            # See if this possible solution is previously known
            known = previously_seen (solutions[I1], possible)
            # If it is known, stop the search
            if known:
                break
    # If the possible solution is not previously known, store it
    if not known:
        solutions.push_back(possible)

cdef bool check_allowed (vector[size_t]& solution, 
                         vector[vector[bool]]& allowed, size_t INDEX):
    """Determine if a structure is compatible with all previous structures"""
    # The value identified by this function starts out as True
    cdef bool answer = True
    # Looping variables
    cdef I, J
    J = solution.size()
    for I in range(J):
        # Check the value in the allowed vector
        answer = allowed[solution[I]][INDEX]
        # If there is an incompatibility, return that
        if not answer:
            return answer
    return answer

cdef void recursive_solution_finder (vector[size_t]& solution, size_t INDEX,
          vector[vector[bool]]& allowed, vector[vector[size_t]]& solutions,
          size_t minNum):
    """Recursively identify all unique acceptable solutions"""
    # If the index indicates that all values have been checked
    if INDEX >= allowed.size():
        # Try to store the solution and be done
        store_solution (solution, solutions, minNum)
        return
    # Determine whether or not this result can be part of a solution
    cdef bool OK = check_allowed (solution, allowed, INDEX)
    # If it can be allowed, include it in a possible solution
    cdef vector[size_t] newSolution
    cdef size_t I, J
    if OK:
        # The number of current values in the solution
        J = solution.size()
        # Set up the new solution to store that many values + 1
        newSolution.reserve(J+1)
        # Store each value currently in the solution in the new solution
        for I in range(J):
            newSolution.push_back(solution[I])
        # Store this index in the new solution, too
        newSolution.push_back(INDEX)
        # call this function on the new solution and the next index
        recursive_solution_finder (newSolution, INDEX+1, allowed, solutions, 
                                   minNum)
    # Always consider the option of NOT including that value in a solution
    recursive_solution_finder (solution, INDEX+1, allowed, solutions, minNum)

cdef void second_level_compatibility (instructions,
          vector[vector[size_t]]& firstLevel, 
          vector[vector[size_t]]& solutions, 
          vector[Position]& positions, dict loops) except +:
    """Do the second round of finding compatible movement combinations"""
    # The first round found all positions that were compatible with a specific
    # position. Now we need to ensure that the solutions reported contain
    # positions that are all compatible with one another.

    # If there were no previous values, end the function
    if firstLevel.size() == 0:
        return
    # Set up the variables that identify the solutions
    # A vector as the starting point for the recursive search. Every answer
    # must include the first position in the vector
    cdef vector[size_t] initial
    initial.push_back(0)
    # The minimum number of Positions in a viable solution
    cdef size_t minNum = instructions["interactions"]
    # Looping variables
    cdef size_t I1, I2, I3, J1, J2, J3
    # The number of solutions from the first function that need to be checked
    J1 = firstLevel.size()
    # Pointers to Position objects
    cdef Position * p1
    cdef Position * p2
    # They will be used for calculating which positions are compatible with
    # one another and which are not
    cdef vector[vector[bool]] allowed
    cdef bool permitted
    # Temporary solutions identified in a given set of values from the first
    # level of compatibility searching
    cdef vector[vector[size_t]] tempSolutions
    # Variables for converting those temporary solutions into a final solution
    cdef vector[size_t] tempSolution
    # Go through the solutions identified in the first level of compatibility
    # searching
    for I1 in range(J1):
        # Get the number of values in this solution
        J2 = firstLevel[I1].size()
        # Set up the vector that indicates whether or not 2 positions are
        # mutually compatible
        allowed.clear()
        allowed.resize(J2)
        for I2 in range(J2):
            allowed[I2].resize(J2)
            for I3 in range(J2):
                allowed[I2][I3] = True
        # Go through every pair of positions not involving the first - based
        # on the first level compatibility information we know that position
        # is compatible with all the other ones here
        for I2 in range(1, J2-1):
            # Get the pointer to the position
            p1 = &(positions[firstLevel[I1][I2]])
            # Go through the subsequent positions
            for I3 in range(I2+1, J2):
                # Get the pointer
                p2 = &(positions[firstLevel[I1][I3]])
                # The two positions are known to have similar movements
                # because they were both compatible with the first position.
                # Now make sure they aren't replicating the same hydrogen
                # bond.
                permitted = p1.distance_compatible(p2)
                # If that is allowed, check the labelling information
                if permitted:
                    permitted = compatible_positions(p1, p2, loops)
                # Store that information
                allowed[I2][I3] = permitted
                allowed[I3][I2] = permitted
        # Make sure the temporary solution vector is empty
        tempSolutions.clear()
        # Use the recursive function to find the possible solutions from this
        # set of positions
        recursive_solution_finder (initial, 1, allowed, tempSolutions, minNum)
        # If there are some solutions to store
        J2 = tempSolutions.size()
        if J2 > 0:
            # Go through the temporary solutions
            for I2 in range(J2):
                # Clear the temporary solution
                tempSolution.clear()
                # Go through the values in the solution
                J3 = tempSolutions[I2].size()
                for I3 in range(J3):
                    # Store the index of the position value
                    tempSolution.push_back(firstLevel[I1][tempSolutions[I2][I3]])
                # Store the temporary solution in the solutions vector
                solutions.push_back(tempSolution)

cdef void find_movements (instructions, vector[Position]& positions,
          vector[vector[vector[size_t]]]& solutions, dict loops) except +:
    """Find the positions that are compatible with one another"""
    # Do the first level of the search
    #a = time.time()
    cdef vector[vector[size_t]] first 
    first_level_compatibility(instructions, first, positions, loops)
    #b = time.time()
    #message = "The first compatibility search took " + format(b-a, '.3f') 
    #message += " seconds and identified " + str(first.size()) + " combinations"
    #print(message)
    #sys.stdout.flush()
    # Do the second level of the search
    #a = time.time()
    cdef vector[vector[size_t]] second 
    second_level_compatibility(instructions, first, second, positions, loops)
    #b = time.time()
    #message = "The second compatibility search took " + format(b-a, '.3f')
    #message += " seconds and identified " + str(second.size()) + " combinations"
    #print(message)
    #sys.stdout.flush()
    # Clear the first vector now that it is no longer needed
    first.clear()
    # If no solutions were identified, end the function
    if second.size() == 0:
        return
    # Find the maximum and minimum number of solutions
    cdef size_t maxSize = 0
    cdef size_t minSize = 1000000000
    cdef size_t I, J, K
    J = second.size()
    for I in range(J):
        K = second[I].size()
        if K > maxSize:
            maxSize = K
        if K < minSize:
            minSize = K
    # Allocate space in the solutions vector appropriately
    solutions.resize(maxSize - minSize + 1)
    # Store the results
    for I in range(J):
        K = maxSize - second[I].size()
        solutions[K].push_back(second[I])
    # Clear the second vector now that it is no longer needed
    second.clear()

def find_some_solutions (instructions, loops, interactions, epitope, 
                         antigens, float n1, float n2, float n3,
                         SolutionStorage solutions):
    """Search for solutions given a specific rotation of the antigen"""
    # Use a PyMatrix to rotate the antigen
    cdef PyMatrix matrix = PyMatrix()
    # Do the X, Y and Z rotations of the antigen
    matrix.specified_rotation (n1, [1.0, 0.0, 0.0])
    cdef PyProtein antigen
    for antigen in antigens:
        antigen.rotate(matrix)
    matrix.specified_rotation (n2, [0.0, 1.0, 0.0])
    for antigen in antigens:
        antigen.rotate(matrix)
    matrix.specified_rotation (n3, [0.0, 0.0, 1.0])
    for antigen in antigens:
        antigen.rotate(matrix)
    # Find the positions
    #a = time.time()
    cdef vector[Position] positions
    find_positions(interactions, epitope, positions)
    #b = time.time()
    #message = "Finding positions took " + format(b-a, '.3f') + " seconds and "
    #message += str(positions.size()) + " positions were identified"
    #print(message)
    #sys.stdout.flush()
    # Get the combinations of positions that can be part of answers
    cdef vector[vector[vector[size_t]]] combinations
    find_movements(instructions, positions, combinations, loops)
    # Store the solutions
    solutions.store_solutions(positions, combinations, n1, n2, n3)    
    # Reverse the movements of the antigen
    matrix.specified_rotation (-n3, [0.0, 0.0, 1.0])
    for antigen in antigens:
        antigen.rotate(matrix)
    matrix.specified_rotation (-n2, [0.0, 1.0, 0.0])
    for antigen in antigens:
        antigen.rotate(matrix)
    matrix.specified_rotation (-n1, [1.0, 0.0, 0.0])
    for antigen in antigens:
        antigen.rotate(matrix)

def find_possible_solutions (instructions, loops, interactions, epitope,
                             antigens, summary, xa, ya, za):
    """Search positions of the antigen to find possible solutions"""
    # Write a message to the summary file saying what is happening
    message = "The preliminary search for promising antigen positions "
    message += "started on " + time_stamp()
    summary.write(message)
    summary.flush()
    # Create a solution storage object
    cdef SolutionStorage solutions = SolutionStorage ()
    # The angles of rotation
    XA = math.radians(xa)
    YA = math.radians(ya)
    ZA = math.radians(za)
    # Time the process
    #a = time.time()
    # Find some solutions
    find_some_solutions (instructions, loops, interactions, epitope, antigens,
                         XA, YA, ZA, solutions)
    # Write a message to the screen saying how long this took
    #b = time.time()
    #message = "Finding possible solutions took " + format(b-a, '.3f')
    #message += " seconds\n" + str(solutions)
    #print(message)
    #sys.stdout.flush()
    # Update the summary file with what was found
    message = "The preliminary search for solutions ended on " + time_stamp()
    message += str(solutions) + "\n"
    summary.write(message)
    summary.flush()
    return solutions
