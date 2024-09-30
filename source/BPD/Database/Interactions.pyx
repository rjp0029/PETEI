# Written in 2018 by Dr. Robert Pantazes. Copyright owned by Auburn
# University.
# This file contains the functions needed to find optimal / near optimal
# binding interactions between proteins, as used in designing binding proteins

# Cython compiler directives
# cython: c_string_type=str
# cython: c_string_encoding=ascii
# cython: language_level=3
# distutils: language=c++

# Include pure Python modules
import os 
import sys
import math
# Include Binding Protein Design methods
from BPD.BindingProteinError import BindingProteinError
from BPD.Instructions import time_stamp
from BPD.Database.Preliminary import database_folder_name
from BPD.Database.Preliminary import summary_file_name
from BPD.General import get_atoms
# Include Protein structure classes
from Proteins.Matrix.PyMatrix cimport PyMatrix
from Proteins.Proteins cimport Matrix
from Proteins.Protein.PyProtein cimport PyProtein
from Proteins.Residue.PyResidue cimport PyResidue
from Proteins.Atom.PyAtom cimport PyAtom
from Proteins.Proteins cimport Atom
# Include C++ classes
from libcpp cimport bool
from libcpp.vector cimport vector

cdef void do_positioning_movement (PyResidue residue, PyMatrix matrix, 
                                   bool subtract = True):
    """Do or reverse an interaction positioning movement"""
    # If the matrix indicates that a movement should be done
    if matrix._thisptr.rows() != matrix._thisptr.columns():
        residue.move(matrix, subtract)
    # Otherwise, do a rotation
    else:
        residue.rotate(matrix)

cdef void reverse_positioning_movements (PyResidue residue, list movements):
    """Reverse a set of interaction positioning movements"""
    # Looping variables
    cdef size_t I, J
    J = len(movements)
    # indicate that any translations should ADD the specified coordinates to
    # the structure
    cdef bool subtract = False
    # A PyMatrix object
    cdef PyMatrix matrix
    # go through the movements backwards
    for I in range(J, 0, -1):
        # Get the PyMatrix
        matrix = movements[I-1]
        # do the positioning movement
        do_positioning_movement (residue, matrix, subtract)

cdef PyMatrix first_motion (PyResidue residue, PyAtom atom):
    """The first in a standard series of atom motions"""
    # Move the atoms so that the first atom is at the origin
    cdef PyMatrix movement = PyMatrix()
    movement.PyAtom_allocation(atom)
    do_positioning_movement(residue, movement)
    # Return the PyMatrix
    return movement

cdef PyMatrix second_motion (PyResidue residue, PyAtom atom):
    """The second in a standard series of atom motions"""
    # Rotate around the Z-axis so that the atom is in the XZ plane with
    # a positive X coordinate. 
    cdef float angle = math.fabs(math.acos(atom[0] /
               math.sqrt(math.pow(atom[0],2) + math.pow(atom[1],2))))
    if atom[1] > 0:
        angle = -angle
    cdef PyMatrix rotation = PyMatrix ()
    rotation.specified_rotation(angle, [0.0, 0.0, 1.0])
    do_positioning_movement(residue, rotation)
    # The matrix that gets returned is one that does the opposite rotation
    rotation.specified_rotation(-angle, [0.0, 0.0, 1.0])
    return rotation

cdef PyMatrix third_motion (PyResidue residue, PyAtom atom):
    """The third in a standard series of atom motions"""
    # Rotate around the Y-axis so that the atom is on the positive Z axis
    cdef float angle = math.fabs(math.acos(atom[2] /
               math.sqrt(math.pow(atom[2],2) + math.pow(atom[0],2))))
    if atom[0] > 0:
        angle = -angle
    cdef PyMatrix rotation = PyMatrix ()
    rotation.specified_rotation (angle, [0.0, 1.0, 0.0])
    do_positioning_movement(residue, rotation)
    rotation.specified_rotation (-angle, [0.0, 1.0, 0.0])
    return rotation

cdef PyMatrix fourth_motion (PyResidue residue, PyAtom atom):
    """The fourth in a series of standard atom motions"""
    return second_motion (residue, atom)

cdef list standard_positioning (PyResidue residue, list atoms):
    """Do the standard positioning of a list of atoms"""
    # This carries out the 4 movements in the standard set
    cdef list results = []
    results.append(first_motion(residue, atoms[0]))
    results.append(second_motion(residue, atoms[1]))
    results.append(third_motion(residue, atoms[1]))
    results.append(fourth_motion(residue, atoms[2]))
    return results

# Many of the interactions involve calculating a plane that goes through 3
# Atoms. These two functions do that.
cdef vector[float] calculate_normal_vector (PyAtom atom1, PyAtom atom2, 
                                            PyAtom atom3):
    """Using three Atoms, calculate a vector normal to their Plane"""
    # Looping variables
    cdef size_t I, J
    J = 3
    # The vector between atom2 and atom1
    cdef vector[float] v1
    v1.resize(J)
    for I in range(J):
        v1[I] = atom2._thisptr[0][I] - atom1._thisptr[0][I]
    # The vector between atom3 and atom1
    cdef vector[float] v2
    v2.resize(J)
    for I in range(J):
        v2[I] = atom3._thisptr[0][I] - atom1._thisptr[0][I]
    # The normal vector from the cross products of vectors 1 and 2
    cdef vector[float] answer
    answer.resize(3)
    answer[0] = v1[1]*v2[2] - v2[1]*v1[2]
    answer[1] = v1[2]*v2[0] - v2[2]*v1[0]
    answer[2] = v1[0]*v2[1] - v2[0]*v1[1]
    return answer

cdef vector[float] calculate_plane (PyAtom atom1, PyAtom atom2, PyAtom atom3):
    """Calculate the plane that goes through the 3 Atoms"""
    # Calculate a normal vector using the 3 atoms
    cdef vector[float] normal = calculate_normal_vector(atom1, atom2, atom3)
    # The equation of the plane is normal[0]*x + normal[1]*y + normal[2]*z =
    # d. This loop calculates D
    cdef float d = 0.0
    cdef size_t I, J
    J = normal.size()
    for I in range(J):
        d += normal[I]*atom1._thisptr[0][I]
    # The information for the plane will be stored using -d (so the equation
    # adds up to 0)
    normal.reserve(J+1)
    normal.push_back(-d)
    return normal

# Some functions that list the coordinates of an atom or the terms involved in
# a plane
cdef str list_coordinates (PyAtom atom):
    """List the coordinates of an Atom"""
    cdef size_t i1 = 0
    cdef size_t i2 = 1
    cdef size_t i3 = 2
    cdef str output = format(atom._thisptr[0][i1], '.3f') + " " \
                    + format(atom._thisptr[0][i2], '.3f') + " " \
                    + format(atom._thisptr[0][i3], '.3f')
    return output

cdef str list_plane (vector[float]& plane):
    """List the equation for a plane"""
    cdef str output = format(plane[0], '.3f') + " " + format(plane[1], '.3f')\
             + " " + format(plane[2], '.3f') + " " + format(plane[3], '.3f')
    return output

# Some of the more detailed rotations need to calculate an angle between 3
# atoms. To do that, I need to:
# 1) Make unit vectors between atoms
# 2) calculate dot product of the two unit vectors
# 3) calculate the angle

cdef vector[float] make_unit_vector(PyAtom atom1, PyAtom atom2):
    """Make a unit vector from the 2 atoms"""
    # Store the vector here
    cdef vector[float] answer
    # Looping variables
    cdef size_t I, J
    J = 3
    # Size the vector appropriately
    answer.resize(J)
    # Store the values in the vector
    for I in range(J):
        answer[I] = atom2._thisptr[0][I] - atom1._thisptr[0][I]
    # Calculate the magnitude of the vector
    cdef float mag = 0.0
    for I in range(J):
        mag += answer[I]**2
    mag = mag**0.5
    # Convert it into a unit vector
    for I in range(J):
        answer[I] /= mag
    return answer

cdef float calculate_angle (PyAtom atom1, PyAtom atom2, PyAtom atom3):
    """Calculate the angle between 3 atoms"""
    # Make 2 unit vectors
    cdef vector[float] v1 = make_unit_vector(atom1, atom2)
    cdef vector[float] v2 = make_unit_vector(atom1, atom3)
    # Calculate the dot product of those vectors
    cdef float dot = 0.0
    cdef size_t I, J
    J = 3
    for I in range(J):
        dot += (v1[I] * v2[I])
    # Calculate the angle. Based on how it will be used, I always want the
    # positive value
    cdef float angle = math.fabs(math.acos(dot))

# The definitions used in the rest of this file are non-standard, but
# necessary. Specifically, residues in the paratope are referred to as
# "takers" during interactions and residues in the epitope are referred to as
# "givers". I'm not using the words paratope or epitope in the functions
# because they could be changed or reversed in the future (this is already the
# second time I'm implementing similar concepts and they're reversing here).
# What matters is that when comparing different interactions, the "taker"
# residue is always in the same protein (i.e. the binding protein) and the
# "giver" residue is always in the same (other) protein (i.e. the antigen)

# Interaction 1 is a N-H in the "taker" protein having a hydrogen bond with an
# O=C in the giving protein
cdef list Interaction_1_Taker (PyResidue residue, list atoms):
    """Position the "taker" residue in Interaction 1"""
    # 2 Atoms are needed - the N and the H. Do the first 3 motions to position
    # them
    cdef list movements = []
    movements.append(first_motion(residue, atoms[0]))
    movements.append(second_motion(residue, atoms[1]))
    movements.append(third_motion(residue, atoms[1]))
    return movements

cdef void Interaction_1_Giver (PyResidue residue, list atoms):
    """Position the "giver" residue in Interaction 1"""
    # 2 Atoms are needed - the O and the C. Start by doing the first three
    # movements on them
    cdef PyMatrix movement = first_motion(residue, atoms[0])
    movement = second_motion(residue, atoms[1])
    movement = third_motion(residue, atoms[1])
    # Rotate the atoms 60 degrees counterclockwise around the Y-axis. This
    # moves the C into the positive XZ plane and will position the O so that
    # the proper angle is formed from the N - O - C. Note that in this case
    # the target angle is 120 degrees, so a 60 rotation is correct.
    movement.specified_rotation(math.radians(60.0), [0.0, 1.0, 0.0])
    do_positioning_movement(residue, movement)
    # Raise the Atoms so that the O will be 2.89 angstroms from the N
    movement.PyAtoms_allocation_for_centering([[0.0, 0.0, 2.89]])
    do_positioning_movement(residue, movement, False)

# Many Interactions involve having specific coordinates for one atom and a
# plane for another. This function does the calculations associated with
# those.
cdef void Spot_and_Plane (list movements, PyAtom atom, str spotName, 
        str planeName, PyResidue rotamer1, str LoopName, size_t StructNum, 
        size_t PosIndex, size_t RotNum, output):
    """Do Interaction steps involving a spot and a plane"""
    # Create a PyMatrix for rotating the rotamer 120 counterclockwise around
    # the Z axis
    cdef PyMatrix rotate = PyMatrix ()
    rotate.specified_rotation (math.radians(120.0), [0.0, 0.0, 1.0])
    # Create two more copies of the rotamer
    cdef PyResidue rotamer2 = rotamer1.duplicate()
    rotamer2.rotate(rotate)
    cdef PyResidue rotamer3 = rotamer2.duplicate()
    rotamer3.rotate(rotate)
    # Position the rotamers around the paratope
    reverse_positioning_movements (rotamer1, movements)
    reverse_positioning_movements (rotamer2, movements)
    reverse_positioning_movements (rotamer3, movements)
    # Get the atoms used for calculating the plane
    cdef PyAtom atom1 = get_atoms(rotamer1, [planeName])[0]
    cdef PyAtom atom2 = get_atoms (rotamer2, [planeName])[0]
    cdef PyAtom atom3 = get_atoms (rotamer3, [planeName])[0]
    # Calculate the plane parameters
    cdef vector[float] plane = calculate_plane(atom1, atom2, atom3)
    # Get the atom that goes in a specific spot
    cdef PyAtom spot = get_atoms (rotamer1, [spotName])[0]
    # Create the string of text describing all of this information
    cdef str text = LoopName + " " + str(StructNum) + " " + str(PosIndex) \
                  + " " + str(RotNum) + " " + list_coordinates(atom) + " " \
                  + list_coordinates(spot) + " " + list_plane(plane) + "\n"
    # Write that text to the output file
    output.write(text)

cdef void Interaction_1 (PyResidue residue, list atoms, str LoopName,
          size_t StructNum, size_t PosIndex, size_t RotNum, dict Rotamers,
          list outputs):
    """Generate and output information for Interaction 1"""
    # Do the positioning of the residue
    cdef list movements = Interaction_1_Taker (residue, atoms)
    # Put it back where it belongs
    reverse_positioning_movements (residue, movements)
    # The rotamer of GLN will be used as the representative structure. 
    cdef PyResidue rotamer1 = Rotamers["GLN"].duplicate()
    # Get the atoms from the rotamer
    cdef list rotAtoms = get_atoms(rotamer1, ['OE1', 'CD'])
    # Position it
    Interaction_1_Giver (rotamer1, rotAtoms)
    # Use the Spot_and_Plane function to do the rest of the calculations
    Spot_and_Plane (movements, atoms[1], 'OE1', 'CD', rotamer1, LoopName, 
                    StructNum, PosIndex, RotNum, outputs[0])

# Interaction 2 is a N-H in the "taker" protein having a hydrogen bond with a
# HO-C in the giving protein
cdef list Interaction_2_Taker (PyResidue residue, list atoms):
    """Position the "taker" residue in Interaction 2"""
    # This is the same as interaction 1
    return Interaction_1_Taker (residue, atoms)

cdef void Interaction_2_Giver (PyResidue residue, list atoms):
    """Position the "giver" residue in Interaction 2"""
    # Three atoms are needed here - O, C, and H. Do the standard positioning
    # on them
    cdef list movements = standard_positioning (residue, atoms)
    # Make a PyMatrix for additional movements
    cdef PyMatrix movement = PyMatrix ()
    # Calculate the angle between the atoms
    cdef float angle = calculate_angle (atoms[0], atoms[1], atoms[2])
    # Rotate the Atoms angle/2 degrees clockwise around the Y-axis to bring 
    # the C and H evenly above the XY plane
    movement.specified_rotation(-angle/2.0, [0.0, 1.0, 0.0])
    do_positioning_movement(residue, movement)
    # Tilt the residue so that a lone-pair on the O is pointing down the Z
    # axis
    movement.specified_rotation (math.radians(-54.25), [1.0, 0.0, 0.0])
    do_positioning_movement(residue, movement)
    # Raise the residue so that the O is 3.00 angstroms above the N
    movement.PyAtoms_allocation_for_centering ([[0.0, 0.0, 3.00]])
    do_positioning_movement(residue, movement, False)

cdef void Interaction_2 (PyResidue residue, list atoms, str LoopName,
          size_t StructNum, size_t PosIndex, size_t RotNum, dict Rotamers,
          list outputs):
    """Generate and output information for Interaction 2"""
    # Position the residue
    cdef list movements = Interaction_2_Taker (residue, atoms)
    # Reverse the position
    reverse_positioning_movements (residue, movements)
    # The SER rotamer is the representative rotamer in this case
    cdef PyResidue rotamer1 = Rotamers["SER"].duplicate()
    # Get the atoms
    cdef list rotAtoms = get_atoms (rotamer1, ['OG', 'CB', 'HG1'])
    # Position the rotamer
    Interaction_2_Giver (rotamer1, rotAtoms)
    # Finish the calculations
    Spot_and_Plane (movements, atoms[1], 'OG', 'CB', rotamer1, LoopName, 
                    StructNum, PosIndex, RotNum, outputs[1])

# Interaction 3 should be a double bonded O in the taker protein having an
# H-bond with an H-N in the giver protein (and Interaction 4 being with an
# H-O)
# However, these interactions are currently the only ones that have different
# properties than all the others. In particular, there is a plane that the H
# should be positioned in and a separate plane for the heavy atom. These are
# the only interactions that do not have specific points involved in moving
# atoms. As a result, they are different to handle than any other type and I
# am skipping them for now. Assuming the algorithm works well (and it should!)
# we should be able to add them in at a later date.

# Interaction 5 is a C-O-H in the taker protein forming an H-bond with a O=C
# in the giver protein.
cdef list Interaction_5_Taker (PyResidue residue, list atoms):
    """Position the "taker" residue in Interaction 5"""
    return Interaction_1_Taker (residue, atoms)

cdef void Interaction_5_Giver (PyResidue residue, list atoms):
    """Position the "giver" residue in Interaction 5"""
    # This is the same as interaction 1, except the motion is 2.79 angstroms
    cdef PyMatrix movement = first_motion (residue, atoms[0])
    movement = second_motion(residue, atoms[1])
    movement = third_motion (residue, atoms[1])
    movement.specified_rotation(math.radians(60.0), [0.0, 1.0, 0.0])
    do_positioning_movement(residue, movement)
    movement.PyAtoms_allocation_for_centering([[0.0, 0.0, 2.79]])
    do_positioning_movement(residue, movement, False)

cdef void Interaction_5 (PyResidue residue, list atoms, str LoopName, 
          size_t StructNum, size_t PosIndex, size_t RotNum, dict Rotamers,
          list outputs):
    """Generate and output information for Interaction 5"""
    # Repeat the steps from Interaction 1
    cdef list movements = Interaction_5_Taker (residue, atoms)
    reverse_positioning_movements (residue, movements)
    cdef PyResidue rotamer1 = Rotamers["GLN"].duplicate()
    cdef list rotAtoms = get_atoms (rotamer1, ['OE1', 'CD'])
    Interaction_5_Giver (rotamer1, rotAtoms)
    Spot_and_Plane (movements, atoms[1], 'OE1', 'CD', rotamer1, LoopName, 
                    StructNum, PosIndex, RotNum, outputs[4])

# Interaction 6 is a C-O-H in the taker protein forming an H-bond with a HO-C
# in the giver protein
cdef list Interaction_6_Taker (PyResidue residue, list atoms):
    """Position the "taker" residue in Interaction 6"""
    return Interaction_1_Taker (residue, atoms)

cdef void Interaction_6_Giver (PyResidue residue, list atoms):
    """Position the "giver" residue in Interaction 6"""
    # This is the same as interaction 2, except the distance is 2.86 angstroms
    cdef list movements = standard_positioning (residue, atoms)
    cdef PyMatrix movement = PyMatrix ()
    cdef float angle = calculate_angle(atoms[0], atoms[1], atoms[2])
    movement.specified_rotation(-angle/2.0, [0.0, 1.0, 0.0])
    do_positioning_movement (residue, movement)
    movement.specified_rotation(math.radians(-54.25), [1.0, 0.0, 0.0])
    movement.PyAtoms_allocation_for_centering ([[0.0, 0.0, 2.86]])
    do_positioning_movement (residue, movement, False)

cdef void Interaction_6 (PyResidue residue, list atoms, str LoopName,
          size_t StructNum, size_t PosIndex, size_t RotNum, dict Rotamers,
          list outputs):
    """Generate and output information for Interaction 6"""
    # Repeat the steps from Interaction 2
    cdef list movements = Interaction_6_Taker (residue, atoms)
    reverse_positioning_movements (residue, movements)
    cdef PyResidue rotamer1 = Rotamers["SER"].duplicate()
    cdef list rotAtoms = get_atoms(rotamer1, ['OG', 'CB', 'HG1'])
    Interaction_6_Giver (rotamer1, rotAtoms)
    Spot_and_Plane (movements, atoms[1], 'OG', 'CB', rotamer1, LoopName, 
                    StructNum, PosIndex, RotNum, outputs[5])

# Interaction 7 is a C-O-H in the taker protein acting as an H-bond acceptor
# with an N-H in the giving protein. In this case, there are two relevant
# orientations for the taker residue. That wasn't the situation in Interaction
# 2 - work through the positions of Atoms in each interaction to figure out
# why!
cdef list Interaction_7_Taker (PyResidue residue, list atoms, int which):
    """Position the "taker" residue in Interaction 7"""
    # The three atoms are O, C, H. Do the standard positioning of them
    cdef list movements = standard_positioning (residue, atoms)
    # The residue needs to be rotated counterclockwise around the Y-axis
    # so that the C and H are evenly below the X-axis
    cdef float angle = calculate_angle (atoms[0], atoms[1], atoms[2])
    angle = (2.0*math.pi - angle)/2.0
    cdef PyMatrix rotate1 = PyMatrix()
    rotate1.specified_rotation (angle, [0.0, 1.0, 0.0])
    do_positioning_movement(residue, rotate1)
    # Store the inverse of that movement
    rotate1.specified_rotation(-angle, [0.0, 1.0, 0.0])
    movements.append(rotate1)
    # Depending on the which integer, rotate 54.25 degrees clockwise or
    # counterclockwise around the X-axis to point a lone-pair of electrons on
    # the O up the Z-axis
    cdef PyMatrix rotate2 = PyMatrix()
    if which == 1:
        rotate2.specified_rotation (math.radians(54.25), [1.0, 0.0, 0.0])
    else:
        rotate2.specified_rotation(math.radians(-54.25), [1.0, 0.0, 0.0])
    do_positioning_movement(residue, rotate2)
    # Store the inverse of that rotation
    if which == 1:
        rotate2.specified_rotation(math.radians(-54.25), [1.0, 0.0, 0.0])
    else:
        rotate2.specified_rotation (math.radians(54.25), [1.0, 0.0, 0.0])
    movements.append(rotate2)
    # Return the movements list
    return movements

cdef void Interaction_7_Giver (PyResidue residue, list atoms):
    """Position the "giver" residue for Interaction 7"""
    # Do the first three movements
    cdef PyMatrix movement = first_motion (residue, atoms[0])
    movement = second_motion(residue, atoms[1])
    movement = third_motion (residue, atoms[1])
    # Rotate 180 counterclockwise around the Y-axis to point the H down the
    # Z-axis
    movement.specified_rotation(math.radians(180.0), [0.0, 1.0, 0.0])
    do_positioning_movement(residue, movement)
    # Raise the atoms so the N is 3.00 angstroms above the O
    movement.PyAtoms_allocation_for_centering([[0.0, 0.0, 3.00]])
    do_positioning_movement(residue, movement, False)

# Another common set of Interaction information are a pair of spots. This
# function does those calculations
cdef void Spot_and_Spot (list movements, str name1, str name2, 
        PyResidue rotamer1, str LoopName, size_t StructNum, size_t PosIndex, 
        size_t RotNum, output):
    """Do Interaction calculations involving 2 spots"""
    # Reverse the positioning of the rotamer
    reverse_positioning_movements (rotamer1, movements)
    # Get the two Atoms
    cdef PyAtom atom1 = get_atoms(rotamer1, [name1])[0]
    cdef PyAtom atom2 = get_atoms(rotamer1, [name2])[0]
    # Create the string of text listing the information
    cdef str text = LoopName + " " + str(StructNum) + " " + str(PosIndex) \
                  + " " + str(RotNum) + " " + list_coordinates(atom1) \
                  + " " + list_coordinates(atom2) + "\n"
    # Write the information to the file
    output.write(text)

cdef void Interaction_7 (PyResidue residue, list atoms, str LoopName,
          size_t StructNum, size_t PosIndex, size_t RotNum, dict Rotamers,
          list outputs):
    """Generate and output the data for Interaction 7"""
    # This is repeated twice, to generate the information for each lone pair
    # on the Oxygen
    cdef int i
    for i in range(2):
        movements = Interaction_7_Taker (residue, atoms, i)
        reverse_positioning_movements (residue, movements)
        # The GLN residue is used as the representative structure
        rotamer = Rotamers["GLN"].duplicate()
        # Get the atoms
        rotAtoms = get_atoms(rotamer, ['NE2', 'HE21'])
        # Position the Rotamer
        Interaction_7_Giver (rotamer, rotAtoms)
        # Do the spot and spot calculations
        Spot_and_Spot (movements, 'HE21', 'NE2', rotamer, LoopName, StructNum,
                       PosIndex, RotNum, outputs[6])

# Interaction 8 is a C-O-H in the taker protein being a hydrogen bond acceptor
# with a H-O in the giving protein
cdef list Interaction_8_Taker (PyResidue residue, list atoms, int which):
    """Position the "taker" residue in Interaction 8"""
    return Interaction_7_Taker (residue, atoms, which)

cdef void Interaction_8_Giver (PyResidue residue, list atoms):
    """Position the "giver" residue in Interaction 8"""
    # Do the same things as in Interaction 7, but the distance should be 2.86
    # angstroms
    cdef PyMatrix movement = first_motion(residue, atoms[0])
    movement = second_motion(residue, atoms[1])
    movement = third_motion(residue, atoms[1])
    movement.specified_rotation(math.radians(180.0), [0.0, 1.0, 0.0])
    do_positioning_movement(residue, movement)
    movement.PyAtoms_allocation_for_centering([[0.0, 0.0, 2.86]])
    do_positioning_movement(residue, movement, False)

cdef void Interaction_8 (PyResidue residue, list atoms, str LoopName, 
          size_t StructNum, size_t PosIndex, size_t RotNum, dict Rotamers,
          list outputs):
    """Generate and output the data for Interaction 8"""
    # This is basically the same as Interaction 7, just using the SER atoms
    # instead
    cdef int i
    for i in range(2):
        movements = Interaction_8_Taker (residue, atoms, i)
        reverse_positioning_movements (residue, movements)
        rotamer = Rotamers["SER"].duplicate()
        rotAtoms = get_atoms(rotamer, ['OG', 'HG1'])
        Interaction_8_Giver (rotamer, rotAtoms)
        Spot_and_Spot (movements, "HG1", "OG", rotamer, LoopName, StructNum,
                       PosIndex, RotNum, outputs[7])

# Interaction 9 is an ARG in the taker protein forming double hydrogen bonds
# with the carboxylic acid of GLU or ASP
cdef list Interaction_9_Taker (PyResidue residue, list atoms):
    """Position the "taker" residue in Interaction 9"""
    # There are 3 atoms - C and 2 Ns. Do the standard positioning of them
    cdef list movements = standard_positioning (residue, atoms)
    # Rotate clockwise around the Y-axis to bring both N's evenly over the
    # x-axis. Start by calculating the angle between the atoms
    cdef float angle = calculate_angle (atoms[0], atoms[1], atoms[2])
    cdef PyMatrix rotate = PyMatrix ()
    rotate.specified_rotation(-angle/2.0, [0.0, 1.0, 0.0])
    do_positioning_movement(residue, rotate)
    # Store the inverse of that rotation
    rotate.specified_rotation(math.radians(angle/2.0), [0.0, 1.0, 0.0])
    movements.append(rotate)
    return movements

cdef void Interaction_9_Giver (PyResidue residue, list atoms):
    """Position the "giver" residue in Interaction 9"""
    # There should be three atoms - the C and 2 Os
    cdef list movements = standard_positioning (residue, atoms)
    # Rotate counterclockwise around the Y-axis to bring both atoms evenly
    # below the x-axis. To do so, start by calculating the angle between the
    # Atoms
    cdef float angle = calculate_angle (atoms[0], atoms[1], atoms[2])
    # Change this to the "remaining" (i.e. outside) angle of the system
    # divided by 2.
    angle = (2.0*math.pi - angle) / 2.0
    # Rotate by that amount to position the atoms
    cdef PyMatrix movement = PyMatrix ()
    movement.specified_rotation(angle, [0.0, 1.0, 0.0])
    do_positioning_movement (residue, movement)
    # Raise the residue so that the C is 3.96 angstroms from the C of the ARG
    movement.PyAtoms_allocation_for_centering([[0.0, 0.0, 3.96]])
    do_positioning_movement (residue, movement, False)

cdef void Interaction_9 (PyResidue residue, list atoms, str LoopName,
          size_t StructNum, size_t PosIndex, size_t RotNum, dict Rotamers,
          list outputs):
    """Generate and output the data for Interaction 9"""
    # Generate the motions to position the residue
    cdef list movements = Interaction_9_Taker (residue, atoms)
    # Put the residue back into the taker protein
    reverse_positioning_movements (residue, movements)
    # Get a GLU rotamer
    cdef PyResidue rotamer = Rotamers["GLU"].duplicate()
    # Get the needed atoms for positioning it
    cdef list rotAtoms = get_atoms(rotamer, ['CD', 'OE1', 'OE2'])
    # Position the rotamer
    Interaction_9_Giver (rotamer, rotAtoms)
    # Reverse it's positioning information
    reverse_positioning_movements (rotamer, movements)
    # Output the coordinate information of the three Atoms
    cdef str text = LoopName + " " + str(StructNum) + " " + str(PosIndex) \
                  + " " + str(RotNum) + " " + list_coordinates(rotAtoms[0]) \
                  + " " + list_coordinates(rotAtoms[1]) + " " \
                  + list_coordinates(rotAtoms[2]) + "\n"
    outputs[8].write(text)

# Interaction 10 is a LYS in the taker protein forming a strong hydrogen bond
# complex with the carboxylic acid of GLU or ASP
cdef list Interaction_10_Taker (PyResidue residue, list atoms):
    """Position the "taker" residue in Interaction 10"""
    # There are 2 atoms - the N and C. Position them so that the N is at the
    # origin and the C is on the NEGATIVE Z-axis
    cdef list movements = []
    movements.append(first_motion(residue, atoms[0]))
    movements.append(second_motion(residue, atoms[1]))
    movements.append(third_motion(residue, atoms[1]))
    cdef PyMatrix rotate = PyMatrix ()
    rotate.specified_rotation(math.radians(180.0), [0.0, 1.0, 0.0])
    do_positioning_movement (residue, rotate)
    # This is the one case where it is OK to not make an inverted version of
    # the rotation - but do so anyway
    rotate.specified_rotation(math.radians(-180.0), [0.0, 1.0, 0.0])
    movements.append(rotate)
    return movements

cdef void Interaction_10_Giver (PyResidue residue, list atoms):
    """Position the "giver" residue in Interaction 10"""
    # Do the same things as in Interaction 9, but now the distance is 3.18
    # angstroms
    cdef list movements = standard_positioning (residue, atoms)
    cdef PyMatrix movement = PyMatrix ()
    cdef float angle = calculate_angle (atoms[0], atoms[1], atoms[2])
    angle = (2.0*math.pi - angle)/2.0
    movement.specified_rotation(angle, [0.0, 1.0, 0.0])
    do_positioning_movement (residue, movement)
    movement.PyAtoms_allocation_for_centering ([[0.0, 0.0, 3.18]])
    do_positioning_movement (residue, movement, False)

cdef void Interaction_10 (PyResidue residue, list atoms, str LoopName,
          size_t StructNum, size_t PosIndex, size_t RotNum, dict Rotamers,
          list outputs):
    """Generate and output the data for Interaction 10"""
    # Do the typical set of preliminary steps
    cdef list movements = Interaction_10_Taker (residue, atoms)
    reverse_positioning_movements (residue, movements)
    cdef PyResidue rotamer = Rotamers["GLU"].duplicate()
    cdef list rotAtoms = get_atoms (rotamer, ["CD", "OE1", "OE2"])
    Interaction_10_Giver (rotamer, rotAtoms)
    # This is a case where a point (for the carbon) and a plane (for the
    # oxygens) is appropriate
    Spot_and_Plane (movements, atoms[0], 'CD', 'OE1', rotamer, LoopName, 
                    StructNum, PosIndex, RotNum, outputs[9])

# Interaction 11 is a GLU/ASP in the taker protein forming a double hydrogen
# bond with the guanidine (?) in a giver protein ARG
cdef list Interaction_11_Taker (PyResidue residue, list atoms):
    """Position the "taker" residue in Interaction 11"""
    return Interaction_9_Taker (residue, atoms)

cdef void Interaction_11_Giver (PyResidue residue, list atoms):
    """Position the "giver" residue In Interaction 11"""
    Interaction_9_Giver (residue, atoms)

cdef void Interaction_11 (PyResidue residue, list atoms, str LoopName,
          size_t StructNum, size_t PosIndex, size_t RotNum, dict Rotamers,
          list outputs):
    """Generate and output the data for Interaction 11"""
    # Do the standard initial steps
    cdef list movements = Interaction_11_Taker (residue, atoms)
    reverse_positioning_movements (residue, movements)
    # Get an ARG residue
    cdef PyResidue rotamer = Rotamers["ARG"].duplicate()
    # Get it's Atoms
    cdef list rotAtoms = get_atoms (rotamer, ['CZ', 'NH1', 'NH2'])
    # Position it and output the atom coordinates
    Interaction_11_Giver (rotamer, atoms)
    reverse_positioning_movements (rotamer, movements)
    cdef str text = LoopName + " " + str(StructNum) + " " + str(PosIndex) \
                  + " " + str(RotNum) + " " + list_coordinates(rotAtoms[0]) \
                  + " " + list_coordinates(rotAtoms[1]) + " " \
                  + list_coordinates(rotAtoms[2]) + "\n"
    outputs[10].write(text)

# Interaction 12 is a GLU/ASP in the taker protein forming a strong hydrogen
# bond complex with a LYS in the giver protein
cdef list Interaction_12_Taker (PyResidue residue, list atoms):
    """Position the "taker" residue in Interaction 12"""
    return Interaction_9_Taker(residue, atoms)

cdef void Interaction_12_Giver (PyResidue residue, list atoms):
    """Position the "giver" residue in Interaction 12"""
    # Position the N at the origin and the C on the positive Z axis, then
    # raise them by 3.18 angstroms
    cdef PyMatrix movement = first_motion(residue, atoms[0])
    movement = second_motion (residue, atoms[1])
    movement = third_motion (residue, atoms[1])
    movement.PyAtoms_allocation_for_centering([[0.0, 0.0, 3.18]])
    do_positioning_movement(residue, movement, False)

cdef void Interaction_12 (PyResidue residue, list atoms, str LoopName, 
          size_t StructNum, size_t PosIndex, size_t RotNum, dict Rotamers,
          list outputs):
    """Generate and output the data for Interaction 12"""
    # Position the residue
    cdef list movements = Interaction_12_Taker (residue, atoms)
    reverse_positioning_movements (residue, movements)
    # Get a LYS rotamer
    cdef PyResidue rotamer = Rotamers["LYS"].duplicate()
    # Get the appropriate atoms
    cdef list rotAtoms = get_atoms (rotamer, ['NZ', 'CE'])
    # Position the rotamer
    Interaction_12_Giver (rotamer, rotAtoms)
    reverse_positioning_movements (rotamer, movements)
    # This is a case where 2 atom positions are needed for the interaction
    # calculation
    Spot_and_Spot (movements, 'NZ', 'CE', rotamer, LoopName, StructNum,
                   PosIndex, RotNum, outputs[11])

# Interaction 13 is an ARG in the taker protein positioned above an aromatic
# ring in the giver protein
cdef list Interaction_13_Taker (PyResidue residue, list atoms):
    """Position the "taker" residue in Interaction 13"""
    # Do the standard positioning of the provided atoms (C, N, N)
    return standard_positioning (residue, atoms)

cdef void Interaction_13_Giver (PyResidue residue, list atoms):
    """Position the "giver" residue in Interaction 13"""
    # Do the standard positioning on the first 3 atoms. This puts the aromatic
    # ring in the XZ plane
    cdef list movements = standard_positioning(residue, atoms)
    # There should be 6 atoms in this list of atoms. Move them so that their
    # center of mass is at the origin
    cdef PyMatrix movement = PyMatrix ()
    movement.PyAtoms_allocation_for_centering(atoms)
    do_positioning_movement (residue, movement)
    # Move the residue 3.75 angstroms along the positive Y axis
    movement.PyAtoms_allocation_for_centering([[0.0, 3.75, 0.0]])
    do_positioning_movement (residue, movement, False)

# Some aromatic interactions need to know where the center of a ring is and
# the plane the ring atoms are located in
cdef void Center_and_Plane (list atoms, str LoopName, size_t StructNum, 
          size_t PosIndex, size_t RotNum, output):
    """Output the center of a ring and the plane it is located in"""
    # Calculate the center of the ring
    cdef PyMatrix center = PyMatrix ()
    center.PyAtoms_allocation_for_centering (atoms)
    # Calculate the plane of the ring
    cdef vector[float] plane = calculate_plane (atoms[0], atoms[1], atoms[2])
    # Create a string of the information
    cdef str text = LoopName + " " + str(StructNum) + " " + str(PosIndex) \
                  + " " + str(RotNum) \
                  + " " + format(center(0, 0), '.3f') \
                  + " " + format(center(0, 1), '.3f') \
                  + " " + format(center(0, 2), '.3f') \
                  + " " + list_plane(plane) + "\n"
    # Write that to the file
    output.write(text)

cdef void Interaction_13 (PyResidue residue, list atoms, str LoopName,
          size_t StructNum, size_t PosIndex, size_t RotNum, dict Rotamers,
          list outputs):
    """Generate and output Interaction 13 data"""
    cdef list movements = Interaction_13_Taker (residue, atoms)
    reverse_positioning_movements (residue, movements)
    # Get a PHE rotamer
    cdef PyResidue rotamer = Rotamers["PHE"].duplicate()
    # Get the atoms of the PHE
    cdef list rotAtoms = get_atoms(rotamer, ['CZ', 'CE1', 'CD1', 'CG', 'CD2', 'CE2'])
    # Position the rotamer
    Interaction_13_Giver (rotamer, rotAtoms)
    reverse_positioning_movements (rotamer, movements)
    # Output the information
    Center_and_Plane (rotAtoms, LoopName, StructNum, PosIndex, RotNum,
                      outputs[12])

# Interaction 14 is a LYS in the taker protein positioned above an aromatic
# ring in the giver protein
cdef list Interaction_14_Taker (PyResidue residue, list atoms):
    """Position the "taker" residue in Interaction 14"""
    # There are 2 atoms. The N and the C. Position the N at the origin and the
    # C on the negative Y axis
    cdef list movements = []
    movements.append(first_motion(residue, atoms[0]))
    movements.append(second_motion(residue, atoms[1]))
    movements.append(third_motion(residue, atoms[1]))
    cdef PyMatrix rotate = PyMatrix ()
    rotate.specified_rotation(math.radians(90.0), [1.0, 0.0, 0.0])
    do_positioning_movement(residue, rotate)
    rotate.specified_rotation (math.radians(-90.0), [1.0, 0.0, 0.0])
    movements.append(rotate)
    return movements

cdef void Interaction_14_Giver (PyResidue residue, list atoms):
    """Position the "giver" residue in Interaction 14"""
    # This is the same as Interaction 13, but the distance is 4.20 angstroms
    # away on the Y axis
    cdef list movements = standard_positioning(residue, atoms)
    cdef PyMatrix movement = PyMatrix()
    movement.PyAtoms_allocation_for_centering(atoms)
    do_positioning_movement (residue, movement)
    movement.PyAtoms_allocation_for_centering([[0.0, 4.20, 0.0]])
    do_positioning_movement (residue, movement, False)

cdef void Interaction_14 (PyResidue residue, list atoms, str LoopName, 
          size_t StructNum, size_t PosIndex, size_t RotNum, dict Rotamers,
          list outputs):
    """Generate and output the data for Interaction 14"""
    cdef list movements = Interaction_14_Taker (residue, atoms)
    reverse_positioning_movements (residue, movements)
    # Get a PHE rotamer and its aromatic atoms
    cdef PyResidue rotamer = Rotamers["PHE"].duplicate()
    cdef list rotAtoms = get_atoms (rotamer, ['CZ','CE1','CD1','CG','CD2','CE2'])
    # Position it
    Interaction_14_Giver (rotamer, rotAtoms)
    reverse_positioning_movements (rotamer, movements)
    # Output the information
    Center_and_Plane (rotAtoms, LoopName, StructNum, PosIndex, RotNum,
                      outputs[13])

# Interaction 15 is an aromatic residue in the taker protein interacting with
# an ARG in the giver protein
cdef list Interaction_15_Taker (PyResidue residue, list atoms):
    """Position the "taker" residue for Interaction 15"""
    # Do the standard positioning to get the aromatic ring in the XZ plane
    cdef list movements = standard_positioning(residue, atoms)
    # Move the ring so its center of mass is at the origin
    cdef PyMatrix movement = PyMatrix ()
    movement.PyAtoms_allocation_for_centering (atoms)
    do_positioning_movement(residue, movement)
    movements.append(movement)
    return movements

cdef void Interaction_15_Giver (PyResidue residue, list atoms):
    """Position the "giver" residue for Interaction 15"""
    # Do the standard positioning of the atoms so the guanidine (?) group of
    # the arginine is in the XZ plane
    cdef list movements = standard_positioning(residue, atoms)
    # Move it out by 3.75 angstroms along the positive Y axis
    cdef PyMatrix movement = PyMatrix ()
    movement.PyAtoms_allocation_for_centering([[0.0, 3.75, 0.0]])
    do_positioning_movement(residue, movement, False)

cdef void Interaction_15 (PyResidue residue, list atoms, str LoopName, 
          size_t StructNum, size_t PosIndex, size_t RotNum, dict Rotamers,
          list outputs):
    """Generate and output the data for Interaction 15"""
    # Position the residue
    cdef list movements = Interaction_15_Taker (residue, atoms)
    reverse_positioning_movements (residue, movements)
    # Get an ARG residue and its atoms
    cdef PyResidue rotamer = Rotamers["ARG"].duplicate()
    cdef list rotAtoms = get_atoms(rotamer, ['CZ', 'NH1', 'NH2'])
    # Position the ARG
    Interaction_15_Giver (rotamer, rotAtoms)
    reverse_positioning_movements (rotamer, movements)
    # Calculate the plane of the three atoms
    cdef vector[float] plane = calculate_plane (rotAtoms[0], rotAtoms[1],
                                                rotAtoms[2])
    # Create the summary text
    cdef str text = LoopName + " " + str(StructNum) + " " + str(PosIndex) \
                  + " " + str(RotNum) + " " + list_coordinates(rotAtoms[0]) \
                  + " " + list_plane(plane) + "\n"
    # Write it to the file
    outputs[14].write(text)

# Interaction 16 is an aromatic residue in the taker protein interaction with
# the cation of LYS in the giver protein
cdef list Interaction_16_Taker (PyResidue residue, list atoms):
    """Position the "taker" residue for Interaction 16"""
    return Interaction_15_Taker (residue, atoms)

cdef void Interaction_16_Giver (PyResidue residue, list atoms):
    """Position the "giver" residue for Interaction 16"""
    # Position the N at the origin, the C on the positive Y axis, then move
    # the atoms 4.20 angstroms along the Y-axis
    cdef PyMatrix movement = first_motion(residue, atoms[0])
    movement = second_motion(residue, atoms[1])
    movement = third_motion (residue, atoms[1])
    movement.specified_rotation (math.radians(-90.0), [1.0, 0.0, 0.0])
    do_positioning_movement (residue, movement)
    movement.PyAtoms_allocation_for_centering([[0.0, 4.20, 0.0]])
    do_positioning_movement (residue, movement, False)

cdef void Interaction_16 (PyResidue residue, list atoms, str LoopName,
          size_t StructNum, size_t PosIndex, size_t RotNum, dict Rotamers,
          list outputs):
    """Generate and output the data for Interaction 16"""
    cdef list movements = Interaction_16_Taker (residue, atoms)
    reverse_positioning_movements (residue, movements)
    # Get a LYS residue and its atoms
    cdef PyResidue rotamer = Rotamers["LYS"].duplicate()
    cdef list rotAtoms = get_atoms(rotamer, ['NZ', 'CE'])
    # Position that rotamer
    Interaction_16_Giver (rotamer, rotAtoms)
    # This is a case where 2 spots are appropriate
    Spot_and_Spot (movements, 'NZ', 'CE', rotamer, LoopName, StructNum,
                   PosIndex, RotNum, outputs[15])

# Interaction 17 is an aromatic residue in the taker protein in a parallel
# displaced pi-pi stacking interaction with an aromatic residue in the giver
# protein
cdef list Interaction_17_Taker (PyResidue residue, list atoms):
    """Position the "taker" residue for Interaction 17"""
    # Do the standard positioning of the atoms
    cdef list movements = standard_positioning (residue, atoms)
    # Move the atoms so that the center of the ring is at the origin. In this
    # case, there are actually 7 atoms, but the first 6 form the ring
    cdef PyMatrix movement = PyMatrix ()
    movement.PyAtoms_allocation_for_centering(atoms[:6])
    do_positioning_movement(residue, movement)
    movements.append(movement)
    # Rotate the 7th atom onto the positive Z axis (it should already be in
    # the plane by virtue of the ring being in the plane
    movements.append(third_motion (residue, atoms[6]))
    return movements

cdef void Interaction_17_Giver (PyResidue residue, list atoms):
    """Position the "giver" residue for Interaction 17"""
    # Do the standard positioning of the atoms
    cdef list movements = standard_positioning (residue, atoms)
    # Move the first 6 atoms so the center of mass of the ring is at the
    # origin
    cdef PyMatrix movement = PyMatrix ()
    movement.PyAtoms_allocation_for_centering (atoms[:6])
    do_positioning_movement(residue, movement)
    # Rotate the 7th atom onto the positive Z axis
    movement = third_motion(residue, atoms[6])
    # Rotate 180 counterclockwise around the Y-axis so that the 7th atom is on
    # the negative Z axis
    movement.specified_rotation(math.radians(180.0), [0.0, 1.0, 0.0])
    do_positioning_movement(residue, movement)
    # Move the ring 3.75 angstroms along the positve Y axis and 1.41 angstroms
    # along the negative Z axis
    movement.PyAtoms_allocation_for_centering([[0.0, 3.75, -1.41]])
    do_positioning_movement(residue, movement, False)

# Another interaction data option is a coordinate for a specific atom in an
# aromatic ring (a "tip" atom) and the equation for the plane the atoms are in
cdef void Tip_and_Plane (list atoms, str tipName, str LoopName, 
          size_t StructNum, size_t PosIndex, size_t RotNum, output):
    """List information for an atom position and the plane of a ring"""
    # Calculate the plane
    cdef vector[float] plane = calculate_plane(atoms[0], atoms[1], atoms[2])
    # Find the atom that is of interest
    atom = None
    for value in atoms:
        if value.name() == tipName:
            atom = value
            break
    # Create the text
    cdef str text = LoopName + " " + str(StructNum) + " " + str(PosIndex) \
                  + " " + str(RotNum) + " " + list_coordinates(atom) \
                  + " " + list_plane (plane) + "\n"
    # Write that information to the file
    output.write(text)

cdef void Interaction_17 (PyResidue residue, list atoms, str LoopName,
          size_t StructNum, size_t PosIndex, size_t RotNum, dict Rotamers,
          list outputs):
    """Generate and output data for Interaction 17"""
    # Position the residue
    cdef list movements = Interaction_17_Taker (residue, atoms)
    reverse_positioning_movements (residue, movements)
    # Get a PHE rotamer
    cdef PyResidue rotamer = Rotamers["PHE"].duplicate()
    # Get the atoms that are needed
    cdef list rotAtoms = get_atoms (rotamer, ['CZ', 'CE1', 'CD1', 'CG', 'CD2',
                                              'CE2', 'CB'])
    # Position the PHE
    Interaction_17_Giver (rotamer, rotAtoms)
    reverse_positioning_movements (rotamer, movements)
    # Use the most recent output function
    Tip_and_Plane (rotAtoms, 'CZ', LoopName, StructNum, PosIndex, RotNum,
                   outputs[16])

# Interaction 18 is an aromatic residue in the taker protein being the top of
# the "T" in a pi-pi interaction with an aromatic residue in the giver protein
cdef list Interaction_18_Taker (PyResidue residue, list atoms):
    """Position the "taker" residue for Interaction 18"""
    # Do the same motions as in Interaction 17
    return Interaction_17_Taker (residue, atoms)

cdef void Interaction_18_Giver (PyResidue residue, list atoms):
    """Position the "giver" residue for Interaction 18"""
    # Do the standard positioning of the atoms
    cdef list movements = standard_positioning (residue, atoms)
    # Move the first 6 atoms so the center of mass of the ring as at the
    # origin
    cdef PyMatrix movement = PyMatrix ()
    movement.PyAtoms_allocation_for_centering(atoms[:6])
    do_positioning_movement (residue, movement)
    # Rotate the 7th atom onto the positive Z axis
    movement = third_motion (residue, atoms[6])
    # Rotate the ring 90 clockwise around the X axis
    movement.specified_rotation (math.radians(-90.0), [1.0, 0.0, 0.0])
    do_positioning_movement (residue, movement)
    # Move the ring 5.0 angstroms along the positive Y axis
    movement.PyAtoms_allocation_for_centering ([[0.0, 5.0, 0.0]])
    do_positioning_movement(residue, movement, False)

# For this interaction, the 'tip' of the rotamer needs to be specified, as
# well as the center of mass of the ring it is in
cdef void Tip_and_Center (list atoms, str tipName, str LoopName, 
          size_t StructNum, size_t PosIndex, size_t RotNum, output):
    """Output information about an Atom and the center of a Ring's spots"""
    # Find the 'tip' atom
    atom = None
    for value in atoms:
        if value.name() == tipName:
            atom = value
            break
    # Calculate the center of the ring
    cdef PyMatrix matrix = PyMatrix ()
    matrix.PyAtoms_allocation_for_centering(atoms)
    # Include the information in a string of text
    cdef str text = LoopName + " " + str(StructNum) + " " + str(PosIndex) \
                  + " " + str(RotNum) + " " + list_coordinates(atom) \
                  + " " + format(matrix(0, 0), '.3f') \
                  + " " + format(matrix(0, 1), '.3f') \
                  + " " + format(matrix(0, 2), '.3f') + "\n"
    # Write the information to the output file
    output.write(text)

cdef void Interaction_18 (PyResidue residue, list atoms, str LoopName,
          size_t StructNum, size_t PosIndex, size_t RotNum, dict Rotamers,
          list outputs):
    """Generate and output data for Interaction 18"""
    cdef list movements = Interaction_18_Taker (residue, atoms)
    reverse_positioning_movements (residue, movements)
    cdef PyResidue rotamer = Rotamers["PHE"].duplicate()
    cdef list rotAtoms = get_atoms (rotamer, ['CZ', 'CE1', 'CD1', 'CG', 'CD2',
                                              'CE2', 'CB'])
    Interaction_18_Giver (rotamer, rotAtoms)
    reverse_positioning_movements (rotamer, movements)
    Tip_and_Center (rotAtoms, 'CZ', LoopName, StructNum, PosIndex, RotNum, 
                    outputs[17])

# Interaction 19 is an aromatic residue in the taker protein being the stem of
# the "T" in a pi-pi interaction with an aromatic residue in the giver protein
cdef list Interaction_19_Taker (PyResidue residue, list atoms):
    """Position the "taker" residue for Interaction 19"""
    # Do the same motion as in Interactions 17 and 18
    return Interaction_17_Taker (residue, atoms)

cdef void Interaction_19_Giver (PyResidue residue, list atoms):
    """Position the "giver" residue for Interaction 19"""
    # This will be the same as Interaction 18, except the ring will be moved
    # 5.0 angstroms down the Z axis instead of along the Y-axis
    cdef list movements = standard_positioning (residue, atoms)
    cdef PyMatrix movement = PyMatrix ()
    movement.PyAtoms_allocation_for_centering (atoms[:6])
    do_positioning_movement (residue, movement)
    movement = third_motion(residue, atoms[6])
    movement.specified_rotation (math.radians(-90.0), [1.0, 0.0, 0.0])
    do_positioning_movement (residue, movement)
    movement.PyAtoms_allocation_for_centering ([[0.0, 0.0, -5.0]])
    do_positioning_movement (residue, movement, False)

cdef void Interaction_19 (PyResidue residue, list atoms, str LoopName, 
          size_t StructNum, size_t PosIndex, size_t RotNum, dict Rotamers,
          list outputs):
    """Generate and output data for Interaction 19"""
    cdef list movements = Interaction_19_Taker (residue, atoms)
    reverse_positioning_movements (residue, movements)
    cdef PyResidue rotamer = Rotamers["PHE"].duplicate()
    cdef list rotAtoms = get_atoms (rotamer, ['CZ', 'CE1', 'CD1', 'CG', 'CD2',
                                              'CE2', 'CB'])
    Interaction_19_Giver(rotamer, rotAtoms)
    reverse_positioning_movements (rotamer, movements)
    # In this case, the Center and Plane information is appropriate
    Center_and_Plane (rotAtoms[:6], LoopName, StructNum, PosIndex, RotNum,
                      outputs[18])

# Histidines can act as cations just like ARG and LYS can. These next 4
# interactions account for HIS - pi interactions
# In interaction 20, a histidine in the taker protein is parallel to an
# aromatic ring in the giver protein
cdef list Interaction_20_Taker (PyResidue residue, list atoms):
    """Position the "taker" residue for Interaction 20"""
    # Position the HIS so it's ring is in the XZ plane, centered at the origin
    return Interaction_15_Taker (residue, atoms)

cdef void Interaction_20_Giver (PyResidue residue, list atoms):
    """Position the "giver" residue for Interaction 20"""
    # Do the standard positioning of the atoms
    cdef list movements = standard_positioning (residue, atoms)
    # Move the atoms in the ring so their center is at the origin
    cdef PyMatrix movement = PyMatrix ()
    movement.PyAtoms_allocation_for_centering (atoms)
    do_positioning_movement (residue, movement)
    # Move the Residue 3.5 angstroms out along the Y-axis
    movement.PyAtoms_allocation_for_centering([[0.0, 3.50, 0.0]])
    do_positioning_movement (residue, movement, False)

cdef void Interaction_20 (PyResidue residue, list atoms, str LoopName, 
          size_t StructNum, size_t PosIndex, size_t RotNum, dict Rotamers,
          list outputs):
    """Generate and output data for Interaction 20"""
    # Position the residue
    cdef list movements = Interaction_20_Taker (residue, atoms)
    reverse_positioning_movements (residue, movements)
    # Get a PHE rotamer and it's ring atoms
    cdef PyResidue rotamer = Rotamers["PHE"].duplicate()
    cdef list rotAtoms = get_atoms (rotamer, ['CZ', 'CE1', 'CD1', 'CG', 'CD2',
                                              'CE2'])
    # Position that rotamer
    Interaction_20_Giver (rotamer, rotAtoms)
    reverse_positioning_movements (rotamer, movements)
    # Store the location of the center of the aromatic ring and the plane the
    # ring is located in
    Center_and_Plane (rotAtoms, LoopName, StructNum, PosIndex, RotNum,
                      outputs[19])

# In Interaction 21, histidine is the stem of a "T" pointing at the center of
# an aromatic ring
cdef list Interaction_21_Taker (PyResidue residue, list atoms):
    """Position the "taker" residue for Interaction 21"""
    # Do the standard positioning of the three atoms, with the "pointing" atom
    # placed at the origin and the other 2 being the ones around it in the
    # ring
    cdef list movements = standard_positioning (residue, atoms)
    # Calculate the angle between the three atoms
    cdef float angle = calculate_angle (atoms[0], atoms[1], atoms[2])
    # Calculate half of the remaining angle in the system
    angle = (2.0*math.pi - angle)/2.0
    # Rotate that many radians counterclockwise around the y-axis
    cdef PyMatrix rotate = PyMatrix ()
    rotate.specified_rotation (angle, [0.0, 1.0, 0.0])
    do_positioning_movement (residue, rotate)
    # Store the opposite of that rotation
    rotate.specified_rotation (-angle, [0.0, 1.0, 0.0])
    movements.append(rotate)
    return movements

cdef void Interaction_21_Giver (PyResidue residue, list atoms):
    """Position the "giver" residue for Interaction 21"""
    # Do the standard positioning of the residue atoms to get the ring into
    # the XZ plane
    cdef list movements = standard_positioning (residue, atoms)
    # Move the ring so that it is centered at the origin
    cdef PyMatrix movement = PyMatrix ()
    movement.PyAtoms_allocation_for_centering(atoms)
    do_positioning_movement (residue, movement)
    # Rotate the ring so that it is in the XY plane
    movement.specified_rotation (math.radians(90.0), [1.0, 0.0, 0.0])
    do_positioning_movement (residue, movement)
    # Raise the ring by 4.25 angstroms
    movement.PyAtoms_allocation_for_centering([[0.0, 0.0, 4.25]])
    do_positioning_movement (residue, movement, False)

cdef void Interaction_21 (PyResidue residue, list atoms, str LoopName, 
          size_t StructNum, size_t PosIndex, size_t RotNum, dict Rotamers,
          list outputs):
    """Generate and output data for Interaction 21"""
    # Get the positioning of the residue
    cdef list movements = Interaction_21_Taker (residue, atoms)
    reverse_positioning_movements (residue, movements)
    # Get a PHE rotamer and its ring atoms
    cdef PyResidue rotamer = Rotamers["PHE"].duplicate()
    cdef list rotAtoms = get_atoms (rotamer, ['CZ', 'CE1', 'CD1', 'CG', 'CD2',
                                              'CE2'])
    # Position the rotamer
    Interaction_21_Giver (rotamer, rotAtoms)
    reverse_positioning_movements (rotamer, movements)
    # Store the needed information
    Center_and_Plane (rotAtoms, LoopName, StructNum, PosIndex, RotNum,
                      outputs[20])

# Interaction 22 is the opposite of Interaction 20 - an aromatic ring in the
# taker protein parallel to a histidine ring in the giver protein
cdef list Interaction_22_Taker (PyResidue residue, list atoms):
    """Position the "taker" residue for Interaction 22"""
    # The aromatic ring should be in the XZ plane centered at the origin.
    # This is the same motions as interaction 20 (and 15)
    return Interaction_20_Taker (residue, atoms)

cdef void Interaction_22_Giver (PyResidue residue, list atoms):
    """Position the "giver" residue for Interaction 22"""
    # These are the same motions as in Interaction 20
    Interaction_20_Giver (residue, atoms)

cdef void Interaction_22 (PyResidue residue, list atoms, str LoopName,
          size_t StructNum, size_t PosIndex, size_t RotNum, dict Rotamers,
          list outputs):
    """Generate and output data for Interaction 22"""
    # Position the taker residue
    cdef list movements = Interaction_22_Taker (residue, atoms)
    reverse_positioning_movements (residue, movements)
    # Get a HIS rotamer and it's ring atoms
    cdef PyResidue rotamer = Rotamers["HIS"].duplicate()
    cdef list rotAtoms = get_atoms (rotamer, ['CG', 'CD2', 'NE2', 'CE1', 'ND1'])
    # Position the rotamer
    Interaction_22_Giver (rotamer, rotAtoms)
    reverse_positioning_movements (rotamer, movements)
    # Store the needed information
    Center_and_Plane (rotAtoms, LoopName, StructNum, PosIndex, RotNum,
                      outputs[21])

# Interaction 23 is the opposite of Interaction 21 - an aromatic ring being
# "pointed" at by a histidine ring
cdef list Interaction_23_Taker (PyResidue residue, list atoms):
    """Position the "taker" residue for Interaction 23"""
    # Once again, we want the aromatic ring in the XZ plane centered at the
    # origin
    return Interaction_20_Taker (residue, atoms)

cdef void Interaction_23_Giver (PyResidue residue, list atoms):
    """Position the "giver" residue for Interaction 23"""
    # Start by pointing the "tip" atom of the Histidine ring up along the Z
    # axis as is done in Interaction 21
    cdef list movements = Interaction_21_Taker (residue, atoms)
    # Rotate around 90 degrees counterclockwise around the X axis
    cdef PyMatrix movement = PyMatrix ()
    movement.specified_rotation (math.radians(90.0), [1.0, 0.0, 0.0])
    do_positioning_movement(residue, movement)
    # Move the residue out 4.25 angstroms along the Y axis
    movement.PyAtoms_allocation_for_centering([[0.0, 4.25, 0.0]])
    do_positioning_movement (residue, movement, False)

cdef void Interaction_23 (PyResidue residue, list atoms, str LoopName,
          size_t StructNum, size_t PosIndex, size_t RotNum, dict Rotamers,
          list outputs):
    """Generate and output data for Interaction 23"""
    # Position the taker residue
    cdef list movements = Interaction_23_Taker (residue, atoms)
    reverse_positioning_movements (residue, movements)
    # Get a Histidine rotamer
    cdef PyResidue rotamer = Rotamers["HIS"].duplicate()
    # Either the CE1 or NE2 can point at the aromatic ring, but whichever is
    # used doesn't change the resulting calculations. I'll use NE2 as the
    # pointer in this case
    cdef list rotAtoms = get_atoms (rotamer, ['NE2', 'CD2', 'CE1'])
    # Position the rotamer
    Interaction_23_Giver (rotamer, rotAtoms)
    reverse_positioning_movements (rotamer, movements)
    # Get the atoms that make up the whole ring of the histidine
    rotAtoms = get_atoms (rotamer, ['NE2', 'CE1', 'ND1', 'CG', 'CD2'])
    # Output Tip and Plane information for the NE2 atom
    Tip_and_Center (rotAtoms, 'NE2', LoopName, StructNum, PosIndex, RotNum,
                    outputs[22])

# The next set of functions are for each individual amino acid and the
# corresponding Interactions they can have
cdef void Every_Residue (PyResidue residue, str L, size_t S, size_t P,
                         dict Rotamers, list outputs):
    """Every Residue participates in certain interactions"""
    # Because the interactions here are backbone interactions, use a rotamer
    # index of 0 to indicate that information
    cdef size_t RI = 0
    # Use the backbone N and H to do Interactions 1 and 2
    atoms = get_atoms (residue, ['N', 'HN'])
    # If that is a list, the calculations can be done
    if isinstance(atoms, list):
        Interaction_1(residue, atoms, L, S, P, RI, Rotamers, outputs)
        Interaction_2(residue, atoms, L, S, P, RI, Rotamers, outputs)
    # In the future, Interactions 3 and 4 (hydrogen bonds involving the double
    # bonded O in the backbone of the peptide) will be added here. Those
    # calculations do not yet exist

cdef void ALA (PyResidue residue, str rotlib, str L, size_t S, size_t P, 
               dict Rotamers, list outputs):
    """The Interactions of ALA amino acids"""
    Every_Residue (residue, L, S, P, Rotamers, outputs)

cdef void CYS (PyResidue residue, str rotlib, str L, size_t S, size_t P,
               dict Rotamers, list outputs):
    """The Interactions of CYS amino acids"""
    Every_Residue (residue, L, S, P, Rotamers, outputs)

cdef void ASP (PyResidue residue, str rotlib, str L, size_t S, size_t P,
               dict Rotamers, list outputs):
    """The Interactions of ASP amino acids"""
    # Start by doing the backbone interactions
    Every_Residue (residue, L, S, P, Rotamers, outputs)
    # Load the rotamers of the residue
    cdef list rotamers = residue.collect_rotamers(rotlib)
    # Get the number of rotamers and loop through them
    cdef size_t I, J
    J = len(rotamers)
    cdef PyResidue rotamer
    for I in range(J):
        rotamer = rotamers[I]
        # The OD1 and OD2 in the ASP can participate in hydrogen bonds.
        # These are interactions 3 and 4, which are not yet implemented
        # ASP can have strong interactions with ARG and LYS residues
        atoms = get_atoms(rotamer, ['CG', 'OD1', 'OD2'])
        Interaction_11 (rotamer, atoms, L, S, P, I+1, Rotamers, outputs)
        Interaction_12 (rotamer, atoms, L, S, P, I+1, Rotamers, outputs)

cdef void GLU (PyResidue residue, str rotlib, str L, size_t S, size_t P,
               dict Rotamers, list outputs):
    """The Interactions of GLU amino acids"""
    # This is effectively identical to the ASP function
    Every_Residue (residue, L, S, P, Rotamers, outputs)
    cdef list rotamers = residue.collect_rotamers(rotlib)
    cdef size_t I, J
    J = len(rotamers)
    cdef PyResidue rotamer
    for I in range(J):
        rotamer = rotamers[I]
        atoms = get_atoms(rotamer, ['CD', 'OE1', 'OE2'])
        Interaction_11 (rotamer, atoms, L, S, P, I+1, Rotamers, outputs)
        Interaction_12 (rotamer, atoms, L, S, P, I+1, Rotamers, outputs)

cdef void PHE (PyResidue residue, str rotlib, str L, size_t S, size_t P, 
               dict Rotamers, list outputs):
    """The Interactions of PHE amino acids"""
    # Do the every residue interactions
    Every_Residue (residue, L, S, P, Rotamers, outputs)
    # Get the rotamers of the residue
    cdef list rotamers = residue.collect_rotamers(rotlib)
    # Loop through the rotamers
    cdef size_t I, J
    J = len(rotamers)
    cdef PyResidue rotamer
    for I in range(J):
        rotamer = rotamers[I]
        # Get the ring atoms in 2 sets. This lets the "top" and "bottom" of
        # the aromatic ring be used for interactions
        atoms1 = get_atoms (rotamer, ['CZ', 'CE1', 'CD1', 'CG', 'CD2', 'CE2'])
        atoms2 = get_atoms (rotamer, ['CZ', 'CD1', 'CE1', 'CG', 'CD2', 'CE2'])
        Interaction_15 (rotamer, atoms1, L, S, P, I+1, Rotamers, outputs)
        Interaction_15 (rotamer, atoms2, L, S, P, I+1, Rotamers, outputs)
        Interaction_16 (rotamer, atoms1, L, S, P, I+1, Rotamers, outputs)
        Interaction_16 (rotamer, atoms2, L, S, P, I+1, Rotamers, outputs)
        Interaction_22 (rotamer, atoms1, L, S, P, I+1, Rotamers, outputs)
        Interaction_22 (rotamer, atoms2, L, S, P, I+1, Rotamers, outputs)
        Interaction_23 (rotamer, atoms1, L, S, P, I+1, Rotamers, outputs)
        Interaction_23 (rotamer, atoms2, L, S, P, I+1, Rotamers, outputs)
        # The atom lists now need to be expanded to include a 7th atom
        atoms1.extend(get_atoms(rotamer, ['CB']))
        atoms2.extend(get_atoms(rotamer, ['CB']))
        Interaction_17 (rotamer, atoms1, L, S, P, I+1, Rotamers, outputs)
        Interaction_17 (rotamer, atoms2, L, S, P, I+1, Rotamers, outputs)
        Interaction_18 (rotamer, atoms1, L, S, P, I+1, Rotamers, outputs)
        Interaction_18 (rotamer, atoms2, L, S, P, I+1, Rotamers, outputs)
        Interaction_19 (rotamer, atoms1, L, S, P, I+1, Rotamers, outputs)
        # Interaction 19 does not need to be duplicated with the second atom
        # set

cdef void GLY (PyResidue residue, str rotlib, str L, size_t S, size_t P,
               dict Rotamers, list outputs):
    """The interactions of glycine amino acids"""
    Every_Residue (residue, L, S, P, Rotamers, outputs)

cdef void HIS (PyResidue residue, str rotlib, str L, size_t S, size_t P, 
               dict Rotamers, list outputs):
    """The interactions of histidine amino acids"""
    # Do the backbone interactions
    Every_Residue (residue, L, S, P, Rotamers, outputs)
    # Get the rotamers of histidine
    cdef list rotamers = residue.collect_rotamers (rotlib)
    # Loop through the rotamers
    cdef size_t I, J
    J = len(rotamers)
    cdef PyResidue rotamer
    for I in range(J):
        rotamer = rotamers[I]
        # The ND1 - HD1 atoms in the HIS can participate in hydrogen bonds.
        # These are interactions 1 and 2
        atoms = get_atoms (rotamer, ['ND1', 'HD1'])
        Interaction_1 (rotamer, atoms, L, S, P, I+1, Rotamers, outputs)
        Interaction_2 (rotamer, atoms, L, S, P, I+1, Rotamers, outputs)
        # Histidine can participate in cation - pi interactions with aromatic
        # residues
        atoms = get_atoms (rotamer, ["NE2", 'CE1', 'ND1', 'CG', 'CD2'])
        Interaction_20 (rotamer, atoms, L, S, P, I+1, Rotamers, outputs)
        atoms = get_atoms (rotamer, ['NE2', 'ND1', 'CE1', 'CG', 'CD2'])
        Interaction_20 (rotamer, atoms, L, S, P, I+1, Rotamers, outputs)
        atoms = get_atoms (rotamer, ['NE2', 'CD2', 'CE1'])
        Interaction_21 (rotamer, atoms, L, S, P, I+1, Rotamers, outputs)
        atoms = get_atoms (rotamer, ['CE1', 'NE2', 'ND1'])
        Interaction_21 (rotamer, atoms, L, S, P, I+1, Rotamers, outputs)

cdef void ILE (PyResidue residue, str rotlib, str L, size_t S, size_t P,
               dict Rotamers, list outputs):
    """The interactions of isoleucine amino acids"""
    Every_Residue (residue, L, S, P, Rotamers, outputs)

cdef void LYS (PyResidue residue, str rotlib, str L, size_t S, size_t P,
               dict Rotamers, list outputs):
    """The interactions of lysine amino acids"""
    # Do the backbone interaction calculations
    Every_Residue (residue, L, S, P, Rotamers, outputs)
    # Get the rotamers of Lysine that are possible at this position
    cdef list rotamers = residue.collect_rotamers (rotlib)
    # Loop through the rotamers
    cdef size_t I, J
    J = len(rotamers)
    cdef PyResidue rotamer
    for I in range(J):
        rotamer = rotamers[I]
        # Do the cation interactions of Lysine first
        atoms = get_atoms (rotamer, ['NZ', 'CE'])
        Interaction_10 (rotamer, atoms, L, S, P, I+1, Rotamers, outputs)
        Interaction_14 (rotamer, atoms, L, S, P, I+1, Rotamers, outputs)
        # Each of the three hydrogens attached to the NZ can participate in
        # hydrogen bonds
        for name in ['HZ1', 'HZ2', 'HZ3']:
            atoms = get_atoms (rotamer, ['NZ', name])
            Interaction_1 (rotamer, atoms, L, S, P, I+1, Rotamers, outputs)
            Interaction_2 (rotamer, atoms, L, S, P, I+1, Rotamers, outputs)

cdef void LEU (PyResidue residue, str rotlib, str L, size_t S, size_t P, 
               dict Rotamers, list outputs):
    """The interactions of leucine amino acids"""
    Every_Residue (residue, L, S, P, Rotamers, outputs)

cdef void MET (PyResidue residue, str rotlib, str L, size_t S, size_t P,
               dict Rotamers, list outputs):
    """The interactions of methionine amino acids"""
    Every_Residue (residue, L, S, P, Rotamers, outputs)

cdef void ASN (PyResidue residue, str rotlib, str L, size_t S, size_t P,
               dict Rotamers, list outputs):
    """The interactions of asparagine amino acids"""
    # Do the backbone interactions
    Every_Residue (residue, L, S, P, Rotamers, outputs)
    # Get the rotamers of ASN and loop through them
    cdef list rotamers = residue.collect_rotamers (rotlib)
    # Loop through the rotamers
    cdef size_t I, J
    J = len(rotamers)
    cdef PyResidue rotamer
    for I in range(J):
        rotamer = rotamers[I]
        # The double-bonded side chain O can participate in interactions 3 and
        # 4, once those interactions are added to the calculations
        # The NH2 group can participate in hydrogen bonds in interactions 1
        # and 2
        for name in ['HD21', 'HD22']:
            atoms = get_atoms (rotamer, ['ND2', name])
            Interaction_1 (rotamer, atoms, L, S, P, I+1, Rotamers, outputs)
            Interaction_2 (rotamer, atoms, L, S, P, I+1, Rotamers, outputs)

cdef void PRO (PyResidue residue, str rotlib, str L, size_t S, size_t P,
               dict Rotamers, list outputs):
    """The interactions of Proline amino acids"""
    # Proline does not contain a backbone H-N, so it does not get handled like
    # every other residue
    # When double-bonded O=C are added to the program, Interactions 3 and 4
    # should be done here on the backbone atoms
    pass

cdef void GLN (PyResidue residue, str rotlib, str L, size_t S, size_t P,
               dict Rotamers, list outputs):
    """The interactions of glutamine amino acids"""
    # Do the backbone interactions
    Every_Residue (residue, L, S, P, Rotamers, outputs)
    # Get the rotamers and loop through them
    cdef list rotamers = residue.collect_rotamers (rotlib)
    cdef size_t I, J
    J = len(rotamers)
    cdef PyResidue rotamer
    for I in range(J):
        rotamer = rotamers[I]
        # The side-chain double bonded O (OE1 attached to CD) can participate
        # in interactions 3 and 4 once those interactions are added to the
        # program. 
        # The side chain NH2 group can participate in hydrogen bonds
        for name in ['HE21', 'HE22']:
            atoms = get_atoms (rotamer, ['NE2', name])
            Interaction_1 (rotamer, atoms, L, S, P, I+1, Rotamers, outputs)
            Interaction_2 (rotamer, atoms, L, S, P, I+1, Rotamers, outputs)

cdef void ARG (PyResidue residue, str rotlib, str L, size_t S, size_t P,
               dict Rotamers, list outputs):
    """The interactions of arginine amino acids"""
    # Do the backbone interactions
    Every_Residue (residue, L, S, P, Rotamers, outputs)
    # Get the rotamers and loop through them
    cdef list rotamers = residue.collect_rotamers (rotlib)
    cdef size_t I, J
    J = len(rotamers)
    cdef PyResidue rotamer
    for I in range(J):
        rotamer = rotamers[I]
        # The ARG guanidinium group can have many different interactions. The
        # relevant atoms in it are: CZ, (NE-HE), (NH1-(HH11,HH12)), 
        # (NH2-(HH21,HH22))
        # All of the single hydrogen bonds that can occur
        for pair in [['NE', 'HE'], ['NH1', 'HH11'], ['NH1', 'HH12'],
                     ['NH2', 'HH21'], ['NH2', 'HH22']]:
            atoms = get_atoms (rotamer, pair)
            Interaction_1 (rotamer, atoms, L, S, P, I+1, Rotamers, outputs)
            Interaction_2 (rotamer, atoms, L, S, P, I+1, Rotamers, outputs)
        # The strong interactions of ARG with GLU or ASP
        for three in [['CZ', 'NH1', 'NH2'], ['CZ', 'NH1', 'NE'], 
                      ['CZ', 'NH2', 'NE']]:
            atoms = get_atoms (rotamer, three)
            Interaction_9 (rotamer, atoms, L, S, P, I+1, Rotamers, outputs)
        # Cation - pi interactions. The aromatic ring can go on either side of
        # the ARG, so two sets of atoms are needed
        for three in [['CZ', 'NH1', 'NH2'], ['CZ', 'NH2', 'NH1']]:
            atoms = get_atoms (rotamer, three)
            Interaction_13 (rotamer, atoms, L, S, P, I+1, Rotamers, outputs)

cdef void SER (PyResidue residue, str rotlib, str L, size_t S, size_t P,
               dict Rotamers, list outputs):
    """The interactions of serine amino acids"""
    # The backbone interactions
    Every_Residue (residue, L, S, P, Rotamers, outputs)
    # Get the rotamers and loop through them
    cdef list rotamers = residue.collect_rotamers (rotlib)
    cdef size_t I, J
    J = len(rotamers)
    cdef PyResidue rotamer
    for I in range(J):
        rotamer = rotamers[I]
        # The side chain atoms that can have interactions are:
        atoms = get_atoms(rotamer, ['OG', 'HG1'])
        # Form hydrogen bonds
        Interaction_5 (rotamer, atoms, L, S, P, I+1, Rotamers, outputs)
        Interaction_6 (rotamer, atoms, L, S, P, I+1, Rotamers, outputs)
        # Accept hydrogen bonds
        atoms = get_atoms (rotamer, ['OG', 'CB', 'HG1'])
        Interaction_7 (rotamer, atoms, L, S, P, I+1, Rotamers, outputs)
        Interaction_8 (rotamer, atoms, L, S, P, I+1, Rotamers, outputs)

cdef void THR (PyResidue residue, str rotlib, str L, size_t S, size_t P,
               dict Rotamers, list outputs):
    """The interactions of threonine amino acids"""
    # This is effectively identical to the SER function
    Every_Residue (residue, L, S, P, Rotamers, outputs)
    cdef list rotamers = residue.collect_rotamers (rotlib)
    cdef size_t I, J
    J = len(rotamers)
    cdef PyResidue rotamer
    for I in range(J):
        rotamer = rotamers[I]
        atoms = get_atoms (rotamer, ['OG1', 'HG1'])
        Interaction_5 (rotamer, atoms, L, S, P, I+1, Rotamers, outputs)
        Interaction_6 (rotamer, atoms, L, S, P, I+1, Rotamers, outputs)
        atoms = get_atoms (rotamer, ['OG1', 'CB', 'HG1'])
        Interaction_7 (rotamer, atoms, L, S, P, I+1, Rotamers, outputs)
        Interaction_8 (rotamer, atoms, L, S, P, I+1, Rotamers, outputs)

cdef void VAL (PyResidue residue, str rotlib, str L, size_t S, size_t P,
               dict Rotamers, list outputs):
    """The interactions of valine amino acids"""
    Every_Residue (residue, L, S, P, Rotamers, outputs)

cdef void TRP (PyResidue residue, str rotlib, str L, size_t S, size_t P,
               dict Rotamers, list outputs):
    """The interactions of tryptophan amino acids"""
    # Do the backbone interactions
    Every_Residue (residue, L, S, P, Rotamers, outputs)
    # Get the rotamers and loop through them
    cdef list rotamers = residue.collect_rotamers (rotlib)
    cdef size_t I, J
    J = len(rotamers)
    cdef PyResidue rotamer
    for I in range(J):
        rotamer = rotamers[I]
        # The side chain NH can form hydrogen bonds
        atoms = get_atoms (rotamer, ['NE1', 'HE1'])
        Interaction_1 (rotamer, atoms, L, S, P, I+1, Rotamers, outputs)
        Interaction_2 (rotamer, atoms, L, S, P, I+1, Rotamers, outputs)
        # The aromatic ring can form a number of interactions as well. The
        # ring atoms are: CH2, CZ3, CE3, CD2, CE2, CZ2. CD2 is attached to CG.
        # Arguably, either CH2 or CZ3 could be the "tip" atom for a T-shaped
        # pi - pi interaction. However, for CZ3 the atom attached to CE2 is
        # NE1, and I don't know if that changes any considerations or not
        # The ring has 2 sides, so get two sets of atoms for considering its
        # interactions
        atoms1 = get_atoms (rotamer, ['CH2', 'CZ3', 'CE3', 'CD2', 'CE2', 'CZ2'])
        atoms2 = get_atoms (rotamer, ['CH2', 'CE3', 'CZ3', 'CD2', 'CE2', 'CZ2'])
        # Do Pi-cation interactions
        Interaction_15 (rotamer, atoms1, L, S, P, I+1, Rotamers, outputs)
        Interaction_15 (rotamer, atoms2, L, S, P, I+1, Rotamers, outputs)
        Interaction_16 (rotamer, atoms1, L, S, P, I+1, Rotamers, outputs)
        Interaction_16 (rotamer, atoms2, L, S, P, I+1, Rotamers, outputs)
        Interaction_22 (rotamer, atoms1, L, S, P, I+1, Rotamers, outputs)
        Interaction_22 (rotamer, atoms2, L, S, P, I+1, Rotamers, outputs)
        Interaction_23 (rotamer, atoms1, L, S, P, I+1, Rotamers, outputs)
        Interaction_23 (rotamer, atoms2, L, S, P, I+1, Rotamers, outputs)
        # Add the 7th atom for pi - pi interactions
        atoms1.extend(get_atoms (rotamer, ['CG']))
        atoms2.extend(get_atoms (rotamer, ['CG']))
        # Do the pi - pi interactions
        Interaction_17 (rotamer, atoms1, L, S, P, I+1, Rotamers, outputs)
        Interaction_17 (rotamer, atoms2, L, S, P, I+1, Rotamers, outputs)
        Interaction_18 (rotamer, atoms1, L, S, P, I+1, Rotamers, outputs)
        Interaction_18 (rotamer, atoms2, L, S, P, I+1, Rotamers, outputs)
        # Interaction 19 is at the "tip" of the ring and therefore only gets
        # calculated once
        Interaction_19 (rotamer, atoms1, L, S, P, I+1, Rotamers, outputs)

cdef void TYR (PyResidue residue, str rotlib, str L, size_t S, size_t P,
               dict Rotamers, list outputs):
    """The interactions of tyrosine residues"""
    # Backbone interactions
    Every_Residue (residue, L, S, P, Rotamers, outputs)
    # Get the rotamers and loop through them
    cdef list rotamers = residue.collect_rotamers (rotlib)
    cdef size_t I, J
    J = len(rotamers)
    cdef PyResidue rotamer
    for I in range(J):
        rotamer = rotamers[I]
        # The side chain OH can participate in hydrogen bonds
        atoms = get_atoms (rotamer, ['OH', 'HH'])
        Interaction_5 (rotamer, atoms, L, S, P, I+1, Rotamers, outputs)
        Interaction_6 (rotamer, atoms, L, S, P, I+1, Rotamers, outputs)
        atoms = get_atoms (rotamer, ['OH', 'CZ', 'HH'])
        Interaction_7 (rotamer, atoms, L, S, P, I+1, Rotamers, outputs)
        Interaction_8 (rotamer, atoms, L, S, P, I+1, Rotamers, outputs)
        # The ring atoms are: CZ, CE1, CD1, CG, CD2, CE2. CG is attached to CB
        # Do the pi-cation interactions
        atoms1 = get_atoms (rotamer, ['CZ', 'CE1', 'CD1', 'CG', 'CD2', 'CE2'])
        atoms2 = get_atoms (rotamer, ['CZ', 'CD1', 'CE1', 'CG', 'CD2', 'CE2'])
        Interaction_15 (rotamer, atoms1, L, S, P, I+1, Rotamers, outputs)
        Interaction_15 (rotamer, atoms2, L, S, P, I+1, Rotamers, outputs)
        Interaction_16 (rotamer, atoms1, L, S, P, I+1, Rotamers, outputs)
        Interaction_16 (rotamer, atoms2, L, S, P, I+1, Rotamers, outputs)
        Interaction_22 (rotamer, atoms1, L, S, P, I+1, Rotamers, outputs)
        Interaction_22 (rotamer, atoms2, L, S, P, I+1, Rotamers, outputs)
        Interaction_23 (rotamer, atoms1, L, S, P, I+1, Rotamers, outputs)
        Interaction_23 (rotamer, atoms2, L, S, P, I+1, Rotamers, outputs)
        # Do the pi-pi interactions
        atoms1.extend(get_atoms(rotamer, ['CB']))
        atoms2.extend(get_atoms(rotamer, ['CB']))
        Interaction_17 (rotamer, atoms1, L, S, P, I+1, Rotamers, outputs)
        Interaction_17 (rotamer, atoms2, L, S, P, I+1, Rotamers, outputs)
        Interaction_18 (rotamer, atoms1, L, S, P, I+1, Rotamers, outputs)
        Interaction_18 (rotamer, atoms2, L, S, P, I+1, Rotamers, outputs)
        # Because of the OH, Tyrosine can't "point" at another aromatic ring

# The final section of code in this file deals with actually finding the
# specific interactions in the binding loops. Start by setting up the
# pre-calculation information
def prepare_for_calculations (instructions):
    """Prepare for doing the Interaction calculations"""
    # Get the name of the summary file
    summaryName = summary_file_name (instructions)
    # Open it for reading
    summary = open(summaryName, "r")
    # Get its lines
    lines = summary.readlines()
    summary.close()
    # Store the number of binding loop structures in this dictionary
    counts = {}
    # Go through the lines
    for line in lines:
        if "usable" in line:
            items = line.split()
            counts[items[0]] = int(items[3])
    # Open the summary file for writing
    summary = open(summaryName, "a")
    # Write a message saying what is happening
    message = "Finding optimal interaction positions started on " + time_stamp()
    summary.write(message)
    summary.flush()
    # Get the name of the database folder
    database = database_folder_name (instructions)
    # the folder where the interaction information will be stored is
    folder = database + "Interactions/"
    # Try to make that folder and raise an error if that fails
    try:
        os.mkdir(folder)
    except OSError:
        error = "Failure to make the Interactions folder. It likely already "
        error += "exists.\n"
        raise BindingProteinError (error)
    # Create the Interaction files
    outputs = []
    for I in range(1, 24):
        fileName = folder + instructions["name"] + "_Interactions_" + str(I)
        outputs.append(open(fileName + ".txt", "w"))
    # Load representative rotamers for 19 of the 20 amino acids
    Rotamers = {}
    for name in ['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE',
            'LYS', 'LEU', 'MET', 'ASN', 'GLN', 'ARG', 'SER', 'THR', 'VAL',
            'TRP', 'TYR']:
        tag = name.lower()
        if name == "HIS":
            tag = "hsd"
        # The file is located at 
        fileName = instructions["rotamers"] + tag + "1.pdb"
        # Store the atoms here
        atoms = []
        # Go through the file
        f = open(fileName, "r")
        for line in f:
            if line.startswith("ATOM"):
                atoms.append(line.strip())
        # Close the file
        f.close()
        # Store a PyResidue
        Rotamers[name] = PyResidue (atoms)
    # Return the identified information
    return summary, counts, outputs, Rotamers

def Find_Interactions (instructions):
    """Find the optimal interaction positions"""
    # Do the preliminary calculations
    summary, counts, outputs, Rotamers = prepare_for_calculations(instructions)
    # Needed variables
    cdef size_t I1, I2, I3, J1, J2, J3
    J1 = len(instructions["loops"])
    cdef PyProtein protein
    cdef PyResidue residue
    cdef str LoopName
    cdef str fileName
    cdef str folder
    cdef str rotlib = instructions["rotamers"]
    cdef str name
    # Go through the binding loops
    for I1 in range(J1):
        LoopName = instructions["loops"][I1][0]
        # Get the number of structures of this binding loop
        J2 = counts[LoopName]
        # The folder for the binding loop's structures
        folder = database_folder_name(instructions) + LoopName + "/Structures/"
        # Go through the structures
        for I2 in range(1, J2+1):
            # The name of the structure file
            fileName = instructions["name"] + "_" + LoopName + "_" \
                     + str(I2) + ".pdb"
            # Load the PyProtein object
            protein = PyProtein (fileName, folder)
            # Calculate the dihedral angles of the Protein
            protein.calculate_dihedrals()
            # The number of Residues in the protein
            J3 = len(protein)
            # Go through the residues of the Protein
            for I3 in range(J3):
                # Get the residue
                residue = protein[I3]
                name = residue.name()
                # Use the appropriate function on the Residue
                if name == "ALA":
                    ALA (residue, rotlib, LoopName, I2, I3, Rotamers, outputs)
                elif name == "CYS":
                    CYS (residue, rotlib, LoopName, I2, I3, Rotamers, outputs)
                elif name == "ASP":
                    ASP (residue, rotlib, LoopName, I2, I3, Rotamers, outputs)
                elif name == "GLU":
                    GLU (residue, rotlib, LoopName, I2, I3, Rotamers, outputs)
                elif name == "PHE":
                    PHE (residue, rotlib, LoopName, I2, I3, Rotamers, outputs)
                elif name == "GLY":
                    GLY (residue, rotlib, LoopName, I2, I3, Rotamers, outputs)
                elif name == "HIS":
                    HIS (residue, rotlib, LoopName, I2, I3, Rotamers, outputs)
                elif name == "ILE":
                    ILE (residue, rotlib, LoopName, I2, I3, Rotamers, outputs)
                elif name == "LYS":
                    LYS (residue, rotlib, LoopName, I2, I3, Rotamers, outputs)
                elif name == "LEU":
                    LEU (residue, rotlib, LoopName, I2, I3, Rotamers, outputs)
                elif name == "MET":
                    MET (residue, rotlib, LoopName, I2, I3, Rotamers, outputs)
                elif name == "ASN":
                    ASN (residue, rotlib, LoopName, I2, I3, Rotamers, outputs)
                elif name == "PRO":
                    PRO (residue, rotlib, LoopName, I2, I3, Rotamers, outputs)
                elif name == "GLN":
                    GLN (residue, rotlib, LoopName, I2, I3, Rotamers, outputs)
                elif name == "ARG":
                    ARG (residue, rotlib, LoopName, I2, I3, Rotamers, outputs)
                elif name == "SER":
                    SER (residue, rotlib, LoopName, I2, I3, Rotamers, outputs)
                elif name == "THR":
                    THR (residue, rotlib, LoopName, I2, I3, Rotamers, outputs)
                elif name == "VAL":
                    VAL (residue, rotlib, LoopName, I2, I3, Rotamers, outputs)
                elif name == "TRP":
                    TRP (residue, rotlib, LoopName, I2, I3, Rotamers, outputs)
                elif name == "TYR":
                    TYR (residue, rotlib, LoopName, I2, I3, Rotamers, outputs)
    # Close the output files
    for output in outputs:
        output.close()
    # Write a message to the summary file saying that database generation is
    # finished
    message = "Database Generation finished on " + time_stamp() + "\n\n"
    summary.write(message)
    summary.flush()
