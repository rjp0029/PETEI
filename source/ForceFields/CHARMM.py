# Written in 2018 by Dr. Robert Pantazes. Copyright owned by Auburn
# University.
# This file contains functions that create CHARMM scripts, run them, and
# collect the results.

# Include pure Python modules
import os
import sys
# Include the error class for problems with force field calculations
from ForceFields.EnergyError import EnergyError
# Include the Protein extension types
from Proteins.Protein.PyProtein import PyProtein

def execute_script (scriptFile, outputFile):
    """Actually run a CHARMM script"""
    # This is the command
    command = "/home/rjp0029/charmm/exec/gnu/charmm < " + scriptFile + " > "
    command += outputFile
    # Run it
    os.system(command)

def validate_proteins (proteins):
    """Make sure that the provided inputs are acceptable"""
    # This is OK if the object is a PyProtein
    if isinstance(proteins, PyProtein):
        return
    # There is also an acceptable option if it is a list (normal behaviour)
    elif isinstance(proteins, list):
        # Go through each object and make sure it is a protein
        for protein in proteins:
            if not isinstance(protein, PyProtein):
                error = "The ForceFields.CHARMM module is only intended to "
                error += "work with PyProtein objects\n"
                raise EnergyError (error)
        # Make sure each protein has a unique name
        if len(proteins) == 1:
            return
        for i in range(len(proteins) - 1):
            for j in range(i+1, len(proteins)):
                if proteins[i].name() == proteins[j].name():
                    error = "The ForceFields.CHARMM module cannot work with "
                    error += "two proteins with the same name\n"
                    raise EnergyError (error)
        return
    # Under any other circumstance, raise an error
    error = "The ForceFields.CHARMM module is only intended to work with "
    error += "PyProtein objects\n"
    raise EnergyError (error)

def introduction (purpose = None):
    """Create an introduction to a CHARMM script"""
    # Store the header information in this string
    header = ""
    # If there is a purpose
    if isinstance(purpose, str):
        header += "* " + purpose.strip() + "\n"
    # Otherwise make a generic statement
    else:
        header += "* CHARMM Calculation of Protein Structures\n"
    # The typical behaviours needed in CHARMM need the warning and bomb levels
    # at -2
    header += "\nwrnl -2\nbomb -2\n\n"
    return header

def load_topology_files (fileList):
    """Load the Topology files in CHARMM"""
    # Store the text that does the loading here
    text = "! Load the Topology Files\n"
    # Go through each file
    for I, fileName in enumerate(fileList):
        # Tell CHARMM to open the file
        text += "open read unit 10 form name " + fileName + " card\n"
        # Read in the topology information
        text += "read rtf card unit 10"
        # Possibly append the information if this is not the first topology
        # file
        if I > 0:
            text += " append"
        # Close the topology file
        text += "\nclose unit 10\n\n"
    return text

def load_parameter_files (fileList):
    """Load the Parameter files in CHARMM"""
    # This is the same as the previous function, but uses "para" instead of
    # "rtf"
    text = "! Load the Parameter Files\n"
    for I, fileName in enumerate(fileList):
        text += "open read unit 10 form name " + fileName + " card\n"
        text += "read para card unit 10"
        if I > 0:
            text += " append"
        text += "\nclose unit 10\n\n"
    return text

def input_proteins (proteins):
    """Create text to load proteins into a CHARMM script"""
    # If the proteins is a protein by itself, put it in a list
    if isinstance(proteins, PyProtein):
        useList = [proteins]
    # The other option from the validate script is that it is already a list
    # of proteins
    else:
        useList = proteins
    # The list of files created by this function
    fileNames = []
    # The text listing the commands to load the proteins
    text = ""
    # Renumber residues and atoms as the proteins are output
    rn = 1
    an = 1
    # Go through the proteins
    for protein in useList:
        # Renumber the residues and atoms sequentially
        rn = protein.renumber_residues(rn)
        an = protein.renumber_atoms (an)
        # Create a name for the protein file
        fileName = "protein_"
        if protein.name() != "_":
            fileName += protein.name().lower() + "_"
        fileName += "for_charmm.pdb"
        # Output the protein to that file using its internal numbering after
        # making sure histidines are listed as HSD
        protein.charmm_his_prep()
        protein.output(fileName, True)
        # Store the file name
        fileNames.append(fileName)
        # Add commands to the CHARMM script
        text += "! Load Protein " + protein.name() + "\n"
        # Open the protein in CHARMM to read in its sequence
        text += "open read unit 10 form name " + fileName + "\n"
        text += "read sequ pdb offi unit 10\nclose unit 10\ngene pr"
        if protein.name() != " ":
            text += protein.name().lower()
        text += " setup\nopen read unit 10 form name " + fileName + "\n"
        text += "read coor pdb unit 10 \nclose unit 10\n\n"
    return fileNames, text

def add_missing_atoms ():
    """Add missing atoms to the CHARMM structures"""
    text = "! Add missing atoms\nic fill preserve\nic param\nic build\nhbuild"
    return text + "\n\n"

def solvation_line ():
    """The line of text to use GBMV implicit solvation in CHARMM"""
    text = """GBMV BETA -20 EPSILON 80 DN 1.0 watr 1.4 GEOM -
    TOL 1e-8 BUFR 0.5 Mem 10 CUTA 20 HSX1 -0.125 HSX2 0.25 -
    ALFRQ 1 EMP 0.25 P4 0.0 P6 8.0 P3 0.70 ONX 1.9 OFFX 2.1 -
    WTYP 2 NPHI 38 SHIFT -0.102 SLOPE 0.9085 CORR 1
"""
    return text

def energy_minimization (Solvation = True):
    """Create text to do an energy minimization in CHARMM"""
    # Store the text here
    text = "! Conduct an energy minimization\n"
    # Add the commands to use Generalized Borne with Molecular Volume
    # integration implicit solvation
    if Solvation:
        text += solvation_line()
    # Add the remainder of the text
    text += "\nnbon nbxm 5\n"
    text += "skip all excl angl bond dihe elec impr urey vdw"
    if Solvation:
        text += " gbener"
    text += "\nmini abnr nstep 5000 nprint 10 -\n"
    text += " tolgrd 0.01 tolenr 0.0001 tolstp 0.00\n\n"
    return text

def energy_calculation (Solvation = True):
    """Create text to do an energy calculation in CHARMM"""
    # Store the text to do the calculation here
    text = "! Commands to do an energy calculation\n"
    # Include the solvation command
    if Solvation:
        text += solvation_line()
    # Tell CHARMM to use option 5 in the non-bonded exclusion list 
    text += "\nnbon nbxm 5\n"
    # Include all appropriate energy terms
    text += "skip all excl angl bond dihe elec impr urey vdw"
    if Solvation:
        text += " gbener"
    text += "\n\nener\n"
    return text

def output_proteins (proteins):
    """Create text that tells CHARMM how to output proteins"""
    # make sure this function definitely uses a list of proteins
    if isinstance(proteins, PyProtein):
        useList = [proteins]
    else:
        useList = proteins
    # Store the text here
    text = ""
    # Store the list of file names generated by this function here
    fileNames = []
    # Go through the proteins
    for protein in useList:
        # Create a name for the protein
        fileName = "protein_"
        if protein.name() != "_":
            fileName += protein.name().lower() + "_"
        fileName += "from_charmm.pdb"
        # Store the file name
        fileNames.append(fileName)
        # Create the text that outputs the final structures of the protein to
        # the file
        text += "! Output Protein " + protein.name() + "\n"
        text += "open write unit 10 name " + fileName + " card\n"
        text += "write coor sele segi pr"
        if protein.name() != ' ':
            text += protein.name().lower()
        text += " end pdb unit 10 card\nclose unit 10\n\n"
    return fileNames, text

def load_proteins (proteins, fileNames):
    """Load the proteins that are output by CHARMM"""
    # Make sure this function definitely uses a list of proteins
    if isinstance(proteins, PyProtein):
        useList = [proteins]
    else:
        useList = proteins
    # Go through the proteins and file names
    for I in range(len(useList)):
        # Have the protein load the file from the current folder
        useList[I].load(fileNames[I], "./")
        # Correct the histidine labelling of the residues
        useList[I].charmm_his_fix()

def clean_up (fileNames):
    """Delete the files created by a CHARMM process"""
    # In the event that a name fails to be deleted, it is OK
    for name in fileNames:
        try:
            os.remove(name)
        except OSError:
            pass

def Missing_Atoms (proteins, topologies, parameters, deleteAll = True):
    """Add missing atoms to the proteins"""
    # Make sure the proteins are OK
    validate_proteins (proteins)
    # Create the introduction
    script = introduction ("Add Missing Atoms to Proteins")
    # Load the topology files
    script += load_topology_files (topologies)
    # Load the parameter files
    script += load_parameter_files (parameters)
    # Load the protein structures
    fileNames, text = input_proteins (proteins)
    script += text
    # Add missing atoms
    script += add_missing_atoms ()
    # Output the protein structures
    outputNames, text = output_proteins (proteins)
    script += text
    # End the script
    script += "stop\n"
    # Write the script to a file
    scriptName = "add_missing_atoms.inp"
    f = open(scriptName, "w")
    f.write(script)
    f.close()
    # The output name for the calculations is
    outputName = "missing_atoms_results.out"
    # Run charmm
    execute_script (scriptName, outputName)
    # Load the protein structures
    load_proteins (proteins, outputNames)
    # Definitely delete the PDB files
    fileNames.extend(outputNames)
    # Possibly delete the CHARMM scripts
    if deleteAll:
        fileNames.append(scriptName)
        fileNames.append(outputName)
    clean_up (fileNames)

def Minimization (proteins, topologies, parameters, deleteAll = True, 
                  Solvation = True):
    """Do an energy minimization of a protein complex"""
    # This is the same as the missing atoms function, except it includes the
    # energy minimization text in the script
    validate_proteins (proteins)
    script = introduction ("Minimize the energy of proteins")
    script += load_topology_files (topologies)
    script += load_parameter_files (parameters)
    fileNames, text = input_proteins (proteins)
    script += text
    script += add_missing_atoms()
    script += energy_minimization (Solvation)
    outputNames, text = output_proteins (proteins)
    script += text + "stop\n"
    scriptName = "energy_minimization.inp"
    f = open(scriptName, "w")
    f.write(script)
    f.close()
    outputName = "energy_minimization_results.out"
    execute_script (scriptName, outputName)
    load_proteins (proteins, outputNames)
    fileNames.extend(outputNames)
    if deleteAll:
        fileNames.append(scriptName)
        fileNames.append(outputName)
    clean_up (fileNames)

def Energy (proteins, topologies, parameters, deleteAll = True, 
            Solvation = True):
    """Calculate the energy of the provided protein structures"""
    # Validate that the provided proteins are acceptable
    validate_proteins (proteins)
    # Start the script
    script = introduction ("Calculate the energy of the proteins")
    # Load the topology, parameter and protein files
    script += load_topology_files (topologies)
    script += load_parameter_files (parameters)
    fileNames, text = input_proteins (proteins)
    script += text
    # Do NOT add missing atoms - this function calculates the energy of the
    # proteins in the form they are provided. 
    script += energy_calculation (Solvation)
    # The energy values are written to this file
    energyFile = "charmm_energy_value.txt"
    # Write the energy to that file and end the script
    script += "\nset tot ?ener\n\nopen write card unit 10 name " + energyFile
    script += "\n\nwrite title unit 10\n* @tot\n*\nclose unit 10\nstop\n"
    # Write the script to a file
    scriptName = "energy_calculation.inp"
    f = open(scriptName, "w")
    f.write(script)
    f.close()
    # The output of CHARMM will go to this file
    outputName = "energy_calculation_details.out"
    # Run CHARMM
    execute_script (scriptName, outputName)
    # Load the calculated energy
    try:
        f = open(energyFile, "r")
        energy = float(f.readline())
        f.close()
    except IOError:
        text = "A CHARMM energy calculation failed to calculate an energy. "
        text += "Please review " + outputName + " to determine why.\n"
        raise EnergyError (text)
    # Delete the files created for CHARMM
    if deleteAll:
        fileNames.append(scriptName)
        fileNames.append(outputName)
        fileNames.append(energyFile)
    clean_up(fileNames)
    # Return the calculated energy
    return energy

def Interaction_Energy (design, target, topologies, parameters, 
                        deleteAll = True, Solvation = True):
    """Calculate an interaction energy of a protein complex"""
    # Calculate the design molecules energy
    designEnergy = Energy(design, topologies, parameters, deleteAll, Solvation)
    # Calculate the target molecules energy
    targetEnergy = Energy(target, topologies, parameters, deleteAll, Solvation)
    # Calculate the complex energy
    group = []
    group.extend(design)
    group.extend(target)
    complexEnergy = Energy(group, topologies, parameters, deleteAll, Solvation)
    # Return the interaction energy
    return complexEnergy - (designEnergy + targetEnergy)
