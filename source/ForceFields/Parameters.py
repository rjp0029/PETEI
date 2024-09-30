# Written in 2018 by Dr. Robert Pantazes. Copyright owned by Auburn
# University.
# This file contains the implementation of the Parameters class, which is a
# container of parameters associated with force fields and energy calculations

# Include the standard python modules
import os
import sys
# Include the energy error class
from ForceFields.EnergyError import EnergyError

# The information is stored in this class on a per-atom basis
class Parameter (object):
    """A collection of Atom force field parameters"""

    def __init__ (self, name, res):
        """The initialization function of the Parameter class"""
        # Initialize the Parameter using its Atom and Residue names
        self.name = name
        self.residue = res

    def set_attribute (self, label, value):
        """Set an attribute value for the Atom"""
        self.__dict__[label] = value

# This function initializes Parameters from CHARMM topology files
def load_CHARMM_topology_files (fileNames, overwrite = False):
    """Load Parameter information from CHARMM topology files"""
    # Store the parameters here
    results = {}
    # If the fileNames are a list
    if isinstance(fileNames, list):
        use = fileNames
    # Or if it is a string
    elif isinstance(fileNames, str):
        use = [fileNames]
    # Otherwise, raise an error
    else:
        error = "The load_CHARMM_topology_files function only works with "
        error += "string or list inputs.\n"
        raise EnergyError (error)
    # Loop through the files
    for fileName in use:
        # If it is not a string, raise an error
        if not isinstance(fileName, str):
            error = "Topology files must be specified as strings.\n"
            raise EnergyError (error)
        # Open the file
        try:
            f = open(fileName, "r")
        except IOError:
            text = "Failure to open " + fileName + "\n"
            raise EnergyError (error)
        # Indicate that information is not yet being stored
        active = False
        # Loop through the file's lines
        for line in f:
            # If it starts with "RESI" or "PRES", it indicates that a Residue
            # or Patched Residue information is available
            if line.upper().startswith(("RESI", "PRES")):
                # Split the line into pieces
                pieces = line.split()
                # The residue name is the first entry
                try:
                    res = pieces[1]
                except IndexError:
                    active = False
                    continue
                # Error check on whether or not the residue is already in the
                # results
                if res in results and not overwrite:
                    error = "Multiple topology entries for " + res + "\n"
                    raise EnergyError (error)
                # If there is not already an entry, create an empty dictionary
                if res not in results:
                    results[res] = {}
                # Indicate that information is actively being stored
                active = True
            # If information is actively being loaded and the entry is about
            # an atom
            elif active and line.upper().startswith("ATOM"):
                # Split the line into pieces
                pieces = line.split()
                # The Atom's name is
                name = pieces[1]
                # Error check on that name
                if not overwrite and name in results[res]:
                    error = "Multiple topology entries for Atom " + name
                    error += "in residue " + res + "\n"
                    raise EnergyError (error)
                # Create a Parameter
                results[res][name] = Parameter(name, res)
                # Set the Atom's chemical type
                results[res][name].set_attribute("chemical_type", pieces[2])
                # Set the Atom's charge
                results[res][name].set_attribute("elec_charge", float(pieces[3]))
        # Close the file
        f.close()
    # Return the generated results
    return results

# This function loads information from CHARMM parameter files
def load_CHARMM_vdw_values (fileNames, overwrite = False):
    """Load CHARMM VDW information from parameter files"""
    # The values will be stored here
    values = {}
    # If the file names are a list
    if isinstance(fileNames, list):
        use = fileNames
    # Or a string
    elif isinstance(fileNames, str):
        use = [fileNames]
    # Otherwise, raise an error
    else:
        error = "The load_CHARMM_vdw_values function only works with string "
        error += "or list inputs.\n"
        raise EnergyError (error)
    # Loop through the files
    for fileName in use:
        # Make sure it is a string
        if not isinstance(fileName, str):
            error = "CHARMM parameter files must be specified as strings.\n"
            raise EnergyError (error)
        # Try to open the file
        try:
            f = open(fileName, "r")
        except IOError:
            error = "Failure to open " + fileName + "\n"
            raise EnergyError (error)
        # Indicate that the non-bonded portion of the file has not yet been
        # reached
        active = False
        # Loop through the file's lines
        for line in f:
            # Strip whitespace and capitalize the line
            text = line.strip().upper()
            # If it is empty, skip it
            if len(text) == 0:
                continue
            # Or if it indicates that the non-bonded section has been reached
            elif text.startswith(("NONBONDED", "NBONDED")):
                active = True
            # If the non-bonded section has been reached and this isn't a line
            # that should be skipped
            elif active and not text.startswith(("!", "CUTNB", "END", "HBOND")):
                # Split the line into pieces
                pieces = text.split()
                # If there are not at least 4 pieces, raise an error
                if len(pieces) < 4:
                    error = "Improperly formatted NONBONDED line:\n" + line
                    raise EnergyError (error + "\n")
                # Get the chemical type of the atom
                kind = pieces[0]
                # If this chemical type is already in the dictionary, raise an
                # error
                if kind in values and not overwrite:
                    error = "Multiple " + kind + " NONBONDED entries identified.\n"
                    raise EnergyError (error)
                # Store the information
                values[kind] = {"vdwE":float(pieces[2]), "vdwR":float(pieces[3])}
                # If there are 1-4 parameters
                if len(pieces) >= 7 and not pieces[4].startswith("!"):
                    values[kind]["vdw14E"] = float(pieces[5])
                    values[kind]["vdw14R"] = float(pieces[6])
                # Otherwise, use the default values
                else:
                    values[kind]["vdw14E"] = values[kind]["vdwE"]
                    values[kind]["vdw14R"] = values[kind]["vdwR"]
        # Close the file
        f.close()
    # Return the values
    return values

# Load Lazaridis-Karplus implicit solvation parameters
def load_LK_Solvation_files (fileNames, overwrite = False):
    """Load Lazaridis-Karplus Implicit Solvation parameters"""
    # The results will be stored here
    results = {}
    # If the file names are a list
    if isinstance(fileNames, list):
        use = fileNames
    # Or a string
    elif isinstance(fileNames, str):
        use = [fileNames]
    # Otherwise, raise an error
    else:
        error = "The load_LK_Solvation_files function only works with "
        error += "string or list inputs.\n"
        raise EnergyError (error)
    # Loop through the files
    for fileName in use:
        # Make sure it is a string
        if not isinstance(fileName, str):
            error = "LK Solvation files must be specified as strings.\n"
            raise EnergyError (error)
        # Try to open the file
        try:
            f = open(fileName, "r")
        except IOError:
            error = "Failure to open " + fileName + "\n"
            raise EnergyError (error)
        # Loop through the contents of the file
        for line in f:
            # Split the line into pieces
            pieces = line.split()
            # If there are not the proper number of pieces or it is the header
            # line, skip it
            if len(pieces) != 6 or pieces[0] == "ATOM":
                continue
            # Get the residue and possibly create a residue entry in the
            # results
            res = pieces[1]
            if res not in results:
                results[res] = {}
            # Get the Atom's name
            name = pieces[0]
            # Possibly raise an error
            if name in results[res] and not overwrite:
                error = "Multiple LK Solvation entries for Atom " + name 
                error += " in Residue " + res + "\n"
                raise EnergyError (error)
            # Create a Parameter
            results[res][name] = Parameter(name, res)
            # Store its attribute information
            results[res][name].set_attribute("lkR", float(pieces[2]))
            results[res][name].set_attribute("lkL", float(pieces[3]))
            results[res][name].set_attribute("lkV", float(pieces[4]))
            results[res][name].set_attribute("lkG", float(pieces[5]))
        # Close the file
        f.close()
    # Return the results
    return results

# Load all energy calculation parameterization values
def load_parameters (topologies, parameters, LKSolvation = None, 
                     overwrite = False):
    """Load all of the non-bonded parameters for energy calculations"""
    # Load the topology results
    results = load_CHARMM_topology_files (topologies, overwrite)
    # Load the parameter results
    vdw = load_CHARMM_vdw_values (parameters, overwrite)
    # Loop through the residues
    for res in results:
        # Loop through the atoms
        for name in results[res]:
            # Get the chemical type of the atom
            kind = results[res][name].chemical_type
            # If that information is missing from the vdw data, raise an error
            if kind not in vdw:
                error = "No VDW parameters were identified for atoms of "
                error += "chemical type: " + kind + "\n"
                error += "Atom " + name + " in Residue " + res + " is one "
                error += "such Atom.\n"
                raise EnergyError (error)
            # Store the VDW parameters of the Atom
            for label in ['vdwE', 'vdwR', 'vdw14E', 'vdw14R']:
                results[res][name].set_attribute(label, vdw[kind][label])
    # If the LK Solvation information should be loaded
    if LKSolvation != None:
        LK = load_LK_Solvation_files(LKSolvation, overwrite)
        # Loop through the residues
        for res in LK:
            # Skip residues not in the CHARMM results
            if res not in results:
                continue
            # And its Atoms
            for atom in LK[res]:
                # Skip atoms not in the CHARMM results
                if atom not in results[res]:
                    continue
                # Store the LK parameters in the results atom
                for label in ['lkR', 'lkL', 'lkV', 'lkG']:
                    results[res][atom].set_attribute(label,LK[res][atom].__dict__[label])
    return results
