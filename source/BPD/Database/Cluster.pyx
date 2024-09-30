# Written in 2018 by Dr. Robert Pantazes. Copyright owned by Auburn
# University.
# This file contains classes and functions useful for clustering protein
# structures into groups with similar backbones.

# Cython compiler directives
# cython: c_string_type=str
# cython: c_string_encoding=ascii
# cython: language_level=3
# distutils: language=c++

# Include pure Python modules
import os
import sys
import math
# Include Binding Protein Design contents
from BPD.Instructions import time_stamp
from BPD.BindingProteinError import BindingProteinError
from BPD.Database.Preliminary import database_folder_name
from BPD.Database.Preliminary import summary_file_name
# Include the PyProtein class
from Proteins.Protein.PyProtein cimport PyProtein
from Proteins.Residue.PyResidue cimport PyResidue
# Include C++ classes
from libcpp.vector cimport vector
from libcpp cimport bool

cdef class Protein:
    """A class for Protein information"""
    
    # The class attributes
    cdef public PyProtein structure
    cdef readonly str name
    cdef readonly vector[float] angles

    def __init__ (self, fileName, folder):
        """Initialize a Protein object"""
        # Load the PyProtein
        self.structure = PyProtein(fileName, folder)
        # Set the name of the protein
        self.name = fileName.split(".")[0]
        # Move the structure so that it is at the origin
        how = self.structure.center(True)
        # Calculate the dihedral angles
        self.structure.calculate_dihedrals(True)
        # Store the sines and cosines of the phi, psi, and omega dihedral
        # angles
        self.angles.clear()
        cdef float value
        # Looping variables
        cdef size_t i, j
        j = len(self.structure)
        # Go through the residues
        cdef PyResidue res
        for i in range(j):
            res = self.structure[i]
            # If the residue has a phi dihedral angle
            if i > 0:
                phi = math.radians(res.phi())
                value = math.sin(phi)
                self.angles.push_back(value)
                value = math.cos(phi)
                self.angles.push_back(value)
            # If the residue has psi and omega dihedral angles
            if i < j-1:
                psi = math.radians(res.psi())
                value = math.sin(psi)
                self.angles.push_back(value)
                value = math.cos(psi)
                self.angles.push_back(value)
                omega = math.radians(res.omega())
                value = math.sin(omega)
                self.angles.push_back(value)
                value = math.cos(omega)
                self.angles.push_back(value)

    def __len__ (self):
        """The number of residues in the Protein"""
        return len(self.structure)
    
    def __str__ (self):
        """The names of the Protein's Residues"""
        cdef str answer = ""
        cdef size_t i, j
        j = len(self.structure)
        for i in range(j):
            if i > 0:
                answer += " "
            answer += self.structure[i].name()
        return answer

cdef class Loop:
    """Information about a Binding Loop"""

    # The attributes of the class
    cdef readonly str name
    cdef readonly str folder
    cdef str pdb_folder

    def __init__ (self, label, instructions):
        """Initialize a Loop object"""
        # Store the attributes of the class
        self.name = label
        path = database_folder_name(instructions) + self.name + "/"
        self.folder = path + "Clustering/"
        self.pdb_folder = path + "PDB/"
        # Create the clustering folder
        try:
            os.mkdir(self.folder)
        except OSError:
            error = "Failure to create the Clustering folder for Binding Loop "
            error += self.name + ". It likely already exists.\n"
            raise BindingProteinError (error)
    
    def load_proteins (self):
        """Load the identified structures for the binding loop"""
        # Get the names here
        names = os.listdir(self.pdb_folder)
        # They will be stored by lengths
        lengths = []
        loops = {}
        # Go through the files
        for name in names:
            # Load a Protein object
            protein = Protein(name, self.pdb_folder)
            # Get it's length
            L = len(protein)
            if L not in lengths:
                lengths.append(L)
                loops[L] = []
            # Store the protein
            loops[L].append(protein)
        # Return the binding loop structures and the observed lengths
        lengths.sort()
        return loops, lengths

cdef dict calculate_RMSDs (list proteins):
    """Calculate the RMSDs between proteins"""
    # Store the RMSD values here
    cdef dict rmsds = {}
    # Looping variables over the number of proteins
    cdef size_t i1, i2, j
    j = len(proteins)
    # Protein objects
    cdef Protein protein1, protein2
    # Initialize the RMSD dictionary
    for i1 in range(j):
        # Get the protein
        protein1 = proteins[i1]
        # Set up it's entry as being a dictionary
        rmsds[protein1.name] = {}
        # It's RMSD with itself is 0
        rmsds[protein1.name][protein1.name] = 0.0
    # Go through all pairs of proteins to calculate the RMSDs
    for i1 in range(j-1):
        protein1 = proteins[i1]
        for i2 in range(i1+1, j):
            protein2 = proteins[i2]
            # Calculate the RMSD
            rmsd = protein1.structure.calculate_RMSD(protein2.structure, True)
            # Store the RMSD information
            rmsds[protein1.name][protein2.name] = rmsd
            rmsds[protein2.name][protein1.name] = rmsd
    return rmsds

cdef vector[float] calculate_cluster_average (list proteins):
    """Calculate the average dihedral angles of the proteins"""
    # Store the values here
    cdef vector[float] average
    # The list of proteins contains Proteins
    cdef Protein protein
    # Set it as the first protein
    protein = proteins[0]
    # Looping variables
    cdef size_t I1, I2, J1, J2
    J2 = protein.angles.size()
    # Reserve that much space in the average vector
    average.reserve(J2)
    # Assign 0s to the values
    for I2 in range(J2):
        average.push_back(0.0)
    # Go through the proteins
    J1 = len(proteins)
    for I1 in range(J1):
        protein = proteins[I1]
        # Go through the angles
        for I2 in range(J2):
            average[I2] += protein.angles[I2]
    # Convert the sums to averages
    cdef float N = J1
    for I2 in range(J2):
        average[I2] /= N
    return average

cdef Protein select_model_protein (list proteins):
    """Select the protein that is most similar to the cluster averages"""
    # Calculate the averages of the cluster phi, psi and omega dihderal angles
    cdef vector[float] average = calculate_cluster_average(proteins)
    # Use these variables to store the best results
    cdef Protein best
    cdef float bestScore
    # Looping variables
    cdef size_t I1, I2, J1, J2
    J1 = len(proteins)
    J2 = average.size()
    # Go through the proteins
    cdef Protein protein
    cdef float score
    for I1 in range(J1):
        protein = proteins[I1]
        score = 0.0
        # Go through the angle values
        for I2 in range(J2):
            score += (average[I2] - protein.angles[I2])**2
        # If this result should be stored
        if I1 == 0 or score < bestScore:
            best = protein
            bestScore = score
    # Return the representative structure
    return best

cdef float score_cluster (Protein model, dict rmsds, list proteins):
    """Find the structure with the worst RMSD with the model"""
    # Looping variables
    cdef size_t I, J
    J = len(proteins)
    # The worst value is
    cdef float worst = 0.0
    # An RMSD value is
    cdef float rmsd
    # Use the model rmsd information
    cdef dict use = rmsds[model.name]
    # Go through the proteins
    for I in range(J):
        rmsd = use[proteins[I].name]
        if rmsd > worst:
            worst = rmsd
    return worst

cdef class Cluster:
    """A group of proteins with similar backbones"""

    # The class attributes
    cdef public size_t index
    cdef public bool active
    cdef public list proteins
    cdef public list scores
    cdef public float rmsd

    def __init__ (self, Protein protein, size_t I, size_t N):
        """Initialize a cluster"""
        # Setup the class's attributes
        self.index = I
        self.active = True
        self.proteins = [protein]
        self.rmsd = 0.0
        self.scores = []
        # The scores originally need to be a list of N values
        cdef size_t i
        for i in range(N):
            self.scores.append(None)

    def combine (self, Cluster other, dict rmsds):
        """Combine this cluster with another cluster"""
        # Make a list of all of the possible proteins
        cdef list combined = []
        combined.extend(self.proteins)
        combined.extend(other.proteins)
        # Deactivate the other cluster
        other.active = False
        # Select the model protein
        cdef Protein model = select_model_protein(combined)
        # Calculate the RMSD score for the cluster
        self.rmsd = score_cluster(model, rmsds, combined)
        # Update the list of proteins in the cluster
        self.proteins = [model]
        for protein in combined:
            if protein.name != model.name:
                self.proteins.append(protein)
        # Set every score this cluster knows to None
        cdef size_t i, j
        j = len(self.scores)
        for i in range(j):
            self.scores[i] = None

    def sort (self, rmsds):
        """Order the Cluster's members by RMSD"""
        # Only do this if there are at least 3 members
        if len(self.proteins) >= 3:
            # Get the proper set of rmsd values
            use = rmsds[self.proteins[0].name]
            # Do the sorting
            done = False
            while not done:
                done = True
                for i in range(1, len(self.proteins) - 1):
                    if use[self.proteins[i].name] > use[self.proteins[i+1].name]:
                        done = False
                        temp = self.proteins[i]
                        self.proteins[i] = self.proteins[i+1]
                        self.proteins[i+1] = temp

    def output (self, rmsds):
        """Summarize the Cluster's information"""
        # Sort the cluster
        self.sort(rmsds)
        # Store the summary here
        result = "Cluster " + str(self.index) + "\n"
        result += "Model: " + self.proteins[0].name + "\n"
        result += "Members: " + str(len(self.proteins)) + "\n"
        result += "Maximum RMSD: " + format(self.rmsd, '.3f') + "\n"
        result += "Structures\n"
        # The RMSD values to use
        use = rmsds[self.proteins[0].name]
        # Go through the proteins
        for protein in self.proteins:
            result += protein.name + " " 
            result += format(use[protein.name], '.3f').rjust(6)
            result += " " + str(protein) + "\n"
        return result

    def __len__ (self):
        """The number of members in the cluster"""
        return len(self.proteins)

cdef list initialize_clusters (proteins):
    """Initialize the Cluster objects"""
    # They are stored here
    cdef list clusters = []
    # Looping variables
    cdef size_t I, J
    J = len(proteins)
    # Make the clusters
    for I in range(J):
        clusters.append(Cluster(proteins[I], I, J))
    return clusters

def compare_clusters (Cluster cluster1, Cluster cluster2, rmsds):
    """Compare 2 clusters to see how well they combine"""
    # If the score for combining the two clusters is not already known
    if cluster1.scores[cluster2.index] == None:
        # make a list of proteins
        proteins = []
        proteins.extend(cluster1.proteins)
        proteins.extend(cluster2.proteins)
        # Identify the model structure
        model = select_model_protein(proteins)
        # Calculate the score
        score = score_cluster(model, rmsds, proteins)
        # Store that information
        cluster1.scores[cluster2.index] = score
        cluster2.scores[cluster1.index] = score
    return cluster1.scores[cluster2.index]

def identify_best_combination (clusters, rmsds):
    """Identify the best pair of clusters to combine"""
    # Store the best values here
    bestPair = None
    bestScore = None
    # Go through all pairs of clusters
    for i in range(len(clusters) - 1):
        # If the cluster is not active, skip it
        if not clusters[i].active:
            continue
        # Go through all subsequent clusters
        for j in range(i+1, len(clusters)):
            # If the cluster isn't active, skip it
            if not clusters[j].active:
                continue
            # get the score for comparing the clusters
            score = compare_clusters(clusters[i], clusters[j], rmsds)
            # Possibly store that information
            if bestScore == None or score < bestScore:
                bestPair = [i, j]
                bestScore = score
    return bestPair, bestScore

def combine_clusters (clusters, bestPair, rmsds):
    """Combine the best cluster combinations"""
    clusters[bestPair[0]].combine(clusters[bestPair[1]], rmsds)
    for i in range(len(clusters)):
        clusters[i].scores[bestPair[0]] = None
        clusters[i].scores[bestPair[1]] = None

def Clustering_Iteration (clusters, rmsds, MAX):
    """Do an iteration of clustering binding loop structures"""
    # Identify the best pair to combine
    bestPair, bestScore = identify_best_combination(clusters, rmsds)
    if bestScore == None or bestScore > MAX:
        return False
    else:
        combine_clusters(clusters, bestPair, rmsds)
        return True

def output_clusters(clusters, fileName, rmsds):
    """Output the clusters to the specified file"""
    # Order the clusters by number of members that they have
    order = []
    for i in range(len(clusters)):
        if clusters[i].active:
            order.append(i)
    done = False
    while not done:
        done = True
        for i in range(len(order) - 1):
            if len(clusters[order[i]]) < len(clusters[order[i+1]]):
                done = False
                temp = order[i]
                order[i] = order[i+1]
                order[i+1] = temp
    # Open the file
    f = open(fileName, "w")
    # Go through the active clusters
    cdef size_t I
    for i in range(len(order)):
        I = i+1
        clusters[order[i]].index = I
        f.write(clusters[order[i]].output(rmsds) + "\n")
    f.close()
    return len(order)

def Clustering (instructions):
    """Cluster the binding loop structures together"""
    # Open the summary file for appending information
    summary = open(summary_file_name(instructions), "a")
    # Write a message
    message = "Clustering the binding loop structures started on "
    message += time_stamp()
    summary.write(message)
    summary.flush()
    # Create the binding loops
    loops = []
    for data in instructions["loops"]:
        loops.append(Loop(data[0], instructions))
    # Get the maximum RMSD value permitted
    MAX = instructions["RMSD"]
    # Go through the binding loops
    for loop in loops:
        # Load the proteins
        proteins, lengths = loop.load_proteins()
        # Go through the lengths
        for L in lengths:
            # Calculate the rmsd values for these proteins
            rmsds = calculate_RMSDs(proteins[L])
            # Initialize the clusters
            clusters = initialize_clusters(proteins[L])
            # Do the clustering
            active = True
            while active:
                active = Clustering_Iteration(clusters, rmsds, MAX)
            # Output the identified clusters
            fileName = loop.folder + loop.name + "_" + str(L) + "_Clusters.txt"
            N = output_clusters(clusters, fileName, rmsds)
            # Include a message in the summary file
            message = loop.name + " structures of " + str(L) + " residues: "
            message += str(len(clusters)) + " structures were grouped into "
            message += str(N) + " clusters\n"
            summary.write(message)
            summary.flush()
    # Update the summary file that the step is finished
    message = "Clustering the binding loop structures ended on "
    message += time_stamp() + "\n"
    summary.write(message)
    summary.close()
