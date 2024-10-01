# PETEI
Proteins Engineered with Targeted Exceptional Interactions (PETEI) is a binding protein design algorithm suitable for designing proteins that bind target molecules using loops, such as antibodies and 10th Type III Fibronectin domains. It was created by the PROTEIN PANT(z) Lab at Auburn University, led by Dr. Robert Pantazes.

# Compiling Cython Code
PETEI requires the use of C++ extension files for Python3, which are compiled with Cython. To do so, use the setup.py script provided in the source folder. From that folder, run the command:
```
python3 setup.py build_ext -i
```
The most recent testing of PETEI by the PROTEIN PANT(z) Lab used Python 3.9.10 and Cython 0.29.25. Cython used gcc version 11.3.0 to compile the extension modules.

# Needed External Resources
PETEI uses parameters from the Chemistry at HARvard Molecular Mechanics (CHARMM) program during parts of its calculations. Information about CHARMM is available here: https://academiccharmm.org/. CHARMM is free for use for academic users, but has its own licensing requirements.
The specific CHARMM files that PETEI needs are top_all36_prot.rtf and par_all36_prot.prm, which are the standard topology and parameter files for the program.

In "Extracellular Peptide-Ligand Dimerization Actuator Receptor Design for Reversible and Spatially Dosed 3D Cell-Material Communication", PETEI used an older version of the Dunbrack Rotamer Library. In this distribution, it has been updated to include rotamers from the 2010 library (DOI: 10.1016/j.str.2011.03.019), which is distributed under the CC BY 4.0 License. PETEI includes the rotamers it needs in the /resources/rotamers/ folder.

# Creating Databases
The distribution of PETEI includes one database, entitled 1TTGu, which can be used to design 10th Type III Fibronectin domains (i.e., monobodies) using the framework of chain A from PDB file 1TTG. This database is an "updated" version of the 1TTG database used in "Extracellular Peptide-Ligand Dimerization Actuator Receptor Design for Reversible and Spatially Dosed 3D Cell-Material Communication." The change is solely due to the update in the rotamers that PETEI uses, and the specific loops in the database match those in the manuscript.

PETEI can be used to design any protein that binds target antigens using modular loops. New databases for design tasks can be created using the Database_Maker.py Python program included with the distribution of PETEI. That program is meant to be run from the main folder of PETEI and will automatically place a new database it creates in the /databases/ folder.

There are two primary things that PETEI needs to create a new database. The first is a PDB-formatted structure of the protein framework for the database that includes all atoms (including hydrogens). Experimentally determined protein structures are often missing some atoms, so they must be added back in. The PROTEIN PANT(z) Lab recommends using CHARMM for this purpose, which is the tool we used.

The second thing PETEI needs is a list of the proteins to search for loops. The information must be listed in a file, and each line of the file must contain two strings: the name of a file to search and the folder in which that file is located. Like the protein framework, the analyzed files must have all atoms present, including hydrogens. 

Additional specific details of how to provide instructions to the Database_Maker.py program are in the /instructions/Database_Instructions.txt file of the distribution.

PETEI uses a default tolerance of 0.225 Angstroms on atom-atom distances for loop placement. If a user wishes to change this, it is in line 236 of /source/BPD/Database/LoopSearch.pyx. PETEI will need to be recompiled each time that parameter is changed.

# Designing Binding Proteins
To design proteins to bind a protein, PETEI needs a structure of the target protein of interest. As with creating the databases, that protein structure must contain all atoms, including hydrogens. 

PETEI designs proteins using the Designer.py program included with the distribution of PETEI. That program is meant to be run from the main folder of PETEI and will automatically place its findings in a /results/ folder. It is likely you will have to initially create this folder. Additional specific information about how to provide instructions to the Designer.py program are in the /instructions/Design_Instructions.txt file of the distribution.

PETEI uses a default tolerance of 0.33 Angstroms on deviations of interaction distances from their ideal geometries. If a user wishes to change this, it is in line 281 of /source/BPD/Design/Prepare.py. Because this is a Python file rather than a Cython file, PETEI does not need to be recompiled after changing this parameter. Note that the numbers of solutions found by PETEI are very sensitive to this threshold and it is not recommended that it be changed.

The three target protein structures used in "Extracellular Peptide-Ligand Dimerization Actuator Receptor Design for Reversible and Spatially Dosed 3D Cell-Material Communication" are included in /structures/. Changes in the available rotamers led to changes in the binding loop database, which means that running this version of PETEI will not exactly duplicate the results of "Extracellular Peptide-Ligand Dimerization Actuator Receptor Design for Reversible and Spatially Dosed 3D Cell-Material Communication". The three monobodies used in studies of the paper were HA3, FLAG5, and MYC7. If this distribution of PETEI is run using the 1TTGu database and the provided target structure, successful solution 323 has 2 / 3 binding loops identical to HA3. Because PETEI requires no more than 1 binding to be the same between solutions, this precludes the exact identification of HA3. FLAG5 exactly corresponds to successful solution 78 when this distribution of PETEI is run using the provided FLAG antigen. Finally, successful solution 138 has 2 / 3 binding loops the same as MYC7 with this distribution of PETEI. 

# Installation and Usage
To download and install PETEI, navigate to your working directory:
```
cd /PATH/TO/WORKING/DIRECTORY/
```
Then, clone the repository: 
```
git clone https://github.com/rjp0029/PETEI.git
```
Next, compile the code:
``` 
cd PETEI/source/ && python3 setup.py build_ext -i
```
Obtain CHARMM topology and paramter files (top_all36_prot.rtf and par_all36_prot.prm) and place them in PETEI/resources/\
Next, make the results folder:
```
cd /PATH/TO/WORKING/DIRECTORY/PETEI/ && mkdir results
```
Finally, create instruction files for your database or design calculations and run the database maker and designer scripts:
```
python3 Database_Maker.py ./instructions/Database_Instructions.txt
python3 Designer.py ./instructions/Design_Instructions.txt
```
