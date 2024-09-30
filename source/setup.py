# A setup script to use Cython and distutils to create Python modules and
# packages for analyzing and designing Protein Structures

from distutils.core import setup, Extension
from Cython.Build import cythonize

# Create the objects for each Extension Module

AtomExt = Extension (name = "Proteins.Atom.PyAtom",
                     sources = ['Proteins/Atom/PyAtom.pyx'],
                     include_dirs = ["Proteins/"],
                     language = "c++")

ResExt = Extension (name = "Proteins.Residue.PyResidue",
                    sources = ["Proteins/Residue/PyResidue.pyx"],
                    include_dirs = ["Proteins/"],
                    language = "c++")

ProtExt = Extension (name = "Proteins.Protein.PyProtein",
                     sources = ["Proteins/Protein/PyProtein.pyx"],
                     include_dirs = ["Proteins/"],
                     language = "c++")

MatrixExt = Extension (name = "Proteins.Matrix.PyMatrix",
                       sources = ["Proteins/Matrix/PyMatrix.pyx"],
                       include_dirs = ["Proteins/"],
                       language = "c++")

StrucExt = Extension (name = "Proteins.PDB.PyStructure",
                      sources = ["Proteins/PDB/PyStructure.pyx"],
                      include_dirs = ["Proteins/"],
                      language = "c++")

PDBExt = Extension (name = "Proteins.PDB.PyPDB",
                    sources = ["Proteins/PDB/PyPDB.pyx"],
                    include_dirs = ["Proteins/"],
                    language = "c++")

BPDGeneralExt = Extension (name = "BPD.General",
                           sources = ["BPD/General.pyx"],
                           include_dirs = ["Proteins/"],
                           language = "c++")

BPDLoopExt = Extension (name = "BPD.Database.LoopSearch",
                        sources = ["BPD/Database/LoopSearch.pyx"],
                        include_dirs = ["Proteins/"],
                        language = "c++")

BPDClusterExt = Extension (name = "BPD.Database.Cluster",
                           sources = ["BPD/Database/Cluster.pyx"],
                           include_dirs = ["Proteins/"],
                           language = "c++")

BPDCompatibleExt = Extension (name = "BPD.Database.Compatible",
                              sources = ["BPD/Database/Compatible.pyx"],
                              include_dirs = ["Proteins/"],
                              language = "c++")

BPDIntExt = Extension (name = "BPD.Database.Interactions",
                       sources = ["BPD/Database/Interactions.pyx"],
                       include_dirs = ["Proteins/"],
                       language = "c++")

BPDDesClsExt = Extension (name = "BPD.Design.PyDesignClasses",
                          sources = ["BPD/Design/PyDesignClasses.pyx"],
                          include_dirs = ["Proteins/", "BPD/Design/"],
                          language = "c++")

BPDDesSrchExt = Extension (name = "BPD.Design.Searching",
                           sources = ["BPD/Design/Searching.pyx"],
                           include_dirs = ["Proteins/", "BPD/Design/"],
                           language = "c++")

BPDAssessExt = Extension (name = "BPD.Design.Assessing",
                          sources = ["BPD/Design/Assessing.pyx"],
                          include_dirs = ["Proteins/", "BPD/Design/"],
                          language = "c++")

setup (name = "PANTZ",
       version = "0.1",
       description = "Protein Analysis and design Toolz",
       author = "Robert Pantazes",
       author_email = "rjp0029@auburn.edu",
       ext_modules = cythonize([AtomExt, ResExt, ProtExt, MatrixExt, StrucExt,
                                PDBExt, BPDGeneralExt, BPDLoopExt,
                                BPDClusterExt, BPDCompatibleExt, BPDIntExt,
                                BPDDesClsExt, BPDDesSrchExt, BPDAssessExt])
       )
