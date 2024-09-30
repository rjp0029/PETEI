/* Written in 2018 by Dr. Robert Pantazes. Copyright owned by Auburn
 * University.
 * This file contains the declarations of functions and classes associated
 * with using Protein Data Bank - formatted information. */

// Use a header guard to prevent this file from being included in a compiled
// program multiple times
#ifndef Proteins_Guard
#define Proteins_Guard 1

// Include the "Text.h" header file
#include "Text.h"
// Include standard C++ files
#include <exception>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <cmath>

// Declare the namespace for the classes
namespace PROT {
    // Define some constant information

    // PDB coordinates have 3 decimal places of significance, so a float is
    // appropriate to store that information. Use a typedef to allow that to
    // change later if needed
    typedef float coor;

    // The names of the 20 standard amino acids and 2 variants of histidine
    // that are commonly used
    const string AANames = "ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU "
                           "MET ASN PRO GLN ARG SER THR VAL TRP TYR "
                           "HSD HSE";
    // Split them into a vector
    const vector<string> AA3 = Text::split(AANames);
    // The corresponding 1 letter codes
    const string AACodes = "A C D E F G H I K L M N P Q R S T V W Y";
    const vector<string> AA1 = Text::split(AACodes);
    // The atoms that appear in the backbones of PDB-formatted amino acid
    // residues
    const string bbAtoms = "N HN CA HA C O HN1 HN2 HT1 HT2 HT3 OT1 OT2";
    const vector<string> BackboneAtoms = Text::split(bbAtoms);
    // The number of coordinates an Atom has
    const size_t AtomCoordinates = 3;
    // The number of characters anticipated in a PDB Atom line
    const size_t AtomStringLength = 81;
    // The dielectric constant of water
    const float CCELEC = 331.843;

    // The information about a single Atom
    class Atom;
    // Corresponding classes for amino acids and proteins
    class Residue;
    class Protein;

    // The Matrix class is used to simplify some of the linear algebra
    // operations that need to happen (primarily related to rotating, moving,
    // and calculating the RMSD of protein structures).
    class Matrix;
    // With pointers to 4 atoms, calculate the dihedral angle between them
    coor calculate_dihedral (const Atom *, const Atom *, const Atom *, 
                             const Atom *);

    // PDB files can contain multiple copies of the same protein. The
    // information about them is compiled in this class:
    class Structure;
    // A class for all of the information in a PDB file
    class PDB;

    // End the namespace
}

// A namespace for functions that validate that information follows PDB
// formatting specifications
namespace CHECK {
    // These functions will all throw errors if there is any problem
    void atom_number (const long);
    void atom_name (const string&);
    void alt_location (const char);
    void residue_name (const string&);
    void protein_name (const char);
    void residue_number (const long);
    void insertion_code (const char);
    void atom_coordinate (const PROT::coor);
    void occupancy (const float);
    void temperature (const float);
    void element (const string&);
    void charge (const string&);
    // Check to see whether or not a string is an amino acid name
    bool is_amino_acid (const string&);
    // Check to see whether or not a string is a backbone atom
    bool is_backbone_atom (const string&);
    // End the namespace
}

// Define the Matrix class
class PROT::Matrix {

    // The private variables of the class
    private:
        // The rows and columns of the matrix
        size_t m_rows, m_columns;
        // The values in the matrix (coor is a typedef in the PROT namespace)
        coor * m_values;

    // Private functions to control certain class behaviors
    private:
        // Assign initial attribute values
        void initialize ();
        // Delete dynamically allocated memory
        void clean_up ();
        // Copy information from another matrix
        void copy (const Matrix *);
        // Error check row and column indices
        void error_check (const size_t, const size_t) const;
        // Convert row, column pairs into a linear index
        size_t index (const size_t, const size_t) const;

    // The public interface of the class
    public:
        // The class destructor
        ~Matrix () {clean_up();}
        // The default constructor
        Matrix () {initialize();}
        // Have "allocation" functions for all the non-default ways that a
        // Matrix can be constructed. First is by dimensions
        void allocate (const size_t, const size_t);
        // Next is from an Atom
        void allocate (const Atom *);
        // Set up a rotation matrix using Rodriguez's Rotation Formula
        void allocate (const coor, const coor []);
        // Set up a rotation matrix using a Residue pointer
        void allocate (const Residue *);
        // Set up a Matrix from a Protein
        void allocate (const Protein *, const bool);
        // Set up a Matrix using another Matrix
        void allocate (const Matrix * other) {copy(other);}
        // Set up a Matrix using a vector of Atoms
        void allocate (const vector<Atom *>&);
        // Now have constructors that use those allocation methods
        Matrix (const size_t i, const size_t j) {initialize(); allocate(i, j);}
        Matrix (const Atom * atom) {initialize(); allocate(atom);}
        Matrix (const coor a, const coor v []) {initialize(); allocate(a, v);}
        Matrix (const Residue * res) {initialize(); allocate(res);}
        Matrix (const Protein *, const bool);
        Matrix (const vector<Atom *>& atoms) {initialize(); allocate(atoms);}
        // Copy construction of a Matrix
        Matrix (const Matrix& other) {initialize(); copy(&other);}
        // Copy assignment
        void operator= (const Matrix& other) {copy(&other);}
        // Access to the information in the matrix
        size_t rows () const;
        size_t columns () const;
        coor operator() (const size_t, const size_t) const;
        // Set a matrix value
        void set (const size_t, const size_t, const coor);
        // Modify a value in the matrix
        void increase (const size_t, const size_t, const coor);
        // Multiply two matrices together
        Matrix product (const Matrix *) const;
        // Transpose a matrix
        Matrix transpose () const;
        // Confirm that a matrix can be used to move atoms
        void move_check () const;
        // Or rotate atoms
        void rotate_check () const;
        // A string representation of the Matrix
        string str () const;

        // An additional set of functions that are used in calculating
        // dihedral angles.
        // Convert a vector into a unit vector
        void make_unit_vector ();
        // Subtract one matrix from another
        Matrix operator- (const Matrix&) const;
        // It doesn't get used, but also have addition
        Matrix operator+ (const Matrix&) const;
        // Calculate the cross product between two VECTOR matrices
        Matrix cross (const Matrix&) const;

        // end the class definition
};

// Define the Atom class
class PROT::Atom {

    // The Residue class is a container of Atoms and needs the ability to
    // directly modify the private variables of the Atom (such as the protein
    // name and residue information)
    friend class PROT::Residue;

    // Because the Protein class modifies Residues, it helps that it be a
    // friend of the Atom class, too
    friend class PROT::Protein;
    
    // Let the PDB class have access, too
    friend class PROT::PDB;

    // The Matrix class can be initialized directly from Atom coordinates, so
    // it is appropriate for it to be able to see the coordinates directly
    friend class PROT::Matrix;

    // The private data of the class
    private:
        // The Atom's coordinates (coor is a typedef, AtomCoordinates is a
        // const int)
        coor m_coors [AtomCoordinates];
        // ATOM or HETATM at the start of the PDB line
        string m_type;
        // The atom's name
        string m_name;
        // The name of the atom's residue
        string m_residue;
        // The element of the atom
        string m_element;
        // The atom's charge (often missing, so stored as a string)
        string m_charge;
        // The atom's number
        long m_number;
        // The residue's number
        long m_residue_number;
        // The occupancy
        float m_occupancy;
        // Temperature factor
        float m_temperature;
        // The 'alternative location' specifier
        char m_alt;
        // The residue's insertion code
        char m_insertion;
        // The protein's name
        char m_protein;

        // The previous attributes include everything that is part of a PDB
        // file. These next attributes are associated with calculating
        // specific energy terms
        // VDW energies
        bool vdw; // Whether or not to calculate VDW energies
        float vdwR, vdwE, vdw14R, vdw14E, vdwPhi;
        // Electrostatics
        bool elec;
        float elec_charge;
        // Generalized-borne implicit solvation
        bool gb;
        float gbRadius;
        // Lazaridis-Karplus implicit solvation
        bool lk;
        float lkR, lkL, lkV, lkG;
        // Whether or not the atom is in a buried residue or a surface residue
        bool m_core;

    // Private functions for controlling class behaviour
    private:
        // Assign default values to the class variables
        void initialize ();
        // Copy information from another instance of the class
        void copy (const Atom&);
        // move or rotate the Atom WITHOUT checking the matrix
        void private_move (const Matrix *, const char);
        void private_rotate (const Matrix *);
        // Private functions used in calculating the energy between 2 Atoms
        int exclusions (const Atom *) const;
        float VDW (const Atom *, const coor, const int) const;
        float ELEC (const Atom *, const coor, const int) const;
        float GB (const Atom *, const coor, const int) const;
        float LK (const Atom *, const coor, const int) const;

    // The public interface of the class
    public:
        // The default constructor
        Atom () {initialize();}
        // The standard constructor works from a string of text
        Atom (const string&);
        // Copy construction and assignment
        Atom (const Atom& other) {copy(other);}
        void operator= (const Atom& other) {copy(other);}
        // Access to Atom information
        size_t size () const {return AtomCoordinates;}
        string type () const {return m_type;}
        string name () const {return m_name;}
        string residue () const {return m_residue;}
        string element () const {return m_element;}
        string charge () const {return m_charge;}
        long number () const {return m_number;}
        long residue_number () const {return m_residue_number;}
        float occupancy () const {return m_occupancy;}
        float temperature () const {return m_temperature;}
        char alternative_location () const {return m_alt;}
        char insertion_code () const {return m_insertion;}
        char protein () const {return m_protein;}
        // Access to the Atom's coordinates
        coor operator[] (const size_t) const;
        coor operator[] (const char) const;
        // Generate a PDB formatted string of the Atom's data
        string str () const;
        // Move the Atom, checking the matrix for errors first
        void move (const Matrix *, const char);
        void rotate (const Matrix *);

        // Calculate the distance between this Atom and another
        PROT::coor calculate_distance (const Atom *, const bool) const;
        // Calculate the pairwise additive energy between this Atom and
        // another
        float calculate_energy (const Atom *) const;

        // Provide access to set the VDW calculation parameters
        void set_calculate_vdw (bool how) {vdw = how;}
        void set_vdw_radius (float n) {vdwR = n;}
        // VDW radius information is sometimes needed. Access it with this
        // function
        float get_vdw_radius () const {return vdwR;}
        void set_vdw_epsilon (float n) {vdwE = n;}
        void set_vdw_1_4_radius (float n) {vdw14R = n;}
        void set_vdw_1_4_epsilon (float n) {vdw14E = n;}
        void set_vdw_softening (float n) {vdwPhi = n;}
        // The electrostatics terms
        void set_calculate_elec (bool how) {elec = how;}
        void set_elec_charge (float n) {elec_charge = n;}
        float get_elec_charge () const {return elec_charge;}
        // Generalized-Born implicit solvation values
        void set_calculate_gb (bool how) {gb = how;}
        void set_gb_radius (float n) {gbRadius = n;}
        // Lazaridis-Karplus implicit solvation values
        void set_calculate_lk (bool how) {lk = how;}
        void set_lk_radius (float n) {lkR = n;}
        void set_lk_lambda (float n) {lkL = n;}
        void set_lk_volume (float n) {lkV = n;}
        void set_lk_gibbs (float n) {lkG = n;}

        // Check to see whether or not the Atom's name is a backbone atom
        bool is_backbone_atom () const {return CHECK::is_backbone_atom(m_name);}
        // And whether or not it is a hydrogen
        bool is_hydrogen () const {return (m_name[0] == 'H');}

        // Calculate whether or not two Atoms are closer than 80% of the sum
        // of their VDW radii
        bool check_vdw_overlap (const Atom *) const;
        // Count the number of VDW overlaps between two Atoms (0 or 1)
        size_t calculate_vdw_overlap (const Atom *) const;
    
    // End the class definition
};

// Define the Residue class
class PROT::Residue {

    // The Protein class is a container of Residues. It needs direct access to
    // the contents of a Residue
    friend class PROT::Protein;
    // A Matrix can be created directly from a Residue, so it should be able
    // to access its Atoms directly
    friend class PROT::Matrix;
    // Because PDB files construct Residues directly, they are a friend of the
    // class
    friend class PROT::PDB;

    // The information stored in the class is private
    private:
        // The Atoms in the Residue
        PROT::Atom * m_atoms;
        // The number of them
        size_t m_count;
        // The residue's name
        string m_name;
        // The residue's number
        long m_number;
        // The residue's internal number in a protein or list of residues
        long m_internal;
        // The residue's insertion code (used when two distinct residues have
        // the same number)
        char m_insertion;
        // The residue's protein
        char m_protein;

        // The following attributes are only meaningful in the context of PDB
        // files
        // Whether or not the Residue is present in the PDB file
        bool m_present;
        // The number of Atoms in the Residue that are absent from a PDB file
        size_t m_missing_atoms;

        // Dihedral angles for amino acids in proteins
        coor m_phi;
        coor m_psi;
        coor m_omega;
        
    // Private functions that control the behaviour of the Residue
    private:
        // Assign default values to the Residue's attributes
        void initialize ();
        // Clean up dynamically allocated memory
        void clean_up ();
        // Copy information from another Residue
        void copy (const Residue *);
        // Private functions to set residue information
        void p_set_number (const long, const char, const bool);
        void p_set_protein (const char);
        // Check that a list of Atoms are acceptable for use in a Residue
        void check_atoms (const vector<PROT::Atom *>&) const;
        // The private rotate and move functions of a Residue that do not
        // error check the matrix
        void p_move (const Matrix *, const char);
        void p_rotate (const Matrix *);

    // The public interface of the Residue class
    public:
        // The class destructor
        ~Residue () {clean_up();}
        // The default class constructor
        Residue () {initialize();}
        // Update or entirely create the Residue from a vector of Atoms
        void load (const vector<PROT::Atom *>&, const bool, const bool);
        // Construct a Residue from a vector of Atoms
        Residue (const vector<PROT::Atom *>&);
        // A Residue can be constructed from a string (its name) and a
        // character (the protein's name). This should only be used in PDB
        // files.
        Residue (const string&, const char);
        // Copy construction of a Residue
        Residue (const Residue& other) {initialize(); copy(&other);}
        Residue (const Residue * other) {initialize(); copy(other);}
        // Copy assignment
        void operator= (const Residue& other) {copy(&other);}
        void operator= (const Residue * other) {copy(other);}
        // Access to the Residue's information
        size_t size () const {return m_count;}
        string name () const {return m_name;}
        long number () const {return m_number;}
        long internal_number () const {return m_internal;}
        char insertion_code () const {return m_insertion;}
        char protein () const {return m_protein;}
        // Access to the Residue's Atoms by number or name
        PROT::Atom * operator[] (const size_t);
        PROT::Atom * operator[] (const string&);
        // Functions that allow for the setting of a Residue's number and
        // protein information
        void set_number (const long, const char, const bool);
        void set_protein (const char);
        // Renumber the atoms in the Residue
        long renumber_atoms (long);
        // The number of the last Atom in the Residue
        long last_atom_number () const;
        // A string of all of the Residue's information
        string str (const bool);
        // Move the Residue, with error checking of the provided matrix
        void move (const Matrix *, const char);
        void move (const Matrix *, const bool);
        void rotate (const Matrix *);
        // Move a Residue so that its center of mass is at the origin
        PROT::Matrix center ();
        // Position a Residue for rotamer packaging
        void position (PROT::Matrix&, const bool);

        // Whether or not the Residue is present in the structural data
        bool is_present () const {return m_present;}
        // The number of Atoms the Residue is missing
        size_t missing_atoms () const {return m_missing_atoms;}
        // Calculate a "completeness" score between 0 and 1 for the Residue
        double score () const;
        // Whether or not the Residue is an amino acid
        bool is_amino_acid () const {return CHECK::is_amino_acid(m_name);}

        // Access to the Residue's dihedral angles
        coor phi () const;
        coor psi () const;
        coor omega () const;
        // A function that corrects histidine residue names
        void charmm_his_fix ();

        // This function creates a duplicated copy of the Residue and all of
        // its atoms.
        PROT::Residue duplicate () const;

        // Calculate the VDW overlaps between the Residue's Atoms and another
        // Atom
        size_t calculate_vdw_overlap (const Atom *, const bool) const;
        // Do the same with a Residue
        size_t calculate_vdw_overlap (const Residue *, const bool, 
                                      const bool) const;
        // Calculate the energy between the Residue and an Atom
        float calculate_energy (const Atom *, const bool) const;
        // Calculate the energy between the Residue and another Residue
        float calculate_energy (const Residue *, const bool, const bool) const;

    // End the Residue class definition
};

// Declare the Protein class
class PROT::Protein {
    
    // A Matrix can be assembled directly from a Protein, so it should have
    // internal access to the Protein's information
    friend class PROT::Matrix;

    // The PDB class is a friend, too
    friend class PROT::PDB;

    // The information stored in the Protein is private
    private:
        // The Protein's name
        char m_name;
        // The Residues the protein holds
        size_t m_count;
        PROT::Residue * m_residues;

    // Private functions that control behaviour of the Protein
    private:
        // Assign default initial values to the variables
        void initialize ();
        // Delete dynamically allocated memory
        void clean_up ();
        // Copy the contents of another Protein into this one
        void copy (const Protein *);
        // Make sure that a set of Residues are acceptable for use in the
        // Protein
        void check_residues (vector<PROT::Residue>&) const;
        // Private functions for moving and rotating the protein
        void p_move (const Matrix *, const char);
        void p_rotate (const Matrix *);
        // Create a vector of pointers to the backbone atoms or all atoms in
        // the protein
        vector<Atom *> select_atoms (const bool) const;

    // The public interface of the Protein class
    public:
        // The class destructor
        ~Protein () {clean_up();}
        // A default constructor
        Protein () {initialize ();}
        // A Protein can be loaded from:
        // 1) a file, with an optional path
        void load (const string&, const string);
        // 2) a vector of strings
        void load (const vector<string>&);
        // 3) a vector of Atoms
        void load (vector<Atom>&);
        // 4) a vector of Residues
        void load (vector<Residue>&);
        // Those same inputs can be used for constructing a Protein
        Protein (const string&, const string);
        Protein (const vector<string>&);
        Protein (vector<Atom>&);
        Protein (vector<Residue>&);
        // Copy construction of a Protein
        Protein (const Protein& other) {initialize(); copy(&other);}
        Protein (const Protein * other) {initialize(); copy(other);}
        // Copy assignment
        void operator= (const Protein& other) {copy(&other);}
        void operator= (const Protein * other) {copy(other);}
        // A formatted string of text listing all of the information in the
        // Protein
        string str (const bool);
        // The number of Atoms in the Protein
        size_t number_of_atoms () const;
        // The number of the last atom in the protein
        long last_atom_number () const;
        // The INTERNAL number of the last residue in the protein
        long last_residue_number () const;
        // The number of residues in the protein
        size_t size () const {return m_count;}
        // The name of the protein
        char name () const {return m_name;}
        // Change the internal numbering of the residues
        long renumber_residues (long);
        // Renumber the atoms in the Protein
        long renumber_atoms (long);
        // Change the Protein's name
        void set_name (const char);
        // Access to the Residues in the Protein
        Residue * operator() (const long, const char, const bool) const;
        // Move and rotate the protein
        void move (const Matrix *, const char);
        void move (const Matrix *, const bool);
        void rotate (const Matrix *);
        // Move a Protein so that its center of mass is at the origin
        PROT::Matrix center (const bool);

        // Calculate a completeness score between 0 and 1 for the Protein.
        // This is used only in the context of PDB files
        double score () const;

        // CHARMM has multiple Histidine residues because of its variable
        // protonation states. These functions convert HIS into our
        // CHARMM-preferred variant (HSD) and back again
        void charmm_his_prep ();
        void charmm_his_fix ();

        // Create a string listing the Protein's sequence
        string get_sequence () const;

        // Calculate the Protein's Residues' dihedral angles
        void calculate_dihedrals (const bool);
        // Calculate the RMSD between two proteins
        coor calculate_RMSD (const Protein *, const bool) const;

        // Calculate the VDW overlaps between an Atom and the Protein's Atoms
        size_t calculate_vdw_overlap (const Atom *) const;
        // Between a Residue and the Protein's Atoms
        size_t calculate_vdw_overlap (const Residue *, const bool) const;
        // Between two Proteins
        size_t calculate_vdw_overlap (const Protein *) const;
        // Calculate energies between a Protein and:
        float calculate_energy (const Atom *) const;
        float calculate_energy (const Residue *, const bool) const;
        float calculate_energy (const Protein *) const;

        // Create a duplicated copy of a protein
        PROT::Protein duplicate () const;

    // End the declaration of the class
};

// The Structure class is used exclusively in the PDB class. It is a container
// of information about Proteins with different chain names but that are
// actually the same molecule
class PROT::Structure {

    // The PDB class needs access to internal information in the Structure
    friend class PROT::PDB;

    // The information stored in the class is private
    private:
        // Pointers to the Proteins that make up the Structure
        vector<Protein *> m_proteins;
        // The names of the structure
        vector<string> m_names;

    // Private functions that control the behaviour of the structure
    private:
        // Copy information from another instance of the class
        void copy (const Structure *);
        void copy (const Structure *, vector<Protein>&);

    // The public interface of the class
    public:
        // The default constructor
        Structure () {m_proteins.clear(); m_names.clear();}
        // Copy construction and assignment
        Structure (const Structure& other) {copy(&other);}
        Structure (const Structure * other) {copy(other);}
        void operator= (const Structure& other) {copy(&other);}
        void operator= (const Structure * other) {copy(other);}
        // Copy construction with a different, but matching, set of Proteins
        // (used in copying PDB files)
        Structure (const Structure&, vector<Protein>&);
        // The number of proteins that are part of the structure
        size_t proteins () const {return m_proteins.size();}
        // Access to a pointer
        PROT::Protein * protein (const size_t);
        // The names of the Structure
        size_t names () const {return m_names.size();}
        string name (const size_t) const;
        // Store a protein in the structure
        void store_protein (Protein *);
        // Store a name in the structure
        void store_name (const string&);
        // Create a string representation summarizing the Structure's
        // information
        string str () const;

    // End the class definition
};

// The PDB class organizes the information from a PDB file into a useable form
class PROT::PDB {

    // The information stored in the class is private
    private:
        // The name of the pdb file
        string m_name;
        // The folder the file is located in
        string m_folder;
        // The contents of that file
        vector<string> m_lines;
        // The type of experiment used to generate the information
        string m_type;
        // The resolution of the experimental data
        double m_resolution;
        // The Proteins in the file
        vector<Protein> m_proteins;
        // The structures they group into
        vector<Structure> m_structures;
        // Whether or not the PDB file is obsolete and shouldn't be used
        bool m_obsolete;

    // Private functions that control behaviour of the class
    private:
        // Copy information from another instance of this class
        void copy (const PDB *);
        // Load the contents of the file
        void load ();
        // Identify the experiment type
        void identify_type();
        // Check to see whether or not the PDB file is obsolete or a
        // theoretical model. In either case, it should not be used
        void check_obsolete ();
        // Identify the experiment's resolution
        void identify_resolution ();
        // Use the SEQRES lines to initialize vectors of Residues that will
        // eventually be turned into Proteins
        void initialize_Residues (vector<vector<Residue> >&);
        // Use the Remark 465 lines to identify missing Residues
        void identify_missing_residues (vector<vector<Residue> >&);
        // Collect the Atoms that provide structural details about the
        // Proteins
        void collect_Atoms (vector<vector<Atom> >&);
        // Integrate the missing Residues into the Residues in a Protein
        void integrate_missing_residues (Protein&, vector<Residue>&);
        // Validate the SEQRES information
        void validate_seqres (Protein *, vector<Residue>&);
        // Use Remark 470 lines to identify atoms that are missing in Residues
        // that are present in the Proteins
        void identify_missing_atoms ();
        // Construct the proteins
        void construct_Proteins ();
        // Create the Structures
        void create_Structures ();

    // The public interface of the PDB class
    public:
        // The class constructor
        PDB (const string&, const string);
        // Copy construction and assignment
        PDB (const PDB& other) {copy(&other);}
        PDB (const PDB * other) {copy(other);}
        void operator= (const PDB& other) {copy(&other);}
        void operator= (const PDB * other) {copy(other);}
        // Access to the class information
        string name () const {return m_name;}
        string folder () const {return m_folder;}
        size_t lines () const {return m_lines.size();}
        string line (const size_t) const;
        string type () const {return m_type;}
        double resolution () const {return m_resolution;}
        size_t proteins () const {return m_proteins.size();}
        Protein * protein (const size_t);
        size_t structures () const {return m_structures.size();}
        Structure * structure (const size_t);
        bool obsolete () const {return m_obsolete;}
        // A string representation of the PDB file's information
        string str () const;

    // End the class definition
};

// Include the header files that implement all of these classes
#include "Atom/Check.h"
#include "Matrix/Matrix.h"
#include "Atom/Atom.h"
#include "Residue/Residue.h"
#include "Protein/Protein.h"
#include "PDB/Structure.h"
#include "PDB/PDB.h"

// Implement the calculate dihedral angle function
PROT::coor PROT::calculate_dihedral (const Atom * atom1, const Atom * atom2,
                                     const Atom * atom3, const Atom * atom4) {
    // Make vectors of the differences between certain atoms
    Matrix f = Matrix(atom1) - Matrix(atom2); f.make_unit_vector();
    Matrix g = Matrix(atom2) - Matrix(atom3); g.make_unit_vector();
    // This is supposed to be 4 - 3
    Matrix h = Matrix(atom4) - Matrix(atom3); h.make_unit_vector();
    // Calculate cross products between the vectors
    Matrix a = f.cross(g); a.make_unit_vector();
    Matrix b = h.cross(g); b.make_unit_vector(); b = b.transpose();
    // Calculate the dot product of a and b
    Matrix c = a.product(&b);
    // Calculate the angle
    coor angle = (180.0/M_PI) * acos(c(0,0));
    // Calculate appropriate values for checking if the sign on the angle
    // needs to be changed
    Matrix check1 = a.cross(b); check1 = check1.transpose();
    Matrix check2 = g.product(&check1);
    if (check2(0,0) > 0) {angle = -angle;}
    return angle;
}


// End the header guard from the start of the file
#endif
