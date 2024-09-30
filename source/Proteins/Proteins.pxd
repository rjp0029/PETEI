# Written in 2018 by Dr. Robert Pantazes. Copyright owned by Auburn
# University.
# This file contains Cython declarations of C++ classes that are used for
# storing and working with PDB-formatted protein information. 

# Cython compiler directives
# cython: c_string_type=str
# cython: c_string_encoding=ascii
# cython: language_level=3
# distutils: language=c++

# Include the C++ string, bool and vector classes
from libcpp.string cimport string
from libcpp cimport bool
from libcpp.vector cimport vector

# Include the contents from the PROT namespace of the Proteins.h header file
cdef extern from "Proteins.h" namespace "PROT":

    # Declare the typedef "coor", which is defined as a float in Proteins.h
    ctypedef float coor

    # Declare the Matrix class
    cdef cppclass Matrix:
        Matrix ()
        void allocate (const size_t, const size_t) except +
        void allocate (const Atom *) except +
        void allocate (const coor, const coor [])
        void allocate (const Residue *) except +
        void allocate (const Protein *, const bool) except +
        void allocate (const Matrix *) except +
        void allocate (const vector[Atom *]&) except +
        Matrix (const size_t, const size_t) except +
        Matrix (const Atom *)
        Matrix (const Residue *) except +
        Matrix (const Protein *, const bool) except +
        Matrix (const coor, const coor [])
        Matrix (const Matrix *)
        Matrix (const vector[Atom *]&) except +
        size_t rows () const
        size_t columns () const
        coor operator() (const size_t, const size_t) except +
        void set (const size_t, const size_t, const coor) except +
        void increase (const size_t, const size_t, const coor) except +
        Matrix product (const Matrix *) except +
        Matrix transpose () const
        void move_check () except +
        void rotate_check () except +
        string str () const

    # Declare the Atom class
    cdef cppclass Atom:
        Atom ()
        Atom (const string&) except +
        Atom (const Atom&)
        size_t size () const
        string type () const
        string name () const
        string residue () const
        string element () const
        string charge () const
        long number () const
        long residue_number () const
        float occupancy () const
        float temperature () const
        char alternative_location () const
        char insertion_code () const
        char protein () const
        coor operator[] (const size_t) except +
        string str () const
        void move (const Matrix *, const char) except +
        void rotate (const Matrix *) except +
        coor calculate_distance (const Atom *, const bool) const
        float calculate_energy (const Atom *) const
        void set_calculate_vdw (bool)
        void set_vdw_radius (float)
        float get_vdw_radius () const
        void set_vdw_epsilon (float)
        void set_vdw_1_4_radius (float)
        void set_vdw_1_4_epsilon (float)
        void set_vdw_softening (float)
        void set_calculate_elec (bool)
        void set_elec_charge (float)
        float get_elec_charge () const
        void set_calculate_gb (bool)
        void set_gb_radius (float)
        void set_calculate_lk (bool)
        void set_lk_radius (float)
        void set_lk_lambda (float)
        void set_lk_volume (float)
        void set_lk_gibbs (float)
        bool is_backbone_atom () const
        bool is_hydrogen () const
        bool check_vdw_overlap (const Atom *) const
        size_t calculate_vdw_overlap (const Atom *) const

    # Declare the Residue class
    cdef cppclass Residue:
        Residue ()
        void load (const vector[Atom *]&, const bool, const bool) except +
        Residue (const vector[Atom *]&) except +
        Residue (const Residue&)
        Residue (const Residue *)
        size_t size () const
        string name () const
        long number () const
        long internal_number () const
        char insertion_code () const
        char protein () const
        Atom * operator[] (const size_t) except +
        Atom * operator[] (const string&) except +
        void set_number (const long, const char, const bool) except +
        void set_protein (const char) except +
        long renumber_atoms (long)
        long last_atom_number () const
        string str (const bool)
        void move (const Matrix *, const bool) except +
        void rotate (const Matrix *) except +
        Matrix center () except +
        void position (Matrix&, const bool) except +
        size_t missing_atoms () const
        bool is_amino_acid () const
        bool is_present () const
        coor phi () except +
        coor psi () except +
        coor omega () except +
        void charmm_his_fix ()
        Residue duplicate () const
        size_t calculate_vdw_overlap (const Atom *, const bool) const
        size_t calculate_vdw_overlap (const Residue *, const bool, const bool) const
        size_t calculate_energy (const Atom *, const bool) const
        size_t calculate_energy (const Residue *, const bool, const bool) const

    # Declare the Protein class
    cdef cppclass Protein:
        Protein ()
        void load (const string&, const string) except +
        void load (const vector[string]&) except +
        void load (vector[Atom]&) except + 
        void load (vector[Residue]&) except +
        Protein (const string&, const string) except +
        Protein (const vector[string]&) except +
        Protein (const vector[Atom]&) except +
        Protein (const vector[Residue]&) except +
        Protein (const Protein&)
        Protein (const Protein *)
        string str (const bool) except +
        size_t number_of_atoms () const
        long last_atom_number () const
        long last_residue_number () const
        size_t size () const 
        char name () const
        long renumber_residues (long)
        long renumber_atoms (long)
        void set_name (const char) except +
        Residue * operator() (const long, const char, const bool) except +
        void move (const Matrix *, const bool)
        void rotate (const Matrix *)
        Matrix center (const bool)
        double score () const
        void charmm_his_prep ()
        void charmm_his_fix ()
        string get_sequence () const
        void calculate_dihedrals (const bool) except +
        coor calculate_RMSD (const Protein *, const bool) except +
        size_t calculate_vdw_overlap (const Atom *) const
        size_t calculate_vdw_overlap (const Residue *, const bool) const
        size_t calculate_vdw_overlap (const Protein *) const
        size_t calculate_energy (const Atom *) const
        size_t calculate_energy (const Residue *, const bool) const
        size_t calculate_energy (const Protein *) const
        Protein duplicate () const

    # Declare the Structure class used in PDB files
    cdef cppclass Structure:
        Structure ()
        size_t proteins () const
        Protein * protein (const size_t) except +
        size_t names () const
        string name (const size_t) except +
        string str () const

    # Declare the PDB class
    cdef cppclass PDB:
        PDB (const string&, const string) except +
        string name () const
        string folder () const
        size_t lines () const
        string line (const size_t) except +
        string type () const
        double resolution () const
        size_t proteins () const
        Protein * protein (const size_t) except +
        size_t structures () const
        Structure * structure (const size_t) except +
        bool obsolete () const
        string str () const
