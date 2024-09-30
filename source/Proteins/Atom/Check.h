/* Written in 2018 by Dr. Robert Pantazes. Copyright owned by Auburn
 * University.
 * This file contains the implementation of methods that check to see whether
 * or not pieces of information match required PDB formatting standards. */

// Use a header guard so that this file is only included in a compiled program
// a single time
#ifndef Check_Guard
#define Check_Guard 1

// This file is intended to be included in the Proteins.h header file directly
// during compilation. Do not include any other header files here.

// Implement the methods of the CHECK namespace
// Atom numbers must be between 1 and 99999
void CHECK::atom_number (const long n) {
    if ((n < 1) || (n > 99999)) {
        stringstream c; c << n;
        string error = "PDB atom numbers must be between 1 and 99,999.\n"
                     + c.str() + " is not acceptable.\n";
        throw logic_error (error);}
}

// Atom names must contain 1-4 characters, each of which is a capital letter
// or digit
void CHECK::atom_name (const string& name) {
    if ((name.size() < 1) || (name.size() > 4)) {
        string error = "PDB atom names must contain between 1 and 4 "
                       "characters.\n" + name + " is not acceptable.\n";
        throw logic_error (error);}
    for (size_t i=0; i<name.size(); ++i) {
        // Atom names can only contain certain characters. Convert the
        // character to a integer and compare it to the proper values
        int j = (int) name[i];
        if (((j < 33) || ((j >= 58) && (j <= 64))) || (j > 90)) {
            string error = "This string contains unacceptable characters for an "
                           "Atom name: " + name + "\n";
            throw logic_error (error);}}
}

// The alternative location information for a residue may be a blank space or
// a capital letter
void CHECK::alt_location (const char L) {
    if ((L != 32) && (!Text::is_upper(L))) {
        string error = "PDB atom alternative location values may only be "
                       "blank spaces or capital letters.\n'";
        error += L;
        error += "' is not acceptable.\n";
        throw logic_error (error);}
}

// Residue names must be exactly 3 capital letters or digits
void CHECK::residue_name (const string& name) {
    // I had to update this because apparently DNA residues don't have to
    // follow the exactly 3 characters requirement
    if ((name.size() < 1) || (name.size() > 3)){
        string error = "PDB residue names must be exactly 3 characters in "
                       "length.\n" + name + " is not acceptable.\n";
        throw logic_error (error);}
    for(size_t i=0; i<name.size(); ++i) {
        if ((!Text::is_digit(name[i])) && (!Text::is_upper(name[i]))) {
            string error = "PDB residue names may only contain digits and "
                           "capital letters.\n"+name+" is not acceptable.\n";
            throw logic_error (error);}}
}

// Protein names may be blank, a digit, or a capital letter
void CHECK::protein_name (const char L) {
    if (((L != 32) && (!Text::is_upper(L))) && (!Text::is_digit(L))) {
        string error = "PDB protein names may only be blank spaces, digits, "
                       "or capital letters.\n'";
        error += L;
        error += "' is not acceptable.\n";
        throw logic_error(error);}
}

// Residue numbers must be between -999 and 9999
void CHECK::residue_number (const long n) {
    if ((n < -999) || (n > 9999)) {
        stringstream c; c << n;
        string error = "PDB residue numbers must be betweeen -999 and 9999.\n"
                     + c.str() + " is not acceptable.\n";
        throw logic_error(error);}
}


// Residue insertion codes may be a blank space or a capital letter
void CHECK::insertion_code (const char L) {
    if ((L != 32) && (!Text::is_upper(L))) {
        string error = "PDB residue insertion codes may only be blank spaces "
                       "or capital letters.\n'";
        error += L;
        error += "' is not acceptable.\n";
        throw logic_error (error);}
}

// Atom coordinates must be between -999.999 and 999.999
void CHECK::atom_coordinate (const PROT::coor n) {
    if ((n <= -999.9995) || (n >= 999.9995)) {
        stringstream c; c << fixed << setprecision(3) << n;
        string error = "PDB atom coordinates must be between -999.999 and "
                       "999.999.\n" + c.str() + " is not acceptable.\n";
        throw logic_error (error);}
}

// Atom occupancies must be between 0.0 and 1.0
void CHECK::occupancy (const float n) {
    if ((n < 0.0) || (n > 1.0)) {
        stringstream c; c << fixed << setprecision(2) << n;
        string error = "PDB atom occupancies must be between 0.0 and 1.0.\n"
                     + c.str() + " is not acceptable.\n";
        throw logic_error (error);}
}

// An Atom's temperature must be between 0 and 1000
void CHECK::temperature (const float n) {
    if ((n < 0.0) || (n > 1000.0)) {
        stringstream c; c << fixed << setprecision(2) << n;
        string error = "PDB atom temperatures must be between 0 and 1000.\n"
                     + c.str() + " is not acceptable.\n";
        throw logic_error (error);}
}

// The Atom's element information may be 0-2 capital letters
void CHECK::element (const string& input) {
    if (input.size() == 0) {return;}
    else if (input.size() > 2) {
        string error = "PDB atom elements may be no more than 2 characters.\n"
                     + input + " is not acceptable.\n";
        throw logic_error (error);}
    for(size_t i=0; i<input.size(); ++i) {
        if (!Text::is_upper(input[i])) {
            string error = "PDB atom elements may only contain capital "
                           "letters.\n" + input + " is not acceptable.\n";
            throw logic_error (error);}}
}

// The Atom's charge may be blank or 2 characters. If it is 2 characters, the
// first must be + or - and the second a digit
void CHECK::charge (const string& input) {
    if (input.size() == 0) {return;}
    else if (input.size() != 2) {
        string error = "PDB atom charge values must be exactly 2 "
                       "characters long.\n" + input + " is not acceptable.\n";
        throw logic_error (error);}
    else if ((input[0] != '-') && (input[0] != '+')) {
        string error = "PDB atom charge values must start with either '+' or "
                       "'-'.\n" + input + " is not acceptable.\n";
        throw logic_error (error);}
    else if (!Text::is_digit(input[1])) {
        string error = "PDB atom charge values must be integer numbers.\n"
                     + input + " is not acceptable.\n";
        throw logic_error (error);}
}

// A function that determines whether or not a Residue is an amino acid
bool CHECK::is_amino_acid (const string& text) {
    // Loop through the 3 amino acid codes
    for(size_t i=0; i<PROT::AA3.size(); ++i) {
        if (text == PROT::AA3[i]) {return true;}}
    // Otherwise, return false
    return false;
}

// Determine whether or not a string represents a backbone atom
bool CHECK::is_backbone_atom (const string& text) {
    // Loop through the backbone atoms
    for(size_t i=0; i<PROT::BackboneAtoms.size(); ++i) {
        // If the text matches the atom, it is a backbone atom
        if (text == PROT::BackboneAtoms[i]) {return true;}}
    return false;
}

// End the header guard from the start of the file
#endif
