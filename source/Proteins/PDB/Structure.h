/* Written in 2018 by Dr. Robert Pantazes. Copyright owned by Auburn
 * University.
 * This file contains the definitions of the methods of the Structure class. A
 * structure is a group of one or more "Proteins" from a PDB file that are
 * identical to one another. */

// Use a header guard to prevent this file from being included in a compiled
// program multiple times
#ifndef Structure_Header_Guard
#define Structure_Header_Guard 1

// This file is intended to be included in a compiled program directly by the
// Proteins.h header file. Do not include any other dependencies here.

// Copy the contents from another Structure into this one
void PROT::Structure::copy (const Structure * other) {
    // copy the protein pointers
    m_proteins.clear(); m_proteins.reserve(other->m_proteins.size());
    for(size_t i=0; i<other->m_proteins.size(); ++i) {
        m_proteins.push_back(other->m_proteins[i]);}
    // Copy the names
    m_names.clear(); m_names.reserve(other->m_names.size());
    for(size_t i=0; i<other->m_names.size(); ++i) {
        m_names.push_back(other->m_names[i]);}
    return;
}

// Copy the contents of another Structure, but create pointers to a different
// set of proteins
void PROT::Structure::copy (const Structure * other, vector<Protein>& prots) {
    // Copy the structure's names
    m_names.clear(); m_names.reserve(other->m_names.size());
    for(size_t i=0; i<other->m_names.size(); ++i) {
        m_names.push_back(other->m_names[i]);}
    // Clear the proteins and reserve the proper number of entries
    m_proteins.clear(); m_proteins.reserve(other->m_proteins.size());
    // Loop through the protein options
    for(size_t i=0; i<other->m_proteins.size(); ++i) {
        // Get the name of the other protein
        char L = other->m_proteins[i]->name();
        // Use this flag to indicate whether or not a protein is found
        bool flag = false;
        // Loop through the provided proteins
        for(size_t j=0; j<prots.size(); ++j) {
            if (prots[j].name() == L) {
                flag = true;
                m_proteins.push_back(&(prots[j]));
                break;}}
        // If the protein wasn't identified, throw an error
        if (!flag) {
            string error = "Failure to find Protein ";
            error += L;
            error += " by the Structure copy function.\n";
            throw logic_error (error);}}
    return;
}

// Access to a pointer in a structure
PROT::Protein * PROT::Structure::protein (const size_t i) {
    if (i >= m_proteins.size()) {
        stringstream c1; c1 << i;
        stringstream c2; c2 << m_proteins.size();
        string error = c1.str() + " is not a valid index in a Structure with "
                     + c2.str() + " proteins.\n";
        throw logic_error (error);}
    return m_proteins[i];
}

// Access to a name of the structure
string PROT::Structure::name (const size_t i) const {
    if (i >= m_names.size()) {
        stringstream c1; c1 << i;
        stringstream c2; c2 << m_names.size();
        string error = c1.str() + " is not a valid index in a Structure with "
                     + c2.str() + " names.\n";
        throw logic_error (error);}
    return m_names[i];
}

// The constructor that uses the non-standard copy function
PROT::Structure::Structure (const Structure& other, vector<Protein>& prots) {
    copy(&other, prots);
    return;
}

// Store a protein pointer or name in the structure
void PROT::Structure::store_protein (Protein * ptr) {m_proteins.push_back(ptr);}
void PROT::Structure::store_name (const string& text) {m_names.push_back(text);}

// Create a string summarizing the Structure's information
string PROT::Structure::str () const {
    // Store the output here
    string output = "Molecule\nName(s):";
    // If there are no listed names
    if (m_names.size() == 0) {output += " Unknown\n";}
    // If there are listed names, list them
    else {
        for(size_t i=0; i<m_names.size(); ++i) {
            output += " " + m_names[i];
            if (i < m_names.size() - 1) {output += ";";}
            else {output += "\n";}}}
    // If there are no proteins
    if (m_proteins.size() == 0) {output += "Chains: None\n";}
    // If there are proteins
    else {
        //Loop through them
        for(size_t i=0; i<m_proteins.size(); ++i) {
            // Calculate a score for the protein
            stringstream c1; 
            c1 << fixed << setprecision(4) << m_proteins[i]->score();
            // List the molecule
            output += "Chain ";
            output += m_proteins[i]->name();
            output += ": Score = " + c1.str() + "\n";}}
    // Add an extra end line character
    output += "\n";
    return output;
}

// End the header guard from the start of the file
#endif
