/* Written in 2018 by Dr. Robert Pantazes. Copyright owned by Auburn
 * University.
 * This file contains the implementation of the Protein class - a core
 * construct intended to hold all information about a Protein from a PDB file.
 * */

// Use a header guard to prevent this file from being loaded into a compiled
// program more than once
#ifndef Protein_Header_Guard
#define Protein_Header_Guard 1

// This file is intended to be included in a compiled program directly from
// the Proteins.h header file. Do not include any additional dependencies
// here.

// Assign default values to the Protein's attributes
void PROT::Protein::initialize () {
    m_name = ' ';
    m_count = 0;
    m_residues = 0;
    return;
}

// Delete dynamically allocated memory
void PROT::Protein::clean_up () {
    if (m_residues != 0) {delete[] m_residues; m_residues = 0;}
    return;
}

// Copy the information from another Protein into this one
void PROT::Protein::copy (const Protein * other) {
    // Delete any existing information
    clean_up();
    // Copy the name
    m_name = other->m_name;
    // Copy the number of Residues
    m_count = other->m_count;
    // If there are Residues, allocate memory and copy them
    if (m_count > 0) {
        m_residues = new Residue [m_count];
        for(size_t i=0; i<m_count; ++i) {m_residues[i] = other->m_residues[i];}}
    return;
}

// Determine whether or not a vector of Residues is acceptable for use in a
// Protein
void PROT::Protein::check_residues (vector<PROT::Residue>& residues) const {
    // If there are no residues, raise an error
    if (residues.size() == 0) {
        string text = "A Protein cannot be loaded from an empty vector of "
                      "Residues.\n";
        throw logic_error (text);}
    // If there is only a single residue, be done 
    else if (residues.size() == 1) {return;}
    // Since there are multiple residues, make sure they all have the same
    // protein character
    for(size_t i=1; i<residues.size(); ++i) {
        if (residues[0].m_protein != residues[i].m_protein) {
            string text = "A Protein cannot be loaded from a vector of "
                          "Residues that belong to different proteins.\n";
            throw logic_error (text);}}
    // Check the residues to see which (if any) of them have only HETATMs
    vector<bool> hetatms; hetatms.reserve(residues.size());
    for(size_t i=0; i<residues.size(); ++i) {
        bool het = true;
        for(size_t j=0; j<residues[i].m_count; ++j) {
            if (residues[i].m_atoms[j].type() != "HETATM") {
                het = false; break;}}
        hetatms.push_back(het);}
    // Ensure every Residue has a unique number / insertion code combination
    for(size_t i=0; i<residues.size() - 1; ++i) {
        // Skip this residue if it only contains hetero atoms
        if (hetatms[i]) {continue;}
        // Loop through all subsequent residues
        for (size_t j=i+1; j<residues.size(); ++j) {
            if (hetatms[j]) {continue;}
            // Check the residues' numbering data
            if ((residues[i].m_number == residues[j].m_number) &&
                (residues[i].m_insertion == residues[j].m_insertion)) {
                string error = "Every Residue in a Protein must have a unique "
                               "number.\n";
                error += residues[i].str() + residues[j].str();
                throw logic_error (error);}}}
    return;
}

// Move and rotate the Protein
void PROT::Protein::p_move (const Matrix * mat, const char how) {
    if (m_count > 0) {
        for(size_t i=0; i<m_count; ++i) {m_residues[i].p_move(mat, how);}}
    return;
}

void PROT::Protein::p_rotate (const Matrix * mat) {
    if (m_count > 0) {
        for(size_t i=0; i<m_count; ++i) {m_residues[i].p_rotate(mat);}}
    return;
}

// Create a vector of atom pointers for either every atom in the protein or
// only backbone atoms from amino acids
vector<PROT::Atom *> PROT::Protein::select_atoms (const bool backboneOnly) const {
    // This is the vector
    vector<Atom *> atoms;
    // Since this is a private function, I know it will only be called when
    // appropriate. So we know there are residues to loop through
    for(size_t i=0; i<m_count; ++i) {
        // Determine if the residue is an amino acid
        bool amino = CHECK::is_amino_acid(m_residues[i].m_name);
        // Loop through the Residue's atoms
        for(size_t j=0; j<m_residues[i].m_count; ++j) {
            // get a pointer to the Atom
            Atom * ptr = &(m_residues[i].m_atoms[j]);
            // If this residue is not an amino acid OR backbone only is False,
            // store it
            if ((!amino) || (!backboneOnly)) {
                atoms.push_back(ptr);}
            // If it is an amino acid and only backbone atoms should be used,
            // store N, CA and C
            else if (((ptr->m_name == "N") || (ptr->m_name == "CA")) ||
                      (ptr->m_name == "C")) {
                atoms.push_back(ptr);}}}
    return atoms;
}

// Load a Protein from a file
void PROT::Protein::load (const string& fileName, const string path = "") {
    // Construct the string of the file to be opened
    string source = path + fileName;
    // Attempt to open that file
    ifstream input; input.open(source.c_str());
    if (!input.is_open()) {
        string text = "Failure to open this Protein file:\nName: " 
                    + fileName + "\nLocation: ";
        if (path == "") {text += "current folder\n";}
        else {text += path + "\n";}
        throw logic_error (text);}
    // Store the contents of the file in this vector
    vector<string> contents;
    // Loop through the file and store its contents
    string line; getline(input, line);
    while(!input.eof()) {
        contents.push_back(line); getline(input, line);}
    // close the file
    input.close();
    // Use the load function that works on a vector of strings
    load (contents);
    return;
}

// Load a Protein from a vector of strings
void PROT::Protein::load (const vector<string>& contents) {
    // Ensure that there are strings
    if (contents.size() == 0) {
        string text = "A Protein cannot be loaded from an empty vector of "
                      "strings.\n";
        throw logic_error (text);}
    // Convert the strings to Atoms
    vector<Atom> atoms; atoms.reserve(contents.size());
    // Loop through the strings
    for(size_t i=0; i<contents.size(); ++i) {
        // Try to store the string as an Atom
        try {atoms.push_back(Atom(contents[i]));}
        // If that fails, just skip the line
        catch (logic_error) {continue;}}
    // Use the load function that works on a vector of Atoms
    load (atoms);
    return;
}

// Load a Protein from a vector of Atoms
void PROT::Protein::load (vector<Atom>& atoms) {
    // Ensure that there are Atoms
    if (atoms.size() == 0) {
        string text = "A Protein cannot be loaded from an empty vector of "
                      "Atoms.\n";
        throw logic_error (text);}
    // Convert the Atoms to Residues. Store pointers to Atoms that belong in
    // the same residue here
    vector<Atom *> res;
    // The labelling information that distinguishes the residues from one
    // another
    long n = atoms[0].m_residue_number;
    char l = atoms[0].m_insertion;
    // Store the residues in this vector
    vector<Residue> residues;
    // Loop through the Atoms
    for(size_t i=0; i<atoms.size(); ++i) {
        // If the Atom's labelling information is different
        if ((atoms[i].m_residue_number != n) ||
            (atoms[i].m_insertion != l)) {
            // Store the residue as a Residue
            residues.push_back(Residue(res));
            // Clear the residue vector
            res.clear();
            // Update the labelling information
            n = atoms[i].m_residue_number;
            l = atoms[i].m_insertion;}
        // Store a pointer to the Atom in the res vector only if it has a
        // blank or 'A' for the alternative location parameter
        Atom * ptr = &(atoms[i]);
        // Store the Atom in the vector only if it has a blank or 'A' for its
        // alternative position character
        if ((ptr->m_alt == ' ') || (ptr->m_alt == 'A')) {
            res.push_back(ptr);}}
    // If the res vector is not empty, store it as a residue
    if (res.size() > 0) {residues.push_back(Residue(res));}
    // Call the load function that works on a vector of Residues
    load(residues);
    return;
}

// Load a Protein from a vector of Residues
void PROT::Protein::load (vector<Residue>& residues) {
    // Check that the residues are acceptable for use in a protein
    check_residues(residues);
    // If the protein is currently empty
    if (m_residues == 0) {
        // Get the number of residues
        m_count = residues.size();
        // Allocate memory
        m_residues = new Residue [m_count];
        // Store the name of the protein
        m_name = residues[0].m_protein;
        // Store the Residues, updating their internal numbers as they are
        // stored
        for(size_t i=0; i<m_count; ++i) {
            m_residues[i] = residues[i];
            m_residues[i].m_internal = i+1;}}
    // If the Protein is not currently empty
    else {
        // Make sure the number of Residues matches
        if (residues.size() != m_count) {
            string text = "An existing Protein can only be loaded from a "
                          "vector of Residues that has the same number of "
                          "Residues as the Protein.\n";
            throw logic_error (text);}
        // Loop through the residues
        for(size_t i=0; i<m_count; ++i) {
            // Make a vector of pointers to the Residue's Atoms
            vector<Atom *> ptrs; ptrs.reserve(residues[i].m_count);
            for(size_t j=0; j<residues[i].m_count; ++j) {
                ptrs.push_back(&(residues[i].m_atoms[j]));}
            // Load these Atoms into the existing Residue
            m_residues[i].load(ptrs);}}
    // Ensure that the Atoms are numbered sequentially
    long n = renumber_atoms(1);
    return;
}

// Construct a Protein using various inputs
PROT::Protein::Protein (const string& fileName, const string path = "") {
    initialize(); load(fileName, path);
    return;
}

PROT::Protein::Protein (const vector<string>& contents) {
    initialize(); load(contents);
    return;
}

PROT::Protein::Protein (vector<Atom>& atoms) {
    initialize(); load(atoms);
    return;
}

PROT::Protein::Protein (vector<Residue>& residues) {
    initialize(); load(residues);
    return;
}

// Create a formatted string of text containing the Protein's contents
string PROT::Protein::str (const bool internal = false) {
    // STore the output here
    string output; output.reserve(number_of_atoms() * AtomStringLength);
    // If there are Residues
    if (m_count > 0) {
        // Loop through the Residues
        for(size_t i=0; i<m_count; ++i) {
            // Add the Residue's information to the Protein string
            output.append(m_residues[i].str(internal));}}
    return output;
}

// Calculate the number of Atoms in a Protein
size_t PROT::Protein::number_of_atoms () const {
    // Store the answer here
    size_t answer = 0;
    // If there are Residues
    if (m_count > 0) {
        for(size_t i=0; i<m_count; ++i) {answer += m_residues[i].m_count;}}
    return answer;
}

// The number of the last Atom in the Protein
long PROT::Protein::last_atom_number () const {
    // If there are residues
    if (m_count > 0) {
        // Return the value from the last residue
        return m_residues[m_count-1].last_atom_number();}
    // Otherwise return 0
    return 0;
}

// The INTERNAL number of the last Residue in the Protein
long PROT::Protein::last_residue_number () const {
    if (m_count > 0) {return m_residues[m_count-1].m_internal;}
    return 0;
}

// Change the internal numbers of the Residues
long PROT::Protein::renumber_residues (long n) {
    // If there are Residues
    if (m_count > 0) {
        // Loop through them
        for(size_t i=0; i<m_count; ++i) {
            m_residues[i].m_internal = n;
            ++n;}}
    return n;
}

// Change the numbers of the Atoms in the Protein
long PROT::Protein::renumber_atoms (long n) {
    if (m_count > 0) {
        for(size_t i=0; i<m_count; ++i) {
            n = m_residues[i].renumber_atoms (n);}}
    return n;
}

// Change the Protein's name
void PROT::Protein::set_name (const char L) {
    // Check that the name is valid
    CHECK::protein_name (L);
    // Modify the protein's attribute
    m_name = L;
    // If there are Residues, change them, too
    if (m_count > 0) {
        for(size_t i=0; i<m_count; ++i) {
            m_residues[i].p_set_protein(L);}}
    return;
}

// Provide access to the Protein's residues
PROT::Residue * PROT::Protein::operator() (const long num,
                                           const char insertion = ' ',
                                           const bool position = false) const {
    // Throw an error if the Protein is empty
    if (m_count == 0) {
        string text = "It is not possible to access a Residue in an empty "
                      "Protein.\n";
        throw logic_error (text);}
    // If the residue is being referenced by position
    if (position) {
        // If the number is in the proper range
        if ((num >= 0) && (num < m_count)) {
            return &(m_residues[num]);}
        else {
            stringstream c1; c1 << num;
            stringstream c2; c2 << m_count;
            string error = c1.str() + " is not a valid Residue index in a "
                           "Protein with " + c2.str() + " Residues.\n";
            throw logic_error (error);}}
    // Otherwise, check the residues for number / insertion code combinations
    // that match the provided values
    else {
        for(size_t i=0; i<m_count; ++i) {
            if ((m_residues[i].m_number == num) &&
                (m_residues[i].m_insertion == insertion)) {
                return &(m_residues[i]);}}
        stringstream c1; c1 << num;
        string error = "The Protein does not contain Residue " + c1.str();
        error += insertion;
        throw logic_error (error + "\n");}
    return 0;
}

// Using a Matrix, move the Protein through space
void PROT::Protein::move (const Matrix * mat, const char how = '-') {
    // Error check the matrix
    mat->move_check();
    // use the private move function
    p_move (mat, how);
    return;
}

void PROT::Protein::move (const Matrix * mat, const bool subtract = true) {
    if (subtract) {move(mat, '-');}
    else {move(mat, '+');}
    return;
}

// Or rotate the protein
void PROT::Protein::rotate (const Matrix * mat) {
    mat->rotate_check ();
    p_rotate(mat);
    return;
}

// Move a Protein so that its center of mass is at the origin
PROT::Matrix PROT::Protein::center (const bool backboneOnly = false) {
    // Create a 1 x Atom Coordinates matrix of zeros
    Matrix matrix (1, AtomCoordinates);
    // If there are no Residues, be done
    if (m_count == 0) {return matrix;}
    // Collect pointers to the relevant atoms
    vector<Atom *> ptrs = select_atoms(backboneOnly);
    if (ptrs.size() == 0) {return matrix;}
    // For each coordinate position
    for(size_t i=0; i<AtomCoordinates; ++i) {
        // Store the average value here
        coor average = 0.0;
        // Loop through the Atoms
        for(size_t j=0; j<ptrs.size(); ++j) {
            average += ptrs[j]->operator[](i);}
        matrix.set(0, i, average / ptrs.size());}
    // Move the Protein using those average coordinates
    p_move(&matrix, '-');
    return matrix;
}

// Calculate a completeness score for the Protein. This is only useful in the
// context of PDB files
double PROT::Protein::score () const {
    // Store the score here
    double value = 0.0;
    // If there are no residues, be done
    if (m_count == 0) {return value;}
    // Loop through the Residues
    for(size_t i=0; i<m_count; ++i) {
        // Add the score for the residue to the score for the protein
        value += m_residues[i].score();}
    // Return the value divided by the number of residues
    return value / ((double) m_count);
}

// Convert HIS residues to HSD and back again for use in CHARMM
void PROT::Protein::charmm_his_prep () {
    // If there are Residues
    if (m_count > 0) {
        // Loop through them
        for(size_t i=0; i<m_count; ++i) {
            // If this Residue is a Histidine
            if (m_residues[i].m_name == "HIS") {
                // Change it to HSD
                m_residues[i].m_name = "HSD";
                // If the Residue has atoms, change their residue names, too
                if (m_residues[i].m_count > 0) {
                    for(size_t j=0; j<m_residues[i].m_count; ++j) {
                        m_residues[i].m_atoms[j].m_residue = "HSD";}}}}}
    return;
}

// Convert HSD residues back to HIS
void PROT::Protein::charmm_his_fix () {
    if (m_count == 0) {return;}
    for(size_t i=0; i<m_count; ++i) {
        if (m_residues[i].m_name != "HSD") {continue;}
        m_residues[i].m_name = "HIS";
        if (m_residues[i].m_count == 0) {continue;}
        for(size_t j=0; j<m_residues[i].m_count; ++j) {
            m_residues[i].m_atoms[j].m_residue = "HIS";}}
    return;
}

// Make a string of the Protein's sequence
string PROT::Protein::get_sequence () const {
    // Store the string here
    string output;
    // Reserve an appropriate amount of space
    output.reserve(m_count * 4);
    // Loop through the protein's residues
    for(size_t i=0; i<m_count; ++i) {
        // Include the residue's name in the output string
        Text::ljust_insert(output, m_residues[i].m_name, 4, ' ');}
    // Strip trailing whitespace
    Text::rstrip(output);
    return output;
}

// Calculate the dihedral angles of the Protein's residues
void PROT::Protein::calculate_dihedrals(const bool omega = false) {
    // Specify strings of the backbone atoms
    string name1 = "N";
    string name2 = "CA";
    string name3 = "C";
    // Loop through the Protein's Residues
    for(size_t i=0; i<m_count; ++i) {
        // If the residue is not an amino acid, skip it
        if (!m_residues[i].is_amino_acid()) {continue;}
        // Get the 3 backbone Atoms
        Atom * N = m_residues[i][name1];
        Atom * CA = m_residues[i][name2];
        Atom * C = m_residues[i][name3];
        // If this is not the first Residue, calculate the phi dihedral angle
        // using the C from the previous residue
        if ((i > 0) && (m_residues[i-1].is_amino_acid())) {
            Atom * Cp = m_residues[i-1][name3];
            // Calculate and store the phi angle
            m_residues[i].m_phi = calculate_dihedral(Cp, N, CA, C);}
        // If this is not the last Residue and the next is an amino acid
        if ((i < m_count-1) && (m_residues[i+1].is_amino_acid())) {
            // Get the Nitrogen in the next Residue
            Atom * Nn = m_residues[i+1][name1];
            // Calculate psi
            m_residues[i].m_psi = calculate_dihedral(N, CA, C, Nn);
            // If the omega angle should be calculated
            if (omega) {
                Atom * Cn = m_residues[i+1][name2];
                m_residues[i].m_omega = calculate_dihedral(CA, C, Nn, Cn);}}}
    return;
}

// Calculate the RMSD between two Proteins
PROT::coor PROT::Protein::calculate_RMSD (const Protein * other, 
                          const bool backboneOnly = false) const {
    // Get the vectors of atoms
    vector<Atom *> atoms1 = select_atoms(backboneOnly);
    vector<Atom *> atoms2 = other->select_atoms(backboneOnly);
    // Raise an error if there are differing atom numbers
    if (atoms1.size() != atoms2.size()) {
        string error = "Cannot calculate the RMSD between two proteins "
                       "that selected differing numbers of atoms.\n";
        throw logic_error (error);}
    else if (atoms1.size() == 0) {
        string error = "Cannot calculate an RMSD value using no Atoms.\n";
        throw logic_error (error);}
    // Store the value here
    coor value = 0;
    // Loop through the atoms
    for(size_t i=0; i<atoms1.size(); ++i) {
        value += atoms1[i]->calculate_distance(atoms2[i], true);}
    return sqrt(value / ((coor) atoms1.size()));
}

// Calculate the number of VDW overlaps between the Protein and an Atom
size_t PROT::Protein::calculate_vdw_overlap (const Atom * atom) const {
    // Store the answer here
    size_t answer = 0;
    // If the atom is a hydrogen, there are no overlaps
    if (atom->m_name[0] == 'H') {return answer;}
    // Loop through the Protein's residues
    for(size_t i=0; i<m_count; ++i) {
        // Use the Residue function for this purpose
        answer += m_residues[i].calculate_vdw_overlap(atom);}
    return answer;
}

// VDW overlaps between the Protein and a Residue
size_t PROT::Protein::calculate_vdw_overlap (const Residue * residue,
                      const bool rotamer = false) const {
    // Store the answer here
    size_t answer = 0;
    // Go through the Protein's residues
    for(size_t i=0; i<m_count; ++i) {
        // Use the Residue's function to do the calculations
        answer += m_residues[i].calculate_vdw_overlap(residue, false, rotamer);}
    return answer;
}

// VDW overlaps between 2 Proteins
size_t PROT::Protein::calculate_vdw_overlap (const Protein * other) const {
    // Store the answer here
    size_t answer = 0;
    // Go through the other Protein's residues
    for(size_t i=0; i<other->m_count; ++i) {
        // Get a pointer to that Residue
        Residue * residue = &(other->m_residues[i]);
        // Use the Residue level function
        answer += calculate_vdw_overlap (residue, false);}
    return answer;
}

// Calculate the energy between a Protein and an Atom
float PROT::Protein::calculate_energy (const Atom * atom) const {
    // Store the energy here
    float energy = 0.0;
    // Go through the Protein's Residues
    for(size_t i=0; i<m_count; ++i) {
        // Use the Residue's function to do the calculation
        energy += m_residues[i].calculate_energy(atom);}
    return energy;
}

// Calculate the energy between a Protein and a Residue
float PROT::Protein::calculate_energy (const Residue * residue, 
                     const bool rotamer = false) const {
    // Store the energy here
    float energy = 0.0;
    // Go through the Protein's Residues
    for (size_t i=0; i<m_count; ++i) {
        // Use the input Residue's energy calculation function. This specific
        // order is necessary to get non-bonded exclusions correct when the
        // input residue is a rotamer that belongs in this protein. The
        // non-bonded exclusions require that the atom doing the calculations
        // be in a rotamer to work correctly. Thus the rotamer residue needs
        // to be the one that is controlling the calculations.
        energy += residue->calculate_energy(&(m_residues[i]), rotamer, false);}
    return energy;
}

// Calculate the energy between a pair of Proteins
float PROT::Protein::calculate_energy (const Protein * other) const {
    // Store the energy here
    float energy = 0.0;
    // Go through the other Protein's residues
    for(size_t i=0; i<other->m_count; ++i) {
        // Get a pointer to the residue
        Residue * residue = &(other->m_residues[i]);
        // Calculate the Protein's energy with that residue
        energy += calculate_energy(residue, false);}
    return energy;
}

// Create a duplicated copy of a Protein
PROT::Protein PROT::Protein::duplicate () const {
    // Create a new protein
    Protein output;
    // Copy the name
    output.m_name = m_name;
    // Copy the number of residues
    output.m_count = m_count;
    // If there are residues
    if (m_count > 0) {
        // Allocate memory
        output.m_residues = new Residue [m_count];
        // Copy the residues
        for(size_t i=0; i<m_count; ++i) {
            output.m_residues[i] = m_residues[i];}}
    return output;
}

// End the header guard from the start of the file
#endif
