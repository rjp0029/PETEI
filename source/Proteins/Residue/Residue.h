/* Written in 2018 by Dr. Robert Pantazes. Copyright owned by Auburn
 * University.
 * This file contains the implementation of the Residue class, where a Residue
 * is most typically an amino acid in a Protein. */

// Use a header guard to prevent this file from being loaded into a compiled
// program more than once
#ifndef Residue_Header_Guard
#define Residue_Header_Guard 1

// This file is intended to be included at the end of the Proteins.h header
// file, and as such it does not include any other C++ files.

// The initialization function of the Residue class
void PROT::Residue::initialize () {
    // Assign default values to all variables
    m_atoms = 0;
    m_count = 0;
    m_name = "N/A";
    m_number = 0;
    m_internal = 0;
    m_insertion = ' ';
    m_protein = ' ';
    m_present = true;
    m_missing_atoms = 0;
    m_phi = -1000.0;
    m_psi = -1000.0;
    m_omega = -1000.0;
    return;
}

// Delete any dynamically allocated memory
void PROT::Residue::clean_up () {
    if (m_atoms != 0) {delete[] m_atoms; m_atoms = 0;}
    return;
}

// Copy the information from another Residue into this one
void PROT::Residue::copy (const Residue * other) {
    // Clean up any existing information
    clean_up();
    // Get the number of Atoms in the other residue
    m_count = other->m_count;
    // If there are Atoms, allocate and copy them
    if (m_count > 0) {
        m_atoms = new PROT::Atom [m_count];
        for(size_t i=0; i<m_count; ++i) {m_atoms[i] = other->m_atoms[i];}}
    // Copy the other attributes
    m_name = other->m_name;
    m_number = other->m_number;
    m_internal = other->m_internal;
    m_insertion = other->m_insertion;
    m_protein = other->m_protein;
    m_present = other->m_present;
    m_missing_atoms = other->m_missing_atoms;
    m_phi = other->m_phi;
    m_psi = other->m_psi;
    m_omega = other->m_omega;
    return;
}

// Set the numbering information of the Residue
void PROT::Residue::p_set_number (const long n, const char i = ' ', 
                                  const bool internal = true) {
    // If this is a modification of the Residue's internal number information,
    // then only the internal number is set
    if (internal) {m_internal = n;}
    // Otherwise, update the number and insertion information
    else {m_number = n; m_insertion = i;}
    return;
}

// Set the protein information of the Residue
void PROT::Residue::p_set_protein (const char L) {
    // Store the value of the Residue
    m_protein = L;
    // If there are Atoms, update them, too
    if(m_count > 0) {
        for(size_t i=0; i<m_count; ++i) {m_atoms[i].m_protein = L;}}
    return;
}

// Determine whether or not a vector of Atom pointers are acceptable for use
// in a Residue
void PROT::Residue::check_atoms (const vector<PROT::Atom *>& atoms) const {
    // If there are no Atoms, throw an error
    if (atoms.size() == 0) {
        string error = "A Residue cannot be assembled from an empty vector "
                       "of Atoms.\n";
        throw logic_error (error);}
    // If there is only a single Atom, there cannot be naming / numbering
    // issues with it
    else if (atoms.size() == 1) {return;}
    // Make sure that every Atom has the same residue and molecule information
    for (size_t i=1; i<atoms.size(); ++i) {
        if (((atoms[0]->m_residue != atoms[i]->m_residue) || 
             (atoms[0]->m_residue_number != atoms[i]->m_residue_number)) ||
            ((atoms[0]->m_insertion != atoms[i]->m_insertion) ||
             (atoms[0]->m_protein != atoms[i]->m_protein))) {
            string error = "All Atoms in a Residue must have the same "
                           "residue and protein labelling information. These "
                           "do not:\n";
            for(size_t j=0; j<atoms.size(); ++j) {
                error += atoms[j]->str();}
            throw logic_error (error);}}
    // Confirm that each ATOM entry has a unique name
    for(size_t i=0; i<atoms.size()-1; ++i) {
        // Skip HETATM entries
        if (atoms[i]->m_type == "HETATM") {continue;}
        // Loop through all subsequent atoms
        for(size_t j=i+1; j<atoms.size(); ++j) {
            // Skip HETATM entries
            if (atoms[j]->m_type == "HETATM") {continue;}
            // Throw an error if the two atoms have the same name
            if (atoms[i]->m_name == atoms[j]->m_name) {
                string error = "Each Atom in a Residue must have a unique name."
                               " These do not:\n";
                for(size_t k=0; k<atoms.size(); ++k) {
                    error += atoms[k]->str();}
                throw logic_error (error);}}}
    return;
}

// The private move and rotate functions of the Residue class
void PROT::Residue::p_move (const Matrix * matrix, const char how) {
    // If there are Atoms, use their private move functions
    if (m_count > 0) {
        for (size_t i=0; i<m_count; ++i){m_atoms[i].private_move(matrix, how);}}
    return;
}

void PROT::Residue::p_rotate (const Matrix * matrix) {
    if (m_count > 0) {
        for(size_t i=0; i<m_count; ++i) {
            m_atoms[i].private_rotate(matrix);}}
    return;
}

// The function that loads a Residue from a vector of Atom pointers
void PROT::Residue::load (const vector<PROT::Atom *>& atoms, 
                          const bool sidechain_only = false,
                          const bool complete_load = false) {
    // First, check that the Atoms are useable
    check_atoms (atoms);
    // If this is a complete load of the Residue's information, extract the
    // meaningful data from the atoms
    if (complete_load) {
        m_name = atoms[0]->m_residue;
        m_number = atoms[0]->m_residue_number;
        m_insertion = atoms[0]->m_insertion;
        m_protein = atoms[0]->m_protein;
        // Set the internal number to 0 - it will have to be externally
        // specified
        m_internal = 0;
        // Because the Residue is being loaded from Atoms, it is definitely
        // present
        m_present = true;}
    // If a complete load is not being done and a side chain load is not being
    // done, make sure the provided atoms have the proper residue name
    if ((!complete_load) && (!sidechain_only)) {
        if (atoms[0]->m_residue != m_name) {
            string error = "It is not permitted to load a "
                         + atoms[0]->m_residue + " in place of a "
                         + m_name + ".\n";
            throw logic_error (error);}}
    // The relevant Atoms have to be collected. Store them here
    vector<Atom> use;
    // If this is an amino acid and only side chain atoms should be patched in
    if ((sidechain_only) && (CHECK::is_amino_acid(m_name))) {
        // Update the name of the residue to match that from the provided
        // atoms
        m_name = atoms[0]->m_residue;
        // Loop through the current atoms to identify the backbone atoms
        for(size_t i=0; i<m_count; ++i) {
            if (CHECK::is_backbone_atom(m_atoms[i].m_name)) {
                // Make sure proper conventions are followed for glycine and
                // proline
                if ((m_name == "GLY") && (m_atoms[i].m_name == "HA")) {
                    m_atoms[i].m_name = "HA2";}
                else if ((m_name != "GLY") && (m_atoms[i].m_name == "HA2")) {
                    m_atoms[i].m_name = "HA";}
                else if (m_name == "PRO") {
                    if (m_atoms[i].m_name == "HN") {continue;}
                    else if (m_atoms[i].m_name == "HT3") {continue;}
                    else if (m_atoms[i].m_name == "HT1") {
                        m_atoms[i].m_name = "HN1";}
                    else if (m_atoms[i].m_name == "HT2") {
                        m_atoms[i].m_name = "HN2";}}
                else if ((m_atoms[i].m_name == "HN1") && (m_name != "PRO")) {
                    m_atoms[i].m_name = "HT1";}
                else if ((m_atoms[i].m_name == "HN2") && (m_name != "PRO")) {
                    m_atoms[i].m_name = "HT2";}
                use.push_back(m_atoms[i]);}}
        // Loop through the provided atoms to identify the non-backbone atoms
        for(size_t i=0; i<atoms.size(); ++i) {
            if (!CHECK::is_backbone_atom(atoms[i]->m_name)) {
                use.push_back(*(atoms[i]));}}}
    // Otherwise, just use all of the atoms
    else {
        for(size_t i=0; i<atoms.size(); ++i) {
            use.push_back(*(atoms[i]));}}
    // Delete existing atoms from the residue
    clean_up();
    // Allocate space for the new set of atoms
    m_count = use.size();
    m_atoms = new Atom [m_count];
    // Store the appropriate atoms
    for(size_t i=0; i<m_count; ++i) {m_atoms[i] = use[i];}
    // Make sure that every atom has the proper labelling information
    for(size_t i=0; i<m_count; ++i) {
        m_atoms[i].m_alt = ' ';
        m_atoms[i].m_residue = m_name;
        m_atoms[i].m_residue_number = m_number;
        m_atoms[i].m_insertion = m_insertion;
        m_atoms[i].m_protein = m_protein;}
    // Ensuring consistent Atom numbering is done outside of this function
    return;
    }

// construct a residue from a vector of atoms
PROT::Residue::Residue (const vector<PROT::Atom *>& atoms) {
    // Initialize the residue
    initialize();
    // Load the atoms
    load(atoms, false, true);
    return;
}

// Construct a Residue from a string name and the protein's character
PROT::Residue::Residue (const string& input, const char L) {
    // Initialize the Residue
    initialize();
    // Set the Residue's name
    m_name = input;
    // Set the Protein's name
    m_protein = L;
    // The remaining attributes can be left as their initialized values
    return;
}

// Access to the Atoms of the Residue
PROT::Atom * PROT::Residue::operator[] (const size_t i) {
    // If there are no atoms, raise an error
    if (m_count == 0) {
        string error = "It is not possible to access an Atom in a Residue "
                       "when the Residue is empty.\n";
        throw logic_error (error);}
    // If the index is in the proper range
    else if ((i >= 0) && (i < m_count)) {return &(m_atoms[i]);}
    else if ((i < 0) && (i >= -m_count)) {return &(m_atoms[m_count + i]);}
    // Otherwise, throw an error
    else {
        stringstream c1; c1 << m_count;
        stringstream c2; c2 << i;
        string error = "This Residue contains " + c1.str() + " Atoms. "
                     + c2.str() + " is not a valid index.\n";
        for(size_t j=0; j<m_count; ++j) {error += m_atoms[j].str();}
        throw logic_error (error);}
    return 0;
}

// Access Atoms using a string
PROT::Atom * PROT::Residue::operator[] (const string& label) {
    // If there are no atoms, raise an error
    if (m_count == 0) {
        string error = "It is not possible to access an Atom in a Residue "
                       "when the Residue is empty.\n";
        throw logic_error (error);}
    // Find the atom
    for(size_t i=0; i<m_count; ++i) {
        if (m_atoms[i].m_name == label) {return &(m_atoms[i]);}}
    // If there was no such atom
    string error = "This Residue does not contain a " + label + " Atom.\n";
    for (size_t i=0; i<m_count; ++i) {error += m_atoms[i].str();}
    throw logic_error (error);
    return 0;
}

// Specify numbering information for a residue, while ensuring that the values
// are acceptable
void PROT::Residue::set_number (const long n, const char i = ' ',
                                const bool internal = true) {
    // Check the number
    CHECK::residue_number (n);
    // Check the insertion code
    CHECK::insertion_code (i);
    // Do the renumbering of the residue
    p_set_number (n, i, internal);
    return;
}

// Set the Residue's protein name
void PROT::Residue::set_protein (const char L) {
    CHECK::protein_name (L);
    p_set_protein (L);
    return;
}

// Renumber the Atoms in the Residue
long PROT::Residue::renumber_atoms (long n) {
    // If there are Atoms
    if (m_count > 0) {
        // Loop through them
        for(size_t i=0; i<m_count; ++i) {
            // modify the atom's number
            m_atoms[i].m_number = n;
            // Increment n
            ++n;}}
    return n;
}

// Get the number of the last Atom in the Residue
long PROT::Residue::last_atom_number () const {
    // If there are atoms
    if (m_count > 0) {return m_atoms[m_count-1].m_number;}
    // Otherwise return 0
    return 0;
}

// Generate a string of the Residue's information properly formatted for
// internal or external use
string PROT::Residue::str (bool internal = false) {
    // Store the output string here
    string output;
    // If there are atoms
    if (m_count > 0) {
        // Allocate an appropriate amount of space in the string
        output.reserve(AtomStringLength * m_count);
        // Loop through the Atoms
        for(size_t i=0; i<m_count; ++i) {
            // If internal numbering should be used
            if (internal) {
                m_atoms[i].m_residue_number = m_internal;
                m_atoms[i].m_insertion = ' ';}
            else {
                m_atoms[i].m_residue_number = m_number;
                m_atoms[i].m_insertion = m_insertion;}
            // Add the Atom's string to the output
            output.append(m_atoms[i].str());}}
    return output;
}

// Move and rotate the Residue
void PROT::Residue::move (const Matrix * matrix, const char how = '-') {
    matrix->move_check();
    p_move(matrix, how);
    return;
}

void PROT::Residue::move (const Matrix * matrix, const bool subtract) {
    if (subtract) {move(matrix, '-');}
    else {move(matrix, '+');}
    return;
}

void PROT::Residue::rotate (const Matrix * matrix) {
    matrix->rotate_check();
    p_rotate (matrix);
    return;
}

// Move a Residue so that its center of mass is at the origin
PROT::Matrix PROT::Residue::center () {
    // Allocate a 1 x AtomCoordinates matrix
    Matrix output (1, AtomCoordinates);
    // If there are atoms
    if (m_count > 0) {
        // Loop through the coordinates
        for(size_t i=0; i<AtomCoordinates; ++i) {
            // Add up the values
            coor n = 0.0;
            for(size_t j=0; j<m_count; ++j) {
                n += m_atoms[j].m_coors[i];}
            // Calculate the average value
            n /= m_count;
            // Store it in the matrix
            output.set(0, i, n);}
        // Use the private move function of the class
        p_move(&output, '-');}
    return output;
}

// Position a Residue in a standard way to facilitate rotamer positioning
void PROT::Residue::position (PROT::Matrix& how, const bool forward = true) {
    // Do one set of operations if the Residue is being assigned the initial
    // position
    if (forward) {
        // This function only works on Amino Acid residues (permit Histidine
        // variations, too)
        if ((!CHECK::is_amino_acid (m_name)) && (m_name != "HSD")){
            // Throw an error
            string error = "The Residue position function only works on amino "
                           "acids, not " + m_name + ".\n";
            throw logic_error (error);}
        // The three Atoms that are used in the positioning
        string names [3] = {"CA", "N", "CB"};
        if (m_name == "GLY") {names[2] = "HA1";}
        Atom * atoms [3] = {0, 0, 0};
        // Find the atoms with those names
        for(size_t i=0; i<3; ++i) {
            for(size_t j=0; j<m_count; ++j) {
                if (m_atoms[j].m_name == names[i]) {
                    atoms[i] = &(m_atoms[j]);
                    break;}}
            // If the Atom isn't present, throw an error
            if (atoms[i] == 0) {
                string error = "This Residue does not have the proper Atoms to "
                               "be positioned:\n";
                for(size_t j=0; j<m_count; ++j) {
                    error += m_atoms[i].str();}
                throw logic_error (error);}}
        // Allocate the how matrix to be the proper dimensions
        how.allocate(1, AtomCoordinates + 3);
        // Create a matrix using the coordinates of the first atom
        Matrix first (atoms[0]);
        // Move the Residue so that first atom is at the origin
        p_move(&first, '-');
        // Store the coordinates in the how matrix
        for(size_t i=0; i<AtomCoordinates; ++i) {how.set(0, i, first(0, i));}
        // Calculate the angle of rotation to rotate the second Atom into the
        // XZ plane
        coor angle1 = fabs(acos(atoms[1]->m_coors[0] /
                           sqrt(pow(atoms[1]->m_coors[0], 2) 
                              + pow(atoms[1]->m_coors[1], 2))));
        if (atoms[1]->m_coors[1] > 0) {angle1 = -angle1;}
        // A unit vector for the Z axis
        coor Z [3] = {0.0, 0.0, 1.0};
        // Create a rotation matrix
        Matrix second (angle1, Z);
        // Rotate the Residue
        p_rotate(&second);
        // Store the angle in the how matrix
        how.set(0, AtomCoordinates, angle1);
        // Rotate the seconde Atom onto the Z-axis by rotating around the
        // y-axis
        coor angle2 = fabs(acos(atoms[1]->m_coors[2] /
                           sqrt(pow(atoms[1]->m_coors[2], 2) 
                              + pow(atoms[1]->m_coors[0], 2))));
        if (atoms[1]->m_coors[0] > 0) {angle2 = -angle2;}
        coor Y [3] = {0.0, 1.0, 0.0};
        Matrix third (angle2, Y);
        p_rotate(&third);
        how.set(0, AtomCoordinates+1, angle2);
        // Rotate the third atom into the XZ plane
        coor angle3 = fabs(acos(atoms[2]->m_coors[0] / 
                           sqrt(pow(atoms[2]->m_coors[0], 2)
                              + pow(atoms[2]->m_coors[1], 2))));
        if (atoms[2]->m_coors[1] > 0) {angle3 = -angle3;}
        Matrix fourth (angle3, Z);
        p_rotate(&fourth);
        how.set(0, AtomCoordinates+2, angle3);
        // The residue is now positioned
    }
    // If this is not the initial positioning, reverse the calculations
    else {
        coor Z [3] = {0.0, 0.0, 1.0}; coor Y [3] = {0.0, 1.0, 0.0};
        Matrix fourth (-how(0, AtomCoordinates+2), Z);
        p_rotate(&fourth);
        Matrix third (-how(0, AtomCoordinates+1), Y);
        p_rotate(&third);
        Matrix second (-how(0, AtomCoordinates), Z);
        p_rotate(&second);
        Matrix first (1, AtomCoordinates);
        for(size_t i=0; i<AtomCoordinates; ++i) {first.set(0, i, how(0, i));}
        p_move(&first, '+');}
    return;
}

// Calculate a score for the Residue that tells how "complete" it is
double PROT::Residue::score () const {
    // If the Residue is missing or has no Atoms, the score is 0
    if ((!m_present) || (m_count == 0)) {return 0.0;}
    // If the Residue is present and has no missing Atoms, the score is 1
    else if (m_missing_atoms == 0) {return 1.0;}
    // Otherwise, calculate the score as follows
    return ((double) m_count / ((double) (m_count + m_missing_atoms)));
}

// Access to the Residue's dihedral angles
PROT::coor PROT::Residue::phi () const {
    // Throw an error if the angle was never assigned
    if (m_phi < -999) {
        string error = "No phi angle available.\n";
        throw logic_error (error);}
    return m_phi;
}

PROT::coor PROT::Residue::psi () const {
    if (m_psi < -999) {
        string error = "No psi angle available.\n";
        throw logic_error (error);}
    return m_psi;
}

PROT::coor PROT::Residue::omega () const {
    if (m_omega < -999) {
        string error = "No omega angle available.\n";
        throw logic_error (error);}
    return m_omega;
}

// If the Residue's name should be HIS, fix that.
void PROT::Residue::charmm_his_fix () {
    if ((m_name == "HSD") || (m_name == "HSE")) {
        m_name = "HIS";
        if (m_count > 0) {
            for(size_t i=0; i<m_count; ++i) {
                m_atoms[i].m_residue = "HIS";}}}
    return;
}

// Create another Residue that has the same contents as this Residue
PROT::Residue PROT::Residue::duplicate () const {
    // Create the other Residue
    Residue other;
    // Copy the contents of this Residue to that Residue
    other.m_count = m_count;
    other.m_name = m_name;
    other.m_number = m_number;
    other.m_internal = m_internal;
    other.m_insertion = m_insertion;
    other.m_protein = m_protein;
    other.m_present = m_present;
    other.m_missing_atoms = m_missing_atoms;
    other.m_phi = m_phi;
    other.m_psi = m_psi;
    other.m_omega = m_omega;
    // If there are Atoms, allocate them
    if (m_count > 0) {
        other.m_atoms = new PROT::Atom [m_count];
        for(size_t i=0; i<m_count; ++i) {
            other.m_atoms[i] = m_atoms[i];}}
    // Return the other residue
    return other;
}

// Calculate the VDW overlap of the Residue's Atoms with another Atom
size_t PROT::Residue::calculate_vdw_overlap (const Atom * atom, 
                      const bool rotamer = false) const {
    // If the other Atom is a hydrogen, there are no overlaps
    if (atom->m_name[0] == 'H') {return 0;}
    // Store the count here
    size_t answer = 0;
    // Go through the Residue's atoms
    for(size_t i=0; i<m_count; ++i) {
        // If this Atom is a Hydrogen, skip it
        if (m_atoms[i].m_name[0] == 'H') {continue;}
        // If the Residue is a rotamer and this is a backbone atom, skip it
        else if ((rotamer) && (m_atoms[i].is_backbone_atom())) {continue;}
        // If the two Atoms have a VDW overlap
        if (m_atoms[i].check_vdw_overlap(atom)) {++answer;}}
    return answer;
}

// Calculate the VDW overlap between 2 Residues
size_t PROT::Residue::calculate_vdw_overlap (const Residue * other,
                      const bool selfRotamer = false,
                      const bool otherRotamer = false) const {
    // Store the answer here
    size_t answer = 0;
    // Go through the other Residue's atoms
    for(size_t i=0; i<other->m_count; ++i) {
        // Get an atom pointer
        Atom * atom = &(other->m_atoms[i]);
        // If the other Residue is a rotamer and the atom is a backbone atom,
        // skip it
        if ((otherRotamer) && (atom->is_backbone_atom())) {continue;}
        // use the atom function
        answer += calculate_vdw_overlap(atom);}
    return answer;
}

// Calculate the energy between a Residue and an Atom
float PROT::Residue::calculate_energy (const Atom * atom, 
                                       const bool rotamer = false) const {
    // Store the energy here
    float energy = 0.0;
    // Go through the Residue's atoms
    for (size_t i=0; i<m_count; ++i) {
        // If this is a rotamer and this atom is a backbone atom, exclude it
        // from the calculations
        if ((rotamer) && (m_atoms[i].is_backbone_atom())) {continue;}
        // Calculate the energy of this atom with the other atom
        energy += m_atoms[i].calculate_energy(atom);}
    return energy;
}

// Calculate the energy between a Residue and another Residue
float PROT::Residue::calculate_energy (const Residue * other,
                     const bool selfRotamer = false, 
                     const bool otherRotamer = false) const {
    // Store the energy here
    float energy = 0.0;
    // Go through the other Residue's atoms
    for(size_t i=0; i<other->m_count; ++i) {
        // Get a pointer to the atom
        Atom * atom = &(other->m_atoms[i]);
        // If the other residue is a rotamer and this is a backbone atom,
        // exclude it from the calculations
        if ((otherRotamer) && (atom->is_backbone_atom())) {continue;}
        // Do the Atom-level energy calculation
        energy += calculate_energy(atom, selfRotamer);}
    return energy;
}

// End the header guard from the start of the file
#endif
