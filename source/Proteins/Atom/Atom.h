/* Written in 2018 by Dr. Robert Pantazes. Copyright owned by Auburn
 * University. 
 * This file contains the definition and implementation of a class to store /
 * manipulate the data about an Atom from a PDB file. */

// Use a header guard to make sure this file is not included in a compiled
// program multiple times
#ifndef Atom_Guard
#define Atom_Guard 1

// This file is intended to be included directly by the Proteins.h header file
// and does not need to include any C++ content itself.

// Assign default values to the Atom's attributes
void PROT::Atom::initialize () {
    // PDB attributes
    for(size_t i=0; i<AtomCoordinates; ++i) {m_coors[i] = 0.0;}
    m_type = "ATOM";
    m_name = "NONE";
    m_residue = "NON";
    m_element = "";
    m_charge = "";
    m_number = 1;
    m_residue_number = 1;
    m_occupancy = 1.0;
    m_temperature = 0.0;
    m_alt = ' ';
    m_insertion = ' ';
    m_protein = ' ';
    // Energy calculation attributes
    vdw = false;
    vdwR = 0.0; 
    vdwE = 0.0;
    vdw14R = 0.0;
    vdw14E = 0.0; 
    vdwPhi = 0.0;
    elec = false;
    elec_charge = 0.0;
    gb = false;
    gbRadius = 0.0;
    lk = false;
    lkR = 0.0;
    lkL = 0.0;
    lkV = 0.0;
    lkG = 0.0;
    m_core = false;
    return;
}

// Copy the information from another Atom into this atom
void PROT::Atom::copy (const Atom& other) {
    // PDB attributes
    for(size_t i=0; i<AtomCoordinates; ++i) {m_coors[i] = other.m_coors[i];}
    m_type = other.m_type;
    m_name = other.m_name;
    m_residue = other.m_residue;
    m_element = other.m_element;
    m_charge = other.m_charge;
    m_number = other.m_number;
    m_residue_number = other.m_residue_number;
    m_occupancy = other.m_occupancy;
    m_temperature = other.m_temperature;
    m_alt = other.m_alt;
    m_insertion = other.m_insertion;
    m_protein = other.m_protein;
    // Energy calculation attributes
    vdw = other.vdw;
    vdwR = other.vdwR;
    vdwE = other.vdwE;
    vdw14R = other.vdw14R;
    vdw14E = other.vdw14E;
    vdwPhi = other.vdwPhi;
    elec = other.elec;
    elec_charge = other.elec_charge;
    gb = other.gb;
    gbRadius = other.gbRadius;
    lk = other.lk;
    lkR = other.lkR;
    lkL = other.lkL;
    lkV = other.lkV;
    lkG = other.lkG;
    m_core = other.m_core;
    return;
}

// Move the Atom without error checking the inputs
void PROT::Atom::private_move (const Matrix * matrix, const char how) {
    // The default behaviour is to subtract the coordinates from the matrix
    if (how == '-') {
        for(size_t i=0; i<AtomCoordinates; ++i) {
            m_coors[i] -= matrix->operator()(0, i);}}
    // If they aren't being subtracted, they should be added
    else {
        for(size_t i=0; i<AtomCoordinates; ++i) {
            m_coors[i] += matrix->operator()(0, i);}}
    return;
}

// Rotate the Atom's position without error checking the matrix
void PROT::Atom::private_rotate (const Matrix * matrix) {
    // Create an array of the new coordinates
    coor newCoors [AtomCoordinates];
    // Calculate the new coordinates
    for(size_t i=0; i<AtomCoordinates; ++i) {
        // Set the new coordinate value to 0
        newCoors[i] = 0;
        // Loop through the current coordinates
        for(size_t j=0; j<AtomCoordinates; ++j) {
            // Add the appropriate product value to the new coordinate
            newCoors[i] += (matrix->operator()(i, j) * m_coors[j]);}}
    // Store the new coordinates as the Atom's coordinates
    for(size_t i=0; i<AtomCoordinates; ++i) {m_coors[i] = newCoors[i];}
    return;
}

// Move and rotate the Atom, doing error checking (using the Matrix class to
// do so) to ensure that the calculations are allowed
void PROT::Atom::move (const Matrix * matrix, const char how = '-') {
    matrix->move_check();
    private_move(matrix, how);
    return;
}

void PROT::Atom::rotate (const Matrix * matrix) {
    matrix->rotate_check ();
    private_rotate(matrix);
    return;
}

// Provide access to the Atom's coordinates using an index
PROT::coor PROT::Atom::operator[] (const size_t i) const {
    // If the coordinate is between 0 and AtomCoordinates
    if ((i >= 0) && (i < AtomCoordinates)) {return m_coors[i];}
    // Permit negative indexing as well
    else if ((i < 0) && (i >= -AtomCoordinates)) {
        return m_coors[i+AtomCoordinates];}
    // Throw an error for anything else
    stringstream c1; c1 << AtomCoordinates;
    stringstream c2; c2 << i;
    string error = "Atoms have " + c1.str() + " coordinates.\n"
                 + c2.str() + " is not an acceptable value.\n";
    throw logic_error (error);
    // Return a generic value so the function compiles correctly
    return 0.0;
}

// Provide access to the Atom's coordinates using a character
PROT::coor PROT::Atom::operator[] (const char L) const {
    // Permit access to the x, y, and z coordinates
    if ((L == 'x') || (L == 'X')) {return m_coors[0];}
    else if ((L == 'y') || (L == 'Y')) {return m_coors[1];}
    else if ((L == 'z') || (L == 'Z')) {return m_coors[2];}
    // throw an error
    string error = "'";
    error += L;
    error += "' is not a recognized Atom coordinate.\n";
    throw logic_error(error);
    // Return a generic value so the function compiles correctly
    return 0.0;
}

// The standard constructor of the Atom class uses a string of text
PROT::Atom::Atom (const string& input) {
    // Assign initial values to the Atom
    initialize ();
    // Create a copy of the string that can be modified and analyzed
    string line = input;
    // Strip the whitespace and convert it to all capital letters
    Text::strip(line); Text::upper(line);
    // Put everything in a try statement
    try {
        // There must be at least 66 columns
        if (line.size () < 66) {
            string error = "This text is too short to contain a PDB atom.\n";
            throw logic_error (error);}
        // Determine the Atom's type
        if (Text::startswith(line, "ATOM  ")) {m_type = "ATOM";}
        else if (Text::startswith(line, "HETATM")) {m_type = "HETATM";}
        else {
            string error = "PDB atom lines must start with ATOM or HETATM.\n";
            throw logic_error (error);}
        // Get the Atom's number
        string value = line.substr(6, 5); Text::strip(value);
        if (Text::is_integer(value)) {
            stringstream c; c << value; c >> m_number;
            CHECK::atom_number (m_number);}
        else {
            string error = "Columns 7-11 do not contain an integer.\n";
            throw logic_error (error);}
        // Get the Atom's name
        m_name = line.substr(12, 4); Text::strip(m_name);
        CHECK::atom_name(m_name);
        // Get the alternative location information
        m_alt = line[16];
        CHECK::alt_location(m_alt);
        // Get the residue's name
        m_residue = line.substr(17, 3); Text::strip(m_residue);
        CHECK::residue_name (m_residue);
        // Get the Protein's name
        m_protein = line[21];
        CHECK::protein_name(m_protein);
        // Get the Residue's number
        value = line.substr(22, 4); Text::strip(value);
        if (Text::is_integer(value)) {
            stringstream c; c << value; c >> m_residue_number;
            CHECK::residue_number (m_residue_number);}
        else {
            string error = "Columns 23-26 do not contain an integer.\n";
            throw logic_error (error);}
        // Get the Residue's insertion code
        m_insertion = line[26];
        CHECK::insertion_code (m_insertion);
        // Get the Atom's coordinates (note that AtomCoordinates is not used
        // here because PDB files ONLY have 3 coordinates)
        for(size_t i=0; i<3; ++i) {
            value = line.substr(30 + (i*8), 8); Text::strip(value);
            if (Text::is_number(value)) {
                stringstream c; c << value; c >> m_coors[i];
                CHECK::atom_coordinate(m_coors[i]);}
            else {
                string error = "PDB atom coordinates have to be numbers.\n";
                throw logic_error (error);}}
        // Get and store the Atom's occupancy
        value = line.substr(54, 6); Text::strip(value);
        if (Text::is_number(value)) {
            stringstream c; c << value; c >> m_occupancy;
            CHECK::occupancy(m_occupancy);}
        else {
            string error = "Columns 55-60 do not contain a number.\n";
            throw logic_error (error);}
        // Get and store the Atom's temperature
        value = line.substr(60, 6); Text::strip(value);
        if (Text::is_number(value)) {
            stringstream c; c << value; c >> m_temperature;
            CHECK::temperature(m_temperature);}
        else {
            string error = "Columns 61-66 do not contain a number.\n";
            throw logic_error (error);}
        // The element and charge information may not be present
        if (line.size() >= 78) {
            m_element = line.substr(76, 2); Text::strip(m_element);
            CHECK::element(m_element);
            if (line.size() >= 80) {
                m_charge = line.substr(78, 2); Text::strip(m_charge);
                CHECK::charge (m_charge);}}
    // Catch and handle any errors that occurred during this process
    } catch (logic_error& e) {
        // Store the initial input and strip it of whitespace
        string data = input; Text::strip(data);
        // Throw a new error with the same message but also the input line
        throw logic_error (e.what() + data + "\n");}
    // End the function definition
    return;
}

// Generate a string of PDB-formatted text containing the Atom's information
string PROT::Atom::str () const {
    // The output is stored here
    string output;
    // Reserve the appropriate number of characters
    output.reserve(AtomStringLength);
    // Include the Atom's attributes with proper formatting and spacing
    // The Atom's type
    Text::ljust_insert (output, m_type, 6, ' ');
    // The Atom's number
    stringstream c1; c1 << m_number;
    Text::rjust_insert(output, c1.str(), 5, ' ');
    // A blank space
    output.push_back(' ');
    // The Atom's name
    if (m_name.size() < 4) {
        output.push_back(' '); Text::ljust_insert(output, m_name, 3, ' ');}
    else {Text::ljust_insert(output, m_name, 4, ' ');}
    // The alternate location character
    output.push_back(m_alt);
    // The Residue's name
    Text::ljust_insert(output, m_residue, 3, ' ');
    // A blank space
    output.push_back (' ');
    // The protein's name
    output.push_back(m_protein);
    // The Residue's number
    stringstream c2; c2 << m_residue_number;
    Text::rjust_insert(output, c2.str(), 4, ' ');
    // The residue's insertion code
    output.push_back(m_insertion);
    // Three blank spaces
    for(size_t i=0; i<3; ++i) {output.push_back(' ');}
    // The Atom's coordinates. As with loading from the PDB file, 3
    // coordinates are explicitly used here
    for(size_t i=0; i<3; ++i) {
        stringstream c3; c3 << fixed << setprecision(3) << m_coors[i];
        Text::rjust_insert(output, c3.str(), 8, ' ');}
    // Occupancy
    stringstream c4; c4 << fixed << setprecision(2) << m_occupancy;
    Text::rjust_insert(output, c4.str(), 6, ' ');
    // Temperature
    stringstream c5; c5 << fixed << setprecision(2) << m_temperature;
    Text::rjust_insert(output, c5.str(), 6, ' ');
    // The Atom's element
    Text::rjust_insert(output, m_element, 12, ' ');
    // The Atom's charge
    Text::rjust_insert(output, m_charge, 2, ' ');
    // Add an end line character to terminate the string
    output.push_back('\n');
    return output;
}

// Calculate the distance between 2 atoms
PROT::coor PROT::Atom::calculate_distance (const Atom * other, 
                                           const bool squared = false) const {
    // The distance will be stored here
    coor value = 0.0;
    // Loop through the coordinates
    for(size_t i=0; i<AtomCoordinates; ++i) {
        // Add the square of the difference between the coordinates
        value += pow(m_coors[i] - other->m_coors[i], 2);}
    // If appropriate, return just the squared sum
    if (squared) {return value;}
    // Calculate the square root of the sum
    return sqrt(value);
}

// Determine whether or not non-bonded exclusions should be applied in the
// calculation of energies between 2 Atoms
int PROT::Atom::exclusions (const Atom * other) const {
    // As currently implemented, this function is ONLY correct if THIS atom is
    // a side-chain atom in a Rotamer. Currently, this is the only instance
    // when intra-protein energies are calculated. If that is ever changed,
    // this function will need to be changed significantly

    // If the two Atoms are in different proteins, they definitely don't need
    // non-bonded exclusions
    if (m_protein != other->m_protein) {return 0;}
    // If the two Atoms are at least 2 Residues apart from one another,
    // non-bonded exclusions aren't needed
    else if ((other->m_residue_number < m_residue_number - 1) || 
             (other->m_residue_number > m_residue_number + 1)) {return 0;}
    // If the two atoms are in the same residue and are both side chain atoms,
    // make sure energy calculations are NOT done
    else if ((other->m_residue_number == m_residue_number) && 
            ((!other->is_backbone_atom()) && (!is_backbone_atom()))) {return -1;}
    // If the other atom is in the previous residue
    else if (other->m_residue_number == m_residue_number - 1) {
        // Proline gets special handling, because its delta carbon is bound to
        // its nitrogen. This opens up the possibility of non-bonded
        // exclusions that aren't available to other amino acids
        if (m_residue == "PRO") {
            // If this is the delta carbon
            if (m_name == "CD") {
                // If the other atom is C
                if (other->m_name == "C") {return 3;}
                else if ((other->m_name == "CA") || (other->m_name == "O")) {
                    return 4;}}
            // If this is the Gamma carbon or a delta carbon hydrogen
            else if ((m_name == "CG") || ((m_name == "HD1") || 
                     (m_name == "HD2"))) {
                // In that case, a non-bonded exclusion of 4 is used with the
                // C from the previous residue
                if (other->m_name == "C") {return 4;}}}
        // Even for prolines, if this atom is the beta carbon (or HA1 for
        // glycines) and the other Atom is C, a non-bonded exclusion of 4 is
        // correct
        if (((m_name == "CB") || (m_name == "HA1")) && (other->m_name == "C")) {
            return 4;}
        // In all other cases where the other atom is in the previous residue,
        // non-bonded exclusions are not appropriate
        return 0;}
    // If the other atom is in the next residue
    else if (other->m_residue_number == m_residue_number + 1) {
        // If this atom is the beta carbon or HA1 and the next atom is a
        // nitrogen, a non-bonded exclusion of 4 is correct
        if (((m_name == "CB") || (m_name == "HA1")) && (other->m_name == "N")) {
            return 4;}
        // In all other cases, no non-bonded exclusion is appropriate
        return 0;}
    // If the function has reached this point, we know that the the two Atoms
    // are in the same residue in the same protein. Because of the delta
    // carbon to nitrogen bond of prolines, they need special handling
    else if (m_residue == "PRO") {
        // Non-proline, N-terminal amino acids that are having proline
        // considered as a rotamer have an "extra" atom. Use the smallest
        // non-bonded exclusion value without question to ensure no energies
        // are calculated with that atom
        if (other->m_name == "HT3") {return 2;}
        // A subsequent portion of this function will calculate the values
        // from atoms off of the alpha carbon. Here, do the calculations for
        // atoms off of the nitrogen
        else if (other->m_name == "N") {
            if (m_name == "CD") {return 2;}
            else if ((m_name == "CG") || 
                    ((m_name == "HD1") || (m_name == "HD2"))) {return 3;}
            // The beta carbon is a 3, going through the alpha carbon. That
            // gets picked up later.
            else if ((m_name == "HG1") || (m_name == "HG2")) {return 4;}}
        // Backbone Atoms bound to the nitrogen
        else if ((((other->m_name == "CA") || (other->m_name == "HN")) ||
                  ((other->m_name == "HN1") || (other->m_name == "HN2"))) ||
                  ((other->m_name == "HT1") || (other->m_name == "HT2"))) {
            if (m_name == "CD") {return 3;}
            // The gamma carbon with the alpha carbon is a 3 via the beta
            // carbon. Get that correct here
            else if ((m_name == "CG") || (other->m_name == "CA")) {return 3;}
            else if ((m_name == "CG") || 
                    ((m_name == "HD1") || (m_name == "HD2"))) {return 4;}}
        // Backbone Atoms bound to those atoms in this residue
        else if ((other->m_name == "C") || 
                ((other->m_name == "HA") || (other->m_name == "HA2"))) {
            if (m_name == "CD") {return 4;}}}
    // Search through the possibilities for within the same residue via
    // alpha-carbon connections. This is also done for prolines, since they
    // still have the same alpha-carbon connections
    // Rather than doing many complicated if statements, arrays of values will
    // be checked for matching residue names. First are the backbone atom
    // names
    string first [1] = {"CA"};
    string second [4] = {"C", "N", "HA", "HA2"};
    string third [9] = {"O", "OT1", "OT2", "HN", "HN1", "HN2", "HT1", "HT2",
                        "HT3"};
    // And now the side chain atom names
    string beta [2] = {"CB", "HA1"};
    string gamma [10] = {"CG", "CG1", "CG2", "HB", "HB1", "HB2", "HB3", "OG", 
                         "OG1", "SG"};
    string delta [17] = {"CD", "CD1", "CD2", "HG", "HG1", "HG2", "HG11", "HG12",
                         "HG13", "HG21", "HG22", "HG23", "OD1", "OD2", "ND1",
                         "ND2", "SD"};
    // Check the Beta Carbon Atoms
    for(size_t i=0; i<2; ++i) {
        if (m_name == beta[i]) {
            for(size_t j=0; j<1; ++j) {
                if (other->m_name == first[j]) {return 2;}}
            for(size_t j=0; j<4; ++j) {
                if (other->m_name == second[j]) {return 3;}}
            for(size_t j=0; j<9; ++j) {
                if (other->m_name == third[j]) {return 4;}}}}
    // Gamma position atoms
    for(size_t i=0; i<10; ++i) {
        if (m_name == gamma[i]) {
            for(size_t j=0; j<1; ++j) {
                if (other->m_name == first[j]) {return 3;}}
            for(size_t j=0; j<4; ++j) {
                if (other->m_name == second[j]) {return 4;}}}}
    // Delta position atoms
    for(size_t i=0; i<17; ++i) {
        if (m_name == delta[i]) {
            for(size_t j=0; j<1; ++j) {
                if (other->m_name == first[j]) {return 4;}}}}
    // In any other circumstance, return 0
    return 0;
}

// Calculate the VDW energy between Atoms
float PROT::Atom::VDW (const Atom * other, const coor dis, 
                       const int skip) const {
    // Check all the circumstances where VDW energy should not be calculated
    if (((vdw == false) || (dis >= 12.0)) || ((skip != 0) && (skip != 4))) {
        return 0.0;}
    // VDW is currently the one energy function where it may be on in some
    // atoms and off in others - in particular, it can be turned off for
    // hydrogens. So check to make sure that the other Atom's VDW is not false
    if (other->vdw == false) {return 0.0;}
    // Three parameters are needed: min radius, epsilon, and a softening term
    // (phi).
    float eps0, rmin; float phi = vdwPhi;
    // If skip is 4, use the 1-4 values
    if (skip == 4) {
        eps0 = sqrt(vdw14E * other->vdw14E);
        rmin = vdw14R + other->vdw14R;}
    // Otherwise use the standard values
    else {
        eps0 = sqrt(vdwE * other->vdwE);
        rmin = vdwR + other->vdwR;}
    // If the distance is greater than the minimum radius value, set phi to 1
    if (dis >= rmin) {phi = 1.0;}
    // Calculate and return the energy
    float energy = phi * eps0 * (pow(rmin/dis, 12) - 2.0*pow(rmin/dis, 6));
    return (-eps0 + phi*eps0) + energy;
}

// Calculate the electrostatics energy between two atoms
float PROT::Atom::ELEC (const Atom * other, const coor dis, 
                        const int skip) const {
    // If electrostatics energy should not be calculated, return 0
    if (((elec == false) || (dis >= 12.0)) || ((skip != 0) && (skip != 4))){
        return 0.0;}
    // Calculate a numerator term (CCELEC is a constant defined in Proteins.h)
    float numerator = elec_charge * other->elec_charge * CCELEC;
    // Calculate and return the energy. 
    float energy = pow(1.0 - pow(dis/12.0, 2), 2) * numerator / dis;
    return energy;
}

// Calculate the Generalized-Born implicit solvation energy between 2 atoms.
// This function comes from:
// B. Dominy and C.L. Brooks, III. Development of a Generalized Born Model
// Parameterization for Proteins and Nucleic Acids. J. Phys. Chem. 103,
// 3765-3773 (1999).
float PROT::Atom::GB (const Atom * other, const coor dis, 
                      const int skip) const {
    // If the value should not be calculated
    if (((gb == false) || (dis >= 12.0)) || ((skip != 0) && (skip != 4))) {
        return 0.0;}
    // Calculate the necessary parameters. CCELEC is a constant defined in the
    // Proteins.h header file
    float Dij = pow(dis, 2) / (4.0 * gbRadius * other->gbRadius);
    float numerator = elec_charge * other->elec_charge;
    float denominator = sqrt(pow(dis, 2)+gbRadius*other->gbRadius*exp(-Dij));
    float energy = -CCELEC * (1.0-(1.0/80.0)) * numerator / denominator;
    return energy;
}

// Calculate the Lazaridis-Karplus implicit solvation energy between 2 atoms
float PROT::Atom::LK (const Atom * other, const coor dis, 
                      const int skip) const {
    // Return 0 if: LK energy should not be calculated, the distance is
    // greater than 9 angstroms (different cutoff from other energy terms),
    // skip has a value of 2 or 3, or either atom is a hydrogen 
    if ((((lk == false) || (dis >= 9.0)) || ((skip != 0) && (skip != 4))) ||
         ((lkR <= 1.34) || (other->lkR <= 1.34))) {return 0.0;}
    // Calculate the appropriate values
    float x01 = (dis - lkR)/lkL;
    float x02 = (dis - other->lkR)/other->lkL;
    float t1 = -0.08979/pow(dis, 2);
    float t2 = exp(-pow(x01,2)) * lkG * other->lkV / lkL;
    float t3 = exp(-pow(x02,2)) * other->lkG * lkV / other->lkL;
    return t1 * (t2 + t3);
}

// Calculate the energy between this atom and another
float PROT::Atom::calculate_energy (const Atom * other) const {
    // Calculate the distance between the Atoms
    PROT::coor dis = calculate_distance (other);
    // Calculate the exclusion information between the Atoms
    int skip = exclusions (other);
    // If energy calculations should not be done with this pair of atoms
    if (skip < 0) {return 0.0;}
    // Calculate the energy terms
    float vdw = VDW (other, dis, skip);
    float elec = ELEC (other, dis, skip);
    float solv = GB(other, dis, skip) + LK(other, dis, skip);
    // The final energy will be stored here
    float energy = 0.0;
    // If this is a buried residue, use one set of weights (weights were
    // determined by Matthew Grisewood using machine learning)
    if (m_core) {energy = 4.0*vdw + 2.0*elec + 16.0*solv;}
    // And if it is a surface residue, use another set
    else {energy = 3.0*vdw + 1.0*elec + 9.0*solv;}
    return energy;
}

// Determine whether or not two Atoms have a VDW overlap
bool PROT::Atom::check_vdw_overlap (const Atom * other) const {
    // Calculate the distance between them
    PROT::coor dis = calculate_distance (other);
    // If that distance is < 80% of the sum of their VDW radii
    return (dis < 0.8*(vdwR + other->vdwR));
}

// Count the number of VDW clashes between 2 atoms (0 or 1)
size_t PROT::Atom::calculate_vdw_overlap (const Atom * other) const {
    // If either Atom is a hydrogen, there is no overlap
    if ((m_name[0] == 'H') || (other->m_name[0] == 'H')) {return 0;}
    // Check if there is an overlap
    if (check_vdw_overlap (other)) {return 1;}
    return 0;
}

// End the header guard if-statement from the start of the file
#endif
