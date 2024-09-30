/* Written in 2018 by Dr. Robert Pantazes. Copyright owned by Auburn
 * University. 
 * This file contains the implementation of the Matrix class. */

// Use a header guard to prevent the file from being included in a compiled
// program multiple times.
#ifndef Matrix_Guard
#define Matrix_Guard 1

// This header file is intended to be included in a compiled C++ program by
// the Proteins.h header file and does not need to include any other files
// itself.

// The initialize method of the Matrix class assigns 0s to the attributes
void PROT::Matrix::initialize () {
    m_rows = 0; m_columns = 0; m_values = 0; return;
}

// The clean-up function makes sure that the dynamically allocated memory is
// deleted
void PROT::Matrix::clean_up () {
    if (m_values != 0) {delete[] m_values; m_values = 0;}
    return;
}

// Copy the information from another Matrix instance
void PROT::Matrix::copy (const Matrix * other) {
    // delete existing information
    clean_up();
    // Get the number of rows and columns from the other matrix
    m_rows = other->m_rows;
    m_columns = other->m_columns;
    // Determine how many entries there are
    size_t n = m_rows * m_columns;
    // If there are values
    if (n > 0) {
        m_values = new coor [n];
        // Copy the values
        for(size_t i=0; i<n; ++i) {m_values[i] = other->m_values[i];}}
    return;
}

// Error check row and column indexes for a matrix
void PROT::Matrix::error_check (const size_t r, const size_t c) const {
    if ((r < 0) || (r >= m_rows)) {
        stringstream c1; c1 << m_rows;
        stringstream c2; c2 << r;
        string error = c2.str() + " is not a valid row index in a Matrix with "
                     + c1.str() + " rows.\n";
        throw logic_error (error);}
    else if ((c < 0) || (c >= m_columns)) {
        stringstream c1; c1 << c;
        stringstream c2; c2 << m_columns;
        string error = c1.str() + " is not a valid column index in a Matrix "
                       "with " + c2.str() + " columns.\n";
        throw logic_error (error);}
    return;
}

size_t PROT::Matrix::index (const size_t r, const size_t c) const {
    // Calculate the value of a linear index from a row, column pair
    return (r * m_columns) + c;
}

// Allocate space for the matrix to be r x c
void PROT::Matrix::allocate (const size_t r, const size_t c) {
    // Error check the number of rows and columns
    if (r < 1) {
        stringstream c1; c1 << r;
        string error = "A Matrix must have a positive number of rows.\n"
                     + c1.str() + " is not acceptable.\n";
        throw logic_error (error);}
    else if (c < 1) {
        stringstream c1; c1 << c;
        string error = "A Matrix must have a positive number of columns.\n"
                     + c1.str() + " is not acceptable.\n";
        throw logic_error (error);}
    // Delete current information
    clean_up();
    // Assign the values and set each entry to 0
    m_rows = r;
    m_columns = c;
    size_t n = m_rows * m_columns;
    m_values = new coor [n];
    for(size_t i=0; i<n; ++i) {m_values[i] = 0;}
}

// Set the Matrix up from an Atom
void PROT::Matrix::allocate (const Atom * atom) {
    // Clean up existing information
    clean_up();
    // Set the Matrix up to be a 1 x 3 matrix of the Atom's coordinates
    m_rows = 1;
    m_columns = AtomCoordinates;
    size_t n = m_rows * m_columns;
    m_values = new coor [n];
    for(size_t i=0; i<n; ++i) {m_values[i] = atom->m_coors[i];}
}

// Set up a Matrix for a rotation calculation using Rodriguez's rotation
// formula
void PROT::Matrix::allocate (const coor angle, const coor vector []) {
    // Set up the Matrix to have appropriate dimensions
    clean_up();
    m_rows = AtomCoordinates;
    m_columns = AtomCoordinates;
    m_values = new coor [m_rows * m_columns];
    // Calculate some trig values that are needed
    coor c = cos(angle);
    coor s = sin(angle);
    coor v = 1 - c;
    // Store values in the matrix
    m_values[0] = c + v * vector[0] * vector[0];
    m_values[1] = -s*vector[2] + v * vector[0] * vector[1];
    m_values[2] = s*vector[1] + v * vector[0] * vector[2];
    m_values[3] = s*vector[2] + v * vector[1] * vector[0];
    m_values[4] = c + v * vector[1] * vector[1];
    m_values[5] = -s*vector[0] + v * vector[1] * vector[2];
    m_values[6] = -s*vector[1] + v * vector[2] * vector[0];
    m_values[7] = s*vector[0] + v * vector[2] * vector[1];
    m_values[8] = c + v * vector[2] * vector[2];
    return;
}

// Set up a Matrix to contain all of the Atoms in a Residue
void PROT::Matrix::allocate (const Residue * res) {
    // Confirm that the residue has atoms
    if (res->m_count == 0) {
        string error = "A Matrix cannot be allocated from an empty Residue.\n";
        throw logic_error (error);}
    // Clean up existing values
    clean_up();
    // The number of rows is the number of atoms
    m_rows = res->m_count;
    // The number of columns is the number of coordinates per atom
    m_columns = AtomCoordinates;
    // Allocate memory
    m_values = new coor [m_rows * m_columns];
    // fill in the values
    for(size_t i=0; i<res->m_count; ++i) {
        for(size_t j=0; j<AtomCoordinates; ++j) {
            m_values[i*m_columns + j] = res->m_atoms[i].m_coors[j];}}
    return;
}

// Set up a Matrix using a Protein
void PROT::Matrix::allocate (const Protein * prot, 
                             const bool backboneOnly = false) {
    // Confirm that the protein has residues
    if ((prot == 0) || (prot->number_of_atoms() == 0)) {
        string error = "A Matrix cannot be allocated from an empty Protein.\n";
        throw logic_error (error);}
    // Create a vector of pointers to the relevant Protein Atoms
    vector<Atom *> ptrs = prot->select_atoms(backboneOnly);
    // If there are no Atoms, throw an error
    if (ptrs.size() == 0) {
        string error = "No Atoms were identified for allocating the Matrix.\n";
        throw logic_error (error);}
    // Delete existing information
    clean_up();
    // Set the number of rows and columns
    m_rows = ptrs.size();
    m_columns = AtomCoordinates;
    // Allocate memory
    m_values = new coor [m_rows * m_columns];
    // Loop through the Atoms
    for(size_t i=0; i<ptrs.size(); ++i) {
        // Loop through the coordinates
        for(size_t j=0; j<AtomCoordinates; ++j) {
            // Store the coordinate of the Atom in the matrix
            m_values[i*m_columns + j] = ptrs[i]->m_coors[j];}}
    return;
}

// Set up a Matrix using a vector of Atoms
void PROT::Matrix::allocate (const vector<Atom *>& atoms) {
    // Delete existing memory
    clean_up();
    // If the vector of Atoms is empty, raise an error
    if (atoms.size() == 0) {
        string error = "A Matrix cannot be allocated from an empty vector of "
                       "Atoms.\n";
        throw logic_error (error);}
    // Set the rows and columns appropriately
    m_rows = atoms.size(); m_columns = AtomCoordinates;
    // Allocate memory
    size_t n = m_rows * m_columns;
    m_values = new coor [n];
    // Set the values
    for(size_t i=0; i<m_rows; ++i) {
        for(size_t j=0; j<AtomCoordinates; ++j) {
            m_values[i*AtomCoordinates + j] = atoms[i]->m_coors[j];}}
    return;
}

// Construct a Matrix from a Protein
PROT::Matrix::Matrix (const Protein * prot, const bool backboneOnly = false) {
    initialize(); allocate(prot, backboneOnly);
    return;
}

// Access to the Matrix's information
size_t PROT::Matrix::rows () const {return m_rows;}
size_t PROT::Matrix::columns () const {return m_columns;}

PROT::coor PROT::Matrix::operator() (const size_t r, const size_t c) const {
    // Error check the indices
    error_check (r, c);
    // Get the proper linear index
    size_t n = index(r, c);
    // Return the value
    return m_values[n];
}

// Specify a value in a matrix
void PROT::Matrix::set (const size_t r, const size_t c, const coor v) {
    error_check (r, c);
    size_t n = index(r, c);
    m_values[n] = v;
    return;
}

// Add a value to an entry in a matrix
void PROT::Matrix::increase (const size_t r, const size_t c, const coor v) {
    error_check (r, c);
    size_t n = index(r, c);
    m_values[n] += v;
    return;
}

// Multiply two matrices together
PROT::Matrix PROT::Matrix::product (const Matrix * other) const {
    // Check that the inner dimensions of both matrices match
    if (m_columns != other->m_rows) {
        stringstream c1; c1 << m_rows;
        stringstream c2; c2 << m_columns;
        stringstream c3; c3 << other->m_rows;
        stringstream c4; c4 << other->m_columns;
        string error = "It is not possible to calculate the product of a "
                     + c1.str() + "x" + c2.str() + " matrix and a "
                     + c3.str() + "x" + c4.str() + " matrix.\n";
        throw logic_error (error);}
    // Create the new matrix
    Matrix output (m_rows, other->m_columns);

    // Calculate each entry of the new matrix
    size_t n = 0;
    for(size_t i=0; i<m_rows; ++i) {
        for(size_t j=0; j<other->m_columns; ++j) {
            for(size_t k=0; k<m_columns; ++k) {
                size_t N1 = index(i, k);
                size_t N2 = other->index(k, j);
                output.m_values[n] += m_values[N1] * other->m_values[N2];}
            ++n;}}
    return output;
}

// Calculate the transpose of the matrix
PROT::Matrix PROT::Matrix::transpose () const {
    // Create a new matrix with reversed values from this matrix
    Matrix output (m_columns, m_rows);
    // Assign the values
    for(size_t i=0; i<m_columns; ++i) {
        for(size_t j=0; j<m_rows; ++j) {
            output.m_values[output.index(i,j)] = m_values[index(j, i)];}}
    return output;
}

// Confirm that the Matrix can be used for moving an atom, residue or protein
void PROT::Matrix::move_check () const {
    if ((m_rows != 1) || (m_columns != AtomCoordinates)) {
        stringstream c1; c1 << AtomCoordinates;
        string error = "Only a 1x" + c1.str() + " Matrix can be used in "
                       "movement calculations.\n";
        throw logic_error (error);}
    return;
}

// Or for rotation
void PROT::Matrix::rotate_check () const {
    if ((m_rows != AtomCoordinates) || (m_columns != AtomCoordinates)) {
        stringstream c1; c1 << AtomCoordinates;
        string error = "Only a " + c1.str() + "x" + c1.str() + " Matrix can "
                       "be used in rotation calculations.\n";
        throw logic_error (error);}
    return;
}

// A string representation of the Matrix
string PROT::Matrix::str () const {
    // Store the output here
    string output = "";
    // If there are no values, be done
    if (m_values == 0) {return output;}
    // Loop through all positions
    for(size_t i=0; i<m_rows; ++i) {
        for(size_t j=0; j<m_columns; ++j) {
            // Use a stringstream to convert the value to a string
            stringstream c1; c1 << fixed <<setprecision(3) 
                                << m_values[i*m_columns + j];
            Text::rjust_insert(output, c1.str(), 12, ' ');}
        output += "\n";}
    return output;
}

// Convert a vector matrix into a unit vector
void PROT::Matrix::make_unit_vector () {
    // The matrix must have either 1 row or 1 column
    if ((m_rows != 1) && (m_columns != 1)) {
        stringstream c1; c1 << m_rows;
        stringstream c2; c2 << m_columns;
        string error = "A " + c1.str() + "x" + c2.str() + " matrix is not a "
                       "vector.\n";
        throw logic_error (error);}
    // Calculate the magnitude of the vector
    coor magnitude = 0.0;
    for(size_t i=0; i<(m_rows*m_columns); ++i) {
        magnitude += pow(m_values[i], 2);}
    magnitude = sqrt(magnitude);
    // Divide each value by the magnitude
    for(size_t i=0; i<(m_rows*m_columns); ++i) {
        m_values[i] = m_values[i]/magnitude;}
    return;
}

// Create a new matrix that is the difference between two existing matrices
PROT::Matrix PROT::Matrix::operator- (const Matrix& other) const {
    // Make sure the dimensions match
    if ((m_rows != other.m_rows) || (m_columns != other.m_columns)) {
        string error = "Matrix dimensions must match for subtraction to be a "
                       "valid operation.\n";
        throw logic_error (error);}
    // Create a new empty matrix
    Matrix output (m_rows, m_columns);
    // Loop through the full range of values
    for(size_t i=0; i<(m_rows*m_columns); ++i) {
        output.m_values[i] = m_values[i] - other.m_values[i];}
    return output;
}

// Add two matrices together
PROT::Matrix PROT::Matrix::operator+ (const Matrix& other) const {
    if ((m_rows != other.m_rows) || (m_columns != other.m_columns)) {
        string error = "Matrix dimensions must match for addition to be a "
                       "valid operation.\n";
        throw logic_error (error);}
    Matrix output (m_rows, m_columns);
    for(size_t i=0; i<(m_rows*m_columns); ++i) {
        output.m_values[i] = m_values[i] + other.m_values[i];}
    return output;
}

// Calculate the cross product between two vector matrices
PROT::Matrix PROT::Matrix::cross (const Matrix& O) const {
    // Use a flag to check the dimensions of both matrices
    bool flag1 = false; bool flag2 = false;
    if (((m_rows == 1) && (m_columns == AtomCoordinates)) ||
        ((m_rows == AtomCoordinates) && (m_columns == 1))) {
        flag1 = true;}
    if (((O.m_rows == 1) && (O.m_columns == AtomCoordinates)) ||
        ((O.m_rows == AtomCoordinates) && (O.m_columns == 1))) {
        flag2 = true;}
    // If either flag is false, the matrix dimensions are not currently
    // permitted
    if ((!flag1) || (!flag2)) {
        stringstream c1; c1 << AtomCoordinates;
        string error = "Currently, cross products can only be calculated "
                       "between 1x" + c1.str() + " vectors.\n";
        throw logic_error (error);}
    // Create a new 1 x AtomCoordinates matrix
    Matrix output (1, AtomCoordinates);
    // Set its values appropriately
    output.m_values[0] = m_values[1]*O.m_values[2] - m_values[2]*O.m_values[1];
    output.m_values[1] = m_values[2]*O.m_values[0] - m_values[0]*O.m_values[2];
    output.m_values[2] = m_values[0]*O.m_values[1] - m_values[1]*O.m_values[0];
    return output;
}

// End the header guard from the start of the file
#endif
