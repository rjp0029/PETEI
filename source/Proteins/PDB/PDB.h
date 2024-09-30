/* Written in 2018 by Dr. Robert Pantazes. Copyright owned by Auburn
 * University.
 * This file contains the definitions of the methods of the PDB class, a class
 * intended to process the contents of PDB files */

// Use a header guard to ensure this file is only loaded into a compiled
// program a single time
#ifndef PDB_Header_Guard
#define PDB_Header_Guard 1

// This file is intended to be included directly by the Proteins.h header
// file. As such, it does not include any additional C++ dependencies itself

// Copy the information from another instance of the PDB class
void PROT::PDB::copy (const PDB * other) {
    // Copy the non-vector variables
    m_name = other->m_name;
    m_folder = other->m_folder;
    m_type = other->m_type;
    m_resolution = other->m_resolution;
    m_obsolete = other->m_obsolete;
    // Copy the lines vector
    m_lines.clear(); m_lines.reserve(other->m_lines.size());
    if (other->m_lines.size() > 0) {
        for(size_t i=0; i<other->m_lines.size(); ++i) {
            m_lines.push_back(other->m_lines[i]);}}
    // Copy the proteins
    m_proteins.clear(); m_proteins.reserve(other->m_proteins.size());
    if (other->m_proteins.size() > 0) {
        for(size_t i=0; i<other->m_proteins.size(); ++i) {
            m_proteins.push_back(other->m_proteins[i]);}}
    // Copy the Structures, making sure the pointers will point to the
    // Proteins in this PDB file
    m_structures.clear(); m_structures.reserve(other->m_structures.size());
    if (other->m_structures.size() > 0) {
        for(size_t i=0; i<other->m_structures.size(); ++i) {
            m_structures.push_back(Structure(other->m_structures[i], m_proteins));}}
    return;
}

// Load the contents of the PDB file into this file
void PROT::PDB::load () {
    // Assemble the name of the file 
    string fileName = m_folder + m_name;
    // Attempt to open that file
    ifstream input; input.open(fileName.c_str());
    // If the file is not open, throw an error
    if(!input.is_open()) {
        string error = "Failure to open:\nFile: " + m_name + "\nLocation: ";
        if ((m_folder == "./") || (m_folder == "")) {
            error += "Current Folder\n";}
        else {error += m_folder + "\n";}
        throw logic_error (error);}
    // The lines will be stored here
    string line;
    // For NMR files, only lines containing the first model will be loaded
    bool model_flag = false;
    // Go through the file's lines
    getline(input, line);
    while(!input.eof()) {
        // Strip whitespace from the line
        Text::strip(line);
        // If the line starts with the word "model"
        if (Text::startswith(line, "MODEL")) {
            // Split the line into pieces
            vector<string> parts = Text::split(line);
            // If there are at least two pieces and the second piece is not
            // "1", set the flag to true. If it is 1, set the flag to false
            if (parts.size() > 1) {
                if (parts[1] == "1") {model_flag = false;}
                else {model_flag = true;}}}
        // If the line should be stored, do so
        if(!model_flag) {m_lines.push_back(line);}
        // Get the next line
        getline(input, line);}
    // Close the input file
    input.close();
    // If no contents were identified, raise an error
    if (m_lines.size() == 0) {
        string error = "No contents were identified in:\nFile: " + m_name
                     + "\nLocation: " + m_folder + "\n";
        throw logic_error (error);}
    // End this function
    return;
}

// Search through a PDB file's information to identify the type of experiment
// that generated the data
void PROT::PDB::identify_type () {
    // Set the type of experiment to a default value
    m_type = "UNKNOWN";
    // This flag tells the program when to stop searching through the lines
    bool seen = false;
    // Loop through the lines
    for(size_t i=0; i<m_lines.size(); ++i) {
        // If the line starts with the appropriate phrase
        if (Text::startswith(m_lines[i], "EXPDTA")) {
            // Extract the phrase, going from the 11th character to the end of
            // the string
            string part = m_lines[i].substr(10);
            // If this is the first time the line has been seen
            if (!seen) {
                // Set the flag to true
                seen = true;
                // Set the type to the part
                m_type = part;}
            // If there are multiple lines of type information, concatenate
            // them appropriately
            else {
                // stick a semi-colon at the end of the current type
                // information if there is not one already
                if (m_type[m_type.size()-1] != ';') {m_type += ";";}
                m_type += " " + part;}}
        // If the experiment data type has already been stored, be done
        else if (seen) {break;}}
    return;
}

// Check to see whether or not a PDB file is obsolete. If it is, we aren't
// interested in using its information
void PROT::PDB::check_obsolete () {
    // Set the obsolete value to false
    m_obsolete = false;
    // Search through the contents of the file
    for(size_t i=0; i<m_lines.size(); ++i) {
        // If the line indicates that the file is obsolete
        if (Text::startswith(m_lines[i], "OBSLTE")) {
            m_obsolete = true; return;}
        // If the line is a REMARK or ATOM line, the search can be finished
        // because OBSLTE entries are required to come BEFORE that information
        else if ((Text::startswith(m_lines[i], "REMARK")) ||
                 (Text::startswith(m_lines[i], "ATOM"))) {break;}}
    // If the function has reached this point, check to see if it is a
    // Theoretical model. If so, it also shouldn't be used
    //m_obsolete = Text::contains(m_type, "THEORETICAL");
    return;
}

// Search through the PDB file's information for the resolution information
// for the experiment
void PROT::PDB::identify_resolution () {
    // Set the value to one that indicates no resolution information is
    // available
    m_resolution = -1.0;
    // Loop through the file's lines
    for(size_t i=0; i<m_lines.size(); ++i) {
        // If the line starts with the appropriate phrase
        if (Text::startswith(m_lines[i], "REMARK   2")) {
            // Split the line into pieces
            vector<string> parts = Text::split(m_lines[i]);
            // If there are the proper number of them
            if (parts.size() == 5) {
                // If the last part indicates a measurement
                if ((parts[4] == "ANGSTROMS.") || (parts[4] == "ANGSTROMS")) {
                    stringstream c1; c1 << parts[3]; c1 >> m_resolution;}
                // Stop the search
                break;}}}
    return;
}

// Use the SEQRES information to set up initial vectors of Residues
void PROT::PDB::initialize_Residues (vector<vector<Residue> >& residues) {
    // Set up the residues so that the vector can directly access contents via
    // character indexing
    residues.clear(); residues.resize(128);
    // This flag indicates when the necessary lines have been observed
    bool seen = false;
    // Loop through the contents of the file
    for(size_t i=0; i<m_lines.size(); ++i) {
        // If the line startswith SEQRES
        if (Text::startswith(m_lines[i], "SEQRES")) {
            seen = true;
            // The chain must be in character 12
            char L = m_lines[i][11];
            // Extract the substring containing the residues (positions 20 -
            // end of string)
            string res_part = m_lines[i].substr(19);
            // strip whitespace
            Text::strip(res_part);
            // Split into pieces
            vector<string> parts = Text::split(res_part);
            // Create and store each of them as a Residue
            for(size_t j=0; j<parts.size(); ++j) {
                residues[(int) L].push_back(Residue(parts[j], L));}}
        // If the sequence residue information has been found, be done
        else if (seen) {break;}}
    return;
}

// Use the REMARK 465 lines to identify missing residues
void PROT::PDB::identify_missing_residues (vector<vector<Residue> >& missing) {
    // Set up the missing vector to be accessible by character indexing
    missing.clear(); missing.resize(128);
    // This flag indicates when the relevant lines have been seen
    bool seen = false;
    // Store the appropriate lines in this vector
    vector<string> remarks;
    // Loop through the contents of the file
    for(size_t i=0; i<m_lines.size(); ++i) {
        // If the line starts with REMARK 465. I considered checking to see if
        // it was a REMARK line and then extracting the number, but this
        // SHOULD be more efficient in most cases.
        if (Text::startswith(m_lines[i], "REMARK 465")) {
            // Indicate that these lines have been seen
            seen = true;
            // If the line is long enough, store it with the header removed
            if (m_lines[i].size() > 10) {
                remarks.push_back(m_lines[i].substr(10));}}
        // Stop searching for lines once the REMARK 465 lines have been
        // collected or an ATOM line is seen
        else if (seen) {break;}
        else if (Text::startswith(m_lines[i], "ATOM")) {break;}}
    // If there were no REMARK 465 lines, be done
    if (remarks.size() == 0) {return;}
    // The standards for the REMARK 465 and 470 lines changed over time.
    // However, there is always an initial set of remarks and then a header
    // line that properly positions the contents. Find that header line
    size_t H = remarks.size();
    for(size_t i=0; i<remarks.size(); ++i) {
        if (Text::endswith(remarks[i], "RES C SSSEQI")) {
            H = i; break;}}
    // If that header line was not found, throw an error
    if (H == remarks.size()) {
        string text = "Failure to find the start of the Missing Residues.\n";
        throw logic_error (text);}
    // Since the ending string is always the same, use math to determine the
    // appropriate columns for the information
    size_t n = remarks[H].size();
    // The column with the insertion character is n - 1
    size_t I = n-1;
    // The residue number can start 5 positions before that
    size_t N = I-5;
    // The chain is two positions before that
    size_t C = N-2;
    // And the name of the residue starts 4 positions before that
    size_t R = C-4;
    // Loop through the subsequent lines
    if (H == remarks.size() - 1) {return;}
    for (size_t i=H+1; i<remarks.size(); ++i) {
        // Get the name of the residues
        string resName = remarks[i].substr(R, 3); Text::strip(resName);
        // Get the name of the protein
        char L = remarks[i][C];
        // Get the number of the residue
        string numStr = remarks[i].substr(N, 5); Text::strip(numStr);
        long resNum = 0;
        stringstream c1; c1 << numStr; c1 >> resNum;
        // If there are enough characters, get the insertion code
        char insert = ' ';
        if (remarks[i].size() > I) {insert = remarks[i][I];}
        // Create a Residue using this information
        int j = (int) L;
        missing[j].push_back(Residue(resName, L));
        missing[j][missing[j].size()-1].m_number = resNum;
        missing[j][missing[j].size()-1].m_insertion = insert;
        missing[j][missing[j].size()-1].m_present = false;}
    // End the function
    return;
}

// Collect the Atoms that make up the Proteins
void PROT::PDB::collect_Atoms (vector<vector<Atom> >& atoms) {
    // Set up the atoms vector to be accessible by character indexing
    atoms.clear(); atoms.resize(128);
    // Loop through the contents of the file
    for(size_t i=0; i<m_lines.size(); ++i) {
        // If the line starts with ATOM or HETATM
        if ((Text::startswith(m_lines[i], "ATOM")) ||
            (Text::startswith(m_lines[i], "HETATM"))) {
            // Try to create an atom from the line
            try {
                Atom atom (m_lines[i]);
                // Convert the Atom's character to an integer
                int n = (int) atom.m_protein;
                // If the Atom's alternative location characteristic is 'A' or
                // ' ', store it
                if (atom.m_alt == 'A') {
                    atom.m_alt = ' '; atoms[n].push_back(atom);}
                else if (atom.m_alt == ' ') {
                    atoms[n].push_back(atom);}}
            // If an error occurs, just ignore it
            catch (logic_error) {continue;}}}
    return;
}

// The identify missing Atoms function is run after the Proteins have been
// assembled. It uses the REMARK 470 lines to store how many atoms are missing
// from various residues
void PROT::PDB::identify_missing_atoms () {
    // This flag indicates whether or not REMARK 470 lines have been seen
    bool seen = false;
    // Store the lines that list the missing Atoms here
    vector<string> remarks;
    // Loop through the contents of the file
    for(size_t i=0; i<m_lines.size(); ++i) {
        // If this line is a REMARK 470 line
        if (Text::startswith(m_lines[i], "REMARK 470")) {
            // Indicate that the lines have been seen
            seen = true;
            // If the line is long enough, store it with the REMARK 470 part
            // removed
            if (m_lines[i].size() > 10) {
                remarks.push_back(m_lines[i].substr(10));}}
        // If the lines have been seen and this wasn't one, stop searching for
        // them
        else if (seen) {break;}
        // Or if an ATOM line is reached
        else if (Text::startswith(m_lines[i], "ATOM")) {break;}}
    // If there are no such remarks, end this function
    if (remarks.size() == 0) {return;}
    // The format for this section changed around 2007. The following code
    // should work correctly regardless of the format, but it takes a little
    // bit of work to correctly identify where the information is located.
    // First, find the line that ends the header information and describes
    // where the subsequent values are stored
    size_t H = remarks.size();
    for(size_t i=0; i<remarks.size(); ++i) {
        // In all cases, the line should end with CSSEQI  ATOMS
        if (Text::endswith(remarks[i], "CSSEQI  ATOMS")) {
            H = i; break;}}
    // If an appropriate line was not identified, throw an error
    if (H >= remarks.size() - 1) {
        string error = "Failure to identify missing atoms in REMARK 470 lines.\n";
        throw logic_error (error);}
    // Get the length of that line
    size_t n = remarks[H].size();
    // The listing of the Atoms starts 4 positions before that
    size_t A = n - 4;
    // The insertion character is 4 positions before that
    size_t I = A - 4;
    // The residue number starts 4 positions before that
    size_t N = I - 4;
    // The chain starts 1 position before that
    size_t C = N - 1;
    // And the Residue starts either 4 or 5 positions before that
    size_t R = C - 5;
    if (remarks[H][C-4] == 'R') {R = C - 4;}
    // Loop through the lines listing missing Atoms
    for(size_t i=H+1; i<remarks.size(); ++i) {
        // Identify the residue's name
        string resName = remarks[i].substr(R, 3); Text::strip(resName);
        // Identify the Protein's name
        char L = remarks[i][C];
        // Extract the Residue's number's string
        string numStr = remarks[i].substr(N, 4); Text::strip(numStr);
        long resNum = 0;
        stringstream c1; c1 << numStr; c1 >> resNum;
        // Get the insertion character
        char insert = remarks[i][I];
        // Get the string of Atoms
        string atomStr = remarks[i].substr(A);
        // Split it into pieces
        vector<string> pieces = Text::split(atomStr);
        // A flag to indicate that the information was used
        bool flag = false;
        // Search the proteins to find the residue
        for(size_t j=0; j<m_proteins.size(); ++j) {
            // If this isn't the right protein, continue the search
            if (m_proteins[j].m_name != L) {continue;}
            // Check the protein's residues
            for(size_t k=0; k<m_proteins[j].m_count; ++k) {
                // Get a pointer to the Residue
                Residue * res = &(m_proteins[j].m_residues[k]);
                // If the Residue's information matches
                if ((res->m_name == resName) &&
                    ((res->m_number == resNum) &&
                     (res->m_insertion == insert))) {
                    // Indicate that the information is being used
                    flag = true;
                    // Add the number of missing atoms to the residue's
                    // information
                    res->m_missing_atoms += pieces.size();
                    break;}}
            // If the for loop reaches this point, it means it checked the
            // protein with the matching name. Break out of the search
            break;}
        // If the information was not used, throw an error
        if (!flag) {
            string error = "Failure to identify the Residue missing these "
                           "Atom:\n" + remarks[i] + "\n";
            throw logic_error(error);}}
    // End the function
    return;
}

// Integrate Missing Residues into the list of Residues that make up a
// Protein. Then store that Protein in the PDB file's internal vector
void PROT::PDB::integrate_missing_residues (Protein& prot,
                                            vector<Residue>& missing) {
    // If there are no missing Residues, just store the protein as it is
    if (missing.size() == 0) {m_proteins.push_back(prot); return;}
    // Since there are missing Residues, extract the Residues from the protein
    // and store it as a vector
    vector<Residue> residues; residues.reserve(prot.m_count + missing.size());
    for(size_t i=0; i<prot.m_count; ++i) {
        residues.push_back(prot.m_residues[i]);}
    // Loop through each missing residue
    for(size_t i=0; i<missing.size(); ++i) {
        // A pointer to the missing residue
        Residue * m = &(missing[i]);
        // The vector of residues needs to be searched for the best place to
        // insert this missing residue. Store the best insertion score and the
        // index of the residue this one should be placed before here
        size_t best_score = residues.size() * 100000;
        size_t best_index = residues.size();
        // A pointer to a "Known" residue, starting with the first
        Residue * k = &(residues[0]);
        // Determine if the missing residue goes before the first residue in
        // the protein
        if ((m->m_number < k->m_number) || 
           ((m->m_number == k->m_number) && (m->m_insertion < k->m_insertion))){
            // If it does, the score is guaranteed to be better than the best
            // score, so just save that information
            best_score = (k->m_number - m->m_number) + 1;
            best_index = 0;}
        // Check to see if the missing residue inserts well into other points
        // in the vector of Residues
        for(size_t j=1; j<residues.size(); ++j) {
            // The known residue is this one
            k = &(residues[j]);
            // And a previous residue is the one before it
            Residue * p = &(residues[j-1]);
            // The missing residue must come after the previous residue, so if
            // it doesn't continue the search
            if ((m->m_number < p->m_number) ||
               ((m->m_number == p->m_number) && 
                (m->m_insertion < p->m_insertion))) {continue;}
            // The missing residue must come before the current residue
            if ((k->m_number < m->m_number) ||
               ((k->m_number == m->m_number) && 
                (k->m_insertion < m->m_insertion))) {continue;}
            // Since the missing residue comes after the previous residue and
            // before the current residue, calculate a score for inserting it
            // here
            // This calculation is weird because I've been able to identify at
            // least one case where units of thousands are used to indicate
            // insertions instead of an insertion character
            size_t score;
            if (k->m_number >= p->m_number) {score = k->m_number - p->m_number;}
            else {score = p->m_number - k->m_number;}
            // If this is the best score found so far
            if (score < best_score) {
                best_score = score;
                best_index = j;}}
        // Calculate a score for inserting the residue at the end of the
        // vector
        k = &(residues[residues.size()-1]);
        if ((m->m_number > k->m_number) ||
           ((m->m_number == k->m_number) && (m->m_insertion > k->m_insertion))){
            size_t score = (m->m_number - k->m_number) + 1;
            if (score < best_score) {
                best_score = score;
                best_index = residues.size();}}
        // Insert the missing residue into the best position identified for
        // it. If that is the end of the vector, just use push_back to store
        // it
        if (best_index == residues.size()) {residues.push_back(missing[i]);}
        // Otherwise, insert the residue at the proper position. To do that, I
        // need an iterator to the start of the vector
        else {
            vector<Residue>::iterator it = residues.begin();
            // Insert the residue
            residues.insert(it + best_index, missing[i]);}}
    // Store the residues as a Protein in the PDB file
    m_proteins.push_back(Protein(residues));
    // End the function
    return;
}

// Validate that identified SEQRES information matches the specified Protein's
// information
void PROT::PDB::validate_seqres (Protein * prot, vector<Residue>& seqres) {
    // There must be at least as many residues in the protein as there are in
    // the sequence residue data
    if (prot->m_count < seqres.size()) {
        string error = "The SEQRES data contains more residues than are in the "
                       "ATOM and Missing Residue lines for chain ";
        error += prot->m_name;
        throw logic_error (error + ".\n");}
    // Create a vector of boolean values indicating whether or not the a
    // protein residue should be compared to the SEQRES data
    vector<bool> use; use.resize(prot->m_count);
    // Go through every protein residue
    for(size_t i=0; i<prot->m_count; ++i) {
        // Initially set the value to false (don't use)
        use[i] = false;
        // As a short-cut for HETATM entries like water
        if ((i > 0) && 
            (prot->m_residues[i].m_name == prot->m_residues[i-1].m_name)) {
            use[i] = use[i-1]; continue;}
        // Search the SEQRES data
        for(size_t j=0; j<seqres.size(); ++j) {
            // If there is at least one SEQRES entry with the same Residue
            // name, use this residue in the comparison
            if (seqres[j].m_name == prot->m_residues[i].m_name) {
                use[i] = true; break;}}}
    // Count up the number of usable residues
    size_t count = 0;
    for(size_t i=0; i<use.size(); ++i) {if(use[i]) {++count;}}
    // This count MUST match the number of SEQRES entries
    if (count != seqres.size()) {
        stringstream c1; c1 << seqres.size();
        stringstream c2; c2 << count;
        string error = "Inconsistency between the SEQRES and Sequence "
                       "information for chain ";
        error += prot->m_name;
        error += ".\nThe SEQRES data lists " + c1.str() + " residues "
                 "but " + c2.str() + " residues have sequence-eligible "
                 "structure information.\n";
        throw logic_error (error);}
    // An index for which residue to use in the Protein
    size_t PI = 0;
    // Loop through the SEQRES values
    for(size_t i=0; i<seqres.size(); ++i) {
        // Find the next protein residue to use
        while (!use[PI]) {++PI;}
        // If the residue information does not match, raise an error
        if (prot->m_residues[PI].m_name != seqres[i].m_name) {
            stringstream c1; c1 << i+1;
            stringstream c2; c2 << PI+1;
            stringstream c3; c3 << prot->m_residues[PI].m_number;
            string resNum = c3.str();
            resNum += prot->m_residues[PI].m_insertion;
            Text::strip(resNum);
            string error = "A Mismatch occurred between the SEQRES and "
                           "identified sequence information in chain ";
            error += prot->m_name;
            error += ":\nThis happend at residue " + c1.str() + " in the "
                     "SEQRES data and residue " + c2.str() + " (" + resNum
                   + ") in the sequence information.\n";
            throw logic_error (error);}
        // Increment the protein index
        ++PI;}
    return;
}

// Construct the Proteins from the information in the PDB file
void PROT::PDB::construct_Proteins () {
    // Collect the information on the SEQRES data
    vector<vector<Residue> > seqres;
    initialize_Residues(seqres);
    // Collect the information on the missing Residues
    vector<vector<Residue> > missing;
    identify_missing_residues(missing);
    // Collect the Atoms that make up the Proteins
    vector<vector<Atom> > atoms;
    collect_Atoms (atoms);
    // Generate proteins with those Atoms
    vector<Protein> prots;
    for(size_t i=0; i<128; ++i) {
        if (atoms[i].size() > 0) {prots.push_back(Protein(atoms[i]));}}
    // If there are no proteins, raise an error
    if (prots.size() == 0) {
        string error = "No Proteins were identified.\n";
        throw logic_error (error);}
    // Loop through those proteins
    for(size_t P=0; P<prots.size(); ++P) {
        // Get the character index of this protein
        int i = (int) prots[P].m_name;
        // Integrate the missing residues into that protein and store it in
        // the PDB file
        integrate_missing_residues (prots[P], missing[i]);
        // Validate the SEQRES information for this protein if there is SEQRES
        // data
        if (seqres[i].size() > 0) {
            validate_seqres (&(m_proteins[m_proteins.size()-1]), seqres[i]);}}
    // Identify atoms that may be missing in those proteins
    identify_missing_atoms ();
    return;
}

// After the Proteins have been constructed, they can be grouped into
// Structures
void PROT::PDB::create_Structures () {
    // Identify the compound line statements and separate them based on the
    // molecule ID
    vector<string> compound;
    vector<vector<string> > compounds;
    bool seen = false;
    // Loop through the lines of the PDB file
    for(size_t i=0; i<m_lines.size(); ++i) {
        // If this is the right sort of line
        if (Text::startswith(m_lines[i], "COMPND")) {
            seen = true;
            // Extract the useful part of the line
            string phrase = m_lines[i].substr(10);
            Text::strip(phrase);
            // If this is a Molecule ID line and there is existing compound
            // information, store it
            if (Text::startswith(phrase, "MOL_ID:")) {
                if (compound.size() > 0) {
                    compounds.push_back(compound);
                    compound.clear();}}
            // Otherwise, store this line in the compound
            else {compound.push_back(phrase);}}
        else if (seen) {break;}}
    if (compound.size() > 0) {compounds.push_back(compound);}
    // If no compound lines were seen, end this function
    if (!seen) {return;}
    // Create a vector of boolean values indicating which Proteins have been
    // assigned to a Structure
    vector<bool> used; used.reserve(m_proteins.size());
    for(size_t i=0; i<m_proteins.size(); ++i) {used[i] = false;}
    // Loop through the compound information
    for(size_t C=0; C<compounds.size(); ++C) {
        // Store a structure in the structures vector
        m_structures.push_back(Structure());
        size_t n = m_structures.size() - 1;
        // Loop through the relevant lines
        for(size_t i=0; i<compounds[C].size(); ++i) {
            // I only care about some of the lines. For instance, if the line
            // lists one or more names of the structure
            if ((Text::startswith(compounds[C][i], "MOLECULE:")) ||
                (Text::startswith(compounds[C][i], "SYNONYM:"))) {
                // Extract the meaningful part of the line and strip it of
                // whitespace
                string text = compounds[C][i].substr(9);
                Text::strip(text);
                // Split the text on commas, which separate the names
                vector<string> parts = Text::split(text, ',');
                // Loop through the parts
                for(size_t j=0; j<parts.size(); ++j) {
                    // Strip the part of whitespace
                    Text::strip(parts[j]);
                    // If the last character is a semi-colon, delete it
                    if (parts[j][parts[j].size()-1] == ';') {
                        parts[j].erase(parts[j].size()-1, 1);}
                    // Store the name
                    m_structures[n].m_names.push_back(parts[j]);}}
            else if (Text::startswith(compounds[C][i], "CHAIN:")) {
                // Strip off the leading portion of the string and whitespace
                string text = compounds[C][i].substr(6);
                Text::strip(text);
                // Split it on commas
                vector<string>parts = Text::split(text, ',');
                // Loop through that list
                for(size_t j=0; j<parts.size(); ++j) {
                    // Strip whitespace from the string
                    Text::strip(parts[j]);
                    // Get the first character
                    char L = parts[j][0];
                    // Use this flag to indicate whether or not that protein
                    // is found
                    bool found = false;
                    for(size_t k=0; k<m_proteins.size(); ++k) {
                        if (m_proteins[k].m_name == L) {
                            found = true;
                            // If this protein has already been used, raise an
                            // error
                            if (used[k]) {
                                string error = "Chain ";
                                error += L;
                                error += "is used in more than one Molecule.\n";
                                throw logic_error (error);}
                            // Indicate that the protein is being used
                            used[k] = true;
                            // Store the pointer in the structure
                            m_structures[n].m_proteins.push_back(&(m_proteins[k]));
                            break;}}
                    // If the protein wasn't identified, raise an error
                    if (!found) {
                        string error = "Chain ";
                        error += L;
                        error += " is listed in the COMPND lines, but no Atoms "
                                 "were identified for it.\n";
                        throw logic_error (error);}}}}}
    // If there are any proteins that have not been used, store them in an
    // unamed structure
    bool flag = false;
    size_t n = m_structures.size();
    for(size_t i=0; i<m_proteins.size(); ++i) {
        if (!used[i]) {
            if (!flag) {
                flag = true;
                m_structures.push_back(Structure());}
            m_structures[n].m_proteins.push_back(&(m_proteins[i]));}}
    // End the function
    return;
}

// The constructor of the PDB class
PROT::PDB::PDB (const string& fileName, const string path = "") {
    // Set the name
    m_name = fileName;
    // Set the folder
    m_folder = path;
    // Load the contents of the file
    load();
    // Use a try statement to add to any error that comes up in the subsequent
    // steps
    try {
        identify_type();
        check_obsolete();
        if (m_obsolete) {return;}
        identify_resolution ();
        construct_Proteins ();
        create_Structures ();
        }
    catch (logic_error& e) {
        string error = e.what();
        error += "This occurred in PDB file " + m_name + "\n";
        throw logic_error (error);}
    return;
}

// Exterior access to a line in the PDB file
string PROT::PDB::line (const size_t i) const {
    if (i >= m_lines.size()) {
        stringstream c1; c1 << m_lines[i].size();
        stringstream c2; c2 << i;
        string error = c1.str() + " lines were loaded from PDB file "
                     + m_name + ". An index of " + c2.str() + " is not "
                       "acceptable to access one of them.\n";
        throw logic_error (error);}
    return m_lines[i];
}

// Access to a Protein in the PDB file
PROT::Protein * PROT::PDB::protein (const size_t i) {
    if (i >= m_proteins.size()) {
        stringstream c1; c1 << m_proteins[i].size();
        stringstream c2; c2 << i;
        string error = c1.str() + " Proteins were identified in PDB file "
                     + m_name + ". An index of " + c2.str() + " is not "
                       "acceptable to access one of them.\n";
        throw logic_error (error);}
    return &(m_proteins[i]);
}

// Access to a Structure in the PDB file
PROT::Structure * PROT::PDB::structure (const size_t i) {
    if (i >= m_structures.size()) {
        stringstream c1; c1 << m_structures[i].proteins();
        stringstream c2; c2 << i;
        string error = c1.str() + " Structures were identified in PDB file "
                     + m_name + ". An index of " + c2.str() + " is not "
                       "accepatble to access one of them.\n";
        throw logic_error (error);}
    return &(m_structures[i]);
}

// A string representation of the PDB file's information
string PROT::PDB::str () const {
    // Store the string here
    string output = "PDB File: " + m_name + "\n";
    // If it is obsolete, write that and be done
    if (m_obsolete) {
        output += "THIS FILE IS NOT APPROPRIATE FOR USE.\n";
        return output;}
    // List the experiment type
    output += "Experiment Type: " + m_type + "\n";
    // List the resolution
    output += "Resolution: ";
    if (m_resolution > 0) {
        stringstream c1; c1 << fixed << setprecision(3) << m_resolution;
        output += c1.str() + " Angstroms.\n";}
    else {output += " Not Applicable.\n";}
    // List the structures in the PDB file
    output += "\n";
    for(size_t i=0; i<m_structures.size(); ++i) {
        output += m_structures[i].str();}
    return output;
}

// End the header guard if statement from the start of the file
#endif
