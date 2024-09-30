/* Written in 2018 by Dr. Robert Pantazes. Copyright owned by Auburn
 * University.
 * This file contains classes that are used for designing binding proteins
 * using a database of protein pieces. */

// Use a header guard to prevent the file from being included in a compiled
// program multiple times
#ifndef Design_Guard
#define Design_Guard 1

// Include the Proteins.h header file, which includes several constants and
// the Atom class
#include "Proteins.h"

// Define a namespace for the classes in this file, which involve interactions
namespace INTER {

    // The Spot class contains coordinates in three dimensional space
    class Spot;
    // The Plane class contains the information that defines a plane
    class Plane;
    // The Interaction class contains the information about an Interaction and
    // can assess whether or not a set of Atoms are compatible with the
    // interaction
    class Interaction;
    // The Position class connects an Antigen / Epitope residue with a
    // specific Interaction
    class Position;
    // The solution class is a container of Positions that are known to be
    // compatible with one another
    class Solution;

    // End the namespace
}

// Define the Spot class
class INTER::Spot {

    // Let the Plane class access the Spot's values
    friend class INTER::Plane;

    // The information stored in the class is private
    private:
        // An array of 3 floating point numbers
        float m_points [3];

    // The public interface of the class
    public:
        // A default constructor
        Spot ();
        // The standard constructor
        Spot (const float, const float, const float);
        // Copy construction and assignment
        Spot (const Spot&);
        void operator= (const Spot&);
        // Access to the class values
        size_t size () const {return 3;}
        float operator[] (const size_t) const;
        // Calculate a Spot that is the difference between this Spot and an
        // Atom
        Spot calculate_movement (const PROT::Atom *) const;
        // Calculate the distance between a Spot and an Atom while accouting
        // for a movement
        float calculate_distance (const PROT::Atom *, const Spot *) const;
        float calculate_distance (const vector<float>&, const Spot *) const;
        string str () const;

    // End the Spot class definition
};

// Implement the methods of the Spot class
INTER::Spot::Spot () {for(size_t i=0; i<3; ++i) {m_points[i] = 0.0;}}

INTER::Spot::Spot (const float n1, const float n2, const float n3) {
    m_points[0] = n1; m_points[1] = n2; m_points[2] = n3;}

INTER::Spot::Spot (const Spot& other) {
    for(size_t i=0; i<3; ++i) {m_points[i] = other.m_points[i];}}

void INTER::Spot::operator= (const Spot& other) {
    for(size_t i=0; i<3; ++i) {m_points[i] = other.m_points[i];}}

// Access to the coordinates in the spot
float INTER::Spot::operator[] (const size_t i) const {
    // Throw an error if the index isn't acceptable
    if (i >= 3) {
        stringstream c1; c1 << i;
        string error = "Spots have 3 coordinates, so an index of " + c1.str()
                     + " is not acceptable\n";
        throw logic_error (error);}
    return m_points[i];
}

// Create a new Spot with the difference between the Atom and this spot
// included
INTER::Spot INTER::Spot::calculate_movement (const PROT::Atom * atom) const {
    // Create the new Spot
    Spot answer = Spot ();
    // Store the values in the Spot
    for(size_t i=0; i<3; ++i) {
        answer.m_points[i] = atom->operator[](i) - m_points[i];}
    return answer;
}

// Calculate the distance from the Spot to an Atom while accoutning for a
// movement
float INTER::Spot::calculate_distance (const PROT::Atom * atom, 
                                       const Spot * move) const {
    // Store the distance here
    float dis = 0.0;
    // Go through the coordinates
    for(size_t i=0; i<3; ++i) {
        // Subtract the movement from the Atom's coordinates
        float d = atom->operator[](i) - move->m_points[i];
        // Add the square of that value minus the value in the Spot
        dis += pow(d - m_points[i], 2);}
    // Take the square root of the distance information
    dis = sqrt(dis);
    return dis;
}

// Calculate the distance from the Spot to a vector of three coordinates
// accounting for a movement
float INTER::Spot::calculate_distance (const vector<float>& coors,
                                       const Spot * move) const {
    // Duplicate what was done in the previous function
    float dis = 0.0;
    for(size_t i=0; i<3; ++i) {
        float d = coors[i] - move->m_points[i];
        dis += pow(d - m_points[i], 2);}
    dis = sqrt(dis);
    return dis;
}

// Create a string summarizing the spot's information
string INTER::Spot::str () const {
    // Store the result here
    string answer = "";
    // Go through the values
    for (size_t i=0; i<3; ++i) {
        // Use a stringstream with proper precision
        stringstream c; c << fixed << setprecision(3) << m_points[i];
        if (i > 0) {answer += " ";}
        answer += c.str();}
    return answer;
}

// Define the Plane class
class INTER::Plane {

    // The information stored in the class is private
    private:
        // A plane in 3 dimensional space has 4 defining values
        float m_values [4];
        // And the distance calculations are facilitate by pre-calculating the
        // magnitude of the plane's coefficients
        float m_mag;

    // A private function to copy the information from another instance of the
    // class
    private:
        void copy (const Plane& other) {m_mag = other.m_mag;
            for(size_t i=0; i<4; ++i) {m_values[i] = other.m_values[i];}}

    // The public interface of the class
    public:
        // The default constructor
        Plane ();
        // The standard constructor
        Plane (const float, const float, const float, const float);
        // Copy construction and assignment
        Plane (const Plane& other) {copy(other);}
        void operator= (const Plane& other) {copy(other);}
        // Access to the Plane's equation values
        size_t size () const {return 4;}
        float operator[] (const size_t) const;
        // The distance from an Atom to the plane, accounting for a movement
        float calculate_distance (const PROT::Atom *, const Spot *) const;

    // End the class definition
};

// Implement the methods of the Plane class
INTER::Plane::Plane () {
    m_mag = 0.0; for(size_t i=0; i<4; ++i) {m_values[i] = 0.0;}}

INTER::Plane::Plane (const float n1, const float n2, const float n3,
                     const float n4) {
    m_values[0] = n1; m_values[1] = n2; m_values[2] = n3; m_values[3] = n4;
    // Note that this should NOT include the fourth value
    m_mag = sqrt(pow(n1,2) + pow(n2, 2) + pow(n3, 2));
}

float INTER::Plane::operator[] (const size_t i) const {
    // Validate the index
    if (i >= 4) {
        stringstream c1; c1 << i;
        string error = "A Plane contains 4 floating point numbers. " 
                     + c1.str() + " is not a valid index to access one.\n";
        throw logic_error (error);}
    // Return the value
    return m_values[i];
}

// Calculate the distance from the Plane to an Atom after the Atom has been
// moved
float INTER::Plane::calculate_distance (const PROT::Atom * atom,
                    const Spot * spot) const {
    // Calculate the absolute value of applying the Plane's coefficients to
    // the coordinates
    float dis = m_values[3];
    for(size_t i=0; i<3; ++i) {
        dis += m_values[i] * (atom->operator[](i) - spot->m_points[i]);}
    // Divide by the magnitude of the plane's mutliplier coefficients
    dis = abs(dis) / m_mag;
    return dis;
}

// Define the Interaction class
class INTER::Interaction {

    // The information stored in the class is private
    private:
        // The interaction knows information about where it goes in the
        // paratope
        string m_loop;
        size_t m_struct;
        size_t m_pos;
        size_t m_rot;
        // The Interaction will contain one or more Spots
        vector<Spot> m_spots;
        // The Interaction may contain a Plane
        Plane m_plane;
        // The Interaction knows what type of Interaction it is
        size_t m_type;
        // The interaction may contain X, Y and Z coordinates. For
        // interactions 1-6 and 10, it is important that they not be placed
        // too near one another, which is what those coordinates are for
        float m_X, m_Y, m_Z;

    // Private functions involved in the class behaviour
    private:
        // Copy the information from another instance of the class
        void copy (const Interaction&);
        // Functions for evaluating each of the possible interaction types
        bool Spot_and_Plane (const vector<PROT::Atom *>&, const float) const;
        bool Double_Spot (const vector<PROT::Atom *>&, const float) const;
        bool Spot_and_Center (const vector<PROT::Atom *>&, const float) const;
        bool Triple_Spot (const vector<PROT::Atom *>&, const float) const;
        bool Center_and_Plane (const vector<PROT::Atom *>&, const float) const;

    // The public interface of the class
    public:
        // A default constructor for the class
        Interaction ();
        // The standard constructor specifies the type of interaction and
        // provides a string containing the information
        Interaction (const size_t, const string&);
        // Copy construction and assignment functions
        Interaction (const Interaction& other) {copy(other);}
        void operator= (const Interaction& other) {copy(other);}
        // Access to the class information
        string loop () const {return m_loop;}
        size_t structure () const {return m_struct;}
        size_t position () const {return m_pos;}
        size_t rotamer () const {return m_rot;}
        // The type of interaction
        size_t kind () const {return m_type;}
        // Get the movement of placing an Atom in the guaranteed spot of the
        // interaction
        Spot movement (const vector<PROT::Atom *>&) const;
        // Evaluate how well a set of Atoms matches the Interactions
        // information
        bool evaluate (const vector<PROT::Atom *>&, const float) const;
        // Return True if the Interaction is of type 1, 2, 5, 6, or 10. Those
        // interactions need some special consideration for certain
        // calculations
        bool is_special () const;
        float X () const {return m_X;}
        float Y () const {return m_Y;}
        float Z () const {return m_Z;}

    // End the class definition
};

// The copy function of the Interaction class
void INTER::Interaction::copy (const Interaction& other) {
    // Copy most values directly
    m_loop = other.m_loop;
    m_struct = other.m_struct;
    m_pos = other.m_pos;
    m_rot = other.m_rot;
    m_plane = other.m_plane;
    m_type = other.m_type;
    m_X = other.m_X;
    m_Y = other.m_Y;
    m_Z = other.m_Z;
    // Copy the spots
    m_spots.clear(); m_spots.reserve(other.m_spots.size());
    for(size_t i=0; i<other.m_spots.size(); ++i) {
        m_spots.push_back(other.m_spots[i]);}
}

// The default constructor of the Interaction class
INTER::Interaction::Interaction () {
    m_loop = "NONE";
    m_struct = 0;
    m_pos = 0;
    m_rot = 0;
    m_plane = Plane ();
    m_type = 0;
    m_X = 0.0;
    m_Y = 0.0;
    m_Z = 0.0;
}

// The standard constructor of the Interaction class
INTER::Interaction::Interaction (const size_t N, const string& line) {
    // Set the type
    m_type = N;
    // Throw an error if the type is not supported
    if (((N == 0) || (N > 23)) || ((N == 3) || (N == 4))) {
        stringstream c; c << N;
        string error = c.str() + " is not a supported Interaction type\n";
        throw logic_error (error);}
    // Set the X, Y and Z coordinates to 0, as they are not used for most
    // interaction types
    m_X = 0.0; m_Y = 0.0; m_Z = 0.0;
    // Put the line into a stringstream
    stringstream text; text << line;
    // Every case starts with the labelling information
    text >> m_loop; text >> m_struct; text >> m_pos; text >> m_rot;
    // There are three combinations of data: a spot and a plane, two spots, or
    // three spots. A spot and a plane are most common, so that will be the
    // "else" condition. Two spots are:
    if ((((N == 7) || (N == 8)) || (N == 12)) || 
        (((N == 16) || (N == 18)) || (N == 23))) {
        // 2 times read in 3 floating point numbers, then store them as spots
        for(size_t i=0; i<2; ++i) {
            float n1, n2, n3;
            text >> n1; text >> n2; text >> n3;
            m_spots.push_back(Spot(n1, n2, n3));}}
    // Three spots are
    else if ((N == 9) || (N == 11)){
        // 3 times read in 3 floating point numbers and make a spot
        for(size_t i=0; i<3; ++i) {
            float n1, n2, n3;
            text >> n1; text >> n2; text >> n3;
            m_spots.push_back(Spot(n1, n2, n3));}}
    // Interactions 1, 2, 5, 6, and 10 contain the X, Y and Z coordinates
    // (before the spot) so that it can be checked that they are not placed
    // too near each other
    else if (((N == 1) || (N == 2)) || (((N == 5) || (N == 6)) || (N == 10))) {
        // Store the X, Y and Z coordinates
        text >> m_X; text >> m_Y; text >> m_Z;
        // Store a Spot
        float n1, n2, n3, n4;
        text >> n1; text >> n2; text >> n3;
        m_spots.push_back(Spot(n1, n2, n3));
        // Store the plane
        text >> n1; text >> n2; text >> n3; text >> n4;
        m_plane = Plane (n1, n2, n3, n4);}
    // The other options all involve a spot and a plane
    else {
        // Store a spot
        float n1, n2, n3, n4;
        text >> n1; text >> n2; text >> n3;
        m_spots.push_back(Spot(n1, n2, n3));
        // Store the plane
        text >> n1; text >> n2; text >> n3; text >> n4;
        m_plane = Plane(n1, n2, n3, n4);}
    // End the function
}

// Get the Spot that governs how far an Atom would have to be moved to fit
// into the first Spot in the Interaction
INTER::Spot INTER::Interaction::movement (const vector<PROT::Atom *>& atoms) const {
    // For a certain set of interactions, calculate the center of mass of a
    // ring
    if ((((m_type == 13) || (m_type == 14)) || (m_type == 19)) ||
        (((m_type == 20) || (m_type == 21)) || (m_type == 22))) {
        // Create a vector to store the coordinates
        vector<float> n; n.resize(3);
        // Loop through the coordinates
        for(size_t i=0; i<3; ++i) {
            n[i] = 0.0;
            // Loop through the atoms
            for(size_t j=0; j<atoms.size(); ++j) {
                // Add the value to the one stored in the vector
                n[i] += atoms[j]->operator[](i);}
            // turn that into an average
            n[i] /= ((float) atoms.size());
            // Store that value minus the value in the first spot
            n[i] = n[i] - m_spots[0][i];}
        // Return a spot made from those three values
        return Spot (n[0], n[1], n[2]);}
    // Otherwise, just use the calculate movement function of the first spot
    // on the first atom
    return m_spots[0].calculate_movement(atoms[0]);
}

// The private functions that determine whether or not a set of Atoms are
// compatible with the interaction
bool INTER::Interaction::Spot_and_Plane (const vector<PROT::Atom *>& atoms,
                                         const float threshold) const {
    // Get the spot for moving the first atom
    Spot how = movement(atoms);
    // For every remaining Atom, calculate its distance from the plane
    for(size_t i=1; i<atoms.size(); ++i) {
        float distance = m_plane.calculate_distance(atoms[i], &how);
        // If the distance is greater than the threshold, the function fails
        if (distance > threshold) {return false;}}
    return true;
}

bool INTER::Interaction::Double_Spot (const vector<PROT::Atom *>& atoms,
                                      const float threshold) const {
    // Get the spot for moving the first atom
    Spot how = movement(atoms);
    // Calculate the distance with the second spot
    float distance = m_spots[1].calculate_distance(atoms[1], &how);
    return (distance <= threshold);
}

bool INTER::Interaction::Spot_and_Center (const vector<PROT::Atom *>& atoms,
                                          const float threshold) const {
    // Get the spot for moving the first atom
    Spot how = movement(atoms);
    // Calculate the center of mass of the ring
    vector<float> n; n.resize(3);
    for(size_t i=0; i<3; ++i) {
        n[i] = 0.0;
        for(size_t j=0; j<atoms.size(); ++j) {
            n[i] += atoms[j]->operator[](i);}
        n[i] /= ((float) atoms.size());}
    // Determine the distance between that vector and the second spot
    float distance = m_spots[1].calculate_distance (n, &how);
    return (distance <= threshold);
}

bool INTER::Interaction::Triple_Spot (const vector<PROT::Atom *>& atoms,
                                      const float threshold) const {
    Spot how = movement(atoms);
    for(size_t i=1; i<3; ++i) {
        float distance = m_spots[i].calculate_distance(atoms[i], &how);
        if (distance > threshold) {return false;}}
    return true;
}

bool INTER::Interaction::Center_and_Plane (const vector<PROT::Atom *>& atoms,
                                           const float threshold) const {
    // Get the movement spot
    Spot how = movement(atoms);
    // For every atom, determine whether or not it is in the plane
    for(size_t i=0; i<atoms.size(); ++i) {
        float distance = m_plane.calculate_distance(atoms[i], &how);
        if (distance > threshold) {return false;}}
    return true;
}

// Evaluate whether or not a set of Atoms is compatible with the Interaction.
bool INTER::Interaction::evaluate (const vector<PROT::Atom *>& atoms,
                                   const float threshold) const {
    // Return the value from the appropriate method. 
    if ((m_type == 9) || (m_type == 11)) {return Triple_Spot(atoms, threshold);}
    else if ((m_type == 18) || (m_type == 23)) {
        return Spot_and_Center (atoms, threshold);}
    else if (((m_type == 7) || (m_type == 8)) || 
             ((m_type == 12) || (m_type == 16))) {
        return Double_Spot(atoms, threshold);}
    else if ((((m_type == 13) || (m_type == 14)) || (m_type == 19)) ||
             (((m_type == 20) || (m_type == 21)) || (m_type == 22))) {
        return Center_and_Plane (atoms, threshold);}
    else if ((((m_type == 1) || (m_type == 2)) || 
              ((m_type == 5) || (m_type == 6))) ||
             (((m_type == 10) || (m_type == 15)) || (m_type == 17))) {
        return Spot_and_Plane (atoms, threshold);}
    // Return false for any other interaction type
    return false;
}

// Determine whether or not the Interaction needs special considerations for
// Position compatibility calculations
bool INTER::Interaction::is_special () const {
    if (((m_type == 1) || (m_type == 2)) ||
       (((m_type == 5) || (m_type == 6)) || (m_type == 10))) {return true;}
    return false;
}

// Define the Position class
class INTER::Position {

    // The information stored in the class is private
    private:
        // A pointer to the Interaction involved
        Interaction * m_inter;
        // A spot associated with the movement of the interaction
        Spot m_spot;
        // A "score" for that spot, which is the sum of the X, Y and Z motions
        float m_score;
        // Information about the antigen residue
        string m_res;
        char m_prot;

    // Private functions that control class behaviour
    private:
        void copy (const Position& other);

    // The public interface of the class
    public:
        // A default constructor
        Position ();
        // The standard constructor
        Position (Interaction *, const vector<PROT::Atom *>&, const string&, 
                  const char);
        // Copy construction and assignment
        Position (const Position& other) {copy(other);}
        void operator= (const Position& other) {copy(other);}
        // Access to the class information
        float score () const {return m_score;}
        string residue () const {return m_res;}
        char protein () const {return m_prot;}
        string loop () const {return m_inter->loop();}
        size_t structure () const {return m_inter->structure();}
        size_t position () const {return m_inter->position();}
        size_t rotamer () const {return m_inter->rotamer();}
        size_t interaction () const {return m_inter->kind();}
        Spot * movement () {return &m_spot;}
        float operator[] (const size_t i) const {return m_spot[i];}
        // Determine whether or not two Positions involve similar enough
        // movements to be mutually compatible
        bool distance_compatible (const Position *) const;
        bool compatible (const Position *, const float) const;
        // A string summarizing the Position's information
        string str () const;

    // End the class definition
};

// Implement the methods of the Position class, starting with the copy
// function
void INTER::Position::copy (const Position& other) {
    m_inter = other.m_inter;
    m_spot = other.m_spot;
    m_score = other.m_score;
    m_res = other.m_res;
    m_prot = other.m_prot;
}

// The default constructor of the Position class
INTER::Position::Position () {
    m_inter = 0;
    m_spot = Spot ();
    m_score = 1.0e15;
    m_res = "NONE";
    m_prot = ' ';
}

// The standard constructor of the Protein class
INTER::Position::Position (Interaction * inter, 
                 const vector<PROT::Atom *>& atoms, const string& res,
                 const char L) {
    // Store the information
    m_inter = inter;
    m_res = res;
    m_prot = L;
    m_spot = m_inter->movement(atoms);
    // Add up the X, Y, and Z coordinates of the spot
    m_score = 0.0;
    for(size_t i=0; i<3; ++i) {m_score += m_spot[i];}
}

// Do a check to see whether or not two Positions are incompatible because the
// starting points of their interactions are too close to each other in the
// binding protein
bool INTER::Position::distance_compatible (const Position * other) const {
    // This only matters if both interactions fall into the "special"
    // consideration category
    if ((m_inter->is_special()) && (other->m_inter->is_special())) {
        // Calculate the distance between the trigger points of the
        // interactions
        float dis = sqrt (pow(m_inter->X() - other->m_inter->X(), 2) +
                          pow(m_inter->Y() - other->m_inter->Y(), 2) +
                          pow(m_inter->Z() - other->m_inter->Z(), 2));
        // This distance threshold is arbitrarily selected and could be
        // refined in the future
        if (dis < 2.0) {return false;}}
    return true;
}

// Determine whether or not two positions are compatible
bool INTER::Position::compatible (const Position * other, 
                                  const float threshold) const {
    // Loop through the 3 coordinates
    for(size_t i=0; i<3; ++i) {
        // Get the absolute value of the difference in the movement value
        float diff = abs(m_spot[i] - other->m_spot[i]);
        // If that differnce is greater than the threshold, they are not
        // compatible
        if (diff > threshold) {return false;}}
    // If both Positions' Interactions are of types 1, 2, 5, 6 or 10, make 
    // sure that their X, Y and Z coordinates are not too close together
    if (!distance_compatible(other)) {return false;}
    // If the function reached this point, there is no reason that the
    // Positions are incompatible with one another
    return true;
}

// A string summarizing the interaction listed in the Position
string INTER::Position::str () const {
    // Store the output here
    string output = "";
    size_t m_type = interaction();
    // Based on the type of interaction, provide a description of what the
    // interaction is
    if (m_type <= 8) {output += "A hydrogen bond";}
    else if (m_type <= 12) {output += "A strong hydrogen bond complex";}
    else if ((m_type <= 16) || (m_type >= 20)) {
        output += "A pi-cation interaction";}
    else {output += "A pi-pi interaction";}
    // List the residue in the paratope
    stringstream c1; c1 << (m_inter->position() + 1);
    stringstream c2; c2 << m_inter->structure();
    output += " between the " + c1.str();
    if (m_inter->position() == 0) {output += "st";}
    else if (m_inter->position() == 1) {output += "nd";}
    else if (m_inter->position() == 2) {output += "rd";}
    else {output += "th";}
    output += " residue in " + m_inter->loop() + " structure " + c2.str()
            + " and antigen residue " + m_res + " in chain ";
    output += m_prot;
    output += "\n";
    return output;
}

// Implement the Solution class
class INTER::Solution {

    // The information stored in the class is private
    private:
        // The Positions that are part of this Solution
        vector<Position> m_positions;
        // The X, Y and Z movements
        float m_X, m_Y, m_Z;
        // The X, Y and Z rotations
        float m_xa, m_ya, m_za;

    // Private functions that help control class behaviour
    private:
        // Copy information from one instance of the class to another
        void copy (const Solution&);

    // The public interface of the class
    public:
        // A default constructor
        Solution ();
        // The standard constructor
        Solution (vector<Position>&, const vector<size_t>&, 
                  const float, const float, const float);
        // Copy construction and assignment
        Solution (const Solution& other) {copy(other);}
        void operator= (const Solution& other) {copy(other);}
        // Access to the class information
        size_t size () const {return m_positions.size();}
        float X () const {return m_X;}
        float Y () const {return m_Y;}
        float Z () const {return m_Z;}
        float XA () const {return m_xa;}
        float YA () const {return m_ya;}
        float ZA () const {return m_za;}
        Position * operator[] (const size_t);

    // End the class definition
};

// Implement the methods of the solution class, starting with the copy method
void INTER::Solution::copy (const Solution& other) {
    m_X = other.m_X; m_Y = other.m_Y; m_Z = other.m_Z;
    m_xa = other.m_xa; m_ya = other.m_ya; m_za = other.m_za;
    m_positions.clear();
    if (other.m_positions.size() > 0) {
        m_positions.reserve(other.m_positions.size());
        for(size_t i=0; i<other.m_positions.size(); ++i) {
            m_positions.push_back(other.m_positions[i]);}}
}

// The default constructor
INTER::Solution::Solution () {
    m_X = 0.0; m_Y = 0.0; m_Z = 0.0;
    m_xa = 0.0; m_ya = 0.0; m_za = 0.0;
}

// The standard constructor
INTER::Solution::Solution (vector<Position>& spots, const vector<size_t>& which,
                           const float n1, const float n2, const float n3) {
    // Store the angle information
    m_xa = n1; m_ya = n2; m_za = n3;
    // Store the proper set of Positions
    m_positions.reserve(which.size());
    for(size_t i=0; i<which.size(); ++i) {
        m_positions.push_back(spots[which[i]]);}
    // Calculate the average X, Y and Z movement information
    m_X = 0.0; m_Y = 0.0; m_Z = 0.0;
    for(size_t i=0; i<m_positions.size(); ++i) {
        m_X += m_positions[i][0];
        m_Y += m_positions[i][1];
        m_Z += m_positions[i][2];}
    float N = (float) m_positions.size();
    if (N > 0) {m_X /= N; m_Y /= N; m_Z /= N;}
}

// Access to the Positions in the solution
INTER::Position * INTER::Solution::operator[] (const size_t i) {
    // Validate the index
    if (i >= m_positions.size()) {
        stringstream c1; c1 << i;
        stringstream c2; c2 << m_positions.size();
        string error = "An index of " + c1.str() + " is not valid to access "
                       "a Position in a Solution with " + c2.str()
                     + " Positions\n";
        throw logic_error (error);}
    // Return the address of hte Position
    return &(m_positions[i]);
}

// End the header guard from the start of the file
#endif
