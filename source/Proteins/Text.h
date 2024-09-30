/* Written in 2018 by Dr. Robert Pantazes. Copyright owned by Auburn
 * University.
 * This file includes functions that provide Python-like string editing
 * abilities. */

// Use a header guard to prevent this file from being included in a compiled
// program multiple times
#ifndef Text_Check
#define Text_Check 1

// Include the necessary C++ files
#include <string>
#include <vector>
#include <cstdlib>

// Use the standard namespace
using namespace std;

// Declare a namespace for the Text manipulation functions
namespace Text {

    // Functions that check properties of characters
    bool is_digit (const char);
    bool is_whitespace (const char);
    bool is_upper (const char);
    bool is_lower (const char);

    // Functions that check properties of a string
    bool startswith (const string&, const string&);
    bool startswith (const string&, const char *);
    bool startswith (const string&, const char);
    bool endswith (const string&, const string&);
    bool endswith (const string&, const char *);
    bool endswith (const string&, const char);
    bool contains (const string&, const string&);
    bool contains (const string&, const char *);
    bool contains (const string&, const char);
    bool is_integer (const string&);
    bool is_number (const string&);

    // Functions for removing whitespace from a string
    void rstrip (string&);
    void lstrip (string&);
    void strip (string&);

    // Functions that modify capitalization of strings
    void upper (string&);
    void lower (string&);
    void capitalize (string&);

    // Functions that modify the packing of a string
    void rjust (string&, const size_t, const char);
    void ljust (string&, const size_t, const char);

    // Functions that insert a string into another string with specified
    // packing
    void rjust_insert (string&, const string&, const size_t, const char);
    void ljust_insert (string&, const string&, const size_t, const char);

    // Functions that split a string on whitespace
    void split (vector<string>&, const string&);
    vector<string> split (const string&);
    // Functions that do the same thing, but using a specified character
    void split (vector<string>&, const string&, const char);
    vector<string> split (const string&, const char);

    // The end of the namespace
}

// Implement the methods of the Text namespace

// The is_digit function checks to see if a character is a digit
bool Text::is_digit (const char L) {return ((L >= 48) && (L <= 57));}

// This function checks to see if a character is whitespace
bool Text::is_whitespace (const char L) {return ((L <= 32) || (L == 127));}

// This function checks to see if a character is an upper case letter
bool Text::is_upper (const char L) {return ((L >= 65) && (L <= 90));}

// This function checks to see if a character is a lower case letter
bool Text::is_lower (const char L) {return ((L >= 97) && (L <= 122));}

// This function checks to see if the start of one string matches the contents
// of a second string
bool Text::startswith (const string& input, const string& has) {
    // If the input is empty, it does not start with the 'has' string
    if (input.size() == 0) {return false;}
    // If the input is shorter than the 'has' string, it cannot start with it
    else if (input.size() < has.size()) {return false;}
    // If the 'has' string is empty, the input starts with it
    else if (has.size() == 0) {return true;}
    // We've established that both strings contain content and the input is at
    // least as long as the has string. Check the characters one at a time
    for (size_t i=0; i<has.size(); ++i) {
        if (input[i] != has[i]) {return false;}}
    return true;
}

// Do the same thing using a character array instead of a string for 'has'
bool Text::startswith (const string& input, const char * has) {
    // If the input is empty, it does not start with has
    if (input.size() == 0) {return false;}
    // Define the null character
    char n = (char) 0;
    // Loop through the input string's characters
    for (size_t i=0; i<input.size(); ++i) {
        // If this is the null character in the 'has' array, the input matched
        // it
        if (has[i] == n) {return true;}
        // If the character's don't match
        else if (input[i] != has[i]) {return false;}}
    // If every character in the input string was checked, only return true if
    // the next has character is the null character
    return (has[input.size()] == n);
}

// Do the same thing for if a string startswith a specific character
bool Text::startswith(const string& input, const char has) {
    return ((input.size() > 0) && (input[0] == has));
}

// Write equivalent functions, checking the end of the input string instead of
// its start
bool Text::endswith (const string& input, const string& has) {
    if (input.size() == 0) {return false;}
    else if (input.size() < has.size()) {return false;}
    else if (has.size() == 0) {return true;}
    size_t I = input.size() - 1; size_t H = has.size() - 1;
    for(size_t i=0; i<has.size(); ++i) {
        if (input[I-i] != has[H-i]) {return false;}}
    return true;
}

bool Text::endswith (const string& input, const char * has) {
    if (input.size() == 0) {return false;}
    char n = (char) 0;
    // Determine the number of characters in the c-string
    size_t H = 0;
    while (has[H] != n) {
        ++H;
        // If there are more than 10000 characters in the string,
        // automatically return false
        if (H >= 10000) {return false;}}
    if (H > input.size()) {return false;}
    else if (H == 0) {return true;}
    size_t I = input.size() - 1; --H;
    for(size_t i=0; i<=H; ++i) {
        if (input[I-i] != has[H-i]) {return false;}}
    return true;
}

bool Text::endswith (const string& input, const char has) {
    return ((input.size() > 0) && (input[input.size()-1] == has));
}

// Check to see whether or not a substring is in the string
bool Text::contains (const string& input, const string& has) {
    size_t n = input.find(has);
    if (n >= input.size()) {return false;}
    return true;
}

bool Text::contains (const string& input, const char * has) {
    size_t n = input.find(has);
    if (n >= input.size()) {return false;}
    return true;
}

bool Text::contains (const string& input, const char has) {
    string check = "";
    check += has;
    size_t n = input.find(check);
    if (n >= input.size()) {return false;}
    return true;
}

// Check to see whether or not a string contains an integer
bool Text::is_integer (const string& input) {
    // If the string is empty, it does not contain an integer
    if (input.size() == 0) {return false;}
    // Loop through the characters of the string
    for (size_t i=0; i<input.size(); ++i) {
        // Get the character
        char c = input[i];
        // Check to see if it is a digit
        if (is_digit(c)) {continue;}
        // If it is the first digit, it may also be '-' or '+'
        else if (((i == 0) && (input.size() > 1)) && ((c == 43) || (c == 45))){
            continue;}
        // Otherwise, this is not an integer
        return false;}
    // If every character was acceptable, return true
    return true;
}

// Determine whether or not a string contains a number
bool Text::is_number (const string& input) {
    // If the string is empty, it does not contain a number
    if (input.size() == 0) {return false;}
    // Scientific notation permits 'e' to be used in a number. This variable
    // tracks when one occurs
    int e = 0;
    // This boolean value tells whether or not a decimal has been seen
    bool decimal = false;
    // Loop through the characters in the string
    for (size_t i=0; i<input.size(); ++i) {
        // If an e has been seen, increment the counter. This matters for
        // tracking whether or not the next character is '+' or '-' near the
        // end of this for loop
        if (e > 0) {++e;}
        // Extract the character
        char c = input[i];
        // If it is a digit, it is fine
        if (is_digit(c)) {continue;}
        // If it is the first character of at least 2, it may be '+' or '-'
        else if (((i == 0) && (input.size() > 1)) && ((c == 43) || (c == 45))){
            continue;}
        // If no 'e' value has been seen yet and there is at least one more
        // subsequent character, a '.' is permitted
        else if (((e == 0) && (c == 46)) && 
                 ((i < input.size()-1) && (decimal == false))){
            decimal = true; continue;}
        // If this is an 'e' or 'E', is not the first character, and there are
        // more characters after it
        else if ((((c == 69) || (c == 101)) &&
                  ((i > 0) && (i < input.size() - 1))) && (e == 0)) {
            e = 1; continue;}
        // If this is the first character after an 'e' (i.e. the e variable ==
        // 2), it may be + or -
        else if (((e == 2) && (i < input.size() - 1)) && ((c==43)||(c==45))){
            continue;}
        // Any other value indicates that this is not a number
        else {return false;}}
    return true;
}

// Remove trailing whitespace from a string
void Text::rstrip (string& input) {
    // If the input is empty, don't do anything
    if (input.size() == 0) {return;}
    // Identify the last non-whitespace character in the string
    size_t last = input.size();
    for(size_t i=last-1; i>=0; --i) {
        if(!is_whitespace(input[i])) {last = i; break;}}
    // If the value was never set, the entire string is whitespace and should
    // be cleared
    if (last == input.size()) {input.clear();}
    // If the last character is not whitespace, don't do anything
    else if (last == input.size() - 1) {return;}
    // Delete the appropriate characters
    input.erase(last+1, input.size() - (last+1));
    return;
}

// Remove leading whitespace from a string
void Text::lstrip (string& input) {
    if (input.size() == 0) {return;}
    size_t first = input.size();
    for(size_t i=0; i<input.size(); ++i) {
        if (!is_whitespace(input[i])) {first = i; break;}}
    if (first == 0) {return;}
    else if (first == input.size()) {input.clear();}
    input.erase(0, first);
    return;
}

// Remove leading and trailing whitespace from a string
void Text::strip (string& input) {lstrip(input); rstrip(input);}

// Convert a string to all upper case letters
void Text::upper (string& input) {
    // if the string is empty, be done
    if (input.size() == 0) {return;}
    // Loop through its characters
    for(size_t i=0; i<input.size(); ++i) {
        // If the character is a lower case letter
        if(is_lower(input[i])) {
            // Cast the character to a short and subtract 32
            short n = ((short) input[i]) - 32;
            // Replace the character in the string
            input.replace(i, 1, 1, (char) n);}}
    return;
}

// Convert a string to all lower case letters
void Text::lower (string& input) {
    if (input.size() == 0) {return;}
    for (size_t i=0; i<input.size(); ++i) {
        if (is_upper(input[i])) {
            short n = ((short) input[i]) + 32;
            input.replace(i, 1, 1, (char) n);}}
    return;
}

// Change the first character in a string to a capital letter
void Text::capitalize(string & input) {
    if ((input.size() > 0) && (is_lower(input[0]))) {
        short n = ((short) input[0]) - 32;
        input.replace(0, 1, 1, (char) n);}
    return;
}

// Modify a string to be at least a certain length by inserting leading
// characters
void Text::rjust (string& input, const size_t n, const char L) {
    if (n > input.size()) {input.insert(0, n-input.size(), L);}
    return;
}

// Do the same by adding the characters to the end of the string
void Text::ljust (string& input, const size_t n, const char L) {
    if (n > input.size()) {input.append(n-input.size(), L);}
    return;
}

// Do similar things, but add the padded strings to the end of another string
// as part of the process
void Text::rjust_insert(string& text, const string& use, const size_t n,
                        const char L) {
    if (n > use.size()) {text.append(n-use.size(), L);}
    text.append(use);
    return;
}

void Text::ljust_insert (string& text, const string& use, const size_t n,
                         const char L) {
    text.append(use);
    if (n > use.size()) {text.append(n-use.size(), L);}
    return;
}

// Split a string into its whitespace-separated substrings
void Text::split (vector<string>& output, const string& input) {
    // Make sure the vector is empty
    output.clear();
    // If the input is empty, be done
    if (input.size() == 0) {return;}
    // The position of the first character in a substring
    size_t first = 0;
    // The number of characters in the substring
    size_t count = 0;
    // Whether or not there a substring has been identified
    bool active = false;
    // Loop through the characters
    for (size_t i=0; i<input.size(); ++i) {
        // If this character is whitespace
        if (is_whitespace(input[i])) {
            // If there is an active substring to store
            if (active) {
                // Store the substring
                output.push_back(input.substr(first, count));
                // Reset the tracking information
                active = false; count = 0;}}
        // If the character is not whitesapce
        else {
            // If this is the first character in a new substring
            if (!active) {active = true; first = i;}
            // Increment the counter
            ++count;}}
    // If there is an unstored active substring going at the end of the string
    if (active) {output.push_back(input.substr(first, count));}
    return;
}

// Do the same, but return the vector as part of the function
vector<string> Text::split (const string& input) {
    vector<string> output;
    split(output, input);
    return output;
}

// Split a string on a specified character instead of on whitespace
void Text::split (vector<string>& output, const string& input, const char L) {
    output.clear();
    if (input.size() == 0) {return;}
    size_t first = 0;
    size_t count = 0;
    bool active = false;
    for (size_t i=0; i<input.size(); ++i) {
        if (input[i] == L) {
            if (active) {
                output.push_back(input.substr(first, count));
                active = false; count = 0;}}
        else {
            if (!active) {active = true; first = i;}
            ++count;}}
    if (active) {output.push_back(input.substr(first, count));}
    return;
}

vector<string> Text::split (const string& input, const char L) {
    vector<string> output;
    split(output, input, L);
    return output;
}

// End the header guard from the start of the file
#endif
