#include "matrix.hpp"

// PARSER of a file:
Circuit Parser::parse(const char* filename) {
    // count the number of lines in a file = p = the number of edges
    ifstream fin(filename);
    assert(fin.is_open());

    string str;
    unsigned p = 0;
    while (!fin.eof())
        while (getline(fin, str))
            p++;

    str.erase();
    fin.close();
    //-----------------------------------------------------

    string* s = new string[p]; // array of strings
    int* q1 = new int[p]; // 1st colomn of nodes
    int* q2 = new int[p]; // 2nd colomn of nodes
    double* r = new double[p]; // colomn of resistances
    double* y = new double[p * p]; // array of conductivities
    double* j = new double[p]; // colomn of current sources
    double* e = new double[p]; // colomn of voltage sources

    // read a file and fill arrays
    ifstream fin1(filename);
    assert(fin1.is_open());
        
    int i = 0;
    while (!fin1.eof()) {
        getline(fin1, s[i]);
        assert(parse_str(s[i], q1[i], q2[i], r[i], j[i], e[i], i));
        i++;
    }

    fin1.close();
    //-------------------------------------------------------------------

    //count a number of nodes q as a max node number 
    unsigned q = q1[0];
    for (int i = 0; i < p; i++) {
        if (q1[i] > q)
            q = q1[i];
        if (q2[i] > q)
            q = q2[i];
    }

    // fill an array of conductivities
    for (int i = 0; i < p; i++)
        for (int j = 0; j < p; j++)
            if (i == j)
                y[j + i * p] = 1 / r[i];
            else
                y[j + i * p] = 0;

    double* a = new double[q * p]; // array of connections
        
    // fill an array of connections
    for (int i = 0; i < q; i++)
        for (int j = 0; j < p; j++)
            a[j + i * p] = 0;
    // accepting the node with the highest number for a zero potential
    for (int i = 0; i < p; i++) {
        if (q1[i] != q)
            a[i + q1[i] * p] = 1;
        if (q2[i] != q)
            a[i + q2[i] * p] = -1;
    }
    //-------------------------------------------------------------------

    // add a line to matrix of connections
    double* a_new = new double[(q + 1) * p];
    for (int j = 0; j < p; j++) {
        double sum = 0;
        for (int i = 0; i < q; i++) {
                sum += a[j + i * p];
                a_new[j + i * p] = a[j + i * p];
        }
        a_new[j + q * p] = -sum;
    }
    //-------------------------------------------------------------------

    Circuit MyCircuit(p, q, a, y, j, e, a_new);

    delete[] s;
    delete[] q1;
    delete[] q2;
    delete[] a;
    delete[] r;
    delete[] y;
    delete[] j;
    delete[] e;
    delete[] a_new;

    return MyCircuit;
}
//-----------------------------------------------------------------------------------------

// PARCER of 1 string:

void Parser::skip_spaces(string& s, int& index) {
    while (s[index] == ' ')
        index++;
}

bool Parser::parce_digit(string& s, int& digit, int& index) {
    const char* DIGIT[] = {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9"};
    for (int i = 0; i < 10; i++)
        if (s[index] == *(DIGIT[i])) {
            digit = atoi(DIGIT[i]);
            index++;
            return true;
        }
    return false;
}

bool Parser::parse_number(string& s, int& num, int& index) {
    int digit;
    skip_spaces(s, index);
    if (!parce_digit(s, digit, index))
        return false;

    num = digit;
    while (parce_digit(s, digit, index)) { // make up the number of digits
        num = num * 10 + digit;
        index++;
    }
    return true;
}

bool Parser::parse_dashes(string& s, int& index) {
    skip_spaces(s, index);
    char DASH = '-';
    if ((s[index] == '-') && (s[index + 1] == '-')) {
        index++;
        index++;
        return true;
    }
    return false; 
}

bool Parser::parse_comma(string& s, int& index) {
    skip_spaces(s, index);
    char COMMA = ',';
    if (s[index] == COMMA) {
        index++;
        return true;
    }
    return false; 
}

bool Parser::parse_real(string& s, double& r, int& index) {
    skip_spaces(s, index);
    const char DOT = '.';
    int digit = 0;

    int signum = 1;
    if (s[index] == '-') {
        signum = -1;
        index++;
    }

    if (!parce_digit(s, digit, index))
        return false;

    r = digit;
    while (parce_digit(s, digit, index)) // make an integer part
        r = r * 10 + digit;

    int i = -1; // if real < 0
    if (s[index] == DOT) {
        index++;
        while (parce_digit(s, digit, index)) { // make a fractional part
            r = r + digit * pow(10, i);
            i--;
        }
    }

    r *= signum;
    return true;
}

bool Parser::parse_scolon(string& s, int& index) {
    skip_spaces(s, index);
    const char SCOLON = ';';
        if (s[index] == SCOLON) {
        index++;
        return true;
    }
    return false;
}

bool Parser::parse_A(string& s, int& index) {
    skip_spaces(s, index);
    if (s[index] == 'A') {
        index++;
        return true;
    }

    return false; 
}

bool Parser::parse_V(string& s, int& index) {
    skip_spaces(s, index);
    if (s[index] == 'V') {
        index++;
        return true;
    }

    return false; 
}

// all these functions are bollean and fill the variables implicitly
bool Parser::parse_str(string& s, int& q1, int& q2, double& r, double& j, double& e, int& line) {

    int index = 0;

    if (!parse_number(s, q1, index)) {
        printf("number1 at index %d in line %d\n\n", index, line);
        return false; // string does not start with a digit
    }
    q1--;
    
    if (!parse_dashes(s, index)) {
        printf("dashes at index %d in line %d\n\n", index, line);
        return false; // "--" missed
    }
    
    if (!parse_number(s, q2, index)) {
        printf("number2 at index %d in line %d\n\n", index, line);
        return false;
    }
    q2--;

    if (!parse_comma(s, index)) {
        printf("comma at index %d in line %d\n\n", index, line);
        return false; // "," missed
    }

    if ((!parse_real(s, r, index)) || (r < 0)) {
        printf("real at index %d in line %d\n\n", index, line);
        return false;
    }

    // in a case when a connection does not have a resistance (its R = 0);
    // we should get rid of infinity, which will multiply incorrectly
    if (r == 0) 
      r = 0.00000000001;

    if (!parse_scolon(s, index)) {
        printf("scolon at index %d in line %d\n\n", index, line);
        return false; // ";" missed
    }

    if (parse_real(s, j, index)) {

        if (parse_A(s, index)) 
            e = 0;
        else if (parse_V(s, index)) {
            e = j;
            j = 0;
        }
        else {
            printf("not A, not V at index %d in line %d\n\n", index, line);
            return false;
        }

    } else {
        e = 0;
        j = 0;
    }

    return true;
}