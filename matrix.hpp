#include <cassert>
#include <cstdio>
#include <utility>
#include <iostream>
#include <fstream>
#include <cstring>
#include <string>
#include <cstdbool>
#include <cmath>

using namespace std;

class Matrix final {
    double** matrix;
    unsigned row, col;
    
    public:
    Matrix();
    Matrix(unsigned, unsigned, double*);
	Matrix(unsigned, unsigned);
	Matrix(const Matrix&);
	Matrix(Matrix&&);
	~Matrix();
	void cleanup();

	unsigned get_row() const;
	unsigned get_col() const;
	double** get_matrix() const;
	
	const double& at(unsigned, unsigned) const;
	double& at(unsigned, unsigned);
	const double* operator[](unsigned) const;
	double* operator[](unsigned);

	Matrix& operator=(const Matrix&);
	Matrix& operator=(Matrix&&);

	Matrix& operator+=(const Matrix&);
	Matrix& operator-=(const Matrix&);
	Matrix& operator*=(double);
	Matrix& transpose();
	Matrix& inverse();

	Matrix& add_line(); 
	void correct();
	void row_decrement();

	void set_from_mult(const Matrix&, const Matrix&);
	Matrix& multEq(const Matrix&);

	void print() const;

	Matrix& kroneker(Matrix&);

	double determinant() const;

	Matrix& fill_I(const Matrix&, const Matrix&, const Matrix&, const Matrix&, const Matrix&);

	void fix_A(int& index);
	void fix_Y(int& index);
	void fix_JEI(int& index);
};
	
bool operator==(const Matrix&, const Matrix&);
bool operator!=(const Matrix&, const Matrix&);

Matrix operator+(const Matrix&, const Matrix&);
Matrix operator-(const Matrix&, const Matrix&);

Matrix operator*(const Matrix, double);
Matrix operator*(double, const Matrix);

Matrix mult(const Matrix&, const Matrix&);


class Circuit final {
	Matrix A; // matrix of connections
	Matrix Y; // matrix of conductivities
	Matrix U0; // matrix of node potentials
	Matrix J; // matrix of current sources
	Matrix E; // matrix of voltage sources
	Matrix I; // matrix of currents
	Matrix A_new; // extended matrix of connections

	public:
	Circuit();
	Circuit(unsigned, unsigned, double*, double*, double*, double*, double*);

	void print();
	void solve();
	void delete0();
};

namespace Parser {
	Circuit parse(const char*);

	void skip_spaces(string&, int&);
	bool parce_digit(string&, int&, int&);
	bool parse_number(string&, int&, int&);
	bool parse_dashes(string&, int&);
	bool parse_comma(string&, int&);
	bool parse_real(string&, double&, int&);
	bool parse_scolon(string&, int&);
	bool parse_A(string&, int&);
	bool parse_V(string&, int&);
	bool parse_str(string&, int&, int&, double&, double&, double&, int&);
}