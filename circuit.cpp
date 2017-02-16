#include "matrix.hpp"

Circuit::Circuit() {}

Circuit::Circuit(unsigned p, unsigned q, double* a, double* y, double* j, double* e, double* a_new) { // allocate memory for data
	A = {q, p, a}; // matrix of connections
	Y = {p, p, y}; // matrix of conductivities
	U0 = {q, 1}; // matrix of node potentials
	J = {p, 1, j}; // matrix of current sources
	E = {p, 1, e}; // matrix of voltage sources
	I = {p, 1}; // matrix of currents
	A_new = {q + 1, p, a_new}; // extended matrix of connections
}

void Circuit::print() {
	assert(A.get_matrix() && Y.get_matrix() && U0.get_matrix() && J.get_matrix() && E.get_matrix() && I.get_matrix() && A_new.get_matrix());
	printf("A:\n");
	A.print();
	printf("Y:\n");
	Y.print();
	printf("J:\n");
	J.print();
	printf("E:\n");
	E.print();
	printf("A_new:\n");
	A_new.print();
	printf("U0:\n");
	U0.print();
	printf("I:\n");
	I.print();
}

void Circuit::solve() {
	assert(A.get_matrix() && Y.get_matrix() && U0.get_matrix() && J.get_matrix() && E.get_matrix() && I.get_matrix() && A_new.get_matrix());

	Matrix A1 = A;

	Matrix First = mult((A * (-1)), mult(Y, E) + J);
    Matrix Second = mult(mult(A, Y), A1.transpose()).inverse();
    U0 = mult(Second, First);
   	U0.correct();
    
	// add a line to matrix of node potentials and to matrix of connections,
	// considering that the full number of nodes is (q + 1)
	U0.add_line();
	printf("Node potentials:\n");
	U0.print();

  	I.fill_I(A_new, Y, U0, J, E);
   	printf("Currents (in the same order as the order of connections in input file) :\n");
  	I.correct();
  	I.print();
  	printf("------------------------------------------------------\n");
  	printf("------------------------------------------------------\n");
}

// delete zero connections
void Circuit::delete0() {
	for (int i = 0; i < Y.get_row(); i++)
		if (Y.get_matrix()[i][i] == 1 / 0.00000000001)
			if ((J.get_matrix()[i][0] == 0) && (E.get_matrix()[i][0] == 0)) {
			  	Y.fix_Y(i);
				J.fix_JEI(i);
				E.fix_JEI(i);
				I.fix_JEI(i);
				U0.fix_JEI(i);
				A_new.fix_A(i);
			}
	A = A_new;
	A.row_decrement();
}

Matrix& Matrix::fill_I(const Matrix& A, const Matrix& Y, const Matrix& U0, const Matrix& J, const Matrix& E) {
	// split a matrix A on the columns
	Matrix** Ai = new Matrix* [A.get_col()]; 
	for (unsigned i = 0; i < A.get_col(); i++)
    	Ai[i] = new Matrix (A.get_row(), 1);

	for (unsigned i = 0; i < row; i++)
    	for (unsigned j = 0; j < Ai[0]->get_row(); j++)
    		Ai[i]->matrix[j][0] = A.matrix[j][i];

	int err = -1;
	// get an answer
	for (unsigned i = 0; i < row; i++) {
    	Matrix One = mult(Ai[i]->transpose(), U0);
    	// in a case when potentials are not equal and resistance is 0
   		// we cannot find a current by formula, we will find them below 
    	// as a sum of the currents in one of the connection's node's -->
    	if (abs(Y.matrix[i][i]) < 1000000)
    		matrix[i][0] = (One.matrix[0][0] + E.matrix[i][0]) * Y.matrix[i][i] + J.matrix[i][0];
    	else {
      		err = i;
      		matrix[i][0] = 0;
    	}      
  	}

	// --> here:
  	if (err > 0) {
   		unsigned j = 0;
    	while (j < A.get_col()) {
      		if (A.matrix[j][err] != 0) {
	        	Matrix Two(1, A.get_col(), A.matrix[j]);
	        	Matrix tmp = *this;
	        	Matrix tmp2 = mult(Two, tmp);
	        	matrix[err][0] = -tmp2[0][0];
	        	break;
      		}
      		j++;
    	}
  	}

	// in a case when potentials are equal the current does not flow;
  	// we take into account the inaccuracy of the multiplication
  	correct();

 	for (unsigned i = 0; i < A.get_col(); i++)
  		Ai[i]->cleanup();

  	delete[] Ai;

  	return *this; 
}