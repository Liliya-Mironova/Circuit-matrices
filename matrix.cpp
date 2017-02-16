#include "matrix.hpp"

// default constructor
Matrix::Matrix() {
  matrix = nullptr;
  row = 0;
  col = 0;
}
  
Matrix::Matrix(unsigned rows, unsigned cols, double* arr) : row(rows), col(cols) {
  matrix = new double* [row];
  matrix[0] = new double [row * col];
  for (unsigned i = 1; i < row; i++)
    matrix[i] = matrix[i-1] + col;
      
  unsigned k = 0;
  for (unsigned i = 0; i < row; i++)
    for (unsigned j = 0; j < col; j++) {
      matrix[i][j] = arr[k];
      k++;
    }
}

Matrix::Matrix(unsigned rows, unsigned cols) : row(rows), col(cols) {
  matrix = new double* [row];
  matrix[0] = new double [row * col];
  for (unsigned i = 1; i < row; i++)
    matrix[i] = matrix[i-1] + col;
      
  for (unsigned i = 0; i < row; i++)
    for (unsigned j = 0; j < col; j++)
      matrix[i][j] = 0;
}

// copy constructor
Matrix::Matrix(const Matrix& other) : row(other.row), col(other.col) {
  matrix = new double* [row];
  matrix[0] = new double [row * col];
  for (unsigned i = 1; i < row; i++)
    matrix[i] = matrix[i-1] + col;
      
  for (unsigned i = 0; i < row; i++)
    for (unsigned j = 0; j < col; j++)
      matrix[i][j] = other.matrix[i][j];
}

// move constructor
Matrix::Matrix(Matrix&& other) {     
  row = other.row;
  col = other.col;
  matrix = other.matrix;

  other.row = 0;
  other.col = 0;
  other.matrix = nullptr;
}

Matrix::~Matrix() {
  if (matrix)
    delete[] matrix[0];
  delete[] matrix;
}
//--------------------------------------------------------------------------
void Matrix::cleanup() {
  if (matrix)
    delete[] matrix[0];
  delete[] matrix;
}

unsigned Matrix::get_row() const {
  return row;
}

unsigned Matrix::get_col() const {
  return col;
}

double** Matrix::get_matrix() const {
  return matrix;
}

// access to the elements
const double& Matrix::at(unsigned i, unsigned j) const {
  return *(matrix[i*col + j]);
}

double& Matrix::at(unsigned i, unsigned j) {
  return *(matrix[i*col + j]);
}

const double* Matrix::operator[](unsigned i) const {
  return *(matrix + i*col);
}

double* Matrix::operator[](unsigned i) {
  return *(matrix + i*col);
}
//---------------------------------------------------------------------------

Matrix& Matrix::operator=(const Matrix& other) {
  if (this != &other) { // compare ptrs, not objects
    cleanup();
    row = other.row;
    col = other.col;

    matrix = new double* [row];
    matrix[0] = new double [row * col];
    for (unsigned i = 1; i < row; i++)
      matrix[i] = matrix[i-1] + col;
        
    for (unsigned i = 0; i < row; i++)
      for (unsigned j = 0; j < col; j++)
        matrix[i][j] = other.matrix[i][j];
  }

  return *this;
}

// move assignment
Matrix& Matrix::operator=(Matrix&& other) {
  if (this != &other) {
    cleanup();

    row = other.row;
    col = other.col;
    matrix = other.matrix;

    other.row = 0;
    other.col = 0;
    other.matrix = nullptr;
  }

  return *this;
}

Matrix& Matrix::operator+=(const Matrix& other) {
  assert((row == other.row)&&(col == other.col));

  for (unsigned i = 0; i < row; i++)
    for (unsigned j = 0; j < col; j++)
      matrix[i][j] += other.matrix[i][j];

  return *this;
}

Matrix& Matrix::operator-=(const Matrix& other) {
  assert((row == other.row)&&(col == other.col));

  for (unsigned i = 0; i < row; i++)
    for (unsigned j = 0; j < col; j++)
      matrix[i][j] -= other.matrix[i][j];

  return *this;
}

Matrix& Matrix::operator*=(double a) {
  for (unsigned i = 0; i < row; i++)
    for (unsigned j = 0; j < col; j++)
      matrix[i][j] *= a;
     
  return *this;
}

Matrix& Matrix::transpose() {
  Matrix tmp(col, row);
  for (unsigned i = 0; i < col; i++)
    for (unsigned j = 0; j < row; j++)
      tmp.matrix[i][j] = matrix[j][i];

  *this = move(tmp);
  return *this;
}

Matrix& Matrix::inverse() {
  assert(row == col);
  assert(determinant() != 0);

  // create an identity matrix
  Matrix E(row, row);   
  for (unsigned i = 0; i < row; i++)
    for (unsigned j = 0; j < row; j++)
      E.matrix[i][i] = 1.0;
      
  // direct course of Gauss method
  double temp;
  for (unsigned k = 0; k < row; k++) {
    temp = matrix[k][k];
   
    for (unsigned j = 0; j < row; j++) {
      matrix[k][j] /= temp;
      E.matrix[k][j] /= temp;
    }
   
    for (unsigned i = k + 1; i < row; i++) {
      temp = matrix[i][k];
   
      for (unsigned j = 0; j < row; j++) {
        matrix[i][j] -= matrix[k][j] * temp;
        E.matrix[i][j] -= E.matrix[k][j] * temp;
      }
    }
  }     

  // reverse course of Gauss method
  for (unsigned k = row - 1; k > 0; k--)
    for (int i = k - 1; i >= 0; i--) {
      temp = matrix[i][k];
      for (unsigned j = 0; j < row; j++) {
        matrix[i][j] -= matrix[k][j] * temp;
        E.matrix[i][j] -= E.matrix[k][j] * temp;
      }
    }
	
  *this = move(E);
  return *this;
}

void Matrix::set_from_mult(const Matrix& a, const Matrix& b) {
  assert(a.col == b.row);

  Matrix tmp (a.row, b.col);
  for (unsigned i = 0; i < a.row; i++)
    for (unsigned j = 0; j < b.col; j++)
      for (unsigned inner = 0; inner < a.col; inner++)
        tmp.matrix[i][j] += a.matrix[i][inner] * b.matrix[inner][j];
  *this = move(tmp);
}

Matrix& Matrix::multEq (const Matrix& other) {
  Matrix tmp (move(*this));
  set_from_mult(tmp, other);
      
  return *this;
}

Matrix mult(const Matrix& a, const Matrix& b) {
  Matrix tmp;
  tmp.set_from_mult(a, b);

  return tmp;
}
//------------------------------------------------------------------------

void Matrix::print() const {
  for (unsigned i = 0; i < row; i++) {
    for (unsigned j = 0; j < col; j++)
      printf("%.3g ", matrix[i][j]);
    printf("\n");
  }
  printf("\n");
}
//------------------------------------------------------------------------
    
Matrix& Matrix::kroneker(Matrix& other) {
  // two-dimensional array of block matrices
  Matrix** mas = new Matrix* [row];
  mas[0] = new Matrix [row * col];
  for (unsigned i = 1; i < row; i++)
    mas[i] = mas[i-1] + col;

  // filling of block matrices
  for (unsigned i = 0; i < row; i++)
    for (unsigned j = 0; j < col; j++) {
      mas[i][j] = other * matrix[i][j];
      other * (1 / matrix[i][j]);
  }

  // filling of the resulting matrix
  Matrix* tmp = new Matrix(row * other.row, col * other.col);
  for (unsigned i = 0; i < tmp->row; i++)
    for (unsigned j = 0; j < tmp->col; j++)
      tmp->matrix[i][j] = mas[i/row][j/col].matrix[i%other.row][j%other.col];

  *this = *tmp;
  delete tmp;
  delete[] mas[0];
  delete[] mas;

  return *this;
}
    
// Gauss's method
double Matrix::determinant() const {
  assert(row == col);

  Matrix tmp = *this;
  double d, det = 1.0, sign = 1.0, max, *t;
  unsigned max_row;

      
  // lead the matrix to upper triangular form
  for (unsigned i = 0; i < row - 1; i++) { // rows

    // find line k having maximum tmp[i][k]
    max = tmp.matrix[i][i];
    max_row = i;

    for (unsigned k = i; k < row; k++) { // rows
      if (abs(tmp.matrix[k][i]) > abs(max)) {
        max = tmp.matrix[k][i];
        max_row = k; 
      }
    }
    if (max_row > i) { // swap lines 
      if (i != 0) {
        t = tmp.matrix[i];
        tmp.matrix[i] = tmp.matrix[max_row];
        tmp.matrix[max_row] = t;
        sign *= -1.0;
        t = nullptr;
      }
      else {
        double* temp = new double [row];
        for (unsigned n = 0; n < row; n++) {
          temp[n] = tmp.matrix[i][n];
          tmp.matrix[i][n] = tmp.matrix[max_row][n];
          tmp.matrix[max_row][n] = temp[n];
        }
        delete[] temp;
        sign *= -1.0;
      }
    }
        
    for (unsigned j = i + 1; j < row; j++) {
      d = tmp.matrix[j][i] / tmp.matrix[i][i];
      for (unsigned l = 0; l < row; l++)
        tmp.matrix[j][l] -= d * tmp.matrix[i][l];
    }
  }
    
  // multiply the diagonal elements
  for (unsigned i = 0; i < row; i++)
    det *= tmp.matrix[i][i];

  det *= sign;
  return det;
}
//------------------------------------------------------------------------------------

bool operator==(const Matrix& a, const Matrix& b) {
  if((a.get_row() == b.get_row())&&(a.get_col() == b.get_col())) {
    for (unsigned i = 0; i < a.get_row(); i++)
      for (unsigned j = 0; j < a.get_col(); j++)
        if (a.at(i,j) != b.at(i,j))
          return false;
    return true;
  }
  return false;
}

bool operator!=(const Matrix& a, const Matrix& b) {
  return (a == b);
}

Matrix operator+(const Matrix& a, const Matrix& b) {
  assert((a.get_row() == b.get_row())&&(a.get_col() == b.get_col()));
  Matrix tmp (a);
  tmp += b;
  return tmp;
}
    
Matrix operator-(const Matrix& a, const Matrix& b) {
  assert((a.get_row() == b.get_row())&&(a.get_col() == b.get_col()));
  Matrix tmp (a);
  tmp -= b;
  return tmp;
}

Matrix operator*(const Matrix m, double a) {
  Matrix tmp (m);
  tmp *= a;
  return tmp;
}

Matrix operator*(double a, const Matrix m) {
  Matrix tmp (m);
  tmp *= a;
  return tmp;
}
//-------------------------------------------------------------

Matrix& Matrix::add_line() { // for the function solve()
  Matrix tmp(row + 1, col);
  for (unsigned i = 0; i < row; i++)
    for (unsigned j = 0; j < col; j++)
      tmp.matrix[i][j] = matrix[i][j];
    
  *this = move(tmp);

  return *this;
}

void Matrix::correct() { // get rid of inaccuracy
  for (unsigned i = 0; i < row; i++)
    for (unsigned j = 0; j < col; j++)
      if (abs(matrix[i][j]) < 0.00001)
        matrix[i][j] = 0.0;
}

void Matrix::row_decrement() { // for the function fix_A()
  row--;
  delete[] matrix[row-1];
}

void Matrix::fix_A(int& index) { // for the function delete0()
  int i = 0, flag1, flag2, minf, maxf;

  // find two nodes with the same potential
  // as nodes in zero connection
  while (i < row) {
    if (matrix[i][index] == -1)
      flag1 = i;

    if (matrix[i][index] == 1)
      flag2 = i;
    i++;
  }

  // delete a colomn of a zero connection
  if (index == col - 1) {} 
  else 
    for (int i = 0; i < row; i++)
      for (int j = index; j < col - 1; j++)
        matrix[i][j] = matrix[i][j+1];
  col--;

  minf = flag1 < flag2 ? flag1 : flag2;
  maxf = flag1 > flag2 ? flag1 : flag2;

  // delete a line with a second node of two with the same potentials
  for (int j = 0; j < col; j++)
    matrix[minf][j] += matrix[maxf][j];

  if (maxf != row - 1)
    for (int i = maxf; i < row - 1; i++)
      for (int j = 0; j < col; j++)
        matrix[i][j] = matrix[i+1][j];

  delete[] matrix[row-1];
  row--;
}

void Matrix::fix_Y(int& index) { // for the function delete0()
  // delete a colomn and a line of a zero resistanse
  for (int i = index; i < row - 1; i++)
    for (int j = index; j < col - 1; j++)
      matrix[i][j] = matrix[i+1][j+1];
  
  delete[] matrix[row-1];
  row--;
  col--;
}

void Matrix::fix_JEI(int& index) { // for the function delete0()
  // delete a current/voltage in a zero connection
  for (int i = index; i < row - 1; i++)
    matrix[i][0] = matrix[i+1][0];
  row--;
}