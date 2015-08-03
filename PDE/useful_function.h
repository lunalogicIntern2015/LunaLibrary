#pragma once

//#include <iostream>
//#include <fstream>
//#include <string>
#include <vector>
//#include <math.h>

#include <PDE/Matrix/TridiagonalMatrix.h>

//#include <ql/quantlib.hpp> // to use the Matrix class
#include <PDE/Matrix/Matrix.h>


double inverse_sinh(double x);
bool   if_positive(double x);


const double pi_2 = 0.398942280401433;  // =1/pow(pi,0.5) 

//! vector is supposed to be implemented as (n,1) matrix
//! -------------------------------------
//!
//!     Matrix to vector operattion 
//! 
//! -------------------------------------
void copy_matrixColumnToVector(Matrix& matrix, Matrix& vectorCol, int index_col);
void copy_matrixColumnToVector_withoutFirstAndLastElment(Matrix& matrix, Matrix& vectorCol, int index_col );

void copy_matrixRowToVector(Matrix& matrix, Matrix& vectorRow, int index_row);
void copy_matrixRowToVector_withoutFirstAndLastElment(Matrix& matrix, Matrix& vectorRow, int index_row);

//! -------------------------------------
//!
//!     vector to Matrix operattion 
//! 
//! -------------------------------------
void add_vectorToMatrixColumn( Matrix& matrix, Matrix& vectorCol, int index_col);
void copy_vectorToMatrixColumn( Matrix& matrix, Matrix& vectorCol, int index_col);
void add_vectorToMatrixColumn_withoutFirstAndLastElment(Matrix& matrix, Matrix& vectorCol, int index_col);
void copy_vectorToMatrixColumn_withoutFirstAndLastElment(Matrix& matrix, Matrix& vectorCol, int index_col);
void add_vectorToMatrixRow(Matrix& matrix, Matrix& vectorRow, int index_row);
void copy_vectorToMatrixRow(Matrix& matrix, Matrix& vectorRow, int index_row);
void add_vectorToMatrixRow_withoutFirstAndLastElment(Matrix& matrix, Matrix& vectorRow, int index_row);
void copy_vectorToMatrixRow_withoutFirstAndLastElment(Matrix& matrix, Matrix& vectorRow, int index_row);

//! -------------------------------------
//!
//!     Matrix to Matrix operation
//! 
//! -------------------------------------
void add_smallMatrix_to_bigMatrix_withOutBoundary(Matrix& m_small, Matrix& m_big);
void copy_smallMatrix_to_bigMatrix_withOutBoundary(Matrix& m_small, Matrix& m_big);



//! make matrix-vector equals to zero
void init_to_zero_colMatrix(Matrix& m);


double max(double x, double y);
double digital(double x, double y);

void print(Matrix& m);
void print(Matrix& m, std::string s);
void print(TridiagonalMatrix& m, std::string s);
void print(const std::vector<double>& m, const std::string s);
void print_MatrixToFile(std::string& file_name, Matrix& m);
void print_result(std::string& file_name, std::vector<double>& x_grid, std::vector<double>& y_grid, Matrix& m);


// Inverse Cumulative Normale
double Inverse_CN(double u);
double NormalDensity(double x);
