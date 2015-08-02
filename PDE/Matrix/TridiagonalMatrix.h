//#pragma once
//
//#include <vector>
//#include <iostream>
//
////#include <ql/quantlib.hpp> // to use the Matrix class
////using namespace QuantLib;
//
//#include "Matrix.h"
//
//class TridiagonalMatrix
//{
//// If faut d'etre private !!! To be changed ...
//public:
//	int sizeDiagonal_;
//	int sizeMatrix_;
//
//	//! TODO: exchang 0 and 2, for it is conter-intuitive
//	//! size = 3*n:
//	//! row[0] = upper(i-j=-1), the last element is added and it equals to 0
//	//! row[1] = diagonal(i-j=0),
//	//! row[2] = lower(i-j=1), the first elemen is added and it equals to 0
//	Matrix tridiagonalMatrix_;  
//
//public:
//	//! constructor & deconstructor
//    TridiagonalMatrix(int sizeMatrix);
//	virtual ~TridiagonalMatrix(){};
//
//	Matrix* get_matrixPointer() {return &tridiagonalMatrix_;}
//	int get_sizeDiagonal()  const {return sizeDiagonal_;}
//	int get_sizeMatrix()    const {return sizeMatrix_;}
//	void multipleColVector(const Matrix& v_ColVectorToMultiple, Matrix& v_result) const; 
//	Matrix solve_linear_equation(Matrix& r) const;
//	void   set_Matrix(Matrix& m) const ;  // OK
//	double get_left_up_element()const{return tridiagonalMatrix_[2][0];}   // OK
//	double get_right_bottom_element()const{return tridiagonalMatrix_[0][sizeMatrix_-1];} // OK
//
//	void print() const;
//};
//


#pragma once

#include <vector>
#include <iostream>
#include <string>

//#include <ql/quantlib.hpp> // to use the Matrix class
//using namespace QuantLib;

#include <PDE/Matrix/Matrix.h>

class TridiagonalMatrix
{
// If faut d'etre private !!! To be changed ...
public:
	int sizeDiagonal_;
	int sizeMatrix_;


	//! size = 3*n:
	//! row[0] = upper(i-j=-1), the last element is added and it equals to 0
	//! row[1] = diagonal(i-j=0),
	//! row[2] = lower(i-j=1), the first elemen is added and it equals to 0
	Matrix tridiagonalMatrix_;  
	std::string special_type_; // only take care of "id"

public:
	//! constructor & deconstructor
    TridiagonalMatrix(int sizeMatrix);
	TridiagonalMatrix(int sizeMatrix, const std::string & special_type);
	virtual ~TridiagonalMatrix(){};

	Matrix* get_matrixPointer() {return &tridiagonalMatrix_;}
	int get_sizeDiagonal() const {return sizeDiagonal_;}
	int get_sizeMatrix()const {return sizeMatrix_;}
	double get_left_up_element()const{return tridiagonalMatrix_[2][0];}   // OK
	double get_right_bottom_element()const{return tridiagonalMatrix_[0][sizeMatrix_-1];} // OK

	TridiagonalMatrix operator*(double s) const;  //! not efficient at all ! 
	void multipleScalar(const double val); 
	void addIdMatrix(); 

	void multipleColVector(const Matrix& v_ColVectorToMultiple, Matrix& v_result) const; 
	void multipleColVector(Matrix& v); 

	void solve_linear_equation(const Matrix& r, Matrix& u) const;  // u is the result matrix
	void down_solve_linear_equation(const Matrix& r, Matrix& u) const;
	void up_solve_linear_equation(const Matrix& r,Matrix& u) const;

	void set_Matrix(Matrix& m);
	void print(std::string name = "") const;
	std::vector<double> eigenvalues();

	bool if_id_TriMatrix() const;
};

