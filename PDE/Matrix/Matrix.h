#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <cmath>

using namespace std;

//class Vector;

class Matrix
{
protected:
	int m_rows;
	int m_cols;
	double** m_matrix;

public:
	//!  wrap operator [i][j], it has the same speed as (,) or direct access ! 
	// When you have a decent compiler and if you judiciously use inlining, the compiler should optimize away the TEMPORARY objects.
	// In other words, the operator[]-approach above will hopefully not be slower than what it would have been if you had directly called Matrix::operator()( row,  col) in the first place. Of course you could have made your life simpler and avoided most of the above work by directly calling Matrix::operator()( row,  col) in the first place. 
	// So you might as well directly call Matrix::operator()( row,  col) in the first place. 
	class Row 
	{
	public:
		Matrix  & m_matrix;
		int     m_row;
		Row(Matrix & matrix, int row): m_matrix(matrix), m_row(row){};
		double& operator [] (int col) {return m_matrix(m_row,col);}
	};
	Row operator [](int row)  
	{
	    return Row(*this,row);
	}

	////! const version (TODO, not very clear for me ... )
	class ConstRow 
	{
	public:
		Matrix const & m_matrix;
		int     m_row;
		ConstRow(Matrix const & matrix, int row): m_matrix(matrix), m_row(row){};
		double operator [] (int col) const {return m_matrix(m_row,col);}
	};
	ConstRow operator [](int row) const
	{
	    return ConstRow(*this,row);
	}


	//! ---- ---- ---- ---- Matrix ---- ---- ---- ---- 
	//! constructor 
	Matrix(int rows, int cols, double val = 0);  // creat a matrix, each element is zeros
	Matrix(double** matrix, int rows, int cols);

	Matrix(int rows, int cols, std::string & matrix_type); // generate special type of matrix ... 
	Matrix(vector<double> & x, bool if_column_matrix);
	//! copy constructor and assignment constructor (no pointer attribute, use the default copy/assignement constructor) 
    Matrix(const Matrix& m);
	Matrix& operator = (const Matrix& m);
	//! destructor
	virtual ~Matrix();


	//! get functions ... 
	int rows() const;
	int cols() const;
	double** get_matrix();

	//! functionalilities
	double& operator()(int i, int j);
	double operator()(int i, int j) const;
    Matrix operator+ (const Matrix& m) const;
	Matrix operator- (const Matrix& m) const;
	void   operator+=(const Matrix& m);
	void   operator-=(const Matrix& m);

	Matrix operator*(const Matrix& m) const;
	//Vector operator*(Vector& m);
	Matrix operator~() const; // transpose a matrix  ... 
    
    //! operation with scalar ... 
	void operator+=( double val );
	void operator-=( double val );
	void operator*=( double val );
	void operator/=( double val );

	void set_element_to_value(double v);

	//! ---- ---- ---- ---- Square Matrix ---- ---- ---- ---- 
	bool   if_square_matrix() const;
	double determinant__(double ** m, int n) const;
    double determinant() const;  
	Matrix inverse() const;
	vector<double> eigenvalues() const;

	//! ---- ---- ---- ---- Symetric Matrix ---- ---- ---- ---- 
	bool   if_symetric_matrix() const;
	Matrix Cholesky(Matrix& m_op) const;

	//! ---- ---- ---- ---- Vector ---- ---- ---- ---- 
	bool   if_row_vector() const;
	bool   if_col_vector() const;
	double inner_product(const Matrix& v) const; // about two column vectors

	//! print
	void print(std::string matrix_name = "");
};


//! TODO: 
// to give a structure of the Matrix --> square matrix --> symetric matrix
// ans symetric matrix can be saved using half memeory ...