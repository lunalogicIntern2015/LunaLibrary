#include "Matrix.h"
#include <stdio.h>
#include <stdlib.h>

//! for inverse matrix
#include "svd.h"
//! for eigenvalue
#include <PDE/NR/code/eigen_unsym.h>
//#include "nr3.h"

const double seille = 0.0000000001;


// ---------------------------------------------------------------------
//
//                          MATRIX 
//   
// ---------------------------------------------------------------------

Matrix :: Matrix (int rows, int cols, double val)
        :  
		   m_rows(rows),
		   m_cols(cols)
{
     m_matrix = new double*[rows];
	 for(int i=0; i<m_rows; ++i)
	 {
	     m_matrix[i] = new double[cols];
		 for(int j=0; j<m_cols; ++j)
		 {
		     m_matrix[i][j] = val;
		 }
	 }
}

Matrix :: Matrix(int rows, int cols, std::string & matrix_type)
        :  
		   m_rows(rows),
		   m_cols(cols)
{
    if(matrix_type == "rand")
	{ 
	     m_matrix = new double*[rows];
		 for(int i=0; i<m_rows; ++i)
		 {
			 m_matrix[i] = new double[cols];
			 for(int j=0; j<m_cols; ++j)
			 {
				 m_matrix[i][j] = (rand()/(double)RAND_MAX - 0.5)*2;  // random number between [-1,1]
			 }
		 }
	}
	else
	{
	    throw ("Error in constructor of Matrix, string matrix_type not defined ... ");
	}
}


Matrix :: Matrix(double**  matrix, int rows, int cols)
			: m_rows(rows),
		      m_cols(cols)
{
	//! create matrix 
	m_matrix = new double*[rows];
	for(int i=0; i<m_rows; ++i)
	{
	    m_matrix[i] = new double[cols];
	}

	//! copy matrix 
    for(int i=0; i<rows; ++i)
	{
	    for(int j=0; j<cols; ++j)
		{
		    m_matrix[i][j] = matrix[i][j];
		}
	}
}


Matrix :: Matrix(vector<double> & x, bool if_column_matrix)
{
	 if(if_column_matrix = true) // col matrix
	 {
	     m_rows = (int)x.size();
		 m_cols = 1;
	 }
	 else  // row matrix
	 {
	     m_rows = 1;
		 m_cols = (int)x.size();	 
	 }

     m_matrix = new double*[m_rows];
	 for(int i=0; i<m_rows; ++i)
	 {
	     m_matrix[i] = new double[m_cols];
		 for(int j=0; j<m_cols; ++j)
		 {
			 if(if_column_matrix == true)
			 {
				m_matrix[i][j] = x[i];
			 }
			 else
			 {
			    m_matrix[i][j] = x[j];
			 }
		 }
	 }
}


Matrix :: Matrix(const Matrix& m)  // old object not existed 
{
	m_rows  = m.m_rows;
	m_cols  = m.m_cols;
    m_matrix = new double*[m_rows];
	for(int i=0; i<m_rows; ++i)
	{
	    m_matrix[i] = new double[m_cols];
		for(int j=0; j<m_cols; ++j)
		{
		    m_matrix[i][j] = m.m_matrix[i][j];
		}
	}
}


Matrix& Matrix :: operator = (const Matrix& m)  // old object existed 
{
    if(this == &m)
	{
	    return *this;
	}
	else
	{
		//! if or not clear the old memeory
		if(m.rows()!= this->rows() || m.cols() != this->cols() )
		{
			//! delete old memory
			for(int i=0; i<m_rows; ++i)
			{
			   delete[] m_matrix[i];
			}
			delete[] m_matrix;

			//! create new memory 
			m_rows  = m.m_rows;
			m_cols  = m.m_cols;
			m_matrix = new double*[m_rows];
			for(int i=0; i<m_rows; ++i)
			{
				m_matrix[i] = new double[m_cols];
				//for(int j=0; j<m_cols; ++j)
				//{
				//	m_matrix[i][j] = m.m_matrix[i][j];
				//}
			}
		}
		//! copy
		for(int i=0; i<m_rows; ++i)
		{
			for(int j=0; j<m_cols; ++j)
			{
				m_matrix[i][j] = m.m_matrix[i][j];
			}
		}

		return *this;
	}
}


Matrix :: ~Matrix()
{     
	 for(int i=0; i<m_rows; ++i)
	 {
	     delete[] m_matrix[i];
	 }
	 delete[] m_matrix;
}


int Matrix :: rows() const 
{
    return m_rows;
}

int Matrix :: cols() const
{
    return m_cols;
}

//vector<vector<double>> Matrix ::get_matrix()
double** Matrix ::get_matrix()
{
    return m_matrix;
}

double& Matrix :: operator()(int i, int j) 
{
	if(0<=i && i<m_rows && 0<=j && j<m_cols)
	{
		return m_matrix[i][j];
	}
	else
	{
		throw ("Error in function double& Matrix :: operator()(int i, int j): index non valide : i,j");
	}
}

double Matrix :: operator()(int i, int j) const
{
	if(0<=i && i<m_rows && 0<=j && j<m_cols)
	{
		return m_matrix[i][j];
	}
	else
	{
		throw ("Error in function double& Matrix :: operator()(int i, int j): index non valide : i,j");
	}
}


Matrix Matrix :: operator+(const Matrix& m) const
{
	if(this->m_rows == m.m_rows && this->m_cols == m.m_cols)
	{
         Matrix m2(*this);
		 for(int i=0; i<m_rows; ++i)
		 {
		     for(int j=0; j<m_cols; ++j)
			 {
			     m2.m_matrix[i][j] += m.m_matrix[i][j];
			 }
		 }
		 return m2;
	}
	else
	{
	    throw ("Error from Matrix Matrix :: operator+(Matrix& m), matrix size dismatch ");
	}
}

void Matrix ::  operator+=(const Matrix& m)
{
	if(this->m_rows == m.m_rows && this->m_cols == m.m_cols)
	{
		 for(int i=0; i<m_rows; ++i)
		 {
		     for(int j=0; j<m_cols; ++j)
			 {
			     m_matrix[i][j] += m.m_matrix[i][j];
			 }
		 }
	}
	else
	{
	    throw ("Error from Matrix Matrix :: operator+(Matrix& m), matrix size dismatch ");
	}
}
	

Matrix Matrix :: operator-( const Matrix& m) const 
{
	if(this->m_rows == m.m_rows && this->m_cols == m.m_cols)
	{
         Matrix m2(*this);
		 for(int i=0; i<m_rows; ++i)
		 {
		     for(int j=0; j<m_cols; ++j)
			 {
			     m2.m_matrix[i][j] -= m.m_matrix[i][j];
			 }
		 }
		 return m2;
	}
	else
	{
	    throw ("Error from void Matrix :: operator+(Matrix m), matrix size dismatch ");
	}
}

void  Matrix :: operator-=(const Matrix& m)
{
	if(this->m_rows == m.m_rows && this->m_cols == m.m_cols)
	{
		for(int i=0; i<m_rows; ++i)
		 {
		     for(int j=0; j<m_cols; ++j)
			 {
			     m_matrix[i][j] -= m.m_matrix[i][j];
			 }
		 }
	}
	else
	{
	    throw ("Error from void Matrix :: operator+(Matrix m), matrix size dismatch ");
	}
}


Matrix Matrix :: operator*(const Matrix& m) const 
{
	if(this->m_cols == m.m_rows)
	{
		 Matrix m2(this->m_rows, m.m_cols);
  		 for(int i=0; i<m2.m_rows; ++i)
		 {
		     for(int j=0; j<m2.m_cols; ++j)
			 {
				 for(int k=0; k<m.m_rows; ++k)
				 {
					m2.m_matrix[i][j] += this->m_matrix[i][k] * m.m_matrix[k][j];   
				 }
			 }
		 }
		 return m2;
	}
	else
	{
	    throw ("Error from Matrix Matrix :: operator*(Matrix& m) matrix size dismatch ");
	}
}



//! operator with scalar 
void Matrix :: operator+=( double val)
{
	 for(int i=0; i<m_rows; ++i)
	 {
	     for(int j=0; j<m_cols; ++j)
		 {
		     m_matrix[i][j] += val;
		 }
	 }
}

void Matrix :: operator-=( double val)
{
	 for(int i=0; i<m_rows; ++i)
	 {
	     for(int j=0; j<m_cols; ++j)
		 {
		     m_matrix[i][j] -= val;
		 }
	 }
}

void Matrix :: operator*=( double val)
{
	 for(int i=0; i<m_rows; ++i)
	 {
	     for(int j=0; j<m_cols; ++j)
		 {
		     m_matrix[i][j] *= val;
		 }
	 }
}

void Matrix :: operator/=( double val)
{
	if (val == 0.0)
	{
	    throw ("Error in function Matrix :: operator/=( double val), devided by zero ! ");
	}
	 for(int i=0; i<m_rows; ++i)
	 {
	     for(int j=0; j<m_cols; ++j)
		 {
		     m_matrix[i][j] /= val;
		 }
	 }
}


void Matrix::set_element_to_value(double v)
{
	 for(int i=0; i<m_rows; ++i)
	 {
	     for(int j=0; j<m_cols; ++j)
		 {
		     m_matrix[i][j] = v;
		 }
	 }
}


//Vector Matrix :: operator*(Vector& vec)
//{
//    if(this->m_cols == vec.size() && vec.is_col_vector())
//	{
//        vector<double> v = vec.get_v();
//	    vector<double> v2(this->m_rows, 0.0);
//		for(int i=0; i<m_rows; ++i)
//		{
//			for(int k=0; k<m_cols; ++k)
//			{
//				v2[i] += m_matrix[i][k] * v[k];
//			}
//		}
//		return Vector(v2,true);   // return a column vector ...
//	}
//	else
//	{
//		throw ("Error from Vector operator*(Vector& m), matrix size dismatch ");
//	}
//}


Matrix Matrix :: operator~() const  //! cannot overload bi-operator but only single-operator ...
{
	 Matrix m2(m_cols, m_rows);
	 for(int i=0; i<m_cols; ++i)
	 {
	     for(int j=0; j<m_rows; ++j)
		 {
			m2.m_matrix[i][j] += this->m_matrix[j][i];
		 }
	 }
	 return m2;
}



void Matrix :: print(std::string matrix_name)
{
	cout << " ---- matrix " << matrix_name << "    ,size = ( " << m_rows << "," << m_cols << " )  ---- " << endl;
	for(int i=0; i<m_rows; ++i)
	{
		cout << "| " ;
		for(int j=0; j<m_cols; ++j)
		{
		    cout << m_matrix[i][j]<< " ";
		}
		cout << "|"<< endl;
	}
}






// ---------------------------------------------------------------------
//
//                         SQUARE MATRIX 
//   
// ---------------------------------------------------------------------
bool  Matrix :: if_square_matrix() const 
{
    if(m_rows == m_cols)
	{
	    return true;
	}
	else
	{
	    return false;
	}
}

double Matrix :: determinant__( double ** m, int n) const 
{
	//if(if_square_matrix() != true)
	//{
	//	throw Exception("Error in Matrix :: determinant__(...): not square matrix.");
	//}

	if(n==1)
	{
		return m[0][0];
	}
	else if(n==2)
	{
		return m[0][0] * m[1][1] - m[1][0] * m[0][1];
	}
	else
	{
		double result    = 0;
		int    sign      = 1;

		//! create sub_matrix: (n-1, n-1)
		double** sub_m   = new double*[n-1];
		for(int i=0; i<n-1; ++i)
		{
			sub_m[i] = new double[n-1];
		}

		//! resursive calculate determinant ...
		//! inverse the iteration, algorithm will be more efficient ... 
		for(int j=0; j<n; ++j)
		{
		   //! initialize sub_matrix
		   int k_sub = 0;

		   for(int k=1; k<n; ++k)  // n-1 rows
		   {
			   k_sub = k-1;

			   bool flag = false;
			   int l_sub = 0;
			   for(int l=0; l<n; ++l)  // n-1 cols 
			   {
				   if(!flag)
				   {
					   l_sub = l;
				   }
				   else
				   {
				       l_sub = l-1;
				   }
			       if(l==j)
				   {
					   flag = true;
					   continue;
				   }
                   sub_m[k_sub][l_sub] = m[k][l];
			   }
		   }
		   result +=  sign*m_matrix[0][j] * determinant__(sub_m, n-1);	   
		   sign   *= (-1);
		}

		//! delete sub_matrix
		for(int i=0; i<n-1; ++i)
		{
			delete[] sub_m[i];
		}
		delete[] sub_m;

		return result;
	}
}

double Matrix  :: determinant() const //! there can be 2 interpolatation: 1. deract, 2 recursive ..
{
	if(if_square_matrix() != true)
	{
		throw ("Error in Matrix :: determinant(): not square matrix.");
	}

    int n = m_rows;
	double** m   = new double*[n];
	for(int i=0; i<n; ++i)
	{
		m[i] = new double[n];
	}

	for(int i=0; i<n; ++i)
	{
	    for(int j=0; j<n; ++j)
		{
		    m[i][j] = m_matrix[i][j];
		}
	}

	double result = determinant__(m,n);

	for(int i=0; i<n; ++i)
	{
		delete[] m[i];
	}
	delete[] m;

	return result;
}

//Matrix Matrix  :: inverse() const //! don't need to copy it again, try to make it more efficient ... 
//{
//	if(if_square_matrix() != true)
//	{
//		throw ("Error in Matrix :: inverse(): not square matrix.");
//	}
//
//    int n = m_rows;
//	double** a  = new double*[n];
//	for(int i=0; i<n; ++i)
//	{
//		a[i] = new double[n]; 
//	}
//
//	for(int i=0; i<n; ++i)
//	{
//		for(int j=0; j<n; ++j)
//		{
//			a[i][j] = m_matrix[i][j];
//		}
//	}
//
//	svd_inverse_matrix(a,n);
//
//	Matrix result(a,n,n);
//	
//	for(int i=0; i<n; ++i)
//	{
//		delete[] a[i];
//	}
//	delete[] a;
//
//	return result;
//}
//
//


Matrix Matrix  :: inverse() const //! don't need to copy it again, try to make it more efficient ... 
{
	if(if_square_matrix() != true)
	{
		throw ("Error in Matrix :: inverse(): not square matrix.");
	}

    int n = m_rows;
	double** a  = new double*[n];
	for(int i=0; i<n; ++i)
	{
		a[i] = new double[n]; 
	}

	for(int i=0; i<n; ++i)
	{
		for(int j=0; j<n; ++j)
		{
			a[i][j] = m_matrix[i][j];
		}
	}

	svd_inverse_matrix(a,n);

	Matrix result(a,n,n);
	
	for(int i=0; i<n; ++i)
	{
		delete[] a[i];
	}
	delete[] a;

	return result;
}


vector<double> Matrix :: eigenvalues() const
{
	int n = m_rows;
    MatDoub aa(n,n,0.0);
	for(int i=0; i<n; ++i)
	{
		for(int j=0; j<n; ++j)
		{
		    aa[i][j] = m_matrix[i][j];
		}
	}

	Unsymmeig eigenvalue_solver(aa);

	std::vector<double> eigenvalue(n,0.0);
	for(int i=0; i<n; ++i)
	{    
		double real = eigenvalue_solver.wri[i].real();
		double comp = eigenvalue_solver.wri[i].imag();
		if(comp != 0.0)
		{
			throw ("Error in function Matrix :: eigenvalues(): eigenvalue is a complex number!");
		}
		eigenvalue[i] = real;
	}
	return eigenvalue;
}



// ---------------------------------------------------------------------
//
//                         SYMETRIC MATRIX 
//   
// ---------------------------------------------------------------------
bool  Matrix  :: if_symetric_matrix() const 
{
	if(if_square_matrix() == true)
	{
		for(int i=0; i<m_rows; ++i)
		{
			for(int j=0; j<i; ++j)
			{
			    if( fabs(m_matrix[i][j] - m_matrix[j][i]) > seille) 
				{
				    return false;
				}
			}
		}
		return true;
	}
	else
	{
	    return false;
	}
}


Matrix Matrix :: Cholesky (Matrix& m_op) const // output a lower triangular matrix ... --? m_op[k][l] : k>=l
{
	//! initialize m_op to a 0-matrix ...
	int n = m_op.rows();
	for(int i=0; i<n; ++i)
	{
		for(int j=0; j<n; ++j)
		{
		     m_op.m_matrix[i][j] = 0.0;
		}
	}
    
	for(int i=0; i<n; ++i)
	{
		for(int j=i; j<n; ++j)
		{
			 if(i==j)
			 {
				 if(m_matrix[i][i] >= 0 ) // ?? not sure if ==0 ???? why should not a[i][i] < 0 ? What's the condition for cholesky decompostiion ??? 
					m_op.m_matrix[i][i] = sqrt(m_matrix[i][i]);
			 }
			 else
			 {
				 double temp = 0.0;
				 for(int k=0; k<i; ++k)
				 {
					 temp += m_op.m_matrix[i][k] * m_op.m_matrix[j][k]; 
				 }
				 //temp + m_op[i][i] * m_op[j][i] = m_matrix[i][j];
				 if(m_op.m_matrix[i][i] == 0.0)
				 {
					 throw ("Error in function Symetric_Matrix :: Cholesky(...), m_op[i][i] == 0.0");
				 }
				 else
				 {
					m_op.m_matrix[j][i] = (m_matrix[i][j] - temp)/m_op.m_matrix[i][i];
				 }
			 }
		}
	}
	throw ("Error in  Matrix :: Cholesky(), the algorith is not tested, it doesn't return a value, BUG!");
}


// ---------------------------------------------------------------------
//
//                         Vector
//   
// ---------------------------------------------------------------------
bool Matrix ::   if_row_vector() const 
{
	if(m_cols ==1)
		return true;
	else 
		return false;
}

bool Matrix ::   if_col_vector() const 
{
	if(m_cols ==1) 
		return true; 
	else 
		return false;
}
double Matrix :: inner_product(const Matrix & v) const// about two vectors, don't care if it is column or row vectors.
{
    if(!this->if_col_vector() || ! v.if_col_vector() || this->rows()!=v.rows() )
	{
		throw ("Error in Matrix :: inner_product(...), Direction_Parameters are not column vectors or matrix size don't match.");
	}
	double sum = 0.0;
	for(int i=0; i<v.rows(); ++i)
	{
	    sum += v[i][0]*(*this)[i][0];
	}
	return sum;
}