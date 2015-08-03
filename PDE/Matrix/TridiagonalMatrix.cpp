//#include "TridiagonalMatrix.h"
//#include <vector>
//
//	TridiagonalMatrix::TridiagonalMatrix(int sizeMatrix)
//		:
//		 sizeDiagonal_(3),
//         sizeMatrix_(sizeMatrix),
//	     tridiagonalMatrix_(3,sizeMatrix_,0)
//	{ 
//	   //sizeDiagonal_ = 3;
//       //sizeMatrix_   = sizeMatrix;
//	   //tridiagonalMatrix_ = Matrix(3,sizeMatrix_,0);
//	}
//	
//	void TridiagonalMatrix:: multipleColVector(const Matrix& v_ColVectorToMultiple, Matrix& v_result) const
//	{
//		if( this->tridiagonalMatrix_.cols()!= v_result.rows() 
//			|| v_result.rows()!=v_ColVectorToMultiple.rows() 
//			|| v_ColVectorToMultiple.cols() != 1
//		  )
//		  std::cout << " Error, wrong size, in function void TridiagonalMatrix:: multipleColVector(...) " << std::endl; // to be changed to throw Exeption
//		
//		for( int nRow=0; nRow<v_result.rows(); nRow++)
//		{
//			v_result[nRow][0] = 0;
//			if(nRow==0)
//			{
//			   v_result[nRow][0] += tridiagonalMatrix_[0][nRow]*v_ColVectorToMultiple[nRow+1][0] ;
//			   v_result[nRow][0] += tridiagonalMatrix_[1][nRow]*v_ColVectorToMultiple[nRow][0];
//			}
//			else if (nRow<v_result.rows()-1)
//			{
//			   v_result[nRow][0] += tridiagonalMatrix_[0][nRow]*v_ColVectorToMultiple[nRow+1][0];
//			   v_result[nRow][0] += tridiagonalMatrix_[1][nRow]*v_ColVectorToMultiple[nRow][0];
//			   v_result[nRow][0] += tridiagonalMatrix_[2][nRow]*v_ColVectorToMultiple[nRow-1][0];
//
//			}
//			else
//			{
//			   v_result[nRow][0] += tridiagonalMatrix_[1][nRow]*v_ColVectorToMultiple[nRow][0];
//			   v_result[nRow][0] += tridiagonalMatrix_[2][nRow]*v_ColVectorToMultiple[nRow-1][0];
//			}
//		}
//	}
//
//	//! A*x=r, A matrix, x,b column vector
//	//! A,b known, and the solution is "x"
//	Matrix TridiagonalMatrix :: solve_linear_equation(Matrix& r) const
//	{
//		if( sizeMatrix_ != r.rows() || !(r.if_col_vector()) )
//		{
//			throw ("Error in Matrix :: inverse(): not square matrix.");
//		}
//
//		Matrix u(sizeMatrix_,1);
//
//		int n = sizeMatrix_;
//		double bet;
//		std::vector<double> gam(n);
//		//if(b[0] == 0.0)
//		if(tridiagonalMatrix_[0][0] == 0.0)
//			throw("Error 1 in function Matrix::solve_linear_equation(): 0 in the tridiagonal");
//		//u[0] = r[0] / (bet=b[0]);
//		u[0][0] = r[0][0] / (bet=tridiagonalMatrix_[1][0]);
//
//		for(int j=1; j<n; ++j)
//		{
//		    //gam[j] = c[j-1]/bet;
//			gam[j] = tridiagonalMatrix_[0][j-1]/bet;
//			//bet    = b[j] - a[j]*gam[j];
//			bet    = tridiagonalMatrix_[1][j] - tridiagonalMatrix_[2][j]*gam[j];
//			if(bet == 0.0)
//			{
//				throw("Error 2 in function Matrix::solve_linear_equation(): 0 in the tridiagonal");
//			}
//			//u[j] = (r[j] - a[j]*u[j-1]) /bet;
//			u[j][0] = (r[j][0] - tridiagonalMatrix_[2][j]*u[j-1][0]) /bet;
//		}
//		for(int j=(n-2); j>=0; --j)
//		{
//		    u[j][0] -= gam[j+1]*u[j+1][0]; 
//		}
//		return u;
//	}
//
//
//
//	void TridiagonalMatrix :: set_Matrix(Matrix& m) const
//	{
//		if(m.rows() != m.cols() || m.rows() != this->get_sizeMatrix())
//			std::cout << "matrix size dismatch in TridiagonalMatrix :: set_Matrix(Matrix& m)" << std::endl;
//
//	    for(int i=0; i<sizeMatrix_; i++)
//		{
//		    for(int j=0; j<sizeMatrix_; j++)
//			{
//				if(i-j == -1) // i<j   f_0
//					m[i][j] = tridiagonalMatrix_[0][i];
//				else if(i-j == 0) // i = j  f_1
//					m[i][j] = tridiagonalMatrix_[1][i];
//				else if(i-j == 1) // i > j  f_2
//					m[i][j] = tridiagonalMatrix_[2][i];
//				else
//					m[i][j] = 0;
//			}
//		}
//	}
//
//	void TridiagonalMatrix:: print() const 
//	{
//		std::cout << "TridiagonalMatrix of size (" << sizeMatrix_ << " , " << sizeMatrix_ << ") :" << std::endl;
//	    for(int i=0; i<sizeMatrix_; i++)
//		{
//		    for(int j=0; j<sizeMatrix_; j++)
//			{
//				if(i-j == -1)
//					std::cout << tridiagonalMatrix_[0][i]<< "	";
//				else if(i-j == 0)
//					std::cout << tridiagonalMatrix_[1][i] << "	";
//				else if(i-j == 1)
//					std::cout << tridiagonalMatrix_[2][i] << "	";
//				else
//					std::cout << 0 << "	";
//			}
//			std::cout << std::endl;
//		}
//	}

#include "TridiagonalMatrix.h"
TridiagonalMatrix::TridiagonalMatrix(int sizeMatrix)
	: sizeDiagonal_(3),
      tridiagonalMatrix_(3,sizeMatrix,0),
	  sizeMatrix_(sizeMatrix)
{ 
   //sizeDiagonal_ = 3;
   //sizeMatrix_   = sizeMatrix;
   //tridiagonalMatrix_ = Matrix(3,sizeMatrix_,0);
}

TridiagonalMatrix::TridiagonalMatrix(int sizeMatrix, const std::string & special_type)
	: sizeDiagonal_(3),
      tridiagonalMatrix_(3,sizeMatrix,0),
	  sizeMatrix_(sizeMatrix),
	  special_type_(special_type) // only take care of "id"
{
	if(special_type == "id")
	{
		for(int i=0; i<sizeMatrix_; ++i)
		{
			tridiagonalMatrix_[1][i] = 1.0;
		}
	}
}


TridiagonalMatrix TridiagonalMatrix :: operator*(double s) const
{
    TridiagonalMatrix tri_m(*this); // copy corrent obj
	tri_m.multipleScalar(s);
	return tri_m;
	
}

void TridiagonalMatrix :: multipleScalar(const double val)
{
	for(int i=0; i<tridiagonalMatrix_.rows(); ++i)
	{
		for(int j=0; j<tridiagonalMatrix_.cols(); ++j)
			tridiagonalMatrix_[i][j] *= val;
	}
}

void TridiagonalMatrix :: addIdMatrix()
{
    	//! row[1] = diagonal(i-j=0),
	for(int i=0; i<tridiagonalMatrix_.cols(); ++i)
	{
	    tridiagonalMatrix_[1][i] += 1;
	}
}

void TridiagonalMatrix:: multipleColVector(const Matrix& v_ColVectorToMultiple, Matrix& v_result) const
{
	if( this->tridiagonalMatrix_.cols()!= v_result.rows() 
		|| v_result.rows()!= v_ColVectorToMultiple.rows() 
		|| v_ColVectorToMultiple.cols() != 1
	  )
	  std::cout << " Error, wrong size, in function void TridiagonalMatrix:: multipleColVector(...) " << std::endl; // to be changed to throw Exeption
	
	for(int nRow=0; nRow<v_result.rows(); nRow++)
	{
		v_result[nRow][0] = 0;
		if(nRow==0)
		{
		   v_result[nRow][0] += tridiagonalMatrix_[0][nRow]*v_ColVectorToMultiple[nRow+1][0] ;
		   v_result[nRow][0] += tridiagonalMatrix_[1][nRow]*v_ColVectorToMultiple[nRow][0];
		}
		else if (nRow<v_result.rows()-1)
		{
		   v_result[nRow][0] += tridiagonalMatrix_[0][nRow]*v_ColVectorToMultiple[nRow+1][0];
		   v_result[nRow][0] += tridiagonalMatrix_[1][nRow]*v_ColVectorToMultiple[nRow][0];
		   v_result[nRow][0] += tridiagonalMatrix_[2][nRow]*v_ColVectorToMultiple[nRow-1][0];

		}
		else
		{
		   v_result[nRow][0] += tridiagonalMatrix_[1][nRow]*v_ColVectorToMultiple[nRow][0];
		   v_result[nRow][0] += tridiagonalMatrix_[2][nRow]*v_ColVectorToMultiple[nRow-1][0];
		}
	}
}

void TridiagonalMatrix:: multipleColVector(Matrix& v) // v is a column vector
{
    Matrix v_ColVectorToMultiple(v);  // a copy of v

	if( this->tridiagonalMatrix_.cols()!= v.rows() 
		|| v.rows()!=v_ColVectorToMultiple.rows() 
		|| v_ColVectorToMultiple.cols() != 1
	  )
	  std::cout << " Error, wrong size, in function void TridiagonalMatrix:: multipleColVector(...) " << std::endl; // to be changed to throw Exeption
	
	for( int nRow=0; nRow<v_ColVectorToMultiple.rows(); nRow++)
	{
		v[nRow][0] = 0;
		if(nRow==0)
		{
		   v[nRow][0] += tridiagonalMatrix_[0][nRow]*v_ColVectorToMultiple[nRow+1][0] ;
		   v[nRow][0] += tridiagonalMatrix_[1][nRow]*v_ColVectorToMultiple[nRow][0];
		}
		else if (nRow<v_ColVectorToMultiple.rows()-1)
		{
		   v[nRow][0] += tridiagonalMatrix_[0][nRow]*v_ColVectorToMultiple[nRow+1][0];
		   v[nRow][0] += tridiagonalMatrix_[1][nRow]*v_ColVectorToMultiple[nRow][0];
		   v[nRow][0] += tridiagonalMatrix_[2][nRow]*v_ColVectorToMultiple[nRow-1][0];

		}
		else
		{
		   v[nRow][0] += tridiagonalMatrix_[1][nRow]*v_ColVectorToMultiple[nRow][0];
		   v[nRow][0] += tridiagonalMatrix_[2][nRow]*v_ColVectorToMultiple[nRow-1][0];
		}
	}
}


	//! A*x=r, A matrix, x,b column vector
//! A,b known, and the solution is "x"
void TridiagonalMatrix :: solve_linear_equation(const Matrix& r, Matrix& u) const
{
	if( sizeMatrix_ != r.rows() || !(r.if_col_vector()) )
	{
		throw ("Error in Matrix :: inverse(): not square matrix.");
	}

	//Matrix u(sizeMatrix_,1);

	int n = sizeMatrix_;
	double bet;
	std::vector<double> gam(n);
	//if(b[0] == 0.0)
	if(tridiagonalMatrix_[1][0] == 0.0)
		throw("Error 1 in function Matrix::solve_linear_equation(...): 0 in the tridiagonal");
	//u[0] = r[0] / (bet=b[0]);
	u[0][0] = r[0][0] / (bet=tridiagonalMatrix_[1][0]);

	for(int j=1; j<n; ++j)
	{
	    //gam[j] = c[j-1]/bet;
		gam[j] = tridiagonalMatrix_[0][j-1]/bet;
		//bet    = b[j] - a[j]*gam[j];
		bet    = tridiagonalMatrix_[1][j] - tridiagonalMatrix_[2][j]*gam[j];
		if(bet == 0.0)
		{
			throw("Error in function Matrix::solve_linear_equation(...): 0 in the tridiagonal");
		}
		//u[j] = (r[j] - a[j]*u[j-1]) /bet;
		u[j][0] = (r[j][0] - tridiagonalMatrix_[2][j]*u[j-1][0]) /bet;
	}
	for(int j=(n-2); j>=0; --j)
	{
	    u[j][0] -= gam[j+1]*u[j+1][0]; 
	}
	//return u;
}

void TridiagonalMatrix :: down_solve_linear_equation(const Matrix& r, Matrix& u) const  // TridiagonalMatrix is a lower tridiagonal matrix
{
	if( sizeMatrix_ != r.rows() || !(r.if_col_vector()) )
	{
		throw ("Error in Matrix :: inverse(): not square matrix.");
	}
    
	int n = sizeMatrix_;
	for(int i=0; i<n; ++i)
	{
	    if(tridiagonalMatrix_[1][i] == 0.0)
			throw ("Error in function Matrix::down_solve_linear_equation(...): 0 in the tridiagonal");
	}
    
	//Matrix u(r);
	u = r;
	for(int i=0; i<n-1; ++i)
	{
	    u[i][0]   /= tridiagonalMatrix_[1][i];
		u[i+1][0] = u[i+1][0] - tridiagonalMatrix_[2][i+1]*u[i][0];
	}
	u[n-1][0]   /= tridiagonalMatrix_[1][n-1];

	//return u;
}


void TridiagonalMatrix :: up_solve_linear_equation(const Matrix& r, Matrix& u) const  // TridiagonalMatrix is a upper tridiagonal matrix
{
	if( sizeMatrix_ != r.rows() || !(r.if_col_vector()) )
	{
		throw ("Error in Matrix :: inverse(): not square matrix.");
	}
    
	int n = sizeMatrix_;
	for(int i=0; i<n; ++i)
	{
	    if(tridiagonalMatrix_[1][i] == 0.0)
			throw ("Error in function Matrix::down_solve_linear_equation(...): 0 in the tridiagonal");
	}
    
	//Matrix u(r);
	u = r;
	for(int i=n-1; i>0; --i)
	{
	    u[i][0]   /= tridiagonalMatrix_[1][i];
		u[i-1][0] = u[i-1][0] - tridiagonalMatrix_[0][i-1]*u[i][0];
	}
	u[0][0]   /= tridiagonalMatrix_[1][0];

	//return u;
}


void TridiagonalMatrix :: set_Matrix(Matrix& m)
{
	if(m.rows() != m.cols() || m.rows() != this->get_sizeMatrix())
		std::cout << "matrix size dismatch in TridiagonalMatrix :: set_Matrix(Matrix& m)" << std::endl;

    for(int i=0; i<sizeMatrix_; i++)
	{
	    for(int j=0; j<sizeMatrix_; j++)
		{
			if(i-j == -1)
				m[i][j] = tridiagonalMatrix_[0][i];
			else if(i-j == 0)
				m[i][j] = tridiagonalMatrix_[1][i];
			else if(i-j == 1)
				m[i][j] = tridiagonalMatrix_[2][i];
			else
				m[i][j] = 0;
		}
	}
}


std::vector<double> TridiagonalMatrix::eigenvalues() 
{
	Matrix m(get_sizeMatrix(),get_sizeMatrix(),0.0);
    set_Matrix(m);
	return m.eigenvalues();
}


void TridiagonalMatrix:: print(std::string name) const 
{
	std::cout << "TridiagonalMatrix " << name << "  of size (" << sizeMatrix_ << " , " << sizeMatrix_ << ") :" << std::endl;
    for(int i=0; i<sizeMatrix_; i++)
	{
	    for(int j=0; j<sizeMatrix_; j++)
		{
			if(i-j == -1)
				std::cout << tridiagonalMatrix_[0][i]<< "	";
			else if(i-j == 0)
				std::cout << tridiagonalMatrix_[1][i] << "	";
			else if(i-j == 1)
				std::cout << tridiagonalMatrix_[2][i] << "	";
			else
				std::cout << 0 << "	";
		}
		std::cout << std::endl;
	}
}

bool TridiagonalMatrix::if_id_TriMatrix() const
{
	if(special_type_ != "id") 
	{
		return false;
	}
	else
	{
		for(int i=0; i<sizeMatrix_; i++)
		{
			if(tridiagonalMatrix_[0][i] !=0 || tridiagonalMatrix_[1][i] !=1 || tridiagonalMatrix_[2][i] !=0)
				return false;
		}
		return true;
	}
};
