#include <PDE/useful_function.h>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>

//#include "TridiagonalMatrix.h"

//#include <ql/quantlib.hpp> // to use the Matrix class

#include <PDE/Matrix/Matrix.h>


double inverse_sinh(double x)
{
    return log(x+sqrt(x*x+1));
}

bool if_positive(double x)
{
    if(x>0)
		return true;
	else
		return false;
}


//! vector is supposed to be implemented as (n,1) matrix
//! --------------------------------------
//! 
//!    vector to Matrix operation 
//! 
//! --------------------------------------
void add_vectorToMatrixColumn(Matrix& matrix, Matrix& vectorCol, int index_col)
{
    if(matrix.rows() != vectorCol.rows() || vectorCol.cols() !=1)
		std::cout << "Error in copy_vectorToMatrixColumn(...), matrix size dismatch! "<< std::endl;
	else
	{
	    for(int i=0; i<matrix.rows(); ++i)
		{
		    matrix[i][index_col] += vectorCol[i][0];
		}
	}   
}
void copy_vectorToMatrixColumn(Matrix& matrix, Matrix& vectorCol, int index_col)
{
    if(matrix.rows() != vectorCol.rows() || vectorCol.cols() !=1)
		std::cout << "Error in copy_vectorToMatrixColumn(...), matrix size dismatch! "<< std::endl;
	else
	{
	    for(int i=0; i<matrix.rows(); ++i)
		{
		    matrix[i][index_col] = vectorCol[i][0];
		}
	}   
}
void add_vectorToMatrixColumn_withoutFirstAndLastElment(Matrix& matrix, Matrix& vectorCol, int index_col)
{
    if(matrix.rows() != vectorCol.rows()+2 || vectorCol.cols() !=1)
		std::cout << "Error incopy_vectorToMatrixColumn_withoutFirstAndLastElment(...), matrix size dismatch! " <<  std::endl;
	else
	{
	    for(int i=0; i<vectorCol.rows(); ++i)
		{
		    matrix[i+1][index_col] += vectorCol[i][0];
		}
	}   
}


void copy_vectorToMatrixColumn_withoutFirstAndLastElment(Matrix& matrix, Matrix& vectorCol, int index_col)
{
    if(matrix.rows() != vectorCol.rows()+2 || vectorCol.cols() !=1)
		std::cout << "Error incopy_vectorToMatrixColumn_withoutFirstAndLastElment(...), matrix size dismatch! " <<  std::endl;
	else
	{
	    for(int i=0; i<vectorCol.rows(); ++i)
		{
		    matrix[i+1][index_col] = vectorCol[i][0];
		}
	}   
}


void add_vectorToMatrixRow(Matrix& matrix, Matrix& vectorRow, int index_row)
{
    if(matrix.cols() != vectorRow.rows() || vectorRow.cols() !=1)
		std::cout << "Error matrix size dismatch! "<< std::endl;
	else
	{
	    for(int i=0; i<matrix.cols(); ++i)
		{
		     matrix[index_row][i] += vectorRow[i][0];
		}
	}
}
void copy_vectorToMatrixRow(Matrix& matrix, Matrix& vectorRow, int index_row)
{
    if(matrix.cols() != vectorRow.rows() || vectorRow.cols() !=1)
		std::cout << "Error matrix size dismatch! "<< std::endl;
	else
	{
	    for(int i=0; i<matrix.cols(); ++i)
		{
		     matrix[index_row][i] = vectorRow[i][0];
		}
	}
}
void add_vectorToMatrixRow_withoutFirstAndLastElment(Matrix& matrix, Matrix& vectorRow, int index_row)
{
    if(matrix.cols() != vectorRow.rows()+2 || vectorRow.cols() !=1)
		std::cout << "Error in copy_vectorToMatrixRow_withoutFirstAndLastElment(...), matrix size dismatch! "<< std::endl;
	else
	{
	    for(int i=0; i<vectorRow.rows(); ++i)
		{
		     matrix[index_row][i+1] += vectorRow[i][0];
		}
	}
}

void copy_vectorToMatrixRow_withoutFirstAndLastElment(Matrix& matrix, Matrix& vectorRow, int index_row)
{
    if(matrix.cols() != vectorRow.rows()+2 || vectorRow.cols() !=1)
		std::cout << "Error in copy_vectorToMatrixRow_withoutFirstAndLastElment(...), matrix size dismatch! "<< std::endl;
	else
	{
	    for(int i=0; i<vectorRow.rows(); ++i)
		{
		     matrix[index_row][i+1] = vectorRow[i][0];
		}
	}
}

//! --------------------------------------
//! 
//!    Matrix to vector operation 
//! 
//! --------------------------------------

void copy_matrixColumnToVector(Matrix& matrix, Matrix& vectorCol, int index_col )
{
    if(matrix.rows() != vectorCol.rows() || vectorCol.cols() !=1)
		std::cout << "Error in copy_matrixColumnToVector(...), matrix size dismatch! " << std::endl;
	else
	{
	    for(int i=0; i<matrix.rows(); ++i)
		{
		    vectorCol[i][0] = matrix[i][index_col];
		}
	}
}

void copy_matrixColumnToVector_withoutFirstAndLastElment(Matrix& matrix, Matrix& vectorCol, int index_col )
{
    if(matrix.rows() != vectorCol.rows()+2 || vectorCol.cols() !=1)
		std::cout << "Error in copy_matrixColumnToVector_withoutFirstAndLastElment(...), matrix size dismatch! " << std::endl;
	else
	{
	    for(int i=0; i<vectorCol.rows(); ++i)
		{
		    vectorCol[i][0] = matrix[i+1][index_col];
		}
	}
}



void copy_matrixRowToVector(Matrix& matrix, Matrix& vectorRow, int index_row)
{
    if(matrix.cols() != vectorRow.rows() || vectorRow.cols() !=1)
		std::cout << "Error in copy_matrixRowToVector(...), matrix size dismatch! "<< std::endl;
	else
	{
	    for(int i=0; i<matrix.cols(); ++i)
		{
		    vectorRow[i][0] = matrix[index_row][i];
		}
	}
}
void copy_matrixRowToVector_withoutFirstAndLastElment(Matrix& matrix, Matrix& vectorRow, int index_row)
{
    if(matrix.cols() != vectorRow.rows()+2 || vectorRow.cols() !=1)
		std::cout << "Error in copy_matrixRowToVector_withoutFirstAndLastElment(...), matrix size dismatch! "<< std::endl;
	else
	{
	    for(int i=0; i<vectorRow.rows(); ++i)
		{
		    vectorRow[i][0] = matrix[index_row][i+1];
		}
	}
}



//! --------------------------------------
//! 
//!    Matrix to Matrix operation 
//! 
//! --------------------------------------

void add_smallMatrix_to_bigMatrix_withOutBoundary(Matrix& m_small, Matrix& m_big)
{
    if(m_small.rows() +2 == m_big.rows() && m_small.cols() +2 == m_big.cols())
	{ 
		for(int i=0; i<m_small.rows(); ++i)
		{
			for(int j=0; j<m_small.cols(); ++j)
			{
			    m_big[i+1][j+1] += m_small[i][j];
			}
		}
	}
	else
	{
		throw ("Error in function add_matrix_withOutBoundary(...), matrix sizes don't match");
	}
}

void copy_smallMatrix_to_bigMatrix_withOutBoundary(Matrix& m_small, Matrix& m_big)
{
    if(m_small.rows() +2 == m_big.rows() && m_small.cols() +2 == m_big.cols())
	{ 
		for(int i=0; i<m_small.rows(); ++i)
		{
			for(int j=0; j<m_small.cols(); ++j)
			{
			    m_big[i+1][j+1] = m_small[i][j];
			}
		}
	}
	else
	{
		throw ("Error in function copy_matrix_withOutBoundary(...), matrix sizes don't match");
	}
}




void init_to_zero_colMatrix(Matrix& m)
{
    if(m.cols() != 1)
		std::cout << "Error in function init_to_zero_vector(...), matrix size dismatch "<< std::endl;

	for(int i=0; i<m.rows(); ++i)
		m[i][0] = 0;
}


double max(double x, double y)
{
   if(x>y)
     return  x;
   else 
     return y;
}

double digital(double x, double y)
{
   if(x>y)
	   return 1;
   else 
	   return 0;
}



void print(Matrix& m)
{
	std::cout << "print matrix m ["<< m.rows() << "," << m.cols() <<  "] " << std::endl;
	for(int i=0; i<m.rows(); i++)
	{
	    for(int j=0; j<m.cols(); j++)
		{
			std::cout << m[i][j]<<"		";
		}
		std::cout << std::endl;
	}
}

void print(Matrix& m, std::string s)
{
	std::ofstream myfile(s.c_str());
    myfile << "print matrix m ["<< m.rows() << "," << m.cols() <<  "] \n" << std::endl;
	for(int i=0; i<m.rows(); i++)
	{
	    for(int j=0; j<m.cols(); j++)
		{
			myfile << m[i][j]<<"		";
		}
		myfile << std::endl;
	}
	myfile.close();
}

void print(TridiagonalMatrix& m, std::string s)
{
	    std::ofstream myfile(s.c_str());
		myfile << "TridiagonalMatrix of size (" << m.sizeMatrix_ << " , " << m.sizeMatrix_ << ") :" << std::endl;
	    for(int i=0; i<m.sizeMatrix_; i++)
		{
		    for(int j=0; j<m.sizeMatrix_; j++)
			{
				if(i-j == -1)
					myfile << m.tridiagonalMatrix_[0][i]<< "	";
				else if(i-j == 0)
					myfile << m.tridiagonalMatrix_[1][i] << "	";
				else if(i-j == 1)
					myfile << m.tridiagonalMatrix_[2][i] << "	";
				else
					myfile << 0 << "	";
			}
			myfile << std::endl;
		}
	   myfile.close();
}


void print(const std::vector<double>& m, const std::string s)
{
	    std::ofstream myfile(s.c_str());
		myfile << "vector<double> of size " << m.size()  << std::endl;
	    for(unsigned int i=0; i<m.size(); i++)
		{
			myfile << exp(m[i]) << std::endl;
		}
	   myfile.close();
}

void print_MatrixToFile(std::string& file_name, Matrix& m)
{
	std::ofstream my_file(file_name.c_str());
    for( int i=0; i<m.rows(); ++i)
	{
		for( int j=0; j<m.cols(); ++j)
		{
		    my_file << m[i][j] <<",";
		}
		my_file << std::endl;
	}
} 

void print_result(std::string& file_name, std::vector<double>& x_grid, std::vector<double>& y_grid, Matrix& m)
{
	std::ofstream my_file(file_name.c_str());

	//! y_grid
	my_file << ",";

	my_file.precision(10);
	for( int j=0; j<m.cols(); ++j)
	{
	    my_file << y_grid[j] <<",";
	}
	my_file << std::endl;

    for( int i=0; i<m.rows(); ++i)
	{
		my_file << x_grid[i]<<",";
		for( int j=0; j<m.cols(); ++j)
		{
		    my_file << m[i][j] <<",";
		}
		my_file << std::endl;
	}
} 

// Inverse Cumulative Normale
double Inverse_CN(double u){
	static double a[4] = {2.50662823884, -18.61500062529, 41.39119773534, -25.44106049637};
	static double b[4] = {-8.47351093090, 23.08336743743, -21.06224101826, 3.13082909833};
	static double c[9] = {0.3374754822726147, 0.9761690190917186, 0.1607979714918209, 0.0276438810333863,
	                      0.0038405729373609, 0.0003951896511919, 0.0000321767881768, 0.0000002888167364,
						  0.0000003960315187};

	double x = u - 0.5;
	double r;

	if(fabs(x)<0.42){
	    double  y = x*x;
		r = x*(((a[3]*y+a[2])*y+a[1])*y+a[0]) / ((((b[3]*y+b[2])*y+b[1])*y+b[0])*y+1.0);
	}else{
	  r=u;
	  if(x>0.0) r = 1.0-u;
	  r = log(-log(r));
	  r = c[0]+r*(c[1]+r*(c[2]+r*(c[3]+r*(c[4]+r*(c[5]+r*(c[6]+r*(c[7]+r*c[8])))))));
	  if(x<0.0) r = -r;
	}
	return r;
}


double NormalDensity(double x){
  return pi_2*exp(-x*x/2);
}

//Cumulative Normale
double CN(double x){
	static double a[5] = {0.319381530, -0.356563782, 1.781477937, -1.821255978, 1.330274429};
    double result;

	if(x<-7.0) result = NormalDensity(x) / sqrt(1.+x*x); 
	else{
	   if(x>7.0) result = 1.0 - CN(-x);
	   else{
	       double tmp = 1.0 / (1.0+0.2316419*fabs(x));
		   result = 1 - NormalDensity(x)* (tmp*(a[0]+tmp*(a[1]+tmp*(a[2]+tmp*(a[3]+tmp*a[4])))));
		   if(x<=0.0) result = 1.0 - result;
	   }
	}
	return result;
}

////! Copy from numerical recipe
////! Newton raphson method to find the root in interval [x1,x2], with tolerance xacc
//double rtnewt(void funcd(const double, double &, double &), const double x1, const double x2, const double xacc)
//{
//	const int JMAX = 20;
//	int j;
//	double df,dx,f,rtn;
//
//	rtn=0.5*(x1+x2);
//	for (j=0;j<JMAX;j++) {
//		funcd(rtn,f,df); // x, f(x), f'(x)
//		dx=f/df;
//		rtn -= dx;
//		if ((x1-rtn)*(rtn-x2) < 0.0)
//			std::cout << ("Jumped out of brackets in rtnewt") << std::endl;
//		if (fabs(dx) < xacc) return rtn;
//	}
//	std::cout << ("Maximum number of iterations exceeded in rtnewt") << std::endl;
//	return 0.0;
//}





//! Numerical recipe Chapiter 3
//! Polynomial interpolation
void polint(std::vector<double> &xa,
			std::vector<double> &ya,
			const double x,
			double &y,
			double &dy)
{
	int i,m,ns=0;
	double den,dif,dift,ho,hp,w;

	int n= (int)xa.size();
	std::vector<double> c(n),d(n);
	dif=fabs(x-xa[0]);
	for (i=0;i<n;i++) {
		if ((dift=fabs(x-xa[i])) < dif) {
			ns=i;
			dif=dift;
		}
		c[i]=ya[i];
		d[i]=ya[i];
	}
	y=ya[ns--];
	for (m=1;m<n;m++) {
		for (i=0;i<n-m;i++) {
			ho=xa[i]-x;
			hp=xa[i+m]-x;
			w=c[i+1]-d[i];
			if ((den=ho-hp) == 0.0) std::cout << "Error in routine polint" << std::endl;
			den=w/den;
			d[i]=hp*den;
			c[i]=ho*den;
		}
		y += (dy=(2*(ns+1) < (n-m) ? c[ns+1] : d[ns--]));
	}
}


//! Numerical recipe Chapiter 3
//! Two-dimensional polynomial interpolation
void polin2(std::vector<double>& x1a,
			std::vector<double>& x2a,
			std::vector<std::vector<double>>& ya,  // les points fix connues
			const double x1,const double x2, double &y, double &dy)   // les points a interpoler
{
	int j,k;

	int m= (int)x1a.size();
	int n= (int)x2a.size();
	std::vector<double> ymtmp(m);
	std::vector<double>  ya_t(n);
	for (j=0;j<m;j++) {
		for (k=0;k<n;k++) 
			ya_t[k]=ya[j][k];
		polint(x2a,ya_t,x2,ymtmp[j],dy);
	}
	polint(x1a,ymtmp,x1,y,dy);
}


