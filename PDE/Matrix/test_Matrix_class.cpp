#include "test_Matrix_class.h"
#include "TridiagonalMatrix.h"
//#include <boost/lambda/lambda.hpp>
#include <iostream>
#include <iterator>
#include <algorithm>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>

using namespace std;

double** create_matrix_memory(int rows, int cols)
{
    double ** m = new double*[rows];
	for(int i=0; i<rows; ++i)
	{
	    m[i] = new double[cols];
		for(int j=0; j<cols; ++j)
		{
		    m[i][j] = 0.0;
		}
	}
	return m;
}

void delete_matrix_memory(double ** v, int rows)
{
	for(int i=0; i<rows; ++i)
	{
	    delete[]  v[i];
	}
	delete[] v;
}

void test_constructor()
{
	cout << " \n\n ---- test_constructor ---- " << endl;
	int rows = 2;
	int cols = 3;

	// vector<vector<double>> v(rows, vector<double>(cols, 0.0));
	double ** v = create_matrix_memory(rows, cols);
	v[0][0] = 1;
	v[0][1] = 2;
	v[0][2] = 3;

	v[1][0] = 4;
	v[1][1] = 5;
	v[1][2] = 6;

    Matrix m(v, rows, cols);
	m.print();

	delete_matrix_memory(v,rows);
}

void test_access()
{
	cout << " \n\n  ---- test_access ---- " << endl;
	int rows = 2;
	int cols = 3;
    Matrix m(rows, cols);

	m(0,0) = 1.1;
	m(1,1) = 3.3;
	m.print();
}

void test_add()
{
	cout << " \n\n  ---- test_add ---- " << endl;
	int rows = 2;
	int cols = 3;

	//vector<vector<double>> v(rows, vector<double>(cols, 0.0));
	double ** v = create_matrix_memory(rows, cols);
	v[0][0] = 1;
	v[0][1] = 2;
	v[0][2] = 3;

	v[1][0] = 4;
	v[1][1] = 5;
	v[1][2] = 6;

    Matrix m(v, rows, cols);
	m.print();
	delete_matrix_memory(v,rows);

	//vector<vector<double>> v2(rows, vector<double>(cols, 0.0));
	double ** v2 = create_matrix_memory(rows, cols);
	v2[0][0] = 6;
	v2[0][1] = 5;
	v2[0][2] = 4;

	v2[1][0] = 3;
	v2[1][1] = 2;
	v2[1][2] = 1;

	Matrix m2(v2, rows, cols);
	m2.print();
	delete_matrix_memory(v2,rows);

	Matrix m3 = m + m2;
	m3.print();
}


void test_moins()
{
	cout << " \n\n  ---- test_moins ---- " << endl;
	int rows = 2;
	int cols = 3;

	//vector<vector<double>> v(rows, vector<double>(cols, 0.0));
	double ** v = create_matrix_memory(rows, cols);
	v[0][0] = 1;
	v[0][1] = 2;
	v[0][2] = 3;

	v[1][0] = 4;
	v[1][1] = 5;
	v[1][2] = 6;

    Matrix m(v, rows, cols);
	m.print();
	delete_matrix_memory(v,rows);

	//vector<vector<double>> v2(rows, vector<double>(cols, 0.0));
    double ** v2 = create_matrix_memory(rows, cols);
	v2[0][0] = 6;
	v2[0][1] = 5;
	v2[0][2] = 4;

	v2[1][0] = 3;
	v2[1][1] = 2;
	v2[1][2] = 1;

	Matrix m2(v2, rows, cols);
	m2.print();
	delete_matrix_memory(v2,rows);

	Matrix m3 = m - m2;
	m3.print();
}


void test_transpose()
{
	cout << " \n\n  ---- test_transpose ---- " << endl;
	int rows = 2;
	int cols = 3;

	//vector<vector<double>> v(rows, vector<double>(cols, 0.0));
	double ** v = create_matrix_memory(rows, cols);
	v[0][0] = 1;
	v[0][1] = 2;
	v[0][2] = 3;

	v[1][0] = 4;
	v[1][1] = 5;
	v[1][2] = 6;

    Matrix m(v, rows, cols);
	m.print();
	delete_matrix_memory(v,rows);

	Matrix m2(~m);
	m2.print();

}


void test_mulply()
{
	cout << " \n\n  ---- test_mulply ---- " << endl;
	int rows = 2;
	int cols = 3;

	//vector<vector<double>> v(rows, vector<double>(cols, 0.0));
	double ** v = create_matrix_memory(rows, cols);
	v[0][0] = 1;
	v[0][1] = 2;
	v[0][2] = 3;

	v[1][0] = 4;
	v[1][1] = 5;
	v[1][2] = 6;

    Matrix m(v, rows, cols);
	m.print();
	delete_matrix_memory(v,rows);

	//vector<vector<double>> v2(cols, vector<double>(rows, 0.0));
	double ** v2 = create_matrix_memory(cols, rows);
	v2[0][0] = 6;
	v2[0][1] = 5;

	v2[1][0] = 4;
	v2[1][0] = 3;

	v2[2][0] = 2;
	v2[2][1] = 1;

	Matrix m2(v2, cols, rows);
	m2.print();
	delete_matrix_memory(v2,cols);

	Matrix m3 = m * m2;
	m3.print();
}


void test_determinant()
{
	cout << " \n\n  ---- test_determinant ---- " << endl;
	//vector<vector<double>> v(3, vector<double>(3, 0.0));
	int rows = 3;
	int cols = 3;

	double ** v = create_matrix_memory(rows, cols);
	v[0][0] = 6;
	v[0][1] = 5;
	v[0][2] = 4;

	v[1][0] = 3;
	v[1][1] = 2;
	v[1][2] = 1;

	v[2][0] = 2;
	v[2][1] = 1;
	v[2][2] = 2;

	Matrix m(v, rows, cols);
	delete_matrix_memory(v,rows);
	m.print();
	cout << "determinant  = " << m.determinant()  << endl;
}

//void test_matrix_multiply_vector()
//{
//    cout << " \n\n  ---- test_matrix_multiply_vector ---- " << endl;
//	int rows = 2;
//	int cols = 3;
//
//	//vector<vector<double>> v(rows, vector<double>(cols, 0.0));
//	double ** v = create_matrix_memory(rows, cols);
//	v[0][0] = 1;
//	v[0][1] = 2;
//	v[0][2] = 3;
//
//	v[1][0] = 4;
//	v[1][1] = 5;
//	v[1][2] = 6;
//
//    Matrix m(v, rows, cols);
//	m.print();
//	delete_matrix_memory(v,rows);
//
//	vector<double> v2(cols, 0.0);
//	v2[0] = 1;
//	v2[1] = 2;
//	v2[2] = 3;
//
//	Vector vec(v2, true); // column vector
//	vec.print();
//
//	Vector vec_result = m * vec;
//	vec_result.print();
//}

//void test_vector_multiply_matrix()
//{
//    cout << " \n\n  ---- test_vector_multiply_matrix ---- " << endl;
//	int rows = 2;
//	int cols = 3;
//
//	//vector<vector<double>> v(rows, vector<double>(cols, 0.0));
//	double ** v = create_matrix_memory(rows, cols);
//	v[0][0] = 1;
//	v[0][1] = 2;
//	v[0][2] = 3;
//
//	v[1][0] = 4;
//	v[1][1] = 5;
//	v[1][2] = 6;
//
//    Matrix m(v, rows, cols);
//	m.print();
//	delete_matrix_memory(v,rows);
//
//	vector<double> v2(rows, 0.0);
//	v2[0] = 1;
//	v2[1] = 2;
//
//	Vector vec(v2, false); // row vector
//	vec.print();
//
//	Vector vec_result = vec * m;
//	vec_result.print();
//}

void test_SVD()
{
	cout << " \n\n  ---- test_SVD ---- " << endl;
	int n = 3; 
	double** a = new double*[n];
	double** v = new double*[n];
	for(int i=0; i<n; ++i)
	{
		a[i] = new double[n]; 
		v[i] = new double[n]; 
	}
	double* w = new double[n];

	a[0][0] = 6.0;
	a[0][1] = 5.0;
	a[0][2] = 4.0;

	a[1][0] = 3.0;
	a[1][1] = 2.0;
	a[1][2] = 1.0;

	a[2][0] = 2.0;
	a[2][1] = 1.0;
	a[2][2] = 2.0;
	
    //svdcmp(double **a, int m, int n, double* w, double **v)
	svd(a, n, n, w, v);

	cout << " ---- a ----" << endl;
	for(int i=0; i<n; ++i)
	{
	    for(int j=0; j<n; ++j)
		{
		     cout << a[i][j] << " ";
		}
		cout << endl;
	}

	cout << " ---- w ----" << endl;
	for(int i=0; i<n; ++i)
	{
		cout << w[i] << " ";
	}
	cout << endl;

	cout << " ---- v ----" << endl;
	for(int i=0; i<n; ++i)
	{
	    for(int j=0; j<n; ++j)
		{
		     cout << v[i][j] << " ";
		}
		cout << endl;
	}
}

	
void test_inverse_matrix_by_SVD()
{
	cout << " \n\n  ---- test_inverse_matrix_by_SVD ---- " << endl;
	int n = 3; 
	double** a = new double*[n];
	double** a2 = new double*[n];
	double** m = new double*[n];
	for(int i=0; i<n; ++i)
	{
		a[i] = new double[n]; 
		a2[i] = new double[n]; 
		m[i] = new double[n]; 
	}

	a[0][0] = 6.0;
	a[0][1] = 5.0;
	a[0][2] = 4.0;

	a[1][0] = 3.0;
	a[1][1] = 2.0;
	a[1][2] = 1.0;

	a[2][0] = 2.0;
	a[2][1] = 1.0;
	a[2][2] = 2.0;

	for(int i=0; i<n; ++i)
	{
		for(int j=0; j<n; ++j)
		{
			a2[i][j] = a[i][j];
			m[i][j] = 0.0;
		}
	}

	svd_inverse_matrix(a,n);

	cout << " ---- inverse of a ----" << endl;
	for(int i=0; i<n; ++i)
	{
	    for(int j=0; j<n; ++j)
		{
		     cout << a[i][j] << " ";
		}
		cout << endl;
	}

	//! multiply a and inverse a
	for(int i=0; i<n; ++i)
	{
		for(int j=0; j<n; ++j)
		{
			for(int k=0; k<n; ++k)
			{
                m[i][j] += a[i][k] * a2[k][j];
			}
		}
	}

	cout << " ---- test the inverse ----" << endl;
	for(int i=0; i<n; ++i)
	{
	    for(int j=0; j<n; ++j)
		{
		     cout << m[i][j] << " ";
		}
		cout << endl;
	}

	
	for(int i=0; i<n; ++i)
	{
		delete[] a[i];
		delete[] a2[i];
	}

	delete[] a;
	delete[] a2;
	delete[] m;
}



void test_Squared_Matrix_inverse()
{
	cout << " \n\n  ---- test_Squared_Matrix_inverse ---- " << endl;
	int rows = 3;
	int cols = 3;
	// vector<vector<double>> v(3, vector<double>(3, 0.0));
	double ** v = create_matrix_memory(rows, cols);
	v[0][0] = 6;
	v[0][1] = 5;
	v[0][2] = 4;

	v[1][0] = 3;
	v[1][1] = 2;
	v[1][2] = 1;

	v[2][0] = 2;
	v[2][1] = 1;
	v[2][2] = 2;

	Matrix m(v, rows, cols);
	delete_matrix_memory(v,rows);
	Matrix m_inverse = m.inverse();
	m_inverse.print();
}

void read_vector(string& file_name, vector<double> & X, vector<double>& Y)
{
	ifstream file;
	string   line;
	string   sub_line;
	file.open(file_name.c_str());

	if(file.is_open())
	{
		while(!file.eof())
		{
			getline(file, line);
			istringstream iss(line); 
			int flag = 1;
			while(iss >> sub_line) 
			{
				double value = atof(sub_line.c_str());
				if(flag == 1)
				{
				    X.push_back(value); 
				}
				else
				{
					Y.push_back(value);
				}
				flag *= -1;
			}
		}    
	}

	file.close();
}


void test_regression()
{
	cout << " \n\n  ---- test_regression ---- " << endl;
	string file_name = "C:\\regression.txt";
	vector<double> X;
	vector<double> Y;
	read_vector(file_name,X,Y);

	int num_basis = 3; // 1,x,x^2
	Matrix m_X( (int)X.size(), (int)num_basis);
	for(int i=0; i<(int)m_X.rows(); ++i)
	{
		for(int j=0; j<m_X.cols(); ++j)
		{
            if(j ==0)
			{
				m_X(i,j) = 1.0;	
			}
			else
			{
			    m_X(i,j) = pow(X[i],(double)j);
			}
		}
	}

	Matrix m_Y((int)Y.size(),1);
	for(int i=0; i<m_Y.rows(); ++i)
	{
	    m_Y(i,0) = Y[i];
	}

	Matrix m_Z (((~m_X)*(m_X)).get_matrix(), m_X.cols(), m_X.cols()) ;
    Matrix m_W = m_Z.inverse();
    Matrix m_L = m_W*(~m_X)*m_Y;

	m_L.print();
}


void test_up_down_resolver()
{
	cout.precision(10);
	//! lower matrix - tridiagonal matrix
     TridiagonalMatrix A(3);

	 A.tridiagonalMatrix_[1][0] = 1.9852;
	 A.tridiagonalMatrix_[1][1] = 2.425;
	 A.tridiagonalMatrix_[1][2] = 3.98751;

	 A.tridiagonalMatrix_[2][0] = 0;
	 A.tridiagonalMatrix_[2][1] = 4.21651;
	 A.tridiagonalMatrix_[2][2] = 5.98441;

	 //A.print("A");

	 Matrix r(3,1,0.0);
	 r[0][0] = 3.14;
	 r[1][0] = 2.732;
	 r[2][0] = 1.414;
	 r.print("r");

	 //Matrix d = A.down_solve_linear_equation(r);
	 Matrix d(r.rows(), r.cols(), 0.0);
	 A.down_solve_linear_equation(r,d);
	 //d.print("d");

	 Matrix result(r);
	 A.multipleColVector(d,result);
	 result.print("result");

	 cout << "\n\n\n"<< endl;


	 //! lower matrix - tridiagonal matrix
     TridiagonalMatrix B(3);

	 B.tridiagonalMatrix_[1][0] = 1.0514;
	 B.tridiagonalMatrix_[1][1] = 2.524;
	 B.tridiagonalMatrix_[1][2] = 3.2222;

	 B.tridiagonalMatrix_[0][0] = 4.00124;
	 B.tridiagonalMatrix_[0][1] = 5.98575;
	 B.tridiagonalMatrix_[0][2] = 0;
	 //B.print("B");

	 r.print("r");
	 //Matrix u = B.up_solve_linear_equation(r);
	 Matrix u(r); 
	 B.up_solve_linear_equation(r,u);
	 //u.print("u");

	 B.multipleColVector(u,result);
	 result.print("result");
}

void test_Matrix_class()
{
	 ////! Passed test ...
	 //test_constructor();
	 //test_add();
	 //test_moins();
	 //test_transpose();
	 //test_mulply();
	 //test_determinant();
	 ////test_matrix_multiply_vector();
	 ////test_vector_multiply_matrix();
	 //test_SVD();
	 //test_inverse_matrix_by_SVD();
	 //test_Squared_Matrix_inverse();
	 //test_regression();
	
   	 //! Test to do ...
	 //test_access();

	 test_up_down_resolver();
}
