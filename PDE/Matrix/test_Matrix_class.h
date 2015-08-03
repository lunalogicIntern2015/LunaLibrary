#include "Matrix.h"
//#include "Squared_Matrix.h"
//#include "Vector.h"
#include <vector>
#include "svd.h"


double** create_matrix_memory(int rows,    int cols);
void     delete_matrix_memory(double ** v, int rows);

void test_constructor();
void test_access();
void test_add();
void test_moins();
void test_transpose();
void test_mulply();
void test_determinant();
//void test_matrix_multiply_vector();
void test_vector_multiply_matrix();
void test_SVD();
void test_Squared_Matrix_inverse();
void read_vector(string& file_name, vector<double> & X, vector<double>& Y);
void test_regression();

void test_up_down_resolver();


//! main test function
void  test_Matrix_class();