#include <iostream>
#include <PDE/Solver/Linear_Equation_Solver.h>
#include <PDE/Matrix/Matrix.h>

using namespace std;
const bool DEBUG_optimize_linear_equation_solver = true;  // huge difference :) 

//! A*U^{n+1} + G = B*U^{n} + F
//! A or B = NULL means they are Id-matrix. 
void Linear_Equation_Solver::resolve(
	 const TridiagonalMatrix& A_tri,  
	 const TridiagonalMatrix& B_tri,
	 const Matrix& F,
	 const Matrix& G,
	 Matrix& U,			//input
	 Matrix& Uplus,
	 Matrix& U_temp) 	//output
{
    //Matrix U_temp(U.rows(), U.cols()); // copy not efficient ... !!! 
	B_tri.multipleColVector(U,U_temp);  
	U_temp += F;    
	U_temp -= G; 
	//if(A_tri.if_id_TriMatrix() == true)
	//{
	//	//A_tri.print();
	//	//getchar();
	//    Uplus = U_temp;	
	//}
	//else
	{
		A_tri.solve_linear_equation(U_temp,Uplus);	
	}
	
}


