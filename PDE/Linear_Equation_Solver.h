#pragma once

#include "useful_function.h"
#include <PDE/Matrix/TridiagonalMatrix.h>
//#include <ql/quantlib.hpp> // to use the Matrix class
#include <PDE/Matrix/Matrix.h>
#include <boost/shared_ptr.hpp>




//! Linear Equation to solve//: 
//#include "Matrix.h"
//  A_tri * U^{n+1} + G^{n+1} = B_tri U^{n} + F^n
// 
//! Attensioin of the following different meaning of pointer:
//! A,B = NULL --> A,B = Id

class Linear_Equation_Solver
{
public:
	virtual ~Linear_Equation_Solver(){};

		//! U = U^n, Uplus=U^{n+1} 
	static void resolve(
		 const TridiagonalMatrix& A_tri, // size = (n,3) 
		 const TridiagonalMatrix& B_tri, // size = (n,3) 
		 const Matrix& F,				   // size = (n,1)
		 const Matrix& G,				   // size = (n,1)
		 Matrix& U,				   // size = (n,1)  input
		 Matrix& Uplus,
		 Matrix& U_temp);	       // size = (n,1)  output
};

typedef boost::shared_ptr<Linear_Equation_Solver> Linear_Equation_Solver_PTR;
typedef boost::shared_ptr<const Linear_Equation_Solver> Linear_Equation_Solver_CONSTPTR;