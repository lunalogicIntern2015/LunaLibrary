#pragma once


//#include <PDE/HestonModelYuan.h>
#include <PDE/BoundaryCondition/BoundaryCondition_D2.h>
#include <PDE/Discretization/Discretization.h>
//#include <PDE/SchemasBuilder.h>
#include <PDE/Solver/Linear_Equation_Solver.h>
//#include "U_boundaryVal_Struct.h"
#include <utility>
#include <vector>
#include <PDE/Matrix/Matrix.h>

#include <PDE/Scheme/Scheme_ADI.h>
#include <boost/shared_ptr.hpp>


//! Solve equation: A_tri*U^{n+1} + G^{n+1} = B_tri*U^n + F^n
class PDE_ADI_Solver
{	
    //! keep reference because of polymorphism 
	Scheme_ADI&       scheme_;  //! A,B,H

    //! To stock the solution
	Matrix U;           //! Matrix D2 size: (matrix_size.first+2, matrix_size.second+2)

   

public:

	PDE_ADI_Solver(Scheme_ADI& scheme);

    virtual ~PDE_ADI_Solver();
	void Initialize_U();             // U = BCI
	void solve_PDE();  
    
	//! copy-coller
	Matrix get_result(){return U;}
};

typedef boost::shared_ptr<PDE_ADI_Solver> PDE_ADI_Solver_PTR;
typedef boost::shared_ptr<const PDE_ADI_Solver> PDE_ADI_Solver_CONSTPTR;