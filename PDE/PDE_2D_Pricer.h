#pragma once


#include <PDE/PDE_2D_Model.h>
#include <PDE/BoundaryCondition/BoundaryCondition_D2.h>
#include <PDE/Scheme/Scheme_ADI.h>
#include <PDE/Solver/PDE_ADI_Solver.h>
#include <boost/shared_ptr.hpp>

//! Cheyette model hiearchy is not well done, so for the moment only treat the cheyetteDD model

//! suppose underlying markovian states: (x_t,y_t).

class PDE_2D_Pricer
{
protected:
	//! why need smart ptr other than obj, because cannot use obj of abstract class, so & or ptr needed...
	Discretization_CONSTPTR		  discretization_;
	PDE_2D_Model_CONSTPTR		  pde_2d_model_;
	Scheme_ADI_PTR			      scheme_; // not const !
	PDE_ADI_Solver_CONSTPTR		  solver_;

	BoundaryCondition_D2_PTR      bc_;     // depend on product! 

public:
	//! constructor
	PDE_2D_Pricer(	 Discretization_CONSTPTR       discretization,
					 PDE_2D_Model_CONSTPTR         pde_2d_model,
					 Scheme_ADI_PTR                scheme,
					 PDE_ADI_Solver_CONSTPTR       solver) // this depends on product too: european or american ...
					:   discretization_(discretization),
					    pde_2d_model_(pde_2d_model),
					    scheme_(scheme),
					    solver_(solver)
						//bc_(NULL), // automatically initialized to NULL
						{}

	virtual ~PDE_2D_Pricer(){}

	//! getter
	Discretization_CONSTPTR	get_discretization() const {return discretization_;}
	PDE_2D_Model_CONSTPTR	get_pde_2d_model()   const {return pde_2d_model_;}
	Scheme_ADI_PTR			get_scheme()         const {return scheme_;} 
	PDE_ADI_Solver_CONSTPTR	get_pde_adi_solver() const {return solver_;}
	BoundaryCondition_D2_PTR get_bc()            const {return bc_;}
};

typedef boost::shared_ptr<PDE_2D_Pricer> PDE_2D_Pricer_PTR;
typedef boost::shared_ptr<const PDE_2D_Pricer> PDE_2D_Pricer_CONSTPTR;
