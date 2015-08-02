#pragma once

#include <PDE/PDE_2D_Pricer.h>
#include <Cheyette/Pricer/PDE_Cheyette_Model.h>
#include <boost/shared_ptr.hpp>

class PDE_CheyetteDD_Pricer: public PDE_2D_Pricer
{
public:
	//! constructor
	PDE_CheyetteDD_Pricer(	Discretization_CONSTPTR        discretization,
							PDE_CheyetteDD_Model_CONSTPTR  pde_cheyette_model,
							Scheme_ADI_PTR                 scheme, // not const
							PDE_ADI_Solver_CONSTPTR        solver)
				:   PDE_2D_Pricer(  discretization,
					                pde_cheyette_model,
									scheme,
									solver)
	{}

	//! taking parameter for discretization
	PDE_CheyetteDD_Pricer(	//! discretization param: uniform discretization
		            const Range& range_t,
					const Range& range_x,
					const Range& range_y,
					size_t sizeDiscretization_t,
					size_t sizeDiscretization_x,
					size_t sizeDiscretization_y,

					PDE_CheyetteDD_Model_CONSTPTR  pde_cheyette_model,
					Scheme_ADI_PTR                 scheme,  // not const
					PDE_ADI_Solver_CONSTPTR        solver)
					:   PDE_2D_Pricer(Discretization_CONSTPTR(new Discretization(range_t, sizeDiscretization_t,
																				 range_x, sizeDiscretization_x, 
																			     range_y, sizeDiscretization_y)),  // uniform discretization
									  pde_cheyette_model,
									  scheme,
									  solver)
	{}

	virtual ~PDE_CheyetteDD_Pricer(){};
};

typedef boost::shared_ptr<PDE_CheyetteDD_Pricer> PDE_CheyetteDD_Pricer_PTR;
typedef boost::shared_ptr<const PDE_CheyetteDD_Pricer> PDE_CheyetteDD_Pricer_CONSTPTR;