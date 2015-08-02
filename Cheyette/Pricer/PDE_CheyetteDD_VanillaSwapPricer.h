#pragma once

#include <Cheyette/Pricer/PDE_Cheyette_Pricer.h>	
#include <LMM/instrument/VanillaSwap.h>  
#include <PDE/BoundaryCondition_D2.h>
#include <boost/shared_ptr.hpp>



class PDE_CheyetteDD_VanillaSwapPricer : public PDE_CheyetteDD_Pricer
{
public:
	PDE_CheyetteDD_VanillaSwapPricer( Discretization_CONSTPTR        discretization,
									  PDE_CheyetteDD_Model_CONSTPTR  pde_cheyette_model,
									  Scheme_ADI_PTR                 scheme,
									  PDE_ADI_Solver_CONSTPTR        solver)
	:   PDE_CheyetteDD_Pricer(discretization,
							  pde_cheyette_model,
							  scheme,
							  solver)
	{}

	virtual ~PDE_CheyetteDD_VanillaSwapPricer(){}

	double swapNPV (VanillaSwap_PTR vanillaSwap, size_t valuationIndex) const;

};


typedef boost::shared_ptr<PDE_CheyetteDD_VanillaSwapPricer> PDE_CheyetteDD_VanillaSwapPricer_PTR;
typedef boost::shared_ptr<const PDE_CheyetteDD_VanillaSwapPricer> PDE_CheyetteDD_VanillaSwapPricer_CONSTPTR;