#pragma once
#include <LMM/Pricer/McLmmPricer/McLmmPricer.h>
#include <LMM/Pricer/Longstaff_Schwartz/LS_BackwardAlgo.h>
#include <LMM/Pricer/Longstaff_Schwartz/LS_ForwardAlgo.h>

class McLmmCallableGenericSwapPricer
{
	McLmm_LS_PTR pMcLmm_LS_;
	LS_BackwardAlgo ls_BackwardAlgo_;
	LS_ForwardAlgo ls_ForwardAlgo_;

public:
	//constructor and destructor
	McLmmCallableGenericSwapPricer(	McLmm_LS_PTR pMcLmm_LS, 
									const LS_BackwardAlgo& ls_BackwardAlgo,
									const LS_ForwardAlgo& ls_ForwardAlgo);
	virtual ~McLmmCallableGenericSwapPricer(){}

	double pricing(	CallableInstrument_CONSTPTR callableInstrument, 
					size_t nbSimulation,
					std::vector<std::vector<double>>& basis_value_on_allPath_buffer);
};

