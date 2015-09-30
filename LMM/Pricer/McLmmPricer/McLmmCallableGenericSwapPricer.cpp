#include <LMM/Pricer/McLmmPricer/McLmmCallableGenericSwapPricer.h>


McLmmCallableGenericSwapPricer::McLmmCallableGenericSwapPricer(	McLmm_LS_PTR pMcLmm_LS, 
																const LS_BackwardAlgo& ls_BackwardAlgo,
																const LS_ForwardAlgo& ls_ForwardAlgo)
																:
																pMcLmm_LS_(pMcLmm_LS),
																ls_BackwardAlgo_(ls_BackwardAlgo),
																ls_ForwardAlgo_(ls_ForwardAlgo)
{
}


double McLmmCallableGenericSwapPricer::pricing(	CallableInstrument_CONSTPTR callableInstrument, 
												size_t nbSimulation,
												std::vector<std::vector<double>>& basis_value_on_allPath_buffer)
{
	//forward simulatoin for backward regression
	pMcLmm_LS_->simulateLMM(nbSimulation);
	const std::vector<McLmm_LS::LMMSimulationResult>&  lmmSimualtionResults_backward = pMcLmm_LS_->lmmSimualtionResults_;

	//Backward Regression
	ls_BackwardAlgo_.do_BackwardAlgo(	callableInstrument,
										lmmSimualtionResults_backward,
										basis_value_on_allPath_buffer);

	//forward simulation for pricing forward
	pMcLmm_LS_->simulateLMM(nbSimulation);
	const std::vector<McLmm_LS::LMMSimulationResult>&  lmmSimualtionResults_forward = pMcLmm_LS_->lmmSimualtionResults_;

	//pricing forward
	std::pair<double, double> result = ls_ForwardAlgo_.do_ForwardAlgo(	callableInstrument, 
																		lmmSimualtionResults_forward,
																		ls_BackwardAlgo_.getRegressions());

	double price = 	result.first;

	return price;
	
}

