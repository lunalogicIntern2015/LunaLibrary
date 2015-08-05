//#include <JBLMM/Longstaff_Schwartz/ExplanatoryVariable.h>
//#include <LMM/pricer/McLmmVanillaSwapPricer.h>
//
//
//
//double EV_Libor::evaluate(	McLmm_PTR mclmm,
//							BermudanSwaption_CONSTPTR bermudanSwaption_PTR,
//							LMM::Index liborIndex)const
//{	
//	VanillaSwap_CONSTPTR vanillaswap= bermudanSwaption_PTR->getVanillaSwap_PTR();
//	return mclmm->get_liborMatrix()(liborIndex, liborIndex);
//}
//
//EV_vanillaswaprate::EV_vanillaswaprate(McLmmVanillaSwapPricer_PTR mcLmmVanillaSwapPricer_PTR)
//	:
//	mcLmmVanillaSwapPricer_PTR_(mcLmmVanillaSwapPricer_PTR)
//{
//}
//
//double EV_vanillaswaprate::evaluate(McLmm_PTR mclmm,
//									BermudanSwaption_CONSTPTR bermudanSwaption_PTR,
//									LMM::Index liborIndex)const
//{
//	//To check the indice liborMatrix, numeraire, bermudanSwaption_PTR, exerciceIndex
//
//	LMM::Index endIndex = bermudanSwaption_PTR->getVanillaSwap_PTR()->get_indexEnd();
//	VanillaSwap_CONSTPTR subVanillaSwap_PTR = bermudanSwaption_PTR->getSubVanillaSwap(liborIndex, endIndex);
//
//	mcLmmVanillaSwapPricer_PTR_->resetMcLmm(mclmm);
//
//	double vanillaSwapRate = mcLmmVanillaSwapPricer_PTR_->swapRate(	liborIndex,
//																	*subVanillaSwap_PTR.get(),
//																	mclmm->get_numeraire(), 
//																	mclmm->get_liborMatrix());
//
//	return vanillaSwapRate;
//}