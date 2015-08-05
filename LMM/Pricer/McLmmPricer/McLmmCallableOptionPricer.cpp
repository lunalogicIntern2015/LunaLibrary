//#include <LMM/Helper/GenericPath.h>
////#include <LMM/Helper/Printer.h>
//#include <LMM/Pricer/McLmmPricer/McLmmCallableOptionPricer.h>
//
//
//McLmmCallableOptionPricer::McLmmCallableOptionPricer(McLmm_PTR mcLmm,
//													 size_t nb_1stFwdSimulation,
//													 size_t nb_2ndFwdSimulation)
//													 : mcLmm_(mcLmm)
//													 , nb_1stFwdSimulation_(nb_1stFwdSimulation)
//													 , liborMatrix_1stFwdSimulation_(nb_1stFwdSimulation)
//													 , numeraire_1stFwdSimulation_(nb_1stFwdSimulation)
//													 , nb_2ndFwdSimulation_(nb_2ndFwdSimulation){}
//
//void McLmmCallableOptionPricer::preparePricerForInstrument(CallableOption_PTR callableOption_ptr)
//{
//	buffer_callableDatesIndex_ = callableOption_ptr->get_callableDatesIndex_Ref();
//}
//
////! YY TODO: not efficient at all! 
//void McLmmCallableOptionPricer::fwdSimulation()      // 1st fwd simulation
//{
//	for(size_t i=0; i<nb_1stFwdSimulation_; ++i)
//	{
//		mcLmm_->simulateLMM();
//		liborMatrix_1stFwdSimulation_[i] = mcLmm_->get_liborMatrix(); // YY TODO: copy-coller too expensive!
//
//		mcLmm_->computeNumeraires();
//		numeraire_1stFwdSimulation_[i] = mcLmm_->get_numeraire();  // YY TODO: copy-coller too expensive!
//	}
//}
//
//
//double McLmmCallableOptionPricer::fwdPricing()
//{
//	double price = 0.0;
//
//	for(size_t i=0; i< nb_2ndFwdSimulation_; ++i)
//	{
//        mcLmm_->simulateLMM();
//		const matrix& liborMatrix = mcLmm_->get_liborMatrix();
//
//		mcLmm_->computeNumeraires();
//		const std::vector<double>& numeraire = mcLmm_->get_numeraire();
//
//		for(size_t itr_callableTimeIndex=0; itr_callableTimeIndex<buffer_callableDatesIndex_.size(); ++itr_callableTimeIndex)
//		{
//			double cv = continuousValue(itr_callableTimeIndex, liborMatrix, numeraire);
//			double iv = intrinsicValue (itr_callableTimeIndex, liborMatrix, numeraire);
//
//			if(iv > cv)
//			{
//				LMM::Index callableTimeIndex = buffer_callableDatesIndex_[itr_callableTimeIndex];
//				price += iv*numeraire[0]/numeraire[callableTimeIndex]; // payoff at time callableTimeIndex, actualized to time 0
//			}
//			//! YY TODO: ???? what happens when never iv>cv ???? need to check 
//		}		
//	}
//	return price /= nb_2ndFwdSimulation_;
//}
//
