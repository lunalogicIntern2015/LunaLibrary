#include <JBLMM/Longstaff_Schwartz/LS_BackwardAlgo.h>


void LS_BackwardAlgo::update_Z( CallableInstrument_CONSTPTR callableInstrument, 
								LMM::Index exerciseLiborDateindex,
								LMM::Index regressionIndex,
								const std::vector<McLmm_LS::LMMSimulationResult>& lmmSimulationResult)  
								// Z_buffer_ = max(iv,df*cv) pathwise
{
	if(lmmSimulationResult.size() != Z_buffer_.size())
		throw("LS_BackwardAlgo do_BackwardAlgo size mismatch ");

	Instrument_CONSTPTR instrument = callableInstrument->getSubInstrument(	exerciseLiborDateindex,  //startIndex
																			callableInstrument->getExerciseTimes().back()); // endIndex 


	for(size_t path_index = 0; path_index < Z_buffer_.size(); ++path_index)
	{
		const std::vector<double>& numeraire = lmmSimulationResult[path_index].get_numeraire();
		//! vi  
		// suppose the last exercise index is the end index of the underlying instrument
		double vi=0.0;
		vi = mcLmmInstrumentPricer_->price_on_oneSimulation(	instrument, 
																exerciseLiborDateindex, 
																lmmSimulationResult[path_index].get_liborMatrix(),
																numeraire);

		//! vc
		double vc;
		double df;
		if(regressionIndex!=exerciseDates_.size()-1)
		{
			LS::Regression& regression = regressions_[regressionIndex];
			vc = regression.getRegressionRepresentation().evaluate_representation(lmmSimulationResult[path_index]);	// E[Z_N|S_{N-1}]

			LMM::Index nextExerciseLiborDateindex = exerciseDates_[regressionIndex+1];
			df = numeraire[exerciseLiborDateindex]/numeraire[nextExerciseLiborDateindex];							// B(T_{N-1}, T_N)
		}
		else
		{
			vc=0.0;
			df=0.0;
		}
		Z_buffer_[path_index] = std::max(vi, df*vc);
	}
}

void LS_BackwardAlgo::do_BackwardAlgo(	CallableInstrument_CONSTPTR callableInstrument, 
										const std::vector<McLmm_LS::LMMSimulationResult>& lmmSimulationResult,
										std::vector<std::vector<double>>& basis_value_on_allPath_buffer)
{	
	// backward: i = N -> 0 
	//! 1. pathwise: max(vi, vc): depends on product
	//  suppose: vc = 0 at last exercise date. 
	//  vi: depends on product: subSwap ... 
	//  Z_i = max(vi(i), vc(i)*df), vc = E(.I.)
	// 
	//! 2. regression: Z_i(T_i) regression on basis (T_{i-1})
	if(exerciseDates_.size() > 1)
	{
		for(size_t i = exerciseDates_.size()-1; i>0;  --i)
		{
			//! regression from t2 -> t1
			LMM::Index t2_exerciseLiborDateindex = exerciseDates_[i];
			LMM::Index t1_exerciseLiborDateindex = exerciseDates_[i-1];

			Z_buffer_.resize(lmmSimulationResult.size());

			//LS::Regression& regression_t2 = regressions_[i];

			clock_t update_Z_startTime = clock();

			update_Z(	callableInstrument, 
						t2_exerciseLiborDateindex,
						i,
						lmmSimulationResult); // at t_2, Z_buffer_ = max(iv,cv) pathwise                                               

			LS::Regression& regression_t1 = regressions_[i-1];
			regression_t1.regress(Z_buffer_, lmmSimulationResult, basis_value_on_allPath_buffer);  // regression: Z_i(t2) regression on basis (t1)
		}
	}
	else
	{
		throw("LS_BackwardAlgo::do_BackwardAlgo(), this is european opion, need to treat it in a more elegant way: other than if-else" );
	}
}


