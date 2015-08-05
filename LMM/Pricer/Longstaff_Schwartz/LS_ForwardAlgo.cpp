#include <LMM/Pricer/Longstaff_Schwartz/LS_ForwardAlgo.h>
#include <LMM/Pricer/McLmmPricer/McLmmPricer.h>


std::pair<double, double> LS_ForwardAlgo::do_ForwardAlgo(	CallableInstrument_CONSTPTR callableInstrument, 
															const std::vector<McLmm_LS::LMMSimulationResult>& lmmSimulationResult,
															std::vector<LS::Regression>& regressions)
{
	// backward: i = N -> 0 
	//! 1. pathwise: max(vi, vc): depends on product
	//  suppose: vc = 0 at last exercise date. 
	//  vi: depends on product: subSwap ... 
	//  Z_i = max(vi(i), vc(i)*df), vc = E(.I.)
	// 
	//! 2. regression: Z_i(T_i) regression on basis (T_{i-1})
	double price = 0.0;
	double squareSum = 0.0;
	size_t nbSimu = lmmSimulationResult.size();
	std::vector<Instrument_CONSTPTR>   instrument_vector;

	//save the subInstrument for different exercice dates. 
	for(size_t i=0; i<exerciseDates_.size();i++)
	{
		LMM::Index liborExerciceIndex = exerciseDates_[i];
		LMM::Index endInstrumentIndex = callableInstrument->getExerciseTimes().back();
		instrument_vector.push_back(callableInstrument->getSubInstrument(	liborExerciceIndex, 
																			endInstrumentIndex));
	}

	for(size_t pathIndex = 0; pathIndex<nbSimu; pathIndex++)
	{
			double Value_on_onePath = do_ForwardAlgo_on_onePath(callableInstrument,lmmSimulationResult[pathIndex],regressions, instrument_vector);
			price += Value_on_onePath;
			squareSum += std::pow(Value_on_onePath, 2);
	}

	price /= nbSimu;	
	double variance = squareSum/nbSimu - std::pow(price, 2);

	return std::pair<double, double>(price, variance);

}


double LS_ForwardAlgo::do_ForwardAlgo_on_onePath(CallableInstrument_CONSTPTR callableInstrument, 
												 const McLmm_LS::LMMSimulationResult& lmmSimulationResult,
												 std::vector<LS::Regression>& regressions,
												 const std::vector<Instrument_CONSTPTR>&  instrument_vector)
{
	const ublas::matrix<double>& libormatrix=lmmSimulationResult.get_liborMatrix();
	const std::vector<double>& numeraire=lmmSimulationResult.get_numeraire();
	double oneIterationPrice = 0.0;


	for(size_t exerciceIndex=0; exerciceIndex<exerciseDates_.size(); exerciceIndex++)
	{

		LMM::Index liborExerciceIndex = exerciseDates_[exerciceIndex];		
		LMM::Index endInstrumentIndex = callableInstrument->getExerciseTimes().back();

		Instrument_CONSTPTR instrument = instrument_vector[exerciceIndex];

		//! vi
		double vi		=	mcLmmInstrumentPricer_->price_on_oneSimulation(	instrument, 
																			liborExerciceIndex,
																			libormatrix,
																			numeraire);
		//! vc
		double vc=-1e100;
		double df=-1e100;
		if(exerciceIndex!=exerciseDates_.size()-1)
		{
			LS::Regression& regression = regressions[exerciceIndex];
			vc = regression.getRegressionRepresentation().evaluate_representation(lmmSimulationResult);

			LMM::Index nextLiborExerciceIndex = exerciseDates_[exerciceIndex+1];
			df = numeraire[liborExerciceIndex]/numeraire[nextLiborExerciceIndex]; 
		}
		else
		{
			vc=0.0;   //suppose vc at the last exercice date is 0.
			df=0.0;
		}
			
		//! df_vc
		double df_vc = df*vc;

		//check if iv>=cv
		if(vi>=df_vc)															
		{
			double Value_on_onePath=numeraire[0]/numeraire[liborExerciceIndex]*vi;			
			oneIterationPrice=Value_on_onePath;
			break;
		}
	}

	return oneIterationPrice;
}
