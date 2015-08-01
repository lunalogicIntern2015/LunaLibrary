#include <JBLMM/Test/JBTests.h>
#include <LMM/Test/Tests.h>

#include <iostream>
#include <fstream> 
#include <string> 
#include <iterator>
#include <algorithm>

#include <boost/math/distributions/inverse_gaussian.hpp>
#include <boost/math/distributions/normal.hpp>

#include <LMM/RNGenerator/McGenerator.h>
#include <LMM/RNGenerator/RNGenerator.h>
#include <LMM/helper/InputFileManager.h>
#include <LMM/numeric/NumericalMethods.h>
#include <LMM/LmmModel/Correlation.h>
#include <LMM/LmmModel/McTerminalLmm.h>
#include <LMM/LmmModel/Lmm.h>
#include <LMM/LmmModel/ConstShifted_HGVolatilityFunction.h>
#include <LMM/LmmModel/LmmSwaptionMarketData.h>
#include <LMM/pricer/McLmmVanillaSwaptionPricer.h>

#include <JBLMM/Element/Rate1.h>  
#include <JBLMM/Element/ConstRate.h>  
#include <JBLMM/Element/LiborRate.h>  
#include <JBLMM/Element/VanillaSwapRate.h>  
#include <JBLMM/Longstaff_Schwartz/EV_Basis/Basis.h>
#include <JBLMM/Longstaff_Schwartz/EV_Basis/Basis_Evaluator.h>
#include <JBLMM/Longstaff_Schwartz/Regression/Regression_LS.h>
#include <JBLMM/Longstaff_Schwartz/Simulation/McLmm_LS.h>
#include <JBLMM/Instrument/CallableInstrument.h>
#include <JBLMM/Instrument/CallableSwap.h>
#include <JBLMM/Instrument/InstrumentFactory.h>
#include <JBLMM/Pricer/McLmmPricer.h>
#include <JBLMM/Pricer/McLmmGenericSwapPricer.h>
#include <JBLMM/Longstaff_Schwartz/LS_BackwardAlgo.h>
#include <JBLMM/Longstaff_Schwartz/LS_ForwardAlgo.h>

void Test_LS_pricing_One_SubSet_basis(	const std::vector<std::vector<size_t>>& basis_subset, 
										double strike,
										LMM::Index  indexStart, 
										LMM::Index  indexEnd,
										const Tenor&	floatingLegTenorType,
										const Tenor&	fixedLegTenorType,
										const std::vector<double>& initLiborValues,
										const std::vector<LMM::Index>& exerciseDates,
										CallableInstrument_CONSTPTR callableGenericSwap,
										const std::vector<McLmm_LS::LMMSimulationResult>&  lmmSimualtionResults_backward,
										const std::vector<McLmm_LS::LMMSimulationResult>&  lmmSimualtionResults_forward,
										const std::vector<LMM::Index>& nbSimulation_vect,
										std::vector<std::vector<double>>& basis_value_on_allPath_buffer,
										ofstream& out,
										ofstream& out_test_time)
{

	assert(indexStart%2==0&&indexEnd%2==0);
	LMMTenorStructure_PTR lmmTenorStructure(new LMMTenorStructure(floatingLegTenorType, indexEnd/2));

	McLmm_PTR mcLmm_for_pricer = getMcLmmExample(lmmTenorStructure, initLiborValues, LmmCalibrationConfig());

	double nominal = 1.0;
	std::vector<std::vector<std::vector<size_t>>> subset;

	LS_BackwardAlgo ls_BackwardAlgo(McLmmPricer_CONSTPTR(new McLmmVanillaSwapPricer(mcLmm_for_pricer)));

	std::vector<LS::Regression>& regression_vect = ls_BackwardAlgo.getRegressions();
	ls_BackwardAlgo.getExerciseDates()=exerciseDates;
	

	size_t fixedFloatingRatio = fixedLegTenorType.ratioTo(floatingLegTenorType);
	for(size_t regressionIndex=0; regressionIndex<exerciseDates.size()-1; regressionIndex++)
	{
		LMM::Index liborIndex	= indexStart + regressionIndex*fixedFloatingRatio;
		LMM::Index paymentIndex = indexStart + regressionIndex*fixedFloatingRatio + 1;
		VanillaSwap vanillaSwap(	strike, 
									liborIndex, 
									indexEnd, 
									floatingLegTenorType, 
									fixedLegTenorType, 
									lmmTenorStructure);	
		std::vector<Basis_CONSTPTR> basis_vect;
		std::vector<Basis_Evaluator_CONSTPTR> basis_evaluator_vect;

		for(size_t basisIndex = 0; basisIndex<basis_subset.size(); basisIndex++)
		{
				basis_vect.push_back(getBasis(basis_subset[basisIndex], 1.0, vanillaSwap, strike, liborIndex));
				basis_evaluator_vect.push_back(getBasisEvaluator(basis_subset[basisIndex], McLmmVanillaSwapPricer(mcLmm_for_pricer)));	
		}

		regression_vect.push_back(LS::Regression(LS::RegressionRepresentation(basis_vect, basis_evaluator_vect)));
	}

	// ----------------------------------------------------------------------
	// 
	//                     LS_forward
	//
	// ---------------------------------------------------------------------
	LS_ForwardAlgo ls_ForwardAlgo(McLmmPricer_CONSTPTR(new McLmmVanillaSwapPricer(mcLmm_for_pricer)));
	ls_ForwardAlgo.getExerciseDates()=exerciseDates;
	// ----------------------------------------------------------------------
	// 
	//                     LS_Algorithme
	//
	// ---------------------------------------------------------------------
	std::vector<double> price_vect;
	std::vector<double> leftValue_vect;
	std::vector<double> rightValue_vect;
	std::vector<double> time_vect;
	std::vector<double> variance_vect;
	std::vector<double> diff_vect;

	// ----------------------------------------------------------------------
	// 
	//               print on file
	//
	// ----------------------------------------------------------------------
	representation_basis_subset(basis_subset,out_test_time); out << endl;

	for (size_t i = 0; i<nbSimulation_vect.size(); i++)
	{
		clock_t startTime = clock();

		size_t nbSimulation = nbSimulation_vect[i];
		   //the same seed for each simulation.
		clock_t do_BackwardAlgo_startTime = clock();
		ls_BackwardAlgo.do_BackwardAlgo(callableGenericSwap,lmmSimualtionResults_backward,basis_value_on_allPath_buffer);			
		
		clock_t forward_simulateLMM_startTime = clock();

		clock_t do_ForwardAlgo_startTime = clock();
		std::pair<double, double> result = ls_ForwardAlgo.do_ForwardAlgo(	callableGenericSwap, 
																			lmmSimualtionResults_forward,
																			ls_BackwardAlgo.getRegressions());
		clock_t endTime = clock();

		double price = result.first;
		double variance = result.second;
		boost::math::normal dist(0.0, 1.0);
		double threshold = 0.99;			// 95% of distribution is below q:
		double q = quantile(dist, threshold);
		double leftValue = price - sqrt(variance/nbSimulation)*q;
		double rightValue = price + sqrt(variance/nbSimulation)*q;



#ifdef TIME_TEST		
		clock_t getMcLmmExample_Time		= do_BackwardAlgo_startTime - startTime;	
		clock_t do_BackwardAlgo_Time		= forward_simulateLMM_startTime - do_BackwardAlgo_startTime;
		clock_t forward_simulateLMM_Time	= do_ForwardAlgo_startTime - forward_simulateLMM_startTime;	
		clock_t do_ForwardAlgo_Time			= endTime - do_ForwardAlgo_startTime;
		
		double getMcLmmExample_Second		= getMcLmmExample_Time / (double) CLOCKS_PER_SEC;
		double do_BackwardAlgo_Second		= do_BackwardAlgo_Time / (double) CLOCKS_PER_SEC;
		double forward_simulateLMM_Second	= forward_simulateLMM_Time / (double) CLOCKS_PER_SEC;
		double do_ForwardAlgo_Second		= do_ForwardAlgo_Time / (double) CLOCKS_PER_SEC;

		cout << "do_BackwardAlgo_Second " << do_BackwardAlgo_Second << endl;
		cout << "forward_simulateLMM_Second " << forward_simulateLMM_Second << endl;
#endif			
		clock_t clockTicksTaken			= endTime - startTime;	
		double timeInSeconds = clockTicksTaken / (double) CLOCKS_PER_SEC;

		representation_basis_subset(basis_subset,out); out << "; ";
		out << price << "; " << leftValue << "; " << rightValue << "; " << sqrt(variance/nbSimulation)*q << "; ";
		out << timeInSeconds << "; " << endl;

		cout << "price " << price << endl; 
		cout << "IC: " << "[ "<<leftValue << ", " << rightValue << "]" << endl;
		cout << "demi_longueur_IC : " << sqrt(variance/nbSimulation)*q << endl;
		cout <<  "time : " << timeInSeconds << endl;
	}
}


void representation_basis(const std::vector<size_t>& basis_vect, ofstream& out)
{
	assert(basis_vect.size()==3);
	size_t swaprate_degree			= basis_vect[0];
	size_t payoff_swaprate_degree	= basis_vect[1];
	size_t liborRate_degree			= basis_vect[2];
	if(swaprate_degree==0&&payoff_swaprate_degree==0&&liborRate_degree==0)
		out << 1 << " ";
	if(swaprate_degree>0&&payoff_swaprate_degree==0&&liborRate_degree==0)
	{
		for(size_t i = 0; i<swaprate_degree; i++)
			out << "a";
		out << " ";
	}
	if(swaprate_degree==0&&payoff_swaprate_degree>0&&liborRate_degree==0)
	{
		for(size_t i = 0; i<payoff_swaprate_degree; i++)
			out << "b";
		out << " ";
	}
	if(swaprate_degree==0&&payoff_swaprate_degree==0&&liborRate_degree>0)
	{
		for(size_t i = 0; i<liborRate_degree; i++)
			out << "c";
		out << " ";
	}
	if(swaprate_degree==1&&payoff_swaprate_degree==1&&liborRate_degree==0)
		out << "ab" << " ";
	if(swaprate_degree==0&&payoff_swaprate_degree==1&&liborRate_degree==1)
		out << "bc" << " ";
	if(swaprate_degree==1&&payoff_swaprate_degree==0&&liborRate_degree==1)
		out << "ac" << " ";
}

void representation_basis_subset(const std::vector<std::vector<size_t>>& basis_subset, ofstream& out)
{
	for(size_t i=0; i<basis_subset.size(); i++)
	{
		representation_basis(basis_subset[i], out);
		out << " ";
	}
}




//void Test_LS_pricing_One_SubSet_basis(	const std::vector<std::vector<size_t>>& basis_subset, 
//										double strike,
//										LMM::Index  indexStart, 
//										LMM::Index  indexEnd,
//										const Tenor&	floatingLegTenorType,
//										const Tenor&	fixedLegTenorType,
//										const std::vector<double>& initLiborValues,
//										const std::vector<LMM::Index>& exerciseDates,
//										CallableInstrument_CONSTPTR callableGenericSwap,
//										const std::vector<McLmm_LS::LMMSimulationResult>&  lmmSimualtionResults_backward,
//										const std::vector<McLmm_LS::LMMSimulationResult>&  lmmSimualtionResults_forward,
//										const std::vector<LMM::Index>& nbSimulation_vect,
//										ofstream& out,
//										ofstream& out_test_time)
//{
//	assert(indexStart%2==0&&indexEnd%2==0);
//	LMMTenorStructure_PTR lmmTenorStructure(new LMMTenorStructure(floatingLegTenorType, indexEnd/2));
//
//	McLmm_PTR mcLmm_for_pricer = getMcLmmExample(lmmTenorStructure, initLiborValues);
//
//	double nominal = 1.0;
//
//	//////! Parameter of h
//	//double a = 0.105;
//	//double b = 0.317;
//	//double c = 0.515;
//	//double d = 0.356;
//	//Shifted_HGVolatilityParam::ABCDParameter abcdParam (a,b,c,d);
//	////Parameter of hg
//	//double g_constParam = 1.0;
//	//double shift_constParam = 0.0;
//	//ConstShifted_HGVolatilityParam_PTR hgParam( new ConstShifted_HGVolatilityParam(lmmTenorStructure,abcdParam,g_constParam,shift_constParam));
//	//
//	////! Correlation 1
//	//size_t nbFactor       = 3; // need to test nbFactor  = 3, and nbFactor = 
//	//size_t correlFullRank = lmmTenorStructure->get_horizon()+1;
//	//size_t correlReducedRank = nbFactor;
//	////!"Check Parameters(): Condition not implemented yet."
//	//std::cout << "checkParams(): ";
//	//CorrelationReductionType::CorrelationReductionType correlReductionType = CorrelationReductionType::PCA;
//	//double correlAlpha = 0.0;
//	//double correlBeta  = 0.1;
//	//Correlation_PTR correlation(new XY_beta_Correlation(correlFullRank,correlReducedRank, correlReductionType,correlAlpha,correlBeta));
//	//correlation->calculate(); // for print.
//	//correlation->print("test_McTerminalLmm_Correlation.csv");
//	////hgVolatilityFunction
//	//ConstShifted_HGVolatilityFunction_PTR hgVolatilityFunction (new ConstShifted_HGVolatilityFunction(lmmTenorStructure, correlation, hgParam)); 
//	//hgVolatilityFunction->print("test_McTerminalLmm_Volatility.csv");
//	////! Dispersion
//	//Dispersion dispersion(hgVolatilityFunction);
//
//	//unsigned long seed = 5033;
//	//RNGenerator_PTR  rnGenerator(new McGenerator(seed));
//
//	////build lmm and mcLmm model
//	//Lmm_PTR shiftedLmm (new Lmm(dispersion));
//

//
//	//vector<std::vector<size_t>> set;
//	////set.push_back(std::vector<size_t>(3,0));
//	//for (size_t i= 0 ; i<3; i++)
//	//{
//	//	//set.push_back(std::vector<size_t>(3,0));
//	//	//set.back()[i]+=1;
//	//	for (size_t j=0; j<=i; j++)
//	//	{
//	//		set.push_back(std::vector<size_t>(3,0));
//	//		set.back()[i]+=1;
//	//		set.back()[j]+=1;
//	//	}
//	//}
//
//	std::vector<std::vector<std::vector<size_t>>> subset;
//	//getAllSubsets(set, subset);
//
//	// ----------------------------------------------------------------------
//	// 
//	//                     LS_Backward
//	//
//	// ---------------------------------------------------------------------
//	
//
//	// ----------------------------------------------------------------------
//	// 
//	//                     LS_Backward
//	//
//	// ----------------------------------------------------------------------
//
//
//	LS_BackwardAlgo ls_BackwardAlgo(McLmmPricer_CONSTPTR(new McLmmVanillaSwapPricer(mcLmm_for_pricer)));
//
//	std::vector<LS::Regression>& regression_vect = ls_BackwardAlgo.getRegressions();
//	ls_BackwardAlgo.getExerciseDates()=exerciseDates;
//	
//
//	size_t fixedFloatingRatio = fixedLegTenorType.ratioTo(floatingLegTenorType);
//	for(size_t regressionIndex=0; regressionIndex<exerciseDates.size()-1; regressionIndex++)
//	{
//		LMM::Index liborIndex	= indexStart + regressionIndex*fixedFloatingRatio;
//		LMM::Index paymentIndex = indexStart + regressionIndex*fixedFloatingRatio + 1;
//		//LMM::Index paymentIndex = indexStart + regressionIndex*fixedFloatingRatio + 1;
//		//std::vector<Basis_CONSTPTR> basis_vect;
//		//basis_vect.push_back(Basis_CONSTPTR(new Basis_CappedFlooredCoupon(	
//		//	EV_CONSTPTR(new EV_VanillaSwapRate( Rate1_CONSTPTR( new VanillaSwapRate(
//		//					VanillaSwap(			strike, 
//		//											indexStart, 
//		//											indexEnd, 
//		//											floatingLegTenorType, 
//		//											fixedLegTenorType, 
//		//											lmmTenorStructure))))), 
//		//	Basis_SinglaEVFunctional::Transformer_CONSTPTR(
//		//						new Basis_CappedFlooredCoupon::CappedFlooredCouponTransformer(
//		//							CappedFlooredCoupon_CONSTPTR(new CappedFlooredCoupon(CappedFlooredCoupon(	
//		//																										paymentIndex,
//		//																										1.0,
//		//																										floatingLegTenorType.YearFraction(), 
//		//																										true,
//		//																										0.0,
//		//																										false, 
//		//																										10e100, 
//		//																										Rate1_CONSTPTR(new LiborRate(liborIndex,Tenor(floatingLegTenorType))),
//		//																										1.0,
//		//																										-strike, 
//		//																										liborIndex))))))));
//	
//		//std::vector<Basis_Evaluator_CONSTPTR> basis_evaluator_vect;
//		//basis_evaluator_vect.push_back(	Basis_Evaluator_CONSTPTR(new Basis_SinglaEVFunctional_Evaluator(
//		//									EV_Evaluator_CONSTPTR(new EV_VanillaSwapRate_Evaluator(McLmmVanillaSwapPricer(mcLmm))))));	
//
//		VanillaSwap vanillaSwap(	strike, 
//									liborIndex, 
//									indexEnd, 
//									floatingLegTenorType, 
//									fixedLegTenorType, 
//									lmmTenorStructure);
//		
//		std::vector<Basis_CONSTPTR> basis_vect;
//		std::vector<Basis_Evaluator_CONSTPTR> basis_evaluator_vect;
//
//		for(size_t basisIndex = 0; basisIndex<basis_subset.size(); basisIndex++)
//		{
//				basis_vect.push_back(getBasis(basis_subset[basisIndex], 1.0, vanillaSwap, strike, liborIndex));
//				basis_evaluator_vect.push_back(getBasisEvaluator(basis_subset[basisIndex], McLmmVanillaSwapPricer(mcLmm_for_pricer)));	
//		}
//
//		regression_vect.push_back(LS::Regression(LS::RegressionRepresentation(basis_vect, basis_evaluator_vect)));
//	}
//
//	// ----------------------------------------------------------------------
//	// 
//	//                     LS_forward
//	//
//	// ---------------------------------------------------------------------
//	LS_ForwardAlgo ls_ForwardAlgo(McLmmPricer_CONSTPTR(new McLmmVanillaSwapPricer(mcLmm_for_pricer)));
//	ls_ForwardAlgo.getExerciseDates()=exerciseDates;
//
//
//
//	// ----------------------------------------------------------------------
//	// 
//	//                     LS_Algorithme
//	//
//	// ---------------------------------------------------------------------
//
//	//std::vector<size_t> nbSimulation_vect;
//	//nbSimulation_vect.push_back(20);
//	//nbSimulation_vect.push_back(300);
//	//nbSimulation_vect.push_back(1000);
//	//nbSimulation_vect.push_back(2000);
//	//nbSimulation_vect.push_back(5000);
//	//nbSimulation_vect.push_back(10000);
//	//nbSimulation_vect.push_back(20000);
//	//nbSimulation_vect.push_back(50000);
//	//nbSimulation_vect.push_back(100000);
//	//nbSimulation_vect.push_back(200000);
//	//nbSimulation_vect.push_back(500000);
//
//	std::vector<double> price_vect;
//	std::vector<double> leftValue_vect;
//	std::vector<double> rightValue_vect;
//	std::vector<double> time_vect;
//	std::vector<double> variance_vect;
//	std::vector<double> diff_vect;
//
//	// ----------------------------------------------------------------------
//	// 
//	//               print on file
//	//
//	// ----------------------------------------------------------------------
//
//	for (size_t i = 0; i<nbSimulation_vect.size(); i++)
//	{
//
//		clock_t startTime = clock();
//
//		size_t nbSimulation = nbSimulation_vect[i];
//
//		clock_t simulateLMM_startTime = clock();
//		   //the same seed for each simulation.
//
//
//		clock_t do_BackwardAlgo_startTime = clock();
//
//		//McLmm_PTR mcLmm = getMcLmmExample(lmmTenorStructure, initLiborValues);
//		//McLmm_LS mcLmm_LS(mcLmm);
//		//mcLmm_LS.simulateLMM(nbSimulation);
//		ls_BackwardAlgo.do_BackwardAlgo(callableGenericSwap,lmmSimualtionResults_backward);
//		
//		clock_t forward_simulateLMM_startTime = clock();
//
//		clock_t do_ForwardAlgo_startTime = clock();
//		std::pair<double, double> result = ls_ForwardAlgo.do_ForwardAlgo(	callableGenericSwap, 
//																			lmmSimualtionResults_forward,
//																			ls_BackwardAlgo.getRegressions());
//
//		clock_t endTime = clock();
//
//#ifdef TIME_TEST		
//		clock_t getMcLmmExample_Time		= simulateLMM_startTime - startTime;
//		clock_t simulateLMM_Time			= do_BackwardAlgo_startTime - simulateLMM_startTime;	
//		clock_t do_BackwardAlgo_Time		= forward_simulateLMM_startTime - do_BackwardAlgo_startTime;
//		clock_t forward_simulateLMM_Time	= do_ForwardAlgo_startTime - forward_simulateLMM_startTime;	
//		clock_t do_ForwardAlgo_Time			= endTime - do_ForwardAlgo_startTime;
//		
//		double getMcLmmExample_Second		= getMcLmmExample_Time / (double) CLOCKS_PER_SEC;
//		double simulateLMM_Second			= simulateLMM_Time / (double) CLOCKS_PER_SEC;
//		double do_BackwardAlgo_Second		= do_BackwardAlgo_Time / (double) CLOCKS_PER_SEC;
//		double forward_simulateLMM_Second	= forward_simulateLMM_Time / (double) CLOCKS_PER_SEC;
//		double do_ForwardAlgo_Second		= do_ForwardAlgo_Time / (double) CLOCKS_PER_SEC;
//
//		cout << "getMcLmmExample_Second " << getMcLmmExample_Second << endl;
//		cout << "simulateLMM_Second " << simulateLMM_Second << endl;
//		cout << "do_BackwardAlgo_Second " << do_BackwardAlgo_Second << endl;
//		cout << "forward_simulateLMM_Second " << forward_simulateLMM_Second << endl;
//		cout << "do_ForwardAlgo_Second " << do_ForwardAlgo_Second << endl;
//
//#endif
//			
//		clock_t clockTicksTaken			= endTime - startTime;	
//		double timeInSeconds = clockTicksTaken / (double) CLOCKS_PER_SEC;
//
//		double price = result.first;
//		double variance = result.second;
//		boost::math::normal dist(0.0, 1.0);
//		double threshold = 0.99;			// 95% of distribution is below q:
//		double q = quantile(dist, threshold);
//		double leftValue = price - sqrt(variance/nbSimulation)*q;
//		double rightValue = price + sqrt(variance/nbSimulation)*q;
//
//
//
//		//price_vect.push_back(result.first);
//		//variance_vect.push_back(result.second);
//		//leftValue_vect.push_back(leftValue);
//		//rightValue_vect.push_back(rightValue);
//		//time_vect.push_back(timeInSeconds);
//		//diff_vect.push_back(sqrt(variance/nbSimulation)*q);
//
//	}
//
//	//save 
//
//
//	//out << "Basis: ;" << endl;  
//	//for(size_t basisIndex = 0; basisIndex<basis_subset.size(); basisIndex++)
//	//{
//	//	for(size_t i=0; i<3; i++)
//	//		out << basis_subset[basisIndex][i] << " ;" ;
//	//	out << endl;
//	//}
//	//out << endl;
//	//out << " 8 basis ;   " << endl;
//	//out << " 0.001 ;  " << "L_i ;" <<  "(L_i)^2 ;" <<"max(L_i-K,0) ;" << "S_i ;" <<  "(S_i)^2 ;";  // << "(S_i)^3 ;" <<  "(S_i)^4 ;" ;
//	//out <<"max(S_i-K, 0) ;" <<"(max(S_i-K, 0))^2 ;"<< endl;
//
//	//size_t nbSimu_size = nbSimulation_vect.size();
//
//	////out << "Basis_Monomial; on Libor Rate: ; L_i" << endl;
//
//	//out << "exerciseDates: ;" ;
//	//for(size_t i = 0; i<exerciseDates.size(); i++)
//	//	out << exerciseDates[i] << " ;" ;
//	//out << endl;
//	//
//	//out << "nbSimulation_vect: ;" ;
//	//for(size_t i = 0; i<nbSimulation_vect.size(); i++)
//	//	out << nbSimulation_vect[i] << " ;" ;
//	//out << endl;
//
//
//	//out << "leftValue_vect: ;" ;
//	//for(size_t i = 0; i<nbSimulation_vect.size(); i++)
//	//	out << leftValue_vect[i] << " ;" ;
//	//out << endl;
//
//	//out << "price_vect: ;" ;
//	//for(size_t i = 0; i<nbSimulation_vect.size(); i++)
//	//	out << price_vect[i] << " ;" ;
//	//out << endl;
//
//	//out << "rightValue_vect: ;" ;
//	//for(size_t i = 0; i<nbSimulation_vect.size(); i++)
//	//	out << rightValue_vect[i] << "; " ;
//	//out << endl;
//
//	//out << "diff_vect: ;" ;
//	//for(size_t i = 0; i<nbSimulation_vect.size(); i++)
//	//	out << diff_vect[i] << " ;" ;
//	//out << endl;
//
//	//out << "variance_vect: ;" ;
//	//for(size_t i = 0; i<nbSimulation_vect.size(); i++)
//	//	out << variance_vect[i] << " ;" ;
//	//out << endl;
//
//	//out << "time_vect: ;" ;
//	//for(size_t i = 0; i<nbSimulation_vect.size(); i++)
//	//	out << time_vect[i] << " ;" ;
//	//out << endl;
//
//
//	//// ----------------------------------------------------------------------
//	//// 
//	////               show on terminal
//	////
//	//// ----------------------------------------------------------------------
//
//	//cout << "Basis: ;" << endl;  
//	//for(size_t basisIndex = 0; basisIndex<basis_subset.size(); basisIndex++)
//	//{
//	//	for(size_t i=0; i<3; i++)
//	//		cout << basis_subset[basisIndex][i] << " ;" ;
//	//	cout << endl;
//	//}
//	//cout << endl;
//	//cout << "exerciseDates: " ;
//	//for(size_t i = 0; i<exerciseDates.size(); i++)
//	//	cout << exerciseDates[i] << " " ;
//	//cout << endl;
//
//	//cout << "nbSimulation_vect: " ;
//	//for(size_t i = 0; i<nbSimulation_vect.size(); i++)
//	//	cout << nbSimulation_vect[i] << " " ;
//	//cout << endl;
//
//	//cout << "leftValue_vect: " ;
//	//for(size_t i = 0; i<nbSimulation_vect.size(); i++)
//	//	cout << leftValue_vect[i] << " " ;
//	//cout << endl;
//
//	//cout << "price_vect: " ;
//	//for(size_t i = 0; i<nbSimulation_vect.size(); i++)
//	//	cout << price_vect[i] << " " ;
//	//cout << endl;
//
//	//cout << "rightValue_vect: " ;
//	//for(size_t i = 0; i<nbSimulation_vect.size(); i++)
//	//	cout << rightValue_vect[i] << " " ;
//	//cout << endl;
//
//	//cout << "diff_vect: " ;
//	//for(size_t i = 0; i<nbSimulation_vect.size(); i++)
//	//	cout << diff_vect[i] << " " ;
//	//cout << endl;
//
//	//cout << "variance_vect: " ;
//	//for(size_t i = 0; i<nbSimulation_vect.size(); i++)
//	//	cout << variance_vect[i] << " " ;
//	//cout << endl;
//
//	//cout << "time_vect: " ;
//	//for(size_t i = 0; i<nbSimulation_vect.size(); i++)
//	//	cout << time_vect[i] << " " ;
//	//cout << endl;
//}