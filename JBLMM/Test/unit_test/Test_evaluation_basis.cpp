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

void Test_evaluation_basis()
{
	LMM::Index  indexStart = 2;		//1Y
	LMM::Index  indexEnd   = 20;		//10Y
	Tenor	floatingLegTenorType = Tenor::_6M;
	Tenor	fixedLegTenorType    = Tenor::_1YR;
	assert(indexStart%2==0&&indexEnd%2==0);
	LMMTenorStructure_PTR lmmTenorStructure( new LMMTenorStructure(floatingLegTenorType, indexEnd/2));

	std::vector<std::string> mkt_file_list = InputFileManager::get_VCUB_FileList();
	const std::string& mkt_data_file = mkt_file_list.back();
	std::string folder_name;   // = "TotalCalib\\" ;  config.use_positive_constraint_=true;
	std::string base_name_file = LMMPATH::get_BaseFileName(mkt_data_file) + "\\";
	folder_name+=base_name_file;
	LMMPATH::reset_Output_SubFolder(folder_name );

	LmmCalibrationConfig config;

	config.floatLegTenor_=floatingLegTenorType;
	config.fixedLegTenor_=fixedLegTenorType;

	config.model_nbYear_		=	indexEnd/2;
	size_t fixedFloatRatio		=	config.fixedLegTenor_.ratioTo(config.floatLegTenor_);
	config.correl_FullRank_		=	fixedFloatRatio*config.model_nbYear_+1;

	LmmSwaptionMarketData_PTR pLmmSwaptionMarketData	=	get_LmmSwaptionMarketData(config, mkt_data_file);
	const std::vector<double>&	initLiborValues			=	pLmmSwaptionMarketData->get_LiborQuotes()->get_InitLibor();

	double strike = 0.0137;

	std::vector<LMM::Index> exerciseDates;
	exerciseDates.push_back(2);
	//exerciseDates.push_back(4);
	//exerciseDates.push_back(6);
	//exerciseDates.push_back(8);
	//exerciseDates.push_back(10);
	//exerciseDates.push_back(12);
	//exerciseDates.push_back(14);
	//exerciseDates.push_back(16);
	//exerciseDates.push_back(18);
	exerciseDates.push_back(20);

	McLmm_PTR mcLmm_for_pricer = getMcLmmExample(lmmTenorStructure, initLiborValues, LmmCalibrationConfig());

	size_t fixedFloatingRatio = fixedLegTenorType.ratioTo(floatingLegTenorType);

	std::vector<std::vector<size_t>> subset;
	for(size_t i = 0; i <= 2; i++) 
	{
		subset.push_back(std::vector<size_t>());
		subset.back().push_back(0);
		subset.back().push_back(0);
		subset.back().push_back(i);							
	}
	//for(size_t j = 1; j <= 2; j++) 
	//{
	//	subset.push_back(std::vector<size_t>());
	//	subset.back().push_back(0);
	//	subset.back().push_back(j);
	//	subset.back().push_back(0);
	//}	
/*	for(size_t k = 1; k <= 2; k++) 
	{
		subset.push_back(std::vector<size_t>());
		subset.back().push_back(k);
		subset.back().push_back(0);
		subset.back().push_back(0);
	}	*/			



	size_t regressionIndex=2;
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

	for(size_t basisIndex = 0; basisIndex<subset.size(); basisIndex++)
	{
			basis_vect.push_back(getBasis(subset[basisIndex], 1.0, vanillaSwap, strike, liborIndex));
			basis_evaluator_vect.push_back(getBasisEvaluator(subset[basisIndex], McLmmVanillaSwapPricer(mcLmm_for_pricer)));	
	}

	LS::Regression rg(LS::RegressionRepresentation(basis_vect, basis_evaluator_vect));

	McLmm_PTR mcLmm = getMcLmmExample(lmmTenorStructure, initLiborValues, LmmCalibrationConfig());
	McLmm_LS mcLmm_LS(mcLmm);
	mcLmm_LS.simulateLMM(1);

	size_t nb = 30000;
	clock_t startTime = clock();
	std::vector<std::vector<double>> vect(nb);
	for(size_t i=0; i<nb; i++)
	{
		//rg.getRegressionRepresentation().evaluate_basis(mcLmm_LS.lmmSimualtionResults_[0]);
		vect[i].resize(3);
		vect[i]=rg.getRegressionRepresentation().getBasis_val_buffer();
		vect[i]=std::vector<double>(3, 1.0);
	}

	clock_t endTime = clock();

	clock_t time = endTime - startTime;
	double time_in_second = time/(double) CLOCKS_PER_SEC;

	cout << "time_in_second  "<< time_in_second << endl;

	const matrix& m = mcLmm_LS.lmmSimualtionResults_[0].get_liborMatrix();
	const std::vector<double>& numeraire = mcLmm_LS.lmmSimualtionResults_[0].get_numeraire();

	std::vector<size_t> basis1;
	basis1.push_back(1);
	basis1.push_back(0);
	basis1.push_back(0);
	Basis_CONSTPTR basisA=getBasis(basis1, 1.0, vanillaSwap, strike, liborIndex);
	Basis_Evaluator_CONSTPTR basis_EvaluatorA = getBasisEvaluator(basis1, McLmmVanillaSwapPricer(mcLmm_for_pricer));

	basis_EvaluatorA->evaluate(basisA,m, numeraire);
	clock_t startTime1 = clock();
	for(size_t i=0; i<1000; i++)
		basis_EvaluatorA->evaluate(basisA,m, numeraire);

	clock_t endTime1 = clock();

	clock_t time1 = endTime1 - startTime1;
	double time1_in_second = time1/(double) CLOCKS_PER_SEC;

	cout << "time1_in_second  "<< time1_in_second << endl;

	std::vector<size_t> basis2;
	basis2.push_back(0);
	basis2.push_back(1);
	basis2.push_back(0);
	Basis_CONSTPTR basisB=getBasis(basis2, 1.0, vanillaSwap, strike, liborIndex);
	Basis_Evaluator_CONSTPTR basis_EvaluatorB = getBasisEvaluator(basis2, McLmmVanillaSwapPricer(mcLmm_for_pricer));

	clock_t startTime2 = clock();
	for(size_t i=0; i<1000; i++)
		basis_EvaluatorB->evaluate(basisB, m, numeraire);

	clock_t endTime2 = clock();

	clock_t time2 = endTime2 - startTime2;
	double time2_in_second = time2/(double) CLOCKS_PER_SEC;

	cout << "time2_in_second  "<< time2_in_second << endl;

	std::vector<size_t> basis3;
	basis3.push_back(0);
	basis3.push_back(0);
	basis3.push_back(1);
	Basis_CONSTPTR basisC=getBasis(basis3, 1.0, vanillaSwap, strike, liborIndex);
	Basis_Evaluator_CONSTPTR basis_EvaluatorC = getBasisEvaluator(basis3, McLmmVanillaSwapPricer(mcLmm_for_pricer));

	
	clock_t startTime3 = clock();
	for(size_t i=0; i<1000; i++)
		basis_EvaluatorC->evaluate(basisC, m, numeraire);

	clock_t endTime3 = clock();

	clock_t time3 = endTime3 - startTime3;
	double time3_in_second = time3/(double) CLOCKS_PER_SEC;

	cout << "time3_in_second  "<< time3_in_second << endl;

	EV_Evaluator_CONSTPTR ev_evaluator_libor(new EV_LiborRate_Evaluator());
	EV_Evaluator_CONSTPTR ev_evaluator_swaprate(new EV_VanillaSwapRate_Evaluator(McLmmVanillaSwapPricer(mcLmm_for_pricer)));

	EV_CONSTPTR ev_swaprate(new EV_VanillaSwapRate( Rate1_CONSTPTR( new VanillaSwapRate(VanillaSwap(vanillaSwap)))));

	EV_CONSTPTR ev_libor(new EV_LiborRate(Rate1_CONSTPTR(new LiborRate(2,Tenor(Tenor::_6M)))));

	clock_t startTime4 = clock();
	for(size_t i=0; i<1000; i++)
		ev_evaluator_swaprate->evaluate(ev_swaprate,m,numeraire);

	clock_t endTime4 = clock();

	clock_t time4 = endTime4 - startTime4;
	double time4_in_second = time4/(double) CLOCKS_PER_SEC;

	cout << "time4_in_second  "<< time4_in_second << endl;



	clock_t startTime5 = clock();
	for(size_t i=0; i<1000; i++)
		ev_evaluator_libor->evaluate(ev_libor,m,numeraire);

	clock_t endTime5 = clock();

	clock_t time5 = endTime5 - startTime5;
	double time5_in_second = time5/(double) CLOCKS_PER_SEC;

	cout << "time5_in_second  "<< time5_in_second << endl;

}