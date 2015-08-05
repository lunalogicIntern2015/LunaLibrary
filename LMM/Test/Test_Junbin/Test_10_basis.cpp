#include <JBLMM/Test/JBTests.h>
#include <LMM/Test/Tests.h>

#include <iostream>
#include <fstream> 
#include <string> 
#include <iterator>
#include <algorithm>

#include <boost/math/distributions/inverse_gaussian.hpp>
#include <boost/math/distributions/normal.hpp>

#include <RNGenerator/McGenerator.h>
#include <RNGenerator/RNGenerator.h>
#include <LMM/Helper/InputFileManager.h>
#include <Numeric/NumericalMethods.h>
#include <LMM/Model/Correlation.h>
#include <LMM/Mc/McTerminalLmm.h>
#include <LMM/Model/Lmm.h>
#include <LMM/Model/ConstShifted_HGVolatilityFunction.h>
#include <LMM/LmmSwaptionMarketData.h>
#include <LMM/Pricer/McLmmPricer/McLmmVanillaSwaptionPricer.h>

#include <Instrument/Rate/Rate1.h>  
#include <Instrument/Rate/ConstRate.h>  
#include <Instrument/Rate/LiborRate.h>  
#include <Instrument/Rate/VanillaSwapRate.h>  
#include <LMM/Pricer/Longstaff_Schwartz/Basis.h>
#include <LMM/Pricer/Longstaff_Schwartz/Basis_Evaluator.h>
#include <LMM/Pricer/Longstaff_Schwartz/Regression_LS.h>
#include <LMM/Pricer/Longstaff_Schwartz/McLmm_LS.h>
#include <Instrument/CallableOption/CallableInstrument.h>
#include <JBInstrument/CallableSwap.h>
#include <JBInstrument/InstrumentFactory.h>
#include <LMM/Pricer/McLmmPricer/McLmmPricer.h>
#include <JBLMM/Pricer/McLmmGenericSwapPricer.h>
#include <LMM/Pricer/Longstaff_Schwartz/LS_BackwardAlgo.h>
#include <LMM/Pricer/Longstaff_Schwartz/LS_ForwardAlgo.h>


//volatility
//correlation
//nb Year
//nb Exercise

void Test_LS_pricing_parameter()
{
	LmmCalibrationConfig config;
	config.model_nbYear_=20;
	config.nbSimulation_=30;

	//config.g = 1.0; config.correl_beta_=0.1;
	//Test_LS_pricing_vol_correl(config);

	config.g = 1.50; config.correl_beta_=0.1;
	Test_LS_pricing_vol_correl(config);

	//config.g = 1.0; config.correl_beta_=0.5;
	//Test_LS_pricing_vol_correl(config);

	//config.correl_beta_=0.1; config.g = 0.5;
	//Test_LS_pricing_vol_correl(config);

	//config.correl_beta_=0.1; config.g = 2.0;
	//Test_LS_pricing_vol_correl(config);
}

void Test_LS_pricing_vol_correl(const LmmCalibrationConfig& config)
{
	LMM::Index  indexStart = 2;		//1Y
	LMM::Index  indexEnd   = 2*config.model_nbYear_;		//10Y
	Tenor	floatingLegTenorType = Tenor::_6M;
	Tenor	fixedLegTenorType    = Tenor::_1YR;
	assert(indexStart%2==0&&indexEnd%2==0);
	LMMTenorStructure_PTR lmmTenorStructure( new LMMTenorStructure(floatingLegTenorType, config.model_nbYear_));

	std::vector<std::string> mkt_file_list = InputFileManager::get_VCUB_FileList();
	const std::string& mkt_data_file = mkt_file_list.back();
	std::string folder_name;   // = "TotalCalib\\" ;  config.use_positive_constraint_=true;
	std::string base_name_file = LMMPATH::get_BaseFileName(mkt_data_file) + "\\";
	folder_name+=base_name_file;
	LMMPATH::reset_Output_SubFolder(folder_name );

	size_t fixedFloatRatio		=	fixedLegTenorType.ratioTo(floatingLegTenorType);

	LmmSwaptionMarketData_PTR pLmmSwaptionMarketData	=	get_LmmSwaptionMarketData(config, mkt_data_file);
	const std::vector<double>&	initLiborValues			=	pLmmSwaptionMarketData->get_LiborQuotes()->get_InitLibor();

	//config.correl_ReducedRank_= 3; config.correl_alpha_ = 0.0 ; config.correl_beta_  = 0.1;
	//QuantLib::Array found_abcd = marketData_LMM_ABCD_calibration(config,pLmmSwaptionMarketData);

	std::vector<LMM::Index> exerciseDates;
	for(size_t i=1 ; i<=config.model_nbYear_; i++)
	{
		exerciseDates.push_back(2*i);
	}
	//exerciseDates.push_back(2);
	//exerciseDates.push_back(4);
	//exerciseDates.push_back(6);
	//exerciseDates.push_back(8);
	//exerciseDates.push_back(10);
	//exerciseDates.push_back(12);
	//exerciseDates.push_back(14);
	//exerciseDates.push_back(16);
	//exerciseDates.push_back(18);
	//exerciseDates.push_back(20);
	//exerciseDates.push_back(22);
	//exerciseDates.push_back(24);
	//exerciseDates.push_back(26);
	//exerciseDates.push_back(28);
	//exerciseDates.push_back(30);
	//exerciseDates.push_back(32);
	//exerciseDates.push_back(34);
	//exerciseDates.push_back(36);
	//exerciseDates.push_back(38);
	//exerciseDates.push_back(40);


	//std::vector<double> initLiborValues;
	//initLiborValues.push_back(0.00304663);
	//initLiborValues.push_back(0.00283432);
	//initLiborValues.push_back(0.00314012);
	//initLiborValues.push_back(0.0037196);
	//initLiborValues.push_back(0.0048919);
	//initLiborValues.push_back(0.00490389);
	//initLiborValues.push_back(0.00745992);
	//initLiborValues.push_back(0.00748785);
	//initLiborValues.push_back(0.0104202);
	//initLiborValues.push_back(0.0104748);
	//initLiborValues.push_back(0.0140121);
	//initLiborValues.push_back(0.0141109);
	//initLiborValues.push_back(0.0173241);
	//initLiborValues.push_back(0.0174755);
	//initLiborValues.push_back(0.0204022);
	//initLiborValues.push_back(0.0206124);
	//initLiborValues.push_back(0.0226241);
	//initLiborValues.push_back(0.022883);
	//initLiborValues.push_back(0.0243152);
	//initLiborValues.push_back(0.0246144);
	//initLiborValues.push_back(0.0253627);

	double strike = 0.0137;

	McLmm_PTR mcLmm_for_pricer = getMcLmmExample(lmmTenorStructure, initLiborValues, config);

	//double nominal = 1.0;
	//GenericSwap_CONSTPTR genericVanillaSwap = InstrumentFactory::createVanillaSwap(	strike, 
		//																				indexStart, 
		//																				indexEnd, 
		//																				floatingLegTenorType, 
		//																				fixedLegTenorType,
		//																				lmmTenorStructure,
		//																				nominal);
	//CallableInstrument_PTR callableGenericSwap(new CallableGenericSwap(genericVanillaSwap, exerciseDates));

	VanillaSwap_CONSTPTR vanillaSwap(new VanillaSwap(	strike, 
														indexStart, 
														indexEnd, 
														floatingLegTenorType, 
														fixedLegTenorType,
														lmmTenorStructure));

	CallableInstrument_PTR callableBermudanSwap(new BermudanVanillaSwaption(vanillaSwap,exerciseDates));

	////for ab, bc, ca
	//vector<std::vector<size_t>> set;
	//set.push_back(std::vector<size_t>());
	//set.back().push_back(1);
	//set.back().push_back(1);
	//set.back().push_back(0);
	//set.push_back(std::vector<size_t>());
	//set.back().push_back(0);
	//set.back().push_back(1);
	//set.back().push_back(1);
	//set.push_back(std::vector<size_t>());
	//set.back().push_back(1);
	//set.back().push_back(0);
	//set.back().push_back(1);
	//for (size_t i= 0 ; i<3; i++)
	//{
	//	//set.push_back(std::vector<size_t>(3,0));
	//	//set.back()[i]+=1;
	//	for (size_t j=0; j<=i; j++)
	//	{
	//		set.push_back(std::vector<size_t>(3,0));
	//		set.back()[i]+=1;
	//		set.back()[j]+=1;
	//	}
	//}

	std::vector<std::vector<std::vector<size_t>>> subset;
	////getAllSubsets(set, subset);

	size_t degre = 2;
	size_t counter = 0;

	for(size_t swaprate_Degre = 0; swaprate_Degre <= degre; swaprate_Degre++) 
	{
		for(size_t pay_off_swaprate_Degre = 0; pay_off_swaprate_Degre <= degre; pay_off_swaprate_Degre++) 
		{
			for(size_t liborRate_Degre = 0; liborRate_Degre <= degre; liborRate_Degre++) 
			{
				subset.push_back(std::vector<std::vector<size_t>>());
				for(size_t i = 0; i <= liborRate_Degre; i++) 
				{
					subset.back().push_back(std::vector<size_t>());
					subset.back().back().push_back(0);
					subset.back().back().push_back(0);
					subset.back().back().push_back(i);							
				}
				for(size_t j = 1; j <= pay_off_swaprate_Degre; j++) 
				{
					subset.back().push_back(std::vector<size_t>());
					subset.back().back().push_back(0);
					subset.back().back().push_back(j);
					subset.back().back().push_back(0);
				}	
				for(size_t k = 1; k <= swaprate_Degre; k++) 
				{
					subset.back().push_back(std::vector<size_t>());
					subset.back().back().push_back(k);
					subset.back().back().push_back(0);
					subset.back().back().push_back(0);
				}						
			}			
		}		
	}

	//for(size_t liborRate_Degre = 2; liborRate_Degre <= degre; liborRate_Degre++) 
	//{
	//		subset.push_back(choosed_basis_set);
	//		for(size_t i = 3; i <= liborRate_Degre; i++) 
	//		{
	//			subset.back().push_back(std::vector<size_t>());
	//			subset.back().back().push_back(0);
	//			subset.back().back().push_back(0);
	//			subset.back().back().push_back(i);							
	//		}
	//}

	std::vector<size_t> nbSimulation_vect;
	//nbSimulation_vect.push_back(20);
	nbSimulation_vect.push_back(config.nbSimulation_);

	std::vector<std::vector<double>> basis_value_on_allPath_buffer(nbSimulation_vect[0]);

	McLmm_PTR mcLmm = getMcLmmExample(lmmTenorStructure, initLiborValues, config);


	VanillaSwaption vanillaswaption(*vanillaSwap.get(), OptionType::OptionType::CALL);
	LmmVanillaSwaptionApproxPricer_Rebonato  LmmVanillaSwaptionApproxPricer_Rebonato(mcLmm->get_lmm());
	double black_vol =  LmmVanillaSwaptionApproxPricer_Rebonato.volBlack(vanillaswaption, initLiborValues);
	cout << "g : " << config.g << " beta: " << config.correl_beta_ <<endl;
	cout << "black_vol : " << black_vol << endl;

	McLmm_LS mcLmm_LS_backward(mcLmm);
	McLmm_LS mcLmm_LS_forward(mcLmm);

	std::stringstream outputFileName_test_time_s; 
	outputFileName_test_time_s << "Test_LS_pricing_with_parameter_time_test"<<".csv";
	std::string outputFileName_test_time = LMMPATH::get_Root_OutputPath() + outputFileName_test_time_s.str();
	ofstream out_test_time;
	out_test_time.open(outputFileName_test_time,  ios::out | ios::app );
	//out.open(outputFileName,  ios::out);
	out_test_time<<endl;
	out_test_time<<endl;
	out_test_time<<endl;


	mcLmm_LS_backward.simulateLMM(nbSimulation_vect[0]);
	const std::vector<McLmm_LS::LMMSimulationResult>&  lmmSimualtionResults_backward = mcLmm_LS_backward.lmmSimualtionResults_;
	mcLmm_LS_backward.write_to_stream(out_test_time);
	mcLmm_LS_forward.simulateLMM(nbSimulation_vect[0]);
	const std::vector<McLmm_LS::LMMSimulationResult>&  lmmSimualtionResults_forward = mcLmm_LS_forward.lmmSimualtionResults_;
	mcLmm_LS_forward.write_to_stream(out_test_time);

	std::stringstream outputFileName_s; 
	outputFileName_s<<"Test_LS_pricing_with_parameter"<<".csv";
	std::string outputFileName = LMMPATH::get_Root_OutputPath() + outputFileName_s.str();
	ofstream out;
	out.open(outputFileName,  ios::out | ios::app );
	//out.open(outputFileName,  ios::out);
	out<<endl;
	out<<endl;
	out<<endl;


	time_t _time;
	struct tm timeInfo;
	char format[32];
	time(&_time);
	localtime_s(&timeInfo, &_time);
	strftime(format, 32, "%Y-%m-%d %H-%M", &timeInfo);

	out << format << endl;
	out<< endl;
	out << "nbYear : " <<  config.model_nbYear_ << endl;
	out << "nbSimulation : " << config.nbSimulation_ << " ;" << endl;
	out<< "g : ;" << config.g << "; " << "correl_beta_ : ;" << config.correl_beta_   << " ;" <<endl;
	out << "black_vol : ; " << black_vol << endl;
	out<< "Combinaison ;" << "prix ;" << "inf IC ;" << "sup IC ;" << "demi longeur IC ;" << "temps ;" <<endl;



	for(size_t i = 0; i < subset.size(); i++) 
	{
		cout << "iteration : " << i << endl;		
		Test_LS_pricing_One_SubSet_basis(	subset[i], 
											strike,
											indexStart, 
											indexEnd,
											floatingLegTenorType,
											fixedLegTenorType,
											initLiborValues,
											exerciseDates,
											callableBermudanSwap,
											lmmSimualtionResults_backward,
											lmmSimualtionResults_forward,
											nbSimulation_vect,
											basis_value_on_allPath_buffer,
											out,
											out_test_time);		
	}

	//VanillaSwap vanillaSwap(	strike, 
	//							2, 
	//							indexEnd, 
	//							floatingLegTenorType, 
	//							fixedLegTenorType, 
	//							lmmTenorStructure);
	//LmmVanillaSwaptionApproxPricer_Rebonato  LmmVanillaSwaptionApproxPricer_Rebonato(mcLmm_for_pricer);
	//	double LmmVanillaSwaptionApproxPricer_Rebonato::volBlack(vanillaswaption, liborsInitValue)
}

void Test_10_basis()
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

	//config.correl_ReducedRank_= 3; config.correl_alpha_ = 0.0 ; config.correl_beta_  = 0.1;
	//QuantLib::Array found_abcd = marketData_LMM_ABCD_calibration(config,pLmmSwaptionMarketData);

	std::vector<LMM::Index> exerciseDates;
	exerciseDates.push_back(2);
	//exerciseDates.push_back(4);
	exerciseDates.push_back(6);
	//exerciseDates.push_back(8);
	exerciseDates.push_back(10);
	//exerciseDates.push_back(12);
	exerciseDates.push_back(14);
	//exerciseDates.push_back(16);
	exerciseDates.push_back(18);
	exerciseDates.push_back(20);


	//std::vector<double> initLiborValues;
	//initLiborValues.push_back(0.00304663);
	//initLiborValues.push_back(0.00283432);
	//initLiborValues.push_back(0.00314012);
	//initLiborValues.push_back(0.0037196);
	//initLiborValues.push_back(0.0048919);
	//initLiborValues.push_back(0.00490389);
	//initLiborValues.push_back(0.00745992);
	//initLiborValues.push_back(0.00748785);
	//initLiborValues.push_back(0.0104202);
	//initLiborValues.push_back(0.0104748);
	//initLiborValues.push_back(0.0140121);
	//initLiborValues.push_back(0.0141109);
	//initLiborValues.push_back(0.0173241);
	//initLiborValues.push_back(0.0174755);
	//initLiborValues.push_back(0.0204022);
	//initLiborValues.push_back(0.0206124);
	//initLiborValues.push_back(0.0226241);
	//initLiborValues.push_back(0.022883);
	//initLiborValues.push_back(0.0243152);
	//initLiborValues.push_back(0.0246144);
	//initLiborValues.push_back(0.0253627);

	double strike = 0.0137;

	McLmm_PTR mcLmm_for_pricer = getMcLmmExample(lmmTenorStructure, initLiborValues, LmmCalibrationConfig());



	VanillaSwap_CONSTPTR vanillaSwap(new VanillaSwap(	strike, 
														indexStart, 
														indexEnd, 
														floatingLegTenorType, 
														fixedLegTenorType,
														lmmTenorStructure));

	CallableInstrument_PTR callableBermudanSwap(new BermudanVanillaSwaption(vanillaSwap,exerciseDates));

	std::vector<std::vector<std::vector<size_t>>> subset;
	get_bisis_subset(subset);

	size_t degre = 2;
	size_t counter = 0;

	std::vector<size_t> nbSimulation_vect;
	//nbSimulation_vect.push_back(20);
	nbSimulation_vect.push_back(10000);

	std::vector<std::vector<double>> basis_value_on_allPath_buffer(nbSimulation_vect[0]);

	McLmm_PTR mcLmm = getMcLmmExample(lmmTenorStructure, initLiborValues, LmmCalibrationConfig());
	McLmm_LS mcLmm_LS(mcLmm);
	mcLmm_LS.simulateLMM(nbSimulation_vect[0]);
	std::vector<McLmm_LS::LMMSimulationResult>  lmmSimualtionResults_backward = mcLmm_LS.lmmSimualtionResults_;
	mcLmm_LS.simulateLMM(nbSimulation_vect[0]);
	std::vector<McLmm_LS::LMMSimulationResult>  lmmSimualtionResults_forward = mcLmm_LS.lmmSimualtionResults_;
	
	std::stringstream outputFileName_s; 
	outputFileName_s<<"Test_LS_price_10_of_basis"<<".csv";
	std::string outputFileName = LMMPATH::get_Root_OutputPath() + outputFileName_s.str();
	ofstream out;
	out.open(outputFileName,  ios::out | ios::app );
	//out.open(outputFileName,  ios::out);
	out<<endl;
	out<<endl;
	out<<endl;

	std::stringstream outputFileName_test_time_s; 
	outputFileName_test_time_s<<"Test_LS_test_10_subset_function_time"<<".csv";
	std::string outputFileName_test_time = LMMPATH::get_Root_OutputPath() + outputFileName_test_time_s.str();
	ofstream out_test_time;
	out_test_time.open(outputFileName_test_time,  ios::out | ios::app );
	//out.open(outputFileName,  ios::out);
	out_test_time<<endl;
	out_test_time<<endl;
	out_test_time<<endl;

	time_t _time;
	struct tm timeInfo;
	char format[32];
	time(&_time);
	localtime_s(&timeInfo, &_time);
	strftime(format, 32, "%Y-%m-%d %H-%M", &timeInfo);

	out << format << endl;
	out<< endl;

	for(size_t i = 0; i < subset.size(); i++) 
	{
			cout << "iteration : " << i << endl;
			Test_LS_pricing_One_SubSet_basis(	subset[i], 
												strike,
												indexStart, 
												indexEnd,
												floatingLegTenorType,
												fixedLegTenorType,
												initLiborValues,
												exerciseDates,
												callableBermudanSwap,
												lmmSimualtionResults_backward,
												lmmSimualtionResults_forward,
												nbSimulation_vect,
												basis_value_on_allPath_buffer,
												out,
												out_test_time);
		}
}

void get_bisis_subset(std::vector<std::vector<std::vector<size_t>>>& subset)
{
		subset.clear();
		subset.push_back(std::vector<std::vector<size_t>>());
		for(size_t i = 0; i <= 1; i++) 
		{
			subset.back().push_back(std::vector<size_t>());
			subset.back().back().push_back(i);
			subset.back().back().push_back(0);
			subset.back().back().push_back(0);							
		}
		for(size_t j = 1; j <= 2; j++) 
		{
			subset.back().push_back(std::vector<size_t>());
			subset.back().back().push_back(0);
			subset.back().back().push_back(j);
			subset.back().back().push_back(0);
		}	
		for(size_t k = 1; k <= 1; k++) 
		{
			subset.back().push_back(std::vector<size_t>());
			subset.back().back().push_back(0);
			subset.back().back().push_back(0);
			subset.back().back().push_back(k);
		}	


		subset.push_back(std::vector<std::vector<size_t>>());
		for(size_t i = 0; i <= 2; i++) 
		{
			subset.back().push_back(std::vector<size_t>());
			subset.back().back().push_back(i);
			subset.back().back().push_back(0);
			subset.back().back().push_back(0);							
		}
		for(size_t j = 1; j <= 1; j++) 
		{
			subset.back().push_back(std::vector<size_t>());
			subset.back().back().push_back(0);
			subset.back().back().push_back(j);
			subset.back().back().push_back(0);
		}	
		for(size_t k = 1; k <= 1; k++) 
		{
			subset.back().push_back(std::vector<size_t>());
			subset.back().back().push_back(0);
			subset.back().back().push_back(0);
			subset.back().back().push_back(k);
		}

		subset.push_back(std::vector<std::vector<size_t>>());
		for(size_t i = 0; i <= 2; i++) 
		{
			subset.back().push_back(std::vector<size_t>());
			subset.back().back().push_back(i);
			subset.back().back().push_back(0);
			subset.back().back().push_back(0);							
		}
		for(size_t j = 1; j <= 2; j++) 
		{
			subset.back().push_back(std::vector<size_t>());
			subset.back().back().push_back(0);
			subset.back().back().push_back(j);
			subset.back().back().push_back(0);
		}	
		for(size_t k = 1; k <= 1; k++) 
		{
			subset.back().push_back(std::vector<size_t>());
			subset.back().back().push_back(0);
			subset.back().back().push_back(0);
			subset.back().back().push_back(k);
		}	

		subset.push_back(std::vector<std::vector<size_t>>());
		for(size_t i = 0; i <= 0; i++) 
		{
			subset.back().push_back(std::vector<size_t>());
			subset.back().back().push_back(i);
			subset.back().back().push_back(0);
			subset.back().back().push_back(0);							
		}
		for(size_t j = 1; j <= 2; j++) 
		{
			subset.back().push_back(std::vector<size_t>());
			subset.back().back().push_back(0);
			subset.back().back().push_back(j);
			subset.back().back().push_back(0);
		}	
		for(size_t k = 1; k <= 1; k++) 
		{
			subset.back().push_back(std::vector<size_t>());
			subset.back().back().push_back(0);
			subset.back().back().push_back(0);
			subset.back().back().push_back(k);
		}	

		subset.push_back(std::vector<std::vector<size_t>>());
		for(size_t i = 0; i <= 1; i++) 
		{
			subset.back().push_back(std::vector<size_t>());
			subset.back().back().push_back(i);
			subset.back().back().push_back(0);
			subset.back().back().push_back(0);							
		}
		for(size_t j = 1; j <= 2; j++) 
		{
			subset.back().push_back(std::vector<size_t>());
			subset.back().back().push_back(0);
			subset.back().back().push_back(j);
			subset.back().back().push_back(0);
		}	
		for(size_t k = 1; k <= 2; k++) 
		{
			subset.back().push_back(std::vector<size_t>());
			subset.back().back().push_back(0);
			subset.back().back().push_back(0);
			subset.back().back().push_back(k);
		}

		subset.push_back(std::vector<std::vector<size_t>>());
		for(size_t i = 0; i <= 2; i++) 
		{
			subset.back().push_back(std::vector<size_t>());
			subset.back().back().push_back(i);
			subset.back().back().push_back(0);
			subset.back().back().push_back(0);							
		}
		for(size_t j = 1; j <= 1; j++) 
		{
			subset.back().push_back(std::vector<size_t>());
			subset.back().back().push_back(0);
			subset.back().back().push_back(j);
			subset.back().back().push_back(0);
		}	
		for(size_t k = 1; k <= 2; k++) 
		{
			subset.back().push_back(std::vector<size_t>());
			subset.back().back().push_back(0);
			subset.back().back().push_back(0);
			subset.back().back().push_back(k);
		}

		subset.push_back(std::vector<std::vector<size_t>>());
		for(size_t i = 0; i <= 2; i++) 
		{
			subset.back().push_back(std::vector<size_t>());
			subset.back().back().push_back(i);
			subset.back().back().push_back(0);
			subset.back().back().push_back(0);							
		}
		for(size_t j = 1; j <= 0; j++) 
		{
			subset.back().push_back(std::vector<size_t>());
			subset.back().back().push_back(0);
			subset.back().back().push_back(j);
			subset.back().back().push_back(0);
		}	
		for(size_t k = 1; k <= 1; k++) 
		{
			subset.back().push_back(std::vector<size_t>());
			subset.back().back().push_back(0);
			subset.back().back().push_back(0);
			subset.back().back().push_back(k);
		}	

		subset.push_back(std::vector<std::vector<size_t>>());
		for(size_t i = 0; i <= 2; i++) 
		{
			subset.back().push_back(std::vector<size_t>());
			subset.back().back().push_back(i);
			subset.back().back().push_back(0);
			subset.back().back().push_back(0);							
		}
		for(size_t j = 1; j <= 2; j++) 
		{
			subset.back().push_back(std::vector<size_t>());
			subset.back().back().push_back(0);
			subset.back().back().push_back(j);
			subset.back().back().push_back(0);
		}	
		for(size_t k = 1; k <= 2; k++) 
		{
			subset.back().push_back(std::vector<size_t>());
			subset.back().back().push_back(0);
			subset.back().back().push_back(0);
			subset.back().back().push_back(k);
		}	

		subset.push_back(std::vector<std::vector<size_t>>());
		for(size_t i = 0; i <= 1; i++) 
		{
			subset.back().push_back(std::vector<size_t>());
			subset.back().back().push_back(i);
			subset.back().back().push_back(0);
			subset.back().back().push_back(0);							
		}
		for(size_t j = 1; j <= 0; j++) 
		{
			subset.back().push_back(std::vector<size_t>());
			subset.back().back().push_back(0);
			subset.back().back().push_back(j);
			subset.back().back().push_back(0);
		}	
		for(size_t k = 1; k <= 2; k++) 
		{
			subset.back().push_back(std::vector<size_t>());
			subset.back().back().push_back(0);
			subset.back().back().push_back(0);
			subset.back().back().push_back(k);
		}	

		subset.push_back(std::vector<std::vector<size_t>>());
		for(size_t i = 0; i <= 0; i++) 
		{
			subset.back().push_back(std::vector<size_t>());
			subset.back().back().push_back(i);
			subset.back().back().push_back(0);
			subset.back().back().push_back(0);							
		}
		for(size_t j = 1; j <= 2; j++) 
		{
			subset.back().push_back(std::vector<size_t>());
			subset.back().back().push_back(0);
			subset.back().back().push_back(j);
			subset.back().back().push_back(0);
		}	
		for(size_t k = 1; k <= 2; k++) 
		{
			subset.back().push_back(std::vector<size_t>());
			subset.back().back().push_back(0);
			subset.back().back().push_back(0);
			subset.back().back().push_back(k);
		}	
}

