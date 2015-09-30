#include <LMM/Test/Test_Junbin/JBTests.h>
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
#include <Instrument/InstrumentFactory.h>
#include <LMM/Pricer/McLmmPricer/McLmmPricer.h>
#include <LMM/Pricer/McLmmPricer/McLmmGenericSwapPricer.h>
#include <LMM/Pricer/Longstaff_Schwartz/LS_BackwardAlgo.h>
#include <LMM/Pricer/Longstaff_Schwartz/LS_ForwardAlgo.h>


// It is an example of McLmm
McLmm_PTR getMcLmmExample(LMMTenorStructure_PTR lmmTenorStructure, const std::vector<double>& initLiborValues, const LmmCalibrationConfig& config)
{
	////! Parameter of h
	double a = 0.105;
	double b = 0.317;
	double c = 0.515;
	double d = 0.356;
	Shifted_HGVolatilityParam::ABCDParameter abcdParam (a,b,c,d);
	//Parameter of hg
	double g_constParam = config.g;
	double shift_constParam = 0.0;
	ConstShifted_HGVolatilityParam_PTR hgParam( new ConstShifted_HGVolatilityParam(lmmTenorStructure,abcdParam,g_constParam,shift_constParam));
	
	//! Correlation 1
	size_t nbFactor       = 3; // need to test nbFactor  = 3, and nbFactor = 
	size_t correlFullRank = lmmTenorStructure->get_horizon()+1;
	size_t correlReducedRank = nbFactor;
	//!"Check Parameters(): Condition not implemented yet."
	std::cout << "checkParams(): ";
	CorrelationReductionType::CorrelationReductionType correlReductionType = CorrelationReductionType::PCA;
	double correlAlpha = config.correl_alpha_;
	double correlBeta  = config.correl_beta_;
	Correlation_PTR correlation(new XY_beta_Correlation(correlFullRank,correlReducedRank, correlReductionType,correlAlpha,correlBeta));
	correlation->calculate(); // for print.
	correlation->print("test_McTerminalLmm_Correlation.csv");
	//hgVolatilityFunction
	ConstShifted_HGVolatilityFunction_PTR hgVolatilityFunction (new ConstShifted_HGVolatilityFunction(lmmTenorStructure, correlation, hgParam)); 
	hgVolatilityFunction->print("test_McTerminalLmm_Volatility.csv");
	//! Dispersion
	Dispersion dispersion(hgVolatilityFunction);

	unsigned long seed = 5033;
	RNGenerator_PTR  rnGenerator(new McGenerator(seed));

	//build lmm and mcLmm model
	Lmm_PTR shiftedLmm (new Lmm(dispersion));

	McLmm_PTR mcLmm(new McTerminalLmm(shiftedLmm, initLiborValues, rnGenerator, MCSchemeType::EULER));
	return mcLmm;
}

void getAllSubsets(const vector<std::vector<size_t>>& set, std::vector<std::vector<std::vector<size_t>>>& subset)
{
    std::vector<std::vector<size_t>> empty;
    subset.push_back( empty );

    for (size_t i = 0; i < set.size(); i++)
    {
        std::vector<std::vector<std::vector<size_t>>> subsetTemp = subset;

        for (size_t j = 0; j < subsetTemp.size(); j++)
            subsetTemp[j].push_back( set[i] );

        for (size_t j = 0; j < subsetTemp.size(); j++)
            subset.push_back( subsetTemp[j] );
    }
}

Basis_CONSTPTR getBasis(const std::vector<size_t>& basis_vect ,
						double coeff, 
						const VanillaSwap& vanillaSwap,
						double strike, 
						LMM::Index liborIndex)
{
	assert(basis_vect.size()==3);
	size_t swaprate_degree			= basis_vect[0];
	size_t payoff_swaprate_degree	= basis_vect[1];
	size_t liborRate_degree			= basis_vect[2];
	if(swaprate_degree==0&&payoff_swaprate_degree==0&&liborRate_degree==0)
		return Basis_CONSTPTR(new Basis_ConstUnity());
	if(swaprate_degree>0&&payoff_swaprate_degree==0&&liborRate_degree==0)
	{
		return 	getSingleEVBasis_swaprate(vanillaSwap, swaprate_degree);
	}
	if(swaprate_degree==0&&payoff_swaprate_degree>0&&liborRate_degree==0)
	{
		std::vector<Basis_CONSTPTR> basis_vect;
		basis_vect.push_back(getSingleEVBasis_swaprateCall(vanillaSwap, strike, liborIndex));	
		return Basis_CONSTPTR(new Basis_Polynomial(basis_vect, coeff, std::vector<size_t>(1,payoff_swaprate_degree)));
	}

	if(swaprate_degree==0&&payoff_swaprate_degree==0&&liborRate_degree>0)
		return getSingleEVBasis_Libor(liborIndex,liborRate_degree);

	if(swaprate_degree==1&&payoff_swaprate_degree==1&&liborRate_degree==0)
	{
		std::vector<Basis_CONSTPTR> basis_vect;
		basis_vect.push_back(getSingleEVBasis_swaprate(vanillaSwap, swaprate_degree));
		basis_vect.push_back(getSingleEVBasis_swaprateCall(vanillaSwap, strike, liborIndex));	
		return Basis_CONSTPTR(new Basis_Polynomial(basis_vect, coeff, std::vector<size_t>(2,swaprate_degree)));
	}

	if(swaprate_degree==0&&payoff_swaprate_degree==1&&liborRate_degree==1)
	{
		std::vector<Basis_CONSTPTR> basis_vect;
		basis_vect.push_back(getSingleEVBasis_swaprateCall(vanillaSwap, strike, liborIndex));	
		basis_vect.push_back(getSingleEVBasis_Libor(liborIndex,liborRate_degree));
		return Basis_CONSTPTR(new Basis_Polynomial(basis_vect, coeff, std::vector<size_t>(2,payoff_swaprate_degree)));
	}
	if(swaprate_degree==1&&payoff_swaprate_degree==0&&liborRate_degree==1)
	{
		std::vector<Basis_CONSTPTR> basis_vect;
		basis_vect.push_back(getSingleEVBasis_swaprate(vanillaSwap, swaprate_degree));	
		basis_vect.push_back(getSingleEVBasis_Libor(liborIndex,liborRate_degree));
		return Basis_CONSTPTR(new Basis_Polynomial(basis_vect, coeff, std::vector<size_t>(2,swaprate_degree)));
	}

	throw("basis_vect mismatch the product");

	//basis_vect.push_back(getSingleEVBasis_swaprate(vanillaSwap));
	//basis_vect.push_back(getSingleEVBasis_swaprateCall(vanillaSwap, strike, liborIndex));
	//basis_vect.push_back(getSingleEVBasis_Libor(liborIndex));

	//return Basis_CONSTPTR(new Basis_Polynomial(basis_vect, coeff, power_vect));
}

Basis_Evaluator_CONSTPTR getBasisEvaluator(const std::vector<size_t>& basis_evaluator_vect ,
										   const McLmmVanillaSwapPricer& mcLmmVanillaSwapPricer)
{

	assert(basis_evaluator_vect.size()==3);
	size_t swaprate_degree			= basis_evaluator_vect[0];
	size_t payoff_swaprate_degree	= basis_evaluator_vect[1];
	size_t liborRate_degree			= basis_evaluator_vect[2];
	if(swaprate_degree==0&&payoff_swaprate_degree==0&&liborRate_degree==0)
		return Basis_Evaluator_CONSTPTR(new Basis_SinglaEVFunctional_Evaluator(
										EV_Evaluator_CONSTPTR(new EV_ConstRate_Evaluator())));
	if(swaprate_degree>0&&payoff_swaprate_degree==0&&liborRate_degree==0)
		return 	Basis_Evaluator_CONSTPTR(	new Basis_SinglaEVFunctional_Evaluator(
											EV_Evaluator_CONSTPTR(new EV_VanillaSwapRate_Evaluator(McLmmVanillaSwapPricer(mcLmmVanillaSwapPricer)))));
	if(swaprate_degree==0&&payoff_swaprate_degree>0&&liborRate_degree==0)
	{
		std::vector<Basis_Evaluator_CONSTPTR> basis_evaluator_vect;
	    basis_evaluator_vect.push_back(	Basis_Evaluator_CONSTPTR(new Basis_SinglaEVFunctional_Evaluator(
										EV_Evaluator_CONSTPTR(new EV_VanillaSwapRate_Evaluator(McLmmVanillaSwapPricer(mcLmmVanillaSwapPricer))))));
		return Basis_Evaluator_CONSTPTR(new Basis_Polynomial_Evaluator(basis_evaluator_vect));
	}


	if(swaprate_degree==0&&payoff_swaprate_degree==0&&liborRate_degree>0)
		return Basis_Evaluator_CONSTPTR(new Basis_SinglaEVFunctional_Evaluator(
										EV_Evaluator_CONSTPTR(new EV_LiborRate_Evaluator())));

	if(swaprate_degree==1&&payoff_swaprate_degree==1&&liborRate_degree==0)
	{
		std::vector<Basis_Evaluator_CONSTPTR> basis_evaluator_vect;
	    basis_evaluator_vect.push_back(	Basis_Evaluator_CONSTPTR(new Basis_SinglaEVFunctional_Evaluator(
										EV_Evaluator_CONSTPTR(new EV_VanillaSwapRate_Evaluator(McLmmVanillaSwapPricer(mcLmmVanillaSwapPricer))))));
		basis_evaluator_vect.push_back(	Basis_Evaluator_CONSTPTR(new Basis_SinglaEVFunctional_Evaluator(
										EV_Evaluator_CONSTPTR(new EV_VanillaSwapRate_Evaluator(McLmmVanillaSwapPricer(mcLmmVanillaSwapPricer))))));
		return Basis_Evaluator_CONSTPTR(new Basis_Polynomial_Evaluator(basis_evaluator_vect));
	}

	if(swaprate_degree==0&&payoff_swaprate_degree==1&&liborRate_degree==1)
	{
		std::vector<Basis_Evaluator_CONSTPTR> basis_evaluator_vect;
	    basis_evaluator_vect.push_back(	Basis_Evaluator_CONSTPTR(new Basis_SinglaEVFunctional_Evaluator(
										EV_Evaluator_CONSTPTR(new EV_VanillaSwapRate_Evaluator(McLmmVanillaSwapPricer(mcLmmVanillaSwapPricer))))));
		basis_evaluator_vect.push_back(	Basis_Evaluator_CONSTPTR(new Basis_SinglaEVFunctional_Evaluator(
										EV_Evaluator_CONSTPTR(new EV_LiborRate_Evaluator()))));
		return Basis_Evaluator_CONSTPTR(new Basis_Polynomial_Evaluator(basis_evaluator_vect));
	}
	if(swaprate_degree==1&&payoff_swaprate_degree==0&&liborRate_degree==1)
	{
		std::vector<Basis_Evaluator_CONSTPTR> basis_evaluator_vect;
	    basis_evaluator_vect.push_back(	Basis_Evaluator_CONSTPTR(new Basis_SinglaEVFunctional_Evaluator(
										EV_Evaluator_CONSTPTR(new EV_VanillaSwapRate_Evaluator(McLmmVanillaSwapPricer(mcLmmVanillaSwapPricer))))));
		basis_evaluator_vect.push_back(	Basis_Evaluator_CONSTPTR(new Basis_SinglaEVFunctional_Evaluator(
										EV_Evaluator_CONSTPTR(new EV_LiborRate_Evaluator()))));
		return Basis_Evaluator_CONSTPTR(new Basis_Polynomial_Evaluator(basis_evaluator_vect));
	}

	throw("basis_evaluator_vect mismatch the basis");

	//std::vector<Basis_Evaluator_CONSTPTR> basis_evaluator_vect;
	//    basis_evaluator_vect.push_back(	Basis_Evaluator_CONSTPTR(new Basis_SinglaEVFunctional_Evaluator(
	//									EV_Evaluator_CONSTPTR(new EV_VanillaSwapRate_Evaluator(McLmmVanillaSwapPricer(mcLmmVanillaSwapPricer))))));
	//    basis_evaluator_vect.push_back(	Basis_Evaluator_CONSTPTR(new Basis_SinglaEVFunctional_Evaluator(
	//									EV_Evaluator_CONSTPTR(new EV_VanillaSwapRate_Evaluator(McLmmVanillaSwapPricer(mcLmmVanillaSwapPricer))))));
	//	basis_evaluator_vect.push_back(	Basis_Evaluator_CONSTPTR(new Basis_SinglaEVFunctional_Evaluator(
	//									EV_Evaluator_CONSTPTR(new EV_LiborRate_Evaluator()))));	


	//return Basis_Evaluator_CONSTPTR(new Basis_Polynomial_Evaluator(basis_evaluator_vect));
}


Basis_CONSTPTR getSingleEVBasis_Libor(LMM::Index liborIndex,  size_t liborRate_degree)
{
	return Basis_CONSTPTR(new Basis_Monomial(	EV_CONSTPTR(new EV_LiborRate(Rate1_CONSTPTR(new LiborRate(liborIndex,Tenor(Tenor::_6M))))), 
												Basis_SinglaEVFunctional::Transformer_CONSTPTR(new Basis_Monomial::MonomialTransformer(1.0,liborRate_degree))));
}

Basis_CONSTPTR getSingleEVBasis_swaprate(const VanillaSwap& vanillaSwap, size_t swaprate_degree)
{
	return Basis_CONSTPTR(new Basis_Monomial(				
			EV_CONSTPTR(new EV_VanillaSwapRate( Rate1_CONSTPTR( new VanillaSwapRate(
							VanillaSwap(vanillaSwap))))), 
			Basis_SinglaEVFunctional::Transformer_CONSTPTR(new Basis_Monomial::MonomialTransformer(1.0,swaprate_degree))));
}

Basis_CONSTPTR getSingleEVBasis_swaprateCall(const VanillaSwap& vanillaSwap, double strike, LMM::Index liborIndex)
{
	return 		Basis_CONSTPTR(new Basis_CappedFlooredCoupon(	
						EV_CONSTPTR(new EV_VanillaSwapRate( Rate1_CONSTPTR( new VanillaSwapRate(VanillaSwap(vanillaSwap))))), 
						Basis_SinglaEVFunctional::Transformer_CONSTPTR(
								new Basis_CappedFlooredCoupon::CappedFlooredCouponTransformer(
									CappedFlooredCoupon_CONSTPTR(new CappedFlooredCoupon(CappedFlooredCoupon(	
																												liborIndex,
																												1.0,
																												1.0, 
																												true,
																												0.0,
																												false, 
																												10e100, 
																												Rate1_CONSTPTR(new Rate1()),
																												10.0,
																												-10.0*strike, 
																												liborIndex)))))));
}


//void Test_LS_pricing_allSubSet_basis()
//{
//	LMM::Index  indexStart = 2;		//1Y
//	LMM::Index  indexEnd   = 20;		//10Y
//	Tenor	floatingLegTenorType = Tenor::_6M;
//	Tenor	fixedLegTenorType    = Tenor::_1YR;
//	assert(indexStart%2==0&&indexEnd%2==0);
//	LMMTenorStructure_PTR lmmTenorStructure( new LMMTenorStructure(floatingLegTenorType, indexEnd/2));
//
//	std::vector<std::string> mkt_file_list = InputFileManager::get_VCUB_FileList();
//	const std::string& mkt_data_file = mkt_file_list.back();
//	std::string folder_name;   // = "TotalCalib\\" ;  config.use_positive_constraint_=true;
//	std::string base_name_file = LMMPATH::get_BaseFileName(mkt_data_file) + "\\";
//	folder_name+=base_name_file;
//	LMMPATH::reset_Output_SubFolder(folder_name );
//
//	LmmCalibrationConfig config;
//
//	config.floatLegTenor_=floatingLegTenorType;
//	config.fixedLegTenor_=fixedLegTenorType;
//
//	config.model_nbYear_		=	indexEnd/2;
//	size_t fixedFloatRatio		=	config.fixedLegTenor_.ratioTo(config.floatLegTenor_);
//	config.correl_FullRank_		=	fixedFloatRatio*config.model_nbYear_+1;
//
//
//	LmmSwaptionMarketData_PTR pLmmSwaptionMarketData	=	get_LmmSwaptionMarketData(config, mkt_data_file);
//	const std::vector<double>&	initLiborValues			=	pLmmSwaptionMarketData->get_LiborQuotes()->get_InitLibor();
//
//	//config.correl_ReducedRank_= 3; config.correl_alpha_ = 0.0 ; config.correl_beta_  = 0.1;
//	//QuantLib::Array found_abcd = marketData_LMM_ABCD_calibration(config,pLmmSwaptionMarketData);
//
//	std::vector<LMM::Index> exerciseDates;
//	exerciseDates.push_back(2);
//	exerciseDates.push_back(4);
//	exerciseDates.push_back(6);
//	exerciseDates.push_back(8);
//	exerciseDates.push_back(10);
//	exerciseDates.push_back(12);
//	exerciseDates.push_back(14);
//	exerciseDates.push_back(16);
//	exerciseDates.push_back(18);
//	exerciseDates.push_back(20);
//
//
//	//std::vector<double> initLiborValues;
//	//initLiborValues.push_back(0.00304663);
//	//initLiborValues.push_back(0.00283432);
//	//initLiborValues.push_back(0.00314012);
//	//initLiborValues.push_back(0.0037196);
//	//initLiborValues.push_back(0.0048919);
//	//initLiborValues.push_back(0.00490389);
//	//initLiborValues.push_back(0.00745992);
//	//initLiborValues.push_back(0.00748785);
//	//initLiborValues.push_back(0.0104202);
//	//initLiborValues.push_back(0.0104748);
//	//initLiborValues.push_back(0.0140121);
//	//initLiborValues.push_back(0.0141109);
//	//initLiborValues.push_back(0.0173241);
//	//initLiborValues.push_back(0.0174755);
//	//initLiborValues.push_back(0.0204022);
//	//initLiborValues.push_back(0.0206124);
//	//initLiborValues.push_back(0.0226241);
//	//initLiborValues.push_back(0.022883);
//	//initLiborValues.push_back(0.0243152);
//	//initLiborValues.push_back(0.0246144);
//	//initLiborValues.push_back(0.0253627);
//
//	double strike = 0.0137;
//
//	McLmm_PTR mcLmm_for_pricer = getMcLmmExample(lmmTenorStructure, initLiborValues,config);
//
//	//double nominal = 1.0;
//	//GenericSwap_CONSTPTR genericVanillaSwap = InstrumentFactory::createVanillaSwap(	strike, 
//	//																				indexStart, 
//	//																				indexEnd, 
//	//																				floatingLegTenorType, 
//	//																				fixedLegTenorType,
//	//																				lmmTenorStructure,
//	//																				nominal);
//	//CallableInstrument_PTR callableGenericSwap(new CallableGenericSwap(genericVanillaSwap, exerciseDates));
//
//
//	VanillaSwap_CONSTPTR vanillaSwap(new VanillaSwap(	strike, 
//														indexStart, 
//														indexEnd, 
//														floatingLegTenorType, 
//														fixedLegTenorType,
//														lmmTenorStructure));
//
//	CallableInstrument_PTR callableBermudanSwap(new BermudanVanillaSwaption(vanillaSwap,exerciseDates));
//
//	////for ab, bc, ca
//	//vector<std::vector<size_t>> set;
//	//set.push_back(std::vector<size_t>());
//	//set.back().push_back(1);
//	//set.back().push_back(1);
//	//set.back().push_back(0);
//	//set.push_back(std::vector<size_t>());
//	//set.back().push_back(0);
//	//set.back().push_back(1);
//	//set.back().push_back(1);
//	//set.push_back(std::vector<size_t>());
//	//set.back().push_back(1);
//	//set.back().push_back(0);
//	//set.back().push_back(1);
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
//	////getAllSubsets(set, subset);
//
//	size_t degre = 2;
//	size_t counter = 0;
//
//	//for 1 + a + a^2 + b + c 
//	//std::vector<std::vector<size_t>> choosed_basis_set;
//	//for(size_t i = 0; i <= 1; i++) 
//	//{
//	//	choosed_basis_set.push_back(std::vector<size_t>());
//	//	choosed_basis_set.back().push_back(0);
//	//	choosed_basis_set.back().push_back(0);
//	//	choosed_basis_set.back().push_back(i);							
//	//}
//	//for(size_t j = 1; j <= 1; j++) 
//	//{
//	//	choosed_basis_set.push_back(std::vector<size_t>());
//	//	choosed_basis_set.back().push_back(0);
//	//	choosed_basis_set.back().push_back(j);
//	//	choosed_basis_set.back().push_back(0);
//	//}	
//	//for(size_t k = 1; k <= 2; k++) 
//	//{
//	//	choosed_basis_set.push_back(std::vector<size_t>());
//	//	choosed_basis_set.back().push_back(k);
//	//	choosed_basis_set.back().push_back(0);
//	//	choosed_basis_set.back().push_back(0);
//	//}
//
//	//for 1 + a + a^2 + b + c + ab + bc +ac
//	//choosed_basis_set.push_back(std::vector<size_t>());
//	//choosed_basis_set.back().push_back(1);
//	//choosed_basis_set.back().push_back(1);
//	//choosed_basis_set.back().push_back(0);
//	//choosed_basis_set.push_back(std::vector<size_t>());
//	//choosed_basis_set.back().push_back(0);
//	//choosed_basis_set.back().push_back(1);
//	//choosed_basis_set.back().push_back(1);
//	//choosed_basis_set.push_back(std::vector<size_t>());
//	//choosed_basis_set.back().push_back(1);
//	//choosed_basis_set.back().push_back(0);
//	//choosed_basis_set.back().push_back(1);
//
//	//choosed_basis_set.push_back(std::vector<size_t>());
//	//choosed_basis_set.back().push_back(0);
//	//choosed_basis_set.back().push_back(3);
//	//choosed_basis_set.back().push_back(0);
//	//choosed_basis_set.push_back(std::vector<size_t>());
//	//choosed_basis_set.back().push_back(0);
//	//choosed_basis_set.back().push_back(4);
//	//choosed_basis_set.back().push_back(0);
//	//choosed_basis_set.push_back(std::vector<size_t>());
//	//choosed_basis_set.back().push_back(3);
//	//choosed_basis_set.back().push_back(0);
//	//choosed_basis_set.back().push_back(0);
//	//choosed_basis_set.push_back(std::vector<size_t>());
//	//choosed_basis_set.back().push_back(4);
//	//choosed_basis_set.back().push_back(0);
//	//choosed_basis_set.back().push_back(0);
//
//	//for (size_t i = 0; i<subset.size(); i++)
//	//{
//	//	subset[i].insert(subset[i].begin(), choosed_basis_set.begin(),choosed_basis_set.end());
//
//	//}
//
//	for(size_t swaprate_Degre = 0; swaprate_Degre <= degre; swaprate_Degre++) 
//	{
//		for(size_t pay_off_swaprate_Degre = 0; pay_off_swaprate_Degre <= degre; pay_off_swaprate_Degre++) 
//		{
//			for(size_t liborRate_Degre = 0; liborRate_Degre <= degre; liborRate_Degre++) 
//			{
//					subset.push_back(std::vector<std::vector<size_t>>());
//					for(size_t i = 0; i <= liborRate_Degre; i++) 
//					{
//						subset.back().push_back(std::vector<size_t>());
//						subset.back().back().push_back(0);
//						subset.back().back().push_back(0);
//						subset.back().back().push_back(i);							
//					}
//					for(size_t j = 1; j <= pay_off_swaprate_Degre; j++) 
//					{
//						subset.back().push_back(std::vector<size_t>());
//						subset.back().back().push_back(0);
//						subset.back().back().push_back(j);
//						subset.back().back().push_back(0);
//					}	
//					for(size_t k = 1; k <= swaprate_Degre; k++) 
//					{
//						subset.back().push_back(std::vector<size_t>());
//						subset.back().back().push_back(k);
//						subset.back().back().push_back(0);
//						subset.back().back().push_back(0);
//					}						
//			}			
//		}		
//	}
//
//	//for(size_t liborRate_Degre = 2; liborRate_Degre <= degre; liborRate_Degre++) 
//	//{
//	//		subset.push_back(choosed_basis_set);
//	//		for(size_t i = 3; i <= liborRate_Degre; i++) 
//	//		{
//	//			subset.back().push_back(std::vector<size_t>());
//	//			subset.back().back().push_back(0);
//	//			subset.back().back().push_back(0);
//	//			subset.back().back().push_back(i);							
//	//		}
//	//}
//
//
//
//
//	std::vector<size_t> nbSimulation_vect;
//	//nbSimulation_vect.push_back(20);
//	nbSimulation_vect.push_back(10000);
//
//	std::vector<std::vector<double>> basis_value_on_allPath_buffer(nbSimulation_vect[0]);
//
//	McLmm_PTR mcLmm = getMcLmmExample(lmmTenorStructure, initLiborValues, LmmCalibrationConfig());
//
//	McLmm_LS mcLmm_LS_backward(mcLmm);
//	McLmm_LS mcLmm_LS_forward(mcLmm);
//
//	mcLmm_LS_backward.simulateLMM(nbSimulation_vect[0]);
//	const std::vector<McLmm_LS::LMMSimulationResult>&  lmmSimualtionResults_backward = mcLmm_LS_backward.lmmSimualtionResults_;
//	mcLmm_LS_forward.simulateLMM(nbSimulation_vect[0]);
//	const std::vector<McLmm_LS::LMMSimulationResult>&  lmmSimualtionResults_forward = mcLmm_LS_forward.lmmSimualtionResults_;
//	
//	std::stringstream outputFileName_s; 
//	outputFileName_s<<"Test_LS_price_allSubSet_of_basis"<<".csv";
//	std::string outputFileName = LMMPATH::get_Root_OutputPath() + outputFileName_s.str();
//	ofstream out;
//	out.open(outputFileName,  ios::out | ios::app );
//	//out.open(outputFileName,  ios::out);
//	out<<endl;
//	out<<endl;
//	out<<endl;
//
//	std::stringstream outputFileName_test_time_s; 
//	outputFileName_test_time_s<<"Test_LS_test_function_time"<<".csv";
//	std::string outputFileName_test_time = LMMPATH::get_Root_OutputPath() + outputFileName_test_time_s.str();
//	ofstream out_test_time;
//	out_test_time.open(outputFileName_test_time,  ios::out | ios::app );
//	//out.open(outputFileName,  ios::out);
//	out_test_time<<endl;
//	out_test_time<<endl;
//	out_test_time<<endl;
//
//	time_t _time;
//	struct tm timeInfo;
//	char format[32];
//	time(&_time);
//	localtime_s(&timeInfo, &_time);
//	strftime(format, 32, "%Y-%m-%d %H-%M", &timeInfo);
//
//	out << format << endl;
//	out<< endl;
//	out<< "g : ;" << config.g << "; " << "correl_beta_ : ;" << config.correl_beta_   << " ;" <<endl;
//	out<< "Combinaison ;" << "prix ;" << "inf IC ;" << "sup IC ;" << "demi longeur IC ;" << "temps ;" <<endl;
//
//	for(size_t i = 0; i < subset.size(); i++) 
//	{
//			cout << "iteration : " << i << endl;
//			Test_LS_pricing_One_SubSet_basis(	subset[i], 
//												strike,
//												indexStart, 
//												indexEnd,
//												floatingLegTenorType,
//												fixedLegTenorType,
//												initLiborValues,
//												exerciseDates,
//												callableBermudanSwap,
//												lmmSimualtionResults_backward,
//												lmmSimualtionResults_forward,
//												nbSimulation_vect,
//												basis_value_on_allPath_buffer,
//												out,
//												out_test_time);
//	}
//
//	//VanillaSwap vanillaSwap(	strike, 
//	//							2, 
//	//							indexEnd, 
//	//							floatingLegTenorType, 
//	//							fixedLegTenorType, 
//	//							lmmTenorStructure);
//	//LmmVanillaSwaptionApproxPricer_Rebonato  LmmVanillaSwaptionApproxPricer_Rebonato(mcLmm_for_pricer);
//	//	double LmmVanillaSwaptionApproxPricer_Rebonato::volBlack(vanillaswaption, liborsInitValue)
//}

//void Test_Longstaff_Schwartz()
//{
//	LMM::Index  indexStart = 2;		//1Y
//	LMM::Index  indexEnd   = 6;		//4Y
//	Tenor	floatingLegTenorType = Tenor::_6M;
//	Tenor	fixedLegTenorType    = Tenor::_1YR;
//	assert(indexStart%2==0&&indexEnd%2==0);
//	LMMTenorStructure_PTR lmmTenorStructure( new LMMTenorStructure(floatingLegTenorType, indexEnd/2));	
//
//	std::vector<LMM::Index> exerciseDates;
//	exerciseDates.push_back(2);
//	exerciseDates.push_back(4);
//	exerciseDates.push_back(6);
//
//	double fwdRate = 0.02;
//	std::vector<double> initLiborValues(lmmTenorStructure->get_horizon()+1, fwdRate);
//	McLmm_PTR mcLmm = getMcLmmExample(lmmTenorStructure, initLiborValues);
//	double strike = fwdRate;
//	VanillaSwap_PTR vanillaSwap_PTR(new VanillaSwap(strike, 
//													indexStart, 
//													indexEnd, 
//													floatingLegTenorType, 
//													fixedLegTenorType, 
//													lmmTenorStructure));
//
//	CallableInstrument_CONSTPTR callableVanillaSwaption_PTR(new BermudanVanillaSwaption(vanillaSwap_PTR,exerciseDates));
//
//
//	// ----------------------------------------------------------------------
//	// 
//	//                     LS_Backward
//	//
//	// ---------------------------------------------------------------------
//	
//    McLmm_LS mcLmm_LS(mcLmm);
//
//	// ----------------------------------------------------------------------
//	// 
//	//                     LS_Backward
//	//
//	// ----------------------------------------------------------------------
//	LS_BackwardAlgo ls_BackwardAlgo(McLmmPricer_CONSTPTR(new McLmmVanillaSwapPricer(mcLmm)));
//
//	std::vector<LS::Regression>& regression_vect = ls_BackwardAlgo.getRegressions();
//	ls_BackwardAlgo.getExerciseDates()=exerciseDates;
//	
//	size_t fixedFloatingRatio = fixedLegTenorType.ratioTo(floatingLegTenorType);
//	for(size_t regressionIndex=0; regressionIndex<2; regressionIndex++)
//	{
//		LMM::Index liborIndex = indexStart + regressionIndex*fixedFloatingRatio;
//		std::vector<Basis_CONSTPTR> basis_vect;
//		basis_vect.push_back(Basis_CONSTPTR(new Basis_Monomial(	EV_CONSTPTR(new EV_LiborRate(Rate1_CONSTPTR(new LiborRate(liborIndex,Tenor(floatingLegTenorType))))), 
//																Basis_SinglaEVFunctional::Transformer_CONSTPTR(new Basis_Monomial::MonomialTransformer(1.0,1)))));
//	
//		std::vector<Basis_Evaluator_CONSTPTR> basis_evaluator_vect;
//		basis_evaluator_vect.push_back(	Basis_Evaluator_CONSTPTR(new Basis_SinglaEVFunctional_Evaluator(
//										EV_Evaluator_CONSTPTR(new EV_LiborRate_Evaluator()))));	
//		
//		regression_vect.push_back(LS::Regression(LS::RegressionRepresentation(basis_vect, basis_evaluator_vect)));
//	}
//
//	// ----------------------------------------------------------------------
//	// 
//	//                     LS_forward
//	//
//	// ---------------------------------------------------------------------
//	LS_ForwardAlgo ls_ForwardAlgo(McLmmPricer_CONSTPTR(new McLmmVanillaSwapPricer(mcLmm)));
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
//	size_t nbSimulation = 50;
//
//	mcLmm_LS.simulateLMM(nbSimulation);
//
//	ls_BackwardAlgo.do_BackwardAlgo(callableVanillaSwaption_PTR,mcLmm_LS.lmmSimualtionResults_);
//
//	mcLmm_LS.simulateLMM(nbSimulation);
//
//	std::pair<double, double> result = ls_ForwardAlgo.do_ForwardAlgo(	callableVanillaSwaption_PTR, 
//																		mcLmm_LS.lmmSimualtionResults_,
//																		ls_BackwardAlgo.getRegressions());
//
//	double price = result.first;
//	cout << price << endl;
//
//}
//
//
//void Test_Longstaff_Schwartz_CallableSwap()
//{
//	LMM::Index  indexStart = 2;		//1Y
//	LMM::Index  indexEnd   = 6;		//4Y
//	Tenor	floatingLegTenorType = Tenor::_6M;
//	Tenor	fixedLegTenorType    = Tenor::_1YR;
//	assert(indexStart%2==0&&indexEnd%2==0);
//	LMMTenorStructure_PTR lmmTenorStructure( new LMMTenorStructure(floatingLegTenorType, indexEnd/2));
//
//	std::vector<LMM::Index> exerciseDates;
//	exerciseDates.push_back(2);
//	exerciseDates.push_back(4);
//	exerciseDates.push_back(6);
//
//	double fwdRate = 0.02;
//	std::vector<double> initLiborValues(lmmTenorStructure->get_horizon()+1, fwdRate);
//	double strike = fwdRate;
//
//	McLmm_PTR mcLmm = getMcLmmExample(lmmTenorStructure, initLiborValues);
//
//	double nominal = 1.0;
//	GenericSwap_CONSTPTR genericVanillaSwap = InstrumentFactory::createVanillaSwap(	strike, 
//																					indexStart, 
//																					indexEnd, 
//																					floatingLegTenorType, 
//																					fixedLegTenorType,
//																					lmmTenorStructure,
//																					nominal);
//
//	CallableInstrument_PTR callableGenericSwap(new CallableGenericSwap(genericVanillaSwap, exerciseDates));
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
//	LS_BackwardAlgo ls_BackwardAlgo(McLmmPricer_CONSTPTR(new McLmmGenericSwapPricer(mcLmm)));
//
//	std::vector<LS::Regression>& regression_vect = ls_BackwardAlgo.getRegressions();
//	ls_BackwardAlgo.getExerciseDates()=exerciseDates;
//	
//
//	size_t fixedFloatingRatio = fixedLegTenorType.ratioTo(floatingLegTenorType);
//	for(size_t regressionIndex=0; regressionIndex<2; regressionIndex++)
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
//		std::vector<Basis_CONSTPTR> basis_vect;
//		basis_vect.push_back(Basis_CONSTPTR(new Basis_ConstUnity()));
//		basis_vect.push_back(Basis_CONSTPTR(new Basis_Monomial(	EV_CONSTPTR(new EV_LiborRate(Rate1_CONSTPTR(new LiborRate(liborIndex,Tenor(floatingLegTenorType))))), 
//																Basis_SinglaEVFunctional::Transformer_CONSTPTR(new Basis_Monomial::MonomialTransformer(1.0,1)))));
//		//basis_vect.push_back(Basis_CONSTPTR(new Basis_Monomial(	EV_CONSTPTR(new EV_LiborRate(Rate1_CONSTPTR(new LiborRate(liborIndex,Tenor(floatingLegTenorType))))), 
//		//														Basis_SinglaEVFunctional::Transformer_CONSTPTR(new Basis_Monomial::MonomialTransformer(1.0,2)))));
//
//		basis_vect.push_back(Basis_CONSTPTR(new Basis_CappedFlooredCoupon(	
//			EV_CONSTPTR(new EV_LiborRate(Rate1_CONSTPTR(new LiborRate(liborIndex,Tenor(floatingLegTenorType))))),
//			Basis_SinglaEVFunctional::Transformer_CONSTPTR(
//								new Basis_CappedFlooredCoupon::CappedFlooredCouponTransformer(
//									CappedFlooredCoupon_CONSTPTR(new CappedFlooredCoupon(CappedFlooredCoupon(	
//																												paymentIndex,
//																												1.0,
//																												1.0, 
//																												true,
//																												0.0,
//																												false, 
//																												10e100, 
//																												Rate1_CONSTPTR(new Rate1()),
//																												1.0,
//																												-strike, 
//																												liborIndex))))))));
//
//		basis_vect.push_back(Basis_CONSTPTR(new Basis_Monomial(				
//			EV_CONSTPTR(new EV_VanillaSwapRate( Rate1_CONSTPTR( new VanillaSwapRate(
//							VanillaSwap(			strike, 
//													liborIndex, 
//													indexEnd, 
//													floatingLegTenorType, 
//													fixedLegTenorType, 
//													lmmTenorStructure))))), 
//			Basis_SinglaEVFunctional::Transformer_CONSTPTR(new Basis_Monomial::MonomialTransformer(1.0,1)))));
//
//		basis_vect.push_back(Basis_CONSTPTR(new Basis_Monomial(				
//			EV_CONSTPTR(new EV_VanillaSwapRate( Rate1_CONSTPTR( new VanillaSwapRate(
//							VanillaSwap(			strike, 
//													liborIndex, 
//													indexEnd, 
//													floatingLegTenorType, 
//													fixedLegTenorType, 
//													lmmTenorStructure))))), 
//			Basis_SinglaEVFunctional::Transformer_CONSTPTR(new Basis_Monomial::MonomialTransformer(1.0,2)))));
//
//		//basis_vect.push_back(Basis_CONSTPTR(new Basis_Monomial(				
//		//	EV_CONSTPTR(new EV_VanillaSwapRate( Rate1_CONSTPTR( new VanillaSwapRate(
//		//					VanillaSwap(			strike, 
//		//											liborIndex, 
//		//											indexEnd, 
//		//											floatingLegTenorType, 
//		//											fixedLegTenorType, 
//		//											lmmTenorStructure))))), 
//		//	Basis_SinglaEVFunctional::Transformer_CONSTPTR(new Basis_Monomial::MonomialTransformer(10000.0,3)))));
//
//		//basis_vect.push_back(Basis_CONSTPTR(new Basis_Monomial(				
//		//	EV_CONSTPTR(new EV_VanillaSwapRate( Rate1_CONSTPTR( new VanillaSwapRate(
//		//					VanillaSwap(			strike, 
//		//											liborIndex, 
//		//											indexEnd, 
//		//											floatingLegTenorType, 
//		//											fixedLegTenorType, 
//		//											lmmTenorStructure))))), 
//		//	Basis_SinglaEVFunctional::Transformer_CONSTPTR(new Basis_Monomial::MonomialTransformer(100000.0,4)))));
//
//
//		basis_vect.push_back(Basis_CONSTPTR(new Basis_CappedFlooredCoupon(	
//			EV_CONSTPTR(new EV_VanillaSwapRate( Rate1_CONSTPTR( new VanillaSwapRate(
//							VanillaSwap(			strike, 
//													liborIndex, 
//													indexEnd, 
//													floatingLegTenorType, 
//													fixedLegTenorType, 
//													lmmTenorStructure))))), 
//			Basis_SinglaEVFunctional::Transformer_CONSTPTR(
//								new Basis_CappedFlooredCoupon::CappedFlooredCouponTransformer(
//									CappedFlooredCoupon_CONSTPTR(new CappedFlooredCoupon(CappedFlooredCoupon(	
//																												paymentIndex,
//																												1.0,
//																												1.0, 
//																												true,
//																												0.0,
//																												false, 
//																												10e100, 
//																												Rate1_CONSTPTR(new Rate1()),
//																												1.0,
//																												-strike, 
//																												liborIndex))))))));
//
//		basis_vect.push_back(Basis_CONSTPTR(new Basis_Composited(
//			
//			Basis_CONSTPTR(new Basis_CappedFlooredCoupon(	
//				EV_CONSTPTR(new EV_VanillaSwapRate( Rate1_CONSTPTR( new VanillaSwapRate(
//							VanillaSwap(			strike, 
//													liborIndex, 
//													indexEnd, 
//													floatingLegTenorType, 
//													fixedLegTenorType, 
//													lmmTenorStructure))))), 
//				Basis_SinglaEVFunctional::Transformer_CONSTPTR(
//								new Basis_CappedFlooredCoupon::CappedFlooredCouponTransformer(
//									CappedFlooredCoupon_CONSTPTR(new CappedFlooredCoupon(CappedFlooredCoupon(	
//																												paymentIndex,
//																												100.0,
//																												floatingLegTenorType.YearFraction(), 
//																												true,
//																												0.0,
//																												false, 
//																												10e100, 
//																												Rate1_CONSTPTR(new Rate1()),
//																												1.0,
//																												-strike, 
//																												liborIndex))))))),
//			Basis_CONSTPTR(new Basis_Monomial(				
//				EV_CONSTPTR(new EV_VanillaSwapRate( Rate1_CONSTPTR( new VanillaSwapRate(
//							VanillaSwap(			strike, 
//													liborIndex, 
//													indexEnd, 
//													floatingLegTenorType, 
//													fixedLegTenorType, 
//													lmmTenorStructure))))), 
//				Basis_SinglaEVFunctional::Transformer_CONSTPTR(new Basis_Monomial::MonomialTransformer(100.0,2)))),
//			Basis_Composited::Compositor::ADD)));
//	
//		std::vector<Basis_Evaluator_CONSTPTR> basis_evaluator_vect;
//		basis_evaluator_vect.push_back(	Basis_Evaluator_CONSTPTR(new Basis_SinglaEVFunctional_Evaluator(
//										EV_Evaluator_CONSTPTR(new EV_ConstRate_Evaluator()))));	
//		basis_evaluator_vect.push_back(	Basis_Evaluator_CONSTPTR(new Basis_SinglaEVFunctional_Evaluator(
//										EV_Evaluator_CONSTPTR(new EV_LiborRate_Evaluator()))));	
//		//basis_evaluator_vect.push_back(	Basis_Evaluator_CONSTPTR(new Basis_SinglaEVFunctional_Evaluator(
//		//								EV_Evaluator_CONSTPTR(new EV_LiborRate_Evaluator()))));	
//		basis_evaluator_vect.push_back(	Basis_Evaluator_CONSTPTR(new Basis_SinglaEVFunctional_Evaluator(
//										EV_Evaluator_CONSTPTR(new EV_LiborRate_Evaluator()))));	
//	    basis_evaluator_vect.push_back(	Basis_Evaluator_CONSTPTR(new Basis_SinglaEVFunctional_Evaluator(
//										EV_Evaluator_CONSTPTR(new EV_VanillaSwapRate_Evaluator(McLmmVanillaSwapPricer(mcLmm))))));
//		basis_evaluator_vect.push_back(	Basis_Evaluator_CONSTPTR(new Basis_SinglaEVFunctional_Evaluator(
//										EV_Evaluator_CONSTPTR(new EV_VanillaSwapRate_Evaluator(McLmmVanillaSwapPricer(mcLmm))))));
//	    //basis_evaluator_vect.push_back(	Basis_Evaluator_CONSTPTR(new Basis_SinglaEVFunctional_Evaluator(
//					//					EV_Evaluator_CONSTPTR(new EV_VanillaSwapRate_Evaluator(McLmmVanillaSwapPricer(mcLmm))))));
//	    //basis_evaluator_vect.push_back(	Basis_Evaluator_CONSTPTR(new Basis_SinglaEVFunctional_Evaluator(
//					//					EV_Evaluator_CONSTPTR(new EV_VanillaSwapRate_Evaluator(McLmmVanillaSwapPricer(mcLmm))))));
//		basis_evaluator_vect.push_back(	Basis_Evaluator_CONSTPTR(new Basis_SinglaEVFunctional_Evaluator(
//										EV_Evaluator_CONSTPTR(new EV_VanillaSwapRate_Evaluator(McLmmVanillaSwapPricer(mcLmm))))));
//		basis_evaluator_vect.push_back(	Basis_Evaluator_CONSTPTR(new Basis_Composited_Function_Evaluator(
//										Basis_Evaluator_CONSTPTR(new Basis_SinglaEVFunctional_Evaluator(
//											EV_Evaluator_CONSTPTR(new EV_VanillaSwapRate_Evaluator(McLmmVanillaSwapPricer(mcLmm))))),
//										Basis_Evaluator_CONSTPTR(new Basis_SinglaEVFunctional_Evaluator(
//											EV_Evaluator_CONSTPTR(new EV_VanillaSwapRate_Evaluator(McLmmVanillaSwapPricer(mcLmm)))))
//										)));
//
//		regression_vect.push_back(LS::Regression(LS::RegressionRepresentation(basis_vect, basis_evaluator_vect)));
//	}
//
//	// ----------------------------------------------------------------------
//	// 
//	//                     LS_forward
//	//
//	// ---------------------------------------------------------------------
//	LS_ForwardAlgo ls_ForwardAlgo(McLmmPricer_CONSTPTR(new McLmmGenericSwapPricer(mcLmm)));
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
//	std::vector<size_t> nbSimulation_vect;
//	//nbSimulation_vect.push_back(20);
//	nbSimulation_vect.push_back(500);
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
//		McLmm_LS mcLmm_LS(mcLmm);
//
//		mcLmm_LS.simulateLMM(nbSimulation);   //the same seed for each simulation.
//
//		ls_BackwardAlgo.do_BackwardAlgo(callableGenericSwap,mcLmm_LS.lmmSimualtionResults_);
//
//		mcLmm_LS.simulateLMM(nbSimulation);
//
//		std::pair<double, double> result = ls_ForwardAlgo.do_ForwardAlgo(	callableGenericSwap, 
//																			mcLmm_LS.lmmSimualtionResults_,
//																			ls_BackwardAlgo.getRegressions());
//
//
//		clock_t endTime = clock();
//		clock_t clockTicksTaken = endTime - startTime;	
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
//		price_vect.push_back(result.first);
//		variance_vect.push_back(result.second);
//		leftValue_vect.push_back(leftValue);
//		rightValue_vect.push_back(rightValue);
//		time_vect.push_back(timeInSeconds);
//		diff_vect.push_back(sqrt(variance/nbSimulation)*q);
//
//	}
//
//	//save 
//	std::stringstream outputFileName_s; 
//	outputFileName_s<<"Test_LS_different_basis_price"<<".csv";
//	std::string outputFileName = LMMPATH::get_Root_OutputPath() + outputFileName_s.str();
//	ofstream out;
//	out.open(outputFileName,  ios::out | ios::app );
//	//out.open(outputFileName,  ios::out);
//	out<<endl;
//	out<<endl;
//	out<<endl;
//	out << " 8 basis ;   " << endl;
//	out << " 0.001 ;  " << "L_i ;" <<  "(L_i)^2 ;" <<"max(L_i-K,0) ;" << "S_i ;" <<  "(S_i)^2 ;";  // << "(S_i)^3 ;" <<  "(S_i)^4 ;" ;
//	out <<"max(S_i-K, 0) ;" <<"(max(S_i-K, 0))^2 ;"<< endl;
//
//	size_t nbSimu_size = nbSimulation_vect.size();
//
//	//out << "Basis_Monomial; on Libor Rate: ; L_i" << endl;
//
//	out << "exerciseDates: ;" ;
//	for(size_t i = 0; i<exerciseDates.size(); i++)
//		out << exerciseDates[i] << " ;" ;
//	out << endl;
//	
//	out << "nbSimulation_vect: ;" ;
//	for(size_t i = 0; i<nbSimulation_vect.size(); i++)
//		out << nbSimulation_vect[i] << " ;" ;
//	out << endl;
//
//
//	out << "leftValue_vect: ;" ;
//	for(size_t i = 0; i<nbSimulation_vect.size(); i++)
//		out << leftValue_vect[i] << " ;" ;
//	out << endl;
//
//	out << "price_vect: ;" ;
//	for(size_t i = 0; i<nbSimulation_vect.size(); i++)
//		out << price_vect[i] << " ;" ;
//	out << endl;
//
//	out << "rightValue_vect: ;" ;
//	for(size_t i = 0; i<nbSimulation_vect.size(); i++)
//		out << rightValue_vect[i] << "; " ;
//	out << endl;
//
//	out << "diff_vect: ;" ;
//	for(size_t i = 0; i<nbSimulation_vect.size(); i++)
//		out << diff_vect[i] << " ;" ;
//	out << endl;
//
//	out << "variance_vect: ;" ;
//	for(size_t i = 0; i<nbSimulation_vect.size(); i++)
//		out << variance_vect[i] << " ;" ;
//	out << endl;
//
//	out << "time_vect: ;" ;
//	for(size_t i = 0; i<nbSimulation_vect.size(); i++)
//		out << time_vect[i] << " ;" ;
//	out << endl;
//
//
//	// ----------------------------------------------------------------------
//	// 
//	//               show on terminal
//	//
//	// ----------------------------------------------------------------------
//	cout << "exerciseDates: " ;
//	for(size_t i = 0; i<exerciseDates.size(); i++)
//		cout << exerciseDates[i] << " " ;
//	cout << endl;
//
//	cout << "nbSimulation_vect: " ;
//	for(size_t i = 0; i<nbSimulation_vect.size(); i++)
//		cout << nbSimulation_vect[i] << " " ;
//	cout << endl;
//
//	cout << "leftValue_vect: " ;
//	for(size_t i = 0; i<nbSimulation_vect.size(); i++)
//		cout << leftValue_vect[i] << " " ;
//	cout << endl;
//
//	cout << "price_vect: " ;
//	for(size_t i = 0; i<nbSimulation_vect.size(); i++)
//		cout << price_vect[i] << " " ;
//	cout << endl;
//
//	cout << "rightValue_vect: " ;
//	for(size_t i = 0; i<nbSimulation_vect.size(); i++)
//		cout << rightValue_vect[i] << " " ;
//	cout << endl;
//
//	cout << "diff_vect: " ;
//	for(size_t i = 0; i<nbSimulation_vect.size(); i++)
//		cout << diff_vect[i] << " " ;
//	cout << endl;
//
//	cout << "variance_vect: " ;
//	for(size_t i = 0; i<nbSimulation_vect.size(); i++)
//		cout << variance_vect[i] << " " ;
//	cout << endl;
//
//	cout << "time_vect: " ;
//	for(size_t i = 0; i<nbSimulation_vect.size(); i++)
//		cout << time_vect[i] << " " ;
//	cout << endl;
//}
//
//void Test_Longstaff_Schwartz_CallableSwap_for_test()
//{
//	LMM::Index  indexStart = 2;		//1Y
//	LMM::Index  indexEnd   = 8;		//4Y
//	Tenor	floatingLegTenorType = Tenor::_6M;
//	Tenor	fixedLegTenorType    = Tenor::_1YR;
//	assert(indexStart%2==0&&indexEnd%2==0);
//	LMMTenorStructure_PTR lmmTenorStructure( new LMMTenorStructure(floatingLegTenorType, indexEnd/2));
//
//	std::vector<LMM::Index> exerciseDates;
//	exerciseDates.push_back(2);
//	exerciseDates.push_back(4);
//	exerciseDates.push_back(6);
//	exerciseDates.push_back(8);
//
//	double fwdRate = 0.02;
//	std::vector<double> initLiborValues(lmmTenorStructure->get_horizon()+1, fwdRate);
//	double strike = fwdRate;
//
//	McLmm_PTR mcLmm = getMcLmmExample(lmmTenorStructure, initLiborValues);
//
//	double nominal = 1.0;
//	GenericSwap_CONSTPTR genericVanillaSwap = InstrumentFactory::createVanillaSwap(	strike, 
//																					indexStart, 
//																					indexEnd, 
//																					floatingLegTenorType, 
//																					fixedLegTenorType,
//																					lmmTenorStructure,
//																					nominal);
//
//	CallableInstrument_PTR callableGenericSwap(new CallableGenericSwap(genericVanillaSwap, exerciseDates));
//	
//    
//
//	// ----------------------------------------------------------------------
//	// 
//	//                     LS_Backward
//	//
//	// ----------------------------------------------------------------------
//	LS_BackwardAlgo ls_BackwardAlgo(McLmmPricer_CONSTPTR(new McLmmGenericSwapPricer(mcLmm)));
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
//		
//		std::vector<Basis_CONSTPTR> basis_vect;
//		basis_vect.push_back(Basis_CONSTPTR(new Basis_ConstUnity()));
//
//		basis_vect.push_back(Basis_CONSTPTR(new Basis_Monomial(				
//			EV_CONSTPTR(new EV_VanillaSwapRate( Rate1_CONSTPTR( new VanillaSwapRate(
//							VanillaSwap(			strike, 
//													liborIndex, 
//													indexEnd, 
//													floatingLegTenorType, 
//													fixedLegTenorType, 
//													lmmTenorStructure))))), 
//			Basis_SinglaEVFunctional::Transformer_CONSTPTR(new Basis_Monomial::MonomialTransformer(1.0,1)))));
//
//		basis_vect.push_back(Basis_CONSTPTR(new Basis_CappedFlooredCoupon(	
//			EV_CONSTPTR(new EV_VanillaSwapRate( Rate1_CONSTPTR( new VanillaSwapRate(
//							VanillaSwap(			strike, 
//													liborIndex, 
//													indexEnd, 
//													floatingLegTenorType, 
//													fixedLegTenorType, 
//													lmmTenorStructure))))),
//			Basis_SinglaEVFunctional::Transformer_CONSTPTR(
//								new Basis_CappedFlooredCoupon::CappedFlooredCouponTransformer(
//									CappedFlooredCoupon_CONSTPTR(new CappedFlooredCoupon(CappedFlooredCoupon(	
//																												paymentIndex,
//																												1.0,
//																												1.0, 
//																												true,
//																												0.0,
//																												false, 
//																												10e100, 
//																												Rate1_CONSTPTR(new Rate1()),
//																												1.0,
//																												-strike, 
//																												liborIndex))))))));
//
//
//
//		basis_vect.push_back(Basis_CONSTPTR(new Basis_Monomial(	EV_CONSTPTR(new EV_LiborRate(Rate1_CONSTPTR(new LiborRate(liborIndex,Tenor(floatingLegTenorType))))), 
//														Basis_SinglaEVFunctional::Transformer_CONSTPTR(new Basis_Monomial::MonomialTransformer(1.0,1)))));
//
//		//basis_vect.push_back(Basis_CONSTPTR(new Basis_Monomial(				
//		//	EV_CONSTPTR(new EV_VanillaSwapRate( Rate1_CONSTPTR( new VanillaSwapRate(
//		//					VanillaSwap(			strike, 
//		//											liborIndex, 
//		//											indexEnd, 
//		//											floatingLegTenorType, 
//		//											fixedLegTenorType, 
//		//											lmmTenorStructure))))), 
//		//	Basis_SinglaEVFunctional::Transformer_CONSTPTR(new Basis_Monomial::MonomialTransformer(1.0,2)))));
//
//		//basis_vect.push_back(Basis_CONSTPTR(new Basis_CappedFlooredCoupon(	
//		//	EV_CONSTPTR(new EV_VanillaSwapRate( Rate1_CONSTPTR( new VanillaSwapRate(
//		//					VanillaSwap(			strike, 
//		//											liborIndex, 
//		//											indexEnd, 
//		//											floatingLegTenorType, 
//		//											fixedLegTenorType, 
//		//											lmmTenorStructure))))), 
//		//	Basis_SinglaEVFunctional::Transformer_CONSTPTR(
//		//						new Basis_CappedFlooredCoupon::CappedFlooredCouponTransformer(
//		//							CappedFlooredCoupon_CONSTPTR(new CappedFlooredCoupon(CappedFlooredCoupon(	
//		//																										paymentIndex,
//		//																										1.0,
//		//																										1.0, 
//		//																										true,
//		//																										0.0,
//		//																										false, 
//		//																										10e100, 
//		//																										Rate1_CONSTPTR(new Rate1()),
//		//																										1.0,
//		//																										-strike, 
//		//																										liborIndex))))))));
//
//		//basis_vect.push_back(Basis_CONSTPTR(new Basis_Composited(
//		//	
//		//	Basis_CONSTPTR(new Basis_CappedFlooredCoupon(	
//		//		EV_CONSTPTR(new EV_VanillaSwapRate( Rate1_CONSTPTR( new VanillaSwapRate(
//		//					VanillaSwap(			strike, 
//		//											liborIndex, 
//		//											indexEnd, 
//		//											floatingLegTenorType, 
//		//											fixedLegTenorType, 
//		//											lmmTenorStructure))))), 
//		//		Basis_SinglaEVFunctional::Transformer_CONSTPTR(
//		//						new Basis_CappedFlooredCoupon::CappedFlooredCouponTransformer(
//		//							CappedFlooredCoupon_CONSTPTR(new CappedFlooredCoupon(CappedFlooredCoupon(	
//		//																										paymentIndex,
//		//																										1.0,
//		//																										1.0, 
//		//																										true,
//		//																										0.0,
//		//																										false, 
//		//																										10e100, 
//		//																										Rate1_CONSTPTR(new Rate1()),
//		//																										1.0,
//		//																										-strike, 
//		//																										liborIndex))))))),
//		//	Basis_CONSTPTR(new Basis_Monomial(				
//		//		EV_CONSTPTR(new EV_VanillaSwapRate( Rate1_CONSTPTR( new VanillaSwapRate(
//		//					VanillaSwap(			strike, 
//		//											liborIndex, 
//		//											indexEnd, 
//		//											floatingLegTenorType, 
//		//											fixedLegTenorType, 
//		//											lmmTenorStructure))))), 
//		//		Basis_SinglaEVFunctional::Transformer_CONSTPTR(new Basis_Monomial::MonomialTransformer(1.0,2)))),
//		//	Basis_Composited::Compositor::ADD)));
//	
//		std::vector<Basis_Evaluator_CONSTPTR> basis_evaluator_vect;
//		basis_evaluator_vect.push_back(	Basis_Evaluator_CONSTPTR(new Basis_SinglaEVFunctional_Evaluator(
//										EV_Evaluator_CONSTPTR(new EV_ConstRate_Evaluator()))));	
//	    basis_evaluator_vect.push_back(	Basis_Evaluator_CONSTPTR(new Basis_SinglaEVFunctional_Evaluator(
//										EV_Evaluator_CONSTPTR(new EV_VanillaSwapRate_Evaluator(McLmmVanillaSwapPricer(mcLmm))))));
//		basis_evaluator_vect.push_back(	Basis_Evaluator_CONSTPTR(new Basis_SinglaEVFunctional_Evaluator(
//										EV_Evaluator_CONSTPTR(new EV_VanillaSwapRate_Evaluator(McLmmVanillaSwapPricer(mcLmm))))));
//		basis_evaluator_vect.push_back(	Basis_Evaluator_CONSTPTR(new Basis_SinglaEVFunctional_Evaluator(
//										EV_Evaluator_CONSTPTR(new EV_LiborRate_Evaluator()))));	
//		//basis_evaluator_vect.push_back(	Basis_Evaluator_CONSTPTR(new Basis_SinglaEVFunctional_Evaluator(
//		//								EV_Evaluator_CONSTPTR(new EV_LiborRate_Evaluator()))));	
//	    //basis_evaluator_vect.push_back(	Basis_Evaluator_CONSTPTR(new Basis_SinglaEVFunctional_Evaluator(
//					//					EV_Evaluator_CONSTPTR(new EV_VanillaSwapRate_Evaluator(McLmmVanillaSwapPricer(mcLmm))))));
//
//		//basis_evaluator_vect.push_back(	Basis_Evaluator_CONSTPTR(new Basis_Composited_Function_Evaluator(
//		//								Basis_Evaluator_CONSTPTR(new Basis_SinglaEVFunctional_Evaluator(
//		//									EV_Evaluator_CONSTPTR(new EV_VanillaSwapRate_Evaluator(McLmmVanillaSwapPricer(mcLmm))))),
//		//								Basis_Evaluator_CONSTPTR(new Basis_SinglaEVFunctional_Evaluator(
//		//									EV_Evaluator_CONSTPTR(new EV_VanillaSwapRate_Evaluator(McLmmVanillaSwapPricer(mcLmm)))))
//		//								)));
//		regression_vect.push_back(LS::Regression(LS::RegressionRepresentation(basis_vect, basis_evaluator_vect)));
//	}
//
//	// ----------------------------------------------------------------------
//	// 
//	//                     LS_forward
//	//
//	// ---------------------------------------------------------------------
//	LS_ForwardAlgo ls_ForwardAlgo(McLmmPricer_CONSTPTR(new McLmmGenericSwapPricer(mcLmm)));
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
//	std::vector<size_t> nbSimulation_vect;
//	nbSimulation_vect.push_back(50);
//	//nbSimulation_vect.push_back(500);
//	//nbSimulation_vect.push_back(1000);
//	//nbSimulation_vect.push_back(2000);
//	//nbSimulation_vect.push_back(5000);
//	//nbSimulation_vect.push_back(10000);
//	//nbSimulation_vect.push_back(20000);
//	//nbSimulation_vect.push_back(50000);
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
//	//save 
//	std::stringstream outputFileName_s; 
//	outputFileName_s<<"Test_LS_allData_forTest"<<".csv";
//	std::string outputFileName = LMMPATH::get_Root_OutputPath() + outputFileName_s.str();
//	ofstream out;
//	//out.open(outputFileName,  ios::out | ios::app );
//	out.open(outputFileName,  ios::out);
//	out<<endl;
//	out<<endl;
//	out<<endl;
//	out << " 7 basis ;   " << endl;
//	out << " 1 ;  " << "L_i ;" << "max(L_i-K,0) ;" << "S_i ;" << "max(S_i-K, 0) ;" << "(S_i)^2 ;" << "(max(S_i-K, 0))^2 ;" <<endl ;
//
//	clock_t startTime = clock();
//
//	size_t nbSimulation = nbSimulation_vect[0];
//
//	McLmm_LS mcLmm_LS(mcLmm);
//
//	mcLmm_LS.simulateLMM(nbSimulation);
//	mcLmm_LS.write_to_stream(out);
//
//	ls_BackwardAlgo.do_BackwardAlgo(callableGenericSwap,mcLmm_LS.lmmSimualtionResults_);
//	ls_BackwardAlgo.write_to_stream(out);
//
//	mcLmm_LS.simulateLMM(nbSimulation);
//	mcLmm_LS.write_to_stream(out);
//
//	std::pair<double, double> result = ls_ForwardAlgo.do_ForwardAlgo(	callableGenericSwap, 
//																		mcLmm_LS.lmmSimualtionResults_,
//																		ls_BackwardAlgo.getRegressions());
//
//	ls_ForwardAlgo.write_to_stream(out);
//
//	clock_t endTime = clock();
//	clock_t clockTicksTaken = endTime - startTime;	
//	double timeInSeconds = clockTicksTaken / (double) CLOCKS_PER_SEC;
//
//	double price = result.first;
//	double variance = result.second;
//
//	boost::math::normal dist(0.0, 1.0);
//	double threshold = 0.99;			// 95% of distribution is below q:
//	double q = quantile(dist, threshold);
//	double leftValue = price - sqrt(variance/nbSimulation)*q;
//	double rightValue = price + sqrt(variance/nbSimulation)*q;
//
//	price_vect.push_back(result.first);
//	variance_vect.push_back(result.second);
//	leftValue_vect.push_back(leftValue);
//	rightValue_vect.push_back(rightValue);
//	time_vect.push_back(timeInSeconds);
//	diff_vect.push_back(sqrt(variance/nbSimulation)*q);
//
//	
//
//
//	size_t nbSimu_size = nbSimulation_vect.size();
//
//	out << "Basis_Monomial; on Libor Rate: ; L_i" << endl;
//
//	out << "exerciseDates: ;" ;
//	for(size_t i = 0; i<exerciseDates.size(); i++)
//		out << exerciseDates[i] << " ;" ;
//	out << endl;
//	
//	out << "nbSimulation_vect: ;" ;
//	for(size_t i = 0; i<nbSimulation_vect.size(); i++)
//		out << nbSimulation_vect[i] << " ;" ;
//	out << endl;
//
//
//	out << "leftValue_vect: ;" ;
//	for(size_t i = 0; i<nbSimulation_vect.size(); i++)
//		out << leftValue_vect[i] << " ;" ;
//	out << endl;
//
//	out << "price_vect: ;" ;
//	for(size_t i = 0; i<nbSimulation_vect.size(); i++)
//		out << price_vect[i] << " ;" ;
//	out << endl;
//
//	out << "rightValue_vect: ;" ;
//	for(size_t i = 0; i<nbSimulation_vect.size(); i++)
//		out << rightValue_vect[i] << "; " ;
//	out << endl;
//
//	out << "diff_vect: ;" ;
//	for(size_t i = 0; i<nbSimulation_vect.size(); i++)
//		out << diff_vect[i] << " ;" ;
//	out << endl;
//
//	out << "variance_vect: ;" ;
//	for(size_t i = 0; i<nbSimulation_vect.size(); i++)
//		out << variance_vect[i] << " ;" ;
//	out << endl;
//
//	out << "time_vect: ;" ;
//	for(size_t i = 0; i<nbSimulation_vect.size(); i++)
//		out << time_vect[i] << " ;" ;
//	out << endl;
//
//
//	// ----------------------------------------------------------------------
//	// 
//	//               show on terminal
//	//
//	// ----------------------------------------------------------------------
//	cout << "exerciseDates: " ;
//	for(size_t i = 0; i<exerciseDates.size(); i++)
//		cout << exerciseDates[i] << " " ;
//	cout << endl;
//
//	cout << "nbSimulation_vect: " ;
//	for(size_t i = 0; i<nbSimulation_vect.size(); i++)
//		cout << nbSimulation_vect[i] << " " ;
//	cout << endl;
//
//	cout << "leftValue_vect: " ;
//	for(size_t i = 0; i<nbSimulation_vect.size(); i++)
//		cout << leftValue_vect[i] << " " ;
//	cout << endl;
//
//	cout << "price_vect: " ;
//	for(size_t i = 0; i<nbSimulation_vect.size(); i++)
//		cout << price_vect[i] << " " ;
//	cout << endl;
//
//	cout << "rightValue_vect: " ;
//	for(size_t i = 0; i<nbSimulation_vect.size(); i++)
//		cout << rightValue_vect[i] << " " ;
//	cout << endl;
//
//	cout << "diff_vect: " ;
//	for(size_t i = 0; i<nbSimulation_vect.size(); i++)
//		cout << diff_vect[i] << " " ;
//	cout << endl;
//
//	cout << "variance_vect: " ;
//	for(size_t i = 0; i<nbSimulation_vect.size(); i++)
//		cout << variance_vect[i] << " " ;
//	cout << endl;
//
//	cout << "time_vect: " ;
//	for(size_t i = 0; i<nbSimulation_vect.size(); i++)
//		cout << time_vect[i] << " " ;
//	cout << endl;
//}
//
//
//void Test_VanillaSwaptionPricing()
//{
//	LMM::Index  indexStart = 2;		//1Y
//	LMM::Index  indexEnd   = 6;		//4Y
//	Tenor	floatingLegTenorType = Tenor::_6M;
//	Tenor	fixedLegTenorType    = Tenor::_1YR;
//	assert(indexStart%2==0&&indexEnd%2==0);
//	LMMTenorStructure_PTR lmmTenorStructure( new LMMTenorStructure(floatingLegTenorType, indexEnd/2));
//
//	double fwdRate = 0.02;
//	std::vector<double> initLiborValues(lmmTenorStructure->get_horizon()+1, fwdRate);
//	McLmm_PTR mcLmm = getMcLmmExample(lmmTenorStructure, initLiborValues);
//	double strike = fwdRate;
//	VanillaSwap_CONSTPTR vanillaSwap_PTR(new VanillaSwap(strike, 
//													indexStart, 
//													indexEnd, 
//													floatingLegTenorType, 
//													fixedLegTenorType, 
//													lmmTenorStructure));
//
//	VanillaSwaption vanillaSwaption(*vanillaSwap_PTR.get() , OptionType::OptionType::CALL);
//
//	McLmmVanillaSwaptionPricer mcLmmVanillaSwaptionPricer(mcLmm);
//
//	std::vector<size_t> nbSimulation_vect;
//	nbSimulation_vect.push_back(500);
//	nbSimulation_vect.push_back(1000);
//	nbSimulation_vect.push_back(2000);
//	nbSimulation_vect.push_back(5000);
//	nbSimulation_vect.push_back(10000);
//	nbSimulation_vect.push_back(20000);
//	nbSimulation_vect.push_back(50000);
//
//	std::vector<double> price_vect;
//	//std::vector<double> leftValue_vect;
//	//std::vector<double> rightValue_vect;
//	std::vector<double> time_vect;
//	//std::vector<double> variance_vect;
//	//std::vector<double> diff_vect;
//
//	for(size_t i = 0; i<nbSimulation_vect.size(); i++)
//	{
//		clock_t startTime = clock();
//
//		size_t nbSimulation = nbSimulation_vect[i];
//
//		//this is the part of LS backward, let the first simulation go
//		mcLmmVanillaSwaptionPricer.price(vanillaSwaption, nbSimulation);
//
//		double price = mcLmmVanillaSwaptionPricer.price(vanillaSwaption, nbSimulation);
//
//		clock_t endTime = clock();
//		clock_t clockTicksTaken = endTime - startTime;	
//		double timeInSeconds = clockTicksTaken / (double) CLOCKS_PER_SEC;
//
//		price_vect.push_back(price);
//		time_vect.push_back(timeInSeconds);
//
//		//double price = result.first;
//		//double variance = result.second;
//
//		//boost::math::normal dist(0.0, 1.0);
//		//double threshold = 0.99;			// 95% of distribution is below q:
//		//double q = quantile(dist, threshold);
//		//double leftValue = price - sqrt(variance/nbSimulation)*q;
//		//double rightValue = price + sqrt(variance/nbSimulation)*q;
//
//		//price_vect.push_back(result.first);
//		//variance_vect.push_back(result.second);
//		//leftValue_vect.push_back(leftValue);
//		//rightValue_vect.push_back(rightValue);
//		//time_vect.push_back(timeInSeconds);
//		//diff_vect.push_back(sqrt(variance/nbSimulation)*q);
//	}
//
//	size_t nbSimu_size = nbSimulation_vect.size();
//
//	cout << "nbSimulation_vect: " ;
//	for(size_t i = 0; i<nbSimulation_vect.size(); i++)
//		cout << nbSimulation_vect[i] << " " ;
//	cout << endl;
//
//	//cout << "leftValue_vect: " ;
//	//for(size_t i = 0; i<nbSimulation_vect.size(); i++)
//	//	cout << leftValue_vect[i] << " " ;
//	//cout << endl;
//
//	cout << "price_vect: " ;
//	for(size_t i = 0; i<nbSimulation_vect.size(); i++)
//		cout << price_vect[i] << " " ;
//	cout << endl;
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
//	cout << "time_vect: " ;
//	for(size_t i = 0; i<nbSimulation_vect.size(); i++)
//		cout << time_vect[i] << " " ;
//	cout << endl;
//
//
//	//save 
//	std::stringstream outputFileName_s; 
//	outputFileName_s<<"Test_Longstaff_Schwartz_CallableSwap"<<".csv";
//	std::string outputFileName = LMMPATH::get_Root_OutputPath() + outputFileName_s.str();
//	ofstream out;
//	out.open(outputFileName,  ios::out | ios::app );
//	out<<endl;
//	out<<endl;
//	out<<endl;
//
//	out << "VanillaSwaptionPricer 1Y*2Y" << endl;
//	
//	out << "nbSimulation_vect: ;" ;
//	for(size_t i = 0; i<nbSimulation_vect.size(); i++)
//		out << nbSimulation_vect[i] << " ;" ;
//	out << endl;
//
//
//	//out << "leftValue_vect: ;" ;
//	//for(size_t i = 0; i<nbSimulation_vect.size(); i++)
//	//	out << leftValue_vect[i] << " ;" ;
//	//out << endl;
//
//	out << "price_vect: ;" ;
//	for(size_t i = 0; i<nbSimulation_vect.size(); i++)
//		out << price_vect[i] << " ;" ;
//	out << endl;
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
//	out << "time_vect: ;" ;
//	for(size_t i = 0; i<nbSimulation_vect.size(); i++)
//		out << time_vect[i] << " ;" ;
//	out << endl;
//}





