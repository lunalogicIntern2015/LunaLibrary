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


void Test_regression_LS()
{
	LMM::Index  indexStart = 2;		//1Y
	LMM::Index  indexEnd   = 6;		//4Y
	Tenor	floatingLegTenorType = Tenor::_6M;
	Tenor	fixedLegTenorType    = Tenor::_1YR;
	assert(indexStart%2==0&&indexEnd%2==0);
	LMMTenorStructure_PTR lmmTenorStructure( new LMMTenorStructure(floatingLegTenorType, indexEnd/2));

	std::vector<LMM::Index> exerciseDates;
	exerciseDates.push_back(2);
	exerciseDates.push_back(4);
	exerciseDates.push_back(6);

	double fwdRate = 0.02;
	std::vector<double> initLiborValues(lmmTenorStructure->get_horizon()+1, fwdRate);
	double strike = fwdRate;

	McLmm_PTR mcLmm = getMcLmmExample(lmmTenorStructure, initLiborValues);

	double nominal = 1.0;
	GenericSwap_CONSTPTR genericVanillaSwap = InstrumentFactory::createVanillaSwap(	strike, 
																					indexStart, 
																					indexEnd, 
																					floatingLegTenorType, 
																					fixedLegTenorType,
																					lmmTenorStructure,
																					nominal);

	CallableInstrument_PTR callableGenericSwap(new CallableGenericSwap(genericVanillaSwap, exerciseDates));

	// ----------------------------------------------------------------------
	// 
	//                     LS_Backward
	//
	// ---------------------------------------------------------------------
	
    McLmm_LS mcLmm_LS(mcLmm);

	// ----------------------------------------------------------------------
	// 
	//                     LS_Backward
	//
	// ----------------------------------------------------------------------
	LS_BackwardAlgo ls_BackwardAlgo(McLmmPricer_CONSTPTR(new McLmmGenericSwapPricer(mcLmm)));

	std::vector<LS::Regression>& regression_vect = ls_BackwardAlgo.getRegressions();
	ls_BackwardAlgo.getExerciseDates()=exerciseDates;
	
	size_t nbRegression = exerciseDates.size()-1;
	size_t fixedFloatingRatio = fixedLegTenorType.ratioTo(floatingLegTenorType);
	for(size_t regressionIndex=0; regressionIndex<nbRegression; regressionIndex++)
	{
		LMM::Index liborIndex	= indexStart + regressionIndex*fixedFloatingRatio;
		LMM::Index paymentIndex = indexStart + regressionIndex*fixedFloatingRatio + 1;
		
		std::vector<Basis_CONSTPTR> basis_vect;
		basis_vect.push_back(Basis_CONSTPTR(new Basis_Monomial(	EV_CONSTPTR(new EV_LiborRate(Rate1_CONSTPTR(new LiborRate(liborIndex,Tenor(floatingLegTenorType))))), 
																Basis_SinglaEVFunctional::Transformer_CONSTPTR(new Basis_Monomial::MonomialTransformer(1.0,1)))));


		basis_vect.push_back(Basis_CONSTPTR(new Basis_CappedFlooredCoupon(	
			EV_CONSTPTR(new EV_LiborRate(Rate1_CONSTPTR(new LiborRate(liborIndex,Tenor(floatingLegTenorType))))),
			Basis_SinglaEVFunctional::Transformer_CONSTPTR(
								new Basis_CappedFlooredCoupon::CappedFlooredCouponTransformer(
									CappedFlooredCoupon_CONSTPTR(new CappedFlooredCoupon(CappedFlooredCoupon(	
																												paymentIndex,
																												1.0,
																												floatingLegTenorType.YearFraction(), 
																												true,
																												0.0,
																												false, 
																												10e100, 
																												Rate1_CONSTPTR(new Rate1()),
																												1.0,
																												-strike, 
																												liborIndex))))))));

		basis_vect.push_back(Basis_CONSTPTR(new Basis_CappedFlooredCoupon(	
			EV_CONSTPTR(new EV_VanillaSwapRate( Rate1_CONSTPTR( new VanillaSwapRate(
							VanillaSwap(			strike, 
													liborIndex, 
													indexEnd, 
													floatingLegTenorType, 
													fixedLegTenorType, 
													lmmTenorStructure))))), 
			Basis_SinglaEVFunctional::Transformer_CONSTPTR(
								new Basis_CappedFlooredCoupon::CappedFlooredCouponTransformer(
									CappedFlooredCoupon_CONSTPTR(new CappedFlooredCoupon(CappedFlooredCoupon(	
																												paymentIndex,
																												1.0,
																												floatingLegTenorType.YearFraction(), 
																												true,
																												0.0,
																												false, 
																												10e100, 
																												Rate1_CONSTPTR(new Rate1()),
																												1.0,
																												-strike, 
																												liborIndex))))))));
	
		std::vector<Basis_Evaluator_CONSTPTR> basis_evaluator_vect;
		basis_evaluator_vect.push_back(	Basis_Evaluator_CONSTPTR(new Basis_SinglaEVFunctional_Evaluator(
										EV_Evaluator_CONSTPTR(new EV_LiborRate_Evaluator()))));	
		basis_evaluator_vect.push_back(	Basis_Evaluator_CONSTPTR(new Basis_SinglaEVFunctional_Evaluator(
								EV_Evaluator_CONSTPTR(new EV_LiborRate_Evaluator()))));	
		basis_evaluator_vect.push_back(	Basis_Evaluator_CONSTPTR(new Basis_SinglaEVFunctional_Evaluator(
									EV_Evaluator_CONSTPTR(new EV_VanillaSwapRate_Evaluator(McLmmVanillaSwapPricer(mcLmm))))));

		regression_vect.push_back(LS::Regression(LS::RegressionRepresentation(basis_vect, basis_evaluator_vect)));
	}


	//save 
	std::stringstream outputFileName_s; 
	outputFileName_s<<"Test_regression_LS"<<".csv";
	std::string outputFileName = LMMPATH::get_Root_OutputPath() + outputFileName_s.str();
	ofstream out;
	//out.open(outputFileName,  ios::out | ios::app );
	out.open(outputFileName,  ios::out);
	out<<endl;
	out<<endl;
	out<<endl;

	size_t nbSimulation = 10;

	mcLmm_LS.simulateLMM(nbSimulation);
	mcLmm_LS.write_to_stream(out);

	ls_BackwardAlgo.do_BackwardAlgo(callableGenericSwap,mcLmm_LS.lmmSimualtionResults_);
	ls_BackwardAlgo.write_to_stream(out);

	const std::vector<LS::Regression>& regressions = ls_BackwardAlgo.getRegressions();
	cout << "regressionCoeffs_: " << endl;
	for(size_t i = 0; i<regressions.size(); i++)
	{
		cout << " " ;
		const std::vector<double>& rgCoef = regressions[i].getRegressionRepresentation().getRegressionCoef();
		for(size_t j = 0; j<rgCoef.size(); j++)
			cout << rgCoef[j] << " " ;
		cout << endl;
	}

	//0.00692709, -1.65404, 3.02145
	//valid getSubInstrument
	//valid price_for_oneSimulation
}