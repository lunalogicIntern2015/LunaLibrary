//#include <JBLMM/Test/JBTests.h>
//
//#include <iostream>
//#include <fstream> 
//#include <string> 
//
//
//#include <LMM/RNGenerator/McGenerator.h>
//#include <LMM/RNGenerator/RNGenerator.h>
//
//#include <LMM/LmmModel/Correlation.h>
//#include <LMM/LmmModel/McTerminalLmm.h>
//
//#include <LMM/LmmModel/Lmm.h>
//#include <LMM/LmmModel/ConstShifted_HGVolatilityFunction.h>
//
//#include <JBLMM/Pricer/McGeneticSwapLMMPricer.h>
//#include <JBLMM/Instrument/InstrumentFactory.h>
//#include <JBLMM/Pricer/GeneticVanillaSwapPricer.h>
//#include <JBLMM/Longstaff_Schwartz/ExplanatoryVariable.h>
//
//#include <time.h>
//
//void Test_EV()
//{
//	double strike = 0.02;
//	LMM::Index  indexStart = 2;    //1Y
//	LMM::Index  indexEnd   = 6;	//3Y
//	assert(indexEnd%2==0);
//	Tenor	floatingLegTenorType = Tenor::_6M;
//	Tenor	fixedLegTenorType    = Tenor::_1YR;
//	Tenor	tenorStruture	=	Tenor::_6M;
//	size_t	horizonYear		=	10;
//	LMMTenorStructure_PTR lmmTenorStructure( new LMMTenorStructure(tenorStruture, indexEnd/2));
//
//	VanillaSwap_CONSTPTR vanillaSwap_PTR(new VanillaSwap(strike, indexStart, indexEnd, floatingLegTenorType, fixedLegTenorType, lmmTenorStructure));
//	LMM::Index liborExercice = 4;
//	std::vector<LMM::Index> exerciceIndex;
//	exerciceIndex.push_back(liborExercice);	
//
//
//	BermudanSwaption_PTR bermudanSwaption_PTR(new BermudanSwaption(vanillaSwap_PTR,exerciceIndex));
//
//
//		//---------------------------Build Lmm and McLmm's structure--------------------------------------
//	//! Parameter of h
//	double a = -0.06;
//	double b = 0.17;
//	double c = 0.54;
//	double d = 0.17;
//	Shifted_HGVolatilityParam::ABCDParameter abcdParam (a,b,c,d);
//	//Parameter of hg
//	double g_constParam = 1.0;
//	double shift_constParam = -0.01;
//	ConstShifted_HGVolatilityParam_PTR hgParam( new ConstShifted_HGVolatilityParam(lmmTenorStructure,abcdParam,g_constParam,shift_constParam));
//	
//	//! Correlation 1
//	size_t nbFactor       = 3; // need to test nbFactor  = 3, and nbFactor = 
//	size_t correlFullRank = lmmTenorStructure->get_horizon()+1;
//	size_t correlReducedRank = nbFactor;
//	//!"Check Parameters(): Condition not implemented yet."
//	std::cout << "checkParams(): ";
//	CorrelationReductionType::CorrelationReductionType correlReductionType = CorrelationReductionType::PCA;
//	double correlAlpha = 0.0;
//	double correlBeta  = 0.1;
//	Correlation_PTR correlation(new XY_beta_Correlation(correlFullRank,correlReducedRank, correlReductionType,correlAlpha,correlBeta));
//	correlation->calculate(); // for print.
//	correlation->print("test_McTerminalLmm_Correlation.csv");
//	//hgVolatilityFunction
//	ConstShifted_HGVolatilityFunction_PTR hgVolatilityFunction (new ConstShifted_HGVolatilityFunction(lmmTenorStructure, correlation, hgParam)); 
//	hgVolatilityFunction->print("test_McTerminalLmm_Volatility.csv");
//	//! Dispersion
//	Dispersion dispersion(hgVolatilityFunction);
//
//	unsigned long seed = 5033;
//	RNGenerator_PTR  rnGenerator(new McGenerator(seed));
//
//	//build lmm and mcLmm model
//	Lmm_PTR shiftedLmm (new Lmm(dispersion));
//	double fwdRate=0.02;
//	std::vector<double> liborsInitValue(lmmTenorStructure->get_horizon()+1, fwdRate);
//	McLmm_PTR mcLmm(new McTerminalLmm(shiftedLmm, liborsInitValue, rnGenerator, MCSchemeType::EULER));
//
//	mcLmm->simulateLMM();
//
//	std::stringstream outputFileName_s; 
//	outputFileName_s<<"Test_EV"<<".csv";
//	std::string outputFileName = LMMPATH::get_Root_OutputPath() + outputFileName_s.str();
//	ofstream out;
//	out.open(outputFileName,  ios::out | ios::app );
//	out<<endl;
//	out<<endl;
//	out<<endl;
//	const matrix& gMatrix = mcLmm->get_liborMatrix();
//	out << "gmatrix: " << endl; 
//	for(size_t i = 0; i<gMatrix.size1(); i++)
//	{
//		for(size_t j = 0; j<gMatrix.size2(); j++)
//			out << gMatrix(i,j) << "; ";
//		out << endl;
//	}
//
//	const std::vector<double>& numeraire = mcLmm->get_numeraire();
//	out << "numeraire: " << endl; 
//	for(size_t i = 0; i<numeraire.size(); i++)
//		out << numeraire[i] << "; ";
//	out << endl;
//
//
//	McLmmVanillaSwapPricer_PTR mcLmmVanillaSwapPricer_PTR(new McLmmVanillaSwapPricer(mcLmm));
//	ExplanatoryVariable::Evaluator_CONSTPTR evaluator_VanillaSwapRate_CONSTPTR(new EV_vanillaswaprate(mcLmmVanillaSwapPricer_PTR));
//	ExplanatoryVariable_CONSTPTR ev_VanillaSwapRate_CONSTPTR(new ExplanatoryVariable(evaluator_VanillaSwapRate_CONSTPTR));
//	double vanillaSwapRate=ev_VanillaSwapRate_CONSTPTR->getEvaluator_PTR()->evaluate(mcLmm, bermudanSwaption_PTR, liborExercice);
//	out << "vanillaSwapRate: ;" << endl; 
//	out << vanillaSwapRate <<";"<< endl; 
//
//	ExplanatoryVariable::Evaluator_CONSTPTR evaluator_libor_CONSTPTR(new EV_Libor());
//	ExplanatoryVariable_CONSTPTR ev_libor_CONSTPTR(new ExplanatoryVariable(evaluator_libor_CONSTPTR));
//	double libor=ev_libor_CONSTPTR->getEvaluator_PTR()->evaluate(mcLmm, bermudanSwaption_PTR, liborExercice);
//	out << "libor: ;" << endl; 
//	out << libor << ";"<<endl; 
//
//}