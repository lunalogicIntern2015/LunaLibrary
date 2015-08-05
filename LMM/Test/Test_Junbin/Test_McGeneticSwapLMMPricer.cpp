//#include "JBLMM/Test/JBTests.h"
//
//#include <iostream> 
//#include <fstream> 
//#include <string> 
//
//#include <RNGenerator/McGenerator.h>
//#include <RNGenerator/RNGenerator.h>
//
//#include <LMM/Model/Correlation.h>
//#include <LMM/Mc/McTerminalLmm.h>
//
//#include <LMM/Model/Lmm.h>
//#include <LMM/Model/ConstShifted_HGVolatilityFunction.h>
//
//#include <JBLMM/Pricer/McLmmGenericSwapPricer.h>
//#include <JBInstrument/InstrumentFactory.h>
//#include <JBLMM/Pricer/GenericVanillaSwapPricer.h>
//
//#include <time.h>
//
//
//void Test_McGeneticSwapLMMPricer()
//{
//
//	//! Parameters
//	double	strike			=	0.02;
//	LMM::Index	indexStart	=	0;
//	LMM::Index	indexEnd	=	20;
//	Tenor	floatingTenor	=	Tenor::_6M;
//	Tenor	fixedTenor		=	Tenor::_12M;
//	Tenor	tenorStruture	=	Tenor::_6M;
//	size_t	horizonYear		=	10;
//	LMMTenorStructure_PTR lmmTenorStructure( new LMMTenorStructure(tenorStruture, horizonYear));
//
//	cout	<<	"strike:                        "	<<	strike												<<	endl;
//	cout	<<	"indexStart:                    "	<<	indexStart											<<	endl;
//	cout	<<	"indexEnd:                      "	<<	indexEnd											<<	endl;
//	cout	<<	"tenorStrutureYearFraction:     "	<<	lmmTenorStructure->get_tenorType().YearFraction()	<<	endl;
//	cout	<<	"floatingVStenorStrutureRatio:  "	<<	floatingTenor.ratioTo(tenorStruture)				<<	endl;
//	cout	<<	"fixedVStenorStrutureRatio:     "	<<	fixedTenor.ratioTo(tenorStruture)					<<	endl;
//
//	double	fwdRate		=	0.02;
//	std::vector<double> liborsInitValue(lmmTenorStructure->get_horizon()+1, fwdRate);
//	cout	<<	"myInitialLibor:  ";
//	for (size_t i = 0; i <liborsInitValue.size(); i++)
//	{
//		cout	<<	liborsInitValue[i]	<<	" ";
//	}
//	cout	<<	 endl;
//	cout	<<	 endl;
//
//
//	//VanillaSwap_Chi_Trang
//	VanillaSwap firstVersionVanillaSwap(strike, indexStart , indexEnd, floatingTenor, fixedTenor, lmmTenorStructure);
//	LmmVanillaSwapPricer myVSP(lmmTenorStructure);
//	double FirstVersionSwapPrice	=	myVSP.swapNPV_Analytical_1(firstVersionVanillaSwap, liborsInitValue);
//
//	//---------------------------Build Lmm and McLmm's structure--------------------------------------
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
//	McLmm_PTR mcLmm(new McTerminalLmm(shiftedLmm, liborsInitValue, rnGenerator, MCSchemeType::EULER));
//
//	//build a McGeneticSwapLMMPricer
//	McLmmGenericSwapPricer_PTR mcGeneticSwapLMMPricer(new McLmmGenericSwapPricer(mcLmm));
//
//	//build the geneticVanillaSwap
//	GenericSwap_CONSTPTR vanillaSwap_Genetic=InstrumentFactory::createVanillaSwap(
//			strike,indexStart,indexEnd,floatingTenor,fixedTenor,lmmTenorStructure,1.0);
//
//	//use Monte Carlo Method
//	size_t nbSimulation=10;
//	double	MonteCarloPrice		=	mcGeneticSwapLMMPricer->swapNPV(vanillaSwap_Genetic, nbSimulation);
//
//	//ordinaryGeneticVanillaSwapPricer
//	GenericVanillaSwapPricer_PTR genericVanillaSwapPricer(new GenericVanillaSwapPricer());
//	double	OrdinaryGeneticVanillaSwapPrice		=	genericVanillaSwapPricer->genericVanillaSwap_Analytical(vanillaSwap_Genetic,liborsInitValue);
//
//	//subVanillaSwap
//	LMM::Index	subIndexStart	=	10;
//	LMM::Index	subIndexEnd		=	16;
//	GenericSwap_CONSTPTR subVanillaSwap_Genetic=InstrumentFactory::createVanillaSwap(
//			strike,subIndexStart,subIndexEnd,floatingTenor,fixedTenor,lmmTenorStructure,1.0);
//	double	subMonteCarloPrice	=	mcGeneticSwapLMMPricer->swapNPV(subVanillaSwap_Genetic, nbSimulation);
//	//subOrdinaryGeneticVanillaSwapPrice
//	double	subOrdinaryGeneticVanillaSwapPrice		=	genericVanillaSwapPricer->genericVanillaSwap_Analytical(subVanillaSwap_Genetic,liborsInitValue);
//	//subFirstVersionVanillaSwapPrice
//	VanillaSwap subFirstVersionVanillaSwap(strike, subIndexStart , subIndexEnd, floatingTenor, fixedTenor, lmmTenorStructure);
//	double subFirstVersionSwapPrice		=	myVSP.swapNPV_Analytical_1(subFirstVersionVanillaSwap, liborsInitValue);
//
//
//	cout	<<	"MonteCarloPrice: "																		<<	MonteCarloPrice													<<	endl;
//	cout	<<	"OrdinaryGeneticVanillaSwapPrice: "														<<	OrdinaryGeneticVanillaSwapPrice									<<	endl;
//	cout	<< "FirstVersionSwapPrice: "																<< FirstVersionSwapPrice											<< endl;
//	cout	<<	"Difference between MonteCarloPrice and OrdinaryGeneticVanillaSwapPrice: "				<<	MonteCarloPrice-OrdinaryGeneticVanillaSwapPrice					<<	endl;	
//	cout	<<	"Difference between OrdinaryGeneticVanillaSwapPrice and FirstVersionSwapPrice: "		<<	OrdinaryGeneticVanillaSwapPrice-FirstVersionSwapPrice			<<	endl;	
//	cout	<<	"subMonteCarloPrice: "																	<<	subMonteCarloPrice												<<	endl;
//	cout	<<	"subOrdinaryGeneticVanillaSwapPrice: "													<<	subOrdinaryGeneticVanillaSwapPrice								<<	endl;
//	cout	<< "subFirstVersionSwapPrice: "																<< subFirstVersionSwapPrice											<< endl;
//	cout	<<	"Difference between subMonteCarloPrice and subOrdinaryGeneticVanillaSwapPrice: "		<<	subMonteCarloPrice-subOrdinaryGeneticVanillaSwapPrice			<<	endl;	
//	cout	<<	"Difference between subOrdinaryGeneticVanillaSwapPrice and subFirstVersionSwapPrice: "	<<	subOrdinaryGeneticVanillaSwapPrice-subFirstVersionSwapPrice		<<	endl;
//
//
//
//	time_t _time;
//	struct tm timeInfo;
//	char format[32];
//	time(&_time);
//	localtime_s(&timeInfo, &_time);
//	strftime(format, 32, "%Y-%m-%d", &timeInfo);
//
//	std::stringstream outputFileName_s; 
//	outputFileName_s<<"TestResult_GeneticVanillaSwap_"<<format<<".csv";
//	std::string outputFileName = LMMPATH::get_Root_OutputPath() + outputFileName_s.str();
//
//
//	ofstream o;
//	o.open(outputFileName,  ios::out | ios::app );
//	o	<<	endl;
//	o	<<	endl;
//	o	<<	endl;
//	o	<<	"strike: "									<<	strike														<<	endl;
//	o	<<	"indexStart: "								<<	indexStart													<<	endl;
//	o	<<	"indexEnd: "								<<	indexEnd													<<	endl;
//	o	<<	"floatingTenorVSLmmStructureTenorRatio: "	<<	floatingTenor.ratioTo(lmmTenorStructure->get_tenorType())	<<	endl;
//	o	<<	"fixedTenorVSLmmStructureTenorRatio: "		<<	fixedTenor.ratioTo(lmmTenorStructure->get_tenorType())		<<	endl;
//	o	<<	"floatingTenorYearFraction: "				<<	floatingTenor.YearFraction()								<<	endl;
//	o	<<	"fixedTenorYearFraction: "					<<	fixedTenor.YearFraction()									<<	endl;
//	o	<<	"horizonYear: "								<<	horizonYear													<<	endl;
//	o	<<	"liborsInitValue: ";
//	for(size_t i=0; i<liborsInitValue.size(); i++)
//	{
//		o	<<	liborsInitValue[i]	<<	" ";
//	}
//	o	<<	endl;
//	o	<<	"nbSimulation: "																		<<	nbSimulation													<<	endl;
//	o	<<	"PRICES: "	<<	endl;
//	o	<<	"MonteCarloPrice: "																		<<	MonteCarloPrice													<<	endl;
//	o	<<	"OrdinaryGeneticVanillaSwapPrice: "														<<	OrdinaryGeneticVanillaSwapPrice									<<	endl;
//	o	<<	"FirstVersionSwapPrice: "																<< FirstVersionSwapPrice											<<	endl;
//	o	<<	"Difference between MonteCarloPrice and OrdinaryGeneticVanillaSwapPrice: "				<<	MonteCarloPrice-OrdinaryGeneticVanillaSwapPrice					<<	endl;	
//	o	<<	"Difference between OrdinaryGeneticVanillaSwapPrice and FirstVersionSwapPrice: "		<<	OrdinaryGeneticVanillaSwapPrice-FirstVersionSwapPrice			<<	endl;
//	o	<<	"SUBPRICES: "<<	endl;
//	o	<<	"subIndexStart: "																		<<	subIndexStart													<<	endl;
//	o	<<	"subIndexEnd: "																			<<	subIndexEnd														<<	endl;
//	o	<<	"subMonteCarloPrice: "																	<<	subMonteCarloPrice												<<	endl;
//	o	<<	"subOrdinaryGeneticVanillaSwapPrice: "													<<	subOrdinaryGeneticVanillaSwapPrice								<<	endl;
//	o	<<	"subFirstVersionSwapPrice: "															<< subFirstVersionSwapPrice											<<	endl;
//	o	<<	"Difference between subMonteCarloPrice and subOrdinaryGeneticVanillaSwapPrice: "		<<	subMonteCarloPrice-subOrdinaryGeneticVanillaSwapPrice			<<	endl;	
//	o	<<	"Difference between subOrdinaryGeneticVanillaSwapPrice and subFirstVersionSwapPrice: "	<<	subOrdinaryGeneticVanillaSwapPrice-subFirstVersionSwapPrice		<<	endl;	
//
//	o.close();
//}
