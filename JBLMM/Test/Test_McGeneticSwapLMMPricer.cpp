#include "JBLMM/Test/JBTests.h"

#include <iostream> 
#include <fstream> 
#include <string> 

#include <LMM/RNGenerator/McGenerator.h>
#include <LMM/RNGenerator/RNGenerator.h>

#include <LMM/LmmModel/Correlation.h>
#include <LMM/LmmModel/McTerminalLmm.h>

#include <LMM/LmmModel/Lmm.h>
#include <LMM/LmmModel/ConstShifted_HGVolatilityFunction.h>

#include <JBLMM/Pricer/McGeneticSwapLMMPricer.h>
#include <JBLMM/Instrument/InstrumentFactory.h>
#include <JBLMM/Pricer/GeneticVanillaSwapPricer.h>


void Test_McGeneticSwapLMMPricer()
{

	//! Parameters
	double	strike			=	0.02;
	LMM::Index	indexStart	=	0;
	LMM::Index	indexEnd	=	20;
	Tenor	floatingTenor	=	Tenor::_6M;
	Tenor	fixedTenor		=	Tenor::_12M;
	Tenor	tenorStruture	=	Tenor::_6M;
	size_t	horizonYear		=	10;
	LMMTenorStructure_PTR lmmTenorStructure( new LMMTenorStructure(tenorStruture, horizonYear));

	cout	<<	"strike:                        "	<<	strike												<<	endl;
	cout	<<	"indexStart:                    "	<<	indexStart											<<	endl;
	cout	<<	"indexEnd:                      "	<<	indexEnd											<<	endl;
	cout	<<	"tenorStrutureYearFraction:     "	<<	lmmTenorStructure->get_tenorType().YearFraction()	<<	endl;
	cout	<<	"floatingVStenorStrutureRatio:  "	<<	floatingTenor.ratioTo(tenorStruture)				<<	endl;
	cout	<<	"fixedVStenorStrutureRatio:     "	<<	fixedTenor.ratioTo(tenorStruture)					<<	endl;

	double	fwdRate		=	0.02;
	std::vector<double> liborsInitValue(lmmTenorStructure->get_horizon()+1, fwdRate);
	cout	<<	"myInitialLibor:  ";
	for (size_t i = 0; i <liborsInitValue.size(); i++)
	{
		cout	<<	liborsInitValue[i]	<<	" ";
	}
	cout	<<	 endl;
	cout	<<	 endl;


	//VanillaSwap_Chi_Trang
	VanillaSwap firstVersionVanillaSwap(strike, indexStart , indexEnd, floatingTenor, fixedTenor, lmmTenorStructure);
	LmmVanillaSwapPricer myVSP(lmmTenorStructure);
	double FirstVersionSwapPrice	=	myVSP.swapNPV_Analytical_1(firstVersionVanillaSwap, liborsInitValue);

	//---------------------------Build Lmm and McLmm's structure--------------------------------------
	//! Parameter of h
	double a = -0.06;
	double b = 0.17;
	double c = 0.54;
	double d = 0.17;
	Shifted_HGVolatilityParam::ABCDParameter abcdParam (a,b,c,d);
	//Parameter of hg
	double g_constParam = 1.0;
	double shift_constParam = -0.01;
	ConstShifted_HGVolatilityParam_PTR hgParam( new ConstShifted_HGVolatilityParam(lmmTenorStructure,abcdParam,g_constParam,shift_constParam));
	
	//! Correlation 1
	size_t nbFactor       = 3; // need to test nbFactor  = 3, and nbFactor = 
	size_t correlFullRank = lmmTenorStructure->get_horizon()+1;
	size_t correlReducedRank = nbFactor;
	//!"Check Parameters(): Condition not implemented yet."
	std::cout << "checkParams(): ";
	CorrelationReductionType::CorrelationReductionType correlReductionType = CorrelationReductionType::PCA;
	double correlAlpha = 0.0;
	double correlBeta  = 0.1;
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
	McLmm_PTR mcLmm(new McTerminalLmm(shiftedLmm, liborsInitValue, rnGenerator, MCSchemeType::EULER));

	//build a McGeneticSwapLMMPricer
	McGeneticSwapLMMPricer_PTR mcGeneticSwapLMMPricer(new McGeneticSwapLMMPricer(mcLmm));

	//build the geneticVanillaSwap
	GeneticSwap_CONSTPTR vanillaSwap_Genetic=InstrumentFactory::createVanillaSwap(
			strike,indexStart,indexEnd,floatingTenor,fixedTenor,lmmTenorStructure,1.0);

	//use Monte Carlo Method
	size_t nbSimulation=1000;
	double	MonteCarloPrice		=	mcGeneticSwapLMMPricer->swapNPV(vanillaSwap_Genetic, nbSimulation);

	//ordinaryGeneticVanillaSwapPricer
	GeneticVanillaSwapPricer_PTR geneticVanillaSwapPricer(new GeneticVanillaSwapPricer());
	double	OrdinaryGeneticVanillaSwapPrice		=	geneticVanillaSwapPricer->geneticVanillaSwap_Analytical(vanillaSwap_Genetic,liborsInitValue);

	//subVanillaSwap
	LMM::Index	subIndexStart	=	10;
	LMM::Index	subIndexEnd		=	16;
	GeneticSwap_CONSTPTR subVanillaSwap_Genetic=InstrumentFactory::createVanillaSwap(
			strike,subIndexStart,subIndexEnd,floatingTenor,fixedTenor,lmmTenorStructure,1.0);
	double	subMonteCarloPrice	=	mcGeneticSwapLMMPricer->swapNPV(subVanillaSwap_Genetic, nbSimulation);
	//subOrdinaryGeneticVanillaSwapPrice
	double	subOrdinaryGeneticVanillaSwapPrice		=	geneticVanillaSwapPricer->geneticVanillaSwap_Analytical(subVanillaSwap_Genetic,liborsInitValue);
	//subFirstVersionVanillaSwapPrice
	VanillaSwap subFirstVersionVanillaSwap(strike, subIndexStart , subIndexEnd, floatingTenor, fixedTenor, lmmTenorStructure);
	double subFirstVersionSwapPrice		=	myVSP.swapNPV_Analytical_1(subFirstVersionVanillaSwap, liborsInitValue);


	cout	<<	"MonteCarloPrice: "																		<<	MonteCarloPrice													<<	endl;
	cout	<<	"OrdinaryGeneticVanillaSwapPrice: "														<<	OrdinaryGeneticVanillaSwapPrice									<<	endl;
	cout	<< "FirstVersionSwapPrice: "																<< FirstVersionSwapPrice											<< endl;
	cout	<<	"Difference between MonteCarloPrice and OrdinaryGeneticVanillaSwapPrice: "				<<	MonteCarloPrice-OrdinaryGeneticVanillaSwapPrice					<<	endl;	
	cout	<<	"Difference between OrdinaryGeneticVanillaSwapPrice and FirstVersionSwapPrice: "		<<	OrdinaryGeneticVanillaSwapPrice-FirstVersionSwapPrice			<<	endl;	
	cout	<<	"subMonteCarloPrice: "																	<<	subMonteCarloPrice												<<	endl;
	cout	<<	"subOrdinaryGeneticVanillaSwapPrice: "													<<	subOrdinaryGeneticVanillaSwapPrice								<<	endl;
	cout	<< "subFirstVersionSwapPrice: "																<< subFirstVersionSwapPrice											<< endl;
	cout	<<	"Difference between subMonteCarloPrice and subOrdinaryGeneticVanillaSwapPrice: "		<<	subMonteCarloPrice-subOrdinaryGeneticVanillaSwapPrice			<<	endl;	
	cout	<<	"Difference between subOrdinaryGeneticVanillaSwapPrice and subFirstVersionSwapPrice: "	<<	subOrdinaryGeneticVanillaSwapPrice-subFirstVersionSwapPrice		<<	endl;

	ofstream o;
	o.open("TestResult_GeneticVanillaSwap_05_06.csv",  ios::out | ios::app );
	o	<<	endl;
	o	<<	endl;
	o	<<	endl;
	o	<<	"strike: "									<<	strike														<<	endl;
	o	<<	"indexStart: "								<<	indexStart													<<	endl;
	o	<<	"indexEnd: "								<<	indexEnd													<<	endl;
	o	<<	"floatingTenorVSLmmStructureTenorRatio: "	<<	floatingTenor.ratioTo(lmmTenorStructure->get_tenorType())	<<	endl;
	o	<<	"fixedTenorVSLmmStructureTenorRatio: "		<<	fixedTenor.ratioTo(lmmTenorStructure->get_tenorType())		<<	endl;
	o	<<	"floatingTenorYearFraction: "				<<	floatingTenor.YearFraction()								<<	endl;
	o	<<	"fixedTenorYearFraction: "					<<	fixedTenor.YearFraction()									<<	endl;
	o	<<	"horizonYear: "								<<	horizonYear													<<	endl;
	o	<<	"liborsInitValue: ";
	for(size_t i=0; i<liborsInitValue.size(); i++)
	{
		o	<<	liborsInitValue[i]	<<	" ";
	}
	o	<<	endl;
	o	<<	"nbSimulation: "																		<<	nbSimulation													<<	endl;
	o	<<	"PRICES: "	<<	endl;
	o	<<	"MonteCarloPrice: "																		<<	MonteCarloPrice													<<	endl;
	o	<<	"OrdinaryGeneticVanillaSwapPrice: "														<<	OrdinaryGeneticVanillaSwapPrice									<<	endl;
	o	<<	"FirstVersionSwapPrice: "																<< FirstVersionSwapPrice											<<	endl;
	o	<<	"Difference between MonteCarloPrice and OrdinaryGeneticVanillaSwapPrice: "				<<	MonteCarloPrice-OrdinaryGeneticVanillaSwapPrice					<<	endl;	
	o	<<	"Difference between OrdinaryGeneticVanillaSwapPrice and FirstVersionSwapPrice: "		<<	OrdinaryGeneticVanillaSwapPrice-FirstVersionSwapPrice			<<	endl;
	o	<<	"SUBPRICES: "<<	endl;
	o	<<	"subIndexStart: "																		<<	subIndexStart													<<	endl;
	o	<<	"subIndexEnd: "																			<<	subIndexEnd														<<	endl;
	o	<<	"subMonteCarloPrice: "																	<<	subMonteCarloPrice												<<	endl;
	o	<<	"subOrdinaryGeneticVanillaSwapPrice: "													<<	subOrdinaryGeneticVanillaSwapPrice								<<	endl;
	o	<<	"subFirstVersionSwapPrice: "															<< subFirstVersionSwapPrice											<<	endl;
	o	<<	"Difference between subMonteCarloPrice and subOrdinaryGeneticVanillaSwapPrice: "		<<	subMonteCarloPrice-subOrdinaryGeneticVanillaSwapPrice			<<	endl;	
	o	<<	"Difference between subOrdinaryGeneticVanillaSwapPrice and subFirstVersionSwapPrice: "	<<	subOrdinaryGeneticVanillaSwapPrice-subFirstVersionSwapPrice		<<	endl;	

	o.close();
}
