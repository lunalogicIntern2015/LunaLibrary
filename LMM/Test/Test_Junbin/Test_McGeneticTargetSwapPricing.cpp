#include <JBLMM/Test/JBTests.h>

#include <iostream> 
#include <fstream> 
#include <string> 

#include <RNGenerator/McGenerator.h>
#include <RNGenerator/RNGenerator.h>

#include <LMM/Model/Correlation.h>
#include <LMM/Mc/McTerminalLmm.h>

#include <LMM/Model/Lmm.h>
#include <LMM/Model/ConstShifted_HGVolatilityFunction.h>

#include <JBLMM/Pricer/McLmmGenericSwapPricer.h>
#include <JBInstrument/InstrumentFactory.h>
#include <JBLMM/Pricer/GenericVanillaSwapPricer.h>
#include <LMM/Pricer/McLmmPricer/McLmmGenericTargetSwapPricer.h>


void Test_McGeneticTargetSwapLMMPricing()
{
	//! Parameter
	double	strike			=	0.02;
	LMM::Index	indexStart	=	0;
	LMM::Index	indexEnd	=	4;
	Tenor	floatingTenor	=	Tenor::_6M;
	Tenor	fixedTenor		=	Tenor::_12M;
	Tenor	tenorType		=	Tenor::_6M;
	size_t	horizonYear		=	2;
	LMMTenorStructure_PTR lmmTenorStructure( new LMMTenorStructure(tenorType, horizonYear));
	double	nominal			=	1.0;
	double	target			=	0.02;

	double	fwdRate		=	0.02;
	std::vector<double> liborsInitValue(lmmTenorStructure->get_horizon()+1, fwdRate);


	//! volatility function
	double a = -0.06;
	double b = 0.17;
	double c = 0.54;
	double d = 0.17;
	Shifted_HGVolatilityParam::ABCDParameter abcdParam (a,b,c,d);

	double g_constParam = 1.0;
	double shift_constParam = -0.01;

	ConstShifted_HGVolatilityParam_PTR hgParam( new ConstShifted_HGVolatilityParam(lmmTenorStructure,abcdParam,g_constParam,shift_constParam));
	
	//! Correlation 1
	size_t nbFactor       = 3; // need to test nbFactor  = 3, and nbFactor = 
	size_t correlFullRank = lmmTenorStructure->get_horizon()+1;
	size_t correlReducedRank = nbFactor;
	CorrelationReductionType::CorrelationReductionType correlReductionType = CorrelationReductionType::PCA;
	double correlAlpha = 0.0;
	double correlBeta  = 0.1;
	Correlation_PTR correlation(new XY_beta_Correlation(correlFullRank,correlReducedRank, correlReductionType,correlAlpha,correlBeta));
	correlation->calculate(); // for print.
	correlation->print("test_McTerminalLmm_Correlation.csv");

	ConstShifted_HGVolatilityFunction_PTR hgVolatilityFunction (new ConstShifted_HGVolatilityFunction(lmmTenorStructure, correlation, hgParam)); 
	hgVolatilityFunction->print("test_McTerminalLmm_Volatility.csv");
	//! Dispersion
	Dispersion dispersion(hgVolatilityFunction);

	unsigned long seed = 503333;
	RNGenerator_PTR  rnGenerator(new McGenerator(seed));

	//build lmm and mcLmm model
	Lmm_PTR shiftedLmm (new Lmm(dispersion));
	McLmm_PTR mcLmm(new McTerminalLmm(shiftedLmm, liborsInitValue, rnGenerator, MCSchemeType::EULER));

	//build a McGeneticSwapLMMPricer
	McLmmGenericSwapPricer_PTR mcGeneticTargetSwapLMMPricer(new McLmmGenericTargetSwapPricer(mcLmm));

	//build the vanillaSwap in the way of GeneticSwap
	GenericSwap_CONSTPTR targetSwap_Genetic=InstrumentFactory::createStandardTARNSwap(
			strike, indexStart, indexEnd, floatingTenor, fixedTenor, lmmTenorStructure, nominal, target);

	//use Monte Carlo Method
	size_t	nbSimulation	=	100;
	double	MonteCarloTargetPrice		=	mcGeneticTargetSwapLMMPricer->swapNPV(targetSwap_Genetic, nbSimulation);

	cout	<<	"MonteCarloTargetPrice: "	<<	MonteCarloTargetPrice	<<	endl;

}
