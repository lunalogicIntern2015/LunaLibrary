#include <LMM/Test/Test_Junbin/JBTests.h>

#include <iostream> 
#include <fstream> 
#include <string> 

#include <RNGenerator/McGenerator.h>
#include <RNGenerator/RNGenerator.h>

#include <LMM/Model/Correlation.h>
#include <LMM/Mc/McTerminalLmm.h>

#include <LMM/Model/Lmm.h>
#include <LMM/Model/ConstShifted_HGVolatilityFunction.h>

#include <LMM/Pricer/McLmmPricer/McLmmGenericSwapPricer.h>
#include <Instrument/InstrumentFactory.h>
#include <LMM/Pricer/GenericVanillaSwapPricer.h>
#include <LMM/Pricer/McLmmPricer/McLmmGenericTargetSwapPricer.h>
#include <Instrument/Coupon/TargetCoupon.h>


void Test_GeneticTargetSwapOneTrajectory()
{
	//Only one trajectory
	const size_t nbSimulation	=	1;

	//! Parameters
	double	strike			=	0.02;
	LMM::Index	indexStart	=	0;
	LMM::Index	indexEnd	=	4;
	Tenor	floatingTenor	=	Tenor::_6M;
	Tenor	fixedTenor		=	Tenor::_12M;
	Tenor	tenorStruture		=	Tenor::_6M;
	size_t	horizonYear		=	2;
	LMMTenorStructure_PTR lmmTenorStructure( new LMMTenorStructure(tenorStruture, horizonYear));
	double	nominal			=	1.0;
	double	target			=	0.019;

	double	fwdRate		=	0.02;
	std::vector<double> liborsInitValue(lmmTenorStructure->get_horizon()+1, fwdRate);

	ofstream out;
	out.open("libormatrix and numeraire.csv", ios::out | ios::app);
	out	<<	endl;
	out	<<	endl;
	out	<<	endl;
	out	<<	"strike:                        "	<<	strike												<<	endl;
	out	<<	"indexStart:                    "	<<	indexStart											<<	endl;
	out	<<	"indexEnd:                      "	<<	indexEnd											<<	endl;
	out	<<	"tenorStrutureYearFraction:     "	<<	lmmTenorStructure->get_tenorType().YearFraction()	<<	endl;
	out	<<	"floatingVStenorStrutureRatio:  "	<<	floatingTenor.ratioTo(tenorStruture)				<<	endl;
	out	<<	"fixedVStenorStrutureRatio:     "	<<	fixedTenor.ratioTo(tenorStruture)					<<	endl;
	out	<<	"target:                        "	<<	target												<<	endl;

	out	<<	"myInitialLibor:  ";
	for (size_t i = 0; i <liborsInitValue.size(); i++)
	{
		out	<<	liborsInitValue[i]	<<	" ";
	}
	out	<<	 endl;


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

	Lmm_PTR shiftedLmm (new Lmm(dispersion));

	//build McLmmStructure
	McLmm_PTR mcLmm(new McTerminalLmm(shiftedLmm, liborsInitValue, rnGenerator, MCSchemeType::EULER));
	//build a McGeneticSwapLMMPricer
	McLmmGenericTargetSwapPricer_PTR mcGeneticTargetSwapLMMPricer(new McLmmGenericTargetSwapPricer(mcLmm));
	
	//build the GvanillaSwap
	GenericSwap_CONSTPTR targetSwap_Genetic=InstrumentFactory::createStandardTARNSwap(
			strike, indexStart, indexEnd, floatingTenor, fixedTenor, lmmTenorStructure, nominal, target);

	//------------------------------------test clone copy----------------------------------------------
	////print a targetCoupon
	//Coupon_CONSTPTR coupon = targetSwap_Genetic->getLeg2()->getLeg()[0];
	//TargetCoupon_CONSTPTR targetCoupon	=	boost::dynamic_pointer_cast<const TargetCoupon>(coupon);
	//targetCoupon->show();

	////test a copy
	//Coupon_PTR targetCouponCopy	=	targetCoupon->clone();
	//targetCouponCopy->show();


	//use Monte Carlo Method
	double	MonteCarloTargetPrice		=	mcGeneticTargetSwapLMMPricer->swapNPV(targetSwap_Genetic, nbSimulation);

	//print liborMatrix and array
	out	<<	"liborMatrix:"	<<	endl;
	for(LMM::Index i=0; i<mcLmm->get_lmm()->get_horizon()+1; i++)
	{
		for(LMM::Index j=0; j<mcLmm->get_lmm()->get_horizon()+1; j++)
		{
			out <<	mcLmm->get_liborMatrix()(i,j) << " ";
		}
		out << endl;
	}
	out	<<	"numeraire: ";
	for(LMM::Index i=0; i<mcLmm->get_lmm()->get_horizon()+1; i++)
	{
		out <<	mcLmm->get_numeraire(i) << " ";
	}
	out << endl;

	out	<<	"MonteCarloTargetPrice: "	<<	MonteCarloTargetPrice	<<	endl;

	out.close();
}
