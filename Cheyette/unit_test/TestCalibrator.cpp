#include "TestCalibrator.h"



void testCalib()
{
//market Data
	std::vector<size_t> a_expiry(1) ; a_expiry[0] = 1 ; //sur la tenor structure ?
	std::vector<size_t> b_tenor(1) ; b_tenor[0] = 1 ;	
	std::vector<double> strikeATM(1) ; strikeATM[0] =  2.78/100 ;  //1.2341/100 ; 0.36/100 ;
	std::vector<double> volQuotes(1) ; volQuotes[0] =  29.41/100 ; //70.97/100 ;50.81/100 ;
	std::vector<double> skew(1) ; skew[0] =  -7.4 ; //-3.8 ;
	double shift = 5.0/10000 ;
	CourbeInput_PTR courbe_PTR_test(createCourbeInput(5));

	double strike = strikeATM[0] ;
	LMM::Index  indexStart = 2 ; 
	LMM::Index  indexEnd = 4 ; 
	Tenor floatingLegTenorType = Tenor::_6M ;
	Tenor fixedLegTenorType = Tenor::_1YR ;
	LMMTenorStructure_PTR lmmTenorStructure(new LMMTenorStructure(Tenor::_6M, 30)) ;

	VanillaSwap vanillaSwap(strike, indexStart, indexEnd, floatingLegTenorType, fixedLegTenorType, lmmTenorStructure); 

	VanillaSwaption_PTR swaptionPTR(new VanillaSwaption(vanillaSwap , OptionType::OptionType::CALL));
	std::vector<VanillaSwaption_PTR> vectSwaptions(1) ;
	vectSwaptions[0] = swaptionPTR ;
	MarketData_PTR mktData(new MarketData(a_expiry, b_tenor, strikeATM, volQuotes, skew, shift, vectSwaptions, "Black")) ;

	
//Cheyette DD model

	std::vector<double> x, sigma_y, m_y ;
	x.push_back(0) ; x.push_back(1) ; 

	double m_param = 0.4 ;
	double k = 0.2 ;
	double sigm = 0.20 ;

	m_y.push_back(m_param) ; 
	sigma_y.push_back(sigm) ;
	
	Piecewiseconst_RR_Function sigma(x, sigma_y) ; 
	Piecewiseconst_RR_Function m(x, m_y) ; 

	CheyetteDD_Model::CheyetteDD_Parameter monStruct(k, sigma, m) ;
	int shiftChoice = 1 ;
	CheyetteDD_Model_PTR modele_test_PTR(new CheyetteDD_Model(courbe_PTR_test, monStruct, shiftChoice)) ;
	modele_test_PTR->show() ;

	CheyetteDD_VanillaSwaptionApproxPricer_PTR approx(
					new CheyetteDD_VanillaSwaptionApproxPricer(modele_test_PTR, vectSwaptions[0]));

//calibrator
	QuantLib::Size maxIterations	= 500 ;
	QuantLib::Real rootEpsilon		= 1/10000 ;
	QuantLib::Real functionEpsilon	= 1/10000 ;

	size_t indexSwaption = 1 ;

	CheyetteBaseCostFunction_PTR cheyetteBaseCostFunction_PTR(
				new CheyetteDD_CostFunctionLevel(mktData, indexSwaption, approx)) ;

	Array sigmaInitiate(1) ; 
	sigmaInitiate[0] = sigm ;
	Array mInitiate(1) ;
	mInitiate[0] = m_param ;
	Array calibrated_sigma(1) ; 
	calibrated_sigma = sigmaInitiate ;
	Array calibrated_m(1) ;
	calibrated_m = mInitiate ;

	CheyetteDD_LocalCalibrator calibrator(	maxIterations,  rootEpsilon,   functionEpsilon,    
											cheyetteBaseCostFunction_PTR,
											sigmaInitiate, mInitiate,
											calibrated_sigma, calibrated_m) ; 

	calibrator.calibrate() ;

}