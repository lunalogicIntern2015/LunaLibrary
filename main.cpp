#include <LMM/Test/Tests.h>
#include <stdio.h>



#include <Cheyette/unit_test/TestFonction.h>
#include <Cheyette/unit_test/Test_CheyetteDD_Model.h>
#include <Cheyette/unit_test/TestIntegrator1D_2D.h>
#include <Cheyette/unit_test/TestMCGeneric.h>
#include <Cheyette/unit_test/TestApproxDD.h>
#include <Cheyette/CheyetteModel/CheyetteDD_Model.h>
#include <LMM/numeric/Integrator1D.h>
#include <Cheyette/unit_test/TestMC.h>
#include <Cheyette/unit_test/TestMC_vs_approx.h>



int main()
{
	//testVanillaSwap() ;
	//testgenericSwap() ;
	//testgenericSwaption() ;

//k, sigma, m
	std::vector<double> x, m_y, sigma_y_low, sigma_y_high ;
	x.push_back(0.) ; x.push_back(10.) ; x.push_back(20.) ; 
	double m_param = 0. ;
	double k = 0.2 ;
	double sigm_low = 0.05 ;
	double sigm_high = 0.2 ;
	m_y.push_back(m_param) ; m_y.push_back(m_param) ;
	sigma_y_low.push_back(sigm_low) ; sigma_y_low.push_back(2. * sigm_low) ;
	sigma_y_high.push_back(sigm_high) ; sigma_y_high.push_back(2. * sigm_high) ;
	int curveChoice = 2 ;
//simu MC	
	std::vector<size_t> nbSimus(3) ;
	nbSimus[0] = 10000 ; nbSimus[1] = 40000 ; nbSimus[2] = 100000 ; //nbSimus[3] = 100000 ;// nbSimus[4] = 200000 ;

	size_t a = 1 ; size_t b = 1 ;
//
	test_approx(0.008, a, b, Tenor::_6M, Tenor::_1YR, curveChoice, 1, x, m_y, sigma_y_low, k, nbSimus) ;
////	test_approx(0.008, a, b, Tenor::_6M, Tenor::_1YR, curveChoice, 1, x, m_y, sigma_y_high, k, nbSimus) ;
////	test_approx_ATM(a, b, Tenor::_6M, Tenor::_1YR, curveChoice, 1, x, m_y, sigma_y_low, k, nbSimus) ;
////	test_approx_ATM(a, b, Tenor::_6M, Tenor::_1YR, curveChoice, 1, x, m_y, sigma_y_high, k, nbSimus) ;
//////
////	m_param = 0.5 ;
////	m_y.clear() ; m_y.push_back(m_param) ; m_y.push_back(m_param) ;
////	test_approx(0.008, a, b, Tenor::_6M, Tenor::_1YR, curveChoice, 1, x, m_y, sigma_y_low, k, nbSimus) ;
////	test_approx(0.008, a, b, Tenor::_6M, Tenor::_1YR, curveChoice, 1, x, m_y, sigma_y_high, k, nbSimus) ;
////	test_approx_ATM(a, b, Tenor::_6M, Tenor::_1YR, curveChoice, 1, x, m_y, sigma_y_low, k, nbSimus) ;
////	test_approx_ATM(a, b, Tenor::_6M, Tenor::_1YR, curveChoice, 1, x, m_y, sigma_y_high, k, nbSimus) ;
//////
////	m_param = 1. ;
////	m_y.clear() ; m_y.push_back(m_param) ; m_y.push_back(m_param) ;
////	test_approx(0.008, a, b, Tenor::_6M, Tenor::_1YR, curveChoice, 1, x, m_y, sigma_y_low, k, nbSimus) ;
////	test_approx(0.008, a, b, Tenor::_6M, Tenor::_1YR, curveChoice, 1, x, m_y, sigma_y_high, k, nbSimus) ;
////	test_approx_ATM(a, b, Tenor::_6M, Tenor::_1YR, curveChoice, 1, x, m_y, sigma_y_low, k, nbSimus) ;
////	test_approx_ATM(a, b, Tenor::_6M, Tenor::_1YR, curveChoice, 1, x, m_y, sigma_y_high, k, nbSimus) ;
//
////5Y 5Y
//	a = 5 ; b = 5 ; 
//	m_param = 0. ;
//	m_y.clear() ; m_y.push_back(m_param) ; m_y.push_back(m_param) ;
//	test_approx(0.008, a, b, Tenor::_6M, Tenor::_1YR, curveChoice, 1, x, m_y, sigma_y_low, k, nbSimus) ;
	//test_approx(0.008, a, b, Tenor::_6M, Tenor::_1YR, curveChoice, 1, x, m_y, sigma_y_high, k, nbSimus) ;
//	test_approx_ATM(a, b, Tenor::_6M, Tenor::_1YR, curveChoice, 1, x, m_y, sigma_y_low, k, nbSimus) ;
	//test_approx_ATM(a, b, Tenor::_6M, Tenor::_1YR, curveChoice, 1, x, m_y, sigma_y_high, k, nbSimus) ;
////
	//m_param = 0.5 ;
	//m_y.clear() ; m_y.push_back(m_param) ; m_y.push_back(m_param) ;
	//test_approx(0.008, a, b, Tenor::_6M, Tenor::_1YR, curveChoice, 1, x, m_y, sigma_y_low, k, nbSimus) ;
	//test_approx(0.008, a, b, Tenor::_6M, Tenor::_1YR, curveChoice, 1, x, m_y, sigma_y_high, k, nbSimus) ;
	//test_approx_ATM(a, b, Tenor::_6M, Tenor::_1YR, curveChoice, 1, x, m_y, sigma_y_low, k, nbSimus) ;
	//test_approx_ATM(a, b, Tenor::_6M, Tenor::_1YR, curveChoice, 1, x, m_y, sigma_y_high, k, nbSimus) ;
////
	//m_param = 1. ;
	//m_y.clear() ; m_y.push_back(m_param) ; m_y.push_back(m_param) ;
	//test_approx(0.008, a, b, Tenor::_6M, Tenor::_1YR, curveChoice, 1, x, m_y, sigma_y_low, k, nbSimus) ;
	////test_approx(0.008, a, b, Tenor::_6M, Tenor::_1YR, curveChoice, 1, x, m_y, sigma_y_high, k, nbSimus) ;
	//test_approx_ATM(a, b, Tenor::_6M, Tenor::_1YR, curveChoice, 1, x, m_y, sigma_y_low, k, nbSimus) ;
	//test_approx_ATM(a, b, Tenor::_6M, Tenor::_1YR, curveChoice, 1, x, m_y, sigma_y_high, k, nbSimus) ;

//10Y 10Y
	//a = 10 ; b = 10 ; 
	//m_param = 0. ;
	//m_y.clear() ; m_y.push_back(m_param) ; m_y.push_back(m_param) ;
	//test_approx(0.008, a, b, Tenor::_6M, Tenor::_1YR, curveChoice, 1, x, m_y, sigma_y_low, k, nbSimus) ;
	//test_approx(0.008, a, b, Tenor::_6M, Tenor::_1YR, curveChoice, 1, x, m_y, sigma_y_high, k, nbSimus) ;
	//test_approx_ATM(a, b, Tenor::_6M, Tenor::_1YR, curveChoice, 1, x, m_y, sigma_y_low, k, nbSimus) ;
	//test_approx_ATM(a, b, Tenor::_6M, Tenor::_1YR, curveChoice, 1, x, m_y, sigma_y_high, k, nbSimus) ;
//
//	m_param = 0.5 ;
//	m_y.clear() ; m_y.push_back(m_param) ; m_y.push_back(m_param) ;
	//test_approx(0.008, a, b, Tenor::_6M, Tenor::_1YR, curveChoice, 1, x, m_y, sigma_y_low, k, nbSimus) ;
	//test_approx(0.008, a, b, Tenor::_6M, Tenor::_1YR, curveChoice, 1, x, m_y, sigma_y_high, k, nbSimus) ;
	//test_approx_ATM(a, b, Tenor::_6M, Tenor::_1YR, curveChoice, 1, x, m_y, sigma_y_low, k, nbSimus) ;
	//test_approx_ATM(a, b, Tenor::_6M, Tenor::_1YR, curveChoice, 1, x, m_y, sigma_y_high, k, nbSimus) ;
////
	//m_param = 1. ;
	//m_y.clear() ; m_y.push_back(m_param) ; m_y.push_back(m_param) ;
	//test_approx(0.008, a, b, Tenor::_6M, Tenor::_1YR, curveChoice, 1, x, m_y, sigma_y_low, k, nbSimus) ;
	//test_approx(0.008, a, b, Tenor::_6M, Tenor::_1YR, curveChoice, 1, x, m_y, sigma_y_high, k, nbSimus) ;
//	test_approx_ATM(a, b, Tenor::_6M, Tenor::_1YR, curveChoice, 1, x, m_y, sigma_y_low, k, nbSimus) ;
//	test_approx_ATM(a, b, Tenor::_6M, Tenor::_1YR, curveChoice, 1, x, m_y, sigma_y_high, k, nbSimus) ;

	//test black implied vol
	//returns vol sigma
	//double bs_call_price = 7.97 ;
	//double fwd = 100 ;
	//double strike = 100 ;
	//double T = 1 ;
	//std::cout << NumericalMethods::Black_impliedVolatility(bs_call_price, fwd, strike, T) << std::endl;
	 
	//double bs_call_price = 5.4529/100. ;
	//double annuity0 = 4.49911 ; //A0
	//double fwd = 0.0201199 ; //s0
	//double strike = 0.008 ;   
	//double T = 5.0 ;
	//std::cout << NumericalMethods::Black_SwaptionImpliedVolatility(bs_call_price, annuity0, fwd, strike, T) ;

	//testgenericSwap() ; 
	//testgenericSwaption() ;

	getchar();
}