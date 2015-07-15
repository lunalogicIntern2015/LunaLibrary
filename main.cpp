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
#include <Cheyette/unit_test/TestCalibrator.h>

#include <LMM\Test\Tests.h>



int main()
{
	//file numero 0, coterminal 2Y
	size_t numFile = 1 ;
	size_t coterminal = 10 ;
	testCalib2(numFile, coterminal) ;   //2004 04 07 : file 27
//	test_CheyetteCalibrationMarketData() ;
	//test_LmmCalibrationMarketData() ;
		//testCalib() ;

	//double bs_call_price = 0.0747381/100 ; //0.0611958/100 ;
	//double annuity0 = 83.86/100  ; //A0
	//double s0 = 2.6342/100 ; //s0
	//double strike = s0 ;   
	//double T = 10.0 ;
	//std::cout << NumericalMethods::Black_SwaptionImpliedVolatility(bs_call_price, annuity0, s0, strike, T) ;

//	size_t dimMatrix = 2 ;
//	std::vector<size_t> expiry_maturity(dimMatrix) ;
//	//expiry_maturity[0] = 1 ; expiry_maturity[1] = 3 ; 
//	expiry_maturity[0] = 0 ; expiry_maturity[1] = 1 ; 
//	//expiry_maturity[2] = 6 ; expiry_maturity[3] = 7 ; expiry_maturity[4] = 10 ;
//	//expiry_maturity[5] = 12 ; //expiry_maturity[7] = 15 ; 
//
////k, sigma, m
//	std::vector<double> x, m_y, sigma_y_low, sigma_y_high ;
//	x.push_back(0.) ; x.push_back(2.) ; //x.push_back(20.) ; x.push_back(30.) ; 
//	double m_param = 0.2 ;
//	double k = 0.2 ;
//	double sigm_low = 5.0/100 ;
//	double sigm_high = 80.0/100 ;
//	m_y.push_back(m_param) ; //m_y.push_back(m_param) ; m_y.push_back(m_param) ;
//	sigma_y_low.push_back(sigm_low) ; //sigma_y_low.push_back(1.3 * sigm_low) ; sigma_y_low.push_back(1.5 * sigm_low) ;
//	sigma_y_high.push_back(sigm_high) ;// sigma_y_high.push_back(1.3 * sigm_high) ; sigma_y_high.push_back(1.5 * sigm_high) ;
//	int curveChoice = 4 ;
//	int shiftChoice = 1 ;
////simu MC	
//	std::vector<size_t> nbSimus(1) ;
//	nbSimus[0] = 5 ; // nbSimus[1] = 400 ; nbSimus[2] = 100 ; //nbSimus[3] = 100000 ;// nbSimus[4] = 200000 ;
//
//	test_approx_ATM(1, 1, Tenor::_6M, Tenor::_1YR, 
//					 curveChoice, shiftChoice, 
//					  x, m_y, sigma_y_low, k, nbSimus);

	
//	//std::cout << " ---------------------------------------------- " << std::endl ;
//
//	//	test_approx_ATM(10, 1, Tenor::_6M, Tenor::_1YR, 
//	//				 curveChoice, shiftChoice, 
//	//				  x, m_y, sigma_y_high, k, nbSimus);
//
////smile
//	double strikeATM = 1.28581/100 ;
//	smile(3, 3, strikeATM, Tenor::_6M, Tenor::_1YR, curveChoice, shiftChoice, x, m_y, sigma_y_high, k) ;
//
//	strikeATM = 2.01199/100 ;
//	smile(5, 5, strikeATM, Tenor::_6M, Tenor::_1YR, curveChoice, shiftChoice, x, m_y, sigma_y_high, k) ;

//
////ATM
	//m_param = 0. ;
	//m_y.clear() ; m_y.push_back(m_param) ; m_y.push_back(m_param) ; m_y.push_back(m_param) ;
	//test_printMatrix_ATM(expiry_maturity, Tenor::_6M, Tenor::_1YR, curveChoice, shiftChoice, x, m_y, sigma_y_low, k, nbSimus) ;
	//test_printMatrix_ATM(expiry_maturity, Tenor::_6M, Tenor::_1YR, curveChoice, shiftChoice, x, m_y, sigma_y_high, k, nbSimus) ;
	//
	//m_param = 0.5 ;
	//m_y.clear() ; m_y.push_back(m_param) ; m_y.push_back(m_param) ; m_y.push_back(m_param) ;
	////test_printMatrix_ATM(expiry_maturity, Tenor::_6M, Tenor::_1YR, curveChoice, shiftChoice, x, m_y, sigma_y_low, k, nbSimus) ;
	//test_printMatrix_ATM(expiry_maturity, Tenor::_6M, Tenor::_1YR, curveChoice, shiftChoice, x, m_y, sigma_y_high, k, nbSimus) ;

	////m_param = 1. ;
	//m_y.clear() ; m_y.push_back(m_param) ; m_y.push_back(m_param) ; m_y.push_back(m_param) ;
	//test_printMatrix_ATM(expiry_maturity, Tenor::_6M, Tenor::_1YR, curveChoice, shiftChoice, x, m_y, sigma_y_low, k, nbSimus) ;
	//test_printMatrix_ATM(expiry_maturity, Tenor::_6M, Tenor::_1YR, curveChoice, shiftChoice, x, m_y, sigma_y_high, k, nbSimus) ;

	//test_printMatrix(0.008, Tenor::_6M, Tenor::_1YR, curveChoice, shiftChoice, x, m_y, sigma_y_low, k, nbSimus) ;
	//m_param = 0.5 ;
	//m_y.clear() ; m_y.push_back(m_param) ; m_y.push_back(m_param) ;
	//test_printMatrix(0.008, Tenor::_6M, Tenor::_1YR, curveChoice, 1, x, m_y, sigma_y_low, k, nbSimus) ;
	//m_param = 1. ;
	//m_y.clear() ; m_y.push_back(m_param) ; m_y.push_back(m_param) ;
	//test_printMatrix(0.008, Tenor::_6M, Tenor::_1YR, curveChoice, 1, x, m_y, sigma_y_low, k, nbSimus) ;
	
	//test black implied vol
	//returns vol sigma
	//double bs_call_price = 7.97 ;
	//double fwd = 100 ;
	//double strike = 100 ;
	//double T = 1 ;
	//std::cout << NumericalMethods::Black_impliedVolatility(bs_call_price, fwd, strike, T) << std::endl; 



	//testgenericSwap() ; 
	//testgenericSwaption() ;

	getchar();
}