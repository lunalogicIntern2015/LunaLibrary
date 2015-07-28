#include <stdio.h>

#include <Cheyette/CheyetteModel/CheyetteDD_Model.h>
#include <LMM/numeric/Integrator1D.h>
#include <Cheyette/unit_test/TestMC.h>
#include <Cheyette/unit_test/TestCalibrator.h>
#include <LMM/Test/Tests.h>
#include <Cheyette/unit_test/TestApproxDD.h>


int main()
{

//test MC swap contre swap courbe + swaption OK
	//TestMCSwapPricer_annuity(1, 1, simus, floatingLegTenor, fixedLegTenor) ;
	//TestMCSwapPricer_annuity(1, 5, simus, floatingLegTenor, fixedLegTenor) ;

	size_t a = 10 ;
	size_t b = 10 ;
	size_t simus = 30000 ;

	Tenor floatingLegTenor	= Tenor::_6M ;
	Tenor fixedLegTenor		= Tenor::_12M ;

	//TestMCSwapPricer(a, b, simus, floatingLegTenor, fixedLegTenor) ;
	//TestMCSwapPricer_annuity(a, b, simus, floatingLegTenor, fixedLegTenor) ;

	MCforward(a, b, simus, floatingLegTenor, fixedLegTenor) ;
	MCannuity(a, b, simus, floatingLegTenor, fixedLegTenor) ; 

	//int curveChoice = 1 ; //courbe plate à 1%
	//int shiftChoice = 1 ;
	//std::vector<double> x(a+b+2) ;
	//std::vector<double> m_y(a+b+1) ;
	//std::vector<double> sigma_y(a+b+1) ;
	//for (size_t i = 0 ; i < a+b+2 ; ++i)
	//{
	//	x[i] = i ;
	//}
	//for (size_t i = 0 ; i < a+b+1 ; ++i)
	//{
	//	m_y[i] = 0. ;
	//	sigma_y[i] = 0.20 ;
	//}
	//double k = 0.02 ;
	//std::vector<size_t> nbSimus(1) ;
	//nbSimus[0] = 5000 ; nbSimus[1] = 10000 ; 
	//nbSimus[0] = 20000 ; 
	//nbSimus[3] = 50000 ; nbSimus[4] = 100000 ;

	//test_approx_ATM(a, b, floatingLegTenor, fixedLegTenor, 
	//				 curveChoice, shiftChoice, x, m_y, sigma_y, k, nbSimus) ;


	//size_t numFile = 6 ;
	//size_t coterminal = 14 ;
	//testCalib(numFile, coterminal) ;   //2004 04 07 : file 27

	//double bs_call_price = 0.0747381/100 ; //0.0611958/100 ;
	//double annuity0 = 83.86/100  ; //A0
	//double s0 = 2.6342/100 ; //s0
	//double strike = s0 ;   
	//double T = 10.0 ;
	//std::cout << NumericalMethods::Black_SwaptionImpliedVolatility(bs_call_price, annuity0, s0, strike, T) ;

	//size_t dimMatrix = 2 ;
	//std::vector<size_t> expiry_maturity(dimMatrix) ;
	//expiry_maturity[0] = 1 ; expiry_maturity[1] = 3 ; 
	//expiry_maturity[0] = 0 ; expiry_maturity[1] = 1 ; 
	//expiry_maturity[2] = 6 ; expiry_maturity[3] = 7 ; expiry_maturity[4] = 10 ;
	//expiry_maturity[5] = 12 ; //expiry_maturity[7] = 15 ; 

//k, sigma, m


//simu MC
	//std::vector<size_t> nbSimus(2) ;
	//nbSimus[0] = 5000 ; nbSimus[1] = 10000 ;
	////std::cout << " ---------------------------------------------- " << std::endl ;


	//	test_approx_ATM(10, 1, Tenor::_6M, Tenor::_1YR, 
	//				 curveChoice, shiftChoice, 
	//				  x, m_y, sigma_y_high, k, nbSimus);

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

	getchar();
}