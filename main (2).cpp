#include <stdio.h>


#include <LMM/Pricer/LmmApproximationPricer/LmmVanillaSwaptionApproxPricer_Rebonato.h>
#include <Numeric/NumericalMethods.h>
//#include <JBLMM/Test/JBTests.h>

#include <Cheyette/unit_test/TestFonction.h>
#include <Cheyette/unit_test/TestIntegrator1D_2D.h>
#include <Cheyette/unit_test/TestApproxDD.h>

#include <Cheyette/Model/CheyetteDD_Model.h>
#include <Numeric/Integrator1D.h>
#include <Cheyette/unit_test/TestMC.h>
#include <Cheyette/unit_test/TestCalibrator.h>
#include <LMM/Test/Tests.h>
#include <Cheyette/unit_test/TestApproxDD.h>


int main()
{

	//std::vector<double> nbSimulation_vect;
	//nbSimulation_vect.push_back(20);
	//nbSimulation_vect.push_back(30000);
	//nbSimulation_vect.push_back(20);
	//nbSimulation_vect.push_back(30000);
	//nbSimulation_vect.push_back(20);
	//nbSimulation_vect.push_back(30000);
	//nbSimulation_vect.push_back(20);
	//nbSimulation_vect.push_back(30000);

	//size_t nb=300000;
	//std::vector<std::vector<double>> basis;
	//clock_t startTime3 = clock();
	//basis.resize(nb);
	//clock_t midTime3 = clock();
	//for(size_t i=0; i<nb; i++)
	//	basis[i]=nbSimulation_vect;

	//clock_t endTime3 = clock();

	//clock_t resize_time = midTime3 - startTime3;
	//clock_t copy_time = endTime3 - midTime3;
	//clock_t time3 = endTime3 - startTime3;
	//double resize_in_second = resize_time/(double) CLOCKS_PER_SEC;
	//double copy_in_second = copy_time/(double) CLOCKS_PER_SEC;
	//double all_in_second = time3/(double) CLOCKS_PER_SEC;
	//cout << "resize_in_second " << resize_in_second << endl;
	//cout << "copy_time " << copy_in_second << endl;
	//cout << "all_in_second" << all_in_second << endl;



	//Test_LS_pricing_parameter();
	//Test_10_basis();
	//Test_evaluation_basis();
	//Test_LS_pricing_allSubSet_basis();
	//Test_Longstaff_Schwartz_CallableSwap();	// use GenericSwap
	//Test_pricing_forward_LS();
	//Test_regression_LS();						//test algo regression backward
	//Test_VanillaSwaptionPricing();			//vanillaSwaptionPricer for comparing
	//Test_Longstaff_Schwartz_CallableSwap();	// use GenericSwap
	//Test_Longstaff_Schwartz_CallableSwap_for_test();  
	//Test_Longstaff_Schwartz();				// use bermudan option
	//Test_RegressionLS();
	//Test_basis();
	//test_beginner();							//test_JB
	//vanillaSwapComparaisonExemple();			//test_JB		
	//Test_McGeneticSwapLMMPricer();				//test_JB
	//Test_McGeneticTargetSwapLMMPricing();		//test_JB
	//Test_GeneticTargetSwapOneTrajectory();	//test_JB
	//JB_test_LmmCalibrationMarketData();		//test_JB

	//test_Integrator1D();
	//createDDModel_xyz() ;
	
	/************    tests JL   ************************/
	//TestFonctionConstanteMorceaux() ;
	//TestInterpolation_RR_Function() ;
	//TestCompositionMultiplicationFonction() ;

	//TestIntegrator1D_Riemann() ;
	//TestIncrementalIntegrator1D_Riemann() ;

	//CheyetteDD_VanillaSwaptionApproxPricer_PTR test = createApproxPricer_PTR() ;

	//test_Derivative_ZC() ;  

	//TestIncrementalIntegrator2D_Riemann() ;    
	//std::cout << "pb lors de 2 appels consecutifs a integrale" << std::endl ;

	//createSwap() ;
	//createDDModel() ;
	//test_y_barre(0) ;
	//test_y_barre(0.5) ;
	//test_y_barre(1.0) ;
	//test_Derivative_ZC() ;
	
	//test_time_average() ; 
	//test_y_bar_cas_limite() ;
	//test_Integrator1D();

	//testSwap() ; 

	//test_incremental_integrale() ;

	//UneTrajectoireEuler() ;

	//TestMCSwapPricer() ; 

 



	/************  tests LMM   *************************/
	//test_Noise();  
	//test_HGVolatility();  
	//test_Functional(); 
	//test_BlackGreek(); 
	//test_McTerminalLmm(); 

	//test_VanillaFRAPricer();  
	//test_VanillaCapletPricer(); 
	//test_VanillaSwapPricer();  
	//test_VanillaSwaptionPricer();    
	//test_VanialSwaptionPricer_MCvsApprox(); 


	//test_SwaptionMarketDataContainer();		
	//test_CalibrationWithATMSwaptionVol();		
	//test_CalibrationShiftWithSwaptionVolSkew();
	
	//test_LmmVirtualCorrelationCalibration();  
	//test_LmmCalibrationConfig(); 
	//test_LmmGnuplotPrinterMatrix();  
	//test_Interpolation();  
	//test_GMatrixMapping(); 
	//test_UpperTriangleVanillaSwaptionQuotes(); 
	//test_LmmSwaptionMarketData(); 
	//test_LmmABCDFinder();	
	//test_LmmVirtualCalibration();	
	//test_LmmCalibrationSensitivity();	
	//test_LmmCalibrationMarketData();
	//test_LmmRegularizedCalibrationMarketData();  

	//test_CalibrationShiftWithSwaptionVolSkew();

	//MCforward_vs_annuity() ;
//test MC swap contre swap courbe + swaption OK
	//TestMCSwapPricer_annuity(1, 1, simus, floatingLegTenor, fixedLegTenor) ;
	//TestMCSwapPricer_annuity(1, 5, simus, floatingLegTenor, fixedLegTenor) ;

	//size_t a = 1 ;
	//size_t b = 9 ;
	//size_t simus = 30000 ;

	//Tenor floatingLegTenor	= Tenor::_6M ;
	//Tenor fixedLegTenor		= Tenor::_12M ;
	//MCforward(a, b, simus, floatingLegTenor, fixedLegTenor) ;
	//MCannuity(a, b, simus, floatingLegTenor, fixedLegTenor) ;

//calibration sur fichier et comparaison MC / approx
/*	size_t numfile = 24 ;
	size_t coterminal = 12 ;
	testCalib(numfile, coterminal) ;*/   //2004 04 07 : file 27


	//swap MC contre courbe
	//////TestMCSwapPricer(a, b, simus, floatingLegTenor, fixedLegTenor) ;
	//////TestMCSwapPricer_annuity(a, b, simus, floatingLegTenor, fixedLegTenor) ;

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