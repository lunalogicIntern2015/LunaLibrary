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

//	std::cout << NumericalMethods::Black_Price_vol2(100, 100, 0.2 * sqrt(2), 2) <<std::endl ;

	//test_P_t_T() ;
	//test_beginner();							//test_JB
	//vanillaSwapComparaisonExemple();			//test_JB		

	//Test_McGeneticSwapLMMPricer();			//test_JB
	//Test_McGeneticTargetSwapLMMPricing();		//test_JB

	//test_Integrator1D();
	//createDDModel_xyz() ;
	
	/************    tests JL   ************************/
	//TestFonctionConstanteMorceaux() ;
	//TestInterpolation_RR_Function() ;
	//TestCompositionMultiplicationFonction() ;

	//TestIntegrator1D_Riemann() ;
	//TestIncrementalIntegrator1D_Riemann() ;
	//TestIncrementalIntegrator2D_Riemann() ;

	//createSwap() ;
	//test_derivatives() ;
	//test_y_barre() ;
	//test_fonction_inverse() ;
	
	//test_swapRate_inverse() ;

	//test_ZC_swapRate_Num_Denom() ;  
	//createCourbeInput() ;
	
	//test_time_average() ; 
	//testgenericSwaption() ;

	test_approx() ;

	//test black implied vol
	//returns vol sigma
	/*double bs_call_price = 7.97 ;
	double fwd = 100 ;
	double strike = 100 ;
	double T = 1 ;*/
	//std::cout << NumericalMethods::Black_impliedVolatility(bs_call_price, fwd, strike, T) << std::endl;
	 
	//double bs_call_price = 0.00197137 ;
	//double annuity0 = 0.970445533548508 ; //A0
	//double fwd = 0.010050167084168 ; //s0
	//double strike = 0.9/100 ;
	//double T = 2.0 ;
	//std::cout << NumericalMethods::Black_SwaptionImpliedVolatility(bs_call_price, annuity0, fwd, strike, T) ;

	//testgenericSwap() ; 
	//testgenericSwaption() ;

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

	getchar();
}