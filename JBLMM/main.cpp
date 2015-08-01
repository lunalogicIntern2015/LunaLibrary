#include <LMM/Test/Tests.h>
#include <stdio.h>

#include <LMM/pricer/LmmVanillaSwaptionApproxPricer_Rebonato.h>
#include <LMM/numeric/NumericalMethods.h>
#include <JBLMM/Test/JBTests.h>

#include <Cheyette/unit_test/TestFonction.h>
#include <Cheyette/unit_test/TestIntegrator1D_2D.h>
#include <Cheyette/unit_test/TestApproxDD.h>
#include <Cheyette/CheyetteModel/CheyetteDD_Model.h>
#include <LMM/numeric/Integrator1D.h>
#include <Cheyette/unit_test/TestMC.h>



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



	Test_LS_pricing_parameter();
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

	getchar();
}