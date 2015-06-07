#include "TestMC.h"


void UneTrajectoireEuler()
{
	unsigned long seed = 47;
	RNGenerator_PTR  rnGenerator(new McGenerator(seed));

	std::vector<double> x, y ;
	x.push_back(0) ; x.push_back(0.5) ; x.push_back(1) ; 
	y.push_back(0.25) ; y.push_back(0.5) ;
	double k(1) ;

	std::vector<double> listeMatu, tauxZC ;
	listeMatu.push_back(0) ;	tauxZC.push_back(0.8/100) ; 
	listeMatu.push_back(1) ;	tauxZC.push_back(0.85/100) ; 
	listeMatu.push_back(2) ;	tauxZC.push_back(0.9/100) ; 
	listeMatu.push_back(3) ;	tauxZC.push_back(0.92/100) ;  
	listeMatu.push_back(4) ;	tauxZC.push_back(0.95/100) ; 
	listeMatu.push_back(5) ;	tauxZC.push_back(1.00/100) ; 
	listeMatu.push_back(10) ;	tauxZC.push_back(1.5/100) ; 
	listeMatu.push_back(15) ;	tauxZC.push_back(2.0/100) ;  
	listeMatu.push_back(20) ;	tauxZC.push_back(2.5/100) ;
	listeMatu.push_back(25) ;	tauxZC.push_back(2.3/100) ;
	CourbeInput_PTR courbe_PTR_test(new CourbeInput(listeMatu, tauxZC));

	Piecewiseconst_RR_Function sigma	= Piecewiseconst_RR_Function(x, y) ; 
	Piecewiseconst_RR_Function m		= Piecewiseconst_RR_Function(x, y) ;
	CheyetteDD_Model::CheyetteDD_Parameter monStruct = CheyetteDD_Model::CheyetteDD_Parameter(k, sigma, m) ;
	CheyetteDD_Model_PTR modele_test_PTR(new CheyetteDD_Model(courbe_PTR_test, monStruct)) ;

	double fwdProbaT = 1 ;
	std::vector<double>			datesOfSimulation ;	
	datesOfSimulation.push_back(0.5) ; datesOfSimulation.push_back(1) ; 
	std::vector<size_t>			discretizationBetweenDates ;
	discretizationBetweenDates.push_back(10) ; discretizationBetweenDates.push_back(5) ;

	MC_Cheyette MC_Cheyette_Test(modele_test_PTR, rnGenerator, fwdProbaT, datesOfSimulation, discretizationBetweenDates) ;


	std::vector<double> x_t = MC_Cheyette_Test.get_x_t_Cheyette_() ;
	std::vector<double> y_t = MC_Cheyette_Test.get_y_t_Cheyette_() ;

	//DEBUG
	for (size_t i = 0 ; i < x_t.size() ; ++i)    
	{			
		std::cout << "x_t " << i << " : " << x_t[i] << ", y_t " << i << " : " << y_t[i] << std::endl ;
	}


}

//double ZCVasicek(double matuT, double mean_rev, double level, double sigma)
//{
//	double A_0_T = (1 - exp(- mean_rev * matuT)) / mean_rev ;
//	double D_0_T = (level - sigma*sigma /(2*mean_rev*mean_rev) ) * (A_0_T - matuT) - sigma*sigma * A_0_T*A_0_T/(4 * mean_rev) ;
//
//}

void TestMCSwapPricer()
{
	unsigned long seed = 47;
	RNGenerator_PTR  rnGenerator(new McGenerator(seed));

	std::vector<double> listeMatu, tauxZC ;
	listeMatu.push_back(0) ;	tauxZC.push_back(0.8/100) ; 
	listeMatu.push_back(1) ;	tauxZC.push_back(0.85/100) ; 
	listeMatu.push_back(2) ;	tauxZC.push_back(0.9/100) ; 
	listeMatu.push_back(3) ;	tauxZC.push_back(0.92/100) ;  
	listeMatu.push_back(4) ;	tauxZC.push_back(0.95/100) ; 
	listeMatu.push_back(5) ;	tauxZC.push_back(1.00/100) ; 
	listeMatu.push_back(10) ;	tauxZC.push_back(1.5/100) ; 
	listeMatu.push_back(15) ;	tauxZC.push_back(2.0/100) ;  
	listeMatu.push_back(20) ;	tauxZC.push_back(2.5/100) ;
	listeMatu.push_back(25) ;	tauxZC.push_back(2.3/100) ;
	CourbeInput_PTR courbe_PTR_test(new CourbeInput(listeMatu, tauxZC));


	std::vector<double> x, y ;
	x.push_back(0) ; x.push_back(1) ; x.push_back(2) ; 
	y.push_back(0.25) ; y.push_back(0.5) ;
	double k(0.25) ;
	Piecewiseconst_RR_Function sigma	= Piecewiseconst_RR_Function(x, y) ; 
	Piecewiseconst_RR_Function m		= Piecewiseconst_RR_Function(x, y) ;
	CheyetteDD_Model::CheyetteDD_Parameter monStruct = CheyetteDD_Model::CheyetteDD_Parameter(k, sigma, m) ;
	CheyetteDD_Model_PTR modele_test_PTR(new CheyetteDD_Model(courbe_PTR_test, monStruct)) ;
	modele_test_PTR->show() ;

	double strike          = 0.04;
	LMM::Index  indexStart = 2; 
	LMM::Index  indexEnd   = 4; 
	Tenor	floatingLegTenorType = Tenor::_6M;
	Tenor	fixedLegTenorType    = Tenor::_1YR;
	LMMTenorStructure_PTR simulationStructure(new LMMTenorStructure(Tenor::_6M , 10) );
	VanillaSwap swap = VanillaSwap(strike, indexStart, indexEnd, floatingLegTenorType, fixedLegTenorType, simulationStructure);
	swap.show() ;
	double fwdProbaT = 1 ;
	std::vector<double>	datesOfSimulation ;
	std::vector<size_t> discretizationBetweenDates ;
	//for(size_t i = 1 ; i < 20 ; ++i)
	//{
	//	datesOfSimulation.push_back(0.5 * i) ;
	//	discretizationBetweenDates.push_back(4) ; 
	//}
	datesOfSimulation.push_back(1) ; datesOfSimulation.push_back(1.5) ; datesOfSimulation.push_back(2) ;
	discretizationBetweenDates.push_back(200) ; discretizationBetweenDates.push_back(100) ; discretizationBetweenDates.push_back(100) ;

	MC_Cheyette_PTR mc_Cheyette_Test_PTR(new MC_Cheyette(modele_test_PTR, 
														rnGenerator, 
														fwdProbaT, 
														datesOfSimulation, 
														discretizationBetweenDates) ) ;
	MC_CheyetteDD_VanillaSwapPricer mc_Cheyette_vanillaSwpaPricer_Test(mc_Cheyette_Test_PTR) ;
	
	std::cout << "voir pour changer la graine a chaque simu" << std::endl ;
	double t_valo = 1 ;
	size_t nbSimu = 1000 ;

	MC_CheyetteDD_ZCPricer mc_CheyetteDD_ZCPricer(mc_Cheyette_Test_PTR) ;
	std::cout << "Prix ZC : " << mc_CheyetteDD_ZCPricer.price(t_valo, 2, nbSimu) << std::endl ;

	std::cout << "Prix swap MC : " << mc_Cheyette_vanillaSwpaPricer_Test.swapNPV(t_valo, swap, nbSimu) << std::endl ;

	VanillaSwaption swaption(swap, OptionType::OptionType::CALL) ;
	MC_CheyetteDD_VanillaSwaptionPricer mc_Cheyette_vanillaSwaptionPricer_Test(mc_Cheyette_Test_PTR) ;
	std::cout << "Prix swaption MC : " << mc_Cheyette_vanillaSwaptionPricer_Test.price(t_valo, swaption, nbSimu) << std::endl ;

}