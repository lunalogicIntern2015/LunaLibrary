#include "TestMC.h"
#include "Cheyette\unit_test\TestApproxDD.h"
//
//void UneTrajectoireEuler()
//{
//	unsigned long seed = 47;
//	RNGenerator_PTR  rnGenerator(new McGenerator(seed));
//
//	std::vector<double> x, y ;
//	x.push_back(0) ; x.push_back(0.5) ; x.push_back(1) ; 
//	y.push_back(0.25) ; y.push_back(0.5) ;
//	double k(1) ;
//
//	std::vector<double> listeMatu, tauxZC ;
//	listeMatu.push_back(0) ;	tauxZC.push_back(0.8/100) ; 
//	listeMatu.push_back(1) ;	tauxZC.push_back(0.85/100) ; 
//	listeMatu.push_back(2) ;	tauxZC.push_back(0.9/100) ; 
//	listeMatu.push_back(3) ;	tauxZC.push_back(0.92/100) ;  
//	listeMatu.push_back(4) ;	tauxZC.push_back(0.95/100) ; 
//	listeMatu.push_back(5) ;	tauxZC.push_back(1.00/100) ; 
//	listeMatu.push_back(10) ;	tauxZC.push_back(1.5/100) ; 
//	listeMatu.push_back(15) ;	tauxZC.push_back(2.0/100) ;  
//	listeMatu.push_back(20) ;	tauxZC.push_back(2.5/100) ;
//	listeMatu.push_back(25) ;	tauxZC.push_back(2.3/100) ;
//	CourbeInput_PTR courbe_PTR_test(new CourbeInput(listeMatu, tauxZC));
//
//	Piecewiseconst_RR_Function sigma	= Piecewiseconst_RR_Function(x, y) ; 
//	Piecewiseconst_RR_Function m		= Piecewiseconst_RR_Function(x, y) ;
//	CheyetteDD_Model::CheyetteDD_Parameter monStruct = CheyetteDD_Model::CheyetteDD_Parameter(k, sigma, m) ;
//	CheyetteDD_Model_PTR modele_test_PTR(new CheyetteDD_Model(courbe_PTR_test, monStruct)) ;
//
//	double fwdProbaT = 1 ;
//	std::vector<double>			datesOfSimulation ;	
//	datesOfSimulation.push_back(0.5) ; datesOfSimulation.push_back(1) ; 
//	std::vector<size_t>			discretizationBetweenDates ;
//	discretizationBetweenDates.push_back(10) ; discretizationBetweenDates.push_back(5) ;
//
//	MC_Cheyette MC_Cheyette_Test(modele_test_PTR, rnGenerator, fwdProbaT, datesOfSimulation, discretizationBetweenDates) ;
//
//
//	std::vector<double> x_t = MC_Cheyette_Test.get_x_t_Cheyette() ;
//	std::vector<double> y_t = MC_Cheyette_Test.get_y_t_Cheyette() ;
//
//	//DEBUG
//	for (size_t i = 0 ; i < x_t.size() ; ++i)    
//	{			
//		std::cout << "x_t " << i << " : " << x_t[i] << ", y_t " << i << " : " << y_t[i] << std::endl ;
//	}
//
//
//}
//
//
void TestMCSwapPricer()
{
	unsigned long seed = 47;
	RNGenerator_PTR  rnGenerator(new McGenerator(seed));

	CourbeInput_PTR courbe_PTR_test(createCourbeInput(3));

	std::vector<double> x, sigma_y, m_y ;
	x.push_back(0) ; x.push_back(5) ; x.push_back(10) ; 
	sigma_y.push_back(0.15) ; sigma_y.push_back(0.20) ;
	m_y.push_back(0.3) ; m_y.push_back(0.1) ;
	double k(0.02) ;
	Piecewiseconst_RR_Function sigma	= Piecewiseconst_RR_Function(x, sigma_y) ; 
	Piecewiseconst_RR_Function m		= Piecewiseconst_RR_Function(x, m_y) ;
	CheyetteDD_Model::CheyetteDD_Parameter monStruct = CheyetteDD_Model::CheyetteDD_Parameter(k, sigma, m) ;
	int shiftChoice = 1 ;
	CheyetteDD_Model_PTR modele_test_PTR(new CheyetteDD_Model(courbe_PTR_test, monStruct, shiftChoice)) ;
	modele_test_PTR->show() ;

	double strike          = - 1000. ;

	//swaption 1Y1Y
	//LMM::Index  indexStart = 2; 
	//LMM::Index  indexEnd   = 4; 

	//swaption 3Y1Y
	//LMM::Index  indexStart = 6; 
	//LMM::Index  indexEnd   = 8; 

	//swaption 1Y3Y
	//LMM::Index  indexStart = 2; 
	//LMM::Index  indexEnd   = 8; 

	//swaption 1Y9Y
	LMM::Index  indexStart = 2; 
	LMM::Index  indexEnd   = 20; 

	Tenor	floatingLegTenorType = Tenor::_6M;
	Tenor	fixedLegTenorType    = Tenor::_1YR;
	LMMTenorStructure_PTR simulationStructure(new LMMTenorStructure(Tenor::_6M , 20) );
	VanillaSwap swap = VanillaSwap(strike, indexStart, indexEnd, floatingLegTenorType, fixedLegTenorType, simulationStructure);
	VanillaSwap_PTR pSwap(new 
		VanillaSwap(strike, indexStart, indexEnd, floatingLegTenorType, fixedLegTenorType, simulationStructure)) ;
	swap.show() ;
	double fwdProbaT = 3 ;
	std::vector<size_t>	indexOfSimulation ;
	std::vector<size_t> discretizationBetweenDates ;
		indexOfSimulation.push_back(0) ;		
		discretizationBetweenDates.push_back(0) ; 
	for(size_t i = 1 ; i <= indexEnd ; ++i)
	{
		indexOfSimulation.push_back(i) ;		//INDEX !!! de simulation 
		discretizationBetweenDates.push_back(20) ; 
	}
	//indexOfSimulation.push_back(1) ; indexOfSimulation.push_back(1.5) ; indexOfSimulation.push_back(2) ;
	//discretizationBetweenDates.push_back(200) ; discretizationBetweenDates.push_back(100) ; discretizationBetweenDates.push_back(100) ;

	MC_CheyetteDD_VanillaSwapPricer_PTR mc_swap_pricer_PTR( 
					new	MC_CheyetteDD_VanillaSwapPricer(modele_test_PTR, 
														rnGenerator, 
														floatingLegTenorType,
														fwdProbaT,
														indexOfSimulation, 
														discretizationBetweenDates) ) ;

	size_t nbSimu = 20000 ;
	

	//float - fixed
	std::cout << "Prix swap MC : " << mc_swap_pricer_PTR->swapNPV(pSwap, nbSimu, floatingLegTenorType, fixedLegTenorType) 
								<< std::endl ;
	double dateStart	= indexStart / 2.;
	double dateEnd		= indexEnd / 2.; 
	double tauxStart	= courbe_PTR_test->get_tauxZC0(dateStart) ;
	double tauxEnd		= courbe_PTR_test->get_tauxZC0(dateEnd) ;
	double ZC_T0 = exp( - dateStart * tauxStart ) ; //exp( - indexStart / 2. * tauxStart ) ;
	double ZC_TN = exp( - dateEnd * tauxEnd ) ;
	double floatLeg = ZC_T0 - ZC_TN ;
	

	std::vector<size_t> fixedLegIndex = pSwap->get_fixedLegPaymentIndexSchedule() ;
	double fixedLeg = 0. ;
	for (size_t i = 0 ; i < fixedLegIndex.size() ; ++i)
	{
		double date	= fixedLegIndex[i] / 2.;
		fixedLeg += exp( - date * courbe_PTR_test->get_tauxZC0(date) ) ;
	}

	std::cout << "Prix swap courbe : " << floatLeg - strike * fixedLeg << std::endl ;


	VanillaSwaption_PTR pSwaption(new VanillaSwaption(swap, OptionType::OptionType::CALL)) ;
	MC_CheyetteDD_VanillaSwaptionPricer_PTR mc_swaption_pricer(
							new MC_CheyetteDD_VanillaSwaptionPricer(modele_test_PTR, 
																			rnGenerator, 
																			floatingLegTenorType,
																			fwdProbaT,
																			indexOfSimulation, 
																			discretizationBetweenDates) ) ;

	std::cout << "Prix swaption MC : " << mc_swaption_pricer->price(pSwaption, nbSimu)[0] << std::endl ;

}