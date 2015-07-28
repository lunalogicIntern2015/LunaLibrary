#include "TestMC.h"
#include "Cheyette\unit_test\TestApproxDD.h"

#include <fstream>
#include <ostream>


//courbe 0 : courbe plate � 1%
//courbe 1 : courbe test 
//courbe 2 : courbe march� interbancaire du 22-06-15
CourbeInput_PTR createCourbeInput(int curveChoice)
{
	switch (curveChoice)
	{
	case 1:{
		std::vector<double> listeMatu, tauxZC ;
		double translation = 0. ;
		listeMatu.push_back(0) ;	tauxZC.push_back(1./100 + translation) ; 
		listeMatu.push_back(1) ;	tauxZC.push_back(1./100 + translation) ; 
		listeMatu.push_back(2) ;	tauxZC.push_back(1./100 + translation) ; 
		listeMatu.push_back(3) ;	tauxZC.push_back(1./100 + translation) ;  
		listeMatu.push_back(4) ;	tauxZC.push_back(1./100 + translation) ; 
		listeMatu.push_back(5) ;	tauxZC.push_back(1./100 + translation) ; 
		listeMatu.push_back(10) ;	tauxZC.push_back(1./100 + translation) ; 
		listeMatu.push_back(15) ;	tauxZC.push_back(1./100 + translation) ;  
		listeMatu.push_back(20) ;	tauxZC.push_back(1./100 + translation) ;
		listeMatu.push_back(25) ;	tauxZC.push_back(1./100 + translation) ;
		CourbeInput_PTR courbe_PTR_test(new CourbeInput(listeMatu, tauxZC));
		return courbe_PTR_test ;
		break ;
		   }
	case 2:{
		std::vector<double> listeMatu, tauxZC ;
		double translation = 0.0 ;
		listeMatu.push_back(0) ;	tauxZC.push_back(0.8/100 + translation) ; 
		listeMatu.push_back(1) ;	tauxZC.push_back(0.85/100 + translation) ; 
		listeMatu.push_back(2) ;	tauxZC.push_back(0.9/100 + translation) ; 
		listeMatu.push_back(3) ;	tauxZC.push_back(0.92/100 + translation) ;  
		listeMatu.push_back(4) ;	tauxZC.push_back(0.95/100 + translation) ; 
		listeMatu.push_back(5) ;	tauxZC.push_back(1.00/100 + translation) ; 
		listeMatu.push_back(10) ;	tauxZC.push_back(1.5/100 + translation) ; 
		listeMatu.push_back(15) ;	tauxZC.push_back(2.0/100 + translation) ;  
		listeMatu.push_back(20) ;	tauxZC.push_back(2.5/100 + translation) ;
		listeMatu.push_back(25) ;	tauxZC.push_back(2.3/100 + translation) ;
		CourbeInput_PTR courbe_PTR_test(new CourbeInput(listeMatu, tauxZC));
		return courbe_PTR_test ;
		break ;
		   }
	case 3:{
		std::vector<double> listeMatu, tauxZC ;
		double translation = 0.0 ;
		listeMatu.push_back(0) ;			tauxZC.push_back(-0.12/100  + translation) ; 
		listeMatu.push_back(0.002777778) ;	tauxZC.push_back(-0.12/100  + translation) ; 
		listeMatu.push_back(0.019444444) ;	tauxZC.push_back(-0.131/100  + translation) ; 
		listeMatu.push_back(0.083333333) ;	tauxZC.push_back(-0.064/100  + translation) ; 
		listeMatu.push_back(0.166666667) ;	tauxZC.push_back(-0.039/100  + translation) ;  
		listeMatu.push_back(0.25) ;			tauxZC.push_back(-0.014/100  + translation) ; 
		listeMatu.push_back(0.5) ;			tauxZC.push_back(0.05/100  + translation) ; 
		listeMatu.push_back(0.583333333) ;	tauxZC.push_back(0.045/100  + translation) ; 
		listeMatu.push_back(0.666666667) ;	tauxZC.push_back(0.045/100  + translation) ;  
		listeMatu.push_back(0.75) ;			tauxZC.push_back(0.05/100  + translation) ;
		listeMatu.push_back(0.833333333) ;	tauxZC.push_back(0.056/100  + translation) ;
		listeMatu.push_back(0.916666667) ;	tauxZC.push_back(0.064/100  + translation) ; 
		listeMatu.push_back(1) ;			tauxZC.push_back(0.073/100  + translation) ; 
		listeMatu.push_back(1.5) ;			tauxZC.push_back(0.098/100  + translation) ; 
		listeMatu.push_back(2 ) ;			tauxZC.push_back(0.135/100  + translation) ;  
		listeMatu.push_back(3 ) ;			tauxZC.push_back(0.24/100  + translation) ; 
		listeMatu.push_back(4 ) ;			tauxZC.push_back(0.382/100  + translation) ; 
		listeMatu.push_back(5 ) ;			tauxZC.push_back(0.544/100  + translation) ; 
		listeMatu.push_back(6 ) ;			tauxZC.push_back(0.701/100  + translation) ;  
		listeMatu.push_back(7 ) ;			tauxZC.push_back(0.847/100  + translation) ;
		listeMatu.push_back(8 ) ;			tauxZC.push_back(0.981/100  + translation) ;
		listeMatu.push_back(9 ) ;			tauxZC.push_back(1.1/100  + translation) ;  
		listeMatu.push_back(10) ;			tauxZC.push_back(1.198/100  + translation) ; 
		listeMatu.push_back(11) ;			tauxZC.push_back(1.286/100  + translation) ; 
		listeMatu.push_back(12) ;			tauxZC.push_back(1.359/100  + translation) ; 
		listeMatu.push_back(15) ;			tauxZC.push_back(1.516/100  + translation) ;  
		listeMatu.push_back(20) ;			tauxZC.push_back(1.629/100  + translation) ;
		listeMatu.push_back(25) ;			tauxZC.push_back(1.657/100  + translation) ;
		listeMatu.push_back(30) ;			tauxZC.push_back(1.666/100  + translation) ; 
		listeMatu.push_back(35) ;			tauxZC.push_back(1.677/100  + translation) ; 
		listeMatu.push_back(40) ;			tauxZC.push_back(1.675/100  + translation) ;  
		listeMatu.push_back(45) ;			tauxZC.push_back(1.656/100  + translation) ;
		listeMatu.push_back(50) ;			tauxZC.push_back(1.634/100  + translation) ;
		CourbeInput_PTR courbe_PTR_test(new CourbeInput(listeMatu, tauxZC));		   
		return courbe_PTR_test ;
		break ;
		   }
	case 4:{
		std::vector<double> listeMatu, tauxZC ;
		double translation = 0.0 ;
		listeMatu.push_back(0) ;			tauxZC.push_back(0.2906/100  + translation) ; 
		listeMatu.push_back(0.25) ;			tauxZC.push_back(0.2906/100  + translation) ; 
		listeMatu.push_back(0.5) ;			tauxZC.push_back(0.3502/100  + translation) ; 
		listeMatu.push_back(14./12.) ;		tauxZC.push_back(0.5522/100  + translation) ; 
		listeMatu.push_back(20./12.) ;		tauxZC.push_back(0.7285/100  + translation) ;  
		listeMatu.push_back(4 ) ;			tauxZC.push_back(1.4524/100  + translation) ; 
		listeMatu.push_back(6 ) ;			tauxZC.push_back(1.9009/100  + translation) ;  
		listeMatu.push_back(8 ) ;			tauxZC.push_back(2.2156/100  + translation) ;
		listeMatu.push_back(10 ) ;			tauxZC.push_back(2.4322/100  + translation) ;
		listeMatu.push_back(12 ) ;			tauxZC.push_back(2.5869/100  + translation) ;
		listeMatu.push_back(15 ) ;			tauxZC.push_back(2.7385/100  + translation) ;
		listeMatu.push_back(20 ) ;			tauxZC.push_back(2.8776/100  + translation) ;
		CourbeInput_PTR courbe_PTR_test(new CourbeInput(listeMatu, tauxZC));		   
		return courbe_PTR_test ;
		break ;
		   }
	case 7:{
		std::vector<double> listeMatu, tauxZC ;
		double translation = 0.0 ;
		listeMatu.push_back(0) ;			tauxZC.push_back(4.7/100  + translation) ; 
		listeMatu.push_back(0.5) ;			tauxZC.push_back(4.7/100  + translation) ; 
		listeMatu.push_back(0.8) ;			tauxZC.push_back(4.623/100  + translation) ; 
		listeMatu.push_back(1) ;			tauxZC.push_back(4.573/100  + translation) ;  
		listeMatu.push_back(1.5 ) ;			tauxZC.push_back(4.437/100  + translation) ; 
		listeMatu.push_back(2 ) ;			tauxZC.push_back(4.353/100  + translation) ;  
		CourbeInput_PTR courbe_PTR_test(new CourbeInput(listeMatu, tauxZC));		   
		return courbe_PTR_test ;
		break ;
		   }
	default :
		throw "courbe non existante" ;
	}
}



//CheyetteDD Model + swaption ATM
CheyetteDD_VanillaSwaptionApproxPricer_PTR creeModeleATM(size_t a, size_t b, const Tenor& floatTenor, const Tenor& fixedTenor)
{

//courbe spot
	int curveChoice = 1 ; //courbe plate � 1%
	int shiftChoice = 1 ;
	std::vector<double> x(a+b+1) ;
	std::vector<double> m_y(a+b) ;
	std::vector<double> sigma_y(a+b) ;
	for (size_t i = 0 ; i <= a+b ; ++i)
	{
		x[i] = i ;
	}
	for (size_t i = 0 ; i < a+b ; ++i)
	{
		m_y[i] = 0. ;
		sigma_y[i] = 0.20 ;
	}
	Piecewiseconst_RR_Function sigma(x, sigma_y) ; 
	Piecewiseconst_RR_Function m(x, m_y) ;
	double k = 0.02 ;
	CheyetteDD_Model::CheyetteDD_Parameter monStruct(k, sigma, m) ;

	CourbeInput_PTR courbe_PTR_test(createCourbeInput(curveChoice));
	CheyetteDD_Model_PTR modele_test_PTR(new CheyetteDD_Model(courbe_PTR_test, monStruct, shiftChoice)) ;

//approx
	double tenor = std::min(floatTenor.YearFraction() , fixedTenor.YearFraction() ) ;
	
	size_t indexStart = size_t(a / tenor) ;
	size_t indexEnd = size_t((a+b) / tenor) ;
	
	double strike =  - 1000. ; //temporaire avant de mettre strike ATM
	LMMTenorStructure_PTR simulationStructure(new LMMTenorStructure(floatTenor, a+b+1) );
	VanillaSwap swap = VanillaSwap(strike, indexStart, indexEnd, floatTenor, fixedTenor, simulationStructure);
	
	VanillaSwaption_PTR swaption_PTR_test(new VanillaSwaption(swap, OptionType::OptionType::CALL)) ;
	CheyetteDD_VanillaSwaptionApproxPricer approx(modele_test_PTR, swaption_PTR_test);

//calcul du strike ATM pour le swap
	double strikeATM = approx.swapRate0() ;
	swap.set_strike(strikeATM) ;
	VanillaSwaption_PTR swaption_PTR_test_ATM(new VanillaSwaption(swap, OptionType::OptionType::CALL)) ;
	swaption_PTR_test_ATM->show() ;
	CheyetteDD_VanillaSwaptionApproxPricer_PTR pApproxATM(
				new CheyetteDD_VanillaSwaptionApproxPricer(modele_test_PTR, swaption_PTR_test_ATM)) ;

	return pApproxATM ;
	
}


//CheyetteDD Model + swaption de strike choisi
CheyetteDD_VanillaSwaptionApproxPricer_PTR creeModele(size_t a, size_t b, 
														const Tenor& floatTenor, const Tenor& fixedTenor,
														double strike)
{

//courbe spot
	int curveChoice = 1 ; //courbe plate � 1%
	int shiftChoice = 1 ;
	std::vector<double> x(a+b+2) ;
	std::vector<double> m_y(a+b+1) ;
	std::vector<double> sigma_y(a+b+1) ;
	for (size_t i = 0 ; i < a+b+2 ; ++i)
	{
		x[i] = i ;
	}
	for (size_t i = 0 ; i < a+b+1 ; ++i)
	{
		m_y[i] = 0. ;
		sigma_y[i] = 0.20 ;
	}
	Piecewiseconst_RR_Function sigma(x, sigma_y) ; 
	Piecewiseconst_RR_Function m(x, m_y) ;
	double k = 0.02 ;
	CheyetteDD_Model::CheyetteDD_Parameter monStruct(k, sigma, m) ;

	CourbeInput_PTR courbe_PTR_test(createCourbeInput(curveChoice));
	CheyetteDD_Model_PTR modele_test_PTR(new CheyetteDD_Model(courbe_PTR_test, monStruct, shiftChoice)) ;

//approx
	double tenor = std::min(floatTenor.YearFraction() , fixedTenor.YearFraction() ) ;
	
	size_t indexStart = size_t(a / tenor) ;
	size_t indexEnd = size_t((a+b) / tenor) ;
	
	LMMTenorStructure_PTR simulationStructure(new LMMTenorStructure(floatTenor, a+b+1) );
	VanillaSwap swap = VanillaSwap(strike, indexStart, indexEnd, floatTenor, fixedTenor, simulationStructure);
	
	VanillaSwaption_PTR swaption_PTR_test(new VanillaSwaption(swap, OptionType::OptionType::CALL)) ;

	CheyetteDD_VanillaSwaptionApproxPricer_PTR pApprox(
				new CheyetteDD_VanillaSwaptionApproxPricer(modele_test_PTR, swaption_PTR_test)) ;

	return pApprox ;	
}


//proba forward QT
void TestMCSwapPricer(size_t a, size_t b, size_t nbSimus, const Tenor& floatTenor, const Tenor& fixedTenor)
{

	double strike          = 0.02 ; //- 1000. ;

	CheyetteDD_VanillaSwaptionApproxPricer_PTR pApproxPricer = creeModele(a, b, floatTenor, fixedTenor, strike) ;
	CheyetteDD_Model_PTR pModel = pApproxPricer->get_CheyetteDD_Model() ;
	VanillaSwaption_PTR pSwaption = pApproxPricer->get_VanillaSwaption() ;
	VanillaSwap_PTR pSwap(new VanillaSwap(pSwaption->getUnderlyingSwap())) ;

	unsigned long seed = 47;
	RNGenerator_PTR  rnGenerator(new McGenerator(seed));

	double tenor = std::min(floatTenor.YearFraction() , fixedTenor.YearFraction() ) ;
	size_t indexStart = size_t(a / tenor) ;
	size_t indexEnd = size_t((a+b) / tenor) ;

	size_t fwdProbaT = a + b ; 
	size_t discretizationBetweenDates = 180 ;
	LMMTenorStructure_PTR simulationStructure(new LMMTenorStructure(floatTenor, a+b+1) );

	MC_CheyetteDD_VanillaSwapPricer_PTR mc_swap_pricer_PTR( 
					new	MC_CheyetteDD_VanillaSwapPricer(pModel, rnGenerator, simulationStructure,
														fwdProbaT, discretizationBetweenDates) ) ;


	//float - fixed
	size_t valuationIndex = 0 ;
	std::cout << "Prix swap MC : " << mc_swap_pricer_PTR->swapNPV(pSwap, nbSimus, valuationIndex, floatTenor, fixedTenor) 
								<< std::endl ;

	CourbeInput_PTR pCourbe = pApproxPricer->get_CheyetteDD_Model()->get_courbeInput_PTR() ;

	double tauxStart	= pCourbe->get_tauxZC0(a) ;
	double tauxEnd		= pCourbe->get_tauxZC0(a + b) ;
	double ZC_T0 = exp( - double(a) * tauxStart ) ; //exp( - indexStart / 2. * tauxStart ) ;
	double ZC_TN = exp( - double((a+b)) * tauxEnd ) ;
	double floatLeg = ZC_T0 - ZC_TN ;
	

	std::vector<size_t> fixedLegIndex = pSwap->get_fixedLegPaymentIndexSchedule() ;
	double fixedLeg = 0. ;
	for (size_t i = 0 ; i < fixedLegIndex.size() ; ++i)
	{
		double date	= fixedLegIndex[i] * tenor ;
		fixedLeg += exp( - date * pCourbe->get_tauxZC0(date) ) ;
	}

	std::cout << "Prix swap courbe : " << floatLeg - strike * fixedLeg << std::endl ;

	MC_CheyetteDD_VanillaSwaptionPricer_PTR mc_swaption_pricer(
							new MC_CheyetteDD_VanillaSwaptionPricer(pModel, rnGenerator, simulationStructure,
														fwdProbaT, discretizationBetweenDates) ) ;

	std::cout << "Prix swaption MC : " << mc_swaption_pricer->price(pSwaption, nbSimus)[0] << std::endl ;

}


void TestMCSwapPricer_annuity(size_t a, size_t b, size_t nbSimus, const Tenor& floatTenor, const Tenor& fixedTenor)
{
	double strike          = 0.02 ;//- 1000. ;

	CheyetteDD_VanillaSwaptionApproxPricer_PTR pApproxPricer = creeModele(a, b, floatTenor, fixedTenor, strike) ;
	CheyetteDD_Model_PTR pModel = pApproxPricer->get_CheyetteDD_Model() ;
	VanillaSwaption_PTR pSwaption = pApproxPricer->get_VanillaSwaption() ;
	VanillaSwap_PTR pSwap(new VanillaSwap(pSwaption->getUnderlyingSwap())) ;

	double S0 = pApproxPricer->swapRate0() ;	
	double r0 = pApproxPricer->get_CheyetteDD_Model()->get_courbeInput_PTR()->get_f_0_t(0) ;

	double sigma = pApproxPricer->get_CheyetteDD_Model()->get_CheyetteDD_Parameter().sigma_(0.5) ;//cste
	double k = pApproxPricer->get_CheyetteDD_Model()->get_CheyetteDD_Parameter().k_ ;
	
	unsigned long seed = 47;
	RNGenerator_PTR  rnGenerator(new McGenerator(seed));
	
	std::vector<double> res(3) ;
	double somme_xi_swap	  = 0. ;
	double somme_xi2_swap	  = 0. ;
	double somme_xi_swaption  = 0. ;
	double somme_xi2_swaption = 0. ;
	
	size_t nbPas = 200 ;
	double dt = double(a)/double(nbPas) ;
	for (size_t simu = 1 ; simu <= nbSimus ; ++simu)    
	{
		if ((simu*10) % nbSimus == 0){std::cout << double(simu)/double(nbSimus)*100 << "%" << std::endl ;}
		std::vector<double>   gaussian_tmp(nbPas);  
		rnGenerator->generate(gaussian_tmp);		// generate Gaussian.
		double t = 0. ;
		double S_t = S0 ;

		for (size_t pasDiscretisation = 1 ; pasDiscretisation <= nbPas ; ++pasDiscretisation)    
		{
			double y_t, inv ;
			if (t == 0.)
			{
				y_t = 0. ;
				inv = 0. ;		
			}
			else 
			{
				y_t = sigma * sigma * r0 * r0 / (2 * k) * (1 - exp(-2 * k *t)) ;
					//pApproxPricer->get_buffer_y_bar_t(t) ;   // ...
				inv = pApproxPricer->inverse(t, S_t) ;
			}

			double derivee = pApproxPricer->swapRate_1stDerivative(t, inv, y_t) ;
			double sigma = pApproxPricer->get_CheyetteDD_Model()->sigma_r(t, inv, y_t) ;

			double S_t_plus_dt = S_t + derivee * sigma * sqrt(dt) * gaussian_tmp[pasDiscretisation - 1] ;
			S_t = S_t_plus_dt ;
			t += dt ;
		}
		double resSwap = S_t - strike ;
		double resSwaption = std::max(S_t - strike, 0.) ; //  (S(T0) - K)+ 
		
		somme_xi_swap += resSwap ;
		somme_xi2_swap += resSwap * resSwap ;

		somme_xi_swaption += resSwaption ;
		somme_xi2_swaption += resSwaption * resSwaption ;
	}

	double A0 = pApproxPricer->swapRateDenominator(0., 0., 0.) ;
//swap
	double mean_x_swap	= A0 * somme_xi_swap / nbSimus; // A(0) * E^QA( (S(T0) - K)+ )
	double mean_x2_swap	= A0 * A0 * somme_xi2_swap / nbSimus; 
	double variance_swap = mean_x2_swap - mean_x_swap * mean_x_swap ;

	std::cout <<"variance : " << variance_swap << std::endl ;

	double IC_left_swap	= mean_x_swap - 2.57*std::sqrt(variance_swap / nbSimus);
	double IC_right_swap = mean_x_swap + 2.57*std::sqrt(variance_swap / nbSimus);

	std::cout   << "SWAP " << std::endl;
	std::cout   << "MC swap : " << mean_x_swap << std::endl;
	std::cout	<< "nb simulations : " << nbSimus << std::endl;
	std::cout   << "99% confidence interval  [" << IC_left_swap << " , " << IC_right_swap	<< "]" << std::endl;
	
//swaption
	double mean_x_swaption	= A0 * somme_xi_swaption / nbSimus; // A(0) * E^QA( (S(T0) - K)+ )
	double mean_x2_swaption	= A0 * A0 * somme_xi2_swaption / nbSimus; 
	double variance_swaption = mean_x2_swaption - mean_x_swaption * mean_x_swaption ;

	std::cout <<"variance : " << variance_swaption << std::endl ;

	double IC_left_swaption	= mean_x_swaption - 2.57*std::sqrt(variance_swaption / nbSimus);
	double IC_right_swaption = mean_x_swaption + 2.57*std::sqrt(variance_swaption / nbSimus);

	std::cout   << "SWAPTION " << std::endl;
	std::cout   << "MC swaption : " << mean_x_swaption << std::endl;
	std::cout	<< "nb simulations : " << nbSimus << std::endl;
	std::cout   << "99% confidence interval  [" << IC_left_swaption << " , " << IC_right_swaption	<< "]" << std::endl;

	//ofstream o;
	//std::stringstream fileName_s ;
	//std::string directory = LMMPATH::get_output_path() ;
	//fileName_s << directory << "TestMC_" << ".csv" ; 
	//std::string fileName = fileName_s.str();

	//o.open(fileName,  ios::out | ios::app );
	//o	<<	endl;
	//o	<<	endl;
	//o   << "MC Phi(" << a << "Y, S_t) : " << mean_x << std::endl;
	//o	<< "nb simulations : " << nbSimus << std::endl;
	//o   << "99% confidence interval  [" << IC_left << " , " << IC_right	<< "]" << std::endl;
	//o	<< "Phi(1Y, s barre) : " << pApproxPricer->Phi(a, pApproxPricer->inverse(a, pApproxPricer->swapRate0())) ;

	//o.close() ;
}


void MCannuity(size_t a, size_t b, size_t nbSimus, const Tenor& floatTenor, const Tenor& fixedTenor)
{
	
	size_t nbPas = 200 ;
	double dt = double(a)/double(nbPas) ;

	CheyetteDD_VanillaSwaptionApproxPricer_PTR pApproxPricer = creeModeleATM(a, b, floatTenor, fixedTenor) ;
	VanillaSwaption_PTR pSwaption = pApproxPricer->get_VanillaSwaption() ;
	double strikeATM = pSwaption->getUnderlyingSwap().get_strike() ;
	std::cout << "strikeATM : " << strikeATM << std::endl ; 

	double S0 = pApproxPricer->swapRate0() ;		//strike et S(0)
	double sigma = pApproxPricer->get_CheyetteDD_Model()->get_CheyetteDD_Parameter().sigma_(0.5) ;//cste
	double k = pApproxPricer->get_CheyetteDD_Model()->get_CheyetteDD_Parameter().k_ ;
	double r0 = pApproxPricer->get_CheyetteDD_Model()->get_courbeInput_PTR()->get_f_0_t(0) ;

	unsigned long seed = 47;
	RNGenerator_PTR  rnGenerator(new McGenerator(seed));
	
	std::vector<double> res(3) ;
	double somme_xi	 = 0. ;
	double somme_xi2 = 0. ;

	for (size_t simu = 1 ; simu <= nbSimus ; ++simu)    
	{
		if ((simu*10) % nbSimus == 0){std::cout << double(simu)/double(nbSimus)*100 << "%" << std::endl ;}
		std::vector<double>   gaussian_tmp(nbPas);  
		rnGenerator->generate(gaussian_tmp);		// generate Gaussian.
		double t = 0. ;
		double S_t = S0 ;

		for (size_t pasDiscretisation = 1 ; pasDiscretisation <= nbPas ; ++pasDiscretisation)    
		{
			double y_t, inv ;
			if (t == 0.)
			{
				y_t = 0. ;
				inv = 0. ;		
			}
			else 
			{
				y_t = sigma * sigma * r0 * r0 / (2 * k) * (1 - exp(-2 * k *t)) ;
					//pApproxPricer->get_buffer_y_bar_t(t) ;   // ...
				inv = pApproxPricer->inverse(t, S_t) ;
			}

			double derivee = pApproxPricer->swapRate_1stDerivative(t, inv, y_t) ;
			double sigma = pApproxPricer->get_CheyetteDD_Model()->sigma_r(t, inv, y_t) ;

			double S_t_plus_dt = S_t + derivee * sigma * sqrt(dt) * gaussian_tmp[pasDiscretisation - 1] ;
			S_t = S_t_plus_dt ;
			t += dt ;
		}
		double res = std::max(S_t - S0, 0.) ; //  (S(T0) - K)+ 
		
		somme_xi += res ;
		somme_xi2 += res * res ;
	}

	double A0 = pApproxPricer->swapRateDenominator(0., 0., 0.) ;
	double mean_x	= A0 * somme_xi / nbSimus; // A(0) * E^QA( (S(T0) - K)+ )
	double mean_x2	= A0 * A0 * somme_xi2 / nbSimus; 
	double variance = mean_x2 - mean_x * mean_x ;

	double IC_left	= mean_x - 2.57*std::sqrt(variance / nbSimus);
	double IC_right = mean_x + 2.57*std::sqrt(variance / nbSimus);

	res[0] = mean_x ;
	res[1] = IC_left ;
	res[2] = IC_right ;

	std::cout   << "annuity : " << pApproxPricer->swapRateDenominator(0., 0., 0.) << std::endl;
	std::cout	<< "swap rate 0 : " << pApproxPricer->swapRate0() << std::endl;

	std::cout   << "MC swaption : " << mean_x << std::endl;
	std::cout	<< "nb simulations : " << nbSimus << std::endl;
	std::cout   << "99% confidence interval  [" << IC_left << " , " << IC_right	<< "]" << std::endl;
	
	//ofstream o;
	//std::stringstream fileName_s ;
	//std::string directory = LMMPATH::get_output_path() ;
	//fileName_s << directory << "TestMC_" << ".csv" ; 
	//std::string fileName = fileName_s.str();

	//o.open(fileName,  ios::out | ios::app );
	//o	<<	endl;
	//o	<<	endl;
	//o   << "MC Phi(" << a << "Y, S_t) : " << mean_x << std::endl;
	//o	<< "nb simulations : " << nbSimus << std::endl;
	//o   << "99% confidence interval  [" << IC_left << " , " << IC_right	<< "]" << std::endl;
	//o	<< "Phi(1Y, s barre) : " << pApproxPricer->Phi(a, pApproxPricer->inverse(a, pApproxPricer->swapRate0())) ;

	//o.close() ;
}



void MCforward(size_t a, size_t b, size_t nbSimus, const Tenor& floatTenor, const Tenor& fixedTenor)
{
	CheyetteDD_VanillaSwaptionApproxPricer_PTR pApproxPricer = creeModeleATM(a, b, floatTenor, fixedTenor) ;
	CheyetteDD_Model_PTR pModel = pApproxPricer->get_CheyetteDD_Model() ;
	VanillaSwaption_PTR pSwaption = pApproxPricer->get_VanillaSwaption() ;
	double strikeATM = pSwaption->getUnderlyingSwap().get_strike() ;
	std::cout << "strikeATM : " << strikeATM << std::endl ; 
	
	double tenor = std::min(floatTenor.YearFraction() , fixedTenor.YearFraction() ) ;
	size_t indexStart = pSwaption->getUnderlyingSwap().get_indexStart() ; //size_t(a / tenor) ;
	size_t indexEnd = pSwaption->getUnderlyingSwap().get_indexEnd() ;   //size_t((a+b) / tenor) ;
	
	double dateStart = pSwaption->getUnderlyingSwap().get_StartDate() ; //size_t(a / tenor) ;
	double dateEnd = pSwaption->getUnderlyingSwap().get_EndDate() ;   //size_t((a+b) / tenor) ;
	std::cout << "date start : " << dateStart << ", dateEnd : " << dateEnd << std::endl ;

//param MC
	size_t fwdProbaT = a + b ; 	
	size_t discretizationBetweenDates = 200 ;
	LMMTenorStructure_PTR simulationStructure(new LMMTenorStructure(floatTenor, a+b+1) );
	
	unsigned long seed = 47;
	RNGenerator_PTR  rnGenerator(new McGenerator(seed));
	MC_CheyetteDD_VanillaSwaptionPricer_PTR mc(new MC_CheyetteDD_VanillaSwaptionPricer(pModel,
															rnGenerator, simulationStructure, fwdProbaT,
															discretizationBetweenDates   ) ) ;

	std::vector<double> vectPrixMC(1), vectICinf(1), vectICsup(1) ; 
	std::vector<double> resMC = mc->price(pSwaption, nbSimus) ;

	std::cout   << "annuity : " << pApproxPricer->swapRateDenominator(0., 0., 0.) << std::endl;
	std::cout	<< "swap rate 0 : " << pApproxPricer->swapRate0() << std::endl;

	//ofstream o;
	//std::stringstream fileName_s ;
	//std::string directory = LMMPATH::get_output_path() ;
	//fileName_s << directory << "TestMC_" << ".csv" ; 
	//std::string fileName = fileName_s.str();

	//o.open(fileName,  ios::out | ios::app );
	//o	<<	endl;
	//o	<<	endl;
	//o   << "MC Phi(" << a << "Y, S_t) : " << mean_x << std::endl;
	//o	<< "nb simulations : " << nbSimus << std::endl;
	//o   << "99% confidence interval  [" << IC_left << " , " << IC_right	<< "]" << std::endl;
	//o	<< "Phi(1Y, s barre) : " << pApproxPricer->Phi(a, pApproxPricer->inverse(a, pApproxPricer->swapRate0())) ;

	//o.close() ;
}






//param 
//swaption aY bY
//floating leg tenor vs fixed leg tenor (pas vraiment un tenor mais frequence des flux)
//size(x) = size(m_y) + 1 = size(sigma_y) + 1
void test_approx(double strike, size_t a, size_t b, Tenor floatingLegTenor, Tenor fixedLegTenor, 
				 int curveChoice, int shiftChoice, 
				 std::vector<double> x, std::vector<double> m_y, std::vector<double> sigma_y, double k, 
				 std::vector<size_t> nbSimus)
{	

	double tenor = std::min(floatingLegTenor.YearFraction() , fixedLegTenor.YearFraction() ) ;
	size_t indexStart = size_t(a / tenor) ;
	size_t indexEnd = size_t((a+b) / tenor) ;

//courbe spot
	CourbeInput_PTR courbe_PTR_test(createCourbeInput(curveChoice));
	
	Piecewiseconst_RR_Function sigma(x, sigma_y) ; 
	Piecewiseconst_RR_Function m(x, m_y) ;
	CheyetteDD_Model::CheyetteDD_Parameter monStruct(k, sigma, m) ;

	CheyetteDD_Model_PTR modele_test_PTR(new CheyetteDD_Model(courbe_PTR_test, monStruct, shiftChoice)) ;
	modele_test_PTR->show() ;

//param MC
	size_t fwdProbaT = a + b +1 ; 	
	size_t discretizationBetweenDates = 200 ;
	LMMTenorStructure_PTR simulationStructure(new LMMTenorStructure(floatingLegTenor, a+b+1) );

	unsigned long seed = 47;
	RNGenerator_PTR  rnGenerator(new McGenerator(seed));

	MC_CheyetteDD_VanillaSwaptionPricer_PTR mc(
					new MC_CheyetteDD_VanillaSwaptionPricer(modele_test_PTR, rnGenerator,
															simulationStructure, fwdProbaT,
															discretizationBetweenDates   ) ) ;

//approx

	VanillaSwap swap = VanillaSwap(strike, indexStart, indexEnd, floatingLegTenor, fixedLegTenor, simulationStructure);
	
	VanillaSwaption_PTR swaption_PTR_test(new VanillaSwaption(swap, OptionType::OptionType::CALL)) ;
	CheyetteDD_VanillaSwaptionApproxPricer approx(modele_test_PTR, swaption_PTR_test);

	std::cout << "prixSwaption" << std::endl ;
	double approxPrice = approx.prixSwaptionApproxPiterbarg() ;
	double b_barre = approx.get_buffer_b_barre_() ;
	std::cout << "approxPrice : " << approxPrice << std::endl ;
	std::cout << "b barre :     " << b_barre << std::endl ;
	std::cout << "  " << std::endl ;

//MC
	//double nominal = 1 ;
	//LMMTenorStructure_CONSTPTR swapStructure(new LMMTenorStructure(floatingLegTenor, a+b+1)) ;
	//GeneticSwap_CONSTPTR genericSwapTest = InstrumentFactory::createVanillaSwap(strike, indexStart, 
	//												indexEnd, floatingLegTenor, fixedLegTenor,
	//												swapStructure, nominal);

	//size_t maturity = indexStart ; //index
	//GeneticSwaption_CONSTPTR 	genericSwaptionTest(new GeneticSwaption(maturity, genericSwapTest)) ;

	size_t nbMC = nbSimus.size() ;
	std::vector<double> vectPrixMC(nbMC), vectICinf(nbMC), vectICsup(nbMC) ; 
	for (size_t i = 0 ; i < nbMC ; ++i)
	{
		//std::vector<double> resMC = mc->price(genericSwaptionTest, nbSimus[i]) ;  //Generic
		std::vector<double> resMC = mc->price(swaption_PTR_test, nbSimus[i]) ;  //Vanilla
		vectPrixMC[i] = resMC[0] ;
		vectICinf[i]  = resMC[1] ;
		vectICsup[i]  = resMC[2] ;
	}
	
	double annuityA0 = approx.swapRateDenominator(0., 0., 0.) ;
	double swapRateS0 = approx.swapRate0() ;

	double blackPrice = NumericalMethods::Black_SwaptionImpliedVolatility(approxPrice, annuityA0, 
																		swapRateS0, strike, a) ;

	//mc->printMC_vs_approx(approxPrice, b_barre, annuityA0, swapRateS0, blackPrice,
	//						a, b, genericSwaptionTest, nbSimus, vectPrixMC, vectICinf, vectICsup) ;

	mc->printMC_vs_approx(approxPrice, b_barre, annuityA0, swapRateS0, blackPrice,
							a, b, swaption_PTR_test, nbSimus, vectPrixMC, vectICinf, vectICsup) ;

}

//param 
//swaption aY bY
//floating leg tenor vs fixed leg tenor (pas vraiment un tenor mais frequence des flux)
//size(x) = size(m_y) + 1 = size(sigma_y) + 1
void test_approx_ATM(size_t a, size_t b, Tenor floatingLegTenor, Tenor fixedLegTenor, 
					 int curveChoice, int shiftChoice, 
					 std::vector<double> x, std::vector<double> m_y, std::vector<double> sigma_y, double k, 
					 std::vector<size_t> nbSimus)
{	

	double tenor = std::min(floatingLegTenor.YearFraction() , fixedLegTenor.YearFraction() ) ;
	size_t indexStart = size_t(a / tenor) ;
	size_t indexEnd = size_t((a+b) / tenor) ;
	size_t fwdProbaT = a + b ; 

//courbe spot
	CourbeInput_PTR courbe_PTR_test(createCourbeInput(curveChoice));
	
	Piecewiseconst_RR_Function sigma(x, sigma_y) ; 
	Piecewiseconst_RR_Function m(x, m_y) ;
	CheyetteDD_Model::CheyetteDD_Parameter monStruct(k, sigma, m) ;

	CheyetteDD_Model_PTR modele_test_PTR(new CheyetteDD_Model(courbe_PTR_test, monStruct, shiftChoice)) ;
	modele_test_PTR->show() ;

//approx
	double strike = 1000. ; //temporaire avant de mettre strike ATM
	LMMTenorStructure_PTR simulationStructure(new LMMTenorStructure(floatingLegTenor, a+b+1) );
	VanillaSwap swap = VanillaSwap(strike, indexStart, indexEnd, floatingLegTenor, fixedLegTenor, simulationStructure);
	
	VanillaSwaption_PTR swaption_PTR_test(new VanillaSwaption(swap, OptionType::OptionType::CALL)) ;
	CheyetteDD_VanillaSwaptionApproxPricer approx(modele_test_PTR, swaption_PTR_test);

//calcul du strike ATM pour le swap
	double strikeATM = approx.swapRate0() ;
	swap.set_strike(strikeATM) ;
	VanillaSwaption_PTR swaption_PTR_test_ATM(new VanillaSwaption(swap, OptionType::OptionType::CALL)) ;
	swaption_PTR_test_ATM->show() ;

//approx
	CheyetteDD_VanillaSwaptionApproxPricer approxATM(modele_test_PTR, swaption_PTR_test_ATM);
	double approxPrice = approxATM.prixSwaptionApproxPiterbarg() ;
	double b_barre = approxATM.get_buffer_b_barre_() ;

	std::cout << "prixSwaption" << std::endl ;
	std::cout << "approxPrice : " << approxPrice << std::endl ;
	std::cout << "b barre :     " << b_barre << std::endl ;
	std::cout << "  " << std::endl ;

//param MC
	size_t	discretizationBetweenDates = 200 ;
	
	unsigned long seed = 47;
	RNGenerator_PTR  rnGenerator(new McGenerator(seed));
	MC_CheyetteDD_VanillaSwaptionPricer_PTR mc(
					new MC_CheyetteDD_VanillaSwaptionPricer(modele_test_PTR, rnGenerator,
															simulationStructure, fwdProbaT,
															discretizationBetweenDates   ) ) ;
	size_t nbMC = nbSimus.size() ;
	std::vector<double> vectPrixMC(nbMC), vectICinf(nbMC), vectICsup(nbMC) ; 
	for (size_t i = 0 ; i < nbMC ; ++i)
	{
		//std::vector<double> resMC = mc->price(genericSwaptionTest, nbSimus[i]) ;
		std::vector<double> resMC = mc->price(swaption_PTR_test_ATM, nbSimus[i]) ;
		vectPrixMC[i] = resMC[0] ;
		vectICinf[i]  = resMC[1] ;
		vectICsup[i]  = resMC[2] ;
	}
	
	double annuityA0 = approxATM.swapRateDenominator(0., 0., 0.) ;
	double swapRateS0 = approxATM.swapRate0() ;

	double volBlack = NumericalMethods::Black_SwaptionImpliedVolatility(approxPrice, annuityA0, 
																		swapRateS0, strikeATM, a) ;

	std::cout << "strike :    " << strikeATM << std::endl ;
	std::cout << "expiry :    " << a << std::endl ;
	std::cout << "annuity :   " << annuityA0 << std::endl ;
	std::cout << "swapRate :  " << swapRateS0 << std::endl ;
	std::cout << "vol Black : " << volBlack << std::endl ;

	//mc->printMC_vs_approx(approxPrice, b_barre, annuityA0, swapRateS0, blackPrice,
	//						a, b, genericSwaptionTest, nbSimus, vectPrixMC, vectICinf, vectICsup) ;
	std::cout << "MonteCarlo proba forward QT, T = " << fwdProbaT << std::endl ;
	mc->printMC_vs_approx(approxPrice, b_barre, annuityA0, swapRateS0, volBlack,
							a, b, swaption_PTR_test_ATM, nbSimus, vectPrixMC, vectICinf, vectICsup) ;
}

//param 
//swaption aY bY
//floating leg tenor vs fixed leg tenor (pas vraiment un tenor mais frequence des flux)
//size(x) = size(m_y) + 1 = size(sigma_y) + 1
void test_approx_ATM(size_t a, size_t b, Tenor floatingLegTenor, Tenor fixedLegTenor, 
					 CheyetteDD_Model_PTR pCheyetteDD_Model,
					 std::vector<size_t> nbSimus, std::ofstream& o)
{	

	double tenor = std::min(floatingLegTenor.YearFraction() , fixedLegTenor.YearFraction() ) ;
	size_t indexStart = size_t(a / tenor) ;
	size_t indexEnd = size_t((a+b) / tenor) ;
	double fwdProbaT = double(b) ;  //size_t vers double ok

//param MC
	size_t	discretizationBetweenDates = 200 ;
	LMMTenorStructure_PTR simulationStructure(new LMMTenorStructure(floatingLegTenor, a+b+1) );

	unsigned long seed = 47;
	RNGenerator_PTR  rnGenerator(new McGenerator(seed));
	MC_CheyetteDD_VanillaSwaptionPricer_PTR mc(
					new MC_CheyetteDD_VanillaSwaptionPricer(pCheyetteDD_Model, rnGenerator,
															simulationStructure, fwdProbaT,
															discretizationBetweenDates   ) ) ;

//approx
	double strike =  - 1000. ; //temporaire avant de mettre strike ATM
	VanillaSwap swap = VanillaSwap(strike, indexStart, indexEnd, floatingLegTenor, fixedLegTenor, simulationStructure);
	
	VanillaSwaption_PTR swaption_PTR_test(new VanillaSwaption(swap, OptionType::OptionType::CALL)) ;
	CheyetteDD_VanillaSwaptionApproxPricer approx(pCheyetteDD_Model, swaption_PTR_test);

//calcul du strike ATM pour le swap
	double strikeATM = approx.swapRate0() ;

	o << "strike ATM pour swaption " << a << "Y" << b << "Y : " << strikeATM << std::endl ;
	std::cout << "strike ATM pour swaption " << a << "Y" << b << "Y : " << strikeATM << std::endl ;
	std::cout << std::endl ;

	swap.set_strike(strikeATM) ;
	VanillaSwaption_PTR swaption_PTR_test_ATM(new VanillaSwaption(swap, OptionType::OptionType::CALL)) ;
	swaption_PTR_test_ATM->show() ;
	CheyetteDD_VanillaSwaptionApproxPricer approxATM(pCheyetteDD_Model, swaption_PTR_test_ATM);

	std::cout << "prixSwaption" << std::endl ;
	double approxPrice = approxATM.prixSwaptionApproxPiterbarg() ;
	double b_barre = approxATM.get_buffer_b_barre_() ;
	std::cout << "approxPrice : " << approxPrice << std::endl ;
	std::cout << "b barre :     " << b_barre << std::endl ;
	std::cout << "  " << std::endl ;

//MC
	size_t nbMC = nbSimus.size() ;
	std::vector<double> vectPrixMC(nbMC), vectICinf(nbMC), vectICsup(nbMC) ; 
	for (size_t i = 0 ; i < nbMC ; ++i)
	{
		std::vector<double> resMC = mc->price(swaption_PTR_test_ATM, nbSimus[i]) ;
		vectPrixMC[i] = resMC[0] ;
		vectICinf[i]  = resMC[1] ;
		vectICsup[i]  = resMC[2] ;
	}
	
	double annuityA0 = approxATM.swapRateDenominator(0., 0., 0.) ;
	double swapRateS0 = approxATM.swapRate0() ;

	double volBlack = NumericalMethods::Black_SwaptionImpliedVolatility(approxPrice, annuityA0, 
																		swapRateS0, strikeATM, a) ;

	std::cout << "strike :    " << strikeATM << std::endl ;
	std::cout << "expiry :    " << a << std::endl ;
	std::cout << "annuity :   " << annuityA0 << std::endl ;
	std::cout << "swapRate :  " << swapRateS0 << std::endl ;
	std::cout << "vol Black : " << volBlack << std::endl ;

	mc->printMC_vs_approx(o, approxPrice, b_barre, annuityA0, swapRateS0, volBlack,
							a, b, swaption_PTR_test_ATM, nbSimus, vectPrixMC, vectICinf, vectICsup) ;
}

/* *********************************/

void test_printMatrix(std::vector<size_t> expiry_maturity, double strike, Tenor floatingLegTenor, Tenor fixedLegTenor, 
					 int curveChoice, int shiftChoice, 
					 std::vector<double> x, std::vector<double> m_y, std::vector<double> sigma_y, double k, 
					 std::vector<size_t> nbSimus)
{	
//courbe spot
	CourbeInput_PTR courbe_PTR_test(createCourbeInput(curveChoice));
	
	Piecewiseconst_RR_Function sigma(x, sigma_y) ; 
	Piecewiseconst_RR_Function m(x, m_y) ;
	CheyetteDD_Model::CheyetteDD_Parameter monStruct(k, sigma, m) ;

	CheyetteDD_Model_PTR modele_test_PTR(new CheyetteDD_Model(courbe_PTR_test, monStruct, shiftChoice)) ;
	modele_test_PTR->show() ;

//swaption grid
	size_t dimMatrix = expiry_maturity.size() ;
	Matrice matriceApprox (dimMatrix, dimMatrix) ;
	Matrice matriceS0 (dimMatrix, dimMatrix) ;
	Matrice matriceBbarre (dimMatrix, dimMatrix) ;
	Matrice matriceVolBlack (dimMatrix, dimMatrix) ;
	Matrice matriceMCprice (dimMatrix, dimMatrix) ;
	Matrice matriceMC_ICinf (dimMatrix, dimMatrix) ;
	Matrice matriceMC_ICsup (dimMatrix, dimMatrix) ;
	Matrice matriceRelativeError (dimMatrix, dimMatrix) ;
	Matrice matriceAbsoluteError (dimMatrix, dimMatrix) ;
	Matrice matriceKtilde (dimMatrix, dimMatrix) ;
	for (size_t index1 = 0 ; index1 < expiry_maturity.size() ; ++index1)
	{
		for (size_t index2 = 0 ; index2 < expiry_maturity.size() ; ++index2)
		{
			test_printElement(strike, index1, index2, expiry_maturity, floatingLegTenor, fixedLegTenor, curveChoice, shiftChoice, 
									modele_test_PTR, nbSimus,
									matriceApprox, matriceS0, matriceBbarre, matriceVolBlack, 
									matriceMCprice, matriceMC_ICinf, matriceMC_ICsup, 
									matriceRelativeError, matriceAbsoluteError, matriceKtilde) ;	
		}	
		std::cout << "index 1 : " << index1 << std::endl ;
	}


	std::stringstream fileName_s ;
	std::string directory = LMMPATH::get_runtime_datapath() ;
	fileName_s << directory << "test.csv" ; 
	std::string fileName = fileName_s.str();

	ofstream o;
	o.open(fileName,  ios::out | ios::app );
	o	<<	endl;
	modele_test_PTR->print(o) ;
	o	<<	endl;
	o	<< "strike choisi : " << strike << endl ;
	switch (shiftChoice)
	{
		case 1 :
			o << "shift 1 : r(t) / r(0)" << endl ;
			break ;
		case 2 :
			o << "shift 2 : r(t) / f(0, t)" << endl ;
			break ;
		case 3 :
			o << "shift 3 : S(t) / S(0)" << endl ;
			break ;
		default:
			std::cout << "invalide shiftChoice" << std::endl ;
			throw "exception" ;
			break;
	}
	
	courbe_PTR_test->print(o) ;
	o << endl ;

	matrixPrint(o, expiry_maturity, "approximation price", matriceApprox) ;
	matrixPrint(o, expiry_maturity, "MC price - nbSimu ", matriceMCprice) ;
	matrixPrint(o, expiry_maturity, "Erreur relative", matriceRelativeError) ;
	matrixPrint(o, expiry_maturity, "Erreur absolue : MCprice - approx", matriceAbsoluteError) ;
	matrixPrint(o, expiry_maturity, "S0 (strike atm)", matriceS0) ;
	matrixPrint(o, expiry_maturity, "b barre", matriceBbarre) ;
	matrixPrint(o, expiry_maturity, "vol Black", matriceVolBlack) ;
	matrixPrint(o, expiry_maturity, "MC IC inf", matriceMC_ICinf) ;
	matrixPrint(o, expiry_maturity, "MC IC sup", matriceMC_ICsup) ;
	matrixPrint(o, expiry_maturity, "K tilde", matriceKtilde) ;

	o.close();
}

void test_printElement(double strike, size_t index1, size_t index2, std::vector<size_t> expiry_maturity, 
						  Tenor floatingLegTenor, Tenor fixedLegTenor, 
						 int curveChoice, int shiftChoice, 
						 CheyetteDD_Model_PTR modele_test_PTR, 
						 std::vector<size_t> nbSimus, 
						 Matrice& matriceApprox, Matrice& matriceS0, Matrice& matriceBbarre, Matrice& matriceVolBlack,
						 Matrice& matriceMCprice, Matrice& matriceMC_ICinf, Matrice& matriceMC_ICsup, 
						 Matrice& matriceRelativeError, Matrice& matriceAbsoluteError, Matrice& matriceKtilde)
{	
	double tenor = std::min(floatingLegTenor.YearFraction() , fixedLegTenor.YearFraction() ) ;
	
	size_t a = expiry_maturity[index1] ;
	size_t b = expiry_maturity[index2] ;

	size_t indexStart = size_t(a / tenor) ;
	size_t indexEnd = size_t((a+b) / tenor) ;
	size_t fwdProbaT = a + b ; 

//param MC
	LMMTenorStructure_PTR simulationStructure(new LMMTenorStructure(floatingLegTenor, a+b+1) );
	size_t	discretizationBetweenDates = 200 ;

	unsigned long seed = 47;
	RNGenerator_PTR  rnGenerator(new McGenerator(seed));
	MC_CheyetteDD_VanillaSwaptionPricer_PTR mc(
					new MC_CheyetteDD_VanillaSwaptionPricer(modele_test_PTR, rnGenerator,
															simulationStructure, fwdProbaT,
															discretizationBetweenDates   ) ) ;

	VanillaSwap swap = VanillaSwap(strike, indexStart, indexEnd, floatingLegTenor, fixedLegTenor, simulationStructure);
	
	VanillaSwaption_PTR swaption_PTR_test(new VanillaSwaption(swap, OptionType::OptionType::CALL)) ;
	CheyetteDD_VanillaSwaptionApproxPricer approx(modele_test_PTR, swaption_PTR_test);

	std::cout << "prixSwaption" << std::endl ;
	double approxPrice = approx.prixSwaptionApproxPiterbarg() ;
	double b_barre = approx.get_buffer_b_barre_() ;
	std::cout << "approxPrice : " << approxPrice << std::endl ;
	std::cout << "b barre :     " << b_barre << std::endl ;
	std::cout << "  " << std::endl ;
	matriceApprox (index1, index2)  = approxPrice ; 
	matriceBbarre (index1, index2)  = b_barre ; 

//MC
	//double nominal = 1 ;
	//LMMTenorStructure_CONSTPTR swapStructure(new LMMTenorStructure(floatingLegTenor, a+b+1)) ;
	//GeneticSwap_CONSTPTR genericSwapTest = InstrumentFactory::createVanillaSwap(strikeATM, indexStart, 
	//												indexEnd, floatingLegTenor, fixedLegTenor,
	//												swapStructure, nominal);

	//size_t maturity = indexStart ; //index
	//GeneticSwaption_CONSTPTR 	genericSwaptionTest(new GeneticSwaption(maturity, genericSwapTest)) ;

	size_t nbMC = nbSimus.size() ;
	std::vector<double> vectPrixMC(nbMC), vectICinf(nbMC), vectICsup(nbMC) ; 
	for (size_t i = 0 ; i < nbMC ; ++i)
	{
		//std::vector<double> resMC = mc->price(genericSwaptionTest, nbSimus[i]) ;
		std::vector<double> resMC = mc->price(swaption_PTR_test, nbSimus[i]) ;
		matriceMCprice (index1, index2)  = resMC[0] ;
		matriceMC_ICinf (index1, index2) = resMC[1] ;
		matriceMC_ICsup (index1, index2) = resMC[2] ;
	}

	double annuityA0 = approx.swapRateDenominator(0., 0., 0.) ;
	double swapRateS0 = approx.swapRate0() ;
	matriceS0(index1, index2) = swapRateS0 ;
	matriceKtilde(index1, index2) = b_barre * strike + (1 - b_barre) * swapRateS0 ;
	double volBlack = NumericalMethods::Black_SwaptionImpliedVolatility(approxPrice, annuityA0, 
																	swapRateS0, strike, a) ;
	matriceVolBlack (index1, index2) = volBlack ;
}

/**********************************/

//param 
//floating leg tenor vs fixed leg tenor (pas vraiment un tenor mais frequence des flux)
//size(x) = size(m_y) + 1 = size(sigma_y) + 1
void test_printMatrix_ATM(std::vector<size_t> expiry_maturity, Tenor floatingLegTenor, Tenor fixedLegTenor, 
					 int curveChoice, int shiftChoice, 
					 std::vector<double> x, std::vector<double> m_y, std::vector<double> sigma_y, double k, 
					 std::vector<size_t> nbSimus)
{	
//courbe spot
	CourbeInput_PTR courbe_PTR_test(createCourbeInput(curveChoice));
	
	Piecewiseconst_RR_Function sigma(x, sigma_y) ; 
	Piecewiseconst_RR_Function m(x, m_y) ;
	CheyetteDD_Model::CheyetteDD_Parameter monStruct(k, sigma, m) ;

	CheyetteDD_Model_PTR modele_test_PTR(new CheyetteDD_Model(courbe_PTR_test, monStruct, shiftChoice)) ;
	modele_test_PTR->show() ;

//swaption grid
	size_t dimMatrix = expiry_maturity.size() ;
	Matrice matriceApprox (dimMatrix, dimMatrix) ;
	Matrice matriceS0 (dimMatrix, dimMatrix) ;
	Matrice matriceBbarre (dimMatrix, dimMatrix) ;
	Matrice matriceVolBlack (dimMatrix, dimMatrix) ;
	Matrice matriceMCprice (dimMatrix, dimMatrix) ;
	Matrice matriceMC_ICinf (dimMatrix, dimMatrix) ;
	Matrice matriceMC_ICsup (dimMatrix, dimMatrix) ;
	Matrice matriceRelativeError (dimMatrix, dimMatrix) ;
	Matrice matriceAbsoluteError (dimMatrix, dimMatrix) ;
	Matrice matriceKtilde (dimMatrix, dimMatrix) ;
	for (size_t index1 = 0 ; index1 < expiry_maturity.size() ; ++index1)
	{
		for (size_t index2 = 0 ; index2 < expiry_maturity.size() ; ++index2)
		{
			test_printElement_ATM(index1, index2, expiry_maturity, floatingLegTenor, fixedLegTenor, curveChoice, shiftChoice, 
									modele_test_PTR, nbSimus,
									matriceApprox, matriceS0, matriceBbarre, matriceVolBlack, 
									matriceMCprice, matriceMC_ICinf, matriceMC_ICsup, 
									matriceRelativeError, matriceAbsoluteError, matriceKtilde) ;	
		}	
		std::cout << "index 1 : " << index1 << std::endl ;
	}


	std::stringstream fileName_s ;
	std::string directory = LMMPATH::get_runtime_datapath() ;
	fileName_s << directory << "test_shift.csv" ; 
	std::string fileName = fileName_s.str();

	ofstream o;
	o.open(fileName,  ios::out | ios::app );
	o	<<	endl;
	modele_test_PTR->print(o) ;
	o	<<	endl;
	switch (shiftChoice)
	{
		case 1 :
			o << "shift 1 : r(t) / r(0)" << endl ;
			break ;
		case 2 :
			o << "shift 2 : r(t) / f(0, t)" << endl ;
			break ;
		case 3 :
			o << "shift 3 : S(t) / S(0)" << endl ;
			break ;
		default:
			std::cout << "invalide shiftChoice" << std::endl ;
			throw "exception" ;
			break;
	}
	
	courbe_PTR_test->print(o) ;
	o << endl ;
	

	matrixPrint(o, expiry_maturity, "approximation price ATM", matriceApprox) ;
	matrixPrint(o, expiry_maturity, "MC price - nbSimu ", matriceMCprice) ;
	matrixPrint(o, expiry_maturity, "Erreur relative", matriceRelativeError) ;
	matrixPrint(o, expiry_maturity, "Erreur absolue : MCprice - approx", matriceAbsoluteError) ;
	matrixPrint(o, expiry_maturity, "S0 (strike atm)", matriceS0) ;
	matrixPrint(o, expiry_maturity, "b barre", matriceBbarre) ;
	matrixPrint(o, expiry_maturity, "vol Black", matriceVolBlack) ;
	matrixPrint(o, expiry_maturity, "MC IC inf", matriceMC_ICinf) ;
	matrixPrint(o, expiry_maturity, "MC IC sup", matriceMC_ICsup) ;
	matrixPrint(o, expiry_maturity, "K tilde", matriceKtilde) ;

	o.close();
}

void matrixPrint(std::ofstream& o, std::vector<size_t> expiry_maturity, std::string chaine, const Matrice& M)
{
	o << endl ;
	o << endl ;
	o << chaine << endl ;
	o << "expiry (a) \\ tenor swap (b) ; " ;
		for (size_t i = 0 ; i < expiry_maturity.size() ; i++)
		{
			o << expiry_maturity[i] << "Y ;" ;
		}
		o	<<	endl;
		MatrixPrintElement<Matrice> matrixPrinter("azzrehrj", M) ;
		for (size_t rowIndex = 0 ; rowIndex < expiry_maturity.size() ; rowIndex++)
		{
			o << expiry_maturity[rowIndex] << "Y ;" ;
			matrixPrinter.print_element2(o, rowIndex) ;
			o << endl ;
		}
}

void test_printElement_ATM(size_t index1, size_t index2, std::vector<size_t>& expiry_maturity, 
						  Tenor floatingLegTenor, Tenor fixedLegTenor, 
						 int curveChoice, int shiftChoice, 
						 CheyetteDD_Model_PTR modele_test_PTR, 
						 std::vector<size_t>& nbSimus, 
						 Matrice& matriceApprox, Matrice& matriceS0, Matrice& matriceBbarre, Matrice& matriceVolBlack,
						 Matrice& matriceMCprice, Matrice& matriceMC_ICinf, Matrice& matriceMC_ICsup, 
						 Matrice& matriceRelativeError, Matrice& matriceAbsoluteError, Matrice& matriceKtilde)
{	
	double tenor = std::min(floatingLegTenor.YearFraction() , fixedLegTenor.YearFraction() ) ;
	
	size_t a = expiry_maturity[index1] ;
	size_t b = expiry_maturity[index2] ;

	size_t indexStart = size_t(a / tenor) ;
	size_t indexEnd = size_t((a+b) / tenor) ;

//approx
	double strike = 1000. ; //temporaire avant de mettre strike ATM
	LMMTenorStructure_PTR simulationStructure(new LMMTenorStructure(floatingLegTenor, a+b+1) );
	VanillaSwap swap = VanillaSwap(strike, indexStart, indexEnd, floatingLegTenor, fixedLegTenor, simulationStructure);
	
	VanillaSwaption_PTR swaption_PTR_test(new VanillaSwaption(swap, OptionType::OptionType::CALL)) ;
	CheyetteDD_VanillaSwaptionApproxPricer approx(modele_test_PTR, swaption_PTR_test);

//calcul du strike ATM pour le swap
	double strikeATM = approx.swapRate0() ;
	swap.set_strike(strikeATM) ;
	VanillaSwaption_PTR swaption_PTR_test_ATM(new VanillaSwaption(swap, OptionType::OptionType::CALL)) ;
	CheyetteDD_VanillaSwaptionApproxPricer approxATM(modele_test_PTR, swaption_PTR_test_ATM);


//param MC
	size_t fwdProbaT = a + b ;
	size_t	discretizationBetweenDates = 200 ;

	unsigned long seed = 47;
	RNGenerator_PTR  rnGenerator(new McGenerator(seed));
	MC_CheyetteDD_VanillaSwaptionPricer_PTR mc(
					new MC_CheyetteDD_VanillaSwaptionPricer(modele_test_PTR, rnGenerator,
															simulationStructure, fwdProbaT,
															discretizationBetweenDates   ) ) ;


	std::cout << "prixSwaption" << std::endl ;
	double approxPrice = approxATM.prixSwaptionApproxPiterbarg() ;
	double b_barre = approxATM.get_buffer_b_barre_() ;
	std::cout << "approxPrice : " << approxPrice << std::endl ;
	std::cout << "b barre :     " << b_barre << std::endl ;
	std::cout << "  " << std::endl ;
	matriceApprox (index1, index2)  = approxPrice ; 
	matriceBbarre (index1, index2)  = b_barre ; 

//MC
	//double nominal = 1 ;
	//LMMTenorStructure_CONSTPTR swapStructure(new LMMTenorStructure(floatingLegTenor, a+b+1)) ;
	//GeneticSwap_CONSTPTR genericSwapTest = InstrumentFactory::createVanillaSwap(strikeATM, indexStart, 
	//												indexEnd, floatingLegTenor, fixedLegTenor,
	//												swapStructure, nominal);

	//size_t maturity = indexStart ; //index
	//GeneticSwaption_CONSTPTR 	genericSwaptionTest(new GeneticSwaption(maturity, genericSwapTest)) ;

	size_t nbMC = nbSimus.size() ;
	std::vector<double> vectPrixMC(nbMC), vectICinf(nbMC), vectICsup(nbMC) ; 
	for (size_t i = 0 ; i < nbMC ; ++i)
	{
		//std::vector<double> resMC = mc->price(genericSwaptionTest_ATM, nbSimus[i]) ;
		std::vector<double> resMC = mc->price(swaption_PTR_test_ATM, nbSimus[i]) ;
		vectPrixMC[i] = resMC[0] ;
		vectICinf[i]  = resMC[1] ;
		vectICsup[i]  = resMC[2] ;
		matriceMCprice (index1, index2)  = resMC[0] ;
		matriceMC_ICinf (index1, index2) = resMC[1] ;
		matriceMC_ICsup (index1, index2) = resMC[2] ;
	}

	matriceAbsoluteError(index1, index2) = matriceMCprice(index1, index2) - matriceApprox(index1, index2) ;
	matriceRelativeError(index1, index2) = matriceAbsoluteError(index1, index2) / matriceMCprice(index1, index2) ;

	double annuityA0 = approxATM.swapRateDenominator(0., 0., 0.) ;
	double swapRateS0 = approxATM.swapRate0() ;
	matriceS0(index1, index2) = swapRateS0 ;

	double volBlack = NumericalMethods::Black_SwaptionImpliedVolatility(approxPrice, annuityA0, 
																	swapRateS0, strikeATM, a) ;
	matriceVolBlack (index1, index2) = volBlack ;
	matriceKtilde (index1, index2) = b_barre * strikeATM + (1 - b_barre) * swapRateS0 ; 
	//mc->printMC_vs_approx(approxPrice, b_barre, annuityA0, swapRateS0, blackPrice,
	//						a, b, genericSwaptionTest, nbSimus, vectPrixMC, vectICinf, vectICsup) ;

}


//swaption aY bY... 
void smile(size_t a, size_t b, double strikeATM, Tenor floatingLegTenor, Tenor fixedLegTenor, 
			int curveChoice, int shiftChoice, 
			std::vector<double> x, std::vector<double> m_y, std::vector<double> sigma_y, double k)
{
	size_t nbStrikes = 17 ;
	std::vector<double> ATMpourcent(nbStrikes) ; 
	ATMpourcent[0] = 0.7 ; 
	ATMpourcent[1] = 0.75 ; 
	ATMpourcent[2] = 0.8 ; 
	ATMpourcent[3] = 0.85 ; 
	ATMpourcent[4] = 0.9 ; 
	ATMpourcent[5] = 0.92 ;
	ATMpourcent[6] = 0.94 ;
	ATMpourcent[7] = 0.96 ;
	ATMpourcent[8] = 0.98 ;
	ATMpourcent[9] = 1.0 ;  
	ATMpourcent[10] = 1.02 ;  
	ATMpourcent[11] = 1.04 ; 
	ATMpourcent[12] = 1.06 ; 
	ATMpourcent[13] = 1.1 ; 
	ATMpourcent[14] = 1.15 ;
	ATMpourcent[15] = 1.2 ;
	ATMpourcent[16] = 1.25 ;

//courbe spot
	CourbeInput_PTR courbe_PTR_test(createCourbeInput(curveChoice));
	
	Piecewiseconst_RR_Function sigma(x, sigma_y) ; 
	Piecewiseconst_RR_Function m(x, m_y) ;
	CheyetteDD_Model::CheyetteDD_Parameter monStruct(k, sigma, m) ;

	CheyetteDD_Model_PTR modele_test_PTR(new CheyetteDD_Model(courbe_PTR_test, monStruct, shiftChoice)) ;
	modele_test_PTR->show() ;
	
	std::stringstream fileName_s ;
	std::string directory = LMMPATH::get_runtime_datapath() ;
	fileName_s << directory << "smile2.csv" ; 
	std::string fileName = fileName_s.str();

	ofstream o;
	o.open(fileName,  ios::out | ios::app );
	o	<<	endl;
	modele_test_PTR->print(o) ;
	switch (shiftChoice)
	{
		case 1 :
			o << "shift 1 : r(t) / r(0)" << endl ;
			break ;
		case 2 :
			o << "shift 2 : r(t) / f(0, t)" << endl ;
			break ;
		case 3 :
			o << "shift 3 : S(t) / S(0)" << endl ;
			break ;
		default:
			std::cout << "invalide shiftChoice" << std::endl ;
			throw "exception" ;
			break;
	}
	
	courbe_PTR_test->print(o) ;
	o << endl ;

	double tenor = std::min(floatingLegTenor.YearFraction() , fixedLegTenor.YearFraction() ) ;
	
	size_t indexStart = size_t(a / tenor) ;
	size_t indexEnd = size_t((a+b) / tenor) ;
	double fwdProbaT = b ;  //size_t vers double ok
	LMMTenorStructure_PTR simulationStructure(new LMMTenorStructure(floatingLegTenor, a+b+1) );
	
	o << "moneyness ; strike ; volBlack ; b_barre ; approxPrice" << endl ;
	for (size_t i = 0 ; i < ATMpourcent.size() ; ++i)
	{
		VanillaSwap swap = VanillaSwap(ATMpourcent[i] * strikeATM, indexStart, indexEnd, floatingLegTenor, fixedLegTenor, simulationStructure);
		VanillaSwaption_PTR swaption_PTR_test(new VanillaSwaption(swap, OptionType::OptionType::CALL)) ;

		CheyetteDD_VanillaSwaptionApproxPricer approx(modele_test_PTR, swaption_PTR_test);

		double approxPrice = approx.prixSwaptionApproxPiterbarg() ;
		double b_barre = approx.get_buffer_b_barre_() ;

		double annuityA0 = approx.swapRateDenominator(0., 0., 0.) ;
		double swapRateS0 = approx.swapRate0() ;
		double volBlack = NumericalMethods::Black_SwaptionImpliedVolatility(approxPrice, annuityA0, 
																		swapRateS0, ATMpourcent[i] * strikeATM, a) ;
		o << ATMpourcent[i] << ";" << ATMpourcent[i] * strikeATM << ";" << volBlack << ";" 
																		<< b_barre << ";" << approxPrice << endl ;
	}

	o.close();
						 
}