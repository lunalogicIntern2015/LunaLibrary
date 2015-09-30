#include "Test_CheyetteDD_Model.h"
#include <iostream>


//courbe 0 : courbe plate à 1%
//courbe 1 : courbe test 
//courbe 2 : courbe marché interbancaire du 22-06-15
CourbeInput_PTR createCourbeInput(int curveChoice)
{
	switch (curveChoice)
	{
	case 1:{
		std::vector<double> listeMatu, tauxZC ;
		double translation = 3./100 ; //0. ;
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



CheyetteDD_Model_PTR creeCheyetteDD_Modele_PTR(size_t xmax, int numCourbe, 
											   double k, double sigmaValue, double mValue)
{

	std::vector<double> x(xmax) ;				//a+b + 1
	std::vector<double> m_y(xmax - 1) ;			//a+b
	std::vector<double> sigma_y(xmax - 1) ;		//a+b
	for (size_t i = 0 ; i <= xmax - 1 ; ++i)
	{
		x[i] = i ;
	}
	for (size_t i = 0 ; i < xmax - 1 ; ++i)
	{
		m_y[i] = mValue ;
		sigma_y[i] = sigmaValue ;
	}
	Piecewiseconst_RR_Function sigma(x, sigma_y) ; 
	Piecewiseconst_RR_Function m(x, m_y) ;

	CourbeInput_PTR courbe_PTR_test(createCourbeInput(numCourbe));

//3 versions du modele Cheyette Displaced Diffusion

	CheyetteDD_Model_PTR modele_test_PTR1(new CheyetteDD_Model_ShortRateVersion(courbe_PTR_test, k, sigma, m)) ;
	CheyetteDD_Model_PTR modele_test_PTR2(new CheyetteDD_Model_ForwardRateVersion(courbe_PTR_test, k, sigma, m)) ;
	
	//VanillaSwaption_PTR pSwaption = createSwaptionTest() ;
	//CheyetteDD_Model_PTR modele_test_PTR3(new CheyetteDD_Model_SwapRateVersion(numCourbe, k, sigma, m, pSwaption)) ;

	return modele_test_PTR1 ;
}

CheyetteDD_Model_PTR creeCheyetteDD_Modele_PTR(size_t xmax, CourbeInput_PTR pCourbeInput, 
											   double k, double sigmaValue, double mValue)
{

	std::vector<double> x(xmax) ;				//a+b + 1
	std::vector<double> m_y(xmax - 1) ;			//a+b
	std::vector<double> sigma_y(xmax - 1) ;		//a+b
	for (size_t i = 0 ; i <= xmax - 1 ; ++i)
	{
		x[i] = i ;
	}
	for (size_t i = 0 ; i < xmax - 1 ; ++i)
	{
		m_y[i] = mValue ;
		sigma_y[i] = sigmaValue ;
	}
	Piecewiseconst_RR_Function sigma(x, sigma_y) ; 
	Piecewiseconst_RR_Function m(x, m_y) ;

//3 versions du modele Cheyette Displaced Diffusion

	CheyetteDD_Model_PTR modele_test_PTR1(new CheyetteDD_Model_ShortRateVersion(pCourbeInput, k, sigma, m)) ;
	CheyetteDD_Model_PTR modele_test_PTR2(new CheyetteDD_Model_ForwardRateVersion(pCourbeInput, k, sigma, m)) ;

	VanillaSwaption_PTR pSwaption = createSwaptionTest() ;
	CheyetteDD_Model_PTR modele_test_PTR3(new CheyetteDD_Model_SwapRateVersion(pCourbeInput, k, sigma, m, pSwaption)) ;
	
	return modele_test_PTR1 ;
}