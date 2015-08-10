#include <LMM/Test/Tests.h>
#include <stdio.h>

#include <LMM/Pricer/LmmApproximationPricer/LmmVanillaSwaptionApproxPricer_Rebonato.h>
#include <Numeric/NumericalMethods.h>

#include <Cheyette/unit_test/TestFonction.h>
#include <Cheyette/unit_test/TestIntegrator1D_2D.h>
#include <Cheyette/unit_test/TestApproxDD.h>
#include <Cheyette/Model/CheyetteDD_Model.h>
#include <Numeric/Integrator1D.h>
#include <Cheyette/unit_test/TestMC.h>



int main()
{

	//lancementAuto() ;

	//size_t fileNumber = 1 ;
	//size_t coterminal = 12 ;
	//printAllResults_calibratedData(fileNumber, coterminal) ;

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
	//size_t numfile = 24 ;
	//size_t coterminal = 12 ;
	//testCalib(numfile, coterminal) ;   //2004 04 07 : file 27


	//swap MC contre courbe
	//////TestMCSwapPricer(a, b, simus, floatingLegTenor, fixedLegTenor) ;
	//////TestMCSwapPricer_annuity(a, b, simus, floatingLegTenor, fixedLegTenor) ;

	//double bs_call_price = 0.0747381/100 ; //0.0611958/100 ;
	//double annuity0 = 83.86/100  ; //A0
	//double s0 = 2.6342/100 ; //s0
	//double strike = s0 ;   
	//double T = 10.0 ;
	//std::cout << NumericalMethods::Black_SwaptionImpliedVolatility(bs_call_price, annuity0, s0, strike, T) ;

	getchar();
}