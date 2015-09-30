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

#include <Cheyette/unit_test/Test_comparaison_MC_approx.h>
#include <Cheyette/unit_test/TestCalibrator.h>
#include <Cheyette/unit_test/Test_LectureFichier.h>
#include <Cheyette/Calibration/CheyetteMarketData_2.h>

int main()
{
//	lancementQualiteApprox() ;

	size_t fileNumber = 3 ; //2008 07 02
//	size_t fileNumber = 14 ; //2011 04 06
//	size_t fileNumber = 37 ; //2014_08_08

	size_t coterminal = 11 ;
	printAllResults_calibratedData(fileNumber, coterminal) ;

	//fileNumber = 14 ;
	//printAllResults_calibratedData(fileNumber, coterminal) ;

	//size_t fileNumber = 37 ;
	//printAllResults_calibratedData(fileNumber, coterminal) ;

	//testRecupData(0) ;
	//testReadFile(0) ;
	//test_MC_approx() ;
	//lancementCalibOneFile() ;

	//std::vector<double> v_sigma(3, 1) ; //3 fois l'élément 1
	//std::cout << v_sigma[0] << v_sigma[1] << v_sigma[2] << std::endl ;

	//test_approx() ;

	//testMCSwap_SwaptionPricer() ;
//	lancementAuto() ;


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


	getchar();
}