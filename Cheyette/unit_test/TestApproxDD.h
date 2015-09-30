#pragma once

#include <iostream>
#include <vector>

#include "TestMC.h"

#include <Instrument/VanillaSwaption.h>
#include <Cheyette/Fonction.h>
#include <Numeric/Integrator1D.h>

#include <Cheyette/Model/CourbeInput.h>	
#include <Cheyette/Pricer/Cheyette_SwaptionPricer_LinearApprox.h>
#include <Cheyette/Model/CheyetteModel.h>
#include <Cheyette/unit_test/Test_CheyetteDD_Model.h>
#include <fstream>

VanillaSwap_PTR createSwapTest() ;
VanillaSwaption_PTR createSwaptionTest();

void testSwap() ; //test swap (dates flux fixes et flottants) 

VanillaSwaption_PTR createSwaption(double strike, LMM::Index  indexStart, LMM::Index  indexEnd, 
									Tenor floatingLegTenorType, Tenor fixedLegTenorType) ;

VanillaSwaption_PTR setSwaptionATM_DD(	const CheyetteDD_Model_PTR modele_test_PTR,
										const Tenor& floatTenor, const Tenor& fixedTenor, 
										const size_t a, const size_t b) ;

//modele
	Cheyette_SwaptionPricer_LinearApprox_PTR createApproxPricer_PTR(size_t xmax, int numCourbe, 
																	double k, double sigmaValue, double mValue) ;

	Cheyette_SwaptionPricer_LinearApprox_PTR createApproxPricer_PTR(size_t xmax, CourbeInput_PTR pCourbeInput, 
																	double k, double sigmaValue, double mValue) ;

//test qualite approximation prix swaption
	void testQualiteApprox(size_t xmax, int numCourbe, double k, double sigmaValue, double mValue) ;
	void lancementQualiteApprox() ;

//autres tests
	void test_y_barre() ; 

	//test des derivees de ZC et swapRateNumerator, swapRateDenominator
	void test_ZC_swapRate_Num_Denom() ;

	void test_swapRate_inverse() ;

	void test_fonction_inverse() ;  //necessite swapRate teste avant

	void test_derivatives() ;

	void test_time_average() ;

