#pragma once

#include <Cheyette\CheyetteModel\CourbeInput.h>	
#include <LMM/instrument/VanillaSwaption.h>
#include <Cheyette/Fonction.h>
#include <Cheyette/Pricer/CheyetteDD_VanillaSwaptionApproxPricer.h>
#include <Cheyette/CheyetteModel/CheyetteDD_Model.h>
#include <Cheyette/CheyetteModel/CourbeInput.h>
#include <LMM/numeric/Integrator1D.h>

#include <iostream>
#include <vector>

CourbeInput_PTR createCourbeInput() ;
CheyetteDD_Model_PTR createCheyetteDD_Model() ;
VanillaSwaption_PTR createSwap();
void testSwap() ; //test swap (dates flux fixes et flottants) 
VanillaSwaption_PTR createSwaption() ;
CheyetteDD_VanillaSwaptionApproxPricer_PTR createApproxPricer_PTR() ;

void test_y_barre() ; 

//test des derivees de ZC et swapRateNumerator, swapRateDenominator
void test_ZC_swapRate_Num_Denom() ;

void test_fonction_inverse() ;  //necessite swapRate teste avant

void test_derivatives() ;

void test_time_average() ;

void test_y_bar_cas_limite() ;
