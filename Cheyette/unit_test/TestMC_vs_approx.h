#pragma once

#include <Cheyette/Pricer/MC_Cheyette.h>
#include <Cheyette/Pricer/MC_CheyetteDD_GenericSwapPricer.h>
#include <Cheyette/Pricer/MC_CheyetteDD_GenericSwaptionPricer.h>
//#include <Cheyette/Pricer/MC_CheyetteDD_VanillaSwaptionPricer.h>
#include <Cheyette/Pricer/MC_CheyetteDD_ZCPricer.h>
#include <Cheyette/Fonction.h>
#include <LMM/RNGenerator/RNGenerator.h>
#include <LMM/RNGenerator/McGenerator.h>

#include <Cheyette\CheyetteModel\CourbeInput.h>	
#include <LMM/instrument/VanillaSwaption.h>
#include <Cheyette/Pricer/CheyetteDD_VanillaSwaptionApproxPricer.h>
#include <Cheyette/CheyetteModel/CheyetteDD_Model.h>
#include <Cheyette/CheyetteModel/CourbeInput.h>
#include <LMM/numeric/Integrator1D.h>

#include <iostream>
#include <vector>
#include <ctime>
#include <string>

void test_approx() ;
