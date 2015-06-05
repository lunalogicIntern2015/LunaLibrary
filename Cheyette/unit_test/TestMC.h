#pragma once

#include <Cheyette/Pricer/MC_Cheyette.h>
#include <Cheyette/Pricer/MC_CheyetteDD_VanillaSwapPricer.h>
#include <Cheyette/Pricer/MC_CheyetteDD_VanillaSwaptionPricer.h>
#include <Cheyette/Pricer/MC_CheyetteDD_ZCPricer.h>
#include <Cheyette/Fonction.h>
#include <LMM/RNGenerator/RNGenerator.h>
#include <LMM/RNGenerator/McGenerator.h>

#include <iostream>
#include <vector>


void UneTrajectoireEuler() ;

void TestMCSwapPricer() ;

//double ZCVasicek(double matuT, double mean_rev, double level, double sigma) ;