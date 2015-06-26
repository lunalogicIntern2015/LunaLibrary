#pragma once

#include <Cheyette/Pricer/MC_Cheyette.h>
#include <Cheyette/Pricer/MC_CheyetteDD_VanillaSwapPricer.h>
#include <Cheyette/Pricer/MC_CheyetteDD_GenericSwapPricer.h>
#include <Cheyette/Pricer/MC_CheyetteDD_GenericSwaptionPricer.h>
//#include <Cheyette/Pricer/MC_CheyetteDD_VanillaSwaptionPricer.h>
#include <Cheyette/Pricer/MC_CheyetteDD_ZCPricer.h>
#include <Cheyette/Fonction.h>
#include <LMM/RNGenerator/RNGenerator.h>
#include <LMM/RNGenerator/McGenerator.h>

#include <JBLMM/Instrument/InstrumentFactory.h>

#include <iostream>
#include <vector>
#include <ctime>
#include <string>

//void print() ;

void testgenericSwap() ;

void testgenericSwaption() ;

void testVanillaSwap() ;

void testVanillaSwaption() ;