#pragma once

#include <Cheyette/Pricer/MC_Cheyette.h>
#include <Cheyette/Pricer/MC_CheyetteDD_GenericSwapPricer.h>
#include <Cheyette/Pricer/MC_CheyetteDD_GenericSwaptionPricer.h>
#include <Cheyette/Pricer/MC_CheyetteDD_VanillaSwaptionPricer.h>

#include <Cheyette/Fonction.h>
#include <Cheyette/Pricer/CheyetteDD_VanillaSwaptionApproxPricer.h>
#include <Cheyette/CheyetteModel/CourbeInput.h>

#include <LMM/numeric/Integrator1D.h>
#include <JBLMM/Instrument/InstrumentFactory.h>

#include <iostream>
#include <vector>
#include <ctime>
#include <string>

CourbeInput_PTR createCourbeInput(int curveChoice) ;

void test_approx(double strike, size_t a, size_t b, Tenor floatingLegTenor, Tenor fixedLegTenor, 
				 int curveChoice, int shiftChoice, 
				 std::vector<double> x, std::vector<double> m_y, std::vector<double> sigma_y, double k, 
				 std::vector<size_t> nbSimus) ;

void test_approx_ATM(size_t a, size_t b, Tenor floatingLegTenor, Tenor fixedLegTenor, 
					 int curveChoice, int shiftChoice, 
					 std::vector<double> x, std::vector<double> m_y, std::vector<double> sigma_y, double k, 
					 std::vector<size_t> nbSimus) ;
