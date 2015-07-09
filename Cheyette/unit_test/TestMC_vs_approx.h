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

#include <LMM/helper/Printer.h>
#include <boost/numeric/ublas/matrix.hpp>

CourbeInput_PTR createCourbeInput(int curveChoice) ;

void test_approx(double strike, size_t a, size_t b, Tenor floatingLegTenor, Tenor fixedLegTenor, 
				 int curveChoice, int shiftChoice, 
				 std::vector<double> x, std::vector<double> m_y, std::vector<double> sigma_y, double k, 
				 std::vector<size_t> nbSimus) ;

void test_approx_ATM(size_t a, size_t b, Tenor floatingLegTenor, Tenor fixedLegTenor, 
					 int curveChoice, int shiftChoice, 
					 std::vector<double> x, std::vector<double> m_y, std::vector<double> sigma_y, double k, 
					 std::vector<size_t> nbSimus) ;

void test_printMatrix_ATM(std::vector<size_t> expiry_maturity, Tenor floatingLegTenor, Tenor fixedLegTenor, 
					 int curveChoice, int shiftChoice, 
					 std::vector<double> x, std::vector<double> m_y, std::vector<double> sigma_y, double k, 
					 std::vector<size_t> nbSimus) ;

void matrixPrint(std::ofstream& o, std::vector<size_t> expiry_maturity, std::string chaine, const Matrice& M) ;

void test_printElement_ATM(size_t index1, size_t index2, std::vector<size_t>& expiry_maturity, 
						  Tenor floatingLegTenor, Tenor fixedLegTenor, 
						 int curveChoice, int shiftChoice, 
						 CheyetteDD_Model_PTR modele_test_PTR, 
						 std::vector<size_t>& nbSimus, 
						 Matrice& matriceApprox, Matrice& matriceS0, Matrice& matriceBbarre, Matrice& matriceVolBlack,
						 Matrice& matriceMCprice, Matrice& matriceMC_ICinf, Matrice& matriceMC_ICsup, 
						 Matrice& matriceRelativeError, Matrice& matriceAbsoluteError, Matrice& matriceKtilde) ;

void test_printMatrix(std::vector<size_t> expiry_maturity, double strike, Tenor floatingLegTenor, Tenor fixedLegTenor, 
					 int curveChoice, int shiftChoice, 
					 std::vector<double> x, std::vector<double> m_y, std::vector<double> sigma_y, double k, 
					 std::vector<size_t> nbSimus) ;

void test_printElement(double strike, size_t index1, size_t index2, std::vector<size_t> expiry_maturity, 
						  Tenor floatingLegTenor, Tenor fixedLegTenor, 
						 int curveChoice, int shiftChoice, 
						 CheyetteDD_Model_PTR modele_test_PTR, 
						 std::vector<size_t> nbSimus, 
						 Matrice& matriceApprox, Matrice& matriceS0, Matrice& matriceBbarre, Matrice& matriceVolBlack,
						 Matrice& matriceMCprice, Matrice& matriceMC_ICinf, Matrice& matriceMC_ICsup, 
						 Matrice& matriceRelativeError, Matrice& matriceAbsoluteError, Matrice& matriceKtilde) ;

void smile(size_t a, size_t b, double strikeATM, Tenor floatingLegTenor, Tenor fixedLegTenor, 
			int curveChoice, int shiftChoice, 
			std::vector<double> x, std::vector<double> m_y, std::vector<double> sigma_y, double k) ;
