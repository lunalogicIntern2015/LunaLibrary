#pragma once

#include <Cheyette/Pricer/MC_Cheyette.h>
#include <Cheyette/Pricer/CheyetteDD_VanillaSwaptionApproxPricer.h>
#include <Cheyette/Pricer/MC_CheyetteDD_VanillaSwapPricer.h>
#include <Cheyette/Pricer/MC_CheyetteDD_VanillaSwaptionPricer.h>
#include <Cheyette/Pricer/MC_CheyetteDD_ZCPricer.h>
#include <Cheyette/Fonction.h>
#include <LMM/RNGenerator/RNGenerator.h>
#include <LMM/RNGenerator/McGenerator.h>

#include <Cheyette/unit_test/TestApproxDD.h>

#include <LMM/helper/Printer.h>

#include <iostream>
#include <vector>

CourbeInput_PTR createCourbeInput(int curveChoice) ;

CheyetteDD_VanillaSwaptionApproxPricer_PTR creeModeleATM(size_t a, size_t b, 
														 const Tenor& floatTenor, const Tenor& fixedTenor) ;

CheyetteDD_VanillaSwaptionApproxPricer_PTR creeModele(size_t a, size_t b, 
														const Tenor& floatTenor, const Tenor& fixedTenor,
														double strike) ;

void TestMCSwapPricer(size_t a, size_t b, size_t nbSimus, const Tenor& floatTenor, const Tenor& fixedTenor) ;
void TestMCSwapPricer_annuity(size_t a, size_t b, size_t nbSimus, const Tenor& floatTenor, const Tenor& fixedTenor) ;


void MCannuity(size_t a, size_t b, size_t nbSimus, const Tenor& floatTenor, const Tenor& fixedTenor) ; 
void MCforward(size_t a, size_t b, size_t nbSimus, const Tenor& floatTenor, const Tenor& fixedTenor) ;


void test_approx(double strike, size_t a, size_t b, Tenor floatingLegTenor, Tenor fixedLegTenor, 
				 int curveChoice, int shiftChoice, 
				 std::vector<double> x, std::vector<double> m_y, std::vector<double> sigma_y, double k, 
				 std::vector<size_t> nbSimus) ;

void test_approx_ATM(size_t a, size_t b, Tenor floatingLegTenor, Tenor fixedLegTenor, 
					 int curveChoice, int shiftChoice, 
					 std::vector<double> x, std::vector<double> m_y, std::vector<double> sigma_y, double k, 
					 std::vector<size_t> nbSimus) ;

void test_approx_ATM(size_t a, size_t b, Tenor floatingLegTenor, Tenor fixedLegTenor, 
					 CheyetteDD_Model_PTR pCheyetteDD_Model,
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

void test_approx_ATM(size_t a, size_t b, Tenor floatingLegTenor, Tenor fixedLegTenor, 
					 CheyetteDD_Model_PTR pCheyetteDD_Model,
					 std::vector<size_t> nbSimus, std::ofstream& o) ;