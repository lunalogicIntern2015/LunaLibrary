#pragma once

#include <Cheyette/Pricer/MC_Cheyette.h>
#include <Cheyette/unit_test/Test_CheyetteDD_Model.h>
#include <Cheyette/unit_test/Test_CheyetteQuad_Model.h>

#include <Cheyette/unit_test/TestApproxDD.h>
#include <Cheyette/unit_test/TestApproxCall.h>

#include <Cheyette/Pricer/Cheyette_SwaptionPricer_Approx.h>
#include <Cheyette/Pricer/MC_Cheyette_VanillaSwapPricer.h>
#include <Cheyette/Pricer/MC_Cheyette_VanillaSwaptionPricer.h>

//#include <Cheyette/Pricer/MC_Cheyette_ZCPricer.h>
#include <Cheyette/Fonction.h>

//#include <Cheyette/unit_test/TestCalibrator.h>

#include <LMM/helper/Printer.h>

#include <iostream>
#include <vector>

MC_Cheyette_VanillaSwapPricer_PTR creeMC_SwapPricer_PTR(CheyetteModel_PTR		cheyetteModel_PTR, 
															LMMTenorStructure_PTR	pTenorStructure,
															size_t					fwdProbaT,
															size_t					discretizationBetweenDates) ;

MC_Cheyette_VanillaSwaptionPricer_PTR creeMC_SwaptionPricer_PTR(CheyetteModel_PTR		cheyetteModel_PTR, 
															LMMTenorStructure_PTR	pTenorStructure,
															size_t					fwdProbaT,
															size_t					discretizationBetweenDates) ;

void testMCSwapPricer() ;

//void testMCSwapPricer_annuity() ;

void testMCSwap_SwaptionPricer() ;

//approx : show / print


//
////MONTE CARLO
//
//std::vector<double> MCforward(size_t nbSimus, CheyetteDD_VanillaSwaptionApproxPricer_PTR pApproxPricer) ;
//
//std::vector<std::vector<double>> MCforward_MultipleStrikes(size_t nbSimus, 
//															CheyetteDD_VanillaSwaptionApproxPricer_PTR pApproxPricer, 
//															double	sigma_ATM) ;
//
//void print_MCforward(size_t nbSimus, CheyetteDD_VanillaSwaptionApproxPricer_PTR pApproxPricer, std::ofstream& o) ;
//
////MC EQUATION 1
//std::vector<double> MCannuity(size_t nbSimus, CheyetteDD_VanillaSwaptionApproxPricer_PTR pApproxPricer) ;
//
//price swaptions pour différents strikes avec un seul jeu de simulations pour les differents strikes
//std::vector<std::vector<double>> MCannuityMultipleStrikes(size_t fileNumber, size_t a, size_t b,
//												  double strikeATM_Bloomberg, size_t nbSimus) ;
//
//std::vector<std::vector<double>> MCannuity_2_MultipleStrikes(  size_t fileNumber, size_t a, size_t b,
//															double strikeATM_Bloomberg, size_t nbSimus) ;
//
//void print_MCannuity(size_t nbSimus, CheyetteDD_VanillaSwaptionApproxPricer_PTR pApproxPricer, std::ofstream& o) ;
//
////MC EQUATION 2
//std::vector<double> MCannuity2(size_t nbSimus, CheyetteDD_VanillaSwaptionApproxPricer_PTR pApproxPricer) ;
//
//void print_MCannuity2(size_t nbSimus, CheyetteDD_VanillaSwaptionApproxPricer_PTR pApproxPricer, std::ofstream& o) ;
//
////MC EQUATION 2
//void MCannuity3(size_t a, size_t b, size_t nbSimus, 
//			   CheyetteDD_VanillaSwaptionApproxPricer_PTR pApproxPricer, std::ofstream& o) ;
//
//void helpPrinter(const std::string& descriptionData, const std::vector<double>& data, std::ofstream& o) ;
//
////1er element du vecteur : strikes
////2eme element du vecteur : prix
////3eme element du vecteur : vol Black
////pour un strike de 2%, rentrer 2
//std::vector<std::vector<double>> smileVect(size_t fileNumber, size_t a, size_t b,
//												  double strikeATM_Bloomberg, size_t nbSimus) ;
//
////resultats pour toutes les swaptions coterminales
////ecrit dans le fichier "Comparaison_Approx_MC.csv"
//
//void test_print_resultats(size_t coterminal, double k, double sigma, double m) ;
//
//
////tests convexite pour Phi et smile
//
//void evaluatePhi(double t_min, double t_max, double S_min, double S_max) ;
//
//void evaluatePhi_lineaire(double t_min, double t_max, double S_min, double S_max) ;
//
//
////prints smile for swaption aY bY in a file
////model calibrated on a given input file
//void smile(size_t a, size_t b, double strikeATM_Bloomberg) ;
//
//void print_smile(size_t a, size_t b, double strikeATM_Bloomberg, 
//				 double price, double annuity0, double swapRate0, std::ofstream& o) ;
//
//
//
//std::vector<double> priceSwaptionsMultipleStrikes(double strikeATM_Bloomberg, size_t nbSimus) ;
//
//
//
//void test_printMatrix_ATM(std::vector<size_t> expiry_maturity, Tenor floatingLegTenor, Tenor fixedLegTenor, 
//					 int curveChoice, int shiftChoice, 
//					 std::vector<double> x, std::vector<double> m_y, std::vector<double> sigma_y, double k, 
//					 std::vector<size_t> nbSimus) ;
//
//void matrixPrint(std::ofstream& o, std::vector<size_t> expiry_maturity, std::string chaine, const Matrice& M) ;
//
//void test_printElement_ATM(size_t index1, size_t index2, std::vector<size_t>& expiry_maturity, 
//						  Tenor floatingLegTenor, Tenor fixedLegTenor, 
//						 int curveChoice, int shiftChoice, 
//						 CheyetteDD_Model_PTR modele_test_PTR, 
//						 std::vector<size_t>& nbSimus, 
//						 Matrice& matriceApprox, Matrice& matriceS0, Matrice& matriceBbarre, Matrice& matriceVolBlack,
//						 Matrice& matriceMCprice, Matrice& matriceMC_ICinf, Matrice& matriceMC_ICsup, 
//						 Matrice& matriceRelativeError, Matrice& matriceAbsoluteError, Matrice& matriceKtilde) ;
//
//void test_printMatrix(std::vector<size_t> expiry_maturity, double strike, Tenor floatingLegTenor, Tenor fixedLegTenor, 
//					 int curveChoice, int shiftChoice, 
//					 std::vector<double> x, std::vector<double> m_y, std::vector<double> sigma_y, double k, 
//					 std::vector<size_t> nbSimus) ;
//
//void test_printElement(double strike, size_t index1, size_t index2, std::vector<size_t> expiry_maturity, 
//						  Tenor floatingLegTenor, Tenor fixedLegTenor, 
//						 int curveChoice, int shiftChoice, 
//						 CheyetteDD_Model_PTR modele_test_PTR, 
//						 std::vector<size_t> nbSimus, 
//						 Matrice& matriceApprox, Matrice& matriceS0, Matrice& matriceBbarre, Matrice& matriceVolBlack,
//						 Matrice& matriceMCprice, Matrice& matriceMC_ICinf, Matrice& matriceMC_ICsup, 
//						 Matrice& matriceRelativeError, Matrice& matriceAbsoluteError, Matrice& matriceKtilde) ;