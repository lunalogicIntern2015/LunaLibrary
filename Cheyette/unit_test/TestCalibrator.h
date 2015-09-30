#pragma once

#include <Cheyette/unit_test/TestMC.h>

#include <Cheyette/Model/CheyetteModel.h>

#include <Cheyette/Model/CheyetteDD_Model.h>
#include <Cheyette/Calibration/CheyetteDD_LocalCalibrator.h>
#include <Cheyette/Calibration/CheyetteQuad_LocalCalibrator.h>

#include <Cheyette/Calibration/Cheyette_CostFunctionLevel.h>
#include <Cheyette/Calibration/Cheyette_CostFunctionSkew.h>
#include <Cheyette/Calibration/MarketDataConvertor.h>
#include <Cheyette/Calibration/CheyetteMarketData.h>
#include <Cheyette/Calibration/CheyetteMarketData_2.h>

//lecture dans fichier VCUB
#include <LMM/helper/InputFileManager.h>
#include <LMM/LmmSwaptionMarketData.h>
//#include <LMM/LmmSwaptionMarketDataFull.h>
#include <LMM/helper/LMMTenorStructure.h>

#include <fstream>

#include <LMM/Test/Tests.h>

void helpPrinter(std::string description, std::vector<double> data, std::ofstream& o) ;

CheyetteMarketData_PTR recoverMarketData(const size_t fileNumber, const size_t coterminal) ;

CheyetteDD_LocalCalibrator_PTR creeCheyetteDD_Calibrator(size_t fileNumber, size_t coterminal) ;

void lancementCalibOneFile() ;

//pour Quad
CheyetteQuad_LocalCalibrator_PTR creeCheyetteQuad_Calibrator(size_t fileNumberMktData, 
															 size_t fileNumberMktData_2, 
															 size_t coterminal) ;

//ecriture dans fichier 
//calibration et smile
CheyetteModel_PTR getCalibratedModel(size_t fileNumber, size_t coterminal) ;

void printAllResults_calibratedData(size_t fileNumber, size_t coterminal) ;

void printSwaptionMultipleStrikesMC(size_t coterminal, Tenor floatTenor, 
									CheyetteModel_PTR calibratedModel,
									VanillaSwaption_PTR pSwaption, 
									std::vector<size_t> nbSimus, 
									std::vector<double> shifts_bp, double annuity0, double swapRate0,
									std::ofstream& o) ;
void lancementAuto() ; 


