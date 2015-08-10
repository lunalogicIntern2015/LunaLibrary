#pragma once

#include <Cheyette/unit_test/TestMC.h>
#include <Cheyette/Model/CheyetteDD_Model.h>
#include <Cheyette/Calibration/CheyetteDD_LocalCalibrator.h>
#include <Cheyette/Calibration/CheyetteDD_CostFunctionLevel.h>
#include <Cheyette/Calibration/CheyetteDD_CostFunctionSkew.h>
#include <Cheyette/Calibration/MarketDataConvertor.h>

//lecture dans fichier VCUB
#include <LMM/helper/InputFileManager.h>
#include <LMM/LmmSwaptionMarketData.h>
#include <LMM/LmmSwaptionMarketDataFull.h>
#include <LMM/helper/LMMTenorStructure.h>

#include <fstream>

#include <LMM/Test/Tests.h>


void printAllResults_calibratedData(size_t fileNumber, size_t coterminal) ;

CheyetteDD_Model_PTR getCalibratedModel(size_t fileNumber, size_t coterminal) ;

void lancementAuto() ; 

//std::vector<double> arrayToVector(QuantLib::Array a) ;

QuantLib::Array vectorToArray(std::vector<double> v) ;

