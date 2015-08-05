#pragma once

#include <Cheyette/unit_test/TestMC.h>
#include <Cheyette/Model/CheyetteDD_Model.h>
#include <Cheyette/Calibration/CheyetteDD_LocalCalibrator.h>
#include <Cheyette/Calibration/CheyetteDD_CostFunctionLevel.h>
#include <Cheyette/Calibration/CheyetteDD_CostFunctionSkew.h>

//lecture dans fichier VCUB
#include <LMM/Helper/InputFileManager.h>
#include <LMM/LmmSwaptionMarketData.h>
#include <LMM/LmmSwaptionMarketDataFull.h>
#include <LMM/Helper/LMMTenorStructure.h>

#include <fstream>

#include <LMM/Test/Tests.h>
#include <Cheyette/Calibration/MarketDataConvertor.h>


void testCalib(size_t fileNumber, size_t coterminal) ;

void generateSmile(size_t model_nbYear, size_t fileNumber, size_t coterminal, std::ofstream& o) ;


//std::vector<double> arrayToVector(QuantLib::Array a) ;

QuantLib::Array vectorToArray(std::vector<double> v) ;
