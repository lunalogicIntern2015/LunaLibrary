#pragma once

#include <Cheyette/unit_test/TestMC.h>
#include <Cheyette/CheyetteModel/CheyetteDD_Model.h>
#include <Cheyette/Calibrator/CheyetteDD_LocalCalibrator.h>
#include <Cheyette/Calibrator/CheyetteDD_CostFunctionLevel.h>
#include <Cheyette/Calibrator/CheyetteDD_CostFunctionSkew.h>

//lecture dans fichier VCUB
#include <LMM/helper/InputFileManager.h>
#include <LMM/LmmModel/LmmSwaptionMarketData.h>
#include <LMM/LmmModel/LmmSwaptionMarketDataFull.h>
#include <LMM/helper/LMMTenorStructure.h>

#include <fstream>

#include <LMM/Test/Tests.h>
#include <Cheyette/Calibrator/MarketDataConvertor.h>


void testCalib(size_t fileNumber, size_t coterminal) ;

void generateSmile(size_t model_nbYear, size_t fileNumber, size_t coterminal, std::ofstream& o) ;


//std::vector<double> arrayToVector(QuantLib::Array a) ;

QuantLib::Array vectorToArray(std::vector<double> v) ;