#pragma once

#include <Cheyette/unit_test/TestMC_vs_approx.h>
#include <Cheyette/CheyetteModel/CheyetteDD_Model.h>
#include <Cheyette/Calibrator/CheyetteDD_LocalCalibrator.h>
#include <Cheyette/Calibrator/CheyetteDD_CostFunctionLevel.h>
#include <Cheyette/Calibrator/CheyetteDD_CostFunctionSkew.h>
#include <Cheyette/Calibrator/CheyetteDD_CalibrationConfig.h>

//lecture dans fichier VCUB
#include <LMM/helper/InputFileManager.h>
#include <LMM/LmmModel/LmmSwaptionMarketData.h>
#include <LMM/Test/Tests.h>

//void testCalib() ;

CoTerminalSwaptionVol_PTR upperMatrixToCoTerminal(LmmSwaptionMarketData_PTR pLmmSwaptionMarketData) ;

void test_CheyetteCalibrationMarketData() ;

void test_calib_Cheyette_OneFile(const std::string& mkt_data_file) ;

void marketData_Cheyette_LocalCalibration(	const CheyetteDD_CalibrationConfig& config,
											LmmSwaptionMarketData_PTR pLmmSwaptionMarketData) ; 

std::vector<double> arrayToVector(QuantLib::Array a) ;

