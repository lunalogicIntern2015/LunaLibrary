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
#include <LMM/helper/LMMTenorStructure.h>

#include <fstream>

#include <LMM/Test/Tests.h>

std::vector<size_t> getVectorExpiry(size_t coterminal) ;
std::vector<size_t> getVectorTenor(size_t coterminal) ;
std::vector<double> getVectorStrike(LmmSwaptionMarketData_PTR pLmmSwaptionMarketData, size_t coterminal) ;
std::vector<double> getVectorVol(LmmSwaptionMarketData_PTR pLmmSwaptionMarketData, size_t coterminal) ;
std::vector<double> getVectorSkew(LmmSwaptionMarketData_PTR pLmmSwaptionMarketData, size_t coterminal) ;
std::vector<VanillaSwaption_PTR> getVectorSwaptions(const Tenor& tenorFloat, const Tenor& tenorFixed,
													const std::vector<size_t>& vectExpiry, 
													const std::vector<size_t>& vectTenor, 
													const std::vector<double>& vectStrike) ;

void testCalib2(size_t fileNumber, size_t coterminal) ;
//
//void testCalib() ;
//
//CoTerminalSwaptionVol_PTR upperMatrixToCoTerminal(LmmSwaptionMarketData_PTR pLmmSwaptionMarketData) ;
//
//void test_CheyetteCalibrationMarketData() ;
//
//void test_calib_Cheyette_OneFile(const std::string& mkt_data_file) ;
//
//void marketData_Cheyette_LocalCalibration(	const CheyetteDD_CalibrationConfig& config,
//											LmmSwaptionMarketData_PTR pLmmSwaptionMarketData) ; 
//
//std::vector<double> arrayToVector(QuantLib::Array a) ;

QuantLib::Array vectorToArray(std::vector<double> v) ;

//CoTerminalSwaptionVol upperMatrixToCoTerminalVol(UpperTriangleVanillaSwaptionQuotes_PTR upperQuotes,
//												  size_t coterminal) ;
//
//CoTerminalSwaptionSkew upperMatrixToCoTerminalSkew(UpperTriangleVanillaSwaptionQuotes_PTR upperQuotes,
//												  size_t coterminal, double shift) ;
//
//CourbeInput upperMatrixZC_To_CourbeInput_yield(LiborQuotes_PTR liborQuotes_PTR) ;