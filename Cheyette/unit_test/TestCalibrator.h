#pragma once

#include <Cheyette/unit_test/TestMC_vs_approx.h>
#include <Cheyette/CheyetteModel/CheyetteDD_Model.h>
#include <Cheyette/Calibrator/CheyetteDD_LocalCalibrator.h>
#include <Cheyette/Calibrator/CheyetteDD_CostFunctionLevel.h>
#include <Cheyette/Calibrator/CheyetteDD_CostFunctionSkew.h>
#include <Cheyette/Calibrator/MarketData.h>

//lecture dans fichier VCUB
#include <LMM/LmmModel/LmmSwaptionMarketData.h>
#include <LMM/Test/Tests.h>

void testCalib() ;