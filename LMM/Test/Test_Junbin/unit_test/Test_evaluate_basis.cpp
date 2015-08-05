//#include <JBLMM/Test/JBTests.h>
//#include <LMM/Test/Tests.h>
//
//#include <iostream>
//#include <fstream> 
//#include <string> 
//#include <iterator>
//#include <algorithm>
//
//#include <boost/math/distributions/inverse_gaussian.hpp>
//#include <boost/math/distributions/normal.hpp>
//
//#include <RNGenerator/McGenerator.h>
//#include <RNGenerator/RNGenerator.h>
//#include <LMM/Helper/InputFileManager.h>
//#include <Numeric/NumericalMethods.h>
//#include <LMM/Model/Correlation.h>
//#include <LMM/Mc/McTerminalLmm.h>
//#include <LMM/Model/Lmm.h>
//#include <LMM/Model/ConstShifted_HGVolatilityFunction.h>
//#include <LMM/LmmSwaptionMarketData.h>
//#include <LMM/Pricer/McLmmPricer/McLmmVanillaSwaptionPricer.h>
//
//#include <Instrument/Rate/Rate1.h>  
//#include <Instrument/Rate/ConstRate.h>  
//#include <Instrument/Rate/LiborRate.h>  
//#include <Instrument/Rate/VanillaSwapRate.h>  
//#include <LMM/Pricer/Longstaff_Schwartz/Basis.h>
//#include <LMM/Pricer/Longstaff_Schwartz/Basis_Evaluator.h>
//#include <LMM/Pricer/Longstaff_Schwartz/Regression_LS.h>
//#include <LMM/Pricer/Longstaff_Schwartz/McLmm_LS.h>
//#include <Instrument/CallableOption/CallableInstrument.h>
//#include <JBInstrument/CallableSwap.h>
//#include <JBInstrument/InstrumentFactory.h>
//#include <LMM/Pricer/McLmmPricer/McLmmPricer.h>
//#include <JBLMM/Pricer/McLmmGenericSwapPricer.h>

//void Test_evaluate_basis_time()
//{
//	std::vector<std::vector<size_t>> subset;
//	for(size_t i = 0; i <= 2; i++) 
//	{
//		subset.push_back(std::vector<size_t>());
//		subset.back().push_back(0);
//		subset.back().push_back(0);
//		subset.back().push_back(i);							
//	}
//	for(size_t j = 1; j <= 2; j++) 
//	{
//		subset.push_back(std::vector<size_t>());
//		subset.back().push_back(0);
//		subset.back().push_back(j);
//		subset.back().push_back(0);
//	}	
//	for(size_t k = 1; k <= 2; k++) 
//	{
//		subset.push_back(std::vector<size_t>());
//		subset.back().push_back(k);
//		subset.back().push_back(0);
//		subset.back().push_back(0);
//	}	
//
//
//}
