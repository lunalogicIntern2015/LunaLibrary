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
//#include <LMM/RNGenerator/McGenerator.h>
//#include <LMM/RNGenerator/RNGenerator.h>
//#include <LMM/helper/InputFileManager.h>
//#include <LMM/numeric/NumericalMethods.h>
//#include <LMM/LmmModel/Correlation.h>
//#include <LMM/LmmModel/McTerminalLmm.h>
//#include <LMM/LmmModel/Lmm.h>
//#include <LMM/LmmModel/ConstShifted_HGVolatilityFunction.h>
//#include <LMM/LmmModel/LmmSwaptionMarketData.h>
//#include <LMM/pricer/McLmmVanillaSwaptionPricer.h>
//
//#include <JBLMM/Element/Rate1.h>  
//#include <JBLMM/Element/ConstRate.h>  
//#include <JBLMM/Element/LiborRate.h>  
//#include <JBLMM/Element/VanillaSwapRate.h>  
//#include <JBLMM/Longstaff_Schwartz/EV_Basis/Basis.h>
//#include <JBLMM/Longstaff_Schwartz/EV_Basis/Basis_Evaluator.h>
//#include <JBLMM/Longstaff_Schwartz/Regression/Regression_LS.h>
//#include <JBLMM/Longstaff_Schwartz/Simulation/McLmm_LS.h>
//#include <JBLMM/Instrument/CallableInstrument.h>
//#include <JBLMM/Instrument/CallableSwap.h>
//#include <JBLMM/Instrument/InstrumentFactory.h>
//#include <JBLMM/Pricer/McLmmPricer.h>
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
