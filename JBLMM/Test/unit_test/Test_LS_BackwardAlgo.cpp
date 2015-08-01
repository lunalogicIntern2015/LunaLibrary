#include <JBLMM/Test/JBTests.h>
#include <LMM/Test/Tests.h>

#include <iostream>
#include <fstream> 
#include <string> 
#include <iterator>
#include <algorithm>


#include <LMM/RNGenerator/McGenerator.h>
#include <LMM/RNGenerator/RNGenerator.h>
#include <LMM/helper/InputFileManager.h>
#include <LMM/numeric/NumericalMethods.h>

#include <LMM/LmmModel/Correlation.h>
#include <LMM/LmmModel/McTerminalLmm.h>

#include <LMM/LmmModel/Lmm.h>
#include <LMM/LmmModel/ConstShifted_HGVolatilityFunction.h>
#include <LMM/LmmModel/LmmSwaptionMarketData.h>

#include <JBLMM/Element/Rate1.h>  
#include <JBLMM/Element/ConstRate.h>  
#include <JBLMM/Element/LiborRate.h>  
#include <JBLMM/Element/VanillaSwapRate.h>  
#include <JBLMM/Longstaff_Schwartz/EV_Basis/Basis.h>
#include <JBLMM/Longstaff_Schwartz/EV_Basis/Basis_Evaluator.h>
#include <JBLMM/Longstaff_Schwartz/Regression/Regression_LS.h>
#include <JBLMM/Longstaff_Schwartz/Simulation/McLmm_LS.h>
#include <JBLMM/Instrument/CallableInstrument.h>
#include <JBLMM/Instrument/CallableSwap.h>
#include <JBLMM/Instrument/InstrumentFactory.h>
#include <JBLMM/Pricer/McLmmPricer.h>
#include <JBLMM/Pricer/McLmmGenericSwapPricer.h>
#include <JBLMM/Longstaff_Schwartz/LS_BackwardAlgo.h>
#include <JBLMM/Longstaff_Schwartz/LS_ForwardAlgo.h>

//void LS_BackwardAlgo()
//{
//	
//}

