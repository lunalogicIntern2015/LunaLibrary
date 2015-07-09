#pragma once

#include <Cheyette/Calibrator/CheyetteBaseCostFunction.h>
#include <Cheyette/Calibrator/CoTerminalSwaptionVol.h>

/*For a setup of your own function, you need to implement a class 
which derives from the CostFunction class and implements the following virtual functions :

*/

class CheyetteDD_CostFunctionLevel : public CheyetteBaseCostFunction
{
private:
	CoTerminalSwaptionVol_CONSTPTR	coTerminalSwaptionVol_PTR_ ;
public:

	CheyetteDD_CostFunctionLevel(	CoTerminalSwaptionVol_CONSTPTR coTerminalSwaptionVol_PTR, 
									size_t indexSwaption,  
									CheyetteDD_VanillaSwaptionApproxPricer_PTR cheyetteApprox_PTR)
		: CheyetteBaseCostFunction(indexSwaption, cheyetteApprox_PTR), coTerminalSwaptionVol_PTR_(coTerminalSwaptionVol_PTR)
	{}

	virtual ~CheyetteDD_CostFunctionLevel()
	{ 
		//� compl�ter ?
	}

	//m fix�, on fait varier sigma
	virtual Disposable<Array> values(const Array& param_sigma) const ;

};

typedef boost::shared_ptr<CheyetteDD_CostFunctionLevel> CheyetteDD_CostFunctionLevel_PTR;
typedef boost::shared_ptr<const CheyetteDD_CostFunctionLevel> CheyetteDD_CostFunctionLevel_CONSTPTR;

