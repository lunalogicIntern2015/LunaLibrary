#pragma once

#include <Cheyette/Calibration/CheyetteBaseCostFunction.h>
#include <Cheyette/Calibration/CoTerminalSwaptionVol.h>

/*For a setup of your own function, you need to implement a class 
which derives from the CostFunction class and implements the following virtual functions :

*/

class CheyetteDD_CostFunctionLevel : public CheyetteBaseCostFunction
{
private:
	CoTerminalSwaptionVol_CONSTPTR	coTerminalSwaptionVol_PTR_ ;
public:

	CheyetteDD_CostFunctionLevel(	ostream& o, CoTerminalSwaptionVol_CONSTPTR coTerminalSwaptionVol_PTR,   
									CheyetteDD_VanillaSwaptionApproxPricer_PTR cheyetteApprox_PTR,
									size_t indexSwaption)
		: CheyetteBaseCostFunction(o, cheyetteApprox_PTR, indexSwaption), coTerminalSwaptionVol_PTR_(coTerminalSwaptionVol_PTR)
	{}

	virtual ~CheyetteDD_CostFunctionLevel(){}

	CoTerminalSwaptionVol_CONSTPTR	getCoTerminalSwaptionVol_PTR() const{return coTerminalSwaptionVol_PTR_ ;}

	//m fixé, on fait varier sigma
	virtual QuantLib::Disposable<QuantLib::Array> values(const QuantLib::Array& param_sigma1D) const ;

};

typedef boost::shared_ptr<CheyetteDD_CostFunctionLevel> CheyetteDD_CostFunctionLevel_PTR ;
typedef boost::shared_ptr<const CheyetteDD_CostFunctionLevel> CheyetteDD_CostFunctionLevel_CONSTPTR;

