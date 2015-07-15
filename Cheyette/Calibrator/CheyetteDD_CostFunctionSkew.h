#pragma once

#include <Cheyette/Calibrator/CheyetteBaseCostFunction.h>
#include <Cheyette/Calibrator/CoTerminalSwaptionSkew.h>

/*For a setup of your own function, you need to implement a class 
which derives from the CostFunction class and implements the following virtual functions :

*/

class CheyetteDD_CostFunctionSkew : public CheyetteBaseCostFunction
{
private:
	CoTerminalSwaptionSkew_CONSTPTR coTerminalSwaptionSkew_PTR_ ;
public:
	CheyetteDD_CostFunctionSkew(	ostream& o, CoTerminalSwaptionSkew_CONSTPTR coTerminalSwaptionSkew_PTR, 
									CheyetteDD_VanillaSwaptionApproxPricer_PTR cheyetteApprox_PTR,
									size_t indexSwaption)
		: CheyetteBaseCostFunction(o, cheyetteApprox_PTR, indexSwaption), coTerminalSwaptionSkew_PTR_(coTerminalSwaptionSkew_PTR)
	{}

	virtual ~CheyetteDD_CostFunctionSkew(){}

	CoTerminalSwaptionSkew_CONSTPTR getCoTerminalSwaptionSkew_PTR() const {return coTerminalSwaptionSkew_PTR_ ;}

	//sigma fixé, on fait varier m pour prendre en compte le skew
	virtual Disposable<Array> values(const Array& param_m1D) const ;

	double volShift(double strike, double shift) const ;  //plus d'index swaption car il n'y a plus qu 'une swaption
};

typedef boost::shared_ptr<CheyetteDD_CostFunctionSkew> CheyetteDD_CostFunctionSkew_PTR;
typedef boost::shared_ptr<const CheyetteDD_CostFunctionSkew> CheyetteDD_CostFunctionSkew_CONSTPTR;

