#pragma once

#include <Cheyette/Calibration/CheyetteBaseCostFunction.h>
#include <Cheyette/Calibration/CheyetteMarketData_2.h>

/*For a setup of your own function, you need to implement a class 
which derives from the CostFunction class and implements the following virtual functions :

*/

class Cheyette_CostFunctionConvexity : public CheyetteBaseCostFunction
{
private:
	size_t indexSwaption_ ;		
	CheyetteMarketData_2_PTR cheyetteMarketData_2_PTR_ ;

	std::vector<double> convexity_ ;
public:

	Cheyette_CostFunctionConvexity(	Cheyette_SwaptionPricer_Approx_PTR cheyetteApprox_PTR,
									CheyetteMarketData_PTR cheyetteMarketData_PTR,
									size_t indexSwaption,
									CheyetteMarketData_2_PTR cheyetteMarketData_2_PTR
									)
		:	CheyetteBaseCostFunction(cheyetteApprox_PTR, cheyetteMarketData_PTR), indexSwaption_(indexSwaption),
			cheyetteMarketData_2_PTR_(cheyetteMarketData_2_PTR)
	{
		calculateConvexity() ; //initialise le vecteur convexity_
	}

	virtual ~Cheyette_CostFunctionConvexity(){}
	void calculateConvexity() ;

	void setIndexSwaption(size_t indexSwaption){indexSwaption_ = indexSwaption ;}

	virtual QuantLib::Disposable<QuantLib::Array> values(const QuantLib::Array& param_c_1D) const ;
	double volShift(double strike, double shift) const ;

};

typedef boost::shared_ptr<Cheyette_CostFunctionConvexity>		Cheyette_CostFunctionConvexity_PTR ;
typedef boost::shared_ptr<const Cheyette_CostFunctionConvexity> Cheyette_CostFunctionConvexity_CONSTPTR;

