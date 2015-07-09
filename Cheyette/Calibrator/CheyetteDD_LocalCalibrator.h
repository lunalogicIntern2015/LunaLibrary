#pragma once

#include <boost/shared_ptr.hpp>

#include <ql/math/optimization/costfunction.hpp>
#include <ql/math/optimization/endcriteria.hpp>
#include <ql/math/optimization/constraint.hpp>
#include <ql/math/optimization/problem.hpp>
#include <ql/math/optimization/simplex.hpp>
#include <ql/math/optimization/bfgs.hpp> 
#include <ql/math/optimization/conjugategradient.hpp>
#include <ql/math/optimization/levenbergmarquardt.hpp>

#include <Cheyette/Calibrator/CheyetteBaseCalibrator.h>
#include <Cheyette/Fonction.h>


class CheyetteDD_LocalCalibrator : public CheyetteBaseCalibrator
{
private:
	//fonctions de coût
	CheyetteDD_CostFunctionLevel_PTR		costFunctionLevel_PTR_ ;
	CheyetteDD_CostFunctionSkew_PTR			costFunctionSkew_PTR_ ;
	//points de depart
	Array	sigmaInitiate_ ;	//param initial pour sigma(t)				
	Array	mInitiate_;			//param initial pour m(t)	

	mutable Array	calibrated_sigma_ ;			//sigma(t) calibre				
	mutable Array	calibrated_m_;				//m(t) calibre

public :
	CheyetteDD_LocalCalibrator(	const QuantLib::Size & maxIterations,       
								const QuantLib::Real & rootEpsilon,        
								const QuantLib::Real & functionEpsilon,    
								Array sigmaInitiate, Array mInitiate,
								Array calibrated_sigma, Array calibrated_m,
								CheyetteDD_CostFunctionLevel_PTR		costFunctionLevel_PTR,
								CheyetteDD_CostFunctionSkew_PTR		costFunctionSkew_PTR)
		: CheyetteBaseCalibrator(maxIterations, rootEpsilon , functionEpsilon),
			sigmaInitiate_(sigmaInitiate), mInitiate_(mInitiate), 
			calibrated_sigma_(calibrated_sigma), calibrated_m_(calibrated_m),
			costFunctionLevel_PTR_(costFunctionLevel_PTR), costFunctionSkew_PTR_(costFunctionSkew_PTR)
	{}

	virtual ~CheyetteDD_LocalCalibrator(){}

	//void minimizeNoConstraint(	QuantLib::Array xInitiate, 
	//							QuantLib::LevenbergMarquardt minimizationSolver, 
	//							CheyetteBaseCostFunction_PTR cheyetteBaseCostFunction_PTR) ;

	void minimizePositiveConstraint(	QuantLib::Array& calibratedArray, 
										QuantLib::LevenbergMarquardt& minimizationSolver, 
										CheyetteBaseCostFunction_PTR cheyetteBaseCostFunction_PTR) ;
	
	virtual std::string Name() const 
	{
		if(isVirtualCalibration_) return "VIRTUAL LocalCalibrator";
		else  return "Local Calibrator";
	}

	//calibration de sigma -> m -> sigma sur une swaption
	void calibrateOneSwaption(size_t indexSwaption) ;

	//boucle sur les differentes swaptions
	void solve() ;

	//void retrieve_calib_info();
	
	//print plus plus without erease
	//void printPlusPlus(const std::string& base_filename) const ;	

};

typedef boost::shared_ptr<CheyetteDD_LocalCalibrator> CheyetteDD_LocalCalibrator_PTR;
typedef boost::shared_ptr<const CheyetteDD_LocalCalibrator> CheyetteDD_LocalCalibrator_CONSTPTR;

