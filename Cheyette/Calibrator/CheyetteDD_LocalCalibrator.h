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
	Array	sigmaInitiate_1D_ ;	//param initial pour sigma(t)				
	Array	mInitiate_1D_ ;			//param initial pour m(t)	

	mutable Array	calibrated_sigma_1D_ ;			//sigma(t) calibre				
	mutable Array	calibrated_m_1D_ ;				//m(t) calibre

public :
	CheyetteDD_LocalCalibrator(	std::ostream& o,
								const QuantLib::Size & maxIterations,       
								const QuantLib::Real & rootEpsilon,        
								const QuantLib::Real & functionEpsilon,    
								Array sigmaInitiate_1D, Array mInitiate_1D,
								Array calibrated_sigma_1D, Array calibrated_m_1D,
								CheyetteDD_CostFunctionLevel_PTR		costFunctionLevel_PTR,
								CheyetteDD_CostFunctionSkew_PTR		costFunctionSkew_PTR)
		: CheyetteBaseCalibrator(o, maxIterations, rootEpsilon , functionEpsilon),
			sigmaInitiate_1D_(sigmaInitiate_1D), mInitiate_1D_(mInitiate_1D), 
			calibrated_sigma_1D_(calibrated_sigma_1D), calibrated_m_1D_(calibrated_m_1D),
			costFunctionLevel_PTR_(costFunctionLevel_PTR), costFunctionSkew_PTR_(costFunctionSkew_PTR)
	{
		assert(sigmaInitiate_1D_.size() ==1) ;
		assert(mInitiate_1D_.size() ==1) ;
	}

	virtual ~CheyetteDD_LocalCalibrator(){}

	//void minimizeNoConstraint(	QuantLib::Array xInitiate, 
	//							QuantLib::LevenbergMarquardt minimizationSolver, 
	//							CheyetteBaseCostFunction_PTR cheyetteBaseCostFunction_PTR) ;

	void minimizePositiveConstraint(QuantLib::Array& calibratedArray1D, 	
									QuantLib::LevenbergMarquardt& minimizationSolver, 
									CheyetteDD_CostFunctionLevel pCostFunctionLevel) ;

	void minimizeBoundaryConstraint(QuantLib::Array& calibratedArray1D, 
									QuantLib::LevenbergMarquardt& minimizationSolver,
									CheyetteDD_CostFunctionSkew pCostFunctionSkew) ;

	//calibration de sigma -> m -> sigma sur une swaption
	void solve() ;

	//void retrieve_calib_info();
	
	//print plus plus without erease
	//void printPlusPlus(const std::string& base_filename) const ;	

};

typedef boost::shared_ptr<CheyetteDD_LocalCalibrator> CheyetteDD_LocalCalibrator_PTR;
typedef boost::shared_ptr<const CheyetteDD_LocalCalibrator> CheyetteDD_LocalCalibrator_CONSTPTR;

