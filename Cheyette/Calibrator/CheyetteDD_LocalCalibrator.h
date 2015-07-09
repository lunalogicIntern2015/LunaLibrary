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
	Array	sigmaInitiate_ ;	//param initial pour sigma(t)				
	Array	mInitiate_;			//param initial pour m(t)	

	mutable Array	calibrated_sigma_ ;			//sigma(t) calibre				
	mutable Array	calibrated_m_;				//m(t) calibre

public :
	CheyetteDD_LocalCalibrator(	const QuantLib::Size & maxIterations,       
								const QuantLib::Real & rootEpsilon,        
								const QuantLib::Real & functionEpsilon,    
								CheyetteBaseCostFunction_PTR & cheyetteBaseCostFunction_PTR,
								Array sigmaInitiate, Array mInitiate,
								Array calibrated_sigma, Array calibrated_m)
		: CheyetteBaseCalibrator(maxIterations, rootEpsilon , functionEpsilon, cheyetteBaseCostFunction_PTR),
			sigmaInitiate_(sigmaInitiate), mInitiate_(mInitiate), 
			calibrated_sigma_(calibrated_sigma), calibrated_m_(calibrated_m)
	{}

	virtual ~CheyetteDD_LocalCalibrator()
	{
		//à compléter
	}

	//void minimizeNoConstraint(	QuantLib::Array xInitiate, 
	//							QuantLib::LevenbergMarquardt minimizationSolver, 
	//							CheyetteBaseCostFunction_PTR cheyetteBaseCostFunction_PTR) ;

	void minimizePositiveConstraint(	QuantLib::Array& calibratedArray, 
										QuantLib::LevenbergMarquardt& minimizationSolver, 
										CheyetteBaseCostFunction_PTR cheyetteBaseCostFunction_PTR) ;
	
	//calibration de sigma -> m -> sigma sur une swaption
	void calibrateOneSwaption(size_t indexSwaption) ;

	//boucle sur les differentes swaptions
	void calibrate() ;

	////! print plus plus without erease
	//void printPlusPlus(const std::string& base_filename) const ;
	//	
	//virtual void retrieve_calib_info();
};

typedef boost::shared_ptr<CheyetteDD_LocalCalibrator> CheyetteDD_LocalCalibrator_PTR;
typedef boost::shared_ptr<const CheyetteDD_LocalCalibrator> CheyetteDD_LocalCalibrator_CONSTPTR;

