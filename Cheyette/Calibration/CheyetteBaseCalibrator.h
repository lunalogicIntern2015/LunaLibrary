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

#include <Cheyette/Calibration/CheyetteBaseCostFunction.h>
#include <Cheyette/Calibration/Cheyette_CostFunctionLevel.h>
#include <Cheyette/Calibration/Cheyette_CostFunctionSkew.h>

#include <ql/math/optimization/costfunction.hpp>


namespace HelperArray
{
	//std::vector<double> arrayToVector(QuantLib::Array a) ;
	QuantLib::Array vectorToArray(const std::vector<double>& v) ;
};

//classe mere. Pour le moment, seul LocalCalibrator en dérive
//possibilité de faire GlobalCalibrator
class CheyetteBaseCalibrator
{
protected :
	QuantLib::EndCriteria					stopCriteria_;			//construit avec maxIterations, rootEpsilon, functionEpsilon
	QuantLib::EndCriteria::Type				endConvergenceType_;
	boost::shared_ptr<QuantLib::Constraint> pConstraint_;

public :
	CheyetteBaseCalibrator( const QuantLib::Size& maxIterations,
							const QuantLib::Real& rootEpsilon,        
							const QuantLib::Real& functionEpsilon) ;    
							
	virtual ~CheyetteBaseCalibrator(){}

	void minimizeNoConstraint(	QuantLib::Array& xInitiate, 
								QuantLib::LevenbergMarquardt& minimizationSolver, 
								CheyetteBaseCostFunction_PTR pCostFunction) const ;

	void minimizePositiveConstraint(QuantLib::Array& calibratedArray1D, 	
									QuantLib::LevenbergMarquardt& minimizationSolver, 
									CheyetteBaseCostFunction_PTR pCostFunction) const ;

	void minimizeBoundaryConstraint(QuantLib::Array& calibratedArray1D, 
									QuantLib::LevenbergMarquardt& minimizationSolver,
									CheyetteBaseCostFunction_PTR pCostFunction) const ;

	void printMinimizationInfo(QuantLib::Problem optimizationProblem, 
							   QuantLib::EndCriteria::Type endConvergenceType) const ;

	virtual void solve() const = 0 ;		//minimise pour une swaption
	virtual void calibrate() const = 0 ;	//calibration du modèle sur le vecteur des swaptions coterminales
};

typedef boost::shared_ptr<CheyetteBaseCalibrator> CheyetteBaseCalibrator_PTR;
typedef boost::shared_ptr<const CheyetteBaseCalibrator> CheyetteBaseCalibrator_CONSTPTR;