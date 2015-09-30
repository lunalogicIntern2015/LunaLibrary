#include "CheyetteBaseCalibrator.h"

#include <fstream>

namespace HelperArray
{
QuantLib::Array vectorToArray(const std::vector<double>& v)
{
	QuantLib::Array a(v.size()) ;
	for (size_t i = 0 ; i < v.size() ; ++i)
	{
		a[i] = v[i] ;
	}
	return a ;
}

//std::vector<double> arrayToVector(QuantLib::Array a)
//{
//	std::vector<double> vect ;
//	vect.reserve(a.size()) ;
//	for (size_t i = 0 ; i < a.size() ; ++i)
//	{
//		vect[i] = a[i] ;
//	}
//	return vect ;
//}

};

CheyetteBaseCalibrator::CheyetteBaseCalibrator( const QuantLib::Size& maxIterations,
												const QuantLib::Real& rootEpsilon,        
												const QuantLib::Real& functionEpsilon)    
	: pConstraint_(new QuantLib::NoConstraint() ) // default constraint is a no constraint
	, stopCriteria_(maxIterations, 100 , rootEpsilon, functionEpsilon, 0.)
{}



void CheyetteBaseCalibrator::minimizeNoConstraint(	QuantLib::Array& xInitiate, 
													QuantLib::LevenbergMarquardt& minimizationSolver, 
													CheyetteBaseCostFunction_PTR pCostFunction) const
{
	QuantLib::Array					start_point = xInitiate ;
	QuantLib::NoConstraint			noConstraint ;
	QuantLib::Problem				optimizationProblem(*pCostFunction, noConstraint, start_point); 
	
	QuantLib::EndCriteria::Type		endConvergenceType = minimizationSolver.minimize(optimizationProblem, stopCriteria_);	

	std :: cout << " Criteria :" << endConvergenceType << std :: endl ;
	std :: cout << " Root :" << optimizationProblem.currentValue() << std :: endl ;
	std :: cout << " Min Function Value :" << optimizationProblem.functionValue () << std :: endl ;
}

//calibratedArray : point de depart de l'algo d'optimisation
//calibratedArray est ensuite mis à jour avec la valeur calibrée
void CheyetteBaseCalibrator::minimizePositiveConstraint(QuantLib::Array& calibratedArray1D, 
														QuantLib::LevenbergMarquardt& minimizationSolver, 
														CheyetteBaseCostFunction_PTR pCostFunction) const
{
	QuantLib::PositiveConstraint	positiveConstraint ;
	QuantLib::Problem				optimizationProblem(*pCostFunction, positiveConstraint, calibratedArray1D);  
	QuantLib::EndCriteria::Type		endConvergenceType = minimizationSolver.minimize(optimizationProblem, stopCriteria_);	

	calibratedArray1D = optimizationProblem.currentValue() ;
//	calibrated_a_[currentSwaptionIndex_]
	printMinimizationInfo(optimizationProblem, endConvergenceType) ;
}

void CheyetteBaseCalibrator::minimizeBoundaryConstraint(QuantLib::Array& calibratedArray1D, 
															QuantLib::LevenbergMarquardt& minimizationSolver, 
															CheyetteBaseCostFunction_PTR pCostFunction) const
{
	QuantLib::BoundaryConstraint bc(0.000000001, 1.) ;
	QuantLib::Problem				optimizationProblem(*pCostFunction, bc, calibratedArray1D);  //pt de depart de l'algo
	QuantLib::EndCriteria::Type		endConvergenceType = minimizationSolver.minimize(optimizationProblem, stopCriteria_);	
	
	calibratedArray1D = optimizationProblem.currentValue() ;

	printMinimizationInfo(optimizationProblem, endConvergenceType) ;
}

void CheyetteBaseCalibrator::printMinimizationInfo(QuantLib::Problem optimizationProblem, 
													   QuantLib::EndCriteria::Type endConvergenceType) const
{
	std::cout << " Criteria :" << endConvergenceType << std::endl ;
	std::cout << " Root :" << optimizationProblem.currentValue() << std::endl ;
	std::cout << " Min Function Value :" << optimizationProblem.functionValue () << std::endl ;
	std::cout << " --------------------------------------------------------------" << std::endl ;
	std::cout << std::endl ;

	//o_ << " Criteria : ;" << endConvergenceType << std::endl ;
	//o_ << " Root : ;" << optimizationProblem.currentValue() << std::endl ;
	//o_ << " Min Function Value : ;" << optimizationProblem.functionValue () << std::endl ;
	//o_ << std::endl ;
}