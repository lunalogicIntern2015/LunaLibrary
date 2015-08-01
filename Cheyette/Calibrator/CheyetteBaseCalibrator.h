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

#include <Cheyette/Calibrator/CheyetteBaseCostFunction.h>
#include <Cheyette/Calibrator/CheyetteDD_CostFunctionLevel.h>
#include <Cheyette/Calibrator/CheyetteDD_CostFunctionSkew.h>

#include <ql/math/optimization/costfunction.hpp>

using namespace QuantLib ;  //pour CostFunction, Real, Array

//classe mere. Pour le moment, seul LocalCalibrator en d�rive
//possibilit� de faire GlobalCalibrator
class CheyetteBaseCalibrator
{
protected :
	QuantLib::EndCriteria			stopCriteria_; 
	//CheyetteBaseCostFunction_PTR	cheyetteBaseCostFunction_PTR_ ;

	bool isVirtualCalibration_;

	std::string base_general_result_info_;
	
	std::ostream& o_ ;

	//! minimizator
	QuantLib::EndCriteria::Type endConvergenceType_;

	bool use_positive_constraint_;
	boost::shared_ptr<Constraint> pConstraint_;


	//! storage information after minimization
	size_t total_number_called_;
	double total_minimization_time_;
	
	// stored relative error for post calibration
	double max_quote_rel_error_;
	std::pair<size_t,size_t> max_quote_rel_error_position_;
	
	//// stored error for post calibration
	double max_quote_abs_error_;
	std::pair<size_t,size_t> max_quote_abs_error_position_;

	double quote_error_l2_, quote_error_l1_, quote_error_lInf_;

	// post calibration calculated errors
	void retrieve_calib_global_error() ;
	
	//virtual void retrieve_calib_info() = 0 ;

	void printAnnexeStopCriteriaLevenbergMarquardt( std::ofstream & stream ) const ;

public :
	CheyetteBaseCalibrator( std::ostream& o, 
							const QuantLib::Size& maxIterations,
							const QuantLib::Real& rootEpsilon,        
							const QuantLib::Real& functionEpsilon) ;    
						//	CheyetteBaseCostFunction_PTR cheyetteBaseCostFunction) ;

	bool isVirtualCalibration() const { return isVirtualCalibration_; }

	virtual void solve() = 0 ;

	//virtual void printPlusPlus(const std::string& base_filename) const = 0 ;

	//activate the positive constraint to be sure that the cellulle calibration do not have negative solution
	void activate_PositiveConstraint() ;

	// methods retrieving information from QuantLib::EndCriteria
	QuantLib::Size maxIterations()                const {return stopCriteria_.maxIterations() ;}
	QuantLib::Size maxStationaryStateIterations() const {return stopCriteria_.maxStationaryStateIterations() ;}
	QuantLib::Real rootEpsilon()                  const {return stopCriteria_.rootEpsilon() ;}
	QuantLib::Real functionEpsilon()              const {return stopCriteria_.functionEpsilon() ;}
	QuantLib::Real gradientNormEpsilon()          const {return stopCriteria_.gradientNormEpsilon() ;}
	const QuantLib::EndCriteria::Type& get_EndConvergenceType() const { return endConvergenceType_ ;}
		
	const double& get_QuoteError_L2() const { return quote_error_l2_; }
	const double& get_QuoteError_L1() const { return quote_error_l1_; }
	const double& get_QuoteError_LInf() const { return quote_error_lInf_; }


	const std::string& get_BaseGeneral_Result_Info() const { return base_general_result_info_;}
};

typedef boost::shared_ptr<CheyetteBaseCalibrator> CheyetteBaseCalibrator_PTR;
typedef boost::shared_ptr<const CheyetteBaseCalibrator> CheyetteBaseCalibrator_CONSTPTR;