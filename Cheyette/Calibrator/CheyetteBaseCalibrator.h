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

class CheyetteBaseCalibrator
{
protected :
	QuantLib::EndCriteria			stopCriteria_; 
	CheyetteBaseCostFunction_PTR	cheyetteBaseCostFunction_ ;

	//QuantLib::EndCriteria::Type endConvergenceType_;
	//bool use_positive_constraint_;
	////! storage information after minimization
	//size_t total_number_called_;
	//double total_minimization_time_;
	//
	//// stored relative error for post calibration
	//double max_quote_rel_error_;
	//std::pair<size_t,size_t> max_quote_rel_error_position_;
	//UpperTriangularDoubleMatrix rel_quote_error_matrix_;// error in (%)
	//
	////// stored error for post calibration
	//double max_quote_abs_error_;
	//std::pair<size_t,size_t> max_quote_abs_error_position_;
	//UpperTriangularDoubleMatrix abs_quote_error_matrix_;// absolute error 
	//double quote_error_l2_, quote_error_l1_, quote_error_lInf_;

public :
	CheyetteBaseCalibrator( const QuantLib::Size& maxIterations,
							const QuantLib::Real& rootEpsilon,        
							const QuantLib::Real& functionEpsilon,    
							CheyetteBaseCostFunction_PTR cheyetteBaseCostFunction) ;

	//virtual void solve() = 0 ;

	//getters
	QuantLib::Size					getMaxIterations()			const {return stopCriteria_.maxIterations() ;}
	QuantLib::Real					getRootEpsilon()			const {return stopCriteria_.rootEpsilon() ;}
	QuantLib::Real					getFunctionEpsilon()		const {return stopCriteria_.functionEpsilon() ;}
	CheyetteBaseCostFunction_PTR	getCheyetteCostFunction()	const {return cheyetteBaseCostFunction_ ;}

	//	QuantLib::Size getMaxStationaryStateIterations() const {return stopCriteria_.maxStationaryStateIterations() ;}
	//	QuantLib::Real getGradientNormEpsilon()          const {return stopCriteria_.gradientNormEpsilon() ;}
	//	const QuantLib::EndCriteria::Type & get_EndConvergenceType() const { return endConvergenceType_ ;}
		
	//const double& get_QuoteError_L2() const { return quote_error_l2_; }
	//const double& get_QuoteError_L1() const { return quote_error_l1_; }
	//const double& get_QuoteError_LInf() const { return quote_error_lInf_; }
	//const std::string& get_BaseGeneral_Result_Info() const { return base_general_result_info_;}
	//// post calibration calculated errors
	//void retreive_calib_global_error() ;
	//virtual void retreive_calib_info() = 0 ;
	//virtual void printPlusPlus(const std::string& base_filename) const = 0 ;
	//void printAnnexeStopCriteriaLevenbergMarquardt( std::ofstream & stream ) const ;
};

typedef boost::shared_ptr<CheyetteBaseCalibrator> CheyetteBaseCalibrator_PTR;
typedef boost::shared_ptr<const CheyetteBaseCalibrator> CheyetteBaseCalibrator_CONSTPTR;