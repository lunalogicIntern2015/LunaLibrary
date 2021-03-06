#pragma once

#include <ql/types.hpp>
#include <ql/math/array.hpp>
#include <ql/math/optimization/costfunction.hpp>

#include <Instrument/VanillaSwaption.h>
#include <LMM/Calibration/GMatrixMapping.h>
#include <LMM/UpperTriangleVanillaSwaptionQuotes.h>
#include <LMM/Calibration/LmmPenalty.h>
#include <LMM/Pricer/LmmApproximationPricer/LmmVanillaSwaptionApproxPricer_Rebonato.h>

#include <LMM/Model/Shifted_HGVolatilityFunction.h>


/*! \class LmmSkewCostFunction
 *
 *
 *
 */

class LmmSkewCostFunction : public QuantLib::CostFunction
{
public:

	//! constructor 
	LmmSkewCostFunction( LmmVanillaSwaptionApproxPricer_Rebonato_PTR lmmRobonato_ptr                // pricer  
				   , LiborQuotes_ConstPTR liborQuotes_ptr
				   , const double& quoted_strike_bump
		           , UpperTriangleVanillaSwaptionQuotes_ConstPTR pUpperTriangleVanillaSkewQuotes // instrument to calibrate 
				   , Shifted_HGVolatilityParam_PTR pShifted_HGVolatilityParam );
		
	QuantLib::Real value(const QuantLib::Array & x) const ; 
	
	//const Array& param_array
	QuantLib::Disposable<QuantLib::Array> values(const QuantLib::Array& x) const ; 

	// Getter!
	// attention, methods returning matrix copies, do not use in a loop
	UpperTriangularDoubleMatrix get_MarketQuotes() const;
	UpperTriangularDoubleMatrix get_CurrentModelQuotes() const;
	Shifted_HGVolatilityParam_PTR get_Shifted_HGVolatilityParam(); 

	//! use only for print calib evolution in the case of virtual calibration test
	//! in this case, we know the "true parameters" which is used for generating market data
	void reset_reference_calib(const QuantLib::Array & true_param); 

private :
	//! pricing 
	void update_SwaptionMdlSkewValues() const ;

	const double quoted_strike_bump_;
	LiborQuotes_ConstPTR pLiborQuotes_;
	UpperTriangleVanillaSwaptionQuotes_ConstPTR pUpperTriangleVanillaSkewQuotes_;	

	mutable LmmVanillaSwaptionApproxPricer_Rebonato_PTR pLmmVanillaSwaptionApproxPricer_Rebonato_;
	mutable UpperTriangularDoubleMatrix    mdl_swaption_skew;
	mutable Shifted_HGVolatilityParam_PTR buffer_Shifted_HGVolatilityParam_;	
	
	// output de calib debugging
	mutable unsigned int nbCalled;
	std::vector<double> buffer_calib_reference;// used for test
	bool breakForPrintOut(unsigned int nbIter) const { return (nbCalled % 5 == 0) ; }
	QuantLib::Array error_calib(const QuantLib::Array & actual_param) const;

	size_t calc_nbSwaptions() const ; 	
};

typedef boost::shared_ptr<LmmSkewCostFunction> LmmSkewCostFunction_PTR;