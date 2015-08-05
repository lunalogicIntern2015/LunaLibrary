//#pragma once
//
//#include <ql/types.hpp>
//#include <ql/math/array.hpp>
//#include <ql/math/optimization/costfunction.hpp>
//
//#include <LMM/Model/VolatilityFunction.h>
//#include <LMM/Pricer/LmmApproximationPricer/LmmVanillaSwaptionApproxPricer.h>
//#include <LMM/Calibration/ATMSwaptionStructuredData.h>
//
//
//
//class RebonatoVolatilityCostFunction: public QuantLib::CostFunction
//{
//public:
//
//	typedef std::vector<std::vector<double> > matrix_;
//
//	RebonatoVolatilityCostFunction
//		(
//		VolatilityParam_PTR pVolatilityParam,
//		ATMSwaptionStructuredData_PTR pATMSwaptionStructuredData_,
//		const LmmApproxVanillaSwaptionPricer& approximation,						 
//		const std::vector<double>& libor_shifts,
//		matrix_ weights,
//		matrix_ weights_maturity,
//		matrix_ weights_tenor,
//		matrix_ weights_maturity_2,
//		matrix_ weights_tenor_2
//		);
//
//	~RebonatoVolatilityCostFunction();
//
//	QuantLib::Real value(const Array & x) const; 
//	QuantLib::Disposable<QuantLib::Array> values(const QuantLib::Array& x) const; 
//
//	//ctntodo to delete , not need seen new structure data are all vectors
//	matrix_ map_ArrayToMatrix(const QuantLib::Array& x)   const ;
//	QuantLib::Array map_MatrixtoArray(const matrix_& mat) const ;
//
//
//private:
//	//! pointer to a vol param in order to hold the actual value of the volatility parameters
//	mutable VolatilityParam_PTR pVolatilityParamBuffer_;
//
//	ATMSwaptionStructuredData_PTR pATMSwaptionStructuredData_;
//
//	mutable LmmApproxVanillaSwaptionPricer approximation_; // Class implementing Rebonato's formula 
//
//	matrix_ weights_; // weights for regulation 
//	matrix_ weights_maturity_; // weights for regulation (1st order derivative with respect to maturity)
//	matrix_ weights_tenor_; // weights for regulation (1st order derivative with respect to tenor)
//	matrix_ weights_maturity_2_; // weights for regulation (2nd order derivative with respect to maturity)
//	matrix_ weights_tenor_2_; // weights for regulation (2nd order derivative with respect to tenor)s
//
//
//
//	//-- Compute regulation term 
//	//-- coefficients ci control regulation terms
//	QuantLib::Real regularisation(const QuantLib::Array& x, QuantLib::Real c1, QuantLib::Real c2, QuantLib::Real c3, QuantLib::Real c4) const;  // TODO: changer nom coef 
//	QuantLib::Real sum_all_weights_regularisation(const matrix_& weights) const;
//	//Real sum_all_derivatives_regularisation(const matrix_& weights, const matrix_& derivatives);
//
//};
//
