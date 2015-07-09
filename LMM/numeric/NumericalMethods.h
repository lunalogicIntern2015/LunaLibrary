#pragma once

#include <vector>
#include <iostream>


#include <boost/math/distributions.hpp>
#include <boost/function.hpp>
#include <boost/bind.hpp>

#include <ql/math/distributions/normaldistribution.hpp>



//TODO: delete const double &, use only double .
namespace NumericalMethods
{
	//! finite difference method for univariate function f
	double df(const double& f_pp, const double& f_mm, const double& dx);
	
	double d1(const double& fwd, const double& strike, const double& vol, const double& T);

	double d2(const double& fwd, const double& strike, const double& vol, const double& T);

	// Compute a derivative's price using Black's formula
	double Black_Price(double fwd, double strike, double vol, double T);

	double swaptionBlack_Price(double annuity0, double fwd, double strike, double vol, double T) ;

	//avec vol non constante
	double Black_Price_vol2(double fwd, double strike, double vol_T, double T);

	double Black_Price_vol2_allStrike(double fwd, double strike, double sqrt_int_sigma2, double T) ;

	double Black_Vega(const double fwd, const double strike, const double vol, const double T);
	double Black_Vega_swaption(const double annuity0, const double fwd, const double strike, const double vol, const double T) ;

	double Black_Volga(const double& fwd, const double& strike, const double& vol, const double& T);

	double Black_impliedVolatility(const double bs_call_price, const double fwd, const double strike, const double T);
	double Black_SwaptionImpliedVolatility(const double bs_call_price, const double annuity0,   							  
										const double  fwd, const double  strike, const double  T) ;

	double linearInterpolation(
							   const double& t, 
		                       const std::vector<double>& maturities,
		                       const std::vector<double>& set_of_points
							   );

	double linearInterpolation2(
							const double& t, 
		                    const std::vector<double>& maturities,
		                    const std::vector<double>& set_of_points
							);

	double vectorProduct(std::vector<double>& v1, std::vector<double>& v2) ;
	
	//spline cubique
	void spline(const std::vector<double>& x, const std::vector<double>& y, 
			double yp1, double ypn, std::vector<double>& y2) ;
	double splineCubique(const std::vector<double>& xa, const std::vector<double>& ya, 
					 const std::vector<double>& y2a, double x) ;
} // end NumericalMethods

