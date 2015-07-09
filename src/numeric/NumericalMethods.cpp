#include <cassert>
#include <LMM/numeric/NumericalMethods.h>

namespace NumericalMethods
{
	
struct BS_function_helper
{
private:
	boost::function<double(double)> _f;
	boost::function<double(double)> _f_derivative;
	double _target;
	
public:
	//! constructor
	BS_function_helper(boost::function<double(double)>& f,
		               boost::function<double(double)>& f_derivative, 
					   double target) 
					   : _f(f), _f_derivative(f_derivative), _target(target){}

	//! evaluation
	boost::math::tuple<double, double> operator()(double x)
	{
		return boost::math::make_tuple(
			_f(x) - _target,
			_f_derivative(x));
	}
};

double df(const double& f_pp, const double& f_mm, const double& dx)
{
	assert(dx>0);
	double df = f_pp - f_mm;
	return df/(2*dx);
}

double d1(const double& fwd, const double& strike, const double& vol, const double& T)
{
	double variance = vol*vol*T; 
	double d1 = (log(fwd/strike) + 0.5*variance)/sqrt(variance);
	return d1;
}

double d2(const double& fwd, const double& strike, const double& vol, const double& T)
{
	double variance = vol*vol*T; 

	double d1 = (log(fwd/strike) + 0.5*variance)/sqrt(variance);
	double d2 = d1 - sqrt(variance);
	return d2;
}


//! without discount factor: r=0
double Black_Price(double fwd, double strike, double vol, double T) 
{
    assert(vol > 0 && T > 0 && fwd >0 && strike >0);

	double variance = vol*vol*T; 

	double d1 = (log(fwd/strike) + 0.5*variance)/sqrt(variance);
	double d2 = d1 - sqrt(variance);

	boost::math::normal_distribution<> nd(0,1); 
	double N1 = cdf(nd,d1);
	double N2 = cdf(nd,d2); 

	return fwd*N1-strike*N2;
}

//swaption Black price
//dS(t) / S(t) = sigma dW(t)^QA
double swaptionBlack_Price(double annuity0, double fwd, double strike, double vol, double T) 
{
    assert(vol > 0 && T > 0 && fwd >0 && strike >0);

	double variance = vol*vol*T; 

	double d1 = (log(fwd/strike) + 0.5*variance)/sqrt(variance);
	double d2 = d1 - sqrt(variance);

	boost::math::normal_distribution<> nd(0,1); 
	double N1 = cdf(nd,d1);
	double N2 = cdf(nd,d2); 

	return annuity0 * (fwd*N1-strike*N2) ;
}

//avec vol non constante
//sqrt_int_sigma2 = sqrt( \int_0^T \sigma^2(u) du )   
double Black_Price_vol2(double fwd, double strike, double sqrt_int_sigma2, double T) 
{
    assert(sqrt_int_sigma2 > 0 && T > 0 && fwd >0 && strike >0);

	double d1 = (log(fwd/strike) + 0.5 * sqrt_int_sigma2 * sqrt_int_sigma2) / sqrt_int_sigma2 ;
	double d2 = d1 - sqrt_int_sigma2;

	boost::math::normal_distribution<> nd(0,1); 
	double N1 = cdf(nd,d1);
	double N2 = cdf(nd,d2); 

	return fwd*N1-strike*N2;
}

//prix Black avec K < 0 
//cas où le sous jacent est martingale sous la proba voulue
//E^{Q annuity}(S(T0) | F_0) = S0 
double Black_Price_vol2_allStrike(double fwd, double strike, double sqrt_int_sigma2, double T) 
{
	std::cout << "SJ : " << fwd << ", strike : " << strike << ", T : " << T << std::endl ;
    assert(sqrt_int_sigma2 > 0 && T > 0 && fwd >0) ; 
	double res  ;
	if (strike > 0)
	{
		res = Black_Price_vol2(fwd, strike, sqrt_int_sigma2, T) ;
	}
	else
	{
		res = fwd - strike ;	
	}

	return res ;
}

double Black_Vega(const double fwd, const double strike, const double vol, const double T)
{
	double variance = vol*vol*T; 

	double d1 = ( log(fwd/strike) + 0.5*variance) / sqrt(variance);

	QuantLib::NormalDistribution ND;
	return fwd*ND(d1)*sqrt(T);
}

double Black_Vega_swaption(const double annuity0, const double fwd, const double strike, const double vol, const double T)
{
	double variance = vol*vol*T; 

	double d1 = ( log(fwd/strike) + 0.5*variance) / sqrt(variance);

	QuantLib::NormalDistribution ND;
	return annuity0 * fwd*ND(d1)*sqrt(T);
}

double Black_Volga(const double& fwd, const double& strike, const double& vol, const double& T) 
{
	
	double variance = vol*vol*T; 
	double d1 = ( log(fwd/strike) + 0.5*variance ) / sqrt(variance);
	double d2 = d1 - sqrt(variance);

	double vega = Black_Vega(fwd,strike,vol,T);
	return vega*(d1*d2/vol);
}



double Black_impliedVolatility(const double  bs_call_price, 							  
							   const double  fwd, 
							   const double  strike,
							   const double  T)
{
	//throw("Error bug, it doesn't work."); // TODO, correct it latter. 
	//-- !! variable _1 represents the SQUARED volatility 
	boost::function<double(double)> f = boost::bind(Black_Price,fwd,strike,_1,T);
	boost::function<double(double)> f_derivative = boost::bind(Black_Vega,fwd,strike,_1,T);
	
	BS_function_helper bs_function_helper(f,f_derivative,bs_call_price);
	double initial_guess = 1 ;
	double min    = 10e-5;
	double max    = 10 ;
    size_t nDigits   = 15;
	boost::uintmax_t nMaxIter  = 100;
	double result_newton_raphson = boost::math::tools::newton_raphson_iterate(bs_function_helper, initial_guess, min, max, nDigits);
	return result_newton_raphson;

}

double Black_SwaptionImpliedVolatility(const double bs_call_price, const double annuity0,   							  
										const double  fwd, const double  strike, const double  T)
{
	boost::function<double(double)> f = boost::bind(swaptionBlack_Price, annuity0, fwd, strike, _1, T);
	boost::function<double(double)> f_derivative = boost::bind(Black_Vega_swaption, annuity0, fwd,strike,_1,T);   
	
	BS_function_helper bs_function_helper(f,f_derivative,bs_call_price);
	double initial_guess = 1 ;
	double min    = 0;
	double max    = 10 ;
    size_t nDigits   = 15;
	boost::uintmax_t nMaxIter  = 100;
	double result_newton_raphson = boost::math::tools::newton_raphson_iterate(bs_function_helper, initial_guess, min, max, nDigits);
	return result_newton_raphson;

}


//-- For now, we assume the first point of set is the T0=0 maturity rate (L(0,T0) or P(0,T0))
double linearInterpolation(const double& t, 
						   const std::vector<double>& maturities,
						   const std::vector<double>& set_of_points)
{
	
	size_t index_maturiy_before_t = 0;

	try{
		if (maturities.size() != set_of_points.size())
		{
			throw std::string("Exception linearInterpolation : vecteurs grid et value doivent avoir meme longueur");
		}
		else
		{
			//-- Search the maturities bounding date t
			for (size_t i = 0; i < maturities.size(); ++i)
			{
				if (t > maturities[i])
					index_maturiy_before_t++;
			}

			double date_prev = maturities[index_maturiy_before_t-1];
			double date_next = maturities[index_maturiy_before_t];

			double point_prev = set_of_points[index_maturiy_before_t-1];
			double point_next = set_of_points[index_maturiy_before_t];

			double coeff_1 = (date_next - t)/(date_next - date_prev);
			double coeff_2 = (t - date_prev)/(date_next - date_prev);

			double interpolatedValue = coeff_1*point_prev + coeff_2*point_next;
			return interpolatedValue;
		}
	}
	catch(std::string const& chaine) {std::cerr << chaine << std::endl;}

}

double linearInterpolation2(const double& t, 
						   const std::vector<double>& maturities,
						   const std::vector<double>& set_of_points)
{
	
	try{
		size_t N = maturities.size() ;
		if (N != set_of_points.size())
		{
			throw std::string("Exception linearInterpolation : vecteurs grid et value doivent avoir meme longueur");
		}
		if (N == 0)
		{
			throw std::string("Exception linearInterpolation : vecteurs doivent etre non vides");
		}
		if (t < maturities[0] || t > maturities[N - 1])
		{
			throw std::string("Exception linearInterpolation : extrapolation !!");
		}
		size_t i=0 ;
		while (t > maturities[i+1] && i < N-2){++i ;}
	
		return ( set_of_points[i] + (set_of_points[i+1] - set_of_points[i])/(maturities[i+1] - maturities[i])*(t - maturities[i]) ) ;
	}
	catch(std::string const& chaine) {std::cerr << chaine << std::endl;}

}

double vectorProduct(std::vector<double>& v1, std::vector<double>& v2)
{
    assert(v1.size() == v2.size());

	double result = 0.0;
	for(size_t i=0; i<v1.size(); ++i)
		result += v1[i] * v2[i];

	return result;
}

/*Given arrays x[1..n] and y[1..n] containing a tabulated function, i.e., yi = f(xi), with
x1 < x2 <...< xN , and given values yp1 and ypn for the first derivative of the interpolating
function at points 1 and n, respectively, this routine returns an array y2[1..n] that contains
the second derivatives of the interpolating function at the tabulated points xi. If yp1 and/or
ypn are equal to 1 × 10^30 or larger, the routine is signaled to set the corresponding boundary
condition for a natural spline, with zero second derivative on that boundary. */
void spline(const std::vector<double>& x, const std::vector<double>& y, 
			double yp1, double ypn, std::vector<double>& y2)
{
	double qn,un;
	size_t n = x.size() ;
	std::vector<double> u(n-1) ;

	if (yp1 > 0.99e30) //The lower boundary condition is set either to be “natural”
		y2[0]=u[0]=0.0;
	else { //or else to have a specified first derivative.
		y2[0] = -0.5;
		u[0]=(3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
	}
	for (size_t i=1; i<n-1; ++i) { 
		//This is the decomposition loop of the tridiagonal algorithm.
		//y2 and u are used for temporary storage of the decomposed factors.
		double sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
		double p=sig*y2[i-1]+2.0;
		y2[i]=(sig-1.0)/p;
		u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
		u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
	}
	if (ypn > 0.99e30) //The upper boundary condition is set either to be “natural”
	{
		qn = 0.0;
		un = 0.0;
	}
	else 
	{				//or else to have a specified first derivative.
		qn=0.5;
		un=(3.0/(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
	}
	y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-1]+1.0);

	for (int k = n-2 ; k>=0 ; k--) //This is the backsubstitution loop of the tridiagonal algorithm.
	{
		y2[k]=y2[k]*y2[k+1]+u[k];   // --  
	}
}

/* Given the arrays xa[1..n] and ya[1..n], which tabulate a function (with the xai’s in order),
and given the array y2a[1..n], which is the output from spline above, and given a value of
x, this routine returns a cubic-spline interpolated value y. */
double splineCubique(const std::vector<double>& xa, const std::vector<double>& ya, 
					 const std::vector<double>& y2a, double x)
{
	//void nrerror(char error_text[]);
	double h,b,a;
	int klo = 1;
	int khi = xa.size() ;
	while (khi-klo > 1) 
	{
		int k=(khi+klo) >> 1;
		if (xa[k] > x) khi=k;
		else klo=k;
	}
	h=xa[khi]-xa[klo];
	if (h == 0.0) throw "Bad xa input to routine splint" ;
	a=(xa[khi]-x)/h;
	b=(x-xa[klo])/h;
	return a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
}
}// end NumericalMethods