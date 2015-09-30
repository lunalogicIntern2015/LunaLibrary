#pragma once

#include <Cheyette/Model/CourbeInput.h>
#include <Instrument/VanillaSwaption.h>
#include <boost/shared_ptr.hpp>

//1D 
class CheyetteModel  
{
protected :

	CourbeInput_PTR					courbeInput_PTR_ ;          // yield y(0,t)  

public :
	CheyetteModel(CourbeInput_PTR courbeInput_PTR) : courbeInput_PTR_(courbeInput_PTR){}
	
	virtual ~CheyetteModel(){};

	virtual double meanReversion(double t, double x_t, double y_t) const = 0; 
	virtual double localVol(double t, double x_t, double y_t) const = 0; 
	virtual double localVol_1stDerivative(double t, double x_t, double y_t) const = 0; 
	virtual double localVol_2ndDerivative(double t, double x_t, double y_t) const = 0; 

	virtual void	setCheyetteModel_Parameter_Skew(const std::vector<double>&) = 0 ;
	virtual void	setCheyetteModel_Parameter_Level(const std::vector<double>&) = 0 ;
	virtual void	setCheyetteModel_Parameter_Convexity(const std::vector<double>&) = 0 ;

	virtual void	setCheyetteModel_Parameter_Skew(double, size_t index) = 0 ;
	virtual void	setCheyetteModel_Parameter_Level(double, size_t index) = 0 ;
	virtual void	setCheyetteModel_Parameter_Convexity(double, size_t index) = 0 ;

	//getter
	CourbeInput_PTR				get_courbeInput_PTR() const{return courbeInput_PTR_ ;}

	virtual void show() const
	{
		std::cout << "---------------------------------------------" << std::endl ;
		std::cout << "--- creation d'un objet Cheyette Model ---" << std::endl ;
	}

	virtual void print(std::ostream& o) const 
	{
		o << "---------------------------------------------" << std::endl ;
		o << "--- creation d'un objet Cheyette Model ---" << std::endl ;	
	}

	//methodes G(t, T), ZC B(t, T) et taux court r(t)
	virtual double G(double t, double T) const = 0 ;  
	double P(double t, double T, double x_t, double y_t) const 
	{
		double g = G(t,T) ;
		double P_0_t = exp( - courbeInput_PTR_->get_tauxZC0(t) * t ) ;
		double P_0_T = exp( - courbeInput_PTR_->get_tauxZC0(T) * T ) ;
		double res = P_0_T/P_0_t * exp(- x_t * g - 0.5 * y_t * g * g) ;
		return res ; 
	}

	double r_t(double t, double x_t) const	
	{
		double f_0_t = courbeInput_PTR_->get_f_0_t(t) ;
		return f_0_t + x_t ; 
	}

	//T2 = T1 + delta
	double libor(double t, double T1, double T2, double x_t, double y_t) const
	{
		double delta = T2 - T1 ;
		double P_t_T1 = P(t, T1, x_t, y_t) ;
		double P_t_T2 = P(t, T2, x_t, y_t) ;
		return 1/delta * (P_t_T1/P_t_T2 - 1.0) ;
	}

	//EDS : drift et diffusion sous Q^T
	double diffusion_x(double t, double x_t, double y_t) const
	{
		return localVol(t, x_t, y_t) ;
	}

	double drift_x_Q(double t, double x_t, double y_t) const                       // risk neutral proba
	{
		double k		= meanReversion(t, x_t, y_t) ;
		return y_t - k * x_t;	  	 
	}

	double drift_x_QT(double t, double T_proba_fwd, double x_t, double y_t) const
	{
		double k		= meanReversion(t, x_t, y_t) ;
		double sigma	= localVol(t, x_t, y_t) ;

		return y_t - k * x_t - G(t, T_proba_fwd) * sigma * sigma ;	  
	}


	double drift_y(double t, double x_t, double y_t) const
	{
		double k		= meanReversion(t, x_t, y_t) ;
		double sigma	= localVol(t, x_t, y_t) ;
		return sigma * sigma - 2 * k * y_t ;
	}

	//pour la version SwapRate de Displaced Diffusion
//sale, code en doublon avec approxPricer
//faudrait creer des classes ApproxPricerDD_linear, ApproxPricerDD_Quad, ...
	double swapRateNumerator(double t, double x_t, double y_t, const VanillaSwaption_PTR pSwaption) const ;
	double swapRateDenominator(double t, double x_t, double y_t, const VanillaSwaption_PTR pSwaption) const ;
	double swapRate(double t, double x_t, double y_t, const VanillaSwaption_PTR pSwaption) const ;
	double swapRate0(const VanillaSwaption_PTR pSwaption) const ;

	double ZC_1stDerivative_on_xt(double t, double T, double x_t, double y_t, const VanillaSwaption_PTR pSwaption) const ;
	double swapRateNumerator_1stDerivative(double t, double x_t, double y_t, const VanillaSwaption_PTR pSwaption) const ;
	double swapRateDenominator_1stDerivative(double t, double x_t, double y_t, const VanillaSwaption_PTR pSwaption) const;

	double swapRate_1stDerivative(double t, double x_t, double y_t, const VanillaSwaption_PTR pSwaption) const;
};

typedef boost::shared_ptr<CheyetteModel>       CheyetteModel_PTR;
typedef boost::shared_ptr<const CheyetteModel> CheyetteModel_CONSTPTR;











