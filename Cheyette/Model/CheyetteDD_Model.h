#pragma once

#include <Cheyette/Model/CheyetteModel.h>
#include <Instrument/VanillaSwaption.h>		//pour SwapRate Version

#include <vector>
#include <iostream>

#include <Cheyette/Fonction.h>
#include <boost/function.hpp>


//#include <Instrument/VanillaSwap.h> 
//#include <LMM/helper/LMMTenorStructure.h>
//
//#include <ql/math/array.hpp>

/*  -----------------------------------------------------------
	Cheyette Model Displaced Diffusion (DD)

	dx(t) = (y(t) - k * x(t) ) dt + \sigma_r(t) dW(t)       //k_ : parametre (constant) de mean reversion)      
	y(t) = ( \sigma^2_r(t) - 2 k y(t) ) dt 

Plus générique :
	\sigma_r(t) = (m(t) R_t(t, x, y) + (1 - m(t)) R_0(t, x, y) ) \sigma(t)  := sigma(t, x) 

	2 parameters: m(t) and sigma(t) assumed to be piecewise constant time-dependent
	m(t) \in [0,1]
   -----------------------------------------------------------*/


class CheyetteDD_Model  : public CheyetteModel
{

protected :
	double						k_;			//mean reversion parameter : Constant
	Piecewiseconst_RR_Function	sigma_ ;	//sigma(t)				
	Piecewiseconst_RR_Function	m_;			//m(t)					

public :
	//Ctor
	CheyetteDD_Model(	const CourbeInput_PTR courbeInput_PTR, 
						const double k,		
						const Piecewiseconst_RR_Function&	sigma,
						const Piecewiseconst_RR_Function&	m) 
		: CheyetteModel(courbeInput_PTR), k_(k), sigma_(sigma), m_(m){}

	//destructor
	virtual ~CheyetteDD_Model(){};

	virtual double meanReversion(double t, double x_t, double y_t) const {return k_;}
	virtual double localVol(double t, double x_t, double y_t) const
	{
		return sigma_(t) * ( m_(t)*R_t_(t,x_t,y_t) + (1-m_(t))*R_0_(t,x_t,y_t) );
	}

	virtual double R_t_(double t, double x_t, double y_t) const = 0;
	virtual double R_0_(double t, double x_t, double y_t) const = 0;

		//setters		//Ok, les asserts sont dans Fonction.h
	virtual void	setCheyetteModel_Parameter_Skew(const std::vector<double>& m_y){m_.sety_(m_y) ;}
	virtual void	setCheyetteModel_Parameter_Level(const std::vector<double>& sigma_y){sigma_.sety_(sigma_y) ;}
	virtual void	setCheyetteModel_Parameter_Skew(double m, size_t index){m_.sety_(m, index) ;}
	virtual void	setCheyetteModel_Parameter_Level(double sigma, size_t index){sigma_.sety_(sigma, index) ;}

	virtual void	setCheyetteModel_Parameter_Convexity(const std::vector<double>& c_y)
						{throw std::string("no convexity parameter in Displaced Diffusion version") ;}
	virtual void	setCheyetteModel_Parameter_Convexity(double c, size_t index)
						{throw std::string("no convexity parameter in Displaced Diffusion version") ;}

		//getters
	double						get_k()		const {return k_;}			
	Piecewiseconst_RR_Function	get_sigma() const {return sigma_ ;}	
	Piecewiseconst_RR_Function	get_m()		const {return m_ ;}			

	double G(double t, double T) const ; 
	void show() const ;
	void print(std::ostream& o) const ;
	
	virtual std::string getModelType() const = 0 ;
};

typedef boost::shared_ptr<CheyetteDD_Model>       CheyetteDD_Model_PTR;
typedef boost::shared_ptr<const CheyetteDD_Model> CheyetteDD_Model_CONSTPTR;


/*  -----------------------------------------------------------
	Cheyette Model Displaced Diffusion (DD) Short Rate Version

	dx(t) = (y(t) - k * x(t) ) dt + \sigma_r(t) dW(t)       //k_ : parametre (constant) de mean reversion)      
	y(t) = ( \sigma^2_r(t) - 2 k y(t) ) dt 

Plus générique :
	\sigma_r(t) = (m(t) r_t(t, x, y) + (1 - m(t)) r_0 ) \sigma(t)  := sigma(t, x) 

	2 parameters: m(t) and sigma(t) assumed to be piecewise constant time-dependent
	m(t) \in [0,1]
   -----------------------------------------------------------*/

//r(t) / r(0)
class CheyetteDD_Model_ShortRateVersion : public CheyetteDD_Model   
{
public: 

	double R_t_(double t, double x_t, double y_t) const {return r_t(t, x_t) ;}
	double R_0_(double t, double x_t, double y_t) const {return r_t(0., 0.) ;}


	CheyetteDD_Model_ShortRateVersion(	const CourbeInput_PTR courbeInput_PTR, 
										const double k, 
										const Piecewiseconst_RR_Function& sigma, 
										const Piecewiseconst_RR_Function& m) 
			: CheyetteDD_Model(courbeInput_PTR, k, sigma, m){}

	//destructor
	virtual ~CheyetteDD_Model_ShortRateVersion(){};

	virtual double localVol_1stDerivative(double t, double x_t, double y_t) const
	{
		//r(t) = x(t) + f(0, t)
		return sigma_(t) * m_(t) * 1. ;	
	} 

	virtual double localVol_2ndDerivative(double t, double x_t, double y_t) const
	{
		return 0. ;	
	} 

	virtual std::string getModelType() const {return "Version DD : Short Rate Version r(t) / r(0)" ;}
};

typedef boost::shared_ptr<CheyetteDD_Model_ShortRateVersion>       CheyetteDD_Model_ShortRateVersion_PTR;
typedef boost::shared_ptr<const CheyetteDD_Model_ShortRateVersion> CheyetteDD_Model_ShortRateVersion_CONSTPTR;


/*  -----------------------------------------------------------
	Cheyette Model Displaced Diffusion (DD) Forward Rate Version

	dx(t) = (y(t) - k * x(t) ) dt + \sigma_r(t) dW(t)       //k_ : parametre (constant) de mean reversion)      
	y(t) = ( \sigma^2_r(t) - 2 k y(t) ) dt 

Plus générique :
	\sigma_r(t) = (m(t) r_t(t, x, y) + (1 - m(t)) f(0,t) ) \sigma(t)  := sigma(t, x) 

	2 parameters: m(t) and sigma(t) assumed to be piecewise constant time-dependent
	m(t) \in [0,1]
   -----------------------------------------------------------*/

//f(t,t) = r(t) / f(0, t)
class CheyetteDD_Model_ForwardRateVersion : public CheyetteDD_Model   
{
public: 
	double R_t_(double t, double x_t, double y_t) const {return r_t(t, x_t) ;}
	double R_0_(double t, double x_t, double y_t) const {return courbeInput_PTR_->get_f_0_t(t) ;}


	CheyetteDD_Model_ForwardRateVersion(const CourbeInput_PTR courbeInput_PTR, 
										const double k, 
										const Piecewiseconst_RR_Function& sigma, 
										const Piecewiseconst_RR_Function& m) 
			: CheyetteDD_Model(courbeInput_PTR, k, sigma, m){}

	//destructor
	virtual ~CheyetteDD_Model_ForwardRateVersion(){};

	virtual double localVol_1stDerivative(double t, double x_t, double y_t) const 
	{
		//r(t) = x(t) + f(0, t)
		return sigma_(t) * m_(t) * 1. ;	
	} 

	virtual double localVol_2ndDerivative(double t, double x_t, double y_t) const
	{
		return 0. ;	
	} 

	virtual std::string getModelType() const {return "Version DD : Forward Rate Version f(t, t) / f(0, t)" ;}
};

typedef boost::shared_ptr<CheyetteDD_Model_ForwardRateVersion>       CheyetteDD_Model_ForwardRateVersion_PTR;
typedef boost::shared_ptr<const CheyetteDD_Model_ForwardRateVersion> CheyetteDD_Model_ForwardRateVersion_CONSTPTR;

/*  -----------------------------------------------------------
	Cheyette Model Displaced Diffusion (DD) Swap Rate Version

	dx(t) = (y(t) - k * x(t) ) dt + \sigma_r(t) dW(t)       //k_ : parametre (constant) de mean reversion)      
	y(t) = ( \sigma^2_r(t) - 2 k y(t) ) dt 

Plus générique :
	\sigma_r(t) = (m(t) S_t(t, x, y) + (1 - m(t)) S_0 ) \sigma(t)  := sigma(t, x) 

	2 parameters: m(t) and sigma(t) assumed to be piecewise constant time-dependent
	m(t) \in [0,1]
   -----------------------------------------------------------*/

//S(t) /S(0) swap rate

/* !!!!!!!

vérifier que swaption passée en parametre de Cheyette DD Model soit la même que la swaption pour le pricer

mettre un attribut approxPricer en parametre et creer une fonction init qui se servira des fonctions du approxPricer 
pour calculer local Vol et les dérivées

!!!!!!! */

class CheyetteDD_Model_SwapRateVersion : public CheyetteDD_Model   
{
private :
	VanillaSwaption_PTR		pSwaption_ ;
public: 
	CheyetteDD_Model_SwapRateVersion(const CourbeInput_PTR courbeInput_PTR, 
									const double k, 
									const Piecewiseconst_RR_Function& sigma, 
									const Piecewiseconst_RR_Function& m, 
									const VanillaSwaption_PTR	pSwaption) 
		: CheyetteDD_Model(courbeInput_PTR, k, sigma, m), pSwaption_(pSwaption){}

	//destructor
	virtual ~CheyetteDD_Model_SwapRateVersion(){};

	virtual double R_t_(double t, double x_t, double y_t) const {return swapRate(t, x_t, y_t, pSwaption_) ;}
	virtual double R_0_(double t, double x_t, double y_t) const {return swapRate0(pSwaption_) ;}

	virtual double localVol_1stDerivative(double t, double x_t, double y_t) const
	{
		return sigma_(t) * m_(t) * swapRate_1stDerivative(t, x_t, y_t, pSwaption_) ;	
	} 

	virtual double localVol_2ndDerivative(double t, double x_t, double y_t) const
	{
		std::cout << "fonction local vol 2nd derivative non definie dans Cheyette Swap rate Version !! " << std::endl ;
		return 0. ;	
	} 
	virtual std::string getModelType() const {return "Version DD : Swap Rate Version S(t) / S(0)" ;}

	VanillaSwaption_PTR	getSwaption() const {return pSwaption_ ;}
	void setSwaption(VanillaSwaption_PTR pSwaption) {pSwaption_ = pSwaption ;}

	void show() const ;
};

typedef boost::shared_ptr<CheyetteDD_Model_SwapRateVersion>       CheyetteDD_Model_SwapRateVersion_PTR;
typedef boost::shared_ptr<const CheyetteDD_Model_SwapRateVersion> CheyetteDD_Model_SwapRateVersion_CONSTPTR;










