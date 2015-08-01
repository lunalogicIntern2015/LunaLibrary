#pragma once

#include <Cheyette/Fonction.h>
#include <Cheyette/CheyetteModel/CourbeInput.h>

#include <vector>
#include <iostream>

#include <boost/function.hpp>
#include <boost/shared_ptr.hpp>

#include <LMM/instrument/VanillaSwap.h> 
#include <LMM/helper/LMMTenorStructure.h>

#include <ql/math/array.hpp>
using namespace QuantLib ;
/*  -----------------------------------------------------------
	Cheyette Model Displaced Diffusion (DD)

	dx(t) = (y(t) - k * x(t) ) dt + \sigma_r(t) dW(t)       //k_ : parametre (constant) de mean reversion)      
	y(t) = ( \sigma^2_r(t) - 2 k y(t) ) dt 

	\sigma_r(t) = (m(t) r(t) + (1 - m(t)) shift(t) ) \sigma(t)  := sigma(t, x)  
	
	\sigma_r(t) = (m(t) shift(t, x, y) + (1 - m(t)) shift(0, x) ) \sigma(t)  := sigma(t, x) 

Plus générique :
	\sigma_r(t) = (m(t) shift1(t, x, y) + (1 - m(t)) shift2(t, x, y) ) \sigma(t)  := sigma(t, x) 

	2 parameters: m(t) and sigma(t) assumed to be piecewise constant time-dependent
	m(t) \in [0,1]
   -----------------------------------------------------------*/


class CheyetteDD_Model  
{
public:

	struct CheyetteDD_Parameter  // tout ce qu-il faudra calibrer
	{
		double						k_;			// add comment please	// c'est une constante à vérifier 
	    Piecewiseconst_RR_Function	sigma_ ;	//sigma(t)				
		Piecewiseconst_RR_Function	m_;			//m(t)					

		//structure constructor
		CheyetteDD_Parameter(const double k, 
			                 const Piecewiseconst_RR_Function& sigma, 
							 const Piecewiseconst_RR_Function& m)
				: k_(k), sigma_(sigma), m_(m) 
				{
					assert(k > 0) ;
					if(sigma_.getx_() != m_.getx_())
					{
						std::cout << "Warning: pillars of m differ from pillars of sigma" << std::endl ;
					}
				}
	};

private:

	CourbeInput_PTR					courbeInput_PTR_ ;          // yield y(0,t)
	mutable CheyetteDD_Parameter	cheyetteDD_Parameter_;
	int								shiftChoice_ ; 
	Boost_R3R_Function_PTR			shift1_;					//ex: r(t)/r(0), r(t)/f(0,t), S(t)/S(0) ou Li(t)/Li(0)
	Boost_R3R_Function_PTR			shift2_;
	Boost_R3R_Function_PTR			derivative_x_shift1_;
	Boost_R3R_Function_PTR			derivative_x_shift2_;

public:

	//constructor
	CheyetteDD_Model(const CourbeInput_PTR& courbeInput_PTR, const CheyetteDD_Parameter& cheyetteParam, int shiftChoice)
		: courbeInput_PTR_(courbeInput_PTR), cheyetteDD_Parameter_(cheyetteParam), shiftChoice_(shiftChoice) 
	{
		switch(shiftChoice) 
		{
			case 1:{		//r(t) / r(0)
				boost::function<double(double, double, double)> f1 
										= boost::bind(&CheyetteDD_Model::shift_rt, this, _1, _2, 0.);
				boost::function<double(double, double, double)> f2 
										= boost::bind(&CheyetteDD_Model::shift_rt, this, 0., 0., 0.);
				boost::function<double(double, double, double)> fp1 
										= boost::bind(&CheyetteDD_Model::shift_rt_1stDerivative, this, 0., 0., 0.);
				boost::function<double(double, double, double)> fp2 
										= boost::bind(&CheyetteDD_Model::shift_r0_1stDerivative, this, 0., 0., 0.);
				Boost_R3R_Function_PTR f1_ptr(new Boost_R3R_Function(f1)) ;
				Boost_R3R_Function_PTR f2_ptr(new Boost_R3R_Function(f2)) ;
				Boost_R3R_Function_PTR fp1_ptr(new Boost_R3R_Function(fp1)) ;
				Boost_R3R_Function_PTR fp2_ptr(new Boost_R3R_Function(fp2)) ;
				setShiftPointer(f1_ptr, f2_ptr, fp1_ptr, fp2_ptr) ;
				break ;
				   }
			case 2 :{	//r(t) / f(0, t)
				boost::function<double(double, double, double)> f1 
										= boost::bind(&CheyetteDD_Model::shift_rt, this, _1, _2, 0.);
				boost::function<double(double, double, double)> f2 
										= boost::bind(&CheyetteDD_Model::shift_f_0_t, this, _1, 0., 0.);
				boost::function<double(double, double, double)> fp1 
										= boost::bind(&CheyetteDD_Model::shift_rt_1stDerivative, this, 0., 0., 0.);
				boost::function<double(double, double, double)> fp2 
										= boost::bind(&CheyetteDD_Model::shift_f_0_t_1stDerivative, this, 0., 0., 0.);
				Boost_R3R_Function_PTR f1_ptr(new Boost_R3R_Function(f1)) ;
				Boost_R3R_Function_PTR f2_ptr(new Boost_R3R_Function(f2)) ;
				Boost_R3R_Function_PTR fp1_ptr(new Boost_R3R_Function(fp1)) ;
				Boost_R3R_Function_PTR fp2_ptr(new Boost_R3R_Function(fp2)) ;
				setShiftPointer(f1_ptr, f2_ptr, fp1_ptr, fp2_ptr) ;
				break ;
					}
			case 3 :{	//S(t) / S(0)
				//shift defini plus tard une fois connue la swaption
				//cf constructeur de CheyetteDD_VanillaSwaptionApproxPricer
				setShiftPointer(NULL, NULL, NULL, NULL) ;				
				break;
					}
			default :
				throw "Displaced Diffusion: shift not defined" ;
		}
	}
	
	//destructor
	virtual ~CheyetteDD_Model(){};

	//getters
	CourbeInput_PTR			get_courbeInput_PTR() const{return courbeInput_PTR_ ;}
	CheyetteDD_Parameter	get_CheyetteDD_Parameter() const{return cheyetteDD_Parameter_;}
	int						get_shiftChoice() const{return shiftChoice_ ;}

	//setters
	void		setCheyetteDD_Parameter_m(const std::vector<double>& m_y){cheyetteDD_Parameter_.m_.sety_(m_y) ;}
	void		setCheyetteDD_Parameter_sigma(const std::vector<double>& sigma_y){cheyetteDD_Parameter_.sigma_.sety_(sigma_y) ;}

	void		setCheyetteDD_Parameter_m(double m, size_t index){cheyetteDD_Parameter_.m_.sety_(m, index) ;}

	void		setCheyetteDD_Parameter_sigma(double sigma, size_t index){cheyetteDD_Parameter_.sigma_.sety_(sigma, index) ;}	

	//pointeurs vers f1, fprime1, f2, fprime2
	void setShiftPointer(	Boost_R3R_Function_PTR f1_ptr, Boost_R3R_Function_PTR f2_ptr, 
							Boost_R3R_Function_PTR fp1_ptr, Boost_R3R_Function_PTR fp2_ptr)
	{
		shift1_ = f1_ptr ;
		shift2_ = f2_ptr ;
		derivative_x_shift1_ = fp1_ptr;
		derivative_x_shift2_ = fp2_ptr;	
	}

	void show() const ;
	void print(std::ostream& o) const ;

	//fonction de vol locale Displaced Diffusion
	double sigma_r( double t,  double x_t, double y_t) const ;
	double sigma_r_t_1stDerivative( double t,  double x_t, double y_t) const ;  //derivee wrt x_t

	//shift functions
	double shift_rt(double t,  double x_t, double y_t) const ;							
	double shift_f_0_t(double t,  double x_t, double y_t) const ;				
	//derivee du shift wrt x_t
	double shift_rt_1stDerivative(double t, double x_t, double y_t) const ;  
	double shift_r0_1stDerivative(double t, double x_t, double y_t) const ;  
	double shift_f_0_t_1stDerivative(double t,  double x_t, double y_t) const ;  

	//fonctions G(t, T), ZC B(t, T), Libor...
	double G(double t, double T) const ;  
	double P(double t, double T, double x_t, double y_t) const ;  //ZC B(0,t) : donné en input sinon stochastique

	double r_t(double t, double x_t) const ;

	double Libor(double t, double T1, double T2, double x_t, double y_t) const ;  

	//EDS : diffusion
	double diffusion_x(double t, double x_t, double y_t) const ;
	
	double drift_y(double t, double x_t, double y_t) const ;

	//EDS : drift sous Q^T
	double drift_x_QT(double t, double T_proba_fwd, double x_t, double y_t) const ;

};

typedef boost::shared_ptr<CheyetteDD_Model>       CheyetteDD_Model_PTR;
typedef boost::shared_ptr<const CheyetteDD_Model> CheyetteDD_Model_CONSTPTR;











