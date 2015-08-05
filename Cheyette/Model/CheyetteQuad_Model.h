#pragma once

#include <Cheyette/Fonction.h>
#include <Cheyette/Model/CourbeInput.h>

#include <vector>
#include <iostream>

#include <boost/function.hpp>
#include <boost/shared_ptr.hpp>

#include <Instrument/VanillaSwap.h> 
#include <LMM/Helper/LMMTenorStructure.h>

#include <ql/math/array.hpp>

/*  -----------------------------------------------------------
	Cheyette Model with Quadratic Volatility

	dx(t) = (y(t) - k * x(t) ) dt + \sigma_r(t) dW(t)       //k_: mean reversion     
	y(t) = ( \sigma^2_r(t) - 2 k y(t) ) dt 

	\sigma_r(t) = a(t) + b(t) x(t) + c(t) x(t)^2
	
	a(t), b(t) and c(t) are piecewise constant
   -----------------------------------------------------------*/


class CheyetteQuad_Model
{
public:

	struct CheyetteQuad_Parameter
	{
		double						k_;			//mean reversion parameter (constant)
	    Piecewiseconst_RR_Function	a_ ;		//a(t)				
		Piecewiseconst_RR_Function	b_;			//b(t)					
		Piecewiseconst_RR_Function	c_;			//c(t)					

		//structure constructor
		CheyetteQuad_Parameter(const double k, 
			                 const Piecewiseconst_RR_Function& a, 
			                 const Piecewiseconst_RR_Function& b, 
							 const Piecewiseconst_RR_Function& c)
			: k_(k), a_(a), b_(b), c_(c) 
			{
				assert(k > 0) ;
				if(a_.getx_() != b_.getx_() && b_.getx_() != c_.getx_())
				{
					std::cout << "Warning: pillars of piecewise constant functions a, b and c are different" << std::endl ;
				}
			}
	};

private:

	CourbeInput_PTR					courbeInput_PTR_ ;             // yield y(0,t)
	mutable CheyetteQuad_Parameter	cheyetteQuad_Parameter_; 

public :
	//constructor
	CheyetteQuad_Model(	const CourbeInput_PTR& courbeInput_PTR, 
						const CheyetteQuad_Parameter& cheyetteParam)
				: courbeInput_PTR_(courbeInput_PTR), cheyetteQuad_Parameter_(cheyetteParam) {}

	virtual ~CheyetteQuad_Model(){}

	//getters
	CourbeInput_PTR			get_courbeInput_PTR() const{return courbeInput_PTR_ ;}
	CheyetteQuad_Parameter	get_CheyetteQuad_Parameter() const{return cheyetteQuad_Parameter_;}

	//setters
	void		setCheyetteQuad_Parameter_a(const std::vector<double>& a_y)
						{cheyetteQuad_Parameter_.a_.sety_(a_y) ;}
	void		setCheyetteQuad_Parameter_b(const std::vector<double>& b_y)
						{cheyetteQuad_Parameter_.b_.sety_(b_y) ;}
	void		setCheyetteQuad_Parameter_c(const std::vector<double>& c_y)
						{cheyetteQuad_Parameter_.c_.sety_(c_y) ;}

	void		setCheyetteQuad_Parameter_a(double a, size_t index)
						{cheyetteQuad_Parameter_.a_.sety_(a, index) ;}
	void		setCheyetteQuad_Parameter_b(double b, size_t index)
						{cheyetteQuad_Parameter_.b_.sety_(b, index) ;}
	void		setCheyetteQuad_Parameter_c(double c, size_t index)
						{cheyetteQuad_Parameter_.c_.sety_(c, index) ;}

	void show() const ;
	void print(std::ostream& o) const ;

	//fonction de vol locale Quadratic Volatility
	double sigma_r( double t,  double x_t) const ;
	double sigma_r_t_1stDerivative( double t,  double x_t) const ;  //derivee wrt x_t
	double sigma_r_t_2ndDerivative( double t,  double x_t) const ;

	//fonctions G(t, T), ZC B(t, T), Libor...
	double G(double t, double T) const ;  
	double P(double t, double T, double x_t, double y_t) const ;  //ZC B(0,t) : donné en input sinon stochastique

	double r_t(double t, double x_t) const ;

	double Libor(double t, double T1, double T2, double x_t, double y_t) const ;  

		//EDS : diffusion
	double diffusion_x(double t, double x_t) const ;
	
	double drift_y(double t, double x_t, double y_t) const ;

	//EDS : drift sous Q^T
	double drift_x_QT(double t, double T_proba_fwd, double x_t, double y_t) const ;

};

typedef boost::shared_ptr<CheyetteQuad_Model>       CheyetteQuad_Model_PTR;
typedef boost::shared_ptr<const CheyetteQuad_Model> CheyetteQuad_Model_CONSTPTR;