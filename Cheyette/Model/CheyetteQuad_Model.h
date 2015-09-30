#pragma once

#include <Cheyette/Fonction.h>
#include <Cheyette/Model/CheyetteModel.h>
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


class CheyetteQuad_Model : public CheyetteModel 
{

private :
	double						k_;			//mean reversion parameter (constant)
	Piecewiseconst_RR_Function	a_ ;		//a(t)				
	Piecewiseconst_RR_Function	b_;			//b(t)					
	Piecewiseconst_RR_Function	c_;			//c(t)					

public :
	//Ctor
	CheyetteQuad_Model(	const CourbeInput_PTR courbeInput_PTR, 
						double k,		
			            const Piecewiseconst_RR_Function& a, 
			            const Piecewiseconst_RR_Function& b, 
						const Piecewiseconst_RR_Function& c) 
		: CheyetteModel(courbeInput_PTR), k_(k), a_(a), b_(b), c_(c)
		{
			assert(k > 0) ;
			if(a_.getx_() != b_.getx_() && b_.getx_() != c_.getx_())
			{
				std::cout << "Warning: pillars of piecewise constant functions a, b and c are different" << std::endl ;
			}
		}


	virtual ~CheyetteQuad_Model(){}

	virtual double meanReversion(double t, double x_t, double y_t) const {return k_;}
	virtual double localVol(double t, double x_t, double y_t) const 
	{
		return a_(t) + b_(t) * x_t + c_(t) * x_t * x_t ;
	}

	virtual double localVol_1stDerivative(double t, double x_t, double y_t) const
	{
		//r(t) = x(t) + f(0, t)
		return b_(t) + 2 * c_(t) * x_t ;
	} 

	virtual double localVol_2ndDerivative(double t, double x_t, double y_t) const
	{
		return 2 * c_(t) ;
	} 

	//getters
	double						get_k() const {return k_;}			
	Piecewiseconst_RR_Function	get_a() const {return a_ ;}	
	Piecewiseconst_RR_Function	get_b() const {return b_ ;}
	Piecewiseconst_RR_Function	get_c() const {return c_ ;}

	//setters
	virtual void		setCheyetteModel_Parameter_Level	(const std::vector<double>& a_y){a_.sety_(a_y) ;}
	virtual void		setCheyetteModel_Parameter_Skew		(const std::vector<double>& b_y){b_.sety_(b_y) ;}
	virtual void		setCheyetteModel_Parameter_Convexity(const std::vector<double>& c_y){c_.sety_(c_y) ;}

	virtual void		setCheyetteModel_Parameter_Level	(double a, size_t index){a_.sety_(a, index) ;}
	virtual void		setCheyetteModel_Parameter_Skew		(double b, size_t index){b_.sety_(b, index) ;}
	virtual void		setCheyetteModel_Parameter_Convexity(double c, size_t index){c_.sety_(c, index) ;}

	virtual void show() const ;
	virtual void print(std::ostream& o) const ;

	double G(double t, double T) const ;  

};

typedef boost::shared_ptr<CheyetteQuad_Model>       CheyetteQuad_Model_PTR;
typedef boost::shared_ptr<const CheyetteQuad_Model> CheyetteQuad_Model_CONSTPTR;