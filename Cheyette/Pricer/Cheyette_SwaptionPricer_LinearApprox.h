#pragma once

#include <Cheyette/Pricer/Cheyette_SwaptionPricer_Approx.h>

class Cheyette_SwaptionPricer_LinearApprox : public Cheyette_SwaptionPricer_Approx
{
private :
	mutable double buffer_b_barre_ ;

public :
	//ctor
	Cheyette_SwaptionPricer_LinearApprox(	const CheyetteModel_PTR pCheyette_Model) 
						: Cheyette_SwaptionPricer_Approx(pCheyette_Model)
	{}
	
	//virtual dtor
	virtual ~Cheyette_SwaptionPricer_LinearApprox(){}
	
	virtual void preCalculateALL(VanillaSwaption_PTR pSwaption) const;
	virtual void preCalculateModelBuffers() const ;


/***************** Linear approx of swap rate volatility *********************
*
*	dS(t) = Phi(t, S_t) dW^Q_A
*	dS(t) \approx Phi(t, s_bar) + DPhi(t, s_bar) (S_t - s_bar) dW^Q_A
*
******************************************************************************/
		virtual double phi(double t, double x_t) const; 
		double dPhi(double t, double x_t) const; 


/****************  parameter averaging  *********************
*
*	dS(t) \approx Phi(t, s_bar) + DPhi(t, s_bar) (S_t - s_bar) dW^Q_A
*	dS(t) = \lambda(t) ( S0 + b(t) (S(t) - S(0) ) dW^Q_A		-> time averaging: b(t) puis b_barre
*
************************************************************/

		double lambda(double t) const;
		double lambda2(double t) const;
		double b(double t) const;

		double f_outer_num(double t) const ;		//lambda^2(u) * b(u)
		double f_outer_denom(double t) const ;		//lambda^2(u)
		
		//retourne b barre du displaced diffusion
		double timeAverage(double t) const ;

		//prix swaption approximé 
		virtual double price(const VanillaSwaption_PTR pSwaption) const ;		

};

typedef boost::shared_ptr<Cheyette_SwaptionPricer_LinearApprox>       Cheyette_SwaptionPricer_LinearApprox_PTR;
typedef boost::shared_ptr<const Cheyette_SwaptionPricer_LinearApprox> Cheyette_SwaptionPricer_LinearApprox_CONSTPTR;