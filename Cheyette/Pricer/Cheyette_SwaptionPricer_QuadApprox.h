#pragma once

#include <Cheyette/Pricer/Cheyette_SwaptionPricer_Approx.h>

class Cheyette_SwaptionPricer_QuadApprox : public Cheyette_SwaptionPricer_Approx
{
private :
	mutable double						buffer_b_barre_ ;
	mutable double						buffer_c_barre_ ;

public :
	//ctor
	Cheyette_SwaptionPricer_QuadApprox(	const CheyetteModel_PTR pCheyette_Model) 
						: Cheyette_SwaptionPricer_Approx(pCheyette_Model){}
	//virtual dtor
	virtual ~Cheyette_SwaptionPricer_QuadApprox(){}
	
	virtual void preCalculateALL(VanillaSwaption_PTR pSwaption) const;
	virtual void preCalculateModelBuffers() const ;

//derivees 3eme par rapport à x_t
	double ZC_3rdDerivative_on_xt(double t, double T, double x_t, double y_t) const;
	double swapRateNumerator_3rdDerivative(double t, double x_t, double y_t) const;
	double swapRateDenominator_3rdDerivative(double t, double x_t, double y_t) const;	
	double swapRate_3rdDerivative(double t, double x_t, double y_t) const;

/**************** Quadratic approx of swap rate volatility ********************
*
*	dS(t) = phi(t, S_t) dW^Q_A
*	dS(t) \approx [ phi(t, s_bar) + dPhi(t, s_bar) (S_t - s_bar) + 1./2. d2Phi(t, s_bar) (S_t - s_bar)^2 ] dW^Q_A
*
******************************************************************************/

	//DL(2) de Phi(t, s)
		double DL2_Phi_t_s(double t, double s) const ;
		double phi(double t, double x_t) const ;
		double dPhi(double t, double x_t) const ;		// d Phi / ds
		double d2Phi(double t, double x_t) const ;	// d2 Phi / ds


/****************  parameter averaging  *********************
*
*	dS(t) \approx [ phi(t, s_bar) + dPhi(t, s_bar) (S_t - s_bar) + 1./2. d2Phi(t, s_bar) (S_t - s_bar)^2 ] dW^Q_A
*	dS(t) = \lambda(t) ( S0 + b(t) (S(t) - S(0)) + 1./2. c(t) (S(t) - S(0))^2 ) dW^Q_A		
*
*	-> time averaging:  b(t) -> b_barre
*						c(t) -> c_barre
************************************************************/

	/* dS(t) = \lambda(t) ( S0 + b_S(t) (S_t - s_bar) + 1./2. c_S(t) (S_t - s_bar)^2 ) dW_t QA */
		double lambda(double t) const ;
		double b_S(double t) const ;
		double c_S(double t) const ;

	//autres fonctions intermédiaires pour les time averaging
		double lambda2(double t) const ;
		double lambda4(double t) const ;
		double f_outer_num_b(double t) const ;
		double f_outer_num_c(double t) const ;

		double f_outer_denom(double t) const ;

	//time averaging of b_S(t) and c_S(t)
		double timeAverage_b_S(double t) const ;
		double timeAverage_c_S(double t) const ;


	double Phi_St(double St) const ;
	double Phi_St_prime(double St) const ;
	double Phi_St_seconde(double St) const ;
	// 1 / Phi(u)
	double UnSurPhi_St(double St) const ;
	// Phi^2(u)
	double Phi2_St(double St) const ;

/**************************  Time change ****************************/

	double timeChange(double t) const ;

/**********************  Call approximation ************************
** dynamique : 
**
**  dS(t) = \Phi(S_t) dW(t)
********************************************************************/

	double omega(double t, double S) const ;
	double omega0(double S) const ;
	double omega1(double S) const ;

	double omega_ATM(double t) const ;  //S = K = strike de la swaption (en attribut)
	double omega0_ATM() const ;
	double omega1_ATM() const ;

	//prix swaption approximé 
	virtual double price(const VanillaSwaption_PTR pSwaption) const ;		

};

typedef boost::shared_ptr<Cheyette_SwaptionPricer_QuadApprox>       Cheyette_SwaptionPricer_QuadApprox_PTR;
typedef boost::shared_ptr<const Cheyette_SwaptionPricer_QuadApprox> Cheyette_SwaptionPricer_QuadApprox_CONSTPTR;