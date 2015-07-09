﻿#pragma once

#include <Cheyette\CheyetteModel\CheyetteDD_Model.h>						
#include <LMM/instrument/VanillaSwaption.h>
#include <LMM/numeric/Integrator1D.h>

#include <LMM/numeric/NumericalMethods.h>  //pour le prix Black

#include <cassert>
#include <vector>
#include <boost/shared_ptr.hpp>
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <Cheyette/Fonction.h>

/**********************************************************

	approx pour Cheyette Displaced Diffusion
			approx cf papier de Yuan / Piterbarg

**********************************************************/

class CheyetteDD_VanillaSwaptionApproxPricer  
{

private:
	
	CheyetteDD_Model_PTR	cheyetteDD_Model_;  
	mutable VanillaSwaption_PTR		swaption_ ;   //mutable pour la calibration (skew avec strike +/- shift)

	//appel fréquent aux éléments suivants -> buffer
	mutable CourbeInput_PTR			buffer_courbeInput_ ;
	mutable VanillaSwap				buffer_UnderlyingSwap_ ;
	mutable double					buffer_T0_ ;
	mutable double					buffer_TN_ ;
	mutable std::vector<size_t>		buffer_fixedLegPaymentIndexSchedule_ ;
	mutable std::vector<size_t>		buffer_floatingLegPaymentIndexSchedule_ ;
	mutable std::vector<double>		buffer_deltaTFixedLeg_ ;
	mutable double					buffer_s0_;		

	mutable Interpolation_RR_Function	buffer_y_bar_;
	mutable double						buffer_b_barre_ ;

public :
	//constructor  
	CheyetteDD_VanillaSwaptionApproxPricer(	const CheyetteDD_Model_PTR& cheyetteDD_Model, 
											const VanillaSwaption_PTR&	swaption); 

	void initialize_buffers() ;

	//destructor
	virtual ~CheyetteDD_VanillaSwaptionApproxPricer(){};

	//getters
	CheyetteDD_Model_PTR	get_CheyetteDD_Model() const {return cheyetteDD_Model_ ;}  //CheyetteDD_Model_CONSTPTR
	VanillaSwaption_PTR			get_VanillaSwaption() const {return swaption_ ;}

	CourbeInput_PTR				get_buffer_courbeInput_() const {return buffer_courbeInput_ ;}
	VanillaSwap					get_buffer_UnderlyingSwap_() const {return buffer_UnderlyingSwap_ ;}
	double						get_buffer_T0_() const {return buffer_T0_ ;}
	double						get_buffer_TN_() const {return buffer_TN_ ;}
	std::vector<size_t>			get_buffer_fixedLegPaymentIndexSchedule_() const {return buffer_fixedLegPaymentIndexSchedule_ ;}
	std::vector<size_t>			get_buffer_floatingLegPaymentIndexSchedule_() const {return buffer_floatingLegPaymentIndexSchedule_ ;}
	std::vector<double>			get_buffer_deltaTFixedLeg_() const {return buffer_deltaTFixedLeg_ ;}
	double						get_buffer_s0_() const {return buffer_s0_ ;}

	Interpolation_RR_Function&	get_buffer_y_bar_() const {return buffer_y_bar_ ;}
	double						get_buffer_y_bar_t(double t) const {return buffer_y_bar_(t) ;}
	double						get_buffer_b_barre_() const {return buffer_b_barre_ ;}

	//setters 
	void setSwaption(const VanillaSwaption_PTR swaption)
	{
		swaption_ = swaption ;
		std::cout << "changement de la swaption et recalcul des buffers de ApproxPricer" << std::endl ;
		initialize_buffers() ; 
	}

	void setStrike(double strike)
	{
		swaption_->getUnderlyingSwap_RefNonConst().set_strike(strike) ;  // for skew calculation: need to bump strike ...
		buffer_UnderlyingSwap_.set_strike(strike) ;
	}

	//pour la calibration, updateVol permet mise à jour des buffers lorsque la vol de CheyetteDD_Model est modifiée
	void updateSigma_calib(const Array& A)
	{
		cheyetteDD_Model_->setCheyetteDD_Parameter_sigma(A) ;
		std::cout << "changement de la swaption et recalcul des buffers de ApproxPricer" << std::endl ;
		initialize_buffers() ; 
	}

	void updateM_calib(const Array& A)
	{
		cheyetteDD_Model_->setCheyetteDD_Parameter_m(A) ;
		std::cout << "changement de la swaption et recalcul des buffers de ApproxPricer" << std::endl ;
		initialize_buffers() ; 
	}
	//calcul de y_barre(t)
	double to_integrate_y_bar(double t) const ;
	void initialize_y_bar(double t, size_t gridSize) const;		
	// ! pendant la calibration, remettre à jour la valeur de y_barre



/**********************  fonctions et dérivées pour ZC, swap rate ************************
** normalement P(t, T, x_t, y_t), A_{0,N}(t, x_t, y_t, swap)
** y_t -> \bar{y}_t : constante
**
** 2 cas f(t, T, x_t) :
**      - si t = 0, x_t = 0 : utilisation de la courbe spot
**      - si t > 0, x_t paramètre pour la fonction inverse 
******************************************************************************************/

	//dérivée du ZC de maturité T évaluée en t
		double ZC_1stDerivative_on_xt(double t, double T, double x_t) const ;
		double ZC_2ndDerivative_on_xt(double t, double T, double x_t) const;

		
	//Numerateur
		double swapRateNumerator(double t, double x_t) const; 
		//derivee 1ère par rapport à x_t
		double swapRateNumerator_1stDerivative(double t, double x_t) const;
		//derivee seconde par rapport à x_t
		double swapRateNumerator_2ndDerivative(double t, double x_t) const;


	//Denominateur
		double swapRateDenominator(double t, double x_t) const;					
		//derivee 1ère par rapport à x_t
		double swapRateDenominator_1stDerivative(double t, double x_t) const;
		//derivee seconde par rapport à x_t
		double swapRateDenominator_2ndDerivative(double t, double x_t) const;


	//swapRate = swapNumerator / swapDenominator
		double swapRate0() const;
		double swapRate(double t, double x_t) const;
		double swapRate_1stDerivative(double t, double x_t) const;
		double swapRate_2ndDerivative(double t, double x_t) const;

	//inverse / Newton Raphson
	//S(t, x_t) = swapRate(t, x_t) = s 
	//retourne le x_t correspondant
		double inverse(double t, double s) const ;

	//SwapRateVolatility \Phi(t, S_t)
	//dS(t) = \Phi(t, S_t) dW_t^{Q Annuite}

		// d\Phi /ds (t, \bar{s}) : vol de S(t, x_t)
		double swapRateVolatility_1stDerivative(double t, double x_t) const;

		//evaluee en t et pour un taux de swap s_bar = s0_
		double calculate_phi_t_s_bar(double t) const; 

		//approximation de Phi(t, s) en s
		//DL autour de s_bar
		//double swapRateVolatility_approx_lineaire(double t, double s) const;


/****************  parameter averaging  *********************
*
*	dS(t) = \lambda(t) ( (1-b(t)) S0 + b(t) S(t) ) dW^Q_A		-> time averaging: b(t) puis b_barre
*
*	dS(t) = (A(t) + B(t) S(t)) dW^Q_A
*	S(0) = S0_ connu
************************************************************/
		double A(double t) const;
		double B(double t) const;

		double b(double t) const;
		double lambda(double t) const;
		double lambda2(double t) const;
		//double v2(double t) const;  //integrale de 0 à t de lambda^2(u) du

		double f_outer_num(double t) const ;		//lambda^2(u) * b(u)
		double f_outer_denom(double t) const ;		//lambda^2(u)
		
		//retourne b barre du displaced diffusion
		double timeAverage(double t) const ;				//, size_t gridSize	

		//prix swaption approximé 
		double prixSwaptionApproxPiterbarg() const ;		//size_t gridSize


};

typedef boost::shared_ptr<CheyetteDD_VanillaSwaptionApproxPricer> CheyetteDD_VanillaSwaptionApproxPricer_PTR;
typedef boost::shared_ptr<const CheyetteDD_VanillaSwaptionApproxPricer> CheyetteDD_VanillaSwaptionApproxPricer_CONSTPTR;

