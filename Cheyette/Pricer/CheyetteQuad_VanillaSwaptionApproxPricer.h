#pragma once

#include <Cheyette/CheyetteModel/CheyetteQuad_Model.h>						
#include <LMM/instrument/VanillaSwaption.h>
#include <LMM/numeric/Integrator1D.h>

#include <LMM/numeric/NumericalMethods.h>  //pour le prix Black

#include <cassert>
#include <vector>
#include <boost/shared_ptr.hpp>
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <Cheyette/Fonction.h>

class CheyetteQuad_VanillaSwaptionApproxPricer
{
private:
	
	CheyetteQuad_Model_PTR			pCheyetteQuad_Model_;  
	mutable VanillaSwaption_PTR		pSwaption_ ;   //mutable pour la calibration (skew avec strike +/- shift)

	//appel fréquent aux éléments suivants -> buffer
	mutable CourbeInput_PTR			pBuffer_courbeInput_ ;
	mutable VanillaSwap				buffer_UnderlyingSwap_ ;
	mutable double					buffer_T0_ ;
	mutable double					buffer_TN_ ;
	mutable std::vector<size_t>		buffer_fixedLegPaymentIndexSchedule_ ;
	mutable std::vector<size_t>		buffer_floatingLegPaymentIndexSchedule_ ;
	mutable std::vector<double>		buffer_deltaTFixedLeg_ ;
	mutable double					buffer_s0_;			

	mutable Interpolation_RR_Function	buffer_y_bar_;
	mutable double						buffer_b_barre_ ;

public:
	//constructor  
	CheyetteQuad_VanillaSwaptionApproxPricer(	const CheyetteQuad_Model_PTR& pCheyetteQuad_Model, 
												const VanillaSwaption_PTR&	pSwaption); 

	void initialize_buffers() ;

	//destructor
	virtual ~CheyetteQuad_VanillaSwaptionApproxPricer(){};

	//getters
	CheyetteQuad_Model_PTR		get_CheyetteQuad_Model() const {return pCheyetteQuad_Model_ ;} 
	VanillaSwaption_PTR			get_VanillaSwaption() const {return pSwaption_ ;}

	CourbeInput_PTR				get_buffer_courbeInput_() const {return pBuffer_courbeInput_ ;}
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
	void setSwaption(double strike, size_t indexStart, size_t indexEnd)
	{
		buffer_UnderlyingSwap_.set_strike(strike) ;
		buffer_UnderlyingSwap_.set_indexStart(indexStart) ; 
		buffer_UnderlyingSwap_.set_indexEnd(indexEnd) ; 

		initialize_buffers() ; 
	}

	void setSwaption(VanillaSwaption_PTR pSwaption)
	{
		pSwaption_ = pSwaption ;
		initialize_buffers() ; 
	}

	void setStrike(double strike)
	{
		pSwaption_->getUnderlyingSwap_RefNonConst().set_strike(strike) ;  // for skew calculation: need to bump strike ...
		buffer_UnderlyingSwap_.set_strike(strike) ;
	}

	//pour la calibration, updateVol permet mise à jour des buffers lorque a(t), b(t) ou c(t) sont modifiés
	void updateA_calib(std::ostream& o, double a, size_t index)
	{
		o << "a : ;" << a ;
		pCheyetteQuad_Model_->setCheyetteQuad_Parameter_a(a, index) ;
		initialize_buffers() ; 
	}

	void updateB_calib(std::ostream& o, double b, size_t index)
	{
		o << "b : ;" << b ;
		pCheyetteQuad_Model_->setCheyetteQuad_Parameter_b(b, index) ;
		initialize_buffers() ; 
	}

	void updateC_calib(std::ostream& o, double c, size_t index)
	{
		o << "c : ;" << c ;
		pCheyetteQuad_Model_->setCheyetteQuad_Parameter_c(c, index) ;
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
		double ZC_1stDerivative_on_xt(double t, double T, double x_t, double y_t) const ;
		double ZC_2ndDerivative_on_xt(double t, double T, double x_t, double y_t) const;
		double ZC_3rdDerivative_on_xt(double t, double T, double x_t, double y_t) const;
		
	//Numerateur (derivee par rapport à x_t)
		double swapRateNumerator(double t, double x_t, double y_t) const; 
		double swapRateNumerator_1stDerivative(double t, double x_t, double y_t) const;
		double swapRateNumerator_2ndDerivative(double t, double x_t, double y_t) const;
		double swapRateNumerator_3rdDerivative(double t, double x_t, double y_t) const;

	//Denominateur (derivee par rapport à x_t)
		double swapRateDenominator(double t, double x_t, double y_t) const;					
		double swapRateDenominator_1stDerivative(double t, double x_t, double y_t) const;
		double swapRateDenominator_2ndDerivative(double t, double x_t, double y_t) const;
		double swapRateDenominator_3rdDerivative(double t, double x_t, double y_t) const;

	//swapRate = swapNumerator / swapDenominator
		double swapRate0() const;
		double swapRate(double t, double x_t, double y_t) const;
	
		double swapRate_1stDerivative(double t, double x_t, double y_t) const;		
		double swapRate_2ndDerivative(double t, double x_t, double y_t) const;
		double swapRate_3rdDerivative(double t, double x_t, double y_t) const;

	//inverse / Newton Raphson
	//S(t, x_t) = swapRate(t, x_t) = s 
	//retourne le x_t correspondant
		double inverse(double t, double s) const ;

	//DL(2) de Phi(t, s)
		double DL2_Phi_t_s(double t, double s) const ;
		double Phi(double t, double x_t) const ;
		double DPhi(double t, double x_t) const ;		// d Phi / ds
		double D2Phi(double t, double x_t) const ;	// d2 Phi / ds


	//
	/* dS(t) = \lambda(t) ( S0 + b_S(t) (S_t - s_bar) + 1./2. c_S(t) (S_t - s_bar)^2 ) dW_t QA */
		double lambda(double t) const ;
		double b_S(double t) const ;
		double c_S(double t) const ;

	//autres fonctions intermédiaires pour les time averaging
		double lambda2(double t) const ;
		double f_outer_num_b(double t) const ;
		double f_outer_num_c(double t) const ;

		double f_outer_denom(double t) const ;

	//time averaging of b_S(t) and c_S(t)
		double timeAverage_b_S(double t) const ;
		double timeAverage_c_S(double t) const ;
};

