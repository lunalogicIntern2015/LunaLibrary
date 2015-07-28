#pragma once

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

	approx prix swaption pour Cheyette Displaced Diffusion

**********************************************************/

class CheyetteDD_VanillaSwaptionApproxPricer  
{

private:
	
	CheyetteDD_Model_PTR	pCheyetteDD_Model_;  
	mutable VanillaSwaption_PTR		pSwaption_ ;   //mutable pour la calibration (skew avec strike +/- shift)

	//appel fréquent aux éléments suivants -> buffer
	mutable CourbeInput_PTR			buffer_pCourbeInput_ ;
	mutable VanillaSwap				buffer_UnderlyingSwap_ ;
	mutable double					buffer_T0_ ;
	mutable double					buffer_TN_ ;
	mutable std::vector<size_t>		buffer_fixedLegPaymentIndexSchedule_ ;
	mutable std::vector<double>		buffer_deltaTFixedLeg_ ;
	mutable double					buffer_s0_;			

	mutable Interpolation_RR_Function	buffer_y_bar_;
	mutable double						buffer_b_barre_ ;

public :
	//constructor  
	CheyetteDD_VanillaSwaptionApproxPricer(	const CheyetteDD_Model_PTR& pCheyetteDD_Model, 
											const VanillaSwaption_PTR&	pSwaption); 

	void initialize_buffers() ;

	//destructor
	virtual ~CheyetteDD_VanillaSwaptionApproxPricer(){};

	//getters
	CheyetteDD_Model_PTR		get_CheyetteDD_Model() const {return pCheyetteDD_Model_ ;}  //CheyetteDD_Model_CONSTPTR
	VanillaSwaption_PTR			get_VanillaSwaption() const {return pSwaption_ ;}

	CourbeInput_PTR				get_buffer_courbeInput_() const {return buffer_pCourbeInput_ ;}
	VanillaSwap					get_buffer_UnderlyingSwap_() const {return buffer_UnderlyingSwap_ ;}
	double						get_buffer_T0_() const {return buffer_T0_ ;}
	double						get_buffer_TN_() const {return buffer_TN_ ;}
	std::vector<size_t>			get_buffer_fixedLegPaymentIndexSchedule_() const {return buffer_fixedLegPaymentIndexSchedule_ ;}
	std::vector<double>			get_buffer_deltaTFixedLeg_() const {return buffer_deltaTFixedLeg_ ;}
	double						get_buffer_s0_() const {return buffer_s0_ ;}

	Interpolation_RR_Function&	get_buffer_y_bar_() const {return buffer_y_bar_ ;}
//	double						get_buffer_y_bar_t(double t) const {return buffer_y_bar_(t) ;}
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

	//pour la calibration, updateVol permet mise à jour des buffers lorsque la vol de CheyetteDD_Model est modifiée
	void updateSigma_calib(std::ostream& o, double a, size_t index)
	{
		o << "sigma : ;" << a ;
		pCheyetteDD_Model_->setCheyetteDD_Parameter_sigma(a, index) ;
		initialize_buffers() ; 
	}

	void updateM_calib(std::ostream& o, double a, size_t index)
	{
		o << "m : ;" << a ;
		pCheyetteDD_Model_->setCheyetteDD_Parameter_m(a, index) ;
		initialize_buffers() ; 
	}
	//calcul de y_barre(t)
	double to_integrate_y_bar(double t) const ;
	void initialize_y_bar(double t, size_t gridSize) const;		
	// ! pendant la calibration, remettre à jour la valeur de y_barre


/**********************  fonctions et dérivées pour ZC, swap rate ************************/

	//dérivée du ZC de maturité T évaluée en t
		double ZC_1stDerivative_on_xt(double t, double T, double x_t, double y_t) const ;
		double ZC_2ndDerivative_on_xt(double t, double T, double x_t, double y_t) const;
		
	//Numerateur
		double swapRateNumerator(double t, double x_t, double y_t) const; 
		//derivee 1ère par rapport à x_t
		double swapRateNumerator_1stDerivative(double t, double x_t, double y_t) const;
		//derivee seconde par rapport à x_t
		double swapRateNumerator_2ndDerivative(double t, double x_t, double y_t) const;


	//Denominateur
		double swapRateDenominator(double t, double x_t, double y_t) const;					
		//derivee 1ère par rapport à x_t
		double swapRateDenominator_1stDerivative(double t, double x_t, double y_t) const;
		//derivee seconde par rapport à x_t
		double swapRateDenominator_2ndDerivative(double t, double x_t, double y_t) const;


	//swapRate = swapNumerator / swapDenominator
		double swapRate0() const;
		double swapRate(double t, double x_t, double y_t) const;
	
		double swapRate_1stDerivative(double t, double x_t, double y_t) const;
		double swapRate0_1stDerivative(double t, double x_t, double y_t) const ;
		
		double swapRate_2ndDerivative(double t, double x_t, double y_t) const;

	//inverse / Newton Raphson
	//S(t, x_t) = swapRate(t, x_t) = s 
	//retourne le x_t correspondant
		double inverse(double t, double s) const ;

/***************** Linear approx of swap rate volatility *********************
*
*	dS(t) = Phi(t, S_t) dW^Q_A
*	dS(t) \approx Phi(t, s_bar) + DPhi(t, s_bar) (S_t - s_bar) dW^Q_A
*
******************************************************************************/
		double Phi(double t, double x_t) const; 
		double DPhi(double t, double x_t) const; 


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
		double prixSwaptionApproxPiterbarg() const ;		//size_t gridSize


};

typedef boost::shared_ptr<CheyetteDD_VanillaSwaptionApproxPricer> CheyetteDD_VanillaSwaptionApproxPricer_PTR;
typedef boost::shared_ptr<const CheyetteDD_VanillaSwaptionApproxPricer> CheyetteDD_VanillaSwaptionApproxPricer_CONSTPTR;

