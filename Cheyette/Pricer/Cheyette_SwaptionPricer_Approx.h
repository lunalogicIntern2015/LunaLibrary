#pragma once

#include <boost/shared_ptr.hpp>

#include <Cheyette/Model/CheyetteModel.h>
#include <Cheyette/Model/CheyetteDD_Model.h>
#include <Instrument/VanillaSwaption.h>
#include <Numeric/Integrator1D.h>
#include <Numeric/NumericalMethods.h>	//prix Black
#include <Cheyette/Fonction.h>			//interpolation RR Function

class Cheyette_SwaptionPricer_Approx 
{
protected :
	CheyetteModel_PTR				pCheyette_Model_;  

	//appel fréquent aux éléments suivants -> buffer
	mutable CourbeInput_PTR			buffer_pCourbeInput_ ;
	mutable VanillaSwap				buffer_UnderlyingSwap_ ;
	mutable double					buffer_T0_ ;
	mutable double					buffer_TN_ ;
	mutable std::vector<size_t>		buffer_fixedLegPaymentIndexSchedule_ ;
	mutable std::vector<double>		buffer_deltaTFixedLeg_ ;
	mutable double					buffer_s0_;			

	mutable Interpolation_RR_Function	buffer_y_bar_;

public :
	Cheyette_SwaptionPricer_Approx(const CheyetteModel_PTR& pCheyette_Model)  ;
	virtual ~Cheyette_SwaptionPricer_Approx(){}

//getters
	CheyetteModel_PTR			get_CheyetteModel()	  const {return pCheyette_Model_ ;} 
	CourbeInput_PTR				get_buffer_courbeInput_() const {return buffer_pCourbeInput_ ;}
	VanillaSwap					get_buffer_UnderlyingSwap_() const {return buffer_UnderlyingSwap_ ;}
	double						get_buffer_T0_() const {return buffer_T0_ ;}
	double						get_buffer_TN_() const {return buffer_TN_ ;}
	std::vector<size_t>			get_buffer_fixedLegPaymentIndexSchedule_() const {return buffer_fixedLegPaymentIndexSchedule_ ;}
	std::vector<double>			get_buffer_deltaTFixedLeg_() const {return buffer_deltaTFixedLeg_ ;}
	double						get_buffer_s0_() const {return buffer_s0_ ;}

	Interpolation_RR_Function&	get_buffer_y_bar_() const {return buffer_y_bar_ ;}

	//initialisation des buffers :
	//void preCalculateSwaptionBuffers(VanillaSwaption_PTR pSwaption) const ;
	virtual void preCalculateModelBuffers() const = 0;
	virtual void preCalculateALL(VanillaSwaption_PTR pSwaption) const;

//TODO voir pourquoi il y a un probleme en mettant preCalculateSwaptionBuffers(pSwaption) et preCalculateModelBuffers
	//dans le corps de preCalculateALL(pSwaption)
	
//! YY Bad idea, this should be in the Model, but not the pricer. 
	//pour la calibration, updateVol permet mise à jour des buffers lorsque la vol de CheyetteDD_Model est modifiée
	//modifie la valeur des parametres (fonctions constantes par morceaux) de Cheyette Model à l'index désiré
	void updateLevel_calib		(double newValue, size_t index) ;
	void updateSkew_calib		(double newValue, size_t index) ;	
	void updateConvexity_calib	(double newValue, size_t index) ;	
	
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
	//retourne le x_t tel que S(t, x_t) = swapRate(t, x_t) = s 
		double inverse(double t, double s) const ;

	//fonctions virtuelles
		virtual double phi(double t, double x_t) const = 0 ; //sera lineaire ou quadratique
		//prix swaption approximé 
		virtual double price(const VanillaSwaption_PTR pSwaption) const = 0 ;	

		//volBlack meme code pour Lin et Quad, price de volBlack doit etre virtual
		double volBlack(const VanillaSwaption_PTR vanillaSwaption) const ;

		//pricing de swaption pour des strikes donnés
		//fait appel à price
		std::vector<std::vector<double>> priceMultipleStrikes(VanillaSwaption_PTR pSwaption, std::vector<double> shifts_bp) ;

		//pricing de swaption pour des strikes correspondant à une standardized moneyness dans [-5, 5]
		//std::vector<std::vector<double>> priceMultipleStrikes(double sigma_ATM) ;
};

typedef boost::shared_ptr<Cheyette_SwaptionPricer_Approx>       Cheyette_SwaptionPricer_Approx_PTR;
typedef boost::shared_ptr<const Cheyette_SwaptionPricer_Approx> Cheyette_SwaptionPricer_Approx_CONSTPTR;
