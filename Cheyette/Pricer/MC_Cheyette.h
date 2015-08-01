#pragma once

#include <iostream>
#include <cassert>
#include <string>
#include <boost/shared_ptr.hpp>

#include <Cheyette/CheyetteModel/CheyetteDD_Model.h>
#include <LMM/LmmModel/McLmm.h>  //pour le namespace MCSchemeType
#include <LMM/RNGenerator/RNGenerator.h>
#include <LMM/RNGenerator/McGenerator.h>
#include <LMM/helper/TenorType.h>
#include <LMM/helper/LMMTenorStructure.h>

#include <JBLMM/Element/Coupon.h>
#include <JBLMM/Element/CouponLeg.h>
#include <JBLMM/Element/CappedFlooredCoupon.h>
#include <JBLMM/Element/LiborRate.h>


/**************  generation de x(t) et y(t) sous la proba forward Q_T  **************
*
*	dx(t) = [ y(t) - k x(t) - G(t, T) sigma_r(t)^2 ] dt + sigma_r(t) dW(t)^QT
*	dy(t) = ( sigma_r(t)^2 - 2 k y(t) ) dt
*
**************************************************************************************/

class MC_Cheyette
{
protected:
	CheyetteDD_Model_PTR		cheyetteDD_Model_ ;
	RNGenerator_PTR				rnGenerator_;

// UNE UNIQUE TENOR STRUCTURE pour tous les flux (même si 3M / 6M) //
	//LMMTenorStructure(const Tenor&  tenorType, int max_nbOfYear);
	LMMTenorStructure_PTR		pTenorStructure_ ;	

	size_t						fwdProbaT_ ;					// define the numeraire

	size_t						discretizationBetweenDates_ ;	//nb de pas de discretisation entre 2 dates

	//! comes from simulation
	mutable std::vector<double> x_t_Cheyette_ ;
	mutable std::vector<double> y_t_Cheyette_ ;
	mutable std::vector<double> numeraires_; // numeraire[i] = P(T_i, fwd_probaT)  , size = N

public:

	MC_Cheyette(	CheyetteDD_Model_PTR		cheyetteDD_Model,
					RNGenerator_PTR				rnGenerator,
					LMMTenorStructure_PTR		pTenorStructure,
					size_t						fwdProbaT,
					size_t						discretizationBetweenDates   )
		:cheyetteDD_Model_(cheyetteDD_Model), rnGenerator_(rnGenerator), pTenorStructure_(pTenorStructure),
		fwdProbaT_(fwdProbaT), discretizationBetweenDates_(discretizationBetweenDates), 
		x_t_Cheyette_(fwdProbaT / pTenorStructure->get_tenorType().YearFraction() + 1), 
		y_t_Cheyette_(fwdProbaT / pTenorStructure->get_tenorType().YearFraction() + 1), 
		numeraires_(fwdProbaT / pTenorStructure->get_tenorType().YearFraction() + 1) 
	{
			assert(pTenorStructure->get_max_nbOfYear() >= fwdProbaT)	;
	} 

	virtual ~MC_Cheyette(){} 
	
	//-- Getters 
	CheyetteDD_Model_PTR		getCheyetteDD_Model() const{return cheyetteDD_Model_ ;}
	RNGenerator_PTR				getRNGenerator() const{return rnGenerator_ ;}
	LMMTenorStructure_PTR		getTenorStructure() const {return pTenorStructure_ ;}
	size_t						getFwdProbaT() const{return fwdProbaT_ ;}
	size_t						getDiscretizationBetweenDates() const{return discretizationBetweenDates_ ;}

	std::vector<double>			get_x_t_Cheyette() const{return x_t_Cheyette_ ;}
	std::vector<double>			get_y_t_Cheyette() const{return y_t_Cheyette_ ;}
	std::vector<double>			getNumeraires() const{return numeraires_ ;}

	//remplit les vecteurs xt et yt aux dates de la tenor structure sous la mesure Q^T (T = fwdProbaT_)
	void simulate_Euler() const;
	void computeNumeraires() const ;

	//! one simulation - pour les produits vanille derives
	double evaluateFloatLeg(	const size_t valuationIndex,
								const std::vector<size_t>& indexFloatLeg,
								const Tenor tenorFloatLeg) const;

	double evaluateFixedLeg(	const size_t valuationIndex,
								const std::vector<size_t>& indexFixedLeg,
								const Tenor tenorFixedLeg, 
								const double fixedRate) const;

	//! one simulation - pour les produits generiques derives
	//double evaluateCouponLeg(	const LMM::Index valuationIndex,
	//							const CouponLeg_CONSTPTR couponLeg,
	//							const Tenor tenorLeg) const;      // !!! tenorLeg peut etre different de TenorStructure !!! 

	//double evaluateCappedFlooredCoupon( CappedFlooredCoupon_CONSTPTR cappedFlooredCoupon, 
	//									double rateValue) const ;
};

typedef boost::shared_ptr<MC_Cheyette> MC_Cheyette_PTR;
typedef boost::shared_ptr<const MC_Cheyette> MC_Cheyette_CONSTPTR;
