#pragma once
#include <Cheyette/Pricer/MC_Cheyette.h>	
#include <Instrument/GenericSwap/GenericSwap.h>

#include <Instrument/Coupon/Coupon.h>
#include <Instrument/Coupon/CouponLeg.h>
#include <Instrument/Coupon/CappedFlooredCoupon.h>
#include <Instrument/Rate/ConstRate.h>
#include <Instrument/Rate/LiborRate.h>

class MC_Cheyette_GenericSwapPricer : public MC_Cheyette
{

public:
	MC_Cheyette_GenericSwapPricer(	CheyetteModel_PTR			cheyetteModel_PTR,
									RNGenerator_PTR				rnGenerator,
									LMMTenorStructure_PTR		pTenorStructure,
									size_t						fwdProbaT,
									size_t						discretizationBetweenDates   )
		:MC_Cheyette(cheyetteModel_PTR, rnGenerator, pTenorStructure, fwdProbaT, discretizationBetweenDates){}

	virtual ~MC_Cheyette_GenericSwapPricer(){}

	//! one simulation - pour les produits generiques derives
	double evaluateCouponLeg(	const LMM::Index valuationIndex,
								const CouponLeg_CONSTPTR couponLeg,
								const Tenor tenorLeg) const;      // !!! tenorLeg peut etre different de TenorStructure !!! 

	double evaluateCappedFlooredCoupon( CappedFlooredCoupon_CONSTPTR cappedFlooredCoupon, 
										double rateValue) const ;

	//! Pricing at time T0=0
	double swapNPV(GenericSwap_CONSTPTR genericSwap, size_t nbSimulation, Tenor tenorLeg1, Tenor tenorLeg2) const ;



};

typedef boost::shared_ptr<MC_Cheyette_GenericSwapPricer>       MC_Cheyette_GenericSwapPricer_PTR;
typedef boost::shared_ptr<const MC_Cheyette_GenericSwapPricer> MC_Cheyette_GenericSwapPricer_CONSTPTR;