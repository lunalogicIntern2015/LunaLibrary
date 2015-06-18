#pragma once
#include <boost/shared_ptr.hpp>
#include <Cheyette/CheyetteModel/CheyetteDD_Model.h>
#include <Cheyette/Pricer/MC_Cheyette.h>	

#include <JBLMM/Instrument/GeneticSwap.h>
#include <JBLMM/Element/Coupon.h>
#include <JBLMM/Element/CappedFlooredCoupon.h>
#include <JBLMM/Element/LiborRate.h>

#include <LMM/numeric/Integrator1D.h>  //pour la fonction closestDate


class MC_CheyetteDD_GenericSwapPricer
{
private:
	MC_Cheyette_PTR mcCheyette_; 

public:
	MC_CheyetteDD_GenericSwapPricer(const MC_Cheyette_PTR& mcCheyette)
		: mcCheyette_(mcCheyette){}

	virtual ~MC_CheyetteDD_GenericSwapPricer(){}

	//gettor
	MC_Cheyette_PTR getMcCheyette()const{return mcCheyette_;}

	//! Pricing at time T0=0
	//double swapRate(const VanillaSwap& vanillaSwap, size_t nbSimulation) const;
	double swapNPV (GeneticSwap_CONSTPTR geneticSwap, size_t nbSimulation) const;

	//! 
	//void resetGeneratorToinitSeed(){mcCheyette_->get_RNGenerator()->resetGeneratorToinitSeed();}

protected: 

	//! one simulation
	virtual double evaluateCouponLeg(	const LMM::Index indexValuationDate,
										const CouponLeg_CONSTPTR couponLeg,
										const std::vector<double>& numeraire, 
										const std::vector<double>& x_t, 
										const std::vector<double>& y_t,
										const Tenor tenor) const;

	double evaluateCappedFlooredCoupon( CappedFlooredCoupon_CONSTPTR cappedFlooredCoupon, 
										double rateValue) const ;

};

typedef boost::shared_ptr<MC_CheyetteDD_GenericSwapPricer>       MC_CheyetteDD_GenericSwapPricer_PTR;
typedef boost::shared_ptr<const MC_CheyetteDD_GenericSwapPricer> MC_CheyetteDD_GenericSwapPricer_CONSTPTR;