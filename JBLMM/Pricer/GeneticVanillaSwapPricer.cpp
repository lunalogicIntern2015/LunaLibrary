#include "JBLMM/Pricer/GeneticVanillaSwapPricer.h"
#include <JBLMM/Element/Coupon.h>
#include <JBLMM/Element/CappedFlooredCoupon.h>
#include <JBLMM/Element/LiborRate.h>

#include <math.h>
#include <memory>
#include <iostream>


//size_t GeneticVanillaSwapPricer::find(double paymentDay, std::vector<double> dateLibor)const
//{
//	assert(paymentDay>=0);
//	for(size_t i=0; i<dateLibor.size();i++)
//	{
//		if(paymentDay<dateLibor[i]){return i-1;}
//	}
//	return dateLibor.size()-1;
//}

double GeneticVanillaSwapPricer::geneticVanillaSwap_Analytical(GeneticSwap_CONSTPTR geneticVanillaSwap, const std::vector<double>& liborsInitValue)const 
{
	std::vector<Coupon_CONSTPTR> floatingLeg	=	geneticVanillaSwap->getLeg1()->getLeg();
	std::vector<Coupon_CONSTPTR> fixedLeg		=	geneticVanillaSwap->getLeg2()->getLeg();

	LMM::Index liborHorizon			=	liborsInitValue.size();
	LMM::Index floatingLegHorizon	=	floatingLeg.size();
	LMM::Index fixedLegHorizon		=	fixedLeg.size();


	//cast from Coupon_CONSTPTR to CappedFlooredCoupon_PTR
	Coupon_PTR c		=	boost::const_pointer_cast<Coupon>(floatingLeg[0]);
	CappedFlooredCoupon_PTR cfc = boost::dynamic_pointer_cast<CappedFlooredCoupon>(c);

	//get deltaT from the first floatingcoupon 
	const double deltaT				=	cfc->getPeriod();

	//calculate zero-coupon
	std::vector<double> ZC(liborHorizon+1);
	ZC[0] = 1.0;
	for(size_t i=1; i<ZC.size(); ++i)
	{
		ZC[i] = ZC[i-1]/(1+deltaT*liborsInitValue[i-1]);
	}


	//calculate the floatingLeg
	Coupon_CONSTPTR firstfloatingCoupon=floatingLeg[0];
	CappedFlooredCoupon_CONSTPTR firstFloatingCappedFlooredCoupon = boost::dynamic_pointer_cast<const CappedFlooredCoupon>(firstfloatingCoupon);
	if(!firstFloatingCappedFlooredCoupon)
		throw("fail to cast CappedFlooredCoupon");

	Rate_CONSTPTR firstFloatingLegRate = firstFloatingCappedFlooredCoupon->getRate();   
	LiborRate_CONSTPTR	firstLiborRate		=	boost::dynamic_pointer_cast<const LiborRate>(firstFloatingLegRate);
	if(!firstLiborRate)
		throw("fail to cast rate to liborRate");

	double floatingIncome=0.0;
	for(size_t i=0; i<floatingLegHorizon;i++)
	{
		Coupon_CONSTPTR coupon	= floatingLeg[i];
		LMM::Index paymentIndex = coupon->getPaymentIndex();
		CappedFlooredCoupon_CONSTPTR cappedFlooredCoupon = boost::dynamic_pointer_cast<const CappedFlooredCoupon>(coupon);

		Rate_CONSTPTR	rate	=	cappedFlooredCoupon->getRate();
		LiborRate_CONSTPTR liborRate = boost::dynamic_pointer_cast<const LiborRate>(rate);		

		LMM::Index indexLibor = liborRate->getFixingTime();
		double tauxActualisation = ZC[paymentIndex];	
		floatingIncome	+=	tauxActualisation*
							cappedFlooredCoupon->getNominal()*
							cappedFlooredCoupon->getPeriod()*
							std::max(	cappedFlooredCoupon->getFloorStrike(), 
										std::min(	cappedFlooredCoupon->getCapStrike(),
													cappedFlooredCoupon->getMultiFactor()*liborsInitValue[indexLibor]
													+cappedFlooredCoupon->getAddFactor()));
	}


	// Caculate the fixedleg
	Coupon_CONSTPTR firstfixedCoupon=fixedLeg[0];
	CappedFlooredCoupon_CONSTPTR firstFixdCappedFlooredCoupon = boost::dynamic_pointer_cast<const CappedFlooredCoupon>(firstfixedCoupon);
	if(!firstFixdCappedFlooredCoupon)
		throw("fail to cast Coupon to CappedFlooredCoupon");

	Rate_CONSTPTR firstFixedLegRate = firstFixdCappedFlooredCoupon->getRate();   
	ConstRate_CONSTPTR	firstConstRate		=	boost::dynamic_pointer_cast<const ConstRate>(firstFixedLegRate);
	if(!firstConstRate)
		throw("fail to cast rate to ConstRate");
	
	double fixedIncome=0.0;
	for(size_t i=0; i<fixedLegHorizon;i++)
	{
		Coupon_CONSTPTR coupon=fixedLeg[i];		
		LMM::Index paymentIndex = coupon->getPaymentIndex();

		CappedFlooredCoupon_CONSTPTR cappedFlooredCoupon = boost::dynamic_pointer_cast<const CappedFlooredCoupon>(coupon);

		Rate_CONSTPTR rate = cappedFlooredCoupon->getRate();   
		ConstRate_CONSTPTR	constRate		=	boost::dynamic_pointer_cast<const ConstRate>(rate);

		double	constRateValue				=	constRate->getConstRateValue();
		double tauxActualisation = ZC[paymentIndex];
		fixedIncome	+=	tauxActualisation
						*cappedFlooredCoupon->getNominal()
						*cappedFlooredCoupon->getPeriod()*
						std::max(	cappedFlooredCoupon->getFloorStrike(), 
									std::min(	cappedFlooredCoupon->getCapStrike(),
																		cappedFlooredCoupon->getMultiFactor()
																		*constRateValue
																		+cappedFlooredCoupon->getAddFactor()));
	}

	return floatingIncome-fixedIncome;

}
