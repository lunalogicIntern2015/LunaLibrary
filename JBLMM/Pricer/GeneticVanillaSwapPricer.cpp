#include "JBLMM/Pricer/GeneticVanillaSwapPricer.h"
#include <JBLMM/Element/Coupon.h>
#include <JBLMM/Element/CappedFlooredCoupon.h>
#include <JBLMM/Element/LiborRate.h>

#include <math.h>
#include <memory>
#include <iostream>


double GeneticVanillaSwapPricer::geneticVanillaSwap_Analytical(GeneticSwap_CONSTPTR geneticVanillaSwap, const std::vector<double>& liborsInitValue)const 
{
	std::vector<Coupon_CONSTPTR> floatingLeg	=	geneticVanillaSwap->getLeg1()->getLeg();
	std::vector<Coupon_CONSTPTR> fixedLeg		=	geneticVanillaSwap->getLeg2()->getLeg();

	//floatingLeg
	Coupon_CONSTPTR firstFloatingCoupon=floatingLeg[0];
	CappedFlooredCoupon_CONSTPTR firstFloatingCappedFlooredCoupon = boost::dynamic_pointer_cast<const CappedFlooredCoupon>(firstFloatingCoupon);
	if(!firstFloatingCappedFlooredCoupon)
		throw("fail to cast cappedFlooredCoupon");

	//calculate zero-coupon        //
	double deltaT	=	firstFloatingCappedFlooredCoupon->getPeriod();
	std::vector<double> ZC(liborsInitValue.size()+1);
	ZC[0] = 1.0;
	for(size_t i=1; i<ZC.size(); ++i)
	{
		ZC[i] = ZC[i-1]/(1+deltaT*liborsInitValue[i-1]);
	}

	Rate_CONSTPTR firstFloatingRate = firstFloatingCappedFlooredCoupon->getRate();   
	LiborRate_CONSTPTR	firstLiborRate		=	boost::dynamic_pointer_cast<const LiborRate>(firstFloatingRate);

	double floatingLegValue = 0.0;
	if(firstLiborRate)
	{
		for(size_t i=0; i<floatingLeg.size(); i++)
		{
			Coupon_CONSTPTR coupon=floatingLeg[i];		
			LMM::Index paymentIndex = coupon->getPaymentIndex();

			CappedFlooredCoupon_CONSTPTR cappedFlooredCoupon = boost::dynamic_pointer_cast<const CappedFlooredCoupon>(coupon);

			Rate_CONSTPTR rate = cappedFlooredCoupon->getRate();   
			LiborRate_CONSTPTR	liborRate		=	boost::dynamic_pointer_cast<const LiborRate>(rate);

			LMM::Index	indexLibor				=	liborRate->getFixingTime();
			double nominal						=	cappedFlooredCoupon->getNominal();
			double period						=	cappedFlooredCoupon->getPeriod();
			double floor						=	cappedFlooredCoupon->getFloorStrike();
			double cap							=	cappedFlooredCoupon->getCapStrike();
			double multiFactor					=	cappedFlooredCoupon->getMultiFactor();
			double addFactor					=	cappedFlooredCoupon->getAddFactor();

			floatingLegValue	+=	ZC[paymentIndex]*nominal*period*
									std::max(floor, std::min(cap, multiFactor*liborsInitValue[indexLibor]+addFactor));
		}
	}
	else
	{
		throw("fail to cast liborRate");
	}

	//fixedLeg
	Coupon_CONSTPTR firstFixedCoupon=fixedLeg[0];
	CappedFlooredCoupon_CONSTPTR firstFixedCappedFlooredCoupon = boost::dynamic_pointer_cast<const CappedFlooredCoupon>(firstFixedCoupon);
	if(!firstFixedCappedFlooredCoupon)
		throw("fail to cast cappedFlooredCoupon");

	Rate_CONSTPTR firstFixedRate = firstFixedCappedFlooredCoupon->getRate();   
	ConstRate_CONSTPTR	firstConstRate		=	boost::dynamic_pointer_cast<const ConstRate>(firstFixedRate);

	double fixedLegValue = 0.0;
	if(firstConstRate)
	{
		for(size_t i=0; i<fixedLeg.size(); i++)
		{
			Coupon_CONSTPTR coupon=fixedLeg[i];		
			LMM::Index paymentIndex = coupon->getPaymentIndex();

			CappedFlooredCoupon_CONSTPTR cappedFlooredCoupon = boost::dynamic_pointer_cast<const CappedFlooredCoupon>(coupon);

			Rate_CONSTPTR rate = cappedFlooredCoupon->getRate();   
			ConstRate_CONSTPTR	constRate		=	boost::dynamic_pointer_cast<const ConstRate>(rate);

			double constRateValue				=	constRate->getConstRateValue();
			double nominal						=	cappedFlooredCoupon->getNominal();
			double period						=	cappedFlooredCoupon->getPeriod();
			double floor						=	cappedFlooredCoupon->getFloorStrike();
			double cap							=	cappedFlooredCoupon->getCapStrike();
			double multiFactor					=	cappedFlooredCoupon->getMultiFactor();
			double addFactor					=	cappedFlooredCoupon->getAddFactor();

			fixedLegValue		+=		ZC[paymentIndex]*nominal*period*
										std::max(floor, std::min(cap, multiFactor*constRateValue+addFactor));
		}
	}
	else
	{
		throw("fail to cast ConstRate");
	}

	return floatingLegValue-fixedLegValue;
}
