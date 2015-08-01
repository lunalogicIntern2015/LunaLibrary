#include "JBLMM/Pricer/McLmmGenericSwapPricer.h"
#include <cassert>
#include <iostream>
#include <math.h>

#include <JBLMM/Element/Coupon.h>
#include <JBLMM/Element/CappedFlooredCoupon.h>
#include <LMM/helper/LMMTenorStructure.h>
#include <JBLMM/Instrument/GenericSwap.h>
#include <JBLMM/Element/CouponLeg.h>
#include <JBLMM/Element/ConstRate.h>
#include <JBLMM/Element/LiborRate.h>
#include <JBLMM/Element/Rate1.h>


//simulation
double McLmmGenericSwapPricer::swapNPV(const GenericSwap_CONSTPTR genericSwap, size_t nbSimulation)const
{
	double result = 0.0;
	LMM::Index indexValuationDate = 0;
	LMMTenorStructure_PTR lmmTenorStructure=getMcLmm()->get_lmm()->get_LMMTenorStructure();

	//MC
	for(size_t itrSimulation=0; itrSimulation<nbSimulation; ++itrSimulation)
	{
		mcLmm_->simulateLMM();  // YY TODO: not efficient at all, don't need to do all the simulation ...   //
		double npv1  = evaluateCouponLeg(indexValuationDate, genericSwap->getLeg1(), mcLmm_->get_numeraire(), mcLmm_->get_liborMatrix(),lmmTenorStructure);
		double npv2  = evaluateCouponLeg(indexValuationDate, genericSwap->getLeg2(), mcLmm_->get_numeraire(), mcLmm_->get_liborMatrix(),lmmTenorStructure);
		double npvSwap		       =	npv1 - npv2;
		result					  +=	npvSwap;
	}
	result   /=nbSimulation; 

	return result;
}


//! for one couponLeg and one simulation
double McLmmGenericSwapPricer::evaluateCouponLeg(	const LMM::Index indexValuationDate,
													const CouponLeg_CONSTPTR couponLeg,
													const std::vector<double>& numeraire, 
													const matrix& liborMatrix,
													LMMTenorStructure_CONSTPTR lmmTenorStructure)const 
													//suppose all the coupons in couponLeg are the same type.
{
	std::vector<Coupon_CONSTPTR> couponList=couponLeg->getLeg();

	Coupon_CONSTPTR firstCoupon=couponList[0];
	CappedFlooredCoupon_CONSTPTR firstCappedFlooredCoupon = boost::dynamic_pointer_cast<const CappedFlooredCoupon>(firstCoupon);
	if(!firstCappedFlooredCoupon)
		throw("fail to cast cappedFlooredCoupon");

	Rate1_CONSTPTR firstRate = firstCappedFlooredCoupon->getRate();   
	// LiborRate_CONSTPTR	firstLiborRate		=	boost::dynamic_pointer_cast<const LiborRate>(firstRate);
	// ConstRate_CONSTPTR	firstConstRate		=	boost::dynamic_pointer_cast<const ConstRate>(firstRate);

	//if(!firstLiborRate&&!firstConstRate)
	//{
	//	throw("fail to cast LiborRate or ConstRate");
	//}

	double price = 0.0;
	if(boost::dynamic_pointer_cast<const LiborRate>(firstRate))
	{  
		LiborRate_CONSTPTR	firstLiborRate		=	boost::dynamic_pointer_cast<const LiborRate>(firstRate);
		if(firstLiborRate->getDuration()!=lmmTenorStructure->get_tenorType())
			throw("the tenor of the libor demanded is not the same than that in the model");

		for(size_t i=0; i<couponList.size(); i++)
		{
			Coupon_CONSTPTR coupon=couponList[i];		
			LMM::Index paymentIndex = coupon->getPaymentIndex();

			if(paymentIndex<=indexValuationDate)continue;

			CappedFlooredCoupon_CONSTPTR cappedFlooredCoupon = boost::dynamic_pointer_cast<const CappedFlooredCoupon>(coupon);

			Rate1_CONSTPTR rate = cappedFlooredCoupon->getRate();   
			LiborRate_CONSTPTR	liborRate		=	boost::dynamic_pointer_cast<const LiborRate>(rate);

			LMM::Index	indexLibor				=	liborRate->getFixingTime();

			//check if the tenor are the same
			
			double payoffFlow = evaluateCappedFlooredCoupon(cappedFlooredCoupon, liborMatrix(indexLibor,indexLibor));
			double numeraireRatio = numeraire[indexValuationDate]/numeraire[paymentIndex];

			price += payoffFlow*numeraireRatio;
		}
	}
	else if(boost::dynamic_pointer_cast<const ConstRate>(firstRate))
	{
		for(size_t i=0; i<couponList.size(); i++)
		{
			Coupon_CONSTPTR coupon=couponList[i];		
			LMM::Index paymentIndex = coupon->getPaymentIndex();

			if(paymentIndex<=indexValuationDate)continue;

			CappedFlooredCoupon_CONSTPTR cappedFlooredCoupon = boost::dynamic_pointer_cast<const CappedFlooredCoupon>(coupon);

			Rate1_CONSTPTR rate = cappedFlooredCoupon->getRate();   
			ConstRate_CONSTPTR	constRate		=	boost::dynamic_pointer_cast<const ConstRate>(rate);

			double	constRateValue				=	constRate->getConstRateValue();
			price+=numeraire[indexValuationDate]/numeraire[paymentIndex]*evaluateCappedFlooredCoupon(cappedFlooredCoupon, constRateValue);
		}
	}
	else
	{
		throw("fail to cast LiborRate or ConstRate");
	}

	return price;
}

double McLmmGenericSwapPricer::evaluateCappedFlooredCoupon(CappedFlooredCoupon_CONSTPTR cappedFlooredCoupon, double rateValue)const
{
	double	nominal			=	cappedFlooredCoupon->getNominal();
	double	period			=	cappedFlooredCoupon->getPeriod();
	bool	ifFloored		=	cappedFlooredCoupon->getIfFloored();
	double	floor			=	cappedFlooredCoupon->getFloorStrike();
	bool	ifCapped		=	cappedFlooredCoupon->getIfCapped();
	double	cap				=	cappedFlooredCoupon->getCapStrike();
	double	multiFactor		=	cappedFlooredCoupon->getMultiFactor();
	double	addFactor		=	cappedFlooredCoupon->getAddFactor();

	double value		=	multiFactor*rateValue+addFactor;
	if(ifCapped)
		value	=	std::min(cap, value);
	if(ifFloored)
		value	=	std::max(floor, value);

	return nominal*period*value;
}

double McLmmGenericSwapPricer::price_on_oneSimulation(	Instrument_CONSTPTR instrument, 
														LMM::Index indexValuationDate,
														const matrix& liborMatrix, 
														const std::vector<double>& numeraire) const
{
	if(instrument==nullptr)
		return 0.0;

	GenericSwap_CONSTPTR genericSwap = boost::dynamic_pointer_cast<const GenericSwap>(instrument);
	LMMTenorStructure_PTR lmmTenorStructure=getMcLmm()->get_lmm()->get_LMMTenorStructure();
	if(genericSwap)
	{
		double result	= 0.0;
		double npv1		= evaluateCouponLeg(indexValuationDate, genericSwap->getLeg1(), numeraire, liborMatrix,lmmTenorStructure);
		double npv2		= evaluateCouponLeg(indexValuationDate, genericSwap->getLeg2(), numeraire, liborMatrix,lmmTenorStructure);
		result	= npv1 - npv2;
		return result;	
	}
	else
	{
		throw("McLmmGenericTargetSwapPricer do not match the instrument");
	}
}



