#include <JBLMM/Pricer/McLmmGenericTargetSwapPricer.h>

#include <JBLMM/Element/Coupon.h>
#include <JBLMM/Element/CappedFlooredCoupon.h>
#include <LMM/helper/LMMTenorStructure.h>
#include <JBLMM/Instrument/GenericSwap.h>
#include <JBLMM/Element/CouponLeg.h>
#include <JBLMM/Element/ConstRate.h>
#include <JBLMM/Element/LiborRate.h>
#include <JBLMM/Element/Rate1.h>
#include <JBLMM/Element/TargetCoupon.h>


double McLmmGenericTargetSwapPricer::swapNPV(GenericSwap_CONSTPTR geneticSwap, size_t nbSimulation)const
{
	LMM::Index indexValuationDate = 0;
	LMMTenorStructure_PTR lmmTenorStructure=getMcLmm()->get_lmm()->get_LMMTenorStructure();

	double result = 0.0;
	for(size_t itrSimulation=0; itrSimulation<nbSimulation; ++itrSimulation)
	{
		mcLmm_->simulateLMM();  // YY TODO: not efficient at all, don't need to do all the simulation ... 
		double npvFloatingLeg_Get  = evaluateCouponLeg(indexValuationDate, geneticSwap->getLeg1(), mcLmm_->get_numeraire(), mcLmm_->get_liborMatrix(),lmmTenorStructure);
		double npvFixedLeg_Give    = evaluateCouponLeg(indexValuationDate, geneticSwap->getLeg2(), mcLmm_->get_numeraire(), mcLmm_->get_liborMatrix(),lmmTenorStructure);
		double npvSwap		       = npvFloatingLeg_Get - npvFixedLeg_Give;
		result				      += npvSwap;
	}
	result   /=nbSimulation; 

	return result;
}

//! for MC
double McLmmGenericTargetSwapPricer::evaluateCouponLeg(	LMM::Index indexValuationDate,
														const CouponLeg_CONSTPTR couponLeg,
														const std::vector<double>& numeraire, 
														const matrix& liborMatrix,
														LMMTenorStructure_PTR lmmTenorStructure)const 
{
	std::vector<Coupon_CONSTPTR> couponList=couponLeg->getLeg();
	
	Coupon_CONSTPTR firstCoupon = couponList[0];
	TargetCoupon_CONSTPTR firstTargetCoupon = boost::dynamic_pointer_cast<const TargetCoupon>(firstCoupon);
	if(!firstTargetCoupon)
		throw("fail to cast targetCoupon");
	double	target	=	firstTargetCoupon->getTarget();
	Rate1_CONSTPTR firstRate = firstTargetCoupon->getRate();   
	LiborRate_CONSTPTR	firstLiborRate		=	boost::dynamic_pointer_cast<const LiborRate>(firstRate);
	ConstRate_CONSTPTR	firstConstRate		=	boost::dynamic_pointer_cast<const ConstRate>(firstRate);
	if(!firstLiborRate&&!firstConstRate)
	{
		throw("fail to cast LiborRate or ConstRate");
	}

	double price = 0.0;
	if(firstLiborRate)
	{
		if(firstLiborRate->getDuration()!=lmmTenorStructure->get_tenorType())
			throw("the tenor of the libor demanded is not the same than that in the model");
		for(size_t i=0; i<couponList.size(); i++)
		{
			Coupon_CONSTPTR coupon=couponList[i];		
			LMM::Index paymentIndex = coupon->getPaymentIndex();

			if(paymentIndex<=indexValuationDate)continue;

			CappedFlooredCoupon_CONSTPTR cappedFlooredCoupon = boost::dynamic_pointer_cast<const CappedFlooredCoupon>(coupon);
			TargetCoupon_CONSTPTR targetCoupon = boost::dynamic_pointer_cast<const TargetCoupon>(coupon);
			Rate1_CONSTPTR rate = targetCoupon->getRate();   
			LiborRate_CONSTPTR	liborRate		=	boost::dynamic_pointer_cast<const LiborRate>(rate);
			LMM::Index	indexLibor				=	liborRate->getFixingTime();

			assert(indexLibor<=paymentIndex);

			if(price<target)
			{	
				double numeraire1		=	numeraire[indexValuationDate];
				double numeraire2		=	numeraire[paymentIndex];
				double numeraireLocal	=	numeraire1/numeraire2;
				double liborValue		=	liborMatrix(indexLibor, indexLibor);
				double result_calculate	=	evaluateCappedFlooredCoupon(cappedFlooredCoupon, liborValue);			
				price				   +=	numeraireLocal*result_calculate;
			}
		}
	}
	else
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

	return price;
}

double McLmmGenericTargetSwapPricer::price_on_oneSimalation(	Instrument_CONSTPTR instrument, 
																LMM::Index indexValuationDate,
																const matrix& liborMatrix, 
																const std::vector<double>& numeraire) const
{
	GenericSwap_CONSTPTR genericTargetSwap = boost::dynamic_pointer_cast<const GenericSwap>(instrument);
	LMMTenorStructure_PTR lmmTenorStructure=getMcLmm()->get_lmm()->get_LMMTenorStructure();
	if(genericTargetSwap)
	{
		double result	= 0.0;
		double npv1		= evaluateCouponLeg(indexValuationDate, genericTargetSwap->getLeg1(), numeraire, liborMatrix,lmmTenorStructure);
		double npv2		= evaluateCouponLeg(indexValuationDate, genericTargetSwap->getLeg2(), numeraire, liborMatrix,lmmTenorStructure);
		result	= npv1 - npv2;
		return result;	
	}
	else
	{
		throw("McLmmGenericTargetSwapPricer do not match the instrument");
	}
}







