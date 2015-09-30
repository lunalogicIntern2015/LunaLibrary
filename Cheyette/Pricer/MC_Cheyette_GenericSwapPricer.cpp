#include "MC_CheyetteDD_GenericSwapPricer.h"

//! one simulation - pour les produits generiques derives
//! for one couponLeg and one simulation
double MC_Cheyette_GenericSwapPricer::evaluateCouponLeg(	const size_t valuationIndex,
															const CouponLeg_CONSTPTR couponLeg,
															const Tenor tenorLeg)  const       
						// !!! tenorLeg peut etre different de TenorStructure !!! 
								// suppose all the coupons in couponLeg are the same type.  
{
	std::vector<Coupon_CONSTPTR> couponList=couponLeg->getLeg();

	Coupon_CONSTPTR firstCoupon=couponList[0];
	CappedFlooredCoupon_CONSTPTR firstCappedFlooredCoupon = boost::dynamic_pointer_cast<const CappedFlooredCoupon>(firstCoupon);
	if(!firstCappedFlooredCoupon){
		std::cout << "fail to cast cappedFlooredCoupon" << std::endl ;
		throw("fail to cast cappedFlooredCoupon");
	}
	Rate1_CONSTPTR firstRate = firstCappedFlooredCoupon->getRate();   

	double price = 0.0;
	if(boost::dynamic_pointer_cast<const LiborRate>(firstRate))
	{  
		LiborRate_CONSTPTR	firstLiborRate		=	boost::dynamic_pointer_cast<const LiborRate>(firstRate);

		size_t posIndexValuationDate	= findIndex(indexValuationDate, indexOfSimulation_) ;
		for(size_t i=0; i<couponList.size(); ++i)
		{
			Coupon_CONSTPTR coupon=couponList[i];		
			size_t paymentIndex = coupon->getPaymentIndex();

			if(paymentIndex<=indexValuationDate)continue;

			CappedFlooredCoupon_CONSTPTR cappedFlooredCoupon = boost::dynamic_pointer_cast<const CappedFlooredCoupon>(coupon);

			Rate1_CONSTPTR rate = cappedFlooredCoupon->getRate();   
			LiborRate_CONSTPTR	liborRate		=	boost::dynamic_pointer_cast<const LiborRate>(rate);

			double tenorYearFrac	= tenorLeg.YearFraction() ;
			size_t fixingIndex		= liborRate->getFixingTime();
			double fixingDate		= fixingIndex * tenorYearFrac ;
			double paymentDate		= paymentIndex * tenorYearFrac ;
			
			size_t pos = findIndex(fixingIndex, indexOfSimulation_) ;

			double libor = cheyetteModel_PTR_->libor(fixingDate, fixingDate, paymentDate, x_t[pos], y_t[pos]) ; 
			
			//retourne nominal * delta * min(max(value, ...)...) 
			double payoffFlow		= evaluateCappedFlooredCoupon(cappedFlooredCoupon, libor) ;
																		
			size_t posPaymentIndex	= findIndex(paymentIndex, indexOfSimulation_) ;
			
			double numeraireRatio	= numeraire[posIndexValuationDate]/numeraire[posPaymentIndex];

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

			Rate1_CONSTPTR rate				= cappedFlooredCoupon->getRate();   
			ConstRate_CONSTPTR	constRate	= boost::dynamic_pointer_cast<const ConstRate>(rate);
			double	constRateValue			= constRate->getConstRateValue();

			size_t posIndexValuationDate	= findIndex(indexValuationDate, indexOfSimulation_) ;
			size_t posPaymentIndex			= findIndex(paymentIndex, indexOfSimulation_) ;

			price+=numeraire[posIndexValuationDate]/numeraire[posPaymentIndex]
										* evaluateCappedFlooredCoupon(cappedFlooredCoupon, constRateValue);
		}
	}
	else
	{
		throw("fail to cast LiborRate or ConstRate");
	}

	return price;
}

double MC_Cheyette_GenericSwapPricer::evaluateCappedFlooredCoupon(CappedFlooredCoupon_CONSTPTR cappedFlooredCoupon, 
														   double rateValue) const
{
	double	nominal			=	cappedFlooredCoupon->getNominal();
	double	delta			=	cappedFlooredCoupon->getDelta();
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

	return nominal* delta *value;
}


//simulation
double MC_Cheyette_GenericSwapPricer::swapNPV(GenericSwap_CONSTPTR genericSwap, size_t nbSimulation, 
												Tenor tenorLeg1, Tenor tenorLeg2) const
{
	double TEST_sommeLeg1 = 0. ;
	double TEST_sommeLeg2 = 0. ;

	double result	= 0. ;
	size_t indexValuationDate = 0 ;

	//MC
	for(size_t itrSimulation=0; itrSimulation<nbSimulation; ++itrSimulation)
	{
		if ((itrSimulation*10) % nbSimulation == 0){std::cout << double(itrSimulation)/double(nbSimulation)*100 << "%" << std::endl ;}
		//simulate_Euler() ;
		//double npv1  = evaluateCouponLeg(indexValuationDate, genericSwap->getLeg1(), tenorLeg1);

		//double npv2  = evaluateCouponLeg(indexValuationDate, genericSwap->getLeg2(), tenorLeg2);			
		//TEST_sommeLeg1 += npv1 ;
		//TEST_sommeLeg2 += npv2 ;
		//result += npv1 - npv2;
	}
	double TEST_meanLeg1 = TEST_sommeLeg1 / nbSimulation ;
	double TEST_meanLeg2 = TEST_sommeLeg2 / nbSimulation ;
	result   /=nbSimulation; 

	std::cout   << "prix MC swap : " << result << std::endl;
	std::cout   << "prix leg1 : " << TEST_meanLeg1 << std::endl;
	std::cout   << "prix leg2 : " << TEST_meanLeg2 << std::endl;
	
	return result;
}


