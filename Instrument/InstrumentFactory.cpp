#include <Instrument/InstrumentFactory.h>
#include <Instrument/GenericSwap/GenericSwap.h>
#include <vector>
#include <iostream>
#include <cassert>
#include <Instrument/Coupon/TargetCoupon.h>


GenericSwap_CONSTPTR InstrumentFactory::createVanillaSwap(	 double strike, 
															 LMM::Index indexStart, 
												 			 LMM::Index indexEnd, 
															 Tenor leg1Tenor,						//
															 Tenor leg2Tenor,						
															 LMMTenorStructure_CONSTPTR swapStructure,
															 double nominal)
	{
		size_t leg1VsLiborTenorTypeRatio	=	leg1Tenor.ratioTo(swapStructure->get_tenorType());		//lmmTenorStructure is a index base.
		size_t leg2VsLiborTenorTypeRatio    =	leg2Tenor.ratioTo(swapStructure->get_tenorType());

		//! check 
		assert(  indexEnd > indexStart ); assert( indexStart >=0 );		
		assert( leg1VsLiborTenorTypeRatio == 1  );	//TODO: necessary for zero-coupon in GeneticVanillaSwap
		assert( (indexEnd - indexStart)%leg2VsLiborTenorTypeRatio ==0  );	

		size_t nbFloatLeg = (indexEnd - indexStart)/leg1VsLiborTenorTypeRatio;
		size_t nbFixedLeg =	(indexEnd - indexStart)/leg2VsLiborTenorTypeRatio;

		//to get all the payment index
		std::vector<size_t> floatingLegPaymentIndexSchedule;
		std::vector<size_t> fixedLegPaymentIndexSchedule;

		for(size_t i=0; i<nbFloatLeg ; ++i)
		{
			floatingLegPaymentIndexSchedule.push_back(indexStart+(i+1)*leg1VsLiborTenorTypeRatio);
		}

		for(size_t i=0; i< nbFixedLeg ; ++i)
		{
			fixedLegPaymentIndexSchedule.push_back(indexStart+(i+1)*leg2VsLiborTenorTypeRatio);
		}

		bool   floatingLeg_ifFloored	= false;
		double floatingLeg_floorStrike	= -1.0e100;
		bool   floatingLeg_ifCapped		= false;
		double floatingLeg_capStrike	= 1.0e100;
		double floatingLeg_multFactor	= 1.0;
		double floatingLeg_addFactor	= 0.0;

		bool   fixedLeg_ifFloored		= false;
		double fixedLeg_floorStrike		= -1.0e100;
		bool   fixedLeg_ifCapped		= false;
		double fixedLeg_capStrike		= 1.0e100;
		double fixedLeg_multFactor		= 1.0;
		double fixedLeg_addFactor		= 0.0;

		LMM::Index  valuationDateIndex = 0;

		//Construction of couponLegs.
		std::vector<Coupon_CONSTPTR> floatingCouponVector;
		std::vector<Coupon_CONSTPTR> fixedCouponVector;

		//floatingLeg
		for(size_t i=0;i<floatingLegPaymentIndexSchedule.size(); i++)	//
		{
			LMM::Index liborFixingindex = floatingLegPaymentIndexSchedule[i]-leg1VsLiborTenorTypeRatio; // TODO: verifier !!!
			Rate1_CONSTPTR liborRate(new LiborRate(liborFixingindex, leg1Tenor));
			LMM::Index couponPaymentIndex=floatingLegPaymentIndexSchedule[i];
			floatingCouponVector.push_back(Coupon_CONSTPTR(new CappedFlooredCoupon(	couponPaymentIndex, 
																					nominal, 
																					leg1Tenor.YearFraction(), 
																					floatingLeg_ifFloored,
																					floatingLeg_floorStrike,
																					floatingLeg_ifCapped,
																					floatingLeg_capStrike,
																					liborRate,
																					floatingLeg_multFactor,
																					floatingLeg_addFactor,
																					valuationDateIndex)));
		}


		for(size_t i=0;i<fixedLegPaymentIndexSchedule.size(); i++)
		{
			LMM::Index couponPaymentIndex=fixedLegPaymentIndexSchedule[i];
			Rate1_CONSTPTR constRate(new ConstRate(strike));
			fixedCouponVector.push_back(Coupon_PTR(new CappedFlooredCoupon(	couponPaymentIndex, 
																			nominal, 
																			leg2Tenor.YearFraction(), 
																			fixedLeg_ifFloored,
																			fixedLeg_floorStrike, 
																			fixedLeg_ifCapped,
																			fixedLeg_capStrike, 
																			constRate, 
																			fixedLeg_multFactor, 
																			fixedLeg_addFactor, 
																			valuationDateIndex)));
		}
		GenericSwap_CONSTPTR vanillaSwap(	new GenericSwap(	CouponLeg_CONSTPTR(new CouponLeg(floatingCouponVector)),
																CouponLeg_CONSTPTR(new CouponLeg(fixedCouponVector))));
		return vanillaSwap;
	}


GenericSwap_CONSTPTR InstrumentFactory::createStandardTARNSwap(	double strike, 
															LMM::Index indexStart, 
												 			LMM::Index indexEnd, 
															Tenor leg1Tenor,		//
															Tenor leg2Tenor,		
															LMMTenorStructure_CONSTPTR swapStructure,
															double nominal,
															double target)
{
	size_t leg1VsLiborTenorTypeRatio = leg1Tenor.ratioTo(swapStructure->get_tenorType()) ;
	size_t leg2VsLiborTenorTypeRatio    = leg2Tenor.ratioTo(swapStructure->get_tenorType())   ;

	//! check 
	assert( indexEnd > indexStart ); assert( indexStart >=0 );
	assert( (indexEnd - indexStart)%leg1VsLiborTenorTypeRatio ==0  );
	assert( (indexEnd - indexStart)%leg2VsLiborTenorTypeRatio ==0  );

	size_t nbLeg1 = (indexEnd - indexStart)/leg1VsLiborTenorTypeRatio;
	size_t nbLeg2 =	(indexEnd - indexStart)/leg2VsLiborTenorTypeRatio;

	//!get all the paymentIndex
	std::vector<size_t> leg1PaymentIndexSchedule;
	std::vector<size_t> leg2PaymentIndexSchedule;

	for(size_t i=0; i<nbLeg1 ; ++i)
	{
		leg1PaymentIndexSchedule.push_back(indexStart+(i+1)*leg1VsLiborTenorTypeRatio);
	}
	for(size_t i=0; i<nbLeg2; ++i)
	{
		leg2PaymentIndexSchedule.push_back(indexStart+(i+1)*leg2VsLiborTenorTypeRatio);
	}

	bool   leg1_ifFloored	= false;
	double leg1_floorStrike	= -1.0e100;
	bool   leg1_ifCapped	= false;
	double leg1_capStrike	= 1.0e100;
	double leg1_multFactor	= 1.0;
	double leg1_addFactor	= 0.0;

	bool   leg2_ifFloored   = false;
	double leg2_floorStrike = -1.0e100;
	bool   leg2_ifCapped    = false;
	double leg2_capStrike   = 1.0e100;
	double leg2_multFactor	= 1.0;
	double leg2_addFactor	= 0.0;

	LMM::Index  valuationDateIndex = 0 ; 
	std::string leg1CouponDependency = "";
	std::string leg2CouponDependency = "";


	//build couponLeg
	std::vector<Coupon_CONSTPTR> leg1CouponVector;
	std::vector<Coupon_CONSTPTR> leg2CouponVector;

	//Leg1
	for(size_t i=0;i<leg1PaymentIndexSchedule.size(); i++)
	{
		
		LMM::Index couponPaymentIndex	=	leg1PaymentIndexSchedule[i];
		LMM::Index indexLibor			=	couponPaymentIndex-leg1VsLiborTenorTypeRatio;
		Rate1_PTR liborRate(new LiborRate(indexLibor,leg1Tenor));
		leg1CouponVector.push_back(Coupon_CONSTPTR(new TargetCoupon(	couponPaymentIndex, 
																		nominal, 
																		leg1Tenor.YearFraction(), 
																		leg1_ifFloored,
																		leg1_floorStrike,
																		leg1_ifCapped,
																		leg1_capStrike,
																		liborRate,
																		leg1_multFactor,
																		leg1_addFactor,
																		valuationDateIndex,
																		target,
																		leg1CouponDependency)));
		leg1CouponDependency += "l"; 
	}


	for(size_t i=0;i<leg2PaymentIndexSchedule.size(); i++)
	{
		LMM::Index couponPaymentIndex=leg2PaymentIndexSchedule[i];
		Rate1_PTR constRate(new ConstRate(strike));
		leg2CouponVector.push_back(Coupon_CONSTPTR(new TargetCoupon(	couponPaymentIndex, 
																		nominal, 
																		leg2Tenor.YearFraction(), 
																		leg2_ifFloored,
																		leg2_floorStrike,
																		leg2_ifCapped,
																		leg2_capStrike,
																		constRate,
																		leg2_multFactor,
																		leg2_addFactor,
																		valuationDateIndex,
																		target,
																		leg2CouponDependency)));
		leg2CouponDependency += "l"; 
	}
	GenericSwap_CONSTPTR targetSwap(	new GenericSwap(	CouponLeg_CONSTPTR(new CouponLeg(leg1CouponVector)),
															CouponLeg_CONSTPTR(new CouponLeg(leg2CouponVector))));
	return targetSwap;
}


GenericSwap_CONSTPTR InstrumentFactory::createGenericSwap(CouponLeg_CONSTPTR Leg1, CouponLeg_CONSTPTR Leg2)
{
	GenericSwap_CONSTPTR geneticSwap(new GenericSwap(Leg1,Leg2));
	return geneticSwap;
}
