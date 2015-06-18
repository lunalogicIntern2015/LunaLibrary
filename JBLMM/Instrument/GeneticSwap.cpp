#include "JBLMM/Instrument/GeneticSwap.h"

GeneticSwap::GeneticSwap(CouponLeg_CONSTPTR leg1, CouponLeg_CONSTPTR leg2)
	:
		leg1_(leg1),
		leg2_(leg2)
{
}

 GeneticSwap_PTR GeneticSwap::getSubGeneticSwap(const size_t indexStart, const size_t indexEnd) const
{
	return GeneticSwap_PTR(new GeneticSwap(getLeg1()->getSubCouponLeg(indexStart,indexEnd), 
		                                   getLeg2()->getSubCouponLeg(indexStart,indexEnd)));
}

 ////INDEX de 1er fixing, pas de flux !!!
 //size_t  GeneticSwap::getIndexStart() const
 //{
	//std::vector<Coupon_CONSTPTR> vectCoupon1 = getLeg1()->getLeg() ;
	//std::vector<Coupon_CONSTPTR> vectCoupon1 = getLeg2()->getLeg() ;
	//size_t firstFixingIndex ;
	//for (size_t i = 0 ; i < vectCoupon1.size() ; ++i)
	//{
	//	vectCoupon1[i]->
	//}
 //}

 ////INDEX de dernier FLUX !!
 //size_t  GeneticSwap::getIndexEnd() const	
 //{
 //
 //}


void GeneticSwap::show() const
{
	std::vector<Coupon_CONSTPTR> couponLeg1 = this->getLeg1()->getLeg() ;
	std::vector<Coupon_CONSTPTR> couponLeg2 = this->getLeg2()->getLeg() ;
//leg 1
	Coupon_CONSTPTR firstCoupon1=couponLeg1[0];
	CappedFlooredCoupon_CONSTPTR firstCappedFlooredCoupon1 = boost::dynamic_pointer_cast<const CappedFlooredCoupon>(firstCoupon1);
	if(!firstCappedFlooredCoupon1)
		throw("fail to cast cappedFlooredCoupon");

	Rate_CONSTPTR firstRate1 = firstCappedFlooredCoupon1->getRate();
	Tenor tenorLeg1 = Tenor::_Non ;										// bad code, sorry 
	double fixedRate1 ;
	if(boost::dynamic_pointer_cast<const LiborRate>(firstRate1))
	{
		LiborRate_CONSTPTR	firstLiborRate		=	boost::dynamic_pointer_cast<const LiborRate>(firstRate1);
		tenorLeg1 = firstLiborRate->getDuration() ;
	}
	else if(boost::dynamic_pointer_cast<const ConstRate>(firstRate1))
	{
		ConstRate_CONSTPTR	constRate		=	boost::dynamic_pointer_cast<const ConstRate>(firstRate1);
		tenorLeg1 = Tenor::_Non ;
		fixedRate1 = constRate->getConstRateValue() ; 
	}
	else
	{
		throw("fail to cast LiborRate or ConstRate");
	}
//leg 2
	Coupon_CONSTPTR firstCoupon2=couponLeg2[0];
	CappedFlooredCoupon_CONSTPTR firstCappedFlooredCoupon2 = boost::dynamic_pointer_cast<const CappedFlooredCoupon>(firstCoupon2);
	if(!firstCappedFlooredCoupon2)
		throw("fail to cast cappedFlooredCoupon");

	Rate_CONSTPTR firstRate2 = firstCappedFlooredCoupon2->getRate();
	Tenor tenorLeg2 = Tenor::_Non ;										// bad code, sorry 
	double fixedRate2 ;
	if(boost::dynamic_pointer_cast<const LiborRate>(firstRate2))
	{
		LiborRate_CONSTPTR	firstLiborRate		=	boost::dynamic_pointer_cast<const LiborRate>(firstRate2);
		tenorLeg2 = firstLiborRate->getDuration() ;
	}
	else if(boost::dynamic_pointer_cast<const ConstRate>(firstRate2))
	{
		ConstRate_CONSTPTR	constRate		=	boost::dynamic_pointer_cast<const ConstRate>(firstRate2);
		tenorLeg2 = Tenor::_Non ;
		fixedRate2 = constRate->getConstRateValue() ; 
	}
	else
	{
		throw("fail to cast LiborRate or ConstRate");
	}

	std::cout << "---------------------------------------------" << std::endl ;
	std::cout << "--- creation d'un objet GeneticSwap --------" << std::endl ;
//leg1	
	if (tenorLeg1 == Tenor::_Non)
	{
		std::cout << "LEG 1 TenorType_ " << tenorLeg1 << "  fixed rate : " << fixedRate1 << std::endl ;
	}
	else
	{
		std::cout << "LEG 1 TenorType_ " << tenorLeg1 << std::endl ;
	}
//leg2
	if (tenorLeg2 == Tenor::_Non)
	{
		std::cout << "LEG 2 TenorType_ " << tenorLeg2 << "  fixed rate : " << fixedRate2 << std::endl ;
	}
	else
	{
		std::cout << "LEG 2 TenorType_ " << tenorLeg2 << std::endl ;
	}

	//std::cout << "indexStart_           " << indexStart_ << std::endl ;
	//std::cout << "indexEnd_             " << indexEnd_ << std::endl ;
	
	std::cout << "indices des flux Leg 1" << std::endl ;
	for (size_t i= 0 ; i < couponLeg1.size() ; ++i)
	{
		std::cout << couponLeg1[i]->getPaymentIndex() << "  " ;
	}
	std::cout << " " << std::endl ;
	std::cout << "indices des flux Leg 2" << std::endl ;
	for (size_t i= 0 ; i < couponLeg2.size() ; ++i)
	{
		std::cout << couponLeg2[i]->getPaymentIndex() << "  " ;
	}
	std::cout << " " << std::endl ;
	//std::cout << "deltaTFloatingLeg_ i" << std::endl ;
	//for (size_t i = 0 ; i < deltaTFloatingLeg_.size() ; ++i)
	//{
	//	std::cout << deltaTFloatingLeg_[i] << "  " ;
	//} 
	//std::cout << " " << std::endl ;
	//std::cout << "deltaTFixedLeg_ i" << std::endl ;
	//for (size_t i = 0 ; i < deltaTFixedLeg_.size() ; ++i)
	//{
	//	std::cout << deltaTFixedLeg_[i] << "  " ;
	//} 
	//std::cout << " " << std::endl ;
	std::cout << "---------------------------------------------" << std::endl ;
}

void GeneticSwap::print(std::ostream& o) const
{
	std::vector<Coupon_CONSTPTR> couponLeg1 = this->getLeg1()->getLeg() ;
	std::vector<Coupon_CONSTPTR> couponLeg2 = this->getLeg2()->getLeg() ;
//leg 1
	Coupon_CONSTPTR firstCoupon1=couponLeg1[0];
	CappedFlooredCoupon_CONSTPTR firstCappedFlooredCoupon1 = boost::dynamic_pointer_cast<const CappedFlooredCoupon>(firstCoupon1);
	if(!firstCappedFlooredCoupon1)
		throw("fail to cast cappedFlooredCoupon");

	Rate_CONSTPTR firstRate1 = firstCappedFlooredCoupon1->getRate();
	Tenor tenorLeg1 = Tenor::_3M ;										// bad code, sorry 
	double fixedRate1 ;
	if(boost::dynamic_pointer_cast<const LiborRate>(firstRate1))
	{
		LiborRate_CONSTPTR	firstLiborRate		=	boost::dynamic_pointer_cast<const LiborRate>(firstRate1);
		tenorLeg1 = firstLiborRate->getDuration() ;
	}
	else if(boost::dynamic_pointer_cast<const ConstRate>(firstRate1))
	{
		ConstRate_CONSTPTR	constRate		=	boost::dynamic_pointer_cast<const ConstRate>(firstRate1);
		tenorLeg1 = Tenor::_Non ;
		fixedRate1 = constRate->getConstRateValue() ; 
	}
	else
	{
		throw("fail to cast LiborRate or ConstRate");
	}
//leg 2
	Coupon_CONSTPTR firstCoupon2=couponLeg2[0];
	CappedFlooredCoupon_CONSTPTR firstCappedFlooredCoupon2 = boost::dynamic_pointer_cast<const CappedFlooredCoupon>(firstCoupon2);
	if(!firstCappedFlooredCoupon2)
		throw("fail to cast cappedFlooredCoupon");

	Rate_CONSTPTR firstRate2 = firstCappedFlooredCoupon2->getRate();
	Tenor tenorLeg2 = Tenor::_3M ;										// bad code, sorry 
	double fixedRate2 ;
	if(boost::dynamic_pointer_cast<const LiborRate>(firstRate2))
	{
		LiborRate_CONSTPTR	firstLiborRate		=	boost::dynamic_pointer_cast<const LiborRate>(firstRate2);
		tenorLeg2 = firstLiborRate->getDuration() ;
	}
	else if(boost::dynamic_pointer_cast<const ConstRate>(firstRate2))
	{
		ConstRate_CONSTPTR	constRate		=	boost::dynamic_pointer_cast<const ConstRate>(firstRate2);
		tenorLeg2 = Tenor::_Non ;
		fixedRate2 = constRate->getConstRateValue() ; 
	}
	else
	{
		throw("fail to cast LiborRate or ConstRate");
	}

	o << "---------------------------------------" << std::endl ;
	o << "------------- GeneticSwap -------------" << std::endl ;
//leg1	
	if (tenorLeg1 == Tenor::_Non)
	{
		o << "LEG 1 TenorType_ " << tenorLeg1 << "  fixed rate : " << fixedRate1 << std::endl ;
	}
	else
	{
		o << "LEG 1 TenorType_ " << tenorLeg1 << std::endl ;
	}
//leg2
	if (tenorLeg2 == Tenor::_Non)
	{
		o << "LEG 2 TenorType_ " << tenorLeg2 << "  fixed rate : " << fixedRate1 << std::endl ;
	}
	else
	{
		o << "LEG 2 TenorType_ " << tenorLeg2 << std::endl ;
	}

	//std::cout << "indexStart_           " << indexStart_ << std::endl ;
	//std::cout << "indexEnd_             " << indexEnd_ << std::endl ;
	
	o << "indices des flux Leg 1" << std::endl ;
	for (size_t i= 0 ; i < couponLeg1.size() ; ++i)
	{
		o << couponLeg1[i]->getPaymentIndex() << "  " ;
	}
	o << " " << std::endl ;
	o << "indices des flux Leg 2" << std::endl ;
	for (size_t i= 0 ; i < couponLeg2.size() ; ++i)
	{
		o << couponLeg2[i]->getPaymentIndex() << "  " ;
	}
	o << " " << std::endl ;
	//std::cout << "deltaTFloatingLeg_ i" << std::endl ;
	//for (size_t i = 0 ; i < deltaTFloatingLeg_.size() ; ++i)
	//{
	//	std::cout << deltaTFloatingLeg_[i] << "  " ;
	//} 
	//std::cout << " " << std::endl ;
	//std::cout << "deltaTFixedLeg_ i" << std::endl ;
	//for (size_t i = 0 ; i < deltaTFixedLeg_.size() ; ++i)
	//{
	//	std::cout << deltaTFixedLeg_[i] << "  " ;
	//} 
	//std::cout << " " << std::endl ;
	o << "---------------------------------------" << std::endl ;
}