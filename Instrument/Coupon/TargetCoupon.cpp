#include <Instrument/Coupon/TargetCoupon.h>
#include <iostream>

TargetCoupon::TargetCoupon(	LMM::Index paymentIndex,
							double nominal,
							double period, 
							bool   ifFloored,
							double floorStrike,
							bool   ifCapped, 
							double capStrike, 
							Rate1_CONSTPTR rate,
							double multiFactor,
							double addFactor, 
							LMM::Index valuationDateIndex,
							double target,
							std::string couponDependency)
	:
	CappedFlooredCoupon(	paymentIndex,
							nominal,
							period, 
							ifFloored,
							floorStrike,
							ifCapped, 
							capStrike, 
							rate,
							multiFactor,
							addFactor, 
							valuationDateIndex),
	target_(target),
	couponDependency_(couponDependency)
{
}

bool TargetCoupon::isSameTarget(Coupon_CONSTPTR coupon)const
{
	TargetCoupon_CONSTPTR targetCoupon = boost::dynamic_pointer_cast<const TargetCoupon>(coupon);
	if(!targetCoupon)
		throw("fail to cast targetCoupon");

	return (getTarget()==targetCoupon->getTarget());
}

void TargetCoupon::show()const
{
	CappedFlooredCoupon::show();
	std::cout	<<	"target:            "	<<	target_						<<	std::endl;
	std::cout	<<	"couponDependency:  "	<<	couponDependency_			<<	std::endl;
}

Coupon_PTR TargetCoupon::clone()const
{
	return Coupon_PTR(new TargetCoupon(*this));
}