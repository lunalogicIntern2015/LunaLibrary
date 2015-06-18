#include "JBLMM/Element/Coupon.h"
#include <iostream>

Coupon::Coupon(const LMM::Index paymentIndex)
	:paymentIndex_(paymentIndex)
{
}

Coupon_PTR Coupon::clone()const
{
	return Coupon_PTR(new Coupon(*this));
}

void Coupon::show()const
{
	std::cout	<<	"paymentIndex:      "	<<		paymentIndex_		<<	std::endl;
}