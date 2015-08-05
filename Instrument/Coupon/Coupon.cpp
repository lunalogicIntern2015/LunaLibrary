#include <Instrument/Coupon/Coupon.h>
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

void Coupon::write_to_stream(std::ostream& out)const
{
	out	<<	"paymentIndex:      "	<<		paymentIndex_		<<	std::endl;
}