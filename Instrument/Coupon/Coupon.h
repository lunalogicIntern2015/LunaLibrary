#pragma once
#include <boost/shared_ptr.hpp>

#include <Instrument/Rate/Rate1.h>
#include <LMM/Helper/Name.h>

class Coupon
{
protected:
	LMM::Index paymentIndex_;  //T_i
public:
	//getter
	LMM::Index getPaymentIndex()const{return paymentIndex_;}
	//constructor, destructor
	Coupon(const LMM::Index paymentIndex);
	virtual ~Coupon(){} 
	//clone
	virtual boost::shared_ptr<Coupon> clone()const;	 
	//print Coupon
	virtual void show()const;

	virtual void write_to_stream(std::ostream& out)const;
};
typedef boost::shared_ptr<Coupon> Coupon_PTR;
typedef boost::shared_ptr<const Coupon> Coupon_CONSTPTR;

