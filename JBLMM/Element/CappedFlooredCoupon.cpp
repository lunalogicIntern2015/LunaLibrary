#include "JBLMM/Element/CappedFlooredCoupon.h"
#include <LMM/helper/Name.h>
#include <iostream>

CappedFlooredCoupon::CappedFlooredCoupon(	LMM::Index paymentIndex,
											double nominal,
											double period,  
											bool ifFloored,
											double floorStrike,
											bool ifCapped,
											double capStrike,
											Rate_CONSTPTR rate,
											double multiFactor,
											double addFactor,
											LMM::Index valuationDateIndex)
										 : 
											Coupon(paymentIndex),
											nominal_(nominal),
											period_(period),
											ifFloored_(ifFloored),
											floorStrike_(floorStrike),
											ifCapped_(ifCapped),
											capStrike_(capStrike),
											rate_(rate),
											multiFactor_(multiFactor),
											addFactor_(addFactor),
											valuationDate_(valuationDateIndex)
{
}

//TODO?latter too long to write the code ... 
CappedFlooredCoupon::CappedFlooredCoupon(const CappedFlooredCoupon& c): Coupon(c.getPaymentIndex()),
																		nominal_(c.getNominal()),
																		period_(c.getPeriod()),
																		ifFloored_(c.getIfFloored()),
																		floorStrike_(c.getFloorStrike()),
																		ifCapped_(c.getIfCapped()),
																		capStrike_(c.getCapStrike()),
																		rate_((c.getRate()->clone())),
																		multiFactor_(c.getMultiFactor()),
																		addFactor_(c.getAddFactor()),
																		valuationDate_(c.getValuationDate())
{		
}


Coupon_PTR CappedFlooredCoupon::clone()const
{
	return Coupon_PTR(new CappedFlooredCoupon(*this));
}

void CappedFlooredCoupon::show()const
{
	Coupon::show();
	std::cout	<<	"nominal:           "	<<	nominal_			<<	std::endl;
	std::cout	<<	"period:            "	<<	period_				<<	std::endl;
	std::cout	<<	"ifFloored:         ";	ifFloored_? std::cout << "Y": std::cout<< "N"	<<	std::endl;
	std::cout	<<	"floorStrike:       "	<<	floorStrike_		<<	std::endl;
	std::cout	<<	"ifCapped:          ";	ifCapped_?  std::cout << "Y": std::cout<< "N"	<<	std::endl;
	std::cout	<<	"capStrike:         "	<<	capStrike_			<<	std::endl;
	std::cout	<<	"rate:              "	<<	std::endl;	
	rate_->show();
	std::cout	<<	"multiFactor:       "	<<	multiFactor_			<<	std::endl;
	std::cout	<<	"addFactor:         "	<<	addFactor_				<<	std::endl;
	std::cout	<<	"valuationDate:     "	<<	valuationDate_			<<	std::endl;
}