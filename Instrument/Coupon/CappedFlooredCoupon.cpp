#include <Instrument/Coupon/CappedFlooredCoupon.h>
#include <LMM/Helper/Name.h>
#include <iostream>

CappedFlooredCoupon::CappedFlooredCoupon(	LMM::Index paymentIndex,
											double nominal,
											double period,  
											bool ifFloored,
											double floorStrike,
											bool ifCapped,
											double capStrike,
											Rate1_CONSTPTR rate,
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


void CappedFlooredCoupon::write_to_stream(std::ostream& out)const
{
	Coupon::write_to_stream(out);
	out	<<	"nominal:          ;"	<<	nominal_			<<	std::endl;
	out	<<	"period:           ;"	<<	period_				<<	std::endl;
	out	<<	"ifFloored:        ;";	ifFloored_? (out << "Y"): (out<< "N");	out <<	std::endl;
	out	<<	"floorStrike:      ;"	<<	floorStrike_		<<	std::endl;
	out	<<	"ifCapped:         ;";	ifCapped_?  (out << "Y"): (out<< "N");	out <<	std::endl;
	out	<<	"capStrike:        ;"	<<	capStrike_			<<	std::endl;
	out	<<	"rate:             ;"	<<	std::endl;	
	rate_->write_to_stream(out);
	out	<<	"multiFactor:      ;"	<<	multiFactor_			<<	std::endl;
	out	<<	"addFactor:        ;"	<<	addFactor_				<<	std::endl;
	out	<<	"valuationDate:    ;"	<<	valuationDate_			<<	std::endl;
}