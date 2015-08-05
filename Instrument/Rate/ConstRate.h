#pragma once
#include <Instrument/Rate/Rate1.h>

class ConstRate : public Rate1
{
	const double constRateValue_;

public:
	//getter
	const double getConstRateValue()const{return constRateValue_;}
	//constructor
	ConstRate(const double constRateValue);
	//clone
	Rate1_PTR clone()const;
	//print	
	virtual void show()const;

	void write_to_stream(std::ostream& out)const;
};
typedef boost::shared_ptr<ConstRate> ConstRate_PTR;
typedef boost::shared_ptr<const ConstRate> ConstRate_CONSTPTR;

