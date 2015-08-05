#include <Instrument/Rate/ConstRate.h>
#include <iostream>

ConstRate::ConstRate(const double constRateValue)
	:
	constRateValue_(constRateValue)
{
}

Rate1_PTR ConstRate::clone()const
{
	return Rate1_PTR(new ConstRate(*this));
}

void ConstRate::show()const
{
	std::cout	<<	"---ConstRate----"	<<	std::endl;
	Rate1::show();
	std::cout	<<	"constRateValue:    "	<<	constRateValue_					<<	std::endl;
	std::cout	<<	"---ConstRateEND----"	<<	std::endl;
}

void ConstRate::write_to_stream(std::ostream& out)const
{
	out	<<	"---ConstRate----"	<<	std::endl;
	out	<<	"constRateValue:    ;"	<<	constRateValue_					<<	std::endl;
	out	<<	"---ConstRateEND----"	<<	std::endl;
}