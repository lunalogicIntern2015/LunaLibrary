#include "JBLMM/Element/LiborRate.h"
#include <LMM/helper/TenorType.h>

LiborRate::LiborRate(LMM::Index fixingTime, const Tenor duration)
	:
		fixingTime_(fixingTime),
		duration_(duration)
{
}


Rate_PTR LiborRate::clone()const
{
	return Rate_PTR(new LiborRate(*this));
}

void LiborRate::show()const
{

	std::cout	<<	"----------------------------LiborRate-------------------------------------"	<<	std::endl;
	Rate1::show();
	std::cout	<<	"fixingTime:        "	<<	fixingTime_					<<	std::endl;
	std::cout	<<	"duration:          "	<<	duration_.NbOfMonth()		<<	" Months"		<<	std::endl;
	std::cout	<<	"----------------------------LiborRateEND----------------------------------"	<<	std::endl;
}

ConstRate::ConstRate(const double constRateValue)
	:
	constRateValue_(constRateValue)
{
}

Rate_PTR ConstRate::clone()const
{
	return Rate_PTR(new ConstRate(*this));
}

void ConstRate::show()const
{
	std::cout	<<	"----------------------------ConstRate-------------------------------------"	<<	std::endl;
	Rate1::show();
	std::cout	<<	"constRateValue:    "	<<	constRateValue_					<<	std::endl;
	std::cout	<<	"----------------------------ConstRateEND----------------------------------"	<<	std::endl;
}







