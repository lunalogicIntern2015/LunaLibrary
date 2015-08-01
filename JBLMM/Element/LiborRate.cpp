#include "JBLMM/Element/LiborRate.h"
#include <LMM/helper/TenorType.h>

LiborRate::LiborRate(LMM::Index fixingTime, const Tenor duration)
	:
		fixingTime_(fixingTime),
		duration_(duration)
{
}


Rate1_PTR LiborRate::clone()const
{
	return Rate1_PTR(new LiborRate(*this));
}

void LiborRate::show()const
{

	std::cout	<<	"----------------------------LiborRate-------------------------------------"	<<	std::endl;
	Rate1::show();
	std::cout	<<	"fixingTime:        "	<<	fixingTime_					<<	std::endl;
	std::cout	<<	"duration:          "	<<	duration_.NbOfMonth()		<<	" Months"		<<	std::endl;
	std::cout	<<	"----------------------------LiborRateEND----------------------------------"	<<	std::endl;
}


void LiborRate::write_to_stream(std::ostream& out)const
{
	out	<<	"----LiborRate----"	<<	std::endl;
	out	<<	"fixingTime:       ; "	<<	fixingTime_					<<	std::endl;
	out	<<	"duration (Month):         ; "	<<	duration_.NbOfMonth()		<<	std::endl;
	out	<<	"----LiborRateEND---"	<<	std::endl;	
}











