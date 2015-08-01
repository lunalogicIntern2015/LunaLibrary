#include "JBLMM/Element/VanillaSwapRate.h"
#include <boost/shared_ptr.hpp>

VanillaSwapRate::VanillaSwapRate(const VanillaSwap& vanillaSwap)
	:
	vanillaSwap_(vanillaSwap)
{
}

Rate1_PTR VanillaSwapRate::clone()const
{
	return Rate1_PTR(new VanillaSwapRate(*this));
}

void VanillaSwapRate::write_to_stream(std::ostream& o)const
{
	o << "--- VanillaSwapRate Info ---" <<std::endl;
	getVanillaSwap().write_to_stream(o);
}


