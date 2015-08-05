#pragma once

#include <boost/shared_ptr.hpp>

#include<Instrument/Rate/Rate1.h>
#include<LMM/Helper/TenorType.h>
#include<Instrument/VanillaSwap.h>


class VanillaSwapRate : public Rate1
{
	const VanillaSwap vanillaSwap_;
public:
	//getter
	const VanillaSwap& getVanillaSwap()const{return vanillaSwap_;}
	//constructor, destructor
	VanillaSwapRate(const VanillaSwap& vanillaSwap);
	virtual ~VanillaSwapRate(){}
	//clone
	Rate1_PTR clone()const;

	void write_to_stream(std::ostream& o)const;

};
typedef boost::shared_ptr<VanillaSwapRate> VanillaSwapRate_PTR;
typedef boost::shared_ptr<const VanillaSwapRate> VanillaSwapRate_CONSTPTR;

