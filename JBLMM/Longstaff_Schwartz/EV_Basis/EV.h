#pragma once
#include <vector>

#include <boost/shared_ptr.hpp>

#include <LMM/LmmModel/McLmm.h>
#include <LMM/pricer/McLmmVanillaSwapPricer.h>

#include <JBLMM/Element/Rate1.h>
#include <JBLMM/Element/LiborRate.h>
#include <JBLMM/Element/VanillaSwapRate.h>
#include <JBLMM/Instrument/BermudanVanillaSwaption.h>

// ----------------------------------------------------------------------
// 
//               Base EV
//
// ----------------------------------------------------------------------

//! containing a Rate 
class EV
{
protected:
	Rate1_CONSTPTR rate_;									
public:
	EV(Rate1_CONSTPTR rate):rate_(rate){}
	Rate1_CONSTPTR getRate() const {return rate_;}			
	virtual ~EV(){}

	virtual void write_to_stream(std::ostream& out)const=0;
};
typedef boost::shared_ptr<EV> EV_PTR;
typedef boost::shared_ptr<const EV> EV_CONSTPTR;



// ----------------------------------------------------------------------
// 
//               Derived EV
//
// ----------------------------------------------------------------------

class EV_ConstRate: public EV
{						
public:
	EV_ConstRate(Rate1_CONSTPTR rate):EV(rate){}
	virtual ~EV_ConstRate(){}

	void write_to_stream(std::ostream& out)const;
};


class EV_LiborRate: public EV
{
public:
	EV_LiborRate(Rate1_CONSTPTR rate):EV(rate){}
	virtual ~EV_LiborRate(){}

	void write_to_stream(std::ostream& out)const;
};


class EV_VanillaSwapRate: public EV
{
public:
	EV_VanillaSwapRate(Rate1_CONSTPTR rate):EV(rate){}
	virtual ~EV_VanillaSwapRate(){}

	void write_to_stream(std::ostream& out)const;
};


