#pragma once
#include <boost/shared_ptr.hpp>
class Instrument
{
public:
	virtual ~Instrument(){}
};
typedef boost::shared_ptr<Instrument> Instrument_PTR;
typedef boost::shared_ptr<const Instrument> Instrument_CONSTPTR;

