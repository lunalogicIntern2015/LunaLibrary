#pragma once

#include <boost/shared_ptr.hpp>


#include <JBLMM/Element/Rate.h>
#include <LMM/helper/TenorType.h>
#include <LMM/helper/Name.h>

class ConstRate : public Rate1
{
	const double constRateValue_;

public:
	//getter
	double getConstRateValue()const{return constRateValue_;}
	//constructor
	ConstRate(const double constRateValue);
	//clone
	Rate_PTR clone()const;
	//print	
	virtual void show()const;
};
typedef boost::shared_ptr<ConstRate> ConstRate_PTR;
typedef boost::shared_ptr<const ConstRate> ConstRate_CONSTPTR;



class LiborRate : public Rate1
{
	LMM::Index fixingTime_; // T_{i}
	const Tenor  duration_; // 3m, 6m c'est le tenor......²
public:
	//gettor
	LMM::Index getFixingTime()const{return fixingTime_;}
	const Tenor& getDuration()const{return duration_;}
	//constructor, destructor
	LiborRate(LMM::Index fixingTime, const Tenor duration);
	virtual ~LiborRate(){}
	//clone
	Rate_PTR clone()const;	//deep copy
	//print
	virtual void show()const;
};
typedef boost::shared_ptr<LiborRate> LiborRate_PTR;
typedef boost::shared_ptr<const LiborRate> LiborRate_CONSTPTR;

