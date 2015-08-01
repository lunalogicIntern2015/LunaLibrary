#pragma once

#include <boost/shared_ptr.hpp>


#include <JBLMM/Element/Rate1.h>
#include <LMM/helper/TenorType.h>
#include <LMM/helper/Name.h>





class LiborRate : public Rate1
{
	LMM::Index fixingTime_; // T_{i}
	const Tenor  duration_; // 3m, 6m 
public:
	//gettor
	LMM::Index getFixingTime()const{return fixingTime_;}
	const Tenor& getDuration()const{return duration_;}
	//constructor, destructor
	LiborRate(LMM::Index fixingTime, const Tenor duration);
	virtual ~LiborRate(){}
	//clone
	Rate1_PTR clone()const;	//deep copy
	//print
	virtual void show()const;
	
	void write_to_stream(std::ostream& out)const;
};
typedef boost::shared_ptr<LiborRate> LiborRate_PTR;
typedef boost::shared_ptr<const LiborRate> LiborRate_CONSTPTR;

