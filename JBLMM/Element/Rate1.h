#pragma once
#include <boost/shared_ptr.hpp>


class Rate1
{
public:
	virtual ~Rate1(){}
	virtual boost::shared_ptr<Rate1> clone()const;
	//print
	virtual void show()const;

	virtual void write_to_stream(std::ostream& o)const;
};
typedef boost::shared_ptr<Rate1> Rate1_PTR;
typedef boost::shared_ptr<const Rate1> Rate1_CONSTPTR;



