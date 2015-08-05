#pragma once
#include <boost/shared_ptr.hpp>

//#define REGRESSION_TIME_TEST
//#define TIME_TEST
#define DEBUG
#define BACKWARD_TIME_TEST
#define FORWARD_TIME_TEST

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



