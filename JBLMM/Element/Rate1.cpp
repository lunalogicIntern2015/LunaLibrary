#include "JBLMM/Element/Rate1.h"
#include <boost/shared_ptr.hpp>

Rate1_PTR Rate1::clone()const
{
	return Rate1_PTR(new Rate1(*this));
}


void Rate1::show()const
{	
}

void Rate1::write_to_stream(std::ostream& o)const
{
}


