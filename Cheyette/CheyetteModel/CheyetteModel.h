#pragma once

#include <boost/shared_ptr.hpp>

class CheyetteModel  
{
	virtual ~CheyetteModel(){};
};

typedef boost::shared_ptr<CheyetteModel>       CheyetteModel_PTR;
typedef boost::shared_ptr<const CheyetteModel> CheyetteModel_CONSTPTR;











