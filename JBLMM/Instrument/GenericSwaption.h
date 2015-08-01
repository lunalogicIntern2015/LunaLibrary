#pragma once

#include <boost/shared_ptr.hpp>

#include <LMM/helper/Name.h>
#include <JBLMM/Instrument/GenericSwap.h>

class GenericSwaption  : public Instrument
{
	LMM::Index maturity_;
	GenericSwap_CONSTPTR genericSwap_;
	//check maturity_ <= firstindex	
	bool check()const;		
public:
	//getter
	LMM::Index getMaturity()const{return maturity_;}
	GenericSwap_CONSTPTR getGenericSwap()const{return genericSwap_;}  
	//constructor, destructor
	GenericSwaption(const LMM::Index maturity, GenericSwap_CONSTPTR genericSwap);   
	virtual ~GenericSwaption(){}					  
};
typedef boost::shared_ptr<GenericSwaption> GenericSwaption_PTR;
typedef boost::shared_ptr<const GenericSwaption> GenericSwaption_CONSTPTR;


