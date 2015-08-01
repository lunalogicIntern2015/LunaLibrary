#pragma once

#include <iostream>
#include <vector>
#include <string>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/shared_ptr.hpp>	

#include<LMM/helper/TenorType.h>
#include<LMM/helper/LMMTenorStructure.h>

#include <JBLMM/Element/CouponLeg.h>
#include <JBLMM/Instrument/Instrument.h>

#include <LMM/helper/LMMTenorStructure.h>


/*! \class VanillaSwap manage leg's INDICES in the simulation structures see <LMM\lmmModel\LMMTenorStructure.h>
 * description 
*/

class GenericSwap : public Instrument
{
	CouponLeg_CONSTPTR leg1_;
	CouponLeg_CONSTPTR leg2_;
	LMMTenorStructure_PTR lmmTenorStructure_;  //???

public:
	//gettor
	CouponLeg_CONSTPTR getLeg1()const{return leg1_;}
	CouponLeg_CONSTPTR getLeg2()const{return leg2_;}
	LMMTenorStructure_PTR getLMMTenorStructure()const{return lmmTenorStructure_;}
	// constructor, destructor
	GenericSwap(CouponLeg_CONSTPTR leg1, CouponLeg_CONSTPTR leg2);
	GenericSwap(CouponLeg_CONSTPTR leg1, CouponLeg_CONSTPTR leg2, LMMTenorStructure_PTR lmmTenorStructure);
	virtual ~GenericSwap(){}
	//get subGeneticSwap
	Instrument_CONSTPTR getSubGenericSwap(const size_t indexStart, const size_t indexEnd) const;
};
typedef boost::shared_ptr<GenericSwap> GenericSwap_PTR;
typedef boost::shared_ptr<const GenericSwap> GenericSwap_CONSTPTR;





