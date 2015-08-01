#pragma once
#include <vector>

#include <boost/shared_ptr.hpp>

#include <LMM/LmmModel/McLmm.h>
#include <LMM/pricer/McLmmVanillaSwapPricer.h>

#include <JBLMM/Element/Rate1.h>
#include <JBLMM/Element/ConstRate.h>
#include <JBLMM/Element/LiborRate.h>
#include <JBLMM/Longstaff_Schwartz/EV_Basis/EV.h>
//! YY: need a parent class of EV, then derived to LMM_EV (main difference is the evaluator)

// ----------------------------------------------------------------------
// 
//               Base EV_Evaluator
//
// ---------------------------------------------------------------------

class EV_Evaluator
{
public:
	virtual double evaluate(EV_CONSTPTR ev, const matrix& liborMatrix, const std::vector<double>& numeraire)const=0;
	virtual ~EV_Evaluator(){}

	virtual void write_to_stream(std::ostream& out)const=0;
};
typedef boost::shared_ptr<EV_Evaluator> EV_Evaluator_PTR;
typedef boost::shared_ptr<const EV_Evaluator> EV_Evaluator_CONSTPTR;



// ----------------------------------------------------------------------
// 
//               Derived EV_Evaluator
//
// ---------------------------------------------------------------------

class EV_ConstRate_Evaluator: public EV_Evaluator
{
public:
	virtual double evaluate(EV_CONSTPTR ev, const matrix& liborMatrix, const std::vector<double>& numeraire)const;
	virtual ~EV_ConstRate_Evaluator(){}

	void write_to_stream(std::ostream& out)const;
};
typedef boost::shared_ptr<EV_ConstRate_Evaluator> EV_ConstRate_Evaluator_PTR;
typedef boost::shared_ptr<const EV_ConstRate_Evaluator> EV_ConstRate_Evaluator_CONSTPTR;


class EV_LiborRate_Evaluator: public EV_Evaluator
{
public:
	virtual double evaluate(EV_CONSTPTR ev, const matrix& liborMatrix, const std::vector<double>& numeraire)const;
	virtual ~EV_LiborRate_Evaluator(){}

	void write_to_stream(std::ostream& out)const;
};
typedef boost::shared_ptr<EV_LiborRate_Evaluator> EV_LiborRate_Evaluator_PTR;
typedef boost::shared_ptr<const EV_LiborRate_Evaluator> EV_LiborRate_Evaluator_CONSTPTR;


class EV_VanillaSwapRate_Evaluator: public EV_Evaluator
{
public:
	McLmmVanillaSwapPricer mcLmmVanillaSwapPricer_;
	EV_VanillaSwapRate_Evaluator(const McLmmVanillaSwapPricer& mcLmmVanillaSwapPricer);

	virtual double evaluate(EV_CONSTPTR ev, const matrix& liborMatrix, const std::vector<double>& numeraire)const;
	virtual ~EV_VanillaSwapRate_Evaluator(){}

	void write_to_stream(std::ostream& out)const;
};
typedef boost::shared_ptr<EV_VanillaSwapRate_Evaluator> EV_VanillaSwapRate_Evaluator_PTR;
typedef boost::shared_ptr<const EV_VanillaSwapRate_Evaluator> EV_VanillaSwapRate_Evaluator_CONSTPTR;


