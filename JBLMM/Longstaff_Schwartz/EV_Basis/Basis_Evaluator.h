#pragma once
#include <vector>

#include <boost/shared_ptr.hpp>

#include <LMM/LmmModel/McLmm.h>
#include <LMM/pricer/McLmmVanillaSwapPricer.h>

#include <JBLMM/Element/Rate1.h>
#include <JBLMM/Element/CappedFlooredCoupon.h>
#include <JBLMM/Instrument/BermudanVanillaSwaption.h>
#include <JBLMM/Longstaff_Schwartz/EV_Basis/EV_Evaluator.h>
#include <JBLMM/Longstaff_Schwartz/EV_Basis/Basis.h>

//!! YY 
//! Basis = f(EV) or a composition of Basis (polynormials )
//! YY taking firstly transformer on double, then on general EV+Evaluator 

// ----------------------------------------------------------------------
// 
//               Base Basis_Evaluator
//
// ---------------------------------------------------------------------

class Basis_Evaluator
{
public:
	virtual ~Basis_Evaluator(){}

	virtual double evaluate(Basis_CONSTPTR basis, const matrix& liborMatrix, const std::vector<double>& numeraire)const=0;

	virtual void write_to_stream(std::ostream& out)const=0;
};
typedef boost::shared_ptr<Basis_Evaluator> Basis_Evaluator_PTR;
typedef boost::shared_ptr<const Basis_Evaluator> Basis_Evaluator_CONSTPTR;


// ----------------------------------------------------------------------
// 
//                Basis_SinglaEVFunctional_Evaluator
//
// ----------------------------------------------------------------------
class Basis_SinglaEVFunctional_Evaluator : public Basis_Evaluator
{
	EV_Evaluator_CONSTPTR ev_evaluator_;
public:
	//constructor et destructor
	Basis_SinglaEVFunctional_Evaluator(EV_Evaluator_CONSTPTR ev_evaluator) : ev_evaluator_(ev_evaluator){}
	virtual ~Basis_SinglaEVFunctional_Evaluator(){}

	virtual double evaluate(Basis_CONSTPTR basis, const matrix& liborMatrix, const std::vector<double>& numeraire)const;

	void write_to_stream(std::ostream& out)const;
};
typedef boost::shared_ptr<Basis_SinglaEVFunctional_Evaluator> Basis_SinglaEVFunctional_Evaluator_PTR;
typedef boost::shared_ptr<const Basis_SinglaEVFunctional_Evaluator> Basis_SinglaEVFunctional_Evaluator_CONSTPTR;


/*
class Basis_ConstUnity_Evaluator : public Basis_SinglaEVFunctional_Evaluator  // consider 1 as a basis 
{
public:
	virtual double evaluate(const matrix& liborMatrix, const std::vector<double>& numeraire)const;
	virtual ~Basis_SinglaEVFunctional_Evaluator(){}
};
typedef boost::shared_ptr<Basis_ConstUnity_Evaluator> Basis_ConstUnity_Evaluator_PTR;
typedef boost::shared_ptr<const Basis_ConstUnity_Evaluator> Basis_ConstUnity_Evaluator_CONSTPTR;



class Basis_CappedFlooredCoupon_Evaluator: public Basis_SinglaEVFunctional_Evaluator
{
public:
	virtual double evaluate(const matrix& liborMatrix, const std::vector<double>& numeraire)const;
	virtual ~Basis_CappedFlooredCoupon_Evaluator(){}
};
typedef boost::shared_ptr<Basis_CappedFlooredCoupon_Evaluator> Basis_CappedFlooredCoupon_Evaluator_PTR;
typedef boost::shared_ptr<const Basis_CappedFlooredCoupon_Evaluator> Basis_CappedFlooredCoupon_Evaluator_CONSTPTR;




class Basis_Monomial_Evaluator : public Basis_SinglaEVFunctional_Evaluator
{
public:
	virtual double evaluate(const matrix& liborMatrix, const std::vector<double>& numeraire)const;
	virtual ~Basis_Monomial_Evaluator(){}
};
typedef boost::shared_ptr<Basis_Monomial_Evaluator> Basis_Monomial_Evaluator_PTR;
typedef boost::shared_ptr<const Basis_Monomial_Evaluator> Basis_Monomial_Evaluator_CONSTPTR;
*/


// ----------------------------------------------------------------------
// 
//                Basis_Composited_Evaluator
//
// ----------------------------------------------------------------------

//! TODO: change basis_Evaluator1_, basis_Evaluator2_ to more user friendly form !!!!! 

class Basis_Composited_Evaluator : public Basis_Evaluator  // composit tow basis by operator: +,-,*,/ 
{
protected:
	Basis_Evaluator_CONSTPTR basis_Evaluator1_;
	Basis_Evaluator_CONSTPTR basis_Evaluator2_;
public:
	Basis_Composited_Evaluator(	Basis_Evaluator_CONSTPTR basis_Evaluator1, Basis_Evaluator_CONSTPTR basis_Evaluator2);

	virtual double evaluate(Basis_CONSTPTR basis, const matrix& liborMatrix, const std::vector<double>& numeraire)const;
	virtual ~Basis_Composited_Evaluator(){}

	void write_to_stream(std::ostream& out)const;
};
typedef boost::shared_ptr<Basis_Composited_Evaluator> Basis_Composited_Evaluator_PTR;
typedef boost::shared_ptr<const Basis_Composited_Evaluator> Basis_Composited_Evaluator_CONSTPTR;



class Basis_Composited_Function_Evaluator : public Basis_Composited_Evaluator  // composit tow basis by operator: °
{

public:
	Basis_Composited_Function_Evaluator(	Basis_Evaluator_CONSTPTR basis_Evaluator1, Basis_Evaluator_CONSTPTR basis_Evaluator2);

	virtual double evaluate(Basis_CONSTPTR basis, const matrix& liborMatrix, const std::vector<double>& numeraire)const;
	virtual ~Basis_Composited_Function_Evaluator(){}

	void write_to_stream(std::ostream& out)const;
};
typedef boost::shared_ptr<Basis_Composited_Evaluator> Basis_Composited_Evaluator_PTR;
typedef boost::shared_ptr<const Basis_Composited_Evaluator> Basis_Composited_Evaluator_CONSTPTR;



// ----------------------------------------------------------------------
// 
//                Basis_Polynomial_Evaluator
//
// ----------------------------------------------------------------------
class Basis_Polynomial_Evaluator : public Basis_Evaluator  // composit tow basis by operator: °
{
	std::vector<Basis_Evaluator_CONSTPTR> basis_evaluator_vect_;
public:
	Basis_Polynomial_Evaluator(const std::vector<Basis_Evaluator_CONSTPTR>& basis_evaluator_vect);

	virtual double evaluate(Basis_CONSTPTR basis, const matrix& liborMatrix, const std::vector<double>& numeraire)const;
	virtual ~Basis_Polynomial_Evaluator(){}

	void write_to_stream(std::ostream& out)const;
};
typedef boost::shared_ptr<Basis_Composited_Evaluator> Basis_Composited_Evaluator_PTR;
typedef boost::shared_ptr<const Basis_Composited_Evaluator> Basis_Composited_Evaluator_CONSTPTR;
