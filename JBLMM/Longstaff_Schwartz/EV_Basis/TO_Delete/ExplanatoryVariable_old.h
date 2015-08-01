//#pragma once
//#include <vector>
//#include <LMM/LmmModel/McLmm.h>
//#include <LMM/pricer/McLmmVanillaSwapPricer.h>
//#include <JBLMM/Element/Rate.h>
//#include <JBLMM/Instrument/BermudanSwaption.h>
//
//class ExplanatoryVariable
//{
//public:
//	class Evaluator
//	{
//		public:
//		//it cannot be static virtual because static relate to class and virtual relate to instance 
//		virtual double evaluate(McLmm_PTR mclmm,
//								BermudanSwaption_CONSTPTR bermudanSwaption_PTR,
//								LMM::Index liborIndex)const=0;
//	};
//	typedef boost::shared_ptr<Evaluator> Evaluator_PTR;
//	typedef boost::shared_ptr<const Evaluator> Evaluator_CONSTPTR;
//
//private:
//	Evaluator_CONSTPTR evaluator_PTR_;
//
//public:
//	//construtor
//	ExplanatoryVariable(Evaluator_CONSTPTR evaluator_PTR):evaluator_PTR_(evaluator_PTR){ }
//	//destructor
//	virtual ~ExplanatoryVariable(){}
//	Evaluator_CONSTPTR getEvaluator_PTR()const{return evaluator_PTR_;}
//
//};
//typedef boost::shared_ptr<ExplanatoryVariable> ExplanatoryVariable_PTR;
//typedef boost::shared_ptr<const ExplanatoryVariable> ExplanatoryVariable_CONSTPTR;
//
//
//typedef boost::numeric::ublas::matrix<double> matrix;
//class EV_Libor: public ExplanatoryVariable::Evaluator
//{
//	
//public:
//	double evaluate(McLmm_PTR mclmm,
//					BermudanSwaption_CONSTPTR bermudanSwaption_PTR,
//					LMM::Index liborIndex)const;
//};
//
//
//class EV_vanillaswaprate: public ExplanatoryVariable::Evaluator
//{
//	
//	mutable McLmmVanillaSwapPricer_PTR mcLmmVanillaSwapPricer_PTR_;
//
//public:
//	EV_vanillaswaprate(McLmmVanillaSwapPricer_PTR mcLmmVanillaSwapPricer_PTR);
//
//	double evaluate(McLmm_PTR mclmm,
//					BermudanSwaption_CONSTPTR bermudanSwaption_PTR,
//					LMM::Index liborIndex)const;
//};
//
//
//
//
//class Compositor
//{
//public:
//	virtual double compose(double x, double y)const=0;
//};
//typedef boost::shared_ptr<Compositor> Compositor_PTR;
//typedef boost::shared_ptr<const Compositor> Compositor_CONSTPTR;
//
//class ComposedEV
//{
//
//protected:
//	ExplanatoryVariable_CONSTPTR ev1_;
//	ExplanatoryVariable_CONSTPTR ev2_;
//	Compositor_PTR compositor_ptr_;
//
//public:
//	ComposedEV(	ExplanatoryVariable_CONSTPTR ev1,
//				ExplanatoryVariable_CONSTPTR ev2,
//				Compositor_PTR compositor_ptr)
//				:
//				ev1_(ev1),
//				ev2_(ev2),
//				compositor_ptr_(compositor_ptr)
//	{	
//	}
//
//	double evaluate(McLmm_PTR mclmm,
//					BermudanSwaption_CONSTPTR bermudanSwaption_PTR,
//					LMM::Index liborIndex)const
//	{
//		double ev_value1 = ev1_->getEvaluator_PTR()->evaluate(mclmm,bermudanSwaption_PTR,liborIndex);
//		double ev_value2 = ev2_->getEvaluator_PTR()->evaluate(mclmm,bermudanSwaption_PTR,liborIndex);
//		double result = compositor_ptr_->compose(ev_value1,ev_value2);
//		return result;
//	}				
//};
//typedef boost::shared_ptr<ComposedEV> ComposedEV_PTR;
//typedef boost::shared_ptr<const ComposedEV> ComposedEV_CONSTPTR;
//
//class MutiplicationCompositor: public Compositor
//{
//public:
//	double compose(double x, double y)const{return x*y;}
//};
//
//class DivisionCompositor: public Compositor
//{
//public:
//	double compose(double x, double y)const{assert(y!=0.0); return x/y;}
//};
//
//class AdditionCompositor: public Compositor
//{
//public:
//	double compose(double x, double y)const{return x+y;}
//};
//
//class SubtractionCompositor: public Compositor
//{
//public:
//	double compose(double x, double y)const{return x-y;}
//};
//
