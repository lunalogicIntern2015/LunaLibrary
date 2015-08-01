//#pragma once
//#include <vector>
//
//#include <JBLMM/Element/Rate.h>
//#include <JBLMM/Element/CappedFlooredCoupon.h>
//#include <JBLMM/Longstaff_Schwartz/ExplanatoryVariable.h>
//
//class Basis
//{
//public: 
//	class Transformer
//	{
//		public:
//		virtual double evaluate_x(double x)const=0;
//	};
//	typedef boost::shared_ptr<Transformer> Transformer_PTR;
//	typedef boost::shared_ptr<const Transformer> Transformer_CONSTPTR;
//
//private:
//	Transformer_PTR transformer_;
//	
//public:
//	//construtor
//	Basis(const Transformer_PTR transformer):transformer_(transformer){}
//	//destructor
//	virtual ~Basis(){}
//	const Transformer_PTR getTransformer_PTR()const{return transformer_;}
//	virtual double evaluateEV(McLmm_PTR mclmm, BermudanSwaption_CONSTPTR bermudanSwaption_ptr, LMM::Index liborIndex )const=0;
//	virtual double transform(double x)const=0;
//};
//typedef boost::shared_ptr<Basis> Basis_PTR;
//typedef boost::shared_ptr<const Basis> Basis_CONSTPTR;
//
//class CappedFlooredCouponTransformer: public Basis::Transformer
//{
//	CappedFlooredCoupon_PTR cappedFlooredCoupon_;
//public:
//	//constructor
//	CappedFlooredCouponTransformer(CappedFlooredCoupon_PTR cappedFlooredCoupon):cappedFlooredCoupon_(cappedFlooredCoupon){}
//	//destructor
//	virtual ~CappedFlooredCouponTransformer(){}
//	//gettor
//	CappedFlooredCoupon_PTR getCappedFlooredCoupon_PTR()const{return cappedFlooredCoupon_;}
//	double evaluate_x(double x)const; 
//};
//
//
//class ConstantTransformer : public Basis::Transformer
//{
//	double constant_;
//public:
//	ConstantTransformer(double constant):constant_(constant){}
//	virtual ~ConstantTransformer(){}
//	double evaluate_x(double x)const{ return constant_;}
//};
//
//
//class Polynomial : public Basis::Transformer
//{
//	const std::vector<size_t>& coef_;
//public:
//	Polynomial(const std::vector<size_t>& coef):coef_(coef){}																		
//	virtual ~Polynomial(){}
//	double evaluate_x(double x)const;
//};
//
//
//
//
//
//class SingleEvBasis: public Basis
//{
//	ExplanatoryVariable_CONSTPTR ev_ptr_;
//	//mutable double ev_value_;
//public: 
//	SingleEvBasis(Transformer_PTR transformer, ExplanatoryVariable_CONSTPTR ev_ptr):Basis(transformer),ev_ptr_(ev_ptr){}
//	virtual ~SingleEvBasis(){}
//	double evaluateEV(McLmm_PTR mclmm, BermudanSwaption_CONSTPTR bermudanSwaption_ptr, LMM::Index liborIndex )const;
//	double transform(double x)const{return getTransformer_PTR()->evaluate_x(x);}
//};
//
//class ComposedEvBasis: public Basis
//{
//	ComposedEV_PTR composedEV_ptr_;
//	//mutable double ev_value_;
//public: 
//	ComposedEvBasis(Transformer_PTR transformer, ComposedEV_PTR composedEV_ptr):Basis(transformer),composedEV_ptr_(composedEV_ptr){}
//	virtual ~ComposedEvBasis(){}
//	double evaluateEV(McLmm_PTR mclmm, BermudanSwaption_CONSTPTR bermudanSwaption_ptr, LMM::Index liborIndex )const;
//	double transform(double x)const{return getTransformer_PTR()->evaluate_x(x);}
//};
//
////class ConstantBasis: public Basis
////{
////	const double constant_;
////public: 
////	ConstantBasis(Transformer_PTR transformer,const double constant):Basis(transformer),constant_(constant){}
////	virtual ~ConstantBasis(){}
////	double transform(double x)const{return constant_;}
////};