#pragma once
#include <vector>

#include <boost/shared_ptr.hpp>

#include <Instrument/Rate/Rate1.h>
#include <Instrument/Rate/ConstRate.h>
#include <Instrument/Coupon/CappedFlooredCoupon.h>
#include <LMM/Pricer/Longstaff_Schwartz/EV.h>
#include <LMM/Pricer/Longstaff_Schwartz/EV_Evaluator.h>


//!! TODO: YY  
//!  think about the efficency problem latter: avoid to evaluate many times the EV in the basis



class Basis 
{
public:
	virtual ~Basis(){};
	virtual void write_to_stream(std::ostream& out)const;
};
typedef boost::shared_ptr<Basis> Basis_PTR;
typedef boost::shared_ptr<const Basis> Basis_CONSTPTR;

 


// ----------------------------------------------------------------------
// 
//                Basis_SinglaEVFunctional
//
// ----------------------------------------------------------------------
class Basis_SinglaEVFunctional : public Basis
{

public:
	class Transformer
	{
		public:
		virtual double transform(double val)const=0;  // val can be ev_val or basis_val ! 
		virtual void write_to_stream(std::ostream& out)const=0;
	};

	typedef boost::shared_ptr<Transformer> Transformer_PTR;
	typedef boost::shared_ptr<const Transformer> Transformer_CONSTPTR;

protected:
	EV_CONSTPTR ev_;
	Transformer_CONSTPTR transformer_;
	
public:
	//construtor
	Basis_SinglaEVFunctional (EV_CONSTPTR ev, Transformer_CONSTPTR transformer): ev_(ev), transformer_(transformer){}
	//getter
	EV_CONSTPTR getEV()const{return ev_;}
	Transformer_CONSTPTR getTransformer()const{return transformer_;}
	virtual ~Basis_SinglaEVFunctional(){}

};
typedef boost::shared_ptr<Basis_SinglaEVFunctional> Basis_SinglaEVFunctional_PTR;
typedef boost::shared_ptr<const Basis_SinglaEVFunctional> Basis_SinglaEVFunctional_CONSTPTR;



class Basis_ConstUnity : public Basis_SinglaEVFunctional  // consider 1 as a basis 
{
	class ConstUnityTransformer : public Transformer
	{
	public:
		virtual double transform(double val)const{return 0.01;}

		void write_to_stream(std::ostream& out)const;
	};

public:
	Basis_ConstUnity() : Basis_SinglaEVFunctional(EV_CONSTPTR(new EV_ConstRate(Rate1_CONSTPTR(new ConstRate(0.01)))), 
													Transformer_CONSTPTR(new ConstUnityTransformer())){} 

	typedef boost::shared_ptr<ConstUnityTransformer> ConstUnityTransformer_PTR;
	typedef boost::shared_ptr<const ConstUnityTransformer> ConstUnityTransformer_CONSTPTR;

	virtual ~Basis_ConstUnity(){};

	void write_to_stream(std::ostream& out)const;
};
typedef boost::shared_ptr<Basis_ConstUnity> Basis_ConstUnity_PTR;
typedef boost::shared_ptr<const Basis_ConstUnity> Basis_ConstUnity_CONSTPTR;



class Basis_CappedFlooredCoupon: public Basis_SinglaEVFunctional
{
public:
	class CappedFlooredCouponTransformer : public Transformer
	{
		CappedFlooredCoupon_CONSTPTR cappedFlooredCoupon_;
	public:
		CappedFlooredCouponTransformer(CappedFlooredCoupon_CONSTPTR cappedFlooredCoupon): cappedFlooredCoupon_(cappedFlooredCoupon){};
		virtual double transform(double val)const;

		void write_to_stream(std::ostream& out)const;
	};

	typedef boost::shared_ptr<CappedFlooredCouponTransformer> CappedFlooredCouponTransformer_PTR;
	typedef boost::shared_ptr<const CappedFlooredCouponTransformer> CappedFlooredCouponTransformer_CONSTPTR;
	
	Basis_CappedFlooredCoupon(EV_CONSTPTR ev, Transformer_CONSTPTR transformer):Basis_SinglaEVFunctional(ev,transformer){}
	virtual ~Basis_CappedFlooredCoupon(){}

	void write_to_stream(std::ostream& out)const;
};
typedef boost::shared_ptr<Basis_CappedFlooredCoupon> Basis_CappedFlooredCoupon_PTR;
typedef boost::shared_ptr<const Basis_CappedFlooredCoupon> Basis_CappedFlooredCoupon_CONSTPTR;



class Basis_Monomial : public Basis_SinglaEVFunctional
{
public:
	class MonomialTransformer : public Transformer
	{
        double coeff_;
		size_t order_;  
	public:
		MonomialTransformer(double coeff, size_t order): coeff_(coeff), order_(order){};
		virtual double transform(double val)const {return coeff_*std::pow(val,order_);}

		void write_to_stream(std::ostream& out)const;
	};

	typedef boost::shared_ptr<MonomialTransformer> MonomialTransformer_PTR;
	typedef boost::shared_ptr<const MonomialTransformer> MonomialTransformer_CONSTPTR;

	Basis_Monomial(EV_CONSTPTR ev, Transformer_CONSTPTR transformer):Basis_SinglaEVFunctional(ev,transformer){}
	virtual ~Basis_Monomial(){}

	void write_to_stream(std::ostream& out)const;
};
typedef boost::shared_ptr<Basis_Monomial> Basis_Monomial_PTR;
typedef boost::shared_ptr<const Basis_Monomial> Basis_Monomial_CONSTPTR;



// ----------------------------------------------------------------------
// 
//                Basis_Composited
//
// ----------------------------------------------------------------------
//class Basis_Composited : public Basis   // composit tow basis by operator: +,-,*,/ 
//{
//
//public:
//	class Transformer
//	{
//	public: 
//		virtual double transform(double val)const=0;  // val can be ev_val or basis_val ! 
//		//virtual void write_to_stream(std::ostream& out)const=0;
//		virtual ~Transformer(){}
//	};
//	typedef boost::shared_ptr<Transformer> Transformer_PTR;	
//	typedef boost::shared_ptr<const Transformer> Transformer_CONSTPTR;
//
//private: 
//
//	Basis_CONSTPTR basis1_;
//	Basis_CONSTPTR basis2_;
//
//	Transformer_CONSTPTR transformer_;
//
//public:
//	Basis_Composited(Basis_CONSTPTR basis1, Basis_CONSTPTR basis2, Transformer_CONSTPTR transformer);
//	virtual ~Basis_Composited(){}
//
//	Basis_CONSTPTR getBasis1() const {return basis1_;}
//	Basis_CONSTPTR getBasis2() const {return basis2_;}
//	Transformer_CONSTPTR getTransformer() const {return transformer_;}
//
//};
//typedef boost::shared_ptr<Basis_Composited> Basis_Composited_PTR;	
//typedef boost::shared_ptr<const Basis_Composited> Basis_Composited_CONSTPTR;




class Basis_Composited : public Basis  // composit tow basis by operator: +,-,*,/ 
{
public:
	enum class Compositor {ADD, SUB, MUL, DIV}; 
public:
	class Transformer
	{
		Compositor compositor_;
	public:
		Transformer(const Compositor& compositor):compositor_(compositor){}
	    virtual double transform(double val1, double val2)const 
		{
			switch(compositor_)
			{
				case Compositor::ADD:
					return val1 + val2;
				case Compositor::SUB:
					return val1 - val2;
				case Compositor::MUL:
					return val1 * val2;
				case Compositor::DIV:
					assert(val2!=0.0);
					return val1 / val2;
				default:
					throw("class Compositor operator not defined! ");
			}			
		}
	};
	typedef boost::shared_ptr<Transformer> Transformer_PTR;	
	typedef boost::shared_ptr<const Transformer> Transformer_CONSTPTR;

private: 

	Basis_CONSTPTR basis1_;
	Basis_CONSTPTR basis2_;

	Transformer_CONSTPTR transformer_;


public:
	Basis_Composited(Basis_CONSTPTR basis1, Basis_CONSTPTR basis2, const Compositor& compositor);
	virtual ~Basis_Composited(){}

	Basis_CONSTPTR getBasis1() const {return basis1_;}
	Basis_CONSTPTR getBasis2() const {return basis2_;}
	Transformer_CONSTPTR getTransformer() const {return transformer_;}

};
typedef boost::shared_ptr<Basis_Composited> Basis_Composited_PTR;	
typedef boost::shared_ptr<const Basis_Composited> Basis_Composited_CONSTPTR;



//Junbin: create the polynomial
class Basis_Polynomial : public Basis 
{
	//TODO latter ... 
public:
	class Transformer
	{
		double coeff_;
        std::vector<size_t> power_vect_;
	public:
		Transformer(double coeff, const std::vector<size_t>& power_vect);
		virtual double transform(const std::vector<double>& val_vect)const;

		virtual ~Transformer(){}

		void write_to_stream(std::ostream& out)const;
	};
	typedef boost::shared_ptr<Transformer> Transformer_PTR;	
	typedef boost::shared_ptr<const Transformer> Transformer_CONSTPTR;

private:
	std::vector<Basis_CONSTPTR> basis_vect_;
	Transformer_CONSTPTR transformer_;

public:
	//contructor and destructor
	Basis_Polynomial(const std::vector<Basis_CONSTPTR>& basis_vect, double coeff , const std::vector<size_t>& power_vect);
	virtual ~Basis_Polynomial(){}

	//getter
	const std::vector<Basis_CONSTPTR>& get_Basis_Vect()const{return basis_vect_;}
	Transformer_CONSTPTR getTransformer()const{return transformer_;}

	void write_to_stream(std::ostream& out)const;
};
typedef boost::shared_ptr<Basis_Polynomial> Basis_Polynomial_PTR;	
typedef boost::shared_ptr<const Basis_Polynomial> Basis_Polynomial_CONSTPTR;