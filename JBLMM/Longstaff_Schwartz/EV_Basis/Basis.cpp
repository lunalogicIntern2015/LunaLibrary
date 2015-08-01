#include <JBLMM/Longstaff_Schwartz/EV_Basis/Basis.h>


void Basis::write_to_stream(std::ostream& out)const
{
}

Basis_Composited::Basis_Composited(Basis_CONSTPTR basis1, Basis_CONSTPTR basis2, const Compositor& compositor)
		:
			basis1_(basis1), 
			basis2_(basis2), 
			transformer_(new Transformer(compositor))
{
}


double Basis_CappedFlooredCoupon::CappedFlooredCouponTransformer::transform(double val)const
{
	if(cappedFlooredCoupon_->getIfFloored())
		return std::max(cappedFlooredCoupon_->getFloorStrike(), cappedFlooredCoupon_->getMultiFactor()*val + cappedFlooredCoupon_->getAddFactor());


	double nominal						=	cappedFlooredCoupon_->getNominal();
	double period						=	cappedFlooredCoupon_->getPeriod();
	double floor						=	cappedFlooredCoupon_->getFloorStrike();
	double cap							=	cappedFlooredCoupon_->getCapStrike();
	double multiFactor					=	cappedFlooredCoupon_->getMultiFactor();
	double addFactor					=	cappedFlooredCoupon_->getAddFactor();

	double result =nominal*period*std::max(floor, std::min(cap, multiFactor*val + addFactor));

	return result;
}



void Basis_ConstUnity::write_to_stream(std::ostream& out)const
{
	out	<<	"---Basis_ConstUnity---"	<<	std::endl;
	ev_->write_to_stream(out);
	transformer_->write_to_stream(out);
}

void Basis_ConstUnity::ConstUnityTransformer::write_to_stream(std::ostream& out)const
{
	out	<<	"---ConstUnityTransformer---"	<<	std::endl;
	out	<<	"constant ; "	<<	1.0  << std::endl;
}

void Basis_CappedFlooredCoupon::write_to_stream(std::ostream& out)const
{
	out	<<	"---Basis_CappedFlooredCoupon---"	<<	std::endl;
	ev_->write_to_stream(out);
	transformer_->write_to_stream(out);
}

void Basis_CappedFlooredCoupon::CappedFlooredCouponTransformer::write_to_stream(std::ostream& out)const
{
	out	<<	"---CappedFlooredCouponTransformer---"	<<	std::endl;
	cappedFlooredCoupon_->write_to_stream(out);
}


void Basis_Monomial::write_to_stream(std::ostream& out)const
{
	out	<<	"---Basis_Monomial---"	<<	std::endl;
	ev_->write_to_stream(out);
	transformer_->write_to_stream(out);
}

void Basis_Monomial::MonomialTransformer::write_to_stream(std::ostream& out)const
{
	out	<<	"---MonomialTransformer---"	<<	std::endl;
	out	<<	"coeff : ;"	<<	coeff_ <<	std::endl;
	out	<<	"order : ;"	<<	order_ <<	std::endl;
}

Basis_Polynomial::Basis_Polynomial(const std::vector<Basis_CONSTPTR>& basis_vect, double coeff , const std::vector<size_t>& power_vect)
	:
	basis_vect_(basis_vect),
	transformer_(new Transformer(coeff,power_vect))
{
}

Basis_Polynomial::Transformer::Transformer(double coeff, const std::vector<size_t>& power_vect)
	:
	coeff_(coeff),
	power_vect_(power_vect)
{
}

double Basis_Polynomial::Transformer::transform(const std::vector<double>& val_vect)const
{
	assert(val_vect.size()==power_vect_.size());
	double result = coeff_;
	for(size_t i = 0; i<val_vect.size(); i++)
	{
		if(power_vect_[i]!=0)
			result *=std::pow(val_vect[i], power_vect_[i]);
	}

	return result;
}

void Basis_Polynomial::write_to_stream(std::ostream& out)const
{
	out	<<	"---Basis_Polynomial---"	<<	std::endl;
	for(size_t i = 0; i<basis_vect_.size(); i++)
		basis_vect_[i]->write_to_stream(out);
}