//#include <JBLMM/Longstaff_Schwartz/Basis.h>
//
//double CappedFlooredCouponTransformer::evaluate_x(double x)const
//{
//	double nominal						=	cappedFlooredCoupon_->getNominal();
//	double period						=	cappedFlooredCoupon_->getPeriod();
//	double floor						=	cappedFlooredCoupon_->getFloorStrike();
//	double cap							=	cappedFlooredCoupon_->getCapStrike();
//	double multiFactor					=	cappedFlooredCoupon_->getMultiFactor();
//	double addFactor					=	cappedFlooredCoupon_->getAddFactor();
//
//	double result =nominal*period*std::max(floor, std::min(cap, multiFactor*x+addFactor));
//
//	return result;
//}
//
//double Polynomial::evaluate_x(double x)const
//{
//	double res=0.0;
//	for(size_t i = 0; i<coef_.size(); i++)
//	{
//		res+=std::pow(x,coef_[i]);
//	}
//	return res;
//}
//
//
////double Polynomial::evaluate_vect(const std::vector<double>& variableVect)const
////{
////	size_t nbVariableInput = variableVect.size();
////	assert(nbVariableInput==nbVariable_);
////	size_t nbMonomial = coefandPower_.size();
////	double result = 0.0;
////	for(size_t indexMonomial=0;indexMonomial<nbMonomial;indexMonomial++)
////	{
////		
////		double monomialValue=coefandPower_[indexMonomial].first;
////		const std::vector<size_t>& powerVect = coefandPower_[indexMonomial].second;
////		for(size_t variableIndex = 0; variableIndex<nbVariable_; variableIndex++)
////		{
////			size_t power = powerVect[variableIndex];
////			double variableValue = variableVect[variableIndex];
////			monomialValue *= pow(variableValue, power);
////		}
////		result += monomialValue;
////	}
////	return result;
////}
//
////bool  Polynomial::checkVariableDim()const
////{
////	for(size_t indexMonomial=0;indexMonomial<coefandPower_.size();indexMonomial++)
////	{
////		if(coefandPower_[indexMonomial].second.size()!=nbVariable_)
////			return false;
////	}
////	return true;
////}
//
//double SingleEvBasis::evaluateEV(McLmm_PTR mclmm, BermudanSwaption_CONSTPTR bermudanSwaption_ptr, LMM::Index liborIndex )const
//{
//	double ev_value=ev_ptr_->getEvaluator_PTR()->evaluate(mclmm,bermudanSwaption_ptr,liborIndex);
//	return ev_value;
//}
//
//double ComposedEvBasis::evaluateEV(McLmm_PTR mclmm, BermudanSwaption_CONSTPTR bermudanSwaption_ptr, LMM::Index liborIndex )const
//{
//	double ev_value=composedEV_ptr_->evaluate(mclmm,bermudanSwaption_ptr,liborIndex);
//	return ev_value;
//}