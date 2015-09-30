#include <Cheyette/Model/CheyetteModel.h>

//en doublon avec la classe Swaption Pricer Approx

//pour la version SwapRate de Displaced Diffusion
double CheyetteModel::swapRateNumerator(double t, double x_t, double y_t, const VanillaSwaption_PTR pSwaption) const 
{
	assert( t >= 0. ) ;

	double buffer_T0_	= pSwaption->getUnderlyingSwap().get_StartDate() ;	//1ere date de fixing
	double buffer_TN_	= pSwaption->getUnderlyingSwap().get_EndDate() ;	//date dernier flux (setté en T_{N-1})

	double ZC_T0, ZC_TN ;
	ZC_T0 = P(t, buffer_T0_, x_t, y_t) ; 
	ZC_TN = P(t, buffer_TN_, x_t, y_t) ;	
	
	return  ZC_T0 - ZC_TN ; 
}

double CheyetteModel::swapRateDenominator(double t, double x_t, double y_t, const VanillaSwaption_PTR pSwaption) const 
{	
	VanillaSwap			buffer_UnderlyingSwap_					= pSwaption->getUnderlyingSwap() ;
	std::vector<size_t> buffer_fixedLegPaymentIndexSchedule_	= buffer_UnderlyingSwap_.get_fixedLegPaymentIndexSchedule() ;
	std::vector<double> buffer_deltaTFixedLeg_					= buffer_UnderlyingSwap_.get_DeltaTFixedLeg() ;

	double fixed_tenor = buffer_UnderlyingSwap_.get_fixedLegTenorType().YearFraction() ;
	double float_tenor = buffer_UnderlyingSwap_.get_floatingLegTenorType().YearFraction() ;
	double tenor_ref = std::min(fixed_tenor, float_tenor) ;  //le plus petit 

	double denom = 0.0;
	for(size_t itr = 0 ; itr < buffer_fixedLegPaymentIndexSchedule_.size() ; ++itr) 
	{
		double dateEchangeFluxFixe = buffer_fixedLegPaymentIndexSchedule_[itr] * tenor_ref ;
		double delta_T = buffer_deltaTFixedLeg_[itr] ;	
		double ZC		= P(t, dateEchangeFluxFixe, x_t, y_t) ;
		denom += delta_T * ZC ;
	}
	
	return denom ;
}

double CheyetteModel::swapRate(double t, double x_t, double y_t, const VanillaSwaption_PTR pSwaption) const
{
	double n      = swapRateNumerator(t, x_t, y_t, pSwaption);
	double d      = swapRateDenominator(t, x_t, y_t, pSwaption);

	assert (d > 0.001) ;
	return n/d;
}

double CheyetteModel::swapRate0(const VanillaSwaption_PTR pSwaption) const
{
	return swapRate(0., 0., 0., pSwaption);
}

double CheyetteModel::ZC_1stDerivative_on_xt(double t, double T, double x_t, double y_t,
											 const VanillaSwaption_PTR pSwaption) const
{
	double buffer_TN_	= pSwaption->getUnderlyingSwap().get_EndDate() ;	//date dernier flux (setté en T_{N-1})
	assert(0 <= T && T <= buffer_TN_) ;
	if (t==0){assert(x_t == 0 && y_t == 0) ;}

	double ZC = P(t, T, x_t, y_t) ;	
	double g = G(t,T) ;
	double res = - g * ZC ;
	return res ;	
}

double CheyetteModel::swapRateNumerator_1stDerivative(double t, double x_t, double y_t,
													  const VanillaSwaption_PTR pSwaption) const
{
	double buffer_T0_	= pSwaption->getUnderlyingSwap().get_StartDate() ;	//1ere date de fixing
	double buffer_TN_	= pSwaption->getUnderlyingSwap().get_EndDate() ;	//date dernier flux (setté en T_{N-1})
	return ZC_1stDerivative_on_xt(t, buffer_T0_, x_t, y_t, pSwaption) 
					- ZC_1stDerivative_on_xt(t, buffer_TN_, x_t, y_t, pSwaption) ;
}

double CheyetteModel::swapRateDenominator_1stDerivative(double t, double x_t, double y_t,
													  const VanillaSwaption_PTR pSwaption) const
{
	double result = 0. ;

	VanillaSwap			buffer_UnderlyingSwap_					= pSwaption->getUnderlyingSwap() ;
	std::vector<size_t> buffer_fixedLegPaymentIndexSchedule_	= buffer_UnderlyingSwap_.get_fixedLegPaymentIndexSchedule() ;
	std::vector<double> buffer_deltaTFixedLeg_					= buffer_UnderlyingSwap_.get_DeltaTFixedLeg() ;

	double fixed_tenor = buffer_UnderlyingSwap_.get_fixedLegTenorType().YearFraction() ;
	double float_tenor = buffer_UnderlyingSwap_.get_floatingLegTenorType().YearFraction() ;
	double tenor_ref = std::min(fixed_tenor, float_tenor) ;  

	for(size_t itr = 0 ; itr < buffer_fixedLegPaymentIndexSchedule_.size() ; ++itr) 
	{
		double dateEchangeFluxFixe = buffer_fixedLegPaymentIndexSchedule_[itr] * tenor_ref ;
		result += buffer_deltaTFixedLeg_[itr] * ZC_1stDerivative_on_xt(t, dateEchangeFluxFixe, x_t, y_t, pSwaption) ;	
	}
	return result;
}

double CheyetteModel::swapRate_1stDerivative(double t, double x_t, double y_t,
											  const VanillaSwaption_PTR pSwaption) const
{
	double n   = swapRateNumerator(t, x_t, y_t, pSwaption);
	double n_1 = swapRateNumerator_1stDerivative(t, x_t, y_t, pSwaption); 

	double d   = swapRateDenominator(t, x_t, y_t, pSwaption);
	double d_1 = swapRateDenominator_1stDerivative(t, x_t, y_t, pSwaption);

	return  (n_1*d - n*d_1)/(d*d) ;	
}
