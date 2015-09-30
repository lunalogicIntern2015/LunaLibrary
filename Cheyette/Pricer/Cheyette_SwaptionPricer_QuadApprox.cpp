#include <Cheyette/Pricer/Cheyette_SwaptionPricer_QuadApprox.h>

const size_t gridSize  = 101 ;	//mettre 11 ou 101 (11 est suffisant)


//void Cheyette_SwaptionPricer_QuadApprox::initialize_buffers() 
//{
//	Cheyette_SwaptionPricer_Approx::initialize_buffers() ;
//	
//	//buffer_b_barre_ = timeAverage(buffer_T0_) ;
//	//buffer_c_barre_ = timeAverage(buffer_T0_) ;
//}

void Cheyette_SwaptionPricer_QuadApprox::preCalculateALL(VanillaSwaption_PTR pSwaption) const
{
	Cheyette_SwaptionPricer_Approx::preCalculateALL(pSwaption) ;
	buffer_b_barre_ = timeAverage_b_S(buffer_T0_) ;
	buffer_c_barre_ = timeAverage_c_S(buffer_T0_) ;
}

void Cheyette_SwaptionPricer_QuadApprox::preCalculateModelBuffers() const
{
	buffer_b_barre_ = timeAverage_b_S(buffer_T0_) ;
	buffer_c_barre_ = timeAverage_c_S(buffer_T0_) ;
}

//derivees 3eme wrt x_t

double Cheyette_SwaptionPricer_QuadApprox::ZC_3rdDerivative_on_xt(double t, double T, double x_t, double y_t) const
{
	assert(0 <= T && T <= buffer_TN_) ;

	double ZC = pCheyette_Model_->P(t, T, x_t, y_t) ;
	double g = pCheyette_Model_->G(t,T) ;
	double res = - g * g * g * ZC ;
	
	return res ;					
}

double Cheyette_SwaptionPricer_QuadApprox::swapRateNumerator_3rdDerivative(double t, double x_t, double y_t) const
{
	return ZC_3rdDerivative_on_xt(t, buffer_T0_, x_t, y_t) - ZC_3rdDerivative_on_xt(t, buffer_TN_, x_t, y_t) ;
}

double Cheyette_SwaptionPricer_QuadApprox::swapRateDenominator_3rdDerivative(double t, double x_t, double y_t) const
{
	double fixed_tenor = buffer_UnderlyingSwap_.get_fixedLegTenorType().YearFraction() ;
	double float_tenor = buffer_UnderlyingSwap_.get_floatingLegTenorType().YearFraction() ;
	double tenor_ref = std::min(fixed_tenor, float_tenor) ;  

	double result = 0. ;
	for(size_t itr = 0; itr < buffer_fixedLegPaymentIndexSchedule_.size(); ++itr) 
	{
		double dateEchangeFluxFixe = buffer_fixedLegPaymentIndexSchedule_[itr] * tenor_ref ;
		result += buffer_deltaTFixedLeg_[itr] * ZC_3rdDerivative_on_xt(t, dateEchangeFluxFixe, x_t, y_t) ;	
	}
	return result;
}

double Cheyette_SwaptionPricer_QuadApprox::swapRate_3rdDerivative(double t, double x_t, double y_t) const
{
	double n	= swapRateNumerator(t, x_t, y_t);
	double n_1	= swapRateNumerator_1stDerivative(t, x_t, y_t); 
	double n_2	= swapRateNumerator_2ndDerivative(t, x_t, y_t); 
	double n_3	= swapRateNumerator_3rdDerivative(t, x_t, y_t); 

	double d	= swapRateDenominator(t, x_t, y_t);
	double d_1	= swapRateDenominator_1stDerivative(t, x_t, y_t);
	double d_2	= swapRateDenominator_2ndDerivative(t, x_t, y_t);
	double d_3	= swapRateDenominator_3rdDerivative(t, x_t, y_t);

	double result = (n_3 * d + n_2 * d_1 - n_1 * d_2 - n * d_3) * d + (n_2 * d - n * d_2) * (1 - 2 * d_1)
						- (n_1 * d - n * d_1) * 2 * d_2 ;
	result /= pow(d, 6);

	return result;
}

/**************** Quadratic approx of swap rate volatility ********************
*
*	dS(t) = phi(t, S_t) dW^Q_A
*	dS(t) \approx [ phi(t, s_bar) + dPhi(t, s_bar) (S_t - s_bar) + 1./2. d2Phi(t, s_bar) (S_t - s_bar)^2 ] dW^Q_A
*
******************************************************************************/

//DL(2) de Phi(t, s) en s_barre = S0
double Cheyette_SwaptionPricer_QuadApprox::DL2_Phi_t_s(double t, double s) const
{
	double inverse_s_bar =	inverse(t, buffer_s0_) ; 
	return phi(t, inverse_s_bar )	+ dPhi(t, inverse_s_bar ) * (s - buffer_s0_) 
									+ 0.5 * d2Phi(t, inverse_s_bar) * pow(s - buffer_s0_, 2) ;
}

double Cheyette_SwaptionPricer_QuadApprox::phi(double t, double x_t) const
{
	//retourne dS(t)/dx(t) sigma_r(t) (t=0, inverse, y_bar)
	double y_bar = buffer_y_bar_(t) ;

	double sigma_r	= pCheyette_Model_->localVol(t, x_t, y_bar) ; 
	double phi_res = swapRate_1stDerivative(t, x_t, y_bar) * sigma_r ;
	
	return phi_res ;

}

// d Phi / ds
double Cheyette_SwaptionPricer_QuadApprox::dPhi(double t, double x_t) const
{
	double y_bar = buffer_y_bar_(t) ;

	double sigma_r	= pCheyette_Model_->localVol(t, x_t, y_bar) ; 
	double res =  sigma_r *	swapRate_2ndDerivative(t, x_t, y_bar) / swapRate_1stDerivative(t, x_t, y_bar) 
					+ pCheyette_Model_->localVol_1stDerivative(t, x_t, y_bar) ;

	return res ; 
}

// d2 Phi / ds
double Cheyette_SwaptionPricer_QuadApprox::d2Phi(double t, double x_t) const
{
	double y_bar = buffer_y_bar_(t) ;

	double sigma_r		= pCheyette_Model_->localVol(t, x_t, 0.) ; 
	double sigma_r_1	= pCheyette_Model_->localVol_1stDerivative(t, x_t, 0.) ;
	double sigma_r_2	= pCheyette_Model_->localVol_2ndDerivative(t, x_t, 0.) ;
	double S_1	= swapRate_1stDerivative(t, x_t, y_bar); 
	double S_2	= swapRate_2ndDerivative(t, x_t, y_bar); 
	double S_3	= swapRate_3rdDerivative(t, x_t, y_bar); 

	double res = sigma_r * ( S_3/(S_1 * S_1) - pow(S_2, 2) / pow(S_1, 3) )
					+ sigma_r_1 * S_2 / (S_1*S_1) 
					+ sigma_r_2 / S_1 ;
		 
	return res ; 
}

/****************  parameter averaging  *********************
*
*	dS(t) \approx [ phi(t, s_bar) + dPhi(t, s_bar) (S_t - s_bar) + 1./2. d2Phi(t, s_bar) (S_t - s_bar)^2 ] dW^Q_A
*	dS(t) = \lambda(t) ( S0 + b(t) (S(t) - S(0)) + 1./2. c(t) (S(t) - S(0))^2 ) dW^Q_A		
*
*	-> time averaging:  b(t) -> b_barre
*						c(t) -> c_barre
************************************************************/

/* dS(t) = \lambda(t) ( S0 + b_S(t) (S_t - s_bar) + 1./2. c_S(t) (S_t - s_bar)^2 ) dW_t QA */

//en prenant s0 = s_bar
//lambda(t, s_bar)
double Cheyette_SwaptionPricer_QuadApprox::lambda(double t) const
{
	double inverse_s_bar =	inverse(t, buffer_s0_) ; 
	return phi(t, inverse_s_bar) / buffer_s0_ ;
}

//b_S (t, s_bar)
double Cheyette_SwaptionPricer_QuadApprox::b_S(double t) const 
{
	double inverse_s_bar =	inverse(t, buffer_s0_) ; 
	return dPhi(t, inverse_s_bar) / lambda(t) ;
}

//c_S (t, s_bar)
double Cheyette_SwaptionPricer_QuadApprox::c_S(double t) const
{
	double inverse_s_bar =	inverse(t, buffer_s0_) ; 
	return d2Phi(t, inverse_s_bar) / lambda(t) ;
}

//autres fonctions intermédiaires pour les time averaging
double Cheyette_SwaptionPricer_QuadApprox::lambda2(double t) const
{
	return lambda(t) * lambda(t) ; 
}

double Cheyette_SwaptionPricer_QuadApprox::lambda4(double t) const
{
	return pow(lambda(t), 4) ; 
}

double Cheyette_SwaptionPricer_QuadApprox::f_outer_num_b(double t) const
{
	return lambda2(t) * b_S(t) ;  //OK (fonction outer seulement ie sans v^2(u)
}
double Cheyette_SwaptionPricer_QuadApprox::f_outer_denom(double t) const
{
	return lambda2(t)  ;		//OK (fonction outer seulement ie sans v^2(u)
}

double Cheyette_SwaptionPricer_QuadApprox::f_outer_num_c(double t) const
{
	return lambda2(t) * c_S(t) ;  //OK (fonction outer seulement ie sans v^2(u)
}

//time averaging of b_S(t) and c_S(t)
double Cheyette_SwaptionPricer_QuadApprox::timeAverage_b_S(double t) const
{
	double gridStart = 0.0;
	double gridEnd = t ; 

//integrale numerateur
	numeric::IncrementalIntegrator2D_Riemann int2D(gridStart, gridEnd, gridSize) ;
	boost::function<double(double)> f_inner = boost::bind(&Cheyette_SwaptionPricer_QuadApprox::lambda2, boost::ref(*this), _1);
	boost::function<double(double)> f_outer = boost::bind(&Cheyette_SwaptionPricer_QuadApprox::f_outer_num_b, boost::ref(*this), _1);	
	double integrale_numerateur = int2D.integrate(f_outer, f_inner) ;
	
//integrale denominateur
	boost::function<double(double)> f_denom = boost::bind(&Cheyette_SwaptionPricer_QuadApprox::f_outer_denom, boost::ref(*this), _1);
	double integrale_denom = int2D.integrate(f_denom, f_inner);

	double b_barre = integrale_numerateur/integrale_denom ;

	return b_barre ;
}

double Cheyette_SwaptionPricer_QuadApprox::timeAverage_c_S(double t) const
{
	double gridStart = 0.0;
	double gridEnd = t ; 

//integrale numerateur
	numeric::IncrementalIntegrator2D_Riemann int2D(gridStart, gridEnd, gridSize) ;
	boost::function<double(double)> f_inner = boost::bind(&Cheyette_SwaptionPricer_QuadApprox::lambda4, boost::ref(*this), _1);
	boost::function<double(double)> f_outer = boost::bind(&Cheyette_SwaptionPricer_QuadApprox::f_outer_num_c, boost::ref(*this), _1);	
	double integrale_numerateur = int2D.integrate(f_outer, f_inner) ;
	
//integrale denominateur
	boost::function<double(double)> f_denom = boost::bind(&Cheyette_SwaptionPricer_QuadApprox::f_outer_denom, boost::ref(*this), _1);
	double integrale_denom = int2D.integrate(f_denom, f_inner);

	double c_barre = integrale_numerateur/integrale_denom ;

	return c_barre ;
}


double Cheyette_SwaptionPricer_QuadApprox::Phi_St(double St) const
{
	return buffer_s0_ + buffer_b_barre_ * (St - buffer_s0_) + 1./2. * buffer_c_barre_ * pow(St - buffer_s0_, 2.) ;
//	return sqrt(St) ;
//	return St ;
}

double Cheyette_SwaptionPricer_QuadApprox::Phi_St_prime(double St) const
{
	return buffer_b_barre_ + buffer_c_barre_ * (St - buffer_s0_) ;
//	return 1. / (2. * sqrt(St)) ;
//	return 1. ;
}

double Cheyette_SwaptionPricer_QuadApprox::Phi_St_seconde(double St) const
{
	return buffer_c_barre_ ;
//	return - 1./ 4. * pow(St, - 3./2.) ;
//	return 0. ;
}

// 1 / Phi(u)
double Cheyette_SwaptionPricer_QuadApprox::UnSurPhi_St(double St) const
{
	return 1. / Phi_St(St) ;
}

// Phi^2(u)
double Cheyette_SwaptionPricer_QuadApprox::Phi2_St(double St) const
{
	return pow(Phi_St(St), 2.) ;
}

/**************************  Time change ****************************/

double Cheyette_SwaptionPricer_QuadApprox::timeChange(double t) const
{
	double start(0.), end(t) ;
	size_t nbPoints(static_cast<size_t>(252 * t)) ;
	numeric::Integrator1D_Riemann inte(start, end, nbPoints);

	boost::function<double(double)> f = boost::bind(&Cheyette_SwaptionPricer_QuadApprox::lambda2, boost::ref(*this), _1);   
	//lambda2 à verifier 
	double integrale = inte.integrate(f) ;

	return integrale ;
}


/**********************  Call approximation ************************
** dynamique : 
**
**  dS(t) = \Phi(S_t) dW(t)
********************************************************************/

//approximation call

double Cheyette_SwaptionPricer_QuadApprox::omega(double t, double S) const 
{
	return pow(t, 1./2.) * omega0(S) + pow(t, 3./2.) * omega1(S) ; 
}

double Cheyette_SwaptionPricer_QuadApprox::omega0(double S) const 
{
	double K = buffer_UnderlyingSwap_.get_strike() ;

	double m = Phi_St_prime(buffer_s0_) ;
	double p = Phi_St(buffer_s0_) - buffer_s0_ * m ;  
	double alpha = p / m ;

	double start(K), end(S) ;
	size_t nbPoints(100) ;			//amelioration : faire dependre nb de points de longueur de l'intervalle
	numeric::Integrator1D_Riemann inte(start, end, nbPoints);
	boost::function<double(double)> f = boost::bind(&Cheyette_SwaptionPricer_QuadApprox::UnSurPhi_St, boost::ref(*this), _1);
	double integrale = inte.integrate(f) ;

	return log((S + alpha) / (K + alpha)) / integrale ; 	
}

double Cheyette_SwaptionPricer_QuadApprox::omega1(double S) const 
{
	double K = buffer_UnderlyingSwap_.get_strike() ;

	double m = Phi_St_prime(buffer_s0_) ;
	double p = Phi_St(buffer_s0_) - buffer_s0_ * m ;  
	double alpha = p / m ;

	double start(K), end(S) ;
	size_t nbPoints(100) ;		//amelioration : faire dependre nb de points de longueur de l'intervalle
	numeric::Integrator1D_Riemann inte(start, end, nbPoints);
	boost::function<double(double)> f = boost::bind(&Cheyette_SwaptionPricer_QuadApprox::UnSurPhi_St, boost::ref(*this), _1);
	double integrale = inte.integrate(f) ;

	double fraction = (S + alpha) * (K + alpha) / ( Phi_St(S) * Phi_St(K) ) ;
	double omega0_S = omega0(S) ;
	double Log = log(omega0_S * sqrt(fraction) );
	
	return - omega0_S / (integrale * integrale) * Log ;
}


double Cheyette_SwaptionPricer_QuadApprox::omega_ATM(double t) const 
{
	double K = buffer_UnderlyingSwap_.get_strike() ;

	std::cout << "omega0 (vs 1) : " << omega0_ATM() << " , omega 1 : " << omega1_ATM() << std::endl ;  
	return pow(t, 1./2.) * omega0_ATM() + pow(t, 3./2.) * omega1_ATM() ; 
}


double Cheyette_SwaptionPricer_QuadApprox::omega0_ATM() const 
{
	double K = buffer_UnderlyingSwap_.get_strike() ;

	double m = Phi_St_prime(buffer_s0_) ;
	double p = Phi_St(buffer_s0_) - buffer_s0_ * m ;  
	double alpha = p / m ;

	return Phi_St(K) /(K + alpha) ;

}

double Cheyette_SwaptionPricer_QuadApprox::omega1_ATM() const 
{
	double K = buffer_UnderlyingSwap_.get_strike() ;
	std::cout << "omega1 ATM, valeur de K ATM : " << K << std::endl ;

	double m = Phi_St_prime(buffer_s0_) ;
	double p = Phi_St(buffer_s0_) - buffer_s0_ * m ;  
	double alpha = p / m ;

	return 1. / 24.  * pow(omega0_ATM(), 3) * 
		( 1. + (pow(K+alpha, 2.) / pow(Phi_St(K), 2.)) * (2. * Phi_St(K) * Phi_St_seconde(K) - pow(Phi_St_prime(K), 2.)) ) ;
}

//y = mx + p du proxy model
double Cheyette_SwaptionPricer_QuadApprox::price(const VanillaSwaption_PTR pSwaption) const 
{
	preCalculateALL(pSwaption) ;

	double K = buffer_UnderlyingSwap_.get_strike() ;
	double T = timeChange(buffer_T0_) ;  // pour le time change, à vérifier, ou buffer_T0_
	double S0 = buffer_s0_ ;

	double m = Phi_St_prime(buffer_s0_) ;
	double p = Phi_St(buffer_s0_) - buffer_s0_ * m ;  
	double alpha = p / m ;

	double omega_value = (S0==K)? omega_ATM(T) : omega(T, S0) ;
	
	double dPlus =  (log((S0 + alpha)/(K + alpha)) + 0.5* pow(omega_value, 2) )/ omega_value ;
	double dMoins = dPlus - m * sqrt(T);

	boost::math::normal_distribution<> nd(0,1); 
	double N1 = cdf(nd,dPlus);
	double N2 = cdf(nd,dMoins); 

	double prixApprox = (S0 + alpha) * N1 - (K + alpha) * N2 ;
	
	double annuite = swapRateDenominator(0., 0., 0.) ;
	
	return annuite * prixApprox ;
}
