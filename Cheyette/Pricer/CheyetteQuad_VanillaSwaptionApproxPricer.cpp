#include "CheyetteQuad_VanillaSwaptionApproxPricer.h"


const size_t gridSize  = 101 ;


CheyetteQuad_VanillaSwaptionApproxPricer::CheyetteQuad_VanillaSwaptionApproxPricer
						(const CheyetteQuad_Model_PTR&	pCheyetteQuad_Model, 
						const VanillaSwaption_PTR&		pSwaption) 
			: pCheyetteQuad_Model_(pCheyetteQuad_Model), pSwaption_(pSwaption)  					
{
	// domaine de validité de l'approximation : //assert( ... ); 

	pBuffer_courbeInput_	 = pCheyetteQuad_Model->get_courbeInput_PTR() ;
	initialize_buffers() ;  //initialise les autres buffers relatifs à la swaption, y barre, b barre 		
}

void CheyetteQuad_VanillaSwaptionApproxPricer::initialize_buffers()
{
	buffer_UnderlyingSwap_					= pSwaption_->getUnderlyingSwap() ;
	buffer_T0_								= pSwaption_->getUnderlyingSwap().get_StartDate() ;		//1ere date de fixing
	buffer_TN_								= pSwaption_->getUnderlyingSwap().get_EndDate() ;			//date dernier flux (setté en T_{N-1})
	buffer_fixedLegPaymentIndexSchedule_	= buffer_UnderlyingSwap_.get_fixedLegPaymentIndexSchedule() ;
	buffer_floatingLegPaymentIndexSchedule_	= buffer_UnderlyingSwap_.get_floatingLegPaymentIndexSchedule() ;
	buffer_deltaTFixedLeg_					= buffer_UnderlyingSwap_.get_DeltaTFixedLeg() ;
	
	initialize_y_bar(buffer_T0_, gridSize);	
	buffer_s0_								= swapRate0();  
	

	std::cout << "constructeur incomplet, initialiser b_barre_S et c_barre_S" << std::endl ;
	//buffer_b_barre_ = timeAverage(buffer_T0_) ;
}

//fonction qui intervient pour le calcul de y bar
//sera integrée sur s entre 0 et t 
double CheyetteQuad_VanillaSwaptionApproxPricer::to_integrate_y_bar(double s) const
{
	double k			= pCheyetteQuad_Model_->get_CheyetteQuad_Parameter().k_ ;
	double sigma_r_0	= pCheyetteQuad_Model_->sigma_r(s, 0.) ;			//sigma_r^0(s) = sigma_r(s, 0) = a(s) 
	return exp(2 * k * s) * sigma_r_0 * sigma_r_0 ; 
}

//integrale jusqu'à t 
// y_bar(t) = exp( - 2 k t ) * integrale_0^t( exp(2 k s) sigma_r^0(s) sigma_r^0(s) ds ) 
void CheyetteQuad_VanillaSwaptionApproxPricer::initialize_y_bar(double t, size_t gridSize) const
{

	std::cout << "calcul de y barre à améliorer, admet forme exacte dans Quadratic Volatility" << std::endl ;

	assert(t > 0);
	double gridStart	= 0 ;
	double gridEnd		= t ;
	double k			= pCheyetteQuad_Model_->get_CheyetteQuad_Parameter().k_ ;

	//constructeur pour le schema d'integration
	numeric::IncrementalIntegrator1D_Riemann integral(gridStart, gridEnd, gridSize);
	//constructeur pour la fonction à intégrer
	boost::function<double(double)> f = 
				boost::bind(&CheyetteQuad_VanillaSwaptionApproxPricer::to_integrate_y_bar, *this, _1);

	integral.vecteur_integrate(f) ;

	std::vector<double> vect_grids = integral.get_grids() ;
	std::vector<double> vect_values = integral.get_values() ;
	
	//multiplication de l'integrale par exp( - 2 k t_i)
	for (size_t i = 0 ; i < vect_values.size() ; ++i)
	{
		//vect_values[i] = exp( - 2. * k * t) * vect_values[i] ;
		vect_values[i] = exp( - 2. * k * vect_grids[i]) * vect_values[i] ;
	}

	buffer_y_bar_.set_grid_(vect_grids) ;
 	buffer_y_bar_.set_value_(vect_values) ;	
}


double CheyetteQuad_VanillaSwaptionApproxPricer::ZC_1stDerivative_on_xt(double t, double T, double x_t, double y_t) const
{
	assert(0 <= T && T <= buffer_TN_) ;

	double ZC = pCheyetteQuad_Model_->P(t, T, x_t, y_t) ;	
	double g = pCheyetteQuad_Model_->G(t,T) ;
	double res = - g * ZC ;

	return res ;	
}

double CheyetteQuad_VanillaSwaptionApproxPricer::ZC_2ndDerivative_on_xt(double t, double T, double x_t, double y_t) const
{
	assert(0 <= T && T <= buffer_TN_) ;

	double ZC = pCheyetteQuad_Model_->P(t, T, x_t, y_t) ;
	double g = pCheyetteQuad_Model_->G(t,T) ;
	double res = g * g * ZC ;
	
	return res ;					
}

double CheyetteQuad_VanillaSwaptionApproxPricer::ZC_3rdDerivative_on_xt(double t, double T, double x_t, double y_t) const
{
	assert(0 <= T && T <= buffer_TN_) ;

	double ZC = pCheyetteQuad_Model_->P(t, T, x_t, y_t) ;
	double g = pCheyetteQuad_Model_->G(t,T) ;
	double res = - g * g * g * ZC ;
	
	return res ;					
}

// Numerator = P(t, T0) - P(t, TN)
//si t = 0 : appel à la courbe spot
//si t > 0 : passage du paramètre x_t pour la fonction inverse
double CheyetteQuad_VanillaSwaptionApproxPricer::swapRateNumerator(double t, double x_t, double y_t) const 
{
	assert( t >= 0. ) ;
	//mono-curve only	
	double ZC_T0 = pCheyetteQuad_Model_->P(t, buffer_T0_, x_t, y_t) ; 
	double ZC_TN = pCheyetteQuad_Model_->P(t, buffer_TN_, x_t, y_t) ;	
	return ZC_T0 - ZC_TN ; 
}

double CheyetteQuad_VanillaSwaptionApproxPricer::swapRateNumerator_1stDerivative(double t, double x_t, double y_t) const
{
	return ZC_1stDerivative_on_xt(t, buffer_T0_, x_t, y_t) - ZC_1stDerivative_on_xt(t, buffer_TN_, x_t, y_t) ;
}

double CheyetteQuad_VanillaSwaptionApproxPricer::swapRateNumerator_2ndDerivative(double t, double x_t, double y_t) const
{
	return ZC_2ndDerivative_on_xt(t, buffer_T0_, x_t, y_t) - ZC_2ndDerivative_on_xt(t, buffer_TN_, x_t, y_t) ;
}

double CheyetteQuad_VanillaSwaptionApproxPricer::swapRateNumerator_3rdDerivative(double t, double x_t, double y_t) const
{
	return ZC_3rdDerivative_on_xt(t, buffer_T0_, x_t, y_t) - ZC_3rdDerivative_on_xt(t, buffer_TN_, x_t, y_t) ;
}

// Denominator = \sum delta_k P(t,T_k) en t = 0
double CheyetteQuad_VanillaSwaptionApproxPricer::swapRateDenominator(double t, double x_t, double y_t) const 
{	
	double fixed_tenor = buffer_UnderlyingSwap_.get_fixedLegTenorType().YearFraction() ;
	double float_tenor = buffer_UnderlyingSwap_.get_floatingLegTenorType().YearFraction() ;
	double tenor_ref = std::min(fixed_tenor, float_tenor) ;  //le plus petit 
	
	double denom = 0.0;
	for(size_t itr = 0 ; itr < buffer_fixedLegPaymentIndexSchedule_.size() ; ++itr) 
	{
		double dateEchangeFluxFixe = buffer_fixedLegPaymentIndexSchedule_[itr] * tenor_ref ;
		double delta_T = buffer_UnderlyingSwap_.get_DeltaTFixedLeg(itr) ;	
		double ZC		= pCheyetteQuad_Model_->P(t, dateEchangeFluxFixe, x_t, y_t) ;
		denom += delta_T * ZC ;
	}
	
	return denom ;
}

double CheyetteQuad_VanillaSwaptionApproxPricer::swapRateDenominator_1stDerivative(double t, double x_t, double y_t) const
{
	double fixed_tenor = buffer_UnderlyingSwap_.get_fixedLegTenorType().YearFraction() ;
	double float_tenor = buffer_UnderlyingSwap_.get_floatingLegTenorType().YearFraction() ;
	double tenor_ref = std::min(fixed_tenor, float_tenor) ;  
	
	double result = 0. ;
	for(size_t itr = 0 ; itr < buffer_fixedLegPaymentIndexSchedule_.size() ; ++itr) 
	{
		double dateEchangeFluxFixe = buffer_fixedLegPaymentIndexSchedule_[itr] * tenor_ref ;
		result += buffer_deltaTFixedLeg_[itr] * ZC_1stDerivative_on_xt(t, dateEchangeFluxFixe, x_t, y_t) ;	
	}
	return result;
}

double CheyetteQuad_VanillaSwaptionApproxPricer::swapRateDenominator_2ndDerivative(double t, double x_t, double y_t) const
{
	double fixed_tenor = buffer_UnderlyingSwap_.get_fixedLegTenorType().YearFraction() ;
	double float_tenor = buffer_UnderlyingSwap_.get_floatingLegTenorType().YearFraction() ;
	double tenor_ref = std::min(fixed_tenor, float_tenor) ;  

	double result = 0. ;
	for(size_t itr = 0; itr < buffer_fixedLegPaymentIndexSchedule_.size(); ++itr) 
	{
		double dateEchangeFluxFixe = buffer_fixedLegPaymentIndexSchedule_[itr] * tenor_ref ;
		result += buffer_deltaTFixedLeg_[itr] * ZC_2ndDerivative_on_xt(t, dateEchangeFluxFixe, x_t, y_t) ;	
	}
	return result;
}

double CheyetteQuad_VanillaSwaptionApproxPricer::swapRateDenominator_3rdDerivative(double t, double x_t, double y_t) const
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

double CheyetteQuad_VanillaSwaptionApproxPricer::swapRate0() const
{
	return swapRate(0., 0., 0.);
}

double CheyetteQuad_VanillaSwaptionApproxPricer::swapRate(double t, double x_t, double y_t) const
{
	double n      = swapRateNumerator(t, x_t, y_t);
	double d      = swapRateDenominator(t, x_t, y_t);

	assert (d > 0.001) ;
	return n/d;
}

double CheyetteQuad_VanillaSwaptionApproxPricer::swapRate_1stDerivative(double t, double x_t, double y_t) const
{
	double n   = swapRateNumerator(t, x_t, y_t);
	double n_1 = swapRateNumerator_1stDerivative(t, x_t, y_t); 

	double d   = swapRateDenominator(t, x_t, y_t);
	double d_1 = swapRateDenominator_1stDerivative(t, x_t, y_t);
	//std::cout << "swap rate 1st derivative : " << (n_1*d - n*d_1)/(d*d) << std::endl ;
	return (n_1*d - n*d_1)/(d*d);
}

double CheyetteQuad_VanillaSwaptionApproxPricer::swapRate_2ndDerivative(double t, double x_t, double y_t) const
{
	double n	= swapRateNumerator(t, x_t, y_t);
	double n_1	= swapRateNumerator_1stDerivative(t, x_t, y_t); 
	double n_2	= swapRateNumerator_2ndDerivative(t, x_t, y_t); 

	double d	= swapRateDenominator(t, x_t, y_t);
	double d_1	= swapRateDenominator_1stDerivative(t, x_t, y_t);
	double d_2	= swapRateDenominator_2ndDerivative(t, x_t, y_t);

	//double result = (n_2 * d - n * d_2) * d*d - (n_1 * d - n * d_1) * 2 * d * d_1 ;
	//result /= (d*d*d*d);

	//simplification par d
	double result = (n_2 * d - n * d_2) * d - (n_1 * d - n * d_1) * 2 * d_1 ;
	result /= (d*d*d);

	return result;
}

double CheyetteQuad_VanillaSwaptionApproxPricer::swapRate_3rdDerivative(double t, double x_t, double y_t) const
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

/*************************     Newton - Raphson        ************************/
struct Newton_Raphson_struct
{
private:
	boost::function<double(double)> f_ ;
	boost::function<double(double)> f_derivative_ ;
	double target_ ;
	
public:
	Newton_Raphson_struct(	const boost::function<double(double)>& f,
							const boost::function<double(double)>& f_derivative, double target) 
				: f_(f), f_derivative_(f_derivative), target_(target){}

	boost::math::tuple<double, double> operator()(double x)
	{
		return boost::math::make_tuple(f_(x) - target_, f_derivative_(x));
	}
};

//S(t, x_t) = swapRate(t, x_t) = s 
//retourne le x_t correspondant
double CheyetteQuad_VanillaSwaptionApproxPricer::inverse(double t, double s) const
{                                                                  //swapRate(double t, double x_t)
	double y_bar = buffer_y_bar_(t) ;
	boost::function<double(double)> f				= 
				boost::bind(&CheyetteQuad_VanillaSwaptionApproxPricer::swapRate, this, t, _1, y_bar);
	boost::function<double(double)> f_derivative	= 
				boost::bind(&CheyetteQuad_VanillaSwaptionApproxPricer::swapRate_1stDerivative, this, t, _1, y_bar);
	
	Newton_Raphson_struct NR_struct(f, f_derivative, s);   //struct qui contient f, f', target = s
	double initial_guess	= 0.0 ;
	double min				= -5.0;
	double max				= 5.0;
    size_t nDigits			= 15;
//	boost::uintmax_t nMaxIter  = 50 ;

	double result_newton_raphson = boost::math::tools::newton_raphson_iterate(NR_struct, initial_guess, min, max, nDigits);
	return result_newton_raphson;
}

//DL(2) de Phi(t, s) en s_barre = S0
double CheyetteQuad_VanillaSwaptionApproxPricer::DL2_Phi_t_s(double t, double s) const
{
	double inverse_s_bar =	inverse(t, buffer_s0_) ; 
	return Phi(t, inverse_s_bar )	+ DPhi(t, inverse_s_bar ) * (s - buffer_s0_) 
									+ 0.5 * D2Phi(t, inverse_s_bar) * pow(s - buffer_s0_, 2) ;
}

double CheyetteQuad_VanillaSwaptionApproxPricer::Phi(double t, double x_t) const
{
	//retourne dS(t)/dx(t) sigma_r(t) (t=0, inverse, y_bar)
	double y_bar = buffer_y_bar_(t) ;

	double sigma_r_t_inv_sbar	= pCheyetteQuad_Model_->sigma_r(t, x_t) ;	
	double phi_res = swapRate_1stDerivative(t, x_t, y_bar) * sigma_r_t_inv_sbar ;
	
	return phi_res ;

}

// d Phi / ds
double CheyetteQuad_VanillaSwaptionApproxPricer::DPhi(double t, double x_t) const
{
	double y_bar = buffer_y_bar_(t) ;

	double sigma_r_t	= pCheyetteQuad_Model_->sigma_r(t, x_t) ;
	double res =  sigma_r_t *	swapRate_2ndDerivative(t, x_t, y_bar) / swapRate_1stDerivative(t, x_t, y_bar) 
					+ pCheyetteQuad_Model_->sigma_r_t_1stDerivative(t, x_t) ;

	return res ; 
}

// d2 Phi / ds
double CheyetteQuad_VanillaSwaptionApproxPricer::D2Phi(double t, double x_t) const
{
	double y_bar = buffer_y_bar_(t) ;

	double sigma_r		= pCheyetteQuad_Model_->sigma_r(t, x_t) ;
	double sigma_r_1	= pCheyetteQuad_Model_->sigma_r_t_1stDerivative(t, x_t) ;
	double sigma_r_2	= pCheyetteQuad_Model_->sigma_r_t_2ndDerivative(t, x_t) ;
	double S_1	= swapRate_1stDerivative(t, x_t, y_bar); 
	double S_2	= swapRate_2ndDerivative(t, x_t, y_bar); 
	double S_3	= swapRate_3rdDerivative(t, x_t, y_bar); 

	double res = sigma_r * ( S_3/(S_1 * S_1) - pow(S_2, 2) / pow(S_1, 3) )
					+ sigma_r_1 * S_2 / (S_1*S_1) 
					+ sigma_r_2 / S_1 ;
		 
	return res ; 
}

/* dS(t) = \lambda(t) ( S0 + b_S(t) (S_t - s_bar) + 1./2. c_S(t) (S_t - s_bar)^2 ) dW_t QA */

//en prenant s0 = s_bar
//lambda(t, s_bar)
double CheyetteQuad_VanillaSwaptionApproxPricer::lambda(double t) const
{
	double inverse_s_bar =	inverse(t, buffer_s0_) ; 
	return Phi(t, inverse_s_bar) / buffer_s0_ ;
}

//b_S (t, s_bar)
double CheyetteQuad_VanillaSwaptionApproxPricer::b_S(double t) const 
{
	double inverse_s_bar =	inverse(t, buffer_s0_) ; 
	return DPhi(t, inverse_s_bar) / lambda(t) ;
}

//c_S (t, s_bar)
double CheyetteQuad_VanillaSwaptionApproxPricer::c_S(double t) const
{
	double inverse_s_bar =	inverse(t, buffer_s0_) ; 
	return D2Phi(t, inverse_s_bar) / lambda(t) ;
}

//autres fonctions intermédiaires pour les time averaging
double CheyetteQuad_VanillaSwaptionApproxPricer::lambda2(double t) const
{
	return lambda(t) * lambda(t) ; 
}

double CheyetteQuad_VanillaSwaptionApproxPricer::lambda4(double t) const
{
	return pow(lambda(t), 4) ; 
}

double CheyetteQuad_VanillaSwaptionApproxPricer::f_outer_num_b(double t) const
{
	return lambda2(t) * b_S(t) ;  //OK (fonction outer seulement ie sans v^2(u)
}
double CheyetteQuad_VanillaSwaptionApproxPricer::f_outer_denom(double t) const
{
	return lambda2(t)  ;		//OK (fonction outer seulement ie sans v^2(u)
}

double CheyetteQuad_VanillaSwaptionApproxPricer::f_outer_num_c(double t) const
{
	return lambda2(t) * c_S(t) ;  //OK (fonction outer seulement ie sans v^2(u)
}

//time averaging of b_S(t) and c_S(t)
double CheyetteQuad_VanillaSwaptionApproxPricer::timeAverage_b_S(double t) const
{
	double gridStart = 0.0;
	double gridEnd = t ; 

//integrale numerateur
	numeric::IncrementalIntegrator2D_Riemann int2D(gridStart, gridEnd, gridSize) ;
	boost::function<double(double)> f_inner = boost::bind(&CheyetteQuad_VanillaSwaptionApproxPricer::lambda2, *this, _1);
	boost::function<double(double)> f_outer = boost::bind(&CheyetteQuad_VanillaSwaptionApproxPricer::f_outer_num_b, *this, _1);	
	double integrale_numerateur = int2D.integrate(f_outer, f_inner) ;
	
//integrale denominateur
	boost::function<double(double)> f_denom = boost::bind(&CheyetteQuad_VanillaSwaptionApproxPricer::f_outer_denom, *this, _1);
	double integrale_denom = int2D.integrate(f_denom, f_inner);

	double b_barre = integrale_numerateur/integrale_denom ;

	return b_barre ;
}

double CheyetteQuad_VanillaSwaptionApproxPricer::timeAverage_c_S(double t) const
{
	double gridStart = 0.0;
	double gridEnd = t ; 

//integrale numerateur
	numeric::IncrementalIntegrator2D_Riemann int2D(gridStart, gridEnd, gridSize) ;
	boost::function<double(double)> f_inner = boost::bind(&CheyetteQuad_VanillaSwaptionApproxPricer::lambda4, *this, _1);
	boost::function<double(double)> f_outer = boost::bind(&CheyetteQuad_VanillaSwaptionApproxPricer::f_outer_num_c, *this, _1);	
	double integrale_numerateur = int2D.integrate(f_outer, f_inner) ;
	
//integrale denominateur
	boost::function<double(double)> f_denom = boost::bind(&CheyetteQuad_VanillaSwaptionApproxPricer::f_outer_denom, *this, _1);
	double integrale_denom = int2D.integrate(f_denom, f_inner);

	double c_barre = integrale_numerateur/integrale_denom ;

	return c_barre ;
}


/**********************  Call approximation ************************
** dynamique : 
**
**  dS(t) = \Phi(S_t) dW(t)
********************************************************************/




