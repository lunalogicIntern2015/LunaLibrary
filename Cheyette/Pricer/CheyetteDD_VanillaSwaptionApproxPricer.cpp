#include <Cheyette\Pricer\CheyetteDD_VanillaSwaptionApproxPricer.h>
#include "LMM/numeric/NumericalMethods.h"

#include <cassert>
#include <vector>
#include <iostream> 
//#include <boost/pointer_cast.hpp>


const size_t gridSize  = 101 ;


CheyetteDD_VanillaSwaptionApproxPricer::CheyetteDD_VanillaSwaptionApproxPricer
						(const CheyetteDD_Model_PTR&	cheyetteDD_Model, 
						const VanillaSwaption_PTR&		swaption) 
						:cheyetteDD_Model_(cheyetteDD_Model), swaption_(swaption)  					
{
	// domaine de validité de l'approximation : //assert( ... ); 

	buffer_courbeInput_	 = cheyetteDD_Model->get_courbeInput_PTR() ;
	initialize_buffers() ;  //initialise les autres buffers relatifs à la swaption, y barre, b barre 		
}

void CheyetteDD_VanillaSwaptionApproxPricer::initialize_buffers()
{
	buffer_UnderlyingSwap_					= swaption_->getUnderlyingSwap() ;
	buffer_T0_								= swaption_->getUnderlyingSwap().get_StartDate() ;		//1ere date de fixing
	buffer_TN_								= swaption_->getUnderlyingSwap().get_EndDate() ;			//date dernier flux (setté en T_{N-1})
	buffer_fixedLegPaymentIndexSchedule_	= buffer_UnderlyingSwap_.get_fixedLegPaymentIndexSchedule() ;
	buffer_floatingLegPaymentIndexSchedule_	= buffer_UnderlyingSwap_.get_floatingLegPaymentIndexSchedule() ;
	buffer_deltaTFixedLeg_					= buffer_UnderlyingSwap_.get_DeltaTFixedLeg() ;
	
	//initialisation du shift de CheyetteDD_Model pas encore initialisé (c'est le cas si shift = S(t) / S(0) ) 
	//shift = S(t) / S(0) dépend de la swaption. Non initialisable dans CheyetteDD_Model
	if (cheyetteDD_Model_->get_shiftChoice() == 3)
	{
		boost::function<double(double, double)> f1 = boost::bind(&CheyetteDD_VanillaSwaptionApproxPricer::swapRate, this, _1, _2);
		boost::function<double(double, double)> f2 = boost::bind(&CheyetteDD_VanillaSwaptionApproxPricer::swapRate, this, 0., 0.);
		boost::function<double(double, double)> fp1 = boost::bind(&CheyetteDD_VanillaSwaptionApproxPricer::swapRate_1stDerivative, this, _1, _2);
		boost::function<double(double, double)> fp2 = boost::bind(&CheyetteDD_VanillaSwaptionApproxPricer::swapRate_1stDerivative, this, 0., 0.);
		Boost_R2R_Function_PTR f1_ptr(new Boost_R2R_Function(f1)) ;
		Boost_R2R_Function_PTR f2_ptr(new Boost_R2R_Function(f2)) ;
		Boost_R2R_Function_PTR fp1_ptr(new Boost_R2R_Function(fp1)) ;
		Boost_R2R_Function_PTR fp2_ptr(new Boost_R2R_Function(fp2)) ;
		cheyetteDD_Model_->setShiftPointer(f1_ptr, f2_ptr, fp1_ptr, fp2_ptr) ;
	}

	initialize_y_bar(buffer_T0_, gridSize);	
	buffer_s0_								= swapRate0();  //necessite y_bar
	
	buffer_b_barre_ = timeAverage(buffer_T0_) ;
}

//fonction qui intervient pour le calcul de y bar
//sera integrée sur s entre 0 et t 
double CheyetteDD_VanillaSwaptionApproxPricer::to_integrate_y_bar(double s) const
{
	double k			= cheyetteDD_Model_->get_CheyetteDD_Parameter().k_ ;
	double sigma_r_0	= cheyetteDD_Model_->sigma_r(s, 0.) ;			//sigma_r^0(t) = sigma_r(t, 0, 0) 
	return exp(2 * k * s) * sigma_r_0 * sigma_r_0 ; 
}

//integrale jusqu'à t 
// y_bar(t) = exp( - 2 k t ) * integrale_0^t( exp(2 k s) sigma_r^0(s) sigma_r^0(s) ds ) 
void CheyetteDD_VanillaSwaptionApproxPricer::initialize_y_bar(double t, size_t gridSize) const
{
	assert(t > 0);
	double gridStart	= 0 ;
	double gridEnd		= t ;
	double k			= cheyetteDD_Model_->get_CheyetteDD_Parameter().k_ ;

	//constructeur pour le schema d'integration
	numeric::IncrementalIntegrator1D_Riemann integral(gridStart, gridEnd, gridSize);
	//constructeur pour la fonction à intégrer
	boost::function<double(double)> f = boost::bind(&CheyetteDD_VanillaSwaptionApproxPricer::to_integrate_y_bar, *this, _1);

	integral.vecteur_integrate(f) ;
	//integral.integrate(f) ;  

	std::vector<double> vect_grids = integral.get_grids() ;
	std::vector<double> vect_values = integral.get_values() ;
	
	//multiplication de l'integrale par exp( - 2 k t_i)
	for (size_t i = 0 ; i < vect_values.size() ; ++i)
	{
		vect_values[i] = exp( - 2. * k * vect_grids[i]) * vect_values[i] ;
	}

	buffer_y_bar_.set_grid_(vect_grids) ;
 	buffer_y_bar_.set_value_(vect_values) ;	
}


double CheyetteDD_VanillaSwaptionApproxPricer::ZC_1stDerivative_on_xt(double t, double T, double x_t) const
{
	assert(0 <= T && T <= buffer_TN_) ;
	double res ;
	
	if (t == 0){     //courbe spot
		double tauxZC = buffer_courbeInput_->get_tauxZC0(T) ;
		double P_0_T  = exp(- tauxZC * T) ;
		double g = cheyetteDD_Model_->G(0,T) ;
		res = - g * P_0_T ;											//- G(0, T) P(0, T) 
	}
	else{			//modèle
		double y_bar_t = buffer_y_bar_(t) ;
		double ZC = cheyetteDD_Model_->P(t, T, x_t, y_bar_t) ;		//approximation
		double g = cheyetteDD_Model_->G(t,T) ;
		res = - g * ZC ;
	}
	return res ;	
}

double CheyetteDD_VanillaSwaptionApproxPricer::ZC_2ndDerivative_on_xt(double t, double T, double x_t) const
{
	assert(0 <= T && T <= buffer_TN_) ;
	double y_bar_t = buffer_y_bar_(t) ;
	double res ;
	
	if (t == 0){
		double tauxZC = buffer_courbeInput_->get_tauxZC0(T) ;
		double P_0_T  = exp(- tauxZC * T) ;
		double g = cheyetteDD_Model_->G(0,T) ;
		res = g * g * P_0_T ;				//- G(0, T)^2 P(0, T) 
	}
	else{			
		double ZC = cheyetteDD_Model_->P(t, T, x_t, y_bar_t) ;
		double g = cheyetteDD_Model_->G(t,T) ;
		res = g * g * ZC ;
	}
	return res ;					
}

// Numerator = P(t, T0) - P(t, TN)
//si t = 0 : appel à la courbe spot
//si t > 0 : passage du paramètre x_t pour la fonction inverse
double CheyetteDD_VanillaSwaptionApproxPricer::swapRateNumerator(double t, double x_t) const 
{
	assert( t >= 0 ) ;
	double y_bar_t ;
	if (this->get_CheyetteDD_Model()->get_shiftChoice() == 3) //cas shift S(t) / S(0)
	{
		y_bar_t = 0. ;
	}
	else
	{
		y_bar_t = buffer_y_bar_(t) ;
	}
	
	double ZC_T0, ZC_TN ;
//les 2 versions sont cohérentes en t = 0 :	
	if (t == 0){  // bad comparision   // prefer to do: if (t==0) { check x_t ==0 }
		ZC_T0 = exp( - buffer_courbeInput_->get_tauxZC0(buffer_T0_) * buffer_T0_) ;
		ZC_TN = exp( - buffer_courbeInput_->get_tauxZC0(buffer_TN_) * buffer_TN_) ; 
	}
	else
	{
		ZC_T0 = cheyetteDD_Model_->P(t, buffer_T0_, x_t, y_bar_t) ; 
		ZC_TN = cheyetteDD_Model_->P(t, buffer_TN_, x_t, y_bar_t) ;	
	}
	return  ZC_T0 - ZC_TN ; 
}

double CheyetteDD_VanillaSwaptionApproxPricer::swapRateNumerator_1stDerivative(double t, double x_t) const
{
	//mono-curve only
	return ZC_1stDerivative_on_xt(t, buffer_T0_, x_t) - ZC_1stDerivative_on_xt(t, buffer_TN_, x_t) ;
}

double CheyetteDD_VanillaSwaptionApproxPricer::swapRateNumerator_2ndDerivative(double t, double x_t) const
{
	return ZC_2ndDerivative_on_xt(t, buffer_T0_, x_t) - ZC_2ndDerivative_on_xt(t, buffer_TN_, x_t) ;
}

// Denominator = \sum delta_k P(t,T_k) en t = 0
double CheyetteDD_VanillaSwaptionApproxPricer::swapRateDenominator(double t, double x_t) const 
{
	double y_bar_t ;
	if (this->get_CheyetteDD_Model()->get_shiftChoice() == 3) //cas shift S(t) / S(0)
	{
		y_bar_t = 0. ;
	}
	else
	{
		y_bar_t = buffer_y_bar_(t) ;
	}
	
	double price = 0.0;

	double fixed_tenor = buffer_UnderlyingSwap_.get_fixedLegTenorType().YearFraction() ;
	double float_tenor = buffer_UnderlyingSwap_.get_floatingLegTenorType().YearFraction() ;
	double tenor_ref = std::min(fixed_tenor, float_tenor) ;  //le plus petit 

//coherence entre les 2 boucles (la 2eme boucle evaluée en t=0 est identique à la 1ere)
	if (t == 0)
	{
		//somme sur tous les flux fixes
		for(size_t itr = 0; itr < buffer_fixedLegPaymentIndexSchedule_.size(); ++itr) 
		{
			//convertit l'indice/le numero du flux en la date de tombée du flux (ex : flux numero 2 survient à date 1Y)
			double dateEchangeFluxFixe = buffer_fixedLegPaymentIndexSchedule_[itr] * tenor_ref ;

			double delta_T = buffer_UnderlyingSwap_.get_DeltaTFixedLeg(itr);	
			double ZC		= exp( - buffer_courbeInput_->get_tauxZC0(dateEchangeFluxFixe) * dateEchangeFluxFixe);
			price += delta_T * ZC ;		
		}
	}else
	{
		for(size_t itr = 0; itr < buffer_fixedLegPaymentIndexSchedule_.size(); ++itr) 
		{
			double dateEchangeFluxFixe = buffer_fixedLegPaymentIndexSchedule_[itr] * tenor_ref ;
			double delta_T = buffer_UnderlyingSwap_.get_DeltaTFixedLeg(itr) ;	
			double ZC		= cheyetteDD_Model_->P(t, dateEchangeFluxFixe, x_t, y_bar_t) ;
			price += delta_T * ZC ;
		}
	}
	return price;
}

double CheyetteDD_VanillaSwaptionApproxPricer::swapRateDenominator_1stDerivative(double t, double x_t) const
{
	double result = 0.;

	double fixed_tenor = buffer_UnderlyingSwap_.get_fixedLegTenorType().YearFraction() ;
	double float_tenor = buffer_UnderlyingSwap_.get_floatingLegTenorType().YearFraction() ;
	double tenor_ref = std::min(fixed_tenor, float_tenor) ;  

	for(size_t itr = 0; itr < buffer_fixedLegPaymentIndexSchedule_.size(); ++itr) 
	{
		double dateEchangeFluxFixe = buffer_fixedLegPaymentIndexSchedule_[itr] * tenor_ref ;
		result += buffer_deltaTFixedLeg_[itr] * ZC_1stDerivative_on_xt(t, dateEchangeFluxFixe, x_t) ;	
	}
	return result;
}

double CheyetteDD_VanillaSwaptionApproxPricer::swapRateDenominator_2ndDerivative(double t, double x_t) const
{
	double result = 0;

	double fixed_tenor = buffer_UnderlyingSwap_.get_fixedLegTenorType().YearFraction() ;
	double float_tenor = buffer_UnderlyingSwap_.get_floatingLegTenorType().YearFraction() ;
	double tenor_ref = std::min(fixed_tenor, float_tenor) ;  

	for(size_t itr = 0; itr < buffer_fixedLegPaymentIndexSchedule_.size(); ++itr) 
	{
		double dateEchangeFluxFixe = buffer_fixedLegPaymentIndexSchedule_[itr] * tenor_ref ;
		result += buffer_deltaTFixedLeg_[itr] * ZC_2ndDerivative_on_xt(t, dateEchangeFluxFixe, x_t) ;	
	}
	return result;
}

double CheyetteDD_VanillaSwaptionApproxPricer::swapRate0() const
{
	return swapRate(0, 0);
}
double CheyetteDD_VanillaSwaptionApproxPricer::swapRate(double t, double x_t) const
{

	double n      = swapRateNumerator(t, x_t);
	double d      = swapRateDenominator(t, x_t);


	assert (d > 0.001) ;
	return n/d;
}

double CheyetteDD_VanillaSwaptionApproxPricer::swapRate_1stDerivative(double t, double x_t) const
{
	double n   = swapRateNumerator(t, x_t);
	double n_1 = swapRateNumerator_1stDerivative(t, x_t); 

	double d   = swapRateDenominator(t, x_t);
	double d_1 = swapRateDenominator_1stDerivative(t, x_t);
	//std::cout << "swap rate 1st derivative : " << (n_1*d - n*d_1)/(d*d) << std::endl ;
	return (n_1*d - n*d_1)/(d*d);
}

double CheyetteDD_VanillaSwaptionApproxPricer::swapRate_2ndDerivative(double t, double x_t) const
{
	double n	= swapRateNumerator(t, x_t);
	double n_1	= swapRateNumerator_1stDerivative(t, x_t); 
	double n_2	= swapRateNumerator_2ndDerivative(t, x_t); 

	double d	= swapRateDenominator(t, x_t);
	double d_1	= swapRateDenominator_1stDerivative(t, x_t);
	double d_2	= swapRateDenominator_2ndDerivative(t, x_t);

	double result = (n_2 * d - n * d_2) * d*d - (n_1 * d - n * d_1) * 2 * d * d_1 ;
	result /= (d*d*d*d);

	return result;
}

/************************  fonction pour inverser     *************************
*************************     Newton - Raphson        ************************/
struct Newton_Raphson_struct
{
private:
	boost::function<double(double)> f_ ;
	boost::function<double(double)> f_derivative_ ;
	double target_ ;
	
public:
	//! constructor
	Newton_Raphson_struct(	const boost::function<double(double)>& f,
							const boost::function<double(double)>& f_derivative, 
							double target) 
				: f_(f), f_derivative_(f_derivative), target_(target){}

	//! evaluation
	//callable function object that accepts one parameter and returns a tuple:
	//the tuple should have two elements containing the evaluation of the function and it's first derivative.
	boost::math::tuple<double, double> operator()(double x)
	{
		return boost::math::make_tuple(
			f_(x) - target_,
			f_derivative_(x));
	}
};

//S(t, x_t) = swapRate(t, x_t) = s 
//retourne le x_t correspondant
double CheyetteDD_VanillaSwaptionApproxPricer::inverse(double t, double s) const
{                                                                  //swapRate(double t, double x_t)
	boost::function<double(double)> f				= 
				boost::bind(&CheyetteDD_VanillaSwaptionApproxPricer::swapRate, this, t, _1);
	boost::function<double(double)> f_derivative	= 
				boost::bind(&CheyetteDD_VanillaSwaptionApproxPricer::swapRate_1stDerivative, this, t, _1);
	
	Newton_Raphson_struct NR_struct(f, f_derivative, s);   //struct qui contient f, f', target = s
	double initial_guess	= 0.0 ;
	double min				= -5.0;
	double max				= 5.0;
    size_t nDigits			= 15;
//	boost::uintmax_t nMaxIter  = 50 ;

	double result_newton_raphson = boost::math::tools::newton_raphson_iterate(NR_struct, initial_guess, min, max, nDigits);
	return result_newton_raphson;
}


//derivee par rapport à x_t
//prend le paramètre x_t
//en paramètre donner f-1(S) = x_t
double CheyetteDD_VanillaSwaptionApproxPricer::swapRateVolatility_1stDerivative(double t, double x_t) const
{
	double sigma_r_t	= cheyetteDD_Model_->sigma_r(t, x_t) ;
	double res =  sigma_r_t *	swapRate_2ndDerivative(t, x_t) / swapRate_1stDerivative(t, x_t) 
					+ cheyetteDD_Model_->sigma_r_t_1stDerivative(t, x_t) ;

	//std::cout << "t : " << t << "swap rate vol 1st der : " << res << std::endl ;
	return res ; 
	
}

//pour un taux de swap s_bar = s0_
double CheyetteDD_VanillaSwaptionApproxPricer::calculate_phi_t_s_bar(double t) const
{
	//retourne dS(t)/dx(t) sigma_r(t) (t=0, inverse, y_bar)

	double inverse_s_bar =	inverse(t, buffer_s0_) ;  

	double sigma_r_t_inv_sbar	= cheyetteDD_Model_->sigma_r(t, inverse_s_bar) ;	
	double phi_res = swapRate_1stDerivative(t, inverse_s_bar) * sigma_r_t_inv_sbar ;
	
	//std::cout << "derivee = " << swapRate_1stDerivative(t, inverse_s_bar) << std::endl ;
	//std::cout << "vol_t = " << phi_res << std::endl ;
	return phi_res ;
}

//s_bar = S0
//double CheyetteDD_VanillaSwaptionApproxPricer::swapRateVolatility_approx_lineaire(double t, double s) const
//{
//	//dérivée à évaluer en x_t = S^{-1}(s)
//	//calcul de l'inverse de S(t, x_t) = s 
//	//retourne x_t
//	double inverse_s_bar =	inverse(t, buffer_s0_) ;  
//	return calculate_phi_t_s_bar(t) + swapRateVolatility_1stDerivative(t, inverse_s_bar) * (s - buffer_s0_) ;
//}


/****************  parameter averaging  *********************
*
*	dS(t) = \lambda(t) ( (1-b(t)) S0 + b(t) S(t) ) dW^Q_A		-> time averaging: b(t) puis b_barre
*
*	dS(t) = (A(t) + B(t) S(t)) dW^Q_A
*	S(0) = S0_ connu
************************************************************/

//A(t) 
double CheyetteDD_VanillaSwaptionApproxPricer::A(double t) const
{
	double inverse_s_bar =	inverse(t, buffer_s0_) ;
	double A_t = calculate_phi_t_s_bar(t) - swapRateVolatility_1stDerivative(t, inverse_s_bar) * buffer_s0_ ;

	//std::cout.precision(15);
	//std::cout << "t =    " << t << std::endl ;
	//std::cout << inverse_s_bar << std::endl ;
	//std::cout << "A(t) = " << A_t << std::endl ;

	//std::cout << "calculate_phi_t_s_bar            : " << calculate_phi_t_s_bar(t) << std::endl ;
	//std::cout << "swapRateVolatility_1stDerivative : " << swapRateVolatility_1stDerivative(t, inverse_s_bar) << std::endl ;
	//std::cout << "buffer_s0_                       : " << buffer_s0_ << std::endl ;
	return A_t ;
}

//à optimiser pour éviter le recalcul entre A(t) et B(t)
//B(t) 
double CheyetteDD_VanillaSwaptionApproxPricer::B(double t) const
{
	double inverse_s_bar =	inverse(t, buffer_s0_) ;
	double B_t = swapRateVolatility_1stDerivative(t, inverse_s_bar) ;

	//std::cout.precision(15);
	//std::cout << "t    = " << t << std::endl ;
	//std::cout << "B(t) = " << B_t << std::endl ;

	//std::cout << "y barre : " << get_buffer_y_bar_t(t) << std::endl ;
	//std::cout << swapRate(t, inverse_s_bar) << std::endl ;
	//std::cout << "  " << std::endl ;
	return B_t ;
}

//lambda(t) 
double CheyetteDD_VanillaSwaptionApproxPricer::lambda(double t) const
{
	return A(t) / buffer_s0_ + B(t) ;
	//return cheyetteDD_Model_->get_CheyetteDD_Parameter().sigma_(t) ;  //test !!! à commenter après
}	
double CheyetteDD_VanillaSwaptionApproxPricer::lambda2(double t) const
{
	double l = lambda(t) ;
	return l * l ;
}	
//b(t) 
double CheyetteDD_VanillaSwaptionApproxPricer::b(double t) const
{
	//return 1 / (1 + A(t)/(B(t) * buffer_s0_)) ;
	//double b_t = B(t) / lambda(t) ;
	double b_t = 1 / (1 + A(t)/(B(t) * buffer_s0_)) ;
	//std::cout << " t : " << t << ", b_t : " << b_t << std::endl ;
	return  b_t ;
}	


double CheyetteDD_VanillaSwaptionApproxPricer::f_outer_num(double t) const
{
	return lambda2(t) * b(t) ;  //OK (fonction outer seulement ie sans v^2(u)
}
double CheyetteDD_VanillaSwaptionApproxPricer::f_outer_denom(double t) const
{
	return lambda2(t)  ;		//OK (fonction outer seulement ie sans v^2(u)
}



//retourne b_barre
//average over [0, t]
//gridSize : nb de points pour integrale Riemann
//gridSize + 1 : pour 100 mettre 101, delta_t = 1/100
double CheyetteDD_VanillaSwaptionApproxPricer::timeAverage(double t) const	
{
	double gridStart = 0.0;
	double gridEnd = t ; 

//integrale numerateur

	numeric::IncrementalIntegrator2D_Riemann int2D(gridStart, gridEnd, gridSize) ;
	boost::function<double(double)> f_inner = boost::bind(&CheyetteDD_VanillaSwaptionApproxPricer::lambda2, *this, _1);
	boost::function<double(double)> f_outer = boost::bind(&CheyetteDD_VanillaSwaptionApproxPricer::f_outer_num, *this, _1);
	
	double integrale_numerateur = int2D.integrate(f_outer, f_inner) ;
	
//integrale denominateur
	boost::function<double(double)> f_denom = boost::bind(&CheyetteDD_VanillaSwaptionApproxPricer::f_outer_denom, *this, _1);
	
	double integrale_denom = int2D.integrate(f_denom, f_inner);
	double b_barre = integrale_numerateur/integrale_denom ;
	//buffer_b_barre_ = b_barre ;
	return b_barre ;
}

double CheyetteDD_VanillaSwaptionApproxPricer::prixSwaptionApproxPiterbarg() const
{
	//calcul de b_barre, time averaging
	double b_barre  = timeAverage(buffer_T0_) ;		//average jusqu'à T0, date d'entrée dans le swap

	//std::cout	<< "strike K : " << buffer_UnderlyingSwap_.get_strike() 
	//			<< ", b_barre : " << b_barre 
	//			<< ", s0 : " << buffer_s0_
	//			<< ", r0 : " << buffer_courbeInput_->get_f_0_t(0) << std::endl ;

	//calcul de la variance (integrale de lambda_t)
	double gridStart = 0.0;
	double gridEnd = buffer_T0_ ;
	numeric::Integrator1D_Riemann integral_Riemann(gridStart, gridEnd, gridSize) ;

	boost::function<double(double)> func1 = boost::bind(&CheyetteDD_VanillaSwaptionApproxPricer::lambda2, this, _1);
	double integrale = integral_Riemann.integrate(func1);
	
//std::cout << "b barre : " << b_barre << ", integrale de lambda_t : " << integrale << std::endl ;

	//prix Black swaption (approximation)
	double annuity0 = swapRateDenominator(0., 0.) ;	//en t = 0 c'est l'annuité(0)
	double K_tilde	= b_barre * buffer_UnderlyingSwap_.get_strike() + (1 - b_barre) * buffer_s0_ ;
	double sigma_sqrt_T = sqrt(integrale * b_barre * b_barre) ;  //b_barre * sqrt(integrale) 

	//prend en compte strikes positifs et négatifs
	//	double Black_Price_vol2(double fwd, double strike, double vol, double T);
	return annuity0 / b_barre * NumericalMethods::Black_Price_vol2_allStrike(buffer_s0_, K_tilde, sigma_sqrt_T, buffer_T0_) ;
	
}


