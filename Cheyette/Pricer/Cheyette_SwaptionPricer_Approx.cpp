#include <Cheyette/Pricer/Cheyette_SwaptionPricer_Approx.h>

const size_t gridSize  = 101 ;	//mettre 11 ou 101 (11 est suffisant)


Cheyette_SwaptionPricer_Approx::Cheyette_SwaptionPricer_Approx(	const CheyetteModel_PTR& pCheyette_Model) 
						:pCheyette_Model_(pCheyette_Model)			
{
	buffer_pCourbeInput_	 = pCheyette_Model_->get_courbeInput_PTR() ;
}

//void Cheyette_SwaptionPricer_Approx::preCalculateSwaptionBuffers(VanillaSwaption_PTR pSwaption) const
//{
//	buffer_UnderlyingSwap_					= pSwaption->getUnderlyingSwap() ;
//	buffer_T0_								= pSwaption->getUnderlyingSwap().get_StartDate() ;	//1ere date de fixing
//	buffer_TN_								= pSwaption->getUnderlyingSwap().get_EndDate() ;	//date dernier flux (setté en T_{N-1})
//	buffer_fixedLegPaymentIndexSchedule_	= buffer_UnderlyingSwap_.get_fixedLegPaymentIndexSchedule() ;
//	buffer_deltaTFixedLeg_					= buffer_UnderlyingSwap_.get_DeltaTFixedLeg() ;
//}
//
void Cheyette_SwaptionPricer_Approx::preCalculateModelBuffers() const 
{
	initialize_y_bar(buffer_T0_, gridSize);	
	//virtual : initialisation de b barre ou b barre et c barre dans les classes dérivées
	buffer_s0_								= swapRate0();     //à reinitialiser en dernier si bug
}

void Cheyette_SwaptionPricer_Approx::preCalculateALL(VanillaSwaption_PTR pSwaption) const
{
	//preCalculateSwaptionBuffers(pSwaption) ;

	//preCalculateModelBuffers() ;	//pour la calibration, recalcul des buffers lorsque les paramètres sont modifiés

	buffer_UnderlyingSwap_					= pSwaption->getUnderlyingSwap() ;
	buffer_T0_								= pSwaption->getUnderlyingSwap().get_StartDate() ;	//1ere date de fixing
	buffer_TN_								= pSwaption->getUnderlyingSwap().get_EndDate() ;	//date dernier flux (setté en T_{N-1})
	buffer_fixedLegPaymentIndexSchedule_	= buffer_UnderlyingSwap_.get_fixedLegPaymentIndexSchedule() ;
	buffer_deltaTFixedLeg_					= buffer_UnderlyingSwap_.get_DeltaTFixedLeg() ;

	initialize_y_bar(buffer_T0_, gridSize);	
	//virtual : initialisation de b barre ou b barre et c barre dans les classes dérivées
	buffer_s0_								= swapRate0();     //à reinitialiser en dernier si bug
}

//pour la calibration, updateVol permet mise à jour des buffers lorsque les parametres sont modifiés
void Cheyette_SwaptionPricer_Approx::updateLevel_calib(double newValue, size_t index)
{
	//o << "level : ;" << newValue ;
	pCheyette_Model_->setCheyetteModel_Parameter_Level(newValue, index) ;
	preCalculateModelBuffers() ; 
}

void Cheyette_SwaptionPricer_Approx::updateSkew_calib(double newValue, size_t index)
{
	//o << "skew : ;" << newValue ;
	pCheyette_Model_->setCheyetteModel_Parameter_Skew(newValue, index) ;
	preCalculateModelBuffers() ; 
}

void Cheyette_SwaptionPricer_Approx::updateConvexity_calib(double newValue, size_t index)
{
	//o << "convexity : ;" << a ;
	pCheyette_Model_->setCheyetteModel_Parameter_Convexity(newValue, index) ;
	preCalculateModelBuffers() ; 
}

//fonction qui intervient pour le calcul de y bar
//sera integrée sur s entre 0 et t 
double Cheyette_SwaptionPricer_Approx::to_integrate_y_bar(double s) const
{
//ne convient que pour DD et Quad, k supposé constant et ne pas dépendre de t

	double k			= pCheyette_Model_->meanReversion(0., 0., 0.) ; //   get_CheyetteDD_Parameter().k_ ;
	double sigma_r_0	= pCheyette_Model_->localVol(s, 0., 0.) ;	//sigma_r(s, 0., 0.) ;	//sigma_r^0(t) = sigma_r(t, 0, 0) 
	return exp(2 * k * s) * sigma_r_0 * sigma_r_0 ; 
}

//integrale jusqu'à t 
// y_bar(t) = exp( - 2 k t ) * integrale_0^t( exp(2 k s) sigma_r^0(s) sigma_r^0(s) ds ) 
void Cheyette_SwaptionPricer_Approx::initialize_y_bar(double t, size_t gridSize) const
{
	assert(t > 0);
	double gridStart	= 0. ;
	double gridEnd		= t ;
	double k			= pCheyette_Model_->meanReversion(0., 0., 0.) ; //pCheyetteDD_Model_->get_CheyetteDD_Parameter().k_ ;

	//constructeur pour le schema d'integration
	numeric::IncrementalIntegrator1D_Riemann integral(gridStart, gridEnd, gridSize);
	//constructeur pour la fonction à intégrer
	boost::function<double(double)> f = 
				boost::bind(&Cheyette_SwaptionPricer_Approx::to_integrate_y_bar, boost::ref(*this), _1);

	integral.vecteur_integrate(f) ;

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


double Cheyette_SwaptionPricer_Approx::ZC_1stDerivative_on_xt(double t, double T, double x_t, double y_t) const
{
	assert(0 <= T && T <= buffer_TN_) ;
	if (t==0){assert(x_t == 0 && y_t == 0) ;}

	double ZC = pCheyette_Model_->P(t, T, x_t, y_t) ;	
	double g = pCheyette_Model_->G(t,T) ;
	double res = - g * ZC ;
	return res ;	
}

double Cheyette_SwaptionPricer_Approx::ZC_2ndDerivative_on_xt(double t, double T, double x_t, double y_t) const
{
	assert(0 <= T && T <= buffer_TN_) ;

	double ZC = pCheyette_Model_->P(t, T, x_t, y_t) ;
	double g = pCheyette_Model_->G(t,T) ;
	double res = g * g * ZC ;
	
	return res ;					
}

// Numerator = P(t, T0) - P(t, TN)
double Cheyette_SwaptionPricer_Approx::swapRateNumerator(double t, double x_t, double y_t) const 
{
	assert( t >= 0. ) ;

	double ZC_T0, ZC_TN ;
	ZC_T0 = pCheyette_Model_->P(t, buffer_T0_, x_t, y_t) ; 
	ZC_TN = pCheyette_Model_->P(t, buffer_TN_, x_t, y_t) ;	
	
	return  ZC_T0 - ZC_TN ; 
}

double Cheyette_SwaptionPricer_Approx::swapRateNumerator_1stDerivative(double t, double x_t, double y_t) const
{
	return ZC_1stDerivative_on_xt(t, buffer_T0_, x_t, y_t) - ZC_1stDerivative_on_xt(t, buffer_TN_, x_t, y_t) ;
}

double Cheyette_SwaptionPricer_Approx::swapRateNumerator_2ndDerivative(double t, double x_t, double y_t) const
{
	return ZC_2ndDerivative_on_xt(t, buffer_T0_, x_t, y_t) - ZC_2ndDerivative_on_xt(t, buffer_TN_, x_t, y_t) ;
}

// Denominator = \sum delta_k P(t,T_k) en t = 0
double Cheyette_SwaptionPricer_Approx::swapRateDenominator(double t, double x_t, double y_t) const 
{	
	double fixed_tenor = buffer_UnderlyingSwap_.get_fixedLegTenorType().YearFraction() ;
	double float_tenor = buffer_UnderlyingSwap_.get_floatingLegTenorType().YearFraction() ;
	double tenor_ref = std::min(fixed_tenor, float_tenor) ;  //le plus petit 

	double denom = 0.0;
	for(size_t itr = 0 ; itr < buffer_fixedLegPaymentIndexSchedule_.size() ; ++itr) 
	{
		double dateEchangeFluxFixe = buffer_fixedLegPaymentIndexSchedule_[itr] * tenor_ref ;
		double delta_T = buffer_deltaTFixedLeg_[itr] ;	
		double ZC		= pCheyette_Model_->P(t, dateEchangeFluxFixe, x_t, y_t) ;
		denom += delta_T * ZC ;
	}
	
	return denom ;
}

double Cheyette_SwaptionPricer_Approx::swapRateDenominator_1stDerivative(double t, double x_t, double y_t) const
{
	double result = 0. ;

	double fixed_tenor = buffer_UnderlyingSwap_.get_fixedLegTenorType().YearFraction() ;
	double float_tenor = buffer_UnderlyingSwap_.get_floatingLegTenorType().YearFraction() ;
	double tenor_ref = std::min(fixed_tenor, float_tenor) ;  

	for(size_t itr = 0 ; itr < buffer_fixedLegPaymentIndexSchedule_.size() ; ++itr) 
	{
		double dateEchangeFluxFixe = buffer_fixedLegPaymentIndexSchedule_[itr] * tenor_ref ;
		result += buffer_deltaTFixedLeg_[itr] * ZC_1stDerivative_on_xt(t, dateEchangeFluxFixe, x_t, y_t) ;	
	}
	return result;
}

double Cheyette_SwaptionPricer_Approx::swapRateDenominator_2ndDerivative(double t, double x_t, double y_t) const
{
	double result = 0. ;

	double fixed_tenor = buffer_UnderlyingSwap_.get_fixedLegTenorType().YearFraction() ;
	double float_tenor = buffer_UnderlyingSwap_.get_floatingLegTenorType().YearFraction() ;
	double tenor_ref = std::min(fixed_tenor, float_tenor) ;  

	for(size_t itr = 0; itr < buffer_fixedLegPaymentIndexSchedule_.size(); ++itr) 
	{
		double dateEchangeFluxFixe = buffer_fixedLegPaymentIndexSchedule_[itr] * tenor_ref ;
		result += buffer_deltaTFixedLeg_[itr] * ZC_2ndDerivative_on_xt(t, dateEchangeFluxFixe, x_t, y_t) ;	
	}
	return result;
}

double Cheyette_SwaptionPricer_Approx::swapRate0() const
{
	return swapRate(0., 0., 0.);
}

double Cheyette_SwaptionPricer_Approx::swapRate(double t, double x_t, double y_t) const
{
	double n      = swapRateNumerator(t, x_t, y_t);
	double d      = swapRateDenominator(t, x_t, y_t);

	assert (d > 0.001) ;
	return n/d;
}

double Cheyette_SwaptionPricer_Approx::swapRate_1stDerivative(double t, double x_t, double y_t) const
{
	double n   = swapRateNumerator(t, x_t, y_t);
	double n_1 = swapRateNumerator_1stDerivative(t, x_t, y_t); 

	double d   = swapRateDenominator(t, x_t, y_t);
	double d_1 = swapRateDenominator_1stDerivative(t, x_t, y_t);

	return  (n_1*d - n*d_1)/(d*d) ;	
}

double Cheyette_SwaptionPricer_Approx::swapRate0_1stDerivative(double t, double x_t, double y_t) const
{
	return 0. ;
}

double Cheyette_SwaptionPricer_Approx::swapRate_2ndDerivative(double t, double x_t, double y_t) const
{
	double n	= swapRateNumerator(t, x_t, y_t);
	double n_1	= swapRateNumerator_1stDerivative(t, x_t, y_t); 
	double n_2	= swapRateNumerator_2ndDerivative(t, x_t, y_t); 

	double d	= swapRateDenominator(t, x_t, y_t);
	double d_1	= swapRateDenominator_1stDerivative(t, x_t, y_t);
	double d_2	= swapRateDenominator_2ndDerivative(t, x_t, y_t);

	double result = (n_2 * d - n * d_2) * d - (n_1 * d - n * d_1) * 2 * d_1 ;
	result /= (d*d*d);

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

	//callable function object that accepts one parameter and returns a tuple:
	//the tuple should have two elements containing the evaluation of the function and it's first derivative.
	boost::math::tuple<double, double> operator()(double x)
	{
		return boost::math::make_tuple(f_(x) - target_, f_derivative_(x));
	}
};

//S(t, x_t) = swapRate(t, x_t) = s 
//retourne le x_t correspondant
double Cheyette_SwaptionPricer_Approx::inverse(double t, double s) const
{                                                                  //swapRate(double t, double x_t)
	double y_bar = buffer_y_bar_(t) ;
	boost::function<double(double)> f				= 
				boost::bind(&Cheyette_SwaptionPricer_Approx::swapRate, this, t, _1, y_bar);
	boost::function<double(double)> f_derivative	= 
				boost::bind(&Cheyette_SwaptionPricer_Approx::swapRate_1stDerivative, this, t, _1, y_bar);
	
	Newton_Raphson_struct NR_struct(f, f_derivative, s);   //struct qui contient f, f', target = s
	double initial_guess	= 0.0 ;
	double min				= -5.0;
	double max				= 5.0;
    size_t nDigits			= 15;
//	boost::uintmax_t nMaxIter  = 50 ;

	double res_newton_raphson = boost::math::tools::newton_raphson_iterate(NR_struct, initial_guess, min, max, nDigits);
	
	return res_newton_raphson;
}



double Cheyette_SwaptionPricer_Approx::volBlack(const VanillaSwaption_PTR pSwaption) const 
{
	double approxPrice = price(pSwaption) ;		//contient	this->preCalculateALL(pSwaption);
	
	//std::cout << "DEBUG vol imp correspondant au prix. Prix : " << approxPrice << std::endl ;

	const VanillaSwap& vanillaSwap = pSwaption->getUnderlyingSwap();

	double T			= pSwaption->getUnderlyingSwap().get_StartDate() ;	//1ere date de fixing
	double strike		= vanillaSwap.get_strike() ;
	double S0			= get_buffer_s0_() ;
	double annuity0		= swapRateDenominator(0., 0., 0.) ;

	double volBlack		= NumericalMethods::Black_impliedVolatility(approxPrice / annuity0, S0, strike, T) ;

	//std::cout << "vol Black : " << volBlack << std::endl ;

	return volBlack;
}


//pricing de swaption pour des strikes donnés
	//res[0] = prixApprox_SwaptionsPlusieursStrikes ;
	//res[1] = strikes ; 
	//res[2] = vol implicite ; 
// la swaption passée en parametre devrait être la swaption ATM
std::vector<std::vector<double>> Cheyette_SwaptionPricer_Approx::priceMultipleStrikes(VanillaSwaption_PTR pSwaption, 
																					std::vector<double> shifts_bp)
{	
	size_t nbShifts = shifts_bp.size() ;
//strikes equivalents
	std::vector<double> strikes(nbShifts) ;	
	std::vector<double> prixApproxPlusieursStrikes(nbShifts) ;
	std::vector<double> volImpBlack(nbShifts) ;
	double strike_0 = pSwaption->get_strike() ;  //strike initial de la swaption en parametre

	for (size_t i = 0 ; i < nbShifts ; ++i)
	{
		double strike = strike_0 + shifts_bp[i] / 10000. ;
		strikes[i] = strike ; 

		pSwaption->getUnderlyingSwap_RefNonConst().set_strike(strikes[i]) ;

		double approx = price(pSwaption) ;	//preCalculateAll() dans price()
		prixApproxPlusieursStrikes[i] = approx ;

		volImpBlack[i] = volBlack(pSwaption) ;   //pas du tout optimal, recalcul du prix approx de la swaption en double
		//surcharger la fonction volBlack
	}
//remise strike ATM
	pSwaption->getUnderlyingSwap_RefNonConst().set_strike(strike_0) ;
	
	std::vector<std::vector<double>> res(3) ;
	res[0] = prixApproxPlusieursStrikes ;
	res[1] = strikes ; 
	res[2] = volImpBlack ;
	return res ;
}

//pricing de swaption pour des strikes correspondant à une standardized moneyness dans [-5, 5]
	//res[0] = prixApprox_SwaptionsPlusieursStrikes ;
	//res[1] = strikes ; 
	//res[2] = moneyness ;
//std::vector<std::vector<double>> Cheyette_SwaptionPricer_Approx::priceMultipleStrikes(double sigma_ATM)
//{
//
//	//standardized moneyness
//	size_t nbMoneyness = 11 ;		//moneyness = -5, -4, ... , 0, 1, ... 5
//	std::vector<double> moneyness(nbMoneyness) ;			
//	moneyness[0] = 5. ; moneyness[1] = 4. ; moneyness[2] = 3. ; moneyness[3] = 2. ; moneyness[4] = 1. ;
//	moneyness[5] = 0. ;
//	moneyness[6] = -1. ; moneyness[7] = -2. ; moneyness[8] = -3. ; moneyness[9] = -4. ; moneyness[10] = -5. ;
//	
//	//strike equivalent pour une standardized moneyness dans [-5 ; 5]
//	std::vector<double> strikes(nbMoneyness) ;	
//	std::vector<double> prixApproxPlusieursStrikes(nbMoneyness) ;
//	double T0 = pSwaption_->getUnderlyingSwap().get_StartDate() ;
//
//	for (size_t i = 0 ; i < nbMoneyness ; ++i)
//	{
//		double strike = buffer_s0_ / exp(sigma_ATM * sqrt(T0) * moneyness[i]) ; //   (strikeATM_Bloomberg + shiftStrike[i])/100. ;
//		strikes[i] = strike ; 
//
//		pSwaption_->getUnderlyingSwap_RefNonConst().set_strike(strikes[i]) ;
//		initialize_buffers() ; 
//
//		double approx = prixSwaptionApprox() ;	
//		prixApproxPlusieursStrikes[i] = approx ;
//	}
//
//	std::vector<double> prixSwaptionsPlusieursStrikes(nbMoneyness) ;
//
//	std::vector<std::vector<double>> res(3) ;
//	res[0] = prixApproxPlusieursStrikes ;
//	res[1] = strikes ; 
//	res[2] = moneyness ; 
//
//	return res ;
//}
