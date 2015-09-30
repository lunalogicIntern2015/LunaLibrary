#include <Cheyette/Pricer/Cheyette_SwaptionPricer_LinearApprox.h>

const size_t gridSize  = 101 ;	//mettre 11 ou 101 (11 est suffisant)


//void Cheyette_SwaptionPricer_LinearApprox::initialize_buffers() 
//{
//	Cheyette_SwaptionPricer_Approx::initialize_buffers() ;
//	buffer_b_barre_ = timeAverage(buffer_T0_) ;
//}

void Cheyette_SwaptionPricer_LinearApprox::preCalculateALL(VanillaSwaption_PTR pSwaption) const
{
	Cheyette_SwaptionPricer_Approx::preCalculateALL(pSwaption) ;
	buffer_b_barre_ = timeAverage(buffer_T0_) ;
}

void Cheyette_SwaptionPricer_LinearApprox::preCalculateModelBuffers() const
{
	buffer_b_barre_ = timeAverage(buffer_T0_) ;
}


/***************** Linear approx of swap rate volatility *********************
*
*	dS(t) = Phi(t, S_t) dW^Q_A
*	dS(t) \approx Phi(t, s_bar) + DPhi(t, s_bar) (S_t - s_bar) dW^Q_A
*
******************************************************************************/

double Cheyette_SwaptionPricer_LinearApprox::phi(double t, double x_t) const
{
	double y_bar = buffer_y_bar_(t) ;

	double sigma_r	= pCheyette_Model_->localVol(t, x_t, y_bar) ;   //sigma_r(t, x_t, y_bar) ;	
	double phi_res = swapRate_1stDerivative(t, x_t, y_bar) * sigma_r ;

	return phi_res ;
}

double Cheyette_SwaptionPricer_LinearApprox::dPhi(double t, double x_t) const
{
	double y_bar = buffer_y_bar_(t) ;

	double sigma_r	= pCheyette_Model_->localVol(t, x_t, y_bar) ; 
	double res		=  sigma_r *	swapRate_2ndDerivative(t, x_t, y_bar) / swapRate_1stDerivative(t, x_t, y_bar) 
					+ pCheyette_Model_->localVol_1stDerivative(t, x_t, y_bar) ;

	return res ; 
}


/****************  parameter averaging  *********************
*
*	dS(t) \approx Phi(t, s_bar) + DPhi(t, s_bar) (S_t - s_bar) dW^Q_A
*	dS(t) = \lambda(t) ( S0 + b(t) (S(t) - S(0) ) dW^Q_A		-> time averaging: b(t) puis b_barre
*
*	Phi = lambda(t) S0
*	DPhi = lambda(t) b(t)
*
************************************************************/

//lambda(t) 
double Cheyette_SwaptionPricer_LinearApprox::lambda(double t) const
{
	double inverse_s_bar =	inverse(t, buffer_s0_) ;
	return phi(t, inverse_s_bar) / buffer_s0_ ;
}	

double Cheyette_SwaptionPricer_LinearApprox::lambda2(double t) const
{
	double l = lambda(t) ;
	double l2 = l * l ;
	return l2 ;
}	

//b(t) 
double Cheyette_SwaptionPricer_LinearApprox::b(double t) const
{
	double inverse_s_bar =	inverse(t, buffer_s0_) ;
	return dPhi(t, inverse_s_bar) / phi(t, inverse_s_bar) * buffer_s0_ ;
}	


double Cheyette_SwaptionPricer_LinearApprox::f_outer_num(double t) const
{
	return lambda2(t) * b(t) ;  //OK (fonction outer seulement ie sans v^2(u)
}
double Cheyette_SwaptionPricer_LinearApprox::f_outer_denom(double t) const
{
	return lambda2(t)  ;		//OK (fonction outer seulement ie sans v^2(u)
}


//average over [0, t]
//gridSize : nb de points pour integrale Riemann
//gridSize + 1 : pour 100 mettre 101, delta_t = 1/100
double Cheyette_SwaptionPricer_LinearApprox::timeAverage(double t) const	
{
	double gridStart = 0.0;
	double gridEnd = t ; 

//integrale numerateur

	numeric::IncrementalIntegrator2D_Riemann int2D(gridStart, gridEnd, gridSize) ;
	boost::function<double(double)> f_inner = boost::bind(&Cheyette_SwaptionPricer_LinearApprox::lambda2, *this, _1);
	boost::function<double(double)> f_outer = boost::bind(&Cheyette_SwaptionPricer_LinearApprox::f_outer_num, *this, _1);
	
	double integrale_numerateur = int2D.integrate(f_outer, f_inner) ;
	
//integrale denominateur
	boost::function<double(double)> f_denom = boost::bind(&Cheyette_SwaptionPricer_LinearApprox::f_outer_denom, *this, _1);
	
	double integrale_denom = int2D.integrate(f_denom, f_inner);
	double b_barre = integrale_numerateur/integrale_denom ;
	//buffer_b_barre_ = b_barre ;
	return b_barre ;
}



double Cheyette_SwaptionPricer_LinearApprox::price(const VanillaSwaption_PTR pSwaption) const
{
	if(boost::dynamic_pointer_cast<CheyetteDD_Model_SwapRateVersion>(pCheyette_Model_))
	{
		CheyetteDD_Model_SwapRateVersion_PTR model	=	
			boost::dynamic_pointer_cast<CheyetteDD_Model_SwapRateVersion>(pCheyette_Model_);
		model->setSwaption(pSwaption) ;
	}
//calcul des buffers
	preCalculateALL(pSwaption) ;
	
	//calcul de la variance (integrale de lambda_t)
	double gridStart = 0.0;
	double gridEnd = buffer_T0_ ;
	numeric::Integrator1D_Riemann integral_Riemann(gridStart, gridEnd, gridSize) ;

	boost::function<double(double)> func1 = boost::bind(&Cheyette_SwaptionPricer_LinearApprox::lambda2, this, _1);
	double integrale = integral_Riemann.integrate(func1);
	
	double annuity0 = swapRateDenominator(0., 0., 0.) ;	
	double K_tilde	= buffer_b_barre_ * buffer_UnderlyingSwap_.get_strike() + (1 - buffer_b_barre_) * buffer_s0_ ;
	double sigma_sqrt_T = sqrt(integrale * buffer_b_barre_ * buffer_b_barre_) ;  //b_barre * sqrt(integrale) 

	//prend en compte strikes positifs et négatifs
	double prixBlack = NumericalMethods::Black_Price_vol2_allStrike(buffer_s0_, K_tilde, sigma_sqrt_T, buffer_T0_) ;

	return annuity0 / buffer_b_barre_ * prixBlack ;
	
}
