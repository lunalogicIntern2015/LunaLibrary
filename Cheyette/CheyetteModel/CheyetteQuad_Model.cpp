#include "CheyetteQuad_Model.h"

void CheyetteQuad_Model::show() const
{
	std::cout << "--------------------------------------------------" << std::endl ;
	std::cout << "--- creation d'un objet Cheyette Quadratic Vol ---" << std::endl ;
	std::cout << "  A piecewise constant " << std::endl ;
	cheyetteQuad_Parameter_.a_.show() ;
	std::cout << "  B piecewise constant " << std::endl ;
	cheyetteQuad_Parameter_.b_.show() ;
	std::cout << "  C piecewise constant " << std::endl ;
	cheyetteQuad_Parameter_.c_.show() ;

	std::cout << "  k constant " << std::endl ;
	std::cout << "       k = " << cheyetteQuad_Parameter_.k_ << std::endl ;
	std::cout << "---------------------------------------------" << std::endl ;

}

void CheyetteQuad_Model::print(std::ostream& o) const 
{
	o << "-------------------------------------------------" << std::endl ;
	o << "------------- Cheyette Quadratic Vol ------------" << std::endl ;
	o << "  A piecewise constant " << std::endl ;
	cheyetteQuad_Parameter_.a_.print(o) ;
	o << "  B piecewise constant " << std::endl ;
	cheyetteQuad_Parameter_.b_.print(o) ;
	o << "  C piecewise constant " << std::endl ;
	cheyetteQuad_Parameter_.c_.print(o) ;
	o << "  k constant " << std::endl ;
	o << "       k = ; " << cheyetteQuad_Parameter_.k_ << std::endl ;
	o << "---------------------------------------------" << std::endl ;
}

	//x = x_t  !!
//fonction de vol locale Quadratic Volatility
double CheyetteQuad_Model::sigma_r( double t,  double x_t) const 
{
	double a_t		= cheyetteQuad_Parameter_.a_(t) ;
	double b_t		= cheyetteQuad_Parameter_.b_(t) ;
	double c_t		= cheyetteQuad_Parameter_.c_(t) ;

	return a_t + b_t * x_t + c_t * x_t * x_t ;
}

double CheyetteQuad_Model::sigma_r_t_1stDerivative( double t,  double x_t) const 
{
	double b_t		= cheyetteQuad_Parameter_.b_(t) ;
	double c_t		= cheyetteQuad_Parameter_.c_(t) ;

	return b_t + 2 * c_t * x_t ;
}

//fonctions G(t, T), ZC B(t, T)...
//k constant, calcul explicite de G(t, T)
double CheyetteQuad_Model::G(double t, double T) const
{
	assert(t >= 0);  
	double k = cheyetteQuad_Parameter_.k_ ;

	return 1/k * (1 - exp(- k * (T-t))) ;
}

double CheyetteQuad_Model::P(double t, double T, double x_t, double y_t) const 
{
	double g = G(t,T) ;
	double P_0_t = exp( - courbeInput_PTR_->get_tauxZC0(t) * t ) ;
	double P_0_T = exp( - courbeInput_PTR_->get_tauxZC0(T) * T ) ;
	return P_0_T/P_0_t * exp(- x_t * g - 0.5 * y_t * g * g) ; 
}

//f_0_t à coder dans la classe courbeInput
double CheyetteQuad_Model::r_t(double t, double x_t) const 
{
	double f_0_t = courbeInput_PTR_->get_f_0_t(t) ;
	return f_0_t + x_t ; 
}

//T2 = T1 + delta
double CheyetteQuad_Model::Libor(double t, double T1, double T2, double x_t, double y_t) const
{
	double delta = T2 - T1 ;
	double P_t_T1 = P(t, T1, x_t, y_t) ;
	double P_t_T2 = P(t, T2, x_t, y_t) ;
	return 1/delta * (P_t_T1/P_t_T2 - 1.0) ;
}



//EDS : drift et diffusion sous Q^T
double CheyetteQuad_Model::diffusion_x(double t, double x_t) const
{
	return sigma_r(t, x_t) ;
}

double CheyetteQuad_Model::drift_x_QT(double t, double T_proba_fwd, double x_t, double y_t) const
{
	double k			= cheyetteQuad_Parameter_.k_ ;
	double sigma_r_t	= sigma_r(t, x_t) ;

	return y_t - k * x_t - G(t, T_proba_fwd) * sigma_r_t * sigma_r_t ;	  
}


double CheyetteQuad_Model::drift_y(double t, double x_t, double y_t) const
{
	double sigma_r_t	=  sigma_r(t, x_t) ;
	double k			= cheyetteQuad_Parameter_.k_ ;
	return sigma_r_t * sigma_r_t - 2 * k * y_t ;
}
