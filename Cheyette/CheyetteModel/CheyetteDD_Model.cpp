#include "CheyetteDD_Model.h"


void CheyetteDD_Model::show() const
{
	std::cout << "---------------------------------------------" << std::endl ;
	std::cout << "--- creation d'un objet Cheyette DD Model ---" << std::endl ;
	std::cout << "  SIGMA piecewise constant " << std::endl ;
	cheyetteDD_Parameter_.sigma_.show() ;
	std::cout << "  M piecewise constant " << std::endl ;
	cheyetteDD_Parameter_.m_.show() ;
	//std::cout << "  SHIFT piecewise constant " << std::endl ;
	//cheyetteDD_Parameter_.shift_.show() ;
	std::cout << "  k constant " << std::endl ;
	std::cout << "       k = " << cheyetteDD_Parameter_.k_ << std::endl ;
	std::cout << "---------------------------------------------" << std::endl ;

}

void CheyetteDD_Model::print(std::ostream& o) const 
{
	o << "---------------------------------------------" << std::endl ;
	o << "------------- Cheyette DD Model ------------" << std::endl ;
	o << "  SIGMA piecewise constant " << std::endl ;
	cheyetteDD_Parameter_.sigma_.print(o) ;
	o << "  M piecewise constant " << std::endl ;
	cheyetteDD_Parameter_.m_.print(o) ;
	o << "  k constant " << std::endl ;
	o << "       k = " << cheyetteDD_Parameter_.k_ << std::endl ;
	o << "---------------------------------------------" << std::endl ;
}


	//x = x_t  !!
//fonction de vol locale Displaced Diffusion
double CheyetteDD_Model::sigma_r( double t,  double x_t) const 
{
	double sigma_t = cheyetteDD_Parameter_.sigma_(t) ;
	double m_t		= cheyetteDD_Parameter_.m_(t) ;
	double shift1	= shift1_->operator()(t, x_t) ;
	double shift2	= shift2_->operator()(t, x_t) ;  //en pratique pointe vers shift1(0, 0)
	return sigma_t * (m_t * shift1 + (1 - m_t) * shift2) ;
}

double CheyetteDD_Model::sigma_r_t_1stDerivative( double t,  double x_t) const 
{
	double sigma_t = cheyetteDD_Parameter_.sigma_(t) ;
	double m_t		= cheyetteDD_Parameter_.m_(t) ;
	return sigma_t * (m_t * derivative_x_shift1_->operator()(t, x_t) 
					+ (1 - m_t) * derivative_x_shift2_->operator()(t, x_t) )  ;
}

//shift functions
double CheyetteDD_Model::shift_rt(double t,  double x_t) const
{
	return x_t + courbeInput_PTR_->get_f_0_t(t) ;  
}	

double CheyetteDD_Model::shift_r0(double t,  double x_t) const 
{
	return courbeInput_PTR_->get_f_0_t(0) ;  
}	

double CheyetteDD_Model::shift_f_0_t(double t,  double x_t) const
{
	return courbeInput_PTR_->get_f_0_t(t) ;  
}				

//derivee du shift wrt x_t
double CheyetteDD_Model::shift_rt_1stDerivative(double t,  double x_t) const
{
	return 1 ;
}
double CheyetteDD_Model::shift_r0_1stDerivative(double t,  double x_t) const
{
	return 0 ;
} 
double CheyetteDD_Model::shift_f_0_t_1stDerivative(double t,  double x_t) const
{
	return 0 ;
}


//fonctions G(t, T), ZC B(t, T)...
//Displaced Diffusion : k constant, calcul explicite de G(t, T)
double CheyetteDD_Model::G(double t, double T) const
{
	assert(t >= 0);  
	double k = cheyetteDD_Parameter_.k_ ;
	//std::cout << " t : " << t << ", T : " << T << ", G(t, T) : " << 1/k * (1 - exp(- k * (T-t))) << std::endl ;
	return 1/k * (1 - exp(- k * (T-t))) ;
}

double CheyetteDD_Model::P(double t, double T, double x_t, double y_t) const 
{
	double g = G(t,T) ;
	double P_0_t = exp( - courbeInput_PTR_->get_tauxZC0(t) * t ) ;
	double P_0_T = exp( - courbeInput_PTR_->get_tauxZC0(T) * T ) ;
	return P_0_T/P_0_t * exp(- x_t * g - 0.5 * y_t * g * g) ; 
}

//f_0_t à coder dans la classe courbeInput
double CheyetteDD_Model::r_t(double t, double x_t) const 
{
	double f_0_t = courbeInput_PTR_->get_f_0_t(t) ;
	return f_0_t + x_t ; 
}

//T2 = T1 + delta
double CheyetteDD_Model::Libor(double t, double T1, double T2, double x_t, double y_t) const
{
	double delta = T2 - T1 ;
	double P_t_T1 = P(t, T1, x_t, y_t) ;
	double P_t_T2 = P(t, T2, x_t, y_t) ;
	return 1/delta * (P_t_T1/P_t_T2 - 1.0) ;
}

//EDS : drift et diffusion sous Q^T
double CheyetteDD_Model::diffusion_x(double t, double x_t) const
{
	return sigma_r(t, x_t) ;
}

double CheyetteDD_Model::drift_x_QT(double t, double T_proba_fwd, double x_t, double y_t) const
{
	double k			= cheyetteDD_Parameter_.k_ ;
//	double r_t			= x_t + courbeInput_PTR_->get_f_0_t(t) ;
	double sigma_r_t	= sigma_r(t, x_t) ;

	return y_t - k * x_t - G(t, T_proba_fwd) * sigma_r_t * sigma_r_t ;	  
}


double CheyetteDD_Model::drift_y(double t, double x_t, double y_t) const
{
	double sigma_r_t	=  sigma_r(t, x_t) ;
	double k			= cheyetteDD_Parameter_.k_ ;
	return sigma_r_t * sigma_r_t - 2 * k * y_t ;
}
