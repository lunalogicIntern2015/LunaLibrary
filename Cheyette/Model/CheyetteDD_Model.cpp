#include "CheyetteDD_Model.h"

//Displaced Diffusion : k constant, calcul explicite de G(t, T)
double CheyetteDD_Model::G(double t, double T) const
{
	assert(t >= 0);  	
	return 1/k_ * (1 - exp(- k_ * (T-t))) ;
}


void CheyetteDD_Model::show() const
{
	CheyetteModel::show() ;

	std::cout << this->getModelType() << std::endl ;

	std::cout << "  SIGMA piecewise constant " << std::endl ;
	sigma_.show() ;
	std::cout << "  M piecewise constant " << std::endl ;
	m_.show() ;
	std::cout << "  k constant " << std::endl ;
	std::cout << "       k = " << k_ << std::endl ;
	std::cout << "  COURBE YIELD " << std::endl ;
	courbeInput_PTR_->show() ;
	std::cout << "---------------------------------------------" << std::endl ;
}

void CheyetteDD_Model_SwapRateVersion::show() const
{
	CheyetteDD_Model::show() ;

	pSwaption_->show() ;
}

void CheyetteDD_Model::print(std::ostream& o) const 
{
	CheyetteModel::print(o) ;

	o << this->getModelType() << std::endl ;

	o << "  SIGMA piecewise constant " << std::endl ;
	sigma_.print(o) ;
	o << "  M piecewise constant " << std::endl ;
	m_.print(o) ;
	o << "  k constant " << std::endl ;
	o << "       k = ; " << k_ << std::endl ;
	o << "---------------------------------------------" << std::endl ;
}


//drift, diffusion mis dans CheyetteModel

	//EDS : drift et diffusion sous Q^T
	//double CheyetteDD_Model::drift_x_Q(double t, double x_t, double y_t) const                       // risk neutral proba
	//{
	//	double k	= cheyetteDD_Parameter_.k_ ;
	//	return y_t - k * x_t;	  	 
	//}
	//double CheyetteDD_Model::drift_x_QT(double t, double T_proba_fwd, double x_t, double y_t) const  // forward proba
	//{
	//	double sigma_r_t	= sigma_r(t, x_t, y_t) ;
	//
	//	return drift_x_Q(t,x_t,y_t) - G(t, T_proba_fwd) * sigma_r_t * sigma_r_t ;	  	  
	//}
	//double CheyetteDD_Model::diffusion_x(double t, double x_t, double y_t) const
	//{
	//	return sigma_r(t, x_t, y_t) ;
	//}

	//double CheyetteDD_Model::drift_y(double t, double x_t, double y_t) const
	//{
	//	double sigma_r_t	=  sigma_r(t, x_t, y_t) ;
	//	double k			= cheyetteDD_Parameter_.k_ ;
	//	return sigma_r_t * sigma_r_t - 2. * k * y_t ;
	//}
