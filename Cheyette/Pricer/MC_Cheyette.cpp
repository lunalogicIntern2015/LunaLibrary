#include <Cheyette\Pricer\MC_Cheyette.h>


void MC_Cheyette::simulate_Euler() const
{
	//initialisation
	x_t_Cheyette_[0] = 0. ;	
	y_t_Cheyette_[0] = 0. ;

	double tenorYearFrac = pTenorStructure_->get_tenorType().YearFraction() ;	
	double x_t = 0. ;
	double y_t = 0. ;
	double t = 0. ;

	size_t fwdProbaIndex = static_cast<size_t>(fwdProbaT_ /tenorYearFrac) ;
	//simulation de x_t, y_t jusqu'à fwdProbaT
	for (size_t index = 1; index <= fwdProbaIndex ; ++index)    
	{ 
		std::vector<double>   gaussian_tmp(discretizationBetweenDates_);  
		rnGenerator_->generate(gaussian_tmp);		// generate Gaussian.
		double dt = tenorYearFrac / discretizationBetweenDates_ ;

		for (size_t pasDiscretisation = 1 ; pasDiscretisation <= discretizationBetweenDates_ ; ++pasDiscretisation)    
		{
			double t_plus_dt = t + dt ;
			double x_t_plus_dt =  x_t	+ cheyetteModel_PTR_->drift_x_QT(t, fwdProbaT_, x_t, y_t) * dt 
							+ cheyetteModel_PTR_->diffusion_x(t, x_t, y_t) * sqrt(dt) * gaussian_tmp[pasDiscretisation - 1] ;
						
			double y_t_plus_dt =  y_t	+ cheyetteModel_PTR_->drift_y(t, x_t, y_t) * dt ;

			x_t = x_t_plus_dt ;
			y_t = y_t_plus_dt ;
			t = t_plus_dt ;
		}
		x_t_Cheyette_[index] = x_t ; 
		y_t_Cheyette_[index] = y_t ;
	}

	computeNumeraires() ;
}



//numeraires_[i] = B(i * tenor.yearFraction()   , fwdProbaT_)
void MC_Cheyette::computeNumeraires() const
{
	double tenorYearFrac = pTenorStructure_->get_tenorType().YearFraction() ;	
	size_t fwdProbaIndex = static_cast<size_t>(fwdProbaT_ /tenorYearFrac) ;
	for (size_t i = 0 ; i <= fwdProbaIndex ; ++i)
	{
		double x_t = x_t_Cheyette_[i] ;
		double y_t = y_t_Cheyette_[i] ;
		numeraires_[i] = cheyetteModel_PTR_->P(i * tenorYearFrac, fwdProbaT_, x_t, y_t);
	}
}




