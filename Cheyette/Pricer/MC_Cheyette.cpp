#include <Cheyette\Pricer\MC_Cheyette.h>

void MC_Cheyette::simulate_Euler() const
{
	
	//initialisation
	x_t_Cheyette_[0] = 0  ;	y_t_Cheyette_[0] = 0;

	double tenorYearFrac = tenorType_.YearFraction() ;
	double dt, t, t_plus_dt, x_t, x_t_plus_dt, y_t, y_t_plus_dt ;	
	x_t = 0 ; y_t = 0 ; t = 0 ;

	for (size_t date = 0; date < indexOfSimulation_.size() ; ++date)    
	{ 
		std::vector<double>   gaussian_tmp(discretizationBetweenDates_[date]);  
		rnGenerator_->generate(gaussian_tmp);		// generate Gaussian.
		if (date == 0)
			{dt = indexOfSimulation_[date] * tenorYearFrac / discretizationBetweenDates_[date] ;}
		else
			{dt = (indexOfSimulation_[date] - indexOfSimulation_[date-1]) * tenorYearFrac / 
																		discretizationBetweenDates_[date] ;}

		t_plus_dt = t ;
		////DEBUG
		//for (size_t i = 0 ; i < gaussian_tmp.size() ; ++i)    
		//{			
		//	std::cout << "gaussian_tmp " << i << " : " << gaussian_tmp[i] << std::endl ;
		//}

		for (size_t pasDiscretisation = 1 ; pasDiscretisation <= discretizationBetweenDates_[date] ; ++pasDiscretisation)    
		{
			t_plus_dt += dt ;

			x_t_plus_dt =  x_t	+ cheyetteDD_Model_->drift_x_QT(t, fwdProbaT_, x_t, y_t) * dt 
								+ cheyetteDD_Model_->diffusion_x(t, x_t) * sqrt(dt) * gaussian_tmp[pasDiscretisation - 1] ;
			y_t_plus_dt =  y_t	+ cheyetteDD_Model_->drift_y(t, x_t, y_t) * dt ;
			x_t = x_t_plus_dt ;
			y_t = y_t_plus_dt ;
			t = t_plus_dt ;
		}
//		std::cout << "date t : " << t << ", valeur de x_t : " << x_t_plus_dt << ", valeur de y_t : " << y_t_plus_dt << std::endl ; 
		x_t_Cheyette_[date] = x_t ; //c'est aussi x_t_plus_dt ; c'est pour gérer cas t= 0 , x(0) = 0
		y_t_Cheyette_[date] = y_t ; //c'est aussi y_t_plus_dt ;
	}

	computeNumeraires() ;
}


//numeraires_[i] = B(indexOfSimulation[i] * tenor.yearFraction()   , fwdProbaT_)
void MC_Cheyette::computeNumeraires() const
{
	for (size_t i = 0; i < indexOfSimulation_.size() ; ++i)
	{
		double x_t = x_t_Cheyette_[i] ;
		double y_t = y_t_Cheyette_[i] ;
		numeraires_[i] = cheyetteDD_Model_->P(indexOfSimulation_[i] * tenorType_.YearFraction(), fwdProbaT_, x_t, y_t);
	}
}