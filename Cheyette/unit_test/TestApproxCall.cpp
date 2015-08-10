#include <Cheyette\unit_test\TestApproxCall.h>

//dS(t) = Phi(S_t) dW_t
//S(0) = S0
std::vector<double> approxMC_call(const size_t nbSimus, const double strike,
								  const double S0, const double T) //, const Boost_RR_Function_PTR& funcPhi)
{
	boost::function<double(double)> f1 = boost::bind(phi, _1);
//	Boost_RR_Function(const boost::function<double(double)>& func):func_(func){}
	Boost_RR_Function_PTR funcPhi (new Boost_RR_Function(f1)) ;

	std::vector<double> res(3) ;
	double somme_xi   = 0.0;
	double somme_xi_2 = 0.0;

	unsigned long seed = 47;
	RNGenerator_PTR  rnGenerator(new McGenerator(seed));

	for(size_t itrSimulation=0; itrSimulation < nbSimus ; ++itrSimulation)
	{
		if ((itrSimulation*10) % nbSimus == 0){std::cout << double(itrSimulation)/double(nbSimus)*100 << "%" << std::endl ;}	

		//simulation de S_t jusqu'à T
		size_t nbStepsPerYear = 250 ; 
		double dt = 1. / static_cast<double>(nbStepsPerYear) ;
		std::vector<double>   gaussian_tmp(static_cast<size_t>(T * nbStepsPerYear));  
		rnGenerator->generate(gaussian_tmp);		// generate Gaussian.

		double t = 0. ;
		double S_t = S0 ;
		for (size_t pasDiscretisation = 1 ; pasDiscretisation <= T * nbStepsPerYear ; ++pasDiscretisation)    
		{
			double t_plus_dt = t + dt ;
			double S_t_plus_dt =  S_t + funcPhi->operator()(S_t) * sqrt(dt) * gaussian_tmp[pasDiscretisation - 1] ;
						
			S_t = S_t_plus_dt ;
			t = t_plus_dt ;
		}
			
		double payoff = std::max(S_t, strike) ;  
		
		somme_xi	+= payoff ;	
		somme_xi_2	+= payoff * payoff ;
	}

	double mean_x	= somme_xi / nbSimus ; 
	double mean_x2	= somme_xi_2 / nbSimus ; 
	double variance = mean_x2 - mean_x * mean_x ;

	double IC_left	= mean_x - 2.57*std::sqrt(variance / nbSimus);
	double IC_right = mean_x + 2.57*std::sqrt(variance / nbSimus);

	res[0] = mean_x ;
	res[1] = IC_left ;
	res[2] = IC_right ;

	std::cout   << "prix MC swaption : " << mean_x << std::endl;
	std::cout	<< "nbSimulations    : " << nbSimus << std::endl;
	std::cout   << "99% confidence interval  [" << IC_left << " , " << IC_right	<< "]" << std::endl;

	return res ;

}


double phi(double St)
{
	return sqrt(St) ;
}