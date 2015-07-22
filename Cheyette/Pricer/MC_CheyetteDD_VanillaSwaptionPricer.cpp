#include <Cheyette/Pricer/MC_CheyetteDD_VanillaSwapPricer.h>
#include <Cheyette/Pricer/MC_CheyetteDD_VanillaSwaptionPricer.h>




//simulation
std::vector<double> MC_CheyetteDD_VanillaSwaptionPricer::price(VanillaSwaption_PTR vanillaSwaption, size_t nbSimulation) const
{
	std::vector<double> res(3) ;
	double somme_xi	= 0. ;
	double somme_xi2	= 0. ;
	size_t indexValuationDate = 0 ;
	VanillaSwap vanillaSwap = vanillaSwaption->getUnderlyingSwap() ;

	Tenor floatingLegTenor = vanillaSwap.get_floatingLegTenorType() ;
	Tenor fixedLegTenor = vanillaSwap.get_fixedLegTenorType() ;
	//MC
	for(size_t itrSimulation=0; itrSimulation<nbSimulation; ++itrSimulation)
	{
		simulate_Euler() ;
		//float
		double npv1  = evaluateFloatLeg(indexValuationDate, 
										vanillaSwap.get_floatingLegPaymentIndexSchedule(), 
										numeraires_, 
										x_t_Cheyette_,
										y_t_Cheyette_,
										floatingLegTenor);
		//fixed
		double npv2  = evaluateFixedLeg(indexValuationDate, 
										vanillaSwap.get_fixedLegPaymentIndexSchedule(), 
										numeraires_, 
										x_t_Cheyette_,
										y_t_Cheyette_,
										fixedLegTenor, 
										vanillaSwap.get_strike());
		double res = std::max(npv1 - npv2, 0.) ;  
		
		somme_xi += res ;
		somme_xi2 += res*res ;
	}
	double mean_x	= somme_xi / nbSimulation; 
	double mean_x2	= somme_xi2 / nbSimulation; 
 
	double variance = mean_x2 - mean_x * mean_x ;

	double IC_left	= mean_x - 2.57*std::sqrt(variance / nbSimulation);
	double IC_right = mean_x + 2.57*std::sqrt(variance / nbSimulation);

	res[0] = mean_x ;
	res[1] = IC_left ;
	res[2] = IC_right ;

	std::cout   << "prix MC swaption : " << mean_x << std::endl;
	std::cout	<< "nbSimulations    : " << nbSimulation << std::endl;
	std::cout   << "99% confidence interval  [" << IC_left << " , " << IC_right	<< "]" << std::endl;

	return res ;
}


void MC_CheyetteDD_VanillaSwaptionPricer::print(VanillaSwaption_PTR vanillaSwaption, 
												std::vector<size_t> nbSimus, 
												std::vector<double> prixMC,
												std::vector<double> IC_inf,
												std::vector<double> IC_sup) const
{
	assert(nbSimus.size() == prixMC.size() );
	assert(IC_inf.size() == IC_sup.size() ) ;
	assert(prixMC.size() == IC_inf.size() ) ; 
	time_t _time;
	struct tm timeInfo;
	char format[32];
 
	time(&_time);
	localtime_s(&timeInfo, &_time);
 
	strftime(format, 32, "%Y-%m-%d %H-%M", &timeInfo);
 
	std::cout << format << std::endl;
	ofstream o;
	std::stringstream fileName_s ;
	std::string directory = LMMPATH::get_runtime_datapath() ;
	fileName_s << directory << "TestMC_GenericSwaption_" << format << ".csv" ; 
	std::string fileName = fileName_s.str();

	o.open(fileName,  ios::out | ios::app );
	o	<<	endl;
	o	<<	endl;
	o	<<	endl;
	cheyetteDD_Model_->print(o) ;

	vanillaSwaption->getUnderlyingSwap().print(o) ;

	for (size_t i = 0 ; i < nbSimus.size() ; ++i)
	{
		o << "nb simulations : ; "	<< nbSimus[i] << " ; prix MC : ; " 
									<< prixMC[i] << " ; IC inf : ; " 
									<< IC_inf[i] << " ; IC sup : ; " 
									<< IC_sup[i] << endl ;
	}
	o.close();
}

void MC_CheyetteDD_VanillaSwaptionPricer::printMC_vs_approx(double approx, double b_barre, 
															double annuityA0, double swapRateS0, double volBlack, 
															double a, double b,
															VanillaSwaption_PTR vanillaSwaption, 
															std::vector<size_t> nbSimus, 
															std::vector<double> prixMC,
															std::vector<double> IC_inf,
															std::vector<double> IC_sup) const 
{
	assert(nbSimus.size() == prixMC.size() );
	assert(IC_inf.size() == IC_sup.size() ) ;
	assert(prixMC.size() == IC_inf.size() ) ; 
	time_t _time;
	struct tm timeInfo;
	char format[32];
 
	time(&_time);
	localtime_s(&timeInfo, &_time);
 
	strftime(format, 32, "%Y-%m-%d %H-%M", &timeInfo);
 
	std::cout << format << std::endl;
	ofstream o;
	std::stringstream fileName_s ;
	std::string directory = LMMPATH::get_output_path() ;
	fileName_s << directory << "Swaption_" << a << "Y" << b << "Y" << format << ".csv" ; 
	std::string fileName = fileName_s.str();

	o.open(fileName,  ios::out | ios::app );
	o	<<	endl;
	o	<<	endl;
	o	<<	endl;
	cheyetteDD_Model_->print(o) ;

	cheyetteDD_Model_->get_courbeInput_PTR()->print(o) ;

	vanillaSwaption->getUnderlyingSwap().print(o) ;

	o	<<	endl;
	o	<< "Prix approximation : ; " << approx << endl ;
	o	<< "b_barre : ; " << b_barre << endl ;
	o	<< "annuity A(0) : ; " << annuityA0 << endl ;
	o	<< "swap rate S(0) : ; " << swapRateS0 << endl ;
	o	<<	endl;

	o	<< "vol Black impli : ; " << volBlack << endl ;
	o	<<	endl;

	o << "prix approximation ; " <<"nb simulations ; " << " prix MC ; " << "IC inf ; " << "IC sup " << endl ;
	for (size_t i = 0 ; i < nbSimus.size() ; ++i)
	{
		o << approx << " ; " << nbSimus[i] << " ; " << prixMC[i] << " ; " << IC_inf[i] << " ; " << IC_sup[i] << endl ;
	}
	o.close();
}

void MC_CheyetteDD_VanillaSwaptionPricer::printMC_vs_approx(std::ofstream& o,
															double approx, double b_barre, 
															double annuityA0, double swapRateS0, double volBlack, 
															double a, double b,
															VanillaSwaption_PTR vanillaSwaption, 
															std::vector<size_t> nbSimus, 
															std::vector<double> prixMC,
															std::vector<double> IC_inf,
															std::vector<double> IC_sup) const 
{
	assert(nbSimus.size() == prixMC.size() );
	assert(IC_inf.size() == IC_sup.size() ) ;
	assert(prixMC.size() == IC_inf.size() ) ; 

	o << "Swaption_" << a << "Y" << b << "Y" << endl ; 
	o	<<	endl;

	o	<< "Prix approximation : ; " << approx << endl ;
	o	<< "b_barre : ; " << b_barre << endl ;
	o	<< "annuity A(0) : ; " << annuityA0 << endl ;
	o	<< "swap rate S(0) : ; " << swapRateS0 << endl ;
	o	<<	endl;

	o	<< "vol Black impli : ; " << volBlack << endl ;
	o	<<	endl;

	o << "prix approximation ; " <<"nb simulations ; " << " prix MC ; " << "IC inf ; " << "IC sup " << endl ;
	for (size_t i = 0 ; i < nbSimus.size() ; ++i)
	{
		o << approx << " ; " << nbSimus[i] << " ; " << prixMC[i] << " ; " << IC_inf[i] << " ; " << IC_sup[i] << endl ;
	}
}
