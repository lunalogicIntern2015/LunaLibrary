#include <Cheyette/Pricer/MC_CheyetteDD_VanillaSwapPricer.h>
#include <Cheyette/Pricer/MC_CheyetteDD_VanillaSwaptionPricer.h>




//simulation
std::vector<double> MC_CheyetteDD_VanillaSwaptionPricer::price(VanillaSwaption_PTR pVanillaSwaption, size_t nbSimulation) const
{
	VanillaSwap vanillaSwap = pVanillaSwaption->getUnderlyingSwap() ;
	Tenor floatingLegTenor = vanillaSwap.get_floatingLegTenorType() ;
	Tenor fixedLegTenor = vanillaSwap.get_fixedLegTenorType() ;

	assert(floatingLegTenor.YearFraction() >= pTenorStructure_->get_tenorType().YearFraction()) ;
	assert(fixedLegTenor.YearFraction() >= pTenorStructure_->get_tenorType().YearFraction()) ;
	assert(vanillaSwap.get_EndDate() <= fwdProbaT_) ;

	std::vector<double> res(3) ;
	double somme_xi   = 0.0;
	double somme_xi_2 = 0.0;

	size_t valuationIndexSwaption = 0 ;
	size_t valuationIndexSwap = pVanillaSwaption->getUnderlyingSwap().get_indexStart() ;
	
	for(size_t itrSimulation=0; itrSimulation<nbSimulation; ++itrSimulation)
	{
		if ((itrSimulation*10) % nbSimulation == 0){std::cout << double(itrSimulation)/double(nbSimulation)*100 << "%" << std::endl ;}	

		simulate_Euler() ;
		double npvFloatingLeg = evaluateFloatLeg(valuationIndexSwap, 
												vanillaSwap.get_floatingLegPaymentIndexSchedule(), floatingLegTenor);
		double npvFixedLeg = evaluateFixedLeg(valuationIndexSwap, vanillaSwap.get_fixedLegPaymentIndexSchedule(), 
										fixedLegTenor, vanillaSwap.get_strike());
		double payoffAtMaturity    = pVanillaSwaption->payoff(npvFloatingLeg,npvFixedLeg);			
		double value = payoffAtMaturity * numeraires_[valuationIndexSwaption]/numeraires_[valuationIndexSwap] ;

		somme_xi					  += value;	
		somme_xi_2 += value*value;
	}

	double mean_x	= somme_xi / nbSimulation; 
	double mean_x2	= somme_xi_2 / nbSimulation; 
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

//pricing de swaption pour des strikes correspondant � une standardized moneyness dans [-5, 5]
	//res[0] = prixSwaptionsPlusieursStrikes ;
	//res[1] = IC_left ; 
	//res[2] = IC_right ; 
	//res[3] = strikes ; 
	//res[4] = moneyness ;
std::vector<std::vector<double>> MC_CheyetteDD_VanillaSwaptionPricer::priceMultipleStrikes(	VanillaSwaption_PTR pVanillaSwaption, 
																							size_t nbSimulation, 
																							double S0, 
																							double sigma_ATM) const
{
	VanillaSwap vanillaSwap = pVanillaSwaption->getUnderlyingSwap() ;
	Tenor floatingLegTenor = vanillaSwap.get_floatingLegTenorType() ;
	Tenor fixedLegTenor = vanillaSwap.get_fixedLegTenorType() ;

	assert(floatingLegTenor.YearFraction() >= pTenorStructure_->get_tenorType().YearFraction()) ;
	assert(fixedLegTenor.YearFraction() >= pTenorStructure_->get_tenorType().YearFraction()) ;
	assert(vanillaSwap.get_EndDate() <= fwdProbaT_) ;

	//standardized moneyness
	size_t nbMoneyness = 11 ;		//moneyness = -5, -4, ... , 0, 1, ... 5
	std::vector<double> moneyness(nbMoneyness) ;			
	moneyness[0] = 5. ; moneyness[1] = 4. ; moneyness[2] = 3. ; moneyness[3] = 2. ; moneyness[4] = 1. ;
	moneyness[5] = 0. ;
	moneyness[6] = -1. ; moneyness[7] = -2. ; moneyness[8] = -3. ; moneyness[9] = -4. ; moneyness[10] = -5. ;
	
	//strike equivalent pour une standardized moneyness dans [-5 ; 5]
	std::vector<double> strikes(nbMoneyness) ;	
	double T0 = pVanillaSwaption->getUnderlyingSwap().get_StartDate() ;
	for (size_t i = 0 ; i < nbMoneyness ; ++i)
	{
		double strike = S0 / exp(sigma_ATM * sqrt(T0) * moneyness[i]) ; //   (strikeATM_Bloomberg + shiftStrike[i])/100. ;
		strikes[i] = strike ; 
	}

	//prix MC et intervalles de confiance pour les differents strikes
	std::vector<double> somme_xi(nbMoneyness, 0.) ;
	std::vector<double> somme_xi2(nbMoneyness, 0.) ;
	
	size_t valuationIndexSwaption = 0 ;
	size_t valuationIndexSwap = pVanillaSwaption->getUnderlyingSwap().get_indexStart() ;

	for(size_t itrSimulation=0; itrSimulation<nbSimulation; ++itrSimulation)
	{
		if ((itrSimulation*10) % nbSimulation == 0){std::cout << double(itrSimulation)/double(nbSimulation)*100 << "%" << std::endl ;}	

		simulate_Euler() ;
		double npvFloatingLeg = evaluateFloatLeg(valuationIndexSwap, 
												vanillaSwap.get_floatingLegPaymentIndexSchedule(), floatingLegTenor);

		for (size_t i = 0 ; i < nbMoneyness ; ++i)
		{
			double npvFixedLeg = evaluateFixedLeg(valuationIndexSwap, vanillaSwap.get_fixedLegPaymentIndexSchedule(), 
											fixedLegTenor, strikes[i] ) ; //vanillaSwap.get_strike());
			
			double payoffAtMaturity    = pVanillaSwaption->payoff(npvFloatingLeg,npvFixedLeg);			
			double value = payoffAtMaturity * numeraires_[valuationIndexSwaption]/numeraires_[valuationIndexSwap] ;

			somme_xi[i] += value ;
			somme_xi2[i] += value * value ;
		}
	}

	std::vector<double> prixSwaptionsPlusieursStrikes(nbMoneyness) ;
	std::vector<double> IC_left(nbMoneyness) ;
	std::vector<double> IC_right(nbMoneyness) ;
	std::vector<double> volBlack(nbMoneyness) ;
	for (size_t i = 0 ; i < nbMoneyness ; ++i)
	{
		double mean_x	= somme_xi[i] / nbSimulation; 
		double mean_x2	= somme_xi2[i] / nbSimulation; 
		double variance = mean_x2 - mean_x * mean_x ;

		prixSwaptionsPlusieursStrikes[i] = mean_x ;
		IC_left[i]	= mean_x - 2.57*std::sqrt(variance / nbSimulation);
		IC_right[i] = mean_x + 2.57*std::sqrt(variance / nbSimulation);
	}	

	std::vector<std::vector<double>> res(5) ;
	res[0] = prixSwaptionsPlusieursStrikes ;
	res[1] = IC_left ; 
	res[2] = IC_right ; 
	res[3] = strikes ; 
	res[4] = moneyness ; 

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
