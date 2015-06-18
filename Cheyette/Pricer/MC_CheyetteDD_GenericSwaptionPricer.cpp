#include "MC_CheyetteDD_GenericSwaptionPricer.h"

//mettre en commun avec MC_CheyetteDD_GenericSwapPricer	
//avec namespace ou autre
static size_t findIndex2(size_t fixingIndex, std::vector<size_t> indexVector)
{
	size_t pos = 0 ;
	bool foundIt(false) ;
	for (size_t i = 0 ; i < indexVector.size() ; ++i)
	{
		if (indexVector[i] == fixingIndex){pos = i ;}
		if (indexVector[i] == fixingIndex){foundIt = true ;}
	}
	if (!foundIt){throw "pb dans FindIndex !!!!" ; }
	return pos ;
}

//simulation
std::vector<double> MC_CheyetteDD_GenericSwaptionPricer::price(GeneticSwaption_CONSTPTR geneticSwaption, size_t nbSimulation) const
{
	std::vector<double> res(3) ;
	double somme_xi	= 0. ;
	double somme_xi2	= 0. ;
	size_t indexValuationDate = 0 ;
	GeneticSwap_CONSTPTR genericSwap = geneticSwaption->getGeneticSwap() ;

	//MC
	for(size_t itrSimulation=0; itrSimulation<nbSimulation; ++itrSimulation)
	{
		mcCheyette_->simulate_Euler() ;
		double npv1  = evaluateCouponLeg(indexValuationDate, 
										genericSwap->getLeg1(), 
										mcCheyette_->getNumeraires(), 
										mcCheyette_->get_x_t_Cheyette(),
										mcCheyette_->get_y_t_Cheyette(),
										mcCheyette_->getTenorType());
		double npv2  = evaluateCouponLeg(indexValuationDate, 
										genericSwap->getLeg2(), 
										mcCheyette_->getNumeraires(), 
										mcCheyette_->get_x_t_Cheyette(),
										mcCheyette_->get_y_t_Cheyette(),
										mcCheyette_->getTenorType());
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

	std::cout   << "prix MC swap : " << mean_x << std::endl;
	std::cout	<< "nbSimulation : " << nbSimulation << std::endl;
	std::cout   << "99% confidence interval  [" << IC_left << " , " << IC_right	<< "]" << std::endl;

	return res ;
}

//! for one couponLeg and one simulation
double MC_CheyetteDD_GenericSwaptionPricer::evaluateCouponLeg(	const size_t indexValuationDate,
															const CouponLeg_CONSTPTR couponLeg,
															const std::vector<double>& numeraire,
															const std::vector<double>& x_t, 
															const std::vector<double>& y_t, 
															const Tenor tenor)  const 
								// suppose all the coupons in couponLeg are the same type.  
{
	std::vector<Coupon_CONSTPTR> couponList=couponLeg->getLeg();

	Coupon_CONSTPTR firstCoupon=couponList[0];
	CappedFlooredCoupon_CONSTPTR firstCappedFlooredCoupon = boost::dynamic_pointer_cast<const CappedFlooredCoupon>(firstCoupon);
	if(!firstCappedFlooredCoupon)
		throw("fail to cast cappedFlooredCoupon");

	Rate_CONSTPTR firstRate = firstCappedFlooredCoupon->getRate();   

	std::vector<size_t> vectIndex = mcCheyette_->getIndexOfSimulation() ;

	double price = 0.0;
	if(boost::dynamic_pointer_cast<const LiborRate>(firstRate))
	{  
		LiborRate_CONSTPTR	firstLiborRate		=	boost::dynamic_pointer_cast<const LiborRate>(firstRate);
		if(firstLiborRate->getDuration() != tenor)
		{
			throw "pas la meme structure de tenor !" ;
		}
		for(size_t i=0; i<couponList.size(); ++i)
		{
			Coupon_CONSTPTR coupon=couponList[i];		
			size_t paymentIndex = coupon->getPaymentIndex();

			if(paymentIndex<=indexValuationDate)continue;

			CappedFlooredCoupon_CONSTPTR cappedFlooredCoupon = boost::dynamic_pointer_cast<const CappedFlooredCoupon>(coupon);

			Rate_CONSTPTR rate = cappedFlooredCoupon->getRate();   
			LiborRate_CONSTPTR	liborRate		=	boost::dynamic_pointer_cast<const LiborRate>(rate);

			double tenorYearFrac	= tenor.YearFraction() ;
			size_t fixingIndex		= liborRate->getFixingTime();
			double fixingDate		= fixingIndex * tenorYearFrac ;
			double paymentDate		= paymentIndex * tenorYearFrac ;
			
			size_t pos = findIndex2(fixingIndex, mcCheyette_->getIndexOfSimulation()) ;

			double libor = mcCheyette_->getCheyetteDD_Model()->Libor(fixingDate, 
																		fixingDate, paymentDate, x_t[pos], y_t[pos]) ; 
			
			//retourne nominal * delta * min(max(value, ...)...) 
			double payoffFlow = evaluateCappedFlooredCoupon(cappedFlooredCoupon, libor) ;
																		
			size_t posIndexValuationDate	= findIndex2(indexValuationDate, vectIndex) ;
			size_t posPaymentIndex			= findIndex2(paymentIndex, vectIndex) ;
			
			double numeraireRatio = numeraire[posIndexValuationDate]/numeraire[posPaymentIndex];

			price += payoffFlow*numeraireRatio;
		}
	}
	else if(boost::dynamic_pointer_cast<const ConstRate>(firstRate))
	{
		for(size_t i=0; i<couponList.size(); i++)
		{
			Coupon_CONSTPTR coupon=couponList[i];		
			LMM::Index paymentIndex = coupon->getPaymentIndex();

			if(paymentIndex<=indexValuationDate)continue;

			CappedFlooredCoupon_CONSTPTR cappedFlooredCoupon = boost::dynamic_pointer_cast<const CappedFlooredCoupon>(coupon);

			Rate_CONSTPTR rate				= cappedFlooredCoupon->getRate();   
			ConstRate_CONSTPTR	constRate	= boost::dynamic_pointer_cast<const ConstRate>(rate);
			double	constRateValue			= constRate->getConstRateValue();

			size_t posIndexValuationDate	= findIndex2(indexValuationDate, vectIndex) ;
			size_t posPaymentIndex			= findIndex2(paymentIndex, vectIndex) ;

			price+=numeraire[posIndexValuationDate]/numeraire[posPaymentIndex]
										* evaluateCappedFlooredCoupon(cappedFlooredCoupon, constRateValue);
		}
	}
	else
	{
		throw("fail to cast LiborRate or ConstRate");
	}

	return price;
}

double MC_CheyetteDD_GenericSwaptionPricer::evaluateCappedFlooredCoupon(CappedFlooredCoupon_CONSTPTR cappedFlooredCoupon, 
														   double rateValue) const
{
	double	nominal			=	cappedFlooredCoupon->getNominal();
	double	delta			=	cappedFlooredCoupon->getDelta();
	bool	ifFloored		=	cappedFlooredCoupon->getIfFloored();
	double	floor			=	cappedFlooredCoupon->getFloorStrike();
	bool	ifCapped		=	cappedFlooredCoupon->getIfCapped();
	double	cap				=	cappedFlooredCoupon->getCapStrike();
	double	multiFactor		=	cappedFlooredCoupon->getMultiFactor();
	double	addFactor		=	cappedFlooredCoupon->getAddFactor();

	double value		=	multiFactor*rateValue+addFactor;
	if(ifCapped)
		value	=	std::min(cap, value);
	if(ifFloored)
		value	=	std::max(floor, value);

	return nominal* delta *value;
}

void MC_CheyetteDD_GenericSwaptionPricer::print(GeneticSwaption_CONSTPTR genericSwaption, 
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
	mcCheyette_->getCheyetteDD_Model()->print(o) ;

	genericSwaption->getGeneticSwap()->print(o) ;

	for (size_t i = 0 ; i < nbSimus.size() ; ++i)
	{
		o << "nb simulations : ; "	<< nbSimus[i] << " ; prix MC : ; " 
									<< prixMC[i] << " ; IC inf : ; " 
									<< IC_inf[i] << " ; IC sup : ; " 
									<< IC_sup[i] << endl ;
	}
	o.close();
}

void MC_CheyetteDD_GenericSwaptionPricer::printMC_vs_approx(double approx, 
															GeneticSwaption_CONSTPTR genericSwaption, 
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
	mcCheyette_->getCheyetteDD_Model()->print(o) ;

	mcCheyette_->getCheyetteDD_Model()->get_courbeInput_PTR()->print(o) ;

	genericSwaption->getGeneticSwap()->print(o) ;

	o	<<	endl;
	o << "Prix approximation : ; " << approx << endl ;
	o	<<	endl;

	o << "prix approximation ; " <<"nb simulations ; " << " prix MC ; " << "IC inf ; " << "IC sup " << endl ;
	for (size_t i = 0 ; i < nbSimus.size() ; ++i)
	{
		o << approx << " ; " << nbSimus[i] << " ; " << prixMC[i] << " ; " << IC_inf[i] << " ; " << IC_sup[i] << endl ;
	}
	o.close();
}