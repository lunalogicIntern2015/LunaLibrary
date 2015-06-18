#include "MC_CheyetteDD_GenericSwapPricer.h"

//mettre en commun avec MC_CheyetteDD_GenericSwapPricer	
//avec namespace ou autre
static size_t findIndex(size_t fixingIndex, std::vector<size_t> indexVector)
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
double MC_CheyetteDD_GenericSwapPricer::swapNPV(GeneticSwap_CONSTPTR geneticSwap, size_t nbSimulation) const
{
//	std::vector<double> simu(nbSimulation) ;
	double result	= 0. ;
	size_t indexValuationDate = 0 ;

	//MC
	for(size_t itrSimulation=0; itrSimulation<nbSimulation; ++itrSimulation)
	{
		mcCheyette_->simulate_Euler() ;
		double npv1  = evaluateCouponLeg(indexValuationDate, 
										geneticSwap->getLeg1(), 
										mcCheyette_->getNumeraires(), 
										mcCheyette_->get_x_t_Cheyette(),
										mcCheyette_->get_y_t_Cheyette(),
										mcCheyette_->getTenorType());
		double npv2  = evaluateCouponLeg(indexValuationDate, 
										geneticSwap->getLeg2(), 
										mcCheyette_->getNumeraires(), 
										mcCheyette_->get_x_t_Cheyette(),
										mcCheyette_->get_y_t_Cheyette(),
										mcCheyette_->getTenorType());
		result += npv1 - npv2;
//		simu[itrSimulation] = result ;
	}
	result   /=nbSimulation; 

	//double variance = 0. ;
	//for (size_t i = 0 ; i < nbSimulation ; ++i)
	//{
	//	variance += pow(simu[i] - result, 2) ;
	//}

	//variance /= nbSimulation ;
	std::cout   << "prix MC swap : " << result << std::endl;
	//std::cout	<< "nbSimulation : " << nbSimulation << std::endl;
	//std::cout   << "99% confidence interval  [" << result - 2.57*std::sqrt(variance / nbSimulation) 
	//											<< " , " 
	//											<< result + 2.57*std::sqrt(variance / nbSimulation) 
	//											<< "]"
	//											<< std::endl;


	return result;
}

//! for one couponLeg and one simulation
double MC_CheyetteDD_GenericSwapPricer::evaluateCouponLeg(	const size_t indexValuationDate,
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
			
			size_t pos = findIndex(fixingIndex, mcCheyette_->getIndexOfSimulation()) ;

			double libor = mcCheyette_->getCheyetteDD_Model()->Libor(fixingDate, 
																		fixingDate, paymentDate, x_t[pos], y_t[pos]) ; 
			
			//retourne nominal * delta * min(max(value, ...)...) 
			double payoffFlow = evaluateCappedFlooredCoupon(cappedFlooredCoupon, libor) ;
																		
			size_t posIndexValuationDate	= findIndex(indexValuationDate, vectIndex) ;
			size_t posPaymentIndex			= findIndex(paymentIndex, vectIndex) ;
			
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

			size_t posIndexValuationDate	= findIndex(indexValuationDate, vectIndex) ;
			size_t posPaymentIndex			= findIndex(paymentIndex, vectIndex) ;

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

double MC_CheyetteDD_GenericSwapPricer::evaluateCappedFlooredCoupon(CappedFlooredCoupon_CONSTPTR cappedFlooredCoupon, 
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
