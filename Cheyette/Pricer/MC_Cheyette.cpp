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
			double x_t_plus_dt =  x_t	+ cheyetteDD_Model_->drift_x_QT(t, fwdProbaT_, x_t, y_t) * dt 
							+ cheyetteDD_Model_->diffusion_x(t, x_t, y_t) * sqrt(dt) * gaussian_tmp[pasDiscretisation - 1] ;
						
			double y_t_plus_dt =  y_t	+ cheyetteDD_Model_->drift_y(t, x_t, y_t) * dt ;

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
		numeraires_[i] = cheyetteDD_Model_->P(i * tenorYearFrac, fwdProbaT_, x_t, y_t);
	}
}


//! one simulation - pour les produits vanille 
double MC_Cheyette::evaluateFloatLeg(	const size_t valuationIndex,
										const std::vector<size_t>& indexFloatLeg,
										const Tenor tenorFloatLeg) const
{
	//verifie que grille TenorStructure plus fine que les flux du swap
	double tenorFloatYearFrac	= tenorFloatLeg.YearFraction() ;
	double tenorStructYearFrac	= pTenorStructure_->get_tenorType().YearFraction() ;
	assert(tenorStructYearFrac <= tenorFloatYearFrac) ;

	double price = 0.0;

	//rapport tenorFloat / tenorStructure  (ex : swap 6M sur TenorStructure 3M) rapport de 2 = nbIndices
	size_t nbIndices = static_cast<size_t>(tenorFloatYearFrac / tenorStructYearFrac) ;

	for(size_t i=0; i<indexFloatLeg.size(); ++i)
	{
		size_t paymentIndex		= indexFloatLeg[i] ;

		if(paymentIndex <= valuationIndex)continue;

		double paymentDate		= paymentIndex * tenorStructYearFrac ;
		size_t fixingIndex		= paymentIndex - nbIndices ;
		double fixingDate		= fixingIndex * tenorStructYearFrac ;	

		//  !!! FAUX !!!
		//double libor = cheyetteDD_Model_->Libor(fixingDate, fixingDate, paymentDate, 
		//										x_t_Cheyette_[fixingIndex], y_t_Cheyette_[fixingIndex]) ; 

		double libor = cheyetteDD_Model_->Libor(valuationIndex * tenorStructYearFrac, fixingDate, paymentDate, 
												x_t_Cheyette_[valuationIndex], y_t_Cheyette_[valuationIndex]) ; 	

		//payoffFlow = nominal * deltaFloat * value  (ici NOMINAL = 1)
		double payoffFlow		= tenorFloatYearFrac * libor ;		//tenorFloat = deltaFloat																			
		double numeraireRatio	= numeraires_[valuationIndex]/numeraires_[paymentIndex];
		price += payoffFlow * numeraireRatio;
	}

//version 2 : mieux mais encore un petit ecart

// A MULTIPLIER PAR LES NUMERAIRES !!!
	//double startDateSwap	= (indexFloatLeg[0] - nbIndices ) * tenorStructYearFrac ;
	//double endDateSwap		= indexFloatLeg[indexFloatLeg.size() - 1] * tenorStructYearFrac ;
	//double price = 1 - cheyetteDD_Model_->P(startDateSwap, endDateSwap, x_t_Cheyette_[valuationIndex], 
	//																	y_t_Cheyette_[valuationIndex]) ;

	return price;
}


double MC_Cheyette::evaluateFixedLeg(	const size_t valuationIndex,
										const std::vector<size_t>& indexFixedLeg,
										const Tenor tenorFixedLeg, 
										const double fixedRate) const
{
	//verifie que grille TenorStructure plus fine que les flux du swap
	double tenorFixedYearFrac	= tenorFixedLeg.YearFraction() ;
	double tenorStructYearFrac	= pTenorStructure_->get_tenorType().YearFraction() ;
	assert(tenorStructYearFrac <= tenorFixedYearFrac) ;

	double price = 0.0;
	for(size_t i=0; i<indexFixedLeg.size(); ++i)
	{
		size_t paymentIndex		= indexFixedLeg[i] ;

		if(paymentIndex <= valuationIndex)continue;

		double paymentDate		= paymentIndex * tenorStructYearFrac ;
			
		//payoffFlow = nominal * deltaFloat * value  (ici NOMINAL = 1)
		double payoffFlow		= tenorFixedYearFrac * fixedRate ;
																					
		double numeraireRatio	= numeraires_[valuationIndex]  / numeraires_[paymentIndex];
		price += payoffFlow*numeraireRatio;
	}

	return price;
}


//! one simulation - pour les produits generiques derives
//! for one couponLeg and one simulation
//double MC_Cheyette::evaluateCouponLeg(	const size_t valuationIndex,
//										const CouponLeg_CONSTPTR couponLeg,
//										const Tenor tenorLeg)  const       
//						// !!! tenorLeg peut etre different de TenorStructure !!! 
//								// suppose all the coupons in couponLeg are the same type.  
//{
//	std::vector<Coupon_CONSTPTR> couponList=couponLeg->getLeg();
//
//	Coupon_CONSTPTR firstCoupon=couponList[0];
//	CappedFlooredCoupon_CONSTPTR firstCappedFlooredCoupon = boost::dynamic_pointer_cast<const CappedFlooredCoupon>(firstCoupon);
//	if(!firstCappedFlooredCoupon){
//		std::cout << "fail to cast cappedFlooredCoupon" << std::endl ;
//		throw("fail to cast cappedFlooredCoupon");
//	}
//	Rate_CONSTPTR firstRate = firstCappedFlooredCoupon->getRate();   
//
//	double price = 0.0;
//	if(boost::dynamic_pointer_cast<const LiborRate>(firstRate))
//	{  
//		LiborRate_CONSTPTR	firstLiborRate		=	boost::dynamic_pointer_cast<const LiborRate>(firstRate);
//
//		size_t posIndexValuationDate	= findIndex(indexValuationDate, indexOfSimulation_) ;
//		for(size_t i=0; i<couponList.size(); ++i)
//		{
//			Coupon_CONSTPTR coupon=couponList[i];		
//			size_t paymentIndex = coupon->getPaymentIndex();
//
//			if(paymentIndex<=indexValuationDate)continue;
//
//			CappedFlooredCoupon_CONSTPTR cappedFlooredCoupon = boost::dynamic_pointer_cast<const CappedFlooredCoupon>(coupon);
//
//			Rate_CONSTPTR rate = cappedFlooredCoupon->getRate();   
//			LiborRate_CONSTPTR	liborRate		=	boost::dynamic_pointer_cast<const LiborRate>(rate);
//
//			double tenorYearFrac	= tenorLeg.YearFraction() ;
//			size_t fixingIndex		= liborRate->getFixingTime();
//			double fixingDate		= fixingIndex * tenorYearFrac ;
//			double paymentDate		= paymentIndex * tenorYearFrac ;
//			
//			size_t pos = findIndex(fixingIndex, indexOfSimulation_) ;
//
//			double libor = cheyetteDD_Model_->Libor(fixingDate, fixingDate, paymentDate, x_t[pos], y_t[pos]) ; 
//			
//			//retourne nominal * delta * min(max(value, ...)...) 
//			double payoffFlow		= evaluateCappedFlooredCoupon(cappedFlooredCoupon, libor) ;
//																		
//			size_t posPaymentIndex	= findIndex(paymentIndex, indexOfSimulation_) ;
//			
//			double numeraireRatio	= numeraire[posIndexValuationDate]/numeraire[posPaymentIndex];
//
//			price += payoffFlow*numeraireRatio;
//		}
//	}
//	else if(boost::dynamic_pointer_cast<const ConstRate>(firstRate))
//	{
//		for(size_t i=0; i<couponList.size(); i++)
//		{
//			Coupon_CONSTPTR coupon=couponList[i];		
//			LMM::Index paymentIndex = coupon->getPaymentIndex();
//
//			if(paymentIndex<=indexValuationDate)continue;
//
//			CappedFlooredCoupon_CONSTPTR cappedFlooredCoupon = boost::dynamic_pointer_cast<const CappedFlooredCoupon>(coupon);
//
//			Rate_CONSTPTR rate				= cappedFlooredCoupon->getRate();   
//			ConstRate_CONSTPTR	constRate	= boost::dynamic_pointer_cast<const ConstRate>(rate);
//			double	constRateValue			= constRate->getConstRateValue();
//
//			size_t posIndexValuationDate	= findIndex(indexValuationDate, indexOfSimulation_) ;
//			size_t posPaymentIndex			= findIndex(paymentIndex, indexOfSimulation_) ;
//
//			price+=numeraire[posIndexValuationDate]/numeraire[posPaymentIndex]
//										* evaluateCappedFlooredCoupon(cappedFlooredCoupon, constRateValue);
//		}
//	}
//	else
//	{
//		throw("fail to cast LiborRate or ConstRate");
//	}
//
//	return price;
//}
//
//double MC_Cheyette::evaluateCappedFlooredCoupon(CappedFlooredCoupon_CONSTPTR cappedFlooredCoupon, 
//														   double rateValue) const
//{
//	double	nominal			=	cappedFlooredCoupon->getNominal();
//	double	delta			=	cappedFlooredCoupon->getDelta();
//	bool	ifFloored		=	cappedFlooredCoupon->getIfFloored();
//	double	floor			=	cappedFlooredCoupon->getFloorStrike();
//	bool	ifCapped		=	cappedFlooredCoupon->getIfCapped();
//	double	cap				=	cappedFlooredCoupon->getCapStrike();
//	double	multiFactor		=	cappedFlooredCoupon->getMultiFactor();
//	double	addFactor		=	cappedFlooredCoupon->getAddFactor();
//
//	double value		=	multiFactor*rateValue+addFactor;
//	if(ifCapped)
//		value	=	std::min(cap, value);
//	if(ifFloored)
//		value	=	std::max(floor, value);
//
//	return nominal* delta *value;
//}