#include <Cheyette\Pricer\MC_Cheyette.h>


size_t MC_Cheyette::findIndex(size_t fixingIndex, std::vector<size_t> indexVector)
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

void MC_Cheyette::simulate_Euler() const
{
	//initialisation
	x_t_Cheyette_[0] = 0. ;	
	y_t_Cheyette_[0] = 0. ;

	double tenorYearFrac = tenorType_.YearFraction() ;	
	double x_t(0.), y_t(0.), t(0.) ;

	assert(indexOfSimulation_[0] == 0) ;

	for (size_t index = 1; index < indexOfSimulation_.size() ; ++index)    
	{ 
		std::vector<double>   gaussian_tmp(discretizationBetweenDates_[index]);  
		rnGenerator_->generate(gaussian_tmp);		// generate Gaussian.
		double dt ;
		if (index == 0)
			{dt = indexOfSimulation_[index] * tenorYearFrac / discretizationBetweenDates_[index] ;}
		else
			{dt = (indexOfSimulation_[index] - indexOfSimulation_[index-1]) * tenorYearFrac / 
																		discretizationBetweenDates_[index] ;}

		////DEBUG
		//for (size_t i = 0 ; i < gaussian_tmp.size() ; ++i)    
		//{			
		//	std::cout << "gaussian_tmp " << i << " : " << gaussian_tmp[i] << std::endl ;
		//}
		double t_plus_dt = t ;
		for (size_t pasDiscretisation = 1 ; pasDiscretisation <= discretizationBetweenDates_[index] ; ++pasDiscretisation)    
		{
			t_plus_dt += dt ;

			double x_t_plus_dt =  x_t	+ cheyetteDD_Model_->drift_x_QT(t, fwdProbaT_, x_t, y_t) * dt 
								+ cheyetteDD_Model_->diffusion_x(t, x_t) * sqrt(dt) * gaussian_tmp[pasDiscretisation - 1] ;
			double y_t_plus_dt =  y_t	+ cheyetteDD_Model_->drift_y(t, x_t, y_t) * dt ;
			x_t = x_t_plus_dt ;
			y_t = y_t_plus_dt ;
			t = t_plus_dt ;
		}
//		std::cout << "date t : " << t << ", valeur de x_t : " << x_t_plus_dt << ", valeur de y_t : " << y_t_plus_dt << std::endl ; 
		x_t_Cheyette_[index] = x_t ; //c'est aussi x_t_plus_dt 
		y_t_Cheyette_[index] = y_t ; //c'est aussi y_t_plus_dt
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

//! one simulation - pour les produits vanille derives
double MC_Cheyette::evaluateFloatLeg(	const size_t indexValuationDate,
										const std::vector<size_t>& indexFloatLeg,
										const std::vector<double>& numeraire, 
										const std::vector<double>& x_t, 
										const std::vector<double>& y_t,
										const Tenor tenor) const
//parametre Tenor represente le tenor du Libor de la floating leg
//SUPPOSE LE MEME QUE CELUI de la LMM structure
{
	double price = 0.0;
	size_t posIndexValuationDate	= findIndex(indexValuationDate, indexOfSimulation_) ;
	for(size_t i=0; i<indexFloatLeg.size(); ++i)
	{
		double tenorYearFrac	= tenor.YearFraction() ;
		size_t paymentIndex		= indexFloatLeg[i] ;
		double paymentDate		= paymentIndex * tenorYearFrac ;
		size_t fixingIndex		= paymentIndex - 1 ;
		double fixingDate		= fixingIndex * tenorYearFrac ;	//HYPOTHESE : tenor Libor = tenor LMM structure

		size_t pos = findIndex(fixingIndex, indexOfSimulation_) ;

		double libor = cheyetteDD_Model_->Libor(fixingDate, fixingDate, paymentDate, x_t[pos], y_t[pos]) ; 
			
		//payoffFlow = nominal * deltaFloat * value 
		//suppose NOMINAL = 1
		double payoffFlow		= tenorYearFrac * libor ;
																		
		size_t posPaymentIndex	= findIndex(paymentIndex, indexOfSimulation_) ;
			
		double numeraireRatio	= numeraire[posIndexValuationDate]/numeraire[posPaymentIndex];

		price += payoffFlow*numeraireRatio;
	}

	return price;
}


double MC_Cheyette::evaluateFixedLeg(	const size_t indexValuationDate,
										const std::vector<size_t>& indexFixedLeg,
										const std::vector<double>& numeraire, 
										const std::vector<double>& x_t, 
										const std::vector<double>& y_t,
										const Tenor tenorFixedLeg, 
										const double fixedRate) const
{

	double price = 0.0;
	size_t posIndexValuationDate	= findIndex(indexValuationDate, indexOfSimulation_) ;
	for(size_t i=0; i<indexFixedLeg.size(); ++i)
	{
		double tenorYearFrac	= tenorFixedLeg.YearFraction() ;
		size_t paymentIndex		= indexFixedLeg[i] ;
		double paymentDate		= paymentIndex * tenorYearFrac ;
			
		//payoffFlow = nominal * deltaFloat * value 
		//suppose NOMINAL = 1
		double payoffFlow		= tenorYearFrac * fixedRate ;
																		
		size_t posPaymentIndex	= findIndex(paymentIndex, indexOfSimulation_) ;
			
		double numeraireRatio	= numeraire[posIndexValuationDate]/numeraire[posPaymentIndex];

		price += payoffFlow*numeraireRatio;
	}

	return price;
}

//! one simulation - pour les produits generiques derives
//! for one couponLeg and one simulation
double MC_Cheyette::evaluateCouponLeg(	const size_t indexValuationDate,
										const CouponLeg_CONSTPTR couponLeg,
										const std::vector<double>& numeraire,
										const std::vector<double>& x_t, 
										const std::vector<double>& y_t, 
										const Tenor tenorLeg)  const       
						// !!! tenorLeg peut etre different de TenorStructure !!! 
								// suppose all the coupons in couponLeg are the same type.  
{
	std::vector<Coupon_CONSTPTR> couponList=couponLeg->getLeg();

	Coupon_CONSTPTR firstCoupon=couponList[0];
	CappedFlooredCoupon_CONSTPTR firstCappedFlooredCoupon = boost::dynamic_pointer_cast<const CappedFlooredCoupon>(firstCoupon);
	if(!firstCappedFlooredCoupon){
		std::cout << "fail to cast cappedFlooredCoupon" << std::endl ;
		throw("fail to cast cappedFlooredCoupon");
	}
	Rate_CONSTPTR firstRate = firstCappedFlooredCoupon->getRate();   

	double price = 0.0;
	if(boost::dynamic_pointer_cast<const LiborRate>(firstRate))
	{  
		LiborRate_CONSTPTR	firstLiborRate		=	boost::dynamic_pointer_cast<const LiborRate>(firstRate);

		size_t posIndexValuationDate	= findIndex(indexValuationDate, indexOfSimulation_) ;
		for(size_t i=0; i<couponList.size(); ++i)
		{
			Coupon_CONSTPTR coupon=couponList[i];		
			size_t paymentIndex = coupon->getPaymentIndex();

			if(paymentIndex<=indexValuationDate)continue;

			CappedFlooredCoupon_CONSTPTR cappedFlooredCoupon = boost::dynamic_pointer_cast<const CappedFlooredCoupon>(coupon);

			Rate_CONSTPTR rate = cappedFlooredCoupon->getRate();   
			LiborRate_CONSTPTR	liborRate		=	boost::dynamic_pointer_cast<const LiborRate>(rate);

			double tenorYearFrac	= tenorLeg.YearFraction() ;
			size_t fixingIndex		= liborRate->getFixingTime();
			double fixingDate		= fixingIndex * tenorYearFrac ;
			double paymentDate		= paymentIndex * tenorYearFrac ;
			
			size_t pos = findIndex(fixingIndex, indexOfSimulation_) ;

			double libor = cheyetteDD_Model_->Libor(fixingDate, fixingDate, paymentDate, x_t[pos], y_t[pos]) ; 
			
			//retourne nominal * delta * min(max(value, ...)...) 
			double payoffFlow		= evaluateCappedFlooredCoupon(cappedFlooredCoupon, libor) ;
																		
			size_t posPaymentIndex	= findIndex(paymentIndex, indexOfSimulation_) ;
			
			double numeraireRatio	= numeraire[posIndexValuationDate]/numeraire[posPaymentIndex];

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

			size_t posIndexValuationDate	= findIndex(indexValuationDate, indexOfSimulation_) ;
			size_t posPaymentIndex			= findIndex(paymentIndex, indexOfSimulation_) ;

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

double MC_Cheyette::evaluateCappedFlooredCoupon(CappedFlooredCoupon_CONSTPTR cappedFlooredCoupon, 
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