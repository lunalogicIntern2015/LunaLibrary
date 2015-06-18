#pragma once
#include <boost/shared_ptr.hpp>
#include <Cheyette/CheyetteModel/CheyetteDD_Model.h>
#include <Cheyette/Pricer/MC_Cheyette.h>	

#include <JBLMM/Instrument/GeneticSwaption.h>
#include <JBLMM/Element/Coupon.h>
#include <JBLMM/Element/CappedFlooredCoupon.h>
#include <JBLMM/Element/LiborRate.h>


class MC_CheyetteDD_GenericSwaptionPricer
{
private:
	MC_Cheyette_PTR mcCheyette_; 

public:
	MC_CheyetteDD_GenericSwaptionPricer(const MC_Cheyette_PTR& mcCheyette)
		: mcCheyette_(mcCheyette){}

	virtual ~MC_CheyetteDD_GenericSwaptionPricer(){}

	//gettor
	MC_Cheyette_PTR getMcCheyette()const{return mcCheyette_;}

	//! Pricing at time T0=0
	std::vector<double> price(GeneticSwaption_CONSTPTR genericSwaption, size_t nbSimulation) const;

	void print(GeneticSwaption_CONSTPTR genericSwaption, 
				std::vector<size_t> nbSimus, 
				std::vector<double> prixMC,
				std::vector<double> IC_inf,
				std::vector<double> IC_sup) const ;

	void printMC_vs_approx(double approx, 
							GeneticSwaption_CONSTPTR genericSwaption, 
							std::vector<size_t> nbSimus, 
							std::vector<double> prixMC,
							std::vector<double> IC_inf,
							std::vector<double> IC_sup) const ;

protected: 

	//! one simulation
	virtual double evaluateCouponLeg(	const LMM::Index indexValuationDate,
										const CouponLeg_CONSTPTR couponLeg,
										const std::vector<double>& numeraire, 
										const std::vector<double>& x_t, 
										const std::vector<double>& y_t,
										const Tenor tenor) const;

	double evaluateCappedFlooredCoupon( CappedFlooredCoupon_CONSTPTR cappedFlooredCoupon, 
										double rateValue) const ;

};

typedef boost::shared_ptr<MC_CheyetteDD_GenericSwaptionPricer>       MC_CheyetteDD_GenericSwaptionPricer_PTR;
typedef boost::shared_ptr<const MC_CheyetteDD_GenericSwaptionPricer> MC_CheyetteDD_GenericSwaptionPricer_CONSTPTR;