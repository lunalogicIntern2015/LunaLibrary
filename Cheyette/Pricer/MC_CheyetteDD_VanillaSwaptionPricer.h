#pragma once
#include <boost/shared_ptr.hpp>
#include <Cheyette/CheyetteModel/CheyetteDD_Model.h>
#include <Cheyette/Pricer/MC_Cheyette.h>	

#include <LMM/instrument/VanillaSwap.h>  
#include <LMM/instrument/VanillaSwaption.h>
#include <LMM/numeric/Integrator1D.h>  //pour la fonction closestDate

class MC_CheyetteDD_VanillaSwaptionPricer
{
private:
	MC_Cheyette_PTR mcCheyette_; 

public:
	MC_CheyetteDD_VanillaSwaptionPricer(const MC_Cheyette_PTR& mcCheyette)
		: mcCheyette_(mcCheyette){}

	virtual ~MC_CheyetteDD_VanillaSwaptionPricer(){}

	//Price at time T0=0
	//---------     payer swaption     ---------
	double price(double t_valo, const VanillaSwaption& vanillaSwaption, size_t nbSimulation)  const ;

	
};


typedef boost::shared_ptr<MC_CheyetteDD_VanillaSwaptionPricer> MC_CheyetteDD_VanillaSwaptionPricer_PTR;
typedef boost::shared_ptr<const MC_CheyetteDD_VanillaSwaptionPricer> MC_CheyetteDD_VanillaSwaptionPricer_CONSTPTR;
