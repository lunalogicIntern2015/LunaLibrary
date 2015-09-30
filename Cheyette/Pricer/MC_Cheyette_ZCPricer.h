#pragma once

#include <boost/shared_ptr.hpp>
#include <Cheyette/Model/CheyetteDD_Model.h>
#include <Cheyette/Pricer/MC_Cheyette.h>	

//#include "GenericSwap.h"
#include <Instrument/VanillaSwap.h>  
#include <Numeric/Integrator1D.h>  //pour la fonction closestDate

class MC_Cheyette_ZCPricer
{
private:
	MC_Cheyette_PTR mcCheyette_; 

public:
	MC_Cheyette_ZCPricer(const MC_Cheyette_PTR& mcCheyette)
		: mcCheyette_(mcCheyette){}

	virtual ~MC_Cheyette_ZCPricer(){}

	//Price at time t, T maturity of the ZC bond
	double price(double t_valo, double T, size_t nbSimulation)  const ;

};


typedef boost::shared_ptr<MC_Cheyette_ZCPricer>			MC_Cheyette_ZCPricer_PTR;
typedef boost::shared_ptr<const MC_Cheyette_ZCPricer>	MC_Cheyette_ZCPricer_CONSTPTR;
