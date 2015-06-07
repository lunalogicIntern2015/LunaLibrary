#pragma once
#include <boost/shared_ptr.hpp>
#include <Cheyette/CheyetteModel/CheyetteDD_Model.h>
#include <Cheyette/Pricer/MC_Cheyette.h>	

//#include "GeneticSwap.h"
#include <LMM/instrument/VanillaSwap.h>  
#include <LMM/numeric/Integrator1D.h>  //pour la fonction closestDate

class MC_CheyetteDD_VanillaSwapPricer
{
private:
	MC_Cheyette_PTR mcCheyette_; 

public:
	MC_CheyetteDD_VanillaSwapPricer(const MC_Cheyette_PTR& mcCheyette)
		: mcCheyette_(mcCheyette){}

	virtual ~MC_CheyetteDD_VanillaSwapPricer(){}

	//Price at time T0=0
	//---------     payer swap     ---------
	double swapNPV(double t_valo, const VanillaSwap& vanillaSwap, size_t nbSimulation)  const ;

	//double swapRate(LMM::Index indexValuationDate,
	//				const VanillaSwap& vanillaSwap,
	//				const std::vector<double>& numeraire, 
	//				const std::vector<double>& xtCheyette, 
	//				const std::vector<double>& yt_Cheyette) const;

};


typedef boost::shared_ptr<MC_CheyetteDD_VanillaSwapPricer> MC_CheyetteDD_VanillaSwapPricer_PTR;
typedef boost::shared_ptr<const MC_CheyetteDD_VanillaSwapPricer> MC_CheyetteDD_VanillaSwapPricer_CONSTPTR;
