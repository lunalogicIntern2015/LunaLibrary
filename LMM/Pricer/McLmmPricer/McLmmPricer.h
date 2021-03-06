#pragma once
#include <vector>

#include <LMM/Helper/Name.h>

#include <Instrument/Instrument.h>
#include <LMM/Pricer/Longstaff_Schwartz/McLmm_LS.h>

class McLmmPricer
{
public:
	virtual double price_on_oneSimulation(	Instrument_CONSTPTR instrument,
											LMM::Index evaluationDateIndex, 
											const matrix& liborMatrix, 
											const std::vector<double>& numeraire)const=0;

	//virtual double price_on_oneSimalation(	Instrument_CONSTPTR instrument,
	//										LMM::Index evaluationDateIndex, 
	//										McLmm_LS::LMMSimulationResult lmmSimulationResult)const
	//{
	//	return price_on_oneSimalation(	instrument,
	//									evaluationDateIndex, 
	//									lmmSimulationResult.get_liborMatrix(), 
	//									lmmSimulationResult.get_numeraire());
	//}

	virtual ~McLmmPricer(){}
};
typedef boost::shared_ptr<McLmmPricer> McLmmPricer_PTR;
typedef boost::shared_ptr<const McLmmPricer> McLmmPricer_CONSTPTR;

