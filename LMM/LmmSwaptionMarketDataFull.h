#pragma once

#include <LMM/LmmSwaptionMarketData.h>

//reads VCUB files with very ITM/OTM vol quotes
//inherits from LmmSwaptionMarketData which only reads files with vol for ATM, ATM+/- 5bp

class LmmSwaptionMarketDataFull : public LmmSwaptionMarketData
{
private:
	//SwaptionMarketData swaptionMarketData_ATMpp_; // atm ++
	//SwaptionMarketData swaptionMarketData_ATMmm_; // atm --
	
	SwaptionMarketData swaptionMarketData_ATMp50_; // atm + 50bp
	SwaptionMarketData swaptionMarketData_ATMp100_; // atm + 100bp
	SwaptionMarketData swaptionMarketData_ATMp200_; // atm + 200bp

	SwaptionMarketData swaptionMarketData_ATMm50_; // atm - 50bp
	SwaptionMarketData swaptionMarketData_ATMm100_; // atm - 100bp
	SwaptionMarketData swaptionMarketData_ATMm200_; // atm - 200bp

public: 
	//ctor
	LmmSwaptionMarketDataFull(const Tenor& fixedtenor, const Tenor& floattenor, const size_t maxNbYear) ;

	UpperTriangleVanillaSwaptionQuotes_PTR get_SwaptionQuotes_ATMp50()const 
	{ return swaptionMarketData_ATMp50_.first;  }
	UpperTriangleVanillaSwaptionQuotes_PTR get_SwaptionQuotes_ATMp100()const 
	{ return swaptionMarketData_ATMp100_.first; }
	UpperTriangleVanillaSwaptionQuotes_PTR get_SwaptionQuotes_ATMp200()const 
	{ return swaptionMarketData_ATMp200_.first; }

	UpperTriangleVanillaSwaptionQuotes_PTR get_SwaptionQuotes_ATMm50()const 
	{ return swaptionMarketData_ATMm50_.first;  }
	UpperTriangleVanillaSwaptionQuotes_PTR get_SwaptionQuotes_ATMm100()const 
	{ return swaptionMarketData_ATMm100_.first; }
	UpperTriangleVanillaSwaptionQuotes_PTR get_SwaptionQuotes_ATMm200()const 
	{ return swaptionMarketData_ATMm200_.first; }

	virtual void parseFromMarketData(const std::string& filename);
	virtual void buildSwationQuotes( const double & strike_bump
								   , const std::vector<double>& mkt_expirities
								   , const std::vector<double>& mkt_tenor
								   , const std::vector<std::vector<double> >& mkt_quote
								   , const std::vector<std::vector<double> >& mkt_strike) ;
} ;

typedef boost::shared_ptr<LmmSwaptionMarketDataFull> LmmSwaptionMarketDataFull_PTR ;
typedef boost::shared_ptr<const LmmSwaptionMarketDataFull> LmmSwaptionMarketDataFull_CONSTPTR ;