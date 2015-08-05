#pragma once

#include <stdlib.h>
#include <vector>
#include <string>
#include <boost/shared_ptr.hpp>
#include <Instrument/VanillaSwaption.h>


#include <LMM/Helper/Name.h>
#include <LMM/Helper/LMMTenorStructure.h>

#include <LMM/UpperTriangleVanillaSwaptionQuotes.h>

/*! \class LmmSwaptionMarketData
 *  This class parse swaption market data from a .csv file. The file has to respect the specified format
 *
 *  Data file are subsituted into blocs, each bloc end by "#####" (NOTE Exactly 5 repeats)
 *
 *  The block storing libor information has to start by "LIBOR_DF"
 *
 *  The blokc storing swaption information has to start by "SwaptionVol , ATM, bump"
 *  The information about swaption's moneyness is store after the tag ATM, a real indicate the bump vale
 *  If bump == 0, this is a ATM swaaption
 *  If not, this is a bumped swaption
 */


//swaptions ATM coterminales (Black vol)
class MarketData
{
private:

	std::string volType_ ; // ou faire class volType "Black" ou "Normal"
	std::vector<size_t> a_expiry_ ;
	std::vector<size_t> b_tenor_ ;
	std::vector<double> strikeATM_ ;
	std::vector<double> volQuotes_ ;	//vol ATM
	std::vector<double> skew_ ;			//skew ATM
	double shift_ ;					//skew = ( vol(strikeATM + shift_) - vol(strikeATM - shift_) ) / (2 * shift_)  
	std::vector<VanillaSwaption_PTR>  vect_swaptions_ ;

public:
	MarketData(	std::vector<size_t> a_expiry, std::vector<size_t> b_tenor, std::vector<double> strikeATM, 
				std::vector<double> volQuotes, std::vector<double> skew, double shift,
				std::vector<VanillaSwaption_PTR>  vect_swaptions,
				std::string volType = "Black")
		:	a_expiry_(a_expiry), b_tenor_(b_tenor), strikeATM_(strikeATM), 
			volQuotes_(volQuotes), skew_(skew), shift_(shift), vect_swaptions_(vect_swaptions), volType_(volType)
	{
		assert(a_expiry_.size() == b_tenor_.size()) ;
		//assert que les swaptions sont bien coterminales
	}

	std::string getVolType()			const {return volType_ ; } 
	std::vector<size_t> get_aExpiry()	const {return a_expiry_ ; }  
	std::vector<size_t> get_bTenor()	const {return b_tenor_ ; } 
	std::vector<double> getStrikeATM()	const {return strikeATM_ ; } 

	std::vector<double> getVolQuotes()	const {return volQuotes_ ; }	//vol ATM
	std::vector<double> getSkew()		const {return skew_ ; }			//skew ATM
	double				getShift()	const {return shift_ ;}
	std::vector<VanillaSwaption_PTR>  getVect_swaptions() const{return vect_swaptions_ ;} 

	//void parseFromMarketData(const std::string& filename) ;
};

typedef boost::shared_ptr<MarketData>       MarketData_PTR;
typedef boost::shared_ptr<const MarketData> MarketData_CONSTPTR;

