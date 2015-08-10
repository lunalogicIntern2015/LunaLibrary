#pragma once

#include <vector>
#include <Instrument/VanillaSwaption.h>


//ATM
class CoTerminalSwaptionQuotes 
{
protected:
	size_t vectorExpiry_ ;				//1Y, 2Y, 3Y
	size_t vectorTenor_ ;				//1Y, 2Y, 3Y
	double strike_ ;				
	VanillaSwaption_PTR swaption_ ; 
	//std::string data_file_name_ ;

	//LMMTenorStructure_PTR lmmTenorStructure_;
	//size_t lastYear_;         // upperTriangleSwaption(i year, j year), i,j \in [1, lastYear_]
	//Tenor  fixedTenor_;
	//Tenor  floatingTenor_; 
	//size_t indexRatio_;       // fixedTenor/floatingTenor

public:	
	CoTerminalSwaptionQuotes::CoTerminalSwaptionQuotes(	size_t vectorExpiry,
														size_t vectorTenor,
														double strike,
														VanillaSwaption_PTR swaption) 
		: vectorExpiry_(vectorExpiry), vectorTenor_(vectorTenor), strike_(strike), swaption_(swaption) 
	{}

//getters
	size_t getVectorExpiry() const {return vectorExpiry_ ;}
	size_t getVectorTenor() const {return vectorTenor_ ;}
	double getStrike() const {return strike_ ;}
	VanillaSwaption_PTR getSwaption() const {return swaption_ ;}

	virtual double get_MinQuote() const = 0 ;
	virtual double get_MaxQuote() const = 0 ;

	//const std::string get_Data_FileName()const { return data_file_name_;}
	//void set_Data_FileName(const std::string& mkt_data_filename) { data_file_name_=mkt_data_filename;}

//calcul du smile de marché
	//static boost::shared_ptr<const std::vector<double>> create_ATMSwaptionImpliedVol
	//	(
	//		LiborQuotes_ConstPTR libor_quotes_ptr,
	//		const Tenor&  fixedTenor,
	//		const Tenor&  floatingTenor,
	//		LmmVanillaSwaptionApproxPricer_Rebonato_PTR black_vol_approx_ptr,
	//		const double& strike_bump = 0
	//	);


	void print(const std::string& filename,  const bool erase_file=true) const;

};

typedef boost::shared_ptr<CoTerminalSwaptionQuotes> CoTerminalSwaptionQuotes_PTR;
typedef boost::shared_ptr<const CoTerminalSwaptionQuotes> CoTerminalSwaptionQuotes_CONSTPTR;