#pragma once

#include <vector>
#include <LMM/instrument/VanillaSwaption.h>


//ATM
class CoTerminalSwaptionQuotes 
{
protected:
	std::vector<size_t> vectorExpiry_ ;				//1Y, 2Y, 3Y
	std::vector<size_t> vectorTenor_ ;				//1Y, 2Y, 3Y
	std::vector<double> strike_ ;				
	
	//std::string data_file_name_ ;

	//LMMTenorStructure_PTR lmmTenorStructure_;
	//size_t lastYear_;         // upperTriangleSwaption(i year, j year), i,j \in [1, lastYear_]
	//Tenor  fixedTenor_;
	//Tenor  floatingTenor_; 
	//size_t indexRatio_;       // fixedTenor/floatingTenor

public:	
	CoTerminalSwaptionQuotes::CoTerminalSwaptionQuotes(	std::vector<size_t> vectorExpiry,
														std::vector<size_t> vectorTenor,
														std::vector<double> strike) 
		: vectorExpiry_(vectorExpiry), vectorTenor_(vectorTenor), strike_(strike) 
	{
		//assert que les swaptions de calibration sont bien coterminales
	}

//getters
	std::vector<size_t> getVectorExpiry() const {return vectorExpiry_ ;}
	std::vector<size_t> getVectorTenor() const {return vectorTenor_ ;}
	std::vector<double> getStrike() const {return strike_ ;}

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

