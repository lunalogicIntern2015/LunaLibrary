#pragma once

#include <Cheyette/Calibrator/CoTerminalSwaptionQuotes.h>


//ATM
class CoTerminalSwaptionVol : public CoTerminalSwaptionQuotes 
{
private:
	double diagonalSwaptionVol_;		// vol implicites stockées par expiry croissantes
													//ex : 1Y 3Y, 2Y 2Y, 3Y 1Y

public:	
	CoTerminalSwaptionVol::CoTerminalSwaptionVol(	double diagonalSwaptionVol,
													size_t vectorExpiry,
													size_t vectorTenor,
													double strike,
													VanillaSwaption_PTR swaption) 
		: CoTerminalSwaptionQuotes(vectorExpiry, vectorTenor, strike, swaption), diagonalSwaptionVol_(diagonalSwaptionVol)
	{}

//getters
	double getDiagonalSwaptionVol() const {return diagonalSwaptionVol_ ;}

	virtual double get_MinQuote() const {return 0 ;}
	virtual double get_MaxQuote() const {return 0 ;}

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

typedef boost::shared_ptr<CoTerminalSwaptionVol> CoTerminalSwaptionVol_PTR;
typedef boost::shared_ptr<const CoTerminalSwaptionVol> CoTerminalSwaptionVol_CONSTPTR;

