#pragma once

#include <Cheyette/Calibrator/CoTerminalSwaptionQuotes.h>


class CoTerminalSwaptionSkew : public CoTerminalSwaptionQuotes
{
private:
	std::vector<double> diagonalSwaptionSkew_;	// skew stocké par expiry croissantes
													//ex : 1Y 3Y, 2Y 2Y, 3Y 1Y			
	double shift_ ;								//shift appliqué au strike pour calculer le skew

	//std::string data_file_name_ ;

public:	
	CoTerminalSwaptionSkew::CoTerminalSwaptionSkew(	std::vector<double> diagonalSwaptionSkew,
														std::vector<size_t> vectorExpiry,
														std::vector<double> strike, 
														double shift) 
		: CoTerminalSwaptionQuotes(vectorExpiry, strike), diagonalSwaptionSkew_(diagonalSwaptionSkew), shift_(shift) 
	{}

//getters
	std::vector<double> getDiagonalSwaptionSkew() const {return diagonalSwaptionSkew_ ;}
	double				getShift() const {return shift_ ;}

	virtual double get_MinQuote() const {return 0 ;}
	virtual double get_MaxQuote() const {return 0 ;}

	//const std::string get_Data_FileName()const { return data_file_name_;}
	//void set_Data_FileName(const std::string& mkt_data_filename) { data_file_name_=mkt_data_filename;}

	void print(const std::string& filename,  const bool erase_file=true) const;

};

typedef boost::shared_ptr<CoTerminalSwaptionSkew> CoTerminalSwaptionSkew_PTR;
typedef boost::shared_ptr<const CoTerminalSwaptionSkew> CoTerminalSwaptionSkew_CONSTPTR;

