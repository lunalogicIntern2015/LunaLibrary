#include "CheyetteMarketData_2.h"


void CheyetteMarketData_2::readFile(size_t fileNumber)
{
//lecture du fichier et recuperation des market data
	std::vector<std::string> mkt_file_list = InputFileManager::get_VCUB_FileList();
	std::string mkt_data_file = mkt_file_list[fileNumber] ;

	std::ifstream instream;
	std::string input_file_fullpath = LMMPATH::get_runtime_datapath() + mkt_data_file;
	instream.open(input_file_fullpath.c_str());

	std::string text_line;
	getline(instream,text_line);		//lecture 1ere ligne delta ; maturity ; ...
	//std::cout << text_line << std::endl ;
	getline(instream,text_line);
	//std::cout << text_line << std::endl ;

	//remplissage de l'indice 0 des vecteurs swaptions avec swaption par defaut
		LMMTenorStructure_PTR simulationStructure(new LMMTenorStructure(Tenor::_6M, 100) );
		VanillaSwap swap = VanillaSwap();
		VanillaSwaption_PTR pSwaption(new VanillaSwaption(swap, OptionType::OptionType::CALL)) ;
	vect_swaptions_convexity_pp_.push_back(pSwaption) ;vect_vol_convexity_pp_.push_back(-100) ;
	vect_swaptions_convexity_mm_.push_back(pSwaption) ;vect_vol_convexity_mm_.push_back(-100) ;

	while (text_line[0] != '#')
	{
		remplirCheyetteMarketData(text_line) ;
		getline(instream,text_line);
		//std::cout << text_line << std::endl ;
	}
}

std::string CheyetteMarketData_2::parcoursChaine(const std::string& chaine, size_t indexStart)
{
	size_t length = 0 ;
	while (chaine[indexStart + length] != ';' && indexStart + length < chaine.size()){++length ;}
	std::string st = chaine.substr(indexStart, length) ;
	//std::cout << st << std::endl ;
	return st ;
}

void CheyetteMarketData_2::remplirCheyetteMarketData(const std::string& chaine)
{
	//delta ; maturity; tenor; vol ; strike 
	size_t indexStart = 0 ;
	std::string sShift = parcoursChaine(chaine, indexStart) ;
	double shift = std::stod(sShift);

	if (shift > 0){shiftConvexity_ = shift ;}  //shift en %

	size_t sizeChaine = sShift.size() ;

	indexStart += sizeChaine + 1 ;
	std::string sMaturity = parcoursChaine(chaine, indexStart) ;
	int dMaturity = std::stoi(sMaturity);
	size_t maturity = static_cast<size_t>(dMaturity) ;
	sizeChaine = sMaturity.size() ;

	indexStart += sizeChaine + 1 ;
	std::string sTenor = parcoursChaine(chaine, indexStart ) ;
	int dTenor = std::stoi(sTenor);
	size_t tenor = static_cast<size_t>(dTenor) ;
	sizeChaine = sTenor.size() ;

	indexStart += sizeChaine + 1 ;
	std::string sVol = parcoursChaine(chaine, indexStart ) ;
	double vol = std::stod(sVol);
	sizeChaine = sVol.size() ;

	indexStart += sizeChaine + 1 ;
	std::string sStrike = parcoursChaine(chaine, indexStart ) ;
	double strike = std::stod(sStrike);
	sizeChaine = sStrike.size() ;

	LMMTenorStructure_PTR simulationStructure(
				new LMMTenorStructure(floatTenor_, maturity + tenor + 1) );
	size_t indexStartSwap = static_cast<size_t>(maturity / floatTenor_.YearFraction())  ;
	size_t indexEndSwap = static_cast<size_t>((maturity + tenor) / floatTenor_.YearFraction()) ;
	VanillaSwap swap = VanillaSwap(strike, indexStartSwap, indexEndSwap, floatTenor_, fixedTenor_, simulationStructure);
	VanillaSwaption_PTR pSwaption(new VanillaSwaption(swap, OptionType::OptionType::CALL)) ;
	
	if (shift > 0){vect_swaptions_convexity_pp_.push_back(pSwaption) ;vect_vol_convexity_pp_.push_back(vol) ;}
	if (shift < 0){vect_swaptions_convexity_mm_.push_back(pSwaption) ;vect_vol_convexity_mm_.push_back(vol) ;}


}

//void CheyetteMarketData_2::calculateSkew()
//{
//	size_t nbQuotes = vect_vol_skew_pp_.size() ;
//	assert(nbQuotes == vect_vol_skew_mm_.size()) ;
//
//	//quotes en %, shift en bp -> shift / 100
//	double shift = shiftSkew_ / 100. ;
//	for (size_t i = 0 ; i < nbQuotes ; ++i)
//	{
//		double skew = (vect_vol_skew_pp_[i] - vect_vol_skew_mm_[i]) / (2 * shift) ;
//		skew_.push_back(skew) ;
//	}
//	skew_.shrink_to_fit() ;
//}



void CheyetteMarketData_2::init(size_t fileNumber)
{
	readFile(fileNumber) ;
	//calculateSkew() ;
	//calculateConvexity() ;
}