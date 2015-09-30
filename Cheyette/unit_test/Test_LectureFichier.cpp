#include "Test_LectureFichier.h"

void testRecupData(size_t fileNumber)
{
	CheyetteMarketData_2 marketData ;
	marketData.init(fileNumber) ;
	
	int i = 0 ;
}


//void testReadFile(size_t fileNumber)
//{
////lecture du fichier et recuperation des market data
//	std::vector<std::string> mkt_file_list = InputFileManager::get_VCUB_FileList();
//	std::string mkt_data_file = mkt_file_list[fileNumber] ;
//
//	std::ifstream instream;
//	std::string input_file_fullpath = LMMPATH::get_runtime_datapath() + mkt_data_file;
//	instream.open(input_file_fullpath.c_str());
//
//	std::string text_line;
//	getline(instream,text_line);		//lecture 1ere ligne delta ; maturity ; ...
//	std::cout << text_line << std::endl ;
//	getline(instream,text_line);
//	std::cout << text_line << std::endl ;
//	while (text_line[0] != '#')
//	{
//		remplirCheyetteMarketData(text_line) ;
//		getline(instream,text_line);
//		std::cout << text_line << std::endl ;
//	}
//}
//
////string (char*) splitting
//	 // char str[] ="- This, a sample string.";
//  //char * pch;
//  //printf ("a ; bcd ; ze ; a",str);
//  //pch = strtok (str," ,.-");
//  //while (pch != NULL)
//  //{
//  //  printf ("%s\n",pch);
//  //  pch = strtok (NULL, ";");
//
//std::string parcoursChaine(const std::string& chaine, size_t indexStart)
//{
//	size_t length = 0 ;
//	while (chaine[indexStart + length] != ';' && indexStart + length < chaine.size()){++length ;}
//	std::string st = chaine.substr(indexStart, length) ;
//	std::cout << st << std::endl ;
//	return st ;
//}
//
//void remplirCheyetteMarketData(const std::string& chaine)
//{
//	//delta ; maturity; tenor; vol ; strike 
//	size_t indexStart = 0 ;
//	std::string sDelta = parcoursChaine(chaine, indexStart) ;
//	double delta = std::stod(sDelta);
//	size_t sizeChaine = sDelta.size() ;
//
//	indexStart += sizeChaine + 1 ;
//	std::string sMaturity = parcoursChaine(chaine, indexStart ) ;
//	double maturity = std::stod(sMaturity);
//	sizeChaine = sMaturity.size() ;
//
//	indexStart += sizeChaine + 1 ;
//	std::string sTenor = parcoursChaine(chaine, indexStart ) ;
//	double tenor = std::stod(sTenor);
//	sizeChaine = sTenor.size() ;
//
//	indexStart += sizeChaine + 1 ;
//	std::string sVol = parcoursChaine(chaine, indexStart ) ;
//	double vol = std::stod(sVol);
//	sizeChaine = sVol.size() ;
//
//	indexStart += sizeChaine + 1 ;
//	std::string sStrike = parcoursChaine(chaine, indexStart ) ;
//	double strike = std::stod(sStrike);
//	sizeChaine = sStrike.size() ;
//
//}
//
