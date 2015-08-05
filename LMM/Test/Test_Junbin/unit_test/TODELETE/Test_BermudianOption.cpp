//#include <JBLMM/Test/JBTests.h>
//
//#include <iostream>
//
//#include <JBLMM/Instrument/BermudanSwaption.h>
//
//void Test_BermudanOption()
//{
//	double strike = 0.02;
//	LMM::Index  indexStart = 12;    //6Y
//	LMM::Index  indexEnd   = 28;	//14Y
//	Tenor	floatingLegTenorType = Tenor::_6M;
//	Tenor	fixedLegTenorType    = Tenor::_1YR;
//	LMMTenorStructure_PTR simulationStructure(new LMMTenorStructure(Tenor::_6M , 15) );
//
//	VanillaSwap_CONSTPTR vanillaSwap_PTR(new VanillaSwap(strike, indexStart, indexEnd, floatingLegTenorType, fixedLegTenorType, simulationStructure));
//	
//	std::vector<LMM::Index> exerciceIndex;
//	exerciceIndex.push_back(14);
//	exerciceIndex.push_back(16);
//	exerciceIndex.push_back(18);
//	BermudanSwaption_PTR bermudanSwaption_PTR(new BermudanSwaption(vanillaSwap_PTR,exerciceIndex));
//
//	std::stringstream outputFileName_s; 
//	outputFileName_s<<"Test_BermudanOption"<<".csv";
//	std::string outputFileName = LMMPATH::get_Root_OutputPath() + outputFileName_s.str();
//	ofstream out;
//	out.open(outputFileName,  ios::out | ios::app );
//
//	bermudanSwaption_PTR->getVanillaSwap_PTR()->show();
//	bermudanSwaption_PTR->getVanillaSwap_PTR()->print(outputFileName);
//}