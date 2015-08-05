//#include <Cheyette/Calibration/MarketData.h>
//
//void MarketData::parseFromMarketData(const std::string& filename)
//{
//
//	std::ifstream instream;
//	std::string input_file_fullpath = LMMPATH::get_runtime_datapath() + filename;
//	instream.open(input_file_fullpath.c_str());
//
//	getGeneralInfo(instream);
//
//	std::vector<double> df_dates_buffer;
//	std::vector<double> df_values_buffer;
//	getDF_MarketData(instream, df_dates_buffer, df_values_buffer);
//	interpolateAndBuildLiborQuotes(df_dates_buffer, df_values_buffer);
//
//	// If a curve of forward rate quote is specified, use this curve instead of interpolation from Discount Factor curve
//	std::string text_line;
//	getline(instream,text_line);
//	if(!isLineEndOfBlock(text_line))
//	{
//		std::vector<double> libor_quote_values;
//		parseLIBORS_Values(text_line,libor_quote_values);
//		//pLiborQuotes_->reset_Libor(libor_quote_values); just for testing if the data taken from the bloomberg forward curve analysis is better, but it is not. 
//	}
//
//
//	for(size_t itSWPMTab=0;itSWPMTab<3;++itSWPMTab)
//	{
//		double strike_bump = -100000000;
//		std::vector<double> atm_swpm_experities_buffer;
//		std::vector<double> atm_swpm_tenor_buffer;
//		std::vector< std::vector<double> > atm_swpm_quote_buffer;
//		std::vector< std::vector<double> > atm_swpm_strike_buffer;
//		getSwaptionMarketData(instream,strike_bump, atm_swpm_experities_buffer, atm_swpm_tenor_buffer, atm_swpm_quote_buffer, atm_swpm_strike_buffer);
//
//		if(itSWPMTab==0)// first block is ATM Swaption VCUB
//		{
//			getSwaptionQuoteFirstColumn(atm_swpm_experities_buffer, atm_swpm_quote_buffer);
//		}
//
//		buildSwationQuotes( strike_bump, atm_swpm_experities_buffer, atm_swpm_tenor_buffer, atm_swpm_quote_buffer, atm_swpm_strike_buffer );
//	}
//
//	swaptionMarketData_ATM_.first->set_Data_FileName(filename);
//
//	buildSwationQuotesVol_skew();
//
//	swaptionMarketData_ATMpp_.first->set_Data_FileName(filename);
//	swaptionMarketData_ATMmm_.first->set_Data_FileName(filename);
//
//	instream.close();
//}