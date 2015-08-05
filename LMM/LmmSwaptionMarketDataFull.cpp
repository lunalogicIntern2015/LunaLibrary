#include <LMM/LmmSwaptionMarketDataFull.h>

#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <cassert>
#include <algorithm>

#include <LMM/Helper/GenericPath.h>
#include <LMM/Helper/BuildVariable.h>
#include <LMM/Helper/Printer.h>

#include <Numeric/Interpolation.h>
#include <LMM/Pricer/LmmAnalyticalPricer/LmmVanillaSwapPricer.h>

/* input file structure:
- ATM quotes
- ATM + 5bp
- ATM - 5bp
- ATM + 50 bp
- ATM + 100 bp
- ATM + 200 bp
- ATM - 50 bp
- ATM - 100 bp
- ATM - 200 bp 
- skew 
*/

LmmSwaptionMarketDataFull::LmmSwaptionMarketDataFull(const Tenor& fixedtenor, 
													 const Tenor& floattenor, const size_t maxNbYear)
	: LmmSwaptionMarketData(fixedtenor, floattenor, maxNbYear) {}

void LmmSwaptionMarketDataFull::parseFromMarketData(const std::string& filename)
{
	mkt_file_name.clear();
	mkt_file_name = filename ;

	std::ifstream instream;
	std::string input_file_fullpath = LMMPATH::get_runtime_datapath() + filename;
	instream.open(input_file_fullpath.c_str());

	getGeneralInfo(instream);

	std::vector<double> df_dates_buffer;
	std::vector<double> df_values_buffer;
	getDF_MarketData(instream, df_dates_buffer, df_values_buffer);
	interpolateAndBuildLiborQuotes(df_dates_buffer, df_values_buffer);

	// If a curve of forward rate quote is specified, use this curve instead of interpolation from Discount Factor curve
	std::string text_line;
	getline(instream,text_line);
	if(!isLineEndOfBlock(text_line))
	{
		std::vector<double> libor_quote_values;
		parseLIBORS_Values(text_line,libor_quote_values);
		//pLiborQuotes_->reset_Libor(libor_quote_values); just for testing if the data taken from the bloomberg forward curve analysis is better, but it is not. 
	}


	for(size_t itSWPMTab=0 ; itSWPMTab<9 ; ++itSWPMTab)
	{
		double strike_bump = -100000000;
		std::vector<double> atm_swpm_experities_buffer;
		std::vector<double> atm_swpm_tenor_buffer;
		std::vector< std::vector<double> > atm_swpm_quote_buffer;
		std::vector< std::vector<double> > atm_swpm_strike_buffer;
		getSwaptionMarketData(instream,strike_bump, atm_swpm_experities_buffer, atm_swpm_tenor_buffer, atm_swpm_quote_buffer, atm_swpm_strike_buffer);

		if(itSWPMTab==0)// first block is ATM Swaption VCUB
		{
			getSwaptionQuoteFirstColumn(atm_swpm_experities_buffer, atm_swpm_quote_buffer);
		}

		buildSwationQuotes( strike_bump, atm_swpm_experities_buffer, atm_swpm_tenor_buffer, atm_swpm_quote_buffer, atm_swpm_strike_buffer );
	}

	swaptionMarketData_ATM_.first->set_Data_FileName(filename);

	buildSwationQuotesVol_skew();

	swaptionMarketData_ATMpp_.first->set_Data_FileName(filename);
	swaptionMarketData_ATMmm_.first->set_Data_FileName(filename);

	instream.close();
}

void LmmSwaptionMarketDataFull::buildSwationQuotes( const double & strike_bump
											   , const std::vector<double>& mkt_expirities
											   , const std::vector<double>& mkt_tenor
											   , const std::vector<std::vector<double> >& mkt_quote
											   , const std::vector<std::vector<double> >& mkt_strike)
{
	size_t fix_float_ratio = fixedtenor_.ratioTo(floattenor_);
	size_t lastLiborIndex = pLMMTenorStructure_->get_horizon();

	size_t last_year = lastLiborIndex / fix_float_ratio ;
	size_t matrix_size = last_year + 1 ;

	UpperTriangularDoubleMatrix strike_rate_matrix(matrix_size,matrix_size);
	numeric::FullMatrix interpolated_black_vol_matrix(matrix_size,matrix_size); // use fully the quotation black vol matrix for interpolation,

	// stored all existing market swaption position (parsed) in teh swaption matrix
	std::vector<std::pair<size_t, size_t> > exist_mkt_swpm_position;

	for(size_t iMktExperity=0;iMktExperity<mkt_expirities.size();++iMktExperity)
	{
		const double mkt_experity_i = mkt_expirities[iMktExperity]; 
		LMM::Index swpm_startIndex  = pLMMTenorStructure_->get_Index(mkt_experity_i);

		if(swpm_startIndex>=fix_float_ratio)
		{
			for(size_t jMktTenor=0;jMktTenor<mkt_tenor.size();++jMktTenor)
			{
				const double mkt_tenor_j = mkt_tenor[jMktTenor];
				const double swpm_end_date = mkt_experity_i + mkt_tenor_j;

				// put the market strike into strike matrix (upper)
				// missed element will be calculated by analytical swap rate formula
				LMM::Index swpm_endIndex   = pLMMTenorStructure_->get_Index(swpm_end_date);
				if( (swpm_endIndex%fix_float_ratio==0) && ( (swpm_endIndex-swpm_startIndex)%fix_float_ratio==0 ) && swpm_endIndex>fix_float_ratio)
				{
					const size_t iRow = swpm_startIndex/fix_float_ratio;
					const size_t jCol = (swpm_endIndex - swpm_startIndex)/fix_float_ratio;
					exist_mkt_swpm_position.push_back( std::pair<size_t,size_t>(iRow,jCol) );

					strike_rate_matrix(iRow,jCol)=mkt_strike[iMktExperity][jMktTenor];					
				}

				// put the market quote into quote matrix (full)
				// missed element will be interpolated
				LMM::Index swpm_tenorIndex   = pLMMTenorStructure_->get_Index(mkt_tenor_j);
				if (swpm_tenorIndex%fix_float_ratio==0)
				{
					const size_t iRow = swpm_startIndex/fix_float_ratio;
					const size_t jCol = swpm_tenorIndex/fix_float_ratio;
				
					interpolated_black_vol_matrix(iRow,jCol) = mkt_quote[iMktExperity][jMktTenor];					
				}
			}	
		}
		else
		{
			if(LMM::WARNLMM ())
				std::cout<<LMM::WARN_MSG <<" SWPM Experity "<<mkt_experity_i<<" does not match in LMMTenorStructure's Swaption Matrix "<<std::endl;
		}
	}

	//getting the missed position in the swaption matrix
	std::vector<size_t > missed_row_indices;
	std::vector<size_t > missed_col_indices;

	for(size_t iExperityRow=1;iExperityRow<last_year;++iExperityRow)
	{
		std::pair<size_t,size_t> indices_pair_first_col(iExperityRow,1);
		if( std::find(exist_mkt_swpm_position.begin(), exist_mkt_swpm_position.end(),indices_pair_first_col)==exist_mkt_swpm_position.end() ) 
		{
			if( std::find(missed_row_indices.begin(), missed_row_indices.end(),iExperityRow)==missed_row_indices.end() )
				missed_row_indices.push_back(iExperityRow);
		}
	}
	std::sort( missed_row_indices.begin(), missed_row_indices.end() );

	for(size_t jTenorCol=1;jTenorCol<last_year;++jTenorCol)
	{
		std::pair<size_t,size_t> indices_pair_first_row(1,jTenorCol);

		if( std::find(exist_mkt_swpm_position.begin(), exist_mkt_swpm_position.end(),indices_pair_first_row)==exist_mkt_swpm_position.end() ) 
		{
			if( std::find(missed_col_indices.begin(), missed_col_indices.end(),jTenorCol)==missed_col_indices.end() )
				missed_col_indices.push_back(jTenorCol);
		}
	}
	std::sort( missed_col_indices.begin(), missed_col_indices.end() );
	
	fill_missing_strike(strike_rate_matrix,strike_bump,missed_row_indices,missed_col_indices);
	
	numeric::Interpolation interpolator;
	interpolator.fullMatrixInterpolate(interpolated_black_vol_matrix,missed_row_indices,missed_col_indices);

	//assigne the upper triangular part of interpolated matrix to a upper black vol matrix
	UpperTriangularDoubleMatrix black_vol_matrix(matrix_size,matrix_size);
	for(size_t iRow=1;iRow<interpolated_black_vol_matrix.size1();++iRow)
	{
		for(size_t jCol=1;jCol<interpolated_black_vol_matrix.size2()-iRow;++jCol)
		{
			black_vol_matrix(iRow,jCol)= interpolated_black_vol_matrix(iRow,jCol);		
		}
	}

	switch(size_t(strike_bump)) 
	{
		case 0:{		
			swaptionMarketData_ATM_.first.reset(new UpperTriangleVanillaSwaptionQuotes( pLMMTenorStructure_
					, last_year, fixedtenor_, floattenor_, strike_rate_matrix, black_vol_matrix)
												);
			swaptionMarketData_ATM_.second=strike_bump;
			break ;
				}
		case 5 :{	
			swaptionMarketData_ATMpp_.first.reset(new UpperTriangleVanillaSwaptionQuotes( pLMMTenorStructure_
					, last_year, fixedtenor_, floattenor_, strike_rate_matrix, black_vol_matrix)
												);
			swaptionMarketData_ATMpp_.second=strike_bump;
			break ;
				}
		case -5 :{	
			swaptionMarketData_ATMmm_.first.reset(new UpperTriangleVanillaSwaptionQuotes( pLMMTenorStructure_
					, last_year, fixedtenor_, floattenor_, strike_rate_matrix, black_vol_matrix)
												);
			swaptionMarketData_ATMmm_.second=strike_bump;
			break;
				}
		case 50 :{	
			std::cout << "lecture ATM + 50bp" << std::endl ;
			swaptionMarketData_ATMp50_.first.reset(new UpperTriangleVanillaSwaptionQuotes( pLMMTenorStructure_
					, last_year, fixedtenor_, floattenor_, strike_rate_matrix, black_vol_matrix)
												);
			swaptionMarketData_ATMp50_.second=strike_bump;
			break ;
				}
		case 100 :{	
			std::cout << "lecture ATM + 100bp" << std::endl ;
			swaptionMarketData_ATMp100_.first.reset(new UpperTriangleVanillaSwaptionQuotes( pLMMTenorStructure_
					, last_year, fixedtenor_, floattenor_, strike_rate_matrix, black_vol_matrix)
												);
			swaptionMarketData_ATMp100_.second=strike_bump;
			break;
				}
		case 200 :{	
			std::cout << "lecture ATM + 200bp" << std::endl ;
			swaptionMarketData_ATMp200_.first.reset(new UpperTriangleVanillaSwaptionQuotes( pLMMTenorStructure_
					, last_year, fixedtenor_, floattenor_, strike_rate_matrix, black_vol_matrix)
												);
			swaptionMarketData_ATMp200_.second=strike_bump;
			break;
				}
		case -50 :{	
			swaptionMarketData_ATMm50_.first.reset(new UpperTriangleVanillaSwaptionQuotes( pLMMTenorStructure_
					, last_year, fixedtenor_, floattenor_, strike_rate_matrix, black_vol_matrix)
												);
			swaptionMarketData_ATMm50_.second=strike_bump;
			break ;
				}
		case -100 :{	
			swaptionMarketData_ATMm100_.first.reset(new UpperTriangleVanillaSwaptionQuotes( pLMMTenorStructure_
					, last_year, fixedtenor_, floattenor_, strike_rate_matrix, black_vol_matrix)
												);
			swaptionMarketData_ATMm100_.second=strike_bump;
			break;
				}
		case -200 :{	
			swaptionMarketData_ATMm200_.first.reset(new UpperTriangleVanillaSwaptionQuotes( pLMMTenorStructure_
					, last_year, fixedtenor_, floattenor_, strike_rate_matrix, black_vol_matrix)
												);
			swaptionMarketData_ATMm200_.second=strike_bump;
			break;
				}
		default :
			throw "Incorrect strike bump in MarketData" ;
	}

}
