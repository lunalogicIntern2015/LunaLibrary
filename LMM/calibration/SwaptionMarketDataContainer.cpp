#include <LMM/calibration/SwaptionMarketDataContainer.h>

#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <cassert>
#include <algorithm>

#include <LMM/helper/GenericPath.h>
#include <LMM/helper/BuildVariable.h>
#include <LMM/helper/Noise.h>
#include <Numeric/Interpolation.h>

SwaptionMarketDataContainer::SwaptionMarketDataContainer()
: swaption_sparse_step_ (-1000000)
, strike_bump_ (-10000000000)
{
}

size_t SwaptionMarketDataContainer::get_nbSWAPTION()         const { return SWAPTION_VECTOR_.size()  ; }
size_t SwaptionMarketDataContainer::get_nbLIBOR()            const { return LIBOR_STARTDATES_.size() ; }
size_t SwaptionMarketDataContainer::get_SwaptionSparseStep() const { return swaption_sparse_step_;     } 
double SwaptionMarketDataContainer::get_StrikeBump()          const{ return strike_bump_;              }

const std::vector<double>&               SwaptionMarketDataContainer::get_LIBOR_INIT()        const { return LIBOR_INIT_        ; }
//const std::vector<double>&               SwaptionMarketDataContainer::get_LIBOR_SHIFT()       const { return LIBOR_SHIFT_       ; }
const std::vector<double>&               SwaptionMarketDataContainer::get_ZC_BOND()           const { return ZC_BOND_           ; }
const std::vector<double>&               SwaptionMarketDataContainer::get_ZC_MATURITIES()     const { return ZC_MATURITIES_     ; }
const std::vector<double>&               SwaptionMarketDataContainer::get_LIBOR_STARTDATES()  const { return LIBOR_STARTDATES_  ; }
const SwaptionMarketDataContainer::RealVector&     SwaptionMarketDataContainer::get_NUMERAIRE()         const { return NUMERAIRE_  ;        }

const SwaptionMarketDataContainer::RealVector&     SwaptionMarketDataContainer::get_STRIKE_VECTOR()   const { return STRIKE_VECTOR_   ; }
const SwaptionMarketDataContainer::RealVector&     SwaptionMarketDataContainer::get_VOLATILITY_VECTOR() const { return VOLATILITY_VECTOR_ ; }
const SwaptionMarketDataContainer::SwaptionVector& SwaptionMarketDataContainer::get_SWAPTION_VECTOR()   const { return SWAPTION_VECTOR_   ; }
const std::vector<size_t>&                       SwaptionMarketDataContainer::get_SWPN_MATURITY_INDICES() const { return SWPN_MATURITY_INDICES_ ; }
const SwaptionMarketDataContainer::VectorOfMatrixIndices& SwaptionMarketDataContainer::get_HGVOL_VECTOR_INDICES()   const { return HGVOL_VECTOR_INDICES_ ; }

const SwaptionMarketDataContainer::RealMatrix    &        SwaptionMarketDataContainer::get_STRIKE_MATRIX()         const { return STRIKE_MATRIX_; }
const SwaptionMarketDataContainer::SwaptionMatrix& SwaptionMarketDataContainer::get_SWPN_MATRIX()           const { return SWPN_MATRIX_; }
const SwaptionMarketDataContainer::RealMatrix    &        SwaptionMarketDataContainer::get_MKT_VOL_MATRIX()        const { return MKT_VOL_MATRIX_; }
const SwaptionMarketDataContainer::MatrixOfMatrixIndices& SwaptionMarketDataContainer::get_HGVOL_NODE_MAPPING()   const { return HGVOL_NODE_MAPPING_ ; }

void SwaptionMarketDataContainer::clear_all_data()
{
	clear_all_LIBOR_ZC()	  ;

	clear_all_SWAPTION_DATA() ;
}

void SwaptionMarketDataContainer::clear_all_LIBOR_ZC()
{
	LIBOR_INIT_.clear()       ;
	ZC_BOND_.clear()          ;
	ZC_MATURITIES_.clear()    ;
	LIBOR_STARTDATES_.clear() ;

	NUMERAIRE_.clear() ;
}

void SwaptionMarketDataContainer::compute_numeraire()
{
	NUMERAIRE_.clear() ;

	for(size_t i=0; i<ZC_BOND_.size(); ++i)
	{
		NUMERAIRE_.push_back( 1/ZC_BOND_[i] ) ;
	}
}

void SwaptionMarketDataContainer::clear_all_SWAPTION_DATA()
{
	swaption_sparse_step_ =-1000000;

	STRIKE_VECTOR_.clear()         ;
	VOLATILITY_VECTOR_.clear()     ;		
	SWAPTION_VECTOR_.clear()       ;
	HGVOL_VECTOR_INDICES_.clear()  ;
	SWPN_MATURITY_INDICES_.clear() ;

	clear_all_Matrix_data() ;
}

void SwaptionMarketDataContainer::refresh_AllSwaptionStrike()
{
	size_t nbSwaption = STRIKE_VECTOR_.size();
	assert(nbSwaption == SWAPTION_VECTOR_.size() );
	for(size_t iSwaption=0;iSwaption<nbSwaption;++iSwaption)
	{
		SWAPTION_VECTOR_[iSwaption].set_strike(STRIKE_VECTOR_[iSwaption]);
	}
}

void SwaptionMarketDataContainer::clear_all_Matrix_data() 
{
	for(size_t iRow=0;iRow<HGVOL_NODE_MAPPING_.size();++iRow)
	{
		HGVOL_NODE_MAPPING_[iRow].clear();
	}
	HGVOL_NODE_MAPPING_.clear();

	for(size_t iRow=0;iRow<SWPN_MATRIX_.size();++iRow)
	{
		SWPN_MATRIX_[iRow].clear();
	}
	SWPN_MATRIX_.clear();

	for(size_t iRow=0;iRow<MKT_VOL_MATRIX_.size();++iRow)
	{
		MKT_VOL_MATRIX_[iRow].clear();
	}
	MKT_VOL_MATRIX_.clear();
}

void SwaptionMarketDataContainer::build_MatrixDataFromVectorData()
{
	//! creation maturities indices that match with Swaptions in the container
	//! creation pair of indices in vol g matrix, each line correspond to a maturity
	for(size_t iRow=0;iRow<SWPN_MATURITY_INDICES_.size();++iRow)
	{
		VectorOfMatrixIndices row_vector_indices;
		SwaptionVector        row_vector_swaption;
		RealVector            row_vector_mkt_vol;
		RealVector            row_vector_strike;

		size_t maturity_swaption_index = SWPN_MATURITY_INDICES_[iRow];

		// search for swaption having lower the maturity 
		for(size_t iSwaption = 0;iSwaption<SWAPTION_VECTOR_.size();++iSwaption)
		{
			const VanillaSwaption & swaption = SWAPTION_VECTOR_[iSwaption];
			const VanillaSwap& swap = swaption.getUnderlyingSwap();

			if(swap.get_indexStart() == maturity_swaption_index )// check if swap  have lower the maturity 
			{
				row_vector_swaption.push_back(swaption);
				row_vector_mkt_vol.push_back(VOLATILITY_VECTOR_[iSwaption]);
				row_vector_strike.push_back(swaption.get_strike());

				size_t end_swaption_index = swap.get_indexEnd();

				assert( (end_swaption_index-maturity_swaption_index)%swaption_sparse_step_ == 0 );
				for(size_t liborIndex =swaption_sparse_step_ ; liborIndex < end_swaption_index ; liborIndex+=swaption_sparse_step_)// not take in account the last libor dependant
				{
					size_t max_time_index = std::min(maturity_swaption_index,liborIndex);
					for(size_t timeIndex =swaption_sparse_step_ ; timeIndex <=  max_time_index ; timeIndex+=swaption_sparse_step_)// not take in account greater than maturity time
					{
						std::pair<size_t,size_t> indices_pair(liborIndex,timeIndex);

						//! push if this index pair is not already added
						bool indices_already_added = false ;
						// search if indices pair already exist in the current row
						if( std::find(row_vector_indices.begin(), row_vector_indices.end(), indices_pair ) !=row_vector_indices.end()  )
						{
							indices_already_added = true;
						}

						// search if indices pair already exist in the precedent rows
						for(size_t i_line=0;i_line<HGVOL_NODE_MAPPING_.size();++i_line)
						{
							const VectorOfMatrixIndices & line_indices = HGVOL_NODE_MAPPING_[i_line];
							if(std::find(line_indices.begin(), line_indices.end(), indices_pair ) !=line_indices.end() )
							{
								indices_already_added = true;
							}
						}

						//! push if this index pair IS NOT already added
						if(!indices_already_added) row_vector_indices.push_back(indices_pair);
					}	
				}
			}
		}

		assert( !row_vector_indices.empty() );// for a stored maturity, a swaption and dependant volatilities has to be found --> not empty

		STRIKE_MATRIX_.push_back(row_vector_strike);
		MKT_VOL_MATRIX_.push_back(row_vector_mkt_vol);
		SWPN_MATRIX_.push_back(row_vector_swaption);
		HGVOL_NODE_MAPPING_.push_back(row_vector_indices);		
	}
}

bool SwaptionMarketDataContainer::is_ATMSwaptionData(const LmmVanillaSwapPricer_PTR swap_pricer_ptr) const 
{
	// todo check the LMMTenorStructure dates coherent with ZC maturities dates
	size_t nbSwaption = SWAPTION_VECTOR_.size();
	for(size_t iSwaption =0;iSwaption<nbSwaption;++iSwaption)
	{
		const double strike    = SWAPTION_VECTOR_[iSwaption].get_strike();
		const double swap_rate = swap_pricer_ptr->swapRate_Analytical(SWAPTION_VECTOR_[iSwaption].getUnderlyingSwap() , this ->get_LIBOR_INIT() ) ;
		if(strike != swap_rate)
		{
			std::cout<<LMM::NOTIF_MSG<<" SwaptionMarketDataContainer "<<iSwaption<<"th swaption strike is not ATM swaption"<<std::endl;
			return false;
		}
	}
	return true;
}

bool SwaptionMarketDataContainer::check_data_consistance() const 
{

	// check if interest rate data vector a coherent in terms of size
	//if (LIBOR_INIT_.size() !=  LIBOR_SHIFT_.size()      ) data_is_consistant = false;
	if (LIBOR_INIT_.size() !=  LIBOR_STARTDATES_.size() ) return false;
	if (ZC_BOND_.size()    !=  ZC_MATURITIES_.size()    ) return false;
	if (ZC_BOND_.size()    !=  NUMERAIRE_.size()        ) return false;
	if (ZC_BOND_.size()    !=  LIBOR_INIT_.size()+1     ) return false;

	// check if swaptions data vector a coherent in terms of size
	size_t nbSwaption = SWAPTION_VECTOR_.size();
	if( STRIKE_VECTOR_.size()         != nbSwaption ) return false;
	if( VOLATILITY_VECTOR_.size()     != nbSwaption ) return false;
	if( HGVOL_VECTOR_INDICES_.size()  != nbSwaption ) return false;
	// check if strike vector is coherent with strike in swaptions
	for(size_t iSwaption =0;iSwaption<nbSwaption;++iSwaption)
	{
		const double  swaption_strike = SWAPTION_VECTOR_[iSwaption].get_strike();
		if(swaption_strike != STRIKE_VECTOR_[iSwaption])
		{
			std::cout<<LMM::ERROR_MSG<<" SwaptionMarketDataContainer "<<iSwaption<<"th swaption strike is not coherent"<<std::endl;
			return false;
		}
	}

	size_t nbSwaptionMaturity = SWPN_MATURITY_INDICES_.size();
	if( STRIKE_MATRIX_.size()      != nbSwaptionMaturity ) return false;
	if( MKT_VOL_MATRIX_.size()     != nbSwaptionMaturity ) return false;
	if( SWPN_MATRIX_.size()        != nbSwaptionMaturity ) return false;
	if( HGVOL_NODE_MAPPING_.size() != nbSwaptionMaturity ) return false;
	
	//Check if every swaption data matrix has good format
	for(size_t iRow=0;iRow<nbSwaptionMaturity;++iRow)
	{
		size_t nbSwaptionOnRow = SWPN_MATRIX_[iRow].size();
		if(STRIKE_MATRIX_[iRow].size()      != nbSwaptionOnRow) return false;
		if(MKT_VOL_MATRIX_[iRow].size()     != nbSwaptionOnRow) return false;
		if(HGVOL_NODE_MAPPING_[iRow].size() != nbSwaptionOnRow) return false;
	}

	return true;
}

void SwaptionMarketDataContainer::print(const std::string& filename) const
{
	std::string path_FileOut = LMMPATH::get_output_path() + filename;

	std::ofstream outputFile  ;  outputFile.open(path_FileOut.c_str() );

	{// printout zero coupon
		std::stringstream zc_dates  ; 
		std::stringstream zc_values ; 
		outputFile << "ZERO COUPON , " << std::endl;
		assert(ZC_MATURITIES_.size() == ZC_BOND_.size() ) ;
		printToStreamZCDate(zc_dates)   ;
		printToStreamZCValue(zc_values) ;
		zc_dates  <<std::endl; 	zc_values <<std::endl; 
		outputFile << zc_dates.str() << zc_values.str();	
	}

	{// printout numeraire
		outputFile << "NUMERAIRE ," << std::endl;
		for (auto & num_bb : NUMERAIRE_) outputFile <<num_bb << ",	";
		outputFile << std::endl<< std::endl;
	}

	{// printout libors init
		std::stringstream libors_dates  ; 
		std::stringstream libors_values ; 
		outputFile << "LIBORS Init, " << std::endl;
		assert(LIBOR_STARTDATES_.size() == LIBOR_INIT_.size() ) ;
		printToStreamLiborDate(libors_dates);
		printToStreamLiborValue(libors_values);
		libors_dates  <<std::endl; 
		libors_values <<std::endl<<std::endl<<std::endl; 
		outputFile << libors_dates.str() << libors_values.str();
	}

	{// Strike Matrix
		outputFile <<"STRIKE, "<<std::endl;
		std::stringstream strike_matrix  ; 
		printToStreamMatrixData(strike_matrix, STRIKE_MATRIX_);

		outputFile<<strike_matrix.str()<<std::endl<<std::endl;
	}

	{// Implied Vol Matrix
		outputFile <<"BLACK VOLs, "<<std::endl;
		std::stringstream black_vol_matrix  ; 
		printToStreamMatrixData(black_vol_matrix, MKT_VOL_MATRIX_);

		outputFile<<black_vol_matrix.str()<<std::endl<<std::endl;
	}

	{// Swaption Indices
		outputFile <<"HG NODE MAPPING, "<<std::endl;
		std::stringstream indices_matrix  ; 
		printToStreamMatrixData(indices_matrix, HGVOL_NODE_MAPPING_);

		outputFile<<indices_matrix.str()<<std::endl;
	}

	outputFile << std::endl;

	outputFile.close();	
}

void SwaptionMarketDataContainer::printToStreamLiborDate(std::stringstream & stream) const
{
	for (size_t iLibor=0;iLibor< LIBOR_STARTDATES_.size();++iLibor)
	{
		stream  <<LIBOR_STARTDATES_[iLibor] << ", ";
	}
}
void SwaptionMarketDataContainer::printToStreamLiborValue(std::stringstream & stream) const
{
	for (size_t iLibor=0;iLibor< LIBOR_INIT_.size();++iLibor)
	{
		stream <<LIBOR_INIT_[iLibor]       << ", ";
	}
}
void SwaptionMarketDataContainer::printToStreamZCDate(std::stringstream & stream) const
{
	for (size_t iZC=0;iZC< ZC_MATURITIES_.size();++iZC)
	{
		stream  <<ZC_MATURITIES_[iZC] << ", ";
	}
}
void SwaptionMarketDataContainer::printToStreamZCValue(std::stringstream & stream) const
{
	for (size_t iZC=0;iZC< ZC_BOND_.size();++iZC)
	{
		stream  <<ZC_BOND_[iZC] << ", ";
	}
}