#include <Cheyette/Calibration/CoTerminalSwaptionVol.h>


//
//void UpperTriangleVanillaSwaptionQuotes::print(const std::string& filename, const bool erase_file) const
//{
//	std::string path_OutPut = LMMPATH::get_output_path() + filename;
//
//	{
//		std::vector<PrintElement_PTR> elements_print;
//		UpperTriangularIndexPairMatrix swap_indices_matrix = get_UpperTriangularSwaptionIndexMatrix();	
//		PrintElement_PTR swapMatrix_print = PrintElement_PTR(new MatrixPrintElement<UpperTriangularIndexPairMatrix>("swaps",  swap_indices_matrix));
//		PrintElement_PTR mapMatrix_print = PrintElement_PTR(new MatrixPrintElement<UpperTriangularIndexPairMatrix>("mapping Indices",  indexMapping_gDelegate_gTransformed_));
//
//		elements_print.push_back(swapMatrix_print);
//		elements_print.push_back(mapMatrix_print);
//
//		Printer printer(path_OutPut, elements_print);
//		printer.print(erase_file);
//	}
//
//
//	{
//		std::ofstream data_stream ;
//		data_stream.open(path_OutPut.c_str(), std::ios::app);
//		data_stream<<std::endl<<"--------"<<std::endl;
//		data_stream.close();
//	}
//
//	{
//		std::vector<PrintElement_PTR> elements_print;
//
//		UpperTriangularDoubleMatrix quote_matrix  = get_UpperTriangularQuoteValues() ;
//		UpperTriangularDoubleMatrix strike_matrix = get_UpperTriangularStrike() ;
//		PrintElement_PTR quoteMatrix_print = PrintElement_PTR(new MatrixPrintElement<UpperTriangularDoubleMatrix>("Quote",  quote_matrix));
//		PrintElement_PTR strikeMatrix_print = PrintElement_PTR(new MatrixPrintElement<UpperTriangularDoubleMatrix>("Strike",  strike_matrix));
//
//		elements_print.push_back(quoteMatrix_print);
//		elements_print.push_back(strikeMatrix_print);
//		Printer printer(path_OutPut, elements_print);
//		printer.print(false);
//	}
//}
//
//
//UpperTriangleVanillaSwaptionQuotes_ConstPTR UpperTriangleVanillaSwaptionQuotes::create_ATMSwaptionImpliedVol
//	(
//	LiborQuotes_ConstPTR libor_quotes_ptr,
//	const Tenor&  fixedTenor,
//	const Tenor&  floatingTenor,
//	LmmVanillaSwaptionApproxPricer_Rebonato_PTR black_vol_approx_ptr,
//	const double&  strike_bump
//	)
//{
//	LMMTenorStructure_PTR pLMMTenorStructure = libor_quotes_ptr->get_LMMTenorStructure_PTR() ;
//	size_t fix_float_ratio = fixedTenor.ratioTo(floatingTenor);
//	size_t lastLiborIndex  = pLMMTenorStructure->get_nbLIBOR()-1;
//	assert(lastLiborIndex%fix_float_ratio == 0);
//	size_t last_year = lastLiborIndex / fix_float_ratio ;
//
//	size_t matrix_size = last_year + 1 ;
//	UpperTriangularDoubleMatrix swap_rate_matrix(matrix_size,matrix_size);
//	UpperTriangularDoubleMatrix black_vol_matrix(matrix_size,matrix_size);
//
//	//first row and first column never used
//	for(size_t k=0;k<matrix_size;++k)
//	{
//		swap_rate_matrix(k,0) = -1000000; swap_rate_matrix(0,k) = -1000000;
//		black_vol_matrix(k,0) = -1000000; black_vol_matrix(0,k) = -1000000;
//	}
//
//
//	const std::vector<double> & init_libor = libor_quotes_ptr->get_InitLibor();
//	LmmVanillaSwapPricer swap_pricer(pLMMTenorStructure);
//	for(size_t iMaturity=1;iMaturity<matrix_size;++iMaturity)
//	{
//		for(size_t jTenor=1;jTenor<matrix_size - iMaturity ;++jTenor)
//		{
//			size_t start_swap_index = iMaturity * fix_float_ratio;
//			size_t end_swap_index   = (iMaturity+jTenor) * fix_float_ratio;	
//
//			double empty_strike = -100000000;
//			VanillaSwap swap_ij(empty_strike,start_swap_index,end_swap_index,floatingTenor,fixedTenor,pLMMTenorStructure);
//			double swap_rate_ij = strike_bump+swap_pricer.swapRate_Analytical(swap_ij,init_libor);
//
//			swap_ij.set_strike(swap_rate_ij);
//			swap_rate_matrix(iMaturity,jTenor) = swap_rate_ij;
//
//			VanillaSwaption swaption_ij(swap_ij, OptionType::OptionType::CALL);
//			double black_vol_ij = black_vol_approx_ptr->volBlack(swaption_ij,init_libor);
//
//			black_vol_matrix(iMaturity,jTenor)=black_vol_ij;
//		}
//	}
//
//	UpperTriangleVanillaSwaptionQuotes_PTR atm_swaption_implied_vol 
//		(new UpperTriangleVanillaSwaptionQuotes(
//		pLMMTenorStructure,
//		last_year,
//		fixedTenor,
//		floatingTenor,
//		swap_rate_matrix,
//		black_vol_matrix
//		) 
//		);
//
//	atm_swaption_implied_vol->set_Data_FileName("VirtuallyCreatedData");
//
//	return atm_swaption_implied_vol;
//}
//
//
//double UpperTriangleVanillaSwaptionQuotes::get_MinQuote() const 
//{
//	double min_quote(1000000.);
//
//	size_t nbRow = upperTriangleVanillaSwaptionQuotes_.size1();
//	size_t nbCol = upperTriangleVanillaSwaptionQuotes_.size2();
//	assert(nbRow == nbCol);
//
//	for(size_t iExpirity = 1; iExpirity<nbRow; ++iExpirity) // row
//	{
//		for(size_t jTenor = 1; jTenor<nbCol-iExpirity; ++jTenor) // col
//		{
//			const double & quote = upperTriangleVanillaSwaptionQuotes_(iExpirity,jTenor).second;
//			assert(quote > 0) ; // vol quote has to be positive , skew quote to do after
//			if( min_quote > quote )  min_quote = quote ;
//		}
//	}
//
//	return min_quote;
//}
//
//double UpperTriangleVanillaSwaptionQuotes::get_MaxQuote() const 
//{
//	double max_quote(-1000000.);
//
//	size_t nbRow = upperTriangleVanillaSwaptionQuotes_.size1();
//	size_t nbCol = upperTriangleVanillaSwaptionQuotes_.size2();
//	assert(nbRow == nbCol);
//
//	for(size_t iExpirity = 1; iExpirity<nbRow; ++iExpirity) // row
//	{
//		for(size_t jTenor = 1; jTenor<nbCol-iExpirity; ++jTenor) // col
//		{
//			const double & quote = upperTriangleVanillaSwaptionQuotes_(iExpirity,jTenor).second;
//			assert(quote > 0) ; // vol quote has to be positive , skew quote to do after
//			if( max_quote < quote )  max_quote = quote ;
//		}
//	}
//
//	return max_quote ;
//}

