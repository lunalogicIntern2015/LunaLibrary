#include <LMM/Calibration/Calibrator_CostFunction/LmmGlobal_gCostFunction.h>

LmmGlobal_gCostFunction::LmmGlobal_gCostFunction( LmmVanillaSwaptionApproxPricer_Rebonato_PTR lmmRobonato_ptr                // pricer  
												 , LiborQuotes_ConstPTR liborQuotes_ptr
												 , UpperTriangleVanillaSwaptionQuotes_ConstPTR pUpperTriangleVanillaSwaptionQuotes // instrument to calibrate 
												 , GMatrixMapping_PTR gMatrixMapping_ptr
												 , Shifted_HGVolatilityParam_PTR pShifted_HGVolatilityParam
												 , LmmPenalty_PTR penalty_ptr)
												 : LmmBaseCostFunction(liborQuotes_ptr , pUpperTriangleVanillaSwaptionQuotes)
												 , pLmmVanillaSwaptionApproxPricer_Rebonato_(lmmRobonato_ptr)
												 , pGMatrixMapping_(gMatrixMapping_ptr)
												 , buffer_Shifted_HGVolatilityParam_(pShifted_HGVolatilityParam)
												 , pLmmPenalty_(penalty_ptr)								 
{
	size_t N = upperTriangleVanillaSwaptionWeights_.size1() - 2 ; // nb Swaption on first row in swaption matrix
	nbSwaptions_ =( N * (N+1) ) /2;
}

void LmmGlobal_gCostFunction::reset_CalibrationParams(const QuantLib::Array & params) const 
{
	pGMatrixMapping_->reset_gDelegate(params); 
	buffer_Shifted_HGVolatilityParam_->reset_g_matrix( pGMatrixMapping_->get_g_Ref() );
	pLmmVanillaSwaptionApproxPricer_Rebonato_->update_VolatilityParam( buffer_Shifted_HGVolatilityParam_ );
}

//const Array& param_array
QuantLib::Disposable<QuantLib::Array> LmmGlobal_gCostFunction::values(const QuantLib::Array& param_array) const
{
	// NOTE : 
	// Costfunction::values has to compute only values, not absolute values, nor squared values
	// see Introduction to Selected Classes of QuantLib II _ Dimitri Reiswich 2010 p.57/148

	this->reset_CalibrationParams(param_array);

	this->fully_UpdateSwaptionMdlValues();

	size_t nbPenalties = calc_nbPenalties();	

	//! diff upperTriangleVanillaSwaptionMktQuotes_ - upperTriangleVanillaSwaptionMdlValues_
	QuantLib::Array diff_cost( nbSwaptions_ + nbPenalties );
	size_t vector_counter=0;

	size_t swaption_matrix_size = upperTriangleVanillaSwaptionMdlValues_.size1();
	const UpperTriangularVanillaSwaptionQuotes & mkt_swaption_quote = pUpperTriangleVanillaSwaptionQuotes_->get_UpperTriangularVanillaSwaptionQuotes() ;

	for(size_t iExperity=1; iExperity<swaption_matrix_size-1 ; ++ iExperity)
	{
		const size_t maxTenorBound =  swaption_matrix_size - iExperity;

		for(size_t jTenor=1;jTenor<maxTenorBound;++jTenor)
		{
			const double diff = upperTriangleVanillaSwaptionMdlValues_(iExperity,jTenor) - mkt_swaption_quote(iExperity,jTenor).second ;
			diff_cost[vector_counter]= upperTriangleVanillaSwaptionWeights_(iExperity,jTenor) * diff;//*diff; // To check with quantLib. --> checked
			++vector_counter;
		}
	}


	// calculate for penalty
	std::vector<double> penalty;
	if(!pLmmPenalty_->isEmpty() )
		penalty=pLmmPenalty_->get_penalty(pGMatrixMapping_);

	assert(penalty.size() == nbPenalties);

	for(size_t iPel=0;iPel<penalty.size();++iPel)
	{
		diff_cost[vector_counter] = penalty[iPel];
		++vector_counter;		
	}

	// check if penalties fully fill the cost vector
	assert(vector_counter == diff_cost.size() ); 

	if( breakForPrintOut() ) 
		std::cout<<std::endl<<LMM::NOTIF_MSG<<" iter."<<nbCalled;

	if(LMM::DEUBGLMM())
	{
		if( breakForPrintOut() ) 
		{
			std::cout<<std::endl;

			if(buffer_calib_reference.empty() )
			{
				std::cout <<"      x = "<<param_array<<std::endl<<std::endl;
			}
			else
			{
				QuantLib::Array error = error_calib(param_array);
				std::cout <<"diff   = "<<error<<std::endl<<std::endl;
			}
			std::cout <<"f(x,p) = "<<diff_cost<<std::endl;			
		}
	}

	++nbCalled;
	return diff_cost;
}

void LmmGlobal_gCostFunction::fully_UpdateSwaptionMdlValues() const // pricing
{
	//! upperTriangle:

	const UpperTriangularVanillaSwaptionQuotes & mkt_swaption_quotes = pUpperTriangleVanillaSwaptionQuotes_->get_UpperTriangularVanillaSwaptionQuotes() ;
	const std::vector<double>& libor_init = pLiborQuotes_->get_InitLibor();

	size_t swaption_matrix_size = upperTriangleVanillaSwaptionMdlValues_.size1();
	for(size_t iExperity=1; iExperity<swaption_matrix_size-1 ; ++ iExperity)
	{
		const size_t maxTenorBound =  swaption_matrix_size - iExperity;

		for(size_t jTenor=1;jTenor<maxTenorBound;++jTenor)
		{
			const VanillaSwaption & swaption = mkt_swaption_quotes(iExperity,jTenor).first;

			const double blackVol_ij = pLmmVanillaSwaptionApproxPricer_Rebonato_->volBlack(swaption,libor_init);

			upperTriangleVanillaSwaptionMdlValues_(iExperity,jTenor)=blackVol_ij;
		}
	}
}

QuantLib::Array LmmGlobal_gCostFunction::error_calib(const QuantLib::Array & actual_param) const
{
	size_t nbParam = actual_param.size();
	QuantLib::Array error_calib(nbParam,0.);
	for(size_t i=0;i<nbParam;++i)
	{
		error_calib[i]  = actual_param[i] - buffer_calib_reference[i] ;
	}
	return error_calib;
}

size_t LmmGlobal_gCostFunction::calc_nbPenalties() const 
{
	return pLmmPenalty_->get_TotalNbPenalty(pGMatrixMapping_);
} 

void LmmGlobal_gCostFunction::print(const std::string& filename) const
{
	std::string path_OutPut = LMMPATH::get_output_path() + filename;

	{
		std::vector<PrintElement_PTR> elements_print;
		PrintElement_PTR swpm_weight_print = PrintElement_PTR(new MatrixPrintElement<UpperTriangularDoubleMatrix>("Swaption Weights",  upperTriangleVanillaSwaptionWeights_));

		elements_print.push_back(swpm_weight_print);

		Printer printer(path_OutPut, elements_print);
		printer.print(); 
	}
}

void LmmGlobal_gCostFunction::reset_TrueParam(const QuantLib::Array & true_param)
{
	const size_t nbParam = true_param.size();
	buffer_calib_reference.resize(nbParam);
	for(size_t i=0;i<nbParam;++i)
	{
		buffer_calib_reference[i]=true_param[i];
	}
}

double LmmGlobal_gCostFunction::calc_PenaltyTimeHomogeneity_L2() const 
{
	return pLmmPenalty_->calc_PenaltyTimeHomogeneity_L2(pGMatrixMapping_);
}
double LmmGlobal_gCostFunction::calc_PenaltyTimeHomogeneity_L1() const 
{
	return pLmmPenalty_->calc_PenaltyTimeHomogeneity_L1(pGMatrixMapping_);
}
double LmmGlobal_gCostFunction::calc_PenaltyTimeHomogeneity_L_INF() const 
{
	return pLmmPenalty_->calc_PenaltyTimeHomogeneity_L_INF(pGMatrixMapping_);
}

double LmmGlobal_gCostFunction::calc_PenaltySmoothMaturity_L2() const 
{
	return pLmmPenalty_->calc_PenaltySmoothMaturity_L2(pGMatrixMapping_);
}
double LmmGlobal_gCostFunction::calc_PenaltySmoothMaturity_L1() const 
{
	return pLmmPenalty_->calc_PenaltySmoothMaturity_L1(pGMatrixMapping_);
}
double LmmGlobal_gCostFunction::calc_PenaltySmoothMaturity_L_INF() const 
{
	return pLmmPenalty_->calc_PenaltySmoothMaturity_L_INF(pGMatrixMapping_);
}