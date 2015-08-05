//#include <LMM/Calibration/Calibrator_CostFunction/CalibrationVolCostFunction.h>
//#include <RNGenerator/McGenerator.h>
//#include <cassert>
//
//#include <boost/random/mersenne_twister.hpp>
//#include <boost/random/uniform_real_distribution.hpp>
//#include <ctime>
//
//
//CalibrationVolCostFunction::CalibrationVolCostFunction
//	(
//	Shifted_HGVolatilityParam_PTR pVolatilityParam,
//	SwaptionMarketDataContainer_PTR pSwaptionMarketDataContainer,
//	const LmmVanillaSwaptionApproxPricer_Rebonato_PTR approxiPricer_ptr,										   
//	//const std::vector<double>& libor_shifts,
//	vector_ swaptionweights,
//	matrix_ weights_maturity,
//	matrix_ weights_tenor,
//	matrix_ weights_maturity_2,
//	matrix_ weights_tenor_2)
//	: use_penalty(false)
//	, penalty_coeff_(0)
//	, nbCalled(0)
//	, pHGVolatilityParamBuffer_(pVolatilityParam)
//	, hg_IndicesVector_()
//	, pSwaptionMarketDataContainer_(pSwaptionMarketDataContainer)
//	, approxiPricer_ptr_(approxiPricer_ptr)
//	, swaption_weights_(swaptionweights)
//	, weights_maturity_(weights_maturity)
//	, weights_tenor_(weights_tenor)
//	, weights_maturity_2_(weights_maturity_2)
//	, weights_tenor_2_(weights_tenor_2)
//{
//	assert( pSwaptionMarketDataContainer_->check_data_consistance() ) ;
//
//	// the default use of mkt swaption is full
//	this->reset_swpn_Vector( pSwaptionMarketDataContainer_->get_SWAPTION_VECTOR() );
//	this->reset_mkt_vol_Vector( pSwaptionMarketDataContainer_->get_VOLATILITY_VECTOR() );
//}
//
//CalibrationVolCostFunction::~CalibrationVolCostFunction(){}
//
//
//void CalibrationVolCostFunction::reset_reference_calib(const QuantLib::Array & true_param)
//{
//	size_t nbParam = true_param.size();
//	reference_calib.resize(nbParam);
//	for(size_t i=0;i<nbParam;++i)
//	{
//		reference_calib[i] = true_param[i];
//	}
//}
//
//QuantLib::Array CalibrationVolCostFunction::error_calib(const QuantLib::Array & actual_param) const
//{
//	size_t nbParam = actual_param.size();
//	QuantLib::Array error_calib(nbParam,0.);
//	for(size_t i=0;i<nbParam;++i)
//	{
//		error_calib[i]  = actual_param[i] - reference_calib[i] ;
//	}
//	return error_calib;
//}
//
//bool CalibrationVolCostFunction::breakForPrintOut(unsigned int nbIter) const 
//{
//	return (nbIter%5 ==0);
//}
//
//void CalibrationVolCostFunction::reset_PenaltyCoeff(const double& coeff)
//{
//	penalty_coeff_ = coeff;
//}
//
//void CalibrationVolCostFunction::reset_nbCalled() 
//{
//	nbCalled=0;
//}
//
//unsigned int CalibrationVolCostFunction::get_nbCalled() const
//{
//	return nbCalled;
//}
//
//void CalibrationVolCostFunction::reset_HG_Indices(const VectorOfMatrixIndices& indices_vector)
//{
//	hg_IndicesVector_ = indices_vector;
//
//}
//
//void CalibrationVolCostFunction::reset_swpn_Vector(const SwaptionVector& swaption_vector)
//{
//	swpn_Vector_ = swaption_vector;
//}
//
//void CalibrationVolCostFunction::reset_mkt_vol_Vector(const std::vector<double>& vol_vector)
//{
//	mkt_vol_Vector_=vol_vector;
//}
//
//void CalibrationVolCostFunction::reset_swaption_weights(const std::vector<double>& weight_vector)
//{
//	swaption_weights_ = weight_vector ;
//}
//
//LmmVanillaSwaptionApproxPricer_Rebonato_PTR CalibrationVolCostFunction::get_LmmVanillaSwaptionApproxPricer()
//{
//	return approxiPricer_ptr_ ;
//}
//
//QuantLib::Real CalibrationVolCostFunction::value(const QuantLib::Array& vol_param) const 
//{
//
//	QuantLib::Array squared_diffs = values(vol_param);
//	//Array weightArray = map_MatrixtoArray(swaption_weights_);//ctntodo uncomment in order to reuse this weight
//
//	QuantLib::Real res = 0;
//	for (size_t i = 0; i < squared_diffs.size(); ++i)
//	{
//		res += squared_diffs[i];
//	}
//
//	res = sqrt(res);
//
//	//res += regularisation(vol_param,1,1,1,1);
//	//res += regularisation(vol_param,0.1,0.05,0.05,0.08);	
//
//	return res;
//}
//
//
//QuantLib::Disposable<QuantLib::Array> CalibrationVolCostFunction::values(const QuantLib::Array& param_array) const
//{
//	//if( !hg_IndicesVector_.empty() )
//	pHGVolatilityParamBuffer_->reset_G_FromArray(param_array,hg_IndicesVector_,true);
//	//else
//	//pHGVolatilityParamBuffer_->reset_FromArray(param_array);
//
//	approxiPricer_ptr_->update_VolatilityParam(pHGVolatilityParamBuffer_);
//
//	size_t nbSwaption = swpn_Vector_.size() ;
//
//	size_t col_calibration;
//	size_t sparse_step    ;
//	this->get_CalibrationType(sparse_step,col_calibration,hg_IndicesVector_);
//
//	const std::vector<double> & libor            = pSwaptionMarketDataContainer_->get_LIBOR_INIT()      ;	
//
//	QuantLib::Array penalty_vector;
//	if(use_penalty)
//	{
//		penalty_vector = this->penalty(sparse_step, col_calibration);
//		penalty_vector*=this->penaltyCoeff();
//	}
//
//	// concatenation of objectives values and penalty values into a longer vector
//	size_t nbValues = nbSwaption + penalty_vector.size();  // YY bug
//	//size_t nbValues = nbSwaption;
//
//	QuantLib::Array squared_diff_vol( nbValues );
//	for(size_t iSwaption=0;iSwaption<nbSwaption;++iSwaption)
//	{
//		const VanillaSwaption  & swaption = swpn_Vector_[iSwaption]; 
//		const double           & marketSwaptionVol= mkt_vol_Vector_[iSwaption];
//		double approxVol = approxiPricer_ptr_->volBlack(swaption,libor);
//
//		squared_diff_vol[iSwaption] = swaption_weights_[iSwaption] * std::abs(approxVol*approxVol - marketSwaptionVol*marketSwaptionVol) ;
//	}
//
//	if(use_penalty)
//	{
//		for(size_t iPenalty=0;iPenalty< penalty_vector.size();++iPenalty)
//		{			
//			squared_diff_vol[nbSwaption+iPenalty] = penalty_vector[iPenalty]*penalty_vector[iPenalty];
//		}		
//	}
//
//	//double average_penalty = 0.0;
//	//for(size_t i=0; i<penalty_vector.size(); ++i)
//	//{
//	//	average_penalty += penalty_vector[i];
//	//}
//	//average_penalty /= penalty_vector.size();
//
//	//if(use_penalty)
//	//{
//	//	for(size_t iSwaption=0;iSwaption<nbSwaption;++iSwaption)
//	//	{
//	//		squared_diff_vol[iSwaption] = swaption_weights_[iSwaption] * average_penalty;
//	//	}	
//	//}
//
//
//
//	if(LMM::DEUBGLMM)
//	{
//		if( breakForPrintOut(nbCalled) ) 
//		{
//			std::cout<<LMM::NOTIF_MSG<<" iter."<<nbCalled  <<"      penalty_coeff = "<<this->penaltyCoeff()<<std::endl<<std::endl;
//			std::cout                <<"     p = "; if(use_penalty) std::cout<<penalty_vector<<std::endl; else std::cout<<" [none] "<<std::endl;
//			if(reference_calib.empty() ){
//				std::cout <<"      x = "<<param_array<<std::endl;
//			}
//			else
//			{
//				QuantLib::Array error = error_calib(param_array);
//				std::cout <<"error  = "<<error<<std::endl;
//			}
//			std::cout <<"f(x,p) = "<<squared_diff_vol<<std::endl;			
//		}
//	}
//	++nbCalled;
//
//	return squared_diff_vol;
//}
//
//
//
//double CalibrationVolCostFunction::penaltyCoeff() const 
//{
//	return penalty_coeff_;
//}
//
//QuantLib::Disposable<QuantLib::Array>  CalibrationVolCostFunction::penalty(size_t sparse_step , size_t jjColIndex) const 
//{
//	assert(sparse_step>1);
//	assert(jjColIndex % sparse_step == 0);
//
//	size_t N = pHGVolatilityParamBuffer_->get_horizon();
//	assert(N%sparse_step == 0);
//
//	bool is_calibration_by_column=(jjColIndex != 0);
//	bool is_not_first_sparse_column=(jjColIndex>sparse_step);
//
//	size_t nbPenalties;
//	if(is_calibration_by_column)
//	{
//		nbPenalties = (N - sparse_step - jjColIndex) / sparse_step + 1;
//	}
//	else 
//	{
//		size_t nb_sparse_diagonal = N/sparse_step - 2; // number of sparse sub diagonal that need to compute the constraint, the subdiagonal having only one element is not counted
//		// nb penalties = sum from 1 to nb_sparse_diagonal : each sub diagonal has (nb_sparse_diagonal_element - 1) penalties
//		nbPenalties = nb_sparse_diagonal*(nb_sparse_diagonal+1)/2;
//	}
//
//	QuantLib::Array constraint_values(nbPenalties,0.);
//
//	// compute penalties only by column
//	if(is_calibration_by_column && is_not_first_sparse_column)
//	{
//		for(size_t iPenalty = 0; iPenalty < nbPenalties ; ++iPenalty )
//		{
//			size_t jColIndex = jjColIndex - sparse_step;
//			size_t iRowRight = jjColIndex + iPenalty*sparse_step;
//			size_t iRowLeft  = iRowRight - sparse_step;
//			constraint_values[iPenalty]= abs(pHGVolatilityParamBuffer_->g(iRowRight,jjColIndex) - pHGVolatilityParamBuffer_->g(iRowLeft,jColIndex));		
//		}	
//	}
//
//	//compute penalties for global sparse indices, fill penalties columns by columns
//	if(!is_calibration_by_column)
//	{
//		size_t iPenalty = 0 ;	
//		for(size_t jjCol = 2*sparse_step;jjCol<N;jjCol+=sparse_step)
//		{
//			size_t nb_constraint_on_column = (N-jjCol) /sparse_step ;//normalize coeff constraint by number of constraints
//
//			size_t jCol = jjCol - sparse_step ;
//			
//			for(size_t iRowRight = jjCol;iRowRight<N;iRowRight+=sparse_step)
//			{
//				size_t iRowLeft  = iRowRight - sparse_step;
//				constraint_values[iPenalty]= abs(pHGVolatilityParamBuffer_->g(iRowRight,jjCol) - pHGVolatilityParamBuffer_->g(iRowLeft,jCol));
//				// YY don't do this // constraint_values[iPenalty]/=nb_constraint_on_column;//normalize constraint by number of constraints
//				++iPenalty;
//			}
//		}	
//		assert(iPenalty == nbPenalties);// after fully filling penalty, check if the number of penalties is correct
//	}
//
//	//// compute penalties for global sparse indices, fill penalties sub diagonal by sub diagonal
//	//if(!is_calibration_by_column)
//	//{
//	//	size_t iPenalty = 0 ;
//	//	size_t nbSubDiagonal = N/sparse_step -2;
//	//	for(size_t i_sparse_diagonal=0;i_sparse_diagonal<nbSubDiagonal;++i_sparse_diagonal)
//	//	{
//	//		size_t nb_constraint_on_diag = nbSubDiagonal - i_sparse_diagonal ;//normalize coeff constraint by number of constraints
//	//		size_t max_jCol = N - (i_sparse_diagonal +2 )*sparse_step;
//	//		for(size_t jCol=sparse_step; jCol<=max_jCol ;jCol+=sparse_step )
//	//		{
//	//			size_t iRow=jCol + sparse_step * i_sparse_diagonal;
//	//			
//	//			size_t jjCol = jCol+sparse_step;
//	//			size_t iiRow = iRow+sparse_step;
//	//			constraint_values[iPenalty]= pHGVolatilityParamBuffer_->g(iiRow,jjCol) - pHGVolatilityParamBuffer_->g(iRow,jCol);	
//	//			//constraint_values[iPenalty]/=nb_constraint_on_diag;//normalize coeff constraint by number of constraints
//	//			++iPenalty;
//	//		}
//	//	}
//	//	assert(iPenalty == nbPenalties);// after fully filling penalty, check if the number of penalties is correct
//	//}
//
//
//	return constraint_values;
//}
//
//void CalibrationVolCostFunction::get_CalibrationType(size_t& sparse_step, size_t& col_index, const VectorOfMatrixIndices& indices) const
//{
//	//Todo
//	size_t N = pHGVolatilityParamBuffer_->get_horizon();
//
//	size_t bottom_margin = N - indices.back().first;
//
//	sparse_step = bottom_margin;
//	assert(N%sparse_step == 0);
//
//	if(indices.size() == 1)// indices is only one element
//	{
//		//! case where the whole vol matrix has only one sparse elements  --> calibration global but the vol matrix is small
//		if(indices.back().first == sparse_step ) col_index=0;
//
//		//! case of calibration by column and this is the last column that has only one element
//		if(indices.back().first > sparse_step  ) col_index= indices[0].second;	
//	}
//	else // indices has many elements
//	{
//		
//		if(indices.back().first ==indices.back().second )
//			col_index=0;// the case where indices are multiple element, but the last element is on the diagonal, this is a global calibration
//		else
//			col_index=indices.back().second ; // this is a by-maturity calibraiton
//	}
//}
//
//QuantLib::Real CalibrationVolCostFunction::regularisation(const QuantLib::Array& x, QuantLib::Real c1, QuantLib::Real c2, QuantLib::Real c3, QuantLib::Real c4) const
//{
//	//-- Map array x to matrix;
//	matrix_ H = map_ArrayToMatrix(x);
//
//	//-- Compute derivatives
//
//	double sum_derivatives_maturity = 0.;
//	for (size_t i = 0; i < H.size(); ++i)
//	{
//		for (size_t j = 0; j < H[i].size()-1; ++j)
//			sum_derivatives_maturity += weights_maturity_[i][j] * abs(H[i][j+1]-H[i][j]);
//	}
//	sum_derivatives_maturity *= (c1/sum_all_weights_regularisation(weights_maturity_));
//
//	double sum_derivatives_tenor = 0;
//	for (size_t i = 0; i < H.size()-1; ++i)
//	{
//		for (size_t j = 0; j < H[i].size(); ++j)
//			sum_derivatives_tenor += weights_tenor_[i][j] * abs(H[i+1][j]-H[i][j]);
//	}
//	sum_derivatives_tenor *= (c2/sum_all_weights_regularisation(weights_tenor_));
//
//	double delta = 1.;
//	double sum_derivatives2_maturity = 0;
//	for (size_t i = 0; i < H.size(); ++i)
//	{
//		for (size_t j = 0; j < H[i].size()-2; ++j)
//			sum_derivatives2_maturity += weights_maturity_2_[i][j] * abs(H[i][j+2]-2*H[i][j+1]+H[i][j])/delta;
//	}
//	sum_derivatives2_maturity *= (c3/sum_all_weights_regularisation(weights_maturity_2_));
//
//	double sum_derivatives2_tenor = 0;
//	for (size_t i = 0; i < H.size()-2; ++i)
//	{
//		for (size_t j = 0; j < H[i].size(); ++j)
//			sum_derivatives2_tenor += weights_tenor_2_[i][j] * abs(H[i+2][j]-2*H[i+1][j]+H[i][j])/delta;
//	}
//	sum_derivatives2_tenor *= (c4/sum_all_weights_regularisation(weights_tenor_2_));
//
//	QuantLib::Real regulator = sum_derivatives_maturity + sum_derivatives_tenor + sum_derivatives2_maturity + sum_derivatives2_tenor;
//
//	return regulator;
//}
//
//
//QuantLib::Real CalibrationVolCostFunction::sum_all_weights_regularisation(const matrix_& weights) const 
//{
//	QuantLib::Real sum = 0.;
//	for (size_t i = 0; i < weights.size(); ++i)
//	{
//		for (size_t j = 0; j < weights[i].size(); ++j)
//			sum += weights[i][j];
//	}
//	return sum;
//}
//
////Real RebonatoVolatilityCostFunction::sum_all_derivatives_regularisation(const matrix_& weights, const matrix_& derivatives)
////{
////	Real sum = 0.;
////	for (size_t i = 0; i < weights.size(); ++i)
////	{
////		for (size_t j = 0; j < weights[0].size(); ++j)
////			sum += weights[i][j] * derivatives[i][j];
////	}
////	return sum;
////}
//
//
//CalibrationVolCostFunction::matrix_ CalibrationVolCostFunction::map_ArrayToMatrix(const QuantLib::Array& x) const
//{
//	matrix_ H;
//	size_t x_size = x.size();
//	size_t nb_alloc = 1;
//	size_t current_index = 0;
//
//	while (current_index < x_size)
//	{
//		std::vector<double> H_row;
//		for (size_t i = 0; i < nb_alloc; ++i)
//		{
//			H_row.push_back(x[current_index]);
//			current_index++;
//		}
//		H.push_back(H_row);
//		nb_alloc++;
//	}
//
//	for (size_t i = 0; i < H.size(); ++i)
//	{
//		size_t current_index = H[i].size();
//		for (size_t j = current_index; j < H.size(); ++j)
//			H[i].push_back(0.);
//	}
//
//	return H;
//}
//
//QuantLib::Array CalibrationVolCostFunction::map_MatrixtoArray(const matrix_& mat) const
//{
//	std::vector<double> tmpRes;
//
//	for (auto vect : mat)
//	{
//		for (auto val : vect)
//		{
//			if (val == 0.)
//				continue;
//
//			tmpRes.push_back(val);
//		}
//	}
//
//	size_t nbVol = tmpRes.size();
//	QuantLib::Array res(nbVol);
//	for (size_t i = 0; i < nbVol; ++i)
//		res[i] = tmpRes[i];
//
//	return res;
//}
