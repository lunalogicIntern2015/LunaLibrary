#include <LMM/Helper/GenericPath.h>
#include <LMM/Helper/Printer.h>
#include <LMM/Mc/McLmm.h>

using namespace std;

namespace MCSchemeType
{
	std::string mcSchemeType2String(MCSchemeType mcSchemeType)
	{
		switch (mcSchemeType)
		{
		case EULER:
			return "EULER";
		case PC:
			return "PC";
		case IPC:
			return "IPC";
		default:
			throw("Error: not valid MCSchemeType");
		}
	}
}

McLmm::McLmm(Lmm_PTR lmm,
			 const std::vector<double>&         liborsInitValue,
			 RNGenerator_PTR                    rnGenerator,
			 MCSchemeType::MCSchemeType			mcSchemeType)
: //Lmm(dispersion,shifts)
  lmm_(lmm)
, nbFactor_(lmm->get_dispersionRef().getNbFactors())
, horizon_(lmm->get_dispersionRef().get_horizon())
, B_(lmm->get_dispersionRef().get_Correlation_PTR()->get_reducedCorrelMatrixB())
, liborsInitValue_(liborsInitValue)
, numeraires_ (horizon_+1, 0.0)
, liborMatrix_(horizon_+1, horizon_+1)
, rnGenerator_(rnGenerator)
, mcSchemeType_(mcSchemeType)
{
	assert(lmm->get_dispersionRef().get_LMMTenorStructure_PTR()->get_tenorDate()[0] == 0.0); 

	size_t horizon = lmm->get_dispersionRef().get_horizon();
	//assert(shifts.size() == horizon+1);
	assert(liborsInitValue.size() == horizon+1);
	
	initLiborMatrix(liborsInitValue);  // Libor at time 0			
}

//! \int_{T_{k-1}}^{T_k} sigma_i(t)  dW_t,   indexTime = k
//! we suppose the correlation is constant: \int_{T_{k-1}}^{T_k} ||sigma_i(t)||^2 dt * <B,G>, where  G is nbFactor indepedent Gaussian.
QuantLib::Real McLmm::computeIntSto(size_t indexTime, 
						  size_t indexLibor, 
						  const std::vector<double>& G)  const // nbFactor D, independent Gaussian
{
	QuantLib::Real varTmp = sqrt(lmm_->get_covarianceTensor(indexTime, indexLibor, indexLibor));  

	QuantLib::Real gaussian = 0;  // 1D correlated gaussian
	for (size_t k = 0; k < nbFactor_; ++k) 
		gaussian += B_[indexLibor][k] * G[k];

	double res = gaussian*varTmp;
	return res;
}


void McLmm::initLiborMatrix(const std::vector<double>& libors_init)
{
	//liborMatrix_ = matrix_(horizon_+1,std::vector<double>(horizon_+2,0.));
	for (size_t i = 0; i <= horizon_; ++i)
		liborMatrix_(i,0) = libors_init[i];
}




//void McLmm::set_shifts(const std::vector<double>& shifts) 
//{
//	assert(shifts.size() == horizon_+1);
//	shifts_ = shifts;
//}


void McLmm::print(const std::string& filename) const
{
	//std::string fileName = "McLmm.csv";
	std::string path = LMMPATH::get_output_path() + filename;

	std::vector<PrintElement_PTR> elements_print;

	//! CorrelationMatrix
	//PrintElement_PTR to_print      = PrintElement_PTR(new MatrixPrintElement<matrix>("CorrelationMatrix", 
	//													  QLMatrix2BoostMatrix(dispersion_.get_Correlation_PTR()->get_reducedCorrelMatrixApprox())));
	//elements_print.push_back(to_print);

	//! covarianceTensor
	for(size_t indexT = 1; indexT<=horizon_; ++indexT)
	{
		PrintElement_PTR to_print      = PrintElement_PTR(new MatrixPrintElement<matrix>("covarianceTensor"+Number2String<size_t>(indexT), lmm_->get_covarianceTensor()[indexT]));
		elements_print.push_back(to_print);
	}
	
	Printer printer(path, elements_print);
	printer.print();
}





////! YY TODO: Need to rewrite copy constructor, not well done !
//McLmm::McLmm (const McLmm & myLmm) 
//{
//	throw ("Error: the copy constructor of McLmm is not implemented!");
//	//this->horizon_ = myLmm.horizon_;
//	//this->deltaT_ = myLmm.deltaT_;
//	//this->tenorDates = myLmm.tenorDates;
//	////this->myVol = myLmm.myVol;
//	//vol_ = myLmm.vol_;
//	//shifts_= myLmm.shifts_;
//	//liborMatrix_ = myLmm.liborMatrix_;
//	//numeraires_ = myLmm.numeraires_;
//	//cumulated_drifts_ = myLmm.cumulated_drifts_;
//	//cumulated_squared_drifts_ = myLmm.cumulated_squared_drifts_;
//	//generator_ = myLmm.generator_;
//}


//McLmm::~McLmm() 
//{
//	//delete [] myVol;
//}


//-- Rk : small differences between simulated libors frm that method and the TestPricing method
//-- TODO : fix that
// the covariance matrix contains vol covariances_t1_t2 computed between dates t1 and t2
//numeric::ublas::vector<Real> McLmm::drifts(
//	                                       const numeric::ublas::vector<Real> & libor_temp, 
//	                                       size_t index_T,
//	                                       Time t1, 
//										   Time t2,
//										   matrix<Real> &covariances_t1_t2
//										   ) 
//{
//	numeric::ublas::vector<Real> drift_temp(horizon_+1,0.);
//	numeric::ublas::vector<Real> res_drifts(horizon_+1,0.);
//
//	Real tmp = 0.0;
//	double tmp2 = 0.0;
//	
//	// We start at horizon_-1 because Libor_horizon_ has no drift
//	// We stop at index_T because Libor_i(T_{index_T}) has no meaning if i < index_T
//	for (size_t index_libor = horizon_-1; index_libor >= index_T; --index_libor) { 
//		cout << libor_temp(index_libor+1) << " ";
//		tmp = libor_temp(index_libor+1)*deltaT_(index_libor+1);
//		drift_temp(index_libor+1) = tmp/(1+tmp);
//
//		for (size_t k = index_libor+1; k < horizon_+1; ++k) { // Inverse loop
//			res_drifts(index_libor) -= drift_temp(k) * covariances_t1_t2(index_libor,k);
//		}
//	}
//	cout << endl;
//		
//	return res_drifts;
//}
//
//Real McLmm::drifts_ipc(size_t index_libor,
//	               const numeric::ublas::vector<Real> & libor_t1, 
//	               const numeric::ublas::vector<Real> & libor_temp, 
//				   size_t index_T,
//				   Time t1, 
//				   Time t2,
//				   matrix<Real> &covariances_t1_t2
//				   ) 
//{
//	numeric::ublas::vector<Real> drift_temp(horizon_+1,0.);
//	Real res_drifts = 0.;
//
//	Real tmp = 0.0;
//	double tmp2 = 0.0;
//
//	tmp = libor_t1(index_libor)*deltaT_(index_libor);  
//	tmp2 = libor_temp(index_libor)*deltaT_(index_libor); 
//	drift_temp(index_libor) = 0.5*(tmp/(1+tmp) + tmp2/(1+tmp2));
//	for (size_t k = index_libor; k < horizon_+1; ++k) {
//		res_drifts -= drift_temp(k) * covariances_t1_t2(index_libor-1,k);
//	}
//	
//	return res_drifts;
//}
//
//numeric::ublas::vector<Real> McLmm::spot_drifts(
//	                                       const numeric::ublas::vector<Real> & libor_temp, 
//	                                       size_t index_T,
//	                                       Time t1, 
//										   Time t2,
//										   matrix<Real> &covariances_t1_t2
//										   ) 
//{
//	numeric::ublas::vector<Real> drift_temp(horizon_+1,0.);
//	numeric::ublas::vector<Real> res_drifts(horizon_+1,0.);
//
//	Real tmp = 0.0;
//	double tmp2 = 0.0;
//
//	for (size_t index_libor = index_T; index_libor < horizon_; ++index_libor) {
//
//		tmp = libor_temp(index_libor+1)*deltaT_(index_libor+1);
//
//		drift_temp(index_libor+1) = tmp/(1+tmp);
//
//		for (size_t k = index_T; k <= index_libor; ++k) {
//			std::cout << "index libor : " << index_libor << " -- current index : " << k << " -- Covariance : " << covariances_t1_t2(index_libor,k) << std::endl;
//			res_drifts(index_libor) -= drift_temp(k) * covariances_t1_t2(index_libor,k);
//		}
//
//		res_drifts(index_libor) *= libor_temp(index_libor);
//	}
//
//	return res_drifts;
//}
//
//
//numeric::ublas::vector<Real> McLmm::diffusions(
//	                                           const numeric::ublas::vector<Real> & var_temp, 
//						                       const numeric::ublas::vector<Real> &G,
//											   Time index_T,
//											   Time t1, 
//											   Time t2) 
//{
//	numeric::ublas::vector<Real> diffusionVect(horizon_+1,0.);
//
//	for (size_t index_libor = horizon_; index_libor >= index_T; --index_libor) {
//
//		for (size_t nb_factors = 0; nb_factors < myVol.getNbFactors(); ++nb_factors) {
//			Real B = myVol.getterB()(index_libor,nb_factors);
//			diffusionVect(index_libor) += B*G(nb_factors);
//		}
//
//		diffusionVect(index_libor) *= sqrt(var_temp(index_libor));
//	}
//
//	return diffusionVect;
//}


// Rk : lmm matrix is supposed to be already filled with Libors at time t=0
//void McLmm::TerminalLmmModel_Euler() 
//{
//
//	size_t nbFactors = myVol.getNbFactors();
//	numeric::ublas::vector<Real> drift_temp(horizon_+1,0.);
//
//	//-- Create the Gaussian vector
//	std::vector<Real> gaussian_tmp(nbFactors,0.);
//	
//	// Compute each Libor rate's value at every tenor date (except at time 0)
//	for (size_t j = 1; j < horizon_+1; ++j) { // time
//
//		//-- Fill the Gaussian vector
//		generator_.pseudoRandomNumberGenerator(gaussian_tmp);
//
//		for (size_t i = horizon_; i >= j; --i) { // Backward simulations
//		 
//			Real diffusion_factor = exp(-0.5 * covariance_tensor[j-1][i][i])
//							* exp(myVol.computeIntSto(i,j,covariance_tensor,gaussian_tmp));
//
//			Real drift_factor = 1.0;
//			
//			// Drift computation if i < horizon_
//			if (i < horizon_) {
//
//				double tmp = liborMatrix_[i+1][j-1]*deltaT_[i+1];
//				double tmp_shifted = (liborMatrix_[i+1][j-1]+shifts_[i+1])*deltaT_[i+1];
//				drift_temp(i+1) = tmp_shifted/(1+tmp);
//
//				Real drift = 0.;
//				for (size_t k = i+1; k < horizon_+1; ++k) 
//					drift -= drift_temp(k) * covariance_tensor[j-1][i][k];
//
//				drift_factor *= exp(drift);
//			}
//
//			//-- Drift accumulation for MC variance computing
//			cumulated_drifts_[i][j] += drift_factor;
//			cumulated_squared_drifts_[i][j] += drift_factor*drift_factor;
//
//			liborMatrix_[i][j] = (liborMatrix_[i][j-1] + shifts_[i] )*drift_factor*diffusion_factor - shifts_[i];
//		}
//	}
//}	
//
//// YY  predictor-corrector
//void McLmm::TerminalLmmModel_Pc() 
//{  
//	std::vector<Real> gaussian_tmp(myVol.getNbFactors(),0.);
//
//	//matrix<Rate> &lmm_pc(lmm); // YY TODO: make if more efficiently latter...
//	std::vector<double> lmm_pc(liborMatrix_.size());
//
//	// Compute each libor rate's value at every tenor date (except at time 0)
//	for (size_t j = 1; j < horizon_+2; ++j) // time
//	{ 
//		numeric::ublas::vector<Real> drift_temp(horizon_+1,0.);
//
//		//-- Fill the gaussian vector
//		generator_.pseudoRandomNumberGenerator(gaussian_tmp);
//
//		//! ---- predictor-part
//		for (size_t i = horizon_; i >= j; --i) // Backward simulations
//		{ 
//			Real diffusion_factor_pc = exp(-0.5 * covariance_tensor[j-1][i][i])
//									* exp(myVol.computeIntSto(i,j,covariance_tensor,gaussian_tmp));
//			Real drift_factor_pc = 1.;
//
//			// Drift computation if i < horizon_
//			if (i < horizon_) 
//			{
//				double tmp = liborMatrix_[i+1][j-1]*deltaT_[i+1];
//				double tmp_shifted = (liborMatrix_[i+1][j-1]+shifts_[i+1])*deltaT_[i+1];
//				drift_temp(i+1) = tmp_shifted/(1+tmp);
//
//				Real drift = 0.0;
//				for (size_t k = i+1; k < horizon_+1; ++k) 
//					drift -= drift_temp(k) * covariance_tensor[j-1][i][k];
//					
//				drift_factor_pc *= exp(drift);
//			}
//
//			lmm_pc[i] = (liborMatrix_[i][j-1]+shifts_[i])*drift_factor_pc*diffusion_factor_pc - shifts_[i];
//		}
//
//		//! ---- corrector-part
//		for (size_t i = horizon_; i >= j; --i) // Backward simulations
//		{ 
//			Real diffusion_factor = exp(-0.5 * covariance_tensor[j-1][i][i])
//									* exp(myVol.computeIntSto(i,j,covariance_tensor,gaussian_tmp));
//			Real drift_factor = 1.;
//
//			// Drift computation if i < horizon_
//			if (i < horizon_) 
//			{
//				double tmp  = lmm_pc[i+1]*deltaT_[i+1];   // t + Delta_t
//				double tmp2 = liborMatrix_[i+1][j-1]*deltaT_[i+1];  // t
//
//				double tmp_shifted = (lmm_pc[i+1]+shifts_[i+1])*deltaT_[i+1];
//				double tmp2_shifted = (liborMatrix_[i+1][j-1]+shifts_[i+1])*deltaT_[i+1];
//
//				drift_temp(i+1) = 0.5*( tmp_shifted/(1+tmp) + tmp2_shifted/(1+tmp2) ); 
//
//				Real drift = 0.0;
//				for (size_t k = i+1; k < horizon_+1; ++k) 
//					drift -= drift_temp(k) * covariance_tensor[j-1][i][k];
//					
//				drift_factor *= exp(drift);
//			}
//
//			liborMatrix_[i][j] = (liborMatrix_[i][j-1]+shifts_[i])*drift_factor*diffusion_factor - shifts_[i];
//		}
//	}
//	//-- TODO : Work directly with the libor matrix from attributes
//	//liborMatrix_ = lmm;
//}	
//
////YY 
//void McLmm::SpotLmmModel_euler() 
//{
//
//	std::vector<Real> gaussian_tmp(myVol.getNbFactors(),0.);
//
//	// Compute each libor value at time t
//	for (size_t j = 1; j < horizon_+2; ++j) {
//
//		numeric::ublas::vector<Real> drift_temp(horizon_+1,0.);
//
//		//-- Fill the gaussian vector
//		generator_.pseudoRandomNumberGenerator(gaussian_tmp);
//
//		for (size_t i = j; i < horizon_+1; ++i) { // Forward simulations
//			
//			Real diffusion_factor = exp(-0.5 * covariance_tensor[j-1][i][i])
//									* exp(myVol.computeIntSto(i,j,covariance_tensor,gaussian_tmp));
//			Real drift_factor = 1.;
//
//			if (i < horizon_) 
//			{
//				double tmp = liborMatrix_[i][j-1]*deltaT_[i];
//				drift_temp(i) = tmp/(1+tmp);
//				Real drift = 0;
//				for (size_t k = j; k <= i; ++k) 
//					drift += drift_temp(k) * covariance_tensor[j-1][i][k];
//					
//				drift_factor *= exp(drift);
//			}
//
//			liborMatrix_[i][j] = (liborMatrix_[i][j-1]+shifts_[i])*drift_factor*diffusion_factor - shifts_[i];
//		}
//	}
//	//-- TODO : Work directly with the libor matrix from attributes
//	//liborMatrix_ = lmm;
//}
//
//
//void McLmm::SpotLmmModel_pc() 
//{
//	std::vector<Real> gaussian_tmp(myVol.getNbFactors(),0.);
//	std::vector<double> lmm_pc(liborMatrix_.size());
//
//	for(size_t index_time = 1; index_time < horizon_+2; ++index_time) // Time loop
//	{
//		numeric::ublas::vector<Real> drift_temp(horizon_+1,0.);
//
//		//-- Fill the gaussian vector
//		generator_.pseudoRandomNumberGenerator(gaussian_tmp);
//
//		// Prediction (Except for the last libor, as it has no drift)
//		for (size_t index_libor = index_time; index_libor < horizon_; ++index_libor) // loop on libors
//		{
//			Real diffusion_factor_pc = exp(-0.5 * covariance_tensor[index_time-1][index_libor][index_libor])
//									* exp(myVol.computeIntSto(index_libor,index_time,covariance_tensor,gaussian_tmp));
//
//			double tmp = liborMatrix_[index_libor][index_time-1] * deltaT_[index_libor];
//			double tmp_shifted = (liborMatrix_[index_libor][index_time-1] + shifts_[index_libor]) * deltaT_[index_libor];
//			drift_temp(index_libor) = tmp_shifted/(1+tmp); 
//
//			Real drift = 0;
//			for (size_t k = index_time; k <= index_libor; ++k) {
//				drift += drift_temp(k) * covariance_tensor[index_time-1][index_libor][k]; 
//			}
//				
//			//lmm_pc[index_libor] *= exp(drift);	
//			Real drift_factor_pc = exp(drift);	
//
//			lmm_pc[index_libor] = (liborMatrix_[index_libor][index_time-1]+shifts_[index_libor])
//									* drift_factor_pc * diffusion_factor_pc - shifts_[index_libor];
//		}
//
//		// Correction		
//		for (size_t index_libor = index_time; index_libor < horizon_+1; ++index_libor) // loop on libors
//		{
//			Real diffusion_factor = exp(-0.5 * covariance_tensor[index_time-1][index_libor][index_libor])
//									* exp(myVol.computeIntSto(index_libor,index_time,covariance_tensor,gaussian_tmp));
//
//			Real drift_factor = 1.;
//
//			// Compute corrected drifts, except for the last libor
//			if (index_libor < horizon_) 
//			{
//				double tmp = liborMatrix_[index_libor][index_time-1] * deltaT_[index_libor];
//				double tmp2 = lmm_pc[index_libor] * deltaT_[index_libor];
//
//				double tmp_shifted = (liborMatrix_[index_libor][index_time-1]+shifts_[index_libor]) * deltaT_[index_libor];
//				double tmp2_shifted = (lmm_pc[index_libor]+shifts_[index_libor]) * deltaT_[index_libor];
//
//				drift_temp(index_libor) = 0.5*(tmp_shifted/(1+tmp) + tmp2_shifted/(1+tmp2)); 
//
//				Real drift = 0;
//				for (size_t k = index_time; k <= index_libor; ++k) {
//					drift += drift_temp(k) * covariance_tensor[index_time-1][index_libor][k]; 
//				}
//		
//				drift_factor *= exp(drift);		
//			}
//
//			liborMatrix_[index_libor][index_time] = (liborMatrix_[index_libor][index_time-1]+shifts_[index_libor])
//											* drift_factor * diffusion_factor - shifts_[index_libor];
//		}
//	}
//}
//
//
//void McLmm::TerminalLmmSimulation(ApproxType approx) 
//{
//	switch(approx)
//	{
//	case euler:
//		TerminalLmmModel_Euler();
//		break;
//	case ipc:
//		TerminalLmmModel_Pc();
//		break;
//	default:
//		break;
//	}
//}
//
//void McLmm::SpotLmmSimulation(ApproxType approx) 
//{
//	switch(approx)
//	{
//	case euler:
//		SpotLmmModel_euler();
//		break;
//	case ipc:
//		SpotLmmModel_pc();
//		break;
//	default:
//		break;
//	}
//}


//Real McLmm::bondPrice_i_j(size_t i, size_t j, const std::vector<std::vector<double>>& lmm) { // P(T_i,T_j)
//Real McLmm::bondPrice_i_j(size_t i, size_t j) { // P(T_i,T_j)
//
//	Real res = 1.;
//	Real per;
//
//	for (size_t k = i; k < j; ++k) {
//		per = this->deltaT_[k];
//		res *= 1/(1 + per*liborMatrix_[k][i]);
//	}
//
//	return  res;
//}



////Getter & Setter
//size_t                             McLmm::get_horizon()     const { return this->horizon_; }
////DeterministicVol McLmm::getVolMod() const { return this->myVol; }
////InstantaneousVolatility*           McLmm::get_volMod() const { return vol_; }
//const std::vector<double>& McLmm::get_tenorDates()  const { return this->tenorDates; }
//const std::vector<double>& McLmm::get_deltaT()      const { return this->deltaT_; }
//const std::vector<double>&         McLmm::get_shifts()      const {return shifts_;}
//const matrix_&                     McLmm::get_liborMatrix() const {return liborMatrix_;}
//const std::vector<double>&         McLmm::getNumeraire()    const {return numeraires_;}
//double			         McLmm::get_numeraire(size_t index) const {return numeraires_[index];}

//Real McLmm::Spot_numeraire(size_t k, const std::vector<std::vector<double>>& lmm) // spot numeraire, B(T_0=0) = 1.0
//{
//	Real res = 1;
//	Time per = tenorDates[k] - tenorDates[k-1]; 
//	for (size_t j = 0; j < k; ++j) {
//		res *= (1 + per*lmm[j][j]); 
//	}
//	return res;
//}

//
////! YY static function.
//void McLmm::drift_variance(size_t nbSimulations,
//						   RNG_Type & generator,
//						   const std::vector<double>& libors_T0,
//						   const std::vector<matrix_>& covTensor)
//{
//	for (size_t nbSim = 0; nbSim < nbSimulations; ++nbSim)
//	{
//		//-- Libor simulation
//		//TerminalLmmModel_Euler();
//		simulateLMM(MCSchemeType::EULER);
//	}
//	
//	for (size_t i = 0; i < cumulated_drifts_.size(); ++i)
//	{
//		for (size_t j = 0; j < cumulated_drifts_[0].size(); ++j)
//		{
//			cumulated_drifts_[i][j] /= nbSimulations;
//			cumulated_drifts_[i][j] *= ((double)nbSimulations)/(nbSimulations-1)
//											*cumulated_drifts_[i][j];
//		}
//	}
//
//	for (size_t i = 0; i < cumulated_drifts_.size(); ++i)
//	{
//		for (size_t j = 0; j < cumulated_drifts_[0].size(); ++j)
//		{
//			cumulated_squared_drifts_[i][j] /= (nbSimulations-1);
//			cumulated_squared_drifts_[i][j] -= cumulated_drifts_[i][j];
//		}
//	}
//
//	cout << "//---- Drift variances" << endl;
//	for (size_t i = 0; i < cumulated_squared_drifts_.size(); ++i)
//	{
//		for (size_t j = 0; j < cumulated_squared_drifts_[0].size(); ++j)
//			cout << cumulated_squared_drifts_[i][j] << " ";
//
//		cout << endl;
//	}
//	cout << endl;
//	
//}
//
//
////! YY static function.
////Adrien: -- TODO : Finish it
//void McLmm::omega_i_variance(size_t omega_index,
//	                       size_t nbSimulations,
//						   std::vector<size_t> fixing_indices,
//		                   RNG_Type & generator,
//						   const std::vector<double>& libors_T0,
//						   const std::vector<matrix_>& covTensor)
//{
//	throw("Error, not implemented!");
//	//matrix<Rate> lmm(horizon_+1,horizon_+2,0.);
//	//for (size_t i = 0; i < horizon_+1; ++i)
//	//	lmm(i,0) = libors_T0[i];
//	//for (size_t index_time = 0; index_time <= omega_index+1; ++index_time)
//	//{   
//	//	//-- Compute annuity with fixing indices
//	//	double annuity = 0.;
//	//	for each (size_t fixing_index in fixing_indices)
//	//		annuity += bondPrice_i_j(index_time,fixing_index+1,lmm);
//	//	// Compute P(T_indextime,T_omegaindex)
//	//}
//}
//
//
//


//// Printer
//void McLmm::printLiborMatrix() 
//{
//	cout << "//--- Libor matrix:" << endl << endl;
//	for (size_t i = 0; i < liborMatrix_.size1(); ++i)
//	{
//		for (size_t j = 0; j < liborMatrix_.size2(); ++j)
//			cout << liborMatrix_(i,j) << " ";
//
//		cout << endl;
//	}
//	cout << endl;
//}
//void McLmm::printNumeraires()
//{
//	cout << "//--- Numeraires:" << endl;
//	for (size_t i = 0; i < numeraires_.size(); ++i)
//	{
//		cout << numeraires_[i] << " ";
//	}
//	cout << endl;
//}


//
//
////--------------------------------------------------------------------------//
////                                 TESTS                                    //
////--------------------------------------------------------------------------//
//
////-- USE MC GENERATOR CLASS INSTEAD
//void McLmm::simulationTest_Terminal_Euler()
//{
//	////-- Compute a gaussian matrix 
//	//matrix<Real> gaussianMatrix(myVol.getNbFactors(),horizon_+1);
//	//numeric::ublas::vector<Real> G(myVol.getNbFactors(),0.);
//	//boost::variate_generator<boost::mt19937, boost::normal_distribution<> >
//	//		generator(boost::mt19937(0.5),
//	//		boost::normal_distribution<>());
//	//for (size_t k = 0; k < horizon_+1; ++k) {
//	//	gen_normal(generator,G);
//	//	for (size_t i = 0; i < myVol.getNbFactors(); ++i)
//	//		gaussianMatrix(i,k) = G(i);
//	//}
//	////-- Compute covariances ( int_Tk^Tk+1 {sigma_i * sigma_j * rho_ij} )
//	//std::vector<std::vector<std::vector<double>>> covTensor(horizon_+1);
//	//covTensor = myVol.computeAllCovariances(horizon_,tenorDates);
//	//std::vector<std::vector<double>> lmm(horizon_+1);
//	//for (size_t i = 0; i < horizon_+1; ++i)
//	//{
//	//	std::vector<double> lmm_column(horizon_+2,0.);
//	//	lmm.push_back(lmm_column);
//	//}
//	//double r = 0.02;
//	//for (size_t i = 0; i < horizon_+1; ++i)
//	//	lmm[i][0] = libor_0(i,r) + shifts_[i];
//	
//
//	//-- Perform one simulation
//	LmmSimulation(euler);
//
//	matrix_ lmm = liborMatrix_;
//	cout << "----> TERMINAL EULER : " << endl << endl;
//	for(size_t i = 0; i < lmm.size(); ++i) 
//	{
//		for (size_t j = 0; j < lmm[0].size(); ++j)
//			cout << lmm[i][j] << " ";
//
//		cout << endl << endl;
//	}
//}
//
//void McLmm::simulationTest_Terminal_Pc()
//{
//	////-- Compute a gaussian matrix 
//	//matrix<Real> gaussianMatrix(myVol.getNbFactors(),horizon_+1);
//	//numeric::ublas::vector<Real> G(myVol.getNbFactors(),0.);
//	//boost::variate_generator<boost::mt19937, boost::normal_distribution<> >
//	//		generator(boost::mt19937(0.5),
//	//		boost::normal_distribution<>());
//	//for (size_t k = 0; k < horizon_+1; ++k) {
//	//	gen_normal(generator,G);
//	//	for (size_t i = 0; i < myVol.getNbFactors(); ++i)
//	//		gaussianMatrix(i,k) = G(i);
//	//}
//	////-- Compute covariances ( int_Tk^Tk+1 {sigma_i * sigma_j * rho_ij} )
//	//std::vector<std::vector<std::vector<double>>> covTensor(horizon_+1);
//	//covTensor = myVol.computeAllCovariances(horizon_,tenorDates);
//	//std::vector<std::vector<double>> lmm(horizon_+1);
//	//for (size_t i = 0; i < horizon_+1; ++i)
//	//{
//	//	std::vector<double> lmm_column(horizon_+2,0.);
//	//	lmm.push_back(lmm_column);
//	//}
//	//double r = 0.02;
//	//for (size_t i = 0; i < horizon_+1; ++i)
//	//	lmm[i][0] = libor_0(i,r)+shifts_[i];
//
//	//-- Perform one simulation
//	LmmSimulation(ipc);
//	matrix_ lmm = liborMatrix_;
//	cout << "----> TERMINAL PREDICTOR-CORRECTOR : " << endl << endl;
//	for(size_t i = 0; i < lmm.size(); ++i) 
//	{
//		for (size_t j = 0; j < lmm[0].size(); ++j)
//			cout << lmm[i][j] << " ";
//
//		cout << endl << endl;
//	}
//}
//
//void McLmm::simulationTest_Spot_Euler()
//{
//	//-- Compute a gaussian matrix 
//	//matrix<Real> gaussianMatrix(myVol.getNbFactors(),horizon_+1);
//	//numeric::ublas::vector<Real> G(myVol.getNbFactors(),0.);
//	//boost::variate_generator<boost::mt19937, boost::normal_distribution<> >
//	//		generator(boost::mt19937(0.5),
//	//		boost::normal_distribution<>());
//	//for (size_t k = 0; k < horizon_+1; ++k) {
//	//	gen_normal(generator,G);
//	//	for (size_t i = 0; i < myVol.getNbFactors(); ++i)
//	//		gaussianMatrix(i,k) = G(i);
//	//}
//	////-- Compute covariances ( int_Tk^Tk+1 {sigma_i * sigma_j * rho_ij} )
//	//std::vector<std::vector<std::vector<double>>> covTensor(horizon_+1);
//	//covTensor = myVol.computeAllCovariances(horizon_,tenorDates);
//	//std::vector<std::vector<double>> lmm(horizon_+1);
//	//for (size_t i = 0; i < horizon_+1; ++i)
//	//{
//	//	std::vector<double> lmm_column(horizon_+2,0.);
//	//	lmm.push_back(lmm_column);
//	//}
//	//double r = 0.02;
//	//for (size_t i = 0; i < horizon_+1; ++i)
//	//	lmm[i][0] = libor_0(i,r)+shifts_[i];
//
//	//-- Perform one simulation
//	LmmSimulation(euler);
//	matrix_ lmm = liborMatrix_;
//	cout << "----> SPOT EULER : " << endl << endl;
//	for(size_t i = 0; i < lmm.size(); ++i) 
//	{
//		for (size_t j = 0; j < lmm[0].size(); ++j)
//			cout << lmm[i][j] << " ";
//
//		cout << endl << endl;
//	}
//}
//
//void McLmm::simulationTest_Spot_Pc()
//{
//	////-- Compute a gaussian matrix 
//	//matrix<Real> gaussianMatrix(myVol.getNbFactors(),horizon_+1);
//	//numeric::ublas::vector<Real> G(myVol.getNbFactors(),0.);
//	//boost::variate_generator<boost::mt19937, boost::normal_distribution<> >
//	//		generator(boost::mt19937(0.5),
//	//		boost::normal_distribution<>());
//	//for (size_t k = 0; k < horizon_+1; ++k) {
//	//	gen_normal(generator,G);
//	//	for (size_t i = 0; i < myVol.getNbFactors(); ++i)
//	//		gaussianMatrix(i,k) = G(i);
//	//}
//	////-- Compute covariances ( int_Tk^Tk+1 {sigma_i * sigma_j * rho_ij} )
//	//std::vector<std::vector<std::vector<double>>> covTensor(horizon_+1);
//	//covTensor = myVol.computeAllCovariances(horizon_,tenorDates);
//	//std::vector<std::vector<double>> lmm(horizon_+1);
//	//for (size_t i = 0; i < horizon_+1; ++i)
//	//{
//	//	std::vector<double> lmm_column(horizon_+2,0.);
//	//	lmm.push_back(lmm_column);
//	//}
//	//
//	//double r = 0.02;
//	//for (size_t i = 0; i < horizon_+1; ++i)
//	//	lmm[i][0] = libor_0(i,r) + shifts_[i];
//
//	//-- Perform one simulation
//	LmmSimulation(ipc);
//	matrix_ lmm = liborMatrix_;
//	cout << "----> SPOT PREDICTOR-CORRECTOR : " << endl << endl;
//	for(size_t i = 0; i < lmm.size(); ++i) 
//	{
//		for (size_t j = 0; j < lmm[0].size(); ++j)
//			cout << lmm[i][j] << " ";
//
//		cout << endl << endl;
//	}
//}


