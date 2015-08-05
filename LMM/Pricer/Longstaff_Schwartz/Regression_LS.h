#pragma once
#include <vector>
#include <time.h>

#include <boost/shared_ptr.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp> 
#include <boost/numeric/ublas/vector_proxy.hpp> 
#include <boost/numeric/ublas/matrix.hpp> 
#include <boost/numeric/ublas/triangular.hpp> 
#include <boost/numeric/ublas/lu.hpp> 

#include <LMM/Helper/Name.h>

#include <Instrument/Rate/Rate1.h>
#include <LMM/Pricer/Longstaff_Schwartz/Basis.h>
#include <LMM/Pricer/Longstaff_Schwartz/McLmm_LS.h>
#include <LMM/Pricer/Longstaff_Schwartz/Basis_Evaluator.h>


namespace ublas = boost::numeric::ublas; 

typedef ublas::matrix<double> Rmatrix;
typedef ublas::vector<double> Rvector;

namespace LS
{

// represent conditional expectation  or VC
class RegressionRepresentation  // regression(path+i) = sum_k coeff_k * basis_k(path_i)
{

	//! basis + evaluator (to be evaluated on eath path)
	std::vector<Basis_CONSTPTR> basis_;                      // size = nbBasis, evaluated on each simulation path 
	std::vector<Basis_Evaluator_CONSTPTR> basisEvaluators_;  // size = nbBasis 

	mutable std::vector<double> basis_val_buffer_; // basis value evaluated on 1 simulation path, it's a buffer

	//! result: coefficient 
	std::vector<double> regressionCoeffs_;

public:
	RegressionRepresentation(	const std::vector<Basis_CONSTPTR>& basis, 
								const std::vector<Basis_Evaluator_CONSTPTR>& basisEvaluators);

	virtual ~RegressionRepresentation(){}

	//getter  !!! not good to copy vector
	const std::vector<double>& getBasis_val_buffer()const{return basis_val_buffer_;}
	const std::vector<double>&  getRegressionCoef()const{return regressionCoeffs_;}

	const std::vector<Basis_CONSTPTR>& getBasis()const{return basis_;}

	//accessor 
	std::vector<double>& getRegressionCoeffs(){return regressionCoeffs_;}

	// evaluate basis values on one simulation path
	void evaluate_basis(const McLmm_LS::LMMSimulationResult& lmmSimulationResult)const;
	void evaluate_basis(const matrix& liborMatrix, const std::vector<double>& numeraire)const; 

	//!! continuous value (cv)
	double evaluate_representation(const McLmm_LS::LMMSimulationResult& lmmSimulationResult)const; // conditional expectation 

	void write_to_stream(std::ostream& out)const;

};

//! for a given date: T_i  ->   T-{i-1}
//! T_i: Z_i regression on T_{i-1}: basis 
class Regression  // one Step regression
{
	//! T_{i-1}: basis  (to be evaluated on each path)
	//! regression result can be saved as coeff in this class 
	LS::RegressionRepresentation regressionRepresentation_; // define basis to regress on, result save on the coeff ! 
	std::vector<std::vector<double>> basis_value_on_allPath_buffer_;   


public:

	Regression(const RegressionRepresentation& regressionRepresentation);

	virtual ~Regression(){}


	const RegressionRepresentation& getRegressionRepresentation()const{return regressionRepresentation_;} 

	//! move to .cpp
	void constructBasisMatrix(	const std::vector<McLmm_LS::LMMSimulationResult>& lmmSimulationResults, 
								std::vector<std::vector<double>>& basis_value_on_allPath_buffer);  // TODO: YY: which is const, which is not ???  


	void regress(	const std::vector<double>& Z, 
					const std::vector<McLmm_LS::LMMSimulationResult>& lmmSimulationResults,
					std::vector<std::vector<double>>& basis_value_on_allPath_buffer);

	//inversion
	bool InvertMatrix (const ublas::matrix<double>& input, ublas::matrix<double>& inverse);

	void write_to_stream(std::ostream& out)const;

	//Junbin: normalization
	void normalizationMatrix(std::vector<std::vector<double>>& m, std::vector<double>& weight_coeff);
};
typedef boost::shared_ptr<Regression> Regression_PTR;
typedef boost::shared_ptr<const Regression> Regression_CONSTPTR;


}