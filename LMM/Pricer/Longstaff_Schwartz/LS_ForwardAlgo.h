#pragma once
#include <vector>
#include <time.h>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp> 
#include <boost/numeric/ublas/vector_proxy.hpp> 
#include <boost/numeric/ublas/matrix.hpp> 
#include <boost/numeric/ublas/triangular.hpp> 
#include <boost/numeric/ublas/lu.hpp> 

#include <LMM/Helper/Name.h>

#include <Instrument/Rate/Rate1.h>   
#include <LMM/Pricer/Longstaff_Schwartz/Basis.h>
#include <LMM/Pricer/Longstaff_Schwartz/Basis_Evaluator.h>
#include <LMM/Pricer/Longstaff_Schwartz/Regression_LS.h>
#include <LMM/Pricer/Longstaff_Schwartz/McLmm_LS.h>
#include <Instrument/CallableOption/CallableInstrument.h>
#include <LMM/Pricer/McLmmPricer/McLmmPricer.h>

namespace ublas = boost::numeric::ublas; 

typedef ublas::matrix<double> Rmatrix;
typedef ublas::vector<double> Rvector;

class LS_ForwardAlgo  // regression(path+i) = sum_k coeff_k * basis_k(path_i)
{
    std::vector<LMM::Index>   exerciseDates_; // i.e regression dates 
	McLmmPricer_CONSTPTR  mcLmmInstrumentPricer_; // price the underlying noncallable instrument ! 

public:
	//constructor
	LS_ForwardAlgo(McLmmPricer_CONSTPTR  mcLmmInstrumentPricer):mcLmmInstrumentPricer_(mcLmmInstrumentPricer){}


	//! move to: .cpp
	// suppose alreday done fwd simulation, and all simulation datas saved! 
	std::pair<double, double> do_ForwardAlgo(	CallableInstrument_CONSTPTR callableInstrument, 
												const std::vector<McLmm_LS::LMMSimulationResult>& lmmSimulationResult,
												std::vector<LS::Regression>& regressions);

	double do_ForwardAlgo_on_onePath(	CallableInstrument_CONSTPTR callableInstrument, 
										const McLmm_LS::LMMSimulationResult& lmmSimulationResult,
										std::vector<LS::Regression>& regressions,
										const std::vector<Instrument_CONSTPTR>&  instrument_vector);

	
	//accessor
	std::vector<LMM::Index>& getExerciseDates(){return exerciseDates_;}
};
typedef boost::shared_ptr<LS_ForwardAlgo> LS_ForwardAlgo_PTR;
typedef boost::shared_ptr<const LS_ForwardAlgo> LS_ForwardAlgo_CONSTPTR;