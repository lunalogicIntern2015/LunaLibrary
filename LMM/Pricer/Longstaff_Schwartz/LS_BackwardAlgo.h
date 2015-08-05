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

#include <Instrument/Rate/Rate1.h>   //??
#include <LMM/Pricer/Longstaff_Schwartz/Basis.h>
#include <LMM/Pricer/Longstaff_Schwartz/Basis_Evaluator.h>
#include <LMM/Pricer/Longstaff_Schwartz/Regression_LS.h>
#include <LMM/Pricer/Longstaff_Schwartz/McLmm_LS.h>
#include <Instrument/CallableOption/CallableInstrument.h>
#include <LMM/Pricer/McLmmPricer/McLmmPricer.h>


namespace ublas = boost::numeric::ublas; 

typedef ublas::matrix<double> Rmatrix;
typedef ublas::vector<double> Rvector;

class LS_BackwardAlgo  // regression(path+i) = sum_k coeff_k * basis_k(path_i)
{
    std::vector<LMM::Index>     exerciseDates_;   // i.e regression dates 
	std::vector<LS::Regression> regressions_;     // size = exerciseDates_.size()-1, regression[i] define basis on exerciseDate[i], exerciseDate[last] don't have regression basis.  
     
	std::vector<double>   Z_buffer_;              // size= nbSimulation
	McLmmPricer_CONSTPTR  mcLmmInstrumentPricer_; // price the underlying noncallable instrument ! 

public:
	LS_BackwardAlgo(McLmmPricer_CONSTPTR  mcLmmInstrumentPricer):mcLmmInstrumentPricer_(mcLmmInstrumentPricer){}

	void update_Z(CallableInstrument_CONSTPTR callableInstrument, 
		          LMM::Index exerciseLiborDateindex,
				  LMM::Index regressionIndex,
				  const std::vector<McLmm_LS::LMMSimulationResult>& lmmSimulationResult);
				  // Z_buffer_ = max(iv,df*cv) pathwise


	//! move to: .cpp
	// suppose alreday done fwd simulation, and all simulation datas saved! 
	void do_BackwardAlgo(	CallableInstrument_CONSTPTR callableInstrument, 
							const std::vector<McLmm_LS::LMMSimulationResult>& lmmSimulationResult,
							std::vector<std::vector<double>>& basis_value_on_allPath_buffer);

	//accessor 
	std::vector<LMM::Index>& getExerciseDates(){return exerciseDates_;}
	std::vector<LS::Regression>& getRegressions(){return regressions_;}
	std::vector<double>&   get_Z_buffer(){return Z_buffer_;}
	McLmmPricer_CONSTPTR& getMcLmmInstrumentPricer(){return mcLmmInstrumentPricer_;}
};
typedef boost::shared_ptr<LS_BackwardAlgo> LS_BackwardAlgo_PTR;
typedef boost::shared_ptr<const LS_BackwardAlgo> LS_BackwardAlgo_CONSTPTR;