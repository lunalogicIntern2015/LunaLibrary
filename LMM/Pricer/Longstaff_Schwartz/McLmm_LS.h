#pragma once
#include <vector>
#include <Map>

#include <LMM/Mc/McLmm.h>


//! TODO: YY not efficient: save all the Libor matrix, need to do optimization latter !!! 
//! TODO: only need to save the datas for the regression dates + not all libor but only the ones needed to evaluate basis !!! 
//          so do a pre-selection of the datas (simulated libor ) ??? 
class McLmm_LS  // forward simulate LMM, and save all the datas for regression. 
{
	 McLmm_PTR mcLmm_;
public:
	class LMMSimulationResult  // Attention: this involves copy-coller, better to use smart-ptr, to avoid this copy-coller ???
	{
		matrix LiborMatrix_;
		std::vector<double> numeraire_;

	public:
		LMMSimulationResult(){}
		LMMSimulationResult(const matrix& LiborMatrix, const std::vector<double>& numeraire)
			:
			LiborMatrix_(LiborMatrix), 
			numeraire_(numeraire)
		{}

		virtual ~ LMMSimulationResult() {}

		const matrix& get_liborMatrix() const{return LiborMatrix_;}
		const std::vector<double>& get_numeraire() const{return numeraire_;}

		void write_to_stream(std::ostream& out)const;
	};

	//! TODO : YY: This is not efficient, try to improve it latter ! 
	std::vector<LMMSimulationResult>  lmmSimualtionResults_;

	//! constructor 
	McLmm_LS(McLmm_PTR mcLmm):mcLmm_(mcLmm){}

	//! destructor 
	virtual ~McLmm_LS(){}; 

	void simulateLMM(size_t nbSimulation);

	void write_to_stream(std::ostream& out)const;
};
typedef boost::shared_ptr<McLmm_LS> McLmm_LS_PTR;
typedef boost::shared_ptr<const McLmm_LS> McLmm_LS_CONSTPTR;
