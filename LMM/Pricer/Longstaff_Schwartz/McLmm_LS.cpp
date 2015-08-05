#include <LMM/Pricer/Longstaff_Schwartz/McLmm_LS.h>

void McLmm_LS::simulateLMM(size_t nbSimulations)
{
	if(lmmSimualtionResults_.size()!=nbSimulations)
	{
		lmmSimualtionResults_.resize(nbSimulations);
	}

	for(size_t i=0; i<nbSimulations; ++i)
	{
		mcLmm_->simulateLMM();

		//!!!  Problem YY: Very bad code here: to much copy-coller !!!!!! 
		//??? need to change result to smart-ptr + need to construct LMMSimulationResult in the lmm class ???
		lmmSimualtionResults_[i] = LMMSimulationResult(mcLmm_->get_liborMatrix(), mcLmm_->get_numeraire()); 
	}
}

void McLmm_LS::write_to_stream(std::ostream& out)const
{
	out << "lmmSimualtionResults_: ;" << endl;
	for(size_t i = 0; i<lmmSimualtionResults_.size(); i++)
	{
		out << "lmmSimualtionResults " << i << endl ;
		lmmSimualtionResults_[i].write_to_stream(out);
	}
	out << endl;	
}

void McLmm_LS::LMMSimulationResult::write_to_stream(std::ostream& out)const
{

	out << "LiborMatrix_: ;" << endl;
	for(size_t i = 0; i<LiborMatrix_.size1(); i++)
	{
		out << " ;" ;
		for(size_t j = 0; j<LiborMatrix_.size2(); j++)
			out << LiborMatrix_(i,j) << " ;" ;
		out << endl;
	}
	out << endl;

	out << "numeraire_: ;" ;
	for(size_t i = 0; i<numeraire_.size(); i++)
		out << numeraire_[i] << " ;" ;
	out << endl;

}