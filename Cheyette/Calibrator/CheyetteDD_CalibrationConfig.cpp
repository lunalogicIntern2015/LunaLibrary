#include "CheyetteDD_CalibrationConfig.h"

#include <LMM/helper/GenericPath.h>

//#include <cassert>

CheyetteDD_CalibrationConfig::CheyetteDD_CalibrationConfig()
	: model_nbYear_(16),
	//, vol_abcd_(0.17,0.37,1.12,0.3)  //// a very "normal" curve in Rebonato.2009 p.13
	//, correl_alpha_(0.2)
	
	shiftChoice_(1),
	test_folder_("calibration"),
	use_local_calib_(false),
	use_positive_constraint_(false)
{
	reset_result();
}

void CheyetteDD_CalibrationConfig::print(const std::string& filename) const 
{
	std::string path_OutPut = LMMPATH::get_output_path() + filename;

	std::ofstream outputstream;
	outputstream.open(path_OutPut.c_str());

	outputstream<<"CheyetteDD_CalibrationConfig parameters holder"<<std::endl;

	outputstream<<"model_nbYear_ , "<<model_nbYear_<<std::endl; 
//	outputstream<<"a,b,c,d, "<<std::endl; 
//	outputstream<<vol_abcd_.a_<<","<<vol_abcd_.b_<<","<<vol_abcd_.c_<<","<<vol_abcd_.d_<<","<<std::endl; 
	outputstream<<"shiftChoice_, " << shiftChoice_ << std::endl;

	if(use_local_calib_) outputstream<<"Use Local Calibraion , " <<std::endl; 
	
	outputstream<<"Calibration Constraint ,"<<std::endl;
	outputstream<<"Positive Constraint ,,";
	if(use_positive_constraint_) outputstream<<"YES , " <<std::endl; else  outputstream<<"NON , " <<std::endl; 
	
	outputstream.close();
}

void CheyetteDD_CalibrationConfig::reset_result() const 
{
	/// errors always have to be positive or null after computing
	/// initiate errors to small negative, if result keep negative, that is a problem when calculating error

	result_quote_error_l2      =-1e-9;
	result_quote_error_l1      =-1e-9;
	result_quote_error_linf    =-1e-9;

	result_pelTime_error_l2    =-1e-9;
	result_pelTime_error_l1    =-1e-9;
	result_pelTime_error_linf  =-1e-9;

	result_pelLibor_error_l2   =-1e-9;
	result_pelLibor_error_l1   =-1e-9;
	result_pelLibor_error_linf =-1e-9;
}