#include "CheyetteDD_CalibrationConfig.h"

#include <LMM/helper/GenericPath.h>

//#include <cassert>

CheyetteDD_CalibrationConfig::CheyetteDD_CalibrationConfig(CheyetteDD_Model::CheyetteDD_Parameter	cheyetteDD_Param,
														   CourbeInput_PTR							courbeInput_PTR)
	: model_nbYear_(3),		//!!! aussi present dans testCalibrator .cpp ligne 97 
	cheyetteDD_Param_(cheyetteDD_Param),
	courbeInput_PTR_(courbeInput_PTR),
	strikeBump_(5./10000),
	shiftChoice_(1),
	test_folder_("calibration"),
	use_local_calib_(true),
	use_positive_constraint_(true)
{
	reset_result();	//initialisation des autres attributs du struct
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

void CheyetteDD_CalibrationConfig::print(const std::string& filename) const 
{
	std::string path_OutPut = LMMPATH::get_output_path() + filename;

	std::ofstream outputstream;
	outputstream.open(path_OutPut.c_str());

	outputstream<<"CheyetteDD_CalibrationConfig parameters holder"<<std::endl;

	outputstream<<"model_nbYear_ , "<<model_nbYear_<<std::endl; 
	outputstream<<"shiftChoice_, " << shiftChoice_ << std::endl;

	if(use_local_calib_) outputstream<<"Use Local Calibraion , " <<std::endl; 
	
	outputstream<<"Calibration Constraint ,"<<std::endl;
	outputstream<<"Positive Constraint ,,";
	if(use_positive_constraint_) outputstream<<"YES , " <<std::endl; else  outputstream<<"NON , " <<std::endl; 
	
	outputstream.close();
}

