#pragma once

#include <Cheyette/CheyetteModel/CheyetteDD_Model.h>
#include <fstream>

// CheyetteDD_CalibrationConfig plays the role of parameters holders for the whole model
 
struct CheyetteDD_CalibrationConfig 
{	
public:
	CheyetteDD_CalibrationConfig();	
	
//param du modele
	CheyetteDD_Model::CheyetteDD_Parameter cheyetteDD_Param_ ;
	CourbeInput_PTR courbeInput_PTR_ ;
	size_t shiftChoice_ ;

	size_t model_nbYear_;  //nb year max, 16 pour calibrer par ex

	std::string test_folder_;

	mutable bool use_local_calib_;
	bool use_positive_constraint_;

	virtual void print(const std::string& filename) const ;	

	// storing result
	void reset_result() const ;

	mutable double result_quote_error_l2;
	mutable double result_quote_error_l1;
	mutable double result_quote_error_linf;

	mutable double result_pelTime_error_l2;
	mutable double result_pelTime_error_l1;
	mutable double result_pelTime_error_linf;

	mutable double result_pelLibor_error_l2;
	mutable double result_pelLibor_error_l1;
	mutable double result_pelLibor_error_linf;
};

