#pragma once

#include <Cheyette/Model/CheyetteDD_Model.h>
#include <fstream>

// CheyetteDD_CalibrationConfig plays the role of parameters holders for the whole model
 
struct CheyetteDD_CalibrationConfig 
{	
public:

//param
	size_t									model_nbYear_;  //nb year max, 16 pour calibrer par ex
	CheyetteDD_Model::CheyetteDD_Parameter	cheyetteDD_Param_ ;
	double									strikeBump_ ;
	CourbeInput_PTR							courbeInput_PTR_ ;  
	size_t									shiftChoice_ ;
	std::string								test_folder_;

	mutable bool							use_local_calib_;
	bool									use_positive_constraint_;
	mutable double result_quote_error_l2;
	mutable double result_quote_error_l1;
	mutable double result_quote_error_linf;

	mutable double result_pelTime_error_l2;
	mutable double result_pelTime_error_l1;
	mutable double result_pelTime_error_linf;

	mutable double result_pelLibor_error_l2;
	mutable double result_pelLibor_error_l1;
	mutable double result_pelLibor_error_linf;

//methodes
	CheyetteDD_CalibrationConfig(	CheyetteDD_Model::CheyetteDD_Parameter	cheyetteDD_Param,
									CourbeInput_PTR							courbeInput_PTR) ;
	void reset_result() const ;
	virtual void print(const std::string& filename) const ;	
};

