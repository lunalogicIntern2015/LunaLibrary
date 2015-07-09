#pragma once

// CheyetteDD_CalibrationConfig plays the role of parameters holders for the whole model
//
// 
struct CheyetteDD_CalibrationConfig 
{	
public:
	
	CheyetteDD_CalibrationConfig();	
	
	//param du modele
//	Shifted_HGVolatilityParam::ABCDParameter vol_abcd_;
	CheyetteDD_Model::CheyetteDD_Parameter cheyetteDD_Param_ ;
	
	CourbeInput_PTR courbeInput_PTR_ ;

	size_t shiftChoice_ ;

	size_t model_nbYear_;  //nb year max, 16 pour calibrer par ex


//autres param fixés
//	double correl_alpha_;
//	double correl_beta_;

	std::string test_folder_;

	mutable bool use_local_calib_;
	bool use_positive_constraint_;

	virtual void print(const std::string& filename) const ;	

	//tranformer le struct cheyette dd param en Array pour la calibration
	QuantLib::Array get_abcdArray() const
	{
		QuantLib::Array abcd(4);
		abcd[0]=vol_abcd_.a_;
		abcd[1]=vol_abcd_.b_;
		abcd[2]=vol_abcd_.c_;
		abcd[3]=vol_abcd_.d_;
		return abcd;
	}

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

