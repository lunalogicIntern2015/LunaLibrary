#pragma once 
#include <iostream>
#include <string>
#include <vector>
using namespace std;

//! usual Mathmatical model

class Model
{
public:
	vector<string> list_modelType_;
	Model()
	{
		list_modelType_.push_back("Heston");
		list_modelType_.push_back("SABR");
	};
	virtual string get_modelType(){return "Abstract!";}// = 0;
	virtual ~Model(){};
};




//class SLV_Heston_Model: public Model
//{
//public:
//	string modelType_;
//
//    double r_;
//	double sigma_;
//	double kappa_;
//	double theta_;
//	double rho_;
//	double lambda_;
//
//	SLV_Heston_Model(double r, double sigma, double kappa, double theta, double rho, double lambda = 0)
//	{	   
//	       modelType_ = "SLV_Heston_";
//
//		   r_      = r;
//		   sigma_  = sigma;
//		   kappa_  = kappa;
//		   theta_  = theta;
//		   rho_    = rho;
//		   lambda_ = lambda;
//		   //! maybe delete this constrain later ...
//		   if(lambda_!=0)
//		   {
//		       throw("Error in constructor of SLV_Heston_Model, lambda should be 0.0");
//		   }
//	}
//	~SLV_Heston_Model(){};
//
//	string get_modelType(){return modelType_;}
//
//	double get_r()      const {return r_;}
//	double get_sigma()  const {return sigma_;}
//	double get_kappa()  const {return kappa_;}
//	double get_theta()  const {return theta_;}
//	double get_rho()    const {return rho_;}
//	double get_lambda() const {return lambda_;}
//};


//class MaModel: public Model
//{
//public:
//	string modelType_;
//
//	double alpha_;
//	double sigma_;
//	double rho_;
//
//	MaModel(double alpha, double sigma, double rho)
//	{	   
//	       modelType_ = "MaModel";
//
//	 alpha_ = alpha;
//	 sigma_ = sigma;
//	 rho_   = rho ;
//
//	}
//	~MaModel(){};
//
//	string get_modelType(){return modelType_;}
//
//	double get_alpha()  const {return alpha_;}
//	double get_sigma()  const {return sigma_;}
//	double get_rho()    const {return rho_;}
//
//};
//
