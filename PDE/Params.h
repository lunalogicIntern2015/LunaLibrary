#pragma once
#include <string>

const double epsilon_compare = 0.0000000001;  //10e-10

//! Fwd PDE
const int if_fwd_PDE = true;


//!-------------------------------------------------------
//!
//!                 HESTON Model
//!
//!-------------------------------------------------------

//! Heston Direction_Parameters
const double Heston_params_kappa = 1.5;
const double Heston_params_theta = 0.04;
const double Heston_params_sigma = 0.3;
const double Heston_params_rho   = -0.5;
const double Heston_params_r     = 0.025;

//! Feller condition: 2*kappa*theta > sigma*sigma

//! Attention: when 2*kappa*theta >= sigma*sigma (condition for CIR)
//! then the vol processus is always positive, which make the Heston SDE valid.

const double Heston_params_num_analytical_fomula_discretization = 1024*4; 

//! Product Direction_Parameters
const double Option_params_K = 100.0; 
const double Option_params_T = 1; 


//! Fwd PDE
const double fwd_PDE_S0 = 100.0;
const double fwd_PDE_v0 = 0.04;

//! Printing
//const std::string DEBUG_output_path = "C:\\Users\\huojie\\Documents\\Visual Studio 2008\\Projects\\PDE\\PDE_D2\\Output\\";
const std::string DEBUG_output_path = "C:\\Users\\yuan LI\\Documents\\Visual Studio 2008\\Projects\\PDE\\PDE_D2\\Output\\";

namespace Direction_Parameters
{
	enum Direction_Parameters
	{
         BC_direction_t, BC_direction_x, BC_direction_y
	};  
}



//!-------------------------------------------------------
//!
//!                 SABR Model
//!
//!-------------------------------------------------------

const double fwd_PDE_SABR_S0   = 0.05;//0.05;
const double fwd_PDE_SABR_v0   = 0.04;

const double SABR_params_alpha = 0.4;
const double SABR_params_beta  = 0.5; // 0.5;
const double SABR_params_rho   = -0.25;



//!-------------------------------------------------------
//!
//!                 SLV SABR Model
//!
//!-------------------------------------------------------

const double fwd_PDE_SLV_SABR_S0   = 0.2; //c'est le taux de change ... 
const double fwd_PDE_SLV_SABR_v0   = 0.04;

const double SLV_SABR_params_alpha = 0.4;
const double SLV_SABR_params_beta  = 0.5;
const double SLV_SABR_params_rho   = -0.25;

const double SLV_SABR_params_r     = 0.04;
