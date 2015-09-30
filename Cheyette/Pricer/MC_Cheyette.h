#pragma once

#include <iostream>
#include <cassert>
#include <string>
#include <boost/shared_ptr.hpp>

#include <Cheyette/Model/CheyetteModel.h>

#include <RNGenerator/RNGenerator.h>
#include <RNGenerator/McGenerator.h>
#include <LMM/Helper/TenorType.h>
#include <LMM/Helper/LMMTenorStructure.h>




/**************  generation de x(t) et y(t) sous la proba forward Q_T  **************
*
*	dx(t) = [ y(t) - k x(t) - G(t, T) sigma_r(t)^2 ] dt + sigma_r(t) dW(t)^QT
*	dy(t) = ( sigma_r(t)^2 - 2 k y(t) ) dt
*
**************************************************************************************/

class MC_Cheyette
{
protected:
	CheyetteModel_PTR			cheyetteModel_PTR_ ;

	RNGenerator_PTR				rnGenerator_;

// UNE UNIQUE TENOR STRUCTURE pour tous les flux (même si 3M / 6M) //
	//LMMTenorStructure(const Tenor&  tenorType, int max_nbOfYear);
	LMMTenorStructure_PTR		pTenorStructure_ ;	

	double						fwdProbaT_ ;					// define the numeraire

	size_t						discretizationBetweenDates_ ;	//nb de pas de discretisation entre 2 dates

	//! comes from simulation
	//x_t, y_t simulés sur la grille TenorStructure
	mutable std::vector<double> x_t_Cheyette_ ;
	mutable std::vector<double> y_t_Cheyette_ ;
	mutable std::vector<double> numeraires_; // numeraire[i] = P(T_i, fwd_probaT)  , size = N

public:

	MC_Cheyette(	CheyetteModel_PTR			cheyetteModel_PTR,
					RNGenerator_PTR				rnGenerator,
					LMMTenorStructure_PTR		pTenorStructure,
					double						fwdProbaT,
					size_t						discretizationBetweenDates   )
		:	cheyetteModel_PTR_(cheyetteModel_PTR)
			, rnGenerator_(rnGenerator)
			, pTenorStructure_(pTenorStructure)
			, fwdProbaT_(fwdProbaT)
			, discretizationBetweenDates_(discretizationBetweenDates)
			, x_t_Cheyette_(static_cast<size_t>(fwdProbaT / pTenorStructure->get_tenorType().YearFraction() + 1))
			, y_t_Cheyette_(static_cast<size_t>(fwdProbaT / pTenorStructure->get_tenorType().YearFraction() + 1))
			, numeraires_(static_cast<size_t>(fwdProbaT / pTenorStructure->get_tenorType().YearFraction() + 1)) 
	{
			assert(pTenorStructure->get_max_nbOfYear() >= fwdProbaT)	;
	} 

	virtual ~MC_Cheyette(){} 
	
	//-- Getters 
	RNGenerator_PTR				getRNGenerator() const{return rnGenerator_ ;}
	LMMTenorStructure_PTR		getTenorStructure() const {return pTenorStructure_ ;}
	double						getFwdProbaT() const{return fwdProbaT_ ;}
	size_t						getDiscretizationBetweenDates() const{return discretizationBetweenDates_ ;}

	std::vector<double>			get_x_t_Cheyette() const{return x_t_Cheyette_ ;}
	std::vector<double>			get_y_t_Cheyette() const{return y_t_Cheyette_ ;}
	std::vector<double>			getNumeraires() const{return numeraires_ ;}

	//remplit les vecteurs xt et yt aux dates de la tenor structure sous la mesure Q^T (T = fwdProbaT_)
	void simulate_Euler() const ;
	void computeNumeraires() const ;

};
typedef boost::shared_ptr<MC_Cheyette> MC_Cheyette_PTR;
typedef boost::shared_ptr<const MC_Cheyette> MC_Cheyette_CONSTPTR;

