#pragma once

#include <iostream>
#include <cassert>
#include <string>
#include <boost/shared_ptr.hpp>

#include <Cheyette/CheyetteModel/CheyetteDD_Model.h>
#include <LMM/LmmModel/McLmm.h>  //pour le namespace MCSchemeType
#include <LMM/RNGenerator/RNGenerator.h>
#include <LMM/RNGenerator/McGenerator.h>
#include <LMM/helper/TenorType.h>

#include <JBLMM/Element/Coupon.h>
#include <JBLMM/Element/CouponLeg.h>
#include <JBLMM/Element/CappedFlooredCoupon.h>
#include <JBLMM/Element/LiborRate.h>



class MC_Cheyette
{
protected:
	CheyetteDD_Model_PTR		cheyetteDD_Model_ ;
	RNGenerator_PTR				rnGenerator_;

// UNE UNIQUE TENOR STRUCTURE pour tous les flux (même si 3M / 6M) //
	Tenor						tenorType_ ;				//tenor de la Tenor Structure / indices where flows occur
	double						fwdProbaT_;					// define the numeraire

	//dates ou on veut recuperer le vecteur de simulation :
	//indexOfSimulation_ : indices, pas dates. (Retour aux dates avec l'attribut tenorType_)
	//le 1er element doit etre 0 (x0 = 0 et y0 = 0)
	std::vector<size_t>			indexOfSimulation_;
	//discretizationBetweenDates_ : nb de simus entre 2 dates : 100 simus, ..., 100 simus de taille N
	//doit etre de meme taille que indexOfSimulation_
	//le 1er élément devrait valoir 0 (0 simulation pour calculer x0 et y0)
	std::vector<size_t>			discretizationBetweenDates_ ;	

	//! comes from simulation
	mutable std::vector<double> x_t_Cheyette_ ;
	mutable std::vector<double> y_t_Cheyette_ ;
	mutable std::vector<double> numeraires_; // numeraire[i] = P(T_i,T_{N+1})  , size = N

public:

	MC_Cheyette(	CheyetteDD_Model_PTR		cheyetteDD_Model,
					RNGenerator_PTR				rnGenerator,
					Tenor						tenorType,
					double						fwdProbaT,
					std::vector<size_t>&		indexOfSimulation,		
					std::vector<size_t>&		discretizationBetweenDates   )
		:cheyetteDD_Model_(cheyetteDD_Model), rnGenerator_(rnGenerator), tenorType_(tenorType), fwdProbaT_(fwdProbaT), 
		indexOfSimulation_(indexOfSimulation), discretizationBetweenDates_(discretizationBetweenDates), 
		x_t_Cheyette_(indexOfSimulation_.size()), y_t_Cheyette_(indexOfSimulation_.size()), numeraires_(indexOfSimulation_.size()) 
	{
			if (indexOfSimulation_.size() != discretizationBetweenDates_.size())
				throw "In MC_Cheyette, vectors datesOfSimulation and discretizationBetweenDates must have the same size" ;		
	} 

	virtual ~MC_Cheyette(){} 
	
	//-- Getters 
	CheyetteDD_Model_PTR		getCheyetteDD_Model() const{return cheyetteDD_Model_ ;}
	RNGenerator_PTR				getRNGenerator() const{return rnGenerator_ ;}
	Tenor						getTenorType() const {return tenorType_ ;}
	double						getFwdProbaT() const{return fwdProbaT_ ;}
	std::vector<size_t>			getIndexOfSimulation() const{return indexOfSimulation_ ;}
	std::vector<size_t>			getDiscretizationBetweenDates() const{return discretizationBetweenDates_ ;}
	std::vector<double>			get_x_t_Cheyette() const{return x_t_Cheyette_ ;}
	std::vector<double>			get_y_t_Cheyette() const{return y_t_Cheyette_ ;}
	std::vector<double>			getNumeraires() const{return numeraires_ ;}

	static size_t findIndex(size_t fixingIndex, std::vector<size_t> indexVector) ; 

	//remplit les vecteurs xt et yt aux dates de la tenor structure sous la mesure Q^T (T = fwdProbaT_)
	void simulate_Euler() const;
	void computeNumeraires() const ;

	//! one simulation - pour les produits vanille derives
	double evaluateFloatLeg(	const size_t indexValuationDate,
								const std::vector<size_t>& indexFloatLeg,
								const std::vector<double>& numeraire, 
								const std::vector<double>& x_t, 
								const std::vector<double>& y_t,
								const Tenor tenorFloatLeg) const;

	double evaluateFixedLeg(	const size_t indexValuationDate,
								const std::vector<size_t>& indexFixedLeg,
								const std::vector<double>& numeraire, 
								const std::vector<double>& x_t, 
								const std::vector<double>& y_t,
								const Tenor tenorFixedLeg, 
								const double fixedRate) const;

	//! one simulation - pour les produits generiques derives
	double evaluateCouponLeg(	const LMM::Index indexValuationDate,
								const CouponLeg_CONSTPTR couponLeg,
								const std::vector<double>& numeraire, 
								const std::vector<double>& x_t, 
								const std::vector<double>& y_t,
								const Tenor tenorLeg) const;      // !!! tenorLeg peut etre different de TenorStructure !!! 

	double evaluateCappedFlooredCoupon( CappedFlooredCoupon_CONSTPTR cappedFlooredCoupon, 
										double rateValue) const ;


};

typedef boost::shared_ptr<MC_Cheyette> MC_Cheyette_PTR;
typedef boost::shared_ptr<const MC_Cheyette> MC_Cheyette_CONSTPTR;
