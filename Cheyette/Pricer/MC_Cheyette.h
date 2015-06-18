#pragma once

#include <iostream>
#include <cassert>
#include <string>

#include <Cheyette/CheyetteModel/CheyetteDD_Model.h>
#include <LMM/LmmModel/McLmm.h>  //pour le namespace MCSchemeType
#include <LMM/RNGenerator/RNGenerator.h>
#include <LMM/RNGenerator/McGenerator.h>
#include <LMM/helper/TenorType.h>
//simulation MC de x_t et y_t sous la mesure forward Q^T (T = fwdProbaT_)





class MC_Cheyette
{
private:
	CheyetteDD_Model_PTR		cheyetteDD_Model_ ;
	RNGenerator_PTR				rnGenerator_;
	Tenor						tenorType_ ;				//tenor de la Tenor Structure / indices where flows occur
	double						fwdProbaT_;					// define the numeraire

	//dates ou on veut recuperer le vecteur de simulation :
	std::vector<size_t>			indexOfSimulation_;				//indices 
	std::vector<size_t>			discretizationBetweenDates_ ;	//nb de simus entre 2 dates : 100 simus, ..., 100 simus de taille N

	//! comes from simulation
	mutable std::vector<double> x_t_Cheyette_;
	mutable std::vector<double> y_t_Cheyette_ ;
	mutable std::vector<double> numeraires_; // numeraire[i] = P(T_i,T_{N+1})  , size = N+1  

public:

	MC_Cheyette(	CheyetteDD_Model_PTR		cheyetteDD_Model,
					RNGenerator_PTR				rnGenerator,
					Tenor						tenorType,
					double						fwdProbaT,
					std::vector<size_t>&		indexOfSimulation,		
					std::vector<size_t>&		discretizationBetweenDates   )
				:cheyetteDD_Model_(cheyetteDD_Model), rnGenerator_(rnGenerator), 
				tenorType_(tenorType), fwdProbaT_(fwdProbaT), 
				indexOfSimulation_(indexOfSimulation), discretizationBetweenDates_(discretizationBetweenDates), 
				x_t_Cheyette_(indexOfSimulation_.size()), y_t_Cheyette_(indexOfSimulation_.size()), numeraires_(indexOfSimulation_.size()) 
	{
			//remplit x_t_Cheyette_ et y_t_Cheyette_ aux dates voulues
			//simulate_Euler() ;
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

	//remplit les vecteurs xt et yt aux dates de la tenor structure sous la mesure Q^T (T = fwdProbaT_)
	void simulate_Euler() const;
	void computeNumeraires() const ;


};

typedef boost::shared_ptr<MC_Cheyette> MC_Cheyette_PTR;
typedef boost::shared_ptr<const MC_Cheyette> MC_Cheyette_CONSTPTR;
