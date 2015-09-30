#pragma once

#include <Cheyette/Calibration/CheyetteBaseCalibrator.h>
#include <Cheyette/Fonction.h>

class CheyetteDD_LocalCalibrator : public CheyetteBaseCalibrator
{
private:
	Cheyette_CostFunctionLevel_PTR		costFunctionLevel_PTR_ ;
	Cheyette_CostFunctionSkew_PTR		costFunctionSkew_PTR_ ;
	
	QuantLib::Array	sigmaInitiate_ /*1D_*/ ;					//param initial pour sigma(t)				
	QuantLib::Array	mInitiate_ /*1D_*/ ;						//param initial pour m(t)	

	mutable QuantLib::Array	calibrated_sigma_ /*1D_*/ ;		//sigma(t) calibre				
	mutable QuantLib::Array	calibrated_m_ /*1D_*/ ;			//m(t) calibre

	mutable size_t currentSwaptionIndex_ ;		//swaption sur laquelle on minimise
	size_t nbIterations_ ;	//nb de boucles pour les calibrations successives de sigma, m

public :
	CheyetteDD_LocalCalibrator(	const QuantLib::Size & maxIterations,       
								const QuantLib::Real & rootEpsilon,        
								const QuantLib::Real & functionEpsilon,    
								QuantLib::Array& sigmaInitiate, QuantLib::Array& mInitiate,
								QuantLib::Array& calibrated_sigma, QuantLib::Array& calibrated_m,
								Cheyette_CostFunctionLevel_PTR&	costFunctionLevel_PTR,
								Cheyette_CostFunctionSkew_PTR&	costFunctionSkew_PTR,
								size_t nbIterations = 2)
		: CheyetteBaseCalibrator(maxIterations, rootEpsilon , functionEpsilon),
			sigmaInitiate_(sigmaInitiate), mInitiate_(mInitiate), 
			calibrated_sigma_(calibrated_sigma), calibrated_m_(calibrated_m),
			costFunctionLevel_PTR_(costFunctionLevel_PTR), costFunctionSkew_PTR_(costFunctionSkew_PTR), 
			currentSwaptionIndex_(1),	//0 ou 1 ??
			nbIterations_(nbIterations)
	{}

	virtual ~CheyetteDD_LocalCalibrator(){}

	Cheyette_CostFunctionLevel_PTR		getCostFunctionLevel_PTR() const {return costFunctionLevel_PTR_ ;}
	Cheyette_CostFunctionSkew_PTR		getCostFunctionSkew_PTR() const {return costFunctionSkew_PTR_ ;}

	//calibration de sigma -> m -> sigma sur une swaption
	virtual void solve() const ;			//minimise pour une swaption	
	virtual void calibrate() const ;		//calibration du modèle sur le vecteur des swaptions coterminales
};

typedef boost::shared_ptr<CheyetteDD_LocalCalibrator> CheyetteDD_LocalCalibrator_PTR;
typedef boost::shared_ptr<const CheyetteDD_LocalCalibrator> CheyetteDD_LocalCalibrator_CONSTPTR;

