#pragma once

#include <Cheyette/Calibration/CheyetteBaseCalibrator.h>
#include <Cheyette/Calibration/Cheyette_CostFunctionConvexity.h>

#include <Cheyette/Fonction.h>

class CheyetteQuad_LocalCalibrator : public CheyetteBaseCalibrator
{
private:
	Cheyette_CostFunctionLevel_PTR		costFunctionLevel_PTR_ ;
	Cheyette_CostFunctionSkew_PTR		costFunctionSkew_PTR_ ;
	Cheyette_CostFunctionConvexity_PTR	costFunctionConvexity_PTR_ ;

	//points de depart
	QuantLib::Array	a_Initiate_ /*1D_*/ ;					//param initial pour a(t)				
	QuantLib::Array	b_Initiate_ /*1D_*/ ;					//param initial pour b(t)	
	QuantLib::Array	c_Initiate_ /*1D_*/ ;					//param initial pour c(t)	

	mutable QuantLib::Array	calibrated_a_ /*1D_*/ ;			
	mutable QuantLib::Array	calibrated_b_ /*1D_*/ ;			
	mutable QuantLib::Array	calibrated_c_ /*1D_*/ ;			

	mutable size_t currentSwaptionIndex_ ;				//swaption sur laquelle on minimise
	size_t nbIterations_ ;								//nb de boucles pour les calibrations successives de sigma, m

public :
	CheyetteQuad_LocalCalibrator(	const QuantLib::Size & maxIterations,       
								const QuantLib::Real & rootEpsilon,        
								const QuantLib::Real & functionEpsilon,    
								QuantLib::Array& a_Initiate, QuantLib::Array& b_Initiate, QuantLib::Array& c_Initiate,
								QuantLib::Array& calibrated_a, QuantLib::Array& calibrated_b, QuantLib::Array& calibrated_c,
								Cheyette_CostFunctionLevel_PTR&	costFunctionLevel_PTR,
								Cheyette_CostFunctionSkew_PTR&	costFunctionSkew_PTR,
								Cheyette_CostFunctionConvexity_PTR&	costFunctionConvexity_PTR,
								size_t nbIterations = 2)
		: CheyetteBaseCalibrator(maxIterations, rootEpsilon , functionEpsilon),
			a_Initiate_(a_Initiate), b_Initiate_(b_Initiate),  c_Initiate_(c_Initiate), 
			calibrated_a_(calibrated_a), calibrated_b_(calibrated_b), calibrated_c_(calibrated_c),
			costFunctionLevel_PTR_(costFunctionLevel_PTR), 
			costFunctionSkew_PTR_(costFunctionSkew_PTR), 
			costFunctionConvexity_PTR_(costFunctionConvexity_PTR), 
			currentSwaptionIndex_(1),	//0 ou 1 ??
			nbIterations_(nbIterations)
	{}

	virtual ~CheyetteQuad_LocalCalibrator(){}

	Cheyette_CostFunctionLevel_PTR		getCostFunctionLevel_PTR()		const {return costFunctionLevel_PTR_ ;}
	Cheyette_CostFunctionSkew_PTR		getCostFunctionSkew_PTR()		const {return costFunctionSkew_PTR_ ;}
	Cheyette_CostFunctionConvexity_PTR	getCostFunctionConvexity_PTR()	const {return costFunctionConvexity_PTR_ ;}

	//calibration de a -> b -> c -> a sur une swaption
	virtual void solve() const ;			//minimise pour une swaption	
	//calibration du modèle sur le vecteur des swaptions coterminales
	virtual void calibrate() const ;		
};

typedef boost::shared_ptr<CheyetteQuad_LocalCalibrator>			CheyetteQuad_LocalCalibrator_PTR;
typedef boost::shared_ptr<const CheyetteQuad_LocalCalibrator>	CheyetteQuad_LocalCalibrator_CONSTPTR;

