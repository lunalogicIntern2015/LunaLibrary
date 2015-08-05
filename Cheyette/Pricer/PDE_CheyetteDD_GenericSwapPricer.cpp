#include <cassert>
#include <iostream>
#include <cmath>

#include <Cheyette/Pricer/PDE_CheyetteDD_GenericSwapPricer.h>
#include <PDE/BoundaryCondition/BoundaryCondition_D2.h>
#include <Instrument/GenericSwap/GenericSwap.h>
//
////simulation
//double PDE_CheyetteDD_GenericSwapPricer::swapNPV (GenericSwap_PTR genericSwap, size_t valuationIndex) const
//{
//	//! construct boundary condition
//	TimeInitial_BoundaryCondition_PTR bc_I = // depend on product
//	BoundaryCondition_D2_CONSTPTR     bc =  BoundaryCondition_D2_Factory::create_BC_D2_space1stDerivative0(bc_I, discretization_);                  // need to construct it here   
//
//
//}
//
