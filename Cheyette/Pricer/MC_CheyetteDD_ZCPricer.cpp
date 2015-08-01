#include <cassert>
#include <iostream>
#include <cmath>
#include <Cheyette\Pricer\MC_CheyetteDD_ZCPricer.h>



//Price at time t, T maturity of the ZC bond
//double MC_CheyetteDD_ZCPricer::price(double t_valo, double T, size_t nbSimulation)  const
//{
//	double x_t, y_t, value(0) ;
//	std::vector<double> x_t_one_sim, y_t_one_sim ;  //one simulation 
//	std::vector<double> dates = mcCheyette_->getDatesOfSimulation_() ;
//
//	size_t index = numeric::findClosestDate(t_valo, dates) ;
//
//	for(size_t itrSimulation=0; itrSimulation<nbSimulation; ++itrSimulation)
//	{
//		mcCheyette_->simulate_Euler() ;			(T)?	<------------------------------------- TODO -------->
//		x_t_one_sim = mcCheyette_->get_x_t_Cheyette_() ;
//		y_t_one_sim = mcCheyette_->get_y_t_Cheyette_() ;
//
//		x_t = x_t_one_sim[index] ;
//		y_t = y_t_one_sim[index] ;
//
//		value += mcCheyette_->getCheyetteDD_Model_()->P(t_valo, T, x_t, y_t) ; 
//	}
//
//	return value / nbSimulation ;
//
//}

