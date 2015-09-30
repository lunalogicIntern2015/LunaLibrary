#include "CheyetteQuad_Model.h"

//Quadratic : k constant, calcul explicite de G(t, T)
double CheyetteQuad_Model::G(double t, double T) const
{
	assert(t >= 0);  
	return 1/k_ * (1 - exp(- k_ * (T-t))) ;
}


void CheyetteQuad_Model::show() const
{
	CheyetteModel::show() ;

	std::cout << "--- Quadratic Vol ---" << std::endl ;
	std::cout << "  A piecewise constant " << std::endl ;
	a_.show() ;
	std::cout << "  B piecewise constant " << std::endl ;
	b_.show() ;
	std::cout << "  C piecewise constant " << std::endl ;
	c_.show() ;

	std::cout << "  k constant " << std::endl ;
	std::cout << "       k = " << k_ << std::endl ;
	std::cout << "---------------------------------------------" << std::endl ;

}

void CheyetteQuad_Model::print(std::ostream& o) const 
{
	CheyetteModel::print(o) ;
	o << "------------- Quadratic Vol ------------" << std::endl ;
	o << "  A piecewise constant " << std::endl ;
	a_.print(o) ;
	o << "  B piecewise constant " << std::endl ;
	b_.print(o) ;
	o << "  C piecewise constant " << std::endl ;
	c_.print(o) ;
	o << "  k constant " << std::endl ;
	o << "       k = ; " << k_ << std::endl ;
	o << "---------------------------------------------" << std::endl ;
}





