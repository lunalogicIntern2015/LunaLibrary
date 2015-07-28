#include <Cheyette\CheyetteModel\CourbeInput.h>


const double epsilon = 1.e-6;  // bad implementation 

CourbeInput::CourbeInput(std::vector<double> listeMatu, std::vector<double> tauxZC)
	: listeMatu_(listeMatu), tauxZC_(tauxZC)
{
	assert(listeMatu_.size() == tauxZC_.size()) ;
}

double CourbeInput::get_tauxZC0(double T) const
{
//interpolation lineaire	
	return NumericalMethods::linearInterpolation2(T, listeMatu_, tauxZC_) ;
//spline cubique (a tester)
	//std::vector<double> vect_y2_derivees_secondes(listeMatu_.size()) ;
	//double yp1 = 2 * pow(10,30) ;//conditions pour spline cubique naturelle
	//double ypn = 2 * pow(10,30) ;
	//NumericalMethods::spline(listeMatu_, tauxZC_, yp1, ypn, vect_y2_derivees_secondes) ;
	//return NumericalMethods::splineCubique(listeMatu_, tauxZC_, vect_y2_derivees_secondes, T) ;
}

double CourbeInput::get_f_0_t(double t) const
{
	//double yieldT = get_tauxZC0(t)*t;
	//double t_bump = t + epsilon;
	//double yieldT_bump = get_tauxZC0(t_bump)*t_bump;

	//return (yieldT_bump-yieldT)/epsilon;  //Yuan

	double derivee_yield = (get_tauxZC0(t + epsilon) - get_tauxZC0(t))/epsilon ;
	return get_tauxZC0(t) + t * derivee_yield ;

}

void CourbeInput::show() const
{
	int N = listeMatu_.size() ;
	std::cout << "listeMatu   |   tauxZC" << std::endl ;
	for (int i = 0 ; i < N ; ++i)
	{
		std::cout << listeMatu_[i] << "  |  " << tauxZC_[i] << std::endl ; 
	}
}

void CourbeInput::print(std::ostream& o) const
{
	int N = listeMatu_.size() ;
	o << "courbe des taux ZC spot (yield)" << std::endl ;
	o << "listeMatu   ;   tauxZC" << std::endl ;
	for (int i = 0 ; i < N ; ++i)
	{
		o << listeMatu_[i] << " ; " << tauxZC_[i] << std::endl ; 
	}
	o	<<	std::endl;
}

void CourbeInput::printHorizontal(std::ostream& o) const
{
	int N = listeMatu_.size() ;
	o << "courbe des taux ZC spot (yield)" << std::endl ;

	o << "listeMatu ; " << std::endl ;
	for (int i = 0 ; i < N ; ++i)
	{
		o << listeMatu_[i] << " ; " ; 
	}
	o << std::endl ;
	o << "tauxZC ; " << std::endl ;
	for (int i = 0 ; i < N ; ++i)
	{
		o << tauxZC_[i] << " ; " ; 
	}
	o	<<	std::endl;
}