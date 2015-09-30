#include "Test_CheyetteQuad_Model.h"

//xmax : matu maximum
CheyetteQuad_Model_PTR creeCheyetteQuad_Modele_PTR(size_t xmax, int numCourbe, 
												   double k, double aValue, double bValue, double cValue)
{
	CourbeInput_PTR pCourbeInput(createCourbeInput(numCourbe)) ;

	CheyetteQuad_Model_PTR modele_test_PTR = creeCheyetteQuad_Modele_PTR(xmax, pCourbeInput, 
																		 k, aValue, bValue, cValue) ;

	return modele_test_PTR ;
}


CheyetteQuad_Model_PTR creeCheyetteQuad_Modele_PTR(size_t xmax, CourbeInput_PTR pCourbeInput, 
												   double k, double aValue, double bValue, double cValue)
{
	std::vector<double> x(xmax) ;				//a+b + 1
	std::vector<double> a_y(xmax - 1) ;			//a+b
	std::vector<double> b_y(xmax - 1) ;			//a+b
	std::vector<double> c_y(xmax - 1) ;			//a+b
	for (size_t i = 0 ; i <= xmax - 1 ; ++i)
	{
		x[i] = i ;
	}
	for (size_t i = 0 ; i < xmax - 1 ; ++i)
	{
		a_y[i] = aValue ;
		b_y[i] = bValue ;
		c_y[i] = cValue ;
	}
	Piecewiseconst_RR_Function aFunc(x, a_y) ; 
	Piecewiseconst_RR_Function bFunc(x, b_y) ; 
	Piecewiseconst_RR_Function cFunc(x, c_y) ; 

	CheyetteQuad_Model_PTR modele_test_PTR(new CheyetteQuad_Model(pCourbeInput, k, aFunc, bFunc, cFunc)) ;

	return modele_test_PTR ;
}
