//approximation du call dans le quadratic vol model

#include <Cheyette/Fonction.h>
#include <RNGenerator/RNGenerator.h>
#include <RNGenerator/McGenerator.h>
#include <Numeric\Integrator1D.h>

#include <Cheyette/Pricer/Cheyette_SwaptionPricer_QuadApprox.h>

#include <Cheyette/unit_test/Test_CheyetteQuad_Model.h>
#include <Cheyette/unit_test/TestApproxDD.h>


Cheyette_SwaptionPricer_QuadApprox_PTR createQuadApproxPricer_PTR(size_t xmax,  int numCourbe, 
																double k, double aValue, double bValue, double cValue) ;

Cheyette_SwaptionPricer_QuadApprox_PTR createQuadApproxPricer_PTR(size_t xmax,  CourbeInput_PTR pCourbeInput, 
																double k, double aValue, double bValue, double cValue) ;

void test_MC_approx() ;

void testCallMC() ;
void testCallApprox() ;

//Euler sur dS(t) = Phi(S_t) dW_t
std::vector<double> approxMC_call(const size_t nbSimus, const double strike, const double S0, const double T) ;


double phi(double St) ;
double phiP(double St) ;   //phi prime
double phiS(double St) ;	 //phi seconde
double phiMoinsUn(double St) ;

VanillaSwaption_PTR setSwaptionATM_Quad(const CheyetteQuad_Model_PTR modele_test_PTR,
										const Tenor& floatTenor, const Tenor& fixedTenor, 
										const size_t a, const size_t b) ;