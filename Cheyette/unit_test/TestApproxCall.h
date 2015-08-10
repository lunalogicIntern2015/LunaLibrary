//approximation du call dans le quadratic vol model

#include <Cheyette\Fonction.h>
#include <RNGenerator/RNGenerator.h>
#include <RNGenerator/McGenerator.h>

//Euler sur dS(t) = Phi(S_t) dW_t

double oneTrajEuler(const double S0, const double T, const Boost_RR_Function_PTR& funcPhi) ;

std::vector<double> approxMC_call(const size_t nbSimus, const double S0, const double T) ;
//, const Boost_RR_Function_PTR& funcPhi) ;

double phi(double St) ;