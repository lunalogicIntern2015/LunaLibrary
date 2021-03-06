#pragma once

#include <vector>
#include "assert.h"
#include <boost/shared_ptr.hpp>
#include <iostream>

#include <Numeric/NumericalMethods.h>
//c'est la courbe spot des taux ZC

//Singleton...
class CourbeInput
{
private:
	std::vector<double> listeMatu_ ;
	std::vector<double> tauxZC_ ;  //yield

public:
	CourbeInput(){} //utile pour la classe test

	CourbeInput(std::vector<double> listeMatu, std::vector<double> tauxZC);
	virtual ~CourbeInput(void){}

	double get_tauxZC0(double T) const ;
	double get_f_0_t(double t) const ;
	void show() const ;
	void print(std::ostream& o) const ;
	void printHorizontal(std::ostream& o) const ;

};

typedef boost::shared_ptr<CourbeInput> CourbeInput_PTR;
typedef boost::shared_ptr<const CourbeInput> CourbeInput_CONSTPTR;