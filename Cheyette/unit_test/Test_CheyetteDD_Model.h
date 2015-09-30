#pragma once

#include <Cheyette\Model\CheyetteDD_Model.h>
#include <Cheyette\unit_test\TestApproxDD.h>

CourbeInput_PTR createCourbeInput(int curveChoice) ;

//pour les tests sur les courbes fictives
CheyetteDD_Model_PTR creeCheyetteDD_Modele_PTR(size_t xmax, int numCourbe, 
											   double k, double sigmaValue, double mValue) ;

//pour les tests sur données de marché issues de VCUB
CheyetteDD_Model_PTR creeCheyetteDD_Modele_PTR(size_t xmax, CourbeInput_PTR pCourbeInput, 
											   double k, double sigmaValue, double mValue) ;
