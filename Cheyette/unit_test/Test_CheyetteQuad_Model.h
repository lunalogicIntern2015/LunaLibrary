#pragma once

#include <Cheyette/Model/CheyetteQuad_Model.h>
#include <Cheyette/unit_test/Test_CheyetteDD_Model.h>   //pour createCourbeInput

CheyetteQuad_Model_PTR creeCheyetteQuad_Modele_PTR(size_t xmax, int numCourbe, 
												   double k, double aValue, double bValue, double cValue) ;

CheyetteQuad_Model_PTR creeCheyetteQuad_Modele_PTR(size_t xmax, CourbeInput_PTR pCourbeInput, 
												   double k, double aValue, double bValue, double cValue) ;

