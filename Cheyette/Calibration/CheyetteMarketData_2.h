#pragma once

#include <stdlib.h>
#include <vector>
#include <string>
#include <boost/shared_ptr.hpp>
#include <Instrument/VanillaSwaption.h>
#include <Cheyette/Model/CourbeInput.h>

#include <LMM/Helper/Name.h>
#include <LMM/Helper/LMMTenorStructure.h>

#include <LMM/UpperTriangleVanillaSwaptionQuotes.h>
#include <LMM/LmmSwaptionMarketData.h>
#include <LMM/helper/InputFileManager.h>

//structure du fichier 
//Shift_bp	Maturity	Tenor	Vol	Strike
//quotes ATM
//quotes + grand shift pour convexity
//quotes - grand shift pour convexity
//#

//cette classe est en complement de CheyetteMarketData
//recupere les donnees de marche pour calibrer la convexite du smile

class CheyetteMarketData_2	
{
private:
	//std::vector<size_t> a_expiry_ ;
	//std::vector<size_t> b_tenor_ ;

	//std::vector<double> strikes_ ;

	double shiftConvexity_ ;

	std::vector<VanillaSwaption_PTR>  vect_swaptions_convexity_pp_ ;
	std::vector<VanillaSwaption_PTR>  vect_swaptions_convexity_mm_ ;  //pour la convexity 

	std::vector<double>  vect_vol_convexity_pp_ ;  //pour la convexity 
	std::vector<double>  vect_vol_convexity_mm_ ;  //pour la convexity 

	Tenor floatTenor_ ;
	Tenor fixedTenor_ ;

public:
	//faire commencer les vecteurs de swaption par une swaption nulle
	CheyetteMarketData_2() 
		: shiftConvexity_(-1000), floatTenor_(Tenor::_6M), fixedTenor_(Tenor::_12M)
	{}

	//std::vector<size_t> get_aExpiry()	const {return a_expiry_ ; }		//maturity swaption aY bY
	//std::vector<size_t> get_bTenor()	const {return b_tenor_ ; }		//tenor swpation

	//std::vector<double> getStrikes()	const {return strikes_ ; } 

	double getShiftConvexity()	const {return shiftConvexity_ ; }	

	std::vector<VanillaSwaption_PTR>  getVect_swaptions_convexity_pp_() const{return vect_swaptions_convexity_pp_ ;}
	std::vector<VanillaSwaption_PTR>  getVect_swaptions_convexity_mm_() const{return vect_swaptions_convexity_mm_ ;}
	
	std::vector<double> getVect_vol_convexity_pp()	const {return vect_vol_convexity_pp_ ; }	
	std::vector<double> getVect_vol_convexity_mm()	const {return vect_vol_convexity_mm_ ; }	

private :
	void readFile(size_t fileNumber) ;
	void remplirCheyetteMarketData(const std::string& chaine) ;
	std::string parcoursChaine(const std::string& chaine, size_t indexStart) ;

public :
	void init(size_t fileNumber) ;

};

typedef boost::shared_ptr<CheyetteMarketData_2>       CheyetteMarketData_2_PTR;
typedef boost::shared_ptr<const CheyetteMarketData_2> CheyetteMarketData_2_CONSTPTR;

