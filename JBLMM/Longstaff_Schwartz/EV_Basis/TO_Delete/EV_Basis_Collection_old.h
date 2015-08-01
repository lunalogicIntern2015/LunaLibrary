//#pragma once
//#include <vector>
//#include <LMM/LmmModel/McLmm.h>
//#include <JBLMM/Element/Rate.h>
//#include <JBLMM/Instrument/BermudanSwaption.h>
//#include <JBLMM/Longstaff_Schwartz/ExplanatoryVariable.h>
//#include <JBLMM/Longstaff_Schwartz/Basis.h>
//
//typedef std::pair<std::vector<size_t>, Basis_CONSTPTR> EV_Basis_pair;
//class EV_Basis_Collection
//{
//
//	std::vector<ExplanatoryVariable_CONSTPTR> ev_vect_; //all EV
//	std::vector<EV_Basis_pair> evBasisPairvect_;		//all EV-Basis pair
//
//	mutable std::vector<double> evValueVect_;		//save all EV values
//	//mutable std::vector<double> basisValueVect_;	//save all Basis values
//
//public:
//	//construtor et destructor
//	EV_Basis_Collection(const std::vector<ExplanatoryVariable_CONSTPTR>& ev_vect,
//						const std::vector<EV_Basis_pair>& evBasisPairvect)
//						:
//						ev_vect_(ev_vect), evBasisPairvect_(evBasisPairvect){}
//	virtual ~EV_Basis_Collection(){};
//	//getters
//	const std::vector<ExplanatoryVariable_CONSTPTR>& getEvVect()const{return ev_vect_;}
//	const std::vector<EV_Basis_pair>& getEvBasisPairVect() const {return evBasisPairvect_;}
//
//	const size_t getMaxNbEV() const {return ev_vect_.size();}
//	const size_t getNbPair() const {return evBasisPairvect_.size();}
//	
//	
//	//evaluate ev values
//	void evaluate(	McLmm_PTR mclmm,
//					BermudanSwaption_CONSTPTR bermudanSwaption_PTR,
//					LMM::Index liborIndex)const;
//
//	//get basis values and add them to basisValueVect
//	void addBasisValueToVector(std::vector<double>& valueVect)const;
//
//private:
//	void getEvValueForBasis(	const std::vector<size_t>& evVect, 
//								std::vector<double>& basisValueVect)const;
//
//};
//typedef boost::shared_ptr<EV_Basis_Collection> EV_Basis_Collection_PTR;
//typedef boost::shared_ptr<const EV_Basis_Collection> EV_Basis_Collection_CONSTPTR;
//
