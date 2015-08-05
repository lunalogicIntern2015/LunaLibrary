//#pragma once
//#include <vector>
//#include <boost/numeric/ublas/matrix.hpp>
//#include <cmath>
//
//#include <JBLMM/Element/Rate.h>
//#include <JBLMM/Element/ZC.h>
//#include <JBLMM/Instrument/BermudanSwaption.h>
//#include <JBLMM/Longstaff_Schwartz/RegressionLS.h>
//#include <JBLMM/Longstaff_Schwartz/Basis.h>
//#include <JBLMM/Longstaff_Schwartz/ExplanatoryVariable.h>
//#include <JBLMM/Longstaff_Schwartz/EV_Basis_Collection.h>
//
//#include <LMM/LmmModel/McLmm.h>
//#include <LMM/LmmModel/Shifted_HGVolatilityFunction.h>
//#include <LMM/pricer/LmmVanillaSwapPricer.h>
//
//namespace ublas = boost::numeric::ublas; 
//typedef boost::numeric::ublas::matrix<double> matrix;
//
// //*  !!the la last exercice time is T_N!!
// //*     T[i]				0     ...   T1   T2   T3    ...                   T_N			\\ nbDates  = T_N
// //* timeline				*---*---*---*----*----*---*----*----*---*----*----*			\\ nbExercice = N
// //* intrinsic Value					T1								T_{N-1}				\\ vi.size = N-1
// //* constiuation Value					T1								T_{N-1}				\\ vc.size = N-1
// //*vector indice						0								 N-2  N-1
//
//
//class LS
//{
//	McLmm_PTR mclmm_;  
//	//std::vector<Basis> basis_;
//	//ExplanatoryVariable ev_;
//	
//public:
//	//constructor
//	LS(McLmm_PTR mclmm):mclmm_(mclmm){}
//
//	//destructor
//	virtual ~LS(){}
//
//	//gettor
//	McLmm_PTR getMcLmm_PTR()const{return mclmm_;}
//	
//	double simulate(	BermudanSwaption_CONSTPTR bermudanSwaption_PTR, 
//						const std::vector<Basis_CONSTPTR>& basisVect,
//						RegressionLS_PTR rgLS,
//						const size_t nbSimu1,
//						const size_t nbSimu2);
//
//
//	//1er simulation: backward et regression
//	void simulateBackward(	BermudanSwaption_CONSTPTR bermudanSwaption_PTR, 
//							const std::vector<Basis_CONSTPTR>& basisVect,
//							const size_t nbSimu1,
//							RegressionLS_PTR rgLS);
//
//
//
//	
//	//2nd simulation: Pricing forward
//	double simulateForward(	BermudanSwaption_CONSTPTR bermudanSwaption_PTR, 
//							const std::vector<Basis_CONSTPTR>& basisVect,
//							const size_t nbSimu2, 
//							const std::vector<ublas::vector<double>>& param);
//
//
//private:
//
//	//suppose the last exercie index is vanillaSwap's endIndex
//	void oneTrajectoryBackward(	BermudanSwaption_CONSTPTR bermudanSwaption_PTR,
//								const std::vector<Basis_CONSTPTR>& basisVect,
//								std::vector<double>& intrinsicValueVectorT,
//								std::vector<double>& continuationValueVectorT, 
//								std::vector<double>& priceVectorT,
//								std::vector<double>& explanaryVariablesVector)const;
//
//	//To improuve
//	double iv(	BermudanSwaption_CONSTPTR bermudanSwaption_PTR,
//				const LMM::Index evaluationIndex)const;
//
//
//	double cv(	BermudanSwaption_CONSTPTR bermudanSwaption_PTR,
//				const size_t exerciceIndex,
//				const std::vector<Basis_CONSTPTR>& basisVect,
//				const std::vector<ublas::vector<double>>& param)const;
//};
//
//typedef boost::shared_ptr<LS> LS_PTR;
//typedef boost::shared_ptr<const LS> LS_CONSTPTR;
//
//
