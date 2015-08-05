//#include <JBLMM/Longstaff_Schwartz/EV_Basis_Collection.h>
//#include <vector>
//
//
//void EV_Basis_Collection::evaluate(	McLmm_PTR mclmm,
//									BermudanSwaption_CONSTPTR bermudanSwaption_PTR,
//									LMM::Index liborIndex)const
//{
//	evValueVect_.clear();
//	size_t nbEv = ev_vect_.size();
//	for(size_t indexEv = 0; indexEv<nbEv; indexEv++ )
//	{
//		ExplanatoryVariable_CONSTPTR ev_ptr = getEvVect()[indexEv];
//		double evValue = ev_ptr->getEvaluator_PTR()->evaluate(	mclmm,
//																bermudanSwaption_PTR,
//																liborIndex);
//		evValueVect_.push_back(evValue);
//	}
//}
//
//void EV_Basis_Collection::addBasisValueToVector(std::vector<double>& basisValueVect)const
//{
//	assert(!evValueVect_.empty());
//	size_t maxNbEv = evValueVect_.size();
//	size_t nbEvBasisPair = evBasisPairvect_.size();
//	for(size_t indexPair = 0; indexPair<nbEvBasisPair;indexPair++ )
//	{
//		Basis_CONSTPTR basis_ptr = evBasisPairvect_[indexPair].second;
//		std::vector<double> evValueVect;
//		getEvValueForBasis(evBasisPairvect_[indexPair].first, evValueVect);
//		double basisValue = basis_ptr->transform(evValueVect[0]);
//		basisValueVect.push_back(basisValue);
//	}
//}
//
//void EV_Basis_Collection::getEvValueForBasis(	const std::vector<size_t>& evVect, 
//												std::vector<double>& evValueVect)const
//{
//	assert(!evValueVect_.empty());
//	size_t maxNbEv = evValueVect_.size();
//	size_t nbEv = evVect.size();
//
//	for(size_t i = 0; i<nbEv; i++)
//	{
//		size_t evValueIndex = evVect[i];
//		assert(evValueIndex<maxNbEv);
//
//		double evValue = evValueVect_[evValueIndex];
//		evValueVect.push_back(evValue);
//	}
//}