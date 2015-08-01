//#include <JBLMM/Longstaff_Schwartz/LS.h>
//#include <LMM/pricer/McLmmVanillaSwapPricer.h>
//
//double LS::simulate(	BermudanSwaption_CONSTPTR bermudanSwaption_PTR, 
//						const std::vector<Basis_CONSTPTR>& basisVect,
//						RegressionLS_PTR rgLS,
//						const size_t nbSimu1,
//						const size_t nbSimu2)
//{
//
//	simulateBackward(	bermudanSwaption_PTR, basisVect, nbSimu1,rgLS);
//
//	double result = simulateForward(bermudanSwaption_PTR, basisVect, nbSimu2, rgLS->getParam());
//
//	return result;
//}
//
//void LS::simulateBackward(	BermudanSwaption_CONSTPTR bermudanSwaption_PTR, 
//							const std::vector<Basis_CONSTPTR>& basisVect,
//							const size_t nbSimu,
//							RegressionLS_PTR rgLS)
//{
//	const std::vector<LMM::Index>& exerciceIndexVector = bermudanSwaption_PTR->getExerciseTimes();
//	size_t nbExercice = exerciceIndexVector.size();
//	size_t regressionNbExercice = nbExercice-1;
//
//	size_t nbBasis = basisVect.size();
//
//	rgLS->setNbBasis(nbBasis);
//	rgLS->setNbExercice(regressionNbExercice);
//	rgLS->setNbSimulation(nbSimu);
//	rgLS->getCV().clear();
//	rgLS->getExplanatoryData().clear();
//	rgLS->clearParam();
//
//	for(size_t i=0; i<nbSimu; i++)
//	{
//		//for each simulation
//		mclmm_->simulateLMM();
//		std::vector<double> intrinsicVectorT(regressionNbExercice);
//		std::vector<double> continuationVectorT(regressionNbExercice);
//		std::vector<double> priceVectorT(nbExercice);
//		
//		std::vector<double> explanaryVariablesVector;
//		const matrix& libormatrix = mclmm_->get_liborMatrix();
//		const std::vector<double>& numeraire = mclmm_->get_numeraire();
//		oneTrajectoryBackward(	bermudanSwaption_PTR, 
//								basisVect, 
//								intrinsicVectorT, 
//								continuationVectorT, 
//								priceVectorT,
//								explanaryVariablesVector);
//
//		rgLS->getExplanatoryData().push_back(explanaryVariablesVector);
//		rgLS->getCV().push_back(continuationVectorT);
//	}
//
//	rgLS->regAll();
//}
//
//double LS::simulateForward(	BermudanSwaption_CONSTPTR bermudanSwaption_PTR, 
//							const std::vector<Basis_CONSTPTR>& basisVect,
//							const size_t nbSimu, 
//							const std::vector<ublas::vector<double>>& param)
//{	
//	std::vector<LMM::Index> exerciceVector=bermudanSwaption_PTR->getExerciseTimes();
//	size_t nbExercice=exerciceVector.size();
//	size_t regressionNbExercice = nbExercice-1;
//
//	double result = 0.0;
//	double variance = 0.0;
//	for(size_t i=0; i<nbSimu; i++)
//	{
//		mclmm_->simulateLMM();
//		const ublas::matrix<double>& libormatrix=mclmm_->get_liborMatrix();
//		//terminal numeraire
//		const std::vector<double>& numeraire=mclmm_->get_numeraire();
//		double oneIterationPrice = 0.0;
//		for(size_t exerciceIndex=0; exerciceIndex<regressionNbExercice; exerciceIndex++)
//		{
//			LMM::Index liborExerciceIndex = exerciceVector[exerciceIndex];
//			double intrinsicValue		=	iv(bermudanSwaption_PTR,liborExerciceIndex);
//			double continuationValue	=	cv(bermudanSwaption_PTR,exerciceIndex,basisVect,param);
//			//if iv>=cv, we exercice the option 
//			if(intrinsicValue>=continuationValue)															
//			{
//				oneIterationPrice=numeraire[0]/numeraire[liborExerciceIndex]*intrinsicValue;			
//				break;
//			}
//		}
//		result+=oneIterationPrice;
//		variance+=oneIterationPrice*oneIterationPrice;
//	}
//	result /= nbSimu;
//	return result;
//}
//
//
//void LS::oneTrajectoryBackward(	BermudanSwaption_CONSTPTR bermudanSwaption_PTR,
//								const std::vector<Basis_CONSTPTR>& basisVect,
//								std::vector<double>& intrinsicValueVectorT,
//								std::vector<double>& continuationValueVectorT, 
//								std::vector<double>& priceVectorT,
//								std::vector<double>& explanaryVariablesVector)const
//{
//	const std::vector<LMM::Index>&  exerciceIndexVect = bermudanSwaption_PTR->getExerciseTimes();
//
//	//check if the last exercice index is the last vanillaswap index
//	LMM::Index endIndex = bermudanSwaption_PTR->getVanillaSwap_PTR()->get_indexEnd();	
//	assert(exerciceIndexVect[exerciceIndexVect.size()-1]==endIndex);
//	assert(continuationValueVectorT.size()==exerciceIndexVect.size()-1);
//	assert(intrinsicValueVectorT.size()==exerciceIndexVect.size()-1);
//	assert(priceVectorT.size()==exerciceIndexVect.size());
//
//	size_t nbExercice = exerciceIndexVect.size();
//	size_t regressionNbExercice = nbExercice-1;
//
//	//save explanaryValues
//	for(size_t exerciceIndex=0; exerciceIndex<regressionNbExercice; exerciceIndex++)
//	{
//		LMM::Index liborExerciceIndex = exerciceIndexVect[exerciceIndex];
//		for(size_t basisIndex=0; basisIndex<basisVect.size(); basisIndex++)
//		{
//			double basisvalue=basisVect[basisIndex]->evaluateEV(mclmm_,bermudanSwaption_PTR, liborExerciceIndex);	
//			double basisValue=basisVect[basisIndex]->transform(basisvalue);
//			explanaryVariablesVector.push_back(basisValue);
//		}
//	}
//
//	//calculate continuation value at each exercice time
//	//backward regression
//	priceVectorT[nbExercice-1] = 0.0;
//	for(int indexExercice=regressionNbExercice-1; indexExercice>=0; --indexExercice)
//	{
//		LMM::Index liborExerciceIndex = exerciceIndexVect[indexExercice];	
//		LMM::Index nextExerciceLiborIndex = exerciceIndexVect[indexExercice+1];
//		const std::vector<double>& numeraire=mclmm_->get_numeraire();
//		double factorAtualization=numeraire[liborExerciceIndex]/numeraire[nextExerciceLiborIndex];
//		continuationValueVectorT[indexExercice]=factorAtualization*priceVectorT[indexExercice+1];
//		intrinsicValueVectorT[indexExercice] = iv(bermudanSwaption_PTR,liborExerciceIndex);
//		priceVectorT[indexExercice] = max(	intrinsicValueVectorT[indexExercice], 
//											continuationValueVectorT[indexExercice]);
//	}
//}
//
//double LS::iv(	BermudanSwaption_CONSTPTR bermudanSwaption_PTR,
//				const LMM::Index evaluationIndex)const
//{
//
//		LMM::Index endIndex = bermudanSwaption_PTR->getVanillaSwap_PTR()->get_indexEnd();
//		VanillaSwap_CONSTPTR subVanillaSwap_PTR = bermudanSwaption_PTR->getSubVanillaSwap(evaluationIndex, endIndex);
//
//		McLmmVanillaSwapPricer_PTR mcLmmVanillaSwapPricer_PTR(new McLmmVanillaSwapPricer(mclmm_));
//		const matrix& libormatrix = mclmm_->get_liborMatrix();
//		const std::vector<double>& numeraire = mclmm_->get_numeraire();
//
//		double floatingPV = mcLmmVanillaSwapPricer_PTR->pvFloatingLeg(	evaluationIndex, 
//																		*subVanillaSwap_PTR.get(), 
//																		numeraire,
//																		libormatrix);
//
//		double fixedPV = mcLmmVanillaSwapPricer_PTR->pvFixedLeg(evaluationIndex, 
//																*subVanillaSwap_PTR.get(),
//																numeraire);
//
//		double swapPriceAtEvaluationDay		=	floatingPV-fixedPV;
//
//		return swapPriceAtEvaluationDay;
//}
//
////evaluate continuation value for pricing forward
//double LS::cv(	BermudanSwaption_CONSTPTR bermudanSwaption_PTR,
//				const size_t exerciceIndex,
//				const std::vector<Basis_CONSTPTR>& basisVect,
//				const std::vector<ublas::vector<double>>& param)const
//{
//	size_t nbBasis = basisVect.size();
//	double res=0.0;
//	size_t basisCounter=0;
//	const std::vector<LMM::Index>& exerciceVect=bermudanSwaption_PTR->getExerciseTimes();
//	LMM::Index liborExerciceIndex = exerciceVect[exerciceIndex];
//	//evaluate ev values
//	for(size_t basisIndex=0; basisIndex<basisVect.size(); basisIndex++)
//	{
//
//	}
//
//	for(size_t basisIndex=0; basisIndex<nbBasis; basisIndex++)
//	{
//		double evvalue=basisVect[basisIndex]->evaluateEV(mclmm_,bermudanSwaption_PTR, liborExerciceIndex);	
//		double basisValue=basisVect[basisIndex]->transform(evvalue);
//		double coefValue = param[exerciceIndex][basisCounter];	
//
//		double elementValue = coefValue*basisValue;
//		res+=elementValue;
//		basisCounter++;
//	}
//	return res;
//}
