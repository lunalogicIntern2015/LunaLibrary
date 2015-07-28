#pragma once

//pour LMMTenorStructure
const int max_nbOfYear = 40 ;

//ici l'indice 0 des vecteurs n'est pas utilisé
//pour avoir : vectExpiry[1] = 1Y, ...
//vectSwaption[1] = swaption 1Y (coterminal - 1)Y, la 1ere swaption de calibration

//longueur des vecteurs : coterminal
//il y a (coterminal - 1) swaptions de calibration

inline std::vector<size_t> getVectorExpiry(size_t coterminal)		//contient des dates, pas des index
{
	std::vector<size_t>	vectExpiry(coterminal) ;
	vectExpiry[0] = 0 ;
	for (size_t i = 1 ; i < coterminal ; ++i)
	{
		vectExpiry[i] = i ;
	}
	return vectExpiry ;
}

inline std::vector<size_t> getVectorTenor(size_t coterminal)		//contient des dates, pas des index
{
	std::vector<size_t>	vectTenor(coterminal) ;
	vectTenor[0] = 0 ;
	for (size_t i = 1 ; i < coterminal ; ++i)
	{
		vectTenor[i] = coterminal - i ;
	}
	return vectTenor ;
}

inline std::vector<double> getVectorStrike(LmmSwaptionMarketData_PTR pLmmSwaptionMarketData, size_t coterminal)
{
	std::vector<double>	vectStrike(coterminal) ;
	vectStrike[0] = 0 ;
	UpperTriangularDoubleMatrix strike = pLmmSwaptionMarketData->get_SwaptionQuotes_ATM()->get_UpperTriangularStrike() ;
	for (size_t i = 1 ; i < coterminal ; ++i)
	{
		vectStrike[i] = strike(i, coterminal - i) ;
	}
	return vectStrike ;
}


inline std::vector<double> getVectorVol(LmmSwaptionMarketData_PTR pLmmSwaptionMarketData, size_t coterminal)
{
	// 1st row and column not used! size = nbYear + 1
	size_t nbLignes		= pLmmSwaptionMarketData->get_SwaptionQuotes_ATM()->size1() ;	//nb year +1
	assert(coterminal <= nbLignes) ;

	std::vector<double>	vectVol(coterminal) ;
	vectVol[0] = 0 ;
	for (size_t i = 1 ; i < coterminal ; ++i)
	{
		UpperTriangularVanillaSwaptionQuotes quotes = 
					pLmmSwaptionMarketData->get_SwaptionQuotes_ATM()->get_UpperTriangularVanillaSwaptionQuotes() ;
		SwaptionQuote qu = quotes(i, coterminal - i) ;			//typedef std::pair<VanillaSwaption, double> SwaptionQuote; 
		vectVol[i] = qu.second ; 
	}
	return vectVol ;
}

inline std::vector<double> getVectorSkew(LmmSwaptionMarketData_PTR pLmmSwaptionMarketData, size_t coterminal)
{
	// 1st row and column not used! size = nbYear + 1
	size_t nbLignes		= pLmmSwaptionMarketData->get_SwaptionQuotes_skew()->size1() ;	//nb year +1
	assert(coterminal <= nbLignes) ;

	std::vector<double>	vectSkew(coterminal) ;
	vectSkew[0] = 0 ;
	for (size_t i = 1 ; i < coterminal ; ++i)
	{
		UpperTriangularVanillaSwaptionQuotes quotes = 
					pLmmSwaptionMarketData->get_SwaptionQuotes_skew()->get_UpperTriangularVanillaSwaptionQuotes() ;
		SwaptionQuote qu = quotes(i, coterminal - i) ;			//typedef std::pair<VanillaSwaption, double> SwaptionQuote; 
		vectSkew[i] = qu.second ; 
	}
	return vectSkew ;
}


//get_VectDiscountFactor() interpole les ZERO COUPON
//on a donc les yield correspondant aux ZC interpolés
//pas interpolation des yields
inline CourbeInput upperMatrixZC_To_CourbeInput_yield(LiborQuotes_PTR liborQuotes_PTR)
{
	std::vector<double> listeMatu(liborQuotes_PTR->get_LMMTenorStructure_PTR()->get_tenorDate()) ; 	
	std::vector<double> ZC(liborQuotes_PTR->get_VectDiscountFactor()) ;
	std::vector<double> tauxZC(ZC.size()) ;

	for (size_t i = 1 ; i < tauxZC.size() ; ++i)
	{
		tauxZC[i] = - 1./listeMatu[i] * log(ZC[i]) ;
	}
	tauxZC[0] = tauxZC[1] ;
	return CourbeInput(listeMatu, tauxZC);
}

inline std::vector<VanillaSwaption_PTR> getVectorSwaptions(const Tenor& tenorFloat, const Tenor& tenorFixed,
												const std::vector<size_t>& vectExpiry, 
												const std::vector<size_t>& vectTenor, 
												const std::vector<double>& vectStrike)
{
	std::vector<VanillaSwaption_PTR> vectSwaptions(vectExpiry.size()) ;
	VanillaSwaption_PTR pSwaptionInit(new VanillaSwaption()) ;			//swaption 0 non utilisée
	vectSwaptions[0] = pSwaptionInit ;

	double float_tenor = tenorFloat.YearFraction() ;
	double fixed_tenor = tenorFixed.YearFraction() ;
	double tenor_ref = std::min(fixed_tenor, float_tenor) ; 

	for (size_t i = 1 ; i < vectSwaptions.size() ; ++i)  
	{
		size_t indexStart		= size_t(vectExpiry[i] / tenor_ref) ;
		size_t indexEnd			= indexStart + size_t(vectTenor[i] / tenor_ref) ;
		
		LMMTenorStructure_PTR pTenorStructure(new LMMTenorStructure(tenorFloat, max_nbOfYear)) ;
		VanillaSwap swap = VanillaSwap(vectStrike[i], indexStart, indexEnd, tenorFloat, tenorFixed, pTenorStructure) ;
		VanillaSwaption_PTR pSwaption(new VanillaSwaption(swap, OptionType::OptionType::CALL)) ;
		vectSwaptions[i] = pSwaption ;				
	}
	return vectSwaptions ;
}
