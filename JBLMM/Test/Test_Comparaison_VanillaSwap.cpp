#include <JBLMM/Test/JBTests.h>

#include <vector>

#include <LMM/helper/LMMTenorStructure.h>

#include <LMM/pricer/LmmVanillaSwapPricer.h>
#include <JBLMM/Instrument/GenericSwap.h>
#include <JBLMM/Pricer/GenericVanillaSwapPricer.h>
#include <JBLMM/Element/Coupon.h>
#include <JBLMM/Element/CappedFlooredCoupon.h>
#include <JBLMM/Instrument/InstrumentFactory.h>
#include <JBLMM/Element/CouponLeg.h>



void vanillaSwapComparaisonExemple()
{

	double strike		=	0.02;
	size_t indexStart	=	3;
	size_t indexEnd		=	19;
	Tenor tenorStruture	=	Tenor::_6M;
	Tenor floatingTenor	=	Tenor::_6M;
	Tenor fixedTenor	=	Tenor::_12M;
	LMMTenorStructure_PTR lmmTenorStructure(new LMMTenorStructure(tenorStruture, 10));
	size_t liborIndex	=	lmmTenorStructure->get_horizon()+1;
	std::vector<double> myInitialLibor(liborIndex);

	for (size_t i = 0; i <myInitialLibor.size(); i++)
	{
		myInitialLibor[i]=0.02;							//+((double)i)*0.01;
	}

	cout	<<	"strike:                        "	<<	strike												<<	endl;
	cout	<<	"indexStart:                    "	<<	indexStart											<<	endl;
	cout	<<	"indexEnd:                      "	<<	indexEnd											<<	endl;
	cout	<<	"tenorStrutureYearFraction:     "	<<	lmmTenorStructure->get_tenorType().YearFraction()	<<	endl;
	cout	<<	"floatingVStenorStrutureRatio:  "	<<	floatingTenor.ratioTo(tenorStruture)				<<	endl;
	cout	<<	"fixedVStenorStrutureRatio:     "	<<	fixedTenor.ratioTo(tenorStruture)					<<	endl;
	cout	<<	"myInitialLibor:  ";
	for (size_t i = 0; i <myInitialLibor.size(); i++)
	{
		cout	<<	myInitialLibor[i]	<<	" ";
	}
	cout	<<	 endl;
	cout	<<	 endl;

	//VanillaSwap_Chi_Trang
	VanillaSwap myVS(strike, indexStart , indexEnd, floatingTenor, fixedTenor, lmmTenorStructure);
	LmmVanillaSwapPricer myVSP(lmmTenorStructure);
	double prix_swap=myVSP.swapNPV_Analytical_1(myVS, myInitialLibor);


	//GeneticSwap test
	//build geneticVanillaSwap
	GenericSwap_CONSTPTR vanillaSwap_Genetic=InstrumentFactory::createVanillaSwap(
		strike,indexStart,indexEnd,floatingTenor,fixedTenor,lmmTenorStructure,1.0);

	GenericVanillaSwapPricer_PTR geneticVanillaSwapPricer(new GenericVanillaSwapPricer());
	double geneticPrice=geneticVanillaSwapPricer->genericVanillaSwap_Analytical(vanillaSwap_Genetic, myInitialLibor);

	cout << "FirstVersionSwapPrice: "	<<	prix_swap				<< endl;
	cout << "GeneticSwapTest:       "	<<	geneticPrice			<< endl;
	cout << "Difference:            "	<<	geneticPrice-prix_swap	<< endl;
}
