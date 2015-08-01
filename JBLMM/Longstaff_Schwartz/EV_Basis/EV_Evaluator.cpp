#include <JBLMM/Longstaff_Schwartz/EV_Basis/EV_Evaluator.h>

double EV_ConstRate_Evaluator::evaluate(EV_CONSTPTR ev, 
										const matrix& liborMatrix, 
										const std::vector<double>& numeraire)const
{
	ConstRate_CONSTPTR constRate = boost::dynamic_pointer_cast<const ConstRate>(ev->getRate());

	if(constRate)
	{
		double constRateValue = constRate->getConstRateValue();
		return constRateValue;
	}
	else
	{
		throw("EV_ConstRate_Evaluator don't match the rate");
	}
}

void EV_ConstRate_Evaluator::write_to_stream(std::ostream& out)const
{
	out	<<	"---EV_ConstRate_Evaluator---"	<<	std::endl;
}

double EV_LiborRate_Evaluator::evaluate(EV_CONSTPTR ev, 
										const matrix& liborMatrix, 
										const std::vector<double>& numeraire)const
{
	LiborRate_CONSTPTR liborRate = boost::dynamic_pointer_cast<const LiborRate>(ev->getRate());

	if(liborRate)
	{
		LMM::Index fixingIndex = liborRate->getFixingTime();
		double libor = liborMatrix(fixingIndex,fixingIndex);
		return libor;
	}
	else
	{
		throw("EV_LiborRate_Evaluator don't match the rate");
	}
}

void EV_LiborRate_Evaluator::write_to_stream(std::ostream& out)const
{
	out	<<	"---EV_LiborRate_Evaluator---"	<<	std::endl;
}


EV_VanillaSwapRate_Evaluator::EV_VanillaSwapRate_Evaluator(const McLmmVanillaSwapPricer& mcLmmVanillaSwapPricer)
	:
	mcLmmVanillaSwapPricer_(mcLmmVanillaSwapPricer)
{
}

double EV_VanillaSwapRate_Evaluator::evaluate(EV_CONSTPTR ev, 
											  const matrix& liborMatrix, 
											  const std::vector<double>& numeraire)const
{
	VanillaSwapRate_CONSTPTR vanillaSwapRate = boost::dynamic_pointer_cast<const VanillaSwapRate>(ev->getRate());
	if(vanillaSwapRate)
	{
		const VanillaSwap& vanillaSwap = vanillaSwapRate->getVanillaSwap();
		double vanillaSwapRate = mcLmmVanillaSwapPricer_.swapRate(vanillaSwap.get_indexStart(),vanillaSwap,numeraire,liborMatrix);
		return vanillaSwapRate;
	}
	else
	{
		throw("EV_VanillaSwapRate_Evaluator don't match the rate");
	}
}

void EV_VanillaSwapRate_Evaluator::write_to_stream(std::ostream& out)const
{
	out	<<	"---EV_VanillaSwapRate_Evaluator---"	<<	std::endl;
}

