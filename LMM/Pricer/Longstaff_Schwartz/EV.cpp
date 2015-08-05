#include <LMM/Pricer/Longstaff_Schwartz/EV.h>

//void EV::write_to_stream(std::ostream& out)const
//{
//}

void EV_ConstRate::write_to_stream(std::ostream& out)const
{
	out	<<	"---EV_ConstRate---"	<<	std::endl;
	rate_->write_to_stream(out);
}

void EV_LiborRate::write_to_stream(std::ostream& out)const
{
	out	<<	"---EV_LiborRate---"	<<	std::endl;
	rate_->write_to_stream(out);
}

void EV_VanillaSwapRate::write_to_stream(std::ostream& out)const
{
	out	<<	"---EV_VanillaSwapRate---"	<<	std::endl;
	rate_->write_to_stream(out);
}