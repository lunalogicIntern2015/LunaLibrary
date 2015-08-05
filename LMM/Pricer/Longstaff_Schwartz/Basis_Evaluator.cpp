#include <LMM/Pricer/Longstaff_Schwartz/Basis_Evaluator.h>


double Basis_SinglaEVFunctional_Evaluator::evaluate(Basis_CONSTPTR basis, 
													const matrix& liborMatrix, 
													const std::vector<double>& numeraire)const
{
	Basis_SinglaEVFunctional_CONSTPTR basis_SingleEVFunctional = boost::dynamic_pointer_cast<const Basis_SinglaEVFunctional>(basis);
	if(basis_SingleEVFunctional)
	{
		double ev_val = ev_evaluator_->evaluate(basis_SingleEVFunctional->getEV(),liborMatrix,numeraire);
		return basis_SingleEVFunctional->getTransformer()->transform(ev_val);
	}
	else
	{
		throw("Basis_SinglaEVFunctional_Evaluator mismatch the basis");
	}
};

void Basis_SinglaEVFunctional_Evaluator::write_to_stream(std::ostream& out)const
{
	out	<<	"---Basis_SinglaEVFunctional_Evaluator---"	<<	std::endl;
	ev_evaluator_->write_to_stream(out);
}

Basis_Composited_Evaluator::Basis_Composited_Evaluator(	Basis_Evaluator_CONSTPTR basis_Evaluator1, 
														Basis_Evaluator_CONSTPTR basis_Evaluator2)
	:
		basis_Evaluator1_(basis_Evaluator1),
		basis_Evaluator2_(basis_Evaluator2)
{
}

double Basis_Composited_Evaluator::evaluate(Basis_CONSTPTR basis, 
											const matrix& liborMatrix, 
											const std::vector<double>& numeraire)const
{
	Basis_Composited_CONSTPTR basis_Composited = boost::dynamic_pointer_cast<const Basis_Composited>(basis);

	if(basis_Composited)
	{
		double basis_val1 = basis_Evaluator1_->evaluate(basis_Composited->getBasis1(),liborMatrix,numeraire);
		double basis_val2 = basis_Evaluator2_->evaluate(basis_Composited->getBasis2(),liborMatrix,numeraire);
		return basis_Composited->getTransformer()->transform(basis_val1,basis_val2);
	}
	else
	{
		throw("Basis_Composited_Evaluator mismatch the basis");
	}
}



void Basis_Composited_Evaluator::write_to_stream(std::ostream& out)const
{
	out	<<	"---Basis_Composited_Evaluator---"	<<	std::endl;
	basis_Evaluator1_->write_to_stream(out);
	basis_Evaluator2_->write_to_stream(out);
}


Basis_Composited_Function_Evaluator::Basis_Composited_Function_Evaluator(	Basis_Evaluator_CONSTPTR basis_Evaluator1, 
																			Basis_Evaluator_CONSTPTR basis_Evaluator2)
	:
		Basis_Composited_Evaluator(basis_Evaluator1,basis_Evaluator2)
{
}

void Basis_Composited_Function_Evaluator::write_to_stream(std::ostream& out)const
{
	out	<<	"---Basis_Composited_Function_Evaluator---"	<<	std::endl;
	basis_Evaluator1_->write_to_stream(out);
	basis_Evaluator2_->write_to_stream(out);
}

double Basis_Composited_Function_Evaluator::evaluate(	Basis_CONSTPTR basis, 
														const matrix& liborMatrix, 
														const std::vector<double>& numeraire)const
{
	Basis_Composited_CONSTPTR basis_Composited = boost::dynamic_pointer_cast<const Basis_Composited>(basis);
	Basis_SinglaEVFunctional_CONSTPTR singleEV_basis2 =boost::dynamic_pointer_cast<const Basis_SinglaEVFunctional>(basis_Composited->getBasis2());

	if(basis_Composited&&singleEV_basis2)
	{
		double basis_val1 = basis_Evaluator1_->evaluate(basis_Composited->getBasis1(),liborMatrix,numeraire);		
		double composedValue = singleEV_basis2->getTransformer()->transform(basis_val1);
		return composedValue;
	}
	else
	{
		throw("Basis_Composited_Function_Evaluator mismatch the basis");
	}
}

Basis_Polynomial_Evaluator::Basis_Polynomial_Evaluator(const std::vector<Basis_Evaluator_CONSTPTR>& basis_evaluator_vect)
	:
	basis_evaluator_vect_(basis_evaluator_vect)
{
}

double Basis_Polynomial_Evaluator::evaluate(Basis_CONSTPTR basis, const matrix& liborMatrix, const std::vector<double>& numeraire)const
{
	Basis_Polynomial_CONSTPTR basis_polynomial = boost::dynamic_pointer_cast<const Basis_Polynomial>(basis);
	std::vector<Basis_CONSTPTR> basis_vect = basis_polynomial->get_Basis_Vect();

	assert(basis_vect.size()==basis_evaluator_vect_.size());
	
	std::vector<double> val_vect(basis_vect.size(),1.0);
	for(size_t i = 0; i<basis_vect.size(); i++)
		val_vect[i] = basis_evaluator_vect_[i]->evaluate(basis_vect[i],liborMatrix,numeraire);

	double result = basis_polynomial->getTransformer()->transform(val_vect);
	return result;
}

void Basis_Polynomial_Evaluator::write_to_stream(std::ostream& out)const
{
	out	<<	"---Basis_Polynomial_Evaluator---"	<<	std::endl;
	for(size_t i = 0; i<basis_evaluator_vect_.size(); i++)
		basis_evaluator_vect_[i]->write_to_stream(out);
}

