#include <JBLMM/Longstaff_Schwartz/Regression/Regression_LS.h>


LS::RegressionRepresentation::RegressionRepresentation(	const std::vector<Basis_CONSTPTR>& basis, 
													const std::vector<Basis_Evaluator_CONSTPTR>& basisEvaluators)
												:
													basis_(basis),
													basisEvaluators_(basisEvaluators), 
													basis_val_buffer_(basis.size()),
													regressionCoeffs_(basis.size())
{
	if(basis.size() != basisEvaluators.size())
		throw ("RegressionResult constructor size mismatch!");
}


void LS::RegressionRepresentation::evaluate_basis(const McLmm_LS::LMMSimulationResult& lmmSimulationResult)const
{
	evaluate_basis(lmmSimulationResult.get_liborMatrix(), lmmSimulationResult.get_numeraire());
}

void LS::RegressionRepresentation::evaluate_basis(const matrix& liborMatrix, const std::vector<double>& numeraire)const 
{
	//std::vector<double> basis_vals(basis_.size());
	for(size_t i=0; i<basis_.size(); ++i)
	{
		basis_val_buffer_[i] = basisEvaluators_[i]->evaluate(basis_[i], liborMatrix, numeraire);  // TODO: YY may not be efficient !!! reevaluate EV !!! 
	}
}

//!! continuous value
double LS::RegressionRepresentation::evaluate_representation(const McLmm_LS::LMMSimulationResult& lmmSimulationResult)const // conditional expectation 
{
	evaluate_basis(lmmSimulationResult);
	double result = 0.0;
	for(size_t i=0; i<basis_.size(); ++i)
	{
		result += basis_val_buffer_[i]*regressionCoeffs_[i];
	}
	return result; 
}

void LS::RegressionRepresentation::write_to_stream(std::ostream& out)const
{
	out	<<	"---RegressionRepresentation---"	<<	std::endl;

	out << "basis_val_buffer_: ;" ;
	for(size_t i = 0; i<basis_val_buffer_.size(); i++)
		out << basis_val_buffer_[i] << " ;" ;
	out << endl;

	out << "basis_vector: ;" << endl;
	for(size_t i = 0; i<basis_.size(); i++)
	{
		out << "basis " << i << endl ;
		basis_[i]->write_to_stream(out);
	}
	out << endl;

	out << "basisEvaluators__vector: ;" << endl;
	for(size_t i = 0; i<basisEvaluators_.size(); i++)
	{
		out << "basisEvaluators " << i << endl ;
		basisEvaluators_[i]->write_to_stream(out);
	}
	out << endl;

	out << "regressionCoeffs_: ;" ;
	for(size_t i = 0; i<regressionCoeffs_.size(); i++)
		out << regressionCoeffs_[i] << " ;" ;
	out << endl;
}

LS::Regression::Regression(const RegressionRepresentation& regressionRepresentation) 
	: 
	regressionRepresentation_(regressionRepresentation)
{
}

void LS::Regression::regress(const std::vector<double>& Z, 
						     const std::vector<McLmm_LS::LMMSimulationResult>& lmmSimulationResults,
							 std::vector<std::vector<double>>& basis_value_on_allPath_buffer)
{
	
	constructBasisMatrix(lmmSimulationResults,basis_value_on_allPath_buffer);

	std::vector<double> coef_weight;

	normalizationMatrix(basis_value_on_allPath_buffer, coef_weight);
	
	//---------------------change vector of vector to Rmatrix-----------------------------------------------//
	size_t nbCol = basis_value_on_allPath_buffer[0].size();
	size_t nbRow= basis_value_on_allPath_buffer.size();
	Rmatrix basis_value_on_allPath_boostMatrix(nbRow,nbCol);
	for(size_t i=0; i<nbRow;i++)
	{
		for(size_t j=0; j<nbCol;j++)
			basis_value_on_allPath_boostMatrix(i,j)=basis_value_on_allPath_buffer[i][j];
	}
	Rvector Z_vector(Z.size());
	std::copy(Z.begin(), Z.end(),Z_vector.begin());
	//---------------------End of change vector of vector to Rmatrix---------------------------------------//

	// do matrix multiplication: Z regression on basis, result save to regressionRepresentation_'s coeff ! 

	//TODO: junbin : copy the matrix multiplication code here ... 
	Rmatrix inverse_TransposeX_multiple_X(nbCol,nbCol); 
	const Rmatrix& transposeX = ublas::trans(basis_value_on_allPath_boostMatrix);
	assert(InvertMatrix(ublas::prod(transposeX, basis_value_on_allPath_boostMatrix), inverse_TransposeX_multiple_X));
	const ublas::vector<double>& matrix_transposeX_multiple_Y = ublas::prod(transposeX,Z_vector);
	const ublas::vector<double>& rg_coef_vect = ublas::prod(inverse_TransposeX_multiple_X, matrix_transposeX_multiple_Y);


	size_t nbRgCoef = rg_coef_vect.size();
	if(regressionRepresentation_.getRegressionCoeffs().size()!=nbRgCoef)
		regressionRepresentation_.getRegressionCoeffs().resize(nbRgCoef);
	assert(coef_weight.size()==nbRgCoef);
	for(size_t i = 0; i<nbRgCoef; i++)
		regressionRepresentation_.getRegressionCoeffs()[i]=rg_coef_vect[i]/coef_weight[i];
}

void LS::Regression::constructBasisMatrix(	const std::vector<McLmm_LS::LMMSimulationResult>& lmmSimulationResults, 
											std::vector<std::vector<double>>& basis_value_on_allPath_buffer)  // Junbin: use a commun vector of vector to save basis value  
{

	size_t nbSimulation = lmmSimulationResults.size();
	if(basis_value_on_allPath_buffer.size()!=nbSimulation)
		basis_value_on_allPath_buffer.resize(nbSimulation); 
	size_t nbBasis = regressionRepresentation_.getBasis().size();
	
	//! construct Basis matrix 
	for(size_t path_index = 0; path_index < lmmSimulationResults.size(); ++path_index)
	{

		regressionRepresentation_.evaluate_basis(lmmSimulationResults[path_index]);
		const std::vector<double>& basis_value_vect =  regressionRepresentation_.getBasis_val_buffer();
		basis_value_on_allPath_buffer[path_index]=basis_value_vect; // YY: bad code, too many copy-coller !!!!!
		//const std::vector<double>& basis_value_vect =  regressionRepresentation_.getBasis_val_buffer();
		//basis_value_on_allPath_buffer_[path_index].resize(nbBasis);
		//basis_value_on_allPath_buffer_[path_index]=basis_value_vect; // YY: bad code, too many copy-coller !!!!!
		//clock_t copy_basis_one_line_endTime = clock();
	}
}

bool LS::Regression::InvertMatrix (const ublas::matrix<double>& input, ublas::matrix<double>& inverse) 
{  
	typedef ublas::permutation_matrix<std::size_t> pmatrix; 
	// create a working copy of the input 
	ublas::matrix<double> A(input); 
	// create a permutation matrix for the LU-factorization 
	pmatrix pm(A.size1()); 
	// perform LU-factorization 
	int res = lu_factorize(A,pm); 
	if( res != 0 )
		return false; 
	// create identity matrix of "inverse" 
	inverse.assign(ublas::identity_matrix<double>(A.size1())); 
	// backsubstitute to get the inverse 
	lu_substitute(A, pm, inverse); 
	return true; 
}

void LS::Regression::write_to_stream(std::ostream& out)const
{
	out	<<	"---Regression---"	<<	std::endl;
	regressionRepresentation_.write_to_stream(out);
	out << endl;
}

//Junbin : normalize matrix and set column averge value at vector
void LS::Regression::normalizationMatrix(std::vector<std::vector<double>>& m , std::vector<double>& moyen_column)
{
	size_t nbRow = m.size();

	assert(nbRow>0);
	assert(m[0].size()>0);

	moyen_column.resize(m[0].size(),0.0);

	for(size_t colIndex=0; colIndex<m[0].size(); colIndex++)
	{
		//double squareSum=0.0;
		for(size_t rowIndex=0; rowIndex<nbRow;rowIndex++)
		{
			moyen_column[colIndex]+=m[rowIndex][colIndex];
			//squareSum+=m[rowIndex][colIndex]*m[rowIndex][colIndex];
		}
		moyen_column[colIndex] = moyen_column[colIndex]/nbRow;
		//double variance	= squareSum/nbRow - mean*mean;
		//double sd = sqrt(variance);
		for(size_t rowIndex=0; rowIndex<nbRow;rowIndex++)
		{
			m[rowIndex][colIndex]/=moyen_column[colIndex];
		}
	}
}