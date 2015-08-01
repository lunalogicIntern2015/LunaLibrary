//#include <JBLMM/Longstaff_Schwartz/RegressionLS.h>
//
//RegressionLS::RegressionLS(	size_t nbExercice,							
//							size_t nbBasis,
//							size_t nbSimu,
//							std::vector<std::vector<double>>& explanatoryData,
//							std::vector<std::vector<double>>& continuationValueData)
//	:
//		nbExercice_(nbExercice),
//		nbBasis_(nbBasis),
//		nbSimu_(nbSimu),
//		explanatoryData_(explanatoryData),
//		cv_(continuationValueData)
//{
//	assert(explanatoryData_.size()==nbSimu&&cv_.size()==nbSimu);
//}
//
//void RegressionLS::reg_OneExerciceTime(const ublas::vector<double>& cvVector, 
//									   const Rmatrix& explicatif)
//{
//	Rmatrix inverseMatrix(explicatif.size2(),explicatif.size2()); 
//	const Rmatrix& transposeX = ublas::trans(explicatif);
//	assert(InvertMatrix(ublas::prod(transposeX, explicatif), inverseMatrix));
//	const ublas::vector<double>& matrixX_tY = ublas::prod(transposeX,cvVector);
//	ublas::vector<double> vect = ublas::prod(inverseMatrix,matrixX_tY);
//	param_.push_back(vect);
//}
//
//void RegressionLS::regAll()
//{
//	//normalizationMatrix(cv_);
//	//normalizationMatrix(explanatoryData_);
//
//	param_.clear();
//	size_t nbSimu = explanatoryData_.size();
//	for(size_t exerciceIndex=0; exerciceIndex<nbExercice_; exerciceIndex++)
//	{
//		ublas::vector<double> cvVector(nbSimu);
//
//		for(size_t iteration=0; iteration<nbSimu; iteration++)
//			cvVector[iteration]=cv_[iteration][exerciceIndex];
//
//		Rmatrix explanary_i(nbSimu_, nbBasis_);
//		for(size_t rowIndex = 0; rowIndex<nbSimu_; rowIndex++)
//		{
//			for(size_t colIndex = 0; colIndex<nbBasis_; colIndex++)
//				explanary_i(rowIndex,colIndex)=explanatoryData_[rowIndex][exerciceIndex*nbBasis_+colIndex];
//		}
//
//		reg_OneExerciceTime(cvVector, explanary_i);
//	}
//}
//
//void RegressionLS::normalizationMatrix(std::vector<std::vector<double>>& m)
//{
//	assert(checkVectSize(m));
//	size_t nbRow = m.size();
//	assert(nbRow>0);
//	for(size_t colIndex=0; colIndex<m[0].size(); colIndex++)
//	{
//
//		double sum=0.0;
//		double squareSum=0.0;
//		for(size_t rowIndex=0; rowIndex<nbRow;rowIndex++)
//		{
//			sum+=m[rowIndex][colIndex];
//			squareSum+=m[rowIndex][colIndex]*m[rowIndex][colIndex];
//		}
//		double mean = sum/nbRow;
//		double variance	= squareSum/nbRow - mean*mean;
//		double sd = sqrt(variance);
//
//		//check if it is a constant column 
//		if(sd==0.0)
//		{
//			for(size_t rowIndex=0; rowIndex<nbRow;rowIndex++)
//			{
//				m[rowIndex][colIndex]/=sum;
//			}
//		}
//		else
//		{
//			for(size_t rowIndex=0; rowIndex<nbRow;rowIndex++)
//			{
//				m[rowIndex][colIndex]=(m[rowIndex][colIndex]-mean)/sd;
//			}
//		}
//	}
//}
//
//bool RegressionLS::checkVectSize(const std::vector<std::vector<double>>& m)const
//{
//	bool res = true;
//	if(!m.empty())
//	{
//		size_t length = m[0].size();
//		for(size_t i=1; i<m.size(); i++)
//		{
//			if(m[i].size()!=length)
//			{
//				res=false;
//				break;
//			}
//		}
//	}
//	return res;
//}
//
//bool RegressionLS::InvertMatrix (const ublas::matrix<double>& input, ublas::matrix<double>& inverse) 
//{  
//	typedef ublas::permutation_matrix<std::size_t> pmatrix; 
//	// create a working copy of the input 
//	ublas::matrix<double> A(input); 
//	// create a permutation matrix for the LU-factorization 
//	pmatrix pm(A.size1()); 
//	// perform LU-factorization 
//	int res = lu_factorize(A,pm); 
//	if( res != 0 )
//		return false; 
//	// create identity matrix of "inverse" 
//	inverse.assign(ublas::identity_matrix<double>(A.size1())); 
//	// backsubstitute to get the inverse 
//	lu_substitute(A, pm, inverse); 
//	return true; 
//}