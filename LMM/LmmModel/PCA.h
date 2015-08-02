#pragma once

//YY sadly this is already implemented in: <ql/math/matrixutilities/pseudosqrt.hpp>
//ctnrefactor move to numeric
#include <ql/math/matrix.hpp>




class PCA
{
public:

	//typedef QuantLib::Matrix Matrx;

	//! YY TODO: Now prefer to keep the algorithm clear, so it is highly inefficient. To improve latter. 
	static QuantLib::Matrix doPCA(const QuantLib::Matrix& originalMatrix, size_t reducedRank, bool normalizeDiagonal);

	static bool checkEigenvalue(QuantLib::Array& eigenvalues, size_t reducedRank);
};
