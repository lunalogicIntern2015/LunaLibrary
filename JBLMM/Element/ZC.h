#pragma once
#include <vector>
#include <boost/numeric/ublas/matrix.hpp>

#include <LMM/helper/LMMTenorStructure.h>
#include <LMM/LmmModel/McLmm.h>

typedef boost::numeric::ublas::matrix<double> matrix;

class ZC
{

	matrix zcMatrix_;

public:
	ZC(matrix liborMatrix, std::vector<double>& initLiborValues, LMMTenorStructure_PTR lmmTenorStructure_PTR)
		:
		zcMatrix_(liborMatrix.size1(),liborMatrix.size2())
	{
		initialize_ZCMatrix(liborMatrix, initLiborValues, lmmTenorStructure_PTR);
	}
	virtual ~ZC(){}
	//getter
	const matrix& getZcMatrix()const{return zcMatrix_;}

private:
	void initialize_ZCMatrix(matrix liborMatrix, std::vector<double>& initLiborValues, LMMTenorStructure_PTR lmmTenorStructure_PTR)
	{
		//initialize forward libor at time 0
		zcMatrix_(0,0)=1.0;
		for(size_t i=1; i<zcMatrix_.size2(); i++)
		{
			size_t liborIndex = i-1;
			double deltaT = lmmTenorStructure_PTR->get_deltaT(liborIndex);
			zcMatrix_(0,i)=zcMatrix_(0,i-1)/(1+deltaT*initLiborValues[liborIndex]);
		}

		for(size_t i=1; i<zcMatrix_.size1(); i++)
		{
			zcMatrix_(i,i)=1.0;
			for(size_t j=i+1; j<zcMatrix_.size2(); j++)
			{
				size_t liborIndex = j-1;
				double deltaT = lmmTenorStructure_PTR->get_deltaT(liborIndex);
				zcMatrix_(i,j)=zcMatrix_(i,j-1)/(1+deltaT*liborMatrix(liborIndex,i));				
			}
		}
	}
};

