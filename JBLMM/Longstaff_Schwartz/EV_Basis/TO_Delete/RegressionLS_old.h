//#pragma once
//#include <vector>
//
//#include <boost/numeric/ublas/matrix.hpp>
//#include <boost/numeric/ublas/io.hpp> 
//#include <boost/numeric/ublas/vector_proxy.hpp> 
//#include <boost/numeric/ublas/matrix.hpp> 
//#include <boost/numeric/ublas/triangular.hpp> 
//#include <boost/numeric/ublas/lu.hpp> 
//
//#include <LMM/helper/Name.h>
//#include <JBLMM/Element/Rate.h>
//
//namespace ublas = boost::numeric::ublas; 
//
//typedef ublas::matrix<double> Rmatrix;
//typedef ublas::vector<double> Rvector;
//class RegressionLS
//{
//
//
//private:
//	size_t nbExercice_;
//	size_t nbBasis_;
//	size_t nbSimu_;
//	
//	mutable std::vector<std::vector<double>> explanatoryData_;	//size: (nbSimu)*(nbExercice_*nbBasis_)  vector of vector TODO
//	mutable std::vector<std::vector<double>> cv_;					//size: (nbSimu)*(nbExercice_)
//	std::vector<ublas::vector<double>> param_;				//size: (nbExercice_)*(nbBasis_)
//
//public:
//	//constructor par défaut
//	RegressionLS(){}
//
//	//don't use. copy a lot.
//	RegressionLS(	size_t nbExercice,
//					size_t nbBasis,
//					size_t nbSimu_,
//					std::vector<std::vector<double>>& explanatoryData,
//					std::vector<std::vector<double>>& continuationValueData);
//
//	//destructor
//	virtual ~RegressionLS(){}
//
//	//getter 
//	const std::vector<ublas::vector<double>>& getParam()const{return param_;}
//
//	//setters
//	void setNbExercice(size_t nbExercice){nbExercice_=nbExercice;}
//	void setNbSimulation(size_t nbSimu){nbSimu_=nbSimu;}
//	void setNbBasis(size_t nbBasis){nbBasis_=nbBasis;}
//
//	//assensors
//	std::vector<std::vector<double>>& getExplanatoryData()const{return explanatoryData_;}
//	std::vector<std::vector<double>>& getCV()const{return cv_;}
//
//	//clear Parameters
//	void clearParam(){param_.clear();}
//
//	//regression at One exercice time 
//	void reg_OneExerciceTime(const ublas::vector<double>& observable, const Rmatrix& explicatif);
//
//
//
//
//	void regAll();
//
//
////private:
//	//normalization
//	void normalizationMatrix(std::vector<std::vector<double>>& m);
//	//check if all the sub vector of the vector m have the
//	bool checkVectSize(const std::vector<std::vector<double>>& m)const;
//
//	//inversion
//	bool InvertMatrix(const ublas::matrix<double>& input, ublas::matrix<double>& inverse);
//
//	//void normalizationVector(Rvector& v)
//	//{
//	//	size_t nbElement = v.size();
//	//	assert(nbElement>0);
//	//	double sum=0.0;
//	//	double squareSum=0.0;
//	//	for(size_t i=0; i<nbElement;i++)
//	//	{
//	//		sum+=v[i];
//	//		squareSum+=v[i]*v[i];
//	//	}
//	//	double mean = sum/nbElement;
//	//	double variance	= squareSum/nbElement - mean*mean;
//	//	double sd = sqrt(variance);
//
//	//	//check if it is a constant column 
//	//	if(sd==0.0)
//	//	{
//	//		for(size_t i=0; i<nbElement;i++)
//	//		{
//	//			v[i]/=sum;
//	//		}
//	//	}
//	//	else
//	//	{
//	//		for(size_t i=0; i<nbElement;i++)
//	//		{
//	//			v[i]=(v[i]-mean)/sd;
//	//		}
//	//	}
//	//}
//};
//typedef boost::shared_ptr<RegressionLS> RegressionLS_PTR;
//typedef boost::shared_ptr<const RegressionLS> RegressionLS_CONSTPTR;