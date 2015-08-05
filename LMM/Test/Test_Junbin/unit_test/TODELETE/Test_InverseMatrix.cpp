//#include <JBLMM/Test/JBTests.h>
//
//#include <iostream>
//
//#include <JBLMM/Longstaff_Schwartz/RegressionLS.h>
//
//void Test_RegressionLS()
//{
//	RegressionLS_PTR regressionLS_PTR(new RegressionLS());
//	Rmatrix matrixToInverse(2,2);
//	Rmatrix matrixInverse(2,2);
//	matrixToInverse(0,0)=0.8;
//	matrixToInverse(1,0)=0.6;
//	matrixToInverse(0,1)=-0.6;
//	matrixToInverse(1,1)=0.8;
//
//
//	//time_t _time;
//	//struct tm timeInfo;
//	//char format[32];
//	//time(&_time);
//	//localtime_s(&timeInfo, &_time);
//	//strftime(format, 32, "%Y-%m-%d", &timeInfo);
//
//	std::stringstream outputFileName_s; 
//	outputFileName_s<<"Test_InverseMatrix"<<".csv";
//	std::string outputFileName = LMMPATH::get_Root_OutputPath() + outputFileName_s.str();
//	ofstream out;
//	out.open(outputFileName,  ios::out | ios::app );
//	out<<endl;
//	out<<endl;
//	out<<endl;
//	//out << "origin matrix: " << endl; 
//	//size_t length = matrixToInverse.size1();
//	//for(size_t i = 0; i<length; i++)
//	//{
//	//	for(size_t j = 0; j<length; j++)
//	//		out << matrixToInverse(i,j) << "; ";
//	//	out << endl;
//	//}
//
//	//regressionLS_PTR->InvertMatrix(matrixToInverse,matrixInverse);
//	//out << "inverse matrix: " << endl; 
//	//for(size_t i = 0; i<length; i++)
//	//{
//	//	for(size_t j = 0; j<length; j++)
//	//		out << matrixInverse(i,j) << "; ";
//	//	out << endl;
//	//}
//
//	//ublas::vector<double> vc(2);
//	//vc[0]=1;
//	//vc[1]=2;
//
//	//Rmatrix ex(2,2);
//	//ex(0,0)=1;
//	//ex(1,0)=2;
//	//ex(0,1)=3;
//	//ex(1,1)=4;
//
// //	regressionLS_PTR->reg_OneExerciceTime(vc,ex);
//	//const Rvector& param= regressionLS_PTR->getParam()[0];
//	//out << "parameters: " << endl; 
//	//for(size_t i = 0; i<2; i++)
//	//	out << param[i] << "; ";
//	//out << endl;
//
//	out << "Test_Normalization: " << endl; 
//	size_t length_vect=3;
//	std::vector<std::vector<double>> vect_vect_(length_vect);
//
//	vect_vect_[0].push_back(1);
//	vect_vect_[0].push_back(3);
//	vect_vect_[0].push_back(23);
//	vect_vect_[1].push_back(6);
//	vect_vect_[1].push_back(2);
//	vect_vect_[1].push_back(47);
//	vect_vect_[2].push_back(4);
//	vect_vect_[2].push_back(17);
//	vect_vect_[2].push_back(19);
//
//	out << "origin matrix: " << endl; 
//	for(size_t i = 0; i<length_vect; i++)
//	{
//		for(size_t j = 0; j<length_vect; j++)
//			out << vect_vect_[i][j] << "; ";
//		out << endl;
//	}
//
//	regressionLS_PTR->normalizationMatrix(vect_vect_);
//	out << "matrix after normalization: " << endl;
//	for(size_t i = 0; i<length_vect; i++)
//	{
//		for(size_t j = 0; j<length_vect; j++)
//			out << vect_vect_[i][j] << "; ";
//		out << endl;
//	}
//
//}
