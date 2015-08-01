//#include <JBLMM/Test/JBTests.h>
//
//#include <iostream>
//
//#include <JBLMM/Longstaff_Schwartz/Basis.h>
//
//
//void Test_basis()
//{
//
//	CoefAndPower coefAndPower(3);
//	
//	coefAndPower[0].first=2.0;
//	coefAndPower[0].second.push_back(2);
//	coefAndPower[0].second.push_back(0);
//	coefAndPower[1].first=2.0;
//	coefAndPower[1].second.push_back(0);
//	coefAndPower[1].second.push_back(2);
//	coefAndPower[2].first=-4.0;
//	coefAndPower[2].second.push_back(1);
//	coefAndPower[2].second.push_back(1);
//	Basis::Transformer_PTR polynomialTransformator(new Polynomial(2, coefAndPower));
//	std::vector<double> inputValues;
//	inputValues.push_back(1.0);
//	inputValues.push_back(2.0);
//	double result_Polynomial = polynomialTransformator->evaluate_vect(inputValues);
//
//	std::stringstream outputFileName_s; 
//	outputFileName_s<<"Test_Basis"<<".csv";
//	std::string outputFileName = LMMPATH::get_Root_OutputPath() + outputFileName_s.str();
//	ofstream out;
//	out.open(outputFileName,  ios::out | ios::app );
//	out<<endl;
//	out<<endl;
//	out<<endl;
//	out<< "Polynomial ; coefficients ; ; powers " <<endl;
//	for(size_t i = 0; i<coefAndPower.size(); i++)
//	{
//		out << "; " << coefAndPower[i].first << "; ; ";
//		for(size_t j = 0; j<coefAndPower[i].second.size(); j++)
//			out << coefAndPower[i].second[j] << "; ";
//		out << endl;
//	}
//	out<< "Polynomial result" << endl;
//	out<< result_Polynomial <<endl;
//}