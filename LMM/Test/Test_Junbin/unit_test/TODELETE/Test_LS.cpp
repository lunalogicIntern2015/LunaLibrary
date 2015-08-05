//#include <JBLMM/Test/JBTests.h>
//#include <LMM/Test/Tests.h>
//
//#include <iostream>
//#include <fstream> 
//#include <string> 
//#include <iterator>
//#include <algorithm>
//
//
//#include <LMM/RNGenerator/McGenerator.h>
//#include <LMM/RNGenerator/RNGenerator.h>
//#include <LMM/helper/InputFileManager.h>
//#include <LMM/numeric/NumericalMethods.h>
//
//#include <LMM/LmmModel/Correlation.h>
//#include <LMM/LmmModel/McTerminalLmm.h>
//
//#include <LMM/LmmModel/Lmm.h>
//#include <LMM/LmmModel/ConstShifted_HGVolatilityFunction.h>
//#include <LMM/LmmModel/LmmSwaptionMarketData.h>
//
//#include <JBLMM/Longstaff_Schwartz/EV_Basis_Collection.h>
//#include <JBLMM/Pricer/McGeneticSwapLMMPricer.h>
//#include <JBLMM/Instrument/InstrumentFactory.h>
//#include <JBLMM/Pricer/GeneticVanillaSwapPricer.h>
//#include <JBLMM/Longstaff_Schwartz/ExplanatoryVariable.h>
//#include <JBLMM/Longstaff_Schwartz/LS.h>
//
//#include <time.h>
//
//void Test_LS()
//{
//	std::vector<std::string> mkt_file_list = InputFileManager::get_VCUB_FileList();
//	std::string folder_name;// = "TotalCalib\\" ;  config.use_positive_constraint_=true;
//	std::string base_name_file = LMMPATH::get_BaseFileName(mkt_file_list[0]) + "\\";
//	folder_name+=base_name_file;
//	LMMPATH::reset_Output_SubFolder(folder_name );
//
//	LmmCalibrationConfig config;  //JB: default parameters set: model_nbYear_=16, correl_FullRank_ = 2*model_nbYear_ + 1  
//									// case fixed/float = 2 typically lmm 6MO/1YR
//
//	config.model_nbYear_=4;		// set model_nbYear_
//	config.correl_FullRank_=2*config.model_nbYear_+1;
//
//	LmmSwaptionMarketData_PTR pLmmSwaptionMarketData=JB_get_LmmSwaptionMarketData(config, mkt_file_list[0]);
//
//
//
//	LMM::Index  indexStart = 2;		//1Y
//	LMM::Index  indexEnd   = 8;		//4Y
//	Tenor	floatingLegTenorType = Tenor::_6M;
//	Tenor	fixedLegTenorType    = Tenor::_1YR;
//	assert(indexStart%2==0&&indexEnd%2==0);
//	double strike = pLmmSwaptionMarketData->get_SwaptionQuotes_ATM()->get_UpperTriangularVanillaSwaptionQuotes()(indexStart/2,(indexEnd-indexStart)/2).first.get_strike();
//	LMMTenorStructure_PTR lmmTenorStructure( new LMMTenorStructure(floatingLegTenorType, indexEnd/2));
//
//	VanillaSwap_CONSTPTR vanillaSwap_PTR(new VanillaSwap(strike, indexStart, indexEnd, floatingLegTenorType, fixedLegTenorType, lmmTenorStructure));
//	LMM::Index liborExercice0 = 2;
//	LMM::Index liborExercice1 = 4;
//	LMM::Index liborExercice2 = 6;
//	LMM::Index liborExercice3 = indexEnd;
//	std::vector<LMM::Index> exerciceIndex;
//	exerciceIndex.push_back(liborExercice0);
//	exerciceIndex.push_back(liborExercice1);	
//	//exerciceIndex.push_back(liborExercice2);
//	exerciceIndex.push_back(liborExercice3);
//	BermudanSwaption_PTR bermudanSwaption_PTR(new BermudanSwaption(vanillaSwap_PTR,exerciceIndex));
//
//	config.correl_ReducedRank_= 3; config.correl_alpha_ = 0.000000001 ; config.correl_beta_  = 0.05;
//	QuantLib::Array found_abcd = marketData_LMM_ABCD_calibration(config,pLmmSwaptionMarketData);
//
//	config.penalty_time_homogeneity_ = 1e-4 ; config.penalty_libor_ = 1e-6 ; config.use_positive_constraint_= true;
//	config.use_local_calib_=false;
//	Shifted_HGVolatilityFunction_PTR constShifted_HGVolatilityFunction_PTR = marketData_LMM_Global_gCalibration(config, pLmmSwaptionMarketData, found_abcd , create_InitCorrelation(config), GMatrixMapping_PTR());
//
//
//		//---------------------------Build Lmm and McLmm's structure--------------------------------------
//	////! Parameter of h
//	//double a = -0.06;
//	//double b = 0.17;
//	//double c = 0.54;
//	//double d = 0.17;
//	//Shifted_HGVolatilityParam::ABCDParameter abcdParam (a,b,c,d);
//	////Parameter of hg
//	//double g_constParam = 1.0;
//	//double shift_constParam = -0.01;
//	//ConstShifted_HGVolatilityParam_PTR hgParam( new ConstShifted_HGVolatilityParam(lmmTenorStructure,abcdParam,g_constParam,shift_constParam));
//	//
//	////! Correlation 1
//	//size_t nbFactor       = 3; // need to test nbFactor  = 3, and nbFactor = 
//	//size_t correlFullRank = lmmTenorStructure->get_horizon()+1;
//	//size_t correlReducedRank = nbFactor;
//	////!"Check Parameters(): Condition not implemented yet."
//	//std::cout << "checkParams(): ";
//	//CorrelationReductionType::CorrelationReductionType correlReductionType = CorrelationReductionType::PCA;
//	//double correlAlpha = 0.0;
//	//double correlBeta  = 0.1;
//	//Correlation_PTR correlation(new XY_beta_Correlation(correlFullRank,correlReducedRank, correlReductionType,correlAlpha,correlBeta));
//	//correlation->calculate(); // for print.
//	//correlation->print("test_McTerminalLmm_Correlation.csv");
//	////hgVolatilityFunction
//	//ConstShifted_HGVolatilityFunction_PTR hgVolatilityFunction (new ConstShifted_HGVolatilityFunction(lmmTenorStructure, correlation, hgParam)); 
//	//hgVolatilityFunction->print("test_McTerminalLmm_Volatility.csv");
//	//! Dispersion
//	Dispersion dispersion(constShifted_HGVolatilityFunction_PTR);
//
//	unsigned long seed = 5033;
//	RNGenerator_PTR  rnGenerator(new McGenerator(seed));
//
//	//build lmm and mcLmm model
//	Lmm_PTR shiftedLmm (new Lmm(dispersion));
//
//	const std::vector<double>& initLiborValues = pLmmSwaptionMarketData->get_LiborQuotes()->get_InitLibor();
//	McLmm_PTR mcLmm(new McTerminalLmm(shiftedLmm, initLiborValues, rnGenerator, MCSchemeType::EULER));
//	cout << endl; 
//	cout << "initLiborValues: ;" << endl; 
//	for(size_t i = 0; i<initLiborValues.size(); i++)
//	{
//			cout << initLiborValues[i] << "; ";	
//	}
//	cout<<endl;
//
//	//
//	//LmmVanillaSwaptionApproxPricer_Rebonato_PTR  lmm_Rebonato_PTR(new LmmVanillaSwaptionApproxPricer_Rebonato(shiftedLmm));
//	//double annuity = lmm_Rebonato_PTR->annuity0(*vanillaSwap_PTR.get(),initLiborValues);
//	//double blackPrice = annuity*NumericalMethods::Black_Price( 0.0427, 0.0427,0.157,1.0);
//
//
//	McLmmVanillaSwapPricer_PTR mcLmmVanillaSwapPricer_PTR(new McLmmVanillaSwapPricer(mcLmm));
//	ExplanatoryVariable::Evaluator_CONSTPTR evaluator_VanillaSwapRate_CONSTPTR(new EV_vanillaswaprate(mcLmmVanillaSwapPricer_PTR));
//	ExplanatoryVariable_CONSTPTR ev_VanillaSwapRate_CONSTPTR(new ExplanatoryVariable(evaluator_VanillaSwapRate_CONSTPTR));
//	//ExplanatoryVariable_CONSTPTR ra[1] = {ev_VanillaSwapRate_CONSTPTR};
//	//const std::vector<ExplanatoryVariable_CONSTPTR> ev_vect(ra,ra+1);
//
//	std::vector<size_t> coef(1);
//	coef[0]=1;
//	Basis::Transformer_PTR polynomialTransformator(new Polynomial(coef));
//	Basis_CONSTPTR basis1_CONSTPTR(new SingleEvBasis(polynomialTransformator,ev_VanillaSwapRate_CONSTPTR));
//
//	ExplanatoryVariable::Evaluator_CONSTPTR evaluator_libor_CONSTPTR(new EV_Libor());
//	ExplanatoryVariable_CONSTPTR ev_Libor_CONSTPTR(new ExplanatoryVariable(evaluator_libor_CONSTPTR));
//	Basis_CONSTPTR basis2_CONSTPTR(new SingleEvBasis(polynomialTransformator,ev_Libor_CONSTPTR));
//
//	std::vector<Basis_CONSTPTR> basisVect(2);
//	basisVect[0]=basis1_CONSTPTR;
//	basisVect[1]=basis2_CONSTPTR;
//
//	//std::vector<size_t> evIndex1(1);
//	//evIndex1[0]=0;
//	//std::vector<size_t> evIndex2(1);
//	//evIndex2[0]=0;
//	//EV_Basis_pair ev_Basis_pair1(evIndex1,basis1_CONSTPTR);
//	//EV_Basis_pair ev_Basis_pair2(evIndex2,basis2_CONSTPTR);
//	//EV_Basis_pair rb[2] = {ev_Basis_pair1,ev_Basis_pair2};
//	//const std::vector<EV_Basis_pair> evBasisPairvect(rb,rb+2);
//	//EV_Basis_Collection_PTR ev_Basis_Collection_PTR(new EV_Basis_Collection(ev_vect,evBasisPairvect));
//
//	LS_PTR ls_PTR(new LS(mcLmm));
//	size_t nbSimu1=50000;
//	size_t nbSimu2=50000;
//	RegressionLS_PTR rgLS(new RegressionLS());
//
//	double price = ls_PTR->simulate(	bermudanSwaption_PTR, 
//										basisVect,
//										rgLS,
//										nbSimu1,
//										nbSimu2);
//	const std::vector<ublas::vector<double>>& paramMatrix = rgLS->getParam();									
//	//const ublas::vector<double>& paramVect = paramMatrix[0];	
//	std::stringstream outputFileName_s; 
//	outputFileName_s<<"Test_LS"<<".csv";
//	std::string outputFileName = LMMPATH::get_Root_OutputPath() + outputFileName_s.str();
//	ofstream out;
//	out.open(outputFileName,  ios::out | ios::app );
//	out<<endl;
//	out<<endl;
//	out<<endl;
//	out << "initLiborValues: ;" << endl; 
//	for(size_t i = 0; i<initLiborValues.size(); i++)
//	{
//			out << initLiborValues[i] << "; ";	
//	}
//	out<<endl;
//	out << "exercice index vector: ;" << endl; 
//	for(size_t i = 0; i<exerciceIndex.size(); i++)
//	{
//			out << exerciceIndex[i] << "; ";	
//	}
//	out << endl;
//	out << "parameters: ;" << endl; 
//	out << " ;";
//	for(size_t i = 0; i<paramMatrix[0].size(); i++)
//	{
//		out << "Beta" << i << ": ;";
//	}
//	out << endl;
//	for(size_t i = 0; i<paramMatrix.size(); i++)
//	{
//		out << "T " << i << ": ;";
//		for(size_t j = 0; j<paramMatrix[0].size(); j++)
//			out << paramMatrix[i][j] << "; ";
//		out << endl;
//	}
//	out <<"price" << endl;
//	out <<price << endl;
//
//	cout << "exercice index vector: ;" << endl; 
//	for(size_t i = 0; i<exerciceIndex.size(); i++)
//	{
//			cout << exerciceIndex[i] << "; ";	
//	}
//	cout << endl;
//	cout << "parameters: " << endl; 
//	for(size_t i = 0; i<paramMatrix.size(); i++)
//	{
//		for(size_t j = 0; j<paramMatrix[0].size(); j++)
//			cout << paramMatrix[i][j] << " ";
//		cout << endl;
//	}
//	cout <<"price: " << price <<endl;
//}