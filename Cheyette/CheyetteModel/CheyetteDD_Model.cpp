#include "CheyetteDD_Model.h"

////initialisation des elements statiques de CheyetteDD_Model
//
//	static const double matu[]={0., 1.,2., 3., 4., 5., 40., 15., 20., 25.};
//	static const double tauxZC[]={0.8/100, 0.85/100, 0.9/100, 0.92/100, 0.95/100, 1.00/100, 1.5/100, 2.0/100, 2.5/100, 2.3/100};
//	static std::vector<double> listeMatu(std::begin(matu),std::end(matu));
//	static std::vector<double> listeTauxZC(std::begin(tauxZC),std::end(tauxZC));
//	CourbeInput_PTR CheyetteDD_Model::courbeInput_PTR_(new CourbeInput(listeMatu, listeTauxZC));
//
//	Boost_R2R_Function CheyetteDD_Model::cheyetteDD_Parameter_.shift_(boost::bind(&CheyetteDD_Model::shift, _1, _1));

	void CheyetteDD_Model::show()
{
	std::cout << "---------------------------------------------" << std::endl ;
	std::cout << "--- creation d'un objet Cheyette DD Model ---" << std::endl ;
	std::cout << "  SIGMA piecewise constant " << std::endl ;
	cheyetteDD_Parameter_.sigma_.show() ;
	std::cout << "  M piecewise constant " << std::endl ;
	cheyetteDD_Parameter_.m_.show() ;
	//std::cout << "  SHIFT piecewise constant " << std::endl ;
	//cheyetteDD_Parameter_.shift_.show() ;
	std::cout << "  k constant " << std::endl ;
	std::cout << "       k = " << cheyetteDD_Parameter_.k_ << std::endl ;
	std::cout << "---------------------------------------------" << std::endl ;
}

//fonction de vol locale Displaced Diffusion
double CheyetteDD_Model::sigma_r(double t, double x_t) const 
{
	double sigma_t = cheyetteDD_Parameter_.sigma_(t) ;
	double m_t		= cheyetteDD_Parameter_.m_(t) ;
	double shift_t	= shift_->operator()(t, x_t) ;
	double shift_0	= shift_->operator()(0, 0) ;  //shift_(0, 0) ??
	return sigma_t * (m_t * shift_t + (1 - m_t) * shift_0) ;
}

double CheyetteDD_Model::sigma_r_t_1stDerivative(double t, double x_t) const 
{
	double sigma_t = cheyetteDD_Parameter_.sigma_(t) ;
	double m_t		= cheyetteDD_Parameter_.m_(t) ;
	return sigma_t * m_t * shift_1stDerivative(t, x_t)  ;
}

double CheyetteDD_Model::shift(double t, double x_t) const
{
	//cas shift = r(t)
	return x_t + courbeInput_PTR_->get_f_0_t(t) ;  
}

double CheyetteDD_Model::shift_1stDerivative(double t, double x_t) const
{
	//si shift = r(t)
	//r(t) = x(t) + f(0,t)
	return 1 ;
}

//fonctions G(t, T), ZC B(t, T)...
//Displaced Diffusion : k constant, calcul explicite de G(t, T)
double CheyetteDD_Model::G(double t, double T) const
{
	assert(t >= 0);  
	double k = cheyetteDD_Parameter_.k_ ;
	return 1/k * (1 - exp(- k * (T-t))) ;
}

double CheyetteDD_Model::P(double t, double T, double x_t, double y_t) const 
{
	double g = G(t,T) ;
	double P_0_t = exp( - courbeInput_PTR_->get_tauxZC0(t) * t ) ;
	double P_0_T = exp( - courbeInput_PTR_->get_tauxZC0(T) * T ) ;
	return P_0_T/P_0_t * exp(- x_t * g - 1/2 * y_t * g * g) ; 
}

//f_0_t � coder dans la classe courbeInput
double CheyetteDD_Model::r_t(double t, double x_t) const 
{
	double f_0_t = courbeInput_PTR_->get_f_0_t(t) ;
	return f_0_t + x_t ; 
}

//T2 = T1 + delta
double CheyetteDD_Model::Libor(double t, double T1, double T2, double x_t, double y_t) const
{
	double delta = T2 - T1 ;
	double P_t_T1 = P(t, T1, x_t, y_t) ;
	double P_t_T2 = P(t, T2, x_t, y_t) ;
	return 1/delta * (P_t_T1/P_t_T2 - 1) ;
}

//EDS : drift et diffusion sous Q^T
double CheyetteDD_Model::diffusion_x(double t, double x_t) const
{
	return sigma_r(t, x_t) ;
}

double CheyetteDD_Model::drift_x_QT(double t, double T_proba_fwd, double x_t, double y_t) const
{
	double k			= cheyetteDD_Parameter_.k_ ;
	double r_t			= x_t + courbeInput_PTR_->get_f_0_t(t) ;
	double sigma_r_t	=  sigma_r(t, x_t) ;

	return y_t - k * x_t - G(t, T_proba_fwd) * sigma_r_t * sigma_r_t ;	    
}


double CheyetteDD_Model::drift_y(double t, double x_t, double y_t) const
{
	double sigma_r_t	=  sigma_r(t, x_t) ;
	double k			= cheyetteDD_Parameter_.k_ ;
	return sigma_r_t * sigma_r_t - 2 * k * y_t ;
}

double CheyetteDD_Model::drift_x_Q(double t, double x_t, double y_t) const
{
	double k			= cheyetteDD_Parameter_.k_ ;

	return y_t - k * x_t ;	    
}

//EDS : drift et diffusion sous Q^Annuity
double CheyetteDD_Model::drift_x_QA(double t, double x_t, double y_t, const VanillaSwap& swap) const
{
	double k			= cheyetteDD_Parameter_.k_ ;
	double r_t			= x_t + courbeInput_PTR_->get_f_0_t(t) ;
	double sigma_r_t	=  sigma_r(t, x_t) ;

	double T_k, drift(y_t) ;
	double fixed_tenor = swap.get_fixedLegTenorType().YearFraction() ;
	double float_tenor = swap.get_floatingLegTenorType().YearFraction() ;
	double tenor_ref = std::min(fixed_tenor, float_tenor) ;   
	std::vector<size_t> index = swap.get_fixedLegPaymentIndexSchedule() ;
	for (size_t i = 0 ; i < index.size() ; ++i)
	{
		T_k = index[i] * tenor_ref ;
		drift += - k * x_t - G(t, T_k) * sigma_r_t * sigma_r_t ; 	
	}
	
	return drift ;	
}

//double CheyetteDD_Model::annuity(double t, double x_t, double y_t, const VanillaSwap& vanillaSwap)
//{
//	//TODO assert // check tenor[0] = 0; c'est � dire le 1er pilier de la tenor structure englobe t = 0 (date de valo) 
//
//	double price = 0.0;
//	const std::vector<size_t>& fixedLegPaymentIndexSchedule  = vanillaSwap.get_fixedLegPaymentIndexSchedule();
//	double dateEchangeFluxFixe, delta_T, ZC ;
//
//	double fixed_tenor = vanillaSwap.get_fixedLegTenorType().YearFraction() ;
//	double float_tenor = vanillaSwap.get_floatingLegTenorType().YearFraction() ;
//	double tenor_ref = std::min(fixed_tenor, float_tenor) ;  //le plus petit 
//
//	if (t == 0)
//	{
//		//somme sur tous les flux fixes
//		for(size_t itr = 0; itr < fixedLegPaymentIndexSchedule.size(); ++itr) 
//		{
//			//convertit l'indice/le numero du flux en la date de tomb�e du flux (ex : flux numero 2 survient � date 1Y)
//			dateEchangeFluxFixe = fixedLegPaymentIndexSchedule[itr] * tenor_ref ;
//
//			delta_T = vanillaSwap.get_DeltaTFixedLeg(itr);	
//			ZC		= exp( - courbeInput_PTR_->get_tauxZC0(dateEchangeFluxFixe) * dateEchangeFluxFixe);
//			price += delta_T * ZC ;		
//		}
//	}else
//	{
//		for(size_t itr = 0; itr < fixedLegPaymentIndexSchedule.size(); ++itr) 
//		{
//			dateEchangeFluxFixe = fixedLegPaymentIndexSchedule[itr] * tenor_ref ;
//			delta_T = vanillaSwap.get_DeltaTFixedLeg(itr) ;	
//			ZC		= P(t, dateEchangeFluxFixe, x_t, y_t) ;
//			price += delta_T * ZC ;			
//	}
//	return price;
//}


//double CheyetteDD_Model::txSwap(const VanillaSwap& vanillaSwap)
//{
//	const std::vector<size_t>& fixedLegPaymentIndexSchedule  = vanillaSwap.get_fixedLegPaymentIndexSchedule();
//
//	double fixed_tenor = vanillaSwap.get_fixedLegTenorType().YearFraction() ;
//	double float_tenor = vanillaSwap.get_floatingLegTenorType().YearFraction() ;
//	double tenor_ref = std::min(fixed_tenor, float_tenor) ;  //le plus petit 
//
//	double T_0 = fixedLegPaymentIndexSchedule[0] * tenor_ref ;
//	double ZC_T0 = exp( - courbeInput_PTR_->get_tauxZC0(T_0) * T_0) ;	
//
//	double T_N = fixedLegPaymentIndexSchedule[fixedLegPaymentIndexSchedule.size() - 1] * tenor_ref ;
//	double ZC_TN = exp( - courbeInput_PTR_->get_tauxZC0(T_N) * T_N)  ;	
//
//	double annuite = annuity(vanillaSwap) ;
//	return (ZC_T0 - ZC_TN)/ annuite ;
//}


//S(t, x_t, y_t) = S(t, x_t, y_bar_t) = S(t, x_t)  :  utile pour pouvoir inverser cette fonction par NR
//double S(double t, double x_t, double y_t, const VanillaSwap& vanillaSwap)
//{
//
//
//}



