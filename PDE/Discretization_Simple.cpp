#include "Discretization_Simple.h"
#include <algorithm>
#include "useful_function.h"
#include "Range.h"
#include "Params.h"
//#include <vector>

bool equal_minimum_distance(double x, double y)
{
    if(fabs(x-y) < 0.000000001)
		return true;
	else 
		return false;
}

//const double epsilon_zero = 0.00000000001;
Discretization_Simple :: Discretization_Simple(  const Range& range,		  
											     int sizeDiscretization,
											     const std::vector<bool>&   if_nonuniform, 
												 const std::vector<double>& nonuniform_center_c,
												 const std::vector<double>& nonuniform_scale_param_c,
												 const std::vector<bool>&   if_shift,
												 const std::vector<double>& shifting_center
								 ) 
								      : 
										range_(range),
										sizeDiscretization_(sizeDiscretization),

										if_nonuniform_(if_nonuniform),
										nonuniform_center_c_(nonuniform_center_c),
										nonuniform_scale_param_c_(nonuniform_scale_param_c),

										if_shift_(if_shift),
										shifting_center_(shifting_center)

{
	if(sizeDiscretization<2)
	  throw ("error: size of discretization, in discretization::discretization(...)" );

	
	discretization_  = discret(sizeDiscretization,0); // basic discretization
    int  basic_discretization_size = discretization_.size();
	if(if_nonuniform.size()>1)
	{
	     if(if_nonuniform.size()!=2)
		 {
			 cout << "if_nonuniform.size() = " << if_nonuniform.size() << endl;
		     throw("Error not implemtent the case of more than 3 dim ! ");
		 }

		 //! to simplity the code suppose discretization_[0] = discretization_temp[0] = 0.0		 
		 std::vector<double> discretization_temp =  discret(sizeDiscretization,1);  // new discretization 
		 int index_i = 0;
		 int index_j = 1;
		 std::vector<int> compter(discretization_.size(),0);

         for(size_t i=0; i<discretization_.size(); ++i)
		 {
			 double  x_i_this = discretization_[i];
			 double  x_i_next = -1; // 
			 if(i!=discretization_.size()-1)
				 x_i_next = discretization_[i+1];
			 else
				 x_i_next = 9999999;

			 for(size_t j=index_j; j<discretization_temp.size(); ++j)
			 {
				 double x_j = discretization_temp[j];
			     if(x_j>x_i_this+0.00000001 && x_j<x_i_next-0.00000001)
				 {
					 compter[i]++;
				 }
				 else if(x_j>=x_i_next-0.00000001 && x_j<=x_i_next+0.00000001)
				 {
					 //do nothing
				 }
				 else 
				 {
				     index_j = j;
					 //index_i ++;
					 break;
				 }
			 }
		 }

		 for(int i=0; i<basic_discretization_size; ++i)
		 {
			 if(i!=discretization_.size()-1)
			 {
				 int    n = compter[i];
				 double x = discretization_[i];
				 double interval = discretization_[i+1] - discretization_[i];
				 if(n!=0)
				 {
					 double step = interval/(n+1);
					 for(int k=0; k<n; ++k)
					 {
						discretization_.push_back(x+(k+1)*step);
					 }
				 }
			 }
			 else
			 {
				 int n = compter[discretization_.size()-1];
				 if(n!=0)
				 {
					 for(int k=0; k<n; ++k)
					 {
						 discretization_.push_back(discretization_temp[discretization_temp.size()-(k+1)]);
					 }
				 }
			 }
		 }
		 //! sorting
		 sort(discretization_.begin(),discretization_.end());  
		 vector<double>::iterator itr = unique(discretization_.begin(),discretization_.end(),equal_minimum_distance);
         discretization_.resize(itr-discretization_.begin());
	}
	//for(int l=0; l<discretization_.size(); ++l)
	//cout << discretization_[l] << "  ";
	//cout << endl;
	//getchar();


	sizeDiscretization_ = discretization_.size();
	initialize_other_tables();
	
	//discretization_ = std::vector<double>(sizeDiscretization,0.0);
 //   //! Uniform discretization
	//if(if_nonuniform_ == false)
 //   {
	//	 uniform_discretisize_vector(range_,discretization_);
	//}
	//else
	//{
	//	//! only for the Range [0,x_max], TODO: implemnet general case ... 
	//	//cout.precision(15);
	//	//cout << "if_uniform_x    range_x.get_leftBorder() = " << abs(range_x.get_leftBorder()) << endl;
 //       //if(abs(range_x.get_leftBorder()) > epsilon_zero)
	//	if(range_.get_leftBorder() != 0.0)
	//		throw("Error in constructor of discretization, non-uniform discretization now only works for Range_x = [0,x_max]");
	//
	//	double K = nonuniform_center_c_;
	//	double c = nonuniform_scale_param_c_; //! usually use: x_nonuniform_scale_param_c = 5.
	//	double max = range_.get_rightBorder();
	//    double epsilon_pace  = ( inverse_sinh((max-K)/c) -inverse_sinh(-K/c) )/(discretization_.size()-1);
	//   
	//    for(int i=0; i<(int)discretization_.size(); ++i)
	//    {
	//	   double epsilon_i = inverse_sinh(-K/c) + i*epsilon_pace;
	//	   discretization_[i] = K + c*sinh(epsilon_i);
	//    }
	//}

	////! shiftting
	//if(if_shift_ == true) // shift for uniform or nonuniform x_grid
	//{
	//	 //! only for the Range [0,x_max], TODO: implemnet general case ... 
	//	 //cout << "if_shift_x  range_x.get_leftBorder() = " << range_x.get_leftBorder() << endl;
 //        //if(abs(range_x.get_leftBorder()) > epsilon_zero)
	//	if(range_.get_leftBorder() > 0.0)
	//		throw("Error in constructor of discretization, shift discretization now only works for Range_x = [0,x_max]");
	//
	//	 double K = shifting_center_;

	//	 if(nonuniform_center_c_!=0) 
	//	 {
	//		 ////! grid shiftting, smoothing the payoff  --> nu bi :) 
	//		 double shifting = 0; //! can be positive or negative ...
	//		 for(unsigned int i=0; i<discretization_.size()-1; ++i)
	//		 {
	//		   if(discretization_[i] > K)
	//		   {
	//			   double dis    = discretization_[i] - K;
	//			   double middle = (discretization_[i] - discretization_[i-1])/2 ;
	//			   shifting = (middle - dis);
	//			   break;
	//		   }
	//		   else if(discretization_[i] ==K)
	//		   {   
	//			   shifting = (discretization_[i] - discretization_[i-1])/2;
	//			   //if(discret_x[i-1]-discret_x[i-2] != discret_x[i+1]- discret_x[i])
	//			   //	   throw("Error, the grid cannot be shifted to make K in the middle of 2 grid :) Not sure if the condition is good :) ");
	//			   break;
	//		   }
	//		 }

	//		 for(unsigned int i=0; i<discretization_.size(); ++i)
	//		 {
	//			 discretization_[i] += shifting;
	//		 }

	//		 //! delete < 0 element, because of the shifting ... 
	//		 sort(discretization_.begin(),discretization_.end());   
	//		 vector<double>::iterator itr = find_if(discretization_.begin(),discretization_.end(), if_positive);
	//		 if(itr != discretization_.end())
	//		 {
	//			discretization_.erase(discretization_.begin(), itr);
	//		 }

	//		 //! add zero
	//		 discretization_.push_back(0); //! for shifting :)
	//           
	//		 //! Fwd PDE
	//		 if(if_fwd_PDE==true)
	//			discretization_.push_back(shifting_center_); //! for shifting :)


	//		 //! delet duplicated grid:
	//		 sort(discretization_.begin(),discretization_.end()); 
	//		 unique(discretization_.begin(),discretization_.end(), minimun_distance...);  //! "unique": delete the consecutive duplicated elements 
	//		 discretization_.resize(...)
	//		 sizeDiscretization_ = (int)(discretization_.size());
	//	 }
	//	 else
	//	 {
	//		 //! dont shift only add shift
	//		 //! Fwd PDE
	//		 if(if_fwd_PDE==true)
	//			discretization_.push_back(shifting_center_); //! for shifting :)

	//		 //! delet duplicated grid:
	//		 sort(discretization_.begin(),discretization_.end()); 
	//		 unique(discretization_.begin(),discretization_.end(), minimum_distance ... );  //! "unique": delete the consecutive duplicated elements 
	//		 discretization_.resize(...)
	//		 sizeDiscretization_ = (int)(discretization_.size());
	//	 }
	//}
	//initialize_other_tables();
}


std::vector<double> Discretization_Simple :: discret(int sizeDiscretization, int index)
{
	std::vector<double> discretization_temp(sizeDiscretization,0);

 //! Uniform discretization
	if(if_nonuniform_[index] == false)
    {
		 uniform_discretisize_vector(range_,discretization_temp);
	}
	else
	{
		//! only for the Range [0,x_max], TODO: implemnet general case ... 
		//cout.precision(15);
		//cout << "if_uniform_x    range_x.get_leftBorder() = " << abs(range_x.get_leftBorder()) << endl;
        //if(abs(range_x.get_leftBorder()) > epsilon_zero)
		if(range_.get_leftBorder() != 0.0)
			throw("Error in constructor of discretization, non-uniform discretization now only works for Range_x = [0,x_max]");
	
		double K = nonuniform_center_c_[index];
		double c = nonuniform_scale_param_c_[index]; //! usually use: x_nonuniform_scale_param_c = 5.
		double max = range_.get_rightBorder();
	    double epsilon_pace  = ( inverse_sinh((max-K)/c) -inverse_sinh(-K/c) )/(discretization_temp.size()-1);
	   
	    for(int i=0; i<(int)discretization_temp.size(); ++i)
	    {
		   double epsilon_i = inverse_sinh(-K/c) + i*epsilon_pace;
		   discretization_temp[i] = K + c*sinh(epsilon_i);
	    }
	}
	cout << "discretization_temp's size = " << discretization_temp.size() << endl;

	//! shiftting
	if(if_shift_[index] == true) // shift for uniform or nonuniform x_grid
	{
		 //! only for the Range [0,x_max], TODO: implemnet general case ... 
		 //cout << "if_shift_x  range_x.get_leftBorder() = " << range_x.get_leftBorder() << endl;
         //if(abs(range_x.get_leftBorder()) > epsilon_zero)
		 if(range_.get_leftBorder() > 0.0)
			throw("Error in constructor of discretization, shift discretization now only works for Range_x = [0,x_max]");
	
		 double K = shifting_center_[index];

		 //if(nonuniform_center_c_[index]!=0) 
		 //{
			 double left_most_interval_size = discretization_temp[1] - discretization_temp[0];
			 ////! grid shiftting, smoothing the payoff  --> nu bi :) 
			 double shifting = 0; //! can be positive or negative ...
			 for(unsigned int i=0; i<discretization_temp.size()-1; ++i)
			 {
			   if(discretization_temp[i] > K)
			   {
				   double dis    = discretization_temp[i] - K;
				   double middle = (discretization_temp[i] - discretization_temp[i-1])/2 ;
				   shifting = (middle - dis);
				   break;
			   }
			   else if(discretization_temp[i] ==K)
			   {   
				   shifting = (discretization_temp[i] - discretization_temp[i-1])/2;
				   std::cout << "discretization_temp[i] = " << discretization_temp[i]  << "  discretization_temp[i-1] = " << discretization_temp[i-1]  << std::endl;
				   //if(discret_x[i-1]-discret_x[i-2] != discret_x[i+1]- discret_x[i])
				   //	   throw("Error, the grid cannot be shifted to make K in the middle of 2 grid :) Not sure if the condition is good :) ");
				   break;
			   }
			 }

			 for(unsigned int i=0; i<discretization_temp.size(); ++i)
			 {
				 discretization_temp[i] += shifting;
			 }

			 //! delete < 0 element, because of the shifting ... 
			 sort(discretization_temp.begin(),discretization_temp.end());   
			 vector<double>::iterator itr = find_if(discretization_temp.begin(),discretization_temp.end(), if_positive);
			 if(itr != discretization_temp.end())
			 {
				discretization_temp.erase(discretization_temp.begin(), itr);
			 }

			 //! add zero
			 discretization_temp.push_back(0); //! for shifting :)

			 //! maybe add a flag for this ! 
			 //! left_most_interval_size ---- make sure the density around 0 after shfiting is the same as before ...
             int n       = (int)(discretization_temp[1]/left_most_interval_size);
			 double step = discretization_temp[1]/(n+1);
			 for(int k=0; k<n; ++k)
			 {
				discretization_temp.push_back(step*(k+1));
				//cout << "step*(k+1)  = " << step*(k+1) << endl;
			 }
			 //cout << "discretization_temp[1] = " << discretization_temp[1] << endl;
			 //getchar();

	           
			 //! Fwd PDE
			 if(if_fwd_PDE==true)
				discretization_temp.push_back(shifting_center_[index]); //! for shifting :)


			 //! delet duplicated grid:
			 sort(discretization_temp.begin(),discretization_temp.end()); 
			 itr = unique(discretization_temp.begin(),discretization_temp.end(),equal_minimum_distance);  //! "unique": delete the consecutive duplicated elements 
             discretization_temp.resize(itr-discretization_temp.begin());
			 //sizeDiscretization_temp = (int)(discretization_temp.size());
		 //}
		 //else
		 //{
			// //! dont shift only add shift
			// //! Fwd PDE
			// if(if_fwd_PDE==true)
			//	discretization_temp.push_back(shifting_center_[index]); //! for shifting :)

			// //! delet duplicated grid:
			// sort(discretization_temp.begin(),discretization_temp.end()); 
			// unique(discretization_temp.begin(),discretization_temp.end(), cretiere);  //! "unique": delete the consecutive duplicated elements 
			// discretization_temp.resize(itr-discretization_temp.begin());
			// //sizeDiscretization_temp = (int)(discretization_temp.size());
		 //}
	}
	else
	{
       		 //! Fwd PDE
			 if(if_fwd_PDE==true)
				discretization_temp.push_back(fwd_PDE_v0); //! for shifting :)



			 //! delet duplicated grid: using "unique"
			 //! 1. unique compare only the consequtive elementes (need sort before calling unique)
			 //! 2. return the end of the new arry (those after the new pointer are not predictable, so need to delete the rest !)
			 sort(discretization_temp.begin(),discretization_temp.end()); 
			 vector<double>::iterator itr = unique(discretization_temp.begin(),discretization_temp.end(), equal_minimum_distance);  
			 discretization_temp.resize( itr - discretization_temp.begin() );  
	}
	cout << "discretization_temp's size = " << discretization_temp.size() << endl;
	return discretization_temp;
}

//Discretization_Simple Discretization_Simple :: create_new_discretization_by_refine(int refine_ratio)   //! refine each interval to refine_ration times, equal-distance;
//{
//    throw("Error in function  Discretization_Simple :: create_new_discretization_by_refine(...), function not implemented.");
//    /*return discretization(  range_t_,   (sizeDiscretization_t_-1)*refine_ratio_time  + 1,
//						    range_x_,   (sizeDiscretization_x_-1)*refine_ratio_space + 1, 
//							if_uniform_x_, x_nonuniform_center_c_ , x_nonuniform_scale_param_c_,
//							if_shift_x_,   x_shifting_center_);*/
//}

void Discretization_Simple :: initialize_other_tables()
{
    //! discretization_delta_  
    discretization_delta_ = std::vector<double>(discretization_.size()-1,0);
    for(unsigned int i=0; i<discretization_delta_.size(); ++i)
    {
        discretization_delta_[i] = discretization_[i+1] - discretization_[i];
    }

    //! discretization_tilde_
    discretization_tilde_= std::vector<double>(discretization_.size()-2,0);
    for(unsigned int i = 0; i<discretization_tilde_.size(); i++)
	 	discretization_tilde_[i] = discretization_[i+1];
}


void Discretization_Simple::uniform_discretisize_vector(Range& range, std::vector<double>& discretization) const
{
    discretization[0] = range.get_leftBorder();
    discretization[discretization.size()-1] = range.get_rightBorder();
    double pas = (range.get_rightBorder() - range.get_leftBorder())/(discretization.size()-1);
    for(unsigned int i=1; i<discretization.size()-1; i++)
    {
	 	 discretization[i] = discretization[i-1] + pas;
    }
}

void Discretization_Simple :: add_discret(std::vector<double> & discret_v)  // add t and reorder :) 
{
	//! insert t
	for(int i=0; i<(int)discret_v.size(); ++i)
	{
		discretization_.push_back(discret_v[i]);
	}
	//! sort / delete duplicated elementes
	sort(discretization_.begin(), discretization_.end());
	vector<double>::iterator itr = unique(discretization_.begin(),discretization_.end(), equal_minimum_distance);
	discretization_.resize( itr - discretization_.begin() );  

	initialize_other_tables();  // not efficient, but hehe :) I don't care ...
}


double Discretization_Simple:: get_discret(int index) const
{
	if(index<0 || index>= (int)discretization_.size())
		throw ("error of the position in time discretization, in discretization::get_discret_t(int i)" );
	return discretization_[index];
}

double Discretization_Simple :: get_discret_tilde(int index) const
{
	if(index<0 || index>= (int)discretization_tilde_.size())
		throw ( "error of the position in time discretization, in discretization::get_discret_x(int i)");
	return discretization_tilde_[index];
}

double Discretization_Simple :: get_discret_delta(int index) const
{
	if(index<0 || index>= (int)discretization_delta_.size())
		throw ( "error of the position in time discretization, in discretization::get_discret_x(int i)");
	return discretization_delta_[index];
}

double Discretization_Simple :: get_discret_pas() const
{
    if(if_uniform() == true)
	{
		return get_discret_delta(0);    
	}
	else
	{
		throw ("Error Discretization_Simple :: get_discret_pas(), works only for uniform discretization without any other modification of grid.");    
	}
}

bool Discretization_Simple :: if_uniform() const
{
	if(if_nonuniform_.size()==1 && if_nonuniform_[0] == false)
		return true;
	else 
		return false;  //! 2 uniform together is considered as nonuniform! 
  
	//if(if_nonuniform_==false && if_shift_==false)
	//{
	//	return true;
	//}
	//else
	//{
	//	return false;
	//}
}


//void Discretization_Simple::refine(int ratio)
//{
//     return discretization(  range_t_,   (sizeDiscretization_t_-1)*refine_ratio_time  + 1,
//						    range_x_,   (sizeDiscretization_x_-1)*refine_ratio_space + 1, 
//							if_uniform_x_, x_nonuniform_center_c_ , x_nonuniform_scale_param_c_,
//							if_shift_x_,   x_shifting_center_);
//
//	int    sizeDiscretization_;
//								  
//	std::vector<double> discretization_;
//	std::vector<double> discretization_tilde_;  // size = x.size()-2 (without the firste and the last)
//	std::vector<double> discretization_delta_;  // discretization_x_delta_[i] = discretization_x_[i+1] - discretization_x_[i]
//}

void Discretization_Simple :: print() const
{
	std::cout << "  discretization: " << std::endl;

	for(int i=0; i<(int)discretization_.size(); ++i)
	{
		std::cout << "discret[" << i << "] = " << discretization_[i]<< std::endl;
	}
	//system("Pause");
}

