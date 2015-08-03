#include "Interpolator.h"
#include <cmath>


template<class T>
inline  const T &MAX( const T &a,  const T &b)
        {return b > a ? (b) : (a);}

template<class T>
inline const T &MIN( const T &a,  const  T &b)
        {return b > a ? (a) : (b);}



Interpolator :: Interpolator()
						:x_(0), y_(0), n_(-1), jsav_(-1), dj_(-1), y2(0){};

Interpolator :: Interpolator(const vector<double>& x_grid, const vector<double>& y_value)
						:x_(x_grid), y_(y_value), n_((int)x_grid.size()), jsav_(0), dj_(MIN(1,(int)pow((double)n_,0.25))), y2((int)x_grid.size(), 0.0)
{
    set_y2();
};

double  Interpolator :: interpolate(double x)
{
	int jl = locate(x);

	int klo=jl, khi=jl+1;
	double y,h,b,a;
	h=x_[khi]-x_[klo];
	if (h == 0.0) throw("Bad input to routine splint");
	a=(x_[khi]-x)/h;
	b=(x-x_[klo])/h;
	y=a*y_[klo]+b*y_[khi]+((a*a*a-a)*y2[klo]
		+(b*b*b-b)*y2[khi])*(h*h)/6.0;
	return y;
}

void Interpolator :: set_y2()
{
	double p,qn,sig,un;

	//int n=y2.size();
	int n = (int)x_.size();
	vector<double> u(n-1);
	//! boundary condtion, second derivative = 0
	//if (yp1 > 0.99e30)
	y2[0]=u[0]=0.0;
	//else {
	//	y2[0] = -0.5;
	//	u[0]=(3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
	//}
	for (int i=1; i<n-1; i++) 
	{
		sig=(x_[i]-x_[i-1])/(x_[i+1]-x_[i-1]);
		p=sig*y2[i-1]+2.0;
		y2[i]=(sig-1.0)/p;
		u[i]=(y_[i+1]-y_[i])/(x_[i+1]-x_[i]) - (y_[i]-y_[i-1])/(x_[i]-x_[i-1]);
		u[i]=(6.0*u[i]/(x_[i+1]-x_[i-1])-sig*u[i-1])/p;
	}
	//! boundary condtion, second derivative = 0
	//if (ypn > 0.99e30)
	qn=un=0.0;
	//else {
	//	qn=0.5;
	//	un=(3.0/(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
	//}
	y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.0);
	for (int k=n-2;k>=0;k--)
		y2[k]=y2[k]*y2[k+1]+u[k];
}


int Interpolator :: locate(double x)
{
	int ju,jm,jl;
	//if (n < 2 || mm < 2 || mm > n) throw("locate size error");
	bool  ascnd;
	if(x_[n_-1] >= x_[0])
	{
		ascnd = true;
	}
	else
	{
	    ascnd = false;
	}

	jl=0;
	ju=n_-1;
	while (ju-jl > 1) {
		jm = (ju+jl) >> 1;
		if (x >= x_[jm] == ascnd)
			jl=jm;
		else
			ju=jm;
	}

	int cor = abs(jl-jsav_) > dj_ ? 0 : 1;
	jsav_ = jl;
	int mm = 2;
	return MAX(0,MIN(n_-mm,jl-((mm-2)>>1)));
}


vector<double> Interpolator :: interpolate(vector<double>& x_interp)
{
	vector<double> y_interp(x_interp.size(),0.0);
    for(unsigned int i=0; i<y_interp.size(); ++i)
	{
	     y_interp[i] = interpolate(x_interp[i]);
	}
	return y_interp;
}


double Interpolator :: calculate_integral_approximation() const
{
	double integral_approximation = 0.0;
    for(int i=0; i<n_-1; ++i)
	{
		double x0		= x_[i];
		double x1		= x_[i+1];
		double constant = -1/24.0*pow(x1-x0,3);

		//! linear part
		double lin = ( y_[i] + y_[i+1] ) * (x1-x0)/2.0; 


		//! second derivative part 
        double curvature = constant*(y2[i+1] + y2[i]);
         
	    integral_approximation += (lin + curvature);
	}
	return integral_approximation; 
}



//! ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- 
//! 
//!				Interpolator_D2
//!
//! ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- 


Interpolator_D2::Interpolator_D2(const vector<double>& x_grid, 
								 const vector<double>& y_grid,
								 const Matrix& value)
								 :	x_grid_(x_grid),
									y_grid_(y_grid),
									grid_value_(value),
									y_interp_v(x_grid.size())     // use default interpolator

{
	//x_interp = Interpolator(x_grid,value);
	vector<double> y_value(y_grid.size(),0.0); 
	for(size_t i=0; i<x_grid.size(); ++i)
	{
		//! copy matrix row to vector
		for(size_t j=0; j<y_grid.size(); ++j)
		{
		    y_value[j] = grid_value_[i][j];
		}
		//! construct interpolator for each rows
	    y_interp_v[i] = Interpolator(y_grid,y_value);
	}
}


Matrix Interpolator_D2::interpolate(const vector<double>& x_grid_interp, const vector<double>& y_grid_interp)
{
	if( y_grid_interp[0] < y_grid_[0] || y_grid_interp.back() > y_grid_.back()
		||x_grid_interp[0] < x_grid_[0] || x_grid_interp.back() > x_grid_.back()
	  )
	{
		throw ("Error in function Interpolator_D2::interpolate, interplated grid in not contained in the original grids ! ");
	}

	//! suppose x_grid_interp and y_grid_interp are increasing ! 
    Matrix m(x_grid_interp.size(), y_grid_interp.size(), 0.0);   

	vector<double> x_value(x_grid_.size(),0.0);

	for(size_t j=0; j<y_grid_interp.size(); ++j)
	{
		//! construct x direction interpolator
		for(size_t i=0; i<x_grid_.size(); ++i)
		{
			x_value[i] = y_interp_v[i].interpolate(y_grid_interp[j]);
		}
		Interpolator x_interp(x_grid_, x_value);

		//! interpolate
		for(size_t i=0; i<x_grid_interp.size(); ++i)
		{
			m[i][j] = x_interp.interpolate(x_grid_interp[i]);
		}
	}
	return m;
}
