//! some Questions:
//! 1. When I need to use the "shfifted_FT"
//! 2. What's the relation of the windows of the time and spectrum : information theory ???
//! 3. what's the analysis of "time-spectrum" ... 

#include "FFT.h"
#include <PDE/NR/code/fourier.h>
using namespace std;

const double pi  = 3.14159265358979;
//const double two_pi = 6.28;


//! FFT transformation for function f in interval [a,b]
//! We want to calculate the Fourier transformation of function f, 
//! but to use the centered windows, we have to calcualte first FFT for f(x-a) then to refind FFT of f(x).


FFT :: FFT(double border, int n, Func f, double shift)
			: m_left(-border), m_right(border), m_n(n), m_N((int)pow(2.0,n)), m_shift(shift),
			  m_x(m_N,0.0),
			  m_f(m_N,0.0),
			  m_xi(m_N,0.0)
			  //m_fft_f(m_N,0.0)
{
	if(border <=0)
		throw ("Error in construtor of class FFT, border should be >0 ! ");

    m_z  = new double[2*m_N];
	set_x();
	set_xi();
	set_f(f);
}

FFT::~FFT()
{
    delete[] m_z;
}

void FFT :: set_x()
{
	double middle    = (m_left + m_right)/2.0;
	double pace      = (m_right - m_left)/m_N;  // in fact use N+1 points, but as the function is periodic, so the (left & right)boundary value is saved at the same points.

	//! set x 
	//! 1st part: middle points --> the right points(not included) 
    m_x[0] = middle;
	for(int i=1; i<m_N/2; ++i)  
	{
	    m_x[i] = m_x[i-1] + pace;
	}
	//! 2nd part: left point --> the middle points (not included), 
	//! left points represent both the value of left point and right point :) ---- because of the periodical feature.
	m_x[m_N/2] = m_left;
	for(int i=m_N/2+1; i<m_N; ++i)
	{
	    m_x[i] = m_x[i-1] + pace;
	}
}

void FFT :: set_xi()
{
	//! 1st part: middle points --> the right points(not included) 
	double middle = (m_right + m_left)/2.0;
	double T = (m_right - m_left);
	double pace_xi   = 1.0/T;
    m_xi[0] = middle;
	for(int i=1; i<m_N/2; ++i)  
	{
	    m_xi[i] = m_xi[i-1] + pace_xi;
	}
	//! 2nd part: left point --> the middle points (not included), 
	//! left points represent both the value of left point and right point :) ---- because of the periodical feature.
	m_xi[m_N/2] = -m_N/2*pace_xi;
	for(int i=m_N/2+1; i<m_N; ++i)
	{
	    m_xi[i] = m_xi[i-1] + pace_xi;
	}
}


void FFT :: set_f(Func f)
{
    //vector<complex<double>> f(x.size(),0.0);
	for(int i=0; i<m_N; ++i)
	{
	    m_f[i] = f(m_x[i]-m_shift);
	}
}

void FFT :: calculate_FFT() // use the function test_func as input function ...
//! fft_result is size of N = 2^n
{
	double T = (m_right-m_left);
	for(int i=0; i<m_N; i++)
	{
		m_z[2*i]   = m_f[i].real();
		m_z[2*i+1] = m_f[i].imag();
	}
	four1(m_z,m_N,-1);
	//!  adjust amplitude ... 
	for(int i=0; i<2*m_N; ++i)
	{
		m_z[i] /= m_N/T;
	}
	//! adjust shift 
	if(m_shift != 0.0)
	{
		for(int i=0; i<m_N; ++i)
		{
			// F(f(x)) = F(f(x-a))*exp(-2*pi*m_shfit*m_xi[i])
            double r_   = m_z[2*i];
			double i_   = m_z[2*i+1];   
			double cos_ = cos(2*pi*m_shift*m_xi[i]);
			double sin_ = sin(2*pi*m_shift*m_xi[i]);		

			m_z[2*i]    = r_*cos_ - i_*sin_;
			m_z[2*i+1]  = r_*sin_ + i_*cos_;
		}
	}
}

void FFT :: calculate_inverse_FFT() //!  size of fft_result = 2*N
{
	four1(m_z,m_N,1);

	double T = (m_right-m_left);
	for(int i=0; i<2*m_N; ++i)
	{
		m_z[i] /= T;
	}
}

void FFT :: print_x_f(string file_name)  // x,f are of size N
{
	ofstream my_file;
	my_file.open(file_name.c_str());

	my_file.precision(10);
    for(int i=0; i<m_N; ++i)
	{
		my_file << m_x[i] << "," << m_f[i].real() << "," << m_f[i].imag() << endl;
	}
	my_file.close();
}

void FFT :: print_fft_result(string file_name)  // fft_result is of size 2*N
{
	ofstream my_file;
	my_file.open(file_name.c_str());
	my_file.precision(10);

	for(int i=0; i<m_N; ++i)
	{
		my_file << m_xi[i] << ","<< m_z[2*i] << "," << m_z[2*i+1] << "," ;
		my_file << sqrt(2*pi)*exp(-2*pi*pi*m_xi[i]*m_xi[i]) << "," << 0 << endl;
	}
	my_file.close();
}

void FFT :: print_inverse_fft_result(string file_name)  // fft_result is of size 2*N
{
	ofstream my_file;
	my_file.open(file_name.c_str());
	my_file.precision(10);

	for(int i=0; i<m_N; ++i)
	{
		my_file << m_x[i] << ","<< m_z[2*i] << "," << m_z[2*i+1] << endl;
	}
	my_file.close();
}
