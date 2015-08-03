#pragma once

#include <PDE/NR/code/nr3.h>
#include <cmath>
#include <vector>
#include <fstream>
#include <utility>
#include <complex>
#include <fstream>
#include <string>

class FFT
{
public:
	double                  m_left;    
	double                  m_right;  //!  function at [left, right] 
	double                  m_shift;  //!  f(x) is shifted to f(x-a), to put f in the center of the window of FFT!

	int                     m_n ;
	int                     m_N;      // N = 2^n 

	std::vector<double>          m_x;      // time discretization   
	std::vector<std::complex<double>> m_f;      // f(x)

	std::vector<double>          m_xi;     // spectrum discretization
	//vector<complex<double>> m_fft_f;  //  fft result of f

	//! working place
	double* m_z;

	typedef std::complex<double> (*Func)(double);
	FFT(double border, int n, Func f, double shift);  //!  the interval must be symetric centered at 0: [-border, border];
	                            //!  if the most signifigent par of the function is far aways from 0, 
	                                 // you need to translate it to the center: F(f(x-a)) = exp{-2*pi*i*a*\xi}*F(f(x)) 
	~FFT();
	void set_x();
	void set_f(Func f);
	void set_xi();     
	void calculate_FFT();
	void calculate_inverse_FFT();
	void get_ordered_fft_result_from_z();  // get ordered result to m_fft_f/m_inverse_fft_f, from z;
	void print_x_f(std::string file_name);
	void print_fft_result(std::string file_name);
	void print_inverse_fft_result(std::string file_name);
};
