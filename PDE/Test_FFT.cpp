#include "FFT.h"
#include <PDE/Params.h>

complex<double> gaussian_f(double x)
{ 
	double y = x;
	return complex<double>(exp(-0.5*y*y),0.0);
}


void test_FFT()
{
	int    n       = 10;    // N = 2^n
	int    N       = (int) pow(2.0,n);
	int    N2      = 2*N;
	double border  = 6;  // for interval [-border, border];
	double shift   = 0.0;

	string file_name_x_f           = DEBUG_output_path + "x(non-ordered)_f.csv";
	string file_name_fft_f         = DEBUG_output_path + "fft_f.csv";
	string file_name_inverse_fft_f = DEBUG_output_path + "inverse_fft_f.csv";

	FFT fft_obj(border, n, gaussian_f, shift); // original input: f(x)
	fft_obj.print_x_f(file_name_x_f);

	fft_obj.calculate_FFT(); // after fft: F(t) 
	fft_obj.print_fft_result(file_name_fft_f);

	fft_obj.calculate_inverse_FFT(); // after inverse_fft, should get the original f(x)
	fft_obj.print_inverse_fft_result(file_name_inverse_fft_f);
}

