#ifndef __KERNELGEN_H_INCLUDED__ 
#define __KERNELGEN_H_INCLUDED__

#include <vector>

class KernelGen
{
public:
	std::vector<double> kernel;
	double gauss_sigma;
	double norm_term;

	KernelGen();

  	KernelGen(double sigma, int lValue);

	std::vector<double> & 
	getKernel();
};

#endif