#include "kernel_gen.h"

#include <cmath>
#include <iostream>

KernelGen::KernelGen() {}

/*
  sigma = 0.85, se maggior allora la campana si allarga
  se usassi una gaussiana asimmetrica centrata?
*/
KernelGen::KernelGen(double sigma, int lValue)
{
  kernel.resize(lValue);
  gauss_sigma = sigma;
  norm_term = 0.0;

  double den = 2.0 * pow(gauss_sigma, 2);
  for (int i = 0; i < lValue; i++)
  {
    kernel[lValue - 1 - i] = exp(- pow(i, 2) / den);
    norm_term += kernel[lValue - 1 - i];
    std::cout << "kernel[lValue - 1 - i]" << kernel[lValue - 1 - i] << std::endl;
  }

  std::cout << "norm_term " << norm_term << std::endl;
}

std::vector<double> & 
KernelGen::getKernel()
{
  return kernel;
}