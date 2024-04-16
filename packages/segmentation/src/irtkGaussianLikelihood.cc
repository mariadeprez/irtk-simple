/*

Date: 22 Aug 2013
Function: Compute Gaussian likelihood

*/

#include <irtkImage.h>

#include <irtkGaussianLikelihood.h>

/// Assign values of mean and variance
void irtkGaussianLikelihood::Initialise(const double &mn, const double &vr)
{
  _mn = mn;
  _vr = vr;
}

/// Evaluate Gaussian likelihood given intensity of voxel
double irtkGaussianLikelihood::Evaluate(irtkRealPixel intensity)
{
  const double pi = 3.14159265359;
  return (1 / sqrt(2*pi*_vr) * exp(- (intensity - _mn) * (intensity - _mn) / (2*_vr)));
}
