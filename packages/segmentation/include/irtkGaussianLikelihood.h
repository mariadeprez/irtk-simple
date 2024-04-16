/*

Date: 22 Aug 2013
Function: Compute Gaussian likelihood

*/

#ifndef _IRTKGAUSSIANLIKELIHOOD_H

#define _IRTKGAUSSIANLIKELIHOOD_H

#include <irtkImage.h>

class irtkGaussianLikelihood: public irtkObject
{

protected:

  /// Class mean
  double _mn;
  
  /// Class variance
  double _vr;
  
public:

  /// Assign values of mean and variance
  void Initialise(const double &mn, const double &vr);
  
  /// Evaluate Gaussian likelihood given intensity of voxel currently pointed by the image pointer
  double Evaluate(irtkRealPixel intensity);
  
};

#endif
  
  
