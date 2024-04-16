/*

Date:     7 July 2016

Function: Myelin segmentation without PV modelling

*/

#ifndef _IRTKMYELINCLASSIFICATION_H

#define _IRTKMYELINCLASSIFICATION_H

#include <irtkImage.h>

#include <irtkProbabilityMaps.h>

#include <irtkGaussianLikelihood.h>

class irtkMyelinClassification: public irtkObject
{

protected:

  /// Input image
  irtkRealImage _image;
  
  /// Mask of ROI
  irtkRealImage _mask;
  
  /// Initial myelin segmentation
  irtkRealImage _init;
  
  /// Myelin segmentation by maximum vote
  irtkRealImage _hard_myelin;
  
  /// Posterior probability maps
  irtkProbabilityMaps _posteriors;
  
  /// MRF prior probability maps
  irtkProbabilityMaps _MRFpriors;
  
  /// Image dimensions
  int _Dx, _Dy, _Dz;
  
  /// Voxel dimensions
  double _dx, _dy, _dz;
  
  /// Total number of voxels
  int _I;  
  
  /// Total number of GMM classes
  int _K;
  
  /// GMM parameters
  double *_mn;
  double *_vr;
  double *_mx;

  /// 2D connectivity matrix for first-order MRFs
  irtkMatrix _M;
  
  /// Parameter of M
  double _m1;  
  bool _mrf1;

  /// Objective function of negative log likelihood to be minimised
  double _f;
  
  /// E-step denominator saved to compute objective function
  irtkRealImage _denom;

public:

  /// Constructor
  irtkMyelinClassification(const irtkRealImage &image, const irtkRealImage &mask, const irtkRealImage &init, double m1);
  
  /// Destructor
  ~irtkMyelinClassification();
  
  /// Construct 2D connectivity matrix for first-order MRFs
  void Construct2DConnectivityMatrix();
  
  /// Initialise posterior probability maps
  void InitialisePosteriors();
  
  /// Initialise MRF prior probability maps
  void InitialisePriors();

  /// Initialise GMM parameters
  void InitialiseParams();

  /// Return sum of neighbour probabilities weighted respectively by inverse Euclidean distances at voxel (x,y,z) of class k 
  irtkRealPixel GetNeighbourhood(int x, int y, int z, int k);
  
  /// Return first-order penalty at voxel (x,y,z) of class k
  irtkRealPixel GetFirstOrderPenalty(int x, int y, int z, int k);
  
  /// EM algorithm
  void EStep();
  void MStep();
  
  /// Return relative decrease of objective function
  double LogLikelihood();
    
  /// Print GMM parameters to screen
  void PrintToScreen();
   
  /// Print GMM parameters to file
  void PrintToFile();
  
  /// Execute one EM iteration and return relative decrease of objective function
  double Iterate();
      
  /// Calculate final myelin segmentation
  void GetMyelinSegment();
  
  /// Write posterior probability maps
  void WritePosteriors();
        
  /// Print progress
  void PrintProgress(int x);  
};

#endif
