/*

Date:     7 July 2016

Function: Myelin segmentation with PV modelling

*/

#ifndef _IRTKMYELINCLASSIFICATIONPV_H

#define _IRTKMYELINCLASSIFICATIONPV_H

#include <irtkImage.h>

#include <irtkProbabilityMaps.h>

#include <irtkGaussianLikelihood.h>

class irtkMyelinClassificationPV: public irtkObject
{

protected:

  /// Input image
  irtkRealImage _image;
  
  /// Mask of ROI
  irtkRealImage _mask;
  
  /// Initial myelin segmentation
  irtkRealImage _init;
  
  /// Soft segmentation as fraction of myelin class
  irtkRealImage _myelin_soft_segment;
  
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
  
  /// Number of PV classes
  int _N;
  
  /// Step size of changing fraction contents of composing tissues in PV voxels
  double _lambda;
  
  /// GMM parameters
  double *_mn;
  double *_vr;
  double *_mx;
 
  /// 3D connectivity tensor for second-order MRFs
  vector<irtkMatrix> _T;

  /// Parameters of T
  double _t1;
  double _t2;
  double _t3;
  bool _mrf2;
      
  /// Objective function of negative log likelihood to be minimised
  double _f;  
  
  /// E-step denominator saved to compute objective function
  irtkRealImage _denom;  
      
public:

  /// Constructor
  irtkMyelinClassificationPV(const irtkRealImage &image, const irtkRealImage &mask, const irtkRealImage &init, double t1, double t2, double t3);
  
  /// Destructor
  ~irtkMyelinClassificationPV();

  /// Construct 3D connectivity tensor for second-order MRFs
  void Construct3DConnectivityTensor(); 
  
  /// Initialise posterior probability maps
  void InitialisePosteriors();

  /// Initialise MRF prior probability maps
  void InitialisePriors();

  /// Initialise GMM parameters
  void InitialiseParams();
  
  /// Return sum of neighbour probabilities weighted respectively by inverse Euclidean distances at voxel (x,y,z) of class k 
  irtkRealPixel GetNeighbourhood(int x, int y, int z, int k);
   
  /// Return second-order penalty at voxel (x,y,z) of class k 
  irtkRealPixel GetSecondOrderPenalty(int x, int y, int z, int k);  
        
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
