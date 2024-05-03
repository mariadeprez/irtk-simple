/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
  Visual Information Processing (VIP), 2011 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

  =========================================================================*/

#ifndef _irtkReconstructionDTI_H

#define _irtkReconstructionDTI_H

#include <irtkReconstruction.h>
//#include <irtkReconstructionb0.h>
#include <irtkSphericalHarmonics.h>


#include <vector>
using namespace std;

/*

  Reconstruction of volume from 2D slices

*/

class irtkReconstructionDTI : public irtkReconstruction
{

 private:
 vector<vector<double> > _directions;
 vector<double> _bvalues;
 irtkRealImage _simulated_signal;
 irtkRealImage _SH_coeffs;
 int _order;
 int _coeffNum;
 irtkMatrix _dirs;
 irtkSphericalHarmonics _sh;
 irtkSphericalHarmonics _sh_vol;

 
 double _motion_sigma;
 
 double _lambdaLB;
 vector<int> _slice_order;
 
 public:
   irtkReconstructionDTI();
   ~irtkReconstructionDTI();
   void GaussianReconstruction4D(int nStacks, vector<irtkRigidTransformation> &stack_transformations, int order);
   void GaussianReconstruction4D2(int nStacks);
   void GaussianReconstruction4D3();
   void Init4D(int nStacks);
   void Init4DGauss(int nStacks, vector<irtkRigidTransformation> &stack_transformations, int order);
   void Init4DGaussWeighted(int nStacks, vector<irtkRigidTransformation> &stack_transformations, int order);
   void RotateDirections(double &dx, double &dy, double &dz, int i);
   void RotateDirections(double &dx, double &dy, double &dz, irtkRigidTransformation &t);
   void CreateSliceDirections(vector< vector<double> >& directions, vector<double>& bvalues);
   void InitSH(irtkMatrix dirs,int order);
   void InitSHT(irtkMatrix dirs,int order);
   void SimulateSlicesDTI();
   void SimulateStacksDTI(vector<irtkRealImage>& stacks, bool simulate_excluded=false);
   void SimulateStacksDTIIntensityMatching(vector<irtkRealImage>& stacks, bool simulate_excluded=false);
   void SuperresolutionDTI(int iter, bool tv = false, double sh_alpha = 5, double regDTI = 0);
   double LaplacianSmoothnessDTI();
   void SaveSHcoeffs(int iteration);
   void SimulateSignal(char *output_name=NULL);
   void SimulateSignal(int iter);
   void SetLambdaLB(double lambda);
   void SetSliceOrder(vector<int> slice_order);
   void PostProcessMotionParameters(irtkRealImage target, bool empty=false);
   void PostProcessMotionParameters2(irtkRealImage target, bool empty=false);
   void PostProcessMotionParametersHS(irtkRealImage target, bool empty=false);
   void PostProcessMotionParametersHS2(irtkRealImage target, bool empty=false);
   void SetMotionSigma(double sigma);
   void SliceToVolumeRegistrationSH();
   void SetSimulatedSignal(irtkRealImage signal);
   void CorrectStackIntensities(vector<irtkRealImage>& stacks);
   
   void NormaliseBiasSH(int iter);
   void NormaliseBiasDTI(int iter, vector<irtkRigidTransformation> &stack_transformations, int order);
   void NormaliseBiasSH2(int iter, vector<irtkRigidTransformation> &stack_transformations, int order);
   
   inline void WriteSimulatedSignal(char * name);
   
   inline void SetSH(irtkRealImage sh);

   double ConsistencyDTI();

   
   friend class ParallelSimulateSlicesDTI;
   friend class ParallelSuperresolutionDTI;
   friend class ParallelSliceToVolumeRegistrationSH;
   friend class ParallelNormaliseBiasDTI;
   friend class ParallelNormaliseBiasDTI2;

};

inline void irtkReconstructionDTI::SetLambdaLB(double lambda)
{
  _lambdaLB=lambda;
}

inline void irtkReconstructionDTI::SetSliceOrder(vector<int> slice_order)
{
  _slice_order=slice_order;
}

inline void irtkReconstructionDTI::SetMotionSigma(double sigma)
{
  _motion_sigma=sigma;
}

inline void irtkReconstructionDTI::SetSimulatedSignal(irtkRealImage signal)
{
  _simulated_signal=signal;
}

inline void irtkReconstructionDTI::WriteSimulatedSignal(char * name)
{
  _simulated_signal.Write(name);
}


inline void irtkReconstructionDTI::SetSH(irtkRealImage sh)
{
  _SH_coeffs=sh;
}



#endif
