/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
  Visual Information Processing (VIP), 2011 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

  =========================================================================*/

#include <irtkReconstruction.h>
#include <irtkReconstructionDTI.h>
#include <irtkResampling.h>
#include <irtkRegistration.h>
#include <irtkImageRigidRegistration.h>
#include <irtkImageRigidRegistrationWithPadding.h>
#include <irtkTransformation.h>



irtkReconstructionDTI::irtkReconstructionDTI():irtkReconstruction()
{
    _recon_type = _interpolate;

}

irtkReconstructionDTI::~irtkReconstructionDTI() 
{ 
  _lambdaLB = 0;
}

void irtkReconstructionDTI::Init4D(int nStacks)
{
  irtkImageAttributes attr = _reconstructed.GetImageAttributes();
  attr._t = nStacks;
  cout<<nStacks<<" stacks."<<endl;
  irtkRealImage recon4D(attr); 
 
  
  for(int t=0; t<recon4D.GetT();t++)
    for(int k=0; k<recon4D.GetZ();k++)
      for(int j=0; j<recon4D.GetY();j++)
        for(int i=0; i<recon4D.GetX();i++)
        {
	  recon4D(i,j,k,t)=_reconstructed(i,j,k);
        }
        
  //recon4D.Write("recon4D.nii.gz");
  _simulated_signal = recon4D;
}

void irtkReconstructionDTI::Init4DGauss(int nStacks, vector<irtkRigidTransformation> &stack_transformations, int order)
{
    unsigned int inputIndex;
    int i, j, k, t, n;
    irtkRealImage slice;
    double scale;
    POINT3D p;
  
  irtkImageAttributes attr = _reconstructed.GetImageAttributes();
  attr._t = nStacks;
  cout<<nStacks<<" stacks."<<endl;
  irtkRealImage recon4D(attr); 
  irtkRealImage weights(attr);

  for (inputIndex = 0; inputIndex < _slices.size(); ++inputIndex) {
      //copy the current slice
      slice = _slices[inputIndex];
      //alias the current bias image
      irtkRealImage& b = _bias[inputIndex];
      //read current scale factor
      scale = _scale[inputIndex];
      //cout<<scale<<" ";
	
        //Distribute slice intensities to the volume
        for (i = 0; i < slice.GetX(); i++)
            for (j = 0; j < slice.GetY(); j++)
                if (slice(i, j, 0) != -1) {
                    
		    //biascorrect and scale the slice
		    //if(origDir==1)
                      slice(i, j, 0) *= exp(-b(i, j, 0)) * scale;

                    //number of volume voxels with non-zero coefficients
                    //for current slice voxel
                    n = _volcoeffs[inputIndex][i][j].size();

                    //add contribution of current slice voxel to all voxel volumes
                    //to which it contributes
                    for (k = 0; k < n; k++) {
                        p = _volcoeffs[inputIndex][i][j][k];
                        recon4D(p.x, p.y, p.z, _stack_index[inputIndex]) += _slice_weight[inputIndex] * p.value * slice(i, j, 0);
                        weights(p.x, p.y, p.z, _stack_index[inputIndex]) += _slice_weight[inputIndex] * p.value;
                    }
                }
      //} //end of loop for origDir
	
    }

    //normalize the volume by proportion of contributing slice voxels
    //for each volume voxe
    recon4D.Write("recon4D.nii.gz");
    recon4D /= weights;
    recon4D.Write("recon4DGaussian.nii.gz");
   
    // rotate directions
    int dirIndex;
    double gx,gy,gz;
    irtkMatrix dirs_xyz(nStacks,3);   
    for (dirIndex = 0; dirIndex < nStacks; dirIndex++)  
    {
      gx=_directions[0][dirIndex+1];
      gy=_directions[1][dirIndex+1];
      gz=_directions[2][dirIndex+1];
      RotateDirections(gx,gy,gz,stack_transformations[dirIndex]);
      dirs_xyz(dirIndex,0) = gx;
      dirs_xyz(dirIndex,1) = gy;
      dirs_xyz(dirIndex,2) = gz;     
    }
    irtkSphericalHarmonics sh;
    sh.InitSHT(dirs_xyz,order);
    _SH_coeffs = sh.Signal2Coeff(recon4D);
    _SH_coeffs.Write("_SH_coeffs.nii.gz");
    SimulateSignal();
}

void irtkReconstructionDTI::Init4DGaussWeighted(int nStacks, vector<irtkRigidTransformation> &stack_transformations, int order)
{
    unsigned int inputIndex;
    int i, j, k, t, n;
    irtkRealImage slice;
    double scale;
    POINT3D p;
  
  irtkImageAttributes attr = _reconstructed.GetImageAttributes();
  attr._t = nStacks;
  cout<<nStacks<<" stacks."<<endl;
  irtkRealImage recon4D(attr); 
  irtkRealImage weights(attr);

  for (inputIndex = 0; inputIndex < _slices.size(); ++inputIndex) {
      //copy the current slice
      slice = _slices[inputIndex];
      //alias the current bias image
      irtkRealImage& b = _bias[inputIndex];
      //read current scale factor
      scale = _scale[inputIndex];
      //cout<<scale<<" ";
	
        //Distribute slice intensities to the volume
        for (i = 0; i < slice.GetX(); i++)
            for (j = 0; j < slice.GetY(); j++)
                if (slice(i, j, 0) != -1) {
                    
		    //biascorrect and scale the slice
		    //if(origDir==1)
                      slice(i, j, 0) *= exp(-b(i, j, 0)) * scale;

                    //number of volume voxels with non-zero coefficients
                    //for current slice voxel
                    n = _volcoeffs[inputIndex][i][j].size();

                    //add contribution of current recon-test4D-gauss-weightedslice voxel to all voxel volumes
                    //to which it contributes
                    for (k = 0; k < n; k++) {
                        p = _volcoeffs[inputIndex][i][j][k];
                        recon4D(p.x, p.y, p.z, _stack_index[inputIndex]) += _slice_weight[inputIndex] * p.value * slice(i, j, 0);
                        weights(p.x, p.y, p.z, _stack_index[inputIndex]) += _slice_weight[inputIndex] * p.value;
                    }
                }
      //} //end of loop for origDir
	
    }

    //normalize the volume by proportion of contributing slice voxels
    //for each volume voxe
    recon4D.Write("recon4D.nii.gz");
    recon4D /= weights;
    recon4D.Write("recon4DGaussian.nii.gz");
   
    // rotate directions
    int dirIndex;
    double gx,gy,gz;
    irtkMatrix dirs_xyz(nStacks,3);   
    for (dirIndex = 0; dirIndex < nStacks; dirIndex++)  
    {
      gx=_directions[0][dirIndex+1];
      gy=_directions[1][dirIndex+1];
      gz=_directions[2][dirIndex+1];
      RotateDirections(gx,gy,gz,stack_transformations[dirIndex]);
      dirs_xyz(dirIndex,0) = gx;
      dirs_xyz(dirIndex,1) = gy;
      dirs_xyz(dirIndex,2) = gz;     
    }
    irtkSphericalHarmonics sh;
    //sh.InitSHT(dirs_xyz,order);
    //_SH_coeffs = sh.Signal2Coeff(recon4D);
    //_SH_coeffs.Write("_SH_coeffs.nii.gz");
    //SimulateSignal();
    
    irtkMatrix SHT = sh.SHbasis(dirs_xyz,order);
    irtkMatrix tSHT = SHT;
    tSHT.Transpose();
 
    attr._t=SHT.Cols();
    irtkRealImage shc(attr);
    irtkVector s(nStacks),c,ev;
    irtkMatrix w(nStacks,nStacks);
    irtkMatrix tmp;
    irtkMatrix tmp1,tmp2;
    irtkRealImage mask = _reconstructed;
    mask=0;
    
    double det;
    for(i=0;i<shc.GetX();i++)
      for(j=0;j<shc.GetY();j++)
        for(k=0;k<shc.GetZ();k++)
	{
	  for(t=0;t<nStacks;t++)
	  {
	    s(t)=recon4D(i,j,k,t);
 	    w(t,t)=weights(i,j,k,t);
	  }

	  tmp = tSHT*w*SHT;
	  tmp.Eigenvalues(tmp1,ev,tmp2);
	  
	  det=1;
	  for(t=0;t<ev.Rows();t++)
	    det*=ev(t);
	  //cerr<<det<<" ";
	  if(det>0.001)
	  {
	    mask(i,j,k)=1;
	    tmp.Invert();
	    c = tmp*tSHT*w*s;
	    for(t=0;t<SHT.Cols();t++)
	    {
	      shc(i,j,k,t)=c(t);
	    }
	  }
	}

    _SH_coeffs = shc;
    _SH_coeffs.Write("_SH_coeffs.nii.gz");
    mask.Write("shmask4D.nii.gz");   
    SimulateSignal();

}


void irtkReconstructionDTI::GaussianReconstruction4D(int nStacks, vector<irtkRigidTransformation> &stack_transformations, int order)
{
    cout << "Gaussian reconstruction ... ";
    unsigned int inputIndex;
    int i, j, k, t, n;
    irtkRealImage slice;
    double scale;
    POINT3D p;
    int dirIndex;
    double gx,gy,gz;
    irtkImageAttributes attr = _reconstructed.GetImageAttributes();
    
    
    //create SH matrix
    irtkMatrix dirs_xyz(nStacks,3);   
    irtkMatrix SHT;
    for (dirIndex = 0; dirIndex < nStacks; dirIndex++)  
    {
      gx=_directions[0][dirIndex+1];
      gy=_directions[1][dirIndex+1];
      gz=_directions[2][dirIndex+1];
      RotateDirections(gx,gy,gz,stack_transformations[dirIndex]);
      dirs_xyz(dirIndex,0) = gx;
      dirs_xyz(dirIndex,1) = gy;
      dirs_xyz(dirIndex,2) = gz;     
    }
    irtkSphericalHarmonics sh;
    sh.InitSHT(dirs_xyz,order);
    SHT = sh.SHbasis(dirs_xyz,order);
    irtkMatrix tSHT = SHT;
    tSHT.Transpose();
    SHT.Print();
    tSHT.Print();
    attr._t=SHT.Cols();
    irtkRealImage shc(attr);
    //shc.Write("shc.nii.gz");



    // Calculate estimate of the signal
    
    attr._t = nStacks;
    cout<<nStacks<<" stacks."<<endl;
    irtkRealImage recon4D(attr);
    irtkRealImage weights(attr);

    for (inputIndex = 0; inputIndex < _slices.size(); ++inputIndex) {
      //copy the current slice
      slice = _slices[inputIndex];
      //alias the current bias image
      irtkRealImage& b = _bias[inputIndex];
      //read current scale factor
      scale = _scale[inputIndex];
      //cout<<scale<<" ";
	
        //Distribute slice intensities to the volume
        for (i = 0; i < slice.GetX(); i++)
            for (j = 0; j < slice.GetY(); j++)
                if (slice(i, j, 0) != -1) {
                    
		    //biascorrect and scale the slice
		    //if(origDir==1)
                      slice(i, j, 0) *= exp(-b(i, j, 0)) * scale;

                    //number of volume voxels with non-zero coefficients
                    //for current slice voxel
                    n = _volcoeffs[inputIndex][i][j].size();

                    //add contribution of current slice voxel to all voxel volumes
                    //to which it contributes
                    for (k = 0; k < n; k++) {
                        p = _volcoeffs[inputIndex][i][j][k];
                        recon4D(p.x, p.y, p.z, _stack_index[inputIndex]) += _slice_weight[inputIndex] * p.value * slice(i, j, 0);
                        weights(p.x, p.y, p.z, _stack_index[inputIndex]) += _slice_weight[inputIndex] * p.value;
                    }
                }
      //} //end of loop for origDir
	
    }

    //normalize the volume by proportion of contributing slice voxels
    //for each volume voxe
   recon4D.Write("recon4D.nii.gz");
   recon4D /= weights;
   _simulated_signal = recon4D;
   irtkRealImage mask = recon4D;
   mask=0;
   irtkRealImage  _test_sh_coeffs = sh.Signal2Coeff(recon4D);
   _test_sh_coeffs.Write("_test_sh_coeffs.nii.gz");
        
    cout << "done." << endl;

    if (_debug)
    {
        recon4D.Write("recon4Dgaussian.nii.gz");
        weights.Write("weights.nii.gz");
    }
    
    //Calculate SH coeffs using the estimated signal and weights
    irtkVector s(nStacks),c,ev;
    irtkMatrix w(nStacks,nStacks);
    irtkMatrix tmp;
    irtkMatrix tmp1,tmp2;
    
    double det;
    //cout<<"SHT:"<<endl;
    //SHT.Print();
    for(i=0;i<shc.GetX();i++)
      for(j=0;j<shc.GetY();j++)
        for(k=0;k<shc.GetZ();k++)
	{
	  for(t=0;t<nStacks;t++)
	  {
	    s(t)=recon4D(i,j,k,t);
 	    w(t,t)=weights(i,j,k,t);
	  }

	  tmp = tSHT*w*SHT;
	  tmp.Eigenvalues(tmp1,ev,tmp2);
	  
	  det=1;
	  for(t=0;t<ev.Rows();t++)
	    det*=ev(t);
	  //cerr<<det<<" ";
	  if(det>0.001)
	  {
	    mask(i,j,k)=1;
	    tmp.Invert();
	    c = tmp*tSHT*w*s;
	    for(t=0;t<SHT.Cols();t++)
	    {
	      shc(i,j,k,t)=c(t);
	    }
	  }
	}

    _SH_coeffs = shc;
    _SH_coeffs.Write("_SH_coeffs.nii.gz");
    mask.Write("shmask4D.nii.gz");
}

void irtkReconstructionDTI::GaussianReconstruction4D2(int nStacks)
{
    cout << "Gaussian reconstruction ... ";
    unsigned int inputIndex;
    int i, j, k, n;
    irtkRealImage slice;
    double scale;
    POINT3D p;
    int dirIndex, origDir;
    double bval,gx,gy,gz,dx,dy,dz,dotp,sigma=0.02,w,tw;

    irtkImageAttributes attr = _reconstructed.GetImageAttributes();
    attr._t = nStacks;
    cout<<nStacks<<" stacks."<<endl;
    irtkRealImage recon4D(attr);
    irtkRealImage weights(attr);

    for (inputIndex = 0; inputIndex < _slices.size(); ++inputIndex) {
      //copy the current slice
      slice = _slices[inputIndex];
      //alias the current bias image
      irtkRealImage& b = _bias[inputIndex];
      //read current scale factor
      scale = _scale[inputIndex];
      cout<<scale<<" ";
        
      //direction for current slice
      dirIndex = _stack_index[inputIndex]+1;
      gx=_directions[0][dirIndex];
      gy=_directions[1][dirIndex];
      gz=_directions[2][dirIndex];
      RotateDirections(gx,gy,gz,inputIndex);
      bval=_bvalues[dirIndex];
	
      for (origDir=1; origDir<_bvalues.size();origDir++)
      {
	//cout<<origDir<<" ";
	
	  dx=_directions[0][origDir];
	  dy=_directions[1][origDir];
	  dz=_directions[2][origDir];
	  
	  dotp = (dx*gx+dy*gy+dz*gz)/sqrt((dx*dx+dy*dy+dz*dz)*(gx*gx+gy*gy+gz*gz));
	  //cout<<"origDir="<<origDir<<": "<<dotp<<"; ";
	  tw = (fabs(dotp)-1)/sigma;
	  w=exp(-tw*tw)/(6.28*sigma);
	  //cout<<"weight = "<<w<<"; ";

	
        //Distribute slice intensities to the volume
        for (i = 0; i < slice.GetX(); i++)
            for (j = 0; j < slice.GetY(); j++)
                if (slice(i, j, 0) != -1) {
                    
		    //biascorrect and scale the slice
		    if(origDir==1)
                      slice(i, j, 0) *= exp(-b(i, j, 0)) * scale;

                    //number of volume voxels with non-zero coefficients
                    //for current slice voxel
                    n = _volcoeffs[inputIndex][i][j].size();

                    //add contribution of current slice voxel to all voxel volumes
                    //to which it contributes
                    for (k = 0; k < n; k++) {
                        p = _volcoeffs[inputIndex][i][j][k];
                        recon4D(p.x, p.y, p.z,origDir-1) += _slice_weight[inputIndex] * w * p.value * slice(i, j, 0);
                        weights(p.x, p.y, p.z,origDir-1) += _slice_weight[inputIndex] * w * p.value;
                    }
                }
      } //end of loop for origDir
      //cout<<endl;
        //end of loop for a slice inputIndex
        //sprintf(buffer, "corslice%i.nii.gz", inputIndex);
	//slice
	
    }

    //normalize the volume by proportion of contributing slice voxels
    //for each volume voxe
   recon4D.Write("recon4D.nii.gz");
   recon4D /= weights;
   _simulated_signal = recon4D;
        
    cout << "done." << endl;

    if (_debug)
    {
        recon4D.Write("recon4Dgaussian.nii.gz");
        weights.Write("weights.nii.gz");
    }

}

void irtkReconstructionDTI::GaussianReconstruction4D3()
{
    cout << "Gaussian reconstruction SH ... ";
    unsigned int inputIndex;
    int i, j, k, n;
    irtkRealImage slice;
    double scale;
    POINT3D p;
    int dirIndex;
    double bval,gx,gy,gz;

    irtkImageAttributes attr = _reconstructed.GetImageAttributes();
    attr._t = _coeffNum;
    cout<<_coeffNum<<" SH coefficients."<<endl;
    irtkRealImage recon4D(attr);
    irtkRealImage weights(attr);

    for (inputIndex = 0; inputIndex < _slices.size(); ++inputIndex) {
      //copy the current slice
      slice = _slices[inputIndex];
      //alias the current bias image
      irtkRealImage& b = _bias[inputIndex];
      //read current scale factor
      scale = _scale[inputIndex];
      //cout<<scale<<" ";
        
      //direction for current slice
      dirIndex = _stack_index[inputIndex]+1;
      gx=_directions[0][dirIndex];
      gy=_directions[1][dirIndex];
      gz=_directions[2][dirIndex];
      RotateDirections(gx,gy,gz,inputIndex);
      bval=_bvalues[dirIndex];
      
      irtkSphericalHarmonics sh;
      irtkMatrix dir(1,3);
      dir(0,0)=gx;
      dir(0,1)=gy;
      dir(0,2)=gz;
      irtkMatrix basis = sh.SHbasis(dir,_order);
      if(basis.Cols() != recon4D.GetT())
      {
	cerr<<"GaussianReconstruction4D3:basis numbers does not match SH coefficients number."<<endl;
	exit(1);
	
      }
	
        //Distribute slice intensities to the volume
        for (i = 0; i < slice.GetX(); i++)
            for (j = 0; j < slice.GetY(); j++)
                if (slice(i, j, 0) != -1) {
                    
		    //biascorrect and scale the slice
                      slice(i, j, 0) *= exp(-b(i, j, 0)) * scale;

                    //number of volume voxels with non-zero coefficients
                    //for current slice voxel
                    n = _volcoeffs[inputIndex][i][j].size();

                    //add contribution of current slice voxel to all voxel volumes
                    //to which it contributes
                    for (k = 0; k < n; k++) {
                      p = _volcoeffs[inputIndex][i][j][k];
                      for(unsigned int l = 0; l < basis.Cols(); l++ )
		      {
			if(l==0)
			{
                          recon4D(p.x, p.y, p.z,l) += basis(0,l) *_slice_weight[inputIndex] * p.value * slice(i, j, 0);
                          weights(p.x, p.y, p.z,l) += basis(0,l) * _slice_weight[inputIndex] * p.value;
			}
		      }
                    }
                }
    //end of loop for a slice inputIndex
    }

    //normalize the volume by proportion of contributing slice voxels
    //for each volume voxe
   recon4D.Write("recon4D.nii.gz");
   recon4D /= weights;
   _SH_coeffs = recon4D;
        
    cout << "done." << endl;

    if (_debug)
    {
        recon4D.Write("recon4DgaussianSH.nii.gz");
        weights.Write("weights.nii.gz");
    }

}

void irtkReconstructionDTI::SimulateStacksDTI(vector<irtkRealImage>& stacks, bool simulate_excluded)
{
    if (_debug)
        cout<<"Simulating stacks."<<endl;
    cout.flush();
  
    unsigned int inputIndex;
    int i, j, k, n;
    irtkRealImage sim;
    POINT3D p;
    double weight;
  
    int z, current_stack;
    z=-1;//this is the z coordinate of the stack
    current_stack=-1; //we need to know when to start a new stack
    
    double threshold = 0.5;
    if (simulate_excluded)
      threshold = -1;

    //_reconstructed.Write("reconstructed.nii.gz");
    for (inputIndex = 0; inputIndex < _slices.size(); inputIndex++) {
      
	cout<<inputIndex<<" ";
	cout.flush();
        // read the current slice
        irtkRealImage& slice = _slices[inputIndex];

        //Calculate simulated slice
        sim.Initialize( slice.GetImageAttributes() );
        sim = 0;
	
	//direction for current slice
        int dirIndex = _stack_index[inputIndex]+1;
        double gx=_directions[0][dirIndex];
        double gy=_directions[1][dirIndex];
        double gz=_directions[2][dirIndex];
        RotateDirections(gx,gy,gz,inputIndex);
        double bval=_bvalues[dirIndex];
	irtkSphericalHarmonics sh;
	irtkMatrix dir(1,3);
	dir(0,0)=gx;
	dir(0,1)=gy;
	dir(0,2)=gz;
	irtkMatrix basis = sh.SHbasis(dir,_order);
	double sim_signal;

	//do not simulate excluded slice
        if(_slice_weight[inputIndex]>threshold)
	{
          for (i = 0; i < slice.GetX(); i++)
            for (j = 0; j < slice.GetY(); j++)
                if (slice(i, j, 0) != -1) {
                    weight=0;
                    n = _volcoeffs[inputIndex][i][j].size();
                    for (k = 0; k < n; k++) {
                        p = _volcoeffs[inputIndex][i][j][k];
			//signal simulated from SH
			sim_signal = 0;
			for(unsigned int l = 0; l < basis.Cols(); l++ )
			  sim_signal += _SH_coeffs(p.x, p.y, p.z,l)*basis(0,l);
			//update slice
                        sim(i, j, 0) += p.value *sim_signal;
                        weight += p.value;
                    }
                    if(weight>0.98)
                        sim(i,j,0)/=weight;
		    else
		      sim(i,j,0)=0;
                }
	}

        if (_stack_index[inputIndex]==current_stack)
            z++;
        else {
            current_stack=_stack_index[inputIndex];
            z=0;
        }
        
        for(i=0;i<sim.GetX();i++)
            for(j=0;j<sim.GetY();j++) {
                stacks[_stack_index[inputIndex]](i,j,z)=sim(i,j,0);
            }
        //end of loop for a slice inputIndex
    }   
}

void irtkReconstructionDTI::SimulateStacksDTIIntensityMatching(vector<irtkRealImage>& stacks, bool simulate_excluded)
{
    if (_debug)
        cout<<"Simulating stacks."<<endl;
    cout.flush();
  
    unsigned int inputIndex;
    int i, j, k, n;
    irtkRealImage sim;
    POINT3D p;
    double weight;
  
    int z, current_stack;
    z=-1;//this is the z coordinate of the stack
    current_stack=-1; //we need to know when to start a new stack
    
    double threshold = 0.5;
    if (simulate_excluded)
      threshold = -1;

    _reconstructed.Write("reconstructed.nii.gz");
    for (inputIndex = 0; inputIndex < _slices.size(); inputIndex++) {
      
	cout<<inputIndex<<" ";
	cout.flush();
        // read the current slice
        irtkRealImage& slice = _slices[inputIndex];
        //read the current bias image
        irtkRealImage& b = _bias[inputIndex];               
        //identify scale factor
        double scale = _scale[inputIndex];


        //Calculate simulated slice
        sim.Initialize( slice.GetImageAttributes() );
        sim = 0;
	irtkRealImage simulatedslice(sim), simulatedsliceint(sim), simulatedweights(sim);
	
	//direction for current slice
        int dirIndex = _stack_index[inputIndex]+1;
        double gx=_directions[0][dirIndex];
        double gy=_directions[1][dirIndex];
        double gz=_directions[2][dirIndex];
        RotateDirections(gx,gy,gz,inputIndex);
        double bval=_bvalues[dirIndex];
	irtkSphericalHarmonics sh;
	irtkMatrix dir(1,3);
	dir(0,0)=gx;
	dir(0,1)=gy;
	dir(0,2)=gz;
	irtkMatrix basis = sh.SHbasis(dir,_order);
	double sim_signal;

	//do not simulate excluded slice
        if(_slice_weight[inputIndex]>threshold)
	{
          for (i = 0; i < slice.GetX(); i++)
            for (j = 0; j < slice.GetY(); j++)
                if (slice(i, j, 0) != -1) {
                    weight=0;
                    n = _volcoeffs[inputIndex][i][j].size();
                    for (k = 0; k < n; k++) {
                        p = _volcoeffs[inputIndex][i][j][k];
			//signal simulated from SH
			sim_signal = 0;
			for(unsigned int l = 0; l < basis.Cols(); l++ )
			  sim_signal += _SH_coeffs(p.x, p.y, p.z,l)*basis(0,l);
			//update slice
                        sim(i, j, 0) += p.value *sim_signal;
                        weight += p.value;
                    }
                    simulatedweights(i,j,0)=weight;
                    if(weight>0.98)
                        sim(i,j,0)/=weight;
		    else
		      sim(i,j,0)=0;
		    simulatedslice(i,j,0)=sim(i,j,0);
		    //intensity parameters
		    if(_intensity_matching_GD)
		    {
		      double a=b(i, j, 0) * scale;
		      if(a>0)
		        sim(i,j,0)/=a;
		      else
			sim(i,j,0)=0;
		    }
		    else
		      sim(i,j,0)/=(exp(-b(i, j, 0)) * scale);
		    simulatedsliceint(i,j,0)=sim(i,j,0);
                }
	}

        if (_stack_index[inputIndex]==current_stack)
            z++;
        else {
            current_stack=_stack_index[inputIndex];
            z=0;
        }
        
        for(i=0;i<sim.GetX();i++)
            for(j=0;j<sim.GetY();j++) {
                stacks[_stack_index[inputIndex]](i,j,z)=sim(i,j,0);
            }
	  
        //end of loop for a slice inputIndex
        
        if (_debug)
	{
          char buffer[256];
	  sprintf(buffer,"simulatedslice%i.nii.gz",inputIndex);
	  simulatedslice.Write(buffer);
	  sprintf(buffer,"simulatedsliceint%i.nii.gz",inputIndex);
	  simulatedsliceint.Write(buffer);
	  sprintf(buffer,"simulatedbias%i.nii.gz",inputIndex);
	  b.Write(buffer);
	  sprintf(buffer,"simulatedweights%i.nii.gz",inputIndex);
	  simulatedweights.Write(buffer);
	}
    }   
}



void irtkReconstructionDTI::CorrectStackIntensities(vector<irtkRealImage>& stacks)
{
  int z, current_stack;
  z=-1;//this is the z coordinate of the stack
  current_stack=-1; //we need to know when to start a new stack
  
  for (int inputIndex = 0; inputIndex < _slices.size(); inputIndex ++)
  {
    if (_stack_index[inputIndex]==current_stack)
      z++;
    else 
    {
      current_stack=_stack_index[inputIndex];
      z=0;
     }
        
     for(int i=0;i<stacks[_stack_index[inputIndex]].GetX();i++)
       for(int j=0;j<stacks[_stack_index[inputIndex]].GetY();j++) 
       {
	 stacks[_stack_index[inputIndex]](i,j,z) *= (exp(-_bias[inputIndex](i, j, 0)) * _scale[inputIndex]);	 
      }
  }
}

class ParallelSimulateSlicesDTI {
    irtkReconstructionDTI *reconstructor;
        
public:
    ParallelSimulateSlicesDTI( irtkReconstructionDTI *_reconstructor ) : 
    reconstructor(_reconstructor) { }

    void operator() (const blocked_range<size_t> &r) const {
        for ( size_t inputIndex = r.begin(); inputIndex != r.end(); ++inputIndex ) {
            //Calculate simulated slice
            reconstructor->_simulated_slices[inputIndex].Initialize( reconstructor->_slices[inputIndex].GetImageAttributes() );
            reconstructor->_simulated_slices[inputIndex] = 0;

            reconstructor->_simulated_weights[inputIndex].Initialize( reconstructor->_slices[inputIndex].GetImageAttributes() );
            reconstructor->_simulated_weights[inputIndex] = 0;

            reconstructor->_simulated_inside[inputIndex].Initialize( reconstructor->_slices[inputIndex].GetImageAttributes() );
            reconstructor->_simulated_inside[inputIndex] = 0;            

            reconstructor->_slice_inside[inputIndex] = false;
	    
	    
	    //direction for current slice
            int dirIndex = reconstructor->_stack_index[inputIndex]+1;
            double gx=reconstructor->_directions[0][dirIndex];
            double gy=reconstructor->_directions[1][dirIndex];
            double gz=reconstructor->_directions[2][dirIndex];
            reconstructor->RotateDirections(gx,gy,gz,inputIndex);
            double bval=reconstructor->_bvalues[dirIndex];
	    irtkSphericalHarmonics sh;
	    irtkMatrix dir(1,3);
	    dir(0,0)=gx;
	    dir(0,1)=gy;
	    dir(0,2)=gz;
	    irtkMatrix basis = sh.SHbasis(dir,reconstructor->_order);
	    if(basis.Cols() != reconstructor->_SH_coeffs.GetT())
	    {
	      cerr<<"ParallelSimulateSlicesDTI:basis numbers does not match SH coefficients number."<<endl;
	      exit(1);
	    }
	    double sim_signal;
            
            POINT3D p;
            for ( unsigned int i = 0; i < reconstructor->_slices[inputIndex].GetX(); i++ )
                for ( unsigned int j = 0; j < reconstructor->_slices[inputIndex].GetY(); j++ )
                    if ( reconstructor->_slices[inputIndex](i, j, 0) != -1 ) {
                        double weight = 0;
                        int n = reconstructor->_volcoeffs[inputIndex][i][j].size();
                        for ( unsigned int k = 0; k < n; k++ ) {
			     //PSF
                            p = reconstructor->_volcoeffs[inputIndex][i][j][k];
			    //signal simulated from SH
			     sim_signal = 0;
			     for(unsigned int l = 0; l < basis.Cols(); l++ )
			       sim_signal += reconstructor->_SH_coeffs(p.x, p.y, p.z,l)*basis(0,l);
			     //update slice
                            reconstructor->_simulated_slices[inputIndex](i, j, 0) += p.value * sim_signal;
                            weight += p.value;
                            if (reconstructor->_mask(p.x, p.y, p.z) == 1) {
                                reconstructor->_simulated_inside[inputIndex](i, j, 0) = 1;
                                reconstructor->_slice_inside[inputIndex] = true;
                            }
                        }                    
                        if( weight > 0 ) {
                            reconstructor->_simulated_slices[inputIndex](i,j,0) /= weight;
                            reconstructor->_simulated_weights[inputIndex](i,j,0) = weight;
                        }
                    }
            
        }
    }
    
    // execute
    void operator() () const {
        task_scheduler_init init(tbb_no_threads);
        parallel_for( blocked_range<size_t>(0, reconstructor->_slices.size() ),
                      *this );
        init.terminate();
    }

};

void irtkReconstructionDTI::SimulateSlicesDTI()
{
    if (_debug)
        cout<<"Simulating slices DTI."<<endl;

    ParallelSimulateSlicesDTI parallelSimulateSlicesDTI( this );
    parallelSimulateSlicesDTI();

    if (_debug)
        cout<<"done."<<endl;    
}

double irtkReconstructionDTI::ConsistencyDTI()
{
  double ssd=0;
  int num = 0;
  double diff;
  for(int index = 0; index< _slices.size(); index++)
  {
    for(int i=0;i<_slices[index].GetX();i++)
      for(int j=0;j<_slices[index].GetY();j++)
	if((_slices[index](i,j,0)>=0)&&(_simulated_inside[index](i,j,0)==1))
	{
	  diff = _slices[index](i,j,0)-_simulated_slices[index](i,j,0);
	  ssd+=diff*diff;
	  num++;
	}
  }
  cout<<"consistency: "<<ssd<<" "<<sqrt(ssd/num)<<endl;
  return ssd;
}

class ParallelSuperresolutionDTI {
    irtkReconstructionDTI* reconstructor;
public:
    irtkRealImage confidence_map;
    irtkRealImage addon;
    
    void operator()( const blocked_range<size_t>& r ) {
        for ( size_t inputIndex = r.begin(); inputIndex < r.end(); ++inputIndex) {
            // read the current slice
            irtkRealImage slice = reconstructor->_slices[inputIndex];
                
            // read the current slice
            irtkRealImage sim = reconstructor->_simulated_slices[inputIndex];
                
            //read the current weight image
            irtkRealImage& w = reconstructor->_weights[inputIndex];
                
            //read the current bias image
            irtkRealImage& b = reconstructor->_bias[inputIndex];
                
            //identify scale factor
            double scale = reconstructor->_scale[inputIndex];
	    if(scale==0)
	    {
	      cerr<<"Scale is 0 for slice "<<inputIndex<<endl;
	      exit(1);
	    }
	    //direction for current slice
            int dirIndex = reconstructor->_stack_index[inputIndex]+1;
            double gx=reconstructor->_directions[0][dirIndex];
            double gy=reconstructor->_directions[1][dirIndex];
            double gz=reconstructor->_directions[2][dirIndex];
            reconstructor->RotateDirections(gx,gy,gz,inputIndex);
            double bval=reconstructor->_bvalues[dirIndex];
	    irtkSphericalHarmonics sh;
	    irtkMatrix dir(1,3);
	    dir(0,0)=gx;
	    dir(0,1)=gy;
	    dir(0,2)=gz;
	    irtkMatrix basis = sh.SHbasis(dir,reconstructor->_order);
	    if(basis.Cols() != reconstructor->_SH_coeffs.GetT())
	    {
	      cerr<<"ParallelSimulateSlicesDTI:basis numbers does not match SH coefficients number."<<endl;
	      exit(1);
	    }

            //Update reconstructed volume using current slice

            //Distribute error to the volume
            POINT3D p;
            for ( int i = 0; i < slice.GetX(); i++)
                for ( int j = 0; j < slice.GetY(); j++)
                    if (slice(i, j, 0) != -1) {
                        //bias correct and scale the slice
		        if(reconstructor->_intensity_matching_GD)
                          sim(i, j, 0) *= b(i, j, 0) / scale;
			else
                          sim(i, j, 0) *= exp(b(i, j, 0)) / scale;
                        
                        if ( reconstructor->_simulated_weights[inputIndex](i,j,0) > 0.5 )
                            slice(i,j,0) -= sim(i,j,0);
                        else
                            slice(i,j,0) = 0;

                        int n = reconstructor->_volcoeffs[inputIndex][i][j].size();
                        for (int k = 0; k < n; k++) {
                            p = reconstructor->_volcoeffs[inputIndex][i][j][k];
			    if(reconstructor->_robust_slices_only)
			    {
			      for(unsigned int l = 0; l < basis.Cols(); l++ )
			      {				
                                addon(p.x, p.y, p.z,l) += p.value * basis(0,l) * slice(i, j, 0) * reconstructor->_slice_weight[inputIndex];
                                confidence_map(p.x, p.y, p.z,l) += p.value *reconstructor->_slice_weight[inputIndex];
			      }
			      
			    }
			    else
			    {
                              for(unsigned int l = 0; l < basis.Cols(); l++ )
			       {
				 if((p.x==20)&&(p.y==50)&&(p.z==40))
				  {
				    //cerr<<inputIndex<<" "<<p.value<<" "<<basis(0,l)<<" "<<slice(i, j, 0)<<" "<<w(i, j, 0)<<" "<<reconstructor->_slice_weight[inputIndex]<<endl;
				    //cerr<<inputIndex<<" "<<exp(b(i, j, 0)) / scale<<" "<<b(i,j,0)<<" "<<scale<<" "<<slice(i, j, 0)<<" "<<w(i, j, 0)<<" "<<reconstructor->_slices[inputIndex](i,j,0)<<" "<<reconstructor->_weights[inputIndex](i,j,0)<<endl;
				  }
				  addon(p.x, p.y, p.z, l) += p.value * basis(0,l) * slice(i, j, 0) * w(i, j, 0) * reconstructor->_slice_weight[inputIndex]* exp(b(i, j, 0)) / scale;
                                confidence_map(p.x, p.y, p.z, l) += p.value * w(i, j, 0) * reconstructor->_slice_weight[inputIndex];
                                //p.value * basis(0,l) * w(i, j, 0) * reconstructor->_slice_weight[inputIndex];
			       }
			    }
                        }
                    }
        } //end of loop for a slice inputIndex
    }
 
    ParallelSuperresolutionDTI( ParallelSuperresolutionDTI& x, split ) :
        reconstructor(x.reconstructor)
    {
        //Clear addon
        addon.Initialize( reconstructor->_SH_coeffs.GetImageAttributes() );
        addon = 0;

        //Clear confidence map
        confidence_map.Initialize( reconstructor->_SH_coeffs.GetImageAttributes() );
        confidence_map = 0;
    }
 
    void join( const ParallelSuperresolutionDTI& y ) {
        addon += y.addon;
        confidence_map += y.confidence_map;
    }
             
    ParallelSuperresolutionDTI( irtkReconstructionDTI *reconstructor ) :
    reconstructor(reconstructor)
    {
        //Clear addon
        addon.Initialize( reconstructor->_SH_coeffs.GetImageAttributes() );
        addon = 0;

        //Clear confidence map
        confidence_map.Initialize( reconstructor->_SH_coeffs.GetImageAttributes() );
        confidence_map = 0;
    }

    // execute
    void operator() () {
        task_scheduler_init init(tbb_no_threads);
        parallel_reduce( blocked_range<size_t>(0,reconstructor->_slices.size()),
                         *this );
        init.terminate();
    }         
};

void irtkReconstructionDTI::SuperresolutionDTI(int iter, bool tv, double sh_alpha, double regDTI)
{
    if (_debug)
        cout << "SuperresolutionDTI " << iter << endl;
    
    int i, j, k, l,t;
    irtkRealImage addon, original;
    

    //Remember current reconstruction for edge-preserving smoothing
    original = _SH_coeffs;

    ParallelSuperresolutionDTI parallelSuperresolutionDTI(this);
    parallelSuperresolutionDTI();
    addon = parallelSuperresolutionDTI.addon;
    _confidence_map = parallelSuperresolutionDTI.confidence_map;
    //_confidence4mask = _confidence_map;
    
    if(_debug) {
        char buffer[256];
        //sprintf(buffer,"confidence-map%i.nii.gz",iter);
        //_confidence_map.Write(buffer);
	_confidence_map.Write("confidence-map-superDTI.nii.gz");
        sprintf(buffer,"addon%i.nii.gz",iter);
        addon.Write(buffer);
    }
/*
    if (!_adaptive) 
        for (i = 0; i < addon.GetX(); i++)
            for (j = 0; j < addon.GetY(); j++)
                for (k = 0; k < addon.GetZ(); k++)
                     for (l = 0; l < addon.GetT(); l++)
                     if (_confidence_map(i, j, k, l) > 0) {
                        // ISSUES if _confidence_map(i, j, k) is too small leading
                        // to bright pixels
                        addon(i, j, k,l) /= _confidence_map(i, j, k, l);
                        //this is to revert to normal (non-adaptive) regularisation
                        _confidence_map(i,j,k,l) = 1;
                    }
                    */

     double alpha;     
     alpha=sh_alpha/_average_volume_weight;
     /*
     if(iter<3)
       alpha = 3;
     else
       if(iter<10)
	 alpha = 1.5;
       else
	 if(iter<30)
	   alpha = 0.5;
	   else
	     alpha=0.25;
	   */
	   
    _SH_coeffs += addon * alpha;// _alpha; //* _average_volume_weight; //_average_volume_weight;
    cout<<"alpha = "<<alpha<<endl;
    
    ////regularisation
    
    if (_lambdaLB > 0)
    {
      cout<<"_lambdaLB = "<<_lambdaLB<<endl;
      irtkRealImage regul(addon);
      regul=0;
      double lb=0;
      double factor;
      for (t = 0; t < regul.GetT(); t++)
        for (i = 0; i < regul.GetX(); i++)
          for (j = 0; j < regul.GetY(); j++)
            for (k = 0; k < regul.GetZ(); k++) 
	    {
	      if(t==0) lb=0;
	      if((t>=1)&&(t<=5)) lb=36.0/1764.0;
	      if((t>=6)&&(t<=14)) lb=400.0/1764.0;
              if((t>=15)&&(t<=27)) lb=1764.0/1764.0;
	      factor = _lambdaLB * lb; 
	      if(factor>1) factor=1;//this is to ensure that there cannot be negative effect
	      regul(i,j,k,t)=factor * _SH_coeffs(i,j,k,t);//alpha * 
	    }
      //regul.Write("regul.nii.gz");
      _SH_coeffs -= regul;
      
    }
    
    //TV on SH basis
    if(tv)
    {
      irtkRealImage o,r;
      
      cout<<"tv="<<tv<<endl;
      for (int steps = 0; steps < _regul_steps; steps ++)
      {
      r=_reconstructed;
      for(t=0;t<_SH_coeffs.GetT();t++)
      {
	o=original.GetRegion(0,0,0,t,original.GetX(),original.GetY(),original.GetZ(),t+1);
	//o.Write("o.nii.gz");
	_reconstructed=_SH_coeffs.GetRegion(0,0,0,t,original.GetX(),original.GetY(),original.GetZ(),t+1);
	//_reconstructed.Write("r.nii.gz");
	
	//L22Regularization(iter, o);
	//AdaptiveRegularization(iter, o);
        // make sure it is unified
        float l = _lambda;
        _lambda = regDTI;
	LaplacianRegularization(iter,t,o);
        _lambda = l;
	
	//o.Write("o1.nii.gz");
	//_reconstructed.Write("r1.nii.gz");
	
	//put it back to _SH_coeffs
        for (i = 0; i < _SH_coeffs.GetX(); i++)
          for (j = 0; j < _SH_coeffs.GetY(); j++)
            for (k = 0; k < _SH_coeffs.GetZ(); k++) 
	    {
	      _SH_coeffs(i,j,k,t)=_reconstructed(i,j,k);
	    }
      }
          _reconstructed=r;
      }
    }

    
  /*  
    //bound the intensities
    for (i = 0; i < _reconstructed.GetX(); i++)
        for (j = 0; j < _reconstructed.GetY(); j++)
            for (k = 0; k < _reconstructed.GetZ(); k++) {
                if (_reconstructed(i, j, k) < _min_intensity * 0.9)
                    _reconstructed(i, j, k) = _min_intensity * 0.9;
                if (_reconstructed(i, j, k) > _max_intensity * 1.1)
                    _reconstructed(i, j, k) = _max_intensity * 1.1;
            }

    //Smooth the reconstructed image
    AdaptiveRegularization(iter, original);
    //Remove the bias in the reconstructed volume compared to previous iteration
    if (_global_bias_correction)
        BiasCorrectVolume(original);
*/
}

void irtkReconstructionDTI::NormaliseBiasSH(int iter)
{
  NormaliseBias(iter, _SH_coeffs);
  /*
  int i,j,k,t;
  irtkRealImage r;
  r=_reconstructed;
  for(t=0;t<_SH_coeffs.GetT();t++)
  {
    _reconstructed=_SH_coeffs.GetRegion(0,0,0,t,_SH_coeffs.GetX(),_SH_coeffs.GetY(),_SH_coeffs.GetZ(),t+1);
    NormaliseBias(iter);
    for (i = 0; i < _SH_coeffs.GetX(); i++)
      for (j = 0; j < _SH_coeffs.GetY(); j++)
	for (k = 0; k < _SH_coeffs.GetZ(); k++) 
	{
	  _SH_coeffs(i,j,k,t)=_reconstructed(i,j,k);
	}
    }
    _reconstructed=r;
    */
}

double irtkReconstructionDTI::LaplacianSmoothnessDTI()
{
    if (_debug)
        cout << "LaplacianSmoothnessDTI "<< endl;
    
    double sum=0;
    irtkRealImage o,r;
    r=_reconstructed;
    for(int t=0;t<_SH_coeffs.GetT();t++)
    {
      o=_SH_coeffs.GetRegion(0,0,0,t,_SH_coeffs.GetX(),_SH_coeffs.GetY(),_SH_coeffs.GetZ(),t+1);
      _reconstructed=_SH_coeffs.GetRegion(0,0,0,t,_SH_coeffs.GetX(),_SH_coeffs.GetY(),_SH_coeffs.GetZ(),t+1);
      sum+=LaplacianSmoothness(o);
      
    }     
    _reconstructed=r;
    return sum;
}





void irtkReconstructionDTI::RotateDirections(double &dx, double &dy, double &dz, int i)
{

  //vector end-point
  double x,y,z;
  //origin
  double ox,oy,oz;
  
  if (_debug)
  {
    //cout<<"Original direction "<<i<<"(dir"<<_stack_index[i]+1<<"): ";
    //cout<<dx<<", "<<dy<<", "<<dz<<". ";
    //cout<<endl;
  }

  //origin
  ox=0;oy=0;oz=0;
  _transformations[i].Transform(ox,oy,oz);
    
  //end-point
  x=dx;
  y=dy;
  z=dz;
  _transformations[i].Transform(x,y,z);
    
  dx=x-ox;
  dy=y-oy;
  dz=z-oz;
    
  if (_debug)
  {
    //cout<<"Rotated direction "<<i<<"(dir"<<_stack_index[i]+1<<"): ";
    //cout<<dx<<", "<<dy<<", "<<dz<<". ";
    //cout<<endl;
  }
}

void irtkReconstructionDTI::RotateDirections(double &dx, double &dy, double &dz, irtkRigidTransformation &t)
{

  //vector end-point
  double x,y,z;
  //origin
  double ox,oy,oz;
  
  if (_debug)
  {
    //cout<<"Original direction "<<i<<"(dir"<<_stack_index[i]+1<<"): ";
    //cout<<dx<<", "<<dy<<", "<<dz<<". ";
    //cout<<endl;
  }

  //origin
  ox=0;oy=0;oz=0;
  t.Transform(ox,oy,oz);
    
  //end-point
  x=dx;
  y=dy;
  z=dz;
  t.Transform(x,y,z);
    
  dx=x-ox;
  dy=y-oy;
  dz=z-oz;
    
  if (_debug)
  {
    //cout<<"Rotated direction "<<i<<"(dir"<<_stack_index[i]+1<<"): ";
    //cout<<dx<<", "<<dy<<", "<<dz<<". ";
    //cout<<endl;
  }
}


void irtkReconstructionDTI::CreateSliceDirections(vector< vector<double> >& directions, vector<double>& bvalues)
{
  _directions = directions;
  _bvalues = bvalues;
  
  cout<<"B-values: ";
  for(uint i=0;i<_bvalues.size(); i++)
    cout<<_bvalues[i]<<" ";
  cout<<endl;

  cout<<"B-vectors: ";
  for(uint j=0; j<_directions.size(); j++)
  {
    for(uint i=0;i<_directions[j].size(); i++)
      cout<<_directions[j][i]<<" ";
    cout<<endl;
  }
  cout<<endl;
}

void irtkReconstructionDTI::InitSH(irtkMatrix dirs, int order)
{
  //irtkSphericalHarmonics sh;
  if(_lambdaLB > 0)
    _sh.InitSHTRegul(dirs,_lambdaLB,order);
  else
    _sh.InitSHT(dirs,order);
  
  _SH_coeffs = _sh.Signal2Coeff(_simulated_signal);
  _SH_coeffs.Write("_SH_coeffs.nii.gz");
  _order = order;
  _dirs = dirs;
  _coeffNum = _sh.NforL(_order);
}

void irtkReconstructionDTI::InitSHT(irtkMatrix dirs, int order)
{
  irtkSphericalHarmonics sh;
  if(_lambdaLB > 0)
    sh.InitSHTRegul(dirs,_lambdaLB,order);
  else
    sh.InitSHT(dirs,order);
  
  //_SH_coeffs = sh.Signal2Coeff(_simulated_signal);
  //_SH_coeffs.Write("_SH_coeffs.nii.gz");
  _order = order;
  _dirs = dirs;
  _coeffNum = sh.NforL(_order);
}

void irtkReconstructionDTI::SimulateSignal()//irtkMatrix dirs, int order)
{
  irtkSphericalHarmonics sh;
  sh.InitSHT(_dirs,_order);//lambda,order);
  _simulated_signal=sh.Coeff2Signal(_SH_coeffs);
  _simulated_signal.Write("_simulated_signal.nii.gz");
  //_SH_coeffs = sh.Signal2Coeff(_simulated_signal);
  //_SH_coeffs.Write("_SH_coeffs.nii.gz");
}

void irtkReconstructionDTI::SimulateSignal(int iter)//irtkMatrix dirs, int order)
{
  irtkSphericalHarmonics sh;
  sh.InitSHT(_dirs,_order);
  irtkRealImage signal=sh.Coeff2Signal(_SH_coeffs);
  char buffer[256];
  sprintf(buffer,"simulated_signal_%i.nii.gz",iter);
  signal.Write(buffer);
}

void irtkReconstructionDTI::SaveSHcoeffs(int iteration)
{
  char buffer[256];
  sprintf(buffer, "shCoeff%i.nii.gz", iteration);
  _SH_coeffs.Write(buffer);
}

void irtkReconstructionDTI::PostProcessMotionParameters(irtkRealImage target, bool empty) //, int packages)
{
  if(_slice_order.size()!=_transformations.size())
  {
    cerr<<"Slice order does not match the number of transformations."<<endl;
    exit(1);
  }
  
  //Reset origin for transformations
  int i;
  irtkGreyImage t = target;
  irtkRigidTransformation offset;
  ResetOrigin(t,offset);
  irtkMatrix mo = offset.GetMatrix();
  irtkMatrix imo = mo;
  imo.Invert();
  
  target.Write("target.nii.gz");
  t.Write("t.nii.gz");
  
  for(i=0;i<_transformations.size();i++)
  {
    irtkMatrix m = _transformations[i].GetMatrix();
    m=imo*m*mo;
    _transformations[i].PutMatrix(m);
  }
  
  //standard deviation for gaussian kernel in number of slices
  double sigma = _motion_sigma;
  int j,iter,par;
  vector<double> kernel;
  for(i=-3*sigma; i<=3*sigma; i++)
  {
    kernel.push_back(exp(-(i/sigma)*(i/sigma)));
    //cout<<kernel[kernel.size()-1]<<" ";
  }
  
  irtkMatrix parameters(6,_transformations.size());
  irtkMatrix weights(6,_transformations.size());
  cout<<"Filling parameters and weights: "<<endl;
  //irtkRealPixel smin, smax;

  ofstream fileOut("motion.txt", ofstream::out | ofstream::app);
  for(i=0;i<_transformations.size();i++)
  {
      parameters(0,i)=_transformations[_slice_order[i]].GetTranslationX();
      parameters(1,i)=_transformations[_slice_order[i]].GetTranslationY();
      parameters(2,i)=_transformations[_slice_order[i]].GetTranslationZ();
      parameters(3,i)=_transformations[_slice_order[i]].GetRotationX();
      parameters(4,i)=_transformations[_slice_order[i]].GetRotationY();
      parameters(5,i)=_transformations[_slice_order[i]].GetRotationZ();
      
      for(j=0;j<6;j++)
      {
	fileOut<<parameters(j,i);
	if(j<5)
	  fileOut<<",";
	else
	  fileOut<<endl;
      }
      
      if (empty)
      {
	
	irtkRealPixel smin, smax;
        _slices[_slice_order[i]].GetMinMax(&smin,&smax);
        if(smax>-1)
	  for(j=0;j<6;j++)
	    weights(j,i)=1;
        else
	  for(j=0;j<6;j++)
	    weights(j,i)=0;
	  
      }
      else
      {	
        if(_slice_inside.size()>0)
        {
          if(_slice_inside[_slice_order[i]])
	    for(j=0;j<6;j++)
	      weights(j,i)=1;
	  else
	    for(j=0;j<6;j++)
	      weights(j,i)=0;
	}
      }
  }
  
  //parameters.Print();
  //weights.Print();
  
  irtkMatrix den(6,_transformations.size());
  irtkMatrix num(6,_transformations.size());
  irtkMatrix kr(6,_transformations.size());
  //irtkMatrix error(6,_transformations.size());
  vector<double> error,tmp;
  double median;
  int packages = 2;
  int dim = _transformations.size()/packages;
  int pack;
  
  cerr<<endl<<endl<<endl<<endl;
  
  for(iter = 0; iter<50;iter++)
  {
    //cerr<<"iter="<<iter<<endl;
    //kernel regression
    for(pack=0;pack<packages;pack++)
    {
      //cerr<<"pack="<<pack<<endl;
      
      for(par=0;par<6;par++)
      {
	//cerr<<"par="<<par<<endl;
        for(i=pack*dim;i<(pack+1)*dim;i++)
	{
	  //cerr<<"i="<<i<<" pack*dim = "<<pack*dim<<" dim = "<<dim<<endl;
	  for(j=-3*sigma;j<=3*sigma;j++)
	    if(((i+j)>=pack*dim) && ((i+j)<(pack+1)*dim))
	    {
	      //cerr<<"j="<<j<<" i+j="<<i+j<<" par="<<par<<" j+3*sigma="<<j+3*sigma<<endl;
	      num(par,i)+=parameters(par,i+j)*kernel[j+3*sigma]*weights(par,i+j);
	      den(par,i)+=kernel[j+3*sigma]*weights(par,i+j);
	    }
	}
      }
    }
    
    for(par=0;par<6;par++)
      for(i=0;i<_transformations.size();i++)
	kr(par,i)=num(par,i)/den(par,i);
    
    //recalculate weights
    for(par=0;par<6;par++)
    {
      error.clear();
      tmp.clear();
     
      for(i=0;i<_transformations.size();i++)
      {
        error.push_back(fabs(parameters(par,i)-kr(par,i)));
        tmp.push_back(fabs(parameters(par,i)-kr(par,i)));
      }
     
      sort(tmp.begin(),tmp.end());
      median = tmp[round(tmp.size()*0.5)];

      for(i=0;i<_transformations.size();i++)
      {
        if(weights(par,i)>0)
        {
          if(error[i]<=median*1.35)
	    weights(par,i)=1;
          else
	    weights(par,i)=median*1.35/error[i];
        }
      }
    }
    
    //test - unify weights
    /*
    for(i=0;i<_transformations.size();i++)
    {
      double w=0;
      for(par=0;par<6;par++)
	//if(weights(par,i)<w)
	  w+=weights(par,i);
      w=w/6;  
      for(par=0;par<6;par++)
	weights(par,i)=w;
    }
   */
  }  
  
  ofstream fileOut2("motion-processed.txt", ofstream::out | ofstream::app);
  for(i=0;i<_transformations.size();i++)
  {
      _transformations[_slice_order[i]].PutTranslationX(kr(0,i));
      _transformations[_slice_order[i]].PutTranslationY(kr(1,i));
      _transformations[_slice_order[i]].PutTranslationZ(kr(2,i));
      _transformations[_slice_order[i]].PutRotationX(kr(3,i));
      _transformations[_slice_order[i]].PutRotationY(kr(4,i));
      _transformations[_slice_order[i]].PutRotationZ(kr(5,i));  

      for(j=0;j<6;j++)
      {
	fileOut2<<kr(j,i);
	if(j<5)
	  fileOut2<<",";
	else
	  fileOut2<<endl;
      }
  }
  weights.Print();
  
  //Put origin back
  //mo.Invert();
  for(i=0;i<_transformations.size();i++)
  {
    irtkMatrix m = _transformations[i].GetMatrix();
    m=mo*m*imo;
    _transformations[i].PutMatrix(m);
  }  
}

void irtkReconstructionDTI::PostProcessMotionParametersHS(irtkRealImage target, bool empty) //, int packages)
{
  if(_slice_order.size()!=_transformations.size())
  {
    cerr<<"Slice order does not match the number of transformations."<<endl;
    exit(1);
  }
  
  //Reset origin for transformations
  int i;
  irtkGreyImage t = target;
  irtkRigidTransformation offset;
  ResetOrigin(t,offset);
  irtkMatrix m;
  irtkMatrix mo = offset.GetMatrix();
  irtkMatrix imo = mo;
  imo.Invert();
  
  target.Write("target.nii.gz");
  t.Write("t.nii.gz");

  
  for(i=0;i<_transformations.size();i++)
  {
    m = _transformations[i].GetMatrix();
    m=imo*m*mo;
    _transformations[i].PutMatrix(m);
  }
  
  //standard deviation for gaussian kernel in number of slices
  double sigma = _motion_sigma;
  int j,iter,par;
  vector<double> kernel;
  for(i=-3*sigma; i<=3*sigma; i++)
  {
    kernel.push_back(exp(-(i/sigma)*(i/sigma)));
    //cout<<kernel[kernel.size()-1]<<" ";
  }
  
  irtkMatrix parameters(6,_transformations.size());
  irtkMatrix weights(6,_transformations.size());
  cout<<"Filling parameters and weights: "<<endl;
  //irtkRealPixel smin, smax;

  ofstream fileOut("motion.txt", ofstream::out | ofstream::app);
  for(i=0;i<_transformations.size();i++)
  {
      parameters(0,i)=_transformations[_slice_order[i]].GetTranslationX();
      parameters(1,i)=_transformations[_slice_order[i]].GetTranslationY();
      parameters(2,i)=_transformations[_slice_order[i]].GetTranslationZ();
      parameters(3,i)=_transformations[_slice_order[i]].GetRotationX();
      parameters(4,i)=_transformations[_slice_order[i]].GetRotationY();
      parameters(5,i)=_transformations[_slice_order[i]].GetRotationZ();
      
      for(j=0;j<6;j++)
      {
	fileOut<<parameters(j,i);
	if(j<5)
	  fileOut<<",";
	else
	  fileOut<<endl;
      }
      
      if (empty)
      {
	
	irtkRealPixel smin, smax;
        _slices[_slice_order[i]].GetMinMax(&smin,&smax);
        if(smax>-1)
	  for(j=0;j<6;j++)
	    weights(j,i)=1;
        else
	  for(j=0;j<6;j++)
	    weights(j,i)=0;
	  
      }
      else
      {	
        if(_slice_inside.size()>0)
        {
          if(_slice_inside[_slice_order[i]])
	    for(j=0;j<6;j++)
	      weights(j,i)=1;
	  else
	    for(j=0;j<6;j++)
	      weights(j,i)=0;
	}
      }
  }
  
  //parameters.Print();
  //weights.Print();
  
  irtkMatrix den(6,_transformations.size());
  irtkMatrix num(6,_transformations.size());
  irtkMatrix kr(6,_transformations.size());
  //irtkMatrix error(6,_transformations.size());
  vector<double> error,tmp;
  double median;
  int packages = 2;
  int dim = _transformations.size()/packages;
  int pack;
  
  cerr<<endl<<endl<<endl<<endl;
  
  for(iter = 0; iter<50;iter++)
  {
    //cerr<<"iter="<<iter<<endl;
    //kernel regression
    for(pack=0;pack<packages;pack++)
    {
      //cerr<<"pack="<<pack<<endl;
      
      for(par=0;par<6;par++)
      {
	//cerr<<"par="<<par<<endl;
        for(i=pack*dim;i<(pack+1)*dim;i++)
	{
	  //cerr<<"i="<<i<<" pack*dim = "<<pack*dim<<" dim = "<<dim<<endl;
	  for(j=-3*sigma;j<=3*sigma;j++)
	    if(((i+j)>=pack*dim) && ((i+j)<(pack+1)*dim))
	    {
	      //cerr<<"j="<<j<<" i+j="<<i+j<<" par="<<par<<" j+3*sigma="<<j+3*sigma<<endl;
	      num(par,i)+=parameters(par,i+j)*kernel[j+3*sigma]*weights(par,i+j);
	      den(par,i)+=kernel[j+3*sigma]*weights(par,i+j);
	    }
	}
      }
    }
    
    for(par=0;par<6;par++)
      for(i=0;i<_transformations.size();i++)
	kr(par,i)=num(par,i)/den(par,i);
    
      cout<<"Iter "<<iter<<": median = ";
    
    //recalculate weights
    //for(par=0;par<6;par++)
    //{
      error.clear();
      tmp.clear();
     
      for(i=0;i<_transformations.size();i++)
      {
	
	///CHange: TRE
	irtkRigidTransformation processed;
        processed.PutTranslationX(kr(0,i));
        processed.PutTranslationY(kr(1,i));
        processed.PutTranslationZ(kr(2,i));
        processed.PutRotationX(kr(3,i));
        processed.PutRotationY(kr(4,i));
        processed.PutRotationZ(kr(5,i));  

	irtkRigidTransformation orig = _transformations[_slice_order[i]];
	
	//need to convert the transformations back to the original coordinate system
	m = orig.GetMatrix();
        m=mo*m*imo;
        orig.PutMatrix(m);

	m = processed.GetMatrix();
        m=mo*m*imo;
        processed.PutMatrix(m);

	irtkRealImage slice = _slices[_slice_order[i]];
	
	double tre=0;
	int n=0;
	double x,y,z,xx,yy,zz,e;
	for(int ii=0;ii<slice.GetX();ii++)
	  for(int jj=0;jj<slice.GetY();jj++)
	    if(slice(ii,jj,0)>-1)
	    {
	      x=ii; y=jj; z=0;
	      slice.ImageToWorld(x,y,z);
	      xx=x;yy=y;zz=z;
	      orig.Transform(x,y,z);
	      processed.Transform(xx,yy,zz);
	      x-=xx;
	      y-=yy;
	      z-=zz;
	      e += sqrt(x*x+y*y+z*z);
	      n++;      
	    }
	if(n>0)
	{
	  e/=n;      
	  error.push_back(e);
	  tmp.push_back(e);
	}
	else
	  error.push_back(-1);
	
	///original
        //error.push_back(fabs(parameters(par,i)-kr(par,i)));
        //tmp.push_back(fabs(parameters(par,i)-kr(par,i)));
      }
     
      sort(tmp.begin(),tmp.end());
      median = tmp[round(tmp.size()*0.5)];
      cout<<median<<" ";

      for(i=0;i<_transformations.size();i++)
      {
        //if(weights(par,i)>0)
	if(error[i]>=0)
        {
          if(error[i]<=median*1.35)
	    for(par=0;par<6;par++)
	      weights(par,i)=1;
          else
	    for(par=0;par<6;par++)
	      weights(par,i)=median*1.35/error[i];
        }
        else
	  for(par=0;par<6;par++)
	      weights(par,i)=0;
      }
    //}
    cout<<endl;
   
  }  
  
  ofstream fileOut2("motion-processed.txt", ofstream::out | ofstream::app);
  ofstream fileOut3("weights.txt", ofstream::out | ofstream::app);
  ofstream fileOut4("outliers.txt", ofstream::out | ofstream::app);
  ofstream fileOut5("empty.txt", ofstream::out | ofstream::app);
  for(i=0;i<_transformations.size();i++)
  {
    //decide whether outlier
    bool outlier  = false;
    for(j=0;j<6;j++)
    {
      if (weights(j,i)<0.5)
	outlier = true;
      fileOut3<<weights(j,i);
      if(j<5)
	fileOut3<<",";
      else
	fileOut3<<endl;
    }
    
    if (outlier)
    {
      if(weights(0,i)>0)
        fileOut4<<_slice_order[i]<<" ";
      else
        fileOut5<<_slice_order[i]<<" ";

      _transformations[_slice_order[i]].PutTranslationX(kr(0,i));
      _transformations[_slice_order[i]].PutTranslationY(kr(1,i));
      _transformations[_slice_order[i]].PutTranslationZ(kr(2,i));
      _transformations[_slice_order[i]].PutRotationX(kr(3,i));
      _transformations[_slice_order[i]].PutRotationY(kr(4,i));
      _transformations[_slice_order[i]].PutRotationZ(kr(5,i));  

      for(j=0;j<6;j++)
      {
	fileOut2<<kr(j,i);
	if(j<5)
	  fileOut2<<",";
	else
	  fileOut2<<endl;
      }
    }
    else
    {
      for(j=0;j<6;j++)
      {
	fileOut2<<parameters(j,i);
	if(j<5)
	  fileOut2<<",";
	else
	  fileOut2<<endl;
      }
    }

  }
  //weights.Print();
  
  //Put origin back
  //mo.Invert();
  for(i=0;i<_transformations.size();i++)
  {
    m = _transformations[i].GetMatrix();
    m=mo*m*imo;
    _transformations[i].PutMatrix(m);
  }
  
}


void irtkReconstructionDTI::PostProcessMotionParametersHS2(irtkRealImage target, bool empty) //, int packages)
{
  if(_slice_order.size()!=_transformations.size())
  {
    cerr<<"Slice order does not match the number of transformations."<<endl;
    exit(1);
  }
  
  //Reset origin for transformations
  int i;
  irtkGreyImage t = target;
  irtkRigidTransformation offset;
  ResetOrigin(t,offset);
  irtkMatrix m;
  irtkMatrix mo = offset.GetMatrix();
  irtkMatrix imo = mo;
  imo.Invert();
  
  target.Write("target.nii.gz");
  t.Write("t.nii.gz");

  
  for(i=0;i<_transformations.size();i++)
  {
    m = _transformations[i].GetMatrix();
    m=imo*m*mo;
    _transformations[i].PutMatrix(m);
  }
  
  //standard deviation for gaussian kernel in number of slices
  double sigma = _motion_sigma;
  int j,iter,par;
  vector<double> kernel;
  for(i=-3*sigma; i<=3*sigma; i++)
  {
    kernel.push_back(exp(-(i/sigma)*(i/sigma)));
    //cout<<kernel[kernel.size()-1]<<" ";
  }
  
  irtkMatrix parameters(6,_transformations.size());
  irtkMatrix weights(6,_transformations.size());
  cout<<"Filling parameters and weights: "<<endl;
  //irtkRealPixel smin, smax;
  
  ofstream fileOut("motion.txt", ofstream::out | ofstream::app);
  for(i=0;i<_transformations.size();i++)
  {
      parameters(0,i)=_transformations[_slice_order[i]].GetTranslationX();
      parameters(1,i)=_transformations[_slice_order[i]].GetTranslationY();
      parameters(2,i)=_transformations[_slice_order[i]].GetTranslationZ();
      parameters(3,i)=_transformations[_slice_order[i]].GetRotationX();
      parameters(4,i)=_transformations[_slice_order[i]].GetRotationY();
      parameters(5,i)=_transformations[_slice_order[i]].GetRotationZ();
      
      for(j=0;j<6;j++)
      {
	fileOut<<parameters(j,i);
	if(j<5)
	  fileOut<<",";
	else
	  fileOut<<endl;
      }
      
      if (empty)
      {
	
	irtkRealPixel smin, smax;
        _slices[_slice_order[i]].GetMinMax(&smin,&smax);
        if(smax>-1)
	  for(j=0;j<6;j++)
	    weights(j,i)=1;
        else
	  for(j=0;j<6;j++)
	    weights(j,i)=0;
	  
      }
      else
      {	
        if(_slice_inside.size()>0)
        {
          if(_slice_inside[_slice_order[i]])
	    for(j=0;j<6;j++)
	      weights(j,i)=1;
	  else
	    for(j=0;j<6;j++)
	      weights(j,i)=0;
	}
      }
  }
  
  //NEW - need to normalize rotation so that the first one is zero
  // this is to deal with cases when rotation is around 180
  
  double offsetrot[6];
  offsetrot[0] = 0;
  offsetrot[1] = 0;
  offsetrot[2] = 0;
  offsetrot[3] = parameters(3,0);
  offsetrot[4] = parameters(4,0);
  offsetrot[5] = parameters(5,0);
  
  /*
  for (i=0; i<_transformations.size(); i++)
  {
    parameters(3,i)-=rx0;
    parameters(4,i)-=ry0;
    parameters(5,i)-=rz0;
  }
  */
  //parameters.Print();
  //weights.Print();
  
  irtkMatrix den(6,_transformations.size());
  irtkMatrix num(6,_transformations.size());
  irtkMatrix kr(6,_transformations.size());
  //irtkMatrix error(6,_transformations.size());
  vector<double> error,tmp;
  double median;
  int packages = 2;
  int dim = _transformations.size()/packages;
  int pack;
  
  cerr<<endl<<endl<<endl<<endl;
  
  for(iter = 0; iter<50;iter++)
  {
    //cerr<<"iter="<<iter<<endl;
    //kernel regression
    for(pack=0;pack<packages;pack++)
    {
      //cerr<<"pack="<<pack<<endl;
      
      for(par=0;par<6;par++)
      {
	//cerr<<"par="<<par<<endl;
        for(i=pack*dim;i<(pack+1)*dim;i++)
	{
	  //cerr<<"i="<<i<<" pack*dim = "<<pack*dim<<" dim = "<<dim<<endl;
	  for(j=-3*sigma;j<=3*sigma;j++)
	    if(((i+j)>=pack*dim) && ((i+j)<(pack+1)*dim))
	    {
	      //cerr<<"j="<<j<<" i+j="<<i+j<<" par="<<par<<" j+3*sigma="<<j+3*sigma<<endl;
	      double value = parameters(par,i+j)-offsetrot[par];
	      if(value>180)
		value-=360;
	      if(value<-180)
		value+=360;
	      num(par,i)+=(value)*kernel[j+3*sigma]*weights(par,i+j);
	      den(par,i)+=kernel[j+3*sigma]*weights(par,i+j);
	    }
	}
      }
    }
    
    for(par=0;par<6;par++)
      for(i=0;i<_transformations.size();i++)
      {
	double value=num(par,i)/den(par,i)+offsetrot[par];
	if(value<-180)
	  value+=360;
	if(value>180)
	  value-=360;
	kr(par,i)=value;
      }
    
      cout<<"Iter "<<iter<<": median = ";
      cout.flush();
    
    //recalculate weights
    //for(par=0;par<6;par++)
    //{
      error.clear();
      tmp.clear();
     
      for(i=0;i<_transformations.size();i++)
      {
	//cout<<i<<" ";
	//cout.flush();
	///CHange: TRE
	irtkRigidTransformation processed;
        processed.PutTranslationX(kr(0,i));
        processed.PutTranslationY(kr(1,i));
        processed.PutTranslationZ(kr(2,i));
        processed.PutRotationX(kr(3,i));
        processed.PutRotationY(kr(4,i));
        processed.PutRotationZ(kr(5,i));  

	irtkRigidTransformation orig = _transformations[_slice_order[i]];
	
	//need to convert the transformations back to the original coordinate system
	m = orig.GetMatrix();
        m=mo*m*imo;
        orig.PutMatrix(m);

	m = processed.GetMatrix();
        m=mo*m*imo;
        processed.PutMatrix(m);

	irtkRealImage slice = _slices[_slice_order[i]];
	
	double tre=0;
	int n=0;
	double x,y,z,xx,yy,zz,e;
	for(int ii=0;ii<_reconstructed.GetX();ii=ii+3)
	  for(int jj=0;jj<_reconstructed.GetY();jj=jj+3)
	    for(int kk=0;kk<_reconstructed.GetZ();kk=kk+3)
	      if(_reconstructed(ii,jj,kk)>-1)
	      {
		//cout<<ii<<" "<<jj<<" "<<kk<<endl;
		//cout.flush();
	        x=ii; y=jj; z=kk;
	        _reconstructed.ImageToWorld(x,y,z);
	        xx=x;yy=y;zz=z;
	        orig.Transform(x,y,z);
	        processed.Transform(xx,yy,zz);
	        x-=xx;
	        y-=yy;
	        z-=zz;
	        e += sqrt(x*x+y*y+z*z);
	        n++;      
	      }
	if(n>0)
	{
	  e/=n;      
	  error.push_back(e);
	  tmp.push_back(e);
	}
	else
	  error.push_back(-1);
	
	///original
        //error.push_back(fabs(parameters(par,i)-kr(par,i)));
        //tmp.push_back(fabs(parameters(par,i)-kr(par,i)));
      }
     
      sort(tmp.begin(),tmp.end());
      median = tmp[round(tmp.size()*0.5)];
      cout<<median<<" ";

      for(i=0;i<_transformations.size();i++)
      {
        //if(weights(par,i)>0)
	if(error[i]>=0)
        {
          if(error[i]<=median*1.35)
	    for(par=0;par<6;par++)
	      weights(par,i)=1;
          else
	    for(par=0;par<6;par++)
	      weights(par,i)=median*1.35/error[i];
        }
        else
	  for(par=0;par<6;par++)
	      weights(par,i)=0;
      }
    //}
    cout<<endl;
   
  }  
  
  ofstream fileOut2("motion-processed.txt", ofstream::out | ofstream::app);
  ofstream fileOut3("weights.txt", ofstream::out | ofstream::app);
  ofstream fileOut4("outliers.txt", ofstream::out | ofstream::app);
  ofstream fileOut5("empty.txt", ofstream::out | ofstream::app);
  for(i=0;i<_transformations.size();i++)
  {
    //decide whether outlier
    bool outlier  = false;
    for(j=0;j<6;j++)
    {
      if (weights(j,i)<0.5)
	outlier = true;
      fileOut3<<weights(j,i);
      if(j<5)
	fileOut3<<",";
      else
	fileOut3<<endl;
    }
    
    if (outlier)
    {
      if(weights(0,i)>0)
        fileOut4<<_slice_order[i]<<" ";
      else
        fileOut5<<_slice_order[i]<<" ";

      _transformations[_slice_order[i]].PutTranslationX(kr(0,i));
      _transformations[_slice_order[i]].PutTranslationY(kr(1,i));
      _transformations[_slice_order[i]].PutTranslationZ(kr(2,i));
      _transformations[_slice_order[i]].PutRotationX(kr(3,i));
      _transformations[_slice_order[i]].PutRotationY(kr(4,i));
      _transformations[_slice_order[i]].PutRotationZ(kr(5,i));  

      for(j=0;j<6;j++)
      {
	fileOut2<<kr(j,i);
	if(j<5)
	  fileOut2<<",";
	else
	  fileOut2<<endl;
      }
    }
    else
    {
      for(j=0;j<6;j++)
      {
	fileOut2<<parameters(j,i);
	if(j<5)
	  fileOut2<<",";
	else
	  fileOut2<<endl;
      }
    }

  }
  //weights.Print();
  
  //Put origin back
  //mo.Invert();
  for(i=0;i<_transformations.size();i++)
  {
    m = _transformations[i].GetMatrix();
    m=mo*m*imo;
    _transformations[i].PutMatrix(m);
  }
  
}


void irtkReconstructionDTI::PostProcessMotionParameters2(irtkRealImage target, bool empty) //, int packages)
{
  if(_slice_order.size()!=_transformations.size())
  {
    cerr<<"Slice order does not match the number of transformations."<<endl;
    exit(1);
  }
  
  //Reset origin for transformations
  int i;
  irtkGreyImage t = target;
  irtkRigidTransformation offset;
  ResetOrigin(t,offset);
  irtkMatrix m;
  irtkMatrix mo = offset.GetMatrix();
  irtkMatrix imo = mo;
  imo.Invert();
  
  target.Write("target.nii.gz");
  t.Write("t.nii.gz");

  
  for(i=0;i<_transformations.size();i++)
  {
    m = _transformations[i].GetMatrix();
    m=imo*m*mo;
    _transformations[i].PutMatrix(m);
  }
  
  //standard deviation for gaussian kernel in number of slices
  double sigma = _motion_sigma;
  int j,iter,par;
  vector<double> kernel;
  for(i=-3*sigma; i<=3*sigma; i++)
  {
    kernel.push_back(exp(-(i/sigma)*(i/sigma)));
    //cout<<kernel[kernel.size()-1]<<" ";
  }
  
  irtkMatrix parameters(6,_transformations.size());
  irtkMatrix weights(6,_transformations.size());
  cout<<"Filling parameters and weights: "<<endl;
  //irtkRealPixel smin, smax;

  ofstream fileOut("motion.txt", ofstream::out | ofstream::app);
  for(i=0;i<_transformations.size();i++)
  {
      parameters(0,i)=_transformations[_slice_order[i]].GetTranslationX();
      parameters(1,i)=_transformations[_slice_order[i]].GetTranslationY();
      parameters(2,i)=_transformations[_slice_order[i]].GetTranslationZ();
      parameters(3,i)=_transformations[_slice_order[i]].GetRotationX();
      parameters(4,i)=_transformations[_slice_order[i]].GetRotationY();
      parameters(5,i)=_transformations[_slice_order[i]].GetRotationZ();
      
      for(j=0;j<6;j++)
      {
	fileOut<<parameters(j,i);
	if(j<5)
	  fileOut<<",";
	else
	  fileOut<<endl;
      }
      
      if (empty)
      {
	
	irtkRealPixel smin, smax;
        _slices[_slice_order[i]].GetMinMax(&smin,&smax);
        if(smax>-1)
	  for(j=0;j<6;j++)
	    weights(j,i)=1;
        else
	  for(j=0;j<6;j++)
	    weights(j,i)=0;
	  
      }
      else
      {	
        if(_slice_inside.size()>0)
        {
          if(_slice_inside[_slice_order[i]])
	    for(j=0;j<6;j++)
	      weights(j,i)=1;
	  else
	    for(j=0;j<6;j++)
	      weights(j,i)=0;
	}
      }
  }
  
  //parameters.Print();
  //weights.Print();
  
  irtkMatrix den(6,_transformations.size());
  irtkMatrix num(6,_transformations.size());
  irtkMatrix kr(6,_transformations.size());
  //irtkMatrix error(6,_transformations.size());
  vector<double> error,tmp;
  double median;
  int packages = 2;
  int dim = _transformations.size()/packages;
  int pack;
  
  cerr<<endl<<endl<<endl<<endl;
  
  for(iter = 0; iter<50;iter++)
  {
    //cerr<<"iter="<<iter<<endl;
    //kernel regression
    for(pack=0;pack<packages;pack++)
    {
      //cerr<<"pack="<<pack<<endl;
      
      for(par=0;par<6;par++)
      {
	//cerr<<"par="<<par<<endl;
        for(i=pack*dim;i<(pack+1)*dim;i++)
	{
	  //cerr<<"i="<<i<<" pack*dim = "<<pack*dim<<" dim = "<<dim<<endl;
	  for(j=-3*sigma;j<=3*sigma;j++)
	    if(((i+j)>=pack*dim) && ((i+j)<(pack+1)*dim))
	    {
	      //cerr<<"j="<<j<<" i+j="<<i+j<<" par="<<par<<" j+3*sigma="<<j+3*sigma<<endl;
	      num(par,i)+=parameters(par,i+j)*kernel[j+3*sigma]*weights(par,i+j);
	      den(par,i)+=kernel[j+3*sigma]*weights(par,i+j);
	    }
	}
      }
    }
    
    for(par=0;par<6;par++)
      for(i=0;i<_transformations.size();i++)
	kr(par,i)=num(par,i)/den(par,i);
    
      cout<<"Iter "<<iter<<": median = ";
      cout.flush();
    
    //recalculate weights
    //for(par=0;par<6;par++)
    //{
      error.clear();
      tmp.clear();
     
      for(i=0;i<_transformations.size();i++)
      {
	//cout<<i<<" ";
	//cout.flush();
	///CHange: TRE
	irtkRigidTransformation processed;
        processed.PutTranslationX(kr(0,i));
        processed.PutTranslationY(kr(1,i));
        processed.PutTranslationZ(kr(2,i));
        processed.PutRotationX(kr(3,i));
        processed.PutRotationY(kr(4,i));
        processed.PutRotationZ(kr(5,i));  

	irtkRigidTransformation orig = _transformations[_slice_order[i]];
	
	//need to convert the transformations back to the original coordinate system
	m = orig.GetMatrix();
        m=mo*m*imo;
        orig.PutMatrix(m);

	m = processed.GetMatrix();
        m=mo*m*imo;
        processed.PutMatrix(m);

	irtkRealImage slice = _slices[_slice_order[i]];
	
	double tre=0;
	int n=0;
	double x,y,z,xx,yy,zz,e;
	for(int ii=0;ii<_reconstructed.GetX();ii=ii+3)
	  for(int jj=0;jj<_reconstructed.GetY();jj=jj+3)
	    for(int kk=0;kk<_reconstructed.GetZ();kk=kk+3)
	      if(_reconstructed(ii,jj,kk)>-1)
	      {
		//cout<<ii<<" "<<jj<<" "<<kk<<endl;
		//cout.flush();
	        x=ii; y=jj; z=kk;
	        _reconstructed.ImageToWorld(x,y,z);
	        xx=x;yy=y;zz=z;
	        orig.Transform(x,y,z);
	        processed.Transform(xx,yy,zz);
	        x-=xx;
	        y-=yy;
	        z-=zz;
	        e += sqrt(x*x+y*y+z*z);
	        n++;      
	      }
	if(n>0)
	{
	  e/=n;      
	  error.push_back(e);
	  tmp.push_back(e);
	}
	else
	  error.push_back(-1);
	
	///original
        //error.push_back(fabs(parameters(par,i)-kr(par,i)));
        //tmp.push_back(fabs(parameters(par,i)-kr(par,i)));
      }
     
      sort(tmp.begin(),tmp.end());
      median = tmp[round(tmp.size()*0.5)];
      cout<<median<<" ";

      for(i=0;i<_transformations.size();i++)
      {
        //if(weights(par,i)>0)
	if(error[i]>=0)
        {
          if(error[i]<=median*1.35)
	    for(par=0;par<6;par++)
	      weights(par,i)=1;
          else
	    for(par=0;par<6;par++)
	      weights(par,i)=median*1.35/error[i];
        }
        else
	  for(par=0;par<6;par++)
	      weights(par,i)=0;
      }
    //}
    cout<<endl;
   
  }  
  
  ofstream fileOut2("motion-processed.txt", ofstream::out | ofstream::app);
  ofstream fileOut3("weights.txt", ofstream::out | ofstream::app);
  ofstream fileOut4("outliers.txt", ofstream::out | ofstream::app);
  ofstream fileOut5("empty.txt", ofstream::out | ofstream::app);
  
  for(i=0;i<_transformations.size();i++)
  {
      fileOut3<<weights(j,i);
    
      if(weights(0,i)<=0)
        fileOut5<<_slice_order[i]<<" ";

      _transformations[_slice_order[i]].PutTranslationX(kr(0,i));
      _transformations[_slice_order[i]].PutTranslationY(kr(1,i));
      _transformations[_slice_order[i]].PutTranslationZ(kr(2,i));
      _transformations[_slice_order[i]].PutRotationX(kr(3,i));
      _transformations[_slice_order[i]].PutRotationY(kr(4,i));
      _transformations[_slice_order[i]].PutRotationZ(kr(5,i));  

      for(j=0;j<6;j++)
      {
	fileOut2<<kr(j,i);
	if(j<5)
	  fileOut2<<",";
	else
	  fileOut2<<endl;
      }

  }
  
  //Put origin back
  for(i=0;i<_transformations.size();i++)
  {
    m = _transformations[i].GetMatrix();
    m=mo*m*imo;
    _transformations[i].PutMatrix(m);
  }
  
}





class ParallelSliceToVolumeRegistrationSH {
public:
    irtkReconstructionDTI *reconstructor;

    ParallelSliceToVolumeRegistrationSH(irtkReconstructionDTI *_reconstructor) : 
    reconstructor(_reconstructor) { }

    void operator() (const blocked_range<size_t> &r) const {

        irtkImageAttributes attr = reconstructor->_reconstructed.GetImageAttributes();
        
        for ( size_t inputIndex = r.begin(); inputIndex != r.end(); ++inputIndex ) {

            irtkImageRigidRegistrationWithPadding registration;
            irtkGreyPixel smin, smax;
            irtkGreyImage target;
            irtkRealImage slice, w, b, t;
            irtkResamplingWithPadding<irtkRealPixel> resampling(attr._dx,attr._dx,attr._dx,-1);         
            irtkReconstruction dummy_reconstruction;
	    char buffer[256];
        
            //target = _slices[inputIndex];
            t = reconstructor->_slices[inputIndex];
            resampling.SetInput(&reconstructor->_slices[inputIndex]);
            resampling.SetOutput(&t);
            resampling.Run();
            target=t;
	    int ii = inputIndex;
	    sprintf(buffer,"target%i.nii.gz",ii);
	    target.Write(buffer);

            target.GetMinMax(&smin, &smax);
        
            if (smax > -1) {
                //put origin to zero
                irtkRigidTransformation offset;
                dummy_reconstruction.ResetOrigin(target,offset);
                irtkMatrix mo = offset.GetMatrix();
                irtkMatrix m = reconstructor->_transformations[inputIndex].GetMatrix();
                m=m*mo;
                reconstructor->_transformations[inputIndex].PutMatrix(m);

                irtkGreyImage source = reconstructor->_simulated_signal.GetRegion(0,0,0,reconstructor->_stack_index[inputIndex],reconstructor->_simulated_signal.GetX(),reconstructor->_simulated_signal.GetY(),reconstructor->_simulated_signal.GetZ(),reconstructor->_stack_index[inputIndex]+1);
		sprintf(buffer,"source%i.nii.gz",ii);
		source.Write(buffer);
                registration.SetInput(&target, &source);
                registration.SetOutput(&reconstructor->_transformations[inputIndex]);
                registration.GuessParameterSliceToVolume();
                registration.SetTargetPadding(-1);
                registration.Run();
                //undo the offset
                mo.Invert();
                m = reconstructor->_transformations[inputIndex].GetMatrix();
                m=m*mo;
                reconstructor->_transformations[inputIndex].PutMatrix(m);
            }      
        }
    }

    // execute
    void operator() () const {
        task_scheduler_init init(tbb_no_threads);
        parallel_for( blocked_range<size_t>(0, reconstructor->_slices.size() ),
                      *this );
        init.terminate();
    }

};

void irtkReconstructionDTI::SliceToVolumeRegistrationSH()
{
    if (_debug)
        cout << "SliceToVolumeRegistration" << endl;
    ParallelSliceToVolumeRegistrationSH registration(this);
    registration();
}


class ParallelNormaliseBiasDTI{
    irtkReconstructionDTI* reconstructor;
public:
    irtkRealImage bias, weights;
    
    void operator()( const blocked_range<size_t>& r ) {
        for ( size_t inputIndex = r.begin(); inputIndex < r.end(); ++inputIndex) {

            if(reconstructor->_debug) {
                cout<<inputIndex<<" ";
            }
            
            // alias the current slice
            irtkRealImage& slice = reconstructor->_slices[inputIndex];
                
            //read the current bias image
            irtkRealImage b = reconstructor->_bias[inputIndex];
                
            //read current scale factor
            double scale = reconstructor->_scale[inputIndex];

            irtkRealPixel *pi = slice.GetPointerToVoxels();
            irtkRealPixel *pb = b.GetPointerToVoxels();
            for(int i = 0; i<slice.GetNumberOfVoxels(); i++) {
                if((*pi>-1)&&(scale>0))
                    *pb -= log(scale);
                pb++;
                pi++;
            }
                
            //Distribute slice intensities to the volume
            POINT3D p;
	    int stackIndex = reconstructor->_stack_index[inputIndex];
	    double sliceWeight = reconstructor->_slice_weight[inputIndex];
            for (int i = 0; i < slice.GetX(); i++)
                for (int j = 0; j < slice.GetY(); j++)
                    if (slice(i, j, 0) != -1) {
                        //number of volume voxels with non-zero coefficients for current slice voxel
                        int n = reconstructor->_volcoeffs[inputIndex][i][j].size();
                        //add contribution of current slice voxel to all voxel volumes
                        //to which it contributes
                        for (int k = 0; k < n; k++) {
                            p = reconstructor->_volcoeffs[inputIndex][i][j][k];
                            bias(p.x, p.y, p.z,stackIndex) += sliceWeight * p.value * b(i, j, 0);
                            weights(p.x, p.y, p.z,stackIndex) += sliceWeight * p.value;
                        }
                    }
            //end of loop for a slice inputIndex                
        }
    }
 
    ParallelNormaliseBiasDTI( ParallelNormaliseBiasDTI& x, split ) :
        reconstructor(x.reconstructor)
    {
        bias.Initialize( reconstructor->_simulated_signal.GetImageAttributes() );
        bias = 0;    
	weights.Initialize( reconstructor->_simulated_signal.GetImageAttributes() );
	weights = 0;
    }
 
    void join( const ParallelNormaliseBiasDTI& y ) {       
        bias += y.bias;
	weights += y.weights;
    }
             
    ParallelNormaliseBiasDTI( irtkReconstructionDTI *reconstructor ) :
    reconstructor(reconstructor)
    {
        bias.Initialize( reconstructor->_simulated_signal.GetImageAttributes() );
        bias = 0;    
	weights.Initialize( reconstructor->_simulated_signal.GetImageAttributes() );
	weights = 0;
    }

    // execute
    void operator() () {
        task_scheduler_init init(tbb_no_threads);
        parallel_reduce( blocked_range<size_t>(0,reconstructor->_slices.size()),
                         *this );
        init.terminate();
    }        
};
/*
class ParallelNormaliseBiasDTI2{
    irtkReconstructionDTI* reconstructor;
    irtkRealImage& bias;
public:
    ParallelNormaliseBiasDTI2( irtkReconstructionDTI *_reconstructor, irtkRealImage &_bias ) :
    reconstructor(_reconstructor), bias(_bias)
    {}
    
    void operator()( const blocked_range<size_t>& r ) {
        for ( size_t inputIndex = r.begin(); inputIndex < r.end(); ++inputIndex) {

            if(reconstructor->_debug) {
                cout<<inputIndex<<" ";
            }
            
            // alias the current slice
            irtkRealImage& slice = reconstructor->_slices[inputIndex];
                
            //read the current bias image
            irtkRealImage b = reconstructor->_bias[inputIndex];
	    b=0;
                
            //read current scale factor
            //double scale = reconstructor->_scale[inputIndex];

            //irtkRealPixel *pi = slice.GetPointerToVoxels();
            //irtkRealPixel *pb = b.GetPointerToVoxels();
            //for(int i = 0; i<slice.GetNumberOfVoxels(); i++) {
            //    if((*pi>-1)&&(scale>0))
            //        *pb -= log(scale);
            //    pb++;
            //    pi++;
            //}
                
            //Distribute slice intensities to the volume
            POINT3D p;
	    int stackIndex = reconstructor->_stack_index[inputIndex];
	    double sliceWeight = reconstructor->_slice_weight[inputIndex];
            for (int i = 0; i < slice.GetX(); i++)
                for (int j = 0; j < slice.GetY(); j++)
                    if (slice(i, j, 0) != -1) {
		        double weight = 0;
                        //number of volume voxels with non-zero coefficients for current slice voxel
                        int n = reconstructor->_volcoeffs[inputIndex][i][j].size();
                        //add contribution of current slice voxel to all voxel volumes
                        //to which it contributes
                        for (int k = 0; k < n; k++) {
                            p = reconstructor->_volcoeffs[inputIndex][i][j][k];
			    b(i,j,0) += p.value * bias(p.x, p.y, p.z,stackIndex);
			    weight += p.value;
                        }
                        if( weight > 0 ) {
			  b(i,j,0)/=weight;
			}
		    reconstructor->_bias[inputIndex] -= b;
                    }
            //end of loop for a slice inputIndex                
        }
    }
 


    // execute
    void operator() () {
        task_scheduler_init init(tbb_no_threads);
        parallel_reduce( blocked_range<size_t>(0,reconstructor->_slices.size()),
                         *this );
        init.terminate();
    }        
};
*/

class ParallelNormaliseBiasDTI2 {
    irtkReconstruction *reconstructor;
    irtkRealImage& bias;    
public:
    ParallelNormaliseBiasDTI2( irtkReconstructionDTI *_reconstructor, irtkRealImage &_bias ) :
    reconstructor(_reconstructor), bias(_bias)
    {}
    
    void operator() (const blocked_range<size_t> &r) const {
      for ( size_t inputIndex = r.begin(); inputIndex < r.end(); ++inputIndex) {
            if(reconstructor->_debug) {
                cout<<inputIndex<<" ";
            }
            int ii = inputIndex;
	     char buffer[256];
	     
            // alias the current slice
            irtkRealImage& slice = reconstructor->_slices[inputIndex];
                
            //read the current bias image
            irtkRealImage b(slice.GetImageAttributes());
    	    //sprintf(buffer,"biaszero%i.nii.gz",ii);
            //b.Write(buffer);

	    POINT3D p;
	    int stackIndex = reconstructor->_stack_index[inputIndex];
	    double sliceWeight = reconstructor->_slice_weight[inputIndex];
            for (int i = 0; i < slice.GetX(); i++)
                for (int j = 0; j < slice.GetY(); j++)
                    if (slice(i, j, 0) != -1) {
		        double weight = 0;
                        //number of volume voxels with non-zero coefficients for current slice voxel
                        int n = reconstructor->_volcoeffs[inputIndex][i][j].size();
                        //add contribution of current slice voxel to all voxel volumes
                        //to which it contributes
                        for (int k = 0; k < n; k++) {
                            p = reconstructor->_volcoeffs[inputIndex][i][j][k];
			    b(i,j,0) += p.value * bias(p.x, p.y, p.z,stackIndex);
			    weight += p.value;
                        }
                        if( weight > 0 ) {
			  b(i,j,0)/=weight;
			}
   

                    }
   	  //sprintf(buffer,"biasold%i.nii.gz",ii);
	  //reconstructor->_bias[inputIndex].Write(buffer);
	  reconstructor->_bias[inputIndex] -= b;
		        
   
       
	
        //sprintf(buffer,"b%i.nii.gz",ii);
        //b.Write(buffer);
        //sprintf(buffer,"biasnew%i.nii.gz",ii);
        //reconstructor->_bias[inputIndex].Write(buffer);

            //end of loop for a slice inputIndex  
      }
    }

    // execute
    void operator() () const {
        task_scheduler_init init(tbb_no_threads);
        parallel_for( blocked_range<size_t>(0, reconstructor->_slices.size()),
                      *this );
        init.terminate();
    }

};

void irtkReconstructionDTI::NormaliseBiasDTI(int iter, vector<irtkRigidTransformation> &stack_transformations, int order)
{
    if(_debug)
        cout << "Normalise Bias ... ";

    ParallelNormaliseBiasDTI parallelNormaliseBiasDTI(this);
    parallelNormaliseBiasDTI();
    irtkRealImage bias = parallelNormaliseBiasDTI.bias;
    irtkRealImage weights = parallelNormaliseBiasDTI.weights;
    
    // normalize the volume by proportion of contributing slice voxels for each volume voxel
    //char buffer[256];
    //sprintf(buffer,"ab%i.nii.gz",iter);
    //bias.Write(buffer);
    //sprintf(buffer,"abw%i.nii.gz",iter);
    //weights.Write(buffer);

    bias /= weights;
    if (_debug) {
        char buffer[256];
        sprintf(buffer,"averagebiasorig%i.nii.gz",iter);
        bias.Write(buffer);
    }
        
    irtkRealImage m = bias;
    m=0;
    irtkRealPixel *pm, *pb;
    pm = m.GetPointerToVoxels();
    pb = bias.GetPointerToVoxels();
    for (int i = 0; i<m.GetNumberOfVoxels();i++) {
        if(*pb!=0)
            *pm =1;
        pm++;
        pb++;
    }
    
    if (_debug) {
        char buffer[256];
        sprintf(buffer,"averagebiasmask%i.nii.gz",iter);
        m.Write(buffer);
    }
    
 
    
    irtkGaussianBlurring<irtkRealPixel> gb(_sigma_bias);
    gb.SetInput(&bias);
    gb.SetOutput(&bias);
    gb.Run();
    gb.SetInput(&m);
    gb.SetOutput(&m);
    gb.Run();
    bias/=m;

    if (_debug) {
        char buffer[256];
        sprintf(buffer,"averagebiasSH%i.nii.gz",iter);
        bias.Write(buffer);
    }

    ParallelNormaliseBiasDTI2 parallelNormaliseBiasDTI2(this,bias);
    parallelNormaliseBiasDTI2();
    

    if(_debug)
        cout << "done." << endl;
 /*
       // rotate directions
    int dirIndex;
    double gx,gy,gz;
    int nStacks = bias.GetT();
    irtkMatrix dirs_xyz(nStacks,3);   
    for (dirIndex = 0; dirIndex < nStacks; dirIndex++)  
    {
      gx=_directions[0][dirIndex+1];
      gy=_directions[1][dirIndex+1];
      gz=_directions[2][dirIndex+1];
      RotateDirections(gx,gy,gz,stack_transformations[dirIndex]);
      dirs_xyz(dirIndex,0) = gx;
      dirs_xyz(dirIndex,1) = gy;
      dirs_xyz(dirIndex,2) = gz;     
    }
    
    irtkSphericalHarmonics sh;
    irtkMatrix SHT = sh.SHbasis(dirs_xyz,order);
    irtkMatrix tSHT = SHT;
    tSHT.Transpose();
 
    irtkRealImage biasSH(_SH_coeffs);
    biasSH=0;
    irtkVector s(nStacks),c,ev;
    irtkMatrix w(nStacks,nStacks);
    irtkMatrix tmp;
    irtkMatrix tmp1,tmp2;
    irtkRealImage mask = _reconstructed;
    mask=0;
    
    double det;
    int i,j,k,t;
    for(i=0;i<biasSH.GetX();i++)
      for(j=0;j<biasSH.GetY();j++)
        for(k=0;k<biasSH.GetZ();k++)
	{
	  for(t=0;t<nStacks;t++)
	  {
	    s(t)=bias(i,j,k,t);
 	    w(t,t)=weights(i,j,k,t);
	  }

	  tmp = tSHT*w*SHT;
	  tmp.Eigenvalues(tmp1,ev,tmp2);
	  
	  det=1;
	  for(t=0;t<ev.Rows();t++)
	    det*=ev(t);
	  //cerr<<det<<" ";
	  if(det>0.001)
	  {
	    mask(i,j,k)=1;
	    tmp.Invert();
	    c = tmp*tSHT*w*s;
	    for(t=0;t<SHT.Cols();t++)
	    {
	      biasSH(i,j,k,t)=c(t);
	    }
	  }
	}
        if (_debug) {
        char buffer[256];
        sprintf(buffer,"averagebiasSHSHorig%i.nii.gz",iter);
        biasSH.Write(buffer);
    }
    
    
    irtkRealImage m = biasSH;
    m=0;
    irtkRealPixel *pm, *pb;
    pm = m.GetPointerToVoxels();
    pb = biasSH.GetPointerToVoxels();
    for (int i = 0; i<m.GetNumberOfVoxels();i++) {
        if(*pb!=0)
            *pm =1;
        pm++;
        pb++;
    }
    
    if (_debug) {
        char buffer[256];
        sprintf(buffer,"averagebiasSHmask%i.nii.gz",iter);
        m.Write(buffer);
    }
    
 
    
    irtkGaussianBlurring<irtkRealPixel> gb(_sigma_bias);
    gb.SetInput(&biasSH);
    gb.SetOutput(&biasSH);
    gb.Run();
    gb.SetInput(&m);
    gb.SetOutput(&m);
    gb.Run();
    biasSH/=m;

    if (_debug) {
        char buffer[256];
        sprintf(buffer,"averagebiasSH%i.nii.gz",iter);
        biasSH.Write(buffer);
        sprintf(buffer,"origSH%i.nii.gz",iter);
        _SH_coeffs.Write(buffer);
    }

    
   
    irtkRealPixel *pi, *pt;
    pi = _SH_coeffs.GetPointerToVoxels();
    pb = biasSH.GetPointerToVoxels();
    irtkRealImage b = biasSH;
    b=0;
    pt = b.GetPointerToVoxels();
    for (int i = 0; i<_SH_coeffs.GetNumberOfVoxels();i++) {
        if(*pi!=0)
	{
            *pi /=exp(-(*pb));
	    *pt = exp(-(*pb));
	}
        pi++;
        pb++;
	pt++;
    }
    char buffer[256];
    sprintf(buffer,"b%i.nii.gz",iter);
    b.Write(buffer);
    */
}


void irtkReconstructionDTI::NormaliseBiasSH2(int iter, vector<irtkRigidTransformation> &stack_transformations, int order)
{
    if(_debug)
        cout << "Normalise Bias ... ";

    // Calculate stack bias
    ParallelNormaliseBiasDTI parallelNormaliseBiasDTI(this);
    parallelNormaliseBiasDTI();
    irtkRealImage bias = parallelNormaliseBiasDTI.bias;
    irtkRealImage weights = parallelNormaliseBiasDTI.weights;
    
    // normalize the volume by proportion of contributing slice voxels for each volume voxel
    //char buffer[256];
    //sprintf(buffer,"ab%i.nii.gz",iter);
    //bias.Write(buffer);
    //sprintf(buffer,"abw%i.nii.gz",iter);
    //weights.Write(buffer);

    bias /= weights;
    if (_debug) {
        char buffer[256];
        sprintf(buffer,"averagebiasorig%i.nii.gz",iter);
        bias.Write(buffer);
    }
    
    ////////////////////////////////////////

    //transfer stack bias to SH space
    // rotate directions
    int dirIndex;
    double gx,gy,gz;
    int nStacks = bias.GetT();
    irtkMatrix dirs_xyz(nStacks,3);   
    for (dirIndex = 0; dirIndex < nStacks; dirIndex++)  
    {
      gx=_directions[0][dirIndex+1];
      gy=_directions[1][dirIndex+1];
      gz=_directions[2][dirIndex+1];
      RotateDirections(gx,gy,gz,stack_transformations[dirIndex]);
      dirs_xyz(dirIndex,0) = gx;
      dirs_xyz(dirIndex,1) = gy;
      dirs_xyz(dirIndex,2) = gz;     
    }
    
    irtkSphericalHarmonics sh;
    irtkMatrix SHT = sh.SHbasis(dirs_xyz,order);
    irtkMatrix tSHT = SHT;
    tSHT.Transpose();
 
    irtkRealImage biasSH(_SH_coeffs);
    biasSH=0;
    irtkRealImage mask(_SH_coeffs);
    mask=0;
    
    int nCoeff = _SH_coeffs.GetT();
    irtkVector c(nCoeff),s,sb,cb,ev;
    irtkMatrix w(nStacks,nStacks),tmp,tmp1,tmp2;
    double b,det;
    int i,j,k,t;
    for(i=0;i<biasSH.GetX();i++)
      for(j=0;j<biasSH.GetY();j++)
        for(k=0;k<biasSH.GetZ();k++)
	  if(_SH_coeffs(i,j,k,0)>0)
	  {
	    //SH coeff
	    for(t=0; t<nCoeff;t++)
	      c(t)=_SH_coeffs(i,j,k,t);
	    //simulate HR signal
	    s = SHT*c;
	    //add bias to it
	    for(t=0; t<nStacks;t++)
	    {
	      b=bias(i,j,k,t);
	      s(t)/=exp(-b);
	      w(t,t)=weights(i,j,k,t);
	    }
	    // calculate SH coefficients with the effect of the bias
	    //and bias itself
	    tmp = tSHT*w*SHT;
	    tmp.Eigenvalues(tmp1,ev,tmp2);
	  
	    det=1;
	    for(t=0;t<ev.Rows();t++)
	      det*=ev(t);
	    if(det>0.001)
	    {
	      tmp.Invert();
	      cb = tmp*tSHT*w*s;
	      for(t=0;t<SHT.Cols();t++)
	        if((cb(t)!=0)&&(c(t)!=0))
	        {
	          biasSH(i,j,k,t)=log(cb(t)/c(t));
		  mask(i,j,k,t)=1;
	        }
	    }

	  }
    
    
    //////////////////////////////////
     
    
    if (_debug) {
        char buffer[256];
        sprintf(buffer,"averagebiasmask%i.nii.gz",iter);
        mask.Write(buffer);
    }
    
    
 
    // blur
    irtkGaussianBlurring<irtkRealPixel> gb(_sigma_bias);
    gb.SetInput(&biasSH);
    gb.SetOutput(&biasSH);
    gb.Run();
    gb.SetInput(&mask);
    gb.SetOutput(&mask);
    gb.Run();
    biasSH/=mask;

    if (_debug) {
        char buffer[256];
        sprintf(buffer,"averagebiasSH%i.nii.gz",iter);
        biasSH.Write(buffer);
    }
    
    //apply to stacks
    irtkRealImage shCoeffs = _SH_coeffs/biasSH;
    
    for(i=0;i<biasSH.GetX();i++)
      for(j=0;j<biasSH.GetY();j++)
        for(k=0;k<biasSH.GetZ();k++)
	  if(_SH_coeffs(i,j,k,0)>0)
	  {
	    //SH coeff
	    for(t=0; t<nCoeff;t++)
	    {
	      c(t)=_SH_coeffs(i,j,k,t);
	      cb(t)=shCoeffs(i,j,k,t);
	    }
	    //simulate HR signal
	    s = SHT*c;
	    sb = SHT*cb;
	    for(t=0; t<nCoeff;t++)
	    {
	      bias(i,j,k,t)=log(sb(t)/s(t));
	    }
	  }
    
    
    //redistribute to slices
    ParallelNormaliseBiasDTI2 parallelNormaliseBiasDTI2(this,bias);
    parallelNormaliseBiasDTI2();
    

    if(_debug)
        cout << "done." << endl;
 /*
    // rotate directions
    int dirIndex;
    double gx,gy,gz;
    int nStacks = bias.GetT();
    irtkMatrix dirs_xyz(nStacks,3);   
    for (dirIndex = 0; dirIndex < nStacks; dirIndex++)  
    {
      gx=_directions[0][dirIndex+1];
      gy=_directions[1][dirIndex+1];
      gz=_directions[2][dirIndex+1];
      RotateDirections(gx,gy,gz,stack_transformations[dirIndex]);
      dirs_xyz(dirIndex,0) = gx;
      dirs_xyz(dirIndex,1) = gy;
      dirs_xyz(dirIndex,2) = gz;     
    }
    
    irtkSphericalHarmonics sh;
    irtkMatrix SHT = sh.SHbasis(dirs_xyz,order);
    irtkMatrix tSHT = SHT;
    tSHT.Transpose();
 
    irtkRealImage biasSH(_SH_coeffs);
    biasSH=0;
    irtkVector s(nStacks),c,ev;
    irtkMatrix w(nStacks,nStacks);
    irtkMatrix tmp;
    irtkMatrix tmp1,tmp2;
    irtkRealImage mask = _reconstructed;
    mask=0;
    
    double det;
    int i,j,k,t;
    for(i=0;i<biasSH.GetX();i++)
      for(j=0;j<biasSH.GetY();j++)
        for(k=0;k<biasSH.GetZ();k++)
	{
	  for(t=0;t<nStacks;t++)
	  {
	    s(t)=bias(i,j,k,t);
 	    w(t,t)=weights(i,j,k,t);
	  }

	  tmp = tSHT*w*SHT;
	  tmp.Eigenvalues(tmp1,ev,tmp2);
	  
	  det=1;
	  for(t=0;t<ev.Rows();t++)
	    det*=ev(t);
	  //cerr<<det<<" ";
	  if(det>0.001)
	  {
	    mask(i,j,k)=1;
	    tmp.Invert();
	    c = tmp*tSHT*w*s;
	    for(t=0;t<SHT.Cols();t++)
	    {
	      biasSH(i,j,k,t)=c(t);
	    }
	  }
	}
        if (_debug) {
        char buffer[256];
        sprintf(buffer,"averagebiasSHSHorig%i.nii.gz",iter);
        biasSH.Write(buffer);
    }
    
    
    irtkRealImage m = biasSH;
    m=0;
    irtkRealPixel *pm, *pb;
    pm = m.GetPointerToVoxels();
    pb = biasSH.GetPointerToVoxels();
    for (int i = 0; i<m.GetNumberOfVoxels();i++) {
        if(*pb!=0)
            *pm =1;
        pm++;
        pb++;
    }
    
    if (_debug) {
        char buffer[256];
        sprintf(buffer,"averagebiasSHmask%i.nii.gz",iter);
        m.Write(buffer);
    }
    
 
    
    irtkGaussianBlurring<irtkRealPixel> gb(_sigma_bias);
    gb.SetInput(&biasSH);
    gb.SetOutput(&biasSH);
    gb.Run();
    gb.SetInput(&m);
    gb.SetOutput(&m);
    gb.Run();
    biasSH/=m;

    if (_debug) {
        char buffer[256];
        sprintf(buffer,"averagebiasSH%i.nii.gz",iter);
        biasSH.Write(buffer);
        sprintf(buffer,"origSH%i.nii.gz",iter);
        _SH_coeffs.Write(buffer);
    }

    
   
    irtkRealPixel *pi, *pt;
    pi = _SH_coeffs.GetPointerToVoxels();
    pb = biasSH.GetPointerToVoxels();
    irtkRealImage b = biasSH;
    b=0;
    pt = b.GetPointerToVoxels();
    for (int i = 0; i<_SH_coeffs.GetNumberOfVoxels();i++) {
        if(*pi!=0)
	{
            *pi /=exp(-(*pb));
	    *pt = exp(-(*pb));
	}
        pi++;
        pb++;
	pt++;
    }
    char buffer[256];
    sprintf(buffer,"b%i.nii.gz",iter);
    b.Write(buffer);
    */
}

