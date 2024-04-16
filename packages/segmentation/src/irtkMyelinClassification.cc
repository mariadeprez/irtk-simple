/*

Date:     7 July 2016

Function: Myelin segmentation without PV modelling

*/

#include <irtkImage.h>

#include <irtkMyelinClassification.h>

#include <irtkGaussianLikelihood.h>

/// Constructor
irtkMyelinClassification::irtkMyelinClassification(const irtkRealImage &image, const irtkRealImage &mask, const irtkRealImage &init, double m1)
{
  cerr << "Setting up input image ... ";
  _image = image;
  cerr << "done" << endl;
  
  cerr << "Setting up ROI mask ... ";
  _mask  = mask;
  cerr << "done" << endl;

  cerr << "Setting up initial myelin segmentation ... ";
  _init = init;
  cerr << "done" << endl;

  _I  = _image.GetNumberOfVoxels();
  _Dx = _image.GetX();
  _Dy = _image.GetY();
  _Dz = _image.GetZ();
  _dx = _image.GetXSize();
  _dy = _image.GetYSize();
  _dz = _image.GetZSize();
  cerr << "Image dimensions = " << _Dx << " " << _Dy << " " << _Dz << endl;
  cerr << "Voxel dimensions = " << _dx << " " << _dy << " " << _dz << endl;

  _K = 2;
  cerr << "Total number of classes = " << _K << endl;  
  _mn = new double [_K];
  _vr = new double [_K];
  _mx = new double [_K];
  
  _f     = 0;
  _denom = mask;  
  _m1    = m1;
  _mrf1  = false; // Flag for first-order MRFs
  
  if (_m1 != 0) {
    _mrf1 = true;
    cerr << "Parameter of M = " << _m1 << endl;    
  }
}


/// Destructor
irtkMyelinClassification::~irtkMyelinClassification()
{
  delete []_mn;
  delete []_vr;
  delete []_mx;
}


/// Construct 2D connectivity matrix for first-order MRFs
void irtkMyelinClassification::Construct2DConnectivityMatrix()
{
  _M = irtkMatrix(_K, _K);
  
  _M.Put(0, 0, 0);   _M.Put(0, 1, _m1);
  _M.Put(1, 0, _m1); _M.Put(1, 1, 0); 
   
  if (_mrf1 == true) {    
    cerr << "Constructing 2D connectivity matrix for first-order MRFs ... done" << endl;  
    _M.Print();     
  }
}


/// Initialise posterior probability maps
void irtkMyelinClassification::InitialisePosteriors()
{
  int i, k;

  for (k=0; k<_K; k++) {
    _posteriors.AddProbabilityMap(_init);
  }
  
  irtkRealPixel *pmask = _mask.GetPointerToVoxels();
  irtkRealPixel *pinit = _init.GetPointerToVoxels();
  _posteriors.MoveToFirstVoxel();
  
  for (i=0; i<_I; i++) {
    if (*pmask == 1) {
      if (*pinit == 1) {
        _posteriors.SetValue(0, 0);
      }
      else {
        _posteriors.SetValue(0, 1);
      }
    }
    pmask++;
    pinit++;
    _posteriors.MoveToNextVoxel();
  }  
}


/// Initialise MRF prior probability maps
void irtkMyelinClassification::InitialisePriors()
{
  _MRFpriors = _posteriors;
  
  if (_mrf1 == true) {
    cerr << "Initialising MRF prior probability maps ... done" << endl;
  } 
}


/// Initialise GMM parameters
void irtkMyelinClassification::InitialiseParams()
{
  int i, k;
  double temp = _K;
  
  double *A = new double [_K];
  double *B = new double [_K]; 
  double *C = new double [_K]; 
  double sum_of_B = 0;
  for (k=0; k<_K; k++) {
    A[k] = 0;
    B[k] = 0;
    C[k] = 0;
  }

  irtkRealPixel *pimage = _image.GetPointerToVoxels();
  irtkRealPixel *pmask  = _mask.GetPointerToVoxels();
  _posteriors.MoveToFirstVoxel();
 
  for (i=0; i<_I; i++) {
    if (*pmask == 1) {
      for (k=0; k<_K; k++) {
        A[k] += _posteriors.GetValue(k) * (*pimage);
        B[k] += _posteriors.GetValue(k);
      }
    }
    pimage++;
    pmask++;
    _posteriors.MoveToNextVoxel();
  }
  
  for (k=0; k<_K; k++) {
    if (B[k] > 0) {
      _mn[k] = A[k] / B[k];
    }
    else {
      cerr << "Dividing by 0 when initialising GMM parameters" << endl;
      exit(1);
    }
  }

  pimage = _image.GetPointerToVoxels();
  pmask  = _mask.GetPointerToVoxels();
  _posteriors.MoveToFirstVoxel();
  
  for (i=0; i<_I; i++) {
    if (*pmask == 1) {
      for (k=0; k<_K; k++) {
        C[k] += _posteriors.GetValue(k) * (*pimage - _mn[k]) * (*pimage - _mn[k]);
      }
    }
    pimage++;
    pmask++;
    _posteriors.MoveToNextVoxel();
  }
  
  for (k=0; k<_K; k++) {
    sum_of_B += B[k];
  }
  for (k=0; k<_K; k++) {
    _vr[k] = C[_K-1] / B[_K-1]; // Variances of all classes are equal to that of myelin
    _mx[k] = B[k] / sum_of_B;   
  }
  
  delete []A;
  delete []B;
  delete []C;
}


/// Return sum of neighbour probabilities weighted respectively by inverse Euclidean distances at voxel (x,y,z) of class k 
irtkRealPixel irtkMyelinClassification::GetNeighbourhood(int x, int y, int z, int k)
{
  int nx, ny, nz;
  double s;
  irtkRealPixel sum = 0;
    
  for (nx=-1; nx<2; nx++) {
    for (ny=-1; ny<2; ny++) {
      for (nz=-1; nz<2; nz++) {
	    if ((nx!=0) || (ny!=0) || (nz!=0)) {
	      s = 1 / sqrt(_dx*_dx*nx*nx + _dy*_dy*ny*ny + _dz*_dz*nz*nz);
     	  sum += _posteriors.GetValue(x+nx, y+ny, z+nz, k) * s;
	    }
      }
    }
  }
  return (sum);
}


/// Return first-order penalty at voxel (x,y,z) of class k
irtkRealPixel irtkMyelinClassification::GetFirstOrderPenalty(int x, int y, int z, int k)
{
  int k1;
  irtkRealPixel U1 = 0;
  
  for (k1=0; k1<_K; k1++) {
    if (_M.Get(k, k1)!=0) {
      U1 += _M.Get(k, k1) * (this->GetNeighbourhood(x, y, z, k1));
    }
  }
  return (U1); 
}


/// EM algorithm
void irtkMyelinClassification::EStep()
{
  int x, y, z, k;
  irtkProbabilityMaps posteriors = _posteriors;

  irtkRealPixel *U1       = new irtkRealPixel [_K];
  irtkRealPixel *MRFnumer = new irtkRealPixel [_K];  
  irtkRealPixel *numer    = new irtkRealPixel [_K];
  irtkRealPixel MRFdenom, denom;
   
  irtkGaussianLikelihood *G = new irtkGaussianLikelihood [_K];
  for (k=0; k<_K; k++) {
    G[k].Initialise(_mn[k], _vr[k]);
  }
    
  for (z=1; z<_Dz-1; z++) {
    for (x=1; x<_Dx-1; x++) {
      for (y=1; y<_Dy-1; y++) {	
        if (_mask.Get(x, y, z) == 1) {
        
          if (_mrf1 == true) {
            MRFdenom = 0;
            for (k=0; k<_K; k++) {
              U1[k]       = this->GetFirstOrderPenalty(x, y, z, k);
              MRFnumer[k] = exp(- U1[k]);
              MRFdenom   += MRFnumer[k];
            }
            if (MRFdenom > 0) {
              for (k=0; k<_K; k++) {
                _MRFpriors.SetValue(x, y, z, k, MRFnumer[k]/MRFdenom);
              }
            }
            else {
              cerr << "Dividing by 0 when computing MRF prior at voxel (" << x << " " << y << " " << z << ")"<< endl;
              for (k=0; k<_K; k++) {
                cerr << "U1[" << k << "] = " << U1[k] << endl;
              }
            }
          }
          
          denom = 0;
          for (k=0; k<_K; k++) {
            if (_mrf1 == true) {
              numer[k] = G[k].Evaluate(_image.Get(x, y, z)) * _MRFpriors.GetValue(x, y, z, k);
            }
            else {
	          numer[k] = G[k].Evaluate(_image.Get(x, y, z)) * _mx[k]; // Pure GMM
	        }
            denom += numer[k];
          }          
          if (denom > 0) {
            _denom.Put(x, y, z, denom);
            for (k=0; k<_K; k++) {
              posteriors.SetValue(x, y, z, k, numer[k]/denom);
            }
          }
          else {
            posteriors.SetValue(x, y, z, 0, 1); // Background tissue
            for (k=1; k<_K; k++) {
              posteriors.SetValue(x, y, z, k, 0);
            }
          }
        } 
      }
    }
    this->PrintProgress(z);
  }
  
  _posteriors = posteriors;
  delete []U1;
  delete []MRFnumer;
  delete []numer;
  delete []G;  
}

 
void irtkMyelinClassification::MStep()
{
  int i, k;
  
  double *A = new double [_K];
  double *B = new double [_K];
  double *C = new double [_K];
  double sum_of_B = 0;
  double sum_of_C = 0;
  for (k=0; k<_K; k++) {
    A[k]  = 0;
    B[k]  = 0;
    C[k]  = 0;
  }
  
  irtkRealPixel *pimage = _image.GetPointerToVoxels();
  irtkRealPixel *pmask  = _mask.GetPointerToVoxels();
  _posteriors.MoveToFirstVoxel();
 
  for (i=0; i<_I; i++) {
    if (*pmask == 1) {
      for (k=0; k<_K; k++) {
        A[k] += _posteriors.GetValue(k) * (*pimage);
        B[k] += _posteriors.GetValue(k);
      }
    }
    pimage++;
    pmask++;
    _posteriors.MoveToNextVoxel();
  }
   
  for (k=0; k<_K; k++) {
    if (B[k] > 0) {
      _mn[k] = A[k] / B[k];
    }
    else {
      cerr << "Dividing by 0 in M-step of class " << k << endl;
      exit(1);
    }
  }
  
  pimage = _image.GetPointerToVoxels();
  pmask  = _mask.GetPointerToVoxels();
  _posteriors.MoveToFirstVoxel();
 
  for (i=0; i<_I; i++) {
    if (*pmask == 1) {
      for (k=0; k<_K; k++) {
        C[k] += _posteriors.GetValue(k) * (*pimage - _mn[k]) * (*pimage - _mn[k]);
      }
    }
    pimage++;
    pmask++;
    _posteriors.MoveToNextVoxel();
  }
  
  for (k=0; k<_K; k++) {
    sum_of_B += B[k];
    sum_of_C += C[k];
  }  
  for (k=0; k<_K; k++) {
    _vr[k] = sum_of_C / sum_of_B;
    _mx[k] = B[k] / sum_of_B;
  }

  delete []A;
  delete []B;
  delete []C;
}


/// Return relative decrease of objective function
double irtkMyelinClassification::LogLikelihood()
{
  int i;
  double f, rel_diff;
  ofstream Log ("param.txt", ios::out | ios::app);

  irtkRealPixel *pdenom = _denom.GetPointerToVoxels();
  irtkRealPixel *pmask  = _mask.GetPointerToVoxels();

  f = 0;
  for (i=0; i<_I; i++) {
    if (*pmask == 1) {
      f += log(*pdenom);
    }
    pdenom++;
    pmask++;
  }
  
  f = - f;
  if (_f == 0) {
    rel_diff = 1;
  }
  else {
    rel_diff = (_f - f) / _f;
  }
  _f = f;
  cerr << "f = " << _f << " rel_diff = " << rel_diff << endl;
  Log << "f = " << _f << " rel_diff = " << rel_diff << endl;

  return (rel_diff);
}


/// Print GMM parameters to screen
void irtkMyelinClassification::PrintToScreen()
{
  int k;
  
  cerr << "Mean: " << endl;
  for (k=0; k<_K; k++) {
    cerr << "class " << k << " " << _mn[k] << endl;
  }
  cerr << "Standard deviation: " << endl;
  for (k=0; k<_K; k++) {
    cerr << "class " << k << " " << sqrt(_vr[k]) << endl;
  }
  cerr << "Mixing coefficient: " << endl;
  for (k=0; k<_K; k++) {
    cerr << "class " << k << " " << _mx[k] << endl;
  }
}

 
/// Print GMM parameters to file
void irtkMyelinClassification::PrintToFile()
{
  int k;
  ofstream Log ("param.txt", ios::out | ios::app);
  
  Log << "Mean: " << endl;
  for (k=0; k<_K; k++) {
    Log << "class " << k << " " << _mn[k] << endl;
  }
  Log << "Standard deviation: " << endl;
  for (k=0; k<_K; k++) {
    Log << "class " << k << " " << sqrt(_vr[k]) << endl;
  }
  Log << "Mixing coefficient: " << endl;
  for (k=0; k<_K; k++) {
    Log << "class " << k << " " << _mx[k] << endl;
  }
}


/// Execute one iteration of EM algorithm and return relative decrease of objective function
double irtkMyelinClassification::Iterate()
{
  cerr << "EStep ... ";
  this->EStep();
  cerr << "done" << endl;
  
  cerr << "MStep ... ";
  this->MStep();
  cerr << "done" << endl;
  
  this->PrintToScreen();
  this->PrintToFile();
  
  cerr << "Computing objective function ... " << endl;
  return (this->LogLikelihood());
}
  

/// Calculate final myelin segmentation
void irtkMyelinClassification::GetMyelinSegment()
{
  int i, k;
  double max; int idx;
  
  _hard_myelin = _mask;

  irtkRealPixel *pmask  = _mask.GetPointerToVoxels();
  irtkRealPixel *pseg   = _hard_myelin.GetPointerToVoxels();
  _posteriors.MoveToFirstVoxel();
  
  for (i=0; i<_I; i++) {
    if (*pmask == 1) {      
      max = _posteriors.GetValue(0);
      idx = 0;     
      for (k=0; k<_K; k++) {
        if (_posteriors.GetValue(k) > max) {
          max = _posteriors.GetValue(k);
          idx = k;
        }
      }
      if (idx == 0) {
        *pseg = 0;
      }
      if (idx == _K-1) {
        *pseg = 1;
      }      
    }
    pmask++;
    pseg++;
    _posteriors.MoveToNextVoxel();
  }
  _hard_myelin.Write("hard-mye.nii.gz");  
}


/// Write posterior probability maps
void irtkMyelinClassification::WritePosteriors()
{
  _posteriors.WriteProbabilityMap(0, "posterior-bkg.nii.gz"); // Background tissue
  _posteriors.WriteProbabilityMap(1, "posterior-mye.nii.gz"); // Myelin
}


/// Print progress
void irtkMyelinClassification::PrintProgress(int x)
{
  cerr << x << "...";
}
