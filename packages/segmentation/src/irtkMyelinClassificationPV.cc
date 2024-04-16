/*

Date:     7 July 2016

Function: Myelin segmentation with PV modelling

*/

#include <irtkImage.h>

#include <irtkMyelinClassificationPV.h>

#include <irtkGaussianLikelihood.h>


/// Constructor
irtkMyelinClassificationPV::irtkMyelinClassificationPV(const irtkRealImage &image, const irtkRealImage &mask, const irtkRealImage &init, double t1, double t2, double t3)
{
  double temp;
  
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

  _N = 1; // Number of PV classes
  _K = _N + 2;
  cerr << "Total number of classes = " << _K << endl;
  
  temp    = _N;
  _lambda = 1 / (temp + 1);
  
  _mn = new double [_K];
  _vr = new double [_K];
  _mx = new double [_K];
 
  _f     = 0;
  _denom = mask;

  _t1 = t1; _t2 = t2; _t3 = t3;
  _mrf2 = false; // Flag for second-order MRFs
  
  if ((_t1 + _t2 + _t3) != 0) {
    _mrf2 = true;
    cerr << "Parameters of T = " << _t1 << " " << _t2 << " " << _t3 << endl;
  }
}


/// Destructor
irtkMyelinClassificationPV::~irtkMyelinClassificationPV()
{
  delete []_mn;
  delete []_vr;
  delete []_mx;
}


/// Construct 3D connectivity tensor for second-order MRFs
void irtkMyelinClassificationPV::Construct3DConnectivityTensor()
{ 
  int k1, k2, n;
  irtkMatrix Tk = irtkMatrix(_K, _K);
 
  // Background tissue
  for (k1=0; k1<_K; k1++) {
    for (k2=0; k2<_K; k2++) {
      Tk.Put(k1, k2, 0);
    }
  }
  for (k1=0; k1<_K; k1++) {
    for (k2=0; k2<_K; k2++) {
      if ((k1==_N+1) && (k2==_N+1)) {
        Tk.Put(k1, k2, _t1);
      }
      if (((k1==0) && (k2==_N+1)) || ((k1==_N+1) && (k2==0))) {
        Tk.Put(k1, k2, _t2);
      }
    }
  }
  _T.push_back(Tk);
  if (_mrf2 == true) {
    cerr << "3D connectivity tensor - background tissue:" << endl;
    Tk.Print();
  }
  
  // PV classes
  for (n=1; n<=_N; n++) {
    for (k1=0; k1<_K; k1++) {
      for (k2=0; k2<_K; k2++) {
        Tk.Put(k1, k2, 0);
      }
    }
    for (k1=0; k1<_K; k1++) {
      for (k2=0; k2<_K; k2++) {
        if ((k1==0) && (k2==0)) {
          Tk.Put(k1, k2, _t3);
        }
        if (((k1==0) && (k2>=1) && (k2<=_N)) || ((k1>=1) && (k1<=_N) && (k2==0))) {
          Tk.Put(k1, k2, 0);
        }
        if ((k1==_N+1) && (k2==_N+1)) {
          Tk.Put(k1, k2, 0);
        }
        if (((k1==_N+1) && (k2>=1) && (k2<=_N)) || ((k1>=1) && (k1<=_N) && (k2==_N+1))) {
          Tk.Put(k1, k2, 0);
        }      
      }
    }
    _T.push_back(Tk);
  if (_mrf2 == true) {
      cerr << "3D connectivity tensor - PV classes:" << endl;
      Tk.Print();
    }
  }
  
  // Myelin
  for (k1=0; k1<_K; k1++) {
    for (k2=0; k2<_K; k2++) {
      Tk.Put(k1, k2, 0);
    }
  }
  for (k1=0; k1<_K; k1++) {
    for (k2=0; k2<_K; k2++) {
      if ((k1==0) && (k2==0)) {
        Tk.Put(k1, k2, _t1);
      }
      if (((k1==0) && (k2==_N+1)) || ((k1==_N+1) && (k2==0))) {
        Tk.Put(k1, k2, _t2);
      }    
    }
  }
  _T.push_back(Tk);
  if (_mrf2 == true) {
    cerr << "3D connectivity tensor - myelin:" << endl;
    Tk.Print();
  }
}


/// Initialise posterior probability maps
void irtkMyelinClassificationPV::InitialisePosteriors()
{
  int i, k, n;
   
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
      for (n=1; n<=_N; n++) {
        _posteriors.SetValue(n, 0); // PV classes
      }
    }
    pmask++;
    pinit++;
    _posteriors.MoveToNextVoxel();
  }  
}


/// Initialise MRF prior probability maps
void irtkMyelinClassificationPV::InitialisePriors()
{
  _MRFpriors = _posteriors;
  
  if (_mrf2 == true) {
    cerr << "Initialising MRF prior probability maps ... done" << endl;
  } 
}


/// Initialise GMM parameters
void irtkMyelinClassificationPV::InitialiseParams()
{
  int i, k, n;
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
  
  if ((B[0]>0) && (B[_N+1]>0)) {
    _mn[0]    = A[0] / B[0];       // Background tissue
    _mn[_N+1] = A[_N+1] / B[_N+1]; // Myelin
  }
  else {
    cerr << "Dividing by 0 when initialising GMM parameters" << endl;
    exit(1);
  }
  for (n=1; n<=_N; n++) {
    _mn[n] = (1 - n * _lambda) * _mn[0] + n * _lambda * _mn[_N+1]; // PV classes
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
    _vr[k] = C[_N+1] / B[_N+1]; // Variances of all classes are equal to that of myelin
    _mx[k] = B[k] / sum_of_B;
  }
  
  delete []A;
  delete []B;
  delete []C;
}


/// Return sum of neighbour probabilities weighted respectively by inverse Euclidean distances at voxel (x,y,z) of class k 
irtkRealPixel irtkMyelinClassificationPV::GetNeighbourhood(int x, int y, int z, int k)
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


/// Return second-order penalty at voxel (x,y,z) of class k 
irtkRealPixel irtkMyelinClassificationPV::GetSecondOrderPenalty(int x, int y, int z, int k)
{
  int k1, k2;
  irtkRealPixel U2 = 0;
  irtkRealPixel tmp;
  
  for (k1=0; k1<_K; k1++) {
    tmp = 0;
    for (k2=0; k2<_K; k2++) {
      if (_T[k].Get(k1, k2)!=0) {
        tmp += _T[k].Get(k1, k2) * (this->GetNeighbourhood(x, y, z, k2));
      }
    }
    U2 += (this->GetNeighbourhood(x, y, z, k1)) * tmp;
  }
  return (U2);
}


/// EM algorithm
void irtkMyelinClassificationPV::EStep()
{
  int x, y, z, k;
  irtkProbabilityMaps posteriors = _posteriors;

  irtkRealPixel *U2       = new irtkRealPixel [_K];
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
        
          if (_mrf2 == true) {
            MRFdenom = 0;
            for (k=0; k<_K; k++) {
              U2[k] = this->GetSecondOrderPenalty(x, y, z, k);      
              MRFnumer[k] = exp(- U2[k]);
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
                cerr << "U2[" << k << "] = " << U2[k] << endl;
              }
              exit(1);
            }
          }
           
          denom = 0;
          for (k=0; k<_K; k++) {
            if (_mrf2 == true) {
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
  delete []U2;
  delete []MRFnumer;
  delete []numer;
  delete []G;  
}


void irtkMyelinClassificationPV::MStep()
{
  int i, k, n;
  
  double *A = new double [_K];
  double *B = new double [_K];
  double *C = new double [_K];
  double *a = new double [2];
  double *b = new double [3];
  double sum_of_B = 0;
  double sum_of_C = 0;  

  a[0] = 0; a[1] = 0;
  b[0] = 0; b[1] = 0; b[2] = 0;
  for (k=0; k<_K; k++) {
    A[k]  = 0;
    B[k]  = 0;
    C[k]  = 0;
  }
  
  irtkMatrix Y = irtkMatrix(2, 2); // Yx = Z solves the means of BKG and MYE as 2 independent variables
  irtkMatrix x = irtkMatrix(2, 1);
  irtkMatrix Z = irtkMatrix(2, 1);
    
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
  
  for (n=1; n<=_N; n++) {
    a[0] += (1 - n * _lambda) * A[n];
    a[1] += n * _lambda * A[n];
    b[0] += (1 - n * _lambda) * (1 - n * _lambda) * B[n];
    b[1] += n * _lambda * (1 - n * _lambda) * B[n];
    b[2] += n * _lambda * n * _lambda * B[n];
  }
  
  Y.Put(0, 0, B[0] + b[0]);
  Y.Put(0, 1, b[1]);
  Y.Put(1, 0, b[1]);
  Y.Put(1, 1, B[_N+1] + b[2]);
  Z.Put(0, 0, A[0] + a[0]);
  Z.Put(1, 0, A[_N+1] + a[1]);
  
  if (Y.Det() > 0) {
    Y.Invert();
    x = Y * Z;
    _mn[0]    = x.Get(0, 0); // Background tissue    
    _mn[_N+1] = x.Get(1, 0); // Myelin
    for (n=1; n<=_N; n++) {
      _mn[n] = (1 - n * _lambda) * _mn[0] + n * _lambda * _mn[_N+1]; // PV classes
    }
  }
  else {
  cerr << "Cannot find inverse matrix" << endl;
  exit(1);
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
  delete []a;
  delete []b;
}


/// Return relative decrease of objective function
double irtkMyelinClassificationPV::LogLikelihood()
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
void irtkMyelinClassificationPV::PrintToScreen()
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
void irtkMyelinClassificationPV::PrintToFile()
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
double irtkMyelinClassificationPV::Iterate()
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
void irtkMyelinClassificationPV::GetMyelinSegment()
{
  int i, k;
  double max; int idx;
  
  _myelin_soft_segment = _mask;

  irtkRealPixel *pimage = _image.GetPointerToVoxels();
  irtkRealPixel *pmask  = _mask.GetPointerToVoxels();
  irtkRealPixel *pfrac  = _myelin_soft_segment.GetPointerToVoxels();
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
        *pfrac = 0;
      }
      if (idx == _N+1) {
        *pfrac = 1;
      }      
      if ((idx!=0) && (idx!=_N+1)) {      
        if (*pimage >= _mn[0]) {
          *pfrac = 0;
        }
        if (*pimage <= _mn[_N+1]) {
          *pfrac = 1;
        }
        if ((*pimage < _mn[0]) && (*pimage > _mn[_N+1])) {
          *pfrac = (_mn[0] - *pimage) / (_mn[0] - _mn[_N+1]);
        }
      }
    }
    pimage++;
    pmask++;
    pfrac++;
    _posteriors.MoveToNextVoxel();
  }
  _myelin_soft_segment.Write("myelin-soft-segment.nii.gz");  
}


/// Write posterior probability maps
void irtkMyelinClassificationPV::WritePosteriors()
{  
  _posteriors.WriteProbabilityMap(0, "posterior-bkg.nii.gz"); // Background tissue
  _posteriors.WriteProbabilityMap(1, "posterior-pv.nii.gz");  // PV
  _posteriors.WriteProbabilityMap(2, "posterior-mye.nii.gz"); // Myelin
}


/// Print progress
void irtkMyelinClassificationPV::PrintProgress(int x)
{
  cerr << x << "...";
}