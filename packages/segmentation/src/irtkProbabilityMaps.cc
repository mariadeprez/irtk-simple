/*

Date: 25 July 2013
Function: Processing a number of tissue probability maps together as a group

*/

#include <irtkImage.h>

#include <vector>

#include <irtkProbabilityMaps.h>

/// Constructor
irtkProbabilityMaps::irtkProbabilityMaps()
{
  _number_of_voxels = 0;
  _number_of_classes = 0;
}

/// Append a new probability map to the vector of probability maps
void irtkProbabilityMaps::AddProbabilityMap(irtkRealImage probabilitymap)
{
  if (_probabilitymaps.size() == 0) {
    _number_of_voxels = probabilitymap.GetNumberOfVoxels();
  }
  else {
    if (_number_of_voxels != probabilitymap.GetNumberOfVoxels()) {
      cerr << "Dimensions of the probability maps do not match" << endl;
      exit(1);
    }
  }
  _probabilitymaps.push_back(probabilitymap);
  _pProbabilitymaps.push_back(probabilitymap.GetPointerToVoxels());
  _number_of_classes = _probabilitymaps.size();
}

/// Create a vector of probability maps
void irtkProbabilityMaps::CreateVectorOfProbabilityMaps(int K, irtkRealImage **ppProbabilitymap)
{
  int k;
  for (k=0; k<K; k++) {
    this->AddProbabilityMap(*(ppProbabilitymap[k]));
  }
}

/// Move pointers in all probability maps to the first voxel
void irtkProbabilityMaps::MoveToFirstVoxel()
{
  int k;
  for (k=0; k<_probabilitymaps.size(); k++) {
    _pProbabilitymaps[k] = _probabilitymaps[k].GetPointerToVoxels();
  }
}

/// Move pointers in all probability maps to the next voxel
void irtkProbabilityMaps::MoveToNextVoxel()
{
  int k;
  for (k=0; k<_probabilitymaps.size(); k++) {
    _pProbabilitymaps[k]++;
  }
}

/// Compute probability map for the background
irtkRealImage irtkProbabilityMaps::ComputeBackground()
{  
  int i; //voxel index
  int k; //class index
  irtkRealImage background;
  irtkRealPixel *pBackground;
  irtkRealPixel min, max;
  
  background = _probabilitymaps[0];
  for (k=1; k<_probabilitymaps.size(); k++) {
    background += _probabilitymaps[k];
  }
  background.GetMinMax(&min, &max);
  
  pBackground = background.GetPointerToVoxels();
  for (i=0; i<_number_of_voxels; i++) {
    *pBackground = max - *pBackground;
    pBackground++;
  }
  return (background);
}

/// Normalise probability maps. Compute background first if necessary
void irtkProbabilityMaps::NormaliseProbabilityMaps()
{
  int i;
  int k;
  irtkRealPixel norm;
    
  this->MoveToFirstVoxel();
  for (i=0; i<_number_of_voxels; i++) {
    norm = 0;
    for (k=0; k<_number_of_classes; k++) {
      norm += *(_pProbabilitymaps[k]);
    }
    if (norm>0) {
      for (k=0; k<_number_of_classes; k++) {
        *(_pProbabilitymaps[k]) /= norm;
      }
    }
    //else {
    //  cerr << "Normalising factor equals 0 at i = " << i << endl;
    //}
    this->MoveToNextVoxel();
  }
}

/// Write probability map of a tissue
void irtkProbabilityMaps::WriteProbabilityMap(int k, char *filename)
{
  _probabilitymaps[k].Write(filename);
}

/// Return probability map of a tissue
irtkRealImage irtkProbabilityMaps::GetProbabilityMap(int k)
{
  return (_probabilitymaps[k]);
}

/// Return sum of probability maps
irtkRealImage irtkProbabilityMaps::SumProbabilityMaps(int k1, int k2)
{
  int i;
  irtkRealImage sum = _probabilitymaps[k1];

  irtkRealPixel *psum = sum.GetPointerToVoxels();
  this->MoveToFirstVoxel();
  for (i=0; i<_number_of_voxels; i++) {
    *psum = *_pProbabilitymaps[k1] + *_pProbabilitymaps[k2];
    psum++;
    this->MoveToNextVoxel();
  }
  return (sum);
}

irtkRealImage irtkProbabilityMaps::SumProbabilityMaps(int k1, int k2, int k3)
{
  int i;
  irtkRealImage sum = _probabilitymaps[k1];

  irtkRealPixel *psum = sum.GetPointerToVoxels();
  this->MoveToFirstVoxel();
  for (i=0; i<_number_of_voxels; i++) {
    *psum = *_pProbabilitymaps[k1] + *_pProbabilitymaps[k2] + *_pProbabilitymaps[k3];
    psum++;
    this->MoveToNextVoxel();
  }
  return (sum);
}

/// Return probability value at voxel currently pointed by the pointer of class k
irtkRealPixel irtkProbabilityMaps::GetValue(int k)
{
  return (*_pProbabilitymaps[k]);
}
  
/// Return probability value at voxel (x,y,z) of class k
irtkRealPixel irtkProbabilityMaps::GetValue(int x, int y, int z, int k)
{
  return (_probabilitymaps[k].Get(x, y, z));
}
  
/// Set probability value at voxel currently pointed by the pointer of class k
void irtkProbabilityMaps::SetValue(int k, irtkRealPixel value)
{
  *_pProbabilitymaps[k] = value;
}
  
/// Set probability value at voxel (x,y,z) of class k
void irtkProbabilityMaps::SetValue(int x, int y, int z, int k, irtkRealPixel value)
{
  _probabilitymaps[k].Put(x, y, z, value);
}
  
/// Return the number of voxels
int irtkProbabilityMaps::GetNumberOfVoxels()
{
  return (_number_of_voxels);
}
  
/// Return the number of classes
int irtkProbabilityMaps::GetNumberOfClasses()
{
  return (_number_of_classes);
}

/// Write probability value of all classes at voxel (x,y,z) to file for debugging
void irtkProbabilityMaps::Debug(int x, int y, int z)
{
  int k;
  ofstream debug ("debug.txt", ios::out | ios::app);
  
  for (k=0; k<_number_of_classes; k++) {
    debug << "tissue " << k << " " << this->GetValue(x, y, z, k) << endl;
  }
}
  
