/*

Date: 25 July 2013
Function: Processing a number of tissue probability maps together as a group

*/

#ifndef _IRTKPROBABILITYMAPS_H

#define _IRTKPROBABILITYMAPS_H

#include <irtkObject.h>

#include <irtkImage.h>

#include <vector>

class irtkProbabilityMaps: public irtkObject
{

protected:

  /// Vector of probability maps
  vector<irtkRealImage> _probabilitymaps;
  
  /// Vector of pointers to the voxels of the probability maps
  vector<irtkRealPixel *> _pProbabilitymaps;

  /// Number of voxels
  int _number_of_voxels;
    
  /// Number of classes
  int _number_of_classes;
  
public:

  /// Constructor
  irtkProbabilityMaps();
  
  /// Append a new probability map to the vector of probability maps
  void AddProbabilityMap(irtkRealImage probabilitymap);
  
  /// Create a vector of probability maps
  void CreateVectorOfProbabilityMaps(int K, irtkRealImage **ppProbabilitymap);
  
  /// Move pointers in all probability maps to the first voxel
  void MoveToFirstVoxel();
  
  /// Move pointers in all probability maps to the next voxel
  void MoveToNextVoxel();  

  /// Compute probability map for the background
  irtkRealImage ComputeBackground();
  
  /// Normalise probability maps. Compute background first if necessary
  void NormaliseProbabilityMaps();
  
  /// Write probability map of a tissue
  void WriteProbabilityMap(int k, char *filename);
  
  /// Return probability map of a tissue
  irtkRealImage GetProbabilityMap(int k);
  
  /// Return sum of probability maps
  irtkRealImage SumProbabilityMaps(int k1, int k2);
  
  irtkRealImage SumProbabilityMaps(int k1, int k2, int k3);

  /// Return probability value at voxel currently pointed by the pointer of class k
  irtkRealPixel GetValue(int k);
  
  /// Return probability value at voxel (x,y,z) of class k
  irtkRealPixel GetValue(int x, int y, int z, int k);
  
  /// Set probability value at voxel currently pointed by the pointer of class k
  void SetValue(int k, irtkRealPixel value);
  
  /// Set probability value at voxel (x,y,z) of class k
  void SetValue(int x, int y, int z, int k, irtkRealPixel value);
  
  /// Return the number of voxels
  int GetNumberOfVoxels();
  
  /// Return the number of classes
  int GetNumberOfClasses();
  
  /// Write probability value of all classes at voxel (x,y,z) to file for debugging
  void Debug(int x, int y, int z);
  
};

#endif
  
  