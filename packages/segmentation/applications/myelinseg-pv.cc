/*

Date:     7 July 2016

Function: Myelin segmentation with PV modelling

*/

#include <irtkImage.h>

#include <irtkGaussianLikelihood.h>

#include <irtkProbabilityMaps.h>

#include <irtkMyelinClassificationPV.h>

#include <irtkMatrix.h>

#include <vector>

void usage()
{
  cerr << "Usage: myelinseg-pv [image] [mask] [initial] [t1] [t2] [t3] [iter]" << endl;  
  exit(1);
}

int main(int argc, char **argv)
{
  int iter, Iterations;
  double t1, t2, t3;

  if (argc < 7) {
    usage();
    exit(1);
  }
  
  /// Read input image
  irtkRealImage image;
  cerr << "Reading input image ... ";
  image.Read(argv[1]);
  cerr << "done" << endl;
  argc--;
  argv++;
  
  /// Read mask of ROI
  irtkRealImage mask;
  cerr << "Reading mask of ROI ... ";
  mask.Read(argv[1]);
  cerr << "done" << endl;
  argc--;
  argv++; 
  
  /// Read initial myelin segmentation
  irtkRealImage init;
  cerr << "Reading initial myelin segmentation ... ";
  init.Read(argv[1]);
  cerr << "done" << endl;
  argc--;
  argv++;

  /// Read MRF parameters
  cerr << "Reading MRF parameters ... ";
  t1 = atof(argv[1]);
  argc--;
  argv++;
  t2 = atof(argv[1]);
  argc--;
  argv++;
  t3 = atof(argv[1]);
  argc--;
  argv++;
  cerr << "done" << endl;
  
  /// Read number of iterations
  cerr << "Reading number of iterations ... ";
  Iterations = atoi(argv[1]);
  argc--;
  argv++;
  cerr << "done" << endl;

  irtkMyelinClassificationPV *myelin = new irtkMyelinClassificationPV(image, mask, init, t1, t2, t3);
  
  cerr << "Initialising posterior probability maps ... ";
  myelin->InitialisePosteriors();
  cerr << "done" << endl;

  myelin->InitialisePriors();
  myelin->Construct3DConnectivityTensor(); 
  
  cerr << "Initialising GMM parameters ... ";
  myelin->InitialiseParams();
  cerr << "done" << endl;  
  myelin->PrintToScreen();
  myelin->PrintToFile();

  double rel_diff;
  iter = 0;

  do {
    cerr << "Iteration " << iter << endl;
    rel_diff = myelin->Iterate();
    iter++;
  } while ((rel_diff > 0.0001) && (iter < Iterations));

  cerr << "Calculating final myelin segmentation ... ";
  myelin->GetMyelinSegment();
  cerr << "done" << endl;
  
  myelin->WritePosteriors();

  delete myelin;
}