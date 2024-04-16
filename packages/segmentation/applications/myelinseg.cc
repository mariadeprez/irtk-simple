/*

Date:     7 July 2016

Function: Myelin segmentation without PV modelling

*/

#include <irtkImage.h>

#include <irtkGaussianLikelihood.h>

#include <irtkProbabilityMaps.h>

#include <irtkMyelinClassification.h>

#include <irtkMatrix.h>

#include <vector>

void usage()
{
  cerr << "Usage: myelinseg [image] [mask] [initial] [m1] [iter]" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int iter, Iterations;
  double m1;
  
  if (argc < 5) {
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
  
  /// Read ROI mask
  irtkRealImage mask;
  cerr << "Reading ROI mask ... ";
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

  /// Read MRF parameter
  cerr << "Reading MRF parameter ... ";
  m1 = atof(argv[1]);
  argc--;
  argv++;
  cerr << "done" << endl;
         
  /// Read number of iterations
  cerr << "Reading number of iterations ... ";
  Iterations = atoi(argv[1]);
  argc--;
  argv++;
  cerr << "done" << endl;
  
  irtkMyelinClassification *myelin = new irtkMyelinClassification(image, mask, init, m1);

  cerr << "Initialising posterior probability maps ... ";  
  myelin->InitialisePosteriors();
  cerr << "done" << endl;
    
  myelin->InitialisePriors();
  myelin->Construct2DConnectivityMatrix();  
  
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