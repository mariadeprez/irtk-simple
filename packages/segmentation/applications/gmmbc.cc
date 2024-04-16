/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id: gmm-par.cc 2 2008-12-23 12:40:14Z dr $
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date: 2008-12-23 12:40:14 +0000 (Tue, 23 Dec 2008) $
  Version   : $Revision: 2 $
  Changes   : $Author: dr $

=========================================================================*/

#include <irtkImage.h>
#include <irtkEMClassification.h>
#include <irtkEMClassificationBiasCorrection.h>

char *output_name;
char *bias_name;
double *mean, *var, *c;
irtkRealImage * background;
// Default parameters
double treshold = 0.001;
int iterations = 50;
int padding    = 0;
double sigma = 30;
bool _gd=false;

void usage()
{
  cerr << "Usage: gmm [image] [output] [n] <-mean mean1 ... meann> <-variance s1 ... sn> <-c c1 ... cn> " << endl;
  cerr << "<-iterations> <-padding> <-background> <-treshold> <-bias filename>" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int i, n, ok;

  if (argc < 3) {
    usage();
    exit(1);
  }

  // Input image
  irtkRealImage image;
  image.Read(argv[1]);
  argc--;
  argv++;

  //output segmentation
  output_name = argv[1];
  argc--;
  argv++;

  // Number of tissues
  n = atoi(argv[1]);
  argc--;
  argv++;

  //Default settings for Gaussian Mixture parameters
  mean = new double[n];
  var  = new double[n];
  c    = new double[n];

  irtkRealPixel min, max;
  image.GetMinMax(&min,&max);

  for (i=0;i<n;i++)  mean[i] = min + i*(max-min)/(double) n;
  for (i=0;i<n;i++)  var[i] = ((max-min)/(double) n)*((max-min)/(double) n)/16;
  for (i=0;i<n;i++)  c[i] = 1.0/(double) n;


  // Parse remaining parameters
  while (argc > 1) {
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-mean") == 0)) {
      argc--;
      argv++;

      for (i=0;i<n;i++) {
        mean[i] = atof(argv[1]);
        argc--;
        argv++;
      }
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-variance") == 0)) {
      argc--;
      argv++;

      for (i=0;i<n;i++) {
        var[i] = atof(argv[1])*atof(argv[1]);
        argc--;
        argv++;
      }
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-c") == 0)) {
      argc--;
      argv++;

      for (i=0;i<n;i++) {
        c[i] = atof(argv[1]);
        argc--;
        argv++;
      }
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-iterations") == 0)) {
      argc--;
      argv++;
      iterations = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-padding") == 0)) {
      argc--;
      argv++;
      padding = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-threshold") == 0)) {
      argc--;
      argv++;
      treshold = atof(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-background") == 0)) {
      argc--;
      argv++;
      background = new irtkRealImage;
      background->Read(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-sigma") == 0)) {
      argc--;
      argv++;
      sigma = atof(argv[1]);
      cerr << "sigma  = " << sigma <<endl;
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-bias") == 0)) {
      argc--;
      argv++;
      bias_name = argv[1];
      cerr << "bias name  = " << sigma <<endl;
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-gd") == 0)) {
      argc--;
      argv++;
      _gd = true;
      cerr << "Bias estimation with gradient descent" <<endl;
      ok = true;
    }

    if (ok == false) {
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  irtkEMClassificationBiasCorrection classification(sigma);
  classification.SetInput(image);
  classification.SetPadding(padding);
  classification.CreateMask();
  classification.InitialiseGMMParameters(n,mean,var,c);

  double rel_diff;
  i=1;
  //treshold = 0.00001;
  do {
    cout << "Iteration = " << i << " / " << iterations << endl;
    if(_gd)
      rel_diff = classification.IterateGMMGD(i,false, false);
    else
      rel_diff = classification.IterateGMM(i,false, false);
    i++;
  } 
  while ((rel_diff>treshold)&&(i<iterations));
  //while (i<iterations);


  classification.WriteGaussianParameters("parameters.txt");

  irtkRealImage segmentation;
  classification.ConstructSegmentationNoBG(segmentation);
  segmentation.Write(output_name);
  
  if(bias_name != NULL)
  {
    classification.WriteBias(bias_name);
  }
}

