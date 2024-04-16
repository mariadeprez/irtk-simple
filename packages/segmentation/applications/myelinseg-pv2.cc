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

#include <irtkMultiChannelImage.h>

#include <string.h>

#include <stdio.h>

using namespace std;

char *output_name;

void usage()
{
  cerr << "Usage: myelinseg-pv2 [image] [mask] [percentile] [t1] [t2] [t3] [iter] " << endl;  
  exit(1);
}

int main(int argc, char **argv)
{
  int iter, Iterations;
  double t1, t2, t3;

  if (argc < 8) {
    usage();
    exit(1);
  }
  
  /// Read input image
  irtkRealImage image;
  cerr << "Reading input image " << argv[1] << " ... ";
  image.Read(argv[1]);
  output_name = "tmp.nii.gz";
  cerr << output_name << endl;
  cerr << "done" << endl;
  argc--;
  argv++;
  
  /// Read mask of ROI
  irtkRealImage mask;
  cerr << "Reading mask of ROI " << argv[1] << " ... ";
  mask.Read(argv[1]);
  cerr << "done" << endl;
  argc--;
  argv++; 

  /// Read Percentile: 
  double percentile;
  cerr << "Reading Percentile " << argv[1] << " ... ";
  percentile = atof(argv[1]);
  argc--;
  argv++;
  
  
  cerr<<"done input"<<endl;
  
  /// Read initial myelin segmentation
  int padding = 0;
  irtkRealImage init_tmp, init;
  irtkMultiChannelImage mch;
  mch.SetPadding(padding);
  mch.AddImage(image);
  mch.SetMask(mask);
  mch.Brainmask();
  mch.Write(0,output_name);

  init_tmp.Read(output_name);
  init.Read(output_name);
  vector<double> sort_arr;
  cerr << "Init generated" << endl;
  int x, y, z;
  int count = 0;
  int thrsh_pos;
	for (z = 0; z < init_tmp.GetZ(); z++) {
		for (y = 0; y < init_tmp.GetY(); y++) {
			for (x = 0; x < init_tmp.GetX(); x++) {
				if (init_tmp.Get(x,y,z) > 0){
					sort_arr.push_back(init_tmp.Get(x,y,z));
					count++;
				}				
			}
		}
	}
  thrsh_pos = percentile*count;
  cerr << "Position at: " << thrsh_pos << endl;
  cerr << "Loop 1 finished" << endl;
  sort(sort_arr.begin(),sort_arr.end());
  float thrsh_num;
  thrsh_num = sort_arr[thrsh_pos];
  cerr << thrsh_num << endl;
  for (int z = 0; z < init_tmp.GetZ(); z++) {
		for (int y = 0; y < init_tmp.GetY(); y++) {
			for (int x = 0; x < init_tmp.GetX(); x++) {
                            if (init_tmp.Get(x,y,z) > 0){
				if (init_tmp.Get(x, y, z) <= thrsh_num) {
                                    init.Put(x, y, z, 1);
				} else {
                                    init.Put(x, y, z, 0);
				}
                            }
			}
		}
	}
  init.Write("init.nii.gz");

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
  } 
  while ((rel_diff > 0.0001) && (iter < Iterations));

  cerr << "Calculating final myelin segmentation ... ";
  myelin->GetMyelinSegment();
  cerr << "done" << endl;
  
  myelin->WritePosteriors();

  delete myelin;
  //if( remove( output_name ) != 0 )
  //  perror( "Error deleting file" );
  //else
  //  puts( "File successfully deleted" );
  //return 0;
}
