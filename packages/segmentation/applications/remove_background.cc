/*=========================================================================

  Library   : packages/segmentation
  Module    : $RCSfile: gmm.cc,v $
  Authors   : Maria Murgasova and Daniel Rueckert
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2006
  Purpose   :
  Date      : $Date: 2007/05/24 15:22:19 $
  Version   : $Revision: 1.1 $
  Changes   : $Locker:  $
              $Log: gmm.cc,v $
              Revision 1.1  2007/05/24 15:22:19  mm3
              added command gmm

              Revision 1.1  2006/10/15 21:15:10  dr
              Imported sources


=========================================================================*/
// queue::push/pop

#include <irtkImage.h>
#include <irtkMeanShift.h>

char *output_name;
int padding = 0;
irtkGreyImage image;
bool ok;
int nBins=256;
double treshold;
bool have_treshold = false;


void usage()
{
  cerr << "Usage: remove_background [image] [output] <-padding> <-nBins> <-threshold>"<<endl;
  exit(1);
}

int main(int argc, char **argv)
{

  if (argc < 3) {
    usage();
    exit(1);
  }

  image.Read(argv[1]);
  argc--;
  argv++;

  output_name=argv[1];
  argc--;
  argv++;

  while (argc > 1){
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-threshold") == 0)){
      argc--;
      argv++;

      treshold = atof(argv[1]);
      have_treshold = true;
      cout << "Threshold = " << treshold << endl;

      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-padding") == 0)){
      argc--;
      argv++;

      padding = atof(argv[1]);
      cout << "Padding = " << padding << endl;

      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-nBins") == 0)){
      argc--;
      argv++;

      nBins = atof(argv[1]);
      cout << "nBins = " << padding << endl;

      argc--;
      argv++;
      ok = true;
    }
    if (ok == false){
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  irtkMeanShift msh(image, padding, nBins);
  msh.GenerateDensity();
  if (have_treshold) msh.SetTreshold(treshold);
  else msh.SetTreshold();
  msh.RemoveBackground();
  msh.Write(output_name);
}



