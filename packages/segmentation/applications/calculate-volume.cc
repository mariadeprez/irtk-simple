/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id: reconstructionb0.cc 1000 2013-10-18 16:49:08Z mm3 $
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2011 onwards
  Date      : $Date: 2013-10-18 17:49:08 +0100 (Fri, 18 Oct 2013) $
  Version   : $Revision: 1000 $
  Changes   : $Author: mm3 $

=========================================================================*/

#include <irtkImage.h>
#include <irtkTransformation.h>
#include <vector>
#include <string>
using namespace std;

//Application to perform reconstruction of volumetric MRI from thick slices.

void usage()
{
  cerr << "Usage: calculate-volume [image] <-mask> <-file> <-padding>\n" << endl;
  cerr << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  
  int ok;
  char buffer[256];
  irtkRealImage image, m;  
  irtkRealImage *mask=NULL;
  char *file_name=NULL;
  double padding=0;
    
  //if not enough arguments print help
  if (argc < 2)
    usage();
  
  image.Read(argv[1]);  
  argc--;
  argv++;
  
  // Parse options.
  while (argc > 1){
    ok = false;
    

    //Read binary mask for final volume
    if ((ok == false) && (strcmp(argv[1], "-mask") == 0)){
      argc--;
      argv++;
      mask= new irtkRealImage(argv[1]);
      ok = true;
      argc--;
      argv++;
    }
    
    if ((ok == false) && (strcmp(argv[1], "-file") == 0)){
      argc--;
      argv++;

      file_name = argv[1];

      argc--;
      argv++;
      ok = true;
    }

    if ((ok == false) && (strcmp(argv[1], "-padding") == 0)){
      argc--;
      argv++;

      padding = atof(argv[1]);

      argc--;
      argv++;
      ok = true;
    }
    
    if (ok == false){
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }
  

  if(mask!=NULL)
  {
    if(image.GetNumberOfVoxels()!=mask->GetNumberOfVoxels())
    {
      cout<<"Give mask on the same grid as images!"<<endl;
      exit(1);
    }   
    else
      m=*mask;
  }
  else
  {
    m=image;
    m=1;
  }
  
  irtkRealPixel *p,*pm;
  p = image.GetPointerToVoxels();
  pm = m.GetPointerToVoxels();
  
  double sum=0, num=0,vol;
  
  for(int i=0;i<image.GetNumberOfVoxels();i++)
  {
    if((*p>padding)&&(*pm==1))
    {
      sum += *p;
      num++;
    }
    p++;
    pm++;
  }
  
  irtkImageAttributes attr = image.GetImageAttributes();
  double voxelvol = attr._dx*attr._dy*attr._dz;
  vol = sum*voxelvol;
  cout<<"vol = "<<vol<<endl;

  if (file_name !=NULL)
  {
    ofstream fileOut(file_name, ofstream::out | ofstream::app);

    if(!fileOut)
    {
      cerr << "Can't open file " << file_name << endl;
      exit(1);
    }

    fileOut.precision(6);  
    
    fileOut<<vol<<endl;
  }
  
  //The end of main()
}  
