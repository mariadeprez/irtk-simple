/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id: enlarge_image.cc 335 2011-05-18 08:57:55Z ws207 $
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2009 onwards
  Date      : $Date: 2011-05-18 09:57:55 +0100 (Wed, 18 May 2011) $
  Version   : $Revision: 335 $
  Changes   : $Author: ws207 $

=========================================================================*/

#include <irtkImage.h>

char *output_name;
irtkRealImage image;


void usage()
{
  cerr << "Usage: split-image [image] [num packages] "<<endl;
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

  int packages=atoi(argv[1]);
  argc--;
  argv++;



  irtkImageAttributes attr = image.GetImageAttributes();
  
  //slices in package
  int pkg_z = attr._z/packages;
  double pkg_dz = attr._dz*packages;
  cout<<"packages: "<<packages<<"; slices: "<<attr._z<<"; slices in package: "<<pkg_z<<endl;
  cout<<"slice thickness "<<attr._dz<<"; slickess thickness in package: "<<pkg_dz<<endl;
  
  attr._z = pkg_z;
  attr._dz = pkg_dz;
  
  //fill values in each stack
  irtkRealImage stack(attr);
  char buffer[256];
  int i,j,k,t,l;
  double x,y,z,sx,sy,sz,ox,oy,oz;
  stack.GetOrigin(ox,oy,oz);
  for(l=0;l<packages;l++)
  {
    cout<<"Stack "<<l<<":"<<endl;
    for(t=0; t<stack.GetT();t++)
      for(k=0; k<stack.GetZ();k++)
        for(j=0; j<stack.GetY();j++)
	  for(i=0; i<stack.GetX();i++)
	    stack.Put(i,j,k,t,image(i,j,k*packages+l,t));
	
     //adjust origin
     
     //original image coordinates
     x=0;y=0;z=l;
     image.ImageToWorld(x,y,z);
     cout<<"image: "<<x<<" "<<y<<" "<<z<<endl;
     //stack coordinates
     sx=0;sy=0;sz=0;
     stack.PutOrigin(ox,oy,oz); //adjust to original value
     stack.ImageToWorld(sx,sy,sz);
     cout<<"stack: "<<sx<<" "<<sy<<" "<<sz<<endl;
     //adjust origin
     cout<<"adjustment needed: "<<x-sx<<" "<<y-sy<<" "<<z-sz<<endl;
     stack.PutOrigin(ox + (x-sx), oy + (y-sy), oz + (z-sz));
     sx=0;sy=0;sz=0;
     stack.ImageToWorld(sx,sy,sz);
     cout<<"adjusted: "<<sx<<" "<<sy<<" "<<sz<<endl;
	
     sprintf(buffer,"stack%i.nii.gz",l);
     stack.Write(buffer);
  }
  cout<<"done.";
}



