/*
 * cut-region.cc
 *
 *  Created on: 27-Feb-2009
 *      Author: mm3
 */

#include <irtkImage.h>

char *input_name = NULL, *mask_name = NULL, *output_name = NULL;
bool ok;

void usage()
{
  cerr << "Usage: cut-region [image] [mask] [cut-image] <-xy>" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  irtkRealImage in, mask, out;
  int x1, x2, y1, y2, z1, z2;
  int treshold=0;
  bool cropZ = true;

  if (argc < 4){
    usage();
  }

  input_name  = argv[1];
  argc--;
  argv++;
  mask_name  = argv[1];
  argc--;
  argv++;
  output_name = argv[1];
  argc--;
  argv++;
  
  // Parse options.
  while (argc > 1){
    ok = false;
    //Debug mode
    if ((ok == false) && (strcmp(argv[1], "-xy") == 0)){
      argc--;
      argv++;
      cropZ=false;
      ok = true;
    }

    if (ok == false){
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }    

  // Read input
  in.Read(input_name);
  mask.Read(mask_name);

  // Default roi
  x1 = 0;
  y1 = 0;
  z1 = 0;
  x2 = in.GetX();
  y2 = in.GetY();
  z2 = in.GetZ();

  int i,j,k;
/////////z coordinate
  int sum=0;
  for(k=in.GetZ()-1;k>=0;k--)
  {
    sum=0;
    for(j=in.GetY()-1;j>=0;j--)
      for(i=in.GetX()-1;i>=0;i--)
      {
        if(mask.Get(i,j,k)>0) sum++;
      }
    if (sum>treshold) break;
  }

  z2=k;
  cerr<<z2<<endl;

  sum=0;
  for(k=0; k<=in.GetZ()-1;k++)
  {
    sum=0;
    for(j=in.GetY()-1;j>=0;j--)
      for(i=in.GetX()-1;i>=0;i--)
      {
        if(mask.Get(i,j,k)>0) sum++;
      }
    if (sum>treshold) break;
  }

  z1=k;
  cerr<<z1<<endl;

/////////y coordinate
  sum=0;
  for(j=in.GetY()-1;j>=0;j--)
  {
    sum=0;
    for(k=in.GetZ()-1;k>=0;k--)
      for(i=in.GetX()-1;i>=0;i--)
      {
        if(mask.Get(i,j,k)>0) sum++;
      }
    if (sum>treshold) break;
  }

  y2=j;
  cerr<<y2<<endl;

  sum=0;
  for(j=0; j<=in.GetY()-1;j++)
  {
    sum=0;
    for(k=in.GetZ()-1;k>=0;k--)
      for(i=in.GetX()-1;i>=0;i--)
      {
        if(mask.Get(i,j,k)>0) sum++;
      }
    if (sum>treshold) break;
  }

  y1=j;
  cerr<<y1<<endl;

/////////x coordinate
  sum=0;
  for(i=in.GetX()-1;i>=0;i--)
  {
    sum=0;
    for(k=in.GetZ()-1;k>=0;k--)
      for(j=in.GetY()-1;j>=0;j--)
      {
        if(mask.Get(i,j,k)>0) sum++;
      }
    if (sum>treshold) break;
  }

  x2=i;
  cerr<<x2<<endl;

  sum=0;
  for(i=0; i<=in.GetX()-1;i++)
  {
    sum=0;
    for(k=in.GetZ()-1;k>=0;k--)
      for(j=in.GetY()-1;j>=0;j--)
      {
        if(mask.Get(i,j,k)>0) sum++;
      }
    if (sum>treshold) break;
  }

  x1=i;
  cerr<<x1<<endl;

  cerr<<x1<<" "<<y1<<" "<<z1<<" "<<x2<<" "<<y2<<" "<<z2<<endl;

  if(cropZ)
    in = in.GetRegion(x1, y1, z1, 0, x2+1, y2+1, z2+1, in.GetT());
  else
    in = in.GetRegion(x1, y1, 0, 0, x2+1, y2+1, in.GetZ(), in.GetT());
  in.Write(output_name);


  return 0;
}
