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

#include <irtkImage.h>
#include <irtkMultiChannelImage.h>

char *output_name;
int label = 0;
irtkRealImage image;


void usage()
{
  cerr << "Usage: extract_label [image] [output] [label]"<<endl;
  exit(1);
}

int main(int argc, char **argv)
{

  if (argc < 4) {
    usage();
    exit(1);
  }

  image.Read(argv[1]);
  argc--;
  argv++;

  output_name=argv[1];
  argc--;
  argv++;

  label=atoi(argv[1]);
  argc--;
  argv++;

   irtkMultiChannelImage mch;
   mch.AddImage(image);
   image = mch.ExtractLabel(0,label);
   image.Write(output_name);
   //mch.Write(0,output_name);


}



