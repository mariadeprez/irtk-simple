/*=========================================================================

  Library   : project
  Module    : $RCSfile: em_classificationMultiCh.cc,v $
  Authors   : (C)opyright Maria Murgasova 2004
              See COPYRIGHT statement in top level directory.
  Purpose   :
  Date      : $Date: 2004/01/09 13:00:00 $

=========================================================================*/

#include <irtkImage.h>
#include <irtkMultiChannelImage.h>

char  *output_name = NULL;


void usage()
{
  cerr << "multiply_images [num of images] [image 1] ... [image n] [out]" << endl;
  cerr << "Options:" <<endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int i;
  int num_of_images;
  irtkMultiChannelImage mchImage;

  if (argc < 4)
  {
    usage();
  }

  num_of_images = atoi(argv[1]);
  argc--;
  argv++;

  for (i=0; i < num_of_images; i++)
  {
    irtkRealImage image;
    image.Read(argv[1]);
    mchImage.AddImage(image);
    cerr << "Image " << i <<" = " << argv[1] <<endl;
    argc--;
    argv++;
  }

  mchImage.SetPadding(0);
  output_name = argv[1];
  cerr << "Output = " << argv[1] <<endl;
  argc--;
  argv++;

  irtkRealImage sum = mchImage.Multiply();
  cerr << "Writing image" << endl;
  sum.Write(output_name);

  return 0;
}

