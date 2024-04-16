/*
 * cut-region.cc
 *
 *  Created on: 27-Feb-2009
 *      Author: mm3
 */

#include <irtkImage.h>
#include <irtkMeanShift.h>

char *input_name = NULL,*output_name = NULL;
int label = 1;
bool ok;

void usage()
{
  cerr << "Usage: lcc-queue [image] [output] <-label label>" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  irtkGreyImage image;
  
  if (argc < 3){
    usage();
  }

  input_name  = argv[1];
  argc--;
  argv++;
  output_name = argv[1];
  argc--;
  argv++;

  while (argc > 1){

    ok = false;

    if ((ok == false) && (strcmp(argv[1], "-label") == 0))
    {
      argc--;
      argv++;

      label = atoi(argv[1]);
      cout << "label = " << label << endl;

      argc--;
      argv++;
      ok = true;
    }

    if (ok == false){
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }



  // Read input
  image.Read(input_name);

  irtkMeanShift msh(image);
  msh.SetOutput(&image);
  msh.Lcc(label);

  image.Write(output_name);

  return 0;
}
