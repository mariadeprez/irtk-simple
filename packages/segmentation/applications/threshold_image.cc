/*=========================================================================

  Library   : project
  Module    : $RCSfile: treshold_image.cc$
  Authors   : (C)opyright Maria Murgasova 2004
              See COPYRIGHT statement in top level directory.
  Purpose   :
  Date      : $Date: 2004/01/09 13:00:00 $

=========================================================================*/

#include <irtkImage.h>

char *input_name = NULL, *output_name = NULL;
char *brainmask_name = NULL;
double treshold1 = 0, treshold2 = 0, tp=-1, value=-1;
bool ok, treshold1_given = false, treshold2_given = false, tp_given=false, value_given=false;

void usage()
{
  cerr << "Usage: treshold_image [in] [out]  <-threshold1> <-threshold2> <-padding> <-value>" << endl;
  exit(1);
}

int main(int argc, char **argv)
{

  if ((argc < 3)||(argc > 11))
  {
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


    if ((ok == false) && (strcmp(argv[1], "-threshold1") == 0))
    {
      argc--;
      argv++;

      treshold1 = atof(argv[1]);
      treshold1_given = true;
      cout << "threshold1 = " << treshold1 << endl;

      argc--;
      argv++;
      ok = true;
    }

    if ((ok == false) && (strcmp(argv[1], "-threshold2") == 0))
    {
      argc--;
      argv++;

      treshold2 = atof(argv[1]);
      treshold2_given = true;
      cout << "threshold2 = " << treshold2 << endl;

      argc--;
      argv++;
      ok = true;
    }


    if ((ok == false) && (strcmp(argv[1], "-padding") == 0))
    {
      argc--;
      argv++;

      tp = atof(argv[1]);
      tp_given = true;
      cout << "padding = " << tp << endl;

      argc--;
      argv++;
      ok = true;
    }

    if ((ok == false) && (strcmp(argv[1], "-value") == 0))
    {
      argc--;
      argv++;

      value = atof(argv[1]);
      value_given = true;
      cout << "value = " << value << endl;

      argc--;
      argv++;
      ok = true;
    }

    if (ok == false){
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }


  // Read input and probability maps
  irtkRealImage input;

  input.Read(input_name);

  int n = input.GetNumberOfVoxels();

  irtkRealPixel *ptr_input;
  int i;
  double temp;
  short min, max;

  if(!value_given) value=tp;

  ptr_input = input.GetPointerToVoxels();
  for (i=0; i<n; i++)
  {
    if ((*ptr_input > tp) || (!tp_given))
    {
      if (treshold1_given)
      {
        if (*ptr_input < treshold1) *ptr_input = value;
      }
      if (treshold2_given)
      {
        if (*ptr_input > treshold2) *ptr_input = value;
      }
    }
    ptr_input++;
  }

  input.Write(output_name);

  return 0;
}
