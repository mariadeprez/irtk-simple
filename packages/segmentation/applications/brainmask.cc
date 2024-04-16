
#include <irtkImage.h>
#include <irtkMultiChannelImage.h>

irtkRealImage image, mask;
int padding = 0;
char * output_name = NULL;

void usage()
{
  cerr << "Usage: brainmask [image] [mask] [output] <-padding> " << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int ok;

  if (argc < 3) {
    usage();
    exit(1);
  }

  // Input image
  image.Read(argv[1]);
  argc--;
  argv++;

  //reference image
  mask.Read(argv[1]);
  argc--;
  argv++;

  //reference image
  output_name = argv[1];
  argc--;
  argv++;

  // Parse remaining parameters
  while (argc > 1){
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-padding") == 0)){
      argc--;
      argv++;
      padding = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
   if (ok == false){
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  } 

   irtkMultiChannelImage mch;
   mch.SetPadding(padding);
   mch.AddImage(image);
   mch.SetMask(mask);
   mch.Brainmask();
   mch.Write(0,output_name);

}

