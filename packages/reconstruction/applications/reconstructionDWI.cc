/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2011 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/
#include <vector>
#include <string>
#include <irtkImage.h>
#include <irtkTransformation.h>
#include <irtkReconstruction.h>
#include <irtkReconstructionDTI.h>
#include <irtkSphericalHarmonics.h>

using namespace std;

//Application to perform reconstruction of volumetric MRI from thick slices.

void usage()
{
  cerr << "Usage: reconstructionDWI [reconstructed] [DWI] [bvalues] [directions] <options>\n" << endl;
  cerr << endl;

  cerr << "\t[reconstructed]         Name for the reconstructed volume. Nifti or Analyze format." << endl;
  cerr << "\t[N]                     Number of stacks." << endl;
  cerr << "\t[stack_1] .. [stack_N]  The input stacks. Nifti or Analyze format." << endl;
  cerr << "\t[bvalues]               List of b-values for each image and direction in a text file." << endl;
  cerr << "\t[directions]            List of directions for each image and direction in a text file." << endl;
  cerr << "\t" << endl;
  cerr << "Options:" << endl;
  cerr << "\t-dofin [dof_1]   .. [dof_N]    The transformations of the input stack to template" << endl;
  cerr << "\t                        in \'dof\' format used in IRTK." <<endl;
  cerr << "\t                        Only rough alignment with correct orienation and " << endl;
  cerr << "\t                        some overlap is needed." << endl;
  cerr << "\t                        Use \'id\' for an identity transformation for at least" << endl;
  cerr << "\t                        one stack. The first stack with \'id\' transformation" << endl;
  cerr << "\t                        will be resampled as template." << endl;
  cerr << "\t-thickness [th_1] .. [th_N] Give slice thickness.[Default: twice voxel size in z direction]"<<endl;
  cerr << "\t-mask [mask]              Binary mask to define the region od interest. [Default: whole image]"<<endl;
  cerr << "\t-packages [num]           Give number of packages used during acquisition for each stack."<<endl;
  cerr << "\t                          The stacks will be split into packages during registration iteration 1"<<endl;
  cerr << "\t                          and then into odd and even slices within each package during "<<endl;
  cerr << "\t                          registration iteration 2. The method will then continue with slice to"<<endl;
  cerr << "\t                          volume approach. [Default: slice to volume registration only]"<<endl;
  cerr << "\t-iterations [iter]        Number of registration-reconstruction iterations. [Default: 9]"<<endl;
  cerr << "\t-sigma [sigma]            Stdev for bias field. [Default: 12mm]"<<endl;
  cerr << "\t-resolution [res]         Isotropic resolution of the volume. [Default: 0.75mm]"<<endl;
  cerr << "\t-multires [levels]        Multiresolution smooting with given number of levels. [Default: 3]"<<endl;
  cerr << "\t-average [average]        Average intensity value for stacks [Default: 700]"<<endl;
  cerr << "\t-delta [delta]            Parameter to define what is an edge. [Default: 150]"<<endl;
  cerr << "\t-lambda [lambda]          Smoothing parameter. [Default: 0.02]"<<endl;
  cerr << "\t-lambdaLB [lambdaLB]          Smoothing parameter. [Default: 0.05]"<<endl;
  cerr << "\t-regul_steps [steps]      Number of regularisation steps. [Default: 1]"<<endl;
  cerr << "\t-lastIter [lambda]        Smoothing parameter for last iteration. [Default: 0.01]"<<endl;
  cerr << "\t-smooth_mask [sigma]      Smooth the mask to reduce artefacts of manual segmentation. [Default: 4mm]"<<endl;
  cerr << "\t-global_bias_correction   Correct the bias in reconstructed image against previous estimation."<<endl;
  cerr << "\t-low_intensity_cutoff     Lower intensity threshold for inclusion of voxels in global bias correction."<<endl;
  cerr << "\t-remove_black_background  Create mask from black background."<<endl;
  cerr << "\t-transformations [folder] Use existing slice-to-volume transformations to initialize the reconstruction."<<endl;
  cerr << "\t-force_exclude [number of slices] [ind1] ... [indN]  Force exclusion of slices with these indices."<<endl;
  cerr << "\t-force_exclude_stacks [number of stacks] [ind1] ... [indN]  Force exclusion of whole stacks with these indices."<<endl;
  cerr << "\t-no_intensity_matching    Switch off intensity matching."<<endl;
  cerr << "\t-no_robust_statistics     Switch off robust statistics."<<endl;
  cerr << "\t-no_robust_statistics_sh  Switch off sh robust statistics."<<endl;
  cerr << "\t-exclude_slices_only      Do not exclude individual voxels."<<endl;
  cerr << "\t-bspline                  Use multi-level bspline interpolation instead of super-resolution."<<endl;
  cerr << "\t-no_motion                Switch off volumetric registration."<<endl;
  cerr << "\t-log_prefix [prefix]      Prefix for the log file."<<endl;
  cerr << "\t-info [filename]          Filename for slice information in\
                                       tab-sparated columns."<<endl;
  cerr << "\t-debug                    Debug mode - save intermediate results."<<endl;
  cerr << "\t-save_transformations     Save slice-wise transformations."<<endl;
  cerr << "\t-no_log                   Do not redirect cout and cerr to log files."<<endl;
  cerr << "\t" << endl;
  cerr << "\t" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  
  //utility variables
  int i, ok;
  char buffer[256];
  irtkRealImage stack; 
  
  //declare variables for input
  /// Name for output volume
  char * output_name = NULL;
  /// Slice stacks
  vector<irtkRealImage> stacks;
  vector<irtkRealImage> corrected_stacks;
  vector<string> stack_files;
  /// Stack transformation
  vector<irtkRigidTransformation> stack_transformations;
  /// user defined transformations
  bool have_stack_transformations = false;
  /// Stack thickness
  vector<double > thickness;
  ///number of stacks
  int nStacks;
  /// number of packages for each stack
  vector<int> packages;

  
   irtkRigidTransformation orient;
  
  // Default values.
  int templateNumber=-1;
  irtkRealImage *mask=NULL;
  int iterations = 3;
  bool debug = false;
  bool save_transformations = false;
  double sigma=20;
  double resolution = 0;
  double lambda = 0.02;
  double delta = 0.001;
  int regul_steps = 0;
  int levels = 3;
  double lastIterLambda = 0.01;
  double regDTI = 0;
  int rec_iterations;
  double averageValue = 700;
  double smooth_mask = 0;
  bool global_bias_correction = false;
  double low_intensity_cutoff = 0.01;
  //folder for slice-to-volume registrations, if given
  char * folder=NULL;
  //flag to remove black background, e.g. when neonatal motion correction is performed
  bool remove_black_background = false;
  //flag to swich the intensity matching on and off
  bool intensity_matching = true;
  bool intensity_matching_sh = true;
  bool rescale_stacks = false;

  //flag to swich the robust statistics on and off
  bool robust_statistics = true;
  bool robust_statistics_sh = true;
  bool robust_slices_only = false;
  //flag to replace super-resolution reconstruction by multilevel B-spline interpolation
  bool bspline = false;
  bool recon_1D = false;
  bool recon_interpolate = false;
  bool sh_only=false;
  bool _reset_robust_statistics = false;
  
  int int_matching_first_iter = 0;
  int sr_sh_iterations = 10;
  double sh_alpha = 5;
  
  bool motion_model_hs = true;
  bool motion = true;
  
  double motion_sigma = 15;
  
  int order = 4;
  double lambdaLB = 0;//0.05;
  
  irtkRealImage average;

  string info_filename = "slice_info.tsv";
  string log_id;
  bool no_log = false;

  //forced exclusion of slices
  int number_of_force_excluded_slices = 0;
  vector<int> force_excluded;

  //forced exclusion of stacks
  int number_of_force_excluded_stacks = 0;
  vector<int> force_excluded_stacks;

  //Create reconstruction object
  irtkReconstructionDTI reconstruction;
    
  //if not enough arguments print help
  if (argc < 3)
    usage();
  
  //read output name
  output_name = argv[1];
  argc--;
  argv++;
  cout<<"Recontructed DWI volume name ... "<<output_name<<endl;

  //read 4D image
  irtkRealImage image4D;
  cout<<"Reading stack ... "<<argv[1]<<endl;
  image4D.Read(argv[1]);
  argc--;
  argv++;
  
  irtkImageAttributes attr = image4D.GetImageAttributes();
  nStacks = attr._t;
  cout<<"Number 0f stacks ... "<<nStacks<<endl;

  // Read stacks 
  for (i=0;i<nStacks;i++)
  {
    //cout<<"Stack "<<i<<endl; cout.flush();
    stacks.push_back(image4D.GetRegion(0,0,0,i,attr._x, attr._y,attr._z,i+1));
    corrected_stacks.push_back(stacks[i]);
    if (debug)
    {
      sprintf(buffer,"stack%i.nii.gz",i);
      stacks[i].Write(buffer);
    }
    
  }
  
  //Read the name of the textfile with directions
  char *textfile = argv[1];
  argc--;
  argv++;
  cout<<"Directions textfile is "<<textfile<<endl;
  
  //read target
  irtkRealImage target(argv[1]);
  cout<<"Target is "<<argv[1]<<endl;
  argc--;
  argv++;

  //read initial orientation
  orient.irtkTransformation::Read(argv[1]);
  cout<<"Initial transformation is "<<argv[1]<<endl;
  argc--;
  argv++;
  if (debug)
    orient.Print();


 
  
  // Parse options.
  while (argc > 1){
    ok = false;
    
    //Read stack transformations
    if ((ok == false) && (strcmp(argv[1], "-dofin") == 0)){
      argc--;
      argv++;
      
      for (i=0;i<nStacks;i++)
      {
        irtkTransformation *transformation;
        cout<<"Reading transformation ... "<<argv[1]<<" ... ";
        cout.flush();
        if (strcmp(argv[1], "id") == 0)
        {
          transformation = new irtkRigidTransformation;
          if ( templateNumber < 0) templateNumber = i;
        }
        else
        {
          transformation = irtkTransformation::New(argv[1]);
        }
        cout<<" done."<<endl;

        argc--;
        argv++;
        irtkRigidTransformation *rigidTransf = dynamic_cast<irtkRigidTransformation*> (transformation);
        stack_transformations.push_back(*rigidTransf);
        delete rigidTransf;
      }
      reconstruction.InvertStackTransformations(stack_transformations);
      have_stack_transformations = true;
    }

    //Read slice thickness
    if ((ok == false) && (strcmp(argv[1], "-thickness") == 0)){
      argc--;
      argv++;
      cout<< "Slice thickness is ";
      double thick = atof(argv[1]);
      argc--;
      argv++;
      for (i=0;i<nStacks;i++)
      {
        thickness.push_back(thick);
	//cout<<thickness[i]<<" ";
       }
       cout<<"."<<endl;
      ok = true;
    }
    
    //Read number of packages for each stack
    if ((ok == false) && (strcmp(argv[1], "-packages") == 0)){
      argc--;
      argv++;
      int pnum = atoi(argv[1]);
      argc--;
      argv++;
      cout<< "Package number is "<<pnum<<endl;
      for (i=0;i<nStacks;i++)
      {
        packages.push_back(pnum);
	//cout<<packages[i]<<" ";        
      }
      //cout<<"size="<<packages.size()<<endl;
      cout<<"."<<endl;
      ok = true;
    }

    //Read binary mask for final volume
    if ((ok == false) && (strcmp(argv[1], "-mask") == 0)){
      argc--;
      argv++;
      mask= new irtkRealImage(argv[1]);
      ok = true;
      argc--;
      argv++;
    }
    
    //Read number of registration-reconstruction iterations
    if ((ok == false) && (strcmp(argv[1], "-iterations") == 0)){
      argc--;
      argv++;
      iterations=atoi(argv[1]);
      ok = true;
      argc--;
      argv++;
    }

    //Read number of registration-reconstruction iterations
    if ((ok == false) && (strcmp(argv[1], "-sr_sh_iterations") == 0)){
      argc--;
      argv++;
      sr_sh_iterations=atoi(argv[1]);
      ok = true;
      argc--;
      argv++;
    }
    
    //Read number of registration-reconstruction iterations
    if ((ok == false) && (strcmp(argv[1], "-intensity_matching_first_iter") == 0)){
      argc--;
      argv++;
      int_matching_first_iter=atoi(argv[1]);
      ok = true;
      argc--;
      argv++;
    }


    //Step for SR SH
    if ((ok == false) && (strcmp(argv[1], "-sh_alpha") == 0)){
      argc--;
      argv++;
      sh_alpha=atof(argv[1]);
      ok = true;
      argc--;
      argv++;
    } 

    //Variance of Gaussian kernel to smooth the bias field.
    if ((ok == false) && (strcmp(argv[1], "-sigma") == 0)){
      argc--;
      argv++;
      sigma=atof(argv[1]);
      ok = true;
      argc--;
      argv++;
    } 
    
    //Variance of Gaussian kernel to smooth the bias field.
    if ((ok == false) && (strcmp(argv[1], "-motion_sigma") == 0)){
      argc--;
      argv++;
      motion_sigma=atof(argv[1]);
      ok = true;
      argc--;
      argv++;
    } 
    
    //Smoothing parameter
    if ((ok == false) && (strcmp(argv[1], "-lambda") == 0)){
      argc--;
      argv++;
      lambda=atof(argv[1]);
      ok = true;
      argc--;
      argv++;
    }
    
    //Smoothing parameter
    if ((ok == false) && (strcmp(argv[1], "-regDTI") == 0)){
      argc--;
      argv++;
      regDTI=atof(argv[1]);
      ok = true;
      argc--;
      argv++;
    }
    
    //Smoothing parameter
    if ((ok == false) && (strcmp(argv[1], "-lambdaLB") == 0)){
      argc--;
      argv++;
      lambdaLB=atof(argv[1]);
      ok = true;
      argc--;
      argv++;
    }
    
    //Smoothing parameter for last iteration
    if ((ok == false) && (strcmp(argv[1], "-lastIter") == 0)){
      argc--;
      argv++;
      lastIterLambda=atof(argv[1]);
      ok = true;
      argc--;
      argv++;
    }
    
    //Parameter to define what is an edge
    if ((ok == false) && (strcmp(argv[1], "-delta") == 0)){
      argc--;
      argv++;
      delta=atof(argv[1]);
      ok = true;
      argc--;
      argv++;
    }
    
    //Parameter to define what is an edge
    if ((ok == false) && (strcmp(argv[1], "-regul_steps") == 0)){
      argc--;
      argv++;
      regul_steps=atoi(argv[1]);
      ok = true;
      argc--;
      argv++;
    }
    
    //Isotropic resolution for the reconstructed volume
    if ((ok == false) && (strcmp(argv[1], "-resolution") == 0)){
      argc--;
      argv++;
      resolution=atof(argv[1]);
      ok = true;
      argc--;
      argv++;
    }

    //Number of resolution levels
    if ((ok == false) && (strcmp(argv[1], "-multires") == 0)){
      argc--;
      argv++;
      levels=atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }

    //Smooth mask to remove effects of manual segmentation
    if ((ok == false) && (strcmp(argv[1], "-smooth_mask") == 0)){
      argc--;
      argv++;
      smooth_mask=atof(argv[1]);
      argc--;
      argv++;
      ok = true;
    }

    //Switch off intensity matching
    if ((ok == false) && (strcmp(argv[1], "-motion_model_hs") == 0)){
      argc--;
      argv++;
      motion_model_hs=true;
      ok = true;
    }
    
    //Switch off intensity matching
    if ((ok == false) && (strcmp(argv[1], "-no_intensity_matching") == 0)){
      argc--;
      argv++;
      intensity_matching=false;
      intensity_matching_sh=false;
      ok = true;
    }
    
    //Switch off intensity matching
    if ((ok == false) && (strcmp(argv[1], "-_reset_robust_statistics") == 0)){
      argc--;
      argv++;
      _reset_robust_statistics=true;
      ok = true;
    }
    
    //Switch off intensity matching
    if ((ok == false) && (strcmp(argv[1], "-no_vol_intensity_matching") == 0)){
      argc--;
      argv++;
      intensity_matching=false;
      ok = true;
    }
    
    //Switch off intensity matching
    if ((ok == false) && (strcmp(argv[1], "-no_sh_intensity_matching") == 0)){
      argc--;
      argv++;
      intensity_matching_sh=false;
      ok = true;
    }
    
    //Switch off robust statistics
    if ((ok == false) && (strcmp(argv[1], "-no_robust_statistics") == 0)){
      argc--;
      argv++;
      robust_statistics=false;
      robust_statistics_sh=false;
      ok = true;
    }
    
    //Switch off robust statistics sh
    if ((ok == false) && (strcmp(argv[1], "-no_robust_statistics_sh") == 0)){
      argc--;
      argv++;
      robust_statistics_sh=false;
      ok = true;
    }
    
    //Switch off robust statistics sh
    if ((ok == false) && (strcmp(argv[1], "-no_motion") == 0)){
      argc--;
      argv++;
      motion=false;
      ok = true;
    }
    
    //Switch off robust statistics
    if ((ok == false) && (strcmp(argv[1], "-exclude_slices_only") == 0)){
      argc--;
      argv++;
      robust_slices_only=true;
      ok = true;
    }
    
    //Use multilevel B-spline interpolation instead of super-resolution
    if ((ok == false) && (strcmp(argv[1], "-bspline") == 0)){
      argc--;
      argv++;
      bspline=true;
      ok = true;
    }

    //Perform bias correction of the reconstructed image agains the GW image in the same motion correction iteration
    if ((ok == false) && (strcmp(argv[1], "-global_bias_correction") == 0)){
      argc--;
      argv++;
      global_bias_correction=true;
      ok = true;
    }

    if ((ok == false) && (strcmp(argv[1], "-low_intensity_cutoff") == 0)){
      argc--;
      argv++;
      low_intensity_cutoff=atof(argv[1]);
      argc--;
      argv++;
      ok = true;
    }

    //Debug mode
    if ((ok == false) && (strcmp(argv[1], "-debug") == 0)){
      argc--;
      argv++;
      debug=true;
      ok = true;
    }
     //Save transformations
    if ((ok == false) && (strcmp(argv[1], "-save_transformations") == 0)){
      argc--;
      argv++;
      save_transformations=true;
      ok = true;
    }
    
    //Prefix for log files
    if ((ok == false) && (strcmp(argv[1], "-log_prefix") == 0)){
      argc--;
      argv++;
      log_id=argv[1];
      ok = true;
      argc--;
      argv++;
    }

    //No log files
    if ((ok == false) && (strcmp(argv[1], "-no_log") == 0)){
      argc--;
      argv++;
      no_log=true;
      ok = true;
    }

    // rescale stacks to avoid error:
    // irtkImageRigidRegistrationWithPadding::Initialize: Dynamic range of source is too large
    if ((ok == false) && (strcmp(argv[1], "-rescale_stacks") == 0)){
      argc--;
      argv++;
      rescale_stacks=true;
      ok = true;
    }       

    // Save slice info
    if ((ok == false) && (strcmp(argv[1], "-info") == 0)) {
        argc--;
        argv++;
        info_filename=argv[1];
        ok = true;
        argc--;
        argv++;
    }
 
    //Read transformations from this folder
    if ((ok == false) && (strcmp(argv[1], "-transformations") == 0)){
      argc--;
      argv++;
      folder=argv[1];
      ok = true;
      argc--;
      argv++;
    }

    //Remove black background
    if ((ok == false) && (strcmp(argv[1], "-remove_black_background") == 0)){
      argc--;
      argv++;
      remove_black_background=true;
      ok = true;
    }

    //Force removal of certain slices
    if ((ok == false) && (strcmp(argv[1], "-force_exclude") == 0)){
      argc--;
      argv++;
      number_of_force_excluded_slices = atoi(argv[1]);
      argc--;
      argv++;

      cout<< number_of_force_excluded_slices<< " force excluded slices: ";
      for (i=0;i<number_of_force_excluded_slices;i++)
      {
        force_excluded.push_back(atoi(argv[1]));
	cout<<force_excluded[i]<<" ";
        argc--;
        argv++;
       }
       cout<<"."<<endl;

      ok = true;
    }


    //Force removal of certain stacks
    if ((ok == false) && (strcmp(argv[1], "-force_exclude_stacks") == 0)){
      argc--;
      argv++;
      number_of_force_excluded_stacks = atoi(argv[1]);
      argc--;
      argv++;

      cout<< number_of_force_excluded_stacks<< " force excluded stacks: ";
      for (i=0;i<number_of_force_excluded_stacks;i++)
      {
        force_excluded_stacks.push_back(atoi(argv[1]));
	cout<<force_excluded_stacks[i]<<" ";
        argc--;
        argv++;
       }
       cout<<"."<<endl;

      ok = true;
    }
    
    //Perform reconstruction only in z direction
    if ((ok == false) && (strcmp(argv[1], "-1D") == 0)){
      argc--;
      argv++;
      recon_1D = true;
      ok = true;
    }

    //Perform interpolation like reconstruction
    if ((ok == false) && (strcmp(argv[1], "-interpolate") == 0)){
      argc--;
      argv++;
      recon_interpolate = true;
      ok = true;
    }

    //SH order
    if ((ok == false) && (strcmp(argv[1], "-order") == 0)){
      argc--;
      argv++;
      order = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }




    if (ok == false){
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  if (rescale_stacks)
  {
      for (i=0;i<nStacks;i++)
          reconstruction.Rescale(stacks[i],1000);
  }
  
  //If transformations were not defined by user, set them to identity
  //if(!have_stack_transformations)
  //{
    stack_transformations.clear();
    for (i=0;i<nStacks;i++)
    {
      //irtkRigidTransformation *rigidTransf = new irtkRigidTransformation;
      stack_transformations.push_back(orient);
      //delete rigidTransf;
    }
    templateNumber = 0;  
    reconstruction.InvertStackTransformations(stack_transformations);
  //}

  //Initialise slice thickness to isotropic if not given by user
  if (thickness.size()==0)
  {
    cout<< "Slice thickness is ";
    for (i=0;i<nStacks;i++)
    {
      double dx,dy,dz;
      stacks[i].GetPixelSize(&dx,&dy,&dz);
      thickness.push_back(dz);
      //cout<<thickness[i]<<" ";
    }
    cout<<thickness[0];
    cout<<"."<<endl;
  }
  
    //this is the number of directions we expect at the bval and bvecs file, including b0
    cout<<"Number of diffusion weighted directions is "<<nStacks<<endl;
    vector<vector<double> > directions(3,vector<double>(nStacks+1,0));
    vector<double> bvalues(nStacks+1,0);

    ifstream in(textfile);
    int coord = 0;
    int dir = 0;
    double num;
    
    cout<<"Reading directions: "<<endl;
    if (in.is_open()) 
    {
      while (!in.eof()) 
      {
	in >> num;
	if ((coord<4)&&(dir<nStacks))
	{
	  if(coord<3)
	    directions[coord][dir+1]=num;
	  else
	    bvalues[dir+1]=num;
	}
	cout<<num<<" ";
	coord++;
	if (coord>=4)
	{
	  coord=0;
	  dir++;
	  cout<<endl;
	}
      }
      in.close();
    } 
    else 
    {
      cout << "Unable to open file " << textfile << endl;
      exit(1);
    }
    


  //Output volume
  irtkRealImage reconstructed;

  //Set debug mode
  if (debug) reconstruction.DebugOn();
  else reconstruction.DebugOff();
  
  //Set force excluded slices
  reconstruction.SetForceExcludedSlices(force_excluded);

  //Set low intensity cutoff for bias estimation
  reconstruction.SetLowIntensityCutoff(low_intensity_cutoff)  ;

  
  // Check whether the template stack can be indentified
  if (templateNumber<0)
  {
    cerr<<"Please identify the template by assigning id transformation."<<endl;
    exit(1);
  }  
  //If no mask was given and flag "remove_black_background" is false, try to create mask from the template image in case it was padded
  if ((mask==NULL)&&(!remove_black_background))
  {
    mask = new irtkRealImage(stacks[templateNumber]);
    *mask = reconstruction.CreateMask(*mask);
  }
  /*
  if (mask !=NULL)
  {
    //first resample the mask to the space of the stack
    //for template stact the transformation is identity
    irtkRealImage m = *mask;
    reconstruction.TransformMask(stacks[templateNumber],m,stack_transformations[templateNumber]);
    //Crop template stack
    reconstruction.CropImage(stacks[templateNumber],m);
    if (debug)
    {
      m.Write("maskTemplate.nii.gz"); 
      stacks[templateNumber].Write("croppedTemplate.nii.gz");
    }
  }
  */
  //Create template volume with isotropic resolution 
  //if resolution==0 it will be determined from in-plane resolution of the image
  
  //resolution = reconstruction.CreateTemplate(stacks[templateNumber],resolution);
  if (resolution==0)
  {
    double dx,dy,dz;
    stacks[0].GetPixelSize(&dx,&dy,&dz);
    resolution = dx;
  }
  
  resolution = reconstruction.CreateTemplate(target,resolution);
  
  
  //Set mask to reconstruction object. 
  reconstruction.SetMaskOrient(mask,orient);   
  //reconstruction.SetMask(mask,0);   
  reconstruction.SetRegulSteps(regul_steps);

  //to redirect output from screen to text files
  
  //to remember cout and cerr buffer
  streambuf* strm_buffer = cout.rdbuf();
  streambuf* strm_buffer_e = cerr.rdbuf();
  //files for registration output
  string name;
  name = log_id+"log-registration.txt";
  ofstream file(name.c_str());
  name = log_id+"log-registration-error.txt";
  ofstream file_e(name.c_str());
  //files for reconstruction output
  name = log_id+"log-reconstruction.txt";
  ofstream file2(name.c_str());
  name = log_id+"log-evaluation.txt";
  ofstream fileEv(name.c_str());
  name = log_id +"log-sh.txt";
  ofstream fileSH(name.c_str());
  
  name = log_id +"log-consistency.csv";
  ofstream fileConsistency(name.c_str());
  
  name = log_id +"log-consistency-ex.csv";
  ofstream fileConsistencyEx(name.c_str());
  
  //set precision
  cout<<setprecision(3);
  cerr<<setprecision(3);

  //perform volumetric registration of the stacks
  //redirect output to files
  if ( ! no_log ) {
      cerr.rdbuf(file_e.rdbuf());
      cout.rdbuf (file.rdbuf());
  }
  
  reconstruction.SetTemplateImage(stacks[0],orient);
  
  if (save_transformations)
      reconstruction.DebugOn();
  
  if (motion)
    reconstruction.StackRegistrations(stacks,stack_transformations);//,templateNumber);

  if (save_transformations)
      reconstruction.DebugOff();

  //if remove_black_background flag is set, create mask from black background of the stacks
  if (remove_black_background)
    reconstruction.CreateMaskFromBlackBackground(stacks,stack_transformations, smooth_mask);
  
  cout<<endl;
  //redirect output back to screen
  if ( ! no_log ) {
      cout.rdbuf (strm_buffer);
      cerr.rdbuf (strm_buffer_e);
  }
  

  average = reconstruction.CreateAverage(stacks,stack_transformations);
  if (debug)
    average.Write("average1.nii.gz");
  reconstruction.SetTemplate(average);
  reconstruction.MaskVolume();
  reconstructed = reconstruction.GetReconstructed();
  if (debug)
    reconstructed.Write("image0.nii.gz");
  //exit(1);

  //Create slices and slice-dependent transformations
  reconstruction.CreateSlicesAndTransformations(stacks,stack_transformations,thickness);
  reconstruction.CreateSliceDirections(directions,bvalues);
  //reconstruction.SaveSlices();
  reconstruction.SetForceExcludedStacks(force_excluded_stacks);
  //calculate slice order
  int z = stacks[0].GetZ();
  int j,k;
  vector<int> slice_order;
  cout<<"Slices in stack: "<<z<<endl;
  for(i=0;i<stacks.size();i++) //go though all directions
  {
      for(k=0;k<z;k=k+2) //even
        slice_order.push_back(i*z+k);
      for(k=1;k<z;k=k+2) //odd
        slice_order.push_back(i*z+k);
  }
    
  if (debug)
  {
    for(i=0;i<slice_order.size();i++)
      cout<<slice_order[i]<<" ";
    cout<<"Number of slices: "<<slice_order.size();
  }
  reconstruction.SetSliceOrder(slice_order);
  //exit(1);
  
  
  //Mask all the slices
  reconstruction.MaskSlices();
  
  //Set sigma for the bias field smoothing
  if (sigma>0)
    reconstruction.SetSigma(sigma);
    //TEST
    //reconstruction.SetSigma(15);
  //else
  //{
    //cerr<<"Please set sigma larger than zero. Current value: "<<sigma<<endl;
    //exit(1);
    //reconstruction.SetSigma(20);
  //}
  
  //Set global bias correction flag
  if (global_bias_correction)
    reconstruction.GlobalBiasCorrectionOn();
  else 
    reconstruction.GlobalBiasCorrectionOff();
    
  //if given read slice-to-volume registrations
  if (folder!=NULL)
    reconstruction.ReadTransformation(folder);
  

    
    
  //Initialise data structures for EM
  reconstruction.InitializeEM();
  
  //CHange: only packages and odd even for now
  //interleaved registration-reconstruction iterations
  //iterations=4;
  
  if(sh_only)
    reconstruction.SetSmoothingParameters(delta, lambda);
  else
    reconstruction.SetSmoothingParameters(delta, lambda*8);
  
  if(recon_1D)
  {
    reconstruction.Set1DRecon();
    if (debug)
        cout<<"Setting 1D"<<endl;
  }
  else
  {
    if (recon_interpolate)
    {
      reconstruction.SetInterpolationRecon();
      if (debug)
        cout<<"Setting interpolate"<<endl;
          
    }
    else
    {
      reconstruction.Set3DRecon();
      if (debug)
        cout<<"Setting 3D"<<endl;
    }

  }
  
  
  for (int iter=0;iter<iterations;iter++)
  {
    //only last iteration!!!
    if(sh_only)
      if(iter<iterations)
	continue;
    //Print iteration number on the screen
      if ( ! no_log ) {
          cout.rdbuf (strm_buffer);
      }
    cout<<"Iteration "<<iter<<". "<<endl;

    //perform slice-to-volume registrations - skip the first iteration 
    if ((iter>0)&&(!sh_only))
    {
        if ( ! no_log ) {
            cerr.rdbuf(file_e.rdbuf());
            cout.rdbuf (file.rdbuf());
        }
      cout<<"Iteration "<<iter<<": "<<endl;
      
      //if((packages.size()>0)&&( (iter<(iterations-1))||(iter<3) ))
      //if((packages.size()>0)&&(iter<=iterations*(levels-1)/levels)&&(iter<(iterations-1)))
      if((iter<=1) &&(packages.size()>0))
      {
	if(iter==1)
          reconstruction.PackageToVolume(corrected_stacks,packages,iter);
	else
	{
	  if(iter==2)
            reconstruction.PackageToVolume(corrected_stacks,packages,iter,true);
	  else
	  {
            if(iter==3)
	      reconstruction.PackageToVolume(corrected_stacks,packages,iter,true,true);
	    else
	    {
	      if(iter>=4)
                reconstruction.PackageToVolume(corrected_stacks,packages,iter,true,true,iter-2);
	      else
	        reconstruction.SliceToVolumeRegistration();
	    }
	  }
	}
      }
      else
      {
        reconstruction.SliceToVolumeRegistration();
        //reconstruction.PackageToVolume(stacks,packages,iter,true,true,8);
	
      }
      
      if (iter>1) 
      {
	if(motion_sigma>0)
	{
	  reconstruction.SetMotionSigma(motion_sigma);
	  if(motion_model_hs)
	    reconstruction.PostProcessMotionParametersHS2(stacks[templateNumber]);
	  else
	    reconstruction.PostProcessMotionParameters2(stacks[templateNumber]);
	}
      }
      
      if (debug)
      reconstruction.SaveTransformations();

      reconstruction.UpdateSlices(stacks,thickness);
      reconstruction.MaskSlices();

      cout<<endl;
      if ( ! no_log ) {
          cerr.rdbuf (strm_buffer_e);
      }
    }

    
    //Write to file
    if ( ! no_log ) {
        cout.rdbuf (file2.rdbuf());
    }
    cout<<endl<<endl<<"Iteration "<<iter<<": "<<endl<<endl;
    
    //Set smoothing parameters 
    //amount of smoothing (given by lambda) is decreased with improving alignment  
//     //delta (to determine edges) stays constant throughout
    if(iter==(iterations-1))
      reconstruction.SetSmoothingParameters(delta,lambda);//lastIterLambda);
    else
    {
      double l=lambda;
      for (i=0;i<levels;i++)
      {
        if (iter==iterations*(levels-i-1)/levels)
          reconstruction.SetSmoothingParameters(delta, l);
        l*=2;
      }
    }
    
    //Use faster reconstruction during iterations and slower for final reconstruction
    if ( iter<(iterations-1) )
      reconstruction.SpeedupOn();
    else 
      reconstruction.SpeedupOff();
    
    if(robust_slices_only)
      reconstruction.ExcludeWholeSlicesOnly();
    
    //Initialise values of weights, scales and bias fields
    reconstruction.InitializeEMValues();
    
    //Calculate matrix of transformation between voxels of slices and volume
    if (bspline)
      reconstruction.CoeffInitBSpline();
    else
    {
      //reconstruction.Set3DRecon();
      reconstruction.CoeffInit();
    }
    
    //Initialize reconstructed image with Gaussian weighted reconstruction
    if (bspline)
      reconstruction.BSplineReconstruction();
    else
      reconstruction.GaussianReconstruction(0);

    //Simulate slices (needs to be done after Gaussian reconstruction)
    reconstruction.SimulateSlices();
        
    //Initialize robust statistics parameters
    reconstruction.InitializeRobustStatistics();
    
    //EStep
    if(robust_statistics)
      reconstruction.EStep();

    //number of reconstruction iterations
    if ( iter==(iterations-1) ) 
    {
      if(robust_statistics)
        rec_iterations = 10;      
      else
        rec_iterations = 0;      
    }
    else 
      rec_iterations = 0;
    
    if ((bspline)&&(!robust_statistics)&&(!intensity_matching))
      rec_iterations=0;
        
    //reconstruction iterations
    if(iter==iterations)
      rec_iterations = 10;
    else
      rec_iterations = 10;

    i=0;
    for (i=0;i<rec_iterations;i++)
    {
      cout<<endl<<"  Reconstruction iteration "<<i<<". "<<endl;
      
      //if (i==15)
	//intensity_matching = true;
      
      if (intensity_matching)
      {
        //calculate bias fields
        if (sigma>0)
          reconstruction.Bias();
        //calculate scales
        reconstruction.Scale();
      }
      
      //Update reconstructed volume
      if (bspline)
        reconstruction.BSplineReconstruction();
      else
        reconstruction.Superresolution(i+1);
      
      
      /*
      if (intensity_matching)
      {
        //if((sigma>0)&&(!global_bias_correction))
          reconstruction.NormaliseBias(i);
      }
      */

      // Simulate slices (needs to be done
      // after the update of the reconstructed volume)
      reconstruction.SimulateSlices();
            
      if(robust_statistics)
        reconstruction.MStep(i+1);
      
      //E-step
      if(robust_statistics)
        reconstruction.EStep();
      
    //Save intermediate reconstructed image
    if (debug)
    {
      reconstructed=reconstruction.GetReconstructed();
      sprintf(buffer,"super%i.nii.gz",i);
      reconstructed.Write(buffer);
    }

      
    }//end of reconstruction iterations
    
    //Mask reconstructed image to ROI given by the mask
    if(!bspline)
      reconstruction.MaskVolume();

    //Save reconstructed image
    //if (debug)
    //{
      reconstructed=reconstruction.GetReconstructed();
      sprintf(buffer,"image%i.nii.gz",iter);
      reconstructed.Write(buffer);
      //reconstruction.SaveConfidenceMap();
    //}

   //Evaluate - write number of included/excluded/outside/zero slices in each iteration in the file
    if ( ! no_log ) {
        cout.rdbuf (fileEv.rdbuf());
    }
   reconstruction.Evaluate(iter);
   cout<<endl;

   if ( ! no_log ) {
       cout.rdbuf (strm_buffer);
   }
   
  }// end of interleaved registration-reconstruction iterations
      //exit(1);

  //save final result
  //reconstruction.RestoreSliceIntensities();
  //reconstruction.ScaleVolume();
  /*
  reconstructed=reconstruction.GetReconstructed();
  reconstructed.Write(output_name); 
  reconstruction.SaveTransformations();
  reconstruction.SaveSlices();

  if ( info_filename.length() > 0 )
      reconstruction.SlicesInfo( info_filename.c_str(),
                                 stack_files );
 
  if(debug)
  {
    reconstruction.SaveWeights();
    reconstruction.SaveBiasFields();
    //reconstruction.SaveConfidenceMap();
    reconstruction.SimulateStacks(stacks);
    for (unsigned int i=0;i<stacks.size();i++)
    {
      sprintf(buffer,"simulated%i.nii.gz",i);
      stacks[i].Write(buffer);
    }
  }
  
*/
  
  if (debug)
  {
  
    corrected_stacks.clear();
    for (unsigned int i=0;i<stacks.size();i++)
    {
      corrected_stacks.push_back(stacks[i]);
    }

    reconstruction.CorrectStackIntensities(corrected_stacks);
    for (unsigned int i=0;i<stacks.size();i++)
    {
      sprintf(buffer,"correctedstack%i.nii.gz",i);
      corrected_stacks[i].Write(buffer);
    }
  }
  
  //exit(1);
  
  //TEST 17/04/19 do not update slices but reset bias and scale
  //reconstruction.UpdateSlices(corrected_stacks, thickness);
  vector<irtkRealImage> old_bias;
  vector<double> old_scale;
  //TEST 2 17/04/19 Do not reset
  reconstruction.ResetBiasAndScale(old_bias,old_scale);
  
  
  
  
  
  
  
  /////////////////////////////////////////////////
  //Beginning of SH part
  
  cerr.rdbuf(fileSH.rdbuf());
  cout.rdbuf (fileSH.rdbuf());

  
  
   reconstruction.SetSmoothingParameters(delta, lastIterLambda);
   reconstruction.SetSigma(sigma);

  
   /*
  if(recon_1D)
    reconstruction.Set1DRecon();
  else
    reconstruction.Set3DRecon();
  
  if(iterations==0)
    reconstruction.InitializeEMValues();
  reconstruction.CoeffInit();
  reconstruction.GaussianReconstruction4D(nStacks);
  //reconstruction.SaveBiasFields();
  //exit(1);
  */
  
  /////////TESTING THIS 
  // Reset intensity matching and robust statistics
  if (_reset_robust_statistics)
  {
    reconstruction.InitializeEMValues();
    reconstruction.InitializeRobustStatistics();
  }
  /////////////TESTING this
  
  reconstruction.Init4D(nStacks);
  //directions
  int nDir = nStacks;
  irtkMatrix dirs_xyz(nDir,3);
  for(int i=0;i<nDir;i++)
    for(int j=0;j<3;j++)
      dirs_xyz(i,j) = directions[j][i+1];
  dirs_xyz.Print();
  
  //SH reconstruction
  reconstruction.SetLambdaLB(lambdaLB);
  //TEST
  //reconstruction.SetLambdaLB(1);//to create order=2 effectively
  
  reconstruction.InitSH(dirs_xyz,order);
  
  //reconstruction.Set3DRecon();
  //reconstruction.CoeffInit();
  reconstruction.SimulateSlicesDTI();
  //reconstruction.SaveSimulatedSlices();

  //TEST
  //reconstruction.SetSigma(10);
  
  for(int sh_iter=1; sh_iter<=1;sh_iter++)
  {
    
  for(int iteration = 0; iteration < sr_sh_iterations; iteration++ )
  {
      /*
      if (iteration == 0)
          reconstruction.SetLambdaLB(4.41);
      if (iteration == int_matching_first_iter)
          reconstruction.SetLambdaLB(0);
      */
     /* 
    if (iteration == 0)
    {
      reconstruction.SetSigma(20);
      reconstruction.SetLambdaLB(2);
      
      reconstruction.SimulateStacksDTIIntensityMatching(stacks,true);
      for (unsigned int i=0;i<stacks.size();i++)
      {
        sprintf(buffer,"simulated%i-%i.nii.gz",i,iteration);
        stacks[i].Write(buffer);
      }
    }

    
    if (iteration == 3)
    {
      reconstruction.SetSigma(20);
      reconstruction.SetLambdaLB(2);
      
      reconstruction.SimulateStacksDTIIntensityMatching(stacks,true);
      for (unsigned int i=0;i<stacks.size();i++)
      {
        sprintf(buffer,"simulated%i-%i.nii.gz",i,iteration);
        stacks[i].Write(buffer);
      }

    }

    if (iteration == 6)
    {
      reconstruction.SetSigma(20);
      reconstruction.SetLambdaLB(2);

      reconstruction.SimulateStacksDTIIntensityMatching(stacks,true);
      for (unsigned int i=0;i<stacks.size();i++)
      {
        sprintf(buffer,"simulated%i-%i.nii.gz",i,iteration);
        stacks[i].Write(buffer);
      }

    }
    
    if (iteration == 9)
    {
      reconstruction.SetSigma(20);
      reconstruction.SetLambdaLB(2);
      
      reconstruction.SimulateStacksDTIIntensityMatching(stacks,true);
      for (unsigned int i=0;i<stacks.size();i++)
      {
        sprintf(buffer,"simulated%i-%i.nii.gz",i,iteration);
        stacks[i].Write(buffer);
      }
    }

*/

    if(debug)
    {
      reconstruction.SimulateSignal(iteration);
    }
    
    
    if ((intensity_matching_sh)&&(iteration>=int_matching_first_iter))
    {
      //calculate bias fields
      if (sigma>0)
	reconstruction.Bias();
        //calculate scales
        reconstruction.Scale();
        //reconstruction.NormaliseBiasDTI(iteration,stack_transformations,order);

      }

    cout<<"regDTI="<<regDTI<<endl;
    if(regDTI>0)
      reconstruction.SuperresolutionDTI(iteration,true,sh_alpha,regDTI);
    else
      reconstruction.SuperresolutionDTI(iteration,false,sh_alpha);
    /*
    if(iteration==1)  
    {
      reconstruction.SaveWeights();
      reconstruction.SaveBiasFields();
      reconstruction.SaveSlices();
      reconstruction.SaveSimulatedSlices();
      exit(1);
    }
    */
    if (intensity_matching_sh)
    {
      //if(sigma>0)
	//reconstruction.NormaliseBiasSH(iteration);
    }
    
    reconstruction.SimulateSlicesDTI();
    cout.rdbuf (fileConsistency.rdbuf());
    cout<<reconstruction.Consistency()<<", ";
    cout.rdbuf (fileConsistencyEx.rdbuf());
    cout<<reconstruction.Consistency(true)<<", ";
    cout.rdbuf (strm_buffer);
    cout<<"SH iter "<<iteration<<": consistency = "<<reconstruction.Consistency()<<endl;
    cout.rdbuf (fileSH.rdbuf());
    cout<<"SH iter "<<iteration<<": consistency = "<<reconstruction.Consistency()<<endl;
    
    if(robust_statistics_sh)
        reconstruction.MStep(i+1);
      
      //E-step
      if(robust_statistics_sh)
        reconstruction.EStep();
    reconstruction.SaveSHcoeffs(iteration);
    /*
    if(robust_statistics)
      reconstruction.MStep(i+1);
    //E-step
    if(robust_statistics)
      reconstruction.EStep();
    */
  
   /* 
  reconstruction.SimulateStacksDTIIntensityMatching(stacks,true);
  for (unsigned int i=0;i<stacks.size();i++)
  {
    sprintf(buffer,"simulated-it%i-%i.nii.gz",iteration,i);
    stacks[i].Write(buffer);
  }
    */
  }
      if(sh_iter==0)
    {
      reconstruction.SliceToVolumeRegistrationSH();
      
      if(motion_model_hs)
	reconstruction.PostProcessMotionParametersHS(stacks[templateNumber]);
      else
	reconstruction.PostProcessMotionParameters(stacks[templateNumber]);

    }

  }
  
  
  ///////////////////
  //end of SH part
 
  reconstruction.SimulateSignal(output_name=NULL);//dirs_xyz,order);
 
  if (debug || save_transformations)
      reconstruction.SaveTransformations();
  
  if (debug)
  {
    reconstruction.SaveSimulatedSlices();
    reconstruction.SaveWeights();
    reconstruction.SaveSlices();
    for (unsigned int i=0;i<stacks.size();i++)
    {
      sprintf(buffer,"origstack%i.nii.gz",i);
      stacks[i].Write(buffer);
    }
  }

  //reconstruction.SimulateStacksDTI(stacks,true);
  double consistency = reconstruction.Consistency();

  
  if ( ! no_log ) {
        cout.rdbuf (fileEv.rdbuf());
  }
  reconstruction.Evaluate(10);
  cout<<endl;
  cout<<"Consistency = "<<consistency<<endl;

   if ( ! no_log ) {
       cout.rdbuf (fileSH.rdbuf());
  }

  if (debug)
  {
    corrected_stacks.clear();
    for (unsigned int i=0;i<stacks.size();i++)
    {
      corrected_stacks.push_back(stacks[i]);
    }
  
    reconstruction.SimulateStacksDTI(stacks,false);
    //  reconstruction.SimulateStacksDTIIntensityMatching(stacks,true);
    for (unsigned int i=0;i<stacks.size();i++)
    {
      sprintf(buffer,"simulated-exc%i.nii.gz",i);
      stacks[i].Write(buffer);
    }
  }

  //TEST 17/04/19 switching of correction after volume
  //reconstruction.AddBiasAndScale(old_bias, old_scale);
    
  reconstruction.SimulateStacksDTIIntensityMatching(stacks,true);
  for (unsigned int i=0;i<stacks.size();i++)
  {
    sprintf(buffer,"simulated%i.nii.gz",i);
    stacks[i].Write(buffer);
  }

  /*
  if (debug)
  {
    reconstruction.SimulateStacksDTI(stacks,true);
    for (unsigned int i=0;i<stacks.size();i++)
    {
      sprintf(buffer,"simulated-nocor%i.nii.gz",i);
      stacks[i].Write(buffer);
    }
  
    reconstruction.CorrectStackIntensities(corrected_stacks);
    for (unsigned int i=0;i<stacks.size();i++)
    {
      sprintf(buffer,"correctedstacksh%i.nii.gz",i);
      corrected_stacks[i].Write(buffer);
    }
  }
  */
  if ( ! no_log ) {
       cout.rdbuf (strm_buffer);
  }
  cout<<"Consistency = "<<consistency<<endl;
  //The end of main()
}  
