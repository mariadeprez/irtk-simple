
# SLICE TO VOLUME RECONSTRUCTION (SVR)

During MRI acquisition fetus moves. Fast acquisition of 2D slices freezes motion in time. Between-slice motion is the corrected with slice-to-volume registration and the 3D volume is reconstructed using super-resolution reconstruction.

add image

[IRTK SIMPLE](https://github.com/mariadeprez/irtk-simple) is a  software package has been created to distribute software for motion correction and super-resolution reconstruction of **structural and diffusion MRI of the fetal brain**. The reconstruction algorithms
were developed and implemented by [Dr Maria Deprez](https://www.kcl.ac.uk/people/maria-deprez) (Kuklisova Murgasova) within [IRTK](https://biomedia.doc.ic.ac.uk/software/irtk/) sofware package.

- [installation](https://github.com/mariadeprez/irtk-simple#Installation): Install and compile on UBUNTU
- [references](https://github.com/mariadeprez/irtk-simple#References): Read detailed description of reconstruction algorithms
- [about irtk-simple] (https://github.com/mariadeprez/irtk-simple#about):


References
----------

The irtk-simple package contains three main algorithms:

1. Motion correction and super-resolution reconstruction of fetal structural MRI

Reconstruction of fetal brain MRI with intensity matching and complete outlier removal. 
M Kuklisova-Murgasova, G Quaghebeur, MA Rutherford, JV Hajnal, JA Schanbel. 
Medical Image Analysis 16(8), 2012. doi: [10.1016/j.media.2012.07.004](https://doi.org/10.1016/j.media.2012.07.004)

About irtk-simple
-----
It is called SIMPLE IRTK - because it is based on a lightweight and easy to compile 
subset of package Image Registration Toolkit (IRTK).

The package is parallelized using Threading Building Blocks (TBB), which is included
to ensure compatibility of the version. 

Installation
------------

## INSTRUCTIONS for UBUNTU:

* Install libraries:

	sudo apt-get install git cmake cmake-curses-gui g++ zlib1g-dev 

* Download package from the Gitlab repository

	git clone https://github.com/mariadeprez/irtk-simple.git
  	
	<br>(you will need to enter your username and the personal token)

* Build

	cd irtk-simple <br>
	mkdir build <br>
	cd build <br>
	cmake .. <br>
	make

* run the executable to display input arguments and options

	irtk-simple/build/bin/reconstruction

The motion correction software was created and is distributed by

Maria Deprez (Kuklisova Murgasova) <br> 
Department of Biomedical Engineering <br>
King's College London <br>

*************************************************************************************
*************************************************************************************		

## ABOUT Image Registration Toolkit (IRTK)

Authors: <br> 
Daniel Rueckert <br>
Julia Schnabel <br>

See the COPYRIGHT file in this directory for more information on the copyright and 
license agreement for the software.

The non-rigid registration algorithms implemented in the software are described by 
the following publications:

J. A. Schnabel, D. Rueckert, M. Quist, J. M. Blackall, A. D. Castellano Smith, 
T. Hartkens, G. P. Penney, W. A. Hall, H. Liu, C. L. Truwit, F. A. Gerritsen, 
D. L. G. Hill, and D. J. Hawkes. A generic framework for non-rigid registration 
based on non-uniform multi-level free-form deformations. In Fourth Int. Conf. on 
Medical Image Computing and Computer-Assisted Intervention (MICCAI '01), pages 
573-581, Utrecht, NL, October 2001

D. Rueckert, L. I. Sonoda, C. Hayes, D. L. G. Hill, M. O. Leach, and D. J. Hawkes. 
Non-rigid registration using free-form deformations: Application to breast MR 
images. IEEE Transactions on Medical Imaging, 18(8):712-721, 1999.

E. R. E. Denton, L. I. Sonoda, D. Rueckert, S. C. Rankin, C. Hayes, M. Leach, D. 
L. G. Hill, and D. J. Hawkes. Comparison and evaluation of rigid and non-rigid 
registration of breast MR images. Journal of Computer Assisted Tomography, 
23:800-805, 1999.


*************************************************************************************
*************************************************************************************

Disclaimer:

This software has been developed for research purposes only, and hence should not be 
used as a diagnostic tool. In no event shall the authors or distributors be liable to 
any direct, indirect, special, incidental, or consequential damages arising of the use 
of this software, its documentation, or any derivatives thereof, even if the authors 
have been advised of the possibility of such damage.


