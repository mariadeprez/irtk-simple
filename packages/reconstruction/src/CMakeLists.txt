SET(RECONSTRUCTION_INCLUDES
../include/irtkBSplineReconstruction.h
../include/irtkReconstruction.h
../include/irtkReconstructionDTI.h
../include/irtkSphericalHarmonics.h
../include/irtkLaplacianSmoothing.h
../include/irtkReconstructionb0.h
)

SET(RECONSTRUCTION_SRCS
irtkBSplineReconstruction.cc
irtkReconstruction.cc
irtkReconstructionDTI.cc
irtkSphericalHarmonics.cc
irtkLaplacianSmoothing.cc
irtkReconstructionb0.cc
)



ADD_LIBRARY(reconstruction++ ${RECONSTRUCTION_INCLUDES} ${RECONSTRUCTION_SRCS})
INSTALL_FILES(/include FILES ${RECONSTRUCTION_INCLUDES})

