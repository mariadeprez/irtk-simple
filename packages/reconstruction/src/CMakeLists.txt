SET(RECONSTRUCTION_INCLUDES
../include/irtkBSplineReconstruction.h
../include/irtkReconstruction.h
../include/irtkReconstructionDTI.h
../include/irtkSphericalHarmonics.h
)

SET(RECONSTRUCTION_SRCS
irtkBSplineReconstruction.cc
irtkReconstruction.cc
irtkReconstructionDTI.cc
irtkSphericalHarmonics.cc
)



ADD_LIBRARY(reconstruction++ ${RECONSTRUCTION_INCLUDES} ${RECONSTRUCTION_SRCS})
INSTALL_FILES(/include FILES ${RECONSTRUCTION_INCLUDES})

