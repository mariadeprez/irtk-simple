SET(RECONSTRUCTION_INCLUDES
../include/irtkBSplineReconstruction.h
../include/irtkReconstruction.h
)

SET(RECONSTRUCTION_SRCS
irtkBSplineReconstruction.cc
irtkReconstruction.cc
)



ADD_LIBRARY(reconstruction++ ${RECONSTRUCTION_INCLUDES} ${RECONSTRUCTION_SRCS})
INSTALL_FILES(/include FILES ${RECONSTRUCTION_INCLUDES})

