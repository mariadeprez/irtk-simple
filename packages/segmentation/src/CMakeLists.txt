SET(SEGMENTATION_INCLUDES
../include/irtkBiasCorrection.h
../include/irtkBiasField.h
../include/irtkBSplineBiasField.h
../include/irtkEMClassification.h
../include/irtkEMClassificationBiasCorrection.h
../include/irtkEMClassificationBiasCorrectionfMRI.h
../include/irtkEMClassificationMultiComp.h
../include/irtkGaussian.h
../include/irtkIntensityMatching.h
../include/irtkMultiChannelImage.h
../include/irtkMeanShift.h
../include/irtkProbabilisticAtlas.h
../include/irtkEMClassification2ndOrderMRF.h
../include/irtkBiasCorrectionMask.h
../include/irtkPolynomialBiasField.h
../include/irtkGaussianLikelihood.h
../include/irtkMyelinClassification.h
../include/irtkMyelinClassificationPV.h 
../include/irtkProbabilityMaps.h
)

SET(SEGMENTATION_SRCS
irtkBiasCorrection.cc
irtkBiasField.cc
irtkBSplineBiasField.cc
irtkEMClassification.cc
irtkEMClassificationBiasCorrection.cc
irtkEMClassificationBiasCorrectionfMRI.cc
irtkEMClassificationMultiComp.cc
irtkGaussian.cc
irtkIntensityMatching.cc
irtkMultiChannelImage.cc
irtkMeanShift.cc
irtkProbabilisticAtlas.cc
irtkEMClassification2ndOrderMRF.cc
irtkBiasCorrectionMask.cc
irtkPolynomialBiasField.cc
irtkGaussianLikelihood.cc 
irtkMyelinClassification.cc 
irtkMyelinClassificationPV.cc 
irtkProbabilityMaps.cc

)



ADD_LIBRARY(segmentation++ ${SEGMENTATION_INCLUDES} ${SEGMENTATION_SRCS})
INSTALL_FILES(/include FILES ${SEGMENTATION_INCLUDES})

