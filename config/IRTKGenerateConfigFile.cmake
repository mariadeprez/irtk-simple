# Generate the ITRKConfig.cmake file in the build tree.  Also configure
# one for installation.  The file tells external projects how to use
# ITRK.

SET(IRTK_INCLUDE_DIRS 
	${PROJECT_SOURCE_DIR}/recipes/include
	${PROJECT_SOURCE_DIR}/common++/include
	${PROJECT_SOURCE_DIR}/geometry++/include
	${PROJECT_SOURCE_DIR}/image++/include
	${PROJECT_SOURCE_DIR}/contrib++/include
	${PROJECT_SOURCE_DIR}/packages/transformation/include
	${PROJECT_SOURCE_DIR}/packages/registration/include
	${PROJECT_SOURCE_DIR}/packages/segmentation/include
	${PROJECT_SOURCE_DIR}/packages/reconstruction/include)

SET(IRTK_LIBRARIES_DIR ${PROJECT_BINARY_DIR}/lib/${CMAKE_CFG_INTDIR})

IF (BUILD_WITH_TBB)
	SET(BUILD_WITH_TBB_FLAG 1)
ELSE(BUILD_WITH_TBB)
	SET(BUILD_WITH_TBB_FLAG 0)
ENDIF(BUILD_WITH_TBB)

IF (ZLIB_FOUND)
	SET(ZLIB_FOUND_FLAG 1)
ELSE(ZLIB_FOUND)
	SET(ZLIB_FOUND_FLAG 0)
ENDIF(ZLIB_FOUND)


IF (BUILD_WITH_NIFTI)
	SET(BUILD_WITH_NIFTI_FLAG 1)
ELSE(BUILD_WITH_NIFTI)
	SET(BUILD_WITH_NIFTI_FLAG 0)
ENDIF(BUILD_WITH_NIFTI)


CONFIGURE_FILE(${CMAKE_CURRENT_LIST_DIR}/IRTKConfig.cmake.in
               ${PROJECT_BINARY_DIR}/IRTKConfig.cmake @ONLY IMMEDIATE)
