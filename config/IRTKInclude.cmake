# Compatible with CMake 2.6.x
GET_FILENAME_COMPONENT(CMAKE_CURRENT_LIST_DIR ${CMAKE_CURRENT_LIST_FILE} PATH)

# Prefer project-own Find<Package>.cmake modules for FIND_PACKAGE().
SET(CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR})


# Option to produce multi-threaded executables using TBB
OPTION(BUILD_WITH_TBB "Build multi-threaded executables using TBB" OFF)

IF (BUILD_WITH_TBB)
  FIND_PACKAGE(TBB REQUIRED)
  IF (TBB_FOUND)
    # Attention: DO NOT define TBB_DEPRECATED by default or before including the
    #            other TBB header files, in particular parallel_for. The deprecated
    #            behavior of parallel_for is to not choose the chunk size (grainsize)
    #            automatically!
    #
    # http://software.intel.com/sites/products/documentation/doclib/tbb_sa/help/tbb_userguide/Automatic_Chunking.htm
    ADD_DEFINITIONS(-DHAS_TBB)
    INCLUDE_DIRECTORIES(${TBB_INCLUDE_DIRS})
    LINK_DIRECTORIES(${TBB_LIBRARY_DIRS})
    LINK_LIBRARIES(${TBB_LIBRARIES})
  ENDIF (TBB_FOUND)
ENDIF(BUILD_WITH_TBB)

INCLUDE(${CMAKE_ROOT}/Modules/FindZLIB.cmake)

IF (ZLIB_FOUND)
  ADD_DEFINITIONS(-DHAS_ZLIB -DHAVE_ZLIB)
  INCLUDE_DIRECTORIES(${ZLIB_INCLUDE_DIR})
  LINK_LIBRARIES(${ZLIB_LIBRARIES})
ENDIF (ZLIB_FOUND)

IF (WIN32)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /TP /Ze /W0")
  ADD_DEFINITIONS(-DvtkCommon_EXPORTS)
ELSE (WIN32)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-deprecated -Wno-write-strings")
ENDIF (WIN32)

ADD_DEFINITIONS(-DIMPERIAL -DANSI -DHAS_CONTRIB -DNO_BOUNDS -DENABLE_UNIX_COMPRESS)

INCLUDE_DIRECTORIES(${IRTK_SOURCE_DIR}/recipes/include)
INCLUDE_DIRECTORIES(${IRTK_SOURCE_DIR}/common++/include)
INCLUDE_DIRECTORIES(${IRTK_SOURCE_DIR}/geometry++/include)
INCLUDE_DIRECTORIES(${IRTK_SOURCE_DIR}/image++/include)
INCLUDE_DIRECTORIES(${IRTK_SOURCE_DIR}/contrib++/include)
INCLUDE_DIRECTORIES(${IRTK_SOURCE_DIR}/packages/transformation/include)
INCLUDE_DIRECTORIES(${IRTK_SOURCE_DIR}/packages/registration/include)
INCLUDE_DIRECTORIES(${IRTK_SOURCE_DIR}/packages/segmentation/include)
INCLUDE_DIRECTORIES(${IRTK_SOURCE_DIR}/packages/reconstruction/include)

LINK_DIRECTORIES(${IRTK_BINARY_DIR}/lib) 

LINK_LIBRARIES(segmentation++ reconstruction++ registration++ transformation++ contrib++
image++ geometry++ common++ recipes) 

# Options to build with nifti, znz and possibly fslio
OPTION(BUILD_WITH_NIFTI "Build using NIFTI support" ON)
IF (BUILD_WITH_NIFTI)
   ADD_DEFINITIONS(-DHAS_NIFTI)
   INCLUDE_DIRECTORIES(${IRTK_SOURCE_DIR}/nifti/niftilib)
   INCLUDE_DIRECTORIES(${IRTK_SOURCE_DIR}/nifti/znzlib)
   LINK_LIBRARIES(znz)
   LINK_LIBRARIES(niftiio)
ENDIF (BUILD_WITH_NIFTI)

