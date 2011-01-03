PROJECT(ITKSUBJECTIVESURFACES)

# CMake >= 2.6 has to be installed (fighting for progress)
CMAKE_MINIMUM_REQUIRED( VERSION 2.6 )

MESSAGE( STATUS "SUBJECTIVE SURFACES" )

# On Visual Studio 8 MS deprecated C. This removes all 1.276E1265 security
# warnings. Copied from ITK CMakeLists.
IF( WIN32 )
  IF( NOT BORLAND )
    IF( NOT CYGWIN )
      IF( NOT MINGW )
          ADD_DEFINITIONS(
            -D_CRT_FAR_MAPPINGS_NO_DEPRECATE
            -D_CRT_IS_WCTYPE_NO_DEPRECATE
            -D_CRT_MANAGED_FP_NO_DEPRECATE
            -D_CRT_NONSTDC_NO_DEPRECATE
            -D_CRT_SECURE_NO_DEPRECATE
            -D_CRT_SECURE_NO_DEPRECATE_GLOBALS
            -D_CRT_SETERRORMODE_BEEP_SLEEP_NO_DEPRECATE
            -D_CRT_TIME_FUNCTIONS_NO_DEPRECATE
            -D_CRT_VCCLRIT_NO_DEPRECATE
            -D_SCL_SECURE_NO_DEPRECATE
            )
      ENDIF( NOT MINGW )
    ENDIF( NOT CYGWIN )
  ENDIF( NOT BORLAND )
ENDIF( WIN32 )


#--------------------------------------------------------------------------
# Define CMAKE_whatever_OUTPUT_DIRECTORY
#--------------------------------------------------------------------------
# all target executables will be
# in ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
SET( CMAKE_RUNTIME_OUTPUT_DIRECTORY
  ${ITKSUBJECTIVESURFACES_BINARY_DIR}/bin
)

# all target static libraries will be in ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}
#SET( CMAKE_ARCHIVE_OUTPUT_DIRECTORY
#  ${ITKSUBJECTIVESURFACES_BINARY_DIR}/lib
#)

# all shared libraries will be in ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
# Note that modules (plugins) are considered as shared libraries
#SET( CMAKE_LIBRARY_OUTPUT_DIRECTORY
#  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/plugins
#   ${SEED_BINARY_DIR}/lib
#)

LINK_DIRECTORIES(
  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
  ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}
  ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
  )


#--------------------------------------------------------------------------
# Define TESTING_DATA_PATH
#--------------------------------------------------------------------------
#SET( TESTING_DATA_PATH
#  ${SEED_SOURCE_DIR}/Testing/Data
#  CACHE PATH "Directory which contains data for testing."
#)



#--------------------------------------------------------------------------
# ITK
#--------------------------------------------------------------------------
FIND_PACKAGE( ITK )
IF( ITK_FOUND )
  INCLUDE( ${ITK_USE_FILE} )
ELSE( ITK_FOUND )
  MESSAGE( FATAL_ERROR "Cannot build ITK project without ITK. Please set ITK_DIR.")
ENDIF( ITK_FOUND )


#--------------------------------------------------------------------------
# VTK
#--------------------------------------------------------------------------
#FIND_PACKAGE( VTK )
#IF( VTK_FOUND )
#  INCLUDE( ${VTK_USE_FILE} )
#ELSE( VTK_FOUND )
#  MESSAGE( FATAL_ERROR "Cannot build project without VTK. Please set VTK_DIR.")
#ENDIF( VTK_FOUND )


#-------------------------------------------------------------------------
# PROJECT STRUCTURE
#-------------------------------------------------------------------------
# directory containing the code
SET( ITKSUBJECTIVESURFACES_CODE_DIR
  ${ITKSUBJECTIVESURFACES_SOURCE_DIR}/Code/Classes
)

INCLUDE_DIRECTORIES( BEFORE
  ${ITKSUBJECTIVESURFACES_CODE_DIR}
)

ADD_SUBDIRECTORY( Exe )

