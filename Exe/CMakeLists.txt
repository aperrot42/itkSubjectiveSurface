SET( CONFOREG_EXE_ITK_SRC
  Subjective2D
  Subjective
  Subjective2DVideo
  Feature2D
  Feature
  Initialize2D
  Initialize
  GeodesicDiffusion
  )

FOREACH( var ${CONFOREG_EXE_ITK_SRC} )
  ADD_EXECUTABLE( ${var} ${var} )
  TARGET_LINK_LIBRARIES( ${var}
  ITKIO
  ITKCommon
  ITKStatistics
  )
ENDFOREACH( var ${CONFOREG_EXE_ITK_SRC} )

#TARGET_LINK_LIBRARIES( VTKClustering
#  VTK
#  )
