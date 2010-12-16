SET( CONFOREG_EXE_ITK_SRC
#  Subjective
  Subjective2D
  Subjective2DVideo
  Subjective
  )

FOREACH( var ${CONFOREG_EXE_ITK_SRC} )
  ADD_EXECUTABLE( ${var} ${var}.cxx )
  TARGET_LINK_LIBRARIES( ${var}
  ITKIO
  ITKCommon
 # ITKReview
  )
ENDFOREACH( var ${CONFOREG_EXE_ITK_SRC} )

#TARGET_LINK_LIBRARIES( VTKClustering
#  VTK
#  )
