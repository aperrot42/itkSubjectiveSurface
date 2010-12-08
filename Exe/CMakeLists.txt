SET( CONFOREG_EXE_ITK_SRC
#seed
#  MedianFiltering
#  ResampleImage
#  BinaryThresholdFiltering
#  DistanceMapCalculation
#  KeepRangeFiltering
#  DilateDistanceMap
#  LocalMaximaExtraction
#  ImageMasking
#  SeedsGrouping
#  HysteresisThresholding
#utilities and interfaces
#  SeedImageToSeedTxt
#  RescaleDoubleToUChar
#  ValidationStatistics
#  RegionsCentroids
#NucleiRegenerate
#  RelabelLabelImage
  Subjective
  Subjective2D
  Subjective2Dopt
  )

FOREACH( var ${CONFOREG_EXE_ITK_SRC} )
  ADD_EXECUTABLE( ${var} ${var} )
  TARGET_LINK_LIBRARIES( ${var}
  ITKIO
  ITKCommon
 # ITKReview
  )
ENDFOREACH( var ${CONFOREG_EXE_ITK_SRC} )

#TARGET_LINK_LIBRARIES( VTKClustering
#  VTK
#  )
