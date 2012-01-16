PROJECT(ONE_TIME_JUNK)

ADD_EXECUTABLE(subfield_slice_rules SubfieldSliceRules.cxx)
ADD_EXECUTABLE(ashs_affine_vol2slice ASHSAffineVolumeToSlice.cxx)
ADD_EXECUTABLE(subfield_mask_slices SubfieldSliceMask.cxx)
ADD_EXECUTABLE(subfield_leveler SubfieldLeveler.cxx)
ADD_EXECUTABLE(unc_mart2sean UNCConvertMartinImageToSeanImage.cxx)
ADD_EXECUTABLE(convert3d ConvertImage3D.cxx)
ADD_EXECUTABLE(randfx RandomEffects.cxx)
ADD_EXECUTABLE(hippov HippoOverlap.cxx)
ADD_EXECUTABLE(tempfile TempFile.c)
ADD_EXECUTABLE(maskmunge MaskMunge.cxx)
ADD_EXECUTABLE(findpeaks FindPeaks.cxx)
ADD_EXECUTABLE(overlaps SegmentationOverlaps.cxx)
ADD_EXECUTABLE(coordinate_images CoordinateImages.cxx)
ADD_EXECUTABLE(pmatlas_ibn_origin_shift OriginShiftHack.cxx)

SET(COMMON_LIBS
  ITKCommon 
  ITKBasicFilters 
  ITKIO 
  ITKNumerics
  ITKStatistics
  ITKVoxBoIO)

TARGET_LINK_LIBRARIES(subfield_slice_rules ${COMMON_LIBS})
TARGET_LINK_LIBRARIES(subfield_mask_slices ${COMMON_LIBS})
TARGET_LINK_LIBRARIES(subfield_leveler ${COMMON_LIBS})
TARGET_LINK_LIBRARIES(unc_mart2sean ${COMMON_LIBS})
TARGET_LINK_LIBRARIES(convert3d ${COMMON_LIBS})
TARGET_LINK_LIBRARIES(randfx ${COMMON_LIBS})
TARGET_LINK_LIBRARIES(hippov ${COMMON_LIBS})
TARGET_LINK_LIBRARIES(overlaps ${COMMON_LIBS})
TARGET_LINK_LIBRARIES(maskmunge ${COMMON_LIBS})
TARGET_LINK_LIBRARIES(findpeaks ${COMMON_LIBS})
TARGET_LINK_LIBRARIES(coordinate_images ${COMMON_LIBS})
TARGET_LINK_LIBRARIES(pmatlas_ibn_origin_shift ${COMMON_LIBS})
TARGET_LINK_LIBRARIES(ashs_affine_vol2slice ${COMMON_LIBS})

INSTALL_TARGETS("/bin" 
  convert3d coordinate_images randfx 
  tempfile maskmunge findpeaks
  overlaps pmatlas_ibn_origin_shift subfield_leveler subfield_mask_slices
  ashs_affine_vol2slice subfield_slice_rules)
