cmake_minimum_required(VERSION 2.8.8)
project(MPASAdaptor)

find_package(ParaView REQUIRED COMPONENTS vtkPVCatalyst vtkPVPythonCatalyst)
include(${PARAVIEW_USE_FILE})

option(BUILD_SHARED_LIBS "Build ParaView using shared libraries" ON)

enable_language(Fortran)
include(FortranCInterface)

FortranCInterface_HEADER(MPASAdaptorAPIMangling.h SYMBOLS
  coprocessor_create_grid 
  coprocessor_register_tracer_data coprocessor_add_tracer_data
  coprocessor_register_data coprocessor_add_data
  mpas_init_paraview_catalyst mpas_registerdata mpas_coprocess mpas_finalize_paraview_catalyst mpas_check_coprocess )

include_directories(${MPASAdaptor_BINARY_DIR})

set(MPASAdaptor_SRCS
  GridUtils.cxx
  MPASAdaptor.cxx
  MPAS_XY_NLayer.cxx
  MPAS_XYZ_1Layer.cxx
  MPAS_XYZ_NLayer.cxx
  MPAS_LonLat_1Layer.cxx
  MPAS_LonLat_NLayer.cxx
)

add_library(MPASAdaptor ${MPASAdaptor_SRCS})
target_link_libraries(
  MPASAdaptor vtkPVCatalyst vtkPVPythonCatalyst
)
