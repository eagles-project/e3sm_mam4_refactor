# Relevant Directories within the HOMME repository
# holding source files needed for building the library
# and executable preqx.

macro(preqx_acc_setup)


INCLUDE_DIRECTORIES ( ${Trilinos_INCLUDE_DIRS} ${Trilinos_TPL_INCLUDE_DIRS} )
LINK_DIRECTORIES    ( ${Trilinos_LIBRARY_DIRS} ${Trilinos_TPL_LIBRARY_DIRS} )

SET(TRILINOS_ZOLTAN_DIR ${HOMME_SOURCE_DIR}/src/zoltan)
SET(TARGET_DIR        ${HOMME_SOURCE_DIR}/src/preqx_acc)
SET(UTILS_SHARE_DIR   ${HOMME_SOURCE_DIR}/utils/csm_share)
SET(SRC_DIR           ${HOMME_SOURCE_DIR}/src)
SET(SRC_SHARE_DIR     ${HOMME_SOURCE_DIR}/src/share)
SET(TEST_SRC_DIR      ${HOMME_SOURCE_DIR}/src/test_src)
SET(UTILS_TIMING_DIR  ${HOMME_BINARY_DIR}/utils/cime/CIME/non_py/src/timing)

# Make INCLUDE_DIRS global so the tests can access it
SET (EXEC_INCLUDE_DIRS ${PIO_INCLUDE_DIRS}
                            ${UTILS_TIMING_DIR} )

SET(TARGET_DIR        ${HOMME_SOURCE_DIR}/src/preqx_acc)
FILE(GLOB TARGET_F90  ${TARGET_DIR}/*.F90 ${SRC_SHARE_DIR}/*.F90 ${SRC_DIR}/preqx/share/*.F90 ${TEST_SRC_DIR}/*.F90)

SET(PREQX_SRCS_F90
  ${TARGET_F90}
  ${SRC_DIR}/checksum_mod.F90 
  ${SRC_DIR}/common_io_mod.F90
  ${SRC_DIR}/common_movie_mod.F90 
  ${SRC_DIR}/interp_movie_mod.F90 
  ${SRC_DIR}/interpolate_driver_mod.F90 
  ${SRC_DIR}/netcdf_io_mod.F90 
  ${SRC_DIR}/pio_io_mod.F90 
  ${SRC_DIR}/prim_main.F90 
  ${SRC_DIR}/prim_movie_mod.F90 
  ${SRC_DIR}/prim_restart_mod.F90 
  ${SRC_DIR}/repro_sum_mod.F90 
  ${SRC_DIR}/restart_io_mod.F90 
  ${SRC_DIR}/surfaces_mod.F90 
  ${SRC_DIR}/test_mod.F90
  ${SRC_DIR}/restart_io_mod.F90 
  ${SRC_DIR}/common_io_mod.F90
  ${SRC_DIR}/pio_io_mod.F90 
  ${SRC_DIR}/repro_sum_mod.F90 
  ${SRC_DIR}/interp_movie_mod.F90 
  ${SRC_DIR}/common_movie_mod.F90 
  ${SRC_DIR}/checksum_mod.F90 
  ${SRC_DIR}/prim_restart_mod.F90 
  ${SRC_DIR}/prim_movie_mod.F90 
  ${SRC_DIR}/surfaces_mod.F90 
  ${SRC_DIR}/interpolate_driver_mod.F90 
  ${SRC_DIR}/netcdf_io_mod.F90 
  ${SRC_DIR}/prim_main.F90 
  ${SRC_DIR}/test_mod.F90
  ${UTILS_SHARE_DIR}/shr_kind_mod.F90 
  ${UTILS_SHARE_DIR}/shr_mpi_mod.F90 
  ${UTILS_SHARE_DIR}/shr_sys_mod.F90 
  ${UTILS_SHARE_DIR}/shr_file_mod.F90 
  ${UTILS_SHARE_DIR}/shr_vmath_mod.F90 
  ${UTILS_SHARE_DIR}/shr_const_mod.F90 
  ${UTILS_SHARE_DIR}/shr_spfn_mod.F90 
  )

IF (${HOMME_USE_TRILINOS})
SET (PREQX_SRCS_ZOLTAN   ${TRILINOS_ZOLTAN_DIR}/zoltan_interface.c  ${TRILINOS_ZOLTAN_DIR}/zoltan_cppinterface.cpp)
SET_SOURCE_FILES_PROPERTIES( ${TRILINOS_ZOLTAN_DIR}/zoltan_cppinterface.cpp PROPERTIES LANGUAGE CXX )
ENDIF()


# If the user specified a file for custom compiler options use those
IF (DEFINED PREQX_CUSTOM_FLAGS_FILE)
  setCustomCompilerFlags(PREQX_CUSTOM_FLAGS_FILE PREQX_SRCS_F90)
ENDIF ()

# Make SRCS global so the tests can access it
SET(EXEC_SOURCES ${PREQX_SRCS_F90} ${PREQX_SRCS_ZOLTAN})

# Set up defaults
IF (NOT PREQX_NP)
  SET (PREQX_NP 4) 
ENDIF ()
IF (NOT PREQX_PLEV)
  SET (PREQX_PLEV 20) 
ENDIF ()
IF (NOT PREQX_USE_PIO)
  SET (PREQX_USE_PIO FALSE) 
ENDIF ()
IF (NOT PREQX_USE_ENERGY)
  SET (PREQX_USE_ENERGY FALSE) 
ENDIF ()
IF (NOT PREQX_NC)
  SET (PREQX_NC 4)
ENDIF ()
IF (NOT QSIZE_D)
  SET (QSIZE_D 4)
ENDIF ()


SET(USE_OPENACC TRUE)
SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OPENACC_Fortran_FLAGS}")
SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OPENACC_Linker_FLAGS}")

# nested openMP would never be used with openACC
SET(NESTED_OPENMP FALSE )
SET(COLUMN_OPENMP FALSE )

endmacro(preqx_acc_setup)

############################################################################
# createTestExec(exec_name exec_type NP PLEV USE_PIO USE_ENERGY)
############################################################################

preqx_acc_setup()
createTestExec(preqx_acc preqx_acc ${PREQX_NP} ${PREQX_NC} ${PREQX_PLEV} 
               ${PREQX_USE_PIO}  ${PREQX_USE_ENERGY} ${QSIZE_D})


