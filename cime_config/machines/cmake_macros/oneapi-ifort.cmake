set(FFLAGS "-convert big_endian -assume byterecl -traceback -assume realloc_lhs -fp-model consistent")
if (compile_threaded)
  string(APPEND FFLAGS " -qopenmp")
endif()
if (NOT DEBUG)
  string(APPEND FFLAGS " -O2")
endif()
if (DEBUG)
  string(APPEND FFLAGS " -O0 -g -check uninit -check bounds -check pointers -fpe0 -check noarg_temp_created")
endif()
set(CFLAGS "-fp-model precise -std=gnu99 -traceback")
if (compile_threaded)
  string(APPEND CFLAGS " -qopenmp")
endif()
if (NOT DEBUG)
  string(APPEND CFLAGS " -O2")
endif()
if (DEBUG)
  string(APPEND CFLAGS " -O0 -g")
endif()
set(CXXFLAGS "-std=c++14 -fp-model precise -traceback")
if (compile_threaded)
  string(APPEND CXXFLAGS " -qopenmp")
endif()
if (NOT DEBUG)
  string(APPEND CXXFLAGS " -O2")
endif()
if (DEBUG)
  string(APPEND CXXFLAGS " -O0 -g")
endif()
set(SUPPORTS_CXX "TRUE")
set(CXX_LINKER "FORTRAN")
set(CXX_LDFLAGS "-cxxlib")
string(APPEND CPPDEFS " -DFORTRANUNDERSCORE -DNO_R16 -DCPRINTEL -DHAVE_SLASHPROC")
set(FC_AUTO_R8 "-r8")
set(FFLAGS_NOOPT "-O0")
set(FIXEDFLAGS "-fixed -132")
set(FREEFLAGS "-free")
set(HAS_F2008_CONTIGUOUS "TRUE")
set(MPIFC "mpif90")
set(MPICC "mpicc")
set(MPICXX "mpicxx")
set(SCC "icc")
set(SCXX "icpc")
set(SFC "ifort")
if (compile_threaded)
  string(APPEND LDFLAGS " -qopenmp")
endif()
execute_process(COMMAND $ENV{NETCDF_PATH}/bin/nf-config --flibs OUTPUT_VARIABLE SHELL_CMD_OUTPUT_BUILD_INTERNAL_IGNORE0 OUTPUT_STRIP_TRAILING_WHITESPACE)
string(APPEND SLIBS " ${SHELL_CMD_OUTPUT_BUILD_INTERNAL_IGNORE0} -Wl,-rpath -Wl,$ENV{NETCDF_PATH}/lib -mkl")
execute_process(COMMAND $ENV{NETCDF_PATH}/bin/nc-config --libs OUTPUT_VARIABLE SHELL_CMD_OUTPUT_BUILD_INTERNAL_IGNORE0 OUTPUT_STRIP_TRAILING_WHITESPACE)
string(APPEND SLIBS " ${SHELL_CMD_OUTPUT_BUILD_INTERNAL_IGNORE0}")
set(NETCDF_PATH "$ENV{NETCDF_PATH}")
set(PNETCDF_PATH "$ENV{PNETCDF_PATH}")