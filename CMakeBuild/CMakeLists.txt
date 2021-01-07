cmake_minimum_required(VERSION 3.17)

project(FDS)
enable_language(Fortran)

find_package(MPI)

if(NOT MPI_Fortran_FOUND)
    message(FATAL_ERROR "Could not find Fortran MPI.  Please set MPI_Fortran_COMPILER to point to the mpifort wrapper.")
endif()

include_directories(${MPI_Fortran_INCLUDE_PATH})

add_compile_options(${MPI_Fortran_COMPILE_FLAGS})
add_compile_options(-cpp)
#add_compile_options(--param=num-threads=24)

add_compile_definitions(GITHASH_PP="scarc")
add_compile_definitions(GITDATE_PP="1.1.1970")
add_compile_definitions(BUILDDATE_PP="4.1.1970")
add_compile_definitions(COMPVER_PP="0.0.0-pre2")

# requested source files
file(GLOB_RECURSE sources CONFIGURE_DEPENDS ../Source/*.f90)

#
# GNU Fortran compiler 
#
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")

  set(CMAKE_Fortran_FLAGS_BASIC "-m64 -std=f2008 -ffpe-summary=none -fall-intrinsics")
   
   # Debug version
   if (CMAKE_BUILD_TYPE MATCHES "Debug") 
     set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_BASIC} -O0 -ggdb -Wall -Wcharacter-truncation -Wno-target-lifetime -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow -frecursive")
     set (FDS "fds_gnu_db")
     message ("Compiling GNU debug version ${FDS} with options: ${CMAKE_Fortran_FLAGS_DEBUG}")
   endif ()

   # Release version
   if (CMAKE_BUILD_TYPE MATCHES "Release") 
     set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_BASIC} -O2")
     set (FDS "fds_gnu")
     message ("Compiling GNU release version ${FDS} with options: ${CMAKE_Fortran_FLAGS_RELEASE}")
   endif ()
endif()

#
# GNU Fortran compiler 
#
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")

   set(CMAKE_Fortran_FLAGS_BASIC "-m64 -mt_mpi -traceback -no-wrap-margin") 

   # Debug version
   if (CMAKE_BUILD_TYPE MATCHES "Debug") 
     set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_BASIC} -check all -warn all -O0 -auto -WB -g -fpe0 -fltconsistency -stand:f08")
     set (FDS "fds_intel_db")
     message ("Compiling Intel debug version ${FDS} with options: ${CMAKE_Fortran_FLAGS_DEBUG}")
   endif ()

   # Release version
   if (CMAKE_BUILD_TYPE MATCHES "Release") 
      set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_BASIC} -O2 -ipo")
      set (FDS "fds_intel")
      message ("Compiling Intel release version ${FDS} with options: ${CMAKE_Fortran_FLAGS_RELEASE}")
   endif ()
endif()


add_executable(${FDS} ${sources})

target_link_libraries(${FDS} ${MPI_Fortran_LIBRARIES})
set_property(TARGET ${FDS} APPEND_STRING PROPERTY LINK_FLAGS "${MPI_Fortran_LINK_FLAGS}")

#message (" CMAKE: MPI_Fortran_LIBRARIES : ${MPI_Fortran_LIBRARIES}")
