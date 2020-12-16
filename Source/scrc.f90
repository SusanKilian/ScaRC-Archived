! ================================================================================================================
!> \brief Scalable Recursive Clustering (ScaRC): Collection of alternative solvers for the FDS pressure equation
! ================================================================================================================
! Use of different directives possible
!  - WITH_MKL                    : use MKL routines PARDISO, CLUSTER_SPARSE_SOLVER, DDOT, DAXPY, DAXPBY, DCOPY, DSCAL
!  - WITH_SCARC_VERBOSE          : print more detailed information about ScaRC iterations and workspace allocation
!  - WITH_SCARC_DEBUG            : print detaild debugging info (only for developing purposes)
!  - WITH_SCARC_AMG              : include algebraic multigrid solver
!  - WITH_SCARC_POSTPROCESSING   : dump environment for separate ScaRC postprocessing program
! ================================================================================================================
#undef WITH_SCARC_DEBUG
#define WITH_SCARC_VERBOSE
#define WITH_SCARC_AMG
#undef WITH_SCARC_POSTPROCESSING

#include "scrc/module_scarc_constants.f90"
#include "scrc/module_scarc_types.f90"
#include "scrc/module_scarc_variables.f90"
#include "scrc/module_scarc_pointers.f90"
#include "scrc/module_scarc_message_services.f90"
#include "scrc/module_scarc_error_handling.f90"
#include "scrc/module_scarc_utilities.f90"
#include "scrc/module_scarc_memory_manager.f90"
#include "scrc/module_scarc_iteration_manager.f90"
#include "scrc/module_scarc_timings.f90"
#include "scrc/module_scarc_stack_administration.f90"
#include "scrc/module_scarc_initialization.f90"
#include "scrc/module_scarc_mpi.f90"
#ifdef WITH_MKL
#include "scrc/module_scarc_mkl_environment.f90"
#endif
#include "scrc/module_scarc_linear_algebra.f90"
#include "scrc/module_scarc_discretization.f90"
#include "scrc/module_scarc_matrix_systems.f90"
#include "scrc/module_scarc_fft_environment.f90"
#include "scrc/module_scarc_gmg_environment.f90"
#ifdef WITH_SCARC_AMG
#include "scrc/module_scarc_amg_environment.f90"
#endif
#include "scrc/module_scarc_mgm_environment.f90"
#include "scrc/module_scarc_methods.f90"


MODULE SCRC

USE PRECISION_PARAMETERS, ONLY: EB
USE GLOBAL_CONSTANTS
USE COMP_FUNCTIONS, ONLY: CURRENT_TIME
USE SCARC_CONSTANTS
USE SCARC_ERROR_HANDLING
USE SCARC_TIMINGS
#ifdef WITH_MKL
USE SCARC_MKL_ENVIRONMENT
#endif
USE SCARC_STACK_ADMINISTRATION, ONLY: N_STACK_TOTAL
USE SCARC_TIMINGS
USE SCARC_MESSAGE_SERVICES
USE SCARC_INITIALIZATION
USE SCARC_DISCRETIZATION
USE SCARC_METHODS

IMPLICIT NONE

PUBLIC :: SCARC_SETUP, SCARC_SOLVER

CONTAINS

! ------------------------------------------------------------------------------------------------
!> \brief Initialize ScaRC structures based on SCARC-input parameters from &PRES namelist
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP
USE SCARC_DISCRETIZATION
REAL(EB) :: TNOW

TNOW = CURRENT_TIME()

! Setup mechanisms for own memory management, different messaging services and CPU-time measurements
 
CALL SCARC_SETUP_MEMORY_MANAGER
CALL SCARC_SETUP_MESSAGE_SERVICES
CALL SCARC_SETUP_TIMINGS

! Parse all ScaRC parameters which have been read in read.f90
 
CALL SCARC_PARSE_INPUT                      ; IF (STOP_STATUS==SETUP_STOP) RETURN

! Setup different basic components of ScaRC solver
 
CALL SCARC_SETUP_LEVELS                     ; IF (STOP_STATUS==SETUP_STOP) RETURN
CALL SCARC_SETUP_TYPES                      ; IF (STOP_STATUS==SETUP_STOP) RETURN
CALL SCARC_SETUP_GRIDS                      ; IF (STOP_STATUS==SETUP_STOP) RETURN
CALL SCARC_SETUP_GLOBALS                    ; IF (STOP_STATUS==SETUP_STOP) RETURN
CALL SCARC_SETUP_NEIGHBORS                  ; IF (STOP_STATUS==SETUP_STOP) RETURN
CALL SCARC_SETUP_FACES                      ; IF (STOP_STATUS==SETUP_STOP) RETURN
CALL SCARC_SETUP_SUBDIVISION                ; IF (STOP_STATUS==SETUP_STOP) RETURN

! Setup wall information according to specified discretization type/method
 
IF (HAS_MULTIPLE_GRIDS) THEN
   CALL SCARC_SET_GRID_TYPE (NSCARC_GRID_STRUCTURED)
   CALL SCARC_SETUP_WALLS                   ; IF (STOP_STATUS==SETUP_STOP) RETURN
   CALL SCARC_SET_GRID_TYPE (NSCARC_GRID_UNSTRUCTURED)
   CALL SCARC_SETUP_WALLS                   ; IF (STOP_STATUS==SETUP_STOP) RETURN
ELSE
   CALL SCARC_SET_GRID_TYPE (TYPE_GRID)
   CALL SCARC_SETUP_WALLS                   ; IF (STOP_STATUS==SETUP_STOP) RETURN
ENDIF

! Setup information for data exchanges and matrix systems
 
CALL SCARC_SETUP_EXCHANGES                  ; IF (STOP_STATUS==SETUP_STOP) RETURN
CALL SCARC_SETUP_SYSTEMS                    ; IF (STOP_STATUS==SETUP_STOP) RETURN

! Setup information for algebraic multigrid if requested

#ifdef WITH_SCARC_AMG
IF (HAS_AMG_LEVELS) &
   CALL SCARC_SETUP_ALGEBRAIC_MULTIGRID     ; IF (STOP_STATUS==SETUP_STOP) RETURN
#endif

! Setup environments for the requested methods 

SELECT_METHOD: SELECT CASE(TYPE_METHOD)
   CASE (NSCARC_METHOD_KRYLOV)
      CALL SCARC_SETUP_KRYLOV_ENVIRONMENT(N_STACK_TOTAL)
    CASE (NSCARC_METHOD_MULTIGRID)
       CALL SCARC_SETUP_MULTIGRID_ENVIRONMENT(N_STACK_TOTAL)
   CASE (NSCARC_METHOD_MGM)
       CALL SCARC_SETUP_MGM_ENVIRONMENT(N_STACK_TOTAL)
#ifdef WITH_MKL
   CASE (NSCARC_METHOD_LU)
       CALL SCARC_SETUP_MKL_ENVIRONMENT(N_STACK_TOTAL)
#endif
END SELECT SELECT_METHOD

! Setup vectors for the requested methods 

CALL SCARC_SETUP_VECTORS                    ; IF (STOP_STATUS==SETUP_STOP) RETURN

! Perform some error statistics for pressure if requested
 
#ifdef WITH_SCARC_POSTPROCESSING
CALL SCARC_SETUP_PRESSURE                   ; IF (STOP_STATUS==SETUP_STOP) RETURN
#endif

CPU(MYID)%SETUP   = CPU(MYID)%SETUP   + CURRENT_TIME() - TNOW
CPU(MYID)%OVERALL = CPU(MYID)%OVERALL + CURRENT_TIME() - TNOW
END SUBROUTINE SCARC_SETUP



! ------------------------------------------------------------------------------------------------------
!> \brief Call of requested ScaRC solver 
! ------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SOLVER(DT_CURRENT)
USE SCARC_ITERATION_MANAGER
REAL (EB), INTENT(IN) :: DT_CURRENT
REAL (EB) :: TNOW

TNOW = CURRENT_TIME()

CALL SCARC_SET_ITERATION_STATE (DT_CURRENT)

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'STARTING SCARC_SOLVER METHOD : TYPE_METHOD, TPI = ', TYPE_METHOD, TOTAL_PRESSURE_ITERATIONS
CALL SCARC_DEBUG_METHOD('DEBUG: STARTING SCARC_SOLVER ',6)
#endif

SELECT_METHOD: SELECT CASE (TYPE_METHOD)

   CASE (NSCARC_METHOD_KRYLOV)
      CALL SCARC_METHOD_KRYLOV (NSCARC_STACK_ROOT, NSCARC_STACK_ZERO, NSCARC_RHS_INHOMOGENEOUS, NLEVEL_MIN)
   
   CASE (NSCARC_METHOD_MULTIGRID)
      CALL SCARC_METHOD_MULTIGRID(NSCARC_STACK_ROOT, NSCARC_STACK_ZERO, NSCARC_RHS_INHOMOGENEOUS, NLEVEL_MIN)
   
   CASE (NSCARC_METHOD_MGM)
      CALL SCARC_METHOD_MGM(NSCARC_STACK_ROOT)

#ifdef WITH_MKL
   CASE (NSCARC_METHOD_LU)
      CALL SCARC_METHOD_MKL(NSCARC_STACK_ROOT, NSCARC_STACK_ZERO, NLEVEL_MIN)
#endif
   
END SELECT SELECT_METHOD

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'LEAVING SCARC_SOLVER METHOD : TYPE_METHOD, TPI = ', TYPE_METHOD, TOTAL_PRESSURE_ITERATIONS
CALL SCARC_DEBUG_METHOD('DEBUG: LEAVING SCARC_SOLVER  ',6)
#endif

IF (STOP_STATUS==SETUP_STOP) RETURN

T_USED(5)=T_USED(5)+CURRENT_TIME()-TNOW
CPU(MYID)%SOLVER =CPU(MYID)%SOLVER+CURRENT_TIME()-TNOW
CPU(MYID)%OVERALL=CPU(MYID)%OVERALL+CURRENT_TIME()-TNOW
END SUBROUTINE SCARC_SOLVER

END MODULE SCRC

