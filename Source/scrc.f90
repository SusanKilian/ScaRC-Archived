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
#define WITH_SCARC_DEBUG
#define WITH_SCARC_VERBOSE
#define WITH_SCARC_AMG
#undef WITH_SCARC_POSTPROCESSING

#include "scrc/module_scarc_constants.f90"
#include "scrc/module_scarc_variables.f90"
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
#include "scrc/module_scarc_pointers.f90"
#include "scrc/module_scarc_pointer_routines.f90"
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

! 
! ---------- Type declarations
! 
TYPE (SCARC_TYPE)       , SAVE, DIMENSION(:), ALLOCATABLE, TARGET :: SCARC       !< Main ScaRC data structure
TYPE (SCARC_STACK_TYPE) , SAVE, DIMENSION(:), ALLOCATABLE, TARGET :: STACK       !< Stack of consecutive solvers

TYPE (SCARC_MEMORY_TYPE), SAVE, TARGET :: MEMORY                   !< Memory administration for ScaRC arrays

TYPE (SCARC_SOLVER_TYPE), SAVE, TARGET :: MAIN_CG                  !< Solver structure for Krylov main solver 
TYPE (SCARC_SOLVER_TYPE), SAVE, TARGET :: MAIN_CG_STRUCTURED       !< Solver structure for structured Krylov main solver 
TYPE (SCARC_SOLVER_TYPE), SAVE, TARGET :: MAIN_CG_UNSTRUCTURED     !< Solver structure for unstructured Krylov main solver 
TYPE (SCARC_SOLVER_TYPE), SAVE, TARGET :: MAIN_GMG                 !< Solver structure for Multigrid main solver 

TYPE (SCARC_SOLVER_TYPE), SAVE, TARGET :: COARSE_KRYLOV            !< Solver structure for Krylov coarse grid solver 

TYPE (SCARC_SOLVER_TYPE), SAVE, TARGET :: PRECON_JAC          !< Solver structure for Jacobi preconditioner 
TYPE (SCARC_SOLVER_TYPE), SAVE, TARGET :: PRECON_SSOR         !< Solver structure for SSOR preconditioner 
TYPE (SCARC_SOLVER_TYPE), SAVE, TARGET :: PRECON_LU           !< Solver structure for ILU preconditioner 
TYPE (SCARC_SOLVER_TYPE), SAVE, TARGET :: PRECON_ILU          !< Solver structure for ILU preconditioner 
TYPE (SCARC_SOLVER_TYPE), SAVE, TARGET :: PRECON_FFT          !< Solver structure for FFT preconditioner 
TYPE (SCARC_SOLVER_TYPE), SAVE, TARGET :: PRECON_FFTO         !< Solver structure for FFTO preconditioner (including overlap)
TYPE (SCARC_SOLVER_TYPE), SAVE, TARGET :: PRECON_MG           !< Solver structure for Multigrid preconditioner 
TYPE (SCARC_SOLVER_TYPE), SAVE, TARGET :: PRECON_MJAC         !< Solver structure for Jacobi preconditioner (matrix version)
TYPE (SCARC_SOLVER_TYPE), SAVE, TARGET :: PRECON_MGS          !< Solver structure for Gauss-Seidel preconditioner (matrix version)
TYPE (SCARC_SOLVER_TYPE), SAVE, TARGET :: PRECON_MSGS         !< Solver structure for Sym. Gauss-Seidel preconditioner (matrix vs.)
TYPE (SCARC_SOLVER_TYPE), SAVE, TARGET :: PRECON_MSOR         !< Solver structure for SOR preconditioner (matrix version)
TYPE (SCARC_SOLVER_TYPE), SAVE, TARGET :: PRECON_MSSOR        !< Solver structure for SSOR preconditioner (matrix version)

TYPE (SCARC_SOLVER_TYPE), SAVE, TARGET :: SMOOTH_JAC          !< Solver structure for Jacobi smoother 
TYPE (SCARC_SOLVER_TYPE), SAVE, TARGET :: SMOOTH_SSOR         !< Solver structure for SSOR smoother 
TYPE (SCARC_SOLVER_TYPE), SAVE, TARGET :: SMOOTH_ILU          !< Solver structure for ILU smoother 
TYPE (SCARC_SOLVER_TYPE), SAVE, TARGET :: SMOOTH_FFT          !< Solver structure for FFT smoother 
TYPE (SCARC_SOLVER_TYPE), SAVE, TARGET :: SMOOTH_FFTO         !< Solver structure for FFTO smoother (including overlap)
TYPE (SCARC_SOLVER_TYPE), SAVE, TARGET :: SMOOTH_MJAC         !< Solver structure for Jacobi smoother (matrix version)
TYPE (SCARC_SOLVER_TYPE), SAVE, TARGET :: SMOOTH_MGS          !< Solver structure for Gauss-Seidel smoother (matrix version)
TYPE (SCARC_SOLVER_TYPE), SAVE, TARGET :: SMOOTH_MSGS         !< Solver structure for Sym. Gauss-Seidel smoother (matrix vs.)
TYPE (SCARC_SOLVER_TYPE), SAVE, TARGET :: SMOOTH_MSOR         !< Solver structure for SOR smoother (matrix version)
TYPE (SCARC_SOLVER_TYPE), SAVE, TARGET :: SMOOTH_MSSOR        !< Solver structure for SSOR smoother (matrix version)
TYPE (SCARC_MESSAGE_TYPE), SAVE, TARGET :: MSG
TYPE (SCARC_CPU_TYPE), SAVE, DIMENSION(:), ALLOCATABLE, TARGET :: CPU

#ifdef WITH_MKL
TYPE (SCARC_SOLVER_TYPE), SAVE, TARGET :: MAIN_LU             !< Solver structure for LU-decomposition main solver 
TYPE (SCARC_SOLVER_TYPE), SAVE, TARGET :: COARSE_CLUSTER      !< Solver structure for CLUSTER_SPARSE_SOLVER coarse grid solver 
TYPE (SCARC_SOLVER_TYPE), SAVE, TARGET :: COARSE_PARDISO      !< Solver structure for PARDISO coarse grid solver
TYPE (SCARC_SOLVER_TYPE), SAVE, TARGET :: PRECON_MKL          !< Solver structure for MKL preconditioner 
TYPE (SCARC_SOLVER_TYPE), SAVE, TARGET :: SMOOTH_MKL          !< Solver structure for MKL smoother 
#endif

TYPE (SCARC_SUBDIVISION_TYPE), SAVE, TARGET :: SUBDIVISION    !< Structure to keep information about subdivision

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

SELECT_METHOD: SELECT CASE (TYPE_METHOD)

   CASE (NSCARC_METHOD_KRYLOV)
      CALL SCARC_METHOD_KRYLOV (NSCARC_STACK_ROOT, NSCARC_STACK_ZERO, NSCARC_RHS_INHOMOGENEOUS, NLEVEL_MIN)
   
   CASE (NSCARC_METHOD_MULTIGRID)
      CALL SCARC_METHOD_MULTIGRID(NSCARC_STACK_ROOT, NSCARC_STACK_ZERO, NSCARC_RHS_INHOMOGENEOUS, NLEVEL_MIN)
   
   CASE (NSCARC_METHOD_MGM)
      CALL SCARC_METHOD_MGM()

#ifdef WITH_MKL
   CASE (NSCARC_METHOD_LU)
      CALL SCARC_METHOD_MKL(NSCARC_STACK_ROOT, NSCARC_STACK_ZERO, NLEVEL_MIN)
#endif
   
END SELECT SELECT_METHOD

IF (STOP_STATUS==SETUP_STOP) RETURN

T_USED(5)=T_USED(5)+CURRENT_TIME()-TNOW
CPU(MYID)%SOLVER =CPU(MYID)%SOLVER+CURRENT_TIME()-TNOW
CPU(MYID)%OVERALL=CPU(MYID)%OVERALL+CURRENT_TIME()-TNOW
END SUBROUTINE SCARC_SOLVER

END MODULE SCRC

