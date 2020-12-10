! ================================================================================================================
!> \brief Scalable Recursive Clustering (ScaRC): Collection of alternative solvers for the FDS pressure equation
! ================================================================================================================
! Use of different directives possible
!  - WITH_MKL                    : use MKL routines PARDISO, CLUSTER_SPARSE_SOLVER, DDOT, DAXPY, DAXPBY, DCOPY, DSCAL
!  - WITH_SCARC_VERBOSE          : print more detailed information about ScaRC iterations and workspace allocation
!  - WITH_SCARC_DEBUG            : print detaild debugging info (only for developing purposes)
!  - WITH_SCARC_POSTPROCESSING   : dump environment for separate ScaRC postprocessing program
! ================================================================================================================
!#undef WITH_MKL
#define WITH_SCARC_DEBUG
#define WITH_SCARC_VERBOSE
#undef WITH_SCARC_POSTPROCESSING

#include "scrc/module_scarc_constants.f90"
#include "scrc/module_scarc_types.f90"
#include "scrc/module_scarc_variables.f90"
#include "scrc/module_scarc_pointers.f90"
#include "scrc/module_scarc_pointer_routines.f90"
#include "scrc/module_scarc_message_services.f90"
#include "scrc/module_scarc_iteration_manager.f90"
#include "scrc/module_scarc_error_manager.f90"
#include "scrc/module_scarc_time_measurement.f90"
#include "scrc/module_scarc_utilities.f90"
#include "scrc/module_scarc_memory_manager.f90"
#include "scrc/module_scarc_mpi.f90"
#include "scrc/module_scarc_stack_operations.f90"
#include "scrc/module_scarc_linear_algebra.f90"
#include "scrc/module_scarc_relaxation_method.f90"
#include "scrc/module_scarc_krylov_method.f90"
#include "scrc/module_scarc_multigrid_method.f90"
#include "scrc/module_scarc_mgm_method.f90"

MODULE SCRC

USE PRECISION_PARAMETERS
USE GLOBAL_CONSTANTS
USE MESH_VARIABLES
USE MEMORY_FUNCTIONS, ONLY: CHKMEMERR
USE COMP_FUNCTIONS, ONLY: CURRENT_TIME, GET_FILE_NUMBER, SHUTDOWN
USE MPI
USE SCARC_CONSTANTS
USE SCARC_TYPES
USE SCARC_VARIABLES
USE SCARC_ITERATION_MANAGER
USE SCARC_MESSAGE_SERVICES
USE SCARC_ERROR_MANAGER
USE SCARC_TIME_MEASUREMENT
USE SCARC_MEMORY_MANAGER
USE SCARC_UTILITIES
USE SCARC_MPI
USE SCARC_STACK_OPERATIONS
USE SCARC_LINEAR_ALGEBRA
USE SCARC_RELAXATION_METHOD
USE SCARC_KRYLOV_METHOD
USE SCARC_MULTIGRID_METHOD

#ifdef WITH_MKL
USE MKL_PARDISO
USE MKL_CLUSTER_SPARSE_SOLVER
#endif

IMPLICIT NONE

PUBLIC :: SCARC_SETUP, SCARC_SOLVER

CONTAINS

! ------------------------------------------------------------------------------------------------
!> \brief Initialize ScaRC structures based on SCARC-input parameters from &PRES namelist
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP
REAL(EB) :: TNOW

TNOW = CURRENT_TIME()

! Setup mechanisms for own memory management, different messaging services and CPU-time measurements
 
CALL SCARC_SETUP_MEMORY_MANAGER
CALL SCARC_SETUP_MESSAGE_SERVICES
CALL SCARC_SETUP_TIME_MEASUREMENTS

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
   CALL SCARC_SETUP_GRID_TYPE (NSCARC_GRID_STRUCTURED)
   CALL SCARC_SETUP_WALLS                   ; IF (STOP_STATUS==SETUP_STOP) RETURN
   CALL SCARC_SETUP_GRID_TYPE (NSCARC_GRID_UNSTRUCTURED)
   CALL SCARC_SETUP_WALLS                   ; IF (STOP_STATUS==SETUP_STOP) RETURN
ELSE
   CALL SCARC_SETUP_GRID_TYPE (TYPE_GRID)
   CALL SCARC_SETUP_WALLS                   ; IF (STOP_STATUS==SETUP_STOP) RETURN
ENDIF

! Setup information for data exchanges, matrix systems, aggregation zones, used methods and vectors
 
CALL SCARC_SETUP_EXCHANGES                  ; IF (STOP_STATUS==SETUP_STOP) RETURN
CALL SCARC_SETUP_SYSTEMS                    ; IF (STOP_STATUS==SETUP_STOP) RETURN

IF (HAS_AMG_LEVELS) &
   CALL SCARC_SETUP_ALGEBRAIC_MULTIGRID     ; IF (STOP_STATUS==SETUP_STOP) RETURN

CALL SCARC_SETUP_METHODS                    ; IF (STOP_STATUS==SETUP_STOP) RETURN
CALL SCARC_SETUP_VECTORS                    ; IF (STOP_STATUS==SETUP_STOP) RETURN

! Perform some error statistics for pressure if corresponding flag is set
 
#ifdef WITH_SCARC_POSTPROCESSING
CALL SCARC_SETUP_PRESSURE                   ; IF (STOP_STATUS==SETUP_STOP) RETURN
#endif

CPU(MYID)%SETUP   = CPU(MYID)%SETUP   + CURRENT_TIME() - TNOW
CPU(MYID)%OVERALL = CPU(MYID)%OVERALL + CURRENT_TIME() - TNOW
END SUBROUTINE SCARC_SETUP


! ----------------------------------------------------------------------------------------------------
!> \brief Determine types of input parameters
! ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_PARSE_INPUT

ITERATE_PRESSURE = .TRUE.  ! Although there is no need to do pressure iterations to drive down 
                           ! velocity error leave it .TRUE. to write out velocity error diagnostics

 
! ------------- Set type of discretization
 
SELECT CASE (TRIM(PRES_METHOD))
   CASE ('SCARC')
      TYPE_GRID = NSCARC_GRID_STRUCTURED
      IS_STRUCTURED   = .TRUE.
      IS_UNSTRUCTURED = .FALSE.
   CASE ('USCARC')
      TYPE_GRID = NSCARC_GRID_UNSTRUCTURED
      IS_STRUCTURED   = .FALSE.
      IS_UNSTRUCTURED = .TRUE.
   CASE DEFAULT
      CALL SCARC_SHUTDOWN(NSCARC_ERROR_PARSE_INPUT, SCARC_GRID, NSCARC_NONE)
END SELECT

 
! ------------ Set type of matrix storage (COMPACT/BANDWISE)
 
SELECT CASE (TRIM(SCARC_MATRIX))
   CASE ('COMPACT')
      TYPE_MATRIX = NSCARC_MATRIX_COMPACT
   CASE ('BANDWISE')
      TYPE_MATRIX = NSCARC_MATRIX_BANDWISE
   CASE DEFAULT
      CALL SCARC_SHUTDOWN(NSCARC_ERROR_PARSE_INPUT, SCARC_MATRIX, NSCARC_NONE)
END SELECT

 
! ------------ Set type of matrix stencil (CONSTANT/VARIABLE)
 
SELECT CASE (TRIM(SCARC_STENCIL))
   CASE ('CONSTANT')
      TYPE_STENCIL = NSCARC_STENCIL_CONSTANT
   CASE ('VARIABLE')
      TYPE_STENCIL = NSCARC_STENCIL_VARIABLE
   CASE DEFAULT
      CALL SCARC_SHUTDOWN(NSCARC_ERROR_PARSE_INPUT, SCARC_STENCIL, NSCARC_NONE)
END SELECT

 
! ------------ Set type of global solver
 
SELECT CASE (TRIM(SCARC_METHOD))

   ! ------------------------- Global Krylov solver ----------------------------------
   CASE ('KRYLOV')

      TYPE_METHOD = NSCARC_METHOD_KRYLOV

      ! Set type of two-level method
      SELECT CASE (TRIM(SCARC_TWOLEVEL))
         CASE ('NONE')
            TYPE_TWOLEVEL = NSCARC_TWOLEVEL_NONE
         CASE ('ADDITIVE')
            TYPE_TWOLEVEL = NSCARC_TWOLEVEL_ADD
         CASE ('MULTIPLICATIVE')
            TYPE_TWOLEVEL = NSCARC_TWOLEVEL_MUL
         CASE ('MULTIPLICATIVE2')
            TYPE_TWOLEVEL = NSCARC_TWOLEVEL_MUL2
         CASE ('COARSE')
            TYPE_TWOLEVEL = NSCARC_TWOLEVEL_COARSE
         CASE ('MACRO')
            TYPE_TWOLEVEL = NSCARC_TWOLEVEL_MACRO
         CASE DEFAULT
            CALL SCARC_SHUTDOWN(NSCARC_ERROR_PARSE_INPUT, SCARC_TWOLEVEL, NSCARC_NONE)
      END SELECT

      ! Set type of interpolation for two-level Krylov method
      SELECT CASE (TRIM(SCARC_KRYLOV_INTERPOL))
         CASE ('NONE')
            TYPE_INTERPOL = NSCARC_UNDEF_INT
         CASE ('CONSTANT')
            TYPE_INTERPOL = NSCARC_INTERPOL_CONSTANT
         CASE ('BILINEAR')
            TYPE_INTERPOL = NSCARC_INTERPOL_BILINEAR
         CASE ('BILINEAR2')
            TYPE_INTERPOL = NSCARC_INTERPOL_BILINEAR2
         CASE ('AMG')
            TYPE_INTERPOL = NSCARC_INTERPOL_AMG
         CASE DEFAULT
            CALL SCARC_SHUTDOWN(NSCARC_ERROR_PARSE_INPUT, SCARC_KRYLOV_INTERPOL, NSCARC_NONE)
      END SELECT

      ! Set type of preconditioner (JACOBI/SSOR/MGS/MSGS/MSOR/MSSOR/ILU/LU/FFT/GMG/PARDISO/CLUSTER)
      SELECT CASE (TRIM(SCARC_PRECON))
         CASE ('JACOBI')                                    ! Jacobi preconditioner
            TYPE_PRECON = NSCARC_RELAX_JAC
         CASE ('SSOR')                                      ! Symmetric SOR preconditioner
            TYPE_PRECON = NSCARC_RELAX_SSOR
         CASE ('MJAC')                                      ! Jacobi preconditioner in matrix form
            TYPE_PRECON = NSCARC_RELAX_MJAC
         CASE ('MGS')                                       ! Gauss-Seidel preconditioner in matrix form
            TYPE_PRECON = NSCARC_RELAX_MGS
         CASE ('MSGS')                                      ! Symmetric Gauss-Seidel preconditioner in matrix form
            TYPE_PRECON = NSCARC_RELAX_MSGS
         CASE ('MSOR')                                      ! SOR preconditioner in matrix form
            TYPE_PRECON = NSCARC_RELAX_MSOR
         CASE ('MSSOR')                                     ! Symmetric SOR preconditioner in matrix form
            TYPE_PRECON = NSCARC_RELAX_MSSOR
         CASE ('ILU')                                       ! ILU preconditioner
            TYPE_PRECON = NSCARC_RELAX_ILU
         CASE ('LU')                                        ! LU preconditioner
            TYPE_PRECON = NSCARC_RELAX_LU
         CASE ('MULTIGRID')                                 ! Multigrid preconditioner
            TYPE_PRECON = NSCARC_RELAX_MULTIGRID
            SELECT CASE (TRIM(SCARC_SMOOTH))
               CASE ('JACOBI')
                  TYPE_SMOOTH = NSCARC_RELAX_JAC
               CASE ('SSOR')
                  TYPE_SMOOTH = NSCARC_RELAX_SSOR
               CASE ('MJAC')
                  TYPE_SMOOTH = NSCARC_RELAX_MJAC
               CASE ('MGS')
                  TYPE_SMOOTH = NSCARC_RELAX_MGS
               CASE ('MSGS')
                  TYPE_SMOOTH = NSCARC_RELAX_MSGS
               CASE ('MSOR')
                  TYPE_SMOOTH = NSCARC_RELAX_MSOR
               CASE ('MSSOR')
                  TYPE_SMOOTH = NSCARC_RELAX_MSSOR
               CASE ('ILU')
                  TYPE_SMOOTH = NSCARC_RELAX_ILU
               CASE ('FFT')
                  IF (IS_UNSTRUCTURED) CALL SCARC_SHUTDOWN(NSCARC_ERROR_FFT_GRID, SCARC_NONE, NSCARC_NONE)
                  TYPE_SMOOTH = NSCARC_RELAX_FFT
               CASE ('FFTO')
                  IF (IS_UNSTRUCTURED) CALL SCARC_SHUTDOWN(NSCARC_ERROR_FFT_GRID, SCARC_NONE, NSCARC_NONE)
                  IF (NMESHES == 1) THEN
                     TYPE_SMOOTH = NSCARC_RELAX_FFT
                  ELSE
                     TYPE_SMOOTH = NSCARC_RELAX_FFTO
                  ENDIF
               CASE ('PARDISO')
#ifdef WITH_MKL
                  TYPE_SMOOTH = NSCARC_RELAX_MKL
#else
                  CALL SCARC_SHUTDOWN(NSCARC_ERROR_MKL_PARDISO, SCARC_NONE, NSCARC_NONE)
#endif

               CASE ('CLUSTER')
#ifdef WITH_MKL
                  TYPE_SMOOTH = NSCARC_RELAX_MKL
#else
                  CALL SCARC_SHUTDOWN(NSCARC_ERROR_MKL_STORAGE, SCARC_NONE, NSCARC_NONE)
#endif
            END SELECT
         CASE ('FFT')                                                ! FFT preconditioner
            IF (IS_UNSTRUCTURED) CALL SCARC_SHUTDOWN(NSCARC_ERROR_FFT_GRID, SCARC_NONE, NSCARC_NONE)
            TYPE_PRECON = NSCARC_RELAX_FFT
         CASE ('FFTO')                                               ! FFT with overlap preconditioner
            IF (IS_UNSTRUCTURED) CALL SCARC_SHUTDOWN(NSCARC_ERROR_FFT_GRID, SCARC_NONE, NSCARC_NONE)
            IF (NMESHES == 1) THEN
               TYPE_PRECON = NSCARC_RELAX_FFT
            ELSE
               TYPE_PRECON = NSCARC_RELAX_FFTO
            ENDIF
         CASE ('PARDISO')                                            ! LU preconditioner based on MKL-PARDISO
#ifdef WITH_MKL
            TYPE_PRECON   = NSCARC_RELAX_MKL
            TYPE_MKL(0)   = NSCARC_MKL_LOCAL
            TYPE_SCOPE(1) = NSCARC_SCOPE_LOCAL
#else
            CALL SCARC_SHUTDOWN(NSCARC_ERROR_MKL_PARDISO, SCARC_NONE, NSCARC_NONE)
#endif
         CASE ('CLUSTER')                            !  LU-preconditioner based on MKL Cluster_Sparse_Solver
#ifdef WITH_MKL
            TYPE_PRECON   = NSCARC_RELAX_MKL
            TYPE_MKL(0)   = NSCARC_MKL_GLOBAL
            TYPE_SCOPE(1) = NSCARC_SCOPE_GLOBAL
#else
            CALL SCARC_SHUTDOWN(NSCARC_ERROR_MKL_CLUSTER, SCARC_NONE, NSCARC_NONE)
#endif
         CASE DEFAULT
            CALL SCARC_SHUTDOWN(NSCARC_ERROR_PARSE_INPUT, SCARC_PRECON, NSCARC_NONE)
      END SELECT

      ! set type scope for preconditioner (GLOBAL/LOCAL)
      SELECT CASE (TRIM(SCARC_PRECON_SCOPE))
         CASE ('GLOBAL')
            TYPE_SCOPE(1) = NSCARC_SCOPE_GLOBAL
         CASE ('LOCAL')
            TYPE_SCOPE(1) = NSCARC_SCOPE_LOCAL
         CASE DEFAULT
            CALL SCARC_SHUTDOWN(NSCARC_ERROR_PARSE_INPUT, SCARC_PRECON_SCOPE, NSCARC_NONE)
      END SELECT

   ! ------------------------- Global geometric multigrid solver -------------------------------
   CASE ('MULTIGRID')

      TYPE_METHOD = NSCARC_METHOD_MULTIGRID

      ! Set type of multigrid method (GEOMETRIC/ALGEBRAIC)
      SELECT CASE (TRIM(SCARC_MULTIGRID))
         CASE ('GEOMETRIC')
            TYPE_MULTIGRID = NSCARC_MULTIGRID_GEOMETRIC
            TYPE_COARSENING = NSCARC_COARSENING_GMG         ! GMG-default, may be overwritten by SCARC_COARSENING 
         CASE ('ALGEBRAIC')
            TYPE_MULTIGRID = NSCARC_MULTIGRID_ALGEBRAIC
            TYPE_COARSENING = NSCARC_COARSENING_CUBIC       ! AMG-default, may be overwritten by SCARC_COARSENING 
         CASE DEFAULT
            CALL SCARC_SHUTDOWN(NSCARC_ERROR_PARSE_INPUT, SCARC_MULTIGRID, NSCARC_NONE)
      END SELECT

      ! Set type of smoother (JACOBI/SGS/SSOR/MSSOR/ILU/PARDISO/CLUSTER)
      SELECT CASE (TRIM(SCARC_SMOOTH))                      ! use same parameters as for preconditioner
         CASE ('JACOBI')                                    ! Jacobi preconditioner
            TYPE_SMOOTH = NSCARC_RELAX_JAC
         CASE ('SSOR')                                      ! SSOR preconditioner
            TYPE_SMOOTH = NSCARC_RELAX_SSOR
         CASE ('MJAC')                                      ! Jacobi preconditioner in matrix form
            TYPE_SMOOTH = NSCARC_RELAX_MJAC
         CASE ('MGS')                                       ! Gauss-Seidel preconditioner in matrix form
            TYPE_SMOOTH = NSCARC_RELAX_MGS
         CASE ('MSGS')                                      ! Symmetric Gauss-Seidel preconditioner in matrix form
            TYPE_SMOOTH = NSCARC_RELAX_MSGS
         CASE ('MSOR')                                      ! SOR preconditioner in matrix form
            TYPE_SMOOTH = NSCARC_RELAX_MSOR
         CASE ('MSSOR')                                     ! Symmetric SOR preconditioner in matrix form
            TYPE_SMOOTH = NSCARC_RELAX_MSSOR
         CASE ('ILU')
            TYPE_SMOOTH = NSCARC_RELAX_ILU
         CASE ('FFT')
            IF (IS_UNSTRUCTURED) CALL SCARC_SHUTDOWN(NSCARC_ERROR_FFT_GRID, SCARC_NONE, NSCARC_NONE)
            TYPE_SMOOTH = NSCARC_RELAX_FFT
         CASE ('FFTO')
            IF (IS_UNSTRUCTURED) CALL SCARC_SHUTDOWN(NSCARC_ERROR_FFT_GRID, SCARC_NONE, NSCARC_NONE)
            IF (NMESHES == 1) THEN
               TYPE_SMOOTH = NSCARC_RELAX_FFT
            ELSE
               TYPE_SMOOTH = NSCARC_RELAX_FFTO
            ENDIF
         CASE ('PARDISO')
#ifdef WITH_MKL
            TYPE_SMOOTH = NSCARC_RELAX_MKL
#else
            CALL SCARC_SHUTDOWN(NSCARC_ERROR_MKL_PARDISO, SCARC_NONE, NSCARC_NONE)
#endif
         CASE ('CLUSTER')
#ifdef WITH_MKL
            TYPE_SMOOTH = NSCARC_RELAX_MKL
#else
            CALL SCARC_SHUTDOWN(NSCARC_ERROR_MKL_CLUSTER, SCARC_NONE, NSCARC_NONE)
#endif
         CASE DEFAULT
            CALL SCARC_SHUTDOWN(NSCARC_ERROR_PARSE_INPUT, SCARC_SMOOTH, NSCARC_NONE)
      END SELECT

      ! set type scope for smoother (GLOBAL/LOCAL)
      SELECT CASE (TRIM(SCARC_SMOOTH_SCOPE))
         CASE ('GLOBAL')
            TYPE_SCOPE(2) = NSCARC_SCOPE_GLOBAL
         CASE ('LOCAL')
            TYPE_SCOPE(2) = NSCARC_SCOPE_LOCAL
         CASE DEFAULT
            CALL SCARC_SHUTDOWN(NSCARC_ERROR_PARSE_INPUT, SCARC_SMOOTH_SCOPE, NSCARC_NONE)
      END SELECT

   ! ------------------------- Global LU-decomposition solver -------------------------------
#ifdef WITH_MKL
   CASE ('MKL')

      TYPE_METHOD  = NSCARC_METHOD_LU

      ! Set type of MKL method (global/local)
      SELECT CASE (TRIM(SCARC_MKL))                  
         CASE ('GLOBAL')
#ifdef WITH_MKL
            TYPE_MKL(0)   = NSCARC_MKL_GLOBAL
#else
            CALL SCARC_SHUTDOWN(NSCARC_ERROR_MKL_CLUSTER, SCARC_NONE, NSCARC_NONE)
#endif
         CASE ('LOCAL')
#ifdef WITH_MKL
            TYPE_MKL(0)   = NSCARC_MKL_LOCAL
#else
            CALL SCARC_SHUTDOWN(NSCARC_ERROR_MKL_PARDISO, SCARC_NONE, NSCARC_NONE)
#endif
         CASE DEFAULT
            CALL SCARC_SHUTDOWN(NSCARC_ERROR_PARSE_INPUT, SCARC_MKL, NSCARC_NONE)
      END SELECT
#endif

   ! ------------------------- McKenny-Greengard-Mayo solver -------------------------
   CASE ('MGM')

      !  Both structured and unstructured discretization are required
      !  The second pass of this method is purely locally based (the Laplace solutions) 

      HAS_MULTIPLE_GRIDS = .TRUE.

      TYPE_METHOD   = NSCARC_METHOD_MGM
      TYPE_SCOPE(1) = NSCARC_SCOPE_LOCAL

      ! set type of MGM interface BCs of Laplace problems
      SELECT CASE (TRIM(SCARC_MGM_BC))
         CASE ('MEAN')
            TYPE_MGM_BC = NSCARC_MGM_BC_MEAN
         CASE ('TRUE')
            TYPE_MGM_BC = NSCARC_MGM_BC_TRUE
         CASE ('EXTRAPOLATION')
            TYPE_MGM_BC = NSCARC_MGM_BC_EXPOL
         CASE ('TAYLOR')
            TYPE_MGM_BC = NSCARC_MGM_BC_TAYLOR
         CASE DEFAULT
            CALL SCARC_SHUTDOWN(NSCARC_ERROR_PARSE_INPUT, SCARC_MGM_BC, NSCARC_NONE)
      END SELECT

      ! set type of MGM interpolation for interface BCs of Laplace problems
      SELECT CASE (TRIM(SCARC_MGM_INTERPOL))
         CASE ('LINEAR')
            TYPE_MGM_INTERPOL = NSCARC_MGM_INTERPOL_LINEAR
         CASE ('SQUARE')
            TYPE_MGM_INTERPOL = NSCARC_MGM_INTERPOL_SQUARE
         CASE DEFAULT
            CALL SCARC_SHUTDOWN(NSCARC_ERROR_PARSE_INPUT, SCARC_MGM_INTERPOL, NSCARC_NONE)
      END SELECT

   CASE DEFAULT
      CALL SCARC_SHUTDOWN(NSCARC_ERROR_PARSE_INPUT, SCARC_METHOD, NSCARC_NONE)

END SELECT

 
! If a multigrid solver is used (either as main solver or as preconditioner)
! set types for multigrid, coarse grid solver and cycling pattern
 
IF (TYPE_METHOD == NSCARC_METHOD_MULTIGRID .OR. TYPE_PRECON == NSCARC_RELAX_MULTIGRID) THEN

   ! Set type of multigrid (GEOMETRIC/ALGEBRAIC with corresponding coarsening strategy)
   SELECT CASE (TRIM(SCARC_MULTIGRID))
      CASE ('GEOMETRIC')
         TYPE_MULTIGRID = NSCARC_MULTIGRID_GEOMETRIC
      CASE ('ALGEBRAIC')
         TYPE_MULTIGRID = NSCARC_MULTIGRID_ALGEBRAIC
      CASE DEFAULT
         CALL SCARC_SHUTDOWN(NSCARC_ERROR_PARSE_INPUT, SCARC_MULTIGRID, NSCARC_NONE)
   END SELECT

   ! Set type of cycling pattern (F/V/W)
   SELECT CASE (TRIM(SCARC_MULTIGRID_CYCLE))
      CASE ('F')
         TYPE_CYCLING = NSCARC_CYCLING_F
      CASE ('V')
         TYPE_CYCLING = NSCARC_CYCLING_V
      CASE ('W')
         TYPE_CYCLING = NSCARC_CYCLING_W
      CASE ('FMG')
         TYPE_CYCLING = NSCARC_CYCLING_FMG
      CASE DEFAULT
         CALL SCARC_SHUTDOWN(NSCARC_ERROR_PARSE_INPUT, SCARC_MULTIGRID_CYCLE, NSCARC_NONE)
   END SELECT

   ! Set type of interpolation 
   SELECT CASE (TRIM(SCARC_MULTIGRID_INTERPOL))
      CASE ('STANDARD')
         TYPE_INTERPOL = NSCARC_INTERPOL_STANDARD
      CASE ('CONSTANT')
         TYPE_INTERPOL = NSCARC_INTERPOL_CONSTANT
      CASE ('BILINEAR')
         TYPE_INTERPOL = NSCARC_INTERPOL_BILINEAR
      CASE ('BILINEAR2')
         TYPE_INTERPOL = NSCARC_INTERPOL_BILINEAR2
      CASE ('CLASSICAL')
         TYPE_INTERPOL = NSCARC_INTERPOL_CLASSICAL
      CASE ('DIRECT')
         TYPE_INTERPOL = NSCARC_INTERPOL_DIRECT
      CASE ('AMG')
         TYPE_INTERPOL = NSCARC_INTERPOL_AMG
      CASE DEFAULT
         CALL SCARC_SHUTDOWN(NSCARC_ERROR_PARSE_INPUT, SCARC_MULTIGRID_INTERPOL, NSCARC_NONE)
   END SELECT

ENDIF

 
! ------------ Set type of coarsening strategy in case of multi-level methods
 
SELECT CASE (TRIM(SCARC_COARSENING))
   CASE ('AGGREGATED')
      TYPE_COARSENING = NSCARC_COARSENING_AGGREGATED
   CASE ('AGGREGATEDS')
      TYPE_COARSENING = NSCARC_COARSENING_AGGREGATED_S
   CASE ('CUBIC')
      TYPE_COARSENING = NSCARC_COARSENING_CUBIC
   CASE ('GMG')
      TYPE_COARSENING = NSCARC_COARSENING_GMG
   CASE DEFAULT
      CALL SCARC_SHUTDOWN(NSCARC_ERROR_PARSE_INPUT, SCARC_COARSENING, NSCARC_NONE)
END SELECT


! Set type of coarse grid solver
SELECT CASE (TRIM(SCARC_COARSE))
   CASE ('ITERATIVE')
      TYPE_COARSE = NSCARC_COARSE_ITERATIVE
   CASE ('DIRECT')
#ifdef WITH_MKL
      TYPE_COARSE   = NSCARC_COARSE_DIRECT
      TYPE_MKL(0)   = NSCARC_MKL_COARSE
#else
      CALL SCARC_SHUTDOWN(NSCARC_ERROR_MKL_CLUSTER, SCARC_NONE, NSCARC_NONE)
#endif
   CASE DEFAULT
      CALL SCARC_SHUTDOWN(NSCARC_ERROR_PARSE_INPUT, SCARC_COARSE, NSCARC_NONE)
END SELECT


 
! Set type of accuracy (ABSOLUTE/RELATIVE)
 
SELECT CASE (TRIM(SCARC_ACCURACY))
   CASE ('ABSOLUTE')
      TYPE_ACCURACY = NSCARC_ACCURACY_ABSOLUTE
   CASE ('RELATIVE')
      TYPE_ACCURACY = NSCARC_ACCURACY_RELATIVE
   CASE DEFAULT
      CALL SCARC_SHUTDOWN(NSCARC_ERROR_PARSE_INPUT, SCARC_ACCURACY, NSCARC_NONE)
END SELECT

 
! Set type of precision for MKL solver (SINGLE/DOUBLE)
 
SELECT CASE (TRIM(SCARC_MKL_PRECISION))
   CASE ('SINGLE')
      TYPE_MKL_PRECISION = NSCARC_PRECISION_SINGLE
   CASE ('DOUBLE')
      TYPE_MKL_PRECISION = NSCARC_PRECISION_DOUBLE
   CASE DEFAULT
      CALL SCARC_SHUTDOWN(NSCARC_ERROR_PARSE_INPUT, SCARC_MKL_PRECISION, NSCARC_NONE)
END SELECT

 
! -------- Define some logical variables - just for notational convenience
 
IS_STRUCTURED   = (TYPE_GRID == NSCARC_GRID_STRUCTURED)
IS_UNSTRUCTURED = (TYPE_GRID == NSCARC_GRID_UNSTRUCTURED)

IS_CG     = (TYPE_METHOD == NSCARC_METHOD_KRYLOV)
IS_CG_MG  = IS_CG .AND. (TYPE_PRECON == NSCARC_RELAX_MULTIGRID) 
IS_CG_GMG = IS_CG .AND. (TYPE_PRECON == NSCARC_RELAX_MULTIGRID) .AND. (TYPE_MULTIGRID == NSCARC_MULTIGRID_GEOMETRIC)
IS_CG_AMG = IS_CG .AND. (TYPE_PRECON == NSCARC_RELAX_MULTIGRID) .AND. (TYPE_MULTIGRID == NSCARC_MULTIGRID_ALGEBRAIC)

IS_MG  = (TYPE_METHOD == NSCARC_METHOD_MULTIGRID)
IS_GMG = IS_MG .AND. (TYPE_MULTIGRID == NSCARC_MULTIGRID_GEOMETRIC)
IS_AMG = IS_MG .AND. (TYPE_MULTIGRID == NSCARC_MULTIGRID_ALGEBRAIC)

IS_FFT = (TYPE_PRECON == NSCARC_RELAX_FFT)  .OR. (TYPE_SMOOTH == NSCARC_RELAX_FFT)
IS_FFTO= (TYPE_PRECON == NSCARC_RELAX_FFTO) .OR. (TYPE_SMOOTH == NSCARC_RELAX_FFTO)
IS_MKL = (TYPE_PRECON >= NSCARC_RELAX_MKL)  .OR. (TYPE_SMOOTH >= NSCARC_RELAX_MKL) 

HAS_TWO_LEVELS      = IS_CG .AND. (TYPE_PRECON /= NSCARC_RELAX_MULTIGRID) .AND. (TYPE_TWOLEVEL > NSCARC_TWOLEVEL_NONE)
HAS_MULTIPLE_LEVELS = IS_MG .OR. IS_CG_MG .OR. HAS_TWO_LEVELS 

IS_CG_ADD    = HAS_TWO_LEVELS .AND. (TYPE_TWOLEVEL == NSCARC_TWOLEVEL_ADD)
IS_CG_MUL    = HAS_TWO_LEVELS .AND. (TYPE_TWOLEVEL == NSCARC_TWOLEVEL_MUL)
IS_CG_MACRO  = HAS_TWO_LEVELS .AND. (TYPE_TWOLEVEL == NSCARC_TWOLEVEL_MACRO)
IS_CG_COARSE = HAS_TWO_LEVELS .AND. (TYPE_TWOLEVEL == NSCARC_TWOLEVEL_COARSE)

HAS_GMG_LEVELS = IS_GMG .OR. IS_CG_GMG .OR. HAS_TWO_LEVELS
HAS_AMG_LEVELS = IS_AMG .OR. IS_CG_AMG 

IS_MGM = TYPE_METHOD == NSCARC_METHOD_MGM


END SUBROUTINE SCARC_PARSE_INPUT


! ------------------------------------------------------------------------------------------------
!> \brief Determine number of grid levels 
! NLEVEL_MIN corresponds to finest grid resolution, NLEVEL_MAX to coarsest resolution
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_LEVELS
#ifdef WITH_MKL
INTEGER :: NL
#endif

SELECT_SCARC_METHOD: SELECT CASE (TYPE_METHOD)

 
   ! ---------- Global data-parallel Krylov method 
 
   CASE (NSCARC_METHOD_KRYLOV)

      SELECT_KRYLOV_PRECON: SELECT CASE (TYPE_PRECON)

#ifdef WITH_MKL
 
         ! Preconditioning by defect correction based on LU-decomposition
         ! If two-level method, also use coarse grid level, otherwise only use single (finest) grid level
         ! Either using a global CLUSTER_SPARSE_SOLVER or local PARDISO solvers from MKL
 
         CASE (NSCARC_RELAX_MKL)

            IF (HAS_TWO_LEVELS) THEN
               CALL SCARC_GET_NUMBER_OF_LEVELS(NSCARC_LEVEL_MULTI)
               IF (TYPE_COARSE == NSCARC_COARSE_DIRECT) TYPE_MKL(NLEVEL_MAX) = NSCARC_MKL_GLOBAL
            ELSE
               CALL SCARC_GET_NUMBER_OF_LEVELS(NSCARC_LEVEL_SINGLE)
            ENDIF

            IF (TYPE_SCOPE(1) == NSCARC_SCOPE_GLOBAL) THEN
               TYPE_MKL(NLEVEL_MIN) = NSCARC_MKL_GLOBAL
            ELSE IF (TYPE_SCOPE(1) == NSCARC_SCOPE_LOCAL) THEN
               TYPE_MKL(NLEVEL_MIN) = NSCARC_MKL_LOCAL
            ENDIF

#endif

         ! Preconditioning by defect correction based on geometric multigrid method,
         ! use specified hierarchy of grid levels
 
         CASE (NSCARC_RELAX_MULTIGRID)
            CALL SCARC_GET_NUMBER_OF_LEVELS(NSCARC_LEVEL_MULTI)
#ifdef WITH_MKL
            IF (TYPE_COARSE == NSCARC_COARSE_DIRECT) TYPE_MKL(NLEVEL_MAX) = NSCARC_MKL_GLOBAL
#endif

 
         ! Preconditioning by defect correction based on local basic iterations (JACOBI/SSOR),
         ! if two-level method, also use coarse grid, otherwise only use single (finest) grid level
 
         CASE DEFAULT
            IF (HAS_TWO_LEVELS) THEN
               CALL SCARC_GET_NUMBER_OF_LEVELS(NSCARC_LEVEL_MULTI)
#ifdef WITH_MKL
               IF (TYPE_COARSE == NSCARC_COARSE_DIRECT) TYPE_MKL(NLEVEL_MAX) = NSCARC_MKL_GLOBAL
#endif
            ELSE
               CALL SCARC_GET_NUMBER_OF_LEVELS(NSCARC_LEVEL_SINGLE)
            ENDIF

      END SELECT SELECT_KRYLOV_PRECON

 
   ! ---------- Global data-parallel Multigrid method 
 
   CASE (NSCARC_METHOD_MULTIGRID)

         ! If not specified by user, determine number of possible grid levels

         CALL SCARC_GET_NUMBER_OF_LEVELS(NSCARC_LEVEL_MULTI)

#ifdef WITH_MKL

         ! In case of smoothing by different MKL solvers, mark levels for the use of MKL,
         ! either by locally acting PARDISO solvers or Globally acting CLUSTER_SPARSE_SOLVER

         IF (TYPE_SMOOTH == NSCARC_RELAX_MKL) THEN

            IF (TYPE_SCOPE(2) == NSCARC_SCOPE_LOCAL) THEN
               TYPE_MKL(NLEVEL_MIN:NLEVEL_MAX-1) = NSCARC_MKL_LOCAL
            ELSE IF (TYPE_SCOPE(2) == NSCARC_SCOPE_GLOBAL) THEN
               TYPE_MKL(NLEVEL_MIN:NLEVEL_MAX-1) = NSCARC_MKL_GLOBAL
            ENDIF

         ENDIF

         IF (TYPE_MKL(0) == NSCARC_MKL_COARSE) TYPE_MKL(NLEVEL_MAX) = NSCARC_MKL_GLOBAL
#endif


 
   ! ---------- Global LU-decomposition 
 
   CASE (NSCARC_METHOD_LU)

      ! Only use single (finest) grid level and mark this level for the use of MKL methods

      CALL SCARC_GET_NUMBER_OF_LEVELS(NSCARC_LEVEL_SINGLE)
      TYPE_MKL(NLEVEL_MIN) = TYPE_MKL(0)

 
   ! ---------- Global McKenney-Greengard-Mayo method - only finest level 
 
   CASE (NSCARC_METHOD_MGM)

      CALL SCARC_GET_NUMBER_OF_LEVELS(NSCARC_LEVEL_SINGLE)


END SELECT SELECT_SCARC_METHOD


#ifdef WITH_MKL
 
! Define MKL related logical short names based on number of levels
 
DO NL = NLEVEL_MIN, NLEVEL_MAX
   IS_MKL_LEVEL(NL) = (TYPE_MKL(0)  == NSCARC_MKL_GLOBAL .AND. NL == NLEVEL_MIN) .OR. &
                      (TYPE_MKL(0)  == NSCARC_MKL_COARSE .AND. NL == NLEVEL_MAX) .OR. &
                      (TYPE_MKL(NL) == NSCARC_MKL_GLOBAL)
ENDDO
#endif

END SUBROUTINE SCARC_SETUP_LEVELS


! ------------------------------------------------------------------------------------------------
!> \brief Setup single level in case of default Krylov method
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_GET_NUMBER_OF_LEVELS(NTYPE)
USE SCARC_POINTERS, ONLY: M
INTEGER, INTENT(IN) :: NTYPE
INTEGER :: KLEVEL(3), KLEVEL_MIN, NM, NLEVEL

SELECT_LEVEL_TYPE: SELECT CASE (NTYPE)

   ! only use finest grid level
   CASE(NSCARC_LEVEL_SINGLE)
   
      NLEVEL     = 1
      NLEVEL_MIN = 1
      NLEVEL_MAX = 1
   
   ! determine maximum number of possible levels based on number of grid cells (based on doubling)
   CASE(NSCARC_LEVEL_MULTI)
   
      NLEVEL = NSCARC_LEVEL_MAX
      KLEVEL = NSCARC_LEVEL_MAX
   
      DO NM=1,NMESHES
         M => MESHES(NM)
         KLEVEL(1)=SCARC_GET_MAX_LEVEL(M%IBAR,1)
         IF (.NOT.TWO_D) KLEVEL(2)=SCARC_GET_MAX_LEVEL(M%JBAR,2)
         KLEVEL(3)=SCARC_GET_MAX_LEVEL(M%KBAR,3)
         KLEVEL_MIN = MINVAL(KLEVEL)
         IF (KLEVEL_MIN<NLEVEL) NLEVEL=KLEVEL_MIN
      ENDDO
   
      NLEVEL_MIN  = 1
      IF (IS_MG .OR. IS_CG_MG) THEN
         IF (SCARC_MULTIGRID_LEVEL /= -1) THEN
            NLEVEL_MAX  = NLEVEL_MIN + SCARC_MULTIGRID_LEVEL - 1
         ELSE
            NLEVEL_MAX  = NLEVEL
         ENDIF
      ELSE IF (HAS_TWO_LEVELS) THEN
         IF (SCARC_COARSE_LEVEL /= -1) THEN
            NLEVEL_MAX  = NLEVEL_MIN + SCARC_COARSE_LEVEL - 1
         ELSE
            NLEVEL_MAX  = NLEVEL
         ENDIF
      ENDIF
   
END SELECT SELECT_LEVEL_TYPE

END SUBROUTINE SCARC_GET_NUMBER_OF_LEVELS


! ------------------------------------------------------------------------------------------------
!> \brief Determine maximum number of possible levels 
! In case of GMG- or 2-Level-method, NC must be divisable by 2 at least one time
! ------------------------------------------------------------------------------------------------
INTEGER FUNCTION SCARC_GET_MAX_LEVEL(NC, IOR0)
INTEGER, INTENT(IN) :: NC, IOR0
INTEGER :: NC0, NL

! Print error message if not divisable by 2
IF (IS_GMG .AND. SCARC_MULTIGRID_LEVEL > 1 .AND. MOD(NC,2)/=0) THEN
   SELECT CASE (ABS(IOR0))
      CASE (1)
         CALL SCARC_SHUTDOWN(NSCARC_ERROR_GRID_NUMBERX, SCARC_NONE, NC)
      CASE (2)
         CALL SCARC_SHUTDOWN(NSCARC_ERROR_GRID_NUMBERY, SCARC_NONE, NC)
      CASE (3)
         CALL SCARC_SHUTDOWN(NSCARC_ERROR_GRID_NUMBERZ, SCARC_NONE, NC)
   END SELECT
ENDIF

! Divide by 2 as often as possible or until user defined max-level is reached
IF (SCARC_MULTIGRID_LEVEL > 1) THEN
   NC0=NC
   DO NL=1,NSCARC_LEVEL_MAX
      NC0=NC0/2
      IF (MOD(NC0,2)/=0) EXIT                ! if no longer divisable by two, leave loop ...
      IF (NL==SCARC_MULTIGRID_LEVEL) EXIT    ! if max possible number of levels reached, leave loop ...
      IF (NC0==1) EXIT                       ! if corresponding power of two has been found, leave loop ...
   ENDDO
   SCARC_GET_MAX_LEVEL=NL
ELSE
   SCARC_GET_MAX_LEVEL=NLEVEL_MIN
ENDIF

RETURN
END FUNCTION SCARC_GET_MAX_LEVEL


! ------------------------------------------------------------------------------------------------
!> \brief Allocate basic ScaRC-structures for all needed levels
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_TYPES
USE SCARC_POINTERS, ONLY: S
USE SCARC_POINTER_ROUTINES, ONLY: SCARC_POINT_TO_MESH
INTEGER :: NM

! Basic information for all requested grid levels
ALLOCATE (SCARC(NMESHES), STAT=IERROR)
CALL CHKMEMERR ('SCARC_SETUP', 'SCARC', IERROR)

! Basic solver stack
ALLOCATE (STACK(NSCARC_STACK_MAX), STAT=IERROR)
CALL CHKMEMERR ('SCARC_SETUP', 'STACK', IERROR)

MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_MESH(NM)

   ! Needed information about other meshes
   ALLOCATE (S%OSCARC(NMESHES), STAT=IERROR)
   CALL CHKMEMERR ('SCARC_SETUP_TYPES', 'OSCARC', IERROR)

   ! Information for single grid levels
   ALLOCATE (S%LEVEL(NLEVEL_MIN:NLEVEL_MAX), STAT=IERROR)
   CALL CHKMEMERR ('SCARC_SETUP_TYPES', 'LEVEL', IERROR)

ENDDO MESHES_LOOP

END SUBROUTINE SCARC_SETUP_TYPES


! ----------------------------------------------------------------------------------------------------
!> \brief Assign handles to currently used grid type
!  This routine assumes, that L already points to the correct level NL of mesh NL and
!  additionally sets the requested discretization type
! ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_GRID_TYPE(NTYPE)
INTEGER, INTENT(IN) :: NTYPE

SELECT CASE (NTYPE)
   CASE (NSCARC_GRID_STRUCTURED)
      PRES_ON_WHOLE_DOMAIN = .TRUE.
      TYPE_GRID = NSCARC_GRID_STRUCTURED
      IS_STRUCTURED   = .TRUE.
      IS_UNSTRUCTURED = .FALSE.
   CASE (NSCARC_GRID_UNSTRUCTURED)
      PRES_ON_WHOLE_DOMAIN = .FALSE.
      TYPE_GRID = NSCARC_GRID_UNSTRUCTURED
      IS_STRUCTURED   = .FALSE.
      IS_UNSTRUCTURED = .TRUE.
END SELECT
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) '========> SETUP_GRID_TYPE : TYPE_GRID, PRES_ON_WHOLE_DOMAIN:', TYPE_GRID, PRES_ON_WHOLE_DOMAIN
#endif

END SUBROUTINE SCARC_SETUP_GRID_TYPE


! -----------------------------------------------------------------------------
!> \brief Setup discretization information
! -----------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_GRIDS
USE SCARC_POINTERS, ONLY: M, S, L, G, XCOR, YCOR, ZCOR, XMID, YMID, ZMID
USE SCARC_POINTER_ROUTINES, ONLY: SCARC_POINT_TO_MESH, SCARC_POINT_TO_GRID
INTEGER :: NL, NM, NC, IX, IY, IZ, IO
INTEGER :: IBAR, JBAR, KBAR

CROUTINE = 'SCARC_SETUP_GRIDS'

! ---------- On all grid levels 
! Specify general mesh related geometry information

MESHES_LOOP1: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_MESH(NM)

   ! Store bounds of mesh in SCARC-structure

   S%XS = M%XS;  S%XF = M%XF
   S%YS = M%YS;  S%YF = M%YF
   S%ZS = M%ZS;  S%ZF = M%ZF

   S%IBAR = M%IBAR;  S%JBAR = M%JBAR;  S%KBAR = M%KBAR
   IBAR   = M%IBAR;  JBAR   = M%JBAR;  KBAR   = M%KBAR

   LEVEL_LOOP1: DO NL = NLEVEL_MIN, NLEVEL_MAX

      CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G

      L%NX = IBAR;  L%NY = JBAR;  L%NZ = KBAR

      L%N_CELLS = L%NX*L%NY*L%NZ

      L%N_WALL_CELLS_EXT = M%N_EXTERNAL_WALL_CELLS
      L%N_WALL_CELLS_INT = M%N_INTERNAL_WALL_CELLS
      L%N_WALL_CELLS     = L%N_WALL_CELLS_EXT + L%N_WALL_CELLS_INT

      ! Get coordination information

      L%DX = (S%XF-S%XS)/REAL(L%NX,EB) ; L%DXH = 0.5_EB*L%DX ; L%DX2 = L%DX
      L%DY = (S%YF-S%YS)/REAL(L%NY,EB) ; L%DYH = 0.5_EB*L%DY ; L%DY2 = L%DY
      L%DZ = (S%ZF-S%ZS)/REAL(L%NZ,EB) ; L%DZH = 0.5_EB*L%DZ ; L%DZ2 = L%DZ

      L%DXI  = 1.0_EB/L%DX;  L%DYI  = 1.0_EB/L%DY;  L%DZI =  1.0_EB/L%DZ
      L%DXI2 = L%DXI**2;     L%DYI2 = L%DYI**2;     L%DZI2 = L%DZI**2

      ! Needed in case of GMG with multiple grid levels

      IBAR=IBAR/2
      IF (.NOT.TWO_D) JBAR=JBAR/2
      KBAR=KBAR/2

      IF (NL == NLEVEL_MIN) THEN

         ! On finest level store information about obstructions

         L%N_OBST = M%N_OBST
         ALLOCATE(L%OBST(L%N_OBST), STAT=IERROR)
         CALL ChkMemErr('SCARC_SETUP_GRIDS','OBST',IERROR)

         DO IO = 1, L%N_OBST
            L%OBST(IO)%I1  = M%OBSTRUCTION(IO)%I1
            L%OBST(IO)%I2  = M%OBSTRUCTION(IO)%I2
            L%OBST(IO)%J1  = M%OBSTRUCTION(IO)%J1
            L%OBST(IO)%J2  = M%OBSTRUCTION(IO)%J2
            L%OBST(IO)%K1  = M%OBSTRUCTION(IO)%K1
            L%OBST(IO)%K2  = M%OBSTRUCTION(IO)%K2
         ENDDO

         ! Point to already existing arrays from main FDS program

         XCOR => M%X ;  YCOR => M%Y ;  ZCOR => M%Z
         XMID => M%XC;  YMID => M%YC;  ZMID => M%ZC

#ifdef WITH_SCARC_DEBUG
         IF (L%NX < 10) THEN
            MSG%CFORM1 = "( E14.6)"  ; WRITE(MSG%CFORM1(2:2),'(I1.1)') L%NX
            MSG%CFORM4 = "( I5)"  ; WRITE(MSG%CFORM4(2:2),'(I1.1)') L%NX
         ELSE
            MSG%CFORM1 = "(  E14.6)" ; WRITE(MSG%CFORM1(2:3),'(I2.2)') L%NX
            MSG%CFORM4 = "(  I5)" ; WRITE(MSG%CFORM4(2:3),'(I2.2)') L%NX
         ENDIF
         IF (L%NX+1 < 10) THEN
            MSG%CFORM2 = "( E14.6)"  ; WRITE(MSG%CFORM2(2:2),'(I1.1)') L%NX+1
         ELSE
            MSG%CFORM2 = "(  E14.6)" ; WRITE(MSG%CFORM2(2:3),'(I2.2)') L%NX+1
         ENDIF
         IF (L%NX+2 < 10) THEN
            MSG%CFORM3 = "( E14.6)"  ; WRITE(MSG%CFORM3(2:2),'(I1.1)') L%NX+2
         ELSE
            MSG%CFORM3 = "(  E14.6)" ; WRITE(MSG%CFORM3(2:3),'(I2.2)') L%NX+2
         ENDIF
#endif

      ELSE

         ! Allocate and compute coordinate information for coarser levels

         CALL SCARC_ALLOCATE_REAL1(L%XCOR, 0, L%NX, NSCARC_INIT_NONE, 'L%XCOR', CROUTINE)
         CALL SCARC_ALLOCATE_REAL1(L%YCOR, 0, L%NY, NSCARC_INIT_NONE, 'L%YCOR', CROUTINE)
         CALL SCARC_ALLOCATE_REAL1(L%ZCOR, 0, L%NZ, NSCARC_INIT_NONE, 'L%ZCOR', CROUTINE)

         DO IX = 0, L%NX
            L%XCOR(IX) = S%XS + IX*L%DX
         ENDDO
         DO IY = 0, L%NY
            L%YCOR(IY) = S%YS + IY*L%DY
         ENDDO
         DO IZ = 0, L%NZ
            L%ZCOR(IZ) = S%ZS + IZ*L%DZ
         ENDDO

         ! Allocate and compute midpoint information for coarser levels

         CALL SCARC_ALLOCATE_REAL1(L%XMID, 0, L%NX+1, NSCARC_INIT_NONE, 'L%XMID', CROUTINE)
         CALL SCARC_ALLOCATE_REAL1(L%YMID, 0, L%NY+1, NSCARC_INIT_NONE, 'L%YMID', CROUTINE)
         CALL SCARC_ALLOCATE_REAL1(L%ZMID, 0, L%NZ+1, NSCARC_INIT_NONE, 'L%ZMID', CROUTINE)

         L%XMID(0) = S%XS - 0.5_EB*L%DX
         DO IX = 1, L%NX
            L%XMID(IX) = 0.5_EB*(L%XCOR(IX-1) + L%XCOR(IX))
         ENDDO
         L%XMID(L%NX+1) = S%XF + 0.5_EB*L%DX

         L%YMID(0) = S%YS - 0.5_EB*L%DY
         DO IY = 1, L%NY
            L%YMID(IY) = 0.5_EB*(L%YCOR(IY-1) + L%YCOR(IY))
         ENDDO
         L%YMID(L%NY+1) = S%YF + 0.5_EB*L%DY

         L%ZMID(0) = S%ZS - 0.5_EB*L%DZ
         DO IZ = 1, L%NZ
            L%ZMID(IZ) = 0.5_EB*(L%ZCOR(IZ-1) + L%ZCOR(IZ))
         ENDDO
         L%ZMID(L%NZ+1) = S%ZF + 0.5_EB*L%DZ

         XCOR => L%XCOR;  YCOR => L%YCOR;  ZCOR => L%ZCOR
         XMID => L%XMID;  YMID => L%YMID;  ZMID => L%ZMID

      ENDIF

      ! Allocate vectors for step sizes in different directions

      CALL SCARC_ALLOCATE_REAL1(L%DXL, 0, L%NX, NSCARC_INIT_ZERO, 'L%DXL', CROUTINE)
      CALL SCARC_ALLOCATE_REAL1(L%DYL, 0, L%NY, NSCARC_INIT_ZERO, 'L%DYL', CROUTINE)
      CALL SCARC_ALLOCATE_REAL1(L%DZL, 0, L%NZ, NSCARC_INIT_ZERO, 'L%DZL', CROUTINE)

      ! Set step sizes between cell midpoints, use interior step sizes for ghost cells as initial values
      ! correct sizes for ghost cells are exchanged later

      DO IX = 1, L%NX-1
         L%DXL(IX) = XMID(IX+1) - XMID(IX)
      ENDDO
      L%DXL(0)    = L%DXL(1)
      L%DXL(L%NX) = L%DXL(L%NX-1)

      DO IY = 1, L%NY-1
         L%DYL(IY) = YMID(IY+1) - YMID(IY)
      ENDDO
      L%DYL(0)    = L%DYL(1)
      L%DYL(L%NY) = L%DYL(L%NY-1)

      DO IZ = 1, L%NZ-1
         L%DZL(IZ) = ZMID(IZ+1) - ZMID(IZ)
      ENDDO
      L%DZL(0)    = L%DZL(1)
      L%DZL(L%NZ) = L%DZL(L%NZ-1)

   ENDDO LEVEL_LOOP1
ENDDO MESHES_LOOP1

 
! ---------------------- On finest grid level -------------------------------------------------
! Allocate several arrays for the administration of discretization related data
 
MESHES_LOOP2: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID(NM, NLEVEL_MIN)               ! sets pointers M, L and G

   ! Set pointers to already existing cell and wall index arrays from main program (on finest level)
   CALL SCARC_SETUP_CELL_INDEX(L, M, NLEVEL_MIN)
   CALL SCARC_SETUP_WALL_INDEX(L, G, M, NLEVEL_MIN)

   ! Allocate and initialize IS_SOLID array which indicates the state of a cell (gasphase/solid)

   CALL SCARC_ALLOCATE_LOG3(L%IS_SOLID, 0, L%NX+1, 0, L%NY+1, 0, L%NZ+1, NSCARC_INIT_TRUE, 'L%IS_SOLID', CROUTINE)
   L%IS_SOLID (1:L%NX, 1:L%NY, 1:L%NZ) = .FALSE.

   ! Identify and mark solid obstruction cells in IS_SOLID-part of the discretization

   DO IZ = 1, L%NZ
      DO IY = 1, L%NY
         DO IX = 1, L%NX
            IF (M%SOLID(M%CELL_INDEX(IX, IY, IZ))) L%IS_SOLID(IX, IY, IZ) = .TRUE.
         ENDDO
      ENDDO
   ENDDO

 
   ! If both discretization types (structured/unstructured) must be administrated (MGM method only):
   ! allocate all arrays which are related to a specific discretization type
 
   IF (HAS_MULTIPLE_GRIDS) THEN
 
      ! ---------- First process structured discretization
 
      CALL SCARC_SETUP_GRID_TYPE(NSCARC_GRID_STRUCTURED)
      CALL SCARC_POINT_TO_GRID (NM, NLEVEL_MIN)

      ! Number of local cells per mesh

      CALL SCARC_ALLOCATE_INT1(G%NC_LOCAL , 1, NMESHES, NSCARC_INIT_ZERO, 'G%NC_LOCAL', CROUTINE)
      CALL SCARC_ALLOCATE_INT1(G%NC_OFFSET, 1, NMESHES, NSCARC_INIT_ZERO, 'G%NC_OFFSET', CROUTINE)

      ! Allocate wall information array

      ALLOCATE(G%WALL(L%N_WALL_CELLS), STAT=IERROR)
      CALL ChkMemErr('SCARC_SETUP_GRIDS','WALL',IERROR)

      ! Allocate and preset cell numbers array

      CALL SCARC_ALLOCATE_INT3(G%CELL_NUMBER, 0, L%NX+1, 0, L%NY+1, 0, L%NZ+1, &
                               NSCARC_INIT_UNDEF, 'CELL_NUMBER', CROUTINE)

      ! Allocate index array which specifies I, J, K components for all degrees of freedom

      NC = L%NX * L%NY * L%NZ
      CALL SCARC_ALLOCATE_INT1(G%ICX, 1, NC, NSCARC_INIT_UNDEF, 'G%ICX', CROUTINE)
      CALL SCARC_ALLOCATE_INT1(G%ICY, 1, NC, NSCARC_INIT_UNDEF, 'G%ICY', CROUTINE)
      CALL SCARC_ALLOCATE_INT1(G%ICZ, 1, NC, NSCARC_INIT_UNDEF, 'G%ICZ', CROUTINE)

      ! Define local cell numbers for Poisson equation

      DO IZ=1,L%NZ
         DO IY=1,L%NY
            DO IX=1,L%NX
               G%NC_LOCAL(NM) = G%NC_LOCAL(NM) + 1
               G%CELL_NUMBER(IX,IY,IZ) = G%NC_LOCAL(NM)
               G%ICX(G%NC_LOCAL(NM)) = IX
               G%ICY(G%NC_LOCAL(NM)) = IY
               G%ICZ(G%NC_LOCAL(NM)) = IZ
            ENDDO
         ENDDO
      ENDDO
      G%NC   = G%NC_LOCAL(NM)
      G%NCE  = G%NC_LOCAL(NM)
      G%NCE2 = G%NC_LOCAL(NM)


#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'SETUP_GRIDS: STRUCTURED: NC, NCE, NCE2:', G%NC, G%NCE, G%NCE2
#endif
 
      ! ---------------- Then process unstructured discretization
 
      CALL SCARC_SETUP_GRID_TYPE(NSCARC_GRID_UNSTRUCTURED)
      CALL SCARC_POINT_TO_GRID (NM, NLEVEL_MIN)

      ! Also allocate and preset cell numbers and state arrays for unstructured discretization

      CALL SCARC_ALLOCATE_INT3(G%CELL_NUMBER, 0, L%NX+1, 0, L%NY+1, 0, L%NZ+1, &
                               NSCARC_INIT_UNDEF, 'G%CELL_NUMBER', CROUTINE)

      ! Allocate index array which specifies I, J, K components for all degrees of freedom

      NC = L%NX * L%NY * L%NZ
      CALL SCARC_ALLOCATE_INT1(G%ICX, 1, NC, NSCARC_INIT_UNDEF, 'G%ICX', CROUTINE)
      CALL SCARC_ALLOCATE_INT1(G%ICY, 1, NC, NSCARC_INIT_UNDEF, 'G%ICY', CROUTINE)
      CALL SCARC_ALLOCATE_INT1(G%ICZ, 1, NC, NSCARC_INIT_UNDEF, 'G%ICZ', CROUTINE)

      ! Number of local cells per mesh

      CALL SCARC_ALLOCATE_INT1(G%NC_LOCAL , 1, NMESHES, NSCARC_INIT_ZERO, 'G%NC_LOCAL', CROUTINE)
      CALL SCARC_ALLOCATE_INT1(G%NC_OFFSET, 1, NMESHES, NSCARC_INIT_ZERO, 'G%NC_OFFSET', CROUTINE)

      ! Allocate wall information array

      ALLOCATE(G%WALL(L%N_WALL_CELLS), STAT=IERROR)
      CALL ChkMemErr('SCARC_SETUP_GRIDS','WALL',IERROR)

      ! Define local cell numbers for Poisson equation

      DO IZ=1,L%NZ
         DO IY=1,L%NY
            DO IX=1,L%NX
               IF (.NOT.L%IS_SOLID(IX,IY,IZ)) THEN
                  G%NC_LOCAL(NM) = G%NC_LOCAL(NM) + 1
                  G%CELL_NUMBER(IX,IY,IZ) = G%NC_LOCAL(NM)
                  G%ICX(G%NC_LOCAL(NM)) = IX
                  G%ICY(G%NC_LOCAL(NM)) = IY
                  G%ICZ(G%NC_LOCAL(NM)) = IZ
               ENDIF
            ENDDO
         ENDDO
      ENDDO
      G%NC   = G%NC_LOCAL(NM)
      G%NCE  = G%NC_LOCAL(NM)
      G%NCE2 = G%NC_LOCAL(NM)

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'SETUP_GRIDS: UNSTRUCTURED: NC, NCE, NCE2:', G%NC, G%NCE, G%NCE2
#endif
 
   ! If only one specified type of discretization must be admistrated:
   ! allocate and preset cell numbers and state arrays for requested type of discretization
 
   ELSE

      ! ---------------- Only process specified type of discretization

      CALL SCARC_SETUP_GRID_TYPE(TYPE_GRID)
      CALL SCARC_POINT_TO_GRID (NM, NLEVEL_MIN)

      ! Also allocate and preset cell numbers and state arrays for unstructured discretization

      CALL SCARC_ALLOCATE_INT3(G%CELL_NUMBER, 0, L%NX+1, 0, L%NY+1, 0, L%NZ+1, &
                               NSCARC_INIT_UNDEF, 'CELL_NUMBER', CROUTINE)

      ! Allocate index array which specifies I, J, K components for all degrees of freedom

      NC = L%NX * L%NY * L%NZ
      CALL SCARC_ALLOCATE_INT1(G%ICX, 1, NC, NSCARC_INIT_UNDEF, 'G%ICX', CROUTINE)
      CALL SCARC_ALLOCATE_INT1(G%ICY, 1, NC, NSCARC_INIT_UNDEF, 'G%ICY', CROUTINE)
      CALL SCARC_ALLOCATE_INT1(G%ICZ, 1, NC, NSCARC_INIT_UNDEF, 'G%ICZ', CROUTINE)

      ! Number of local cells per mesh

      CALL SCARC_ALLOCATE_INT1(G%NC_LOCAL , 1, NMESHES, NSCARC_INIT_ZERO, 'G%NC_LOCAL', CROUTINE)
      CALL SCARC_ALLOCATE_INT1(G%NC_OFFSET, 1, NMESHES, NSCARC_INIT_ZERO, 'G%NC_OFFSET', CROUTINE)

      ! Allocate wall information array

      ALLOCATE(G%WALL(L%N_WALL_CELLS), STAT=IERROR)
      CALL ChkMemErr('SCARC_SETUP_GRIDS','WALL',IERROR)

      ! Define local cell numbers for Poisson equation

      DO IZ=1,L%NZ
         DO IY=1,L%NY
            DO IX=1,L%NX
               IF (IS_STRUCTURED .OR. .NOT.L%IS_SOLID(IX,IY,IZ)) THEN

                  G%NC_LOCAL(NM) = G%NC_LOCAL(NM) + 1
                  G%CELL_NUMBER(IX, IY, IZ) = G%NC_LOCAL(NM)

                  G%ICX(G%NC_LOCAL(NM)) = IX
                  G%ICY(G%NC_LOCAL(NM)) = IY
                  G%ICZ(G%NC_LOCAL(NM)) = IZ

               ENDIF
            ENDDO
         ENDDO
      ENDDO
      G%NC   = G%NC_LOCAL(NM)
      G%NCE  = G%NC_LOCAL(NM)
      G%NCE2 = G%NC_LOCAL(NM)

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'SETUP_GRIDS: ', TYPE_GRID,' : NC, NCE, NCE2:', G%NC, G%NCE, G%NCE2
#endif
   ENDIF

ENDDO MESHES_LOOP2

END SUBROUTINE SCARC_SETUP_GRIDS


! -----------------------------------------------------------------------------
!> \brief Setup discretization information on coarser levels
! -----------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_GRID_LEVEL(NL)
USE SCARC_POINTERS, ONLY: LF, LC, GC
USE SCARC_POINTER_ROUTINES, ONLY: SCARC_POINT_TO_MULTIGRID
INTEGER, INTENT(IN) :: NL
INTEGER :: NM, NC, IXF, IYF, IZF, IX, IY, IZ, NSTEP

CROUTINE = 'SCARC_SETUP_GRID_LEVEL'

MESHES_LOOP1: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_MULTIGRID (NM, NLEVEL_MIN, NL)

   CALL SCARC_ALLOCATE_LOG3(LC%IS_SOLID, 0, LC%NX+1, 0, LC%NY+1, 0, LC%NZ+1, NSCARC_INIT_FALSE, 'LC%IS_SOLID', CROUTINE)
   LC%IS_SOLID (1:LC%NX, 1:LC%NY, 1:LC%NZ)  = .FALSE.

   NSTEP = 2**(NL - NLEVEL_MIN)

   SELECT CASE(TYPE_GRID)

 
      ! Get cell numberings for coarser grid in case of structured discretization
 
      CASE (NSCARC_GRID_STRUCTURED)

         CALL SCARC_ALLOCATE_INT1(GC%NC_LOCAL , 1, NMESHES, NSCARC_INIT_ZERO, 'GC%NC_LOCAL', CROUTINE)
         CALL SCARC_ALLOCATE_INT1(GC%NC_OFFSET, 1, NMESHES, NSCARC_INIT_ZERO, 'GC%NC_OFFSET', CROUTINE)

         CALL SCARC_ALLOCATE_INT3(GC%CELL_NUMBER, 0, LC%NX+1, 0, LC%NY+1, 0, LC%NZ+1, &
                                  NSCARC_INIT_UNDEF, 'GC%CELL_NUMBER', CROUTINE)

         NC = LC%NX * LC%NY * LC%NZ
         CALL SCARC_ALLOCATE_INT1(GC%ICX , 1, NC, NSCARC_INIT_UNDEF, 'GC%ICX', CROUTINE)
         CALL SCARC_ALLOCATE_INT1(GC%ICY , 1, NC, NSCARC_INIT_UNDEF, 'GC%ICY', CROUTINE)
         CALL SCARC_ALLOCATE_INT1(GC%ICZ , 1, NC, NSCARC_INIT_UNDEF, 'GC%ICZ', CROUTINE)

         GC%CELL_NUMBER = NSCARC_UNDEF_INT

         DO IZ = 1, LC%NZ
            IZF = (IZ-1)*NSTEP + 1
            DO IY = 1, LC%NY
               IYF = (IY-1)*NSTEP + 1
               DO IX = 1, LC%NX

                  IXF = (IX-1)*NSTEP + 1
                  LC%IS_SOLID(IX,IY,IZ) = LF%IS_SOLID(IXF, IYF, IZF)

                  GC%NC_LOCAL(NM) = GC%NC_LOCAL(NM) + 1
                  GC%CELL_NUMBER(IX,IY,IZ) = GC%NC_LOCAL(NM)

                  GC%ICX(GC%NC_LOCAL(NM)) = IX
                  GC%ICY(GC%NC_LOCAL(NM)) = IY
                  GC%ICZ(GC%NC_LOCAL(NM)) = IZ

               ENDDO
            ENDDO
         ENDDO

         GC%NC = GC%NC_LOCAL(NM)

 
      ! Get cell numberings for coarser grid in case of unstructured discretization
 
      CASE (NSCARC_GRID_UNSTRUCTURED)

         CALL SCARC_ALLOCATE_INT1(GC%NC_LOCAL , 1, NMESHES, NSCARC_INIT_ZERO, 'GC%NC_LOCAL', CROUTINE)
         CALL SCARC_ALLOCATE_INT1(GC%NC_OFFSET, 1, NMESHES, NSCARC_INIT_ZERO, 'GC%NC_OFFSET', CROUTINE)

         CALL SCARC_ALLOCATE_INT3(GC%CELL_NUMBER, 0, LC%NX+1, 0, LC%NY+1, 0, LC%NZ+1, &
                                  NSCARC_INIT_UNDEF, 'GC%CELL_NUMBER', CROUTINE)

         NC = LC%NX * LC%NY * LC%NZ
         CALL SCARC_ALLOCATE_INT1(GC%ICX , 1, NC, NSCARC_INIT_UNDEF, 'GC%ICX', CROUTINE)
         CALL SCARC_ALLOCATE_INT1(GC%ICY , 1, NC, NSCARC_INIT_UNDEF, 'GC%ICY', CROUTINE)
         CALL SCARC_ALLOCATE_INT1(GC%ICZ , 1, NC, NSCARC_INIT_UNDEF, 'GC%ICZ', CROUTINE)

         GC%CELL_NUMBER = NSCARC_UNDEF_INT

         DO IZ = 1, LC%NZ
            IZF = (IZ-1)*NSTEP + 1
            DO IY = 1, LC%NY
               IYF = (IY-1)*NSTEP + 1
               DO IX = 1, LC%NX

                  IXF = (IX-1)*NSTEP + 1
                  LC%IS_SOLID(IX,IY,IZ) = LF%IS_SOLID(IXF, IYF, IZF)

                  IF (.NOT.LF%IS_SOLID(IXF, IYF, IZF)) THEN

                     GC%NC_LOCAL(NM) = GC%NC_LOCAL(NM) + 1
                     GC%CELL_NUMBER(IX,IY,IZ) = GC%NC_LOCAL(NM)

                     GC%ICX(GC%NC_LOCAL(NM)) = IX
                     GC%ICY(GC%NC_LOCAL(NM)) = IY
                     GC%ICZ(GC%NC_LOCAL(NM)) = IZ

                  ENDIF

               ENDDO
            ENDDO
         ENDDO

         GC%NC = GC%NC_LOCAL(NM)

   END SELECT

ENDDO MESHES_LOOP1

END SUBROUTINE SCARC_SETUP_GRID_LEVEL


! ----------------------------------------------------------------------------------------------------
!> \brief Setup structures related to mesh faces on finest grid level
!   - get dimensions for each of the 6 faces of a mesh
!   - get grid width vector along face
!   - get information for adjacent neighbors
!   - allocate pointer arrays for data exchanges with neighbors
! ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_FACES
USE SCARC_POINTERS, ONLY: M, S, L, LC, F, OL
USE SCARC_POINTER_ROUTINES, ONLY: SCARC_POINT_TO_GRID, SCARC_POINT_TO_OTHER_GRID
INTEGER :: NL, NM, NOM
INTEGER :: IFACE, IOR0, JOR0, INBR, IWG, ICW
LOGICAL :: IS_KNOWN(-3:3)
INTEGER :: FACE_NEIGHBORS(-3:3, NSCARC_MAX_FACE_NEIGHBORS)
INTEGER :: MESH_NEIGHBORS(6*NSCARC_MAX_FACE_NEIGHBORS)
INTEGER :: N_FACE_NEIGHBORS(-3:3)
INTEGER :: N_MESH_NEIGHBORS

CROUTINE = 'SCARC_SETUP_FACES'

MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID(NM, NLEVEL_MIN)             ! consider only finest grid level

   ! Allocate FACE arrays on different grid levels

   ALLOCATE(L%FACE(-3:3), STAT=IERROR)
   CALL ChkMemErr('SCARC_SETUP_FACES','FACE',IERROR)

   IF (NLEVEL_MAX > NLEVEL_MIN) THEN
      DO NL = NLEVEL_MIN+1, NLEVEL_MAX
         LC => SCARC(NM)%LEVEL(NL)
         ALLOCATE(LC%FACE(-3:3), STAT=IERROR)
         CALL ChkMemErr('SCARC_SETUP_FACES','FACE',IERROR)
      ENDDO
   ENDIF

   FACE_NEIGHBORS = -1
   MESH_NEIGHBORS = -1

   N_FACE_NEIGHBORS = 0
   N_MESH_NEIGHBORS = 0

   CALL SCARC_SETUP_FACE_BASICS(L)
 
   ! Store first wall cell number for each face
 
   ICW = 1
   FACE_ORDER_LOOP: DO IFACE = 1, 6
      IOR0 = FACE_ORIENTATION(IFACE)
      F => L%FACE(IOR0)
      F%NCW0 = ICW
      ICW = ICW + F%NCW
   ENDDO FACE_ORDER_LOOP

   ! Loop over external wall cells:
   ! store basic data and determine number of adajacent neighbors to each face
 
   EXTERNAL_WALL_CELLS_LOOP: DO IWG = 1, L%N_WALL_CELLS_EXT

      NOM  = M%EXTERNAL_WALL(IWG)%NOM
      IOR0 = M%WALL(IWG)%ONE_D%IOR

      IF (NOM /= 0) THEN
         IS_KNOWN = .FALSE.
         DO JOR0 = -3, 3                                              ! neighbor already known?
            IF (JOR0 == 0) CYCLE
            DO INBR = 1, N_FACE_NEIGHBORS(JOR0)
               IF (FACE_NEIGHBORS(JOR0, INBR) == NOM) THEN
                  IS_KNOWN(JOR0) = .TRUE.
                  EXIT
               ENDIF
            ENDDO
         ENDDO
         IF (.NOT.IS_KNOWN(IOR0)) THEN
            N_FACE_NEIGHBORS(IOR0) = N_FACE_NEIGHBORS(IOR0) + 1       ! increase neighbor counter for face
            FACE_NEIGHBORS(IOR0, N_FACE_NEIGHBORS(IOR0)) = NOM        ! store number of neighbor for face
         ENDIF
         IF (.NOT.ANY(IS_KNOWN)) THEN
            N_MESH_NEIGHBORS = N_MESH_NEIGHBORS + 1                   ! increase neighbor counter for mesh
            MESH_NEIGHBORS(N_FACE_NEIGHBORS(IOR0)) = NOM              ! store number of neighbor for mesh
         ENDIF
      ENDIF

   ENDDO EXTERNAL_WALL_CELLS_LOOP

 
   ! Allocate array which stores numbers of all neighboring meshes
 
   IF (N_MESH_NEIGHBORS /= 0) &
      CALL SCARC_ALLOCATE_INT1(S%NEIGHBORS, 1, N_MESH_NEIGHBORS, NSCARC_INIT_UNDEF, 'S%NEIGHBORS', CROUTINE)
   S%N_NEIGHBORS = N_MESH_NEIGHBORS

   NEIGHBORS_OF_FACE_LOOP: DO IOR0 = -3, 3

      IF (IOR0 == 0) CYCLE NEIGHBORS_OF_FACE_LOOP

      ! If there are neighbors at face IOR0 store information about them

      F => L%FACE(IOR0)
      IF (N_FACE_NEIGHBORS(IOR0) /= 0) THEN

         ! Allocate array for storing the numbers of the single neighbors

         F%N_NEIGHBORS = N_FACE_NEIGHBORS(IOR0)
         CALL SCARC_ALLOCATE_INT1(F%NEIGHBORS, 1, N_FACE_NEIGHBORS(IOR0), NSCARC_INIT_NONE, 'F%NEIGHBORS', CROUTINE)

         ! Store every neighbor and allocate corresponding administration arrays on finest level

         DO INBR = 1, N_FACE_NEIGHBORS(IOR0)

            NOM = FACE_NEIGHBORS(IOR0, INBR)
            F%NEIGHBORS(INBR) = NOM                          ! store NOM as a neighbor of that face and if
            CALL SCARC_STORE_NEIGHBOR(NM, NOM)               ! not already done also as mesh neighbor itself

            CALL SCARC_POINT_TO_OTHER_GRID(NM, NOM, NLEVEL_MIN)
            IF (.NOT.ALLOCATED(OL%FACE)) THEN
               ALLOCATE(OL%FACE(-3:3), STAT=IERROR)
               CALL ChkMemErr('SCARC_SETUP_FACES','OL%FACE',IERROR)
            ENDIF

         ENDDO

      ENDIF
   ENDDO NEIGHBORS_OF_FACE_LOOP

ENDDO MESHES_LOOP

END SUBROUTINE SCARC_SETUP_FACES


! ------------------------------------------------------------------------------------------------
!> \brief Setup subdivision information 
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_SUBDIVISION
USE SCARC_POINTERS, ONLY: SUB
INTEGER, ALLOCATABLE, DIMENSION(:) :: BUFFER_INT
INTEGER, ALLOCATABLE, DIMENSION(:) :: COUNTS_NBR   
INTEGER, ALLOCATABLE, DIMENSION(:) :: DISPLS_NBR    
INTEGER :: N, NM, INBR, IP, MAX_NBR

CROUTINE = 'SCARC_SETUP_SUBDIVISION'

! Determine number of neighbors for each mesh and make them available in a global array
SUB => SUBDIVISION

CALL SCARC_ALLOCATE_INT1 (SUB%N_NEIGHBORS, 1, NMESHES, NSCARC_INIT_ZERO, 'SUB%N_NEIGHBORS', CROUTINE)
DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   SUB%N_NEIGHBORS(NM) = SCARC(NM)%N_NEIGHBORS
ENDDO

IF (N_MPI_PROCESSES > 1) &
   CALL MPI_ALLGATHERV(MPI_IN_PLACE,1,MPI_INTEGER,SUB%N_NEIGHBORS,COUNTS,DISPLS, MPI_INTEGER,MPI_COMM_WORLD,IERROR)
SUB%N_NEIGHBORS_TOTAL = SUM(SUB%N_NEIGHBORS)

CALL SCARC_ALLOCATE_INT1 (COUNTS_NBR, 0, N_MPI_PROCESSES-1, NSCARC_INIT_ZERO, 'COUNTS_NBR', CROUTINE)
CALL SCARC_ALLOCATE_INT1 (DISPLS_NBR, 0, N_MPI_PROCESSES-1, NSCARC_INIT_ZERO, 'DISPLS_NBR', CROUTINE)

DO N = 0, N_MPI_PROCESSES - 1
   DO NM = 1, NMESHES
      IF (PROCESS(NM) == N) COUNTS_NBR(N) = COUNTS_NBR(N) + SUB%N_NEIGHBORS(NM)
   ENDDO
ENDDO
DO N = 1, N_MPI_PROCESSES -1
   DISPLS_NBR(N) = COUNTS_NBR(N-1) + DISPLS_NBR(N-1)
ENDDO

CALL SCARC_ALLOCATE_INT1 (BUFFER_INT, 1, SUB%N_NEIGHBORS_TOTAL, NSCARC_INIT_ZERO, 'BUFFER_INT', CROUTINE)
IP = 1
DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   DO INBR = 1, SCARC(NM)%N_NEIGHBORS
      BUFFER_INT(DISPLS_NBR(PROCESS(NM)) + IP) = SCARC(NM)%NEIGHBORS(INBR)
      IP = IP + 1
   ENDDO
ENDDO

IF (N_MPI_PROCESSES > 1) &
   CALL MPI_ALLGATHERV(MPI_IN_PLACE,1,MPI_INTEGER,BUFFER_INT,COUNTS_NBR,DISPLS_NBR, MPI_INTEGER,MPI_COMM_WORLD,IERROR)
MAX_NBR = MAXVAL(SUB%N_NEIGHBORS)

CALL SCARC_ALLOCATE_INT2 (SUB%NEIGHBORS, 1, MAX_NBR, 1, NMESHES,  NSCARC_INIT_ZERO, 'SUB%NEIGHBORS', CROUTINE)

DO NM = 1, NMESHES
   DO INBR = 1, SUB%N_NEIGHBORS(NM)
      SUB%NEIGHBORS(INBR, NM) = BUFFER_INT(DISPLS_NBR(PROCESS(NM)) + INBR)
   ENDDO
ENDDO

CALL SCARC_DEALLOCATE_INT1(COUNTS_NBR, 'COUNTS_NBR', CROUTINE)
CALL SCARC_DEALLOCATE_INT1(DISPLS_NBR, 'DISPLS_NBR', CROUTINE)
CALL SCARC_DEALLOCATE_INT1(BUFFER_INT, 'BUFFER_INT', CROUTINE)

END SUBROUTINE SCARC_SETUP_SUBDIVISION


! ----------------------------------------------------------------------------------------------------
!> \brief Determine basic data for single faces (orientation, dimensions, numbers)
! ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_FACE_BASICS(L)
USE SCARC_POINTERS, ONLY: F
TYPE (SCARC_LEVEL_TYPE), POINTER, INTENT(IN) :: L
INTEGER:: IOR0

FACES_OF_MESH_LOOP: DO IOR0 = -3, 3

   IF (IOR0 == 0) CYCLE FACES_OF_MESH_LOOP

   F => L%FACE(IOR0)
   
   SELECT CASE (ABS(IOR0))

      ! ---------- Faces in x-direction
 
      CASE (1)

         F%NOP =  L%NX                           ! number of cells between opposite mesh faces

         F%NX  =  1                              ! number of cells in x-direction
         F%NY  =  L%NY                           ! number of cells in y-direction
         F%NZ  =  L%NZ                           ! number of cells in z-direction

         F%NCW =  L%NY*L%NZ                      ! number of wall cells at that face
         F%DH  => L%DXL                          ! step size vector between opposite mesh faces

         F%INCR_BOUNDARY  = L%DXI2               ! contribution due to boundary condition 
         F%SCAL_DIRICHLET = -2.0_EB * L%DXI2

         F%INCRY =  0
         F%INCRZ =  0
         IF (IOR0 > 0) THEN
            F%SCAL_NEUMANN = L%DXI
            F%INCR_FACE    = 2.0_EB/(F%DH(0)*(F%DH(0)+F%DH(1)))
            F%INCR_INSIDE  = 2.0_EB/(F%DH(1)*(F%DH(1)+F%DH(2)))
            F%INCRX =  1                           ! offset to next internal cell in that direction
         ELSE
            F%SCAL_NEUMANN = -L%DXI
            F%INCR_FACE    = 2.0_EB/(F%DH(F%NOP)  *(F%DH(F%NOP-1)+F%DH(F%NOP)))
            F%INCR_INSIDE  = 2.0_EB/(F%DH(F%NOP-1)*(F%DH(F%NOP-2)+F%DH(F%NOP-1)))
            F%INCRX = -1
         ENDIF     
         IF (TWO_D) THEN
            F%INCRS = (/ F%NY, 0, 0, 0, 0,  0, -F%NY /)
         ELSE
            F%INCRS = (/ F%NY, 1, 0, 0, 0, -1, -F%NY /)
         ENDIF
 
      ! ---------- Faces in y-direction
 
      CASE (2)

         F%NOP =  L%NY                   ! dito

         F%NX  =  L%NX
         F%NY  =  1
         F%NZ  =  L%NZ

         F%NCW =  L%NX*L%NZ
         F%DH  => L%DYL

         F%INCR_BOUNDARY  = L%DYI2
         F%SCAL_DIRICHLET = -2.0_EB * L%DYI2

         F%INCRX =  0                           ! offset to next internal cell in that direction
         F%INCRY =  0
         F%INCRZ =  0
         IF (IOR0>0) THEN
            F%SCAL_NEUMANN = L%DYI
            IF (.NOT.TWO_D) THEN
               F%INCR_FACE   = 2.0_EB/(F%DH(0)*(F%DH(0)+F%DH(1)))
               F%INCR_INSIDE = 2.0_EB/(F%DH(1)*(F%DH(1)+F%DH(2)))
               F%INCRY = 1
            ENDIF
         ELSE
            F%SCAL_NEUMANN = -L%DYI
            IF (.NOT.TWO_D) THEN
               F%INCR_FACE   =  2.0_EB/(F%DH(F%NOP)  *(F%DH(F%NOP-1)+F%DH(F%NOP)))
               F%INCR_INSIDE =  2.0_EB/(F%DH(F%NOP-1)*(F%DH(F%NOP-2)+F%DH(F%NOP-1)))
               F%INCRY = -1
            ENDIF
         ENDIF
         IF (TWO_D) THEN
            F%INCRS = (/ 0   , 0, 0, 0,  0, 0, 0     /)             ! special case, not used
         ELSE
            F%INCRS = (/ F%NX, 0, 1, 0, -1, 0, -F%NX /)
         ENDIF

      ! ---------- Faces in z-direction
 
      CASE (3)

         F%NOP =  L%NZ                   ! dito

         F%NX  =  L%NX
         F%NY  =  L%NY
         F%NZ  =  1

         F%NCW =  L%NX*L%NY
         F%DH  => L%DZL

         F%NX  = L%NX
         F%INCR_BOUNDARY  = L%DZI2
         F%SCAL_DIRICHLET = -2.0_EB * L%DZI2

         F%INCRX =  0
         F%INCRY =  0
         IF (IOR0>0) THEN
            F%SCAL_NEUMANN = L%DZI
            F%INCR_FACE    = 2.0_EB/(F%DH(0)*(F%DH(0)+F%DH(1)))
            F%INCR_INSIDE  = 2.0_EB/(F%DH(1)*(F%DH(1)+F%DH(2)))
            F%INCRZ =  1
         ELSE
            F%SCAL_NEUMANN = -L%DZI
            F%INCR_FACE    =  2.0_EB/(F%DH(F%NOP)  *(F%DH(F%NOP-1)+F%DH(F%NOP)))
            F%INCR_INSIDE  =  2.0_EB/(F%DH(F%NOP-1)*(F%DH(F%NOP-2)+F%DH(F%NOP-1)))
            F%INCRZ = -1
         ENDIF
         IF (TWO_D) THEN
            F%INCRS = (/ 0, 0   , 1, 0, -1,     0, 0 /)
         ELSE
            F%INCRS = (/ 0, F%NX, 1, 0, -1, -F%NX, 0 /)
         ENDIF

   END SELECT

ENDDO FACES_OF_MESH_LOOP

END SUBROUTINE SCARC_SETUP_FACE_BASICS


! ----------------------------------------------------------------------------------------------------
!> \brief Setup wall related structures and boundary conditions
! ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_WALLS
USE SCARC_POINTERS, ONLY: M, L, LF, LC, FF, FC, OL, OLF, OLC, G, GC, GF, OG, OGC, OGF, GWC, MWC, EWC
USE SCARC_POINTER_ROUTINES, ONLY: SCARC_POINT_TO_GRID, SCARC_POINT_TO_OTHER_GRID, &  
                                  SCARC_POINT_TO_MULTIGRID, SCARC_POINT_TO_OTHER_MULTIGRID
INTEGER :: NL, NM, NOM
INTEGER :: IREFINE, IFACE, IOR0, JOR0, INBR, IWG, IWC, ICW, IW
LOGICAL :: IS_KNOWN(-3:3), IS_DIRIC, IS_OPEN

CROUTINE = 'SCARC_SETUP_WALLS'
 
! -------- Get dimensionings for wall cells
 
MESHES_LOOP1: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID (NM, NLEVEL_MIN)

 
   ! First loop over external wall cells:
   ! Determine number of adajacent neighbors to each face with corresponding number of IW's
   ! Store neighbors, orientation and number of couplings for a single wall cell
 
   EXTERNAL_WALL_CELLS_LOOP1: DO IWG = 1, L%N_WALL_CELLS_EXT

      MWC => M%WALL(IWG)
      EWC => M%EXTERNAL_WALL(IWG)

      NOM  =  EWC%NOM
      IOR0 =  MWC%ONE_D%IOR

      GWC => G%WALL(IWG)
      GWC%NOM  = NOM                                    ! store number of neighbor in wall cell
      GWC%IOR  = IOR0                                   ! store orientation of that cell

      IF (NOM /= 0) THEN

         IS_KNOWN = .FALSE.
         DO JOR0 = -3, 3
            IF (JOR0 == 0) CYCLE
            DO INBR = 1, L%FACE(JOR0)%N_NEIGHBORS
               IF (L%FACE(JOR0)%NEIGHBORS(INBR) == NOM) THEN
                  IS_KNOWN(JOR0) = .TRUE.
                  EXIT
               ENDIF
            ENDDO
         ENDDO

         G%NCE  = G%NCE  + 1                                                ! increase number of extended grid cells
         IF (HAS_AMG_LEVELS) G%NCE2 = G%NCE2 + 2                            ! increase number of extended grid cells type2
         CALL SCARC_POINT_TO_OTHER_GRID(NM, NOM, NLEVEL_MIN)

         IF (ANY(IS_KNOWN)) OG%NCG = OG%NCG + 1                             ! increase counter for local ghost cells
         IF (OL%GHOST_FIRSTW(IOR0) == 0) OL%GHOST_FIRSTW(IOR0) = OG%NCG     ! save first ghost cell for -IOR0
         IF (OL%GHOST_FIRSTE(IOR0) == 0) OL%GHOST_FIRSTE(IOR0) = OG%NCG     ! save first extended cell for -IOR0
         OL%GHOST_LASTW(IOR0) = OG%NCG                                     
         OL%GHOST_LASTE(IOR0) = OG%NCG                                     

      ENDIF

   ENDDO EXTERNAL_WALL_CELLS_LOOP1
   IF (HAS_AMG_LEVELS) G%ICE2 = G%NCE                                       ! initialize counter for second layer ghost cells

 
   ! Then process internal wall cells
 
   INTERNAL_WALL_CELLS_LOOP1: DO IWG = L%N_WALL_CELLS_EXT+1, L%N_WALL_CELLS_EXT+L%N_WALL_CELLS_INT

      MWC => M%WALL(IWG)
      GWC => G%WALL(IWG)

      GWC%IOR  = MWC%ONE_D%IOR
      GWC%NOM  = 0

      GWC%BTYPE = NEUMANN
      GWC%BOUNDARY_TYPE = M%WALL(IWG)%BOUNDARY_TYPE

      GWC%IXG =  MWC%ONE_D%II                        ! ghost cell indices
      GWC%IYG =  MWC%ONE_D%JJ
      GWC%IZG =  MWC%ONE_D%KK

      GWC%IXW =  MWC%ONE_D%IIG                       ! (internal) wall cell indices
      GWC%IYW =  MWC%ONE_D%JJG
      GWC%IZW =  MWC%ONE_D%KKG

   ENDDO INTERNAL_WALL_CELLS_LOOP1

 
   ! Allocate corresponding pointer arrays for data exchanges with neighbors
 
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'TYPE_GRID =', TYPE_GRID,': SETTING ICE_TO_IWG in length ', G%NCE
#endif
   IF (G%NCE > G%NC) THEN
      G%N_FINE = G%NCE
      CALL SCARC_ALLOCATE_INT1 (G%ICE_TO_IWG, G%NC+1, G%NCE, NSCARC_INIT_ZERO, 'G%ICE_TO_IWG', CROUTINE)
      CALL SCARC_ALLOCATE_INT1 (G%ICE_TO_ICN, G%NC+1, G%NCE, NSCARC_INIT_ZERO, 'G%ICE_TO_IWG', CROUTINE)
   ENDIF

   FACE_NEIGHBORS_LOOP: DO IOR0 = -3, 3

      IF (IOR0 == 0) CYCLE FACE_NEIGHBORS_LOOP
      DO INBR = 1, L%FACE(IOR0)%N_NEIGHBORS

         NOM = L%FACE(IOR0)%NEIGHBORS(INBR)
         CALL SCARC_POINT_TO_OTHER_GRID(NM, NOM, NLEVEL_MIN)

         CALL SCARC_ALLOCATE_INT1 (OG%ICG_TO_IWG, 1, OG%NCG, NSCARC_INIT_ZERO, 'OG%ICG_TO_IWG', CROUTINE)

         IF (HAS_AMG_LEVELS) THEN
            CALL SCARC_ALLOCATE_INT2 (OG%ICG_TO_ICW, 1, OG%NCG, 1, 2, NSCARC_INIT_ZERO, 'OG%ICG_TO_ICW', CROUTINE)
            CALL SCARC_ALLOCATE_INT2 (OG%ICG_TO_ICE, 1, OG%NCG, 1, 2, NSCARC_INIT_ZERO, 'OG%ICG_TO_ICE', CROUTINE)
         ELSE
            CALL SCARC_ALLOCATE_INT2 (OG%ICG_TO_ICW, 1, OG%NCG, 1, 1, NSCARC_INIT_ZERO, 'OG%ICG_TO_ICW', CROUTINE)
            CALL SCARC_ALLOCATE_INT2 (OG%ICG_TO_ICE, 1, OG%NCG, 1, 1, NSCARC_INIT_ZERO, 'OG%ICG_TO_ICE', CROUTINE)
         ENDIF

      ENDDO

   ENDDO FACE_NEIGHBORS_LOOP

 
   ! Second loop over external wall cells:
   ! Store detailed coordinate and cell data and get type of boundary condition
 
   G%ICE = G%NC
   WALL_CELLS_LOOP2: DO IWG = 1, L%N_WALL_CELLS_EXT

      IOR0 = G%WALL(IWG)%IOR
      NOM  = G%WALL(IWG)%NOM

      MWC => M%WALL(IWG)
      EWC => M%EXTERNAL_WALL(IWG)

 
      ! Preset ScaRC's boundary type indicator BTYPE
      ! INTERNAL  : the global Poisson problem is solved, so no BC's along mesh interfaces are needed
      ! DIRICHLET : - in the structured case face-wise BC-settings are used ccording to original FFT-solver
      !               (this also allows to use FFT as local preconditioner)
      !             - in the unstructured case Dirichlet BC's are only used for open boundary cells
      ! NEUMANN   : is used for the rest
 
      IS_DIRIC = MWC%PRESSURE_BC_INDEX == DIRICHLET
      IS_OPEN  = MWC%BOUNDARY_TYPE     == OPEN_BOUNDARY

      GWC => G%WALL(IWG)

      IF (EWC%NOM /= 0) THEN
         GWC%BTYPE = INTERNAL
      ELSE IF ((IS_STRUCTURED .AND. IS_DIRIC) .OR. (IS_UNSTRUCTURED .AND. IS_OPEN)) THEN
         GWC%BTYPE = DIRICHLET
         G%N_DIRIC = G%N_DIRIC + 1
      ELSE
         GWC%BTYPE = NEUMANN
         G%N_NEUMANN = G%N_NEUMANN + 1
      ENDIF

      GWC%BOUNDARY_TYPE = MWC%BOUNDARY_TYPE

      GWC%IXG = MWC%ONE_D%II                                 ! ghost cell indices
      GWC%IYG = MWC%ONE_D%JJ
      GWC%IZG = MWC%ONE_D%KK

      GWC%IXW = MWC%ONE_D%IIG                                ! (internal) wall cell indices
      GWC%IYW = MWC%ONE_D%JJG
      GWC%IZW = MWC%ONE_D%KKG

      ! If there exists a neighbor for that wall cell, setup corresponding neighborship information
      IF (NOM /= 0) CALL SCARC_SETUP_WALL_NEIGHBOR(G, OG, &
                                                   EWC%IIO_MIN, EWC%IIO_MAX, &
                                                   EWC%JJO_MIN, EWC%JJO_MAX, &
                                                   EWC%KKO_MIN, EWC%KKO_MAX, &
                                                   IWG, NM, NOM, NLEVEL_MIN)

   ENDDO WALL_CELLS_LOOP2

ENDDO MESHES_LOOP1

 
! Set dimensions on finest level for requested type(s) of discretization
! and mapping from local to global cell numbering
 
CALL SCARC_SETUP_DIMENSIONS(NLEVEL_MIN)

 
! -------- For multi-level variants get discretization information and dimensions on coarser levels
 
DO NL = NLEVEL_MIN+1, NLEVEL_MAX
   CALL SCARC_SETUP_GRID_LEVEL(NL)
   CALL SCARC_SETUP_DIMENSIONS(NL)
ENDDO

 
! -------- Check whether there are no Dirichlet BC's available - TODO: Check !!!
 
MESH_INT = 0                            
RANK_INT = 0

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   CALL SCARC_POINT_TO_GRID (NM, NLEVEL_MIN)
   MESH_INT(NM) = G%N_DIRIC   
   RANK_INT = RANK_INT + MESH_INT(NM)
ENDDO

IF (N_MPI_PROCESSES>1) &
   CALL MPI_ALLREDUCE(MPI_IN_PLACE, RANK_INT, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, IERROR)
N_DIRIC_GLOBAL(NLEVEL_MIN) = RANK_INT

IS_PURE_NEUMANN = N_DIRIC_GLOBAL(NLEVEL_MIN) == 0 .AND. &
                  (TYPE_PRECON /= NSCARC_RELAX_FFT .OR. TYPE_PRECON /= NSCARC_RELAX_FFTO)


 
! -------- Only for multi-level variants 
! (twolevel-CG or GMG method as main solver or preconditioner):
! Determine WALL, FACE and OSCARC types for coarser levels
 
MULTI_LEVEL_IF: IF (HAS_GMG_LEVELS) THEN

   MESHES_LOOP3: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

      IREFINE=1
      MULTI_LEVELS_LOOP: DO NL = NLEVEL_MIN+1, NLEVEL_MAX

         CALL SCARC_POINT_TO_MULTIGRID(NM, NL-1, NL)

         IREFINE=IREFINE*2
         IF (IS_GMG) CALL SCARC_CHECK_DIVISIBILITY(GF%NCE-GF%NC, 'GF%NCE')

         ! Initialize counts for overlapping and wall cells
 
         GC%NCE = GC%NC + (GF%NCE-GF%NC)/2
         GC%ICE = GC%NC

         LC%N_WALL_CELLS_EXT = SCARC_COUNT_EXTERNAL_WALL_CELLS(LF, LC, GF)
         LC%N_WALL_CELLS_INT = SCARC_COUNT_INTERNAL_WALL_CELLS(LF, LC, GC)

         LC%N_WALL_CELLS = LC%N_WALL_CELLS_EXT + LC%N_WALL_CELLS_INT

         SELECT CASE(TYPE_GRID)
            CASE (NSCARC_GRID_STRUCTURED)
               GC%NW = LC%N_WALL_CELLS_EXT
            CASE (NSCARC_GRID_UNSTRUCTURED)
               GC%NW = LC%N_WALL_CELLS_EXT + LC%N_WALL_CELLS_INT
         END SELECT

         ALLOCATE(GC%WALL(LC%N_WALL_CELLS), STAT=IERROR)
         CALL ChkMemErr('SCARC_SETUP_NEIGHBORS','WALL',IERROR)

         ! First allocate administrative mapping arrays for own mesh if there is an overlap
 
         IF (GC%NCE > GC%NC) THEN
            CALL SCARC_ALLOCATE_INT1 (GC%ICE_TO_IWG , GC%NC+1, GC%NCE, NSCARC_INIT_ZERO, 'GC%ICE_TO_IWG', CROUTINE)
            CALL SCARC_ALLOCATE_INT1 (GC%ICE_TO_ICN , GC%NC+1, GC%NCE, NSCARC_INIT_ZERO, 'GC%ICE_TO_ICN', CROUTINE)
         ENDIF

         ! Setup basic face information for coarser mesh
 
         CALL SCARC_SETUP_FACE_BASICS(LC)

         IWC = 1
         IWG = 1
         ICW = 1
         FACES_LOOP: DO IFACE = 1, 6

            IOR0 = FACE_ORIENTATION(IFACE)

            FF => LF%FACE(IOR0)
            FC => LC%FACE(IOR0)

            ! initialize FACE type for coarser mesh
            FC%NCW0 = ICW
            FC%N_NEIGHBORS = FF%N_NEIGHBORS

            IF (FC%N_NEIGHBORS /= 0) &
               CALL SCARC_ALLOCATE_INT1(FC%NEIGHBORS, 1, FC%N_NEIGHBORS, NSCARC_INIT_NONE, 'FC%FACE_NEIGHBORS', CROUTINE)
            DO INBR= 1, FC%N_NEIGHBORS
               FC%NEIGHBORS(INBR) = FF%NEIGHBORS(INBR)
            ENDDO
            
            FC%NCW = FC%NX * FC%NY * FC%NZ                                ! get number of wall cells for that face
            ICW = ICW + FC%NCW                                            ! increase global wall cell counter

            ! Get related data and pointer structures for every mesh neighbor

            IF (LF%FACE(IOR0)%N_NEIGHBORS /= 0) THEN
               DO INBR = 1, LF%FACE(IOR0)%N_NEIGHBORS

                  NOM = LF%FACE(IOR0)%NEIGHBORS(INBR)
                  CALL SCARC_POINT_TO_OTHER_MULTIGRID(NM, NOM, NL-1, NL)

                  IF (IS_GMG) THEN
                     CALL SCARC_CHECK_DIVISIBILITY(OLF%N_WALL_CELLS_LOCAL, 'OLF%N_WALL_CELLS_LOCAL')
                     CALL SCARC_CHECK_DIVISIBILITY(OGF%NCG, 'OGF%NCG')
                  ENDIF

                  IF (.NOT.TWO_D) THEN
                     OLC%N_WALL_CELLS_LOCAL = OLF%N_WALL_CELLS_LOCAL/4
                     OGC%NCG = OGF%NCG/4
                  ELSE
                     OLC%N_WALL_CELLS_LOCAL = OLF%N_WALL_CELLS_LOCAL/2
                     OGC%NCG = OGF%NCG/2
                  ENDIF

                  CALL SCARC_SETUP_EXCHANGE_DIMENSIONS(GC, OGC, NOM, IREFINE)

                  CALL SCARC_ALLOCATE_INT1(OGC%ICG_TO_IWG, 1, OGC%NCG, NSCARC_INIT_ZERO, 'OGC%ICG_TO_IWG', CROUTINE)
                  CALL SCARC_ALLOCATE_INT2(OGC%ICG_TO_ICW, 1, OGC%NCG, 1, 1, NSCARC_INIT_ZERO, 'OGC%ICG_TO_ICW', CROUTINE)
                  CALL SCARC_ALLOCATE_INT2(OGC%ICG_TO_ICE, 1, OGC%NCG, 1, 1, NSCARC_INIT_ZERO, 'OGC%ICG_TO_ICE', CROUTINE)

               ENDDO
            ENDIF

            ! Setup complete wall information on coarser mesh level

            CALL SCARC_SETUP_WALL_LEVEL(LF, LC, GF, GC, IOR0, IWC, IREFINE, NM, NL)

         ENDDO FACES_LOOP

         CALL SCARC_SETUP_CELL_INDEX (LC, M, NL)
         CALL SCARC_SETUP_WALL_COORDS(LC, GC)
         CALL SCARC_SETUP_WALL_INDEX (LC, GC, M, NL)

         ! Setup order in which ghost cells are processed during data exchanges

         WALLCELLS_LOOP: DO IW = 1, LC%N_WALL_CELLS
         
           NOM = GC%WALL(IW)%NOM
           IF (NOM /= 0) THEN
              IOR0 = GC%WALL(IW)%IOR
              CALL SCARC_POINT_TO_OTHER_MULTIGRID(NM, NOM, NL-1, NL)
              OGC%ICG2 = OGC%ICG2 + 1
              IF (OLC%GHOST_FIRSTW(IOR0) == 0) OLC%GHOST_FIRSTW(IOR0) = OGC%ICG2
              IF (OLC%GHOST_FIRSTE(IOR0) == 0) OLC%GHOST_FIRSTE(IOR0) = OGC%ICG2
              OLC%GHOST_LASTW(IOR0) = OGC%ICG2 
              OLC%GHOST_LASTE(IOR0) = OGC%ICG2 
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'NL, IOR0, NOM, FIRSTW, LASTW, FIRSTE, LASTE:', NL, IOR0, NOM, &
                      OLC%GHOST_FIRSTW(IOR0),OLC%GHOST_LASTW(IOR0),OLC%GHOST_FIRSTE(IOR0),OLC%GHOST_LASTE(IOR0)
#endif
           ENDIF
         
         ENDDO WALLCELLS_LOOP

      ENDDO MULTI_LEVELS_LOOP
   ENDDO MESHES_LOOP3
ENDIF MULTI_LEVEL_IF


! Correct boundary types for cells adjacent to obstructions on ghost cells

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   CALL SCARC_POINT_TO_GRID(NM, NLEVEL_MIN)                  ! sets level and grid pointers L and G
   IF (IS_UNSTRUCTURED) THEN
      CALL SCARC_IDENTIFY_INTERNAL_NEUMANNS(L, G)
      IF (.NOT.HAS_AMG_LEVELS) THEN
         DO NL = NLEVEL_MIN+1, NLEVEL_MAX
            CALL SCARC_IDENTIFY_INTERNAL_NEUMANNS(L, G)
         ENDDO
      ENDIF
   ENDIF
ENDDO

! Debug FACE, WALL and DISCRET structures - only if directive SCARC_DEBUG is set

#ifdef WITH_SCARC_DEBUG
CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_STACK, NLEVEL_MIN, 'STACK')
CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_FACE , NLEVEL_MIN, 'FACE')
CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_WALL , NLEVEL_MIN, 'WALL')
CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_GRID , NLEVEL_MIN, 'DISCRET')
#endif

END SUBROUTINE SCARC_SETUP_WALLS


! -----------------------------------------------------------------------------------------
!> \brief Store all neighbors of a mesh
! -----------------------------------------------------------------------------------------
SUBROUTINE SCARC_STORE_NEIGHBOR(NM, NOM)
INTEGER, INTENT(IN) :: NM, NOM
INTEGER :: INBR
DO INBR = 1, SCARC(NM)%N_NEIGHBORS
   IF (SCARC(NM)%NEIGHBORS(INBR) == NSCARC_UNDEF_INT) EXIT      ! not found, to be stored
   IF (SCARC(NM)%NEIGHBORS(INBR) == NOM) RETURN                 ! nothing to do, already stored
ENDDO
SCARC(NM)%NEIGHBORS(INBR) = NOM
RETURN
END SUBROUTINE SCARC_STORE_NEIGHBOR


! -----------------------------------------------------------------------------------------
!> \brief Setup cells indexing array on coarser grid levels in case of MG method
! -----------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_CELL_INDEX(L, M, NL)
USE SCARC_POINTERS, ONLY: OB
TYPE (MESH_TYPE), POINTER, INTENT(IN) :: M
TYPE (SCARC_LEVEL_TYPE), POINTER, INTENT(IN) :: L
INTEGER, INTENT(IN) :: NL
INTEGER :: I, J, K, NOBST

CROUTINE = 'SCARC_SETUP_CELL_INDEX'

! If finest level, the corresponding CELL_INDEX array is already available by surrounding routines
! on coarser levels, it must still be computed

IF (NL == NLEVEL_MIN) THEN

   L%CELL_INDEX_PTR => M%CELL_INDEX

ELSE

   CALL SCARC_ALLOCATE_INT3(L%CELL_INDEX, 0, L%NX+1, 0, L%NY+1, 0, L%NZ+1, NSCARC_INIT_ZERO, 'L%CELL_INDEX', CROUTINE)
   L%CELL_INDEX_PTR => L%CELL_INDEX
   L%N_CELL_INDEX = 0

   ! Preset it for all grid cells
 
   DO K=0,L%NZ+1
      DO J=0,L%NY+1
         DO I=0,1
            IF (L%CELL_INDEX(I,J,K)==0) THEN
               L%N_CELL_INDEX = L%N_CELL_INDEX + 1
               L%CELL_INDEX(I,J,K) = L%N_CELL_INDEX
            ENDIF
         ENDDO
         DO I=L%NX,L%NX+1
            IF (L%CELL_INDEX(I,J,K)==0) THEN
               L%N_CELL_INDEX = L%N_CELL_INDEX + 1
               L%CELL_INDEX(I,J,K) = L%N_CELL_INDEX
            ENDIF
         ENDDO
      ENDDO
   ENDDO

   DO K=0,L%NZ+1
      DO I=0,L%NX+1
         DO J=0,1
            IF (L%CELL_INDEX(I,J,K)==0) THEN
               L%N_CELL_INDEX = L%N_CELL_INDEX + 1
               L%CELL_INDEX(I,J,K) = L%N_CELL_INDEX
            ENDIF
         ENDDO
         DO J=L%NY,L%NY+1
            IF (L%CELL_INDEX(I,J,K)==0) THEN
               L%N_CELL_INDEX = L%N_CELL_INDEX + 1
               L%CELL_INDEX(I,J,K) = L%N_CELL_INDEX
            ENDIF
         ENDDO
      ENDDO
   ENDDO

   DO J=0,L%NY+1
      DO I=0,L%NX+1
         DO K=0,1
            IF (L%CELL_INDEX(I,J,K)==0) THEN
               L%N_CELL_INDEX = L%N_CELL_INDEX + 1
               L%CELL_INDEX(I,J,K) = L%N_CELL_INDEX
            ENDIF
         ENDDO
         DO K=L%NZ,L%NZ+1
            IF (L%CELL_INDEX(I,J,K)==0) THEN
               L%N_CELL_INDEX = L%N_CELL_INDEX + 1
               L%CELL_INDEX(I,J,K) = L%N_CELL_INDEX
            ENDIF
         ENDDO
      ENDDO
   ENDDO

 
   ! Consider cells in obstructions
 
   DO NOBST=1,L%N_OBST
      OB => L%OBST(NOBST)
      DO K=OB%K1,OB%K2+1
         DO J=OB%J1,OB%J2+1
            DO I=OB%I1,OB%I2+1
               IF (L%CELL_INDEX(I,J,K)==0) THEN
                  L%N_CELL_INDEX = L%N_CELL_INDEX + 1
                  L%CELL_INDEX(I,J,K) = L%N_CELL_INDEX
               ENDIF
            ENDDO
         ENDDO
      ENDDO
   ENDDO

ENDIF

END SUBROUTINE SCARC_SETUP_CELL_INDEX


! -----------------------------------------------------------------------------------------
!> \brief Setup wall cells indexing array on coarser grid levels
! -----------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_WALL_INDEX(L, G, M, NL)
TYPE (MESH_TYPE), POINTER, INTENT(IN) :: M
TYPE (SCARC_LEVEL_TYPE), POINTER, INTENT(IN) :: L
TYPE (SCARC_GRID_TYPE),  POINTER, INTENT(IN) :: G
INTEGER, INTENT(IN) :: NL
INTEGER :: I, J, K, ICG, IW, IOR0

CROUTINE = 'SCARC_SETUP_WALL_INDEX'

! if on finest level, the array WALL_INDEX is already available by surrounding routines
IF (NL == NLEVEL_MIN) THEN

   L%WALL_INDEX_PTR => M%WALL_INDEX

   ! if on coarser levels, it must still be computed
ELSE

   CALL SCARC_ALLOCATE_INT2(L%WALL_INDEX, 1, L%N_CELL_INDEX, -3, 3, NSCARC_INIT_ZERO, 'L%WALL_INDEX', CROUTINE)
   L%WALL_INDEX_PTR => L%WALL_INDEX

   DO IW = 1, L%N_WALL_CELLS_EXT

      I = G%WALL(IW)%IXW
      J = G%WALL(IW)%IYW
      K = G%WALL(IW)%IZW

      IOR0 = G%WALL(IW)%IOR
      ICG  = L%CELL_INDEX(I,J,K)

      L%WALL_INDEX(ICG,-IOR0) = IW

   ENDDO

ENDIF

END SUBROUTINE SCARC_SETUP_WALL_INDEX


! -------------------------------------------------------------------------------------------------
!> \brief Setup all necessary information for a wall cell with neighbor in case of MG method
! Number of obstructions on coarse level is the same as on fine level
! TODO: Only works for special cases which run for GMG, must still be extended!!
! -------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_WALL_COORDS(L, G)
USE SCARC_POINTERS, ONLY: OB, GWC
TYPE (SCARC_LEVEL_TYPE), POINTER, INTENT(IN) :: L
TYPE (SCARC_GRID_TYPE), POINTER, INTENT(IN)  :: G
INTEGER :: IC, IO, IWC
INTEGER :: I, J, K

IWC = L%N_WALL_CELLS_EXT + 1
DO IO = 1, L%N_OBST

   OB => L%OBST(IO)

   ! Analyze IOR = 1

   I = OB%I1
   DO K = OB%K1+1, OB%K2
      DO J = OB%J1+1, OB%J2
         IC = G%CELL_NUMBER(I, J, K)
         IF (IC == NSCARC_UNDEF_INT .OR. IC > G%NC) CYCLE
         GWC => G%WALL(IWC)
         GWC%IXW = I+1; GWC%IYW = J; GWC%IZW = K
         GWC%IXG = I  ; GWC%IYG = J; GWC%IZG = K
         GWC%IOR = 1
         GWC%BTYPE = NEUMANN
         GWC%BOUNDARY_TYPE = SOLID_BOUNDARY
         IWC = IWC + 1
      ENDDO
   ENDDO

   ! Analyze IOR = -1
   I = OB%I2
   DO K = OB%K1+1, OB%K2
      DO J = OB%J1+1, OB%J2
         IC = G%CELL_NUMBER(I+1, J, K)
         IF (IC == NSCARC_UNDEF_INT .OR. IC > G%NC) CYCLE
         GWC => G%WALL(IWC)
         GWC%IXW = I  ; GWC%IYW = J; GWC%IZW = K
         GWC%IXG = I+1; GWC%IYG = J; GWC%IZG = K
         GWC%IOR =-1
         GWC%BTYPE = NEUMANN
         GWC%BOUNDARY_TYPE = SOLID_BOUNDARY
         IWC = IWC + 1
      ENDDO
   ENDDO

   ! Analyze IOR = 2

   J = OB%J1
   DO K = OB%K1+1, OB%K2
      DO I = OB%I1+1, OB%I2
         IC = G%CELL_NUMBER(I, J, K)
         IF (IC == NSCARC_UNDEF_INT .OR. IC > G%NC) CYCLE
         GWC => G%WALL(IWC)
         GWC%IXW = I; GWC%IYW = J+1; GWC%IZW = K
         GWC%IXG = I; GWC%IYG = J  ; GWC%IZG = K
         GWC%IOR = 2
         GWC%BTYPE = NEUMANN
         GWC%BOUNDARY_TYPE = SOLID_BOUNDARY
         IWC = IWC + 1
      ENDDO
   ENDDO

   ! Analyze IOR = -2

   J = OB%J2
   DO K = OB%K1+1, OB%K2
      DO I = OB%I1+1, OB%I2
         IC = G%CELL_NUMBER(I, J+1, K)
         IF (IC == NSCARC_UNDEF_INT .OR. IC > G%NC) CYCLE
         GWC => G%WALL(IWC)
         GWC%IXW = I; GWC%IYW = J  ; GWC%IZW = K
         GWC%IXG = I; GWC%IYG = J+1; GWC%IZG = K
         GWC%IOR =-2
         GWC%BTYPE = NEUMANN
         GWC%BOUNDARY_TYPE = SOLID_BOUNDARY
         IWC = IWC + 1
      ENDDO
   ENDDO

   ! Analyze IOR = 3

   K = OB%K1
   DO J = OB%J1+1, OB%J2
      DO I = OB%I1+1, OB%I2
         IC = G%CELL_NUMBER(I, J, K)
         IF (IC == NSCARC_UNDEF_INT .OR. IC > G%NC) CYCLE
         GWC => G%WALL(IWC)
         GWC%IXW = I; GWC%IYW = J; GWC%IZW = K+1
         GWC%IXG = I; GWC%IYG = J; GWC%IZG = K
         GWC%IOR = 3
         GWC%BTYPE = NEUMANN
         GWC%BOUNDARY_TYPE = SOLID_BOUNDARY
         IWC = IWC + 1
      ENDDO
   ENDDO

   ! Analyze IOR = -3

   K = OB%K2
   DO J = OB%J1+1, OB%J2
      DO I = OB%I1+1, OB%I2
         IC = G%CELL_NUMBER(I, J, K+1)
         IF (IC == NSCARC_UNDEF_INT .OR. IC > G%NC) CYCLE
         GWC => G%WALL(IWC)
         GWC%IXW = I; GWC%IYW = J; GWC%IZW = K
         GWC%IXG = I; GWC%IYG = J; GWC%IZG = K+1
         GWC%IOR =-3
         GWC%BTYPE = NEUMANN
         GWC%BOUNDARY_TYPE = SOLID_BOUNDARY
         IWC = IWC + 1
      ENDDO
   ENDDO

ENDDO

END SUBROUTINE SCARC_SETUP_WALL_COORDS


! -------------------------------------------------------------------------------------------------
!> \brief Correct boundary type array related to internal obstructions on ghost cells
! -------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_IDENTIFY_INTERNAL_NEUMANNS(L, G) 
USE SCARC_POINTERS, ONLY: GWC
TYPE (SCARC_LEVEL_TYPE), POINTER, INTENT(IN) :: L
TYPE (SCARC_GRID_TYPE),  POINTER, INTENT(IN) :: G
INTEGER :: IWG
INTEGER :: IX, IY, IZ, IOR0, BTYPE0

DO IWG = 1, L%N_WALL_CELLS_EXT

   GWC => G%WALL(IWG)
   IF (GWC%NOM == 0) CYCLE                    ! TODO: equal or not equal ??

   IX = GWC%IXW
   IY = GWC%IYW
   IZ = GWC%IZW

   IOR0   = GWC%IOR
   BTYPE0 = GWC%BTYPE

   ! ICG = L%CELL_INDEX_PTR(IX, IY, IZ)
   ! IWG = L%WALL_INDEX_PTR(ICG, IOR0)
   ! IF (GWC%BOUNDARY_TYPE /= INTERPOLATED_BOUNDARY) GWC%BTYPE=NEUMANN
   ! IF (L%IS_SOLID(IX, IY, IZ)) GWC%BTYPE=NEUMANN

   IF (GWC%BOUNDARY_TYPE == SOLID_BOUNDARY) GWC%BTYPE=NEUMANN

ENDDO

END SUBROUTINE SCARC_IDENTIFY_INTERNAL_NEUMANNS


! -------------------------------------------------------------------------------------------------
!> \brief Setup all necessary information for a wall cell with neighbor
! -------------------------------------------------------------------------------------------------
INTEGER FUNCTION SCARC_COUNT_EXTERNAL_WALL_CELLS(LF, LC, GF)
TYPE (SCARC_LEVEL_TYPE), POINTER, INTENT(IN) :: LF, LC
TYPE (SCARC_GRID_TYPE),  POINTER, INTENT(IN) :: GF
INTEGER :: IXC, IYC, IZC
INTEGER :: IXF, IYF, IZF
INTEGER :: IWC, ICF(4)=0, IWF(4)=0, IOR0

ICF = 0
IWC = 0
IWF = 0

IF (TWO_D) THEN
   IYC = 1
   IYF = 1

   ! IOR = 1

   IOR0 = 1
   DO IZC = 1, LC%NZ
      IZF = 2*IZC - 1
      ICF(1) = LF%CELL_INDEX_PTR(1  , IYF  , IZF  )
      ICF(2) = LF%CELL_INDEX_PTR(1  , IYF  , IZF+1)
      IF (IS_EXTERNAL_WALLCELL(LF, GF, IOR0, ICF, 2)) IWC = IWC + 1
   ENDDO

   ! IOR = -1

   IOR0 = -1
   DO IZC = 1, LC%NZ
      IZF = 2*IZC - 1
      ICF(1) = LF%CELL_INDEX_PTR(LF%NX, IYF  , IZF  )
      ICF(2) = LF%CELL_INDEX_PTR(LF%NX, IYF  , IZF+1)
      IF (IS_EXTERNAL_WALLCELL(LF, GF, IOR0, ICF, 2)) IWC = IWC + 1
   ENDDO

   ! IOR = 2

   IOR0 = 2
   DO IZC = 1, LC%NZ
      IZF = 2*IZC - 1
      DO IXC = 1, LC%NX
         IXF = 2*IXC - 1
         ICF(1) = LF%CELL_INDEX_PTR(IXF  , IYF, IZF  )
         ICF(2) = LF%CELL_INDEX_PTR(IXF+1, IYF, IZF  )
         ICF(3) = LF%CELL_INDEX_PTR(IXF  , IYF, IZF+1)
         ICF(4) = LF%CELL_INDEX_PTR(IXF+1, IYF, IZF+1)
         IF (IS_EXTERNAL_WALLCELL(LF, GF, IOR0, ICF, 4)) IWC = IWC + 1
      ENDDO
   ENDDO

   ! IOR = -2

   IOR0 = -2
   IXF  = 1
   DO IZC = 1, LC%NZ
      IZF = 2*IZC - 1
      DO IXC = 1, LC%NX
         IXF = 2*IXC - 1
         ICF(1) = LF%CELL_INDEX_PTR(IXF    , IYF, IZF  )
         ICF(2) = LF%CELL_INDEX_PTR(IXF+1  , IYF, IZF  )
         ICF(3) = LF%CELL_INDEX_PTR(IXF    , IYF, IZF+1)
         ICF(4) = LF%CELL_INDEX_PTR(IXF+1  , IYF, IZF+1)
         IF (IS_EXTERNAL_WALLCELL(LF, GF, IOR0, ICF, 4)) IWC = IWC + 1
      ENDDO
   ENDDO

   ! IOR = 3

   IOR0 = 3
   DO IXC = 1, LC%NX
      IXF = 2*IXC - 1
      ICF(1) = LF%CELL_INDEX_PTR(IXF    , IYF  , 1)
      ICF(2) = LF%CELL_INDEX_PTR(IXF+1  , IYF  , 1)
      IF (IS_EXTERNAL_WALLCELL(LF, GF, IOR0, ICF, 2)) IWC = IWC + 1
   ENDDO

   ! IOR = -3

   IOR0 = -3
   DO IXC = 1, LC%NX
      IXF = 2*IXC - 1
      ICF(1) = LF%CELL_INDEX_PTR(IXF  , IYF  , LF%NZ)
      ICF(2) = LF%CELL_INDEX_PTR(IXF+1, IYF  , LF%NZ)
      IF (IS_EXTERNAL_WALLCELL(LF, GF, IOR0, ICF, 2)) IWC = IWC + 1
   ENDDO

ELSE

   ! IOR = 1

   IOR0 = 1
   DO IZC = 1, LC%NZ
      IZF = 2*IZC - 1
      DO IYC = 1, LC%NY
         IYF = 2*IYC - 1
         ICF(1) = LF%CELL_INDEX_PTR(1  , IYF  , IZF  )
         ICF(2) = LF%CELL_INDEX_PTR(1  , IYF+1, IZF  )
         ICF(3) = LF%CELL_INDEX_PTR(1  , IYF  , IZF+1)
         ICF(4) = LF%CELL_INDEX_PTR(1  , IYF+1, IZF+1)
         IF (IS_EXTERNAL_WALLCELL(LF, GF, IOR0, ICF, 4)) IWC = IWC + 1
      ENDDO
   ENDDO

   ! IOR = -1

   IOR0 = -1
   DO IZC = 1, LC%NZ
      IZF = 2*IZC - 1
      DO IYC = 1, LC%NY
         IYF = 2*IYC - 1
         ICF(1) = LF%CELL_INDEX_PTR(LF%NX, IYF  , IZF  )
         ICF(2) = LF%CELL_INDEX_PTR(LF%NX, IYF+1, IZF  )
         ICF(3) = LF%CELL_INDEX_PTR(LF%NX, IYF  , IZF+1)
         ICF(4) = LF%CELL_INDEX_PTR(LF%NX, IYF+1, IZF+1)
         IF (IS_EXTERNAL_WALLCELL(LF, GF, IOR0, ICF, 4)) IWC = IWC + 1
      ENDDO
   ENDDO

   ! IOR = 2

   IOR0 = 2
   DO IZC = 1, LC%NZ
      IZF = 2*IZC - 1
      DO IXC = 1, LC%NX
         IXF = 2*IXC - 1
         ICF(1) = LF%CELL_INDEX_PTR(IXF  , 1, IZF  )
         ICF(2) = LF%CELL_INDEX_PTR(IXF+1, 1, IZF  )
         ICF(3) = LF%CELL_INDEX_PTR(IXF  , 1, IZF+1)
         ICF(4) = LF%CELL_INDEX_PTR(IXF+1, 1, IZF+1)
         IF (IS_EXTERNAL_WALLCELL(LF, GF, IOR0, ICF, 4)) IWC = IWC + 1
      ENDDO
   ENDDO

   ! IOR = -2

   IOR0 = -2
   IXF  = 1
   DO IZC = 1, LC%NZ
      IZF = 2*IZC - 1
      DO IXC = 1, LC%NX
         IXF = 2*IXC - 1
         ICF(1) = LF%CELL_INDEX_PTR(IXF    , LF%NY, IZF  )
         ICF(2) = LF%CELL_INDEX_PTR(IXF+1  , LF%NY, IZF  )
         ICF(3) = LF%CELL_INDEX_PTR(IXF    , LF%NY, IZF+1)
         ICF(4) = LF%CELL_INDEX_PTR(IXF+1  , LF%NY, IZF+1)
         IF (IS_EXTERNAL_WALLCELL(LF, GF, IOR0, ICF, 4)) IWC = IWC + 1
      ENDDO
   ENDDO

   ! IOR = 3

   IOR0 = 3
   DO IYC = 1, LC%NY
      IYF = 2*IYC - 1
      DO IXC = 1, LC%NX
         IXF = 2*IXC - 1
         ICF(1) = LF%CELL_INDEX_PTR(IXF    , IYF  , 1)
         ICF(2) = LF%CELL_INDEX_PTR(IXF+1  , IYF  , 1)
         ICF(3) = LF%CELL_INDEX_PTR(IXF    , IYF+1, 1)
         ICF(4) = LF%CELL_INDEX_PTR(IXF+1  , IYF+1, 1)
         IF (IS_EXTERNAL_WALLCELL(LF, GF, IOR0, ICF, 4)) IWC = IWC + 1
      ENDDO
   ENDDO

   ! IOR = -3

   IOR0 = -3
   DO IYC = 1, LC%NY
      IYF = 2*IYC - 1
      DO IXC = 1, LC%NX
         IXF = 2*IXC - 1
         ICF(1) = LF%CELL_INDEX_PTR(IXF  , IYF  , LF%NZ)
         ICF(2) = LF%CELL_INDEX_PTR(IXF+1, IYF  , LF%NZ)
         ICF(3) = LF%CELL_INDEX_PTR(IXF  , IYF+1, LF%NZ)
         ICF(4) = LF%CELL_INDEX_PTR(IXF+1, IYF+1, LF%NZ)
         IF (IS_EXTERNAL_WALLCELL(LF, GF, IOR0, ICF, 4)) IWC = IWC + 1
      ENDDO
   ENDDO

ENDIF

SCARC_COUNT_EXTERNAL_WALL_CELLS = IWC
END FUNCTION SCARC_COUNT_EXTERNAL_WALL_CELLS


! -------------------------------------------------------------------------------------------------
!> \brief Setup all necessary information for a wall cell with neighbor
! -------------------------------------------------------------------------------------------------
INTEGER FUNCTION SCARC_COUNT_INTERNAL_WALL_CELLS(LF, LC, GC)
USE SCARC_POINTERS, ONLY: OBF, OBC
TYPE (SCARC_LEVEL_TYPE), POINTER, INTENT(IN) :: LF, LC
TYPE (SCARC_GRID_TYPE),  POINTER, INTENT(IN) :: GC
INTEGER :: IWC, NW_INT
INTEGER :: IC, IO
INTEGER :: I, J, K

LC%N_OBST = LF%N_OBST                   ! Number of obstructions is the same on all levels

ALLOCATE(LC%OBST(LC%N_OBST), STAT=IERROR)
CALL ChkMemErr('SCARC_COUNT_INTERNAL_WALL_CELLS','OBST',IERROR)

NW_INT = 0
IWC = LC%N_WALL_CELLS_EXT + 1

DO IO = 1, LF%N_OBST

   OBF => LF%OBST(IO)
   OBC => LC%OBST(IO)

   OBC%I1 = (OBF%I1+1)/2
   OBC%I2 =  OBF%I2/2

   IF (TWO_D) THEN
      OBC%J1 = 0
      OBC%J2 = 1
   ELSE
      OBC%J1 = (OBF%J1+1)/2
      OBC%J2 =  OBF%J2/2
   ENDIF

   OBC%K1 = (OBF%K1+1)/2
   OBC%K2 =  OBF%K2/2

   ! Analyze IOR = 1

   I = OBC%I1
   DO K = OBC%K1+1, OBC%K2
      DO J = OBC%J1+1, OBC%J2
         IC = GC%CELL_NUMBER(I, J, K)
         IF (IC == NSCARC_UNDEF_INT) CYCLE
         NW_INT = NW_INT + 1
      ENDDO
   ENDDO

   ! Analyze IOR = -1

   I = OBC%I2
   DO K = OBC%K1+1, OBC%K2
      DO J = OBC%J1+1, OBC%J2
         IC = GC%CELL_NUMBER(I+1, J, K)
         IF (IC == NSCARC_UNDEF_INT) CYCLE
         NW_INT = NW_INT + 1
      ENDDO
   ENDDO

   ! Analyze IOR = 2

   J = OBC%J1
   DO K = OBC%K1+1, OBC%K2
      DO I = OBC%I1+1, OBC%I2
         IC = GC%CELL_NUMBER(I, J, K)
         IF (IC == NSCARC_UNDEF_INT) CYCLE
         NW_INT = NW_INT + 1
      ENDDO
   ENDDO

   ! Analyze IOR = -2

   J = OBC%J2
   DO K = OBC%K1+1, OBC%K2
      DO I = OBC%I1+1, OBC%I2
         IC = GC%CELL_NUMBER(I, J+1, K)
         IF (IC == NSCARC_UNDEF_INT) CYCLE
         NW_INT = NW_INT + 1
      ENDDO
   ENDDO

   ! Analyze IOR = 3

   K = OBC%K1
   DO J = OBC%J1+1, OBC%J2
      DO I = OBC%I1+1, OBC%I2
         IC = GC%CELL_NUMBER(I, J, K)
         IF (IC == NSCARC_UNDEF_INT) CYCLE
         NW_INT = NW_INT + 1
      ENDDO
   ENDDO

   ! Analyze IOR = -3

   K = OBC%K2
   DO J = OBC%J1+1, OBC%J2
      DO I = OBC%I1+1, OBC%I2
         IC = GC%CELL_NUMBER(I, J, K+1)
         IF (IC == NSCARC_UNDEF_INT) CYCLE
         NW_INT = NW_INT + 1
      ENDDO
   ENDDO

ENDDO

SCARC_COUNT_INTERNAL_WALL_CELLS = NW_INT
END FUNCTION SCARC_COUNT_INTERNAL_WALL_CELLS


! -------------------------------------------------------------------------------------------------
!> \brief Count external wall cells on specified face if mesh
! -------------------------------------------------------------------------------------------------
LOGICAL FUNCTION IS_EXTERNAL_WALLCELL(L, G, IOR0, ICF, NCNT)
TYPE (SCARC_LEVEL_TYPE), POINTER, INTENT(IN) :: L
TYPE (SCARC_GRID_TYPE),  POINTER, INTENT(IN) :: G
INTEGER, INTENT(IN) :: IOR0, NCNT
INTEGER, DIMENSION(:), INTENT(IN) :: ICF
INTEGER :: I, IWF_LAST, IWF(4)=0
REAL(EB) :: BSUM

IS_EXTERNAL_WALLCELL = .FALSE.

DO I = 1, NCNT
   IWF(I) = L%WALL_INDEX_PTR(ICF(I), -IOR0)
ENDDO

BSUM = 0.0_EB
IWF_LAST = 0

DO I = 1, NCNT
   IF (IWF(I)>0) THEN
      BSUM = BSUM + REAL(G%WALL(IWF(I))%BTYPE,EB)
      IWF_LAST = IWF(I)
   ENDIF
ENDDO

IF (IWF_LAST == 0) RETURN
IF (ABS(BSUM/REAL(NCNT,EB) - REAL(G%WALL(IWF_LAST)%BTYPE,EB)) < 1E-12) THEN
   IS_EXTERNAL_WALLCELL = .TRUE.
   RETURN
ELSE
   CALL SCARC_SHUTDOWN(NSCARC_ERROR_BOUNDARY_SUM, SCARC_NONE, IOR0)
ENDIF

END FUNCTION IS_EXTERNAL_WALLCELL


! -------------------------------------------------------------------------------------------------
!> \brief Setup all necessary information for a wall cell with neighbor
! -------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_WALL_NEIGHBOR(G, OG, NX1, NX2, NY1, NY2, NZ1, NZ2, IWG, NM, NOM, NL)
USE SCARC_POINTERS, ONLY: GWC
USE SCARC_POINTER_ROUTINES, ONLY: SCARC_POINT_TO_OTHER_GRID
TYPE (SCARC_GRID_TYPE),  POINTER, INTENT(IN) :: G, OG
INTEGER, INTENT(IN) :: NX1, NX2, NY1, NY2, NZ1, NZ2
INTEGER, INTENT(IN) :: IWG, NM, NOM, NL
INTEGER :: NOMX, NOMY, NOMZ
INTEGER :: ICG, ICE, IX, IY, IZ, JL, IXW, IYW, IZW, IXG, IYG, IZG

CALL SCARC_POINT_TO_OTHER_GRID (NM, NOM, NL)

ICE  = G%ICE
ICG  = OG%ICG

! set neighboring coordinates
GWC => G%WALL(IWG)

GWC%IXN(1) = NX1
GWC%IXN(2) = NX2
GWC%IYN(1) = NY1
GWC%IYN(2) = NY2
GWC%IZN(1) = NZ1
GWC%IZN(2) = NZ2

NOMX = MESHES(NOM)%IBAR
NOMY = MESHES(NOM)%JBAR
NOMZ = MESHES(NOM)%KBAR

IF (NL > 1) THEN
   DO JL = 2, NL
      NOMX = MESHES(NOM)%IBAR/NL
      IF (.NOT.TWO_D) NOMY = MESHES(NOM)%JBAR/NL
      NOMZ = MESHES(NOM)%KBAR/NL
   ENDDO
ENDIF

! store information about overlapped cells and set mapping arrays
DO IZ = NZ1, NZ2
   DO IY = NY1, NY2
      DO IX = NX1, NX2

         ICG  = ICG  + 1
         ICE  = ICE  + 1

         GWC%ICE = ICE                                         ! number of extended grid cell
         GWC%ICG = ICG                                         ! number of ghost grid cell

         G%ICE_TO_IWG(ICE) = IWG                               ! map extended cell to global wall cell
         
         IXG = G%WALL(IWG)%IXG
         IYG = G%WALL(IWG)%IYG
         IZG = G%WALL(IWG)%IZG

         IXW = G%WALL(IWG)%IXW
         IYW = G%WALL(IWG)%IYW
         IZW = G%WALL(IWG)%IZW

         G%CELL_NUMBER(IXG, IYG, IZG) = ICE

         OG%ICG_TO_IWG(ICG)    = IWG                              ! map ghost cell to global wall cell
         OG%ICG_TO_ICW(ICG, 1) = G%CELL_NUMBER(IXW, IYW, IZW)     ! get cell number of adjacent internal cell

      ENDDO
   ENDDO
ENDDO

G%ICE  = ICE                                                   ! store extended cell counter
OG%ICG = ICG                                                   ! store ghost cell counter
OG%ICG_TO_IWG(ICG) = IWG                                       ! map local wall cell to global wall cell
OG%ICG_TO_ICE(ICG, 1) = ICE                                    ! map local wall cell to global wall cell

END SUBROUTINE SCARC_SETUP_WALL_NEIGHBOR


! ----------------------------------------------------------------------------------------------------
!> \brief Check divisibility by 2 of a given number of elements (in one grid direction)
! ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_CHECK_DIVISIBILITY(NN, CDIR)
INTEGER, INTENT(IN) :: NN
CHARACTER(*) , INTENT(IN) :: CDIR
IF (MOD(NN,2) /= 0) CALL SCARC_SHUTDOWN(NSCARC_ERROR_GRID_NUMBER, CDIR, NSCARC_NONE)
END SUBROUTINE SCARC_CHECK_DIVISIBILITY


! ----------------------------------------------------------------------------------------------------
!> \brief Set wall cell information on coarse level
! ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_WALL_LEVEL(LF, LC, GF, GC, IOR0, IWC, IREFINE, NM, NL)
USE SCARC_POINTERS, ONLY: FF, FC, WF, WC, OGC
USE SCARC_POINTER_ROUTINES, ONLY: SCARC_POINT_TO_OTHER_MULTIGRID
TYPE (SCARC_LEVEL_TYPE), POINTER, INTENT(IN) :: LF, LC
TYPE (SCARC_GRID_TYPE),  POINTER, INTENT(IN) :: GF, GC
INTEGER, INTENT(INOUT) :: IWC
INTEGER, INTENT(IN) :: IOR0, IREFINE, NM, NL
INTEGER :: IWF(4) , IBCF(4), NOMF(4)
INTEGER :: IX,  IY,  IZ, I
INTEGER :: NX1, NY1, NZ1
INTEGER :: NX2, NY2, NZ2
INTEGER :: IX1, IY1, IZ1
INTEGER :: IX2, IY2, IZ2
INTEGER :: IDIFF, JDIFF, KDIFF

WC => GC%WALL
WF => GF%WALL

! set coordinate dimensions for correspoding face

SELECT CASE (ABS(IOR0))
   CASE (1)
      IF (IOR0 > 0) THEN                                           ! set dimensions for wall cell counting
         NX1 = 0; NX2 = 0
      ELSE
         NX1 = LC%NX+1; NX2 = LC%NX+1
      ENDIF
      NY1 = 1;  NY2 = LC%NY
      NZ1 = 1;  NZ2 = LC%NZ
   CASE (2)
      NX1 = 1; NX2 = LC%NX
      IF (IOR0 > 0) THEN
         NY1 = 0; NY2 = 0
      ELSE
         NY1 = LC%NY+1; NY2 = LC%NY+1
      ENDIF
      NZ1 = 1; NZ2 = LC%NZ
   CASE (3)
      NX1 = 1; NX2 = LC%NX
      NY1 = 1; NY2 = LC%NY
      IF (IOR0 > 0) THEN
         NZ1 = 0; NZ2 = 0
      ELSE
         NZ1 =LC%NZ+1; NZ2 =LC%NZ+1
      ENDIF
END SELECT

 
! Loop over all wall cells of face IOR0
 
DO IZ = NZ1, NZ2
   DO IY = NY1, NY2
      DO IX = NX1, NX2

         ! Set orientation of neiboring face, indices of ghost and adjacent cell for coarse IW

         WC(IWC)%IOR = IOR0

         FF => LF%FACE(IOR0)
         FC => LC%FACE(IOR0)

         SELECT CASE (IOR0)
            CASE (1)
               WC(IWC)%ICW = (IZ-1)*LC%NX*LC%NY + (IY-1)*LC%NX + IX + 1
            CASE (-1)
               WC(IWC)%ICW = (IZ-1)*LC%NX*LC%NY + (IY-1)*LC%NX + IX - 1
            CASE (2)
               WC(IWC)%ICW = (IZ-1)*LC%NX*LC%NY +  IY   *LC%NX + IX
            CASE (-2)
               WC(IWC)%ICW = (IZ-1)*LC%NX*LC%NY + (IY-2)*LC%NX + IX
            CASE (3)
               WC(IWC)%ICW =  IZ   *LC%NX*LC%NY + (IY-1)*LC%NX + IX
            CASE (-3)
               WC(IWC)%ICW = (IZ-2)*LC%NX*LC%NY + (IY-1)*LC%NX + IX
         END SELECT

         WC(IWC)%IOR = IOR0

         WC(IWC)%IXG = IX
         WC(IWC)%IYG = IY
         WC(IWC)%IZG = IZ

         SELECT CASE (IOR0)
            CASE (1)
               WC(IWC)%IXW = IX+1
               WC(IWC)%IYW = IY
               WC(IWC)%IZW = IZ
            CASE (-1)
               WC(IWC)%IXW = IX-1
               WC(IWC)%IYW = IY
               WC(IWC)%IZW = IZ
            CASE (2)
               WC(IWC)%IXW = IX
               WC(IWC)%IYW = IY+1
               WC(IWC)%IZW = IZ
            CASE (-2)
               WC(IWC)%IXW = IX
               WC(IWC)%IYW = IY-1
               WC(IWC)%IZW = IZ
            CASE (3)
               WC(IWC)%IXW = IX
               WC(IWC)%IYW = IY
               WC(IWC)%IZW = IZ+1
            CASE (-3)
               WC(IWC)%IXW = IX
               WC(IWC)%IYW = IY
               WC(IWC)%IZW = IZ-1
         END SELECT

         ! ------------------------------------------------------------
         !  2D-version
         ! ------------------------------------------------------------
         IF (TWO_D) THEN

            ! determine fine IW's, which must be merged to one coarse IW

            SELECT CASE (ABS(IOR0))
               CASE ( 1)
                  IWF(1) = FF%NCW0 + 2*(IZ-1)
               CASE ( 2)
                  IWF(1) = FF%NCW0 + 2*(IZ-1)*LF%NX + 2*(IX - 1)
               CASE ( 3)
                  IWF(1) = FF%NCW0 + 2*(IX-1)
            END SELECT
            IWF(2) = IWF(1)+1

            ! set fine cell neighbors (they must be the same for all fine IW's)

            NOMF(1) = WF(IWF(1))%NOM
            NOMF(2) = WF(IWF(2))%NOM
            IF (NOMF(1) /= NOMF(2)) CALL SCARC_SHUTDOWN(NSCARC_ERROR_NEIGHBOR_TYPE, SCARC_NONE, NOMF(1))

            WC(IWC)%NOM = NOMF(1)

            ! set corresponding pressure_bc_index on coarser level

            IBCF(1) = WF(IWF(1))%BTYPE
            IBCF(2) = WF(IWF(2))%BTYPE
            IF (IBCF(1) == INTERNAL .OR. IBCF(2) == INTERNAL) THEN
               WC(IWC)%BTYPE = INTERNAL
            ELSE IF (IBCF(1) == DIRICHLET .OR. IBCF(2) == DIRICHLET) THEN
               WC(IWC)%BTYPE = DIRICHLET
            ELSE
               WC(IWC)%BTYPE = NEUMANN
            ENDIF

            ! set corresponding pressure_bc_index on coarser level

            IBCF(1) = WF(IWF(1))%BOUNDARY_TYPE
            IBCF(2) = WF(IWF(2))%BOUNDARY_TYPE
            IF (IBCF(1)==NULL_BOUNDARY .OR. IBCF(2)==NULL_BOUNDARY) THEN
               WC(IWC)%BOUNDARY_TYPE = NULL_BOUNDARY
            ELSE IF (IBCF(1)==SOLID_BOUNDARY .OR. IBCF(2)==SOLID_BOUNDARY) THEN
               WC(IWC)%BOUNDARY_TYPE = SOLID_BOUNDARY
            ELSE IF (IBCF(1)==OPEN_BOUNDARY .OR. IBCF(2)==OPEN_BOUNDARY) THEN
               WC(IWC)%BOUNDARY_TYPE = OPEN_BOUNDARY
            ELSE IF (IBCF(1)==INTERPOLATED_BOUNDARY .OR. IBCF(2)==INTERPOLATED_BOUNDARY) THEN
               WC(IWC)%BOUNDARY_TYPE = INTERPOLATED_BOUNDARY
            ELSE IF (IBCF(1)==MIRROR_BOUNDARY .OR. IBCF(2)==MIRROR_BOUNDARY) THEN
               WC(IWC)%BOUNDARY_TYPE = MIRROR_BOUNDARY
            ELSE
               CALL SCARC_SHUTDOWN(NSCARC_ERROR_BOUNDARY_TYPE, SCARC_NONE, IBCF(1))
            ENDIF

            ! in case of an internal boundary set neighboring wall cells

            IF (NOMF(1) > 0) THEN

               CALL SCARC_POINT_TO_OTHER_MULTIGRID (NM, NOMF(1), NL-1, NL)

               IY1 = 1
               IY2 = 1
               SELECT CASE (ABS(IOR0))
                  CASE (1)
                     KDIFF = WF(IWF(2))%IZN(1) - WF(IWF(1))%IZN(1)
                     IF (KDIFF == 1) THEN
                        IZ1 = WF(IWF(2))%IZN(2)/2
                        IZ2 = IZ1
                     ELSE IF (KDIFF == 2) THEN
                        IZ1 = WF(IWF(1))%IZN(2)/2
                        IZ2 = WF(IWF(2))%IZN(2)/2
                     ELSE IF (KDIFF == 0) THEN
                        IZ1 = (WF(IWF(1))%IZN(2)+1)/2
                        IZ2 = IZ1
                     ELSE
                        CALL SCARC_SHUTDOWN(NSCARC_ERROR_GRID_RESOLUTION, SCARC_NONE, IOR0)
                     ENDIF
                  CASE (3)
                     IDIFF = WF(IWF(2))%IXN(1) - WF(IWF(1))%IXN(1)
                     IF (IDIFF == 1) THEN
                        IX1 = WF(IWF(2))%IXN(2)/2
                        IX2 = IX1
                     ELSE IF (IDIFF == 2) THEN
                        IX1 = WF(IWF(1))%IXN(2)/2
                        IX2 = WF(IWF(2))%IXN(2)/2
                     ELSE IF (IDIFF == 0) THEN
                        IX1 = (WF(IWF(1))%IXN(2)+1)/2
                        IX2 = IX1
                     ELSE
                        CALL SCARC_SHUTDOWN(NSCARC_ERROR_GRID_RESOLUTION, SCARC_NONE, IOR0)
                     ENDIF
               END SELECT

               SELECT CASE (IOR0)
                  CASE (1)
                     IX1 = MESHES(NOMF(1))%IBAR/IREFINE
                     IX2 = IX1
                  CASE (-1)
                     IX1 = 1
                     IX2 = 1
                  CASE (3)
                     IZ1 = MESHES(NOMF(1))%KBAR/IREFINE
                     IZ2 = IZ1
                  CASE (-3)
                     IZ1 = 1
                     IZ2 = 1
               END SELECT

               WC(IWC)%IXN(1) = IX1
               WC(IWC)%IYN(1) = 1
               WC(IWC)%IZN(1) = IZ1
               WC(IWC)%IXN(2) = IX2
               WC(IWC)%IYN(2) = 1
               WC(IWC)%IZN(2) = IZ2

                
               ! Allocate and specify ICN and ICE arrays for OC
                
               CALL SCARC_SETUP_WALL_NEIGHBOR(GC, OGC, IX1, IX2, 1, 1, IZ1, IZ2, IWC, NM, NOMF(1), NL)

            ENDIF

         ! ------------------------------------------------------------
         ! 3D-version
         ! ------------------------------------------------------------
         ELSE

            ! determine fine IW's, which must be merged to one coarse IW

            SELECT CASE (ABS(IOR0))
               CASE (1)
                  IWF(1) = FF%NCW0 + (2*IZ-2)*LF%NY + 2*IY - 2
                  IWF(3) = FF%NCW0 + (2*IZ-1)*LF%NY + 2*IY - 2
               CASE (2)
                  IWF(1) = FF%NCW0 + (2*IZ-2)*LF%NX + 2*IX - 2
                  IWF(3) = FF%NCW0 + (2*IZ-1)*LF%NX + 2*IX - 2
               CASE (3)
                  IWF(1) = FF%NCW0 + (2*IY-2)*LF%NX + 2*IX - 2
                  IWF(3) = FF%NCW0 + (2*IY-1)*LF%NX + 2*IX - 2
            END SELECT
            IWF(2) = IWF(1)+1
            IWF(4) = IWF(3)+1

            ! set fine cell neighbors (they must be the same for all fine IW's)

            DO I=1,4
               NOMF(I) = WF(IWF(I))%NOM
            ENDDO

            IF (NOMF(1) /= NOMF(2) .OR. NOMF(1) /= NOMF(3) .OR. NOMF(1) /= NOMF(4)) &
               CALL SCARC_SHUTDOWN(NSCARC_ERROR_NEIGHBOR_TYPE, SCARC_NONE, IOR0)
            WC(IWC)%NOM = NOMF(1)

            ! set corresponding pressure_bc_index on coarser level

            DO I=1,4
               IBCF(I) = WF(IWF(I))%BTYPE
            ENDDO
            IF (IBCF(1)==INTERNAL.OR.IBCF(2)==INTERNAL.OR.&
               IBCF(3)==INTERNAL.OR.IBCF(4)==INTERNAL) THEN
               WC(IWC)%BTYPE =INTERNAL
            ELSE IF (IBCF(1)==DIRICHLET.OR.IBCF(2)==DIRICHLET.OR.&
               IBCF(3)==DIRICHLET.OR.IBCF(4)==DIRICHLET) THEN
               WC(IWC)%BTYPE =DIRICHLET
            ELSE
               WC(IWC)%BTYPE =NEUMANN
            ENDIF

            ! set corresponding pressure_bc_index on coarser level

            DO I=1,4
               IBCF(I) = WF(IWF(I))%BOUNDARY_TYPE
            ENDDO
            IF (IBCF(1)==NULL_BOUNDARY.OR.IBCF(2)==NULL_BOUNDARY.OR.&
               IBCF(3)==NULL_BOUNDARY.OR.IBCF(4)==NULL_BOUNDARY) THEN
               WC(IWC)%BOUNDARY_TYPE = NULL_BOUNDARY
            ELSE IF (IBCF(1)==SOLID_BOUNDARY.OR.IBCF(2)==SOLID_BOUNDARY.OR.&
               IBCF(3)==SOLID_BOUNDARY.OR.IBCF(4)==SOLID_BOUNDARY) THEN
               WC(IWC)%BOUNDARY_TYPE = SOLID_BOUNDARY
            ELSE IF (IBCF(1)==OPEN_BOUNDARY.OR.IBCF(2)==OPEN_BOUNDARY.OR.&
               IBCF(3)==OPEN_BOUNDARY.OR.IBCF(4)==OPEN_BOUNDARY) THEN
               WC(IWC)%BOUNDARY_TYPE = OPEN_BOUNDARY
            ELSE IF (IBCF(1)==INTERPOLATED_BOUNDARY.OR.IBCF(2)==INTERPOLATED_BOUNDARY.OR.&
               IBCF(3)==INTERPOLATED_BOUNDARY.OR.IBCF(4)==INTERPOLATED_BOUNDARY) THEN
               WC(IWC)%BOUNDARY_TYPE = INTERPOLATED_BOUNDARY
            ELSE IF (IBCF(1)==MIRROR_BOUNDARY.OR.IBCF(2)==MIRROR_BOUNDARY.OR.&
               IBCF(3)==MIRROR_BOUNDARY.OR.IBCF(4)==MIRROR_BOUNDARY) THEN
               WC(IWC)%BOUNDARY_TYPE = MIRROR_BOUNDARY
            ELSE
               CALL SCARC_SHUTDOWN(NSCARC_ERROR_BOUNDARY_TYPE, SCARC_NONE, -999)
            ENDIF

            ! in case of an internal boundary set WALL(10:15,IWC)

            IF (NOMF(1) > 0) THEN

               SELECT CASE (ABS(IOR0))

                  CASE (1)
                     JDIFF = WF(IWF(2))%IYN(1) - WF(IWF(1))%IYN(1)
                     KDIFF = WF(IWF(3))%IZN(1) - WF(IWF(1))%IZN(1)
                     IF (JDIFF==1 .AND. KDIFF==1) THEN
                        IY1 = WF(IWF(2))%IYN(2)/2
                        IY2 = IY1
                        IZ1 = WF(IWF(3))%IZN(2)/2
                        IZ2 = IZ1
                     ELSE IF (JDIFF==2 .AND. KDIFF==2) THEN
                        IY1 = WF(IWF(1))%IYN(2)/2
                        IY2 = WF(IWF(2))%IYN(2)/2
                        IZ1 = WF(IWF(1))%IZN(2)/2
                        IZ2 = WF(IWF(3))%IZN(2)/2
                     ELSE IF (JDIFF==0 .AND. KDIFF==0) THEN
                        IY1 = WF(IWF(1))%IYN(1)/2
                        IY2 = IY1
                        IZ1 = WF(IWF(1))%IZN(1)/2
                        IZ2 = IZ1
                     ELSE
                        CALL SCARC_SHUTDOWN(NSCARC_ERROR_GRID_RESOLUTION, SCARC_NONE, IOR0)
                     ENDIF

                  CASE (2)
                     IDIFF = WF(IWF(2))%IXN(1) - WF(IWF(1))%IXN(1)
                     KDIFF = WF(IWF(3))%IZN(1) - WF(IWF(1))%IZN(1)
                     IF (IDIFF==1 .AND. KDIFF==1) THEN
                        IX1 = WF(IWF(2))%IXN(2)/2
                        IX2 = IX1
                        IZ1 = WF(IWF(3))%IZN(2)/2
                        IZ2 = IZ1
                     ELSE IF (IDIFF==2 .AND. KDIFF==2) THEN
                        IX1 = WF(IWF(1))%IXN(2)/2
                        IX2 = WF(IWF(2))%IXN(2)/2
                        IZ1 = WF(IWF(1))%IZN(2)/2
                        IZ2 = WF(IWF(3))%IZN(2)/2
                     ELSE IF (IDIFF==0 .AND. KDIFF==0) THEN
                        IX1 = WF(IWF(1))%IXN(2)/2
                        IX2 = IX1
                        IZ1 = WF(IWF(1))%IZN(2)/2
                        IZ2 = IZ1
                     ELSE
                        CALL SCARC_SHUTDOWN(NSCARC_ERROR_GRID_RESOLUTION, SCARC_NONE, IOR0)
                     ENDIF

                  CASE (3)
                     IDIFF = WF(IWF(2))%IXN(1) - WF(IWF(1))%IXN(1)
                     JDIFF = WF(IWF(3))%IYN(1) - WF(IWF(1))%IYN(1)
                     IF (IDIFF==1 .AND. JDIFF==1) THEN
                        IX1 = WF(IWF(2))%IXN(2)/2
                        IX2 = IX1
                        IY1 = WF(IWF(3))%IYN(2)/2
                        IY2 = IY1
                     ELSE IF (IDIFF==2 .AND. JDIFF==2) THEN
                        IX1 = WF(IWF(1))%IXN(2)/2
                        IX2 = WF(IWF(2))%IXN(2)/2
                        IY1 = WF(IWF(1))%IYN(2)/2
                        IY2 = WF(IWF(3))%IYN(2)/2
                     ELSE IF (IDIFF==0 .AND. JDIFF==0) THEN
                        IX1 = WF(IWF(2))%IXN(2)/2
                        IX2 = IX1
                        IY1 = WF(IWF(3))%IYN(2)/2
                        IY2 = IY1
                     ELSE
                        CALL SCARC_SHUTDOWN(NSCARC_ERROR_GRID_RESOLUTION, SCARC_NONE, IOR0)
                     ENDIF

               END SELECT

               SELECT CASE (IOR0)
                  CASE (1)
                     IX1 = MESHES(NOMF(1))%IBAR/IREFINE
                     IX2 = IX1
                  CASE (-1)
                     IX1 = 1
                     IX2 = IX1
                  CASE (2)
                     IY1 = MESHES(NOMF(1))%JBAR/IREFINE
                     IY2 = IY1
                  CASE (-2)
                     IY1 = 1
                     IY2 = IY1
                  CASE (3)
                     IZ1 = MESHES(NOMF(1))%KBAR/IREFINE
                     IZ2 = IZ1
                  CASE (-3)
                     IZ1 = 1
                  IZ2 = IZ1
               END SELECT

               WC(IWC)%IXN(1) = IX1
               WC(IWC)%IYN(1) = IY1
               WC(IWC)%IZN(1) = IZ1
               WC(IWC)%IXN(2) = IX2
               WC(IWC)%IYN(2) = IY2
               WC(IWC)%IZN(2) = IZ2

               CALL SCARC_SETUP_WALL_NEIGHBOR(GC, OGC, IX1, IX2, IY1, IY2, IZ1, IZ2, IWC, NM, NOMF(1), NL)

            ENDIF
         ENDIF
         IWC = IWC + 1
      ENDDO
   ENDDO
ENDDO

END SUBROUTINE SCARC_SETUP_WALL_LEVEL


! -------------------------------------------------------------------------------------------------
!> \brief Setup dimensions for data exchanges
! -------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_EXCHANGE_DIMENSIONS(G, OG, NOM, IREFINE)
USE SCARC_POINTERS, ONLY: GWC
TYPE (SCARC_GRID_TYPE), POINTER, INTENT(IN) :: G, OG
INTEGER, INTENT(IN) :: NOM, IREFINE
INTEGER :: IMIN, IMAX, JMIN, JMAX, KMIN, KMAX, IW
LOGICAL :: FOUND

IMIN=0
IMAX=MESHES(NOM)%IBAR/IREFINE+1
IF (TWO_D) THEN
   JMIN=0
   JMAX=2
ELSE
   JMIN=0
   JMAX=MESHES(NOM)%JBAR/IREFINE+1
ENDIF
KMIN=0
KMAX=MESHES(NOM)%KBAR/IREFINE+1

FOUND = .FALSE.
SEARCH_LOOP: DO IW=1, OG%NCG

   GWC => G%WALL(IW)
   IF (GWC%NOM/=NOM) CYCLE SEARCH_LOOP
   FOUND = .TRUE.

   SELECT CASE (GWC%IOR)
      CASE ( 1)
         IMIN=MAX(IMIN,GWC%IXN(1)-1)
      CASE (-1)
         IMAX=MIN(IMAX,GWC%IXN(2)+1)
      CASE ( 2)
         JMIN=MAX(JMIN,GWC%IYN(1)-1)
      CASE (-2)
         JMAX=MIN(JMAX,GWC%IYN(2)+1)
      CASE ( 3)
         KMIN=MAX(KMIN,GWC%IZN(1)-1)
      CASE (-3)
         KMAX=MIN(KMAX,GWC%IZN(2)+1)
   END SELECT
ENDDO SEARCH_LOOP

N_EXCHANGES = N_EXCHANGES+1

END SUBROUTINE SCARC_SETUP_EXCHANGE_DIMENSIONS


! ----------------------------------------------------------------------------------------------------
!> \brief Allocate several global structures for data exchange
! ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_GLOBALS
INTEGER :: NM, NP

CROUTINE = 'SCARC_SETUP_GLOBALS'

IF (N_MPI_PROCESSES > 1) THEN

   ! Allocate and preset counter and displacement vector for global data exchanges
   CALL SCARC_ALLOCATE_INT1 (COUNTS, 0, N_MPI_PROCESSES-1, NSCARC_INIT_ZERO, 'COUNTS', CROUTINE)
   CALL SCARC_ALLOCATE_INT1 (DISPLS, 0, N_MPI_PROCESSES-1, NSCARC_INIT_ZERO, 'DISPLS', CROUTINE)

   ! Get number of data to send per process
   DO NP = 0, N_MPI_PROCESSES-1
      DO NM = 1, NMESHES
         IF (PROCESS(NM)==NP) COUNTS(NP) = COUNTS(NP) + 1
      ENDDO
   ENDDO

   ! Get displacements on communication vector for all meshes
   DO NP = 1, N_MPI_PROCESSES-1
      DISPLS(NP) = COUNTS(NP-1) + DISPLS(NP-1)
   ENDDO

ENDIF

CALL SCARC_ALLOCATE_INT1 (MESH_INT , 1, NMESHES, NSCARC_INIT_ZERO, 'MESH_INT', CROUTINE)
CALL SCARC_ALLOCATE_REAL1(MESH_REAL, 1, NMESHES, NSCARC_INIT_ZERO, 'MESH_REAL', CROUTINE)

CALL SCARC_SETUP_DIMENSIONS(NLEVEL_MIN)

END SUBROUTINE SCARC_SETUP_GLOBALS


! ----------------------------------------------------------------------------------------------------------
!> \brief Get information about global numbers of unknowns for unstructured discretization
! ----------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_DIMENSIONS(NL)
USE SCARC_POINTERS, ONLY: G
USE SCARC_POINTER_ROUTINES, ONLY: SCARC_POINT_TO_GRID
INTEGER, INTENT(IN) :: NL
INTEGER :: NM, NM2

! Preset communication array MESH_INT with local numbers of cells for all meshes depending on type of discretization
MESHES_LOOP1: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G
   MESH_INT(NM) = G%NC_LOCAL(NM)

!   DO NM2 = LOWER_MESH_INDEX, UPPER_MESH_INDEX
!      MESH_INT(NM2) = G%NC_LOCAL(NM2)
!   ENDDO

ENDDO MESHES_LOOP1


! Broadcast number of local mesh cells on level NL to all and build global sum
IF (N_MPI_PROCESSES > 1) &
   CALL MPI_ALLGATHERV(MPI_IN_PLACE,1,MPI_INTEGER,MESH_INT,COUNTS,DISPLS, MPI_INTEGER,MPI_COMM_WORLD,IERROR)
NC_GLOBAL(NL) = SUM(MESH_INT(1:NMESHES))

! Store information on local and global cells numbers on data structure of corresponding discretization type
MESHES_LOOP2: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G

   G%NC_LOCAL(1:NMESHES) = MESH_INT(1:NMESHES)
   G%NC_GLOBAL = SUM(MESH_INT(1:NMESHES))

   ! compute offset between local grid numberings
   IF (NMESHES > 1) THEN
      DO NM2=2,NMESHES
         G%NC_OFFSET(NM2) = G%NC_OFFSET(NM2-1) + G%NC_LOCAL(NM2-1)
      ENDDO
   ENDIF

ENDDO MESHES_LOOP2

IF (NL == NLEVEL_MIN) THEN
   DO NM = 1, NMESHES
      SCARC(NM)%NC = MESH_INT(NM)
   ENDDO
ENDIF

END SUBROUTINE SCARC_SETUP_DIMENSIONS


! ----------------------------------------------------------------------------------------------------
!> \brief Assign handles to currently used grid type
!  This routine assumes, that L already points to the correct level NL of mesh NL and
!  additionally sets the requested discretization type
! ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_SYSTEM_TYPE(NGRID, NMATRIX)
INTEGER, INTENT(IN) :: NGRID, NMATRIX

SELECT CASE (NMATRIX)
   CASE (NSCARC_MATRIX_POISSON)
      IS_POISSON = .TRUE.
      IS_LAPLACE = .FALSE.
   CASE (NSCARC_MATRIX_LAPLACE)
      IS_POISSON = .FALSE.
      IS_LAPLACE = .TRUE.
END SELECT

SELECT CASE (NGRID)
   CASE (NSCARC_GRID_STRUCTURED)
      PRES_ON_WHOLE_DOMAIN = .TRUE.
      TYPE_GRID = NSCARC_GRID_STRUCTURED
      IS_STRUCTURED   = .TRUE.
      IS_UNSTRUCTURED = .FALSE.
   CASE (NSCARC_GRID_UNSTRUCTURED)
      PRES_ON_WHOLE_DOMAIN = .FALSE.
      TYPE_GRID = NSCARC_GRID_UNSTRUCTURED
      IS_STRUCTURED   = .FALSE.
      IS_UNSTRUCTURED = .TRUE.
END SELECT
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) '========> SETUP_GRID_TYPE : TYPE_GRID, PRES_ON_WHOLE_DOMAIN:', TYPE_GRID, PRES_ON_WHOLE_DOMAIN
#endif

END SUBROUTINE SCARC_SETUP_SYSTEM_TYPE

! ----------------------------------------------------------------------------------------------------
!> \brief Setup system of equations (Poisson matrix + BC's) for different variants of ScaRC
! Define matrix stencils and initialize matrices and boundary conditions on all needed levels
! ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_SYSTEMS
INTEGER :: NM, NL, TYPE_SCOPE_BACKUP, TYPE_GRID_BACKUP
  
CROUTINE = 'SCARC_SETUP_SYSTEMS'

! ------ Setup sizes for system matrices
  
SELECT_SCARC_METHOD_SIZES: SELECT CASE (TYPE_METHOD)

   ! -------- Global Krylov method

   CASE (NSCARC_METHOD_KRYLOV)
   
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'TTTTT: SETUP_SYSTEM:A: TYPE_SCOPE:', TYPE_SCOPE(0)
#endif
      CALL SCARC_SETUP_GRID_TYPE (TYPE_GRID)                      ! process specified discretization type
      CALL SCARC_SETUP_POISSON_SIZES (NLEVEL_MIN)                  ! setup sizes on finest level
   
      IF (HAS_TWO_LEVELS .AND. .NOT.HAS_AMG_LEVELS) &
         CALL SCARC_SETUP_POISSON_SIZES (NLEVEL_MAX)               ! twolevel-precon: also setup size for coarse level
   
      IF (IS_CG_GMG) THEN                                                   
         DO NL=NLEVEL_MIN+1, NLEVEL_MAX
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'TTTTT: SETUP_SYSTEM:B: TYPE_SCOPE:', TYPE_SCOPE(0)
#endif
            CALL SCARC_SETUP_POISSON_SIZES (NL)                    ! GMG-precon: also setup size for all other levels
         ENDDO
      ENDIF
   
   ! -------- Global Multigrid method

   CASE (NSCARC_METHOD_MULTIGRID)
   
      CALL SCARC_SETUP_GRID_TYPE (TYPE_GRID)                      ! process specified discretization type
      SELECT CASE (TYPE_MULTIGRID)
         CASE (NSCARC_MULTIGRID_GEOMETRIC)                                   
            DO NL=NLEVEL_MIN, NLEVEL_MAX
               CALL SCARC_SETUP_POISSON_SIZES (NL)                 ! GMG: setup size for all levels
            ENDDO
         CASE (NSCARC_MULTIGRID_ALGEBRAIC)                                   
            CALL SCARC_SETUP_POISSON_SIZES (NLEVEL_MIN)            ! AMG: setup sizes only on finest level
      END SELECT
   
   ! -------- Global MGM method - currently just proof of concept

   CASE (NSCARC_METHOD_MGM)
   
      CALL SCARC_SETUP_GRID_TYPE (NSCARC_GRID_STRUCTURED)         ! First process structured discretization
      CALL SCARC_SETUP_POISSON_SIZES (NLEVEL_MIN)        
   
      CALL SCARC_SETUP_GRID_TYPE (NSCARC_GRID_UNSTRUCTURED)       ! Then process unstructured discretization
      IF (SCARC_MGM_CHECK_LAPLACE .OR. SCARC_MGM_INIT_EXACT) &
         CALL SCARC_SETUP_POISSON_SIZES (NLEVEL_MIN)              ! ... for global Poisson matrix
      CALL SCARC_SETUP_MGM_LAPLACE_SIZES (NLEVEL_MIN)             ! ... for local Laplace matrices
   
END SELECT SELECT_SCARC_METHOD_SIZES


  
! ------ Assemble system matrices on requested grid levels and set boundary conditions
  
MESHES_POISSON_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   SELECT_SCARC_METHOD: SELECT CASE (TYPE_METHOD)

      ! ---------- Krylov method (CG) as main solver, different preconditioners possible

      CASE (NSCARC_METHOD_KRYLOV)

         ! For all different possible Krylov variants, first setup Poisson matrix on finest level including BC's 

         CALL SCARC_SETUP_POISSON (NM, NLEVEL_MIN)
         CALL SCARC_SETUP_BOUNDARY(NM, NLEVEL_MIN)

         ! Depending on the requested preconditioner, also assemble the Poisson matrix with BC's on specific coarser levels

         SELECT_KRYLOV_PRECON: SELECT CASE (TYPE_PRECON)

            ! In case of multigrid as preconditioner:
            ! only build higher level structures in case of geometric multigrid (algebraic variant is done elsewhere)

            CASE (NSCARC_RELAX_MULTIGRID)

               IF (IS_CG_GMG) THEN
                  DO NL = NLEVEL_MIN+1, NLEVEL_MAX
                     CALL SCARC_SETUP_POISSON (NM, NL)
                     CALL SCARC_SETUP_BOUNDARY(NM, NL)
                  ENDDO
               ENDIF

#ifdef WITH_MKL
            ! In case of LU-decomposition as preconditioner
            ! locally acting: PARDISO from MKL as preconditioners on fine level with possible coarse grid correction

            CASE (NSCARC_RELAX_MKL)

               IF (TYPE_SCOPE(1) == NSCARC_SCOPE_LOCAL .AND. HAS_TWO_LEVELS .AND. .NOT.HAS_AMG_LEVELS) THEN
                  CALL SCARC_SETUP_POISSON (NM, NLEVEL_MAX)
                  CALL SCARC_SETUP_BOUNDARY(NM, NLEVEL_MAX)
               ENDIF
#endif

            ! in case of default preconditioners (JACOBI/SSOR/FFT/...):
            ! if there is an additional coarse grid correction which is NOT AMG-based, 
            ! then also assemble matrix on coarse grid level

            CASE DEFAULT
   
               IF (HAS_TWO_LEVELS .AND. .NOT.HAS_AMG_LEVELS) THEN
                  CALL SCARC_SETUP_POISSON (NM, NLEVEL_MAX)
                  CALL SCARC_SETUP_BOUNDARY(NM, NLEVEL_MAX)
               ENDIF

         END SELECT SELECT_KRYLOV_PRECON


      ! ---------- Multigrid as main solver

      CASE (NSCARC_METHOD_MULTIGRID)

         ! For all different possible multigrid-variants, first setup Poisson matrix on finest level including BC's 

         CALL SCARC_SETUP_POISSON (NM, NLEVEL_MIN)
         CALL SCARC_SETUP_BOUNDARY(NM, NLEVEL_MIN)

         ! On case of a  geometric multigrid, assemble standard n-point-matrix hierarchy on all coarser levels, too
         ! Note: in case of an algebraic multigrid, this will be done in a separate routine later

         IF (TYPE_MULTIGRID == NSCARC_MULTIGRID_GEOMETRIC) THEN
            DO NL = NLEVEL_MIN + 1, NLEVEL_MAX
               CALL SCARC_SETUP_POISSON (NM, NL)
               CALL SCARC_SETUP_BOUNDARY(NM, NL)
            ENDDO
         ENDIF


      ! ---------- McKenny-Greengard-Mayo method:
      ! Solving for the structured and unstructured Poisson matrix
      ! Assemble both, the structured and unstructured Poisson matrix
      ! temporarily they will be stored separately in matrices AC and ACU due to the different
      ! settings along internal boundary cells,
      ! in the medium term, a toggle mechanism will be implemented which only switches the corresponding
      ! entries while keeping the entries which are the same for both discretization types

      CASE (NSCARC_METHOD_MGM)
   
         TYPE_SCOPE(0) = NSCARC_SCOPE_GLOBAL

         ! Then assemble structured matrix with inhomogeneous boundary conditions

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) '========================= MGM: STRUCTURED'
#endif
         CALL SCARC_SETUP_GRID_TYPE (NSCARC_GRID_STRUCTURED)
         CALL SCARC_SETUP_POISSON (NM, NLEVEL_MIN)
         CALL SCARC_SETUP_BOUNDARY(NM, NLEVEL_MIN)

         ! First assemble unstructured matrix with homogeneous Dirichlet boundary conditions

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) '========================= MGM: UNSTRUCTURED'
#endif
         CALL SCARC_SETUP_GRID_TYPE (NSCARC_GRID_UNSTRUCTURED)
         CALL SCARC_SETUP_POISSON (NM, NLEVEL_MIN)
         CALL SCARC_SETUP_BOUNDARY(NM, NLEVEL_MIN)

         TYPE_SCOPE(0) = NSCARC_SCOPE_LOCAL
         !CALL SCARC_SETUP_LAPLACE (NM, NLEVEL_MIN)
         CALL SCARC_SETUP_MGM_LAPLACE (NM, NLEVEL_MIN)
         CALL SCARC_SETUP_MGM_BOUNDARY(NM, NLEVEL_MIN) 

#ifdef WITH_MKL
         TYPE_MKL(NLEVEL_MIN) = NSCARC_MKL_LOCAL
#endif

         TYPE_SCOPE(0) = NSCARC_SCOPE_GLOBAL


   END SELECT SELECT_SCARC_METHOD

ENDDO MESHES_POISSON_LOOP

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) '========================= SUSI: CELL_MAPPING', TYPE_SCOPE(0), TYPE_MATVEC, IS_MGM
#endif

! Setup mappings for the global numbering of vectors and the Poisson matrix (compact storage technique only)
 
IF (TYPE_MATRIX == NSCARC_MATRIX_COMPACT) THEN
   IF (IS_MGM) THEN
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) '========================= SUSI: CELL_MAPPING, UNSTRUCTURED'
#endif
      TYPE_SCOPE = NSCARC_SCOPE_GLOBAL
      !TYPE_SCOPE = NSCARC_SCOPE_LOCAL
      CALL SCARC_SETUP_GRID_TYPE (NSCARC_GRID_UNSTRUCTURED)
      CALL SCARC_SETUP_GLOBAL_CELL_MAPPING(NLEVEL_MIN)
      CALL SCARC_SETUP_GLOBAL_POISSON_COLUMNS(NLEVEL_MIN)
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) '========================= SUSI: CELL_MAPPING, STRUCTURED'
#endif
      TYPE_SCOPE = NSCARC_SCOPE_GLOBAL
      CALL SCARC_SETUP_GRID_TYPE (NSCARC_GRID_STRUCTURED)
      CALL SCARC_SETUP_GLOBAL_CELL_MAPPING(NLEVEL_MIN)
      CALL SCARC_SETUP_GLOBAL_POISSON_COLUMNS(NLEVEL_MIN)
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) '========================= SUSI: CELL_MAPPING, FINISH'
#endif
   ELSE
      CALL SCARC_SETUP_GLOBAL_CELL_MAPPING(NLEVEL_MIN)
      CALL SCARC_SETUP_GLOBAL_POISSON_COLUMNS(NLEVEL_MIN)
   ENDIF
ENDIF
 
! If there is more than one mesh, exchange matrix values in overlapping parts
! This must be done for all multilevel methods at least at the finest grid level
! Furthermore also at all higher levels except for the AMG method,
! in this case it will be done later in routine SETUP_ALGEBRAIC_MULTIGRID

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) '========================= SUSI: OVERLAPS'
#endif
IF (SCARC_GET_MATRIX_TYPE(NLEVEL_MIN) == NSCARC_MATRIX_COMPACT) CALL SCARC_SETUP_GLOBAL_POISSON_OVERLAPS(NLEVEL_MIN)
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) '========================= SUSI: OVERLAPS'
#endif

 MULTI_LEVEL_IF: IF (HAS_MULTIPLE_LEVELS .AND. .NOT.HAS_AMG_LEVELS) THEN

   DO NL = NLEVEL_MIN+1, NLEVEL_MAX
      IF (SCARC_GET_MATRIX_TYPE(NL) /= NSCARC_MATRIX_COMPACT) CYCLE
      CALL SCARC_SETUP_GLOBAL_CELL_MAPPING(NL)
      CALL SCARC_SETUP_GLOBAL_POISSON_COLUMNS(NL)
      CALL SCARC_SETUP_GLOBAL_POISSON_OVERLAPS(NL)
   ENDDO 

ENDIF MULTI_LEVEL_IF

  
! ------ IF MKL-solver is used on specific levels, then setup symmetric Poisson matrix there
  
#ifdef WITH_MKL
IF (IS_MGM) THEN
   TYPE_SCOPE_BACKUP = TYPE_SCOPE(0)
   TYPE_SCOPE(0) = NSCARC_SCOPE_LOCAL
   CALL SCARC_SETUP_GRID_TYPE (NSCARC_GRID_UNSTRUCTURED)        
ENDIF

MESHES_MKL_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   IF (TYPE_MKL(NLEVEL_MIN) /= NSCARC_MKL_NONE) CALL SCARC_SETUP_POISSON_MKL(NM, NLEVEL_MIN)
   IF (HAS_GMG_LEVELS) THEN
      DO NL = NLEVEL_MIN+1, NLEVEL_MAX
         IF (TYPE_MKL(NL) /= NSCARC_MKL_NONE) CALL SCARC_SETUP_POISSON_MKL(NM, NL)
      ENDDO
   ENDIF
ENDDO MESHES_MKL_LOOP

IF (IS_MGM) THEN
   TYPE_SCOPE(0) = TYPE_SCOPE_BACKUP 
   CALL SCARC_SETUP_GRID_TYPE (TYPE_GRID_BACKUP)        
ENDIF

#endif


! Debug matrix and wall structures - only if directive SCARC_DEBUG is set

#ifdef WITH_SCARC_DEBUG
CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_WALL  , NLEVEL_MIN, 'WALL')
CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_FACE  , NLEVEL_MIN, 'FACE_AFTER_SYSTEM')
#endif

END SUBROUTINE SCARC_SETUP_SYSTEMS


! -------------------------------------------------------------------------------------------
!> \brief Setup mapping from local to global cell numbering
! -------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_GLOBAL_CELL_MAPPING(NL)
USE SCARC_POINTERS, ONLY : G
USE SCARC_POINTER_ROUTINES, ONLY: SCARC_POINT_TO_GRID
INTEGER, INTENT(IN) :: NL
INTEGER :: NM, NOM, IC, IW, ICE, ICN

!IF (NMESHES == 1) RETURN
CROUTINE = 'SCARC_SETUP_GLOBAL_CELL_MAPPING'

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'SETUP_GLOBAL_CELL_MAPPING:1:', TYPE_GRID, TYPE_SCOPE(0), TYPE_MATVEC
#endif
!IF (IS_MGM .AND. IS_UNSTRUCTURED) RETURN
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'SETUP_GLOBAL_CELL_MAPPING:2'
#endif

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID (NM, NL)      

   IF (HAS_AMG_LEVELS) THEN
      CALL SCARC_ALLOCATE_INT1 (G%LOCAL_TO_GLOBAL, 1, G%NCE2, NSCARC_INIT_ZERO, 'G%LOCAL_TO_GLOBAL', CROUTINE)  
   ELSE
      CALL SCARC_ALLOCATE_INT1 (G%LOCAL_TO_GLOBAL, 1, G%NCE , NSCARC_INIT_ZERO, 'G%LOCAL_TO_GLOBAL', CROUTINE) 
   ENDIF

   DO IC = 1, G%NC
      G%LOCAL_TO_GLOBAL(IC) = IC + G%NC_OFFSET(NM)
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'LOCAL_TO_GLOBAL(', IC,')=', G%LOCAL_TO_GLOBAL(IC), G%NCE
#endif
   ENDDO

   DO IW = 1, G%NW
      NOM = G%WALL(IW)%NOM
      IF (NOM == 0) CYCLE
      ICE = G%WALL(IW)%ICE
      ICN = G%ICE_TO_ICN(ICE)
      G%LOCAL_TO_GLOBAL(ICE) = ICN + G%NC_OFFSET(NOM)
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'LOCAL_TO_GLOBAL(', ICE,')=', G%LOCAL_TO_GLOBAL(ICE), IW, ICE, ICN, NOM
#endif
   ENDDO

ENDDO

END SUBROUTINE SCARC_SETUP_GLOBAL_CELL_MAPPING


! -------------------------------------------------------------------------------------------
!> \brief Get global numberings for compact column vector of Poisson matrix 
! -------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_GLOBAL_POISSON_COLUMNS(NL)
USE SCARC_POINTERS, ONLY : G, A
USE SCARC_POINTER_ROUTINES, ONLY: SCARC_POINT_TO_GRID, SCARC_POINT_TO_CMATRIX
INTEGER, INTENT(IN) :: NL
INTEGER :: NM, IC, ICOL, JC

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'GLOBAL_POISSON_COLUMNS: TYPE_SCOPE(0)=', TYPE_SCOPE(0)
#endif

IF (NMESHES == 1 .OR. TYPE_SCOPE(0) == NSCARC_SCOPE_LOCAL) THEN
   DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
      CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G
      A => SCARC_POINT_TO_CMATRIX (G, NSCARC_MATRIX_POISSON)
      A%COLG = A%COL
   ENDDO
ELSE
   DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
      CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G
      A => SCARC_POINT_TO_CMATRIX (G, NSCARC_MATRIX_POISSON)
      DO IC = 1, G%NC
         DO ICOL = A%ROW(IC), A%ROW(IC+1)-1
            JC = A%COL(ICOL)
            A%COLG(ICOL) = G%LOCAL_TO_GLOBAL(JC)
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'A%COLG(', ICOL,')=', IC, ICOL, JC, A%COLG(ICOL)
#endif
         ENDDO
      ENDDO
   ENDDO
ENDIF

END SUBROUTINE SCARC_SETUP_GLOBAL_POISSON_COLUMNS


! -------------------------------------------------------------------------------------------
!> \brief Make Poisson matrix global by exchanging adjacent overlaps
! -------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_GLOBAL_POISSON_OVERLAPS(NL)
USE SCARC_POINTERS, ONLY: S, G, OG, A, OA
USE SCARC_POINTER_ROUTINES, ONLY: SCARC_POINT_TO_GRID, SCARC_POINT_TO_OTHER_GRID, &
                                  SCARC_POINT_TO_CMATRIX, SCARC_POINT_TO_OTHER_CMATRIX
INTEGER, INTENT(IN) :: NL
INTEGER :: NM, INBR, NOM

IF (NMESHES == 1 .OR. TYPE_SCOPE(0) == NSCARC_SCOPE_LOCAL) RETURN

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'SETUP_GLOBAL_POISSON_OVERLAPS:A:', TYPE_SCOPE(0), TYPE_MATVEC
#endif
CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_MATRIX_COLS,  NSCARC_MATRIX_POISSON, NL)
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'SETUP_GLOBAL_POISSON_OVERLAPS:B:', TYPE_SCOPE(0), TYPE_MATVEC
#endif
CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_MATRIX_COLSG, NSCARC_MATRIX_POISSON, NL)
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'SETUP_GLOBAL_POISSON_OVERLAPS:C:', TYPE_SCOPE(0), TYPE_MATVEC
#endif
CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_MATRIX_VALS,  NSCARC_MATRIX_POISSON, NL)

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'SETUP_GLOBAL_POISSON_OVERLAPS:D:', TYPE_SCOPE(0), TYPE_MATVEC
#endif
CALL SCARC_EXTRACT_MATRIX_OVERLAPS(NSCARC_MATRIX_POISSON, 1, NL)

MESHES_FINE_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID (NM, NL)    
   A => SCARC_POINT_TO_CMATRIX (G, NSCARC_MATRIX_POISSON)
   CALL SCARC_REDUCE_CMATRIX(A, 'G%POISSON', CROUTINE)

   OMESHES_FINE_LOOP: DO INBR = 1, S%N_NEIGHBORS
      NOM = S%NEIGHBORS(INBR)
      CALL SCARC_POINT_TO_OTHER_GRID(NM, NOM, NL)
      OA => SCARC_POINT_TO_OTHER_CMATRIX (OG, NSCARC_MATRIX_POISSON)
      CALL SCARC_REDUCE_CMATRIX(OA, 'OG%POISSON', CROUTINE)
   ENDDO OMESHES_FINE_LOOP

ENDDO MESHES_FINE_LOOP
    
END SUBROUTINE SCARC_SETUP_GLOBAL_POISSON_OVERLAPS


! ------------------------------------------------------------------------------------------------
!> \brief Check if specified cell is within a given mesh
! ------------------------------------------------------------------------------------------------
LOGICAL FUNCTION SCARC_CELL_WITHIN_MESH(G, NM, IC)
TYPE (SCARC_GRID_TYPE), POINTER, INTENT(IN) :: G
INTEGER, INTENT(IN) :: NM, IC
INTEGER :: IC_START, IC_STOP

SCARC_CELL_WITHIN_MESH = .FALSE.
IC_START = G%NC_OFFSET(NM) + 1
IF (NM < NMESHES) THEN
   IC_STOP  = G%NC_OFFSET(NM+1)
ELSE
   IC_STOP  = G%NC_GLOBAL
ENDIF
IF (IC_START <=  IC .AND. IC <= IC_STOP) SCARC_CELL_WITHIN_MESH = .TRUE.
RETURN

END FUNCTION SCARC_CELL_WITHIN_MESH


! ------------------------------------------------------------------------------------------------
!> \brief Allocate Poisson matrix for the usual 5-point-stencil (2D) or 7-point-stencil (3D)
! Compact storage technique (POISSON)
!    Compression technique to store sparse matrices, non-zero entries are stored
!    in a 1D-vector B(.), row after row,
!    Each row starts with its diagonal entry followed by the other non-zero entries
!    In order to identify each element, pointer arrays ROW and COL are needed,
!    ROW points to the several diagonal entries in vector B(.),
!    COL points to the columns which non-zero entries in the matrix stencil
! Bandwise storage technique (POISSONB)
!    explanation to come ...
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_POISSON (NM, NL)
USE SCARC_POINTERS, ONLY: S, L, G, OG, A, AB, OA, OAB
USE SCARC_POINTER_ROUTINES, ONLY: SCARC_POINT_TO_GRID,    SCARC_POINT_TO_OTHER_GRID, &
                                  SCARC_POINT_TO_CMATRIX, SCARC_POINT_TO_OTHER_CMATRIX, &
                                  SCARC_POINT_TO_BMATRIX, SCARC_POINT_TO_OTHER_BMATRIX
INTEGER, INTENT(IN) :: NM, NL
INTEGER :: IX, IY, IZ, IC, IP, INBR, NOM

CROUTINE = 'SCARC_SETUP_POISSON'
 
! Compute single matrix entries and corresponding row and column pointers
! Along internal boundaries use placeholders for the neighboring matrix entries
! which will be communicated in a following step
 
SELECT_STORAGE_TYPE: SELECT CASE (SCARC_GET_MATRIX_TYPE(NL))

 
   ! ---------- COMPACT Storage technique
 
   CASE (NSCARC_MATRIX_COMPACT)
   
      ! Allocate main matrix on non-overlapping part of mesh

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'SETUP_POISSON: TYPE_SCOPE:', TYPE_SCOPE(0)
#endif
      CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G
      A => SCARC_POINT_TO_CMATRIX (G, NSCARC_MATRIX_POISSON)
      CALL SCARC_ALLOCATE_CMATRIX (A, NL, NSCARC_PRECISION_DOUBLE, NSCARC_MATRIX_FULL, 'G%POISSON', CROUTINE)

      ! For every neighbor allocate small matrix on overlapping part of mesh

      DO INBR = 1, SCARC(NM)%N_NEIGHBORS
         NOM = S%NEIGHBORS(INBR)
         CALL SCARC_POINT_TO_OTHER_GRID(NM, NOM, NL)
         OA => SCARC_POINT_TO_OTHER_CMATRIX (OG, NSCARC_MATRIX_POISSON)
         CALL SCARC_ALLOCATE_CMATRIX(OA, NL, NSCARC_PRECISION_DOUBLE, NSCARC_MATRIX_FULL, 'OG%POISSON', CROUTINE)
      ENDDO

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'SETUP_POISSON:B: TYPE_SCOPE:', TYPE_SCOPE(0)
#endif
      IP = 1
      DO IZ = 1, L%NZ
         DO IY = 1, L%NY
            DO IX = 1, L%NX
   
               IF (IS_UNSTRUCTURED .AND. L%IS_SOLID(IX, IY, IZ)) CYCLE
               IC = G%CELL_NUMBER(IX, IY, IZ)

               ! Main diagonal 

               CALL SCARC_SETUP_MAINDIAG (IC, IX, IY, IZ, IP)
   
               ! Lower subdiagonals

               IF (IS_VALID_DIRECTION(IX, IY, IZ,  3)) CALL SCARC_SETUP_SUBDIAG(IC, IX, IY, IZ, IX  , IY  , IZ-1, IP,  3)
               IF (IS_VALID_DIRECTION(IX, IY, IZ,  2)) CALL SCARC_SETUP_SUBDIAG(IC, IX, IY, IZ, IX  , IY-1, IZ  , IP,  2)
               IF (IS_VALID_DIRECTION(IX, IY, IZ,  1)) CALL SCARC_SETUP_SUBDIAG(IC, IX, IY, IZ, IX-1, IY  , IZ  , IP,  1)
   
               ! Upper subdiagonals

               IF (IS_VALID_DIRECTION(IX, IY, IZ, -1)) CALL SCARC_SETUP_SUBDIAG(IC, IX, IY, IZ, IX+1, IY  , IZ  , IP, -1)
               IF (IS_VALID_DIRECTION(IX, IY, IZ, -2)) CALL SCARC_SETUP_SUBDIAG(IC, IX, IY, IZ, IX  , IY+1, IZ  , IP, -2)
               IF (IS_VALID_DIRECTION(IX, IY, IZ, -3)) CALL SCARC_SETUP_SUBDIAG(IC, IX, IY, IZ, IX  , IY  , IZ+1, IP, -3)
   
            ENDDO
         ENDDO
      ENDDO
   
      A%ROW(G%NC+1) = IP
      A%N_VAL = IP
   
      CALL SCARC_GET_MATRIX_STENCIL_MAX(A, G%NC)

#ifdef WITH_SCARC_DEBUG
CALL SCARC_DEBUG_CMATRIX (A, 'POISSON', 'SETUP_POISSON: NO BDRY')
#endif
 
   ! ---------- bandwise storage technique
 
   CASE (NSCARC_MATRIX_BANDWISE)
   
      ! Allocate main matrix on non-overlapping part of mesh

      CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G
      AB => SCARC_POINT_TO_BMATRIX (G, NSCARC_MATRIX_POISSON)
      CALL SCARC_ALLOCATE_BMATRIX(AB, NL, 'G%POISSONB', CROUTINE)
   
      ! For every neighbor allocate little matrix on overlapping part of mesh

      DO INBR = 1, SCARC(NM)%N_NEIGHBORS
         NOM = S%NEIGHBORS(INBR)
         CALL SCARC_POINT_TO_OTHER_GRID(NM, NOM, NL)
         OAB => SCARC_POINT_TO_BMATRIX (G, NSCARC_MATRIX_POISSON)
         CALL SCARC_ALLOCATE_BMATRIX(OAB, NL, 'OG%POISSONB', CROUTINE)
      ENDDO
   
      IP  = 1
      DO IZ = 1, L%NZ
         DO IY = 1, L%NY
            DO IX = 1, L%NX
   
               IF (IS_UNSTRUCTURED .AND. L%IS_SOLID(IX, IY, IZ)) CYCLE
               IC = G%CELL_NUMBER(IX, IY, IZ)
   
               ! Lower subdiagonals

               IF (IS_VALID_DIRECTION(IX, IY, IZ,  3)) CALL SCARC_SETUP_SUBDIAGB(IC, IX, IY, IZ, IX  , IY  , IZ-1,  3)
               IF (IS_VALID_DIRECTION(IX, IY, IZ,  2)) CALL SCARC_SETUP_SUBDIAGB(IC, IX, IY, IZ, IX  , IY-1, IZ  ,  2)
               IF (IS_VALID_DIRECTION(IX, IY, IZ,  1)) CALL SCARC_SETUP_SUBDIAGB(IC, IX, IY, IZ, IX-1, IY  , IZ  ,  1)
   
               ! Main diagonal

               CALL SCARC_SETUP_MAINDIAGB (IC, IX, IY, IZ)

               ! Upper subdiagonals

               IF (IS_VALID_DIRECTION(IX, IY, IZ, -1)) CALL SCARC_SETUP_SUBDIAGB(IC, IX, IY, IZ, IX+1, IY  , IZ  , -1)
               IF (IS_VALID_DIRECTION(IX, IY, IZ, -2)) CALL SCARC_SETUP_SUBDIAGB(IC, IX, IY, IZ, IX  , IY+1, IZ  , -2)
               IF (IS_VALID_DIRECTION(IX, IY, IZ, -3)) CALL SCARC_SETUP_SUBDIAGB(IC, IX, IY, IZ, IX  , IY  , IZ+1, -3)
   
            ENDDO
         ENDDO
      ENDDO
   
END SELECT SELECT_STORAGE_TYPE

END SUBROUTINE SCARC_SETUP_POISSON


! ---------------------------------------------------------------------------------------------------------------
!> \brief Setup local Laplace matrices 
! ---------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_LAPLACE (NM, NL)
USE SCARC_POINTERS, ONLY: L, G, A
USE SCARC_POINTER_ROUTINES, ONLY: SCARC_POINT_TO_GRID, SCARC_POINT_TO_CMATRIX
INTEGER, INTENT(IN) :: NM, NL
INTEGER :: IX, IY, IZ, IC, IP

CROUTINE = 'SCARC_SETUP_LAPLACE'
 
! Allocate main matrix on non-overlapping part of mesh

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'SETUP_POISSON: TYPE_SCOPE:', TYPE_SCOPE(0)
#endif
CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G
A => SCARC_POINT_TO_CMATRIX (G, NSCARC_MATRIX_LAPLACE)
CALL SCARC_ALLOCATE_CMATRIX (A, NL, NSCARC_PRECISION_DOUBLE, NSCARC_MATRIX_FULL, 'G%POISSON', CROUTINE)

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'SETUP_POISSON:B: TYPE_SCOPE:', TYPE_SCOPE(0)
#endif
IP = 1
DO IZ = 1, L%NZ
   DO IY = 1, L%NY
      DO IX = 1, L%NX

         IF (IS_UNSTRUCTURED .AND. L%IS_SOLID(IX, IY, IZ)) CYCLE
         IC = G%CELL_NUMBER(IX, IY, IZ)

         ! Main diagonal 

         CALL SCARC_SETUP_MAINDIAG (IC, IX, IY, IZ, IP)

         ! Lower subdiagonals

         IF (IS_VALID_DIRECTION(IX, IY, IZ,  3)) CALL SCARC_SETUP_SUBDIAG(IC, IX, IY, IZ, IX  , IY  , IZ-1, IP,  3)
         IF (IS_VALID_DIRECTION(IX, IY, IZ,  2)) CALL SCARC_SETUP_SUBDIAG(IC, IX, IY, IZ, IX  , IY-1, IZ  , IP,  2)
         IF (IS_VALID_DIRECTION(IX, IY, IZ,  1)) CALL SCARC_SETUP_SUBDIAG(IC, IX, IY, IZ, IX-1, IY  , IZ  , IP,  1)

         ! Upper subdiagonals

         IF (IS_VALID_DIRECTION(IX, IY, IZ, -1)) CALL SCARC_SETUP_SUBDIAG(IC, IX, IY, IZ, IX+1, IY  , IZ  , IP, -1)
         IF (IS_VALID_DIRECTION(IX, IY, IZ, -2)) CALL SCARC_SETUP_SUBDIAG(IC, IX, IY, IZ, IX  , IY+1, IZ  , IP, -2)
         IF (IS_VALID_DIRECTION(IX, IY, IZ, -3)) CALL SCARC_SETUP_SUBDIAG(IC, IX, IY, IZ, IX  , IY  , IZ+1, IP, -3)

      ENDDO
   ENDDO
ENDDO
   
A%ROW(G%NC+1) = IP
A%N_VAL = IP
   
CALL SCARC_GET_MATRIX_STENCIL_MAX(A, G%NC)

#ifdef WITH_SCARC_DEBUG
CALL SCARC_DEBUG_CMATRIX (A, 'LAPLACE', 'SETUP_LAPLACE: NO BDRY')
#endif
 
END SUBROUTINE SCARC_SETUP_LAPLACE

! ------------------------------------------------------------------------------------------------
!> \brief Set main diagonal entry for Poisson matrix in compact storage technique
! These values correspond to the full matrix of the global problem
! In case of an equidistant grid, we get the usual 5-point (2d) and 7-point (3d) stencil
! If two meshes with different step sizes meet, we get a weighted stencil along internal wall cells
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_MAINDIAG (IC, IX, IY, IZ, IP)
USE SCARC_POINTERS, ONLY: L, A
INTEGER, INTENT(IN) :: IC, IX, IY, IZ
INTEGER, INTENT(INOUT) :: IP

A%VAL(IP) = - 2.0_EB/(L%DXL(IX-1)*L%DXL(IX))
IF (.NOT.TWO_D) A%VAL(IP) = A%VAL(IP) - 2.0_EB/(L%DYL(IY-1)*L%DYL(IY))
A%VAL(IP) = A%VAL(IP) - 2.0_EB/(L%DZL(IZ-1)*L%DZL(IZ))

A%ROW(IC) = IP
A%COL(IP) = IC

A%STENCIL(0) = A%VAL(IP)

IP = IP + 1
END SUBROUTINE SCARC_SETUP_MAINDIAG


! ------------------------------------------------------------------------------------------------
!> \brief Set subdigonal entries for Poisson matrix in compact storage technique on specified face
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_SUBDIAG (IC, IX1, IY1, IZ1, IX2, IY2, IZ2, IP, IOR0)
USE SCARC_POINTERS, ONLY: L, F, G, A
INTEGER, INTENT(IN) :: IC, IX1, IY1, IZ1, IX2, IY2, IZ2, IOR0
INTEGER, INTENT(INOUT) :: IP
INTEGER :: IW
LOGICAL :: IS_INTERNAL_CELL

! Decide wheter cell is interior or exterior cell

F => L%FACE(IOR0)
SELECT CASE (IOR0)
   CASE ( 1)
      IS_INTERNAL_CELL = IX1 > 1
   CASE (-1)
      IS_INTERNAL_CELL = IX1 < F%NOP
   CASE ( 2)
      IS_INTERNAL_CELL = IY1 > 1
   CASE (-2)
      IS_INTERNAL_CELL = IY1 < F%NOP
   CASE ( 3)
      IS_INTERNAL_CELL = IZ1 > 1
   CASE (-3)
      IS_INTERNAL_CELL = IZ1 < F%NOP
END SELECT

! If IC is an internal cell of the mesh, compute usual matrix contribution for corresponding subdiagonal
IF (IS_INTERNAL_CELL) THEN

   IF (IS_STRUCTURED .OR. .NOT.L%IS_SOLID(IX2, IY2, IZ2)) THEN
      A%VAL(IP) = A%VAL(IP) + F%INCR_INSIDE
      A%COL(IP) = G%CELL_NUMBER(IX2, IY2, IZ2)
      A%STENCIL(-IOR0) = A%VAL(IP)
      IP = IP + 1
   ELSE
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'IX1, IY1, IZ1, IX2, IY2, IZ2, L%IS_SOLID(IX2, IY2, IZ2):', &
                       IX1, IY1, IZ1, IX2, IY2, IZ2, L%IS_SOLID(IX2, IY2, IZ2)
#endif
      CALL SCARC_SHUTDOWN(NSCARC_ERROR_MATRIX_SUBDIAG, SCARC_NONE, NSCARC_NONE)
   ENDIF

! If IC is a boundary cell of the mesh, compute matrix contribution only if there is a neighbor for that cell

ELSE IF (TYPE_SCOPE(0) == NSCARC_SCOPE_GLOBAL .AND. L%FACE(IOR0)%N_NEIGHBORS /= 0) THEN

   IW = SCARC_ASSIGN_SUBDIAG_TYPE (IC, IOR0)           ! get IW of a possibly suitable neighbor at face IOR0
   IF (IW > 0) then                                    ! if available, build corresponding subdiagonal entry
      A%VAL(IP) = A%VAL(IP) + F%INCR_FACE
      A%COL(IP) = G%WALL(IW)%ICE                       ! store its extended number in matrix column pointers
      A%STENCIL(-IOR0) = A%VAL(IP)
      IP = IP + 1
   ENDIF

ENDIF

END SUBROUTINE SCARC_SETUP_SUBDIAG



! ------------------------------------------------------------------------------------------------
!> \brief Determine if cell has a neighbor and, if yes, return corresponding wall cell index
! ------------------------------------------------------------------------------------------------
INTEGER FUNCTION SCARC_ASSIGN_SUBDIAG_TYPE (IC, IOR0)
USE SCARC_POINTERS, ONLY: L, G, F, GWC
INTEGER, INTENT(IN) :: IC, IOR0
INTEGER :: IXW, IYW, IZW
INTEGER :: IXG, IYG, IZG
INTEGER :: IW

SCARC_ASSIGN_SUBDIAG_TYPE = -1

F => L%FACE(IOR0)
SEARCH_WALL_CELLS_LOOP: DO IW = F%NCW0, F%NCW0 + F%NCW - 1

   GWC => G%WALL(IW)

   IF (GWC%NOM == 0) CYCLE
   IXW = GWC%IXW
   IYW = GWC%IYW
   IZW = GWC%IZW

   IF (G%CELL_NUMBER(IXW, IYW, IZW) /= IC) CYCLE
   IXG = GWC%IXG
   IYG = GWC%IYG
   IZG = GWC%IZG

   IF (IS_UNSTRUCTURED.AND.L%IS_SOLID(IXG, IYG, IZG)) RETURN
   SCARC_ASSIGN_SUBDIAG_TYPE = IW
   RETURN

ENDDO SEARCH_WALL_CELLS_LOOP

END FUNCTION SCARC_ASSIGN_SUBDIAG_TYPE


! ------------------------------------------------------------------------------------------------
!> \brief Set main diagonal entry for Poisson matrix in bandwise storage technique
! These values correspond to the full matrix of the global problem
! In case of an equidistant grid, we get the usual 5-point (2d) and 7-point (3d) stencil
! If two meshes with different step sizes meet, we get a weighted stencil along internal wall cells
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_MAINDIAGB (IC, IX, IY, IZ)
USE SCARC_POINTERS, ONLY: L, G, AB
INTEGER, INTENT(IN)  :: IC, IX, IY, IZ
INTEGER :: ID

AB => G%POISSONB
ID = AB%POS(0)               ! get column vector corresponding to matrix diagonal

AB%VAL(IC, ID) = AB%VAL(IC, ID) - 2.0_EB/(L%DXL(IX-1)*L%DXL(IX))
IF (.NOT.TWO_D)  AB%VAL(IC, ID) = AB%VAL(IC, ID) - 2.0_EB/(L%DYL(IY-1)*L%DYL(IY))
AB%VAL(IC, ID) = AB%VAL(IC, ID) - 2.0_EB/(L%DZL(IZ-1)*L%DZL(IZ))

AB%STENCIL(0) = AB%VAL(IC, ID)

END SUBROUTINE SCARC_SETUP_MAINDIAGB


! ------------------------------------------------------------------------------------------------
!> \brief Set subdigonal entries for Poisson matrix in bandwise storage technique on specified face
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_SUBDIAGB (IC, IX1, IY1, IZ1, IX2, IY2, IZ2, IOR0)
USE SCARC_POINTERS, ONLY: L, F, G, AB
INTEGER, INTENT(IN) :: IC, IX1, IY1, IZ1, IX2, IY2, IZ2, IOR0
INTEGER :: IW, ID
LOGICAL  :: IS_INTERNAL_CELL

F => L%FACE(IOR0)

! Decide wheter cell is interior or exterior cell
AB => G%POISSONB
ID = AB%POS(IOR0)                                
SELECT CASE (IOR0)
   CASE ( 1)
      IS_INTERNAL_CELL = IX1 > 1
   CASE (-1)
      IS_INTERNAL_CELL = IX1 < F%NOP
   CASE ( 2)
      IS_INTERNAL_CELL = IY1 > 1
   CASE (-2)
      IS_INTERNAL_CELL = IY1 < F%NOP
   CASE ( 3)
      IS_INTERNAL_CELL = IZ1 > 1
   CASE (-3)
      IS_INTERNAL_CELL = IZ1 < F%NOP
END SELECT

 
! If IC is an internal cell of the mesh, compute usual matrix contribution for corresponding subdiagonal
 
IF (IS_INTERNAL_CELL) THEN

   IF (IS_STRUCTURED .OR. .NOT.L%IS_SOLID(IX2, IY2, IZ2)) THEN
      AB%VAL(IC, ID)   = AB%VAL(IC, ID) + F%INCR_INSIDE
      AB%STENCIL(IOR0) = AB%VAL(IC, ID)
   ELSE
      CALL SCARC_SHUTDOWN(NSCARC_ERROR_MATRIX_SUBDIAG, SCARC_NONE, NSCARC_NONE)
   ENDIF

 
! If IC is a boundary cell of the mesh, compute matrix contribution only if there is a neighbor for that cell
 
!ELSE IF (L%FACE(IOR0)%N_NEIGHBORS /= 0) THEN
ELSE IF (L%FACE(IOR0)%N_NEIGHBORS == 123456) THEN       ! CAUTION: TO FIX AGAIN, ONLY FOR TESTING, IMPOSSIBLE CONDITION

   IW = SCARC_ASSIGN_SUBDIAG_TYPE (IC, IOR0)            ! get IW of a possibly suitable neighbor at face IOR0
   IF (IW /= 0) THEN
      AB%VAL(IC, ID)   = AB%VAL(IC, ID) + F%INCR_FACE
      AB%STENCIL(IOR0) = AB%VAL(IC, ID)
   ENDIF

ENDIF

END SUBROUTINE SCARC_SETUP_SUBDIAGB


! ------------------------------------------------------------------------------------------------
!> \brief Get maximum stencil size in specified matrix 
! This is known to be 7 for the 3D-Poisson matrix on finest level
! In algebraic multigrid-method this size results only in the course and can be much larger
! (required for dimensioning the coarse-level matrices)
! If NTYPE == 0, only internal matrix part is considered, if NTYPE == 1, also the overlap
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_GET_MATRIX_STENCIL_MAX (A, NLEN)
TYPE (SCARC_CMATRIX_TYPE), INTENT(INOUT) :: A
INTEGER, INTENT(IN) :: NLEN
INTEGER :: IC

A%N_STENCIL_MAX = 0
DO IC = 1, NLEN
   A%N_STENCIL_MAX = MAX(A%N_STENCIL_MAX, A%ROW(IC+1)-A%ROW(IC)+1)
ENDDO
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'GET_STENCIL_MAX:', A%N_STENCIL_MAX
#endif

END SUBROUTINE SCARC_GET_MATRIX_STENCIL_MAX


#ifdef WITH_MKL
! ------------------------------------------------------------------------------------------------
!> \brief Setup symmetric version of Poisson matrix for MKL solver in double precision
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_POISSON_MKL (NM, NL)
USE SCARC_POINTERS, ONLY: G, A, AS
USE SCARC_POINTER_ROUTINES, ONLY: SCARC_POINT_TO_GRID, SCARC_POINT_TO_CMATRIX
INTEGER, INTENT(IN) :: NM, NL
INTEGER :: IC, JC, JC0, ICS, JCS, JCG
INTEGER :: ICOL, JCOL, IAS
INTEGER :: ISYM, JSYM, NSYM
REAL(EB) :: VAL = 0.0_EB, VALS = 0.0_EB, DIFF
LOGICAL  :: BSYM, BCHECK_SYMMETRY = .FALSE.
INTEGER, DIMENSION(:), ALLOCATABLE :: ICOL_AUX, IC_AUX
INTEGER, POINTER, DIMENSION(:) :: ACOLG, ASCOLG

CROUTINE = 'SCARC_SETUP_POISSON_MKL'

CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G
A  => SCARC_POINT_TO_CMATRIX (G, NSCARC_MATRIX_POISSON)
AS => SCARC_POINT_TO_CMATRIX (G, NSCARC_MATRIX_POISSON_SYM)

IF (NMESHES == 1 .OR. TYPE_SCOPE(0) == NSCARC_SCOPE_LOCAL) THEN
   ACOLG  => A%COL
ELSE
   ACOLG  => A%COLG
ENDIF

#ifdef WITH_SCARC_DEBUG
CALL SCARC_DEBUG_CMATRIX (A, 'POISSON', 'SETUP_MATRIX_MKL: BEGIN')
WRITE(MSG%LU_DEBUG,*) 'TYPE_SCOPE(',0,')=', TYPE_SCOPE(0), NMESHES
WRITE(MSG%LU_DEBUG,*) 'TYPE_MKL(',NL,')=', TYPE_MKL(NL)
WRITE(MSG%LU_DEBUG,*) 'IS_MKL_LEVEL(',NL,') =', IS_MKL_LEVEL(NL)
WRITE(MSG%LU_DEBUG,*) 'ACOLG:', ACOLG
#endif
  
! ---------- Store only symmetric parts of matrix (diagonal and upper part)
  
IF (SCARC_MKL_MTYPE == 'SYMMETRIC') THEN

   IF (BCHECK_SYMMETRY) THEN
      ! First check whether symmetry of system matrix is guaranteed
      DO IC = 1, G%NC
         COLUMN_LOOP: DO ICOL = A%ROW(IC)+1, A%ROW(IC+1)-1
            ICS = ACOLG(ICOL)
            VAL = A%VAL(ICOL)
            IF (ICS > IC .AND. ICS <= G%NC) THEN
               BSYM = .FALSE.
               DO JCOL = A%ROW(ICS)+1, A%ROW(ICS+1)-1
                  JCS = ACOLG(JCOL)
                  IF (JCS == IC) THEN
                     VALS = A%VAL(JCOL)
                     DIFF = ABS(VAL-VALS)
                     IF (ABS(VAL - VALS) < 1E-6) THEN
                        BSYM=.TRUE.
                        CYCLE COLUMN_LOOP
                     ENDIF
                  ENDIF
               ENDDO
               IF (.NOT.BSYM) CALL SCARC_SHUTDOWN(NSCARC_ERROR_MATRIX_SYMMETRY, SCARC_NONE, NM)
            ENDIF
         ENDDO COLUMN_LOOP
      ENDDO
   ENDIF

 
   ! Compute number of entries in symmetric matrix
 
   AS%N_VAL = 0
   DO IC = 1, G%NC
      DO ICOL = A%ROW(IC), A%ROW(IC+1)-1
         IF (TYPE_MKL(NL) == NSCARC_MKL_LOCAL) THEN
            JC = ACOLG(ICOL)
            IF (JC >= IC .AND. JC <= G%NC) AS%N_VAL = AS%N_VAL+1
         ELSE IF (TYPE_MKL(NL) == NSCARC_MKL_GLOBAL) THEN
            IF (NL == NLEVEL_MIN) THEN
               JCG = G%LOCAL_TO_GLOBAL(ACOLG(ICOL))
            ELSE
               JCG = ACOLG(ICOL)
            ENDIF
            IF (JCG >= IC + G%NC_OFFSET(NM)) AS%N_VAL = AS%N_VAL+1
         ELSE
            CALL SCARC_SHUTDOWN(NSCARC_ERROR_MATRIX_SETUP, SCARC_NONE, TYPE_MKL(NL))
         ENDIF
      ENDDO
   ENDDO

ELSE
   AS%N_VAL = A%N_VAL
ENDIF

! Allocate storage for symmetric matrix and its column and row pointers
  
CALL SCARC_GET_MATRIX_STENCIL_MAX(A, G%NC)
AS%N_ROW = G%NC + 1
AS%N_VAL = A%N_STENCIL_MAX * G%NC
CALL SCARC_ALLOCATE_CMATRIX (AS, NL, TYPE_MKL_PRECISION, NSCARC_MATRIX_FULL, 'G%AS', CROUTINE)

IF (NMESHES == 1 .OR. TYPE_SCOPE(0) == NSCARC_SCOPE_LOCAL) THEN
   ASCOLG  => AS%COL
ELSE
   ASCOLG  => AS%COLG
ENDIF

! If global MKL method is used, also allocate auxiliary space for computation of global numbering

IF (IS_MKL_LEVEL(NL)) THEN
   CALL SCARC_ALLOCATE_INT1(ICOL_AUX, 1, A%N_STENCIL_MAX, NSCARC_HUGE_INT, 'ICOL_AUX', CROUTINE)
   CALL SCARC_ALLOCATE_INT1(IC_AUX  , 1, A%N_STENCIL_MAX, NSCARC_HUGE_INT, 'IC_AUX', CROUTINE)
ENDIF
  
! Subtract symmetric matrix part from usual system matrix
  
IAS = 1
DO IC = 1, AS%N_ROW - 1
   AS%ROW(IC) = IAS

   TYPE_MKL_SELECT: SELECT CASE (TYPE_MKL(NL)) 

      ! Blockwise use of local MKL solvers - no global numbering required

      CASE(NSCARC_MKL_LOCAL) 

         DO ICOL = A%ROW(IC), A%ROW(IC+1)-1
            JC = A%COL(ICOL)
            IF (JC >= IC .AND. JC <= G%NC) THEN
               AS%COL(IAS) = A%COL(ICOL)
               ASCOLG(IAS) = A%COL(ICOL)
               SELECT CASE (TYPE_MKL_PRECISION)
                  CASE (NSCARC_PRECISION_DOUBLE)
                     AS%VAL(IAS) = A%VAL(ICOL)
                  CASE (NSCARC_PRECISION_SINGLE)
                     AS%VAL_FB(IAS) = REAL(A%VAL(ICOL),FB)
                  END SELECT
               IAS = IAS + 1
            ENDIF
         ENDDO
         AS%ROW(IC+1) = IAS

      ! Global use of MKL solver - get global numbering of matrix elements

      CASE(NSCARC_MKL_GLOBAL) 

         ! Store indices of all diagonal and upper-diagonal entries

         ICOL_AUX = 0
         IC_AUX   = NSCARC_HUGE_INT
         ISYM = 1
         JC0 = ACOLG(A%ROW(IC))
         DO ICOL = A%ROW(IC), A%ROW(IC+1)-1
             JC = ACOLG(ICOL)
            IF (SCARC_MKL_MTYPE == 'SYMMETRIC') THEN
               IF (JC >= JC0) THEN
                  ICOL_AUX(ISYM) = ICOL
                  IC_AUX(ISYM) = JC
                  ISYM  = ISYM  + 1
               ENDIF
            ELSE
               ICOL_AUX(ISYM) = ICOL
               IC_AUX(ISYM) = JC
               ISYM  = ISYM  + 1
            ENDIF
         ENDDO
         AS%ROW(IC+1) = IAS

         NSYM = ISYM - 1
         JSYM = 1

         ! Sort them in increasing order (for the use of Cluster_Sparse_Solver and PARDISO functionality)

         SORT_LOOP: DO WHILE (JSYM <= NSYM)
            DO ISYM = 1, NSYM
               JC = IC_AUX(ISYM)
               IF (JC == NSCARC_HUGE_INT) CYCLE
               IF (JC <= MINVAL(ABS(IC_AUX(1:NSYM)))) THEN
                  ICOL = ICOL_AUX(ISYM)
                  SELECT CASE (TYPE_MKL_PRECISION)
                     CASE (NSCARC_PRECISION_DOUBLE)
                        AS%VAL(IAS) = A%VAL(ICOL)
                     CASE (NSCARC_PRECISION_SINGLE)
                        AS%VAL_FB(IAS) = REAL(A%VAL(ICOL), FB)
                  END SELECT
                  AS%COL(IAS) = ASCOLG(ICOL)
                  IC_AUX(ISYM) = NSCARC_HUGE_INT            ! mark entry as already used
                  IAS  = IAS  + 1
               ENDIF
            ENDDO
            JSYM = JSYM + 1
         ENDDO SORT_LOOP

   END SELECT TYPE_MKL_SELECT
ENDDO

AS%ROW(AS%N_ROW) = IAS

IF (IS_MKL_LEVEL(NL)) THEN
   CALL SCARC_DEALLOCATE_INT1 (ICOL_AUX, 'COL_AUX', CROUTINE)
   CALL SCARC_DEALLOCATE_INT1 (IC_AUX,  'IC_AUX', CROUTINE)
ENDIF

CALL SCARC_REDUCE_CMATRIX (AS, 'AS', CROUTINE)

#ifdef WITH_SCARC_DEBUG
CALL SCARC_DEBUG_CMATRIX(AS, 'AS', 'SETUP_MATRIX_MKL: END')
#endif
END SUBROUTINE SCARC_SETUP_POISSON_MKL
#endif


! ------------------------------------------------------------------------------------------------
!> \brief Insert correct boundary conditions into system matrix
!
! If A is a pure Neumann matrix, get neighboring cell indices of communicated stencil legs for 
! condensed system, also save values and column indices of last matrix row of last mesh
!
! Set correct boundary conditions for system matrix
! Take care of whether the structured or unstructured discretization is used
!
! If there are no Dirichlet BC's transform sytem into condensed one by replacing the
! matrix entries in last column and row by the stored ones (zeros and one at diaonal position)
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_BOUNDARY (NM, NL)
USE SCARC_POINTERS, ONLY: L, G, F, GWC, A, AB, ACO, ABCO
USE SCARC_POINTER_ROUTINES, ONLY: SCARC_POINT_TO_GRID
INTEGER, INTENT(IN) :: NM, NL
INTEGER :: I, J, K, IOR0, IW, IC, NOM, IP, ICO, ICOL

CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G

SELECT CASE (SCARC_GET_MATRIX_TYPE(NL))

   ! ---------- Matrix in compact storage technique
 
   CASE (NSCARC_MATRIX_COMPACT)

      A => G%POISSON

      ! Setup condensing if there are no Dirichlet BC's 

      IF (IS_PURE_NEUMANN) CALL SCARC_SETUP_CMATRIX_CONDENSED(NM)

      ! Set correct boundary conditions 

      DO IW = 1, G%NW

         GWC => G%WALL(IW)
         IOR0 = GWC%IOR
         IF (TWO_D .AND. ABS(IOR0) == 2) CYCLE       

         F  => L%FACE(IOR0)

         I = GWC%IXW
         J = GWC%IYW
         K = GWC%IZW

         IF (IS_UNSTRUCTURED .AND. L%IS_SOLID(I, J, K)) CYCLE

         NOM = GWC%NOM
         IC  = G%CELL_NUMBER(I, J, K)
         GWC%ICW = IC

         ! SPD-matrix with mixture of Dirichlet and Neumann BC's according to BTYPE

         IF (N_DIRIC_GLOBAL(NLEVEL_MIN) > 0) THEN

            IP = A%ROW(IC)
            SELECT CASE (GWC%BTYPE)
               CASE (DIRICHLET)
                  A%VAL(IP) = A%VAL(IP) - F%INCR_BOUNDARY
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,'(A,6I6,E14.6)') 'B :DIRICHLET: IW, I, J, K, NOM, IC, A%VAL:', IW, I, J, K, NOM, IC, A%VAL(IP)
#endif
               CASE (NEUMANN)
                  A%VAL(IP) = A%VAL(IP) + F%INCR_BOUNDARY
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,'(A,6I6,E14.6)') 'B :NEUMANN  : IW, I, J, K, NOM, IC, A%VAL:', IW, I, J, K, NOM, IC, A%VAL(IP)
#endif
            END SELECT

         ! Purely Neumann matrix

         ELSE IF (GWC%BTYPE == NEUMANN) THEN
            IP = A%ROW(IC)
            A%VAL(IP) = A%VAL(IP) + F%INCR_BOUNDARY
         ENDIF

      ENDDO 

      ! Transform into condensed system, if there are no Dirichlet BC's 

      IF (IS_PURE_NEUMANN) THEN
         DO ICO = 1, A%N_CONDENSED
            ACO => A%CONDENSED(ICO)
            DO ICOL = 1, ACO%N_COL
               IP = ACO%PTR(ICOL)
               A%VAL(IP) = ACO%VAL2(ICOL)
            ENDDO
         ENDDO
      ENDIF 

#ifdef WITH_SCARC_DEBUG
      CALL SCARC_DEBUG_CMATRIX(A, 'POISSON', 'POISSON WITH BDRY')
#endif

 
   ! ---------- Matrix in bandwise storage technique
 
   CASE (NSCARC_MATRIX_BANDWISE)

      ! Preset matrix switch if no Dirichlet BC's available

      AB => G%POISSONB
      IF (IS_PURE_NEUMANN) CALL SCARC_SETUP_BMATRIX_CONDENSED(NM)

      ! Set right boundary conditions 

      DO IW = 1, G%NW

         GWC => G%WALL(IW)
         IOR0 = GWC%IOR
         IF (TWO_D .AND. ABS(IOR0) == 2) CYCLE     

         F  => L%FACE(IOR0)

         I = GWC%IXW
         J = GWC%IYW
         K = GWC%IZW

         IF (IS_UNSTRUCTURED .AND. L%IS_SOLID(I, J, K)) CYCLE

         NOM  = GWC%NOM
         GWC%ICW =G%CELL_NUMBER(I, J, K)
         IC = G%CELL_NUMBER(I, J, K)

         ! SPD-matrix with mixture of Dirichlet and Neumann BC's according to the SETTING of BTYPE

         IF (N_DIRIC_GLOBAL(NLEVEL_MIN) > 0) THEN

            SELECT CASE (GWC%BTYPE)
               CASE (DIRICHLET)
                  AB%VAL(IC, AB%POS(0)) = AB%VAL(IC, AB%POS(0)) - F%INCR_BOUNDARY
               CASE (NEUMANN)
                  AB%VAL(IC, AB%POS(0)) = AB%VAL(IC, AB%POS(0)) + F%INCR_BOUNDARY
            END SELECT

         ! Purely Neumann matrix

         ELSE
            IF (GWC%BTYPE == NEUMANN) AB%VAL(IC, AB%POS(0)) = AB%VAL(IC, AB%POS(0)) + F%INCR_BOUNDARY
         ENDIF

      ENDDO 
   
      ! Transform into condensed system, if there are no Dirichlet BC's 

      IF (IS_PURE_NEUMANN) THEN
         DO ICO = 1, AB%N_CONDENSED
            ABCO => AB%CONDENSED(ICO)
            IF (ICO == 1) THEN
               AB%VAL(ABCO%ICO, 1:AB%N_STENCIL) = ABCO%VAL2(1:AB%N_STENCIL)
            ELSE
               IP = AB%POS(ABCO%IOR0)
               AB%VAL(ABCO%ICO, IP) = ABCO%VAL2(IP)
            ENDIF
         ENDDO
      ENDIF 
 
END SELECT 

END SUBROUTINE SCARC_SETUP_BOUNDARY


! ------------------------------------------------------------------------------------------------
!> \brief Setup condensed system for compact matrix storage technique
! Define switch entries for toggle between original and condensed values
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_CMATRIX_CONDENSED (NM)
USE SCARC_POINTERS, ONLY: L, G, A, ACO, GWC
INTEGER, INTENT(IN) :: NM
INTEGER :: ICO = 0, NC, NOM, IP, IC, JC, ICE, ICN, ICOL, IOR0, IW, I, J, K

A => G%POISSON
LAST_CELL_IN_LAST_MESH_IF: IF (NM == NMESHES) THEN

   NC = G%NC_LOCAL(NMESHES)
   IP = A%ROW(NC)

   ! Store column indices and values of diagonal and all off-diagonal entries in last row
   ! index '1' corresponds to main diagonal entry

   ICO = ICO + 1
   ACO => A%CONDENSED(ICO)

   ICOL = 1
   ACO%PTR(ICOL)  = IP
   ACO%COL(ICOL)  = A%COL(IP)
   ACO%VAL1(ICOL) = A%VAL(IP)
   ACO%VAL2(ICOL) = 1.0_EB

   DO IP = A%ROW(NC)+1, A%ROW(NC+1)-1
      ICOL = ICOL + 1
      ACO%PTR(ICOL)  = IP
      ACO%COL(ICOL)  = A%COL(IP)
      ACO%VAL1(ICOL) = A%VAL(IP)
      ACO%VAL2(ICOL) = 0.0_EB
   ENDDO
   ACO%N_COL = ICOL                                ! number of stored columns

 
   ! Within last mesh: check which other cells have a connection to the last cell;
   ! in each corresponding matrix row store the column index and value of just that matrix entry
   ! for each direction only one value has to be stored
 
   JC = NC - 1
   DO IP = A%ROW(JC)+1, A%ROW(JC+1)-1
      IF (A%COL(IP) == NC) THEN
         ICO = ICO + 1
         ACO => A%CONDENSED(ICO)
         ACO%PTR(1)  = IP
         ACO%COL(1)  = JC
         ACO%VAL1(1) = A%VAL(IP)                     ! store original value of system matrix
         ACO%VAL2(1) = 0.0_EB                        ! store new value of condensed system matrix
         ACO%N_COL   = 1
         EXIT
      ENDIF
   ENDDO

   JC = NC - L%NX
   DO IP = A%ROW(JC)+1, A%ROW(JC+1)-1
      IF (A%COL(IP) == NC) THEN
         ICO = ICO + 1
         ACO => A%CONDENSED(ICO)
         ACO%PTR(1)  = IP
         ACO%COL(1)  = JC
         ACO%VAL1(1) = A%VAL(IP)                     ! store original value of system matrix
         ACO%VAL2(1) = 0.0_EB                        ! store new value of condensed system matrix
         ACO%N_COL   = 1
         EXIT
      ENDIF
   ENDDO

   IF (.NOT.TWO_D) THEN
      JC = NC - L%NX * L%NY
      DO IP = A%ROW(JC)+1, A%ROW(JC+1)-1
         IF (A%COL(IP) == NC) THEN
            ICO = ICO + 1
            ACO => A%CONDENSED(ICO)
            ACO%PTR(1)  = IP
            ACO%COL(1)  = JC
            ACO%VAL1(1) = A%VAL(IP)                  ! store original value of system matrix
            ACO%VAL2(1) = 0.0_EB                     ! store new value of condensed system matrix
            ACO%N_COL   = 1
            EXIT
         ENDIF
      ENDDO
   ENDIF

ENDIF LAST_CELL_IN_LAST_MESH_IF

 
! Cycle boundary cells to check if there is a periodic communication partner whose stencil is coupled
! with the last cell of last mesh;
! this can be a cell on the opposite side of the own mesh or on a different mesh
! if such a cell exists, store corresponding matrix entry
 
DO IW = 1, G%NW

   GWC => G%WALL(IW)

   IOR0 = GWC%IOR
   IF (TWO_D .AND. ABS(IOR0) == 2) CYCLE

   I = GWC%IXW
   J = GWC%IYW
   K = GWC%IZW

   IF (IS_UNSTRUCTURED .AND. L%IS_SOLID(I, J, K)) CYCLE

   NOM = GWC%NOM
   IC  = G%CELL_NUMBER(I, J, K)
   GWC%ICW = IC

   IF (NOM == NMESHES) THEN

      ICE = GWC%ICE                               ! adjacent ghost cell number
      ICN = G%ICE_TO_ICN(ICE)                     ! get column index of neighboring offdiagonal matrix entry
      IF (ICN /= SCARC(NMESHES)%NC) CYCLE         ! if no relation to last cell in last mesh, cycle

      DO IP = A%ROW(IC)+1, A%ROW(IC+1)-1
         IF (A%COL(IP) == ICE) THEN
            ICO = ICO + 1
            ACO => A%CONDENSED(ICO)
            ACO%PTR(1)  = IP
            ACO%COL(1)  = ICN
            ACO%VAL1(1) = A%VAL(IP)
            ACO%VAL2(1) = 0.0_EB
            ACO%N_COL   = 1
            EXIT
         ENDIF
      ENDDO

   ENDIF 
ENDDO 

A%N_CONDENSED = ICO

END SUBROUTINE SCARC_SETUP_CMATRIX_CONDENSED


! ------------------------------------------------------------------------------------------------
!> \brief Setup condensed system for bandwise matrix storage technique
! Define switch entries for toggle between original and condensed values
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_BMATRIX_CONDENSED (NM)
USE SCARC_POINTERS, ONLY: L, G, AB, ABCO, GWC
INTEGER, INTENT(IN) :: NM
INTEGER :: ICO = 0, NC, NOM, IOR0, IC, JC, ICE, ICN, IW, I, J, K

AB => G%POISSONB
LAST_CELL_IN_LAST_MESH_BANDWISE_IF: IF (NM == NMESHES) THEN

   NC = G%NC_LOCAL(NMESHES)

   ! Store column indices and values of diagonal and all off-diagonal entries in last row
   ! index '1' corresponds to main diagonal entry
   ICO = ICO + 1
   ABCO => AB%CONDENSED(ICO)

   ABCO%IOR0 = 0
   ABCO%ICO  = NC
   ABCO%VAL1(1:AB%N_STENCIL) = AB%VAL(NC, 1:AB%N_STENCIL)
   ABCO%VAL2(1:AB%N_STENCIL) = 0.0_EB
   ABCO%VAL2(AB%POS(0)) = 1.0_EB

   ! Within last mesh: check which other cells have a connection to the last cell;
   ! in each corresponding matrix row store the column index and value of just that matrix entry
   ! for each direction only one value has to be stored
 
   JC = NC - 1
   DO IOR0 = -3, 3
      IF (JC + AB%OFFSET(IOR0) == NC) THEN
         ICO = ICO + 1
         ABCO => AB%CONDENSED(ICO)
         ABCO%IOR0 = IOR0
         ABCO%ICO  = JC
         ABCO%VAL1(1:AB%N_STENCIL) = AB%VAL(JC, 1:AB%N_STENCIL)
         ABCO%VAL2(1:AB%N_STENCIL) = AB%VAL(JC, 1:AB%N_STENCIL)
         ABCO%VAL2(AB%POS(ABCO%IOR0)) = 0.0_EB
         EXIT
      ENDIF
   ENDDO

   JC = NC - L%NX
   DO IOR0 = -3, 3
      IF (JC + AB%OFFSET(IOR0) == NC) THEN
         ICO = ICO + 1
         ABCO => AB%CONDENSED(ICO)
         ABCO%IOR0 = IOR0
         ABCO%ICO  = JC
         ABCO%VAL1(1:AB%N_STENCIL) = AB%VAL(JC, 1:AB%N_STENCIL)
         ABCO%VAL2(1:AB%N_STENCIL) = AB%VAL(JC, 1:AB%N_STENCIL)
         ABCO%VAL2(AB%POS(ABCO%IOR0)) = 0.0_EB
         EXIT
      ENDIF
   ENDDO

   IF (.NOT.TWO_D) THEN
      JC = NC - L%NX * L%NY
      DO IOR0 = -3, 3
         IF (JC + AB%OFFSET(IOR0) == NC) THEN
            ICO = ICO + 1
            ABCO => AB%CONDENSED(ICO)
            ABCO%IOR0 = IOR0
            ABCO%ICO  = JC
            ABCO%VAL1(1:AB%N_STENCIL) = AB%VAL(JC, 1:AB%N_STENCIL)
            ABCO%VAL2(1:AB%N_STENCIL) = AB%VAL(JC, 1:AB%N_STENCIL)
            ABCO%VAL2(AB%POS(ABCO%IOR0)) = 0.0_EB
            EXIT
         ENDIF
      ENDDO
   ENDIF

ENDIF LAST_CELL_IN_LAST_MESH_BANDWISE_IF

 
! Cycle boundary cells to check if there is a periodic communication partner whose stencil is coupled
! with the last cell of last mesh;
! this can be a cell on the opposite side of the own mesh or a cell on a different mesh
! if such a cell exists, store corresponding matrix entry
 
DO IW = 1, G%NW

   GWC => G%WALL(IW)

   IOR0 = GWC%IOR
   IF (TWO_D .AND. ABS(IOR0) == 2) CYCLE

   I    = GWC%IXW
   J    = GWC%IYW
   K    = GWC%IZW

   IF (IS_UNSTRUCTURED .AND. L%IS_SOLID(I, J, K)) CYCLE

   NOM = GWC%NOM
   IC  = G%CELL_NUMBER(I, J, K)
   GWC%ICW = IC

   IF (NOM == NMESHES) THEN
      ICE = GWC%ICE                               ! adjacent ghost cell number
      ICN = G%ICE_TO_ICN(ICE)                     ! get column index of neighboring offdiagonal matrix entry
      IF (ICN /= SCARC(NMESHES)%NC) CYCLE         ! if no relation to last cell in last mesh, cycle
      ICO = ICO + 1
      ABCO => AB%CONDENSED(ICO)
      ABCO%IOR0 = IOR0
      ABCO%ICO  = IC
      ABCO%VAL1(1:AB%N_STENCIL) = AB%VAL(IC, 1:AB%N_STENCIL)
      ABCO%VAL2(1:AB%N_STENCIL) = AB%VAL(IC, 1:AB%N_STENCIL)
      ABCO%VAL2(AB%POS(ABCO%IOR0)) = 0.0_EB
      EXIT
   ENDIF 
ENDDO 

AB%N_CONDENSED = ICO

END SUBROUTINE SCARC_SETUP_BMATRIX_CONDENSED


! ------------------------------------------------------------------------------------------------
!> \brief Setup condensed system in case of periodic or pure Neumann boundary conditions
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_SYSTEM_CONDENSED (NV, NL, ITYPE)
USE SCARC_POINTERS, ONLY: L, G, OG, F, OL, VC, A, ACO, AB, ABCO
USE SCARC_POINTER_ROUTINES, ONLY: SCARC_POINT_TO_GRID, SCARC_POINT_TO_OTHER_GRID, SCARC_POINT_TO_VECTOR
INTEGER, INTENT(IN) :: NV, NL, ITYPE
INTEGER :: NM, NOM, IFACE, ICN, ICE, ICW, JC, NC, ICO, IOR0, IP, ICG, INBR

IF (N_DIRIC_GLOBAL(NLEVEL_MIN) > 0 .OR. &
    TYPE_PRECON == NSCARC_RELAX_FFT .OR. TYPE_PRECON == NSCARC_RELAX_FFTO) RETURN

 
! In last mesh:  subtract B*RHS(end) for internal legs of stencil
 
MESH_REAL = 0.0_EB
IF (UPPER_MESH_INDEX == NMESHES) THEN

   CALL SCARC_POINT_TO_GRID (NMESHES, NL)

   NC =  G%NC_LOCAL(NMESHES)
   VC => SCARC_POINT_TO_VECTOR(NMESHES, NL, NV)

   ! Process last column entries of all rows except of last one
   ! for those rows only one matrix entry was stored, namely that one which connects to the last cell
 
   SELECT CASE (SCARC_GET_MATRIX_TYPE(NL))

      CASE (NSCARC_MATRIX_COMPACT)
         A => G%POISSON
         DO ICO = 2, A%N_CONDENSED
            ACO => A%CONDENSED(ICO)
            JC = ACO%COL(1)
            IF (JC < NC) VC(JC) = VC(JC) - ACO%VAL1(1)*VC(NC)
         ENDDO

      CASE (NSCARC_MATRIX_BANDWISE)
         AB => G%POISSONB
         DO ICO = 2, AB%N_CONDENSED
            ABCO => AB%CONDENSED(ICO)
            IP = AB%POS(ABCO%IOR0)
            JC = ABCO%ICO
            IF (JC < NC) VC(JC) = VC(JC) - ABCO%VAL1(IP)*VC(NC)
        ENDDO

   END SELECT

   MESH_REAL(NMESHES) = VC(NC)     ! store last entry of RHS
   VC(NC) = 0.0_EB                 ! set last entry of last mesh to zero

ENDIF

IF (ITYPE == 0) RETURN

 
! Broadcast last RHS-value of last cell in last mesh to all meshes
 
IF (N_MPI_PROCESSES > 1) &
   CALL MPI_ALLGATHER(MPI_IN_PLACE, 1, MPI_DOUBLE_PRECISION, MESH_REAL, 1, MPI_DOUBLE_PRECISION,&
                      MPI_COMM_WORLD, IERROR)

DO NM = 1, NMESHES
   SCARC(NM)%RHS_END = MESH_REAL(NMESHES)
ENDDO

 
! Only in case of periodic BC's:
! Subtract B*RHS(end) for corresponding entries of all periodic communication partners
 
DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G

   SNODE = PROCESS(NM)
   RNODE = PROCESS(NMESHES)

   IF (.NOT. ARE_NEIGHBORS(NM, NMESHES)) CYCLE

   CALL SCARC_POINT_TO_OTHER_GRID(NM, NMESHES, NL)
   VC => SCARC_POINT_TO_VECTOR (NM, NL, NV)

 
   ! Subtract B*RHS(end) at corresponding positions
 
   DO IFACE = 1, 6                                         ! check if this face has connection to last cell

      IOR0 = FACE_ORIENTATION(IFACE)
      F => L%FACE(IOR0)

      DO INBR = 1, F%N_NEIGHBORS

         NOM = F%NEIGHBORS(INBR)
         IF (NOM /= NMESHES) CYCLE                         ! only check for common matrix entries with last mesh
         CALL SCARC_POINT_TO_OTHER_GRID(NM, NOM, NL)

         DO ICG = OL%GHOST_FIRSTW(IOR0), OL%GHOST_LASTW(IOR0)

            ICW = OG%ICG_TO_ICW(ICG, 1)
            ICE = OG%ICG_TO_ICE(ICG, 1)
            ICN = G%ICE_TO_ICN(ICE)                        ! get column index of neighboring offdiagonal matrix entry

            IF (ICN /= SCARC(NMESHES)%NC) CYCLE            ! if no relation to last cell in last mesh, cycle

            VC(ICW) = VC(ICW) - F%INCR_FACE * SCARC(NM)%RHS_END

         ENDDO

      ENDDO
   ENDDO
ENDDO

END SUBROUTINE SCARC_SETUP_SYSTEM_CONDENSED


! --------------------------------------------------------------------------------------------------------
!> \brief Extract overlapping matrix parts after data exchange with neighbors and add them to main matrix
! --------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_EXTRACT_MATRIX_OVERLAPS (NMATRIX, NTYPE, NL)
USE SCARC_POINTERS, ONLY : G, F, OL, OG, A, OA
USE SCARC_POINTER_ROUTINES, ONLY: SCARC_POINT_TO_GRID, SCARC_POINT_TO_OTHER_GRID, &
                                  SCARC_POINT_TO_CMATRIX, SCARC_POINT_TO_OTHER_CMATRIX
INTEGER, INTENT(IN) :: NL, NMATRIX, NTYPE
INTEGER :: NM, IFACE, NOM, IOR0, ICG, ICE, IP, ICOL, INBR, ICN, ICE1, ICE2

MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID (NM, NL)                                 
   A => SCARC_POINT_TO_CMATRIX(G, NMATRIX)

   IP = A%ROW(G%NC+1)
   FACES_LOOP: DO IFACE = 1, 6               

      IOR0 = FACE_ORIENTATION(IFACE)
      F => SCARC(NM)%LEVEL(NLEVEL_MIN)%FACE(IOR0)
   
      DO INBR = 1, F%N_NEIGHBORS

         NOM = F%NEIGHBORS(INBR)
         CALL SCARC_POINT_TO_OTHER_GRID(NM, NOM, NL)
         OA => SCARC_POINT_TO_OTHER_CMATRIX(OG, NMATRIX)

         ICOL = 1
         DO ICG = OL%GHOST_FIRSTE(IOR0), OL%GHOST_LASTE(IOR0)
  
            ICE = OG%ICG_TO_ICE(ICG, 1)
            A%ROW(ICE) = IP 

            IF (NTYPE == 1) THEN
               ICOL = OA%ROW(ICG)
               ICN = ABS(OA%COLG(ICOL))
               A%COL(IP)  = ICE
               A%COLG(IP) = ICN
               A%VAL(IP) = OA%VAL(ICOL)
               IP = IP + 1
               DO ICOL = OA%ROW(ICG)+1, OA%ROW(ICG+1)-1
                  ICN = OA%COLG(ICOL)
                  IF (SCARC_CELL_WITHIN_MESH(G, NM, ICN)) THEN
                     A%COL(IP) = ABS(OA%COLG(ICOL)) - G%NC_OFFSET(NM)     
                  ELSE
                     A%COL(IP) = -ABS(OA%COLG(ICOL))
                     IF (ICG == OL%GHOST_FIRSTE(IOR0)) THEN
                        ICE2 = OG%ICG_TO_ICE(ICG+1, 1)
                        IF (G%LOCAL_TO_GLOBAL(ICE2) == ICN) A%COL(IP) = ICE2
                     ELSE IF (ICG == OL%GHOST_LASTW(IOR0)) THEN
                        ICE1 = OG%ICG_TO_ICE(ICG-1, 1)
                        IF (G%LOCAL_TO_GLOBAL(ICE1) == ICN) A%COL(IP) = ICE1
                     ELSE
                        ICE1 = OG%ICG_TO_ICE(ICG-1, 1)
                        ICE2 = OG%ICG_TO_ICE(ICG+1, 1)
                        IF (G%LOCAL_TO_GLOBAL(ICE1) == ICN) A%COL(IP) = ICE1
                        IF (G%LOCAL_TO_GLOBAL(ICE2) == ICN) A%COL(IP) = ICE2
                     ENDIF
                  ENDIF
                  A%COLG(IP) = ABS(OA%COLG(ICOL))      
                  A%VAL(IP)  = OA%VAL(ICOL)
                  IP = IP + 1
               ENDDO
            ELSE
               DO ICOL = OA%ROW(ICG), OA%ROW(ICG+1)-1
                  A%COL(IP) = -OA%COL(ICOL)   
                  A%COLG(IP) = ABS(OA%COLG(ICOL))      
                  A%VAL(IP) = OA%VAL(ICOL)
                  IP = IP + 1
               ENDDO
            ENDIF
         ENDDO

         A%ROW(ICE+1) = IP 
         A%N_ROW = ICE + 1
         A%N_VAL = IP - 1

      ENDDO
   ENDDO FACES_LOOP

#ifdef WITH_SCARC_DEBUG
CALL SCARC_DEBUG_CMATRIX (A, 'A', 'AFTER EXTRACT_MATRIX_OVERLAPS')
#endif

ENDDO MESHES_LOOP

END SUBROUTINE SCARC_EXTRACT_MATRIX_OVERLAPS


! ------------------------------------------------------------------------------------------------------
!> \brief Extract diagonal of Poisson matrix and store it in a separate vector for further use
! ------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_EXTRACT_MATRIX_DIAGONAL(NL)
USE SCARC_POINTERS, ONLY: G, A
USE SCARC_POINTER_ROUTINES, ONLY: SCARC_POINT_TO_GRID, SCARC_POINT_TO_CMATRIX
INTEGER, INTENT(IN) :: NL
INTEGER :: NM, IC, JC, ICOL

CROUTINE = 'SCARC_EXTRACT_MATRIX_DIAGONAL'

MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G
   A => SCARC_POINT_TO_CMATRIX (G, NSCARC_MATRIX_POISSON)

   CALL SCARC_ALLOCATE_REAL1 (G%DIAG, 1, G%NCE2, NSCARC_INIT_ZERO, 'G%DIAG', CROUTINE)
   DO IC = 1, G%NC
      DO ICOL = A%ROW(IC), A%ROW(IC+1) - 1
         JC = A%COL(ICOL)
         IF (JC == IC) G%DIAG(IC) = A%VAL(ICOL)
      ENDDO
   ENDDO

ENDDO MESHES_LOOP

! If there are multiple meshes exchange diagonal matrix on overlapping parts
IF (NMESHES > 1) CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_MATRIX_DIAGS, NSCARC_NONE, NL)

END SUBROUTINE SCARC_EXTRACT_MATRIX_DIAGONAL


! -------------------------------------------------------------------------------------------
!> \brief Extract overlapping zone information (including second layers)
! -------------------------------------------------------------------------------------------
SUBROUTINE SCARC_EXTRACT_ZONE_OVERLAPS(NL)
USE SCARC_POINTERS, ONLY : GC, GF, OLF, OGF, F
USE SCARC_POINTER_ROUTINES, ONLY: SCARC_POINT_TO_MULTIGRID, SCARC_POINT_TO_OTHER_MULTIGRID
INTEGER, INTENT(IN) :: NL
INTEGER :: NM, INBR, IOR0, NOM, IZ, ICG, ICE, ICE2, IFOUND, IZL_CURRENT
INTEGER :: IZL1, IZL2, IZG1, IZG2, IFACE

CROUTINE = 'SCARC_EXTRACT_ZONE_OVERLAPS'

MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_MULTIGRID (NM, NL, NL+1)    

   ! Clear overlapping parts and fill them with the recently exchanged data
   GF%ZONES_LOCAL (GF%NC+1: GF%NCE2) = 0
   GF%ZONES_GLOBAL(GF%NC+1: GF%NCE2) = 0

   IZL_CURRENT = GF%N_ZONES + 1

   DO IFACE = 1, 6                                        

      IOR0 = FACE_ORIENTATION(IFACE)
      F => SCARC(NM)%LEVEL(NLEVEL_MIN)%FACE(IOR0)

      DO INBR = 1, F%N_NEIGHBORS

         NOM = F%NEIGHBORS(INBR)
         CALL SCARC_POINT_TO_OTHER_MULTIGRID(NM, NOM, NL, NL+1)

         IZ = 0
         DO ICG = OLF%GHOST_FIRSTE(IOR0), OLF%GHOST_LASTE(IOR0)

            ICE  = OGF%ICG_TO_ICE(ICG, 1)

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,'(A,6I6)') 'ZONE_OVERLAPS: NM, INBR, NOM, IOR0, ICG, ICE:', NM, INBR, NOM, IOR0, ICG, ICE
#endif

            IZG1 = OGF%ICG_TO_GZONE(ICG)
            IFOUND = FINDLOC (GF%ZONES_GLOBAL, VALUE = IZG1, DIM = 1)
            IF (IFOUND == 0) THEN
               IZ = IZ + 1
               GF%N_ZONES = GF%N_ZONES + 1
               IZL1 = GF%N_ZONES
            ELSE
               IZL1 = GF%ZONES_LOCAL(IFOUND)
            ENDIF
            GF%ZONES_LOCAL(ICE)  = IZL1
            GF%ZONES_GLOBAL(ICE) = IZG1
            GC%LOCAL_TO_GLOBAL(IZL1) = IZG1

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,'(A,9I6)') 'A: NM, ICG, IFOUND, N_ZONES IZL1, ICE,  LOCAL(ICE),  GLOBAL(ICE) :', &
          NM, ICG, IFOUND, GF%N_ZONES, IZL1, ICE, GF%ZONES_LOCAL(ICE), GF%ZONES_GLOBAL(ICE), GC%LOCAL_TO_GLOBAL(IZL1)
#endif
 
            IF (NL /= NLEVEL_MIN) CYCLE

            ICE2 = OGF%ICG_TO_ICE(ICG, 2)
            IZG2 = OGF%ICG_TO_GZONE(ICG + OGF%NCG)
            IFOUND = FINDLOC (GF%ZONES_GLOBAL, VALUE = IZG2, DIM = 1)
            IF (IFOUND == 0) THEN
               GF%N_ZONES = GF%N_ZONES + 1
               IZL2 = GF%N_ZONES
            ELSE
               IZL2 = GF%ZONES_LOCAL(IFOUND)
            ENDIF
            GF%ZONES_LOCAL(ICE2)  = IZL2
            GF%ZONES_GLOBAL(ICE2) = IZG2
            GC%LOCAL_TO_GLOBAL(IZL2) = IZG2

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,'(A,9I6)') 'B: NM, ICG, IFOUND, N_ZONES IZL2, ICE2, LOCAL(ICE2), GLOBAL(ICE2):', &
          NM, ICG, IFOUND, GF%N_ZONES, IZL2, ICE2, GF%ZONES_LOCAL(ICE2), GF%ZONES_GLOBAL(ICE2), GC%LOCAL_TO_GLOBAL(IZL2)
#endif

         ENDDO
      ENDDO
   ENDDO 

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) '===================== EXTRACT_ZONE_OVERLAPS: NM=',NM
WRITE(MSG%LU_DEBUG,*) 'LOCAL_TO_GLOBAL'
WRITE(MSG%LU_DEBUG,'(8I6)') GC%LOCAL_TO_GLOBAL
CALL SCARC_DEBUG_ZONES(GF, -1, 1, 'AFTER EXTRACT_ZONES')
CALL SCARC_DEBUG_ZONES(GF, -1, 2, 'AFTER EXTRAXT_ZONES')
#endif

ENDDO MESHES_LOOP

END SUBROUTINE SCARC_EXTRACT_ZONE_OVERLAPS


! -------------------------------------------------------------------------------------------
!> \brief Setup pointers for overlapping zones for a pair of grid levels
! -------------------------------------------------------------------------------------------
SUBROUTINE SCARC_EXTRACT_ZONE_POINTERS(NL)
USE SCARC_POINTERS, ONLY : F, GF, GC, OLF, OLC, OGF, OGC
USE SCARC_POINTER_ROUTINES, ONLY: SCARC_POINT_TO_MULTIGRID, SCARC_POINT_TO_OTHER_MULTIGRID
INTEGER, INTENT(IN) :: NL
INTEGER :: INBR, IOR0, IZ, ICW1, ICW2, ICE1, ICE2, IZL1, IZL2, ICG, IZW, IZE, IFACE
INTEGER :: NM, NOM, NCGE_TOTAL = 0

CROUTINE = 'SCARC_EXTRACT_ZONE_POINTERS'

MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_MULTIGRID (NM, NL, NL+1)    
   
   DO IFACE = 1, 6                                        

      IOR0 = FACE_ORIENTATION(IFACE)
      F => SCARC(NM)%LEVEL(NLEVEL_MIN)%FACE(IOR0)

      DO INBR = 1, F%N_NEIGHBORS

         NOM = F%NEIGHBORS(INBR)
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) '============== ZONE_POINTERS: PROCESSING IFACE =', IFACE,' INBR, NOM =', INBR, NOM, NL
#endif
         CALL SCARC_POINT_TO_OTHER_MULTIGRID(NM, NOM, NL, NL+1)

         CALL SCARC_ALLOCATE_INT1(OGF%ICG_TO_IZONE, 1, 2*OGF%NCG, NSCARC_INIT_ZERO, 'OGF%ICG_TO_IZONE', CROUTINE)
         CALL SCARC_ALLOCATE_INT1(OGF%ICG_TO_EZONE, 1, 2*OGF%NCG, NSCARC_INIT_ZERO, 'OGF%ICG_TO_EZONE', CROUTINE)

         IZ  = 0
         INTERNAL_ZONES_LOOP: DO ICG = OLF%GHOST_FIRSTW(IOR0), OLF%GHOST_LASTW(IOR0)

            ICW1 = OGF%ICG_TO_ICW(ICG, 1)
            IZL1 = GF%ZONES_LOCAL(ICW1)
            IF (FINDLOC(OGF%ICG_TO_IZONE, VALUE = IZL1, DIM = 1) == 0) THEN
               IZ = IZ + 1
               OGF%ICG_TO_IZONE(IZ) = IZL1
            ENDIF

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,'(A, 8I6)') 'NM, INBR, NOM, IOR0, ICG, ICW1, IZL1, IZ :', NM, INBR, NOM, IOR0, ICG, ICW1, IZL1, IZ
#endif
            IF (NL /= NLEVEL_MIN) CYCLE

            ICW2 = OGF%ICG_TO_ICW(ICG, 2)
            IZL2 = GF%ZONES_LOCAL(ICW2)
            IF (FINDLOC(OGF%ICG_TO_IZONE, VALUE = IZL2, DIM = 1) == 0) THEN
               IZ = IZ + 1
               OGF%ICG_TO_IZONE(IZ) = IZL2
            ENDIF

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,'(A, 8I6)') 'NM, ,NBR, NOM, IOR0, ICG, ICW2, IZL2, IZ :', NM, INBR, NOM, IOR0, ICG, ICW2, IZL2, IZ
#endif
         ENDDO INTERNAL_ZONES_LOOP
         OGF%NCGI = IZ
         CALL SCARC_REDUCE_INT1(OGF%ICG_TO_IZONE, OGF%NCGI, 'OGF%ICG_TO_IZONE', CROUTINE)

         !First allocate in fine cell related length

         CALL SCARC_ALLOCATE_INT2(OGC%ICG_TO_ICW, 1, OGF%NCGI, 1, 1, NSCARC_INIT_ZERO, 'OGF%ICG_TO_ICW', CROUTINE)
         
         IZ = 0
         DO ICG = 1, OGF%NCGI
            IZW = OGF%ICG_TO_IZONE(ICG)
            IF (FINDLOC(OGC%ICG_TO_ICW(1:OGF%NCGI,1), VALUE = IZW, DIM = 1) == 0) THEN
               IZ = IZ + 1
               OGC%ICG_TO_ICW(IZ, 1) = IZW
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'NM=',NM,': OGC%ICG_TO_ICW(',IZ,',1)=',IZW
#endif
            ENDIF
         ENDDO
         OGC%NCG  = IZ
         OGC%NCGI = IZ

         OLC%GHOST_FIRSTW(IOR0) = 1
         OLC%GHOST_LASTW(IOR0)  = IZ

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'NM=',NM,': GHOST_FIRSTW(',IOR0,')=', OLC%GHOST_FIRSTW(IOR0)
WRITE(MSG%LU_DEBUG,*) 'NM=',NM,': GHOST_LASTW(',IOR0,') =', OLC%GHOST_LASTW(IOR0)
#endif
         ! Then reduce to zone related length 

         CALL SCARC_REDUCE_INT2(OGC%ICG_TO_ICW, OGC%NCGI, 1, 'OGC%ICG_TO_ICW', CROUTINE)


         IZ  = 0 
         EXTERNAL_ZONES_LOOP: DO ICG = OLF%GHOST_FIRSTE(IOR0), OLF%GHOST_LASTE(IOR0)

            ICE1 = OGF%ICG_TO_ICE(ICG, 1)
            IZL1 = GF%ZONES_LOCAL(ICE1)
            IF (FINDLOC(OGF%ICG_TO_EZONE, VALUE = IZL1, DIM = 1) == 0) THEN
               IZ = IZ + 1
               OGF%ICG_TO_EZONE(IZ) = IZL1
            ENDIF

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,'(A, 8I6)') 'NM, INBR, NOM, IOR0, ICG, ICE1, IZL1, IZ :', NM, INBR, NOM, IOR0, ICG, ICE1, IZL1, IZ
#endif
            IF (NL /= NLEVEL_MIN) CYCLE

            ICE2 = OGF%ICG_TO_ICE(ICG, 2)
            IZL2 = GF%ZONES_LOCAL(ICE2)
            IF (FINDLOC(OGF%ICG_TO_EZONE, VALUE = IZL2, DIM = 1) == 0) THEN
               IZ = IZ + 1
               OGF%ICG_TO_EZONE(IZ) = IZL2
            ENDIF

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,'(A, 8I6)') 'NM, INBR, NOM, IOR0, ICG, ICE2, IZL2, IZ :', NM, INBR, NOM, IOR0, ICG, ICE2, IZL2, IZ
#endif
         ENDDO EXTERNAL_ZONES_LOOP
         OGF%NCGE = IZ
         CALL SCARC_REDUCE_INT1(OGF%ICG_TO_EZONE, OGF%NCGE, 'OGF%ICG_TO_EZONE', CROUTINE)

         ! First allocate in fine cell related length

         CALL SCARC_ALLOCATE_INT2 (OGC%ICG_TO_ICE, 1, 2*OGF%NCGE, 1, 1, NSCARC_INIT_ZERO, 'OGC%ICG_TO_ICE', CROUTINE)
         
         IZ = 0
         DO ICG = 1, OGF%NCGE
            IZE = OGF%ICG_TO_EZONE(ICG)
            IF (FINDLOC(OGC%ICG_TO_ICE(1:OGF%NCGE,1), VALUE = IZE, DIM = 1) == 0) THEN
               IZ = IZ + 1
               OGC%ICG_TO_ICE(IZ, 1) = IZE
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'NM=',NM,': OGC%ICG_TO_ICE(',IZ,',1)=',IZE
#endif
            ENDIF
         ENDDO
         OGC%NCGE = IZ
         NCGE_TOTAL = NCGE_TOTAL + OGC%NCGE

         OLC%GHOST_FIRSTE(IOR0) = 1
         OLC%GHOST_LASTE(IOR0)  = IZ

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'NM=',NM,': GHOST_FIRSTE(',IOR0,')=', OLC%GHOST_FIRSTE(IOR0)
WRITE(MSG%LU_DEBUG,*) 'NM=',NM,': GHOST_LASTE(',IOR0,')= ', OLC%GHOST_LASTE(IOR0)
#endif
         ! Then reduce to zone related length 

         CALL SCARC_REDUCE_INT2(OGC%ICG_TO_ICE, OGC%NCGE, 1, 'OGC%ICG_TO_ICE', CROUTINE)

         GC%N_STENCIL_MAX = 25                  ! TODO: ONLY TEMPORARILY
         OGC%NLEN_BUFFER_LAYER1  = MAX(OGC%NCGI, OGC%NCGE)
         OGC%NLEN_BUFFER_LAYER2  = OGC%NLEN_BUFFER_LAYER1 * 2
         OGC%NLEN_BUFFER_LAYER4  = OGC%NLEN_BUFFER_LAYER1 * 4
         OGC%NLEN_BUFFER_STENCIL = OGC%NLEN_BUFFER_LAYER1 * GC%N_STENCIL_MAX
         OGC%NLEN_BUFFER_FULL    = OGC%NLEN_BUFFER_LAYER1 * GC%N_STENCIL_MAX * 2


#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) '===================== EXTRACT_ZONE_POINTERS: NM', NM
WRITE(MSG%LU_DEBUG,'(A,3I6)') 'ZONE_POINTERS INBR, NOM, IOR0:', INBR, NOM, IOR0
WRITE(MSG%LU_DEBUG,'(A,2I6)') 'EXCHANGE LENGTH WITH NEIGHBOR ',NOM, OGC%NLEN_BUFFER_LAYER1
WRITE(MSG%LU_DEBUG,*) '------------------------------------------------'
WRITE(MSG%LU_DEBUG,'(A,I6)') 'ICG_TO_IZONE(.): OGF%NCGI', OGF%NCGI
WRITE(MSG%LU_DEBUG,'(8I6)')  (OGF%ICG_TO_IZONE(IZ), IZ=1, OGF%NCGI)
WRITE(MSG%LU_DEBUG,*) '------------------------------------------------'
WRITE(MSG%LU_DEBUG,'(A,I6)') 'OGF:ICG_TO_EZONE(.):  OGF%NCGE', OGF%NCGE
WRITE(MSG%LU_DEBUG,'(8I6)') (OGF%ICG_TO_EZONE(IZ), IZ=1, OGF%NCGE)
WRITE(MSG%LU_DEBUG,*) '------------------------------------------------'
WRITE(MSG%LU_DEBUG,'(A,I6)') 'OGC%ICG_TO_ICW: OGC%NCGI', OGC%NCGI
WRITE(MSG%LU_DEBUG,'(8I6)')  OGC%ICG_TO_ICW(1:OGC%NCGI,1)
WRITE(MSG%LU_DEBUG,*) '------------------------------------------------'
WRITE(MSG%LU_DEBUG,'(A,I6)') 'OGC%ICG_TO_ICE: OGC%NCGE', OGC%NCGE
WRITE(MSG%LU_DEBUG,'(8I6)')  OGC%ICG_TO_ICE(1:OGC%NCGE,1)
WRITE(MSG%LU_DEBUG,*) '------------------------------------------------'
WRITE(MSG%LU_DEBUG,*) 'NCGE_TOTAL = ', NCGE_TOTAL
WRITE(MSG%LU_DEBUG,*) 'OGC%NLEN_BUFFER_LAYER1 = ', OGC%NLEN_BUFFER_LAYER1
WRITE(MSG%LU_DEBUG,*) 'OGC%NLEN_BUFFER_LAYER2 = ', OGC%NLEN_BUFFER_LAYER2
WRITE(MSG%LU_DEBUG,*) 'OGC%NLEN_BUFFER_LAYER4 = ', OGC%NLEN_BUFFER_LAYER4
WRITE(MSG%LU_DEBUG,*) 'OGC%NLEN_BUFFER_STENCIL= ', OGC%NLEN_BUFFER_STENCIL
WRITE(MSG%LU_DEBUG,*) 'OGC%NLEN_BUFFER_FULL   = ', OGC%NLEN_BUFFER_FULL
#endif

      ENDDO 
   ENDDO 
ENDDO MESHES_LOOP

END SUBROUTINE SCARC_EXTRACT_ZONE_POINTERS

! -------------------------------------------------------------------------------------------
!> \brief Identify cells on second layer
! -------------------------------------------------------------------------------------------
SUBROUTINE SCARC_IDENTIFY_LAYER2(NL)
USE SCARC_POINTERS, ONLY : S, A, G, OL, OG
USE SCARC_POINTER_ROUTINES, ONLY: SCARC_POINT_TO_GRID, SCARC_POINT_TO_OTHER_GRID, SCARC_POINT_TO_CMATRIX
INTEGER, INTENT(IN) :: NL
INTEGER :: NM, NOM, ICW, ICOL, JC, JCG, INBR, IOR0, ICG, IS

CROUTINE = 'SCARC_IDENTIFY_LAYER2'
 
MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID (NM, NL)    
   A => SCARC_POINT_TO_CMATRIX (G, NSCARC_MATRIX_POISSON)

   NEIGHBORS_LOOP: DO INBR = 1, SCARC(NM)%N_NEIGHBORS

      NOM = S%NEIGHBORS(INBR)
      CALL SCARC_POINT_TO_OTHER_GRID(NM, NOM, NL)
    
      CALL SCARC_ALLOCATE_INT1(OG%ICG_TO_ELAYER2, 1, 2*OG%NCG, NSCARC_INIT_ZERO, 'OG%ICG_TO_ELAYER2', CROUTINE)

      IS = 1
      FACE_LOOP: DO IOR0 = -3, 3
         IF (OL%GHOST_LASTW(IOR0) == 0) CYCLE
         GHOST_CELL_LOOP: DO ICG = OL%GHOST_FIRSTW(IOR0), OL%GHOST_LASTW(IOR0)
            ICW = OG%ICG_TO_ICW(ICG, 1)
            DO ICOL = A%ROW(ICW), A%ROW(ICW+1) - 1
               JC = A%COL(ICOL)
               IF (JC == 0) THEN
                  JCG = A%COLG(ICOL)
                  IF (FINDLOC (OG%ICG_TO_ELAYER2(1:2*OG%NCG), VALUE = JCG, DIM = 1) == 0) THEN 
                  OG%ICG_TO_ELAYER2(IS) = JCG
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,'(A,3I6)') '----------> JCG, IS, ICG_TO_ELAYER2 : ', JCG, IS, OG%ICG_TO_ELAYER2(IS)
#endif
                  IS = IS + 1
                  ENDIF
               ENDIF
            ENDDO
         ENDDO GHOST_CELL_LOOP
      ENDDO FACE_LOOP
      OL%N_LAYER2 = IS - 1

   ENDDO NEIGHBORS_LOOP
ENDDO MESHES_LOOP

END SUBROUTINE SCARC_IDENTIFY_LAYER2


! ----------------------------------------------------------------------------------------------------
!> \brief Setup all needed structures for a specified global ScaRC solver
! ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_METHODS
INTEGER :: NSTACK

SELECT_METHOD: SELECT CASE(TYPE_METHOD)

 
   ! ------------------ Global Krylov method -------------------------------------
 
   CASE (NSCARC_METHOD_KRYLOV)

      CALL SCARC_SETUP_KRYLOV_ENVIRONMENT()


    ! ------------------ Global Multigrid method -------------------------------------
     
    CASE (NSCARC_METHOD_MULTIGRID)

       CALL SCARC_SETUP_MULTIGRID_ENVIRONMENT()

   ! ---------------- McKenny-Greengard-Mayo method (MGM) --------------------
 
   CASE (NSCARC_METHOD_MGM)

      ! Allocate velocity vectors along internal obstructions for the setting of internal BC's

      CALL SCARC_SETUP_MGM(NLEVEL_MIN, NLEVEL_MIN)

      ! ------- First part of method: Setup CG solver for inhomogeneous problem on structured discretization
      !         Use FFT-preconditioning by default

      WRITE(*,*) 'CAUTION: TODO: PRECON SCARC MGM!'
      TYPE_PRECON = NSCARC_RELAX_FFT
      TYPE_PRECON = NSCARC_RELAX_SSOR
      CALL SCARC_SETUP_GRID_TYPE(NSCARC_GRID_STRUCTURED)

      NSTACK = NSCARC_STACK_ROOT
      STACK(NSTACK)%SOLVER => MAIN_CG_STRUCTURED
      CALL SCARC_SETUP_KRYLOV(NSCARC_SOLVER_MAIN, NSCARC_SCOPE_GLOBAL, NSCARC_STAGE_ONE, NSTACK, NLEVEL_MIN, NLEVEL_MIN)

      NSTACK = NSTACK + 1
      STACK(NSTACK)%SOLVER => PRECON_FFT
      CALL SCARC_SETUP_PRECON(NSTACK, NSCARC_SCOPE_LOCAL)
      CALL SCARC_SETUP_FFT(NLEVEL_MIN, NLEVEL_MIN)


      ! ------- Second part of method: Setup CG solver for homogeneous problem on unstructured discretization
      !         Only working for compact matrix storage technique (because of the unstructured grid)
      !         Use LU-preconditioning by default

      WRITE(*,*) 'CAUTION: TODO: PRECON USCARC MGM!'
      TYPE_PRECON = NSCARC_RELAX_SSOR
      TYPE_MATRIX = NSCARC_MATRIX_COMPACT
      CALL SCARC_SETUP_GRID_TYPE(NSCARC_GRID_UNSTRUCTURED)

      NSTACK = NSTACK + 1
      STACK(NSTACK)%SOLVER => MAIN_CG_UNSTRUCTURED
      CALL SCARC_SETUP_KRYLOV(NSCARC_SOLVER_MAIN, NSCARC_SCOPE_GLOBAL, NSCARC_STAGE_ONE, NSTACK, NLEVEL_MIN, NLEVEL_MIN)

      NSTACK = NSTACK + 1
      IF (TYPE_PRECON == NSCARC_RELAX_MKL .AND. TYPE_MATRIX == NSCARC_MATRIX_COMPACT) THEN
#ifdef WITH_MKL
         STACK(NSTACK)%SOLVER => PRECON_MKL
         CALL SCARC_SETUP_PRECON(NSTACK, NSCARC_SCOPE_LOCAL)
         CALL SCARC_SETUP_PARDISO(NLEVEL_MIN, NLEVEL_MIN)          ! use global PARDISO from MKL
#else
         WRITE(*,*) 'MGM-method: MKL-preconditioning required, MKL library not available, using LU preconditioning'
         STACK(NSTACK)%SOLVER => PRECON_SSOR
         CALL SCARC_SETUP_PRECON(NSTACK, NSCARC_SCOPE_LOCAL)
         !CALL SCARC_SETUP_LU(NLEVEL_MIN, NLEVEL_MIN)
#endif

     ELSE
         STACK(NSTACK)%SOLVER => PRECON_SSOR
         CALL SCARC_SETUP_PRECON(NSTACK, NSCARC_SCOPE_LOCAL)
         !CALL SCARC_SETUP_LU(NLEVEL_MIN, NLEVEL_MAX)
     ENDIF



   ! ------------------ MKL method -------------------------------------
 
#ifdef WITH_MKL
   CASE (NSCARC_METHOD_LU)

      NSTACK = NSCARC_STACK_ROOT
      STACK(NSTACK)%SOLVER => MAIN_LU

      CALL SCARC_SETUP_MKL(NSCARC_SOLVER_MAIN, NSCARC_SCOPE_GLOBAL, NSCARC_STAGE_ONE, NSTACK, NLEVEL_MIN, NLEVEL_MIN)

      ! In the multi-mesh case use CLUSTER_SPARSE_SOLVER, else PARDISO solver (only on finest grid level)

      IF (NMESHES > 1) THEN
         CALL SCARC_SETUP_CLUSTER(NLEVEL_MIN, NLEVEL_MIN)
      ELSE
         CALL SCARC_SETUP_PARDISO(NLEVEL_MIN, NLEVEL_MIN)
      ENDIF
#endif

END SELECT SELECT_METHOD

! Store total number of stack entries (used solvers)

N_STACK_TOTAL = NSTACK

END SUBROUTINE SCARC_SETUP_METHODS


! ----------------------------------------------------------------------------------------------------
!> \brief Allocate and initialize vectors for Krylov method
! ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_VECTORS()
USE SCARC_POINTERS, ONLY: G, SV, ST
USE SCARC_POINTER_ROUTINES, ONLY: SCARC_POINT_TO_GRID
INTEGER :: NM, NSTACK, NL

CROUTINE = 'SCARC_SETUP_VECTORS'

! If multiple grid types are used (currently only in MGM method) take care that the vectors are
! allocated in the longest necessary length which corresponds to the structured discretization.
! The related workspaces are also used for possible shorter instances in the unstructured discretization
IF (HAS_MULTIPLE_GRIDS) CALL SCARC_SETUP_GRID_TYPE(NSCARC_GRID_STRUCTURED)   

DO NSTACK = 1, N_STACK_TOTAL

   SV  => STACK(NSTACK)%SOLVER

   DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
      DO NL = SV%TYPE_LEVEL(1), SV%TYPE_LEVEL(2)

         CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G
         ST => SCARC(NM)%LEVEL(NL)%STAGE(SV%TYPE_STAGE)

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'SETUP_VECTORS: NC, NCE =', G%NC, G%NCE, SV%X, SV%B, SV%D
#endif
         IF (SV%X /= NSCARC_UNDEF_INT) CALL SCARC_ALLOCATE_REAL1(ST%X, 1, G%NCE, NSCARC_INIT_ZERO, 'ST%X', CROUTINE)
         IF (SV%B /= NSCARC_UNDEF_INT) CALL SCARC_ALLOCATE_REAL1(ST%B, 1, G%NCE, NSCARC_INIT_ZERO, 'ST%B', CROUTINE)
         IF (SV%D /= NSCARC_UNDEF_INT) CALL SCARC_ALLOCATE_REAL1(ST%D, 1, G%NCE, NSCARC_INIT_ZERO, 'ST%D', CROUTINE)
         IF (SV%R /= NSCARC_UNDEF_INT) CALL SCARC_ALLOCATE_REAL1(ST%R, 1, G%NCE, NSCARC_INIT_ZERO, 'ST%R', CROUTINE)
         IF (SV%V /= NSCARC_UNDEF_INT) CALL SCARC_ALLOCATE_REAL1(ST%V, 1, G%NCE, NSCARC_INIT_ZERO, 'ST%V', CROUTINE)
         IF (SV%Y /= NSCARC_UNDEF_INT) CALL SCARC_ALLOCATE_REAL1(ST%Y, 1, G%NCE, NSCARC_INIT_ZERO, 'ST%Y', CROUTINE)
         IF (SV%Z /= NSCARC_UNDEF_INT) CALL SCARC_ALLOCATE_REAL1(ST%Z, 1, G%NCE, NSCARC_INIT_ZERO, 'ST%Z', CROUTINE)

#ifdef WITH_SCARC_DEBUG
         IF (SV%E /= NSCARC_UNDEF_INT) CALL SCARC_ALLOCATE_REAL1(ST%E, 1, G%NCE, NSCARC_INIT_ZERO, 'ST%E', CROUTINE)
#endif

         IF (TYPE_MKL_PRECISION == NSCARC_PRECISION_SINGLE) THEN
            IF (SV%X_FB /= NSCARC_UNDEF_INT) CALL SCARC_ALLOCATE_REAL1_FB(ST%X_FB, 1, G%NCE, NSCARC_INIT_ZERO, 'ST%X_FB', CROUTINE)
            IF (SV%B_FB /= NSCARC_UNDEF_INT) CALL SCARC_ALLOCATE_REAL1_FB(ST%B_FB, 1, G%NCE, NSCARC_INIT_ZERO, 'ST%B_FB', CROUTINE)
            IF (SV%R_FB /= NSCARC_UNDEF_INT) CALL SCARC_ALLOCATE_REAL1_FB(ST%R_FB, 1, G%NCE, NSCARC_INIT_ZERO, 'ST%R_FB', CROUTINE)
            IF (SV%V_FB /= NSCARC_UNDEF_INT) CALL SCARC_ALLOCATE_REAL1_FB(ST%V_FB, 1, G%NCE, NSCARC_INIT_ZERO, 'ST%V_FB', CROUTINE)
         ENDIF

      ENDDO
   ENDDO
ENDDO

END SUBROUTINE SCARC_SETUP_VECTORS




#ifdef WITH_MKL
! ----------------------------------------------------------------------------------------------------
!> \brief Allocate and initialize vectors for LU-solvers (based on MKL)
! ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_MKL(NSOLVER, NSCOPE, NSTAGE, NSTACK, NLMIN, NLMAX)
USE SCARC_POINTERS, ONLY: SV
INTEGER, INTENT(IN) :: NSOLVER, NSCOPE, NSTAGE, NSTACK, NLMIN, NLMAX

 
! Basic setup of stack information and types for MKL
 
CALL SCARC_SETUP_STACK(NSTACK)

SV  => STACK(NSTACK)%SOLVER
SELECT CASE (NSOLVER)
   CASE (NSCARC_SOLVER_MAIN)
      SV%CNAME = 'SCARC_MAIN_LU'
   CASE (NSCARC_SOLVER_COARSE)
      SV%CNAME = 'SCARC_COARSE_LU'
   CASE DEFAULT
      CALL SCARC_SHUTDOWN(NSCARC_ERROR_PARSE_INPUT, SCARC_NONE, NSOLVER)
END SELECT

! Reset types for LU-decomposition method
 
SV%TYPE_METHOD        = NSCARC_METHOD_LU
SV%TYPE_SOLVER        = NSOLVER
SV%TYPE_SCOPE(0)      = NSCOPE
SV%TYPE_STAGE         = NSTAGE
SV%TYPE_GRID          = TYPE_GRID
SV%TYPE_LEVEL(0)      = NLMAX-NLMIN+1
SV%TYPE_LEVEL(1)      = NLMIN
SV%TYPE_LEVEL(2)      = NLMAX
SV%TYPE_MATRIX        = NSCARC_MATRIX_COMPACT    
SV%TYPE_MKL_PRECISION = TYPE_MKL_PRECISION

 
! Point to solution vectors (in corresponding scope)
 
CALL SCARC_SETUP_REFERENCES(.TRUE.,.TRUE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.TRUE., NSTACK)

END SUBROUTINE SCARC_SETUP_MKL
#endif


! ----------------------------------------------------------------------------------------------------
!> \brief Allocate and initialize vectors for additive or multiplicative coarse grid
! (corresponding to Schwarz domain decomposition method)
! ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_INTERPOLATION(NSTAGE, NLMIN, NLMAX)
USE SCARC_POINTERS, ONLY: G, ST
USE SCARC_POINTER_ROUTINES, ONLY: SCARC_POINT_TO_GRID
INTEGER, INTENT(IN) :: NSTAGE, NLMIN, NLMAX
INTEGER :: NM, NL

CROUTINE = 'SCARC_SETUP_INTERPOLATION'

MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   LEVEL_LOOP: DO NL = NLMIN, NLMAX

      CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G
      ST => SCARC(NM)%LEVEL(NL)%STAGE(NSTAGE)

      CALL SCARC_ALLOCATE_REAL1(ST%X, 1, G%NCE, NSCARC_INIT_ZERO, 'ST%X', CROUTINE)
      CALL SCARC_ALLOCATE_REAL1(ST%B, 1, G%NCE, NSCARC_INIT_ZERO, 'ST%B', CROUTINE)
      CALL SCARC_ALLOCATE_REAL1(ST%V, 1, G%NCE, NSCARC_INIT_ZERO, 'ST%Q', CROUTINE)
      CALL SCARC_ALLOCATE_REAL1(ST%R, 1, G%NCE, NSCARC_INIT_ZERO, 'ST%W', CROUTINE)
      CALL SCARC_ALLOCATE_REAL1(ST%Y, 1, G%NCE, NSCARC_INIT_ZERO, 'ST%Y', CROUTINE)
      CALL SCARC_ALLOCATE_REAL1(ST%Z, 1, G%NCE, NSCARC_INIT_ZERO, 'ST%Z', CROUTINE)

   ENDDO LEVEL_LOOP
ENDDO MESHES_LOOP

END SUBROUTINE SCARC_SETUP_INTERPOLATION


! ----------------------------------------------------------------------------------------------------
!> \brief Setup data structures for the use of blockwise FFT methods as preconditioners
! New here: Perform own initialization of FFT based on H2CZIS/H3CZIS and use own SAVE and WORK arrays
! ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_FFT(NLMIN, NLMAX)
USE SCARC_POINTERS, ONLY: M, S, L, FFT
USE SCARC_POINTER_ROUTINES, ONLY: SCARC_POINT_TO_GRID
USE POIS, ONLY: H2CZIS, H3CZIS
INTEGER, INTENT(IN) :: NLMIN, NLMAX
INTEGER :: NM, NL

CROUTINE = 'SCARC_SETUP_FFT'
 
! Allocate working space for FFT routine
 
MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   LEVEL_LOOP: DO NL = NLMIN, NLMAX

      CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G
      FFT => SCARC(NM)%LEVEL(NL)%FFT

      FFT%LBC = M%LBC
      FFT%MBC = M%MBC
      FFT%NBC = M%NBC

      FFT%XS = S%XS
      FFT%XF = S%XF
      FFT%YS = S%YS
      FFT%YF = S%YF
      FFT%ZS = S%ZS
      FFT%ZF = S%ZF

      FFT%IBAR = L%NX
      FFT%JBAR = L%NY
      FFT%KBAR = L%NZ

      FFT%ITRN = L%NX+1
      IF (TWO_D) THEN
         FFT%JTRN = 1
      ELSE
         FFT%JTRN = L%NY+1
      ENDIF
      FFT%KTRN = L%NZ+1

      FFT%LSAVE = (FFT%ITRN+1)*FFT%JTRN*FFT%KTRN+7*FFT%ITRN+5*FFT%JTRN+6*FFT%KTRN+56
      FFT%LWORK = (FFT%ITRN+1)*FFT%JTRN*FFT%KTRN

      CALL SCARC_ALLOCATE_REAL1(FFT%SAVE1, -3, FFT%LSAVE, NSCARC_INIT_ZERO, 'FFT%SAVE1', CROUTINE)
      CALL SCARC_ALLOCATE_REAL1(FFT%WORK ,  1, FFT%LWORK, NSCARC_INIT_ZERO, 'FFT%WORK', CROUTINE)

      ! Allocate stretching vector (set to 1)
      CALL SCARC_ALLOCATE_REAL1(FFT%HX, 0, FFT%ITRN, NSCARC_INIT_ONE, 'FFT%HX', CROUTINE)

      ! Allocate RHS vector for FFT routine
      IF (L%NY == 1) THEN
         CALL SCARC_ALLOCATE_REAL3(FFT%PRHS, 1, FFT%ITRN, 1,        1, 1, FFT%KTRN, NSCARC_INIT_ZERO, 'FFT%PRHS', CROUTINE)
      ELSE
         CALL SCARC_ALLOCATE_REAL3(FFT%PRHS, 1, FFT%ITRN, 1, FFT%JTRN, 1, FFT%KTRN, NSCARC_INIT_ZERO, 'FFT%PRHS', CROUTINE)
      ENDIF

      ! Allocate boundary data vector for XS
      IF (L%NZ>1) THEN
         IF (L%NY  > 1) CALL SCARC_ALLOCATE_REAL2(FFT%BXS, 1, FFT%JTRN, 1, FFT%KTRN, NSCARC_INIT_ZERO, 'FFT%BXS', CROUTINE)
         IF (L%NY == 1) CALL SCARC_ALLOCATE_REAL2(FFT%BXS, 1,        1, 1, FFT%KTRN, NSCARC_INIT_ZERO, 'FFT%BXS', CROUTINE)
      ELSE
         CALL SCARC_ALLOCATE_REAL2(FFT%BXS, 1, FFT%JTRN, 1, 1, NSCARC_INIT_ZERO, 'FFT%BXS', CROUTINE)
      ENDIF

      ! Allocate boundary data vector for XF
      IF (L%NZ>1) THEN
         IF (L%NY  > 1) CALL SCARC_ALLOCATE_REAL2(FFT%BXF, 1, FFT%JTRN, 1, FFT%KTRN, NSCARC_INIT_ZERO, 'FFT%BXF', CROUTINE)
         IF (L%NY == 1) CALL SCARC_ALLOCATE_REAL2(FFT%BXF, 1,        1, 1, FFT%KTRN, NSCARC_INIT_ZERO, 'FFT%BXF', CROUTINE)
      ELSE
         CALL SCARC_ALLOCATE_REAL2(FFT%BXF, 1, FFT%JTRN, 1, 1, NSCARC_INIT_ZERO, 'FFT%BXF', CROUTINE)
      ENDIF

      ! Allocate boundary data vector for YS
      IF (L%NZ > 1) THEN
         CALL SCARC_ALLOCATE_REAL2(FFT%BYS, 1, FFT%ITRN,1, FFT%KTRN, NSCARC_INIT_ZERO, 'FFT%BYS', CROUTINE)
      ELSE
         CALL SCARC_ALLOCATE_REAL2(FFT%BYS, 1, FFT%ITRN,1,        1, NSCARC_INIT_ZERO, 'FFT%BYS', CROUTINE)
      ENDIF

      ! Allocate boundary data vector for YF
      IF (L%NZ > 1) THEN
         CALL SCARC_ALLOCATE_REAL2(FFT%BYF, 1, FFT%ITRN,1, FFT%KTRN, NSCARC_INIT_ZERO, 'FFT%BYF', CROUTINE)
      ELSE
         CALL SCARC_ALLOCATE_REAL2(FFT%BYF, 1, FFT%ITRN,1,        1, NSCARC_INIT_ZERO, 'FFT%BYF', CROUTINE)
      ENDIF

      ! Allocate boundary data vector for ZS
      IF (L%NY  > 1) CALL SCARC_ALLOCATE_REAL2(FFT%BZS, 1, FFT%ITRN, 1, FFT%JTRN, NSCARC_INIT_ZERO, 'FFT%BZS', CROUTINE)
      IF (L%NY == 1) CALL SCARC_ALLOCATE_REAL2(FFT%BZS, 1, FFT%ITRN, 1,        1, NSCARC_INIT_ZERO, 'FFT%BZS', CROUTINE)

      ! Allocate boundary data vector for ZF
      IF (L%NY  >1)  CALL SCARC_ALLOCATE_REAL2(FFT%BZF, 1, FFT%ITRN, 1, FFT%JTRN, NSCARC_INIT_ZERO, 'FFT%BZF', CROUTINE)
      IF (L%NY == 1) CALL SCARC_ALLOCATE_REAL2(FFT%BZF, 1, FFT%ITRN, 1,        1, NSCARC_INIT_ZERO, 'FFT%BZF', CROUTINE)

      IF (TWO_D) THEN
         CALL H2CZIS(FFT%XS,FFT%XF,FFT%IBAR,FFT%LBC,&
                     FFT%ZS,FFT%ZF,FFT%KBAR,FFT%NBC,&
                     FFT%HX,FFT%XLM,FFT%ITRN,IERROR,FFT%SAVE1)
      ELSE
         CALL H3CZIS(FFT%XS,FFT%XF,FFT%IBAR,FFT%LBC,&
                     FFT%YS,FFT%YF,FFT%JBAR,FFT%MBC,&
                     FFT%ZS,FFT%ZF,FFT%KBAR,FFT%NBC,&
                     FFT%HX,FFT%XLM,FFT%ITRN,FFT%JTRN,IERROR,FFT%SAVE1)
      ENDIF

   ENDDO LEVEL_LOOP
ENDDO MESHES_LOOP

END SUBROUTINE SCARC_SETUP_FFT


! ----------------------------------------------------------------------------------------------------
!> \brief Allocate and initialize vectors for blockwise FFT methods with overlap
! New here: Perform own initialization of FFT based on H2CZIS/H3CZIS and use own SAVE and WORK arrays
! ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_FFTO(NLMIN, NLMAX)
USE SCARC_POINTERS, ONLY: M, S, L, FFT
USE SCARC_POINTER_ROUTINES, ONLY: SCARC_POINT_TO_GRID
USE POIS, ONLY: H2CZIS, H3CZIS
INTEGER, INTENT(IN) :: NLMIN, NLMAX
INTEGER :: NM, NL

CROUTINE = 'SCARC_SEtUP_FFTO'

! Allocate working space for FFT routine
 
MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   LEVEL_LOOP: DO NL = NLMIN, NLMAX

      CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G
      FFT => SCARC(NM)%LEVEL(NL)%FFT

      FFT%LBC = M%LBC
      FFT%MBC = M%MBC
      FFT%NBC = M%NBC

      IF (NM == 1) THEN
         FFT%IBAR = S%IBAR+1
         FFT%XS = S%XS
         FFT%XF = S%XF + L%DX
      ELSE IF (NM == NMESHES) THEN
         FFT%IBAR = S%IBAR+1
         FFT%XS = S%XS - L%DX
         FFT%XF = S%XF
      ELSE 
         FFT%IBAR = S%IBAR+2
         FFT%XS = S%XS - L%DX
         FFT%XF = S%XF + L%DX
      ENDIF

      FFT%JBAR = S%JBAR
      FFT%KBAR = S%KBAR

      FFT%YS = S%YS
      FFT%YF = S%YF
      FFT%ZS = S%ZS
      FFT%ZF = S%ZF

      FFT%ITRN = FFT%IBAR+1
      IF (TWO_D) THEN
         FFT%JTRN = 1
      ELSE
         FFT%JTRN = L%NY+1
      ENDIF
      FFT%KTRN = L%NZ+1

      FFT%LSAVE = (FFT%ITRN+1)*FFT%JTRN*FFT%KTRN+7*FFT%ITRN+5*FFT%JTRN+6*FFT%KTRN+56
      FFT%LWORK = (FFT%ITRN+1)*FFT%JTRN*FFT%KTRN

      CALL SCARC_ALLOCATE_REAL1(FFT%SAVE1, -3, FFT%LSAVE, NSCARC_INIT_ZERO, 'FFT%SAVE1', CROUTINE)
      CALL SCARC_ALLOCATE_REAL1(FFT%WORK ,  1, FFT%LWORK, NSCARC_INIT_ZERO, 'FFT%WORK', CROUTINE)

      ! Allocate stretching vector (set to 1)

      CALL SCARC_ALLOCATE_REAL1(FFT%HX, 0, FFT%ITRN, NSCARC_INIT_ONE, 'FFT', CROUTINE)

      ! Allocate RHS vector for FFT routine

      IF (L%NY == 1) THEN
         CALL SCARC_ALLOCATE_REAL3(FFT%PRHS, 1, FFT%ITRN, 1,        1, 1, FFT%KTRN, NSCARC_INIT_ZERO, 'FFT%PRHS', CROUTINE)
      ELSE
         CALL SCARC_ALLOCATE_REAL3(FFT%PRHS, 1, FFT%ITRN, 1, FFT%JTRN, 1, FFT%KTRN, NSCARC_INIT_ZERO, 'FFT%PRHS', CROUTINE)
      ENDIF

      ! Allocate boundary data vector for XS

      IF (L%NZ>1) THEN
         IF (L%NY  > 1) CALL SCARC_ALLOCATE_REAL2(FFT%BXS, 1, FFT%JTRN, 1, FFT%KTRN, NSCARC_INIT_ZERO, 'FFT%BXS', CROUTINE)
         IF (L%NY == 1) CALL SCARC_ALLOCATE_REAL2(FFT%BXS, 1,        1, 1, FFT%KTRN, NSCARC_INIT_ZERO, 'FFT%BXS', CROUTINE)
      ELSE
         CALL SCARC_ALLOCATE_REAL2(FFT%BXS, 1, FFT%JTRN, 1, 1, NSCARC_INIT_ZERO, 'FFT%BXS', CROUTINE)
      ENDIF

      ! Allocate boundary data vector for XF

      IF (L%NZ>1) THEN
         IF (L%NY  > 1) CALL SCARC_ALLOCATE_REAL2(FFT%BXF, 1, FFT%JTRN, 1, FFT%KTRN, NSCARC_INIT_ZERO, 'FFT%BXF', CROUTINE)
         IF (L%NY == 1) CALL SCARC_ALLOCATE_REAL2(FFT%BXF, 1,        1, 1, FFT%KTRN, NSCARC_INIT_ZERO, 'FFT%BXF', CROUTINE)
      ELSE
         CALL SCARC_ALLOCATE_REAL2(FFT%BXF, 1, FFT%JTRN, 1, 1, NSCARC_INIT_ZERO, 'FFT%BXF', CROUTINE)
      ENDIF

      ! Allocate boundary data vector for YS

      IF (L%NZ > 1) THEN
         CALL SCARC_ALLOCATE_REAL2(FFT%BYS, 1, FFT%ITRN, 1, FFT%KTRN, NSCARC_INIT_ZERO, 'FFT%BYS', CROUTINE)
      ELSE
         CALL SCARC_ALLOCATE_REAL2(FFT%BYS, 1, FFT%ITRN, 1,        1, NSCARC_INIT_ZERO, 'FFT%BYS', CROUTINE)
      ENDIF

      ! Allocate boundary data vector for YF

      IF (L%NZ > 1) THEN
         CALL SCARC_ALLOCATE_REAL2(FFT%BYF, 1, FFT%ITRN,1, FFT%KTRN, NSCARC_INIT_ZERO, 'FFT%BYF', CROUTINE)
      ELSE
         CALL SCARC_ALLOCATE_REAL2(FFT%BYF, 1, FFT%ITRN,1,        1, NSCARC_INIT_ZERO, 'FFT%BYF', CROUTINE)
      ENDIF

      ! Allocate boundary data vector for ZS

      IF (L%NY  > 1) CALL SCARC_ALLOCATE_REAL2(FFT%BZS, 1, FFT%ITRN, 1, FFT%JTRN, NSCARC_INIT_ZERO, 'FFT%BZS', CROUTINE)
      IF (L%NY == 1) CALL SCARC_ALLOCATE_REAL2(FFT%BZS, 1, FFT%ITRN, 1,        1, NSCARC_INIT_ZERO, 'FFT%BZS', CROUTINE)

      ! Allocate boundary data vector for ZF

      IF (L%NY  >1)  CALL SCARC_ALLOCATE_REAL2(FFT%BZF, 1, FFT%ITRN, 1, FFT%JTRN, NSCARC_INIT_ZERO, 'FFT%BZF', CROUTINE)
      IF (L%NY == 1) CALL SCARC_ALLOCATE_REAL2(FFT%BZF, 1, FFT%ITRN, 1,        1, NSCARC_INIT_ZERO, 'FFT%BZF', CROUTINE)

      IF (TWO_D) THEN
         CALL H2CZIS(FFT%XS,FFT%XF,FFT%IBAR,FFT%LBC,&
                     FFT%ZS,FFT%ZF,FFT%KBAR,FFT%NBC,&
                     FFT%HX,FFT%XLM,FFT%ITRN,IERROR,FFT%SAVE1)
      ELSE
         CALL H3CZIS(FFT%XS,FFT%XF,FFT%IBAR,FFT%LBC,&
                     FFT%YS,FFT%YF,FFT%JBAR,FFT%MBC,&
                     FFT%ZS,FFT%ZF,FFT%KBAR,FFT%NBC,&
                     FFT%HX,FFT%XLM,FFT%ITRN,FFT%JTRN,IERROR,FFT%SAVE1)
      ENDIF

   ENDDO LEVEL_LOOP
ENDDO MESHES_LOOP

END SUBROUTINE SCARC_SETUP_FFTO




#ifdef WITH_MKL
! ------------------------------------------------------------------------------------------------
!> \brief Initialize CLUSTER_SPARSE_SOLVER from MKL-library
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_CLUSTER(NLMIN, NLMAX)
USE SCARC_POINTERS, ONLY: L, G, AS, MKL
USE SCARC_POINTER_ROUTINES, ONLY: SCARC_POINT_TO_GRID
INTEGER, INTENT(IN) :: NLMIN, NLMAX
INTEGER :: NM, NL, I 
REAL (EB) :: TNOW
REAL (EB) :: DUMMY(1)=0.0_EB
REAL (FB) :: DUMMY_FB(1)=0.0_FB

TNOW = CURRENT_TIME()
CROUTINE = 'SCARC_SETUP_CLUSTER'

MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   LEVEL_LOOP: DO NL = NLMIN, NLMAX

      CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G

      MKL => L%MKL
      AS  => G%POISSON_SYM

      ! Allocate workspace for parameters and pointers needed in MKL-routine
 
      CALL SCARC_ALLOCATE_INT1(MKL%IPARM, 1, 64, NSCARC_INIT_ZERO, 'MKL%IPARM', CROUTINE)

      IF (.NOT.ALLOCATED(MKL%CT)) THEN
         ALLOCATE(MKL%CT(64), STAT=IERROR)
         CALL CHKMEMERR ('SCARC', 'CT', IERROR)
         DO I=1,64
            MKL%CT(I)%DUMMY = 0
         ENDDO
      ENDIF

      ! Define corresponding parameters
      ! Note: IPARM-vectory is allocate from 1:64, not from 0:63
 
      MKL%NRHS   =  1         ! one right hand side
      MKL%MAXFCT =  1         ! one matrix
      MKL%MNUM   =  1         ! number of matrix to be factorized
      MKL%ERROR  =  0         ! initialize error flag
      MKL%MSGLVL =  0         ! do not print statistical information
      IF (SCARC_MKL_MTYPE == 'SYMMETRIC') THEN
         MKL%MTYPE  = -2      ! Matrix type real and symmetric indefinite
      ELSE
         MKL%MTYPE  = 11      ! Matrix type real and non-symmetric
      ENDIF
      MKL%IPARM(1)  =  1      ! no solver default
      IF (N_MPI_PROCESSES > 4) THEN 
         MKL%IPARM(2) =10     ! 10 = MPI Parallel fill-in reordering from METIS. If 3 = OpenMP parallel reordering in Master
      ELSE                    ! Note IPARM(2)=10 has a bug which has been fixed from Intel MKL 2018 update 2 onwards.
         MKL%IPARM(2) = 3
      ENDIF
      MKL%IPARM(4)  = 0       ! no iterative-direct algorithm
      MKL%IPARM(5)  = 0       ! no user fill-in reducing permutation
      MKL%IPARM(6)  = 0       ! =0 solution on the first n components of x
      MKL%IPARM(8)  = 2       ! numbers of iterative refinement steps
      MKL%IPARM(10) = 13      ! perturb the pivot elements with 1E-13
      MKL%IPARM(11) = 1       ! use nonsymmetric permutation and scaling MPS  !!!!! was 1
      MKL%IPARM(13) = 1       ! maximum weighted matching algorithm is switched-off
                              !(default for symmetric). Try iparm(13) = 1 in case of inappropriate accuracy
      MKL%IPARM(14) = 0       ! Output: number of perturbed pivots
      MKL%IPARM(18) = 0       ! Output: number of nonzeros in the factor LU
      MKL%IPARM(19) = 0       ! Output: Mflops for LU factorization
      MKL%IPARM(20) = 0       ! Output: Numbers of CG Iterations
      MKL%IPARM(21) = 1       ! 1x1 diagonal pivoting for symmetric indefinite matrices.
      MKL%IPARM(24) = 0
      MKL%IPARM(27) = 1       ! Check matrix
      MKL%IPARM(40) = 2       ! Matrix, solution and rhs provided in distributed assembled matrix input format.

      MKL%IPARM(41) = G%NC_OFFSET(NM) + 1                      ! first global cell number for mesh NM
      MKL%IPARM(42) = G%NC_OFFSET(NM) + G%NC_LOCAL(NM)         ! last global cell number for mesh NM
      !MKL%IPARM(39) = 2                                       ! provide matrix in distributed format
      !MKL%IPARM(40) = G%NC_OFFSET(NM)+1                       ! first global cell number for mesh NM
      !MKL%IPARM(41) = G%NC_OFFSET(NM)+G%NC_LOCAL(NM)         ! last global cell number for mesh NM
 
      ! First perform only reordering and symbolic factorization
      ! Then perform only factorization
 
      IF (TYPE_MKL_PRECISION == NSCARC_PRECISION_SINGLE) THEN

         MKL%IPARM(28) = 1         ! single precision
         MKL%PHASE = 11
         CALL CLUSTER_SPARSE_SOLVER_S(MKL%CT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, G%NC_GLOBAL, &
                                      AS%VAL_FB, AS%ROW, AS%COL, MKL%PERM, MKL%NRHS, MKL%IPARM, &
                                      MKL%MSGLVL, DUMMY_FB, DUMMY_FB, MPI_COMM_WORLD, MKL%ERROR)
         MKL%PHASE = 22
         CALL CLUSTER_SPARSE_SOLVER_S(MKL%CT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, G%NC_GLOBAL, &
                                      AS%VAL_FB, AS%ROW, AS%COL, MKL%PERM, MKL%NRHS, MKL%IPARM, &
                                      MKL%MSGLVL, DUMMY_FB, DUMMY_FB, MPI_COMM_WORLD, MKL%ERROR)

         IF (MKL%ERROR /= 0) THEN
            WRITE(*,*) 'ERROR in MKL SETUP, MKL%ERROR=', MKL%ERROR
            CALL MPI_FINALIZE(IERROR)
            STOP
         ENDIF

      ELSE

         MKL%IPARM(28) = 0         ! double precision
         MKL%PHASE = 11
         CALL CLUSTER_SPARSE_SOLVER_D(MKL%CT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, G%NC_GLOBAL, &
                                      AS%VAL, AS%ROW, AS%COL, MKL%PERM, MKL%NRHS, MKL%IPARM, &
                                      MKL%MSGLVL, DUMMY, DUMMY, MPI_COMM_WORLD, MKL%ERROR)
         MKL%PHASE = 22
         CALL CLUSTER_SPARSE_SOLVER_D(MKL%CT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, G%NC_GLOBAL, &
                                      AS%VAL, AS%ROW, AS%COL, MKL%PERM, MKL%NRHS, MKL%IPARM, &
                                      MKL%MSGLVL, DUMMY, DUMMY, MPI_COMM_WORLD, MKL%ERROR)
         IF (MKL%ERROR /= 0) THEN
            WRITE(*,*) 'ERROR in MKL SETUP, MKL%ERROR=', MKL%ERROR
            CALL MPI_FINALIZE(IERROR)
            STOP
         ENDIF

      ENDIF

   ENDDO LEVEL_LOOP
ENDDO MESHES_LOOP

END SUBROUTINE SCARC_SETUP_CLUSTER


! ------------------------------------------------------------------------------------------------
!> \brief Initialize PARDISO solver from MKL-library
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_PARDISO(NLMIN, NLMAX)
USE SCARC_POINTERS, ONLY: L, G, AS, MKL
USE SCARC_POINTER_ROUTINES, ONLY: SCARC_POINT_TO_GRID
INTEGER, INTENT(IN) :: NLMIN, NLMAX
INTEGER :: NM, NL, I, IDUMMY(1)=0
REAL (EB) :: TNOW
REAL (EB) :: DUMMY(1)=0.0_EB
REAL (FB) :: DUMMY_FB(1)=0.0_FB

TNOW = CURRENT_TIME()
CROUTINE = 'SCARC_SETUP_PARDISO'

MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   LEVEL_LOOP: DO NL = NLMIN, NLMAX

      CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G
      MKL => L%MKL
      AS  => G%POISSON_SYM

      ! Allocate workspace for parameters nnd pointers eeded in MKL-routine
 
      CALL SCARC_ALLOCATE_INT1(MKL%IPARM, 1, 64, NSCARC_INIT_ZERO, 'MKL%IPARM', CROUTINE)

      IF (.NOT.ALLOCATED(MKL%PT)) THEN
         ALLOCATE(MKL%PT(64), STAT=IERROR)
         CALL CHKMEMERR ('SCARC', 'PT', IERROR)
         DO I=1,64
            MKL%PT(I)%DUMMY = 0
         ENDDO
      ENDIF

      ! Define corresponding parameters
      ! Note: IPARM-vectory is allocate from 1:64, not from 0:63
 
      MKL%NRHS   = 1
      MKL%MAXFCT = 1
      MKL%MNUM   = 1

      MKL%IPARM(1)  =  1      ! no solver default
      MKL%IPARM(4)  =  0      ! factorization computed as required by phase
      MKL%IPARM(5)  =  0      ! user permutation ignored
      MKL%IPARM(6)  =  0      ! write solution on x
      MKL%IPARM(8)  =  2      ! numbers of iterative refinement steps
      MKL%IPARM(10) = 13      ! perturb the pivot elements with 1E-13
      MKL%IPARM(11) =  0      ! disable scaling (default for SPD)
      MKL%IPARM(13) =  0      ! disable matching
      MKL%IPARM(18) = -1      ! Output: number of nonzeros in the factor LU
      MKL%IPARM(19) = -1      ! Output: number of floating points operations
      MKL%IPARM(20) =  1      ! Output: Numbers of CG Iterations
      MKL%IPARM(27) =  1      ! use matrix checker
      MKL%IPARM(37) =  0      ! matrix storage in COMPACT-format
      MKL%IPARM(40) = 2       ! Matrix, solution and rhs provided in distributed assembled matrix input format.

      MKL%ERROR  =  0         ! initialize error flag
      MKL%MSGLVL =  0         ! do not print statistical information
      MKL%MTYPE  = -2         ! Matrix type real non-symmetric

      ! First perform only reordering and symbolic factorization
      ! Then perform only factorization
 

      IF (TYPE_MKL_PRECISION == NSCARC_PRECISION_SINGLE) THEN
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'SETUP_PARDISO SINGLE: G%NC=',G%NC
CALL SCARC_DEBUG_CMATRIX(AS, 'AS','PARDISO SETUP')
#endif
         MKL%IPARM(28) = 1         ! single precision
         MKL%PHASE = 11
         CALL PARDISO_S(MKL%PT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, G%NC, &
                        AS%VAL_FB, AS%ROW, AS%COL, IDUMMY, &
                        MKL%NRHS, MKL%IPARM, MKL%MSGLVL, DUMMY_FB, DUMMY_FB, MKL%ERROR)
         MKL%PHASE = 22
         CALL PARDISO_S(MKL%PT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, G%NC, &
                        AS%VAL_FB, AS%ROW, AS%COL, IDUMMY, &
                        MKL%NRHS, MKL%IPARM, MKL%MSGLVL, DUMMY_FB, DUMMY_FB, MKL%ERROR)
      ELSE
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'SETUP_PARDISO DOUBLE: G%NC=',G%NC
CALL SCARC_DEBUG_CMATRIX(AS, 'AS','PARDISO SETUP')
#endif
         MKL%IPARM(28) = 0         ! double precision
         MKL%PHASE = 11
         CALL PARDISO_D(MKL%PT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, G%NC, &
                        AS%VAL, AS%ROW, AS%COL, IDUMMY, &
                        MKL%NRHS, MKL%IPARM, MKL%MSGLVL, DUMMY, DUMMY, MKL%ERROR)
         MKL%PHASE = 22
         CALL PARDISO_D(MKL%PT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, G%NC, &
                        AS%VAL, AS%ROW, AS%COL, IDUMMY, &
                        MKL%NRHS, MKL%IPARM, MKL%MSGLVL, DUMMY, DUMMY, MKL%ERROR)
      ENDIF

   ENDDO LEVEL_LOOP
ENDDO MESHES_LOOP

END SUBROUTINE SCARC_SETUP_PARDISO
#endif

! ------------------------------------------------------------------------------------------------
!> \brief Define sizes for system matrix A (including extended regions related to overlaps)
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_POISSON_SIZES(NL)
USE SCARC_POINTERS, ONLY: S, L, G, OG, A, OA, AB, OAB
USE SCARC_POINTER_ROUTINES, ONLY: SCARC_POINT_TO_GRID, SCARC_POINT_TO_OTHER_GRID, &
                                  SCARC_POINT_TO_CMATRIX, SCARC_POINT_TO_OTHER_CMATRIX, &
                                  SCARC_POINT_TO_BMATRIX, SCARC_POINT_TO_OTHER_BMATRIX
INTEGER, INTENT(IN) :: NL
INTEGER :: NM, NOM, INBR

MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   
   CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G
   
   SELECT_MATRIX_TYPE: SELECT CASE (SCARC_GET_MATRIX_TYPE(NL))
   
 
      ! -------- Matrix in compact storage technique
 
      CASE (NSCARC_MATRIX_COMPACT)

         A => SCARC_POINT_TO_CMATRIX (G, NSCARC_MATRIX_POISSON)

         ! Assign IOR settings to corresponding positions in stencil

         IF (TWO_D) THEN
            A%N_STENCIL = 5
            A%POS(-3:3) = (/1,0,2,3,4,0,5/)     
         ELSE
            A%N_STENCIL = 7
            A%POS(-3:3) = (/1,2,3,4,5,6,7/)
         ENDIF

         A%N_VAL  = G%NCE * A%N_STENCIL
         A%N_ROW  = G%NCE + 1

         ! Allocate matrices on overlapping parts for later data exchanges with neighbors

         DO INBR = 1, SCARC(NM)%N_NEIGHBORS
            NOM = S%NEIGHBORS(INBR)
            CALL SCARC_POINT_TO_OTHER_GRID(NM, NOM, NL)
            OA => SCARC_POINT_TO_OTHER_CMATRIX(OG, NSCARC_MATRIX_POISSON)
            OA%N_STENCIL = A%N_STENCIL
            OA%N_VAL = 4 * OG%NCG * A%N_STENCIL            ! TODO: CHECK LENGTH
            OA%N_ROW = OG%NCG + 1
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'POISSON_SIZES: INBR, NVAL, NROW:', INBR, OA%N_VAL, OA%N_ROW
#endif
         ENDDO

 
      ! -------- Matrix in bandwise storage technique
 
      CASE (NSCARC_MATRIX_BANDWISE)

         AB => SCARC_POINT_TO_BMATRIX (G, NSCARC_MATRIX_POISSON)

         IF (TWO_D) THEN
   
            AB%N_STENCIL   = 5                      ! 5-point Laplacian
            AB%POS(-3:3)   = (/5,0,4,3,2,0,1/)      ! assignment of IOR settings to columns in matrix array
   
            AB%OFFSET( 3)  = -L%NX                  ! lower z
            AB%OFFSET( 1)  = -1                     ! lower x
            AB%OFFSET( 0)  =  0                     ! diag
            AB%OFFSET(-1)  =  1                     ! upper x
            AB%OFFSET(-3)  =  L%NX                  ! upper z
   
            AB%SOURCE( 3)   =  1                    ! lower z
            AB%SOURCE( 1)   =  1                    ! lower x
            AB%SOURCE( 0)   =  1                    ! diag
            AB%SOURCE(-1)   =  2                    ! upper x
            AB%SOURCE(-3)   =  L%NX+1               ! upper z
   
            AB%TARGET( 3)   =  L%NX+1               ! lower z
            AB%TARGET( 1)   =  2                    ! lower x
            AB%TARGET( 0)   =  1                    ! diag
            AB%TARGET(-1)   =  1                    ! upper x
            AB%TARGET(-3)   =  1                    ! upper z
   
            AB%LENGTH( 3)  =  G%NC - L%NX           ! lower z
            AB%LENGTH( 1)  =  G%NC - 1              ! lower x
            AB%LENGTH( 0)  =  G%NC                  ! diag
            AB%LENGTH(-1)  =  G%NC - 1              ! upper x
            AB%LENGTH(-3)  =  G%NC - L%NX           ! upper z
   
         ELSE
   
            AB%N_STENCIL   = 7                      ! 7-point Laplacian
            AB%POS(-3:3)   = (/7,6,5,4,3,2,1/)      ! assignment of IOR settings to columns in matrix array
   
            AB%OFFSET( 3)  = -L%NX*L%NY             ! lower z
            AB%OFFSET( 2)  = -L%NX                  ! lower y
            AB%OFFSET( 1)  = -1                     ! lower x
            AB%OFFSET( 0)  =  0                     ! diag
            AB%OFFSET(-1)  =  1                     ! upper x
            AB%OFFSET(-2)  =  L%NX                  ! upper y
            AB%OFFSET(-3)  =  L%NX*L%NY             ! upper z
   
            AB%SOURCE( 3)  =  1                     ! lower z
            AB%SOURCE( 2)  =  1                     ! lower y
            AB%SOURCE( 1)  =  1                     ! lower x
            AB%SOURCE( 0)  =  1                     ! diag
            AB%SOURCE(-1)  =  2                     ! upper x
            AB%SOURCE(-2)  =  L%NX+1                ! upper y
            AB%SOURCE(-3)  =  L%NX*L%NY+1           ! upper z
   
            AB%TARGET( 3)  =  L%NX*L%NY+1           ! lower z
            AB%TARGET( 2)  =  L%NX+1                ! lower y
            AB%TARGET( 1)  =  2                     ! lower x
            AB%TARGET( 0)  =  1                     ! diag
            AB%TARGET(-1)  =  1                     ! upper x
            AB%TARGET(-2)  =  1                     ! upper y
            AB%TARGET(-3)  =  1                     ! upper z
   
            AB%LENGTH( 3)  =  G%NC - L%NX*L%NY      ! lower z
            AB%LENGTH( 2)  =  G%NC - L%NX           ! lower y
            AB%LENGTH( 1)  =  G%NC - 1              ! lower x
            AB%LENGTH( 0)  =  G%NC                  ! diag
            AB%LENGTH(-1)  =  G%NC - 1              ! upper x
            AB%LENGTH(-2)  =  G%NC - L%NX           ! upper y
            AB%LENGTH(-3)  =  G%NC - L%NX*L%NY      ! upper z
   
         ENDIF

         AB%N_VAL  = G%NC * AB%N_STENCIL
         AB%N_DIAG = G%NC

         ! Determine sizes of overlapping parts for later communication with corresponding neighbors

         DO INBR = 1, SCARC(NM)%N_NEIGHBORS
            NOM = S%NEIGHBORS(INBR)
            CALL SCARC_POINT_TO_OTHER_GRID(NM, NOM, NL)
            OAB => SCARC_POINT_TO_OTHER_BMATRIX(OG, NSCARC_MATRIX_POISSON)
            OAB%N_STENCIL = AB%N_STENCIL
            OAB%N_VAL     = OG%NCG * AB%N_STENCIL
            OAB%N_DIAG    = OG%NCG 
         ENDDO

   END SELECT SELECT_MATRIX_TYPE
   
ENDDO MESHES_LOOP
   
 
! -------- Exchange matrix sizes in case of a multi-mesh geometry
 
IF (NMESHES > 1) CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_MATRIX_SIZES, NSCARC_MATRIX_POISSON, NL)

END SUBROUTINE SCARC_SETUP_POISSON_SIZES


! ------------------------------------------------------------------------------------------------------
!> \brief Basic call of ScaRC solver routine
! ------------------------------------------------------------------------------------------------------
!SUBROUTINE SCARC_SOLVER(DT_CURRENT, PRES_ITE, TOTAL_PRES_ITE)
SUBROUTINE SCARC_SOLVER(DT_CURRENT)
USE VELO, ONLY: NO_FLUX
!INTEGER, INTENT(IN) :: PRES_ITE, TOTAL_PRES_ITE
REAL (EB), INTENT(IN) :: DT_CURRENT
REAL (EB) :: TNOW
INTEGER :: IBAR, KBAR

TNOW = CURRENT_TIME()

DT  = DT_CURRENT
DTI = 1.0_EB/DT_CURRENT

ITE_PRES = ITE_PRES + 1
ITE_GLOBAL = ICYC

IBAR = MESHES(MYID+1)%IBAR
KBAR = MESHES(MYID+1)%KBAR

#ifdef WITH_SCARC_DEBUG
CALL SCARC_DEBUG_METHOD('STARTING SCARC',3)        
WRITE(MSG%LU_DEBUG,*) 'TYPE_GRID   =', TYPE_GRID  
WRITE(MSG%LU_DEBUG,*) 'TYPE_METHOD =', TYPE_METHOD
WRITE(MSG%LU_DEBUG,*) 'PREDICTOR =', PREDICTOR
WRITE(MSG%LU_DEBUG,*) 'PRES_ON_WHOLE_DOMAIN =', PRES_ON_WHOLE_DOMAIN
#endif

SELECT_METHOD: SELECT CASE (TYPE_METHOD)

   ! ---------------- Krylov method (CG) -------------------------------------
   CASE (NSCARC_METHOD_KRYLOV)
 
      CALL SCARC_METHOD_KRYLOV (NSCARC_STACK_ROOT, NSCARC_STACK_ZERO, NSCARC_RHS_INHOMOGENEOUS, NLEVEL_MIN)
   
 
   ! ---------------- Multigrid method ---------------------------------------
 
   CASE (NSCARC_METHOD_MULTIGRID)
   
      CALL SCARC_METHOD_MULTIGRID(NSCARC_STACK_ROOT, NSCARC_STACK_ZERO, NSCARC_RHS_INHOMOGENEOUS, NLEVEL_MIN)
   

   ! ---------------- McKenny-Greengard-Mayo method (MGM) --------------------
 
   CASE (NSCARC_METHOD_MGM)
   
      CALL SCARC_METHOD_MGM

   ! ---------------- MKL method ---------------------------------------------
 
#ifdef WITH_MKL
   CASE (NSCARC_METHOD_LU)
   
      SELECT_MKL: SELECT CASE (TYPE_MKL(0))
         CASE (NSCARC_MKL_GLOBAL)
            CALL SCARC_METHOD_CLUSTER(NSCARC_STACK_ROOT, NSCARC_STACK_ZERO, NLEVEL_MIN)
         CASE (NSCARC_MKL_LOCAL)
            CALL SCARC_METHOD_PARDISO(NSCARC_STACK_ROOT, NSCARC_STACK_ZERO, NLEVEL_MIN)
      END SELECT SELECT_MKL
#endif
   
END SELECT SELECT_METHOD

#ifdef WITH_SCARC_DEBUG
CALL SCARC_DEBUG_METHOD('LEAVING SCARC',3)      
#endif

IF (STOP_STATUS==SETUP_STOP) RETURN

T_USED(5)=T_USED(5)+CURRENT_TIME()-TNOW
CPU(MYID)%SOLVER =CPU(MYID)%SOLVER+CURRENT_TIME()-TNOW
CPU(MYID)%OVERALL=CPU(MYID)%OVERALL+CURRENT_TIME()-TNOW
END SUBROUTINE SCARC_SOLVER


! ------------------------------------------------------------------------------------------------------
!> \brief Perform preceding FFT method to improve start solution for ScaRC
! ------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_METHOD_FFT
USE MESH_POINTERS
USE POIS, ONLY: H2CZSS, H3CZSS
REAL(EB), POINTER, DIMENSION(:,:,:) :: HP
INTEGER :: NM, I, J, K
LOGICAL :: WITH_BDRY = .FALSE.

MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   !CALL POINT_TO_MESH(NM)
   
   IF (PREDICTOR) THEN
      HP => H
   ELSE
      HP => HS
   ENDIF
   
   ! Call the Poisson solver
 
   IF (.NOT.TWO_D) CALL H3CZSS(BXS,BXF,BYS,BYF,BZS,BZF,ITRN,JTRN,PRHS,POIS_PTB,SAVE1,WORK,HX)
   IF (TWO_D .AND. .NOT.CYLINDRICAL) CALL H2CZSS(BXS,BXF,BZS,BZF,ITRN,PRHS,POIS_PTB,SAVE1,WORK,HX)
   
   DO K=1,KBAR
     DO J=1,JBAR
        DO I=1,IBAR
            HP(I,J,K) = PRHS(I,J,K)
        ENDDO
      ENDDO
   ENDDO
   
   ! Apply boundary conditions to H
 
   IF (WITH_BDRY) THEN
      DO K=1,KBAR
         DO J=1,JBAR
            IF (LBC==3 .OR. LBC==4)             HP(0,J,K)    = HP(1,J,K)    - DXI*BXS(J,K)
            IF (LBC==3 .OR. LBC==2 .OR. LBC==6) HP(IBP1,J,K) = HP(IBAR,J,K) + DXI*BXF(J,K)
            IF (LBC==1 .OR. LBC==2)             HP(0,J,K)    =-HP(1,J,K)    + 2._EB*BXS(J,K)
            IF (LBC==1 .OR. LBC==4 .OR. LBC==5) HP(IBP1,J,K) =-HP(IBAR,J,K) + 2._EB*BXF(J,K)
            IF (LBC==5 .OR. LBC==6)             HP(0,J,K)    = HP(1,J,K)
            IF (LBC==0) THEN
               HP(0,J,K) = HP(IBAR,J,K)
               HP(IBP1,J,K) = HP(1,J,K)
            ENDIF
         ENDDO
      ENDDO
      
      DO K=1,KBAR
         DO I=1,IBAR
            IF (MBC==3 .OR. MBC==4) HP(I,0,K)    = HP(I,1,K)    - DETA*BYS(I,K)
            IF (MBC==3 .OR. MBC==2) HP(I,JBP1,K) = HP(I,JBAR,K) + DETA*BYF(I,K)
            IF (MBC==1 .OR. MBC==2) HP(I,0,K)    =-HP(I,1,K)    + 2._EB*BYS(I,K)
            IF (MBC==1 .OR. MBC==4) HP(I,JBP1,K) =-HP(I,JBAR,K) + 2._EB*BYF(I,K)
            IF (MBC==0) THEN
               HP(I,0,K) = HP(I,JBAR,K)
               HP(I,JBP1,K) = HP(I,1,K)
            ENDIF
         ENDDO
      ENDDO
      
      DO J=1,JBAR
         DO I=1,IBAR
            IF (NBC==3 .OR. NBC==4)  HP(I,J,0)    = HP(I,J,1)    - DZETA*BZS(I,J)
            IF (NBC==3 .OR. NBC==2)  HP(I,J,KBP1) = HP(I,J,KBAR) + DZETA*BZF(I,J)
            IF (NBC==1 .OR. NBC==2)  HP(I,J,0)    =-HP(I,J,1)    + 2._EB*BZS(I,J)
            IF (NBC==1 .OR. NBC==4)  HP(I,J,KBP1) =-HP(I,J,KBAR) + 2._EB*BZF(I,J)
            IF (NBC==0) THEN
               HP(I,J,0) = HP(I,J,KBAR)
               HP(I,J,KBP1) = HP(I,J,1)
            ENDIF
         ENDDO
      ENDDO
   ENDIF

ENDDO MESHES_LOOP

END SUBROUTINE SCARC_METHOD_FFT




#ifdef WITH_MKL
! ------------------------------------------------------------------------------------------------
!> \brief Perform global Pardiso-method based on MKL
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_METHOD_CLUSTER(NSTACK, NPARENT, NLEVEL)
USE SCARC_POINTERS, ONLY: L, G, MKL, V1, V2, AS, V1_FB, V2_FB
USE SCARC_POINTER_ROUTINES, ONLY: SCARC_POINT_TO_GRID, SCARC_POINT_TO_VECTOR, SCARC_POINT_TO_VECTOR_FB, &
                                  SCARC_POINT_TO_CMATRIX
INTEGER, INTENT(IN) :: NSTACK, NPARENT, NLEVEL
INTEGER ::  NM, NS, NP, NL
REAL (EB) :: TNOW

NS = NSTACK
NP = NPARENT
NL = NLEVEL

TNOW = CURRENT_TIME()

CALL SCARC_SETUP_SOLVER(NS, NP)
CALL SCARC_SETUP_WORKSPACE(NS, NL, NSCARC_RHS_INHOMOGENEOUS)

MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G
   AS => SCARC_POINT_TO_CMATRIX(G, NSCARC_MATRIX_POISSON_SYM)

   V1 => SCARC_POINT_TO_VECTOR (NM, NL, B)
   V2 => SCARC_POINT_TO_VECTOR (NM, NL, X)

   MKL => L%MKL
   MKL%PHASE  = 33                                ! only solving

#ifdef WITH_SCARC_DEBUG2
WRITE(MSG%LU_DEBUG,*) 'G%NC_GLOBAL=', G%NC_GLOBAL
WRITE(MSG%LU_DEBUG,*) 'CLUSTER, PRE, V1:'
WRITE(MSG%LU_DEBUG,'(6E14.6)') V1
WRITE(MSG%LU_DEBUG,*) 'CLUSTER, PRE, V2:'
WRITE(MSG%LU_DEBUG,'(6E14.6)') V2
CALL SCARC_DEBUG_CMATRIX(AS, 'AS','CLUSTER')
#endif

   IF (TYPE_MKL_PRECISION == NSCARC_PRECISION_SINGLE) THEN

      V1_FB => SCARC_POINT_TO_VECTOR_FB (NM, NL, B)
      V2_FB => SCARC_POINT_TO_VECTOR_FB (NM, NL, X)

      V1_FB = REAL(V1, FB)
      V2_FB = 0.0_FB
      CALL CLUSTER_SPARSE_SOLVER_S(MKL%CT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, G%NC_GLOBAL, &
                                   AS%VAL_FB, AS%ROW, AS%COL, MKL%PERM, MKL%NRHS, MKL%IPARM, &
                                   MKL%MSGLVL, V1_FB, V2_FB, MPI_COMM_WORLD, MKL%ERROR)
      V2 = REAL(V2_FB, EB)

   ELSE

      CALL CLUSTER_SPARSE_SOLVER_D(MKL%CT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, G%NC_GLOBAL, &
                                   AS%VAL, AS%ROW, AS%COL, MKL%PERM, MKL%NRHS, MKL%IPARM, &
                                   MKL%MSGLVL, V1, V2, MPI_COMM_WORLD, MKL%ERROR)
   ENDIF

   IF (MKL%ERROR /= 0) CALL SCARC_SHUTDOWN(NSCARC_ERROR_MKL_INTERNAL, SCARC_NONE, MKL%ERROR)

#ifdef WITH_SCARC_DEBUG2
WRITE(MSG%LU_DEBUG,*) 'CLUSTER, POST, V1:'
WRITE(MSG%LU_DEBUG,'(2E14.6)') V1
WRITE(MSG%LU_DEBUG,*) 'CLUSTER, POST, V2:'
WRITE(MSG%LU_DEBUG,'(2E14.6)') V2
#endif
ENDDO MESHES_LOOP

CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_VECTOR_PLAIN, X, NL)

IF (TYPE_SOLVER == NSCARC_SOLVER_MAIN) THEN
   CALL SCARC_UPDATE_MAINCELLS (NLEVEL_MIN)
   CALL SCARC_UPDATE_GHOSTCELLS(NLEVEL_MIN)
ENDIF

CALL SCARC_RELEASE_SOLVER(NS, NP)

END SUBROUTINE SCARC_METHOD_CLUSTER
#endif


#ifdef WITH_MKL
! ------------------------------------------------------------------------------------------------
!> \brief Perform global Pardiso-method based on MKL
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_METHOD_PARDISO(NSTACK, NPARENT, NLEVEL)
USE SCARC_POINTERS, ONLY: L, G, MKL, AS, V1, V2, V1_FB, V2_FB
USE SCARC_POINTER_ROUTINES, ONLY: SCARC_POINT_TO_GRID, SCARC_POINT_TO_VECTOR, SCARC_POINT_TO_VECTOR_FB, &
                                  SCARC_POINT_TO_CMATRIX
INTEGER, INTENT(IN) :: NSTACK, NPARENT, NLEVEL
INTEGER ::  NM, NS, NP, NL
REAL (EB) :: TNOW

TNOW = CURRENT_TIME()

NS = NSTACK
NP = NPARENT
NL = NLEVEL

CALL SCARC_SETUP_SOLVER(NS, NP)
CALL SCARC_SETUP_WORKSPACE(NS, NL, NSCARC_RHS_INHOMOGENEOUS)

MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G
   AS => SCARC_POINT_TO_CMATRIX(G, NSCARC_MATRIX_POISSON_SYM)

   V1 => SCARC_POINT_TO_VECTOR (NM, NL, B)
   V2 => SCARC_POINT_TO_VECTOR (NM, NL, X)

   MKL => L%MKL
   MKL%PHASE  = 33         ! only solving

   IF (TYPE_MKL_PRECISION == NSCARC_PRECISION_SINGLE) THEN

#ifdef WITH_SCARC_DEBUG2
WRITE(MSG%LU_DEBUG,*) 'PARDISO SINGLE, PRE, V1:', G%NC, SIZE(V1)
WRITE(MSG%LU_DEBUG,'(6E14.6)') V1
WRITE(MSG%LU_DEBUG,*) 'PARDISO SINGLE, PRE, V2:', G%NC, SIZE(V2)
WRITE(MSG%LU_DEBUG,'(6E14.6)') V2
!CALL SCARC_DEBUG_CMATRIX(AS, 'AS','PARDISO')
#endif

      V1_FB => SCARC_POINT_TO_VECTOR_FB (NM, NL, B)
      V2_FB => SCARC_POINT_TO_VECTOR_FB (NM, NL, X)

      V1_FB = REAL(V1, FB)
      V2_FB = 0.0_FB
      CALL PARDISO_S(MKL%PT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, G%NC, &
                     AS%VAL_FB, AS%ROW, AS%COL, MKL%PERM, MKL%NRHS, MKL%IPARM, &
                     MKL%MSGLVL, V1_FB, V2_FB, MKL%ERROR)

      V2 = REAL(V2_FB, EB)

   ELSE

#ifdef WITH_SCARC_DEBUG2
WRITE(MSG%LU_DEBUG,*) 'PARDISO DOUBLE, PRE, V1:', G%NC, SIZE(V1)
WRITE(MSG%LU_DEBUG,'(6E14.6)') V1
WRITE(MSG%LU_DEBUG,*) 'PARDISO DOUBLE, PRE, V2:', G%NC, SIZE(V2)
WRITE(MSG%LU_DEBUG,'(6E14.6)') V2
!CALL SCARC_DEBUG_CMATRIX(AS, 'AS','PARDISO')
#endif

      V2 = 0.0_EB
      CALL PARDISO_D(MKL%PT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, G%NC, &
                     AS%VAL, AS%ROW, AS%COL, MKL%PERM, MKL%NRHS, MKL%IPARM, &
                     MKL%MSGLVL, V1, V2, MKL%ERROR)
   ENDIF

   IF (MKL%ERROR /= 0) CALL SCARC_SHUTDOWN(NSCARC_ERROR_MKL_INTERNAL, SCARC_NONE, MKL%ERROR)

#ifdef WITH_SCARC_DEBUG2
WRITE(MSG%LU_DEBUG,*) 'PARDISO, POST, V1:'
WRITE(MSG%LU_DEBUG,'(2E14.6)') V1
WRITE(MSG%LU_DEBUG,*) 'PARDISO, POST, V2:'
WRITE(MSG%LU_DEBUG,'(2E14.6)') V2
#endif
ENDDO MESHES_LOOP

IF (TYPE_SOLVER == NSCARC_SOLVER_MAIN) THEN
   CALL SCARC_UPDATE_MAINCELLS (NLEVEL_MIN)
   CALL SCARC_UPDATE_GHOSTCELLS(NLEVEL_MIN)
ENDIF

CALL SCARC_RELEASE_SOLVER(NSTACK, NPARENT)

END SUBROUTINE SCARC_METHOD_PARDISO
#endif


! ------------------------------------------------------------------------------------------------
!> \brief Perform requested coarse grid solver (iterative/direct)
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_METHOD_COARSE(NSTACK, NPARENT, NLEVEL)
INTEGER, INTENT(IN) :: NSTACK, NPARENT, NLEVEL

SELECT CASE (TYPE_COARSE)

   CASE (NSCARC_COARSE_ITERATIVE)
      CALL SCARC_METHOD_KRYLOV (NSTACK, NPARENT, NSCARC_RHS_DEFECT, NLEVEL)

   CASE (NSCARC_COARSE_DIRECT)
#ifdef WITH_MKL
      !IF (STACK(NPARENT)%SOLVER%TYPE_SCOPE(0) == NSCARC_SCOPE_GLOBAL .AND. NMESHES > 1) THEN
      IF (NMESHES > 1) THEN
         CALL SCARC_METHOD_CLUSTER (NSTACK, NPARENT, NLEVEL)
      ELSE
         CALL SCARC_METHOD_PARDISO (NSTACK, NPARENT, NLEVEL)
      ENDIF
#else
      CALL SCARC_SHUTDOWN(NSCARC_ERROR_DIRECT_NOMKL, SCARC_NONE, NLEVEL)
#endif

END SELECT

END SUBROUTINE SCARC_METHOD_COARSE




! ------------------------------------------------------------------------------------------------
!> \brief  Setup environement for current solver 
! i.e. set pointers to used vectors related to current position in stack
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_SOLVER(NS, NP)
USE SCARC_POINTERS, ONLY: SV, SVP
INTEGER, INTENT(IN) :: NS, NP                          ! references to current stack and parent
 
! If not first solver in stack, store last iteration parameters of parent solver NP
 
IF (NP > 0) THEN
   SVP => STACK(NP)%SOLVER
   SVP%ITE   = ITE
   SVP%RES   = RES
   SVP%RESIN = RESIN
   SVP%CAPPA = CAPPA
ENDIF
 
! Set new environment for solver on stack position NS
 
SV => STACK(NS)%SOLVER

CNAME = SV%CNAME
ITE   = 0
NIT   = SV%NIT
EPS   = SV%EPS
OMEGA = SV%OMEGA
RESIN = SV%RESIN
CAPPA = -1.0

TYPE_PARENT = NP

TYPE_ACCURACY       = SV%TYPE_ACCURACY
TYPE_COARSE         = SV%TYPE_COARSE
TYPE_COARSENING     = SV%TYPE_COARSENING
TYPE_CYCLING        = SV%TYPE_CYCLING
TYPE_GRID           = SV%TYPE_GRID
TYPE_EXCHANGE       = SV%TYPE_EXCHANGE
TYPE_INTERPOL       = SV%TYPE_INTERPOL
TYPE_LEVEL          = SV%TYPE_LEVEL
TYPE_MATRIX         = SV%TYPE_MATRIX
TYPE_METHOD         = SV%TYPE_METHOD
TYPE_MKL            = SV%TYPE_MKL
TYPE_MKL_PRECISION  = SV%TYPE_MKL_PRECISION
TYPE_MULTIGRID      = SV%TYPE_MULTIGRID
TYPE_PARENT         = SV%TYPE_PARENT
TYPE_PRECON         = SV%TYPE_PRECON
TYPE_RELAX          = SV%TYPE_RELAX
TYPE_SCOPE          = SV%TYPE_SCOPE
TYPE_SMOOTH         = SV%TYPE_SMOOTH
TYPE_SOLVER         = SV%TYPE_SOLVER
TYPE_STAGE          = SV%TYPE_STAGE
TYPE_STENCIL        = SV%TYPE_STENCIL
TYPE_TWOLEVEL       = SV%TYPE_TWOLEVEL
TYPE_VECTOR         = SV%TYPE_VECTOR

X = SV%X
B = SV%B
D = SV%D
R = SV%R
V = SV%V
Y = SV%Y
Z = SV%Z

#ifdef WITH_SCARC_DEBUG
E = SV%E
#endif

IF (TYPE_SOLVER == NSCARC_SOLVER_MAIN) ITE_TOTAL = 0

END SUBROUTINE SCARC_SETUP_SOLVER


! ------------------------------------------------------------------------------------------------
!> \brief Reset environment of calling routine when leaving solver
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_RELEASE_SOLVER(NS, NP)
USE SCARC_POINTERS, ONLY: SV, SVP
INTEGER, INTENT(IN)  :: NS, NP                            ! references to current stack and parent

SV  => STACK(NS)%SOLVER

! Store convergence information of preceding solver for FDS dump routine
 
IF (TYPE_SOLVER == NSCARC_SOLVER_MAIN) THEN
   SCARC_CAPPA      = CAPPA
   SCARC_RESIDUAL   = RES
   SCARC_ITERATIONS = ITE
   !CALL SCARC_CLOSE_CSV_FILE()
ENDIF

SV%RESIN = RESIN
SV%RES   = RES
SV%ITE   = ITE
SV%CAPPA = CAPPA

! If not first solver in stack, reset environment of parent (calling) routine
 
IF (NP > 0) THEN

   SVP => STACK(NP)%SOLVER

   ITE   = SVP%ITE
   NIT   = SVP%NIT
   EPS   = SVP%EPS
   RESIN = SVP%RESIN
   RES   = SVP%RES
   OMEGA = SVP%OMEGA
   CAPPA = SVP%CAPPA

   TYPE_ACCURACY       = SVP%TYPE_ACCURACY
   TYPE_COARSE         = SVP%TYPE_COARSE
   TYPE_COARSENING     = SVP%TYPE_COARSENING
   TYPE_CYCLING        = SVP%TYPE_CYCLING
   TYPE_GRID           = SVP%TYPE_GRID
   TYPE_EXCHANGE       = SVP%TYPE_EXCHANGE
   TYPE_INTERPOL       = SVP%TYPE_INTERPOL
   TYPE_LEVEL          = SVP%TYPE_LEVEL
   TYPE_MATRIX         = SVP%TYPE_MATRIX
   TYPE_METHOD         = SVP%TYPE_METHOD
   TYPE_MKL            = SVP%TYPE_MKL
   TYPE_MKL_PRECISION  = SVP%TYPE_MKL_PRECISION
   TYPE_MULTIGRID      = SVP%TYPE_MULTIGRID
   TYPE_PARENT         = SVP%TYPE_PARENT
   TYPE_PRECON         = SVP%TYPE_PRECON
   TYPE_RELAX          = SVP%TYPE_RELAX
   TYPE_SCOPE          = SVP%TYPE_SCOPE
   TYPE_SMOOTH         = SVP%TYPE_SMOOTH
   TYPE_SOLVER         = SVP%TYPE_SOLVER
   TYPE_STAGE          = SVP%TYPE_STAGE
   TYPE_STENCIL        = SVP%TYPE_STENCIL
   TYPE_TWOLEVEL       = SVP%TYPE_TWOLEVEL
   TYPE_VECTOR         = SVP%TYPE_VECTOR

   X = SVP%X
   B = SVP%B
   D = SVP%D
   R = SVP%R
   V = SVP%V
   Y = SVP%Y
   Z = SVP%Z

#ifdef WITH_SCARC_DEBUG
   E = SVP%E
#endif

ENDIF

END SUBROUTINE SCARC_RELEASE_SOLVER


! ----------------------------------------------------------------------------------------------------
!> \brief Set initial solution corresponding to boundary data in BXS, BXF, ...
! ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_WORKSPACE(NS, NL, NRHS)
USE SCARC_POINTERS, ONLY: M, L, F, G, SV, ST, STP, GWC, PRHS, HP
USE SCARC_POINTER_ROUTINES, ONLY: SCARC_POINT_TO_GRID
#ifdef WITH_SCARC_POSTPROCESSING
USE SCARC_POINTERS, ONLY: PR
#endif
INTEGER, INTENT(IN) :: NS, NL, NRHS
REAL(EB) :: VAL
INTEGER  :: NM, IW, IW1, IW2, IOR0, I, J, K, IC
LOGICAL  :: BFIRST_WORKSPACE = .FALSE.

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'STARTING SETUP_WORKSPACE ', NS, NL, NRHS
#endif

SV  => STACK(NS)%SOLVER

SELECT_SOLVER_TYPE: SELECT CASE (SV%TYPE_SOLVER)

         ! ---------- If used as main solver use values from pressure-routine as initialization
 
         CASE (NSCARC_SOLVER_MAIN)
         
            MAIN_MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
         
               CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G
               ST => SCARC(NM)%LEVEL(NL)%STAGE(SV%TYPE_STAGE)
         
               PRHS => M%PRHS
               IF (PREDICTOR) THEN
                  HP => M%H
               ELSE
                  HP => M%HS
               ENDIF


#ifdef WITH_SCARC_POSTPROCESSING
               PR => L%PRESSURE
               IF (PREDICTOR) THEN
                  PR%H_OLD = PR%H_NEW
               ELSE
                  PR%HS_OLD = PR%HS_NEW
               ENDIF
               PR%B_OLD = ST%B
#endif
         
               ! Get right hand side (PRHS from pres.f90) and initial vector (H or HS from last time step)

               SELECT_RHS_TYPE: SELECT CASE (NRHS)
         
                  ! Solve original problem with inhomegeneous boundary conditions
                  CASE (NSCARC_RHS_INHOMOGENEOUS)
            
                     IF (IS_MGM) BFIRST_WORKSPACE = .TRUE.

                     !$OMP PARALLEL DO PRIVATE(IC) SCHEDULE(STATIC)
                     DO IC = 1, G%NC
                        ST%X(IC) = HP(G%ICX(IC), G%ICY(IC), G%ICZ(IC))        ! use last iterate as initial solution
                        ST%B(IC) = PRHS(G%ICX(IC), G%ICY(IC), G%ICZ(IC))      ! get new RHS from surrounding code
                     ENDDO                         
                     !$OMP END PARALLEL DO
                     ST%X = 0.0_EB                      ! CAUTION - ONLY TEMPORARILY - TODO
            
                     !!$OMP PARALLEL 
                     MAIN_INHOMOGENEOUS_LOOP: DO IOR0 = -3, 3, 1 
            
                        IF (IOR0 == 0) CYCLE
                        F => SCARC(NM)%LEVEL(NL)%FACE(IOR0)
                        
                        IW1 = F%NCW0
                        IW2 = F%NCW0 + F%NCW - 1
            
                        !!$OMP DO PRIVATE(IW, GWC, I, J, K, IC, VAL) SCHEDULE(STATIC)
                        FACE_INHOMOGENEOUS_LOOP: DO IW = IW1, IW2
            
                           GWC => G%WALL(IW)
               
                           I = GWC%IXW
                           J = GWC%IYW
                           K = GWC%IZW
               
                           IF (TWO_D .AND. J /= 1) CALL SCARC_SHUTDOWN(NSCARC_ERROR_GRID_INDEX, SCARC_NONE, J)
               
                           IF (IS_UNSTRUCTURED .AND. L%IS_SOLID(I, J, K)) CYCLE
               
                           IC = G%CELL_NUMBER(I,J,K)
               
                           ! ---------- Dirichlet BC's:
                           ! these are based on the SETTING in BTYPE
                           ! in the structured case this corresponds to the face-wise SETTING according to the FFT
                           ! (this allows to use local FFT's as preconditioners)
                           ! in the unstructured case only open boundary cells lead to Dirichlet BC's

                           IF_DIRICHLET: IF (GWC%BTYPE == DIRICHLET) THEN
               
                              SELECT CASE (IOR0)
                              CASE (1)
                                 VAL =  M%BXS(J,K)
                              CASE (-1)
                                 VAL =  M%BXF(J,K)
                              CASE (2)
                                 VAL =  M%BYS(I,K)
                              CASE (-2)
                                 VAL =  M%BYF(I,K)
                              CASE (3)
                                 VAL =  M%BZS(I,J)
                              CASE (-3)
                                 VAL =  M%BZF(I,J)
                              END SELECT
               
                              ST%B(IC) = ST%B(IC) + F%SCAL_DIRICHLET * VAL
#ifdef WITH_SCARC_DEBUG2
WRITE(MSG%LU_DEBUG,'(A, 5I6,2E14.6)') 'SETUP_WORKSPACE: DIRICHLET: IW, I, J, K, IC, VAL, B(IC):', &
                                             IW, I, J, K, IC, VAL, ST%B(IC)
#endif
               
                           ENDIF IF_DIRICHLET
               
                           ! ---------- Neumann BC's:
                           ! Note for the unstructured case only:
                           ! Here, the matrix also contains Neumann BC's for those cells which have a
                           ! PRESSURE_BC_INDEX == DIRICHLET but are NOT open; these cells must be excluded below,
                           ! because BXS, BXF, ... contain the Dirichlet information from pres.f90 there;
                           ! excluding them corresponds to a homogeneous Neumann condition for these cells

                           IF_NEUMANN: IF (GWC%BTYPE == NEUMANN) THEN
               
                              IF (IS_UNSTRUCTURED .AND. M%WALL(IW)%PRESSURE_BC_INDEX /= NEUMANN) CYCLE
               
                              SELECT CASE (IOR0)
                              CASE (1)
                                 VAL =  M%BXS(J,K)
                              CASE (-1)
                                 VAL =  M%BXF(J,K)
                              CASE (2)
                                 VAL =  M%BYS(I,K)
                              CASE (-2)
                                 VAL =  M%BYF(I,K)
                              CASE (3)
                                 VAL =  M%BZS(I,J)
                              CASE (-3)
                                 VAL =  M%BZF(I,J)
                              END SELECT
               
                              ST%B(IC) = ST%B(IC) + F%SCAL_NEUMANN * VAL

#ifdef WITH_SCARC_DEBUG2
WRITE(MSG%LU_DEBUG,'(A, 5I6,2E14.6)') 'SETUP_WORKSPACE: NEUMANN  : IW, I, J, K, IC, VAL, B(IC):', &
                                             IW, I, J, K, IC, VAL, ST%B(IC)
#endif
               
                           ENDIF IF_NEUMANN
            
                        ENDDO FACE_INHOMOGENEOUS_LOOP
                        !!$OMP END DO

                     ENDDO MAIN_INHOMOGENEOUS_LOOP
                     !!$OMP END PARALLEL 
         

                  ! Solve problem with homegeneous boundary conditions (MGM only)

                  CASE (NSCARC_RHS_HOMOGENEOUS)
         
                     ST%B = 0.0_EB                                    ! set RHS to zero
                     ST%X = 0.0_EB                                    ! use zero as initial vector
                     ST%V = 0.0_EB                                 
                     ST%D = 0.0_EB                                
                     ST%R = 0.0_EB                               
                     ST%Y = 0.0_EB                              
                     ST%Z = 0.0_EB                             
  
                     IF (NMESHES > 1) CALL SCARC_SETUP_MGM_INTERFACES(NM, NL)         ! setup BC's along mesh interfaces
                     CALL SCARC_SETUP_MGM_OBSTRUCTIONS                                ! setup BC's along internal obstructions

#ifdef WITH_SCARC_DEBUG
                     CALL SCARC_DEBUG_LEVEL_MESH(ST%B, 'RHS second pass MGM', NSCARC_GRID_UNSTRUCTURED, NM, NL)
#endif
                  BFIRST_WORKSPACE = .FALSE.
         
               END SELECT SELECT_RHS_TYPE
 

#ifdef WITH_SCARC_POSTPROCESSING
               PR%B_NEW = ST%B
#endif
         
            ENDDO MAIN_MESHES_LOOP
            
            ! In case of a Krylov method clear overlapping parts of auxiliary vectors

            IF (IS_CG.OR.IS_MGM.OR.HAS_TWO_LEVELS) THEN
               DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
                  L  => SCARC(NM)%LEVEL(NL)
                  ST => SCARC(NM)%LEVEL(NL)%STAGE(SV%TYPE_STAGE)
                  ST%D(1:G%NCE) = 0.0_EB
                  ST%R(1:G%NCE) = 0.0_EB
                  ST%V(1:G%NCE) = 0.0_EB
                  ST%Y(1:G%NCE) = 0.0_EB
                  ST%Z(1:G%NCE) = 0.0_EB
               ENDDO
            ENDIF
         
            ! In case of a multigrid method as main solver clear
            ! overlapping parts of auxiliary vectors and coarse grid solver vectors

            IF (IS_GMG) THEN
               DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
                  L  => SCARC(NM)%LEVEL(NL)
                  ST => SCARC(NM)%LEVEL(NL)%STAGE(SV%TYPE_STAGE)
                  ST%V(1:G%NCE) = 0.0_EB
                  ST%Z(1:G%NCE) = 0.0_EB
               ENDDO
            ENDIF
         
            ! In case of pure Neumann or periodic BCs, broadcast RHS(end) from last mesh
            ! to all and store it on all meshes

            IF (N_DIRIC_GLOBAL(NLEVEL_MIN) == 0) THEN
               IF (UPPER_MESH_INDEX == NMESHES) THEN
                  L  => SCARC(NMESHES)%LEVEL(NL)
                  ST => SCARC(NMESHES)%LEVEL(NL)%STAGE(SV%TYPE_STAGE)
                  MESH_REAL = ST%B(G%NC)
               ELSE
                  MESH_REAL = 0.0_EB
               ENDIF
               IF (N_MPI_PROCESSES > 1) &
                  CALL MPI_ALLGATHER(MPI_IN_PLACE, 1, MPI_DOUBLE_PRECISION, MESH_REAL, 1, MPI_DOUBLE_PRECISION,&
                                     MPI_COMM_WORLD, IERROR)
               DO NM = 1, NMESHES
                  SCARC(NM)%RHS_END = MESH_REAL(NMESHES)
               ENDDO
            ENDIF
         
 
         ! ---------- If MG is used as Krylov preconditioner, vector R of main Krylov is the new RHS for MG
 
         CASE (NSCARC_SOLVER_PRECON)
         
            IF (IS_CG_MG) THEN
               DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
                  L   => SCARC(NM)%LEVEL(NL)
                  ST  => SCARC(NM)%LEVEL(NL)%STAGE(SV%TYPE_STAGE)              ! current stage
                  STP => SCARC(NM)%LEVEL(NL)%STAGE(NSCARC_STAGE_ONE)           ! parent stage
                  ST%X(1:G%NCE) = 0.0_EB
                  ST%B(1:G%NCE) = STP%R(1:G%NCE)                                                
                  ST%V(1:G%NCE) = 0.0_EB
                  ST%Z(1:G%NCE) = 0.0_EB
               ENDDO
            ENDIF
         
 
         ! ---------- If used as coarse grid solver start with zero initialization
 
         CASE (NSCARC_SOLVER_COARSE)
         
            DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
               ST => SCARC(NM)%LEVEL(NL)%STAGE(SV%TYPE_STAGE)
               ST%X = 0.0_EB
               ST%D = 0.0_EB
               ST%R = 0.0_EB
               ST%V = 0.0_EB
               ST%Y = 0.0_EB
               ST%Z = 0.0_EB
            ENDDO
         
END SELECT SELECT_SOLVER_TYPE

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'LEAVING SETUP_WORKSPACE ', NS, NL, NRHS
#endif
END SUBROUTINE SCARC_SETUP_WORKSPACE




! ------------------------------------------------------------------------------------------------
!> \brief Perform restriction from finer to coarser grid level
!    - 'VF' corresponds to vector on fine   grid
!    - 'VC' corresponds to vector on coarse grid
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_RESTRICTION (NVB, NVC, NLF, NLC)
USE SCARC_POINTERS, ONLY: LC, GF, GC, VF, VC, R
USE SCARC_POINTER_ROUTINES, ONLY: SCARC_POINT_TO_MULTIGRID, SCARC_POINT_TO_VECTOR, SCARC_POINT_TO_CMATRIX
INTEGER, INTENT(IN) :: NVB, NVC, NLF, NLC
REAL(EB) :: DSUM
INTEGER :: NM
INTEGER :: IXF, IYF, IZF, ICF(8)=0, ICFB(-2:2,-2:2)=0
INTEGER :: IXC, IYC, IZC, ICC, IC, ICOL

#ifdef WITH_SCARC_DEBUG
CALL SCARC_DEBUG_LEVEL (NVB, 'RESTRICTION INIT: NVB', NLF)
CALL SCARC_DEBUG_LEVEL (NVC, 'RESTRICTION INIT: NVC', NLC)
#endif
IF (IS_GMG .OR. IS_CG_GMG .OR. IS_CG_ADD .OR. IS_CG_MUL) THEN
 
! ---------- Twolevel-CG or Geometric multigrid (as main solver or preconditioner) 
 
   IF (HAS_MULTIPLE_LEVELS) THEN
   
      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   
         CALL SCARC_POINT_TO_MULTIGRID(NM, NLF, NLC)
   
         VF => SCARC_POINT_TO_VECTOR(NM, NLF, NVB)
         VC => SCARC_POINT_TO_VECTOR(NM, NLC, NVC)
   
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'SIZE(VC)=', SIZE(VC)
WRITE(MSG%LU_DEBUG,*) 'GC%NC=', GC%NC
WRITE(MSG%LU_DEBUG,*) 'GC%NCE=', GC%NCE
#endif
         IF (TWO_D) THEN
   
            SELECT_INTERPOL: SELECT CASE (TYPE_INTERPOL)
   
               ! ---------- Constant Interpolation
   
               CASE (NSCARC_INTERPOL_CONSTANT)
      
                  !$OMP PARALLEL DO PRIVATE(IXC, IZC, IXF, IZF, ICC, ICF) SCHEDULE(STATIC)
                  DO IZC = 1, LC%NZ
                     DO IXC = 1, LC%NX
      
                        IF (IS_UNSTRUCTURED .AND. LC%IS_SOLID(IXC, 1, IZC)) CYCLE
      
                        IXF = 2*IXC
                        IZF = 2*IZC
      
                        ICC = GC%CELL_NUMBER(IXC, 1, IZC)
      
                        ICF(1) = GF%CELL_NUMBER(IXF-1, 1, IZF-1)
                        ICF(2) = GF%CELL_NUMBER(IXF-1, 1, IZF  )
                        ICF(3) = GF%CELL_NUMBER(IXF  , 1, IZF-1)
                        ICF(4) = GF%CELL_NUMBER(IXF  , 1, IZF  )
      
                        VC(ICC) = 0.25_EB * (  VF(ICF(1)) &
                                             + VF(ICF(2)) &
                                             + VF(ICF(3)) &
                                             + VF(ICF(4)) )
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,'(A,7I6, 5E14.6)') 'REST: IXC, IZC, ICC, ICF(1:4), VC(ICC), 0.25*VF(ICF):', &
 IXC, IZC, ICC, ICF(1), ICF(2), ICF(3), ICF(4), &
 VC(ICC), VF(ICF(1)), VF(ICF(2)), VF(ICF(3)), VF(ICF(4))
#endif
                     ENDDO
                  ENDDO
                  !$OMP END PARALLEL DO
      
               ! ---------- Bilinear Interpolation
   
               CASE (NSCARC_INTERPOL_BILINEAR)
      
                  VC=0.0_EB
      
                  !$OMP PARALLEL DO PRIVATE(IXC, IYC, IZC, IXF, IYF, IZF, ICC, ICFB) SCHEDULE(STATIC)
                  DO IZC = 1, LC%NZ
                     DO IXC = 1, LC%NX
      
                        IF (IS_UNSTRUCTURED .AND. LC%IS_SOLID(IXC, 1, IZC)) CYCLE
      
                        IXF = 2*IXC
                        IZF = 2*IZC
      
                        ICC = GC%CELL_NUMBER(IXC, 1, IZC)
      
                        ICFB(-2,-2) = GF%CELL_NUMBER(IXF-2, 1, IZF-2)
                        ICFB(-1,-2) = GF%CELL_NUMBER(IXF-1, 1, IZF-2)
                        ICFB( 1,-2) = GF%CELL_NUMBER(IXF  , 1, IZF-2)
                        ICFB( 2,-2) = GF%CELL_NUMBER(IXF+1, 1, IZF-2)
      
                        ICFB(-2,-1) = GF%CELL_NUMBER(IXF-2, 1, IZF-1)
                        ICFB(-1,-1) = GF%CELL_NUMBER(IXF-1, 1, IZF-1)
                        ICFB( 1,-1) = GF%CELL_NUMBER(IXF  , 1, IZF-1)
                        ICFB( 2,-1) = GF%CELL_NUMBER(IXF+1, 1, IZF-1)
      
                        ICFB(-2, 1) = GF%CELL_NUMBER(IXF-2, 1, IZF)
                        ICFB(-1, 1) = GF%CELL_NUMBER(IXF-1, 1, IZF)
                        ICFB( 1, 1) = GF%CELL_NUMBER(IXF  , 1, IZF)
                        ICFB( 2, 1) = GF%CELL_NUMBER(IXF+1, 1, IZF)
      
                        ICFB(-2, 2) = GF%CELL_NUMBER(IXF-2, 1, IZF+1)
                        ICFB(-1, 2) = GF%CELL_NUMBER(IXF-1, 1, IZF+1)
                        ICFB( 1, 2) = GF%CELL_NUMBER(IXF  , 1, IZF+1)
                        ICFB( 2, 2) = GF%CELL_NUMBER(IXF+1, 1, IZF+1)
      
                        IF (IXC==1.AND.IZC==1) THEN
                           VC(ICC) = SCALR*( &
                              W4 *VF(ICFB(-1, 2)) + W3 *VF(ICFB(1, 2)) + W1*VF(ICFB(2, 2)) + &
                              W12*VF(ICFB(-1, 1)) + W9 *VF(ICFB(1, 1)) + W3*VF(ICFB(2, 1)) + &
                              W16*VF(ICFB(-1,-1)) + W12*VF(ICFB(1,-1)) + W4*VF(ICFB(2,-1)) )
                        ELSE IF (IXC==LC%NX.AND.IZC==  1) THEN
                           VC(ICC) = SCALR*( &
                              W1 *VF(ICFB(-2, 2)) + W3 *VF(ICFB(-1, 2)) + W4 *VF(ICFB(1, 2)) + &
                              W3 *VF(ICFB(-2, 1)) + W9 *VF(ICFB(-1, 1)) + W12*VF(ICFB(1, 1)) + &
                              W4 *VF(ICFB(-2,-1)) + W12*VF(ICFB(-1,-1)) + W16*VF(ICFB(1,-1)) )
                        ELSE IF (IXC==  1.AND.IZC==LC%NZ) THEN
                           VC(ICC) = SCALR*( &
                              W16*VF(ICFB(-1, 1)) + W12*VF(ICFB(1, 1)) + W4*VF(ICFB(2, 1)) + &
                              W12*VF(ICFB(-1,-1)) + W9 *VF(ICFB(1,-1)) + W3*VF(ICFB(2,-1)) + &
                              W4 *VF(ICFB(-1,-2)) + W3 *VF(ICFB(1,-2)) + W1*VF(ICFB(2,-2)) )
                        ELSE IF (IXC==LC%NX.AND.IZC==LC%NZ) THEN
                           VC(ICC) = SCALR*( &
                              W4 *VF(ICFB(-2, 1)) + W12*VF(ICFB(-1, 1)) + W16*VF(ICFB(1, 1)) + &
                              W3 *VF(ICFB(-2,-1)) + W9 *VF(ICFB(-1,-1)) + W12*VF(ICFB(1,-1)) + &
                              W1 *VF(ICFB(-2,-2)) + W3 *VF(ICFB(-1,-2)) + W4 *VF(ICFB(1,-2)) )
                        ELSE IF (IZC==  1) THEN
                           VC(ICC) = SCALR*( &
                              W1*VF(ICFB(-2, 2)) + W3 *VF(ICFB(-1, 2)) + W3 *VF(ICFB(1, 2)) + W1*VF(ICFB(2, 2)) + &
                              W3*VF(ICFB(-2, 1)) + W9 *VF(ICFB(-1, 1)) + W9 *VF(ICFB(1, 1)) + W3*VF(ICFB(2, 1)) + &
                              W4*VF(ICFB(-2,-1)) + W12*VF(ICFB(-1,-1)) + W12*VF(ICFB(1,-1)) + W4*VF(ICFB(2,-1)) )
                        ELSE IF (IZC==LC%NZ) THEN
                           VC(ICC) = SCALR*( &
                              W4*VF(ICFB(-2, 1)) + W12*VF(ICFB(-1, 1)) + W12*VF(ICFB(1, 1)) + W4*VF(ICFB(2, 1)) + &
                              W3*VF(ICFB(-2,-1)) + W9 *VF(ICFB(-1,-1)) + W9 *VF(ICFB(1,-1)) + W3*VF(ICFB(2,-1)) + &
                              W1*VF(ICFB(-2,-2)) + W3 *VF(ICFB(-1,-2)) + W3 *VF(ICFB(1,-2)) + W1*VF(ICFB(2,-2)) )
                        ELSE IF (IXC==  1) THEN
                           VC(ICC) = SCALR*( &
                              W4 *VF(ICFB(-1, 2)) + W3*VF(ICFB(1, 2)) + W1*VF(ICFB(2, 2)) +&
                              W12*VF(ICFB(-1, 1)) + W9*VF(ICFB(1, 1)) + W3*VF(ICFB(2, 1)) +&
                              W12*VF(ICFB(-1,-1)) + W9*VF(ICFB(1,-1)) + W3*VF(ICFB(2,-1)) +&
                              W4 *VF(ICFB(-1,-2)) + W3*VF(ICFB(1,-2)) + W1*VF(ICFB(2,-2)) )
                        ELSE IF (IXC==LC%NX) THEN
                           VC(ICC) = SCALR*( &
                              W1*VF(ICFB(-2, 2)) + W3*VF(ICFB(-1, 2)) + W4 *VF(ICFB(1, 2)) + &
                              W3*VF(ICFB(-2, 1)) + W9*VF(ICFB(-1, 1)) + W12*VF(ICFB(1, 1)) +&
                              W3*VF(ICFB(-2,-1)) + W9*VF(ICFB(-1,-1)) + W12*VF(ICFB(1,-1)) +&
                              W1*VF(ICFB(-2,-2)) + W3*VF(ICFB(-1,-2)) + W4 *VF(ICFB(1,-2)) )
                        ELSE
                           VC(ICC) = SCALR*( &
                              W1*VF(ICFB(-2,-2)) + W3*VF(ICFB(-1,-2)) + W3*VF(ICFB(1,-2)) + W1*VF(ICFB(2,-2)) +&
                              W3*VF(ICFB(-2,-1)) + W9*VF(ICFB(-1,-1)) + W9*VF(ICFB(1,-1)) + W3*VF(ICFB(2,-1)) +&
                              W3*VF(ICFB(-2, 1)) + W9*VF(ICFB(-1, 1)) + W9*VF(ICFB(1, 1)) + W3*VF(ICFB(2, 1)) +&
                              W1*VF(ICFB(-2, 2)) + W3*VF(ICFB(-1, 2)) + W3*VF(ICFB(1, 2)) + W1*VF(ICFB(2, 2)) )
                        ENDIF
                     ENDDO
                  ENDDO
                  !$OMP END PARALLEL DO
      
            END SELECT SELECT_INTERPOL
   
         ! ---------- Constant Interpolation (Note: 3D-bilinear case is still missing)
   
         ELSE
   
            !$OMP PARALLEL DO PRIVATE(IXC, IYC, IZC, IXF, IYF, IZF, ICC, ICF) SCHEDULE(STATIC)
            DO IZC = 1, LC%NZ
               DO IYC = 1, LC%NY
                  DO IXC = 1, LC%NX
   
                     IF (IS_UNSTRUCTURED .AND. LC%IS_SOLID(IXC, IYC, IZC)) CYCLE
   
                     IXF = 2*IXC
                     IYF = 2*IYC
                     IZF = 2*IZC
   
                     ICC = GC%CELL_NUMBER(IXC, IYC, IZC)
   
                     ICF(1) = GF%CELL_NUMBER(IXF-1, IYF-1, IZF-1)
                     ICF(2) = GF%CELL_NUMBER(IXF-1, IYF-1, IZF  )
                     ICF(3) = GF%CELL_NUMBER(IXF-1, IYF  , IZF-1)
                     ICF(4) = GF%CELL_NUMBER(IXF-1, IYF  , IZF  )
                     ICF(5) = GF%CELL_NUMBER(IXF  , IYF-1, IZF-1)
                     ICF(6) = GF%CELL_NUMBER(IXF  , IYF-1, IZF  )
                     ICF(7) = GF%CELL_NUMBER(IXF  , IYF  , IZF-1)
                     ICF(8) = GF%CELL_NUMBER(IXF  , IYF  , IZF  )
   
                     VC(ICC) = 0.125_EB * (  VF(ICF(1)) &
                                           + VF(ICF(2)) &
                                           + VF(ICF(3)) &
                                           + VF(ICF(4)) &
                                           + VF(ICF(5)) &
                                           + VF(ICF(6)) &
                                           + VF(ICF(7)) &
                                           + VF(ICF(8)) )
   
                  ENDDO
               ENDDO
            ENDDO
            !$OMP END PARALLEL DO
   
         ENDIF
      ENDDO
   ENDIF
   
! ---------- Use restriction based on smoothed aggregation method

ELSE

   IF (NMESHES > 1) CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_VECTOR_PLAIN, NVB, NLF)

#ifdef WITH_SCARC_DEBUG
CALL SCARC_DEBUG_LEVEL (NVB, 'RESTRICTION INIT2: NVB', NLF)
#endif
   DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

      CALL SCARC_POINT_TO_MULTIGRID(NM, NLF, NLC)

      VF => SCARC_POINT_TO_VECTOR(NM, NLF, NVB)
      VC => SCARC_POINT_TO_VECTOR(NM, NLC, NVC)
   
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'SIZE(VC)=', SIZE(VC)
WRITE(MSG%LU_DEBUG,*) 'GC%NC=', GC%NC
WRITE(MSG%LU_DEBUG,*) 'GC%NCE=', GC%NCE
WRITE(MSG%LU_DEBUG,*) 'GF%N_COARSE=', GF%N_COARSE
#endif
      R => SCARC_POINT_TO_CMATRIX(GF, NSCARC_MATRIX_RESTRICTION)

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'RESTRICTION_AMG: NM, NLF, NLC, N_COARSE:', NM, NLF, NLC, GF%N_COARSE
#endif

      DO IC = 1, GF%N_COARSE
         DSUM = 0.0_EB
         DO ICOL = R%ROW(IC), R%ROW(IC+1)-1                            
            DSUM =  DSUM + R%VAL(ICOL) * VF(R%COLG(ICOL))
         ENDDO
         VC(IC) = DSUM
      ENDDO

   ENDDO

ENDIF
#ifdef WITH_SCARC_DEBUG
CALL SCARC_DEBUG_LEVEL (NVC, 'RESTRICTION EXIT: NVC', NLC)
#endif

END SUBROUTINE SCARC_RESTRICTION


! ------------------------------------------------------------------------------------------------
!> \brief Perform prolongation from coarser to finer grid level
!    - 'VC' corresponds to coarser grid
!    - 'VF' corresponds to finer   grid
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_PROLONGATION (NVC, NVB, NLC, NLF)
USE SCARC_POINTERS, ONLY: LC, GF, GC, VF, VC, P
USE SCARC_POINTER_ROUTINES, ONLY: SCARC_POINT_TO_MULTIGRID, SCARC_POINT_TO_VECTOR, SCARC_POINT_TO_CMATRIX
INTEGER, INTENT(IN) :: NVC, NVB, NLC, NLF
REAL(EB) :: DSUM
INTEGER :: NM, I
INTEGER :: IXF, IYF, IZF, ICF(8)=0, ICFB(-1:1,-1:1)=0
INTEGER :: IXC, IYC, IZC, ICC, IC, ICOL

#ifdef WITH_SCARC_DEBUG
CALL SCARC_DEBUG_LEVEL (NVB, 'PROLONGATION INIT: NVB', NLF)
CALL SCARC_DEBUG_LEVEL (NVC, 'PROLONGATION INIT: NVC', NLC)
#endif
IF (IS_GMG .OR. IS_CG_GMG .OR. IS_CG_ADD .OR. IS_CG_MUL) THEN
 
! ------------------ Twolevel CG or Geometric Multigrid 
 
   IF (HAS_MULTIPLE_LEVELS) THEN
   
      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   
         CALL SCARC_POINT_TO_MULTIGRID(NM, NLF, NLC)
   
         VC => SCARC_POINT_TO_VECTOR(NM, NLC, NVC)
         VF => SCARC_POINT_TO_VECTOR(NM, NLF, NVB)
   
         IF (TWO_D) THEN
   
            SELECT_INTERPOL: SELECT CASE (TYPE_INTERPOL)
   
               CASE (NSCARC_INTERPOL_CONSTANT)
      
                  !!$OMP PARALLEL DO PRIVATE(IXC, IZC, IXF, IZF, ICC, ICF) SCHEDULE(STATIC)
                  DO IZC = 1, LC%NZ
                     DO IXC = 1, LC%NX
      
                        IF (IS_UNSTRUCTURED .AND. LC%IS_SOLID(IXC, 1, IZC)) CYCLE
      
                        IXF = 2*IXC
                        IZF = 2*IZC
      
                        ICC = GC%CELL_NUMBER(IXC, 1, IZC)
      
                        ICF(1) = GF%CELL_NUMBER(IXF-1, 1, IZF-1)
                        ICF(2) = GF%CELL_NUMBER(IXF-1, 1, IZF  )
                        ICF(3) = GF%CELL_NUMBER(IXF  , 1, IZF-1)
                        ICF(4) = GF%CELL_NUMBER(IXF  , 1, IZF  )
      
                        DO I = 1, 4
                           VF(ICF(I)) = VC(ICC)
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,'(A,4I6, 2E14.6)') 'PROL: IXC, IZC, ICC, ICF(I), VF(ICF(I)),VC(ICC):', &
   IXC, IZC, ICC, ICF(I), VF(ICF(I)), VC(ICC)
#endif
                        ENDDO
                     ENDDO
                  ENDDO
                  !!$OMP END PARALLEL DO 
      
               CASE (NSCARC_INTERPOL_BILINEAR)
      
                  !!$OMP PARALLEL DO PRIVATE(IXC, IZC, IXF, IZF, ICC, ICFB) SCHEDULE(STATIC)
                  DO IZC = 1, LC%NZ
                     DO IXC = 1, LC%NX
      
                        IF (IS_UNSTRUCTURED .AND. LC%IS_SOLID(IXC, 1, IZC)) CYCLE
      
                        IXF = 2*IXC
                        IZF = 2*IZC
      
                        ICC = GC%CELL_NUMBER(IXC, 1, IZC)
      
                        ICFB(-1,-1) = GF%CELL_NUMBER(IXF-1, 1, IZF-1)
                        ICFB(-1, 1) = GF%CELL_NUMBER(IXF-1, 1, IZF  )
                        ICFB( 1,-1) = GF%CELL_NUMBER(IXF  , 1, IZF-1)
                        ICFB( 1, 1) = GF%CELL_NUMBER(IXF  , 1, IZF  )
      
                        IF (IXC==1.AND.IZC==1) THEN
                           VF(ICFB(-1,-1)) = VC(ICC)
                           VF(ICFB(-1, 1)) = SCALP*(W12*VC(ICC)+W4*VC(ICC+LC%NX))
                           VF(ICFB( 1,-1)) = SCALP*(W12*VC(ICC)+W4*VC(ICC+1))
                           VF(ICFB( 1, 1)) = SCALP*(W9 *VC(ICC)+W3*VC(ICC+1)+W3*VC(ICC+LC%NX)+W1*VC(ICC+LC%NX+1))
                        ELSE IF (IXC==1 .AND. IZC==LC%NZ) THEN
                           VF(ICFB(-1,-1)) = SCALP*(W12*VC(ICC)+W4*VC(ICC-LC%NX))
                           VF(ICFB(-1, 1)) = VC(ICC)
                           VF(ICFB( 1,-1)) = SCALP*(W9 *VC(ICC)+W3*VC(ICC+1)+W3*VC(ICC-LC%NX)+W1*VC(ICC-LC%NX+1))
                           VF(ICFB( 1, 1)) = SCALP*(W12*VC(ICC)+W4*VC(ICC+1))
                        ELSE IF (IXC==LC%NX .AND. IZC==1) THEN
                           VF(ICFB(-1,-1)) = SCALP*(W12*VC(ICC)+W4*VC(ICC-1))
                           VF(ICFB(-1, 1)) = SCALP*(W9 *VC(ICC)+W3*VC(ICC-1)+W3*VC(ICC+LC%NX)+W1*VC(ICC+LC%NX-1))
                           VF(ICFB( 1,-1)) = VC(ICC)
                           VF(ICFB( 1, 1)) = SCALP*(W12*VC(ICC)+W4*VC(ICC+LC%NX))
                        ELSE IF (IXC==LC%NX .AND. IZC==LC%NZ) THEN
                           VF(ICFB(-1,-1)) = SCALP*(W9 *VC(ICC)+W3*VC(ICC-1)+W3*VC(ICC-LC%NX)+W1*VC(ICC-LC%NX-1))
                           VF(ICFB(-1, 1)) = SCALP*(W12*VC(ICC)+W4*VC(ICC-1))
                           VF(ICFB( 1,-1)) = SCALP*(W12*VC(ICC)+W4*VC(ICC-LC%NX))
                           VF(ICFB( 1, 1)) = VC(ICC)
                        ELSE IF (IZC==1) THEN
                           VF(ICFB(-1,-1)) = SCALP*(W12*VC(ICC)+W4*VC(ICC-1))
                           VF(ICFB(-1, 1)) = SCALP*(W9 *VC(ICC)+W3*VC(ICC-1)+W3*VC(ICC+LC%NX)+W1*VC(ICC+LC%NX-1))
                           VF(ICFB( 1,-1)) = SCALP*(W12*VC(ICC)+W4*VC(ICC+1))
                           VF(ICFB( 1, 1)) = SCALP*(W9 *VC(ICC)+W3*VC(ICC+1)+W3*VC(ICC+LC%NX)+W1*VC(ICC+LC%NX+1))
                        ELSE IF (IZC==LC%NZ) THEN
                           VF(ICFB(-1,-1)) = SCALP*(W9 *VC(ICC)+W3*VC(ICC-1)+W3*VC(ICC-LC%NX)+W1*VC(ICC-LC%NX-1))
                           VF(ICFB(-1, 1)) = SCALP*(W12*VC(ICC)+W4*VC(ICC-1))
                           VF(ICFB( 1,-1)) = SCALP*(W9 *VC(ICC)+W3*VC(ICC+1)+W3*VC(ICC-LC%NX)+W1*VC(ICC-LC%NX+1))
                           VF(ICFB( 1, 1)) = SCALP*(W12*VC(ICC)+W4*VC(ICC+1))
                        ELSE IF (IXC==1) THEN
                           VF(ICFB(-1,-1)) = SCALP*(W12*VC(ICC)+W4*VC(ICC-LC%NX))
                           VF(ICFB(-1, 1)) = SCALP*(W12*VC(ICC)+W4*VC(ICC+LC%NX))
                           VF(ICFB( 1,-1)) = SCALP*(W9 *VC(ICC)+W3*VC(ICC+1)+W3*VC(ICC-LC%NX)+W1*VC(ICC-LC%NX+1))
                           VF(ICFB( 1, 1)) = SCALP*(W9 *VC(ICC)+W3*VC(ICC+1)+W3*VC(ICC+LC%NX)+W1*VC(ICC+LC%NX+1))
                        ELSE IF (IXC==LC%NX) THEN
                           VF(ICFB(-1,-1)) = SCALP*(W9 *VC(ICC)+W3*VC(ICC-1)+W3*VC(ICC-LC%NX)+W1*VC(ICC-LC%NX-1))
                           VF(ICFB(-1, 1)) = SCALP*(W9 *VC(ICC)+W3*VC(ICC-1)+W3*VC(ICC+LC%NX)+W1*VC(ICC+LC%NX-1))
                           VF(ICFB( 1,-1)) = SCALP*(W12*VC(ICC)+W4*VC(ICC-LC%NX))
                           VF(ICFB( 1, 1)) = SCALP*(W12*VC(ICC)+W4*VC(ICC+LC%NX))
                        ELSE
                           VF(ICFB(-1,-1)) = SCALP*(W9*VC(ICC)+W3*VC(ICC-1)+W3*VC(ICC-LC%NX)+W1*VC(ICC-LC%NX-1))
                           VF(ICFB(-1, 1)) = SCALP*(W9*VC(ICC)+W3*VC(ICC-1)+W3*VC(ICC+LC%NX)+W1*VC(ICC+LC%NX-1))
                           VF(ICFB( 1,-1)) = SCALP*(W9*VC(ICC)+W3*VC(ICC+1)+W3*VC(ICC-LC%NX)+W1*VC(ICC-LC%NX+1))
                           VF(ICFB( 1, 1)) = SCALP*(W9*VC(ICC)+W3*VC(ICC+1)+W3*VC(ICC+LC%NX)+W1*VC(ICC+LC%NX+1))
                        ENDIF
                     ENDDO
                  ENDDO
                  !!$OMP END PARALLEL DO 
      
            END SELECT SELECT_INTERPOL
   
         ELSE
   
            ! ---------- Constant Interpolation (Note: 3D-bilinear case is still missing)
   
            !!$OMP PARALLEL DO PRIVATE(IXC, IYC, IZC, IXF, IYF, IZF, ICC, ICF) SCHEDULE(STATIC)
            DO IZC = 1, LC%NZ
               DO IYC = 1, LC%NY
                  DO IXC = 1, LC%NX
   
                     IF (IS_UNSTRUCTURED .AND. LC%IS_SOLID(IXC, IYC, IZC)) CYCLE
   
                     IXF = 2*IXC
                     IYF = 2*IYC
                     IZF = 2*IZC
   
                     ICC = GC%CELL_NUMBER(IXC, IYC, IZC)
   
                     ICF(1) = GF%CELL_NUMBER(IXF-1, IYF-1, IZF-1)
                     ICF(2) = GF%CELL_NUMBER(IXF-1, IYF-1, IZF  )
                     ICF(3) = GF%CELL_NUMBER(IXF-1, IYF  , IZF-1)
                     ICF(4) = GF%CELL_NUMBER(IXF-1, IYF  , IZF  )
                     ICF(5) = GF%CELL_NUMBER(IXF  , IYF-1, IZF-1)
                     ICF(6) = GF%CELL_NUMBER(IXF  , IYF-1, IZF  )
                     ICF(7) = GF%CELL_NUMBER(IXF  , IYF  , IZF-1)
                     ICF(8) = GF%CELL_NUMBER(IXF  , IYF  , IZF  )
   
                     DO I = 1, 8
                        VF(ICF(I)) = VC(ICC)
                     ENDDO
   
                  ENDDO
               ENDDO
            ENDDO
            !!$OMP END PARALLEL DO 
   
         ENDIF
      ENDDO
   ENDIF
   
! ---------- Use Prolongation matrix based on smoothed aggregation method
ELSE

   DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

      CALL SCARC_POINT_TO_MULTIGRID(NM, NLF, NLC)

      VC => SCARC_POINT_TO_VECTOR(NM, NLC, NVC)
      VF => SCARC_POINT_TO_VECTOR(NM, NLF, NVB)

      P => SCARC_POINT_TO_CMATRIX(GF, NSCARC_MATRIX_PROLONGATION)

      !DO IC = 1, GF%N_FINE
      DO IC = 1, GF%NC
         DSUM = 0.0_EB
#ifdef WITH_SCARC_DEBUG2
WRITE(MSG%LU_DEBUG,'(A,30I6)') 'PROL: IC, P%COLS:', IC, (P%COL(ICOL), ICOL=P%ROW(IC), P%ROW(IC+1)-1)
#endif
         DO ICOL = P%ROW(IC), P%ROW(IC+1)-1                            
#ifdef WITH_SCARC_DEBUG2
WRITE(MSG%LU_DEBUG,*) 'PROL: IC, ICOL, P%COL(ICOL):', IC, ICOL, P%COL(ICOL)
#endif
            DSUM = DSUM + P%VAL(ICOL) * VC(P%COL(ICOL))
         ENDDO
         VF(IC) = DSUM
      ENDDO

   ENDDO

   IF (NMESHES > 1) CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_VECTOR_PLAIN, NVB, NLF)

ENDIF
#ifdef WITH_SCARC_DEBUG
CALL SCARC_DEBUG_LEVEL (NVC, 'PROLONGATION EXIT: NVC', NLC)
#endif

END SUBROUTINE SCARC_PROLONGATION


! --------------------------------------------------------------------------------------------------------
!> \brief Copy final solution from GMG (as preconditioner) to corresponding vector of CG (as main solver)
! --------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_UPDATE_PRECONDITIONER(NL)
INTEGER, INTENT(IN) :: NL
INTEGER :: NM
DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   SCARC(NM)%LEVEL(NL)%STAGE(NSCARC_STAGE_ONE)%V = SCARC(NM)%LEVEL(NL)%STAGE(NSCARC_STAGE_TWO)%X
ENDDO
END SUBROUTINE SCARC_UPDATE_PRECONDITIONER


! --------------------------------------------------------------------------------------------------------
!> \brief Finalize data for pressure vector (predictor/corrector) when local ScaRC solver has finished
! --------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_UPDATE_MAINCELLS(NL)
USE SCARC_POINTERS, ONLY: M, G, L, ST, HP
USE SCARC_POINTER_ROUTINES, ONLY: SCARC_POINT_TO_GRID
INTEGER, INTENT(IN) :: NL
INTEGER :: NM, IC 
#ifdef WITH_SCARC_DEBUG
INTEGER :: I, K
#endif

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G
   ST  => L%STAGE(NSCARC_STAGE_ONE)

   IF (PREDICTOR) THEN
      HP => M%H
   ELSE
      HP => M%HS
   ENDIF

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'UPDATE_MAIN_CELLS:1: HP'
WRITE(MSG%LU_DEBUG,MSG%CFORM3) ((HP(I,1,K), I=0, L%NX+1), K=L%NZ+1,0,-1)
#endif

   HP = 0.0_EB
   !!$OMP PARALLEL DO PRIVATE(IC) SCHEDULE(STATIC)
   DO IC = 1, G%NC
      HP (G%ICX(IC), G%ICY(IC), G%ICZ(IC)) = ST%X(IC)
#ifdef WITH_SCARC_DEBUG2
      WRITE(MSG%LU_DEBUG,'(A, 4I6, E14.6)') 'UPDATE_MAIN_CELLS: IC, IX, IY, IZ, HP(IC):', &
                                             IC, G%ICX(IC), G%ICY(IC), G%ICZ(IC), ST%X(IC)
#endif
   ENDDO
   !!$OMP END PARALLEL DO 

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'UPDATE_MAIN_CELLS:2: HP'
WRITE(MSG%LU_DEBUG,MSG%CFORM3) ((HP(I,1,K), I=0,L%NX+1), K=L%NZ+1,0,-1)
#endif

#ifdef WITH_SCARC_DEBUG
   CALL SCARC_DEBUG_VECTOR3_BIG (HP, NM, 'HP: UPDATE_MAIN_CELLS')
#endif
#ifdef WITH_SCARC_VERBOSE2
   CALL SCARC_VERBOSE_PRESSURE (HP, NM, 'main')
#endif
ENDDO

END SUBROUTINE SCARC_UPDATE_MAINCELLS


! ------------------------------------------------------------------------------------------------
!> \brief Set correct boundary values at external and internal boundaries
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_UPDATE_GHOSTCELLS(NL)
USE SCARC_POINTERS, ONLY: M, L, G, GWC, HP
USE SCARC_POINTER_ROUTINES, ONLY: SCARC_POINT_TO_GRID
INTEGER, INTENT(IN) :: NL
INTEGER :: NM, IW, IOR0, IXG, IYG, IZG, IXW, IYW, IZW 
#ifdef WITH_SCARC_DEBUG
INTEGER :: I, K
#endif

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G

   IF (PREDICTOR) THEN
      HP => M%H
   ELSE
      HP => M%HS
   ENDIF

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'UPDATE_GHOST_CELLS:1: HP'
WRITE(MSG%LU_DEBUG,MSG%CFORM3) ((HP(I,1,K), I=0, L%NX+1), K=L%NZ+1,0,-1)
#endif
   ! Compute ghost cell values
 
   !!$OMP PARALLEL DO SHARED(HP, M, L, G) PRIVATE(IW, IXG, IYG, IZG, IXW, IYW, IZW, IOR0, GWC) SCHEDULE(STATIC)
   WALL_CELLS_LOOP: DO IW = 1, L%N_WALL_CELLS_EXT

      GWC => G%WALL(IW)

      IXG = GWC%IXG
      IYG = GWC%IYG
      IZG = GWC%IZG

      IXW = GWC%IXW
      IYW = GWC%IYW
      IZW = GWC%IZW

      IOR0 = GWC%IOR

      SELECT CASE (IOR0)
         CASE ( 1)
            IF (GWC%BTYPE==DIRICHLET) THEN
               HP(IXG,IYW,IZW) = -HP(IXW,IYW,IZW) + 2.0_EB * M%BXS(IYW,IZW)
            ELSE IF (GWC%BTYPE==NEUMANN) THEN
               HP(IXG,IYW,IZW) =  HP(IXW,IYW,IZW) - L%DX *M%BXS(IYW,IZW)
            ENDIF
         CASE (-1)
            IF (GWC%BTYPE==DIRICHLET) THEN
               HP(IXG,IYW,IZW) = -HP(IXW,IYW,IZW) + 2.0_EB * M%BXF(IYW,IZW)
            ELSE IF (GWC%BTYPE==NEUMANN) THEN
               HP(IXG,IYW,IZW) =  HP(IXW,IYW,IZW) + L%DX *M%BXF(IYW,IZW)
            ENDIF
         CASE ( 2)
            IF (GWC%BTYPE==DIRICHLET) THEN
               HP(IXW,IYG,IZW) = -HP(IXW,IYW,IZW) + 2.0_EB * M%BYS(IXW,IZW)
            ELSE IF (GWC%BTYPE==NEUMANN) THEN
               HP(IXW,IYG,IZW) =  HP(IXW,IYW,IZW) - L%DY *M%BYS(IXW,IZW)
            ENDIF
         CASE (-2)
            IF (GWC%BTYPE==DIRICHLET) THEN
               HP(IXW,IYG,IZW) = -HP(IXW,IYW,IZW) + 2.0_EB * M%BYF(IXW,IZW)
            ELSE IF (GWC%BTYPE==NEUMANN) THEN
               HP(IXW,IYG,IZW) =  HP(IXW,IYW,IZW) + L%DY *M%BYF(IXW,IZW)
            ENDIF
         CASE ( 3)
            IF (GWC%BTYPE==DIRICHLET) THEN
               HP(IXW,IYW,IZG) = -HP(IXW,IYW,IZW) + 2.0_EB * M%BZS(IXW,IYW)
            ELSE IF (GWC%BTYPE==NEUMANN) THEN
               HP(IXW,IYW,IZG) =  HP(IXW,IYW,IZW) - L%DZ *M%BZS(IXW,IYW)
            ENDIF
         CASE (-3)
            IF (GWC%BTYPE==DIRICHLET) THEN
               HP(IXW,IYW,IZG) = -HP(IXW,IYW,IZW) + 2.0_EB * M%BZF(IXW,IYW)
            ELSE IF (GWC%BTYPE==NEUMANN) THEN
               HP(IXW,IYW,IZG) =  HP(IXW,IYW,IZW) + L%DZ *M%BZF(IXW,IYW)
            ENDIF
      END SELECT
#ifdef WITH_SCARC_DEBUG2
      WRITE(MSG%LU_DEBUG,'(A, 5I6, E14.6)') 'UPDATE_GHOST_CELLS: IW, IOR0, IXW, IYW, IZG, HP:',&
                                             IW, IOR0, IXW, IYW, IZG, HP(IXW, IYW, IZG)
#endif

   ENDDO WALL_CELLS_LOOP
   !!$OMP END PARALLEL DO

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'UPDATE_GHOST_CELLS:2: HP'
WRITE(MSG%LU_DEBUG,MSG%CFORM3) ((HP(I,1,K), I=0, L%NX+1), K=L%NZ+1,0,-1)
#endif
#ifdef WITH_SCARC_DEBUG
   CALL SCARC_DEBUG_VECTOR3_BIG (HP, NM, 'HP: UPDATE_GHOST_CELLS')
#endif
#ifdef WITH_SCARC_VERBOSE2
   CALL SCARC_VERBOSE_PRESSURE (HP, NM, 'h')
#endif

ENDDO

! Perform data exchange to achieve consistency of ghost values along internal boundaries
! Note: this is most probably no longer necessary because MESH_EXCHANGE(5) is used after the call of ScaRC

CALL SCARC_EXCHANGE(NSCARC_EXCHANGE_PRESSURE, NSCARC_NONE, NL)

#ifdef WITH_SCARC_DEBUG
DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G
   IF (PREDICTOR) THEN
         HP => M%H
   ELSE
      HP => M%HS
   ENDIF
   CALL SCARC_DEBUG_VECTOR3_BIG (HP, NM, 'HP: UPDATE_GHOST_CELLS - AFTER EXCHANGE')
ENDDO
#endif
   
END SUBROUTINE SCARC_UPDATE_GHOSTCELLS
   
   

! ------------------------------------------------------------------------------------------------
!> \brief Determine if cell should be considered during packing of zone numbers
! ------------------------------------------------------------------------------------------------
LOGICAL FUNCTION SCARC_FORBIDDEN_ZONE(SEND_BUFFER_INT, IZ, ICG1, ICG2)
INTEGER, DIMENSION(:), INTENT(IN) :: SEND_BUFFER_INT
INTEGER, INTENT(IN) :: IZ, ICG1, ICG2
INTEGER :: LL, ICG

SCARC_FORBIDDEN_ZONE = .FALSE.
LL = 5
DO ICG = ICG1, ICG2
   IF (SEND_BUFFER_INT(LL) == IZ) THEN
      SCARC_FORBIDDEN_ZONE = .TRUE.
      RETURN
   ENDIF
   LL = LL + 4
ENDDO
END FUNCTION SCARC_FORBIDDEN_ZONE


! ------------------------------------------------------------------------------------------------
!> \brief Check if difference of two values is less than a given tolerance
! ------------------------------------------------------------------------------------------------
LOGICAL FUNCTION MATCH (VAL1, VAL2)
REAL (EB), INTENT(IN) :: VAL1, VAL2
REAL (EB) :: TOL
TOL = 1.0E-10_EB
MATCH = .FALSE.
IF (Abs(VAL1-VAL2) <= TOL) MATCH = .TRUE.
RETURN
END FUNCTION MATCH


! ----------------------------------------------------------------------------------------------------
!> \brief Filter out mean value
! ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_FILTER_MEANVALUE(NV, NL)
USE SCARC_POINTERS, ONLY: L, G, VC
USE SCARC_POINTER_ROUTINES, ONLY: SCARC_POINT_TO_GRID, SCARC_POINT_TO_VECTOR
INTEGER, INTENT(IN) :: NV, NL
INTEGER :: NM, IC, I, J, K

MESH_REAL = 0.0_EB
DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G
   VC => SCARC_POINT_TO_VECTOR(NM, NL, NV)
   DO IC = 1, G%NC
      MESH_REAL(NM) = MESH_REAL(NM) + VC(IC)
   ENDDO
ENDDO

IF (N_MPI_PROCESSES > 1) &
   CALL MPI_ALLGATHERV(MPI_IN_PLACE,1,MPI_INTEGER,MESH_REAL,COUNTS,DISPLS,&
                       MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,IERROR)

GLOBAL_REAL = SUM(MESH_REAL(1:NMESHES))/REAL(NC_GLOBAL(NL))

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G
   VC  => SCARC_POINT_TO_VECTOR(NM, NL, NV)
   DO K = 1, L%NZ
      DO J = 1, L%NY
         DO I = 1, L%NX
            IF (IS_UNSTRUCTURED .AND. L%IS_SOLID(I, J, K)) CYCLE
            IC = G%CELL_NUMBER(I,J,K)
            VC(IC) = VC(IC) - GLOBAL_REAL
         ENDDO
      ENDDO
   ENDDO
ENDDO

END SUBROUTINE SCARC_FILTER_MEANVALUE


! ------------------------------------------------------------------------------------------------
!> \brief Restore last cell of last mesh
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_RESTORE_LAST_CELL (XX, NL)
USE SCARC_POINTERS, ONLY: S, VC
USE SCARC_POINTER_ROUTINES, ONLY: SCARC_POINT_TO_VECTOR
INTEGER, INTENT(IN) :: XX, NL

IF (UPPER_MESH_INDEX /= NMESHES .OR. TYPE_RELAX == NSCARC_RELAX_FFT) RETURN
S => SCARC(UPPER_MESH_INDEX)

VC => SCARC_POINT_TO_VECTOR (UPPER_MESH_INDEX, NL, XX)
VC(S%NC) = S%RHS_END

END SUBROUTINE SCARC_RESTORE_LAST_CELL





! ====================================================================================================
! ====================================================================================================
! Bundle of routines for the Smoothed Algebraic Multigrid Method
! ====================================================================================================
! ====================================================================================================

! ----------------------------------------------------------------------------------------------------
!> \brief Setup algebraic multigrid structures
! Allocate needed workspace for hierarchy of system matrices, prolongation, restriction, etc.
! Note: all used pointers end with either 'F' or 'C' where:
!     'F' corresponds to fine   level NL
!     'C' corresponds to coarse level NL+1
! Determine mesh hierarchy based on smoothed aggregation
! Compute QR-decomposition of nullspace vector in order to determine tentative prolongator 
! Set nullspace for next level and perform Jacobi relaxation to get the final prolongator
! If the maximum allowed level is not yet reached, set dimensions for next coarser level, 
! define its nullspace and perform relaxation to define the respective Prolongation matrix
! Define Poisson matrix on coarser level by Galerkin approach: A_coarse = R * A_fine * P
! ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_ALGEBRAIC_MULTIGRID
INTEGER :: NL
LOGICAL :: FURTHER_COARSENING_REQUIRED

IF (.NOT.HAS_AMG_LEVELS) RETURN
FURTHER_COARSENING_REQUIRED = .TRUE.

NL = NLEVEL_MIN
!CALL  SCARC_PYTHON_MATRIX(NL, 'A')

COARSENING_LOOP: DO WHILE (FURTHER_COARSENING_REQUIRED)

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) '========================================================'
WRITE(MSG%LU_DEBUG,*) ' ALGEBRAIC MULTIGRID : LEVEL ', NL
WRITE(MSG%LU_DEBUG,*) '========================================================'
#endif

   ! Determine the aggregation order among the meshes
   CALL SCARC_SETUP_AGGREGATION_ORDER                    

   ! Extract matrix diagonal from Poisson matrix A, determine strength of connection matrix and store inverted matrix diagonal 
   CALL SCARC_EXTRACT_MATRIX_DIAGONAL(NL)           
   CALL SCARC_SETUP_CONNECTION(NL)                  
   CALL SCARC_INVERT_MATRIX_DIAGONAL(NL)            

   ! Apply smoothed aggregation heuristic to specify aggregation zones and improve near null space by Jacobi relaxation step
   CALL SCARC_SETUP_AGGREGATION_ZONES(NL)           
   CALL SCARC_RELAX_NULLSPACE(NL)                   

   ! Setup final aggregation Zones matrix Z and Prolongation matrix P based on QR-decomposition
   CALL SCARC_SETUP_ZONE_OPERATOR(NL)               

      ! Setup restriction and prolongation matrices for GMG-like coarsening
   IF (TYPE_COARSENING == NSCARC_COARSENING_GMG) THEN
      CALL SCARC_SETUP_TRANSFER_GMG(NL)

   ! Setup Prolongation matrix P based on QR-decomposition, near nullspace on coarser level and corresponding Restriction matrix R
   ELSE
      CALL SCARC_SETUP_PROLONGATION_AMG(NL)
      CALL SCARC_SETUP_NULLSPACE_COARSE(NL)
      CALL SCARC_SETUP_RESTRICTION(NL)
   ENDIF

   ! First setup A*P matrix to finally build the Galerkin matrix R*A*P
   CALL SCARC_SETUP_POISSON_PROL(NL)             
   CALL SCARC_SETUP_GALERKIN(NL)             

   ! Remove workspace which is no longer used and get the next coarsening round on the wa
   CALL SCARC_CLEAN_WORKSPACE_AMG(NL)
      
   NL = NL + 1
   IF (NL == NLEVEL_MAX) FURTHER_COARSENING_REQUIRED = .FALSE.

ENDDO COARSENING_LOOP

END SUBROUTINE SCARC_SETUP_ALGEBRAIC_MULTIGRID


! ------------------------------------------------------------------------------------------------------
!> \brief  Setup order in which aggregation is performed over mesh decomposition
! ------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_AGGREGATION_ORDER
USE SCARC_POINTERS, ONLY: SUB, S
LOGICAL, ALLOCATABLE, DIMENSION(:) :: NOT_AGGREGATED
INTEGER :: NM, NOM, INBR, ICYCLE

CROUTINE = 'SCARC_SETUP_AGGREGATION_ORDER'
SUB => SUBDIVISION

CALL SCARC_ALLOCATE_INT2(SUB%ORDER, 1, NMESHES, 1, NMESHES, NSCARC_INIT_NONE, 'SUB%ORDER', CROUTINE)
CALL SCARC_ALLOCATE_LOG1(NOT_AGGREGATED, 1, NMESHES, NSCARC_INIT_TRUE, 'NOT_AGGREGATED', CROUTINE)

ICYCLE = 1
DO WHILE (ANY(NOT_AGGREGATED)) 
   SUB%ORDER(1:NMESHES, ICYCLE) = NSCARC_ORDER_UNASSIGNED
   DO NM = 1, NMESHES
      S => SCARC(NM)
      IF (NOT_AGGREGATED(NM) .AND. SUB%ORDER(NM, ICYCLE) /= NSCARC_ORDER_LOCKED) THEN
         SUB%ORDER(NM, ICYCLE) = NSCARC_ORDER_ACTIVE
         DO INBR = 1, SUB%N_NEIGHBORS(NM)
            NOM = SUB%NEIGHBORS(INBR, NM)
            SUB%ORDER(NOM, ICYCLE) = NSCARC_ORDER_LOCKED
         ENDDO
         NOT_AGGREGATED(NM) = .FALSE.
      ENDIF
   ENDDO
   ICYCLE = ICYCLE + 1
ENDDO

SUB%N_CYCLES = ICYCLE - 1

CALL SCARC_DEALLOCATE_LOG1 (NOT_AGGREGATED, 'NOT_AGGREGATED', CROUTINE)

END SUBROUTINE SCARC_SETUP_AGGREGATION_ORDER


! ------------------------------------------------------------------------------------------------------
!> \brief  Invert matrix diagonal which is already stored in DIAG-vector (reuse workspace)
! Scale each matrix element with inverse of diagonal and approximate spectral radius (currently disabled) 
! ------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_INVERT_MATRIX_DIAGONAL(NL)
USE SCARC_POINTERS, ONLY: G
USE SCARC_POINTER_ROUTINES, ONLY: SCARC_POINT_TO_GRID
INTEGER, INTENT(IN) :: NL
INTEGER :: NM, IC

MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G
   DO IC = 1, G%NCE
      G%DIAG(IC) = 1.0_EB/G%DIAG(IC)
   ENDDO
ENDDO MESHES_LOOP

END SUBROUTINE SCARC_INVERT_MATRIX_DIAGONAL


! ------------------------------------------------------------------------------------------------------
!> \brief  Compute a strength of connection matrix based on symmetric smoothed aggregation heuristic. 
! A nonzero connection A[i,j] is considered strong if:
!
!     abs(A[i,j]) >= theta * sqrt( abs(A[i,i]) * abs(A[j,j]) )
!
! The strength matrix S corresponds to the set of nonzero entries of A that are strong connections
! based on a strength of connection tolerance theta
! ------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_CONNECTION(NL)
USE SCARC_POINTERS, ONLY: G, A, S, C, OG, OA, OC
USE SCARC_POINTER_ROUTINES, ONLY: SCARC_POINT_TO_GRID, SCARC_POINT_TO_OTHER_GRID, &
                                  SCARC_POINT_TO_CMATRIX, SCARC_POINT_TO_OTHER_CMATRIX
INTEGER, INTENT(IN) :: NL
REAL(EB):: VAL, EPS, SCAL, CVAL_MAX, THETA
INTEGER :: NM, NOM, IC, JC, ICOL, IZONE, INBR

IF (TYPE_COARSENING == NSCARC_COARSENING_CUBIC) RETURN

CROUTINE = 'SCARC_SETUP_CONNECTION'

MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   
   CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G

   A => SCARC_POINT_TO_CMATRIX (G, NSCARC_MATRIX_POISSON)
   C => SCARC_POINT_TO_CMATRIX (G, NSCARC_MATRIX_CONNECTION)

#ifdef WITH_SCARC_DEBUG
   CALL SCARC_DEBUG_CMATRIX(A, 'A','POISSON')
   WRITE(MSG%LU_DEBUG,*) 'CONNECTION: DIAG:', SIZE(G%DIAG)
   WRITE(MSG%LU_DEBUG,'(8E14.6)') G%DIAG
#endif

   ! Allocate workspace for strength of connection matrix (use same size as Poisson matrix)
   C%N_VAL = A%N_VAL                         
   C%N_ROW = A%N_ROW
   CALL SCARC_ALLOCATE_CMATRIX(C, NL, NSCARC_PRECISION_DOUBLE, NSCARC_MATRIX_LIGHT, 'G%CONNECTION', CROUTINE)

   IF (NL == NLEVEL_MIN) THEN
      THETA = 0.10E+0_EB
   ELSE
      THETA = SCARC_MULTIGRID_THETA
   ENDIF
   THETA = SCARC_MULTIGRID_THETA
   
   DO INBR = 1, SCARC(NM)%N_NEIGHBORS

      NOM = S%NEIGHBORS(INBR)
      CALL SCARC_POINT_TO_OTHER_GRID(NM, NOM, NL)

      OA => SCARC_POINT_TO_OTHER_CMATRIX(OG, NSCARC_MATRIX_POISSON)
      OC => SCARC_POINT_TO_OTHER_CMATRIX(OG, NSCARC_MATRIX_CONNECTION)

      OC%N_VAL = 2*OA%N_VAL                   ! use double layers
      OC%N_ROW = OA%N_ROW           
      CALL SCARC_ALLOCATE_CMATRIX(OC, NL, NSCARC_PRECISION_DOUBLE, NSCARC_MATRIX_LIGHT, 'OG%CONNECTION', CROUTINE)
      
   ENDDO 

   ! Check strength-of-connection criterion  
   IZONE = 1
   C%ROW(1) = 1
   DO IC = 1, G%NC
   
      EPS = THETA**2 * ABS(G%DIAG(IC))                  ! EPS = theta**2 * A_ii
      DO ICOL = A%ROW(IC), A%ROW(IC+1) - 1
         JC  = A%COL(ICOL)
         IF (JC == 0) CYCLE                                             ! omit second layer
         VAL = A%VAL(ICOL)                                              ! VAL = A_ij
   
         ! Always add the diagonal: |A_ii|  >= THETA * sqrt(|A_ii| * |A_ii|)     true!
         IF (IC == JC) THEN
            C%COL(IZONE) = JC
            C%VAL(IZONE) = VAL
            IZONE = IZONE + 1

         ! Check subdiagonal entry: |A_ij|  >= THETA * sqrt(|A_ii| * |A_jj|)     ??
         ELSE IF (VAL**2 >= EPS * ABS(G%DIAG(JC))) THEN
            C%COL(IZONE) = JC
            C%VAL(IZONE) = VAL
            IZONE = IZONE + 1
#ifdef WITH_SCARC_VERBOSE
         ELSE
            WRITE(MSG%LU_VERBOSE,'(I6,A,4I6,6E14.6)') MYID+1,': CONNECTION: NO NEIGHBORS ', &
                                                      IC, ICOL, JC, IZONE, THETA, EPS, &
                                                      G%DIAG(IC), G%DIAG(JC), EPS*ABS(G%DIAG(JC)), VAL**2
#endif
         ENDIF

      ENDDO
      C%ROW(IC+1) = IZONE
   ENDDO
   
#ifdef WITH_SCARC_DEBUG
   CALL SCARC_DEBUG_CMATRIX(C, 'CONNECTION', 'STRENGTH OF CNNECTION - PASS 1')
#endif
   
   DO IC = 1, G%NC
      CVAL_MAX = 0.0_EB
      DO ICOL = C%ROW(IC), C%ROW(IC+1) - 1
         CVAL_MAX = MAX(ABS(C%VAL(ICOL)), CVAL_MAX)
      ENDDO
      SCAL = 1.0_EB/CVAL_MAX
      DO ICOL = C%ROW(IC), C%ROW(IC+1) - 1
         C%VAL(ICOL) = ABS(C%VAL(ICOL))*SCAL
      ENDDO
   ENDDO
   
#ifdef WITH_SCARC_DEBUG
   CALL SCARC_DEBUG_CMATRIX(C, 'C','CONNECTION')
#endif
   
ENDDO MESHES_LOOP
   
! If there are multiple meshes, exchange strength matrix on overlapping parts
IF (NMESHES > 1) THEN
   CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_MATRIX_COLS,  NSCARC_MATRIX_CONNECTION, NL)
   CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_MATRIX_VALS,  NSCARC_MATRIX_CONNECTION, NL)
ENDIF

END SUBROUTINE SCARC_SETUP_CONNECTION
 

! ------------------------------------------------------------------------------------------------------
!> \brief Setup aggregation zones for Smoothed Aggregation Algebraic Multigrid Method
! ------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_AGGREGATION_ZONES(NL)
USE SCARC_POINTERS, ONLY: SUB, C, CF, G, LF, LC, GF, GC
USE SCARC_POINTER_ROUTINES, ONLY: SCARC_POINT_TO_GRID, SCARC_POINT_TO_CMATRIX, &
                                  SCARC_POINT_TO_MULTIGRID
INTEGER, INTENT(IN) :: NL
INTEGER :: NM, NM2, ICYCLE, IC, IZL

CROUTINE = 'SCARC_SETUP_AGGREGATION_ZONES'

SUB => SUBDIVISION
MESH_INT = -1

! Allocate workspaces for coarse points, global and local aggregation zones 
MESHES_ALLOCATION_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G

   CALL SCARC_ALLOCATE_INT1 (G%ZONE_CENTERS,  1, G%NCE,  NSCARC_INIT_ZERO, 'G%ZONE_CENTERS', CROUTINE)
   CALL SCARC_ALLOCATE_INT1 (G%ZONES_GLOBAL, 1, G%NCE2, NSCARC_INIT_ZERO, 'G%ZONES_GLOBAL', CROUTINE)
   CALL SCARC_ALLOCATE_INT1 (G%ZONES_LOCAL,  1, G%NCE2, NSCARC_INIT_ZERO, 'G%ZONES_LOCAL', CROUTINE)

ENDDO MESHES_ALLOCATION_LOOP


COARSENING_TYPE_SELECT: SELECT CASE (TYPE_COARSENING)

 
   ! ---- Default aggregation procedure for SAMG
 
   CASE (NSCARC_COARSENING_AGGREGATED)

   WRITE(*,*) 'COARSENING_AGGREGATED'
      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
         CALL SCARC_POINT_TO_GRID(NM, NL)
         C => SCARC_POINT_TO_CMATRIX (G, NSCARC_MATRIX_CONNECTION)
         CALL SCARC_SETUP_COARSENING_AGGREGATION(G, C)
         MESH_INT (NM) = G%N_ZONES
      ENDDO
      
      ! Exchange overlapping information of active meshes

      IF (NMESHES > 1)  CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_ZONE_TYPES, NSCARC_NONE, NL)
   
   ! ---- Staggered aggregation procedure for SAMG
 
   CASE (NSCARC_COARSENING_AGGREGATED_S)

      CYCLES_LOOP1: DO ICYCLE = 1, SUB%N_CYCLES
   
         ! First aggregate on active meshes

         DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   
            IF (SUB%ORDER(NM, ICYCLE) == NSCARC_ORDER_ACTIVE) THEN
               CALL SCARC_POINT_TO_GRID(NM, NL)
               C => SCARC_POINT_TO_CMATRIX (G, NSCARC_MATRIX_CONNECTION)
               CALL SCARC_SETUP_COARSENING_AGGREGATION(G, C)
               MESH_INT (NM) = G%N_ZONES
            ENDIF
   
         ENDDO
      
         ! Exchange overlapping information of active meshes

         IF (NMESHES > 1)  CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_ZONE_TYPES, NSCARC_NONE, NL)
   
         ! Then aggregate on passive meshes (taking into account overlapping aggregate information of active meshes)

         DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
      
            IF (SUB%ORDER (NM, ICYCLE) /= NSCARC_ORDER_ACTIVE) THEN
               CALL SCARC_POINT_TO_GRID(NM, NL)
               C => SCARC_POINT_TO_CMATRIX (G, NSCARC_MATRIX_CONNECTION)
               CALL SCARC_SETUP_COARSENING_AGGREGATION(G, C)
               MESH_INT (NM) = G%N_ZONES
            ENDIF
   
         ENDDO
   
         ! Exchange overlapping information of passive meshes

         IF (NMESHES > 1)  CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_ZONE_TYPES, NSCARC_NONE, NL)
   
      ENDDO CYCLES_LOOP1
   
 
   ! ---- GMG-like aggregation procedure 
   !      In case of even cell numbers this process corresponds to the usual GMG coarsening
   !      in case of uneven cell number in a coordinate direction, on patch with 3 cells is used, the rest with patches of 2
 
   CASE (NSCARC_COARSENING_CUBIC)

   WRITE(*,*) 'COARSENING_CUBIC'
      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
         CALL SCARC_POINT_TO_MULTIGRID(NM, NL, NL+1)
         CALL SCARC_SETUP_COARSENING_CUBIC(LF, LC, GF, GC)
         MESH_INT(NM) = GF%N_ZONES
      ENDDO
      
      ! Exchange overlapping information 

      IF (NMESHES > 1)  CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_ZONE_TYPES, NSCARC_NONE, NL)
   
   CASE (NSCARC_COARSENING_GMG)

   WRITE(*,*) 'COARSENING_GMG'
      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
         CALL SCARC_POINT_TO_MULTIGRID(NM, NL, NL+1)
         CALL SCARC_SETUP_COARSENING_GMG(LF, LC, GF, GC)
         MESH_INT(NM) = GF%N_ZONES
      ENDDO
      
END SELECT COARSENING_TYPE_SELECT


! Broadcast number of zones of all meshes

IF (N_MPI_PROCESSES>1) &
   CALL MPI_ALLGATHERV(MPI_IN_PLACE, 1, MPI_INTEGER, MESH_INT, COUNTS, DISPLS, MPI_INTEGER, MPI_COMM_WORLD, IERROR)
      

! Prepare grid dimensions of coarse grid level
 
DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_MULTIGRID(NM, NL, NL+1)         

   ! Setup grid dimensions on coarse level 
   GC%NC_LOCAL(1:NMESHES) = MESH_INT(1:NMESHES)
   GC%NC_GLOBAL = SUM(MESH_INT(1:NMESHES))
   GC%NC  = GC%NC_LOCAL(NM)
   GC%NCE = GC%NC_LOCAL(NM)
   IF (NMESHES > 1) THEN
      DO NM2 = 2, NMESHES
         GC%NC_OFFSET(NM2) = GC%NC_OFFSET(NM2-1) + GC%NC_LOCAL(NM2-1)
      ENDDO
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) '============== NM=',NM
WRITE(MSG%LU_DEBUG,*) 'GC%NC_LOCAL(1:NMESHES) =', GC%NC_LOCAL(1:NMESHES)
WRITE(MSG%LU_DEBUG,*) 'GC%NC_GLOBAL ', GC%NC_GLOBAL
WRITE(MSG%LU_DEBUG,*) 'GC%NC ', GC%NC
WRITE(MSG%LU_DEBUG,*) 'GC%NCE', GC%NCE
WRITE(MSG%LU_DEBUG,*) 'GC%NC_OFFSET(1:NMESHES) ', GC%NC_OFFSET(1:NMESHES)
WRITE(MSG%LU_DEBUG,*) 'GF%NCE', GF%NCE
#endif
   ENDIF                   

   ! Setup mapping from local zones to global zones

   CALL SCARC_ALLOCATE_INT1(GC%LOCAL_TO_GLOBAL, 1, GF%NCE, NSCARC_INIT_ZERO, 'G%LOCAL_TO_GLOBAL', CROUTINE)
   DO IZL = 1, GC%NC
      GC%LOCAL_TO_GLOBAL(IZL) = IZL + GC%NC_OFFSET(NM)
#ifdef WITH_SCARC_DEBUG2
WRITE(MSG%LU_DEBUG,*) 'NM=',NM,':GC%LOCAL_TO_GLOBAL(',IZL,')=',GC%LOCAL_TO_GLOBAL(IZL)
#endif
   ENDDO

   DO IC = 1, GF%NC
      GF%ZONES_GLOBAL(IC) = GF%ZONES_LOCAL(IC) + GC%NC_OFFSET(NM)
#ifdef WITH_SCARC_DEBUG2
WRITE(MSG%LU_DEBUG,*) 'NM=',NM,':GF%ZONES_GLOBAL(',IC,')=',GF%ZONES_GLOBAL(IC)
#endif
   ENDDO

ENDDO

! Exchange zones information between meshes

IF (NMESHES > 1)  THEN
   CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_CELL_NEIGHBORS, NSCARC_NONE, NL)
   CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_ZONE_NEIGHBORS, NSCARC_NONE, NL)
   CALL SCARC_EXTRACT_ZONE_OVERLAPS(NL)
   CALL SCARC_EXTRACT_ZONE_POINTERS(NL)
ENDIF

! Determine final grid dimensions on coarser level and reduce zone arrays to correct length

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_MULTIGRID(NM, NL, NL+1)         

   GC%NCE  = GF%N_ZONES
   GC%NCE2 = GF%N_ZONES

   CALL SCARC_REDUCE_INT1(GC%LOCAL_TO_GLOBAL, GC%NCE2, 'GC%LOCAL_TO_GLOBAL', CROUTINE)
   CALL SCARC_REDUCE_INT1(GC%ZONE_CENTERS, GC%NCE2, 'GC%ZONE_CENTERS', CROUTINE)

   GF%N_COARSE = GF%N_ZONES

   CALL SCARC_REDUCE_INT1(GF%ZONES_LOCAL,  GF%NCE2, 'GC%LOCAL_TO_GLOBAL', CROUTINE)
   CALL SCARC_REDUCE_INT1(GF%ZONES_GLOBAL, GF%NCE2, 'GC%LOCAL_TO_GLOBAL', CROUTINE)

   CF => SCARC_POINT_TO_CMATRIX (GF, NSCARC_MATRIX_CONNECTION)
   CALL SCARC_DEALLOCATE_CMATRIX(CF, 'STRENGTH OF CONNECTION', CROUTINE)

#ifdef WITH_SCARC_DEBUG
   WRITE(MSG%LU_DEBUG,*) 'NM=',NM,': ================== END OF SETUP AGGREGATION_ZONES '
   WRITE(MSG%LU_DEBUG,*) 'GC%NC =', GC%NC
   WRITE(MSG%LU_DEBUG,*) 'GC%NCE =', GC%NCE
   WRITE(MSG%LU_DEBUG,*) 'GC%NCE2=', GC%NCE2
   WRITE(MSG%LU_DEBUG,*) 'GC%NC_GLOBAL =', GC%NC_GLOBAL
   WRITE(MSG%LU_DEBUG,*) 'GC%NC_LOCAL(1:NMESHES) =', GC%NC_LOCAL(1:NMESHES)
   WRITE(MSG%LU_DEBUG,*) 'GC%NC_OFFSET(1:NMESHES) =', GC%NC_OFFSET(1:NMESHES)
   WRITE(MSG%LU_DEBUG,*) 'GF%ZONES_LOCAL  ='
   WRITE(MSG%LU_DEBUG,'(8I6)') GF%ZONES_LOCAL
   WRITE(MSG%LU_DEBUG,*) 'GF%ZONES_GLOBAL  ='
   WRITE(MSG%LU_DEBUG,'(8I6)') GF%ZONES_GLOBAL
   WRITE(MSG%LU_DEBUG,*) 'GC%LOCAL_TO_GLOBAL  ='
   WRITE(MSG%LU_DEBUG,'(8I6)') GC%LOCAL_TO_GLOBAL
#endif
#ifdef WITH_SCARC_VERBOSE
   CALL SCARC_VERBOSE_BLENDER_ZONES(NM, NL)
#endif

ENDDO

END SUBROUTINE SCARC_SETUP_AGGREGATION_ZONES


! ------------------------------------------------------------------------------------------------------
!> \brief  Standard aggregation prodecure based on strength of connection matrix
! ------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_COARSENING_AGGREGATION(G, C)
TYPE (SCARC_CMATRIX_TYPE), POINTER, INTENT(IN) :: C
TYPE (SCARC_GRID_TYPE), POINTER, INTENT(IN) :: G
INTEGER :: IC, ICOL, JC, IZONE, JZONE
LOGICAL :: HAS_NEIGHBORS, HAS_AGGREGATED_NEIGHBORS

CROUTINE = 'SCARC_SETUP_COARSENING_AGGREGATION'

! 
! Pass 1 of aggregation:  Setup aggregation zones on internal cells of active mesh
! 
G%ZONES_LOCAL = 0
G%ZONES_GLOBAL = 0
G%ZONE_CENTERS = 0

PASS1_LOOP: DO IC = 1, G%NC

   IF (G%ZONES_LOCAL(IC) /= 0) CYCLE                           ! has cell already been aggregated?

   HAS_NEIGHBORS = .FALSE.
   HAS_AGGREGATED_NEIGHBORS = .FALSE.

   DO ICOL = C%ROW(IC), C%ROW(IC+1)-1                          ! are all neighbors free (not already aggregated)?
      JC = C%COL(ICOL)
      IF (JC /= 0 .AND. IC /= JC .AND. JC <= G%NC) THEN        ! only consider internal cells here
         HAS_NEIGHBORS = .TRUE.
         IF (G%ZONES_LOCAL(JC) /= 0) THEN
            HAS_AGGREGATED_NEIGHBORS = .TRUE.
            EXIT
         ENDIF
      ENDIF
   ENDDO

   IF (.NOT. HAS_NEIGHBORS) THEN                               ! do not aggregate isolated cells
      G%ZONES_LOCAL(IC) = NSCARC_HUGE_INT
   ELSE IF (.NOT. HAS_AGGREGATED_NEIGHBORS) THEN               ! build aggregate of this cell and its neighbors
      G%N_ZONES = G%N_ZONES + 1
      G%ZONES_LOCAL(IC) = G%N_ZONES
      G%ZONE_CENTERS(G%N_ZONES) = IC                
      DO ICOL = C%ROW(IC), C%ROW(IC+1)-1 
         JC = C%COL(ICOL)
         IF (JC /= 0 .AND. JC <= G%NC) G%ZONES_LOCAL(C%COL(ICOL)) = G%N_ZONES
      ENDDO
   ENDIF

#ifdef WITH_SCARC_DEBUG2
CALL SCARC_DEBUG_ZONES(G, IC, 1, 'AFTER ACTIVE PASS1')
#endif

ENDDO PASS1_LOOP


! 
! Pass 2 of Aggregation:  Add unaggregated nodes to neighboring aggregate
! 
PASS2_LOOP: DO IC = 1, G%NC

   IF (G%ZONES_LOCAL(IC) /= 0) CYCLE            
   DO ICOL = C%ROW(IC), C%ROW(IC+1)-1
      JC = C%COL(ICOL)
      JZONE = G%ZONES_LOCAL(JC)
      IF (JZONE > 0) THEN
         IF (JC >= G%NC .OR. G%ZONE_CENTERS(JZONE)>0) THEN
            G%ZONES_LOCAL(IC) = -JZONE
            EXIT
         ENDIF
      ENDIF
   ENDDO
ENDDO PASS2_LOOP
!G%N_ZONES = G%N_ZONES - 1                         !TODO: check
      
#ifdef WITH_SCARC_DEBUG
CALL SCARC_DEBUG_ZONES(G, -1, 1, 'AFTER ACTIVE PASS2')
#endif


! 
! Pass 3 of Aggregation:  Process remaining nodes which have not been aggregated yet
! 
PASS3_LOOP: DO IC = 1, G%NC

   IZONE = G%ZONES_LOCAL(IC)

   ! cell IC has not been aggregated
   IF (IZONE /= 0) THEN
      IF (IZONE > 0) THEN
         G%ZONES_LOCAL(IC) = IZONE 
      ELSE IF (IZONE == NSCARC_HUGE_INT ) THEN
         G%ZONES_LOCAL(IC) = -1
      ELSE
         G%ZONES_LOCAL(IC) = -IZONE 
      ENDIF
      CYCLE PASS3_LOOP
   ENDIF

   G%ZONES_LOCAL(IC) = G%N_ZONES
   G%ZONE_CENTERS(G%N_ZONES) = IC

   DO ICOL = C%ROW(IC), C%ROW(IC+1)-1
      JC = C%COL(ICOL)
      IF (JC <= G%NC .AND. G%ZONES_LOCAL(JC) == 0) G%ZONES_LOCAL(JC) = G%N_ZONES
   ENDDO
   G%N_ZONES = G%N_ZONES + 1

ENDDO PASS3_LOOP

IF (MINVAL(G%ZONES_LOCAL) < 0) THEN
   WRITE(*,*) MYID+1, ':CAUTION: CELL ',MINLOC(G%ZONES_LOCAL),' HAS NOT BEEN AGGREGATED DURING AGGREGATION'
#ifdef WITH_SCARC_VERBOSE
   WRITE(*,*) MYID+1, ':CAUTION: CELL ',MINLOC(G%ZONES_LOCAL),' HAS NOT BEEN AGGREGATED DURING AGGREGATION'
   WRITE(MSG%LU_VERBOSE,*) 'G%ZONES_LOCAL:'
   WRITE(MSG%LU_VERBOSE,'(8I12)') G%ZONES_LOCAL(1:G%NCE)
#endif
   CALL MPI_FINALIZE(IERROR)
   STOP
ENDIF
      
#ifdef WITH_SCARC_DEBUG
CALL SCARC_DEBUG_ZONES(G, -1, 1, 'AFTER ACTIVE PASS3')
#endif

END SUBROUTINE SCARC_SETUP_COARSENING_AGGREGATION


! ------------------------------------------------------------------------------------------------------
!> \brief Selfdefined geometric motivated aggregation procedure using cubic zones
! ------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_COARSENING_CUBIC(LF, LC, GF, GC)
TYPE (SCARC_LEVEL_TYPE), POINTER, INTENT(IN) :: LF, LC
TYPE (SCARC_GRID_TYPE), POINTER, INTENT(IN) :: GF, GC
INTEGER :: NXM, NYM, NZM, NXD, NYD, NZD, NXI, NYI, NZI
INTEGER :: IX, IY, IZ, IXZ, IYZ, IZZ, IXP, IYP, IZP, IX0, IY0, IZ0, IC
INTEGER, DIMENSION(:), ALLOCATABLE :: OFFX, OFFY, OFFZ
LOGICAL :: BFIRST

CROUTINE = 'SCARC_SETUP_COARSENING_CUBIC'

NXM = MOD(LF%NX,2)
NXD = LF%NX/2

IF (TWO_D) THEN
   NYM = 0
   NYD = 1
ELSE
   NYM = MOD(LF%NY,2)
   NYD = LF%NY/2
ENDIF

NZM = MOD(LF%NZ,2)
NZD = LF%NZ/2

! Temporarily - to prevent failure of following algorithm

IF ((LF%NX < 4) .OR. (.NOT.TWO_D .AND. LF%NY < 4) .OR. (LF%NZ < 4)) THEN 
   WRITE(*,*) 'Grid dimensions too small für GMG-like aggregation'
   CALL MPI_FINALIZE(IERROR)
   STOP
ENDIF

CALL SCARC_ALLOCATE_INT1 (OFFX, 1, NXD, NSCARC_INIT_ZERO, 'OFFX', CROUTINE)
CALL SCARC_ALLOCATE_INT1 (OFFY, 1, NYD, NSCARC_INIT_ZERO, 'OFFY', CROUTINE)
CALL SCARC_ALLOCATE_INT1 (OFFZ, 1, NZD, NSCARC_INIT_ZERO, 'OFFZ', CROUTINE)

! If even number of cells in x-direction, use patch length of 2 as in default GMG
! else insert a patch of 3 after the first quarter of cells in x-direction

IF (NXM == 0) THEN
   OFFX = 2
ELSE
   NXI = MAX(NXD/4,1)
   DO IX = 1, NXI
      OFFX(IX) = 2
   ENDDO
   OFFX(NXI+1) = 3
   DO IX = NXI+2, NXD
      OFFX(IX) = 2
   ENDDO
ENDIF

! If even number of cells in y-direction, use patch length of 2 as in default GMG
! else insert a patch of 3 after the first third of cells in x-direction

IF (TWO_D) THEN
   OFFY = 0
ELSE
   IF (NYM == 0) THEN
      OFFY = 2
   ELSE
      NYI = MAX(NYD/3,1)
      DO IY = 1, NYI
         OFFY(IY) = 2
      ENDDO
      OFFY(NYI+1) = 3
      DO IY = NYI+2, NYD
         OFFY(IY) = 2
      ENDDO
   ENDIF
ENDIF

! If even number of cells in x-direction, use patch length of 2 as in default GMG
! else insert a patch of 3 after the first half of cells in x-direction
! the idea is to use different portions in the different coordinate direction to prevent local concentrations

IF (NZM == 0) THEN
   OFFZ = 2
ELSE
   NZI = MAX(NZD/2,1)
   DO IZ = 1, NZI
      OFFZ(IZ) = 2
   ENDDO
   OFFZ(NZI+1) = 3
   DO IZ = NZI+2, NZD
      OFFZ(IZ) = 2
   ENDDO
ENDIF

LC%NX = NXD
LC%NY = NYD
LC%NZ = NZD

CALL SCARC_ALLOCATE_INT3 (GC%CELL_NUMBER, 0, LC%NX+1, 0, LC%NY+1, 0, LC%NZ+1, NSCARC_INIT_UNDEF, 'GC%CELL_NUMBER', CROUTINE)
CALL SCARC_ALLOCATE_LOG3 (LC%IS_SOLID, 0, LC%NX+1, 0, LC%NY+1, 0, LC%NZ+1, NSCARC_INIT_TRUE, 'LC%IS_SOLID', CROUTINE)

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'OFFX:'
WRITE(MSG%LU_DEBUG,'(16I6)') OFFX
WRITE(MSG%LU_DEBUG,*) 'OFFY:'
WRITE(MSG%LU_DEBUG,'(16I6)') OFFY
WRITE(MSG%LU_DEBUG,*) 'OFFZ:'
WRITE(MSG%LU_DEBUG,'(16I6)') OFFZ
WRITE(MSG%LU_DEBUG,*) 'NXD=',NXD
WRITE(MSG%LU_DEBUG,*) 'NYD=',NYD
WRITE(MSG%LU_DEBUG,*) 'NZD=',NZD
#endif

GF%ZONES_LOCAL = 0
GF%ZONES_GLOBAL = 0
GF%ZONE_CENTERS = 0
DIMENSION_IF: IF (TWO_D) THEN

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'TWO_D'
#endif
   IZ0 = 1
   DO IZ = 1, NZD
      IX0 = 1
      DO IX = 1, NXD

         BFIRST = .TRUE.
         DO IZZ = 0, OFFZ(IZ)-1
            DO IXZ = 0, OFFX(IX)-1
               IXP = IX0 + IXZ
               IZP = IZ0 + IZZ
               IF (IS_UNSTRUCTURED .AND. LF%IS_SOLID(IXP, 1, IZP)) CYCLE
               IC = GF%CELL_NUMBER(IXP, 1, IZP)
               IF (BFIRST) THEN
                  GF%N_ZONES = GF%N_ZONES + 1
                  BFIRST = .FALSE. 
                  GF%ZONE_CENTERS(GF%N_ZONES) = IC
                  GC%CELL_NUMBER(IX, 1, IZ) = GF%N_ZONES
                  LC%IS_SOLID(IX, 1, IZ) = .FALSE.
               ENDIF
               GF%ZONES_LOCAL(IC) = GF%N_ZONES
            ENDDO
         ENDDO
         IX0 = IX0 + OFFX(IX)
      ENDDO
      IZ0 = IZ0 + OFFZ(IZ)
   ENDDO

ELSE

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'THREE_D'
#endif
   IZ0 = 1
   DO IZ = 1, NZD
      IY0 = 1
      DO IY = 1, NYD
         IX0 = 1
         DO IX = 1, NXD

            BFIRST = .TRUE.
            DO IZZ = 0, OFFZ(IZ)-1
               DO IYZ = 0, OFFY(IY)-1
                  DO IXZ = 0, OFFX(IX)-1
                     IXP = IX0 + IXZ
                     IYP = IY0 + IYZ
                     IZP = IZ0 + IZZ
                     IF (IS_UNSTRUCTURED .AND. LF%IS_SOLID(IXP, IYP, IZP)) CYCLE
                     IC = GF%CELL_NUMBER(IXP, IYP, IZP)
                     IF (BFIRST) THEN
                        GF%N_ZONES = GF%N_ZONES + 1
                        BFIRST = .FALSE. 
                        GF%ZONE_CENTERS(GF%N_ZONES) = IC
                        GC%CELL_NUMBER(IX, IY, IZ) = GF%N_ZONES
                        LC%IS_SOLID(IX, IY, IZ) = .FALSE.
                     ENDIF
                     IC = GF%CELL_NUMBER(IXP, IYP, IZP)
                     GF%ZONES_LOCAL(IC) = GF%N_ZONES
                  ENDDO
               ENDDO
            ENDDO
            IX0 = IX0 + OFFX(IX)
         ENDDO
         IY0 = IY0 + OFFY(IY)
      ENDDO
      IZ0 = IZ0 + OFFZ(IZ)
   ENDDO
ENDIF DIMENSION_IF

CALL SCARC_DEALLOCATE_INT1 (OFFX, 'OFFX', CROUTINE)
CALL SCARC_DEALLOCATE_INT1 (OFFY, 'OFFY', CROUTINE)
CALL SCARC_DEALLOCATE_INT1 (OFFZ, 'OFFZ', CROUTINE)

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'LC%NX=',LC%NX
WRITE(MSG%LU_DEBUG,*) 'LC%NY=',LC%NY
WRITE(MSG%LU_DEBUG,*) 'LC%NZ=',LC%NZ
WRITE(MSG%LU_DEBUG,*) 'LC%IS_SOLID:'
DO IZ = 1, LC%NZ
   WRITE(MSG%LU_DEBUG,*)
   DO IY = 1, LC%NY
      WRITE(MSG%LU_DEBUG,*) (LC%IS_SOLID(IX, IY, IZ), IX=1, LC%NX)
   ENDDO
ENDDO
WRITE(MSG%LU_DEBUG,*) 'GC%CELL_NUMBER:'
DO IZ = 1, LC%NZ
   WRITE(MSG%LU_DEBUG,*)
   DO IY = 1, LC%NY
      WRITE(MSG%LU_DEBUG,'(8I4)') (GC%CELL_NUMBER(IX, IY, IZ), IX=1, LC%NX)
   ENDDO
ENDDO
CALL SCARC_DEBUG_ZONES(GF, -1, 1, 'AFTER ACTIVE PASS3')
#endif

END SUBROUTINE SCARC_SETUP_COARSENING_CUBIC

! ------------------------------------------------------------------------------------------------------
!> \brief Selfdefined geometric motivated aggregation procedure using cubic zones
! ------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_COARSENING_GMG(LF, LC, GF, GC)
TYPE (SCARC_LEVEL_TYPE), POINTER, INTENT(IN) :: LF, LC
TYPE (SCARC_GRID_TYPE),  POINTER, INTENT(IN) :: GF, GC
INTEGER :: MODX, MODY, MODZ
INTEGER :: RELX, RELY, RELZ
INTEGER :: IX, IY, IZ, IXZ, IYZ, IZZ, IXP, IYP, IZP, IX0, IY0, IZ0, IC
INTEGER, DIMENSION(:), ALLOCATABLE :: OFFX, OFFY, OFFZ
LOGICAL :: BFIRST

CROUTINE = 'SCARC_SETUP_COARSENING_GMG'

MODX = MOD(LF%NX,2)
LC%NX = FLOOR(REAL(LF%NX/2),EB) + MODX

IF (TWO_D) THEN
   MODY = 0
   LC%NY = 1
ELSE
   MODY = MOD(LF%NY,2)
   LC%NY = FLOOR(REAL(LF%NY/2),EB) + MODY
ENDIF

MODZ = MOD(LF%NZ,2)
LC%NZ = FLOOR(REAL(LF%NZ/2),EB) + MODZ

CALL SCARC_ALLOCATE_INT1 (OFFX, 1, LC%NX, NSCARC_INIT_ZERO, 'OFFX', CROUTINE)
CALL SCARC_ALLOCATE_INT1 (OFFY, 1, LC%NY, NSCARC_INIT_ZERO, 'OFFY', CROUTINE)
CALL SCARC_ALLOCATE_INT1 (OFFZ, 1, LC%NZ, NSCARC_INIT_ZERO, 'OFFZ', CROUTINE)

OFFX = 2 ;  IF (MODX /= 0) OFFX(LC%NX) = 1
IF (TWO_D) THEN
   OFFY = 0 ;  LC%NY = 1
ELSE
   OFFY = 2 ;  IF (MODY /= 0) OFFY(LC%NY) = 1
ENDIF
OFFZ = 2 ;  IF (MODZ /= 0) OFFZ(LC%NZ) = 1

RELX = CEILING(REAL(LF%NX/LC%NX),EB)
RELY = CEILING(REAL(LF%NY/LC%NY),EB)
RELZ = CEILING(REAL(LF%NZ/LC%NZ),EB)

CALL SCARC_ALLOCATE_INT3 (GC%CELL_NUMBER, 0, LC%NX+1, 0, LC%NY+1, 0, LC%NZ+1, NSCARC_INIT_UNDEF, 'GC%CELL_NUMBER', CROUTINE)
CALL SCARC_ALLOCATE_LOG3 (LC%IS_SOLID,    0, LC%NX+1, 0, LC%NY+1, 0, LC%NZ+1, NSCARC_INIT_TRUE,  'LC%IS_SOLID', CROUTINE)

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'MODX, RELX, LC%NX=', MODX, RELX, LC%NX
WRITE(MSG%LU_DEBUG,*) 'MODY, RELY, LC%NY=', MODY, RELY, LC%NY
WRITE(MSG%LU_DEBUG,*) 'MODZ, RELZ, LC%NZ=', MODZ, RELZ, LC%NZ
#endif

GF%ZONES_LOCAL  = 0
GF%ZONES_GLOBAL = 0
GF%ZONE_CENTERS = 0

DIMENSION_IF: IF (TWO_D) THEN

   IZ0 = 1
   DO IZ = 1, LC%NZ
      IX0 = 1
      DO IX = 1, LC%NX

         BFIRST = .TRUE.
         DO IZZ = 0, OFFZ(IZ)-1
            DO IXZ = 0, OFFX(IX)-1
               IXP = IX0 + IXZ
               IZP = IZ0 + IZZ
               IF (IS_UNSTRUCTURED .AND. LF%IS_SOLID(IXP, 1, IZP)) CYCLE
               IC = GF%CELL_NUMBER(IXP, 1, IZP)
               IF (BFIRST) THEN
                  GF%N_ZONES = GF%N_ZONES + 1
                  BFIRST = .FALSE. 
                  GF%ZONE_CENTERS(GF%N_ZONES) = IC
                  GC%CELL_NUMBER(IX, 1, IZ) = GF%N_ZONES
                  LC%IS_SOLID(IX, 1, IZ) = .FALSE.
               ENDIF
               GF%ZONES_LOCAL(IC) = GF%N_ZONES
            ENDDO
         ENDDO
         IX0 = IX0 + OFFX(IX)
      ENDDO
      IZ0 = IZ0 + OFFZ(IZ)
   ENDDO

ELSE

   IZ0 = 1
   DO IZ = 1, LC%NZ
      IY0 = 1
      DO IY = 1, LC%NY
         IX0 = 1
         DO IX = 1, LC%NX

            BFIRST = .TRUE.
            DO IZZ = 0, OFFZ(IZ)-1
               DO IYZ = 0, OFFY(IY)-1
                  DO IXZ = 0, OFFX(IX)-1
                     IXP = IX0 + IXZ
                     IYP = IY0 + IYZ
                     IZP = IZ0 + IZZ
                     IF (IS_UNSTRUCTURED .AND. LF%IS_SOLID(IXP, IYP, IZP)) CYCLE
                     IC = GF%CELL_NUMBER(IXP, IYP, IZP)
                     IF (BFIRST) THEN
                        GF%N_ZONES = GF%N_ZONES + 1
                        BFIRST = .FALSE. 
                        GF%ZONE_CENTERS(GF%N_ZONES) = IC
                        GC%CELL_NUMBER(IX, IY, IZ) = GF%N_ZONES
                        LC%IS_SOLID(IX, IY, IZ) = .FALSE.
                     ENDIF
                     IC = GF%CELL_NUMBER(IXP, IYP, IZP)
                     GF%ZONES_LOCAL(IC) = GF%N_ZONES
                  ENDDO
               ENDDO
            ENDDO
            IX0 = IX0 + OFFX(IX)
         ENDDO
         IY0 = IY0 + OFFY(IY)
      ENDDO
      IZ0 = IZ0 + OFFZ(IZ)
   ENDDO

ENDIF DIMENSION_IF

CALL SCARC_DEALLOCATE_INT1 (OFFX, 'OFFX', CROUTINE)
CALL SCARC_DEALLOCATE_INT1 (OFFY, 'OFFY', CROUTINE)
CALL SCARC_DEALLOCATE_INT1 (OFFZ, 'OFFZ', CROUTINE)

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'LC%NX=',LC%NX
WRITE(MSG%LU_DEBUG,*) 'LC%NY=',LC%NY
WRITE(MSG%LU_DEBUG,*) 'LC%NZ=',LC%NZ
WRITE(MSG%LU_DEBUG,*) 'LC%IS_SOLID:'
DO IZ = 1, LC%NZ
   WRITE(MSG%LU_DEBUG,*)
   DO IY = 1, LC%NY
      WRITE(MSG%LU_DEBUG,*) (LC%IS_SOLID(IX, IY, IZ), IX=1, LC%NX)
   ENDDO
ENDDO
WRITE(MSG%LU_DEBUG,*) 'GC%CELL_NUMBER:'
DO IZ = 1, LC%NZ
   WRITE(MSG%LU_DEBUG,*)
   DO IY = 1, LC%NY
      WRITE(MSG%LU_DEBUG,'(8I4)') (GC%CELL_NUMBER(IX, IY, IZ), IX=1, LC%NX)
   ENDDO
ENDDO
CALL SCARC_DEBUG_ZONES(GF, -1, 1, 'AFTER ACTIVE PASS3')
#endif

END SUBROUTINE SCARC_SETUP_COARSENING_GMG

! ------------------------------------------------------------------------------------------------------
!> \brief Remove workspace on specified grid level which will no longer be needed after matrix setup
! ------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_CLEAN_WORKSPACE_SYSTEM(NL)
USE SCARC_POINTER_ROUTINES, ONLY: SCARC_POINT_TO_GRID
INTEGER, INTENT(IN) :: NL
INTEGER:: NM

! TODO: deallocate arrays which are no longer used
DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   CALL SCARC_POINT_TO_GRID (NM, NL)
ENDDO

END SUBROUTINE SCARC_CLEAN_WORKSPACE_SYSTEM


! ------------------------------------------------------------------------------------------------------
!> \brief Remove workspace on specified grid level which will no longer be needed in SAMG method
! ------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_CLEAN_WORKSPACE_AMG(NL)
USE SCARC_POINTERS, ONLY: G
USE SCARC_POINTER_ROUTINES, ONLY: SCARC_POINT_TO_GRID
INTEGER, INTENT(IN) :: NL
INTEGER:: NM

CROUTINE = 'SCARC_CLEAN_WORKSPACE_AMG'

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   CALL SCARC_POINT_TO_GRID (NM, NL)
   IF (ALLOCATED(G%ZONE_CENTERS)) CALL SCARC_DEALLOCATE_INT1 (G%ZONE_CENTERS, 'G%ZONE_CENTERS', CROUTINE)
   IF (ALLOCATED(G%AUX1)) CALL SCARC_DEALLOCATE_REAL1 (G%AUX1, 'G%AUX1', CROUTINE)
   IF (ALLOCATED(G%AUX2)) CALL SCARC_DEALLOCATE_REAL1 (G%AUX2, 'G%AUX2', CROUTINE)
   IF (ALLOCATED(G%RR)) CALL SCARC_DEALLOCATE_REAL1 (G%RR, 'G%RR', CROUTINE)
   IF (ALLOCATED(G%QQ)) CALL SCARC_DEALLOCATE_REAL1 (G%QQ, 'G%QQ', CROUTINE)
ENDDO

END SUBROUTINE SCARC_CLEAN_WORKSPACE_AMG


! ------------------------------------------------------------------------------------------------------
!> \brief Perform relaxation of nullspac
! Perform AMG Jacobi :.. x = x - omega D^{-1} (Ax-b)
! Near-null space vector is given in vector G%NULLSPACE --> corresponds to x
! vector b is zero
! ------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_RELAX_NULLSPACE(NL)
USE SCARC_POINTERS, ONLY: L, G, A, MG
USE SCARC_POINTER_ROUTINES, ONLY: SCARC_POINT_TO_GRID, SCARC_POINT_TO_CMATRIX
INTEGER, INTENT(IN) :: NL
INTEGER :: IC, ICOL, NM, JC, JCG

CROUTINE = 'SCARC_RELAX_NULLSPACE'

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   CALL SCARC_POINT_TO_GRID (NM, NL) 
   CALL SCARC_ALLOCATE_REAL1 (G%AUX1, 1, G%NCE2+1, NSCARC_INIT_ZERO, 'G%AUX1', CROUTINE)
   CALL SCARC_ALLOCATE_REAL1 (G%AUX2, 1, G%NCE2+1, NSCARC_INIT_ZERO, 'G%AUX2', CROUTINE)
ENDDO

! For coarser levels exchange numbers and values of second layer cells with are needed for nullspace computation

IF (SCARC_MULTIGRID_RELAXING .AND. NL > NLEVEL_MIN) THEN
   CALL SCARC_IDENTIFY_LAYER2(NL)
   CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_LAYER2_NUMS, NSCARC_NONE, NL)
   CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_LAYER2_VALS, NSCARC_NONE, NL)
ENDIF

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID (NM, NL)         
   A => SCARC_POINT_TO_CMATRIX (G, NSCARC_MATRIX_POISSON)

   MG => L%MG
   MG%OMEGA = 1.0_EB/MG%APPROX_SPECTRAL_RADIUS
   !MG%OMEGA = 4.0_EB/(3.0_EB*MG%APPROX_SPECTRAL_RADIUS)

#ifdef WITH_SCARC_DEBUG
IF (SCARC_MULTIGRID_RELAXING) THEN
WRITE(MSG%LU_DEBUG,*) 'Using A', G%NCE2, NL
CALL SCARC_DEBUG_CMATRIX(A, 'A','A IN NULLSPACE')
WRITE(MSG%LU_DEBUG,*) 'Using LOCAL_TO_GLOBAL:', G%NCE2, NL
WRITE(MSG%LU_DEBUG,'(8I6)') G%ZONES_LOCAL
WRITE(MSG%LU_DEBUG,*) 'Using ZONES_GLOBAL:', G%NCE2, NL
WRITE(MSG%LU_DEBUG,'(8I6)') G%ZONES_GLOBAL
IF (NL > NLEVEL_MIN .AND. TYPE_COARSENING == NSCARC_COARSENING_AGGREGATED) THEN
   WRITE(MSG%LU_DEBUG,*) 'Using LAYER2_NUMS:', G%NCE2, NL, L%N_LAYER2_TOTAL
   WRITE(MSG%LU_DEBUG,'(8I6)') G%ELAYER2_NUMS(1: L%N_LAYER2_TOTAL)
ENDIF
ENDIF
#endif

   ! Compute defect to near-null-space vector: d = Ax - b, in this case
   !    'x' corresponds to nullspace vector consisting of only '1'-entries 
   !    'b' corresponds to zero vector 

   ! On finest level NULLSPACE vector is preset with 1, so matrix entries can simply be added
   IF (NL == NLEVEL_MIN) THEN

      CALL SCARC_ALLOCATE_REAL1 (G%NULLSPACE, 1, G%NCE2, NSCARC_INIT_ONE, 'G%NULLSPACE', CROUTINE)
      FINE_CELLS_LOOP: DO IC = 1, G%NC
         G%AUX2(IC) = 0.0_EB
         DO ICOL = A%ROW(IC), A%ROW(IC+1)-1                          
            G%AUX2(IC) = G%AUX2(IC) + A%VAL(ICOL)
         ENDDO
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,'(A,1I6, 2E14.6)') &
   'FINE LEVEL: IC, A%VAL(ICOL), AUX2:', IC, A%VAL(ICOL), G%AUX2(IC)
#endif
      ENDDO FINE_CELLS_LOOP

   ! on coarser levels NULLSPACE vector was set in preceding coarsening loop to R-vector from QR-decomposition 
   ELSE

      G%AUX1 = G%NULLSPACE
      COARSE_CELLS_LOOP: DO IC = 1, G%NC
   
         G%AUX2(IC) = 0.0_EB
         DO ICOL = A%ROW(IC), A%ROW(IC+1)-1                          
            JC = A%COL(ICOL)
            IF (JC /= 0) THEN
               G%AUX2(IC) =  G%AUX2(IC) + A%VAL(ICOL) * G%AUX1(JC)
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,'(A,3I6, E14.6)') &
   'RELAXING NULLSPACE:A: IC, ICOL, JC, AUX2:', IC, ICOL, JC, G%AUX2(IC)
#endif
            ELSE 
               JCG = A%COLG(ICOL)
               JC = FINDLOC (G%LOCAL_TO_GLOBAL(1:G%NCE2), VALUE = JCG, DIM = 1)
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,'(A,3I6,A,I6)') &
   'RELAXING NULLSPACE:B: IC, ICOL, JCG     :', IC, ICOL, JCG, ' Searching, found ', JC
#endif
               IF (JC /= 0) THEN
                  G%AUX2(IC) = G%AUX2(IC) + A%VAL(ICOL) * G%AUX1(JC)
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,'(A,3I6, E14.6, I6, 2E14.6)') &
   'RELAXING NULLSPACE:C: IC, ICOL, JC, AUX2:', IC, ICOL, JC, G%AUX2(IC), JCG, A%VAL(ICOL), G%AUX1(JC)
#endif
               ELSE
#ifdef WITH_SCARC_VERBOSE
WRITE(MSG%LU_VERBOSE,*) 'RELAX_NULLSPACE: STENCIL FOR IC = ', IC,': GLOBAL LEG CELL ', JCG, ' NOT FOUND!'
#endif
                  JC = FINDLOC (G%ELAYER2_NUMS(1: L%N_LAYER2_TOTAL), VALUE = JCG, DIM = 1)
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'RELAX_NULLSPACE: STENCIL FOR IC = ', IC,': GLOBAL LEG CELL ', JCG, ' NOT FOUND, but found ', JC
#endif
                  IF (JC /= 0) G%AUX2(IC) = G%AUX2(IC) + A%VAL(ICOL) * G%ELAYER2_VALS(JC)
#ifdef WITH_SCARC_DEBUG
                  IF (JC /= 0) WRITE(MSG%LU_DEBUG,'(A,3I6, E14.6, I6, 2E14.6)') &
   'RELAXING NULLSPACE:D: IC, ICOL, JC, AUX2:', IC, ICOL, JC, G%AUX2(IC), JCG, A%VAL(ICOL), G%ELAYER2_VALS(JC)
#endif
               ENDIF
            ENDIF
         ENDDO
   
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) '======================= IC, G%AUX2(IC) =', IC, G%AUX2(IC)
#endif
      ENDDO COARSE_CELLS_LOOP

   ENDIF
   
   ! Scale it by parameter omega and inverse of diagonal:   d = omega D^{-1} d

   DO IC = 1, G%NC
      G%AUX2(IC) = MG%OMEGA * G%DIAG(IC) * G%AUX2(IC) 
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'RELAXING NULLSPACE: IC, DIAG, AUX2:', IC, G%DIAG(IC), G%AUX2(IC)
#endif
   ENDDO
   
  ! Get new iterate:   x = x - d

#ifdef WITH_MKL
  CALL DAXPBY(G%NC, -1.0_EB, G%AUX2, 1, 1.0_EB, G%NULLSPACE, 1)
#else
  CALL SCARC_DAXPY_CONSTANT_DOUBLE(G%NC, -1.0_EB, G%AUX2, 1.0_EB, G%NULLSPACE)
#endif
   
#ifdef WITH_SCARC_DEBUG
   WRITE(MSG%LU_DEBUG,*) 'NM=',NM,': ======================================================='
   WRITE(MSG%LU_DEBUG,*) 'OMEGA=',MG%OMEGA
   WRITE(MSG%LU_DEBUG,*) 'RELAX_NULLSPACE: AUX1: '
   WRITE(MSG%LU_DEBUG,'(6E14.6)') G%AUX1(1: G%NC)
   WRITE(MSG%LU_DEBUG,*) '---------------------------------'
   WRITE(MSG%LU_DEBUG,'(6E14.6)') G%AUX1(G%NC+1:G%NCE2)
   WRITE(MSG%LU_DEBUG,*) '======================================================='
   WRITE(MSG%LU_DEBUG,*) 'RELAX_NULLSPACE: AUX2: '
   WRITE(MSG%LU_DEBUG,'(6E14.6)') G%AUX2(1: G%NC)
   WRITE(MSG%LU_DEBUG,*) '---------------------------------'
   WRITE(MSG%LU_DEBUG,'(6E14.6)') G%AUX2(G%NC+1:G%NCE2)
   WRITE(MSG%LU_DEBUG,*) '======================================================='
   WRITE(MSG%LU_DEBUG,*) 'RELAX_NULLSPACE: NULLSPACE: '
   WRITE(MSG%LU_DEBUG,'(6E14.6)') G%NULLSPACE(1: G%NC)
   WRITE(MSG%LU_DEBUG,*) '---------------------------------'
   WRITE(MSG%LU_DEBUG,'(6E14.6)') G%NULLSPACE(G%NC+1:G%NCE2)
#endif

ENDDO

IF (NMESHES > 1) CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_NULLSPACE, NSCARC_NONE, NL)

END SUBROUTINE SCARC_RELAX_NULLSPACE


! ------------------------------------------------------------------------------------------------------
!> \brief Setup basic structure of Prolongation matrix
! This concerns the setting of the number of rows and the column pointers
! The values are still missing and are set in a later step
! ------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_ZONE_OPERATOR(NL)
USE SCARC_POINTERS, ONLY: GF, GC, AF, ZF
USE SCARC_POINTER_ROUTINES, ONLY: SCARC_POINT_TO_MULTIGRID, SCARC_POINT_TO_CMATRIX
INTEGER, INTENT(IN) :: NL
INTEGER :: NM, ICC, ICCL, ICCG, IC, IP, N_ROW, N_VAL
LOGICAL :: IS_INCLUDED 

CROUTINE = 'SCARC_SETUP_ZONE_OPERATOR'

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_MULTIGRID (NM, NL, NL+1)                       ! Sets grid pointer GF and GC

   AF => SCARC_POINT_TO_CMATRIX (GF, NSCARC_MATRIX_POISSON)
   ZF => SCARC_POINT_TO_CMATRIX (GF, NSCARC_MATRIX_ZONES)

   ! First use very conservative bounds for the size of the zones operator matrix 
   ! reduce it later once the real size is known

   ZF%N_VAL = AF%N_VAL                  
   ZF%N_ROW = AF%N_VAL                  
   CALL SCARC_ALLOCATE_CMATRIX(ZF, NL, NSCARC_PRECISION_DOUBLE, NSCARC_MATRIX_MINIMAL, 'GF%ZONES', CROUTINE)
  
   ! Again conservative upper bound for length - to be reduced later 

   CALL SCARC_ALLOCATE_INT1(GF%ZONES_LOCAL, 1, GF%NCE2, NSCARC_INIT_ZERO, 'GF%ZONES_LOCAL', CROUTINE)

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'NM=',NM,': =================== SETUP_ZONE_OPERATOR:'
WRITE(MSG%LU_DEBUG,*) 'ZF%ROW:', ZF%N_ROW, ZF%N_VAL, SIZE(ZF%ROW), SIZE(ZF%VAL)
WRITE(MSG%LU_DEBUG,*) 'ZONES_LOCAL:'
WRITE(MSG%LU_DEBUG,'(8I6)') GF%ZONES_LOCAL(1:GF%NCE2)
WRITE(MSG%LU_DEBUG,*) 'ZONES_GLOBAL:'
WRITE(MSG%LU_DEBUG,'(8I6)') GF%ZONES_GLOBAL(1:GF%NCE2)
WRITE(MSG%LU_DEBUG,*) '=================== SETUP_ZONE_OPERATOR: FINE'
WRITE(MSG%LU_DEBUG,*) 'SETUP_ZONE_OPERATOR: GF%N_FINE:', GF%N_FINE
WRITE(MSG%LU_DEBUG,*) 'SETUP_ZONE_OPERATOR: GF%N_COARSE:', GF%N_COARSE
WRITE(MSG%LU_DEBUG,*) 'SETUP_ZONE_OPERATOR: GF%NC:', GF%NC
WRITE(MSG%LU_DEBUG,*) 'SETUP_ZONE_OPERATOR: GF%NCE:', GF%NCE
WRITE(MSG%LU_DEBUG,*) 'SETUP_ZONE_OPERATOR: GF%NCE2:', GF%NCE2
WRITE(MSG%LU_DEBUG,*) '=================== SETUP_ZONE_OPERATOR: COARSE'
WRITE(MSG%LU_DEBUG,*) 'SETUP_ZONE_OPERATOR: GC%NC_LOCAL:', GC%NC_LOCAL
WRITE(MSG%LU_DEBUG,*) 'SETUP_ZONE_OPERATOR: GC%NC:', GC%NC
WRITE(MSG%LU_DEBUG,*) 'SETUP_ZONE_OPERATOR: GC%NCE:', GC%NCE
WRITE(MSG%LU_DEBUG,*) 'SETUP_ZONE_OPERATOR: GC%NCE2:', GC%NCE2
WRITE(MSG%LU_DEBUG,*) 'SETUP_ZONE_OPERATOR: GC%NC_OFFSET:', GC%NC_OFFSET
WRITE(MSG%LU_DEBUG,*) 'SETUP_ZONE_OPERATOR: GC%LOCAL_TO_GLOBAL:'
WRITE(MSG%LU_DEBUG,'(8I6)') GC%LOCAL_TO_GLOBAL(1:GC%NCE2)
#endif

   ! Based on global zone numbers determine local zone numbers within mesh

   IP = 1
   ICC = 1
   ZF%ROW(ICC) = 1
   DO ICCL = 1, GC%NCE2
      ICCG = GC%LOCAL_TO_GLOBAL(ICCL)
      IS_INCLUDED = .FALSE.
      DO IC = 1, GF%NCE2
         IF (GF%ZONES_LOCAL(IC) /= ICCL) CYCLE
         IS_INCLUDED = .TRUE.
         ZF%COL(IP)  = IC
         IP = IP + 1
      ENDDO
      IF (IS_INCLUDED) THEN
         ICC = ICC + 1
         ZF%ROW(ICC) = IP
      ENDIF
   ENDDO

   N_ROW = ICC
   N_VAL = IP - 1

   ZF%N_ROW=N_ROW
   ZF%N_VAL=N_VAL

   CALL SCARC_REDUCE_CMATRIX(ZF, 'ZF', CROUTINE)

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) '------------- NM=',NM
WRITE(MSG%LU_DEBUG,*) 'GF%NCE=',GF%NCE
WRITE(MSG%LU_DEBUG,*) 'GF%N_COARSE=',GF%N_COARSE
WRITE(MSG%LU_DEBUG,*) 'GC%NCE=',GC%NCE
WRITE(MSG%LU_DEBUG,*) 'ZF%N_ROW=',ZF%N_ROW
WRITE(MSG%LU_DEBUG,*) 'ZF%N_VAL=',ZF%N_VAL
CALL SCARC_DEBUG_CMATRIX (ZF, 'ZONES','AFTER SETUP AGGREGATION ZONES 2')
#endif

ENDDO

END SUBROUTINE SCARC_SETUP_ZONE_OPERATOR


! ------------------------------------------------------------------------------------------------------
!> \brief Determine tentative prolongator for current level by computing QR-decomposition of smoothed 
! nullspace vector and set nullspace for next level
! Compute the tentative prolongator, T, which is a tentative interpolation
! matrix from the coarse-grid to the fine-grid.  T exactly interpolates  B_fine = T B_coarse.
! ------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_PROLONGATION_AMG(NL)
USE SCARC_POINTERS, ONLY:  S, L, G, OG, A, OA, P, OP, GF, OGF, GC, PF, OPF, Z, MG
USE SCARC_POINTER_ROUTINES, ONLY: SCARC_POINT_TO_GRID, SCARC_POINT_TO_OTHER_GRID, SCARC_POINT_TO_MULTIGRID, &
                                  SCARC_POINT_TO_CMATRIX, SCARC_POINT_TO_OTHER_CMATRIX
INTEGER, INTENT(IN) :: NL
REAL(EB):: DSUM, SCAL !, TOL = 1.0E-12_EB
INTEGER :: NM, NOM, IC, JC, ICC, ICC0, ICOL, ICCOL, JCCOL, IP0, IP, JCC, IQ, INBR, NLEN

CROUTINE = 'SCARC_SETUP_PROLONGATION'

! Allocate several workspaces (with conservative bounds which will be reduced later)
!    - Prolongation matrix on internal part of mesh 
!    - for every neighbor small Prolongation matrix for corresponding overlap
!    - vectors Q and R for QR-decomposition of aggregation zones operator
! Initialize QR-decomposition
 
DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G

   A => SCARC_POINT_TO_CMATRIX(G, NSCARC_MATRIX_POISSON)
   P => SCARC_POINT_TO_CMATRIX(G, NSCARC_MATRIX_PROLONGATION)
   Z => SCARC_POINT_TO_CMATRIX(G, NSCARC_MATRIX_ZONES)

   P%N_VAL = A%N_VAL + 1        
   P%N_ROW = G%NCE + 1                  
   CALL SCARC_ALLOCATE_CMATRIX(P, NL, NSCARC_PRECISION_DOUBLE, NSCARC_MATRIX_FULL, 'G%PROLONGATION', CROUTINE)

   DO INBR = 1, S%N_NEIGHBORS

      NOM = S%NEIGHBORS(INBR)
      CALL SCARC_POINT_TO_OTHER_GRID(NM, NOM, NL)

      OA => SCARC_POINT_TO_OTHER_CMATRIX(OG, NSCARC_MATRIX_POISSON)
      OP => SCARC_POINT_TO_OTHER_CMATRIX(OG, NSCARC_MATRIX_PROLONGATION)

      OP%N_VAL = OA%N_VAL + 1              ! TODO : CHECK : MUCH TOO BIG !!!
      OP%N_ROW = G%NCE + 1                 ! TODO : CHECK : MUCH TOO BIG !!!
      CALL SCARC_ALLOCATE_CMATRIX(OP, NL, NSCARC_PRECISION_DOUBLE, NSCARC_MATRIX_FULL, 'OG%PROLONGATION', CROUTINE)

   ENDDO

   CALL SCARC_ALLOCATE_REAL1(G%QQ, 1, G%NCE2+1, NSCARC_INIT_ZERO, 'G%QQ', CROUTINE)
   CALL SCARC_ALLOCATE_REAL1(G%RR, 1, G%NCE2+1, NSCARC_INIT_ZERO, 'G%RR', CROUTINE)     ! TODO check length!

 
   ! Copy blocks into Q according to aggregation zones and compute norms for single ZONES
   ! In each cell corresponding to a single zone, store square-sum of entries
 
   IQ = 1
   G%AUX1 = 0.0_EB
   DO ICC = 1, Z%N_ROW-1
      DSUM = 0.0_EB
      DO ICOL = Z%ROW(ICC), Z%ROW(ICC+1)-1
         IC = Z%COL(ICOL)
         G%QQ(IQ) = G%NULLSPACE(IC)
         DSUM = DSUM + G%QQ(IQ)**2
         IQ = IQ + 1
      ENDDO
      DO ICOL = Z%ROW(ICC), Z%ROW(ICC+1)-1
         G%AUX1(Z%COL(ICOL)) = DSUM
      ENDDO
      G%RR(ICC) = DSUM
   ENDDO

   CALL SCARC_REDUCE_REAL1(G%QQ, IQ, 'G%QQ', CROUTINE)
   CALL SCARC_REDUCE_REAL1(G%RR, Z%N_ROW, 'G%RR', CROUTINE)

#ifdef WITH_SCARC_DEBUG
   WRITE(MSG%LU_DEBUG,*) 'NM=',NM,': ================================ PART 0 :'
   WRITE(MSG%LU_DEBUG,*) 'NM=',NM,':Z%N_ROW:', Z%N_ROW
   WRITE(MSG%LU_DEBUG,*) 'NM=',NM,':G%NULLSPACE:'
   WRITE(MSG%LU_DEBUG,'(6E14.6)') G%NULLSPACE
   WRITE(MSG%LU_DEBUG,*) 'NM=',NM,':G%AUX2:'
   WRITE(MSG%LU_DEBUG,'(6E14.6)') G%AUX2
   WRITE(MSG%LU_DEBUG,*) 'NM=',NM,':G%QQ:'
   WRITE(MSG%LU_DEBUG,'(6E14.6)') G%QQ
   WRITE(MSG%LU_DEBUG,*) 'NM=',NM,':G%AUX1:'
   WRITE(MSG%LU_DEBUG,'(6E14.6)') G%AUX1
   WRITE(MSG%LU_DEBUG,*) 'NM=',NM,':G%RR:'
   WRITE(MSG%LU_DEBUG,'(6E14.6)') G%RR
#endif

ENDDO

 
! Exchange sums of nullspace entries within single aggregation zones
 
IF (NMESHES > 1) CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_AUXILIARY, NSCARC_NONE, NL)
   
#ifdef WITH_SCARC_DEBUG
   WRITE(MSG%LU_DEBUG,*) 'NM=',NM,':================================ PART 1 :'
   WRITE(MSG%LU_DEBUG,*) 'NM=',NM,':G%QQ:'
   WRITE(MSG%LU_DEBUG,'(6E14.6)') G%QQ
   WRITE(MSG%LU_DEBUG,*) 'NM=',NM,':G%AUX1:'
   WRITE(MSG%LU_DEBUG,'(6E14.6)') G%AUX1
#endif

 
! Build norms over single zones and scale Q-entries by norms
 
DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID (NM, NL)                          
   Z => SCARC_POINT_TO_CMATRIX(G, NSCARC_MATRIX_ZONES)

   DO ICC = 1, Z%N_ROW-1
      ICOL = Z%ROW(ICC)
      IC = Z%COL(ICOL)
      G%RR(ICC) = SQRT(G%AUX1(IC))
   ENDDO

   IQ = 1
   DO ICC = 1, Z%N_ROW-1
      DO ICOL = Z%ROW(ICC), Z%ROW(ICC+1)-1
         IC = Z%COL(ICOL)
         G%QQ(IQ) = G%QQ(IQ)/G%RR(ICC)
         IQ = IQ + 1
      ENDDO
   ENDDO

ENDDO
   
#ifdef WITH_SCARC_DEBUG
   WRITE(MSG%LU_DEBUG,*) 'NM=',NM,':================================ PART 2 :'
   WRITE(MSG%LU_DEBUG,*) 'NM=',NM,':G%QQ:'
   WRITE(MSG%LU_DEBUG,'(6E14.6)') G%QQ
#endif

! ----------- Relax prolongator:
! Smooth the tentative prolongator, so that it's accuracy is greatly improved for algebraically smooth error.
! Compute:                P =  P - A_Dinv * P   
! with:                   A_Dinv = 4/3 * 1/rho * D^{-1} A   
 
! First step: Compute P_0: = A_Dinv * Q
 
DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID (NM, NL)                                ! Sets grid pointer G

   A => SCARC_POINT_TO_CMATRIX (G, NSCARC_MATRIX_POISSON)
   P => SCARC_POINT_TO_CMATRIX (G, NSCARC_MATRIX_PROLONGATION)
   Z => SCARC_POINT_TO_CMATRIX (G, NSCARC_MATRIX_ZONES)

   MG => L%MG
   MG%OMEGA = 4.0_EB/3.0_EB                                         ! currently used default
   IF (SCARC_MULTIGRID_RELAXING) THEN
      SCAL = MG%OMEGA/MG%APPROX_SPECTRAL_RADIUS                     ! for testing purposes rho is set to 2 currently
   ELSE
      SCAL = 0.0_EB
   ENDIF
   
   IP = 1
   IP0 = IP
   P%ROW(1) = IP
   DO IC = 1, G%NC
      DO ICC = 1, Z%N_ROW-1
   
         DSUM = 0.0_EB
         DO ICCOL = Z%ROW(ICC), Z%ROW(ICC+1)-1
            JC = Z%COL(ICCOL)
            ICOL = SCARC_MATCH_MATRIX_COLUMN(A, IC, JC)
            IF (ICOL /= -1) THEN
               ICC0 = ICC
               DSUM = DSUM - SCAL * G%DIAG(IC) * A%VAL(ICOL) * G%QQ(ICCOL)
            ENDIF
         ENDDO
   
         IF (ABS(DSUM) /= 0.0_EB) THEN
            P%VAL(IP) = DSUM
            P%COL(IP) = ICC
            IP = IP + 1
         ENDIF
   
      ENDDO
 
      ! take care that at least one entry per fine cell is generated
      IF (IP == IP0) THEN
         P%VAL(IP) = 0.0_EB
         P%COL(IP) = ICC0
         IP = IP + 1
      ENDIF
      IP0 = IP

      P%ROW(IC+1) = IP
   ENDDO
 
   P%N_VAL = IP - 1

ENDDO
   
 
! Second step: Compute P: = P - P_0
 
DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G
   P => SCARC_POINT_TO_CMATRIX(G, NSCARC_MATRIX_PROLONGATION)
   Z => SCARC_POINT_TO_CMATRIX(G, NSCARC_MATRIX_ZONES)

#ifdef WITH_SCARC_DEBUG
   WRITE(MSG%LU_DEBUG,*) 'P%N_VAL=',P%N_VAL
   WRITE(MSG%LU_DEBUG,*) 'P%N_ROW=',P%N_ROW
   CALL SCARC_DEBUG_CMATRIX(Z, 'ZONES','AFTER RESORT PROL ')
   CALL SCARC_DEBUG_CMATRIX(P, 'PROLONGATION','AFTER RESORT PROL ')
#endif

   DO ICC = 1, Z%N_ROW-1
      DO ICCOL = Z%ROW(ICC), Z%ROW(ICC+1) - 1
         IC = Z%COL(ICCOL)

         IF (IC > G%NC) CYCLE
         DO JCCOL = P%ROW(IC), P%ROW(IC+1) - 1
            JCC = P%COL(JCCOL)
            IF (JCC == ICC) THEN
               P%VAL(JCCOL) = P%VAL(JCCOL) + G%QQ(ICCOL)
               CYCLE
            ENDIF
         ENDDO
      ENDDO
   ENDDO
   
#ifdef WITH_SCARC_DEBUG
   CALL SCARC_DEBUG_CMATRIX(P, 'G%PROLONGATION','SETUP_PROLONGATION: AFTER RELAX STEP, BEFORE EXCHANGE ')
#endif
ENDDO

! Determine global columns array for Prolongation matrix
 
DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_MULTIGRID (NM, NL, NL+1)                
   PF => SCARC_POINT_TO_CMATRIX(GF, NSCARC_MATRIX_PROLONGATION)

   DO IC = 1, GF%NC
      DO JCCOL = PF%ROW(IC), PF%ROW(IC+1) - 1
         JCC = PF%COL(JCCOL)
         PF%COLG(JCCOL) = GC%LOCAL_TO_GLOBAL(JCC)
      ENDDO
   ENDDO
   
#ifdef WITH_SCARC_DEBUG
   WRITE(MSG%LU_DEBUG,*) '----------- NM =',NM,': NL=',NL
   CALL SCARC_DEBUG_CMATRIX(PF, 'G%PROLONGATION','SETUP_PROLONGATION: AFTER LAST EXCHANGE')
#endif

ENDDO

 
! Exchange resulting columns and values of Prolongation matrix and extract exchanged data from 
! overlapping parts with single neighbors and attach them to main matrix
 
IF (NMESHES > 1) THEN
   CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_MATRIX_COLS,  NSCARC_MATRIX_PROLONGATION, NL)
   CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_MATRIX_COLSG, NSCARC_MATRIX_PROLONGATION, NL)
   CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_MATRIX_VALS,  NSCARC_MATRIX_PROLONGATION, NL)
   CALL SCARC_EXTRACT_MATRIX_OVERLAPS(NSCARC_MATRIX_PROLONGATION, 0, NL)
ENDIF
   
DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_MULTIGRID (NM, NL, NL+1) 
   PF => SCARC_POINT_TO_CMATRIX(GF, NSCARC_MATRIX_PROLONGATION)

   NLEN = MAX(GC%NCE2, 4 * (GC%NCE2 - GC%NC + 2) + 50)           ! TODO check length
   CALL SCARC_ALLOCATE_INT1 (GC%CELLS_LOCAL,  1, NLEN, NSCARC_INIT_ZERO, 'GC%CELLS_LOCAL',  CROUTINE)
   CALL SCARC_ALLOCATE_INT1 (GC%CELLS_GLOBAL, 1, NLEN, NSCARC_INIT_ZERO, 'GC%CELLS_GLOBAL', CROUTINE)

   CALL SCARC_GET_CELL_DEPENDENCIES_GALERKIN(GC, GF, PF, NLEN)

   CALL SCARC_REDUCE_INT1(GC%CELLS_LOCAL, GC%NC_GALERKIN, 'GC%CELLS_LOCAL', CROUTINE)
   CALL SCARC_REDUCE_INT1(GC%CELLS_GLOBAL, GC%NC_GALERKIN, 'GC%CELLS_GLOBAL', CROUTINE)

   CALL SCARC_REDUCE_CMATRIX(PF, 'P%PROLONGATION', CROUTINE)
   DO INBR = 1, S%N_NEIGHBORS
      NOM = S%NEIGHBORS(INBR)
      CALL SCARC_POINT_TO_OTHER_GRID(NM, NOM, NL)
      OPF => SCARC_POINT_TO_OTHER_CMATRIX(OGF, NSCARC_MATRIX_PROLONGATION)
      CALL SCARC_REDUCE_CMATRIX(OPF, 'OP%PROLONGATION', CROUTINE)
   ENDDO

#ifdef WITH_SCARC_DEBUG
   WRITE(MSG%LU_DEBUG,*) '======================= LEVEL ', NL
   CALL SCARC_DEBUG_CMATRIX(PF, 'GF%PROLONGATION','SETUP_PROLONGATION: FINAL')
#endif

ENDDO

END SUBROUTINE SCARC_SETUP_PROLONGATION_AMG


! -------------------------------------------------------------------------------------------
!> \brief Setup restriction and prolongation matrices in case of GMG-like coarsening
! -------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_TRANSFER_GMG(NL)
USE SCARC_POINTERS, ONLY:  S, LC, LF, GF, GC, AF, RF, ZF, PF, OPF, OGF
USE SCARC_POINTER_ROUTINES, ONLY: SCARC_POINT_TO_OTHER_GRID, SCARC_POINT_TO_MULTIGRID, &
                                  SCARC_POINT_TO_CMATRIX, SCARC_POINT_TO_OTHER_CMATRIX
INTEGER, INTENT(IN) :: NL
INTEGER :: NM, IC, JC, ICC, ICF, ICCOL, IRCOL, IP, JCC, JCCOL, NOM, INBR, NLEN
INTEGER :: IS, IXC, IZC, IXF, IZF, IOFFX, IOFFZ, ICC0(4)
INTEGER :: STENCIL(16) = 0
INTEGER :: RN, RE, RS, RW

CROUTINE = 'SCARC_SETUP_PROLONGATION_GMG'

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_MULTIGRID (NM, NL, NL + 1)                                   ! Sets grid pointer G

   AF => SCARC_POINT_TO_CMATRIX(GF, NSCARC_MATRIX_POISSON)
   ZF => SCARC_POINT_TO_CMATRIX(GF, NSCARC_MATRIX_ZONES)
   RF => SCARC_POINT_TO_CMATRIX(GF, NSCARC_MATRIX_RESTRICTION)
   PF => SCARC_POINT_TO_CMATRIX(GF, NSCARC_MATRIX_PROLONGATION)

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'TRANSFER_GMG: GF%NC=', GF%NC
WRITE(MSG%LU_DEBUG,*) 'TRANSFER_GMG: GC%NC=', GC%NC
WRITE(MSG%LU_DEBUG,*) 'TRANSFER_GMG: TYPE_INTERPOL=', TYPE_INTERPOL
#endif

   SELECT CASE (TYPE_INTERPOL)

      CASE (NSCARC_INTERPOL_CONSTANT)

         RF%N_VAL = GF%NCE2
         RF%N_ROW = GC%NCE2 + 1
         CALL SCARC_ALLOCATE_CMATRIX(RF, NL, NSCARC_PRECISION_DOUBLE, NSCARC_MATRIX_FULL, 'GF%RESTRICTION', CROUTINE)

         PF%N_VAL = AF%N_VAL + 1        
         PF%N_ROW = GF%NCE2 + 1
         CALL SCARC_ALLOCATE_CMATRIX(PF, NL, NSCARC_PRECISION_DOUBLE, NSCARC_MATRIX_FULL, 'GF%PROLONGATION', CROUTINE)

         IP = 1
         RF%ROW(1) = IP
         DO ICC = 1, ZF%N_ROW-1
            DO ICCOL = ZF%ROW(ICC), ZF%ROW(ICC+1)-1
               JC = ZF%COL(ICCOL)
               RF%COL(IP) = ZF%COL(ICCOL)
               RF%COLG(IP) = GF%LOCAL_TO_GLOBAL(ZF%COL(ICCOL))
               RF%VAL(IP) = 1.00_EB
               IP = IP + 1
            ENDDO
            RF%ROW(ICC + 1) = IP
            !WRITE(*,*) ICC,':SUM(RF)=',SUM(RF%VAL(RF%ROW(ICC):RF%ROW(ICC+1)-1))
         ENDDO
      
         IP = 1
         PF%ROW(1) = IP
         DO IC = 1, GF%NC
            DO ICC = 1, ZF%N_ROW -1
               COLUMN_LOOP: DO ICCOL = ZF%ROW(ICC), ZF%ROW(ICC+1)-1
                  IF (ZF%COL(ICCOL) == IC) THEN
                     PF%COL(IP) = ICC
                     PF%COLG(IP) = GC%LOCAL_TO_GLOBAL(ICC)
                     PF%VAL(IP) = 0.25_EB
                     IP = IP + 1
                     EXIT COLUMN_LOOP
                  ENDIF
               ENDDO COLUMN_LOOP
            ENDDO
            PF%ROW(IC + 1) = IP
            !WRITE(*,*) ICC,':SUM(PF)=',SUM(PF%VAL(PF%ROW(IC):PF%ROW(IC+1)-1))
         ENDDO

      CASE (NSCARC_INTERPOL_BILINEAR)

         RF%N_VAL = 16*GF%NC
         RF%N_ROW = GC%NC + 1
         CALL SCARC_ALLOCATE_CMATRIX(RF, NL, NSCARC_PRECISION_DOUBLE, NSCARC_MATRIX_FULL, 'GF%RESTRICTION', CROUTINE)

         PF%N_VAL = AF%N_VAL + 1        
         PF%N_ROW = GF%NC + 1
         CALL SCARC_ALLOCATE_CMATRIX(PF, NL, NSCARC_PRECISION_DOUBLE, NSCARC_MATRIX_FULL, 'GF%PROLONGATION', CROUTINE)

         IP = 1
         RF%ROW(1) = IP
         DO IZC = 1, LC%NZ
            DO IXC = 1, LC%NX

               ICC = (IZC - 1) * LC%NX + IXC

               IXF = 2 * IXC
               IZF = 2 * IZC
               
               IF (IXC == 1 .AND. IZC == 1) THEN
                  RN = 1; RE = 1; RS = 0; RW = 0
               ELSE IF (IXC == 1 .AND. IZC == LC%NZ) THEN
                  RN = 0; RE = 1; RS = 1; RW = 0
               ELSE IF (IXC == LC%NX .AND. IZC == 1) THEN
                  RN = 1; RE = 0; RS = 0; RW = 1
               ELSE IF (IXC == LC%NX .AND. IZC == LC%NZ) THEN
                  RN = 0; RE = 0; RS = 1; RW = 1
               ELSE IF (IXC == 1) THEN
                  RN = 1; RE = 1; RS = 1; RW = 0
               ELSE IF (IXC == LC%NX) THEN
                  RN = 1; RE = 0; RS = 1; RW = 1
               ELSE IF (IZC == 1) THEN
                  RN = 1; RE = 1; RS = 0; RW = 1
               ELSE IF (IZC == LC%NZ) THEN
                  RN = 0; RE = 1; RS = 1; RW = 1
               ELSE
                  RN = 1; RE = 1; RS = 1; RW = 1
               ENDIF

               STENCIL(13:16) =  (/    RN *RW,    RN *(2+RW),    RN *(2+RE),    RN *RE /)
               STENCIL( 9:12) =  (/ (2+RN)*RW, (2+RN)*(2+RW), (2+RN)*(2+RE), (2+RN)*RE /)
               STENCIL( 5: 8) =  (/ (2+RS)*RW, (2+RS)*(2+RW), (2+RS)*(2+RE), (2+RS)*RE /)
               STENCIL( 1: 4) =  (/    RS *RW,    RS *(2+RW),    RS *(2+RE),    RS *RE /)

               IS = 1
               DO IOFFZ = -2, 1
                  DO IOFFX = -2, 1
                     CALL PROCESS_FINE_CELL (IXF + IOFFX, 1, IZF + IOFFZ, IP, STENCIL(IS))
                     IS = IS + 1
                  ENDDO
               ENDDO

               RF%ROW(ICC + 1) = IP
               !WRITE(*,*) ICC,':SUM(RF)=',SUM(RF%VAL(RF%ROW(ICC):RF%ROW(ICC+1)-1))


            ENDDO
         ENDDO

         IP = 1
         PF%ROW(1) = IP
         DO IC = 1, GF%NC
            DO ICC = 1, GC%NC 
               COLUMN_LOOP2: DO IRCOL = RF%ROW(ICC), RF%ROW(ICC+1)-1
                  IF (RF%COL(IRCOL) == IC) THEN
                     PF%COL(IP) = ICC
                     PF%COLG(IP) = GC%LOCAL_TO_GLOBAL(ICC)
                     PF%VAL(IP) = RF%VAL(IRCOL)
                     IP = IP + 1
                     EXIT COLUMN_LOOP2
                  ENDIF
               ENDDO COLUMN_LOOP2
            ENDDO
            PF%ROW(IC + 1) = IP
            !WRITE(*,*) ICC,':SUM(PF)=',SUM(PF%VAL(PF%ROW(IC):PF%ROW(IC+1)-1))
         ENDDO

         PF%VAL = PF%VAL/16.0_EB
         IF (TWO_D) THEN
            RF%VAL = RF%VAL/4.0_EB
         ELSE
            RF%VAL = RF%VAL/2.0_EB
         ENDIF

      CASE (NSCARC_INTERPOL_BILINEAR2)

         RF%N_VAL = 16*GF%NC
         RF%N_ROW = GC%NC + 1
         CALL SCARC_ALLOCATE_CMATRIX(RF, NL, NSCARC_PRECISION_DOUBLE, NSCARC_MATRIX_FULL, 'GF%RESTRICTION', CROUTINE)

         PF%N_VAL = 4* GF%NC
         PF%N_ROW = GF%NC + 1
         CALL SCARC_ALLOCATE_CMATRIX(PF, NL, NSCARC_PRECISION_DOUBLE, NSCARC_MATRIX_FULL, 'GF%PROLONGATION', CROUTINE)

         IP = 1
         PF%ROW(1) = IP
         DO IZF = 1, LF%NZ
            DO IXF = 1, LF%NX

               ICF = GF%CELL_NUMBER(IXF, 1, IZF)

               !IXC = CEILING(REAL((IXF+1)/2),EB) 
               !IZC = CEILING(REAL((IZF+1)/2),EB)
               IXC = MOD(IXF,2) + 1
               IZC = MOD(IZF,2) + 1

               ICC0(1) = (IZC - 1) * LC%NX + IXC
               ICC0(2) = (IZC - 1) * LC%NX + IXC + 1
               ICC0(3) = IZC * LC%NX + IXC
               ICC0(4) = IZC * LC%NX + IXC + 1

#ifdef WITH_SCARC_DEBUG
 WRITE(MSG%LU_DEBUG,*) 'IXF, IZF, ICF, IXC, IZC, ICC0(1:4):', IXF, IZF, ICF, IXC, IZC, ICC0(1:4)
#endif
 WRITE(*,*) 'IXF, IZF, ICF, IXC, IZC, ICC0(1:4):', IXF, IZF, ICF, IXC, IZC, ICC0(1:4)
               
 IF (IXC /= 0 .AND. IZC /= 0) CALL PROCESS_COARSE_CELL (ICC0(1), IP, 9)
 IF (IXC /= 0 .AND. IZC /= 0) CALL PROCESS_COARSE_CELL (ICC0(2), IP, 3)
 IF (IXC /= 0 .AND. IZC /= 0) CALL PROCESS_COARSE_CELL (ICC0(3), IP, 3)
 IF (IXC /= 0 .AND. IZC /= 0) CALL PROCESS_COARSE_CELL (ICC0(4), IP, 1)

               PF%ROW(ICF + 1) = IP

            ENDDO
         ENDDO

         IP = 1
         PF%ROW(1) = IP
         DO IC = 1, GF%NC
            DO ICC = 1, GC%NC 
               COLUMN_LOOP3: DO IRCOL = RF%ROW(ICC), RF%ROW(ICC+1)-1
                  IF (RF%COL(IRCOL) == IC) THEN
                     PF%COL(IP) = ICC
                     PF%COLG(IP) = GC%LOCAL_TO_GLOBAL(ICC)
                     PF%VAL(IP) = RF%VAL(IRCOL)
                     IP = IP + 1
                     EXIT COLUMN_LOOP3
                  ENDIF
               ENDDO COLUMN_LOOP3
            ENDDO
            PF%ROW(IC + 1) = IP
            WRITE(*,*) ICC,':SUM(PF)=',SUM(PF%VAL(PF%ROW(IC):PF%ROW(IC+1)-1))
         ENDDO

         PF%VAL = PF%VAL/16.0_EB
         IF (TWO_D) THEN
            RF%VAL = RF%VAL/4.0_EB
         ELSE
            RF%VAL = RF%VAL/2.0_EB
         ENDIF

   END SELECT

#ifdef WITH_SCARC_DEBUG
   CALL SCARC_DEBUG_CMATRIX(RF, 'RESTRICTION','IN TRANSFER GMG')
   CALL SCARC_DEBUG_CMATRIX(PF, 'PROLONGATION','IN TRANSFER GMG')
#endif

ENDDO

! Determine global columns array for Prolongation matrix
 
DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_MULTIGRID (NM, NL, NL+1)                
   PF => SCARC_POINT_TO_CMATRIX(GF, NSCARC_MATRIX_PROLONGATION)

   DO IC = 1, GF%NC
      DO JCCOL = PF%ROW(IC), PF%ROW(IC+1) - 1
         JCC = PF%COL(JCCOL)
         PF%COLG(JCCOL) = GC%LOCAL_TO_GLOBAL(JCC)
      ENDDO
   ENDDO
   
#ifdef WITH_SCARC_DEBUG
   WRITE(MSG%LU_DEBUG,*) '----------- NM =',NM,': NL=',NL
   CALL SCARC_DEBUG_CMATRIX(PF, 'G%PROLONGATION','SETUP_PROLONGATION: AFTER LAST EXCHANGE')
#endif

ENDDO

 
! Exchange resulting columns and values of Prolongation matrix and extract exchanged data from 
! overlapping parts with single neighbors and attach them to main matrix
 
IF (NMESHES > 1 .AND. TYPE_COARSENING /= NSCARC_COARSENING_GMG) THEN
   CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_MATRIX_COLS,  NSCARC_MATRIX_PROLONGATION, NL)
   CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_MATRIX_COLSG, NSCARC_MATRIX_PROLONGATION, NL)
   CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_MATRIX_VALS,  NSCARC_MATRIX_PROLONGATION, NL)
   CALL SCARC_EXTRACT_MATRIX_OVERLAPS(NSCARC_MATRIX_PROLONGATION, 0, NL)
ENDIF
   
DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_MULTIGRID (NM, NL, NL+1) 
   PF => SCARC_POINT_TO_CMATRIX(GF, NSCARC_MATRIX_PROLONGATION)

   NLEN = MAX(GC%NCE2, 4 * (GC%NCE2 - GC%NC + 2) + 50)           ! TODO check length
   CALL SCARC_ALLOCATE_INT1 (GC%CELLS_LOCAL,  1, NLEN, NSCARC_INIT_ZERO, 'GC%CELLS_LOCAL',  CROUTINE)
   CALL SCARC_ALLOCATE_INT1 (GC%CELLS_GLOBAL, 1, NLEN, NSCARC_INIT_ZERO, 'GC%CELLS_GLOBAL', CROUTINE)

   CALL SCARC_GET_CELL_DEPENDENCIES_GALERKIN(GC, GF, PF, NLEN)

   CALL SCARC_REDUCE_INT1(GC%CELLS_LOCAL, GC%NC_GALERKIN, 'GC%CELLS_LOCAL', CROUTINE)
   CALL SCARC_REDUCE_INT1(GC%CELLS_GLOBAL, GC%NC_GALERKIN, 'GC%CELLS_GLOBAL', CROUTINE)

   CALL SCARC_REDUCE_CMATRIX(PF, 'P%PROLONGATION', CROUTINE)
   IF (TYPE_COARSENING /= NSCARC_COARSENING_GMG) THEN
      DO INBR = 1, S%N_NEIGHBORS
         NOM = S%NEIGHBORS(INBR)
         CALL SCARC_POINT_TO_OTHER_GRID(NM, NOM, NL)
         OPF => SCARC_POINT_TO_OTHER_CMATRIX(OGF, NSCARC_MATRIX_PROLONGATION)
         CALL SCARC_REDUCE_CMATRIX(OPF, 'OP%PROLONGATION', CROUTINE)
      ENDDO
   ENDIF

#ifdef WITH_SCARC_DEBUG
   WRITE(MSG%LU_DEBUG,*) '======================= LEVEL ', NL
   CALL SCARC_DEBUG_CMATRIX(PF, 'GF%PROLONGATION','SETUP_PROLONGATION: FINAL')
#endif

ENDDO

END SUBROUTINE SCARC_SETUP_TRANSFER_GMG

SUBROUTINE PROCESS_FINE_CELL(IXF, IYF, IZF, IP, VAL)
USE SCARC_POINTERS, ONLY: RF, GF
INTEGER, INTENT(IN) :: IXF, IYF, IZF, VAL
INTEGER, INTENT(INOUT) :: IP
IF (VAL == 0) RETURN
RF%COL(IP) = ICF (IXF, IYF, IZF)
RF%COLG(IP) = GF%LOCAL_TO_GLOBAL(RF%COL(IP))
RF%VAL(IP) = VAL
IP = IP + 1
END SUBROUTINE PROCESS_FINE_CELL

SUBROUTINE PROCESS_COARSE_CELL(ICC, IP, VAL)
USE SCARC_POINTERS, ONLY: PF, GC
INTEGER, INTENT(IN) :: ICC, VAL
INTEGER, INTENT(INOUT) :: IP
IF (VAL == 0) RETURN
PF%COL(IP) = ICC
PF%COLG(IP) = GC%LOCAL_TO_GLOBAL(ICC)
PF%VAL(IP) = VAL
IP = IP + 1
END SUBROUTINE PROCESS_COARSE_CELL

INTEGER FUNCTION ICF (IXF, IYF, IZF)
USE SCARC_POINTERS, ONLY : LF
INTEGER, INTENT(IN) :: IXF, IYF, IZF
ICF = (IZF - 1) * LF%NX * LF%NY + (IYF - 1) * LF%NX + IXF
RETURN
END FUNCTION

INTEGER FUNCTION ICC (IXC, IYC, IZC)
USE SCARC_POINTERS, ONLY : LC
INTEGER, INTENT(IN) :: IXC, IYC, IZC
ICC = (IZC - 1) * LC%NX * LC%NY + (IYC - 1) * LC%NX + IXC
RETURN
END FUNCTION

! -------------------------------------------------------------------------------------------
!> \brief Determine on which overlapping global coarse cells are given mesh depends (also considering diagonal connections)
! -------------------------------------------------------------------------------------------
SUBROUTINE SCARC_GET_CELL_DEPENDENCIES_GALERKIN(GC, GF, PF, NLEN)
TYPE (SCARC_GRID_TYPE), POINTER, INTENT(IN) :: GC, GF
TYPE (SCARC_CMATRIX_TYPE), POINTER, INTENT(IN) :: PF
INTEGER, INTENT(IN) :: NLEN
INTEGER :: IZL, IZG, IP, ICOL, IC, IFOUND1, IFOUND2, IFOUND3

IP = 1
PROLONGATION_CELLS_LOOP: DO IC = 1, GF%NCE

   ! Check if zone number used in given row of Prolongation matrix is already accounted for
   DO ICOL = PF%ROW(IC), PF%ROW(IC+1) - 1

      IZL = PF%COL(ICOL)
      IZG = PF%COLG(ICOL)

      IFOUND1 = FINDLOC (GC%CELLS_GLOBAL(1:NLEN),  VALUE = IZG, DIM = 1)
      IFOUND2 = FINDLOC (GF%ZONES_GLOBAL(1:GF%NC), VALUE = IZG, DIM = 1)
      IFOUND3 = FINDLOC (GC%LOCAL_TO_GLOBAL(1:GC%NCE2), VALUE = IZG, DIM = 1)
      IF (IFOUND1 <= 0 .AND. IFOUND2 <= 0) THEN  
         GC%CELLS_LOCAL(IP)  = IFOUND3
         GC%CELLS_GLOBAL(IP) = IZG
         IP = IP + 1
      ENDIF
   ENDDO

ENDDO PROLONGATION_CELLS_LOOP
GC%NC_GALERKIN = IP - 1

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'GET_CELL_DEPENDENCIES: GC%NC=',GC%NC
WRITE(MSG%LU_DEBUG,*) 'GET_CELL_DEPENDENCIES: GC%NCE=',GC%NCE
WRITE(MSG%LU_DEBUG,*) 'GET_CELL_DEPENDENCIES: GC%NCE2=',GC%NCE2
WRITE(MSG%LU_DEBUG,*) 'GET_CELL_DEPENDENCIES: GC%NC_GALERKIN=', GC%NC_GALERKIN
WRITE(MSG%LU_DEBUG,*) 'GET_CELL_DEPENDENCIES: GC%CELLS_LOCAL'
WRITE(MSG%LU_DEBUG,'(8I6)') GC%CELLS_LOCAL(1:NLEN)
WRITE(MSG%LU_DEBUG,*) 'GET_CELL_DEPENDENCIES: GC%CELLS_GLOBAL'
WRITE(MSG%LU_DEBUG,'(8I6)') GC%CELLS_GLOBAL(1:NLEN)
#endif

END SUBROUTINE SCARC_GET_CELL_DEPENDENCIES_GALERKIN


! ------------------------------------------------------------------------------------------------------
!> \brief Define nullspace for next coarser level, if coarsest level isn't reached yet
! ------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_NULLSPACE_COARSE(NL)
USE SCARC_POINTERS, ONLY:  GC, GF, ZF
USE SCARC_POINTER_ROUTINES, ONLY: SCARC_POINT_TO_MULTIGRID, SCARC_POINT_TO_CMATRIX
INTEGER, INTENT(IN) :: NL
INTEGER :: NM

CROUTINE = 'SCARC_SETUP_NULLSPACE_COARSE'

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_MULTIGRID(NM, NL, NL+1)                   ! Sets pointers GC and GF
   ZF => SCARC_POINT_TO_CMATRIX (GF, NSCARC_MATRIX_ZONES)

   IF (NL < NLEVEL_MAX) THEN
      GC%N_FINE = GC%NC_LOCAL(NM)
      CALL SCARC_ALLOCATE_REAL1(GC%NULLSPACE, 1, GC%NCE, NSCARC_INIT_ZERO, 'GC%NULLSPACE', CROUTINE)
      GC%NULLSPACE(1:GC%NCE) = GF%RR(1:GC%NCE)
      CALL SCARC_REDUCE_INT1(GC%LOCAL_TO_GLOBAL, GC%NCE, 'GC%LOCAL_TO_GLOBAL', CROUTINE)
   ENDIF

#ifdef WITH_SCARC_DEBUG
   WRITE(MSG%LU_DEBUG,*) 'NM=',NM,'============== SETUP_NULLSPACE_COARSE ================='
   WRITE(MSG%LU_DEBUG,*) 'NM=',NM,': GF%NULLSPACE:'
   WRITE(MSG%LU_DEBUG,'(6E14.6)') GF%NULLSPACE
   WRITE(MSG%LU_DEBUG,*) 'NM=',NM,': GF%RR'
   WRITE(MSG%LU_DEBUG,'(6E14.6)') GF%RR
   WRITE(MSG%LU_DEBUG,*) 'NM=',NM,': GF%QR:'
   WRITE(MSG%LU_DEBUG,'(6E14.6)') GF%QQ
   WRITE(MSG%LU_DEBUG,*) 'NM=',NM,': GC%LOCAL_TO_GLOBAL:'
   WRITE(MSG%LU_DEBUG,'(8I6)') GC%LOCAL_TO_GLOBAL
   WRITE(MSG%LU_DEBUG,*) 'NM=',NM,': GC%N_FINE:', GC%N_FINE
   WRITE(MSG%LU_DEBUG,*) 'NM=',NM,': GC%NC:', GC%NC
   WRITE(MSG%LU_DEBUG,*) 'NM=',NM,': GC%NCE:', GC%NCE
   WRITE(MSG%LU_DEBUG,*) '==============================================='
#endif

ENDDO
   
END SUBROUTINE SCARC_SETUP_NULLSPACE_COARSE


! ------------------------------------------------------------------------------------------------------
!> \brief Determine which columns of system matrix are involved in multiplication with tentative prolongator
! ------------------------------------------------------------------------------------------------------
INTEGER FUNCTION SCARC_MATCH_MATRIX_COLUMN(A, IC, JC)
TYPE (SCARC_CMATRIX_TYPE), POINTER :: A
INTEGER, INTENT(IN) :: IC, JC
INTEGER :: ICOL
SCARC_MATCH_MATRIX_COLUMN = -1
DO ICOL = A%ROW(IC), A%ROW(IC+1)-1
   IF (A%COL(ICOL) == JC) THEN
      SCARC_MATCH_MATRIX_COLUMN = ICOL
      RETURN 
   ENDIF
ENDDO
END FUNCTION SCARC_MATCH_MATRIX_COLUMN


! ------------------------------------------------------------------------------------------------------
!> \brief Setup Restriction matrix: Build transpose of Prolongation matrix
! Compute the Restriction matrix, R, which interpolates from the fine-grid to the coarse-grid.
! ------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_RESTRICTION(NL)
USE SCARC_POINTERS, ONLY: GC, GF, RF, PF
USE SCARC_POINTER_ROUTINES, ONLY: SCARC_POINT_TO_MULTIGRID, SCARC_POINT_TO_CMATRIX
INTEGER, INTENT(IN) :: NL
INTEGER :: NM, IC, IRROW, IRCOL, IPCOL, IPC, ICCL, ICCG, IFOUND 
LOGICAL :: IS_INCLUDED

CROUTINE = 'SCARC_SETUP_RESTRICTION'

! Allocate Restriction matrix R on internal mesh part
MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_MULTIGRID (NM, NL, NL+1)            ! Sets pointers GF and GC to fine and coarse level

   PF => SCARC_POINT_TO_CMATRIX (GF, NSCARC_MATRIX_PROLONGATION)
   RF => SCARC_POINT_TO_CMATRIX (GF, NSCARC_MATRIX_RESTRICTION)

   RF%N_VAL = PF%N_VAL + 100
   RF%N_ROW = PF%N_ROW + 100

   CALL SCARC_ALLOCATE_CMATRIX(RF, NL, NSCARC_PRECISION_DOUBLE, NSCARC_MATRIX_FULL, 'GF%RESTRICTION', CROUTINE)

   IRROW = 1                             ! counter of current row of Restriction matrix - corresponds to coarse cells
   IRCOL = 1                             ! counter of current column of Restriction matrix
   RF%ROW(IRROW) = IRCOL

   LOCAL_COARSE_CELLS_LOOP: DO ICCL = 1, GC%NCE

      ICCG = GC%LOCAL_TO_GLOBAL(ICCL)                     ! corresponding global coarse cell

      IFOUND = -1
      IFOUND = FINDLOC(PF%COLG, VALUE = ICCG, DIM=1)
      IF (IFOUND == -1) CYCLE

      IS_INCLUDED = .FALSE.

      FINE_CELLS_LOOP: DO IC = 1, GF%NCE                  ! counter of fine cell (including overlaps)

         ROW_LOOP: DO IPCOL = PF%ROW(IC), PF%ROW(IC+1)-1
            IPC = PF%COLG(IPCOL)
            IF (IPC == ICCG) THEN
               IS_INCLUDED = .TRUE.
               RF%VAL(IRCOL) = PF%VAL(IPCOL)
               RF%COLG(IRCOL) = IC
               IRCOL = IRCOL + 1
               EXIT ROW_LOOP
            ENDIF
         ENDDO ROW_LOOP
      ENDDO FINE_CELLS_LOOP

      IF (IS_INCLUDED) THEN 
         RF%ROW(IRROW+1) = IRCOL
         IRROW = IRROW + 1
      ENDIF

      RF%N_ROW = IRROW 
      RF%N_VAL = IRCOL - 1

   ENDDO LOCAL_COARSE_CELLS_LOOP

   CALL SCARC_REDUCE_CMATRIX (RF, 'GF%RESTRICTION', CROUTINE)

#ifdef WITH_SCARC_DEBUG
   WRITE(MSG%LU_DEBUG,*) '--------------> NM=',NM
   CALL SCARC_DEBUG_CMATRIX (RF, 'GF%RESTRICTION','AFTER SETUP_RESTRICTION')
#endif

ENDDO MESHES_LOOP

END SUBROUTINE SCARC_SETUP_RESTRICTION


! ------------------------------------------------------------------------------------------------------
!> \brief Find matching column index during matrix-matrix multiplication of compact matrices
! ------------------------------------------------------------------------------------------------------
INTEGER FUNCTION SCARC_FIND_MATCHING_COLUMN(P, JC, ICCG)
TYPE (SCARC_CMATRIX_TYPE), POINTER :: P
INTEGER, INTENT(IN) :: JC, ICCG
INTEGER :: IPCOL

SCARC_FIND_MATCHING_COLUMN = -1
IF (JC == 0) RETURN
DO IPCOL = P%ROW(JC), P%ROW(JC+1)-1
   IF (P%COLG(IPCOL) == ICCG) THEN
      SCARC_FIND_MATCHING_COLUMN = IPCOL
      RETURN 
   ENDIF
ENDDO

END FUNCTION SCARC_FIND_MATCHING_COLUMN


! ------------------------------------------------------------------------------------------------------
!> \brief Find matching components to multiply row of Poisson matrix with column of Prolongation matrix
! ------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_MULTIPLY_POISSON_PROL(A, P, PP, ICC, ICC0, IP, IC)
TYPE (SCARC_CMATRIX_TYPE), POINTER, INTENT(IN) :: A, P, PP
INTEGER, INTENT(IN) :: ICC, IC
INTEGER, INTENT(INOUT) :: IP, ICC0
REAL(EB) :: DSUM, TOL = 1E-12_EB
INTEGER :: IACOL, IPCOL, JC

DSUM = 0.0_EB
DO IACOL = A%ROW(IC), A%ROW(IC+1)-1
   JC = A%COL(IACOL)
#ifdef WITH_SCARC_DEBUG2
WRITE(MSG%LU_DEBUG,*) 'MULTIPLY_POISSON_PROL: IC, ICC, IP, JC:', IC, ICC, IP, JC
#endif
   IF (JC == 0) CYCLE
   IPCOL = SCARC_FIND_MATCHING_COLUMN(P, JC, ICC)
#ifdef WITH_SCARC_DEBUG2
WRITE(MSG%LU_DEBUG,*) 'MULTIPLY_POISSON_PROL: IPCOL =', IPCOL
#endif
   IF (JC < 0 .OR. IPCOL <= 0) CYCLE
   ICC0 = ICC
   DSUM = DSUM + A%VAL(IACOL) * P%VAL(IPCOL)
#ifdef WITH_SCARC_DEBUG2
WRITE(MSG%LU_DEBUG,*) 'MULTIPLY_POISSON_PROL: DSUM, A%VAL, P%VAL:', DSUM, A%VAL(IACOL), P%VAL(IPCOL)
#endif
ENDDO

IF (ABS(DSUM) > TOL) THEN
   PP%COL(IP)  = ICC
   PP%COLG(IP) = ICC
   PP%VAL(IP)  = DSUM
   IP = IP + 1
ENDIF

END SUBROUTINE SCARC_MULTIPLY_POISSON_PROL


! ------------------------------------------------------------------------------------------------------
!> \brief Perform matrix multiplication between fine Poisson matrix and Prolongation matrix 
! ------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_POISSON_PROL(NL)
USE SCARC_POINTERS, ONLY: GC, GF, AF, PF, PPF, OA, OPP, OGF
USE SCARC_POINTER_ROUTINES, ONLY: SCARC_POINT_TO_OTHER_GRID, SCARC_POINT_TO_MULTIGRID, &
                                  SCARC_POINT_TO_CMATRIX, SCARC_POINT_TO_OTHER_CMATRIX
INTEGER, INTENT(IN) :: NL
INTEGER  :: NM, NOM, IC, IP, IP0, INBR, ICC, ICC0 = -1
REAL(EB) :: TNOW, TSUM = 0.0_EB

CROUTINE = 'SCARC_SETUP_POISSON_PROL'

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_MULTIGRID (NM, NL, NL+1)

   AF  => SCARC_POINT_TO_CMATRIX (GF, NSCARC_MATRIX_POISSON)          
   PF  => SCARC_POINT_TO_CMATRIX (GF, NSCARC_MATRIX_PROLONGATION)     
   PPF => SCARC_POINT_TO_CMATRIX (GF, NSCARC_MATRIX_POISSON_PROL)         

   !PPF%N_ROW = AF%N_ROW
   PPF%N_ROW = GF%NCE+1
   PPF%N_VAL = PPF%N_ROW*30            ! TODO: only temporarily
   CALL SCARC_ALLOCATE_CMATRIX(PPF, NL, NSCARC_PRECISION_DOUBLE, NSCARC_MATRIX_FULL, 'GF%PPF', CROUTINE)

   DO INBR = 1, SCARC(NM)%N_NEIGHBORS

      NOM = SCARC(NM)%NEIGHBORS(INBR)
      CALL SCARC_POINT_TO_OTHER_GRID(NM, NOM, NL)

      OA  => SCARC_POINT_TO_OTHER_CMATRIX(OGF, NSCARC_MATRIX_POISSON)
      OPP => SCARC_POINT_TO_OTHER_CMATRIX(OGF, NSCARC_MATRIX_POISSON_PROL)

      OPP%N_VAL = AF%N_VAL              ! TODO : CHECK : MUCH TOO BIG !!!
      OPP%N_ROW = GF%NCE + 1            
      CALL SCARC_ALLOCATE_CMATRIX(OPP, NL, NSCARC_PRECISION_DOUBLE, NSCARC_MATRIX_FULL, 'OGF%PP', CROUTINE)

   ENDDO

   IP = 1
   IP0 = IP
   PPF%ROW(1) = IP

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'START OF PPF LOOP'
WRITE(MSG%LU_DEBUG,*) '  GF%NC=',GF%NC
WRITE(MSG%LU_DEBUG,*) '  GC%NC=',GC%NC
WRITE(MSG%LU_DEBUG,*) '  GC%NCE2=',GC%NCE2
WRITE(MSG%LU_DEBUG,*) '  GC%NC_GALERKIN=',GC%NC_GALERKIN
WRITE(MSG%LU_DEBUG,*) '  GC%LOCAL_TO_GLOBAL=',GC%LOCAL_TO_GLOBAL
WRITE(MSG%LU_DEBUG,*) '  GC%CELLS_GLOBAL=',GC%CELLS_GLOBAL(1:GC%NC_GALERKIN)
#endif

   FINE_CELLS_LOOP: DO IC = 1, GF%NC

      ! TODO: Better time measurement!
      TNOW = CURRENT_TIME()
      IF (MYID == 0 .AND. MOD(IC,1000) == 0) WRITE(*,*) 'ScaRC-AMG-Setup: Processing cell ',IC,' of ',GF%NC

      INTERNAL_COARSE_CELLS_LOOP: DO ICC = 1, GC%NC
         CALL SCARC_MULTIPLY_POISSON_PROL(AF, PF, PPF, GC%LOCAL_TO_GLOBAL(ICC), ICC0, IP, IC)
      ENDDO INTERNAL_COARSE_CELLS_LOOP

      EXTERNAL_COARSE_CELLS_LOOP: DO ICC = 1, GC%NC_GALERKIN
         CALL SCARC_MULTIPLY_POISSON_PROL(AF, PF, PPF, GC%CELLS_GLOBAL(ICC), ICC0, IP, IC)
      ENDDO EXTERNAL_COARSE_CELLS_LOOP

      ! take care that at least one entry per fine cell is generated
      IF (IP == IP0) THEN
         PPF%VAL(IP)  = 0.0_EB
         PPF%COL(IP)  = ICC0
         PPF%COLG(IP) = ICC0
         IP = IP + 1
      ENDIF
      IP0 = IP

      PPF%ROW(IC+1) = IP

CPU(MYID)%AMG =CPU(MYID)%AMG+CURRENT_TIME()-TNOW
TSUM = TSUM + CPU(MYID)%AMG

   ENDDO FINE_CELLS_LOOP

#ifdef WITH_SCARC_DEBUG
   WRITE(MSG%LU_DEBUG,*) '--------------> NM=',NM
   CALL SCARC_DEBUG_CMATRIX (PPF, 'PPF-FINE','AFTER MULTIPLY')
#endif

ENDDO

! Exchange overlapping parts of Prolongation matrix and extract exchanged data
IF (NMESHES > 1) THEN
   CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_MATRIX_COLS,  NSCARC_MATRIX_POISSON_PROL, NL)
   CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_MATRIX_COLSG, NSCARC_MATRIX_POISSON_PROL, NL)
   CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_MATRIX_VALS,  NSCARC_MATRIX_POISSON_PROL, NL)
   CALL SCARC_EXTRACT_MATRIX_OVERLAPS(NSCARC_MATRIX_POISSON_PROL, 0, NL)
ENDIF

! Reduce workspace for Prolongation matrix to really needed size
DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_MULTIGRID (NM, NL, NL+1)
   PPF => SCARC_POINT_TO_CMATRIX (GF, NSCARC_MATRIX_POISSON_PROL)         

#ifdef WITH_SCARC_DEBUG
   WRITE(MSG%LU_DEBUG,*) '--------------> NM=',NM
   CALL SCARC_DEBUG_CMATRIX (PPF, 'AP-FINE','END SETUP_POISSON_PROL')
#endif

   CALL SCARC_REDUCE_CMATRIX (PPF, 'GF%POISSON-PROL', CROUTINE)
   DO INBR = 1, SCARC(NM)%N_NEIGHBORS
      NOM = SCARC(NM)%NEIGHBORS(INBR)
      CALL SCARC_POINT_TO_OTHER_GRID(NM, NOM, NL)
      OPP => SCARC_POINT_TO_OTHER_CMATRIX(OGF, NSCARC_MATRIX_POISSON_PROL)
      CALL SCARC_REDUCE_CMATRIX(OPP, 'OGF%POISSON_PROL', CROUTINE)
   ENDDO

ENDDO

END SUBROUTINE SCARC_SETUP_POISSON_PROL


! ------------------------------------------------------------------------------------------------------
!> \brief Find matching components to multiply row of Restriction matrix with column of Poisson-Prol matrix
! ------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_MULTIPLY_GALERKIN(PPF, RF, AC, ICCL, JCCL, JCCG, JCCG0, IP)
USE SCARC_POINTERS, ONLY: GF
TYPE (SCARC_CMATRIX_TYPE), POINTER, INTENT(IN) :: PPF, RF, AC
INTEGER, INTENT(IN) :: ICCL, JCCL, JCCG
INTEGER, INTENT(INOUT) :: IP, JCCG0
INTEGER :: IAPCOL, IRCOL, JC
REAL(EB) :: DSUM, TOL = 1E-12_EB

DSUM = 0.0_EB
DO IRCOL = RF%ROW(ICCL), RF%ROW(ICCL+1)-1
   JC = RF%COLG(IRCOL)
   IF (TYPE_COARSENING == NSCARC_COARSENING_GMG .AND. JC > GF%NC) RETURN
   IAPCOL = SCARC_FIND_MATCHING_COLUMN(PPF, JC, JCCG) 
   IF (IAPCOL > 0) THEN
      JCCG0 = JCCG
      DSUM = DSUM + RF%VAL(IRCOL) * PPF%VAL(IAPCOL)
   ENDIF
ENDDO

IF (ABS(DSUM) > TOL) THEN
   AC%COL(IP)  = JCCL
   AC%COLG(IP) = JCCG
   AC%VAL(IP) = DSUM
   IP = IP + 1
ENDIF

END SUBROUTINE SCARC_MULTIPLY_GALERKIN


! ------------------------------------------------------------------------------------------------------
!> \brief Setup Galerkin matrix on coarser grid level (AMG only)
! Note: Matrix POISPROL corresponds to POISSON x PROLONGATION  ~ AP
! ------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_GALERKIN(NL)
USE SCARC_POINTERS, ONLY: GF, GC, PPF, RF, AC, OAC, OGC
USE SCARC_POINTER_ROUTINES, ONLY: SCARC_POINT_TO_MULTIGRID, SCARC_POINT_TO_OTHER_MULTIGRID, &
                                  SCARC_POINT_TO_CMATRIX, SCARC_POINT_TO_OTHER_CMATRIX
INTEGER, INTENT(IN) :: NL
INTEGER  :: NM, IP, IP0, INBR, NOM, NLEN, ICCL, ICCG, JCC, JCCL, JCCG, JCCG0 = -1
REAL(EB) :: TNOW, TSUM = 0.0_EB

CROUTINE = 'SCARC_SETUP_GALERKIN'

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_MULTIGRID (NM, NL, NL+1)

   PPF => SCARC_POINT_TO_CMATRIX (GF, NSCARC_MATRIX_POISSON_PROL)         
   RF  => SCARC_POINT_TO_CMATRIX (GF, NSCARC_MATRIX_RESTRICTION)         
   AC  => SCARC_POINT_TO_CMATRIX (GC, NSCARC_MATRIX_POISSON)         

#ifdef WITH_SCARC_DEBUG
   WRITE(MSG%LU_DEBUG,*) '--------------> NM=',NM
   CALL SCARC_DEBUG_CMATRIX (RF, 'RESTRICTION-FINE', 'START OF SETUP_GALERKIN')
   WRITE(MSG%LU_DEBUG,*) '--------------> NM=',NM
   CALL SCARC_DEBUG_CMATRIX (PPF, 'PPF-FINE', 'START OF SETUP_GALERKIN')
#endif

   IF (.NOT.ALLOCATED (AC%VAL)) THEN
      AC%N_ROW = GC%NCE+1
      AC%N_VAL = AC%N_ROW**2             ! only temporarily TODO TOO BIG
      CALL SCARC_ALLOCATE_CMATRIX(AC, NL+1, NSCARC_PRECISION_DOUBLE, NSCARC_MATRIX_FULL, 'GC%POISSON', CROUTINE)
   ENDIF

   DO INBR = 1, SCARC(NM)%N_NEIGHBORS

      NOM = SCARC(NM)%NEIGHBORS(INBR)
      CALL SCARC_POINT_TO_OTHER_MULTIGRID(NM, NOM, NL, NL+1)

      OAC => SCARC_POINT_TO_OTHER_CMATRIX(OGC, NSCARC_MATRIX_POISSON)

      OAC%N_VAL = AC%N_VAL              ! TODO : CHECK : MUCH TOO BIG !!!
      OAC%N_ROW = GC%NCE2 + 1           ! TODO : CHECK : MUCH TOO BIG !!!
      CALL SCARC_ALLOCATE_CMATRIX(OAC, NL, NSCARC_PRECISION_DOUBLE, NSCARC_MATRIX_FULL, 'OAC%POISSON', CROUTINE)

   ENDDO

   CALL SCARC_DEALLOCATE_INT1(GC%CELLS_LOCAL,  'GC%CELLS_LOCAL',  CROUTINE)
   CALL SCARC_DEALLOCATE_INT1(GC%CELLS_GLOBAL, 'GC%CELLS_GLOBAL', CROUTINE)

   NLEN = 4 * (GC%NCE2 - GC%NC + 2)
#ifdef WITH_SCARC_DEBUG
   WRITE(MSG%LU_DEBUG,*) 'NL=',NL,': GC%NC=',GC%NC, ': GC%NCE=',GC%NCE,': GC%NCE2=',GC%NCE2,': NLEN=',NLEN
#endif

   CALL SCARC_ALLOCATE_INT1 (GC%CELLS_LOCAL, 1, NLEN, NSCARC_INIT_ZERO, 'GC%CELLS_LOCAL', CROUTINE)
   CALL SCARC_ALLOCATE_INT1 (GC%CELLS_GLOBAL, 1, NLEN, NSCARC_INIT_ZERO, 'GC%CELLS_GLOBAL', CROUTINE)

   CALL SCARC_GET_CELL_DEPENDENCIES_GALERKIN(GC, GF, PPF, NLEN)

   CALL SCARC_REDUCE_INT1(GC%CELLS_LOCAL, GC%NC_GALERKIN, 'GC%CELLS_LOCAL', CROUTINE)
   CALL SCARC_REDUCE_INT1(GC%CELLS_GLOBAL, GC%NC_GALERKIN, 'GC%CELLS_GLOBAL', CROUTINE)

ENDDO

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_MULTIGRID (NM, NL, NL+1)

   PPF => SCARC_POINT_TO_CMATRIX (GF, NSCARC_MATRIX_POISSON_PROL)         
   RF  => SCARC_POINT_TO_CMATRIX (GF, NSCARC_MATRIX_RESTRICTION)         
   AC  => SCARC_POINT_TO_CMATRIX (GC, NSCARC_MATRIX_POISSON)         

   IP = 1
   IP0 = IP
   AC%ROW(1) = IP

   LOCAL_COARSE_CELLS_LOOP: DO ICCL = 1, GC%NC

TNOW = CURRENT_TIME()
      ICCG = GC%LOCAL_TO_GLOBAL(ICCL)               ! corresponding global coarse cell

      INTERNAL_COARSE_CELLS_LOOP: DO JCCL = 1, GC%NC
         JCCG = GC%LOCAL_TO_GLOBAL(JCCL)
         CALL SCARC_MULTIPLY_GALERKIN(PPF, RF, AC, ICCL, JCCL, JCCG, JCCG0, IP)
      ENDDO INTERNAL_COARSE_CELLS_LOOP

      !EXTERNAL_COARSE_CELLS_LOOP: DO JCCL = GC%NC + 1, GC%NCE2
      EXTERNAL_COARSE_CELLS_LOOP: DO JCC = 1, GC%NC_GALERKIN
         JCCL = GC%CELLS_LOCAL(JCC)
         JCCG = GC%CELLS_GLOBAL(JCC)
         CALL SCARC_MULTIPLY_GALERKIN(PPF, RF, AC, ICCL, JCCL, JCCG, JCCG0, IP)
      ENDDO EXTERNAL_COARSE_CELLS_LOOP

      ! take care that at least one entry per fine cell is generated
      IF (IP == IP0) THEN
         AC%VAL(IP)  = 0.0_EB
         AC%COL(IP)  = JCCG0
         AC%COLG(IP) = JCCG0
         IP = IP + 1
      ENDIF
      IP0 = IP

      AC%ROW(ICCL + 1) = IP

! TODO: better time measurement
CPU(MYID)%AMG =CPU(MYID)%AMG+CURRENT_TIME()-TNOW
TSUM = TSUM + CPU(MYID)%AMG

   ENDDO LOCAL_COARSE_CELLS_LOOP
   AC%N_ROW = ICCL 

#ifdef WITH_SCARC_DEBUG
   WRITE(MSG%LU_DEBUG,*) '--------------> NM=',NM
   CALL SCARC_DEBUG_CMATRIX (AC, 'POISSON-COARSE','END OF RAP')
#endif

   CALL SCARC_REDUCE_CMATRIX (AC, 'POISSON-COARSE', CROUTINE)
   CALL SCARC_GET_MATRIX_STENCIL_MAX(AC, GC%NC)

ENDDO

CALL SCARC_RESORT_MATRIX_ROWS(NL+1)             

MESH_INT = 0                            
RANK_INT = 0

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   CALL SCARC_POINT_TO_MULTIGRID (NM, NL, NL+1)
   AC => SCARC_POINT_TO_CMATRIX (GC, NSCARC_MATRIX_POISSON)         
   MESH_INT(NM) = AC%N_STENCIL_MAX
   RANK_INT = MAX(RANK_INT, MESH_INT(NM))
ENDDO

IF (N_MPI_PROCESSES>1) &
   CALL MPI_ALLREDUCE(MPI_IN_PLACE, RANK_INT, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, IERROR)

#ifdef WITH_MKL
IF (TYPE_MKL(NL+1) == NSCARC_MKL_LOCAL .OR. TYPE_MKL(NL+1) == NSCARC_MKL_GLOBAL) THEN
   DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
      CALL SCARC_SETUP_POISSON_MKL(NM, NL+1)
   ENDDO
ENDIF
#endif

END SUBROUTINE SCARC_SETUP_GALERKIN


! ------------------------------------------------------------------------------------------------------
!> \brief Compute entry of Poisson times Prolongation matrix at specified position
! This consists of a summation over the entries:    P(:,ICC)*A(IC,:) 
! Thus, it must be checked, if - for a given entry of A in row IC - the Prolongation matrix
! has a corresponding non-zero value
! ------------------------------------------------------------------------------------------------------
DOUBLE PRECISION FUNCTION SCARC_VALUE_RAP(IC, ICC)
USE SCARC_POINTERS, ONLY : GC, AF, PF
INTEGER, INTENT(IN) :: IC, ICC
INTEGER :: JC, IA, IP, JCC
REAL(EB) :: DSUM

DSUM = 0.0_EB
DO IA = AF%ROW(IC), AF%ROW(IC+1) - 1
   JC = AF%COL(IA)
   IF (JC < 0) CYCLE
   DO IP = PF%ROW(JC), PF%ROW(JC+1) -1
      JCC = PF%COL(IP) 
      IF (JCC == GC%LOCAL_TO_GLOBAL(ICC)) THEN
         DSUM = DSUM + AF%VAL(IA)*PF%VAL(IP)
         CYCLE
      ENDIF
   ENDDO
ENDDO
SCARC_VALUE_RAP = DSUM

#ifdef WITH_SCARC_DEBUG2
WRITE(MSG%LU_DEBUG,*) 'VALUE_AP: RETURN DSUM=', DSUM
#endif

END FUNCTION SCARC_VALUE_RAP


! ------------------------------------------------------------------------------------------------------
!> \brief Resort matrix entries such that diagonal entry comes first (compact storage technique only)
! ------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_RESORT_MATRIX_ROWS(NL)
USE SCARC_POINTERS, ONLY: G, A
USE SCARC_POINTER_ROUTINES, ONLY: SCARC_POINT_TO_GRID, SCARC_POINT_TO_CMATRIX
INTEGER, INTENT(IN) :: NL
INTEGER, ALLOCATABLE, DIMENSION(:) :: COL_AUX, COLG_AUX
REAL(EB), ALLOCATABLE, DIMENSION(:) :: VAL_AUX
INTEGER:: NM, NCOL, ICOL, JCOL, KCOL, IC
LOGICAL :: COLG_IS_DEFINED = .FALSE.

CROUTINE = 'SCARC_RESORT_MATRIX_ROWS'

! TODO: use correct length of COL
DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   CALL SCARC_POINT_TO_GRID (NM, NL)
   A => SCARC_POINT_TO_CMATRIX (G, NSCARC_MATRIX_POISSON)

   CALL SCARC_ALLOCATE_INT1  (COL_AUX,  1, A%N_STENCIL_MAX, NSCARC_INIT_ZERO, 'COL_AUX',  CROUTINE)
   CALL SCARC_ALLOCATE_INT1  (COLG_AUX, 1, A%N_STENCIL_MAX, NSCARC_INIT_ZERO, 'COLG_AUX', CROUTINE)
   CALL SCARC_ALLOCATE_REAL1 (VAL_AUX,  1, A%N_STENCIL_MAX, NSCARC_INIT_ZERO, 'VAL_AUX',  CROUTINE)

   IF (ALLOCATED(A%COLG)) COLG_IS_DEFINED = .TRUE.

   !DO IC = 1, G%NCE2
   DO IC = 1, G%NC
      COL_AUX = 0
      JCOL = 1
      DO ICOL = A%ROW(IC), A%ROW(IC+1)-1
         COL_AUX(JCOL) = A%COL(ICOL)
         IF (COLG_IS_DEFINED)  COLG_AUX(JCOL) = A%COLG(ICOL)
         VAL_AUX(JCOL) = A%VAL(ICOL)
         JCOL = JCOL + 1
      ENDDO
      NCOL = JCOL - 1

      ! Find column index of diagonal element
      JCOL = 0
      DO WHILE (JCOL <= NCOL)
        JCOL = JCOL + 1
        IF (COL_AUX(JCOL) == IC) EXIT
      ENDDO

      ! Store corresponding index and value in first matrix element of that row
      ICOL = A%ROW(IC)
      A%COL(ICOL)  = COL_AUX(JCOL)
      A%COLG(ICOL) = COLG_AUX(JCOL)
      A%VAL(ICOL)  = VAL_AUX(JCOL)

      COL_AUX(JCOL) = 99999999
      COLG_AUX(JCOL) = 99999999

      IF (COLG_IS_DEFINED) THEN
         JCOL = MINLOC(COLG_AUX(1:NCOL), DIM=1)
      ELSE
         JCOL = MINLOC(COL_AUX(1:NCOL), DIM=1)
      ENDIF
      KCOL = 1
      ICOL = ICOL + 1
      DO WHILE (KCOL < NCOL)
         A%COLG(ICOL) = COLG_AUX(JCOL)
         A%COL(ICOL)  = COL_AUX(JCOL)
         A%VAL(ICOL)  = VAL_AUX(JCOL)
         IF (COLG_IS_DEFINED) THEN
            COLG_AUX(JCOL) = 99999999
            JCOL = MINLOC(COLG_AUX(1:NCOL), DIM=1)
         ELSE
            COL_AUX(JCOL) = 99999999
            JCOL = MINLOC(COL_AUX(1:NCOL), DIM=1)
         ENDIF
         KCOL = KCOL + 1
         ICOL = ICOL + 1
      ENDDO

      IF (ICOL /= A%ROW(IC+1)) WRITE(*,*) 'ERROR IN RESORT_MATRIX_ROWS'

   ENDDO

#ifdef WITH_SCARC_DEBUG
   CALL SCARC_DEBUG_CMATRIX (A, 'A','AFTER RESORT_MATRIX_ROWS')
#endif

   CALL SCARC_DEALLOCATE_INT1  (COL_AUX,  'COL_AUX',  CROUTINE)
   CALL SCARC_DEALLOCATE_INT1  (COLG_AUX, 'COLG_AUX', CROUTINE)
   CALL SCARC_DEALLOCATE_REAL1 (VAL_AUX,  'VAL_AUX',  CROUTINE)

ENDDO

END SUBROUTINE SCARC_RESORT_MATRIX_ROWS



END MODULE SCRC


