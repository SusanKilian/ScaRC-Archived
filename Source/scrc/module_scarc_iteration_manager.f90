! ================================================================================================================
!> \brief Iteration parameters and handles for single variants of ScaRC
! ================================================================================================================
MODULE SCARC_ITERATION_MANAGER

USE PRECISION_PARAMETERS, ONLY: EB
USE SCARC_VARIABLES
USE SCARC_MESSAGE_SERVICES
  
REAL(EB) :: DT                                  !< TS width 
REAL(EB) :: DTI                                 !< Inverse of TS width 
REAL(EB) :: OMEGA                               !< Relaxation parameter for current solver
REAL(EB) :: EPS                                 !< Requested accuracy for current solver
REAL(EB) :: RES                                 !< Current residual of current solver
REAL(EB) :: RESIN = -1.0_EB                     !< Initial residual of current solver
REAL(EB) :: CAPPA = -1.0_EB                     !< Convergence rate of current solver

REAL(EB) :: VELOCITY_ERROR_GLOBAL 

INTEGER :: NIT        = 0                       !< Maximum number of iterations in current solver
INTEGER :: ITE        = 0                       !< Current number of iterations in current solver
INTEGER :: ITE_CG     = 0                       !< Statistical information about number of Krylov iterations
INTEGER :: ITE_MG     = 0                       !< Statistical information about number of multigrid iterations
INTEGER :: ITE_LU     = 0                       !< Statistical information about number of LU iterations
INTEGER :: ITE_PRES   = 0                       !< Statistical information about number of pressure iterations
INTEGER :: ITE_TOTAL  = 0                       !< Statistical information about number of total iterations
INTEGER :: ITE_SMOOTH = 0                       !< Statistical information about number of smoothing iterations
INTEGER :: ITE_COARSE = 0                       !< Statistical information about number of coarse grid iterations
INTEGER :: ITE_GLOBAL = 0                       !< Statistical information about number of global iterations

INTEGER  :: X                                   !< Handle for solution 1D-vector 
INTEGER  :: B                                   !< Handle for right hand side one dimensional-vector 
INTEGER  :: D                                   !< Handle for auxiliary one-dimensional vector
INTEGER  :: R                                   !< Handle for auxiliary one-dimensional vector
INTEGER  :: V                                   !< Handle for auxiliary one-dimensional vector
INTEGER  :: Y                                   !< Handle for auxiliary one-dimensional vector
INTEGER  :: Z                                   !< Handle for auxiliary one-dimensional vector

#ifdef WITH_SCARC_DEBUG
INTEGER  :: E                                   !< Handle for one-dimensional error vector (debugging only)
#endif

TYPE SCARC_SOLVER_TYPE

   CHARACTER(30) :: CNAME = 'NONE'                             !< Name of current solver

   ! Types of different solver components
   INTEGER :: TYPE_ACCURACY      = NSCARC_ACCURACY_ABSOLUTE    !< Type of requested accuracy
   INTEGER :: TYPE_COARSE        = NSCARC_COARSE_DIRECT        !< Type of coarse grid solver for multilevel methods
   INTEGER :: TYPE_COARSENING    = NSCARC_COARSENING_CUBIC     !< Type of grid coarsening 
   INTEGER :: TYPE_CYCLING       = NSCARC_CYCLING_V            !< Type of cycling for multigrid method
   INTEGER :: TYPE_GRID          = NSCARC_GRID_STRUCTURED      !< Type of discretization
   INTEGER :: TYPE_EXCHANGE      = NSCARC_UNDEF_INT            !< Type of data exchange
   INTEGER :: TYPE_INTERPOL      = NSCARC_INTERPOL_CONSTANT    !< Type of interpolation method
   INTEGER :: TYPE_LEVEL(0:2)    = NSCARC_UNDEF_INT            !< Type of levels
   INTEGER :: TYPE_MATRIX        = NSCARC_MATRIX_COMPACT       !< Type of storage for matrix
   INTEGER :: TYPE_METHOD        = NSCARC_METHOD_KRYLOV        !< Type of ScaRC method
   INTEGER :: TYPE_MKL(0:10)     = NSCARC_UNDEF_INT            !< Type of MKL for single levels
   INTEGER :: TYPE_MKL_PRECISION = NSCARC_PRECISION_DOUBLE     !< Type of precision for MKL solver
   INTEGER :: TYPE_MULTIGRID     = NSCARC_MULTIGRID_GEOMETRIC  !< Type of multigrid method
   INTEGER :: TYPE_PARENT        = NSCARC_UNDEF_INT            !< Type of parent (calling) solver
   INTEGER :: TYPE_PRECON        = NSCARC_UNDEF_INT            !< Type of preconditioner for iterative solver
   INTEGER :: TYPE_RELAX         = NSCARC_UNDEF_INT            !< Type of preconditioner for iterative solver
   INTEGER :: TYPE_SCOPE(0:2)    = NSCARC_SCOPE_GLOBAL         !< Type of solver scopes
   INTEGER :: TYPE_SMOOTH        = NSCARC_UNDEF_INT            !< Type of smoother for multigrid method
   INTEGER :: TYPE_SOLVER        = NSCARC_SOLVER_MAIN          !< Type of surrounding solver stage
   INTEGER :: TYPE_STAGE         = NSCARC_STAGE_ONE            !< Type of surrounding solver stage
   INTEGER :: TYPE_STENCIL       = NSCARC_STENCIL_CONSTANT     !< Type of storage for matrix
   INTEGER :: TYPE_TWOLEVEL      = NSCARC_TWOLEVEL_NONE        !< Type of two-level method
   INTEGER :: TYPE_VECTOR        = NSCARC_UNDEF_INT            !< Type of vector to point to

   ! References to different vectors which are needed for the current solver
   INTEGER :: X = NSCARC_UNDEF_INT                             !< Reference to local X-vector, double precision
   INTEGER :: B = NSCARC_UNDEF_INT                             !< Reference to local B-vector, double precision
   INTEGER :: D = NSCARC_UNDEF_INT                             !< Reference to local D-vector, double precision
   INTEGER :: R = NSCARC_UNDEF_INT                             !< Reference to local R-vector, double precision
   INTEGER :: V = NSCARC_UNDEF_INT                             !< Reference to local V-vector, double precision
   INTEGER :: Y = NSCARC_UNDEF_INT                             !< Reference to local Y-vector, double precision
   INTEGER :: Z = NSCARC_UNDEF_INT                             !< Reference to local Z-vector, double precision

   INTEGER :: X_FB = NSCARC_UNDEF_INT                          !< Reference to local X-vector, single precision
   INTEGER :: B_FB = NSCARC_UNDEF_INT                          !< Reference to local B-vector, single precision
   INTEGER :: R_FB = NSCARC_UNDEF_INT                          !< Reference to local R-vector, single precision
   INTEGER :: V_FB = NSCARC_UNDEF_INT                          !< Reference to local V-vector, single precision

#ifdef WITH_SCARC_DEBUG
   INTEGER :: E = NSCARC_UNDEF_INT                             !< Reference to local E-vector, double precision
#endif

   ! Converegence requirements for current solver
   INTEGER  :: NIT   = NSCARC_UNDEF_INT                        !< Maximum iteration number
   INTEGER  :: ITE   = NSCARC_UNDEF_INT                        !< Current iteration number
   REAL(EB) :: EPS   = NSCARC_UNDEF_REAL_EB                    !< Required accuracy
   REAL(EB) :: RES   = NSCARC_UNDEF_REAL_EB                    !< Current residual
   REAL(EB) :: RESIN = NSCARC_UNDEF_REAL_EB                    !< Initial residual
   REAL(EB) :: ERR   = NSCARC_UNDEF_REAL_EB                    !< Initial residual
   REAL(EB) :: OMEGA = NSCARC_UNDEF_REAL_EB                    !< Relaxation parameter
   REAL(EB) :: CAPPA = NSCARC_UNDEF_REAL_EB                    !< Convergence rate

END TYPE SCARC_SOLVER_TYPE


CONTAINS

! ------------------------------------------------------------------------------------------------
!> \brief Set current iteration state
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SET_ITERATION_STATE (DT_CURRENT)
REAL(EB), INTENT(IN) :: DT_CURRENT

DT  = DT_CURRENT
DTI = 1.0_EB/DT_CURRENT

ITE_PRES = ITE_PRES + 1
ITE_GLOBAL = ICYC

END SUBROUTINE SCARC_SET_ITERATION_STATE


! ------------------------------------------------------------------------------------------------
!> \brief Check if solver converges or diverges and print out residual information
! ------------------------------------------------------------------------------------------------
INTEGER FUNCTION SCARC_CONVERGENCE_STATE(ISM, NS, NL)
INTEGER, INTENT(IN) :: NL, NS, ISM
INTEGER :: NSTATE

NSTATE = NSCARC_STATE_PROCEED

SELECT CASE (TYPE_ACCURACY)
   CASE (NSCARC_ACCURACY_RELATIVE)
      IF (RES <= RESIN*EPS .OR. RES <= NSCARC_THRESHOLD_CONVERGENCE) NSTATE = NSCARC_STATE_CONV
   CASE (NSCARC_ACCURACY_ABSOLUTE)
      IF (RES <= EPS .AND. RES <= RESIN) THEN
         IF (ITE == 0) THEN
            NSTATE = NSCARC_STATE_CONV_INITIAL
         ELSE
            NSTATE = NSCARC_STATE_CONV
         ENDIF
         NIT = 0
      ENDIF
END SELECT
IF (RES > NSCARC_THRESHOLD_DIVGERGENCE) NSTATE = NSCARC_STATE_DIVG

SCARC_CONVERGENCE_STATE = NSTATE

IF (HAS_CSV_DUMP) CALL SCARC_DUMP_CSV(ISM, NS, NL)

#ifdef WITH_SCARC_VERBOSE2
IF (TYPE_SOLVER == NSCARC_SOLVER_MAIN) &
   WRITE(MSG%LU_VERBOSE,1100) STACK(NS)%SOLVER%CNAME, NL, ITE, RES
1100 FORMAT (A30,': Level=',I4,': Iteration = ',I8,': Residual =',E14.6)
#endif

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG, 1000) STACK(NS)%SOLVER%CNAME, NL, ITE, RES
1000 FORMAT (A30,': Level=',I4,': Iteration = ',I8,': Residual =',e25.16)
#endif

END FUNCTION SCARC_CONVERGENCE_STATE


! ------------------------------------------------------------------------------------------------
!> \brief Compute convergence rate and print out residual information for final loop
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_CONVERGENCE_RATE(NSTATE, NS, NL)
INTEGER, INTENT(IN) :: NSTATE, NS, NL

IF (NSTATE == NSCARC_STATE_DIVG) THEN
   ITE   = - 1
   CAPPA = 1.0_EB
ELSE
   IF (NSTATE == NSCARC_STATE_CONV_INITIAL) THEN
      ITE= 0
   ELSE IF (NSTATE == NSCARC_STATE_CONV) THEN
      ITE= ITE
   ELSE
      ITE= ITE-1
   ENDIF
   IF (RESIN >= TWO_EPSILON_EB) THEN
      IF (ITE== 0) THEN
         CAPPA = 0.0_EB
      ELSE
         IF (NSTATE == NSCARC_STATE_CONV_INITIAL) THEN
            CAPPA = 0.0E0
         ELSE
            CAPPA = (RES/RESIN) ** (1.0_EB/ITE)
         ENDIF
      ENDIF
   ELSE
      CAPPA = 0.0_EB
   ENDIF
ENDIF

CALL SCARC_DUMP_CSV(0, NS, NL)

#ifdef WITH_SCARC_VERBOSE2
IF (TYPE_SOLVER == NSCARC_SOLVER_MAIN) WRITE(MSG%LU_VERBOSE,2000) STACK(NS)%SOLVER%CNAME, ITE, CAPPA
#endif

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,2000) STACK(NS)%SOLVER%CNAME, ITE, CAPPA
#endif

#if defined (WITH_SCARC_DEBUG) || defined (WITH_SCARC_VERBOSE)
2000 FORMAT (A30,': Iterations: ',i6,':   Convergence Rate =',E14.6,/)
#endif

END SUBROUTINE SCARC_CONVERGENCE_RATE


! ------------------------------------------------------------------------------------------------
!> \brief Increase corresponding iteration count (just for visualization of convergence behavior)
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_INCREASE_ITERATION_COUNTS(ITE0)
INTEGER, INTENT(IN) :: ITE0

SELECT CASE (TYPE_SOLVER)
   CASE (NSCARC_SOLVER_MAIN)
      SELECT CASE (TYPE_METHOD)
         CASE (NSCARC_METHOD_KRYLOV)
            ITE_CG = ITE0
         CASE (NSCARC_METHOD_MULTIGRID)
            ITE_MG = ITE0
         CASE (NSCARC_METHOD_LU)
            ITE_LU = ITE0
      END SELECT
   CASE (NSCARC_SOLVER_PRECON)
      ITE_MG = ITE0
   CASE (NSCARC_SOLVER_SMOOTH)
      ITE_SMOOTH = ITE0
   CASE (NSCARC_SOLVER_COARSE)
      ITE_COARSE = ITE0
END SELECT
ITE_TOTAL = ITE_TOTAL + 1

END SUBROUTINE SCARC_INCREASE_ITERATION_COUNTS

! ----------------------------------------------------------------------------------------------------
!> \brief Dump residual information
! ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_DUMP_CSV(ISM, NS, NL)
INTEGER, INTENT(IN) :: ISM, NS, NL

IF (.NOT.HAS_CSV_DUMP .OR. MYID /= 0) RETURN
IF (ITE_TOTAL == 0 .AND. TYPE_SOLVER /= NSCARC_SOLVER_MAIN) RETURN
IF (TYPE_SOLVER == NSCARC_SOLVER_COARSE) RETURN
WRITE(MSG%LU_STAT,1000) ITE_PRES, NS, ITE_TOTAL, ITE_CG, ITE_MG, NL, ITE_SMOOTH, ISM, ITE_COARSE, ITE_LU, RES, CAPPA

1000 FORMAT(10(I8,','), E14.6,',',E14.6)
END SUBROUTINE SCARC_DUMP_CSV


END MODULE SCARC_ITERATION_MANAGER

