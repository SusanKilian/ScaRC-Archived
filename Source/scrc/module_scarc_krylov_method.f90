MODULE SCARC_KRYLOV_METHOD
  
USE GLOBAL_CONSTANTS
USE PRECISION_PARAMETERS, ONLY: EB, FB
USE MEMORY_FUNCTIONS, ONLY: CHKMEMERR
USE MPI
USE SCARC_CONSTANTS
USE SCARC_VARIABLES
USE SCARC_TYPES
USE SCARC_MESSAGE_SERVICES, ONLY: MSG
USE SCARC_TIME_MEASUREMENT, ONLY: CPU
USE SCARC_ERROR_MANAGER
USE SCARC_RELAXATION_SOLVERS


CONTAINS


! -----------------------------------------------------------------------------------------------------------------
!> \brief  Setup environment for Krylov methods
! -----------------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_KRYLOV_ENVIRONMENT()

NSTACK = NSCARC_STACK_ROOT
STACK(NSTACK)%SOLVER => MAIN_CG
CALL SCARC_SETUP_KRYLOV(NSCARC_SOLVER_MAIN, NSCARC_SCOPE_GLOBAL, NSCARC_STAGE_ONE, NSTACK, NLEVEL_MIN, NLEVEL_MIN)

 
! Setup preconditioner for Krylov solver
 
NSTACK = NSTACK + 1
SELECT_KRYLOV_PRECON: SELECT CASE (TYPE_PRECON)

   ! Jacobi-preconditioning (acting locally by default)

   CASE (NSCARC_RELAX_JAC)
      STACK(NSTACK)%SOLVER => PRECON_JAC
      CALL SCARC_SETUP_PRECON(NSTACK, NSCARC_SCOPE_LOCAL)

   ! SSOR-preconditioning (acting locally by default)

   CASE (NSCARC_RELAX_SSOR)
      STACK(NSTACK)%SOLVER => PRECON_SSOR
      CALL SCARC_SETUP_PRECON(NSTACK, NSCARC_SCOPE_LOCAL)

   ! JACOBI-preconditioning in matrix form (acting locally by default)

   CASE (NSCARC_RELAX_MJAC)
      STACK(NSTACK)%SOLVER => PRECON_MJAC
      CALL SCARC_SETUP_PRECON(NSTACK, NSCARC_SCOPE_LOCAL)
      CALL SCARC_SETUP_MJAC(NLEVEL_MIN, NLEVEL_MAX)

   ! GS-preconditioning in matrix form (acting locally by default)

   CASE (NSCARC_RELAX_MGS)
      STACK(NSTACK)%SOLVER => PRECON_MGS
      CALL SCARC_SETUP_PRECON(NSTACK, NSCARC_SCOPE_LOCAL)
      CALL SCARC_SETUP_MGS(NLEVEL_MIN, NLEVEL_MAX)

   ! SGS-preconditioning in matrix form (acting locally by default)

   CASE (NSCARC_RELAX_MSGS)
      STACK(NSTACK)%SOLVER => PRECON_MSGS
      CALL SCARC_SETUP_PRECON(NSTACK, NSCARC_SCOPE_LOCAL)
      CALL SCARC_SETUP_MSGS(NLEVEL_MIN, NLEVEL_MAX)

   ! SOR-preconditioning in matrix form (acting locally by default)

   CASE (NSCARC_RELAX_MSOR)
      STACK(NSTACK)%SOLVER => PRECON_MSOR
      CALL SCARC_SETUP_PRECON(NSTACK, NSCARC_SCOPE_LOCAL)
      CALL SCARC_SETUP_MSOR(NLEVEL_MIN, NLEVEL_MAX, NSTACK)

   ! SSOR-preconditioning in matrix form (acting locally by default)

   CASE (NSCARC_RELAX_MSSOR)
      STACK(NSTACK)%SOLVER => PRECON_MSSOR
      CALL SCARC_SETUP_PRECON(NSTACK, NSCARC_SCOPE_LOCAL)
      CALL SCARC_SETUP_MSSOR(NLEVEL_MIN, NLEVEL_MAX, NSTACK)

   ! LU-preconditioning in matrix form (acting locally by default)

   CASE (NSCARC_RELAX_LU)
      STACK(NSTACK)%SOLVER => PRECON_LU
      CALL SCARC_SETUP_PRECON(NSTACK, NSCARC_SCOPE_LOCAL)
      CALL SCARC_SETUP_LU(NLEVEL_MIN, NLEVEL_MAX)

   ! ILU(0)-preconditioning in matrix form (acting locally by default)

   CASE (NSCARC_RELAX_ILU)
      STACK(NSTACK)%SOLVER => PRECON_ILU
      CALL SCARC_SETUP_PRECON(NSTACK, NSCARC_SCOPE_LOCAL)
      CALL SCARC_SETUP_ILU(NLEVEL_MIN, NLEVEL_MAX)

   ! FFT-preconditioning (acting locally by default)

   CASE (NSCARC_RELAX_FFT)
      STACK(NSTACK)%SOLVER => PRECON_FFT
      CALL SCARC_SETUP_PRECON(NSTACK, NSCARC_SCOPE_LOCAL)
      CALL SCARC_SETUP_FFT(NLEVEL_MIN, NLEVEL_MIN)

   ! FFT-preconditioning (acting locally by default)

   CASE (NSCARC_RELAX_FFTO)
      STACK(NSTACK)%SOLVER => PRECON_FFT
      CALL SCARC_SETUP_PRECON(NSTACK, NSCARC_SCOPE_LOCAL)
      CALL SCARC_SETUP_FFTO(NLEVEL_MIN, NLEVEL_MIN)

#ifdef WITH_MKL
   ! LU-preconditioning based on MKL (either locally or globally acting depending on user specification)

   CASE (NSCARC_RELAX_MKL)
      STACK(NSTACK)%SOLVER => PRECON_MKL

      SELECT CASE(TYPE_SCOPE(1))

         ! Globally acting - call global CLUSTER_SPARSE_SOLVER from MKL

         CASE (NSCARC_SCOPE_GLOBAL)
            CALL SCARC_SETUP_PRECON(NSTACK, NSCARC_SCOPE_GLOBAL)
            CALL SCARC_SETUP_CLUSTER(NLEVEL_MIN, NLEVEL_MIN)

         ! locally acting - call global PARDISO solver from MKL

         CASE (NSCARC_SCOPE_LOCAL)
            CALL SCARC_SETUP_PRECON(NSTACK, NSCARC_SCOPE_LOCAL)
            CALL SCARC_SETUP_PARDISO(NLEVEL_MIN, NLEVEL_MIN)

      END SELECT
#endif

 
   ! Preconditioning by Geometric multigrid,
   ! either locally or Globally acting, depending on user specification stored in TYPE_SCOPE(1)
 
   CASE (NSCARC_RELAX_MULTIGRID)

      STACK(NSTACK)%SOLVER => PRECON_MG
      CALL SCARC_SETUP_PRECON(NSTACK, TYPE_SCOPE(1))
      CALL SCARC_SETUP_MULTIGRID(NSCARC_SOLVER_PRECON, TYPE_SCOPE(1), NSCARC_STAGE_TWO, NSTACK, &
                                 NLEVEL_MIN, NLEVEL_MAX)

      NSTACK = NSTACK + 1
      SELECT CASE (TYPE_SMOOTH)

         ! Jacobi-smoothing (acting locally by default)

         CASE (NSCARC_RELAX_JAC)
            STACK(NSTACK)%SOLVER => SMOOTH_JAC
            CALL SCARC_SETUP_SMOOTH(NSTACK, NSCARC_SCOPE_LOCAL)

         ! SSOR-smoothing (acting locally by default)

         CASE (NSCARC_RELAX_SSOR)
            STACK(NSTACK)%SOLVER => SMOOTH_SSOR
            CALL SCARC_SETUP_SMOOTH(NSTACK, NSCARC_SCOPE_LOCAL)

         ! Jacobi-preconditioning in matrix form (acting locally by default)

         CASE (NSCARC_RELAX_MJAC)
            STACK(NSTACK)%SOLVER => SMOOTH_MJAC
            CALL SCARC_SETUP_SMOOTH(NSTACK, NSCARC_SCOPE_LOCAL)
            CALL SCARC_SETUP_MJAC(NLEVEL_MIN, NLEVEL_MAX-1)

         ! GS-preconditioning in matrix form (acting locally by default)

         CASE (NSCARC_RELAX_MGS)
            STACK(NSTACK)%SOLVER => SMOOTH_MGS
            CALL SCARC_SETUP_SMOOTH(NSTACK, NSCARC_SCOPE_LOCAL)
            CALL SCARC_SETUP_MGS(NLEVEL_MIN, NLEVEL_MAX-1)

         ! SGS-preconditioning in matrix form (acting locally by default)

         CASE (NSCARC_RELAX_MSGS)
            STACK(NSTACK)%SOLVER => SMOOTH_MSGS
            CALL SCARC_SETUP_SMOOTH(NSTACK, NSCARC_SCOPE_LOCAL)
            CALL SCARC_SETUP_MSGS(NLEVEL_MIN, NLEVEL_MAX-1)

         ! SOR-preconditioning in matrix form (acting locally by default)

         CASE (NSCARC_RELAX_MSOR)
            STACK(NSTACK)%SOLVER => SMOOTH_MSOR
            CALL SCARC_SETUP_SMOOTH(NSTACK, NSCARC_SCOPE_LOCAL)
            CALL SCARC_SETUP_MSOR(NLEVEL_MIN, NLEVEL_MAX-1, NSTACK)

         ! SSOR-preconditioning in matrix form (acting locally by default)

         CASE (NSCARC_RELAX_MSSOR)
            STACK(NSTACK)%SOLVER => SMOOTH_MSSOR
            CALL SCARC_SETUP_SMOOTH(NSTACK, NSCARC_SCOPE_LOCAL)
            CALL SCARC_SETUP_MSSOR(NLEVEL_MIN, NLEVEL_MAX-1, NSTACK)

         ! FFT-smoothing (acting locally by default)

         CASE (NSCARC_RELAX_FFT)
            STACK(NSTACK)%SOLVER => SMOOTH_FFT
            CALL SCARC_SETUP_SMOOTH(NSTACK, NSCARC_SCOPE_LOCAL)
            CALL SCARC_SETUP_FFT(NLEVEL_MIN, NLEVEL_MAX-1)

         ! FFTO-smoothing (acting locally by default)

         CASE (NSCARC_RELAX_FFTO)
            STACK(NSTACK)%SOLVER => SMOOTH_FFTO
            CALL SCARC_SETUP_SMOOTH(NSTACK, NSCARC_SCOPE_LOCAL)
            CALL SCARC_SETUP_FFTO(NLEVEL_MIN, NLEVEL_MAX-1)
#ifdef WITH_MKL
         ! LU-smoothing (acting locally by default)

         CASE (NSCARC_RELAX_MKL)
            STACK(NSTACK)%SOLVER => SMOOTH_MKL
            CALL SCARC_SETUP_SMOOTH(NSTACK, NSCARC_SCOPE_LOCAL)
            CALL SCARC_SETUP_PARDISO(NLEVEL_MIN, NLEVEL_MIN)
#endif
      END SELECT

      ! Coarse grid solver (same scope of action as calling GMG)

      NSTACK = NSTACK + 1
      CALL SCARC_SETUP_COARSE_SOLVER(NSCARC_STAGE_TWO, TYPE_SCOPE(1), NSTACK, NLEVEL_MAX, NLEVEL_MAX)

END SELECT SELECT_KRYLOV_PRECON

 
! If two-level Krylov, allocate intermediate structures for interpolation and workspace for global coarse solver
 
IF (HAS_TWO_LEVELS) THEN

   IF (.NOT.IS_CG_AMG) CALL SCARC_SETUP_INTERPOLATION(NSCARC_STAGE_ONE, NLEVEL_MIN+1, NLEVEL_MAX)

   NSTACK = NSTACK + 1
   CALL SCARC_SETUP_COARSE_SOLVER(NSCARC_STAGE_ONE, NSCARC_SCOPE_GLOBAL, NSTACK, NLEVEL_MAX, NLEVEL_MAX)

ENDIF

END SUBROUTINE SCARC_SETUP_KRYLOV_ENVIRONMENT

! ----------------------------------------------------------------------------------------------------
!> \brief Allocate and initialize vectors for Krylov method
! ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_KRYLOV(NSOLVER, NSCOPE, NSTAGE, NSTACK, NLMIN, NLMAX)
USE SCARC_POINTERS, ONLY: SV
INTEGER, INTENT(IN) :: NSOLVER, NSCOPE, NSTAGE, NSTACK, NLMIN, NLMAX

! Basic setup of stack information and types for Krylov method
 
CALL SCARC_SETUP_STACK(NSTACK)

SV  => STACK(NSTACK)%SOLVER
SV%TYPE_METHOD   = NSCARC_METHOD_KRYLOV
SV%TYPE_SOLVER   = NSOLVER
SV%TYPE_SCOPE(0) = NSCOPE
SV%TYPE_STAGE    = NSTAGE
SV%TYPE_LEVEL(1) = NLMIN
SV%TYPE_LEVEL(2) = NLMAX
SV%TYPE_MATRIX   = TYPE_MATRIX

! Preset iteration parameters for Krylov method

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'SETUP_KRYLOV=', SCARC_KRYLOV_ITERATIONS, NSOLVER, NSCOPE, NSTAGE, NSTACK
#endif
SELECT CASE(NSOLVER)

   ! -------------- Krylov method is used as main solver
   CASE (NSCARC_SOLVER_MAIN)
   
      SV%CNAME = 'SCARC_MAIN_KRYLOV'
   
      SV%EPS = SCARC_KRYLOV_ACCURACY
      SV%NIT = SCARC_KRYLOV_ITERATIONS
   
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'SV%NIT=', SV%NIT,SCARC_KRYLOV_ITERATIONS
#endif
      SV%TYPE_RELAX    = TYPE_PRECON
      SV%TYPE_TWOLEVEL = TYPE_TWOLEVEL
   
   ! -------------- Krylov method is used as coarse grid solver solver
   CASE (NSCARC_SOLVER_COARSE)
   
      SV%CNAME = 'SCARC_COARSE_KRYLOV'
   
      SV%EPS = SCARC_COARSE_ACCURACY
      SV%NIT = SCARC_COARSE_ITERATIONS
   
      SV%TYPE_RELAX    = NSCARC_RELAX_SSOR             ! only use SSOR-preconditioning for coarse solver
      SV%TYPE_TWOLEVEL = NSCARC_TWOLEVEL_NONE          ! only use one level for coarse solver
   
   ! -------------- Otherwise: print error message
   CASE DEFAULT
   
      CALL SCARC_SHUTDOWN(NSCARC_ERROR_PARSE_INPUT, SCARC_NONE, NSOLVER)
   
END SELECT

 
! Point to solution vectors (in corresponding scope)
 
CALL SCARC_SETUP_REFERENCES(.TRUE.,.TRUE.,.TRUE.,.TRUE.,.TRUE.,.TRUE.,.TRUE., NSTACK)

END SUBROUTINE SCARC_SETUP_KRYLOV

! ------------------------------------------------------------------------------------------------
!> \brief Perform global conjugate gradient method based on global Possion-matrix
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_METHOD_KRYLOV(NSTACK, NPARENT, NRHS, NLEVEL)
USE SCARC_ITERATION_ENVIRONMENT
INTEGER, INTENT(IN) :: NSTACK, NPARENT, NRHS, NLEVEL
INTEGER :: NSTATE, NS, NP, NL, NG
REAL (EB) :: ALPHA, BETA, SIGMA, SIGMA0=0.0_EB
REAL (EB) :: TNOW, TNOWI

TNOW = CURRENT_TIME()
ITE_CG = 0

! Get current and parent stack position, and current level
TYPE_MATVEC = NSCARC_MATVEC_GLOBAL
NS  = NSTACK
NP  = NPARENT
NL  = NLEVEL
NG = TYPE_GRID


#ifdef WITH_SCARC_POSTPROCESSING
IF (ICYC == 1) THEN
   CALL SCARC_DUMP_SYSTEM(NS, NSCARC_DUMP_MESH)
   CALL SCARC_DUMP_SYSTEM(NS, NSCARC_DUMP_A)
ENDIF
#endif

! ---------- Initialization:
!   - Get parameters for current scope (note: NL denotes the finest level)
!   - Get right hand side vector and clear solution vectors

CALL SCARC_SETUP_SOLVER(NS, NP)
TYPE_GRID = NG

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) '====================== BEGIN KRYLOV METHOD '
WRITE(MSG%LU_DEBUG,*) 'NSTACK  =', NSTACK
WRITE(MSG%LU_DEBUG,*) 'NPARENT =', NPARENT
WRITE(MSG%LU_DEBUG,*) 'NRHS    =', NRHS
WRITE(MSG%LU_DEBUG,*) 'NLEVEL  =', NLEVEL
WRITE(MSG%LU_DEBUG,*) 'TYPE_GRID  =', TYPE_GRID
WRITE(MSG%LU_DEBUG,*) 'TYPE_MATVEC  =', TYPE_MATVEC
WRITE(MSG%LU_DEBUG,*) 'TYPE_PRECON  =', TYPE_PRECON
WRITE(MSG%LU_DEBUG,*) 'ITYPE  =', STACK(NS)%SOLVER%TYPE_RELAX
WRITE(MSG%LU_DEBUG,*) 'PRES_ON_WHOLE_DOMAIN =', PRES_ON_WHOLE_DOMAIN
WRITE(MSG%LU_DEBUG,*) 'IS_STRUCTURED = ', IS_STRUCTURED
WRITE(MSG%LU_DEBUG,*) 'IS_UNSTRUCTURED = ', IS_UNSTRUCTURED
WRITE(MSG%LU_DEBUG,*) 'IS_LAPLACE = ', IS_LAPLACE
WRITE(MSG%LU_DEBUG,*) 'IS_POISSON = ', IS_POISSON
#endif

CALL SCARC_SETUP_WORKSPACE(NS, NL, NRHS)

!CALL SCARC_PRESET_VECTOR(B, NL)

! In case of pure Neumann boundary conditions setup condensed system

IF (N_DIRIC_GLOBAL(NLEVEL_MIN) == 0) THEN
   CALL SCARC_VECTOR_INIT (X, 0.0_EB, NL)                    
   CALL SCARC_FILTER_MEANVALUE(B, NL)                       
   CALL SCARC_SETUP_SYSTEM_CONDENSED (B, NL, 1)            
ENDIF

#ifdef WITH_SCARC_POSTPROCESSING
CALL SCARC_DUMP_SYSTEM(NS, NSCARC_DUMP_B)
#endif

#ifdef WITH_SCARC_DEBUG
CALL SCARC_DEBUG_LEVEL (X, 'CG-METHOD: X INIT0 ', NL)
CALL SCARC_DEBUG_LEVEL (B, 'CG-METHOD: B INIT0 ', NL)
CALL SCARC_DEBUG_METHOD('BEGIN OF KRYLOV METHOD ',7)                     
#endif

! Compute initial residual 

IF (IS_MGM .AND. NSTACK == 3) THEN
   TYPE_MATVEC = NSCARC_MATVEC_LOCAL
ELSE
   TYPE_MATVEC = NSCARC_MATVEC_GLOBAL
ENDIF
CALL SCARC_MATVEC_PRODUCT (X, R, NL)                         !  r^0 := A*x^0
CALL SCARC_VECTOR_SUM     (B, R, -1.0_EB, 1.0_EB, NL)        !  r^0 := r^0 - b     corresponds to  A*x^0 - b

RES    = SCARC_L2NORM (R, NL)                                !  res   := ||r^0||
RESIN  = RES                                                 !  resin := res
NSTATE = SCARC_CONVERGENCE_STATE (0, NS, NL)                 !  res < tolerance ?

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'SUSI: KRYLOV: NSTACK:', NSTACK, ': TYPE_MATVEC=', TYPE_MATVEC
CALL SCARC_DEBUG_LEVEL (X, 'CG-METHOD: X INIT1 ', NL)
CALL SCARC_DEBUG_LEVEL (B, 'CG-METHOD: B INIT1 ', NL)
#endif

! Perform initial preconditioning

IF (NSTATE /= NSCARC_STATE_CONV_INITIAL) THEN                !  if no convergence yet, call intial preconditioner
   CALL SCARC_PRECONDITIONER(NS, NS, NL)                     !  v^0 := Precon(r^0)
   SIGMA0 = SCARC_SCALAR_PRODUCT(R, V, NL)                   !  SIGMA0 := (r^0,v^0)
#ifdef WITH_SCARC_DEBUG
CALL SCARC_DEBUG_LEVEL (R, 'CG-METHOD: R INIT1 ', NL)
CALL SCARC_DEBUG_LEVEL (V, 'CG-METHOD: V INIT1 ', NL)
WRITE(MSG%LU_DEBUG,*) 'SIGMA0=', SIGMA0
#endif
   CALL SCARC_VECTOR_COPY (V, D, -1.0_EB, NL)                !  d^0 := -v^0
ENDIF

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'RESIN, RES, ITE, SIGMA0:', RESIN, RES, ITE, SIGMA0
CALL SCARC_DEBUG_LEVEL (D, 'CG-METHOD: D INIT1 ', NL)
#endif


! ---------- Perform conjugate gradient looping

CG_LOOP: DO ITE = 1, NIT

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) '========================> CG : ITE =', ITE
#endif
   TNOWI = CURRENT_TIME()
   CALL SCARC_INCREASE_ITERATION_COUNTS(ITE)

   !TYPE_MATVEC = NSCARC_MATVEC_LOCAL
   CALL SCARC_MATVEC_PRODUCT (D, Y, NL)                      !  y^k := A*d^k

   ALPHA = SCARC_SCALAR_PRODUCT (D, Y, NL)                   !  alpha := (d^k,y^k)     corresponds to   (d^k,A*d^k)

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'ALPHA, SIGMA0=', ALPHA, SIGMA0
CALL SCARC_DEBUG_LEVEL (Y, 'CG-METHOD: Y AFTER MAT-VEC ', NL)
#endif

   ALPHA = SIGMA0/ALPHA                                      !  alpha := (r^k,v^k)/(d^k,A*d^k)

   CALL SCARC_VECTOR_SUM (D, X, ALPHA, 1.0_EB, NL)           !  x^{k+1} := x^k + alpha * d^k
   CALL SCARC_VECTOR_SUM (Y, R, ALPHA, 1.0_EB, NL)           !  r^{k+1} := r^k + alpha * y^k   ~  r^k + alpha * A * d^k

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'ITE, ITE_CG=', ITE, ITE_CG
CALL SCARC_DEBUG_LEVEL (X, 'CG-METHOD: X ITE ', NL)
CALL SCARC_DEBUG_LEVEL (Y, 'CG-METHOD: Y ITE ', NL)
CALL SCARC_DEBUG_LEVEL (R, 'CG-METHOD: R ITE ', NL)
WRITE(MSG%LU_DEBUG,*) '======================> CG : ITE2 =', ITE
#endif

   RES = SCARC_L2NORM (R, NL)                                !  res := ||r^{k+1}||
   NSTATE = SCARC_CONVERGENCE_STATE (0, NS, NL)              !  res < tolerance ??
   IF (NSTATE /= NSCARC_STATE_PROCEED) EXIT CG_LOOP

   CALL SCARC_PRECONDITIONER(NS, NS, NL)                     !  v^{k+1} := Precon(r^{k+1})

   SIGMA  = SCARC_SCALAR_PRODUCT (R, V, NL)                  !  sigma := (r^{k+1},v^{k+1})
   BETA   = SIGMA/SIGMA0                                     !  beta  := (r^{k+1},v^{k+1})/(r^k,v^k)
   SIGMA0 = SIGMA                                            !  save last sigma

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) '======================> CG : ITE3 =', ITE
CALL SCARC_DEBUG_LEVEL (V, 'CG-METHOD: V ITE ', NL)
CALL SCARC_DEBUG_LEVEL (D, 'CG-METHOD: D ITE ', NL)
#endif

   CALL SCARC_VECTOR_SUM (V, D, -1.0_EB, BETA, NL)           !  d^{k+1} := -v^{k+1} + beta * d^{k+1}

   CPU(MYID)%ITERATION=MAX(CPU(MYID)%ITERATION,CURRENT_TIME()-TNOWI)

ENDDO CG_LOOP

! ---------- Determine convergence rate and print corresponding information
! In case of CG as main solver:
!   - Transfer ScaRC solution vector X to FDS pressure vector
!   - Set ghost cell values along external boundaries
!   - Exchange values along internal boundaries

CALL SCARC_CONVERGENCE_RATE(NSTATE, NS, NL)

IF (N_DIRIC_GLOBAL(NLEVEL_MIN) == 0) THEN
   CALL SCARC_RESTORE_LAST_CELL(X, NL)
   CALL SCARC_FILTER_MEANVALUE(X, NL)
ENDIF

#ifdef WITH_SCARC_DEBUG
CALL SCARC_DEBUG_LEVEL (X, 'CG-METHOD: X FINAL', NL)
WRITE(MSG%LU_DEBUG,*) '=======================>> CG : END =', ITE
#endif

IF (TYPE_SOLVER == NSCARC_SOLVER_MAIN .AND. .NOT.IS_MGM) THEN
   CALL SCARC_UPDATE_MAINCELLS(NLEVEL_MIN)
   CALL SCARC_UPDATE_GHOSTCELLS(NLEVEL_MIN)
#ifdef WITH_SCARC_POSTPROCESSING
   CALL SCARC_PRESSURE_DIFFERENCE(NLEVEL_MIN)
   CALL SCARC_DUMP_SYSTEM(NS, NSCARC_DUMP_X)
#endif
ENDIF

#ifdef WITH_SCARC_DEBUG
CALL SCARC_DEBUG_METHOD('END OF KRYLOV METHOD ',6)                     
#endif

CALL SCARC_RELEASE_SOLVER(NS, NP)

END SUBROUTINE SCARC_METHOD_KRYLOV

END MODULE SCARC_KRYLOV_METHOD
