MODULE SCARC_MULTIGRID_METHOD
  
USE GLOBAL_CONSTANTS
USE PRECISION_PARAMETERS, ONLY: EB, FB
USE MEMORY_FUNCTIONS, ONLY: CHKMEMERR
USE COMP_FUNCTIONS, ONLY: CURRENT_TIME
USE MPI
USE SCARC_CONSTANTS
USE SCARC_VARIABLES
USE SCARC_TYPES
USE SCARC_UTILITIES
USE SCARC_MESSAGE_SERVICES, ONLY: MSG
USE SCARC_TIME_MEASUREMENT, ONLY: CPU
USE SCARC_STACK_OPERATIONS
USE SCARC_ITERATION_MANAGER
USE SCARC_LINEAR_ALGEBRA
USE SCARC_ITERATION_MANAGER
USE SCARC_RELAXATION_METHOD

CONTAINS

SUBROUTINE SCARC_SETUP_MULTIGRID_ENVIRONMENT()

NSTACK = NSCARC_STACK_ROOT
STACK(NSTACK)%SOLVER => MAIN_GMG
CALL SCARC_SETUP_MULTIGRID(NSCARC_SOLVER_MAIN, NSCARC_SCOPE_GLOBAL, NSCARC_STAGE_ONE, NSTACK, NLEVEL_MIN, NLEVEL_MAX)

NSTACK = NSTACK + 1
SELECT CASE(TYPE_SMOOTH)

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
      CALL SCARC_SETUP_MJAC(NLEVEL_MIN, NLEVEL_MAX)

    ! GS-preconditioning in matrix form (acting locally by default)

   CASE (NSCARC_RELAX_MGS)
      STACK(NSTACK)%SOLVER => SMOOTH_MGS
      CALL SCARC_SETUP_SMOOTH(NSTACK, NSCARC_SCOPE_LOCAL)
      CALL SCARC_SETUP_MGS(NLEVEL_MIN, NLEVEL_MAX)

   ! SGS-preconditioning in matrix form (acting locally by default)

   CASE (NSCARC_RELAX_MSGS)
      STACK(NSTACK)%SOLVER => SMOOTH_MSGS
      CALL SCARC_SETUP_SMOOTH(NSTACK, NSCARC_SCOPE_LOCAL)
      CALL SCARC_SETUP_MSGS(NLEVEL_MIN, NLEVEL_MAX)

   ! SOR-preconditioning in matrix form (acting locally by default)

   CASE (NSCARC_RELAX_MSOR)
      STACK(NSTACK)%SOLVER => SMOOTH_MSOR
      CALL SCARC_SETUP_SMOOTH(NSTACK, NSCARC_SCOPE_LOCAL)
      CALL SCARC_SETUP_MSOR(NLEVEL_MIN, NLEVEL_MAX, NSTACK)

   ! SSOR-preconditioning in matrix form (acting locally by default)

   CASE (NSCARC_RELAX_MSSOR)
      STACK(NSTACK)%SOLVER => SMOOTH_MSSOR
      CALL SCARC_SETUP_SMOOTH(NSTACK, NSCARC_SCOPE_LOCAL)
      CALL SCARC_SETUP_MSSOR(NLEVEL_MIN, NLEVEL_MAX, NSTACK)

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
   ! Smoothing by LU-decomposition

   CASE (NSCARC_RELAX_MKL)
      CALL SCARC_SETUP_SMOOTH(NSTACK, TYPE_SCOPE(2))

      SELECT CASE(TYPE_SCOPE(2))

         ! Globally acting - call global CLUSTER_SPARSE_SOLVER on MKL

         CASE (NSCARC_SCOPE_GLOBAL)
            STACK(NSTACK)%SOLVER => SMOOTH_MKL
            CALL SCARC_SETUP_CLUSTER(NLEVEL_MIN, NLEVEL_MIN)

         ! Locally acting - call local PARDISO solvers based on MKL

         CASE (NSCARC_SCOPE_LOCAL)
            STACK(NSTACK)%SOLVER => SMOOTH_MKL
            CALL SCARC_SETUP_PARDISO(NLEVEL_MIN, NLEVEL_MIN)

      END SELECT
#endif

END SELECT

! Globally acting coarse grid solver

NSTACK = NSTACK + 1
CALL SCARC_SETUP_COARSE_SOLVER(NSCARC_STAGE_ONE, NSCARC_SCOPE_GLOBAL, NSTACK, NLEVEL_MAX, NLEVEL_MAX)

END SUBROUTINE SCARC_SETUP_MULTIGRID_ENVIRONMENT

 
! ----------------------------------------------------------------------------------------------------
!> \brief Allocate and initialize vectors for MKL-methods
! ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_COARSE_SOLVER(NSTAGE, NSCOPE, NSTACK, NLMIN, NLMAX)
INTEGER, INTENT(IN)    :: NSCOPE, NSTAGE, NLMIN, NLMAX
INTEGER, INTENT(INOUT) :: NSTACK

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'SETUP COARSE_SOLVER: START'
#endif
SELECT_COARSE: SELECT CASE (TYPE_COARSE)

   ! -------------- CG-method is used as iterative coarse grid solver
   CASE (NSCARC_COARSE_ITERATIVE)

      ! initialize current stack position as CG-method
      STACK(NSTACK)%SOLVER => COARSE_KRYLOV
      CALL SCARC_SETUP_KRYLOV(NSCARC_SOLVER_COARSE, NSCOPE, NSTAGE, NSTACK, NLMIN, NLMAX)

      ! and next stack position as its SSOR-preconditioner
      NSTACK = NSTACK + 1
      TYPE_PRECON = NSCARC_RELAX_SSOR
      STACK(NSTACK)%SOLVER => PRECON_SSOR
      CALL SCARC_SETUP_PRECON(NSTACK, NSCOPE)

   ! -------------- LU-decomposition (from MKL) is used as direct coarse grid solver
#ifdef WITH_MKL 
   CASE (NSCARC_COARSE_DIRECT)

      ! Global scope in the multi-mesh case:
      ! initialize current stack position as global CLUSTER_SPARSE_SOLVER
      !IF (NSCOPE == NSCARC_SCOPE_GLOBAL .AND. NMESHES > 1) THEN
      IF (NMESHES > 1) THEN
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'SETUP COARSE_SOLVER: CLUSTER'
#endif
         STACK(NSTACK)%SOLVER => COARSE_CLUSTER
         CALL SCARC_SETUP_MKL(NSCARC_SOLVER_COARSE, NSCOPE, NSTAGE, NSTACK, NLMIN, NLMAX)
         CALL SCARC_SETUP_CLUSTER(NLMIN, NLMAX)

      ! Local scope:
      ! initialize current stack position as PARDISO solver
      ELSE
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'SETUP COARSE_SOLVER: PARDISO'
#endif
         STACK(NSTACK)%SOLVER => COARSE_PARDISO
         CALL SCARC_SETUP_MKL(NSCARC_SOLVER_COARSE, NSCOPE, NSTAGE, NSTACK, NLMIN, NLMAX)
         CALL SCARC_SETUP_PARDISO(NLMIN, NLMAX)
      ENDIF
#endif

   ! -------------- Otherwise: print error message
   CASE DEFAULT
      CALL SCARC_SHUTDOWN(NSCARC_ERROR_PARSE_INPUT, SCARC_NONE, TYPE_COARSE)

END SELECT SELECT_COARSE
END SUBROUTINE SCARC_SETUP_COARSE_SOLVER

! ------------------------------------------------------------------------------------------------
!> \brief Perform geometric multigrid method based on global possion-matrix
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_METHOD_MULTIGRID(NSTACK, NPARENT, NRHS, NLEVEL)
INTEGER, INTENT(IN) :: NSTACK, NPARENT, NRHS, NLEVEL
INTEGER :: NS, NP, NL
INTEGER :: NSTATE, ICYCLE
REAL (EB) :: TNOW, TNOW_COARSE

TNOW = CURRENT_TIME()
ITE_MG = 0

! Store current and parent stack position and current level

TYPE_MATVEC = NSCARC_MATVEC_GLOBAL
NS = NSTACK
NP = NPARENT
NL = NLEVEL

 
! ---------- Initialization:
!   - Save SETTING (in case that subsequent solvers with different SETTING are called)
!   - Define parameters for current scope (note: NL denotes the finest level)
!   - Initialize solution, right hand side vector
  
CALL SCARC_SETUP_SOLVER(NS, NP)
CALL SCARC_SETUP_WORKSPACE(NS, NL, NRHS)

  
! ---------- Compute initial defect:  
!            RESIN := || B - A*X ||
!   - Initialize cycle counts for MG-iteration
!   - Perform initial matrix-vector product on finest level
!   - calculate norm of initial residual on finest level
  
#ifdef WITH_SCARC_DEBUG
!CALL SCARC_PRESET_VECTOR(B, NL)
CALL SCARC_DEBUG_LEVEL (X, 'MG INIT: X', NL)
CALL SCARC_DEBUG_LEVEL (B, 'MG INIT: B', NL)
#endif

CALL SCARC_MATVEC_PRODUCT (X, V, NL)                                  !  V := A*X
CALL SCARC_VECTOR_SUM (B, V, 1.0_EB, -1.0_EB, NL)                     !  V := B - V

#ifdef WITH_SCARC_DEBUG
CALL SCARC_DEBUG_LEVEL (V, 'MG INIT: V', NL)
#endif

RES    = SCARC_L2NORM (V, NL)                                         !  RESIN := ||V||
RESIN  = RES
NSTATE = SCARC_CONVERGENCE_STATE (0, NS, NL)                          !  RES < TOL already ??

IF (TYPE_SOLVER == NSCARC_SOLVER_PRECON .AND. RESIN <= 1E-6_EB) THEN
   CALL SCARC_VECTOR_SUM (V, X, 1.0_EB, 1.0_EB, NL)                    !  x := omega * v + x
   CALL SCARC_UPDATE_PRECONDITIONER(NLEVEL_MIN)
   CALL SCARC_RELEASE_SOLVER(NS, NP)
   RETURN
ENDIF

ICYCLE = SCARC_CYCLING_CONTROL(NSCARC_CYCLING_SETUP, NL)

  
! ---------- Perform multigrid-looping (start each iteration on finest level)
  
MULTIGRID_LOOP: DO ITE = 1, NIT

   CALL SCARC_INCREASE_ITERATION_COUNTS(ITE)

   NL = NLEVEL_MIN
   ICYCLE = SCARC_CYCLING_CONTROL(NSCARC_CYCLING_RESET, NL)

   CYCLE_LOOP: DO WHILE (ICYCLE /= NSCARC_CYCLING_EXIT)

      ! Presmoothing  (smoothing/restriction till coarsest level is reached)
      ! initial and final residual are passed via vector V by default
 
      PRESMOOTHING_LOOP: DO WHILE (NL < NLEVEL_MAX)
         !IF (ITE /= 1) CALL SCARC_SMOOTHER (NSCARC_CYCLING_PRESMOOTH, NS+1, NS, NL)         ! D_fine   := Smooth(defect)
         CALL SCARC_SMOOTHER (NSCARC_CYCLING_PRESMOOTH, NS+1, NS, NL)         ! D_fine   := Smooth(defect)
         CALL SCARC_RESTRICTION (V, B, NL, NL+1)                              ! B_coarse := Rest(D_fine)
#ifdef WITH_SCARC_DEBUG
CALL SCARC_DEBUG_LEVEL (V, 'MG PRE: V', NL)
CALL SCARC_DEBUG_LEVEL (B, 'MG PRE: B', NL+1)
#endif
         CALL SCARC_VECTOR_CLEAR (X, NL+1)                                    ! use zero initial guess on coarse level
         NL = NL + 1                                                          ! set coarser level
      ENDDO PRESMOOTHING_LOOP

 
      ! Coarse grid solver
 
      TNOW_COARSE = CURRENT_TIME()
      CALL SCARC_METHOD_COARSE(NS+2, NS, NLEVEL_MAX)                          ! X_coarse := exact_sol(.)
      CPU(MYID)%COARSE =CPU(MYID)%COARSE+CURRENT_TIME()-TNOW_COARSE

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) '==================> AFTER SCARC_METHOD_COARSE'
CALL SCARC_DEBUG_LEVEL (X, 'MG COA: X', NLEVEL_MAX)
CALL SCARC_DEBUG_LEVEL (B, 'MG COA: B', NLEVEL_MAX)
#endif
 
      ! Postsmoothing (smoothing/restriction till finest level is reached again)
 
      POSTSMOOTHING_LOOP: DO WHILE (NL > NLEVEL_MIN)
         NL=NL-1
         CALL SCARC_PROLONGATION (X, V, NL+1, NL)                             ! V_fine := Prol(X_coarse)
         CALL SCARC_VECTOR_SUM (V, X, 1.0_EB, 1.0_EB, NL)                     ! X_fine := V_fine + X_fine

#ifdef WITH_SCARC_DEBUG
CALL SCARC_DEBUG_LEVEL (V, 'MG after PROL: V', NL)
CALL SCARC_DEBUG_LEVEL (X, 'MG new X', NL)
CALL SCARC_DEBUG_LEVEL (B, 'MG before POST: B', NL)
#endif
         CALL SCARC_SMOOTHER (NSCARC_CYCLING_POSTSMOOTH, NS+1, NS, NL)        ! V_fine := Smooth(defect)
#ifdef WITH_SCARC_DEBUG
CALL SCARC_DEBUG_LEVEL (V, 'MG POST: V', NL)
CALL SCARC_DEBUG_LEVEL (X, 'MG POST: X', NL)
#endif
         ICYCLE = SCARC_CYCLING_CONTROL(NSCARC_CYCLING_PROCEED, NL)           ! perform requested cycle
         IF (ICYCLE /= NSCARC_CYCLING_POSTSMOOTH) CYCLE CYCLE_LOOP
      ENDDO POSTSMOOTHING_LOOP

   ENDDO CYCLE_LOOP

   IF (NL /= NLEVEL_MIN) CALL SCARC_SHUTDOWN(NSCARC_ERROR_MULTIGRID_LEVEL, SCARC_NONE, NL)

 
   ! Compute norm of new residual on finest level and  leave loop correspondingly
 
   CALL SCARC_MATVEC_PRODUCT (X, V, NL)                                       ! V := A*X
   CALL SCARC_VECTOR_SUM (B, V, 1.0_EB, -1.0_EB, NL)                          ! V := F - V

#ifdef WITH_SCARC_DEBUG
CALL SCARC_DEBUG_LEVEL (X, 'MG ITE X', NL)
CALL SCARC_DEBUG_LEVEL (V, 'MG RES V', NL)
#endif
   RES = SCARC_L2NORM (V, NL)                                                 ! RES := ||V||
   NSTATE = SCARC_CONVERGENCE_STATE(0, NS, NL)                                ! convergence ?
   IF (NSTATE /= NSCARC_STATE_PROCEED) EXIT MULTIGRID_LOOP

ENDDO MULTIGRID_LOOP

  
! ---------- Determine convergence rate and print corresponding information:
! In case of MG as main solver:
!   - Transfer ScaRC solution vector X to FDS pressure vector
!   - Set ghost cell values along external boundaries
!   - Exchange values along internal boundaries (consistency!)
  
CALL SCARC_CONVERGENCE_RATE(NSTATE, NS, NL)

#ifdef WITH_SCARC_DEBUG
CALL SCARC_DEBUG_LEVEL (X, 'MG method: FINAL X', NL)
#endif

SELECT CASE (TYPE_SOLVER)
   CASE (NSCARC_SOLVER_MAIN)
      CALL SCARC_UPDATE_MAINCELLS(NLEVEL_MIN)
      CALL SCARC_UPDATE_GHOSTCELLS(NLEVEL_MIN)
#ifdef WITH_SCARC_POSTPROCESSING
      CALL SCARC_PRESSURE_DIFFERENCE(NLEVEL_MIN)
#endif
   CASE (NSCARC_SOLVER_PRECON)
      CALL SCARC_UPDATE_PRECONDITIONER(NLEVEL_MIN)
END SELECT

CALL SCARC_RELEASE_SOLVER(NS, NP)

END SUBROUTINE SCARC_METHOD_MULTIGRID


! ------------------------------------------------------------------------------------------------
!> \brief Control multigrid cycling (F/V/W)
! Note: NLEVEL_MIN corresponds to finest level, NLEVEL_MAX to coarsest level
! ------------------------------------------------------------------------------------------------
INTEGER FUNCTION SCARC_CYCLING_CONTROL(NTYPE, NL)
USE SCARC_POINTERS, ONLY: MG
INTEGER, INTENT(IN) :: NTYPE, NL
INTEGER :: NM, IL, ICYCLE

SELECT CASE (NTYPE)

   ! Initialize cycle counts at beginning of multigrid method
 
   CASE (NSCARC_CYCLING_SETUP)
   
      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   
         MG => SCARC(NM)%LEVEL(NL)%MG
         MG%CYCLING(2)=1
   
         DO IL = NLEVEL_MIN+1, NLEVEL_MAX - 1
            MG => SCARC(NM)%LEVEL(IL)%MG
            IF (TYPE_CYCLING==NSCARC_CYCLING_F) THEN
               MG%CYCLING(2)=2
            ELSE
               MG%CYCLING(2)=TYPE_CYCLING
            ENDIF
         ENDDO
      ENDDO
   
      ICYCLE = NSCARC_CYCLING_NEXT
   
 
   ! Reset cycle counts at beginning of each new multigrid iteration
 
   CASE (NSCARC_CYCLING_RESET)
   
      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
         DO IL = NLEVEL_MIN, NLEVEL_MAX
            MG => SCARC(NM)%LEVEL(IL)%MG
            MG%CYCLING(1)=MG%CYCLING(2)
         ENDDO
      ENDDO
      ICYCLE = NSCARC_CYCLING_NEXT
   
 
   ! Determine where to proceed with cycling
 
   CASE (NSCARC_CYCLING_PROCEED)
   
      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   
         MG => SCARC(NM)%LEVEL(NL)%MG
         MG%CYCLING(1)=MG%CYCLING(1)-1
   
         IF (MG%CYCLING(1)==0) THEN
            IF (TYPE_CYCLING==NSCARC_CYCLING_F) THEN
               MG%CYCLING(1)=1
            ELSE
               MG%CYCLING(1)=MG%CYCLING(2)
            ENDIF
            IF (NL == NLEVEL_MIN) THEN
               ICYCLE = NSCARC_CYCLING_EXIT
            ELSE
               ICYCLE = NSCARC_CYCLING_POSTSMOOTH
            ENDIF
         ELSE
            IF (NL == NLEVEL_MIN) THEN
               ICYCLE = NSCARC_CYCLING_EXIT
            ELSE
               ICYCLE = NSCARC_CYCLING_NEXT
            ENDIF
         ENDIF
   
      ENDDO
   
END SELECT

SCARC_CYCLING_CONTROL = ICYCLE
RETURN

END FUNCTION SCARC_CYCLING_CONTROL

END MODULE SCARC_MULTIGRID_METHOD
