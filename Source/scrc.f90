! ================================================================================================================
! Use of different directives possible
!  - WITH_MKL                : use MKL routines PARDISO, CLUSTER_SPARSE_SOLVER, DDOT, DAXPY, DAXPBY, DCOPY, DSCAL
!  - WITH_SCARC_VERBOSE      : print more detailed information about ScaRC iterations and workspace allocation
!  - WITH_SCARC_DEBUG        : print detaild debugging info (only for debugging purposes)
!  - WITH_SCARC_MGM          : experimental tests for McKeeney-Greengard-Mayo method
!  - WITH_SCARC_STANDALONE   : dump environment for ScaRC standalone version
! ================================================================================================================
!#undef WITH_MKL
#define WITH_SCARC_VERBOSE
#define WITH_SCARC_DEBUG
#undef WITH_SCARC_STANDALONE
#undef WITH_SCARC_MGM


! ================================================================================================================
! MODULE 'SCARC_GLOBAL_CONSTANTS': 
! Global SCARC parameters
! ================================================================================================================
MODULE SCARC_GLOBAL_CONSTANTS
USE PRECISION_PARAMETERS
IMPLICIT NONE

INTEGER :: ISUSI = 0

INTEGER, PARAMETER :: NSCARC_ACCURACY_ABSOLUTE       =  1, &    ! absolute accuracy must be reached
                      NSCARC_ACCURACY_RELATIVE       =  2       ! relative accuracy must be reached
                   
INTEGER, PARAMETER :: NSCARC_BUFFER_FULL             =  1, &    ! exchange buffer of type integer
                      NSCARC_BUFFER_LAYER1           =  2, &    ! exchange length of one layer
                      NSCARC_BUFFER_LAYER2           =  3, &    ! exchange length of two STENCIL
                      NSCARC_BUFFER_LAYER4           =  4, &    ! exchange length of four STENCIL
                      NSCARC_BUFFER_STENCIL          =  5, &    ! exchange length of stencil size STENCIL
                      NSCARC_BUFFER_BASIC            =  6       ! basic constant setup length

INTEGER, PARAMETER :: NSCARC_COARSE_ITERATIVE        =  1, &    ! iterative coarse grid solution
                      NSCARC_COARSE_DIRECT           =  2, &    ! direct coarse grid solution
                      NSCARC_COARSE_LEVEL            =  3, &    ! Coarse grid on specified grid level
                      NSCARC_COARSE_MACRO            =  4       ! Coarse grid on macro decomposition

INTEGER, PARAMETER :: NSCARC_COARSENING_GMG          =  1, &    ! GMG coarsening for multigrid method (AMG only)
                      NSCARC_COARSENING_SAMG         = 2       ! SAMG coarsening for multigrid method (AMG only)

INTEGER, PARAMETER :: NSCARC_COUPLING_MAX            = 10       ! maximum of possible couplings in stencil

INTEGER, PARAMETER :: NSCARC_CYCLING_F               =  0, &    ! F-cycle for mg-method
                      NSCARC_CYCLING_V               =  1, &    ! V-cycle for mg-method
                      NSCARC_CYCLING_W               =  2, &    ! W-cycle for mg-method
                      NSCARC_CYCLING_FMG             =  3, &    ! Full multigrid grid cycle for mg-method
                      NSCARC_CYCLING_SETUP           =  4, &    ! initialize cycle
                      NSCARC_CYCLING_RESET           =  5, &    ! reset cycle counts
                      NSCARC_CYCLING_PROCEED         =  6, &    ! proceed cycle counts
                      NSCARC_CYCLING_NEXT            =  7, &    ! perform next cycling loop
                      NSCARC_CYCLING_EXIT            =  8, &    ! exit cycling loop
                      NSCARC_CYCLING_PRESMOOTH       = -1, &    ! presmoothing cycle
                      NSCARC_CYCLING_POSTSMOOTH      =  1       ! postsmoothing cycle

INTEGER, PARAMETER :: NSCARC_DATA_INT                =  1, &    ! data type integer
                      NSCARC_DATA_REAL               =  2       ! data type real

INTEGER, PARAMETER :: NSCARC_DEBUG_BDRY              =  1, &    ! show pressure_bc_index
                      NSCARC_DEBUG_CELL_TYPE         =  2, &    ! show type of cells 
                      NSCARC_DEBUG_FACE              =  3, &    ! show face information
                      NSCARC_DEBUG_GRID              =  4, &    ! show grid information
                      NSCARC_DEBUG_CMATRIX            =  5, &    ! show matrix
                      NSCARC_DEBUG_MATRIX_SYM        =  6, &    ! show symmetric matrix
                      NSCARC_DEBUG_MEASURE           =  7, &    ! show measure of cells
                      NSCARC_DEBUG_RESTRICTION       =  8, &    ! show pressure quantities
                      NSCARC_DEBUG_PRESSURE          =  9, &    ! show pressure quantities
                      NSCARC_DEBUG_PROLONGATION      = 10, &    ! show pressure quantities
                      NSCARC_DEBUG_STACK             = 11, &    ! show matrix
                      NSCARC_DEBUG_CONNECTION        = 12, &    ! show strength of connection matrix
                      NSCARC_DEBUG_WALL              = 13       ! show wall information

#ifdef WITH_SCARC_STANDALONE
INTEGER, PARAMETER :: NSCARC_DUMP_A                  =  1, &    ! dump matrix A
                      NSCARC_DUMP_B                  =  2, &    ! dump right hand side B
                      NSCARC_DUMP_X                  =  3, &    ! dump solution X
                      NSCARC_DUMP_MESH               =  4       ! dump solution mesh information
#endif

INTEGER, PARAMETER :: NSCARC_ERROR_AGGREGATION       =  1, &    ! SAMG-aggregation failed
                      NSCARC_ERROR_BICG_DISABLED     =  2, &    ! BICG temporarily disabled
                      NSCARC_ERROR_BOUNDARY_SUM      =  3, &    ! wrong sum of elements along boundary
                      NSCARC_ERROR_BOUNDARY_TYPE     =  4, &    ! wrong boundary type
                      NSCARC_ERROR_DIRECT_NOMKL      =  5, &    ! wrong boundary type
                      NSCARC_ERROR_EXCHANGE_RECV     =  6, &    ! wrong receive exchange structure
                      NSCARC_ERROR_EXCHANGE_SEND     =  7, &    ! wrong send exchange structure
                      NSCARC_ERROR_FFT_GRID          =  8, &    ! wrong unstructured discretization for FFT
                      NSCARC_ERROR_GRID_INDEX        =  9, &    ! error with grid index
                      NSCARC_ERROR_GRID_NUMBER       = 10, &    ! cell number not divisable by two
                      NSCARC_ERROR_GRID_NUMBERX      = 11, &    ! cell number not divisable by two
                      NSCARC_ERROR_GRID_NUMBERY      = 12, &    ! cell number not divisable by two
                      NSCARC_ERROR_GRID_NUMBERZ      = 13, &    ! cell number not divisable by two
                      NSCARC_ERROR_GRID_RESOLUTION   = 14, &    ! error with grid resolution
                      NSCARC_ERROR_NEIGHBOR_NUMBER   = 15, &    ! wrong neighbor number
                      NSCARC_ERROR_NEIGHBOR_TYPE     = 16, &    ! wrong neighbor type
                      NSCARC_ERROR_MATRIX_ALLOCATION = 17, &    ! error in matrix allocation
                      NSCARC_ERROR_MATRIX_COPY       = 18, &    ! subdiagonal missing
                      NSCARC_ERROR_MATRIX_SETUP      = 19, &    ! error in matrix setup
                      NSCARC_ERROR_MATRIX_SIZE       = 20, &    ! error in matrix size
                      NSCARC_ERROR_MATRIX_SUBDIAG    = 21, &    ! subdiagonal missing
                      NSCARC_ERROR_MATRIX_SYMMETRY   = 22, &    ! matrix not symmetric
                      NSCARC_ERROR_MKL_CLUSTER       = 23, &    ! cluster_sparse_solver not available
                      NSCARC_ERROR_MKL_INTERNAL      = 24, &    ! pardiso solver not available
                      NSCARC_ERROR_MKL_PARDISO       = 25, &    ! pardiso solver not available
                      NSCARC_ERROR_MKL_STORAGE       = 26, &    ! wrong storage scheme for MKL-solvers
                      NSCARC_ERROR_MULTIGRID_LEVEL   = 27, &    ! wrong multigrid level
                      NSCARC_ERROR_PARSE_INPUT       = 28, &    ! wrong input parameter
                      NSCARC_ERROR_POC_STOP          = 29, &    ! temporary stop, only proof of concept
                      NSCARC_ERROR_STACK_MESSAGE     = 30, &    ! error with stack message
                      NSCARC_ERROR_STACK_SOLVER      = 31, &    ! error in solver stack
                      NSCARC_ERROR_STENCIL           = 32, &    ! error in matrix stencil
                      NSCARC_ERROR_RHS_SETUP         = 33, &    ! error with rhs setup
                      NSCARC_ERROR_VECTOR_LENGTH     = 34       ! error in vector length

INTEGER, PARAMETER :: NSCARC_EXCHANGE_BASIC_SIZES    =  1, &    ! exchange neighboring cell numbers on overlap
                      NSCARC_EXCHANGE_CELL_NEIGHBORS =  2, &    ! exchange initial information
                      NSCARC_EXCHANGE_CELL_NUMBERS   =  3, &    ! exchange neighboring grid widths on overlap
                      NSCARC_EXCHANGE_CELL_SIZES     =  4, &    ! exchange neighboring grid widths on overlap
                      NSCARC_EXCHANGE_AUX            = 19, &    ! exchange sum of nullspace entries
                      NSCARC_EXCHANGE_NULLSPACE      =  8, &    ! exchange sum of nullspace entries
                      NSCARC_EXCHANGE_POISSON_COLS   =  9, &    ! exchange columns of neighboring matrix on overlap
                      NSCARC_EXCHANGE_POISSON_COLSG  = 10, &    ! exchange columns of neighboring matrix on overlap
                      NSCARC_EXCHANGE_POISSON_DIAGS  = 11, &    ! exchange size of neighboring matrix 
                      NSCARC_EXCHANGE_POISSON_SIZES  = 12, &    ! exchange size of neighboring matrix 
                      NSCARC_EXCHANGE_POISSON_VALS   = 13, &    ! exchange values of neighboring matrix on overlap
                      NSCARC_EXCHANGE_PRESSURE       = 14, &    ! exchange vector values along internal boundaries
                      NSCARC_EXCHANGE_VECTOR_MEAN    = 15, &    ! exchange vector and build mean values with own data
                      NSCARC_EXCHANGE_VECTOR_PLAIN   = 16, &    ! exchange vector and use neighboring data
                      NSCARC_EXCHANGE_ZONE_NEIGHBORS = 17, &    ! exchange number of aggregation zones (AMG only)
                      NSCARC_EXCHANGE_ZONE_TYPES     = 18       ! exchange aggregation zones types (AMG only)

INTEGER, PARAMETER :: NSCARC_GRID_STRUCTURED         =  1, &    ! structured discretization
                      NSCARC_GRID_UNSTRUCTURED       =  2       ! unstructured discretization

INTEGER, PARAMETER :: NSCARC_HUGE_INT                = -999999999    ! undefined integer value
REAL(EB), PARAMETER:: NSCARC_HUGE_REAL               = -999999999_EB ! undefined integer value

INTEGER, PARAMETER :: NSCARC_INIT_UNDEF              =-999, &   ! initialize allocated array as undefined
                      NSCARC_INIT_NONE               =  -2, &   ! do not initialize allocated arrays
                      NSCARC_INIT_MINUS              =  -1, &   ! initialize allocated array with minus one
                      NSCARC_INIT_ZERO               =   0, &   ! initialize allocated array with zero
                      NSCARC_INIT_ONE                =   1, &   ! initialize allocated array with one
                      NSCARC_INIT_TRUE               =   2, &   ! initialize allocated array with .TRUE.
                      NSCARC_INIT_FALSE              =   3      ! initialize allocated array with .FALSE.

INTEGER, PARAMETER :: NSCARC_INTERPOL_BILINEAR       =  1, &    ! bilinear interpolation
                      NSCARC_INTERPOL_CONSTANT       =  2, &    ! constant interpolation
                      NSCARC_INTERPOL_CLASSICAL      =  3, &    ! classical interpolation
                      NSCARC_INTERPOL_DIRECT         =  4, &    ! direct interpolation
                      NSCARC_INTERPOL_STANDARD       =  5       ! standard interpolation

INTEGER, PARAMETER :: NSCARC_KRYLOV_CG               =  1, &    ! CG   as Krylov solver
                      NSCARC_KRYLOV_BICG             =  3, &    ! BICG as Krylov solver
                      NSCARC_KRYLOV_MAIN             =  1, &    ! Krylov solver as main solver
                      NSCARC_KRYLOV_COARSE           =  2       ! Krylov solver as coarse grid solver

INTEGER, PARAMETER :: NSCARC_LEVEL_MIN               =  0, &    ! minimum multigrid level
                      NSCARC_LEVEL_MAX               = 10, &    ! maximum multigrid level
                      NSCARC_LEVEL_SINGLE            =  1, &    ! only one grid level needed
                      NSCARC_LEVEL_MULTI             =  2, &    ! multiple grid levels needed for GMG meethod
                      NSCARC_LEVEL_AMG               =  3       ! number of grid levels specified by AMG method

INTEGER, PARAMETER :: NSCARC_MATRIX_BANDWISE         =  1, &    ! matrix in bandwise storage 
                      NSCARC_MATRIX_COMPACT          =  2, &    ! matrix in compact storage 
                      NSCARC_MATRIX_CONDENSED        =  3, &    ! matrix condensing for Neumann problems
                      NSCARC_MATRIX_CONNECTION       =  4, &    ! strength of connection matrix (AMG only)
                      NSCARC_MATRIX_FULL             =  5, &    ! full matrix allocation
                      NSCARC_MATRIX_GENERAL          =  6, &    ! general matrix (not necessarily symmetric)
                      NSCARC_MATRIX_LIGHT            =  7, &    ! reduced matrix allocation
                      NSCARC_MATRIX_MINIMAL          =  8, &    ! reduced matrix allocation
                      NSCARC_MATRIX_POISSON          =  9, &    ! Poisson matrix in compact storage (default)
                      NSCARC_MATRIX_POISSON_PROL     = 10, &    ! Poisson matrix times prolongation matrix (AMG only)
                      NSCARC_MATRIX_POISSON_SYM      = 11, &    ! symmetric Poisson matrix in compact storage (default)
                      NSCARC_MATRIX_PROLONGATION     = 12, &    ! Prolongation matrix (AMG only)
                      NSCARC_MATRIX_RESTRICTION      = 13, &    ! Restriction matrix (AMG only)
                      NSCARC_MATRIX_ZONES            = 14       ! aggregation zones matrix (AMG only)

INTEGER, PARAMETER :: NSCARC_MATVEC_GLOBAL           =  1, &    ! perform matrix vector product globally
                      NSCARC_MATVEC_LOCAL            =  2, &    ! perform matrix vector product locally
                      NSCARC_MATVEC_MEAN             =  3       ! perform matrix vector product with averaging overlaps

INTEGER, PARAMETER :: NSCARC_MAX_FACE_NEIGHBORS      = 10       ! max number neighbors per mesh face
INTEGER, PARAMETER :: NSCARC_MAX_STENCIL             =  7       ! max stencil size
INTEGER, PARAMETER :: NSCARC_MAX_BUFFER0             = 10       ! max stencil size

INTEGER, PARAMETER :: NSCARC_METHOD_KRYLOV           =  1, &    ! Krylov   -method as global solver
                      NSCARC_METHOD_MULTIGRID        =  2, &    ! multigrid-method as global solver
                      NSCARC_METHOD_LU               =  3, &    ! LU-decomposition based on MKL solver
                      NSCARC_METHOD_FFT              =  4, &    ! FFT method based on CRAYFISHPAK
                      NSCARC_METHOD_MGM              =  5       ! MGM-method based on LU

INTEGER, PARAMETER :: NSCARC_MKL_LOCAL               =  1, &    ! local LU-decompositions (PARDISO)
                      NSCARC_MKL_GLOBAL              =  2, &    ! global LU-decomposition(CLUSTER_SPARSE_SOLVER)
                      NSCARC_MKL_COARSE              =  3       ! LU-decomposition on coarse grid level

INTEGER, PARAMETER :: NSCARC_MULTIGRID_GEOMETRIC     =  1, &    ! geometric multigrid
                      NSCARC_MULTIGRID_ALGEBRAIC     =  2, &    ! algebraic multigrid
                      NSCARC_MULTIGRID_MAIN          =  1, &    ! multigrid solver as main solver
                      NSCARC_MULTIGRID_PRECON        =  2       ! multigrid solver as preconditioner

INTEGER, PARAMETER :: NSCARC_ORDER_ACTIVE            =  1, &    ! be next in line for aggregation
                      NSCARC_ORDER_LOCKED            = -1, &    ! locked for aggregation
                      NSCARC_ORDER_UNASSIGNED        =  0       ! still not defined for aggregation

INTEGER, PARAMETER :: NSCARC_PRECISION_SINGLE        =  1, &    ! single precision data 
                      NSCARC_PRECISION_DOUBLE        =  2       ! double precision data

INTEGER, PARAMETER :: NSCARC_RELAX_JAC               =  1, &    ! preconditioning by local JACOBI-methods
                      NSCARC_RELAX_SSOR              =  2, &    ! preconditioning by local SSOR-methods
                      NSCARC_RELAX_MJAC              =  3, &    ! preconditioning by local MJAC-methods (matrix form)
                      NSCARC_RELAX_MGS               =  4, &    ! preconditioning by local MGS-methods (matrix form)
                      NSCARC_RELAX_MSGS              =  5, &    ! preconditioning by local MSGS-methods (matrix form)
                      NSCARC_RELAX_MSOR              =  6, &    ! preconditioning by local MSOR-methods (matrix form)
                      NSCARC_RELAX_MSSOR             =  7, &    ! preconditioning by local MSSOR-methods (matrix form)
                      NSCARC_RELAX_ILU               =  8, &    ! preconditioning by local ILU-decompositions (own)
                      NSCARC_RELAX_FFT               =  9, &    ! preconditioning by local FFT-methods
                      NSCARC_RELAX_FFTO              = 10, &    ! preconditioning by local FFT-methods
                      NSCARC_RELAX_GMG               = 11, &    ! preconditioning by local GMG-methods
                      NSCARC_RELAX_AMG               = 12, &    ! preconditioning by local AMG-methods
                      NSCARC_RELAX_MKL               = 13       ! preconditioning by local LU-decompositions (MKL)

INTEGER, PARAMETER :: NSCARC_RHS_HOMOGENEOUS         =  1, &    ! homogeneous boundary conditions for global problem
                      NSCARC_RHS_INHOMOGENEOUS       =  2, &    ! inhomogeneous boundary conditions for global problem
                      NSCARC_RHS_DEFECT              =  3       ! RHS of precond. is set to defect of main iteration

INTEGER, PARAMETER :: NSCARC_SCOPE_GLOBAL            =  0, &    ! scope of defect correction is global
                      NSCARC_SCOPE_LOCAL             =  1       ! scope of defect correction is local

INTEGER, PARAMETER :: NSCARC_SIZE_MATRIX             =  1, &    ! size of system matrix for compact system
                      NSCARC_SIZE_TRANSFER           =  2       ! size of transfer matrices for compact system

INTEGER, PARAMETER :: NSCARC_SOLVER_MAIN             =  1, &    ! Main solver
                      NSCARC_SOLVER_PRECON           =  2, &    ! Preconditioner
                      NSCARC_SOLVER_SMOOTH           =  3, &    ! Smoother
                      NSCARC_SOLVER_COARSE           =  4       ! Coarse grid solver

INTEGER, PARAMETER :: NSCARC_STACK_ZERO              =   0, &   ! zero stack
                      NSCARC_STACK_ROOT              =   1, &   ! root stack
                      NSCARC_STACK_MAX               =  10, &   ! maximum stack
                      NSCARC_STACK_NOPARENT          = -99      ! no stack information

INTEGER, PARAMETER :: NSCARC_STAGE_ONE               =  1, &    ! primary scope for solution vectors
                      NSCARC_STAGE_TWO               =  2       ! secondary scope for solution vectors

INTEGER, PARAMETER :: NSCARC_STATE_PROCEED           =  0, &    ! proceed loop
                      NSCARC_STATE_CONV0             =  1, &    ! check convergence already for initial residual
                      NSCARC_STATE_CONV              =  2, &    ! check convergence for residual
                      NSCARC_STATE_DIVG              =  3       ! check divergence for residual

INTEGER, PARAMETER :: NSCARC_STENCIL_CONSTANT        =  1, &    ! constant matrix entries
                      NSCARC_STENCIL_VARIABLE        =  2, &    ! variable matrix entries, own version of scarc_daypyv
                      NSCARC_STENCIL_CONSTANT_LOCAL  =  3, &    ! constant-matrix entries with local corrections
                      NSCARC_STENCIL_VARIABLE_LOCAL  =  4       ! variable matrix entries with local matvec

REAL(EB), PARAMETER:: NSCARC_THRESHOLD_CONVERGENCE   = 1.0E-15_EB, & ! threshold for convergence
                      NSCARC_THRESHOLD_DIVGERGENCE   = 1.0E+15_EB    ! threshold for divergence

INTEGER, PARAMETER :: NSCARC_TWOLEVEL_NONE           =  0, &    ! no two levels, only one level
                      NSCARC_TWOLEVEL_ADD            =  1, &    ! additive 2-level method
                      NSCARC_TWOLEVEL_MUL            =  2, &    ! multiplicative 2-level method
                      NSCARC_TWOLEVEL_MUL2           =  3, &    ! multiplicative 2-level method
                      NSCARC_TWOLEVEL_COARSE         =  4, &    ! only coarse grid
                      NSCARC_TWOLEVEL_MACRO          =  5, &    ! 2. level is macro solver 
                      NSCARC_TWOLEVEL_SAMG           =  6       ! Coarse level based on SAMG-aggregation

INTEGER, PARAMETER :: NSCARC_VECTOR_ONE_X            =  1, &    ! selection parameter for 1-stage vector X
                      NSCARC_VECTOR_ONE_B            =  2, &    !       ...                              F
                      NSCARC_VECTOR_ONE_D            =  3, &    !       ...                              G
                      NSCARC_VECTOR_ONE_R            =  4, &    !       ...                              D
                      NSCARC_VECTOR_ONE_V            =  5, &    !       ...                              R
                      NSCARC_VECTOR_ONE_Y            =  6, &    !       ...                              Y
                      NSCARC_VECTOR_ONE_Z            =  7, &    !       ...                              Z
                      NSCARC_VECTOR_ONE_E            =  8, &    !       ...                              Z
                      NSCARC_VECTOR_TWO_X            =  9, &    ! selection parameter for 2-stage vector X
                      NSCARC_VECTOR_TWO_B            = 10, &    !       ...                              F
                      NSCARC_VECTOR_TWO_D            = 11, &    !       ...                              G
                      NSCARC_VECTOR_TWO_R            = 12, &    !       ...                              D
                      NSCARC_VECTOR_TWO_V            = 13, &    !       ...                              R
                      NSCARC_VECTOR_TWO_Y            = 14, &    !       ...                              Y
                      NSCARC_VECTOR_TWO_Z            = 15, &    !       ...                              Z
                      NSCARC_VECTOR_TWO_E            = 16, &    !       ...                              Z
                      NSCARC_VECTOR_H                = 17, &    ! selection parameter for vector H
                      NSCARC_VECTOR_HS               = 18       !       ...                      HS

INTEGER, PARAMETER :: NSCARC_ZONES_INITIAL           =  1, &    ! basic initial setup of zones information
                      NSCARC_ZONES_FILLUP            =  2, &    ! only filling up zones information
                      NSCARC_ZONES_COMPLETE          =  3, &    ! complete setting of zones information
                      NSCARC_ZONES_SEND              =  4, &    ! complete setting of zones information
                      NSCARC_ZONES_RECEIVE           =  5       ! complete setting of zones information
                   
INTEGER, PARAMETER :: NSCARC_UNDEF_INT               = -1, &    ! undefined integer value
                      NSCARC_ZERO_INT                =  0, &    ! zero integer value
                      NSCARC_ONE_INT                 =  1       ! one integer value

REAL(EB), PARAMETER:: NSCARC_UNDEF_REAL_EB           = -1.0_EB, &    ! undefined real value
                      NSCARC_ZERO_REAL_EB            =  0.0_EB, &    ! zero real value
                      NSCARC_ONE_REAL_EB             =  1.0_EB       ! one real value

REAL(EB), PARAMETER:: NSCARC_UNDEF_REAL_FB           = -1.0_FB, &    ! undefined real value
                      NSCARC_ZERO_REAL_FB            =  0.0_FB, &    ! zero real value
                      NSCARC_ONE_REAL_FB             =  1.0_FB       ! one real value

INTEGER, PARAMETER :: NSCARC_NONE                    = -123456789    ! dummy integer value for requests
CHARACTER(40), PARAMETER :: SCARC_NONE = 'NONE'                      ! dummy character value for requests

END MODULE SCARC_GLOBAL_CONSTANTS


! ================================================================================================================
! MODULE 'SCARC_VARIABLES': 
! Global SCARC variables 
! ================================================================================================================
MODULE SCARC_VARIABLES
USE PRECISION_PARAMETERS
USE SCARC_GLOBAL_CONSTANTS
IMPLICIT NONE

!
! ---------- Basic definitions
!
CHARACTER(40) :: SCARC_GRID              = 'STRUCTURED'     ! Type of discretization (STRUCTURED/UNSTRUCTURED)
CHARACTER(40) :: SCARC_METHOD            = 'NONE'           ! Requested solver method (KRYLOV/MULTIGRID)
CHARACTER(40) :: SCARC_MATRIX            = 'NONE'           ! Type of matrix storage (COMPACT/BANDWISE)
CHARACTER(40) :: SCARC_STENCIL           = 'VARIABLE'       ! Type of matrix stencil (CONSTANT/VARIABLE)
CHARACTER(40) :: SCARC_TWOLEVEL          = 'NONE'           ! Type of two-level method (NONE/ADDITIVE/MULTIPLICATIVE)

!
! ---------- General iteration parameters
!
CHARACTER(40) :: SCARC_ACCURACY          = 'ABSOLUTE'       ! Accuracy type (ABSOLUTE/RELATIVE)
REAL (EB)     :: SCARC_CAPPA             =  0.0_EB          ! Convergence rate of selected ScarC solver
INTEGER       :: SCARC_ITERATIONS        =  0               ! Number of iterations of selected ScaRC solver
REAL (EB)     :: SCARC_RESIDUAL          =  0.0_EB          ! Residual of global selected solver

!
! ---------- Parameters for coarse grid method
!
CHARACTER(40) :: SCARC_COARSE            = 'DIRECT'         ! Type of coarse grid solver (ITERATIVE/DIRECT)
REAL (EB)     :: SCARC_COARSE_ACCURACY   = 1.E-14_EB        ! Requested accuracy for iterative solver
INTEGER       :: SCARC_COARSE_ITERATIONS = 100              ! Max number of iterations for iterative solver
INTEGER       :: SCARC_COARSE_LEVEL      =  1               ! Coarse grid level for twolevel-Krylov method (default min level)
REAL (EB)     :: SCARC_COARSE_OMEGA      = 0.80E+0_EB       ! Relaxation parameter


!
! ---------- Parameters for Krylov-type methods
!
CHARACTER(40) :: SCARC_KRYLOV            = 'CG'             ! Type of Krylov-method (CG/BICG)
REAL (EB)     :: SCARC_KRYLOV_ACCURACY   = 1.E-8_EB         ! Requested accuracy for convergence
CHARACTER(40) :: SCARC_KRYLOV_INTERPOL   = 'CONSTANT'       ! twolevel-interpolation (CONSTANT/BILINEAR)
INTEGER       :: SCARC_KRYLOV_ITERATIONS = 1000             ! Max number of iterations

!
! ---------- Parameters for multigrid-type methods
!
CHARACTER(40) :: SCARC_MULTIGRID            = 'GEOMETRIC'   ! Type of MG-method (GEOMETRIC/ALGEBRAIC)
REAL (EB)     :: SCARC_MULTIGRID_ACCURACY   = 1.E-8_EB      ! Requested accuracy for convergence
CHARACTER(40) :: SCARC_MULTIGRID_COARSENING = 'GMG'         ! Coarsening strategy  (GMG/SAMG)
CHARACTER(3)  :: SCARC_MULTIGRID_CYCLE      = 'V'           ! Cycling type  (F/V/W)
CHARACTER(40) :: SCARC_MULTIGRID_INTERPOL   = 'CONSTANT'    ! Interpolation strategy (CONSTANT/BILINEAR)
INTEGER       :: SCARC_MULTIGRID_ITERATIONS = 100           ! Max number of iterations
INTEGER       :: SCARC_MULTIGRID_PRESMOOTH  = 4             ! Number of presmoothing iterations
INTEGER       :: SCARC_MULTIGRID_POSTSMOOTH = 4             ! Number of postsmoothing iterations
INTEGER       :: SCARC_MULTIGRID_LEVEL      = -1            ! User defined number of MG-levels (optionally, otherwise maximum)

!
! ---------- Parameters for smoothing method (used in multigrids-methods)
!
CHARACTER(40) :: SCARC_SMOOTH            = 'SSOR'           ! Smoother for MG (JACOBI/SSOR)
REAL (EB)     :: SCARC_SMOOTH_ACCURACY   = 1.E-8_EB         ! Requested accuracy for convergence
INTEGER       :: SCARC_SMOOTH_ITERATIONS = 4                ! Max number of iterations
REAL (EB)     :: SCARC_SMOOTH_OMEGA      = 0.80E+0_EB       ! Relaxation parameter
CHARACTER(40) :: SCARC_SMOOTH_SCOPE      = 'GLOBAL'         ! Scope of action (LOCAL/GLOBAL)

!
! ---------- Parameters for preconditioning method (used in Krylov-methods)
!
CHARACTER(40) :: SCARC_PRECON            = 'NONE'           ! Preconditioner for CG/BICG (JACOBI/SSOR/FFT/PARDISO/MG)
REAL (EB)     :: SCARC_PRECON_ACCURACY   = 1.E-10_EB        ! Requested accuracy for convergence
INTEGER       :: SCARC_PRECON_ITERATIONS = 100              ! Max number of iterations
REAL (EB)     :: SCARC_PRECON_OMEGA      = 1.50E+0_EB       ! Relaxation parameter
CHARACTER(40) :: SCARC_PRECON_SCOPE      = 'LOCAL'          ! Scope of action (LOCAL/GLOBAL)

!
! ---------- Parameter for MKL solver
!
CHARACTER(40) :: SCARC_MKL            = 'GLOBAL'            ! Type of MKL solver (LOCAL:Pardiso/GLOBAL:Cluster_Sparse_solver)
CHARACTER(40) :: SCARC_MKL_MTYPE      = 'SYMMETRIC'         ! Type of MKL matrix (SYMMETRIC/UNSYMMETRIC)
CHARACTER(6)  :: SCARC_MKL_PRECISION  = 'DOUBLE'            ! Single/double precision for MKL solver

!
! ---------- Dump out of error information and error handling
!
LOGICAL :: SCARC_ERROR_FILE = .FALSE.                       ! Print ScaRC statistics into chid_scarc.csv (TRUE/FALSE)
INTEGER :: IERROR = 0                                       ! general error flag - used at different positions

#ifdef WITH_SCARC_STANDALONE
LOGICAL :: SCARC_DUMP = .TRUE.                             ! Dump out several arrays for standalone use of ScaRC
#endif

! 
! ---------- Logical indicators for different methods and mechanisms
! 
LOGICAL :: IS_STRUCTURED      = .FALSE.                     ! is structured discretization?
LOGICAL :: IS_UNSTRUCTURED    = .FALSE.                     ! is unstructured discretization?
LOGICAL :: IS_PURE_NEUMANN    = .FALSE.                     ! is pure Neumann system?
LOGICAL :: IS_CG              = .FALSE.                     ! is Krylov-method?
LOGICAL :: IS_CG_ADD          = .FALSE.                     ! is additive twolevel-Krylov-method?
LOGICAL :: IS_CG_AMG          = .FALSE.                     ! is Krylov-method with AMG-preconditioning?
LOGICAL :: IS_CG_COARSE       = .FALSE.                     ! is only coarse grid preconditiner?
LOGICAL :: IS_CG_GMG          = .FALSE.                     ! is Krylov-method with GMG-preconditioning?
LOGICAL :: IS_CG_MUL          = .FALSE.                     ! is multiplicative Twolevel-Krylov-method ?
LOGICAL :: IS_CG_MACRO        = .FALSE.                     ! is macro coarse grid solver
LOGICAL :: IS_CG_SAMG         = .FALSE.                     ! is Krylov method with additive SAMG-coarse grid?
LOGICAL :: IS_MG              = .FALSE.                     ! is Multigrid-method?
LOGICAL :: IS_GMG             = .FALSE.                     ! is Geometric Multigrid-method?
LOGICAL :: IS_AMG             = .FALSE.                     ! is Algebraic Multigrid-method?
LOGICAL :: IS_FFT             = .FALSE.                     ! is FFT-method?
LOGICAL :: IS_FFTO            = .FALSE.                     ! is FFTO-method?
LOGICAL :: IS_MGM             = .FALSE.                     ! is McKeeney-Greengard-Mayo method?
LOGICAL :: IS_MKL             = .FALSE.                     ! is MKL-method?
LOGICAL :: IS_MKL_LEVEL(10)   = .FALSE.                     ! is level-dependent MKL method?
 
LOGICAL :: HAS_CSV_DUMP       = .FALSE.                     ! has CSV-file to be dumped out
LOGICAL :: HAS_TWO_LEVELS     = .FALSE.                     ! has two grid levels?
LOGICAL :: HAS_MULTI_LEVELS   = .FALSE.                     ! has multiple grid levels?
LOGICAL :: HAS_MULTI_GRIDS    = .FALSE.                     ! has multiple discretization types?

! 
! ---------- Globally used types for description of different solvers
! 
INTEGER :: TYPE_ACCURACY           = NSCARC_ACCURACY_ABSOLUTE    ! default type of requested accuracy
INTEGER :: TYPE_AGGREGATE          = NSCARC_ZONES_INITIAL        ! default aggregation type for exchange
INTEGER :: TYPE_COARSE             = NSCARC_COARSE_DIRECT        ! default type of coarse grid solver for multigrid method
INTEGER :: TYPE_COARSENING         = NSCARC_COARSENING_GMG       ! default GMG-coarsening in case of twolevel method
INTEGER :: TYPE_CYCLING            = NSCARC_CYCLING_V            ! default type of cycling for multigrid method
INTEGER :: TYPE_GRID               = NSCARC_GRID_STRUCTURED      ! default type of discretization (structured/unstructured)
INTEGER :: TYPE_EXCHANGE           = NSCARC_UNDEF_INT            ! no default type of data exchange
INTEGER :: TYPE_EXCHANGE_MATRIX    = NSCARC_MATRIX_POISSON       ! default matrix type for exchange
INTEGER :: TYPE_INTERPOL           = NSCARC_INTERPOL_CONSTANT    ! default type of interpolation method
INTEGER :: TYPE_KRYLOV             = NSCARC_KRYLOV_CG            ! default type of Krylov method (CG/BICG)
INTEGER :: TYPE_LEVEL(0:2)         = NSCARC_UNDEF_INT            ! default type of levels
INTEGER :: TYPE_MATRIX             = NSCARC_MATRIX_COMPACT       ! default type of storage for matrix
INTEGER :: TYPE_MATVEC             = NSCARC_MATVEC_GLOBAL        ! default type of matrix-vector multiplication
INTEGER :: TYPE_METHOD             = NSCARC_METHOD_KRYLOV        ! default type of ScaRC method
INTEGER :: TYPE_MKL(0:10)          = NSCARC_UNDEF_INT            ! no default type of MKL for single levels
INTEGER :: TYPE_MKL_PRECISION      = NSCARC_PRECISION_DOUBLE     ! default double precision MKL solver
INTEGER :: TYPE_MULTIGRID          = NSCARC_MULTIGRID_GEOMETRIC  ! default type of multigrid method (GMG/AMG)
INTEGER :: TYPE_PARENT             = NSCARC_UNDEF_INT            ! no default type of parent (calling) solver
INTEGER :: TYPE_PRECON             = NSCARC_UNDEF_INT            ! no default type of preconditioner for iterative solver
INTEGER :: TYPE_RELAX              = NSCARC_UNDEF_INT            ! no default type of preconditioner for iterative solver
INTEGER :: TYPE_SCOPE(0:2)         = NSCARC_SCOPE_LOCAL          ! default scope types are local
INTEGER :: TYPE_SMOOTH             = NSCARC_UNDEF_INT            ! no default type of smoother for multigrid method
INTEGER :: TYPE_SOLVER             = NSCARC_SOLVER_MAIN          ! default type of surrounding solver stage
INTEGER :: TYPE_STAGE              = NSCARC_STAGE_ONE            ! default type of surrounding solver stage
INTEGER :: TYPE_STENCIL            = NSCARC_STENCIL_VARIABLE     ! default type of storage for matrix
INTEGER :: TYPE_TWOLEVEL           = NSCARC_TWOLEVEL_NONE        ! default type of two-level method
INTEGER :: TYPE_VECTOR             = NSCARC_UNDEF_INT            ! no default type of vector to point to

! 
! ---------- Globally used parameters
! 
INTEGER :: NLEVEL_MIN, NLEVEL_MAX                           ! minimum and maximum number of multigrid levels
INTEGER :: NC_GLOBAL(20) = 0                                ! number of global cells
INTEGER :: N_DIRIC_GLOBAL(20) = 0                           ! global number of Dirichlet BCs
INTEGER :: N_STACK_TOTAL                                    ! maximum number of used solvers in stack

INTEGER :: N_REQ, N_EXCHANGES, TAG                           ! Information for data exchange
INTEGER :: SNODE, RNODE                                     ! Process identifier for data exchange

INTEGER,  ALLOCATABLE, DIMENSION (:)  :: REQ                ! Request array for data exchange
INTEGER,  ALLOCATABLE, DIMENSION (:)  :: COUNTS             ! Counter array for data exchange
INTEGER,  ALLOCATABLE, DIMENSION (:)  :: DISPLS             ! Displacement array for data exchange
INTEGER,  ALLOCATABLE, DIMENSION (:)  :: MESH_INT           ! Local integer data array for data exchange
REAL(EB), ALLOCATABLE, DIMENSION (:)  :: MESH_REAL          ! Local real data array for data exchange

INTEGER  :: GLOBAL_INT,  RANK_INT
REAL(EB) :: GLOBAL_REAL, RANK_REAL

INTEGER :: FACE_ORIENTATION(6) = (/1,-1,2,-2,3,-3/)         ! Coordinate direction related order of mesh faces

CHARACTER(10) :: CFORM_INT, CFORM_REAL
CHARACTER(60) :: CNAME


! 
! ---------- Grid transfer information
! 
REAL(EB), PARAMETER :: SCALR  = 0.015625_EB                 ! some scaling parameters
REAL(EB), PARAMETER :: SCALP  = 0.0625_EB                    
REAL(EB), PARAMETER :: W1     =  1.0_EB                     ! some weighting parameters
REAL(EB), PARAMETER :: W3     =  3.0_EB                      
REAL(EB), PARAMETER :: W4     =  4.0_EB                     
REAL(EB), PARAMETER :: W9     =  9.0_EB                    
REAL(EB), PARAMETER :: W12    = 12.0_EB                   
REAL(EB), PARAMETER :: W16    = 16.0_EB                  


END MODULE SCARC_VARIABLES


! ================================================================================================================
! MODULE 'SCARC_ITERATION_ENVIRONMENT': 
! Iteration parameters and handles
! ================================================================================================================
MODULE SCARC_ITERATION_ENVIRONMENT
USE PRECISION_PARAMETERS
  
REAL(EB) :: DT, DTI
REAL(EB) :: OMEGA, EPS, RES, RESIN = -1.0_EB, CAPPA = -1.0_EB, ERR
INTEGER :: X, B, D, R, V, Y, Z
#ifdef WITH_SCARC_DEBUG
INTEGER :: E
#endif
INTEGER :: NIT, ITE
INTEGER :: ITE_PRES=0, ITE_TOTAL=0, ITE_CG=0, ITE_MG=0, ITE_LU=0, ITE_SMOOTH=0, ITE_COARSE=0, ITE_GLOBAL = 0

END MODULE SCARC_ITERATION_ENVIRONMENT


! ================================================================================================================
! MODULE 'SCARC_TYPES': 
! Collection of different types for ScaRC
! ================================================================================================================
MODULE SCARC_TYPES
USE PRECISION_PARAMETERS
USE SCARC_GLOBAL_CONSTANTS
#ifdef WITH_MKL
USE MKL_PARDISO
USE MKL_CLUSTER_SPARSE_SOLVER
#endif
IMPLICIT NONE

!
! ---------- Messaging and debugging mechanisms
!
TYPE SCARC_MESSAGE_TYPE

   CHARACTER(100) :: TEXT
   CHARACTER(60)  :: FILE_DEBUG, FILE_STAT, FILE_DUMP, FILE_VERBOSE, FILE_CPU
   CHARACTER(60)  :: FILE_CNN, FILE_CNN1, FILE_CNN2, FILE_CNN3
   CHARACTER(120) :: FILE_SCARC

   INTEGER :: LU_DEBUG = 0, LU_STAT = 0, LU_DUMP = 0, LU_VERBOSE = 0, LU_SCARC = 0, LU_CPU = 0
   INTEGER :: LU_CNN = 0, LU_CNN1 = 0, LU_CNN2 = 0, LU_CNN3 = 0

END TYPE SCARC_MESSAGE_TYPE

!
! ---------- Information about the mesh decomposition
!
TYPE SCARC_SUBDIVISION_TYPE

   INTEGER, ALLOCATABLE, DIMENSION (:,:) :: ORDER
   INTEGER, ALLOCATABLE, DIMENSION (:,:) :: NEIGHBORS
   INTEGER, ALLOCATABLE, DIMENSION (:)   :: N_NEIGHBORS
   INTEGER :: N_NEIGHBORS_TOTAL, N_CYCLES

END TYPE SCARC_SUBDIVISION_TYPE

!
! ---------- Face information related to wall cells and neighbors
!
TYPE SCARC_FACE_TYPE

   REAL(EB), POINTER, DIMENSION(:) :: DH              ! step size vector between adjacent faces
   REAL(EB) :: SCAL_FACE                              ! increments for matrix subdiagonal for cells right face
   REAL(EB) :: SCAL_INSIDE                            ! increments for matrix subdiagonal for cells between opposite faces
   REAL(EB) :: SCAL_BOUNDARY                          ! increments for boundary conditions
   REAL(EB) :: SCAL_DIRICHLET                         ! scaling factors for Dirichlet BC's
   REAL(EB) :: SCAL_NEUMANN                           ! scaling factors for Neumann BC's

   INTEGER, ALLOCATABLE, DIMENSION(:) :: NEIGHBORS    ! adjacent neighbors at that face
   INTEGER  :: N_NEIGHBORS = 0                        ! number of adjacent neighbors 

   INTEGER  :: NOP = 0                                ! number of cells between opposite faces
   INTEGER  :: NX, NY, NZ                             ! cells in different directions on that face
   INTEGER  :: NCW0 = 0, NCW = 0                      ! number of first wall cell and total number of wall cells
   INTEGER  :: NOFFX = 0, NOFFY = 0, NOFFZ = 0        ! offsets to next internal cell in that face direction

END TYPE SCARC_FACE_TYPE


!
! ---------- Wall information related to neighbors and BC's
!
TYPE SCARC_WALL_TYPE

   ! different properties of wall cell
   INTEGER :: BTYPE = 0                               ! BTYPE of wall cell (Dirichlet/Neumann/Internal)
   INTEGER :: BOUNDARY_TYPE = 0                       ! boundary type of wall cell (Solid/Interpolated/Open))
   INTEGER :: IOR = 0                                 ! orientation of wall cell
   INTEGER :: NOM = 0                                 ! adjacent neighbor at wall cell

   INTEGER :: ICW = NSCARC_UNDEF_INT                  ! internal wall cell for IW
   INTEGER :: IXG, IYG, IZG                           ! x-, y- and z-indices of ghost cells
   INTEGER :: IXW, IYW, IZW                           ! x-, y- and z-indices of (internal) wall cells
   INTEGER :: IXN(2), IYN(2), IZN(2)                  ! x-, y- and z-indices of neighboring cells
   INTEGER :: ICE, ICG

END TYPE SCARC_WALL_TYPE


!
! ---------- Obstruction information
!
TYPE SCARC_OBST_TYPE
   INTEGER :: I1, I2, J1, J2, K1, K2                         ! cell indices of obstructions
END TYPE SCARC_OBST_TYPE


!
! ---------- Compact matrix entries which will be exchanged during generation of condensed system
!
TYPE SCARC_MATRIX_COMPACT_CONDENSED_TYPE

   REAL(EB) :: VAL1(NSCARC_MAX_STENCIL) = 0.0_EB         ! original values (double precision)
   REAL(EB) :: VAL2(NSCARC_MAX_STENCIL) = 0.0_EB         ! condensed values (double precision)

   INTEGER :: COL(NSCARC_MAX_STENCIL) = 0                ! column pointers
   INTEGER :: PTR(NSCARC_MAX_STENCIL) = 0                ! storage pointer
   INTEGER :: N_COL

END TYPE SCARC_MATRIX_COMPACT_CONDENSED_TYPE


!
! ---------- BANDWISE matrix entries which will exchanged during generation of condensed system
!
TYPE SCARC_MATRIX_BANDWISE_CONDENSED_TYPE

   REAL(EB) :: VAL1(NSCARC_MAX_STENCIL) = 0.0_EB         ! original values (double precision)
   REAL(EB) :: VAL2(NSCARC_MAX_STENCIL) = 0.0_EB         ! condensed values (double precision)
   INTEGER :: IOR0 = 0                                   ! position pointer
   INTEGER :: ICO = 0                                    ! cell pointer

END TYPE SCARC_MATRIX_BANDWISE_CONDENSED_TYPE


!
! ---------- Compact sparse row (COMPACT) storage technique for matrices
! Is based on three arrays:
!    - non-zero matrix values
!    - corresponding columns pointers
!    - row pointers
!
TYPE SCARC_MATRIX_COMPACT_TYPE

   REAL(EB), ALLOCATABLE, DIMENSION (:) :: VAL                ! values of matrix (real precision)
   REAL(EB), ALLOCATABLE, DIMENSION (:) :: ILU                ! ILU-decomposition
   REAL(EB), ALLOCATABLE, DIMENSION (:) :: RELAX              ! workspace for relaxation
   REAL(EB), DIMENSION (-3:3)           :: STENCIL            ! store basic stencil information in single precision

   REAL(FB), ALLOCATABLE, DIMENSION (:) :: VAL_FB             ! values of matrix (single precision)
   REAL(FB), ALLOCATABLE, DIMENSION (:) :: RELAX_FB           ! workspace for relaxation
   REAL(FB), DIMENSION (-3:3)           :: STENCIL_FB         ! store basic stencil information in single precision

   TYPE (SCARC_MATRIX_COMPACT_CONDENSED_TYPE) :: CONDENSED(NSCARC_MAX_STENCIL)

   INTEGER,  ALLOCATABLE, DIMENSION (:) :: ROW                ! row pointer
   INTEGER,  ALLOCATABLE, DIMENSION (:) :: COL, COLG    ! column pointers and global column pointers

   INTEGER :: POS(-3:3) = 0                                   ! Position of IOR's in STENCIL
   INTEGER :: N_CONDENSED = 0                                 ! number of condensed entries in matrix
   INTEGER :: N_VAL = 0                                       ! number of matrix values
   INTEGER :: N_ROW = 0                                       ! number of matrix rows
   INTEGER :: N_STENCIL = 0                                   ! number of points in matrix stencil
   INTEGER :: N_STENCIL_MAX = 0                               ! max stencil size (AMG only)
   INTEGER :: NTYPE = 0                                       ! matrix type

   CHARACTER(40) :: CNAME                                     ! Name of matrix

END TYPE SCARC_MATRIX_COMPACT_TYPE

! 
! ---------- bandwise storage technique for matrices
! The entries are stored one diagonal after the other
! Missing entries of subdiagonals are filled with zero
! Is based on two arrays:
!      - non-zero matrix entries diagonal-wise
!      - the offsets from the main diagonal
!
TYPE SCARC_MATRIX_BANDWISE_TYPE

   REAL(EB), ALLOCATABLE, DIMENSION (:)   :: AUX           ! auxiliary vector (double precision)
   REAL(EB), ALLOCATABLE, DIMENSION (:,:) :: VAL           ! values of matrix (double precision)
   REAL(EB), ALLOCATABLE, DIMENSION (:,:) :: RELAX         ! workspace for relaxation
   REAL(EB), ALLOCATABLE, DIMENSION (:)   :: RELAXD        ! workspace for relaxation - only for diagonal scaling
   REAL(EB), DIMENSION (-3:3)             :: STENCIL       ! store basic stencil information (double precision)

   REAL(FB), ALLOCATABLE, DIMENSION (:)   :: AUX_FB        ! auxiliary vector (single precision)
   REAL(FB), ALLOCATABLE, DIMENSION (:,:) :: VAL_FB        ! values of matrix (single precision)
   REAL(FB), ALLOCATABLE, DIMENSION (:,:) :: RELAX_FB      ! workspace for relaxation
   REAL(FB), DIMENSION (-3:3)             :: STENCIL_FB    ! store basic stencil information (single precision)

   TYPE (SCARC_MATRIX_BANDWISE_CONDENSED_TYPE) :: CONDENSED(NSCARC_MAX_STENCIL)

   CHARACTER(40) :: CNAME                                  ! Name of matrix

   INTEGER,  DIMENSION (-3:3) :: OFFSET                    ! offset pointers
   INTEGER,  DIMENSION (-3:3) :: LENGTH                    ! Relevant diagonal length 
   INTEGER,  DIMENSION (-3:3) :: SOURCE                    ! Source address in corresponding diagonal
   INTEGER,  DIMENSION (-3:3) :: TARGET                    ! Target address in corresponding diagonal

   INTEGER :: POS(-3:3) = 0                                ! position of IOR's in STENCIL and in matrix storage array
   INTEGER :: N_STENCIL = 0                                ! number of points in matrix stencil
   INTEGER :: N_CONDENSED = 0                              ! number of condensed entries in matrix
   INTEGER :: N_VAL = 0                                    ! number of matrix values in general and symmetric cass
   INTEGER :: N_DIAG = 0                                   ! length of main diagonal

END TYPE SCARC_MATRIX_BANDWISE_TYPE

! 
! ---------- Pressure information
! 
TYPE SCARC_PRESSURE_TYPE

   REAL(EB), ALLOCATABLE, DIMENSION (:)       :: B_OLD, B_NEW
   REAL(EB), ALLOCATABLE, DIMENSION (:, :, :) :: H_NEW, HS_NEW
   REAL(EB), ALLOCATABLE, DIMENSION (:, :, :) :: H_OLD, HS_OLD
   REAL(EB), ALLOCATABLE, DIMENSION (:, :, :) :: H_DIFF
   REAL(EB) :: DIFF_H = 0.0_EB, DIFF_HS = 0.0_EB

END TYPE SCARC_PRESSURE_TYPE

! 
! ---------- Workspace for FFT preconditioners
! 
TYPE SCARC_FFT_TYPE

   REAL(EB), ALLOCATABLE, DIMENSION (:)       :: SAVE1, WORK, HX
   REAL(EB), ALLOCATABLE, DIMENSION (:, :)    :: BXS, BXF, BYS, BYF, BZS, BZF
   REAL(EB), ALLOCATABLE, DIMENSION (:, :, :) :: PRHS
   REAL(EB) :: XS, XF, YS, YF, ZS, ZF
   REAL(EB) :: POIS_PTB = 0.0_EB, XLM = 0.0_EB

   INTEGER :: LSAVE, LWORK
   INTEGER :: LBC, MBC, NBC
   INTEGER :: ITRN, JTRN, KTRN
   INTEGER :: IBAR, JBAR, KBAR

END TYPE SCARC_FFT_TYPE

#ifdef WITH_MKL
! 
! ---------- MKL information
! 
TYPE SCARC_MKL_TYPE

   CHARACTER(40) :: CNAME                                  ! Name of matrix

   INTEGER, ALLOCATABLE :: IPARM(:)
   INTEGER :: MAXFCT, MNUM, MTYPE, PHASE, NRHS, ERROR, MSGLVL
   INTEGER :: PERM(1)

   TYPE (MKL_PARDISO_HANDLE),               ALLOCATABLE :: PT_H(:), PT(:)
   TYPE (MKL_CLUSTER_SPARSE_SOLVER_HANDLE), ALLOCATABLE :: CT_H(:), CT(:)

END TYPE SCARC_MKL_TYPE
#endif

! 
! ---------- Different scopes for solution, rhs and auxiliary vectors of different solvers
! 
TYPE SCARC_STAGE_TYPE

   REAL (EB), ALLOCATABLE, DIMENSION (:) :: X          ! solution vector double precision
   REAL (EB), ALLOCATABLE, DIMENSION (:) :: B          ! right hand side vector double precision
   REAL (EB), ALLOCATABLE, DIMENSION (:) :: R          ! residual vector double precision
   REAL (EB), ALLOCATABLE, DIMENSION (:) :: D          ! auxiliary vector double precision
   REAL (EB), ALLOCATABLE, DIMENSION (:) :: V          ! auxiliary vector double precision
   REAL (EB), ALLOCATABLE, DIMENSION (:) :: Y          ! auxiliary vector double precision
   REAL (EB), ALLOCATABLE, DIMENSION (:) :: Z          ! auxiliary vector double precision

   REAL (FB), ALLOCATABLE, DIMENSION (:) :: X_FB       ! solution vector vector single precision
   REAL (FB), ALLOCATABLE, DIMENSION (:) :: B_FB       ! right hand side vector single precision
   REAL (FB), ALLOCATABLE, DIMENSION (:) :: R_FB       ! residual vector single precision
   REAL (FB), ALLOCATABLE, DIMENSION (:) :: V_FB       ! auxiliary vector single precision

#ifdef WITH_SCARC_DEBUG
   REAL (EB), ALLOCATABLE, DIMENSION (:) :: E          ! error vector double precision
#endif

END TYPE SCARC_STAGE_TYPE

! 
! ---------- Multigrid type - to be extended for algebraic multigrid
! 
TYPE SCARC_MULTIGRID_TYPE

   REAL(EB) :: APPROX_SPECTRAL_RADIUS = 2.0_EB             ! Relaxation parameter (AMG only)
   REAL(EB) :: AMG_TOL = 0.25_EB                           ! Tolerance for coarsening
   REAL(EB) :: OMEGA = 1.0_EB                              ! Relaxation parameter

   INTEGER :: CYCLING(2) = 0                               ! Counter for multigrid cycling
   INTEGER :: N_PRESMOOTH, N_POSTSMOOTH                    ! Number of pre- and post-processing steps

END TYPE SCARC_MULTIGRID_TYPE

! 
! ---------- Information related to discretization type (structured/unstructured)
! 
TYPE SCARC_GRID_TYPE

   ! Basic boundary and interface information
   TYPE (SCARC_WALL_TYPE), ALLOCATABLE, DIMENSION(:) :: WALL   ! wall information

   ! Matrices in different storage types
   TYPE (SCARC_MATRIX_BANDWISE_TYPE):: POISSONB                ! Poisson matrix in bandwise storage technique
   TYPE (SCARC_MATRIX_COMPACT_TYPE) :: POISSONC                ! Poisson matrix in compact storage technique (default)
   TYPE (SCARC_MATRIX_COMPACT_TYPE) :: GALERKIN                ! Galerkin matrix (AMG only)
#ifdef WITH_MKL
   TYPE (SCARC_MATRIX_COMPACT_TYPE) :: POISSONC_SYM            ! symmetric part of compact Poisson matrix (only for MKL)
   TYPE (SCARC_MATRIX_COMPACT_TYPE) :: GALERKIN_SYM            ! Galerkin matrix symmetric version (AMG only)
#endif

   TYPE (SCARC_MATRIX_COMPACT_TYPE) :: PROLONGATION            ! prolongation matrix
   TYPE (SCARC_MATRIX_COMPACT_TYPE) :: RESTRICTION             ! restriction matrix
   TYPE (SCARC_MATRIX_COMPACT_TYPE) :: CONNECTION              ! strength of connection matrix
   TYPE (SCARC_MATRIX_COMPACT_TYPE) :: ZONES                   ! aggregation zones matrix
   TYPE (SCARC_MATRIX_COMPACT_TYPE) :: AP                      ! A*P matrix - only temporarily (AMG only)

   REAL(EB), ALLOCATABLE, DIMENSION (:) :: MEASURES            ! measure for grid coarsening (AMG only)
   REAL(EB), ALLOCATABLE, DIMENSION (:) :: NULLSPACE           ! nullspace vector
   REAL(EB), ALLOCATABLE, DIMENSION (:) :: AUX1, AUX2          ! auxiliary vectors
   REAL(EB), ALLOCATABLE, DIMENSION (:) :: DIAG                ! matrix diagonal, possible inverted (AMG only)
   REAL(EB), ALLOCATABLE, DIMENSION (:) :: QQ, RR              ! workspace for QR-decompostion (AMG only)

   INTEGER,  ALLOCATABLE, DIMENSION (:) :: ZONES_GLOBAL        ! global zone numbers (grid coarsening of AMG only)
   INTEGER,  ALLOCATABLE, DIMENSION (:) :: ZONES_LOCAL         ! local  zone numbers (grid coarsening of AMG only)
   INTEGER,  ALLOCATABLE, DIMENSION (:) :: CELLTYPES           ! celltype for grid coarsening (AMG only)
   INTEGER,  ALLOCATABLE, DIMENSION (:) :: CPOINTS             ! ZONES for grid coarsening (AMG only)
   INTEGER,  ALLOCATABLE, DIMENSION (:) :: ORDER               ! search order for aggregation
   INTEGER,  ALLOCATABLE, DIMENSION (:) :: LOCAL_TO_GLOBAL     ! mapping from local to global numbering 

   ! Pointer arrays for data exchange with neighbors
   INTEGER, ALLOCATABLE, DIMENSION (:)   :: ICE_TO_IWG         ! mapping from ICE to IWG
   INTEGER, ALLOCATABLE, DIMENSION (:)   :: ICE_TO_ICN         ! mapping from ICE to ICN
   INTEGER, ALLOCATABLE, DIMENSION (:)   :: ICE_TO_ECELL       ! mapping from ICE to global neighboring cell
   INTEGER, ALLOCATABLE, DIMENSION (:)   :: ICE_TO_ICELL       ! mapping from ICE to local neighboring cell
   INTEGER, ALLOCATABLE, DIMENSION (:)   :: ICE_TO_GCELL       ! mapping from ICE to global neighboring cell
   INTEGER, ALLOCATABLE, DIMENSION (:)   :: ICE_TO_OCELL       ! mapping from ICE to local neighboring cell
   INTEGER, ALLOCATABLE, DIMENSION (:)   :: ICE_TO_EZONE       ! mapping from ICE to internal own zone
   INTEGER, ALLOCATABLE, DIMENSION (:)   :: ICE_TO_IZONE       ! mapping from ICE to external own zone
   INTEGER, ALLOCATABLE, DIMENSION (:)   :: ICE_TO_GZONE       ! mapping from ICE to global neighboring zone
   INTEGER, ALLOCATABLE, DIMENSION (:)   :: ICE_TO_OZONE       ! mapping from ICE to local neighboring zone

   INTEGER, ALLOCATABLE, DIMENSION (:,:) :: ICG_TO_ICW         ! mapping from ICG to ICW
   INTEGER, ALLOCATABLE, DIMENSION (:,:) :: ICG_TO_ICE         ! mapping from ICG to ICE 
   INTEGER, ALLOCATABLE, DIMENSION (:)   :: ICG_TO_IWG         ! mapping from ICG to IWG
   INTEGER, ALLOCATABLE, DIMENSION (:)   :: ICG_TO_ECELL       ! mapping from ICG to global neighboring cell
   INTEGER, ALLOCATABLE, DIMENSION (:)   :: ICG_TO_ICELL       ! mapping from ICG to local neighboring cell
   INTEGER, ALLOCATABLE, DIMENSION (:)   :: ICG_TO_GCELL       ! mapping from ICG to global neighboring cell
   INTEGER, ALLOCATABLE, DIMENSION (:)   :: ICG_TO_OCELL       ! mapping from ICG to local neighboring cell
   INTEGER, ALLOCATABLE, DIMENSION (:)   :: ICG_TO_IZONE       ! mapping from ICG to internal own zone
   INTEGER, ALLOCATABLE, DIMENSION (:)   :: ICG_TO_EZONE       ! mapping from ICG to external own zone
   INTEGER, ALLOCATABLE, DIMENSION (:)   :: ICG_TO_GZONE       ! mapping from ICG to global neighboring zone
   INTEGER, ALLOCATABLE, DIMENSION (:)   :: ICG_TO_OZONE       ! mapping from ICG to local neighboring zone

   ! Assignment of cell coordinates
   INTEGER, ALLOCATABLE, DIMENSION (:)   :: ICX                ! I-coordinate of cell IC
   INTEGER, ALLOCATABLE, DIMENSION (:)   :: ICY                ! J-coordinate of cell IC
   INTEGER, ALLOCATABLE, DIMENSION (:)   :: ICZ                ! J-coordinate of cell IC

   ! Number and state of single cells
   INTEGER, ALLOCATABLE, DIMENSION (:,:,:) :: CELL_NUMBER      ! numbering of single cells

   ! Cell numbers of all meshes and offsets between meshes
   INTEGER, ALLOCATABLE, DIMENSION (:) :: NC_LOCAL             ! number of cells in local meshes
   INTEGER, ALLOCATABLE, DIMENSION (:) :: NC_OFFSET            ! offset in cell numbering between meshes
   INTEGER :: NC_GLOBAL = NSCARC_ZERO_INT                      ! global number of cells in all meshes

   ! Local numbers of internal, extended and ghost cells
   INTEGER :: NC   = NSCARC_ZERO_INT                           ! number of cells needed for matrix
   INTEGER :: NW   = NSCARC_ZERO_INT                           ! number of wall cells
   INTEGER :: NCG  = NSCARC_ZERO_INT                           ! number of ghost cells 
   INTEGER :: NCGI = NSCARC_ZERO_INT                           ! number of ghost zones in internal direction
   INTEGER :: NCGE = NSCARC_ZERO_INT                           ! number of ghost zones in external direction
   INTEGER :: NZG  = NSCARC_ZERO_INT                           ! number of ghost zones 
   INTEGER :: NCE  = NSCARC_ZERO_INT                           ! number of extended cells
   INTEGER :: NCE2 = NSCARC_ZERO_INT                           ! number of extended cells, second layer
   INTEGER :: ICE2 = NSCARC_ZERO_INT                           ! counter for extended cells, second layer

   ! Number of Dirichlet and Neumann boundary cells
   INTEGER :: N_DIRIC   = NSCARC_ZERO_INT                      ! number of Dirichlet BCs
   INTEGER :: N_NEUMANN = NSCARC_ZERO_INT                      ! number of Neumann BCs

   INTEGER :: N_FINE   = NSCARC_ZERO_INT                       ! Number of fine cells (AMG only)
   INTEGER :: N_COARSE = NSCARC_ZERO_INT                       ! Number of coarse cells (AMG only)
   INTEGER :: N_ZONES  = NSCARC_ZERO_INT                       ! Number of zones (AMG only)
   INTEGER :: N_STENCIL_MAX  = 20                              ! Max stencil size (AMG only)

   ! Pointer variables and arrays for data exchange with neighbors
   INTEGER :: ICG = 0, ICG2 = 0, ICE = 0                       ! ghost cell and extended cell pointers

   INTEGER :: NLEN_BUFFER_LAYER1                               ! Length for single layer length exchange on that level
   INTEGER :: NLEN_BUFFER_LAYER2                               ! Length for double layer length exchange on that level
   INTEGER :: NLEN_BUFFER_LAYER4                               ! Length for fourfold length exchange on that level
   INTEGER :: NLEN_BUFFER_STENCIL                              ! Length for stencil layer length exchange on that level
   INTEGER :: NLEN_BUFFER_FULL                                 ! Length for full length exchange on that level

END TYPE SCARC_GRID_TYPE


#ifdef WITH_SCARC_MGM
! 
! ---------- McKenney-Greengard-Mayo method
! 
TYPE SCARC_MGM_TYPE

   REAL(EB), ALLOCATABLE, DIMENSION (:)   :: H1, H2       ! H pressure vectors of different parts
   REAL(EB), ALLOCATABLE, DIMENSION (:)   :: US, VS, WS   ! Velocity components along internal BC's
   REAL(EB), ALLOCATABLE, DIMENSION (:,:) :: A, IA        ! Lower part of LU-decomposition for unstructured AC
   REAL(EB), ALLOCATABLE, DIMENSION (:,:) :: L, IL        ! Lower part of LU-decomposition for unstructured AC
   REAL(EB), ALLOCATABLE, DIMENSION (:,:) :: U, IU        ! Upper part of LU-decomposition for unstructured AC
   REAL(EB), ALLOCATABLE, DIMENSION (:)   :: B, X, Y      ! Right hand side, solution vectors

   INTEGER,  ALLOCATABLE, DIMENSION (:)   :: PERM         ! Permutation vector for reordering of matrix rows
   INTEGER :: NW1, NW2, NWI, NWE                          ! Range of IW's with non-zero B-values

END TYPE SCARC_MGM_TYPE
#endif


! 
! ---------- Collection of grid level related information on single mesh
! 
TYPE SCARC_LEVEL_TYPE

   ! Administrative structures for different components based on given grid level
   TYPE (SCARC_FACE_TYPE), ALLOCATABLE, DIMENSION(:) :: FACE   ! Face information
   TYPE (SCARC_OBST_TYPE), ALLOCATABLE, DIMENSION(:) :: OBST   ! Obstruction information
   TYPE (SCARC_STAGE_TYPE), DIMENSION(2) :: STAGE              ! Hierarchy of solvers and related working vectors
   TYPE (SCARC_GRID_TYPE)      :: STRUCTURED, UNSTRUCTURED     ! Structured and unstructured grid information
   TYPE (SCARC_MULTIGRID_TYPE) :: MG                           ! Multigrid method information
   TYPE (SCARC_FFT_TYPE)       :: FFT                          ! FFT preconditioner based on CRAYFISHPAK
   TYPE (SCARC_PRESSURE_TYPE)  :: PRES                         ! Multigrid method information
#ifdef WITH_SCARC_MGM
   TYPE (SCARC_MGM_TYPE)       :: MGM                          ! McKenney-Greengar-Mayo method 
#endif
#ifdef WITH_MKL
   TYPE (SCARC_MKL_TYPE)       :: MKL                          ! MKL preconditioner based on Intel MKL
#endif

   ! Coordinate information
   REAL(EB), ALLOCATABLE, DIMENSION (:) :: XCOR, YCOR, ZCOR    ! Coordinate vectors in x-, y- and z-direction
   REAL(EB), ALLOCATABLE, DIMENSION (:) :: XMID, YMID, ZMID    ! Midpoint vectors in x-, y- and z-direction
   REAL(EB), ALLOCATABLE, DIMENSION (:) :: DXL, DYL, DZL       ! Step size vectors in x-, y- and z-direction
   REAL(EB) :: DX , DY , DZ                                    ! Step sizes in x-, y- and z-direction
   REAL(EB) :: DXI, DYI, DZI                                   ! Inversed of step sizes in x-, y- and z-direction
   REAL(EB) :: DXI2, DYI2, DZI2                                ! Squared and inversed step sizes in x-, y- and z-direction

   ! Cell and wall index information:
   !    - on the finest level, the original arrays from FDS are used
   !    - separate arrays will only be allocated for coarser levels
   !    - to address them on all levels, corresponding pointers are used
   INTEGER, POINTER, DIMENSION (:,:,:)     :: CELL_INDEX_PTR   ! Pointer to CELL_INDEX
   INTEGER, POINTER, DIMENSION (:,:)       :: WALL_INDEX_PTR   ! Pointer to WALL_INDEX
   INTEGER, ALLOCATABLE, DIMENSION (:,:,:) :: CELL_INDEX       ! Cell index list (only allocated for coarser levels)
   INTEGER, ALLOCATABLE, DIMENSION (:,:)   :: WALL_INDEX       ! Wall index list (only allocated for coarser levels)

   LOGICAL, ALLOCATABLE, DIMENSION (:,:,:) :: IS_SOLID         ! State of single cells (.TRUE. if solid/.FALSE. otherwise)

   ! Orientation and different cell related lengths
   INTEGER :: GHOST_FIRSTW(-3:3) = NSCARC_ZERO_INT             ! First internal ghost cell numbers for all faces
   INTEGER :: GHOST_LASTW(-3:3)  = NSCARC_ZERO_INT             ! Last internal ghost cell numbers for all faces
   INTEGER :: GHOST_FIRSTE(-3:3) = NSCARC_ZERO_INT             ! First external ghost cell numbers for all faces
   INTEGER :: GHOST_LASTE(-3:3)  = NSCARC_ZERO_INT             ! Last external ghost cell numbers for all faces
   INTEGER :: NX, NY, NZ                                       ! Number of grid cells in x-, y- and z-direction

   ! Number of discretizations and obstructions
   INTEGER :: N_DISCRET                                        ! Number of discretization types used
   INTEGER :: N_OBST                                           ! Number of obstructions
   INTEGER :: N_CELL_INDEX                                     ! Number of entries in CELL_INDEX array
   INTEGER :: N_CELLS                                          ! Number of cells in structured discretization
   INTEGER :: N_GHOST_ZONES                                    ! Number of adjacent ghost zones

   ! Different wall related lengths
   INTEGER :: N_WALL_CELLS                                     ! Number of wall cells
   INTEGER :: N_WALL_CELLS_EXT                                 ! Number of external wall cells
   INTEGER :: N_WALL_CELLS_INT                                 ! Number of internal wall cells
   INTEGER :: N_WALL_CELLS_LOCAL                               ! Number of local wall cells

END TYPE SCARC_LEVEL_TYPE

! 
! ---------- Sample sequence of used solvers in stack
! 
TYPE SCARC_SOLVER_TYPE

   CHARACTER(30) :: CNAME = 'NONE'                          ! Name of current solver

   ! Types of different solver components
   INTEGER :: TYPE_ACCURACY      = NSCARC_ACCURACY_ABSOLUTE    ! Default type of requested accuracy
   INTEGER :: TYPE_COARSE        = NSCARC_COARSE_DIRECT        ! Default type of coarse grid solver for multigrid method
   INTEGER :: TYPE_COARSENING    = NSCARC_UNDEF_INT            ! No default type of coarsening algorithm for AMG
   INTEGER :: TYPE_CYCLING       = NSCARC_CYCLING_V            ! Default type of cycling for multigrid method
   INTEGER :: TYPE_GRID          = NSCARC_GRID_STRUCTURED      ! Default type of discretization (structured/unstructured)
   INTEGER :: TYPE_EXCHANGE      = NSCARC_UNDEF_INT            ! No default type of data exchange
   INTEGER :: TYPE_INTERPOL      = NSCARC_INTERPOL_CONSTANT    ! Default type of interpolation method
   INTEGER :: TYPE_KRYLOV        = NSCARC_KRYLOV_CG            ! Default type of Krylov method (CG/BICG)
   INTEGER :: TYPE_LEVEL(0:2)    = NSCARC_UNDEF_INT            ! Default type of levels
   INTEGER :: TYPE_MATRIX        = NSCARC_MATRIX_COMPACT       ! Default type of storage for matrix
   INTEGER :: TYPE_METHOD        = NSCARC_METHOD_KRYLOV        ! Default type of ScaRC method
   INTEGER :: TYPE_MKL(0:10)     = NSCARC_UNDEF_INT            ! No default type of MKL for single levels
   INTEGER :: TYPE_MKL_PRECISION = NSCARC_PRECISION_DOUBLE     ! Default double precision MKL solver
   INTEGER :: TYPE_MULTIGRID     = NSCARC_MULTIGRID_GEOMETRIC  ! Default type of multigrid method (GMG/AMG)
   INTEGER :: TYPE_PARENT        = NSCARC_UNDEF_INT            ! No default type of parent (calling) solver
   INTEGER :: TYPE_PRECON        = NSCARC_UNDEF_INT            ! No default type of preconditioner for iterative solver
   INTEGER :: TYPE_RELAX         = NSCARC_UNDEF_INT            ! No default type of preconditioner for iterative solver
   INTEGER :: TYPE_SCOPE(0:2)    = NSCARC_SCOPE_LOCAL          ! Default scope types are local
   INTEGER :: TYPE_SMOOTH        = NSCARC_UNDEF_INT            ! No default type of smoother for multigrid method
   INTEGER :: TYPE_SOLVER        = NSCARC_SOLVER_MAIN          ! Default type of surrounding solver stage
   INTEGER :: TYPE_STAGE         = NSCARC_STAGE_ONE            ! Default type of surrounding solver stage
   INTEGER :: TYPE_STENCIL       = NSCARC_STENCIL_CONSTANT     ! Default type of storage for matrix
   INTEGER :: TYPE_TWOLEVEL      = NSCARC_TWOLEVEL_NONE        ! Default type of two-level method
   INTEGER :: TYPE_VECTOR        = NSCARC_UNDEF_INT            ! No default type of vector to point to


   ! References to different vectors which are needed for the current solver
   INTEGER :: X = NSCARC_UNDEF_INT                 ! Reference to local X-vector, double precision
   INTEGER :: B = NSCARC_UNDEF_INT                 ! Reference to local B-vector, double precision
   INTEGER :: D = NSCARC_UNDEF_INT                 ! Reference to local D-vector, double precision
   INTEGER :: R = NSCARC_UNDEF_INT                 ! Reference to local R-vector, double precision
   INTEGER :: V = NSCARC_UNDEF_INT                 ! Reference to local V-vector, double precision
   INTEGER :: Y = NSCARC_UNDEF_INT                 ! Reference to local Y-vector, double precision
   INTEGER :: Z = NSCARC_UNDEF_INT                 ! Reference to local Z-vector, double precision

   INTEGER :: X_FB = NSCARC_UNDEF_INT              ! Reference to local X-vector, single precision
   INTEGER :: B_FB = NSCARC_UNDEF_INT              ! Reference to local B-vector, single precision
   INTEGER :: R_FB = NSCARC_UNDEF_INT              ! Reference to local R-vector, single precision
   INTEGER :: V_FB = NSCARC_UNDEF_INT              ! Reference to local V-vector, single precision

#ifdef WITH_SCARC_DEBUG
   INTEGER :: E = NSCARC_UNDEF_INT                 ! Reference to local E-vector, double precision
#endif

   ! Converegence requirements for current solver
   INTEGER  :: NIT   = NSCARC_UNDEF_INT            ! Maximum iteration number
   INTEGER  :: ITE   = NSCARC_UNDEF_INT            ! Current iteration number
   REAL(EB) :: EPS   = NSCARC_UNDEF_REAL_EB        ! Required accuracy
   REAL(EB) :: RES   = NSCARC_UNDEF_REAL_EB        ! Current residual
   REAL(EB) :: RESIN = NSCARC_UNDEF_REAL_EB        ! Initial residual
   REAL(EB) :: ERR   = NSCARC_UNDEF_REAL_EB        ! Initial residual
   REAL(EB) :: OMEGA = NSCARC_UNDEF_REAL_EB        ! Relaxation parameter
   REAL(EB) :: CAPPA = NSCARC_UNDEF_REAL_EB        ! Convergence rate

END TYPE SCARC_SOLVER_TYPE

! 
! ---------- Stack type
! 
TYPE SCARC_STACK_TYPE
   TYPE (SCARC_SOLVER_TYPE), POINTER :: SOLVER     ! Type of current solver
END TYPE SCARC_STACK_TYPE

! 
! ---------- Administration other mesh data needed for the coupling of adjacent neighbors
! 
TYPE SCARC_OSCARC_TYPE

   TYPE (SCARC_LEVEL_TYPE), ALLOCATABLE, DIMENSION(:) :: LEVEL      ! Level related information

   REAL(EB) :: SEND_BUFFER_REAL0(1:NSCARC_MAX_BUFFER0)              ! Constant length send buffer for setup
   REAL(EB) :: RECV_BUFFER_REAL0(1:NSCARC_MAX_BUFFER0)              ! Constant length receive buffer for setup
   REAL(EB), ALLOCATABLE, DIMENSION (:) :: SEND_BUFFER_REAL         ! Real send buffer 
   REAL(EB), ALLOCATABLE, DIMENSION (:) :: RECV_BUFFER_REAL         ! Real receive buffer 

   INTEGER :: SEND_BUFFER_INT0(NSCARC_MAX_BUFFER0)                  ! Constant length send buffer for setup
   INTEGER :: RECV_BUFFER_INT0(NSCARC_MAX_BUFFER0)                  ! Constant length receive buffer for setup
   INTEGER,  ALLOCATABLE, DIMENSION (:) :: SEND_BUFFER_INT          ! Integer send buffer 
   INTEGER,  ALLOCATABLE, DIMENSION (:) :: RECV_BUFFER_INT          ! Integer receive buffer 

   INTEGER :: NLEN_MAX_BUFFER_LAYER1 = NSCARC_HUGE_INT              ! Length for single layer length exchange
   INTEGER :: NLEN_MAX_BUFFER_LAYER2 = NSCARC_HUGE_INT              ! Length for double layer length exchange
   INTEGER :: NLEN_MAX_BUFFER_LAYER4 = NSCARC_HUGE_INT              ! Length for double layer length exchange
   INTEGER :: NLEN_MAX_BUFFER_STENCIL = NSCARC_HUGE_INT             ! Length for stencil layer length exchange
   INTEGER :: NLEN_MAX_BUFFER_FULL = NSCARC_HUGE_INT                ! Length for max length exchange

END TYPE SCARC_OSCARC_TYPE

! 
! ---------- Measurement of CPU times
! 
TYPE SCARC_CPU_TYPE
   REAL(EB) :: OVERALL          = 0.0_EB                      ! complete time for ScaRC
   REAL(EB) :: SETUP            = 0.0_EB                      ! time for setup of requested ScaRC solver
   REAL(EB) :: SOLVER           = 0.0_EB                      ! time for solver 
   REAL(EB) :: ITERATION        = 0.0_EB                      ! time for krylov solver
   REAL(EB) :: MATVEC_PRODUCT   = 0.0_EB                      ! time for matrix vector multiplication
   REAL(EB) :: SCALAR_PRODUCT   = 0.0_EB                      ! time for scalar product
   REAL(EB) :: L2NORM           = 0.0_EB                      ! time for l2-norm
   REAL(EB) :: RELAXATION       = 0.0_EB                      ! time for relaxation
   REAL(EB) :: SMOOTHER         = 0.0_EB                      ! time for smoothing
   REAL(EB) :: COARSE           = 0.0_EB                      ! time for coarse grid solver
   REAL(EB) :: EXCHANGE         = 0.0_EB                      ! time for data exchange
   REAL(EB) :: BUFFER_PACKING   = 0.0_EB                      ! time for data exchange
   REAL(EB) :: BUFFER_UNPACKING = 0.0_EB                      ! time for data exchange
   INTEGER  :: N_TIMER          = 12                          ! total number of timers
END TYPE SCARC_CPU_TYPE

!
! ---------- Basic administration type for ScaRC-method
!
TYPE SCARC_TYPE

   TYPE (SCARC_OSCARC_TYPE), ALLOCATABLE, DIMENSION(:) :: OSCARC   ! ScaRC type on other mesh
   TYPE (SCARC_LEVEL_TYPE) , ALLOCATABLE, DIMENSION(:) :: LEVEL    ! Level related information

   REAL(EB) :: XS, XF, YS, YF, ZS, ZF                              ! x-, y- and z-bounds (corresponding to main prg)
   REAL(EB) :: RHS_END = 0.0_EB                                    ! Very last RHS entry, needed for matrix condensing

   INTEGER, ALLOCATABLE, DIMENSION(:) :: NEIGHBORS                 ! List of adjacent neighbors of whole mesh
   INTEGER :: N_NEIGHBORS = 0                                      ! Number of adjacent neighbors of whole mesh
   INTEGER :: NC = 0                                               ! Total number of cells on that mesh
   INTEGER :: IBAR, JBAR, KBAR                                     ! Number of cells (corresponding to main prg)

END TYPE SCARC_TYPE

END MODULE SCARC_TYPES


! ================================================================================================================
! MODULE 'SCARC_POINTERS' : 
! Collection of different pointers to specify the different meshes, grid levels, discretizations and matrices
! ================================================================================================================
MODULE SCARC_POINTERS
USE MESH_VARIABLES
USE SCARC_TYPES
IMPLICIT NONE

TYPE (MESH_TYPE),  POINTER :: M=>NULL()
TYPE (OMESH_TYPE), POINTER :: OM=>NULL()
TYPE (WALL_TYPE),  POINTER :: MWC=>NULL()
TYPE (EXTERNAL_WALL_TYPE), POINTER :: EWC=>NULL()

TYPE (SCARC_TYPE),          POINTER :: S=>NULL()
TYPE (SCARC_OSCARC_TYPE),   POINTER :: OS=>NULL()
TYPE (SCARC_LEVEL_TYPE),    POINTER :: L=>NULL(), LF=>NULL(), LC=>NULL(), OL=>NULL(), OLF=>NULL(), OLC=>NULL()
TYPE (SCARC_GRID_TYPE),     POINTER :: G=>NULL(), GF=>NULL(), GC=>NULL(), OG=>NULL(), OGF=>NULL(), OGC=>NULL()
TYPE (SCARC_FACE_TYPE),     POINTER :: F=>NULL(), FF=>NULL(), FC=>NULL()
TYPE (SCARC_OBST_TYPE),     POINTER :: OB=>NULL()

TYPE (SCARC_SUBDIVISION_TYPE), POINTER :: SUB=>NULL()

TYPE (SCARC_WALL_TYPE), POINTER :: GWC=>NULL()
TYPE (SCARC_WALL_TYPE), DIMENSION(:), POINTER :: W=>NULL(), WF=>NULL(), WC=>NULL()

TYPE (SCARC_SOLVER_TYPE), POINTER :: SV=>NULL(), SVP=>NULL()
TYPE (SCARC_STAGE_TYPE),  POINTER :: ST=>NULL(), STP=>NULL()

TYPE (SCARC_FFT_TYPE), POINTER :: FFT=>NULL()
#ifdef WITH_SCARC_MGM
TYPE (SCARC_MGM_TYPE), POINTER :: MGM=>NULL()
#endif

TYPE (SCARC_MATRIX_BANDWISE_TYPE), POINTER :: AB=>NULL(), OAB=>NULL()

TYPE (SCARC_MATRIX_COMPACT_TYPE), POINTER :: A=>NULL(), OA=>NULL(), AF=>NULL(), AC=>NULL(), OAF=>NULL(), OAC=>NULL()
TYPE (SCARC_MATRIX_COMPACT_TYPE), POINTER :: P=>NULL(), OP=>NULL(), PF=>NULL(), PC=>NULL(), OPF=>NULL(), OPC=>NULL()
TYPE (SCARC_MATRIX_COMPACT_TYPE), POINTER :: R=>NULL(), OR=>NULL(), RF=>NULL(), RC=>NULL(), ORF=>NULL(), ORC=>NULL()
TYPE (SCARC_MATRIX_COMPACT_TYPE), POINTER :: C=>NULL(), OC=>NULL(), CF=>NULL(), CC=>NULL(), OCF=>NULL(), OCC=>NULL()
TYPE (SCARC_MATRIX_COMPACT_TYPE), POINTER :: Z=>NULL(), OZ=>NULL(), ZF=>NULL(), ZC=>NULL(), OZF=>NULL(), OZC=>NULL()
TYPE (SCARC_MATRIX_COMPACT_TYPE), POINTER :: AP=>NULL(), APF=>NULL(), APC=>NULL(), OAP=>NULL()

TYPE (SCARC_MATRIX_COMPACT_CONDENSED_TYPE),  POINTER :: ACO =>NULL()
TYPE (SCARC_MATRIX_BANDWISE_CONDENSED_TYPE), POINTER :: ABCO=>NULL()

TYPE (SCARC_MULTIGRID_TYPE), POINTER :: MG =>NULL()

REAL(EB), POINTER, DIMENSION(:) :: XCOR=>NULL(), YCOR=>NULL(), ZCOR=>NULL()
REAL(EB), POINTER, DIMENSION(:) :: XMID=>NULL(), YMID=>NULL(), ZMID=>NULL()

REAL(EB), POINTER, DIMENSION(:) :: VC=>NULL(), VF=>NULL(), V1=>NULL(), V2=>NULL()

REAL(EB), POINTER, DIMENSION(:,:,:) :: HP=>NULL(), HVC=>NULL(), PRHS=>NULL()
REAL(EB), POINTER, DIMENSION(:,:,:) :: UU=>NULL(), VV=>NULL(), WW=>NULL()

REAL(EB), POINTER, DIMENSION(:) ::  RECV_BUFFER_REAL
INTEGER,  POINTER, DIMENSION(:) ::  RECV_BUFFER_INT

INTEGER, POINTER :: I_MIN_R, I_MAX_R, I_MIN_S, I_MAX_S
INTEGER, POINTER :: J_MIN_R, J_MAX_R, J_MIN_S, J_MAX_S
INTEGER, POINTER :: K_MIN_R, K_MAX_R, K_MIN_S, K_MAX_S

INTEGER, POINTER, DIMENSION(:) :: IOR_R, IOR_S
INTEGER, POINTER, DIMENSION(:) :: IIO_R, IIO_S
INTEGER, POINTER, DIMENSION(:) :: JJO_R, JJO_S
INTEGER, POINTER, DIMENSION(:) :: KKO_R, KKO_S

#ifdef WITH_MKL
TYPE (SCARC_MKL_TYPE), POINTER :: MKL=>NULL()
TYPE (SCARC_MATRIX_COMPACT_TYPE), POINTER :: AS=>NULL(), ACS=>NULL(), AFS=>NULL()
TYPE (SCARC_MATRIX_COMPACT_TYPE), POINTER :: OAS=>NULL(), OACS=>NULL(), OAFS=>NULL()
REAL(FB), DIMENSION(:), POINTER :: V1_FB=>NULL(), V2_FB=>NULL()
#endif

END MODULE SCARC_POINTERS


! ================================================================================================================
! ROUTINE 'SCRC'
! ================================================================================================================
MODULE SCRC

USE PRECISION_PARAMETERS
USE GLOBAL_CONSTANTS
USE MESH_VARIABLES
USE MEMORY_FUNCTIONS, ONLY: CHKMEMERR
USE COMP_FUNCTIONS, ONLY: CURRENT_TIME, GET_FILE_NUMBER, SHUTDOWN
USE MPI
USE SCARC_GLOBAL_CONSTANTS
USE SCARC_VARIABLES
USE SCARC_TYPES
USE SCARC_POINTERS

#ifdef WITH_MKL
USE MKL_PARDISO
USE MKL_CLUSTER_SPARSE_SOLVER
#endif

IMPLICIT NONE

! 
! ---------- Public variables
! 
PUBLIC :: SCARC_SETUP, SCARC_METHOD, SCARC_SOLVER, SCARC_GRID, SCARC_TWOLEVEL      
PUBLIC :: SCARC_ACCURACY, SCARC_CAPPA, SCARC_ITERATIONS, SCARC_RESIDUAL, SCARC_MATRIX
PUBLIC :: SCARC_DUMP_TIMERS, SCARC_ERROR_FILE, SCARC_MKL_PRECISION

PUBLIC :: SCARC_KRYLOV, SCARC_KRYLOV_ACCURACY, SCARC_KRYLOV_ITERATIONS, SCARC_KRYLOV_INTERPOL
PUBLIC :: SCARC_COARSE, SCARC_COARSE_ACCURACY, SCARC_COARSE_ITERATIONS, SCARC_COARSE_OMEGA, SCARC_COARSE_LEVEL     
PUBLIC :: SCARC_PRECON, SCARC_PRECON_ACCURACY, SCARC_PRECON_ITERATIONS, SCARC_PRECON_OMEGA, SCARC_PRECON_SCOPE 
PUBLIC :: SCARC_SMOOTH, SCARC_SMOOTH_ACCURACY, SCARC_SMOOTH_ITERATIONS, SCARC_SMOOTH_OMEGA, SCARC_SMOOTH_SCOPE 

PUBLIC :: SCARC_MULTIGRID, SCARC_MULTIGRID_ACCURACY, SCARC_MULTIGRID_ITERATIONS, SCARC_MULTIGRID_INTERPOL
PUBLIC :: SCARC_MULTIGRID_CYCLE, SCARC_MULTIGRID_COARSENING, SCARC_MULTIGRID_LEVEL
PUBLIC :: SCARC_MULTIGRID_PRESMOOTH, SCARC_MULTIGRID_POSTSMOOTH


! 
! ---------- Type declarations
! 
TYPE (SCARC_TYPE)        , SAVE, DIMENSION(:), ALLOCATABLE, TARGET :: SCARC
TYPE (SCARC_STACK_TYPE)  , SAVE, DIMENSION(:), ALLOCATABLE, TARGET :: STACK
TYPE (SCARC_CPU_TYPE)    , SAVE, DIMENSION(:), ALLOCATABLE, TARGET :: CPU

TYPE (SCARC_SOLVER_TYPE) , SAVE, TARGET :: MAIN_CG, MAIN_CG_STRUCTURED, MAIN_CG_UNSTRUCTURED
TYPE (SCARC_SOLVER_TYPE) , SAVE, TARGET :: MAIN_GMG
TYPE (SCARC_SOLVER_TYPE) , SAVE, TARGET :: COARSE_KRYLOV
#ifdef WITH_MKL
TYPE (SCARC_SOLVER_TYPE) , SAVE, TARGET :: MAIN_LU
TYPE (SCARC_SOLVER_TYPE) , SAVE, TARGET :: COARSE_CLUSTER, COARSE_PARDISO
#endif
TYPE (SCARC_SOLVER_TYPE) , SAVE, TARGET :: PRECON_JAC, PRECON_SSOR, PRECON_ILU, PRECON_FFT, PRECON_FFTO, PRECON_GMG
TYPE (SCARC_SOLVER_TYPE) , SAVE, TARGET :: PRECON_MJAC, PRECON_MGS, PRECON_MSGS, PRECON_MSOR, PRECON_MSSOR
TYPE (SCARC_SOLVER_TYPE) , SAVE, TARGET :: SMOOTH_JACOBI, SMOOTH_SSOR, SMOOTH_ILU, SMOOTH_FFT, SMOOTH_FFTO
TYPE (SCARC_SOLVER_TYPE) , SAVE, TARGET :: SMOOTH_MJAC, SMOOTH_MGS, SMOOTH_MSGS, SMOOTH_MSOR, SMOOTH_MSSOR
#ifdef WITH_MKL
TYPE (SCARC_SOLVER_TYPE) , SAVE, TARGET :: PRECON_MKL, SMOOTH_MKL
#endif
TYPE (SCARC_MESSAGE_TYPE), SAVE :: MSG

TYPE (SCARC_SUBDIVISION_TYPE), SAVE, TARGET :: SUBDIVISION

CONTAINS


! ------------------------------------------------------------------------------------------------
! Initialize ScaRC structures based on SCARC-input parameters from &PRES namelist
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP
REAL(EB) :: TNOW

TNOW = CURRENT_TIME()

!
! Setup messaging/debugging mechanisms and time measurements
!
CALL SCARC_SETUP_MESSAGING
CALL SCARC_SETUP_TIMING

!
! Parse all ScaRC parameters which have been read in read.f90
!
CALL SCARC_PARSE_INPUT                      ; IF (STOP_STATUS==SETUP_STOP) RETURN

!
! Setup different basic components of ScaRC solver
!
CALL SCARC_SETUP_LEVELS                     ; IF (STOP_STATUS==SETUP_STOP) RETURN
CALL SCARC_SETUP_TYPES                      ; IF (STOP_STATUS==SETUP_STOP) RETURN
CALL SCARC_SETUP_GRIDS                      ; IF (STOP_STATUS==SETUP_STOP) RETURN
CALL SCARC_SETUP_GLOBALS                    ; IF (STOP_STATUS==SETUP_STOP) RETURN
CALL SCARC_SETUP_NEIGHBORS                  ; IF (STOP_STATUS==SETUP_STOP) RETURN
CALL SCARC_SETUP_FACES                      ; IF (STOP_STATUS==SETUP_STOP) RETURN
CALL SCARC_SETUP_SUBDIVISION                ; IF (STOP_STATUS==SETUP_STOP) RETURN

!
! Setup wall information according to specified discretization type/method
!
IF (HAS_MULTI_GRIDS) THEN
   CALL SCARC_ASSIGN_GRID_TYPE (NSCARC_GRID_STRUCTURED)
   CALL SCARC_SETUP_WALLS                   ; IF (STOP_STATUS==SETUP_STOP) RETURN
   CALL SCARC_ASSIGN_GRID_TYPE (NSCARC_GRID_UNSTRUCTURED)
   CALL SCARC_SETUP_WALLS                   ; IF (STOP_STATUS==SETUP_STOP) RETURN
ELSE
   CALL SCARC_ASSIGN_GRID_TYPE (TYPE_GRID)
   CALL SCARC_SETUP_WALLS                   ; IF (STOP_STATUS==SETUP_STOP) RETURN
ENDIF

!
! Setup information for data exchanges, matrix systems, aggregation zones, used methods and vectors
!
CALL SCARC_SETUP_EXCHANGES                  ; IF (STOP_STATUS==SETUP_STOP) RETURN
CALL SCARC_SETUP_SYSTEMS                    ; IF (STOP_STATUS==SETUP_STOP) RETURN

IF (IS_AMG .OR. IS_CG_AMG .OR. IS_CG_SAMG) THEN
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'SETUP_SAMG'
#endif
   CALL SCARC_SETUP_SAMG                    ; IF (STOP_STATUS==SETUP_STOP) RETURN
ENDIF

CALL SCARC_SETUP_METHODS                    ; IF (STOP_STATUS==SETUP_STOP) RETURN
CALL SCARC_SETUP_VECTORS                    ; IF (STOP_STATUS==SETUP_STOP) RETURN

!
! Perform some error statistics for pressure if corresponding flag is set
!
#ifdef WITH_SCARC_STANDALONE
CALL SCARC_SETUP_PRESSURE                   ; IF (STOP_STATUS==SETUP_STOP) RETURN
#endif


CPU(MYID)%SETUP   = CPU(MYID)%SETUP   + CURRENT_TIME() - TNOW
CPU(MYID)%OVERALL = CPU(MYID)%OVERALL + CURRENT_TIME() - TNOW
END SUBROUTINE SCARC_SETUP


! ------------------------------------------------------------------------------------------------
! Setup time measurements
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_TIMING
ALLOCATE (CPU(0:N_MPI_PROCESSES-1), STAT = IERROR)
CALL CHKMEMERR ('SCARC_SETUP_TIMING', 'CPU', IERROR)
END SUBROUTINE SCARC_SETUP_TIMING


! ------------------------------------------------------------------------------------------------
! Setup debug file if requested
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_MESSAGING
#if defined(WITH_SCARC_VERBOSE) || defined(WITH_SCARC_DEBUG)
INTEGER :: NM, LASTID
#endif

IF (SCARC_ERROR_FILE) HAS_CSV_DUMP = .TRUE.

IF (NMESHES == 1) THEN
   CFORM_INT  = "(6I4)"
   CFORM_REAL = "(6E12.4)"
ELSE IF (NMESHES == 2) THEN
   IF (MYID == 0) THEN
      CFORM_INT = "(3I4)"
      CFORM_REAL = "(3E12.4)"
   ELSE
      CFORM_INT = "(3I4)"
      CFORM_REAL = "(3E12.4)"
   ENDIF
ELSE IF (NMESHES == 4) THEN
   IF (MYID == 0 .OR. MYID == 2) THEN
      CFORM_INT = "(7I4)"
      CFORM_REAL = "(7E12.4)"
   ELSE IF (MYID == 1 .OR. MYID == 3) THEN
      CFORM_INT = "(6I4)"
      CFORM_REAL = "(6E12.4)"
   ENDIF
ENDIF

!
! If requested, open file for CSV-information about convergence of different solvers
!
IF (HAS_CSV_DUMP) THEN
   IF (MYID == 0) THEN
      WRITE (MSG%FILE_STAT, '(A,A)') TRIM(CHID),'_scarc.csv'
      MSG%LU_STAT = GET_FILE_NUMBER()
      OPEN (MSG%LU_STAT, FILE=MSG%FILE_STAT)
      WRITE(MSG%LU_STAT,*) '  #Pres,   Stack,  #ScaRC,     #CG,     #MG,   Level, #Smooth, SmoType, ', &
                           '#Coarse,     #LU,    Residual,   Cappa'
   ENDIF
ENDIF

#ifdef WITH_SCARC_VERBOSE
!
! If requested, open file for log-information
!
LASTID = -NSCARC_HUGE_INT
DO NM=LOWER_MESH_INDEX, UPPER_MESH_INDEX
   IF (MYID == LASTID) CYCLE
   WRITE (MSG%FILE_VERBOSE, '(A,A,i3.3)') TRIM(CHID),'.log',MYID+1
   MSG%LU_VERBOSE = GET_FILE_NUMBER()
   OPEN (MSG%LU_VERBOSE, FILE=MSG%FILE_VERBOSE, ACTION = 'readwrite')
   LASTID = MYID
ENDDO
#endif

#ifdef WITH_SCARC_DEBUG
!
! If requested, open file for debug messages
!
LASTID = -NSCARC_HUGE_INT
DO NM=LOWER_MESH_INDEX, UPPER_MESH_INDEX
   IF (MYID == LASTID) CYCLE
   WRITE (MSG%FILE_DEBUG, '(A,A,i3.3)') TRIM(CHID),'.debug',MYID+1
   MSG%LU_DEBUG = GET_FILE_NUMBER()
   OPEN (MSG%LU_DEBUG, FILE=MSG%FILE_DEBUG, ACTION = 'readwrite')
   LASTID = MYID
ENDDO
#endif

END SUBROUTINE SCARC_SETUP_MESSAGING


! ------------------------------------------------------------------------------------------------
! Shutdown ScaRC with error message
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SHUTDOWN(NERROR, CPARAM, NPARAM)
CHARACTER(*), INTENT(IN) :: CPARAM
INTEGER, INTENT(IN) :: NERROR, NPARAM
CHARACTER(80) :: CERROR

!
! Assign error message according to specified error
!
SELECT CASE (NERROR)
   CASE (NSCARC_ERROR_PARSE_INPUT)
      CERROR = 'Wrong input parameter'
   CASE (NSCARC_ERROR_MKL_INTERNAL)
      CERROR = 'The following MKL error was detected'
   CASE (NSCARC_ERROR_MKL_PARDISO)
      CERROR = 'MKL Library compile flag not defined, Pardiso solver not available'
   CASE (NSCARC_ERROR_MKL_CLUSTER)
      CERROR = 'MKL Library compile flag not defined, Cluster_Sparse_Solver not available'
   CASE (NSCARC_ERROR_MKL_STORAGE)
      CERROR = 'Wrong matrix storage scheme for MKL solvers, only COMPACT storage available'
   CASE (NSCARC_ERROR_GRID_NUMBER)
      CERROR = 'Number not divisable by 2'
   CASE (NSCARC_ERROR_GRID_NUMBERX)
      CERROR = 'Number of cells not divisable by 2 in x-direction, NC'
   CASE (NSCARC_ERROR_GRID_NUMBERY)
      CERROR = 'Number of cells not divisable by 2 in y-direction, NC'
   CASE (NSCARC_ERROR_GRID_NUMBERZ)
      CERROR = 'Number of cells not divisable by 2 in z-direction, NC'
   CASE (NSCARC_ERROR_BOUNDARY_SUM)
      CERROR = 'Wrong boundary sum for IOR'
   CASE (NSCARC_ERROR_DIRECT_NOMKL)
      CERROR = 'Direct coarse grid solver is only working in combination with MKL'
   CASE (NSCARC_ERROR_BOUNDARY_TYPE)
      CERROR = 'Wrong boundary type'
   CASE (NSCARC_ERROR_NEIGHBOR_TYPE)
      CERROR = 'Wrong neighbor'
   CASE (NSCARC_ERROR_NEIGHBOR_NUMBER)
      CERROR = 'More than 20 neighbors along one face not allowed'
   CASE (NSCARC_ERROR_MATRIX_ALLOCATION)
      CERROR = 'Wrong specifier during allocation or deallocation of  matrix'
   CASE (NSCARC_ERROR_MATRIX_SUBDIAG)
      CERROR = 'Subdiagonal missing for system matrix'
   CASE (NSCARC_ERROR_MATRIX_SYMMETRY)
      CERROR = 'Matrix not symmetric for mesh'
   CASE (NSCARC_ERROR_MATRIX_SETUP)
      CERROR = 'Matrix setup failed for level type'
   CASE (NSCARC_ERROR_MATRIX_SIZE)
      CERROR = 'Matrix reducing failed because new length is too big for matrix'
   CASE (NSCARC_ERROR_MATRIX_COPY)
      CERROR = 'Matrix copy failed due to too already existing array'
   CASE (NSCARC_ERROR_STENCIL)
      CERROR = 'Wrong type for matrix stencil - only constant or variable allowed'
   CASE (NSCARC_ERROR_STACK_SOLVER)
      CERROR = 'Wrong number of solvers in stack'
   CASE (NSCARC_ERROR_STACK_MESSAGE)
      CERROR = 'Too many messages in calling stack'
   CASE (NSCARC_ERROR_MULTIGRID_LEVEL)
      CERROR = 'Wrong level for multigrid method'
   CASE (NSCARC_ERROR_GRID_RESOLUTION)
      CERROR = 'Wrong grid resolution at IOR'
   CASE (NSCARC_ERROR_GRID_INDEX)
      CERROR = 'Wrong index for J'
   CASE (NSCARC_ERROR_RHS_SETUP)
      CERROR = 'Wrong level for presetting RHS'
   CASE (NSCARC_ERROR_AGGREGATION)
      CERROR = 'Error during aggregation process in algebraic multigrid'
   CASE (NSCARC_ERROR_VECTOR_LENGTH)
      CERROR = 'Inconsistent length for vector allocation'
   CASE (NSCARC_ERROR_POC_STOP)
      CERROR = 'Only one call of solver due to proof of concept'
   CASE (NSCARC_ERROR_BICG_DISABLED)
      CERROR = 'Krylov solver BICG temporarily disabled'
   CASE (NSCARC_ERROR_EXCHANGE_RECV)
      CERROR = 'Wrong receive exchange structure'
   CASE (NSCARC_ERROR_EXCHANGE_SEND)
      CERROR = 'Wrong send exchange structure'
END SELECT

!
! Specify more detailed information if available
!
IF (CPARAM /= SCARC_NONE) THEN
   IF (MYID == 0) WRITE(LU_ERR,1000)  CERROR, CPARAM, TRIM(CHID)
ELSE IF (NPARAM /= NSCARC_NONE) THEN
   IF (MYID == 0) WRITE(LU_ERR,2000)  CERROR, NPARAM, TRIM(CHID)
ELSE
   IF (MYID == 0) WRITE(LU_ERR,3000)  CERROR, TRIM(CHID)
ENDIF

!
! Also print verbose message if enabled
!
#ifdef WITH_SCARC_VERBOSE
IF (CPARAM /= SCARC_NONE) THEN
   WRITE(MSG%LU_VERBOSE,1000)  CERROR, CPARAM, TRIM(CHID)
ELSE IF (NPARAM /= NSCARC_NONE) THEN
   WRITE(MSG%LU_VERBOSE,2000)  CERROR, NPARAM, TRIM(CHID)
ELSE
   WRITE(MSG%LU_VERBOSE,3000)  CERROR, TRIM(CHID)
ENDIF
CLOSE(MSG%LU_VERBOSE)
#endif

#ifdef WITH_SCARC_DEBUG
CLOSE(MSG%LU_DEBUG)
#endif

STOP_STATUS = SETUP_STOP
RETURN

1000 FORMAT('Stop in ScaRC-solver: ', A,' : ',   A, ' (CHID: ',A,')' )
2000 FORMAT('Stop in ScaRC-solver: ', A,' : ', I12, ' (CHID: ',A,')' )
3000 FORMAT('Stop in ScaRC-solver: ', A, ' (CHID: ',A,')' )
END SUBROUTINE SCARC_SHUTDOWN


! ----------------------------------------------------------------------------------------------------
! Determine types of input parameters
! ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_PARSE_INPUT

ITERATE_PRESSURE = .TRUE.  ! Although there is no need to do pressure iterations to drive down 
                           ! velocity error leave it .TRUE. to write out velocity error diagnostics

!
! ------------- Set type of discretization
!
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

!
! ------------ Set type of matrix storage (COMPACT/BANDWISE)
!
SELECT CASE (TRIM(SCARC_MATRIX))
   CASE ('COMPACT')
      TYPE_MATRIX = NSCARC_MATRIX_COMPACT
   CASE ('BANDWISE')
      TYPE_MATRIX = NSCARC_MATRIX_BANDWISE
   CASE DEFAULT
      CALL SCARC_SHUTDOWN(NSCARC_ERROR_PARSE_INPUT, SCARC_MATRIX, NSCARC_NONE)
END SELECT

!
! ------------ Set type of matrix stencil (CONSTANT/VARIABLE)
!
SELECT CASE (TRIM(SCARC_STENCIL))
   CASE ('CONSTANT')
      TYPE_STENCIL = NSCARC_STENCIL_CONSTANT
   CASE ('VARIABLE')
      TYPE_STENCIL = NSCARC_STENCIL_VARIABLE
   CASE ('CONSTANT_LOCAL')
      TYPE_STENCIL = NSCARC_STENCIL_CONSTANT_LOCAL
   CASE ('VARIABLE_LOCAL')
      TYPE_STENCIL = NSCARC_STENCIL_VARIABLE_LOCAL
   CASE DEFAULT
      CALL SCARC_SHUTDOWN(NSCARC_ERROR_PARSE_INPUT, SCARC_STENCIL, NSCARC_NONE)
END SELECT

!
! ------------ Set type of global solver
!
SELECT CASE (TRIM(SCARC_METHOD))

   ! ------------------------- Global Krylov solver ----------------------------------
   CASE ('KRYLOV')

      TYPE_METHOD = NSCARC_METHOD_KRYLOV

      ! Set type of Krylov-method (CG/BICG)
      SELECT CASE (TRIM(SCARC_KRYLOV))
         CASE ('CG')
            TYPE_KRYLOV = NSCARC_KRYLOV_CG
         CASE ('BICG')
            CALL SCARC_SHUTDOWN(NSCARC_ERROR_BICG_DISABLED, SCARC_KRYLOV, NSCARC_NONE)    ! only temporarily
            TYPE_KRYLOV = NSCARC_KRYLOV_BICG
         CASE DEFAULT
            CALL SCARC_SHUTDOWN(NSCARC_ERROR_PARSE_INPUT, SCARC_KRYLOV, NSCARC_NONE)
      END SELECT

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
         CASE ('SAMG')
            TYPE_TWOLEVEL = NSCARC_TWOLEVEL_SAMG
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
         CASE DEFAULT
            CALL SCARC_SHUTDOWN(NSCARC_ERROR_PARSE_INPUT, SCARC_KRYLOV_INTERPOL, NSCARC_NONE)
      END SELECT

      ! Set type of preconditioner (JACOBI/SSOR/MGS/MSGS/MSOR/MSSOR/ILU/FFT/GMG/PARDISO/CLUSTER)
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
         CASE ('GMG')                                       ! Geometric multigrid preconditioner
            TYPE_PRECON = NSCARC_RELAX_GMG
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
            TYPE_SCOPE(0) = NSCARC_SCOPE_LOCAL
#else
            CALL SCARC_SHUTDOWN(NSCARC_ERROR_MKL_PARDISO, SCARC_NONE, NSCARC_NONE)
#endif
         CASE ('CLUSTER')                            !  LU-preconditioner based on MKL Cluster_Sparse_Solver
#ifdef WITH_MKL
            TYPE_PRECON   = NSCARC_RELAX_MKL
            TYPE_MKL(0)   = NSCARC_MKL_GLOBAL
            TYPE_SCOPE(0) = NSCARC_SCOPE_GLOBAL
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
         CASE ('ALGEBRAIC')
            TYPE_MULTIGRID = NSCARC_MULTIGRID_ALGEBRAIC
         CASE DEFAULT
            CALL SCARC_SHUTDOWN(NSCARC_ERROR_PARSE_INPUT, SCARC_MULTIGRID, NSCARC_NONE)
      END SELECT

      ! Set type of smoother (JACOBI/SGS/SSOR/MSSOR/ILU/PARDISO/CLUSTER)
      SELECT CASE (TRIM(SCARC_SMOOTH))                        ! use same parameters as for preconditioner
         CASE ('JACOBI')
            TYPE_SMOOTH = NSCARC_RELAX_JAC
         CASE ('SSOR')
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

#ifdef WITH_SCARC_MGM
   ! ------------------------- McKenny-Greengard-Mayo solver -------------------------
   CASE ('MGM')

      ! Just preset some values for proof of concept
      HAS_MULTI_GRIDS = .TRUE.

      TYPE_METHOD   = NSCARC_METHOD_MGM
      TYPE_KRYLOV   = NSCARC_KRYLOV_CG
      TYPE_TWOLEVEL = NSCARC_TWOLEVEL_NONE
      TYPE_PRECON   = NSCARC_RELAX_FFT
      TYPE_SCOPE(1) = NSCARC_SCOPE_LOCAL
#endif

   CASE DEFAULT
      CALL SCARC_SHUTDOWN(NSCARC_ERROR_PARSE_INPUT, SCARC_METHOD, NSCARC_NONE)

END SELECT

!
! If a multigrid solver is used (either as main solver or as preconditioner)
! set types for multigrid, coarse grid solver and cycling pattern
!
IF (TYPE_METHOD == NSCARC_METHOD_MULTIGRID .OR. TYPE_PRECON == NSCARC_RELAX_GMG .OR. TYPE_PRECON == NSCARC_RELAX_AMG) THEN

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

   ! Set type of coarsening strategy (AMG only)
   SELECT CASE (TRIM(SCARC_MULTIGRID_COARSENING))
      CASE ('GMG')
         TYPE_COARSENING = NSCARC_COARSENING_GMG
      CASE ('SAMG')
         TYPE_COARSENING = NSCARC_COARSENING_SAMG
      CASE DEFAULT
         CALL SCARC_SHUTDOWN(NSCARC_ERROR_PARSE_INPUT, SCARC_MULTIGRID_COARSENING, NSCARC_NONE)
   END SELECT

   ! Set type of interpolation 
   SELECT CASE (TRIM(SCARC_MULTIGRID_INTERPOL))
      CASE ('STANDARD')
         TYPE_INTERPOL = NSCARC_INTERPOL_STANDARD
      CASE ('CONSTANT')
         TYPE_INTERPOL = NSCARC_INTERPOL_CONSTANT
      CASE ('BILINEAR')
         TYPE_INTERPOL = NSCARC_INTERPOL_BILINEAR
      CASE ('CLASSICAL')
         TYPE_INTERPOL = NSCARC_INTERPOL_CLASSICAL
      CASE ('DIRECT')
         TYPE_INTERPOL = NSCARC_INTERPOL_DIRECT
      CASE DEFAULT
         CALL SCARC_SHUTDOWN(NSCARC_ERROR_PARSE_INPUT, SCARC_MULTIGRID_INTERPOL, NSCARC_NONE)
   END SELECT

ENDIF


! Set type of coarse grid solver
SELECT CASE (TRIM(SCARC_COARSE))
   CASE ('ITERATIVE')
      TYPE_COARSE = NSCARC_COARSE_ITERATIVE
      TYPE_KRYLOV = NSCARC_KRYLOV_CG
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


!
! Set type of accuracy (ABSOLUTE/RELATIVE)
!
SELECT CASE (TRIM(SCARC_ACCURACY))
   CASE ('ABSOLUTE')
      TYPE_ACCURACY = NSCARC_ACCURACY_ABSOLUTE
   CASE ('RELATIVE')
      TYPE_ACCURACY = NSCARC_ACCURACY_RELATIVE
   CASE DEFAULT
      CALL SCARC_SHUTDOWN(NSCARC_ERROR_PARSE_INPUT, SCARC_ACCURACY, NSCARC_NONE)
END SELECT

!
! Set type of precision for MKL solver (SINGLE/DOUBLE)
!
SELECT CASE (TRIM(SCARC_MKL_PRECISION))
   CASE ('SINGLE')
      TYPE_MKL_PRECISION = NSCARC_PRECISION_SINGLE
   CASE ('DOUBLE')
      TYPE_MKL_PRECISION = NSCARC_PRECISION_DOUBLE
   CASE DEFAULT
      CALL SCARC_SHUTDOWN(NSCARC_ERROR_PARSE_INPUT, SCARC_MKL_PRECISION, NSCARC_NONE)
END SELECT

!
! -------- define some logical variables - just for notational convenience
!
IS_STRUCTURED   = (TYPE_GRID == NSCARC_GRID_STRUCTURED)
IS_UNSTRUCTURED = (TYPE_GRID == NSCARC_GRID_UNSTRUCTURED)

IS_CG     = (TYPE_METHOD == NSCARC_METHOD_KRYLOV)
IS_CG_GMG = (TYPE_PRECON == NSCARC_RELAX_GMG) .AND. (TYPE_MULTIGRID == NSCARC_MULTIGRID_GEOMETRIC) .AND. IS_CG
IS_CG_AMG = (TYPE_PRECON == NSCARC_RELAX_AMG) .AND. (TYPE_MULTIGRID == NSCARC_MULTIGRID_ALGEBRAIC) .AND. IS_CG

IS_MG  = (TYPE_METHOD == NSCARC_METHOD_MULTIGRID)
IS_GMG = (TYPE_MULTIGRID == NSCARC_MULTIGRID_GEOMETRIC) .AND. IS_MG
IS_AMG = (TYPE_MULTIGRID == NSCARC_MULTIGRID_ALGEBRAIC) .AND. IS_MG

HAS_TWO_LEVELS   = IS_CG .AND. (TYPE_PRECON /= NSCARC_RELAX_GMG) .AND. (TYPE_TWOLEVEL > NSCARC_TWOLEVEL_NONE)
HAS_MULTI_LEVELS = IS_MG .OR. IS_CG_GMG .OR. HAS_TWO_LEVELS 

IS_CG_ADD    = HAS_TWO_LEVELS .AND. (TYPE_TWOLEVEL == NSCARC_TWOLEVEL_ADD)
IS_CG_MUL    = HAS_TWO_LEVELS .AND. (TYPE_TWOLEVEL == NSCARC_TWOLEVEL_MUL)
IS_CG_COARSE = HAS_TWO_LEVELS .AND. (TYPE_TWOLEVEL == NSCARC_TWOLEVEL_COARSE)
IS_CG_MACRO  = HAS_TWO_LEVELS .AND. (TYPE_TWOLEVEL == NSCARC_TWOLEVEL_MACRO)
IS_CG_SAMG   = HAS_TWO_LEVELS .AND. (TYPE_TWOLEVEL == NSCARC_TWOLEVEL_SAMG)

IS_FFT = (TYPE_PRECON == NSCARC_RELAX_FFT)  .OR. (TYPE_SMOOTH == NSCARC_RELAX_FFT)
IS_FFTO= (TYPE_PRECON == NSCARC_RELAX_FFTO) .OR. (TYPE_SMOOTH == NSCARC_RELAX_FFTO)
IS_MKL = (TYPE_PRECON >= NSCARC_RELAX_MKL)  .OR. &
         (TYPE_SMOOTH >= NSCARC_RELAX_MKL)  .OR. &
         (HAS_MULTI_LEVELS .AND. TYPE_COARSE == NSCARC_COARSE_DIRECT)

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'TYPE_TWOLEVEL=',TYPE_TWOLEVEL
WRITE(MSG%LU_DEBUG,*) 'IS_CG_SAMG=',IS_CG_SAMG
#endif

#ifdef WITH_SCARC_MGM
IS_MGM = TYPE_METHOD == NSCARC_METHOD_MGM
#endif

END SUBROUTINE SCARC_PARSE_INPUT


! ------------------------------------------------------------------------------------------------
! Determine number of grid levels  (1 for CG/BICG-method, NLEVEL for MG-method)
! Note: NLEVEL_MIN corresponds to finest grid resolution, NLEVEL_MAX to coarsest resolution
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_LEVELS
#ifdef WITH_MKL
INTEGER :: NL
#endif

SELECT_SCARC_METHOD: SELECT CASE (TYPE_METHOD)

   !
   ! ----------------------- Global data-parallel Krylov method --------------------------------------------
   !
   CASE (NSCARC_METHOD_KRYLOV)

      SELECT_KRYLOV_PRECON: SELECT CASE (TYPE_PRECON)

#ifdef WITH_MKL
         !
         ! Preconditioning by defect correction based on LU-decomposition
         ! if two-level method, also use coarse grid level, otherwise only use single (finest) grid level
         !
         CASE (NSCARC_RELAX_MKL)

            IF (HAS_TWO_LEVELS) THEN
               CALL SCARC_GET_NUMBER_OF_LEVELS(NSCARC_LEVEL_MULTI)
               IF (TYPE_COARSE == NSCARC_COARSE_DIRECT) TYPE_MKL(NLEVEL_MAX) = NSCARC_MKL_GLOBAL
            ELSE
               CALL SCARC_GET_NUMBER_OF_LEVELS(NSCARC_LEVEL_SINGLE)
            ENDIF

            ! Either using a global CLUSTER_SPARSE_SOLVER or local PARDISO solvers from MKL
            SELECT CASE (TYPE_SCOPE(1))
               CASE(NSCARC_SCOPE_GLOBAL)
                  TYPE_MKL(NLEVEL_MIN) = NSCARC_MKL_GLOBAL
               CASE(NSCARC_SCOPE_LOCAL)
                  TYPE_MKL(NLEVEL_MIN) = NSCARC_MKL_LOCAL
            END SELECT

#endif

         !
         ! Preconditioning by defect correction based on geometric multigrid method,
         ! use specified hierarchy of grid levels
         !
         CASE (NSCARC_RELAX_GMG)
            CALL SCARC_GET_NUMBER_OF_LEVELS(NSCARC_LEVEL_MULTI)
#ifdef WITH_MKL
            IF (TYPE_COARSE == NSCARC_COARSE_DIRECT) TYPE_MKL(NLEVEL_MAX) = NSCARC_MKL_GLOBAL
#endif

         !
         ! Preconditioning by defect correction based on local basic iterations (JACOBI/SSOR),
         ! if two-level method, also use coarse grid, otherwise only use single (finest) grid level
         !
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

   !
   ! ----------------------- Global data-parallel Multigrid method -------------------------------------
   !
   CASE (NSCARC_METHOD_MULTIGRID)

         !
         ! If not specified by user, determine number of possible grid levels
         !
         CALL SCARC_GET_NUMBER_OF_LEVELS(NSCARC_LEVEL_MULTI)

#ifdef WITH_MKL
         ! in case of smoothing by different MKL solvers, mark levels for the use of MKL,
         ! either by locally acting PARDISO solvers or Globally acting CLUSTER_SPARSE_SOLVER
         IF (TYPE_SMOOTH == NSCARC_RELAX_MKL) THEN

            MKL_SCOPE_SELECT: SELECT CASE (TYPE_SCOPE(2))
            CASE (NSCARC_SCOPE_LOCAL)
               DO NL = NLEVEL_MIN, NLEVEL_MAX-1
                  TYPE_MKL(NL) = NSCARC_MKL_LOCAL
               ENDDO
            CASE (NSCARC_SCOPE_GLOBAL)
               DO NL = NLEVEL_MIN, NLEVEL_MAX-1
                  TYPE_MKL(NL) = NSCARC_MKL_GLOBAL
               ENDDO
            END SELECT MKL_SCOPE_SELECT

         ENDIF

         IF (TYPE_MKL(0) == NSCARC_MKL_COARSE) TYPE_MKL(NLEVEL_MAX) = NSCARC_MKL_GLOBAL
#endif


   !
   ! ----------------------- Global LU-decomposition -----------------------------------------
   !
   CASE (NSCARC_METHOD_LU)

      ! Only use single (finest) grid level
      CALL SCARC_GET_NUMBER_OF_LEVELS(NSCARC_LEVEL_SINGLE)

      ! and mark this level for the use of MKL methods
      SELECT_MKL_SCOPE: SELECT CASE (TYPE_MKL(0))
         CASE (NSCARC_MKL_LOCAL)
            TYPE_MKL(NLEVEL_MIN) = NSCARC_MKL_LOCAL
         CASE (NSCARC_MKL_GLOBAL)
            TYPE_MKL(NLEVEL_MIN) = NSCARC_MKL_GLOBAL
      END SELECT SELECT_MKL_SCOPE

#ifdef WITH_SCARC_MGM
   !
   ! ----------------------- Global McKenney-Greengard-Mayo method - only finest level --------------------
   !
   CASE (NSCARC_METHOD_MGM)

      CALL SCARC_GET_NUMBER_OF_LEVELS(NSCARC_LEVEL_SINGLE)

#endif

END SELECT SELECT_SCARC_METHOD


#ifdef WITH_MKL
!
! Define MKL related logical short names based on number of levels
!
DO NL = NLEVEL_MIN, NLEVEL_MAX
   IS_MKL_LEVEL(NL) = (TYPE_MKL(0)  == NSCARC_MKL_GLOBAL .AND. NL == NLEVEL_MIN) .OR. &
                      (TYPE_MKL(0)  == NSCARC_MKL_COARSE .AND. NL == NLEVEL_MAX) .OR. &
                      (TYPE_MKL(NL) == NSCARC_MKL_GLOBAL)
ENDDO
#endif

END SUBROUTINE SCARC_SETUP_LEVELS


! ------------------------------------------------------------------------------------------------
! Setup single level in case of default Krylov method
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
      IF (IS_GMG.OR.IS_AMG.OR.IS_CG_GMG.OR.IS_CG_AMG.OR.IS_CG_SAMG) THEN

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
   
   ! determine maximum number of levels for AMG method; if nothing has been specified, use default
   CASE(NSCARC_LEVEL_AMG)
      NLEVEL_MIN = 1
      IF (SCARC_MULTIGRID_LEVEL /= -1) THEN
         NLEVEL_MAX  = SCARC_MULTIGRID_LEVEL
      ELSE
         NLEVEL_MAX  = NSCARC_LEVEL_MAX
      ENDIF
      NLEVEL = NLEVEL_MAX
   
END SELECT SELECT_LEVEL_TYPE

END SUBROUTINE SCARC_GET_NUMBER_OF_LEVELS


! ------------------------------------------------------------------------------------------------
! Determine maximum number of possible levels on direction IOR0 of mesh NM
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
! Allocate basic ScaRC-structures for all needed levels
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_TYPES
USE SCARC_POINTERS, ONLY: S
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
! This routine assumes, that L already points to the correct level NL of mesh NL and
! additionally sets the requested discretization type
! ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_ASSIGN_GRID_TYPE(NTYPE)
INTEGER, INTENT(IN) :: NTYPE

SELECT CASE (NTYPE)
   CASE (NSCARC_GRID_STRUCTURED)
      TYPE_GRID = NSCARC_GRID_STRUCTURED
      IS_STRUCTURED   = .TRUE.
      IS_UNSTRUCTURED = .FALSE.
   CASE (NSCARC_GRID_UNSTRUCTURED)
      TYPE_GRID = NSCARC_GRID_UNSTRUCTURED
      IS_STRUCTURED   = .FALSE.
      IS_UNSTRUCTURED = .TRUE.
END SELECT

END SUBROUTINE SCARC_ASSIGN_GRID_TYPE


! -----------------------------------------------------------------------------
! Setup discretization information
! -----------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_GRIDS
USE SCARC_POINTERS, ONLY: M, S, L, G, XCOR, YCOR, ZCOR, XMID, YMID, ZMID
INTEGER :: NL, NM, NC, IX, IY, IZ, IO
INTEGER :: IBAR, JBAR, KBAR

!
! -------------------- On all grid levels --------------------------------------
! Specify general mesh related geometry information
!

MESHES_LOOP1: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_MESH(NM)

   ! store bounds of mesh in SCARC-structure
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

      ! get coordination information
      L%DX = (S%XF-S%XS)/REAL(L%NX,EB)
      L%DY = (S%YF-S%YS)/REAL(L%NY,EB)
      L%DZ = (S%ZF-S%ZS)/REAL(L%NZ,EB)

      L%DXI  = 1.0_EB/L%DX;  L%DYI  = 1.0_EB/L%DY;  L%DZI =  1.0_EB/L%DZ
      L%DXI2 = L%DXI**2;     L%DYI2 = L%DYI**2;     L%DZI2 = L%DZI**2

      ! needed in case of GMG with multiple grid levels
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

      ELSE

         ! Allocate and compute coordinate information for coarser levels
         CALL SCARC_ALLOCATE_REAL1(L%XCOR, 0, L%NX, NSCARC_INIT_NONE, 'L%XCOR')
         CALL SCARC_ALLOCATE_REAL1(L%YCOR, 0, L%NY, NSCARC_INIT_NONE, 'L%YCOR')
         CALL SCARC_ALLOCATE_REAL1(L%ZCOR, 0, L%NZ, NSCARC_INIT_NONE, 'L%ZCOR')

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
         CALL SCARC_ALLOCATE_REAL1(L%XMID, 0, L%NX+1, NSCARC_INIT_NONE, 'L%XMID')
         CALL SCARC_ALLOCATE_REAL1(L%YMID, 0, L%NY+1, NSCARC_INIT_NONE, 'L%YMID')
         CALL SCARC_ALLOCATE_REAL1(L%ZMID, 0, L%NZ+1, NSCARC_INIT_NONE, 'L%ZMID')

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
      CALL SCARC_ALLOCATE_REAL1(L%DXL, 0, L%NX, NSCARC_INIT_ZERO, 'L%DXL')
      CALL SCARC_ALLOCATE_REAL1(L%DYL, 0, L%NY, NSCARC_INIT_ZERO, 'L%DYL')
      CALL SCARC_ALLOCATE_REAL1(L%DZL, 0, L%NZ, NSCARC_INIT_ZERO, 'L%DZL')

      ! set step sizes between cell midpoints, use interior step sizes for ghost cells as initial values
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

!
! ---------------------- On finest grid level -------------------------------------------------
! Allocate several arrays for the administration of discretization related data
!
MESHES_LOOP2: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   !CALL SCARC_POINT_TO_LEVEL(NM, NLEVEL_MIN)
   CALL SCARC_POINT_TO_GRID(NM, NLEVEL_MIN)

   ! Set pointers to already existing CELL_INDEX WALL_INDEX arrays from main program (on finest level)
   CALL SCARC_SETUP_CELL_INDEX(NM, NLEVEL_MIN)
   CALL SCARC_SETUP_WALL_INDEX(NM, NLEVEL_MIN)

   ! Allocate and initialize IS_SOLID array which indicates the state of a cell (gasphase/solid)
   CALL SCARC_ALLOCATE_LOG3(L%IS_SOLID, 0, L%NX+1, 0, L%NY+1, 0, L%NZ+1, NSCARC_INIT_TRUE, 'L%IS_SOLID')
   L%IS_SOLID (1:L%NX, 1:L%NY, 1:L%NZ) = .FALSE.

   ! Identify and mark solid obstruction cells in IS_SOLID-part of the discretization
   DO IZ = 1, L%NZ
      DO IY = 1, L%NY
         DO IX = 1, L%NX
            IF (M%SOLID(M%CELL_INDEX(IX, IY, IZ))) L%IS_SOLID(IX, IY, IZ) = .TRUE.
         ENDDO
      ENDDO
   ENDDO

   !
   ! If both discretization types (structured/unstructured) must be administrated (MGM method only):
   ! Allocate all arrays which are related to a specific discretization type
   !
   IF (HAS_MULTI_GRIDS) THEN

      !
      ! ---------------- First process structured discretization
      !
      CALL SCARC_ASSIGN_GRID_TYPE(NSCARC_GRID_STRUCTURED)
      CALL SCARC_POINT_TO_GRID (NM, NLEVEL_MIN)

      ! Number of local cells per mesh
      CALL SCARC_ALLOCATE_INT1(G%NC_LOCAL , 1, NMESHES, NSCARC_INIT_ZERO, 'G%NC_LOCAL')
      CALL SCARC_ALLOCATE_INT1(G%NC_OFFSET, 1, NMESHES, NSCARC_INIT_ZERO, 'G%NC_OFFSET')

      ! Allocate wall information array
      ALLOCATE(G%WALL(L%N_WALL_CELLS), STAT=IERROR)
      CALL ChkMemErr('SCARC_SETUP_GRIDS','WALL',IERROR)

      ! Allocate and preset cell numbers array
      CALL SCARC_ALLOCATE_INT3(G%CELL_NUMBER, 0, L%NX+1, 0, L%NY+1, 0, L%NZ+1, &
                               NSCARC_INIT_UNDEF, 'CELL_NUMBER')

      ! Define local cell numbers for Poisson equation
      DO IZ=1,L%NZ
         DO IY=1,L%NY
            DO IX=1,L%NX
               G%NC_LOCAL(NM) = G%NC_LOCAL(NM) + 1
               G%CELL_NUMBER(IX,IY,IZ) = G%NC_LOCAL(NM)
            ENDDO
         ENDDO
      ENDDO
      G%NC   = G%NC_LOCAL(NM)
      G%NCE  = G%NC_LOCAL(NM)
      G%NCE2 = G%NC_LOCAL(NM)


      !
      ! ---------------- Then process unstructured discretization
      !
      CALL SCARC_ASSIGN_GRID_TYPE(NSCARC_GRID_UNSTRUCTURED)
      CALL SCARC_POINT_TO_GRID (NM, NLEVEL_MIN)

      ! Also allocate and preset cell numbers and state arrays for unstructured discretization
      CALL SCARC_ALLOCATE_INT3(G%CELL_NUMBER, 0, L%NX+1, 0, L%NY+1, 0, L%NZ+1, &
                               NSCARC_INIT_UNDEF, 'G%CELL_NUMBER')

      ! Number of local cells per mesh
      CALL SCARC_ALLOCATE_INT1(G%NC_LOCAL , 1, NMESHES, NSCARC_INIT_ZERO, 'G%NC_LOCAL')
      CALL SCARC_ALLOCATE_INT1(G%NC_OFFSET, 1, NMESHES, NSCARC_INIT_ZERO, 'G%NC_OFFSET')

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
               ENDIF
            ENDDO
         ENDDO
      ENDDO
      G%NC   = G%NC_LOCAL(NM)
      G%NCE  = G%NC_LOCAL(NM)
      G%NCE2 = G%NC_LOCAL(NM)

   !
   ! If only one specified type of discretization must be admistrated:
   ! Allocate and preset cell numbers and state arrays for requested type of discretization
   !
   ELSE

      ! ---------------- Only process specified type of discretization
      CALL SCARC_ASSIGN_GRID_TYPE(TYPE_GRID)
      CALL SCARC_POINT_TO_GRID (NM, NLEVEL_MIN)

      ! Also allocate and preset cell numbers and state arrays for unstructured discretization
      CALL SCARC_ALLOCATE_INT3(G%CELL_NUMBER, 0, L%NX+1, 0, L%NY+1, 0, L%NZ+1, &
                               NSCARC_INIT_UNDEF, 'CELL_NUMBER')

      ! Allocate index array which specifies I, J, K components for all degrees of freedom
      NC = L%NX * L%NY * L%NZ
      CALL SCARC_ALLOCATE_INT1(G%ICX, 1, NC, NSCARC_INIT_UNDEF, 'G%ICX')
      CALL SCARC_ALLOCATE_INT1(G%ICY, 1, NC, NSCARC_INIT_UNDEF, 'G%ICY')
      CALL SCARC_ALLOCATE_INT1(G%ICZ, 1, NC, NSCARC_INIT_UNDEF, 'G%ICZ')

      ! Number of local cells per mesh
      CALL SCARC_ALLOCATE_INT1(G%NC_LOCAL , 1, NMESHES, NSCARC_INIT_ZERO, 'G%NC_LOCAL')
      CALL SCARC_ALLOCATE_INT1(G%NC_OFFSET, 1, NMESHES, NSCARC_INIT_ZERO, 'G%NC_OFFSET')

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

   ENDIF

ENDDO MESHES_LOOP2

END SUBROUTINE SCARC_SETUP_GRIDS


! -----------------------------------------------------------------------------
! Setup discretization information on coarser levels
! -----------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_GRID_LEVEL(NL)
USE SCARC_POINTERS, ONLY: LF, LC, GC
INTEGER, INTENT(IN) :: NL
INTEGER :: NM, NC, IXF, IYF, IZF, IX, IY, IZ, NSTEP

MESHES_LOOP1: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   LF => SCARC(NM)%LEVEL(NLEVEL_MIN)
   LC => SCARC(NM)%LEVEL(NL)

   CALL SCARC_ALLOCATE_LOG3(LC%IS_SOLID, 0, LC%NX+1, 0, LC%NY+1, 0, LC%NZ+1, NSCARC_INIT_FALSE, 'LC%IS_SOLID')
   LC%IS_SOLID (1:LC%NX, 1:LC%NY, 1:LC%NZ)  = .FALSE.

   NSTEP = 2**(NL - NLEVEL_MIN)

   SELECT CASE(TYPE_GRID)

      !
      ! Get cell numberings for coarser grid in case of structured discretization
      !
      CASE (NSCARC_GRID_STRUCTURED)

         GC => LC%STRUCTURED

         CALL SCARC_ALLOCATE_INT1(GC%NC_LOCAL , 1, NMESHES, NSCARC_INIT_ZERO, 'GC%NC_LOCAL')
         CALL SCARC_ALLOCATE_INT1(GC%NC_OFFSET, 1, NMESHES, NSCARC_INIT_ZERO, 'GC%NC_OFFSET')

         CALL SCARC_ALLOCATE_INT3(GC%CELL_NUMBER, 0, LC%NX+1, 0, LC%NY+1, 0, LC%NZ+1, &
                                  NSCARC_INIT_UNDEF, 'GC%CELL_NUMBER')

         NC = LC%NX * LC%NY * LC%NZ

         CALL SCARC_ALLOCATE_INT1(GC%ICX , 1, NC, NSCARC_INIT_UNDEF, 'GC%ICX')
         CALL SCARC_ALLOCATE_INT1(GC%ICY , 1, NC, NSCARC_INIT_UNDEF, 'GC%ICY')
         CALL SCARC_ALLOCATE_INT1(GC%ICZ , 1, NC, NSCARC_INIT_UNDEF, 'GC%ICZ')

         GC%CELL_NUMBER = NSCARC_UNDEF_INT

         DO IZ = 1, LC%NZ
            IZF = (IZ-1)*NSTEP + 1
            DO IY = 1, LC%NY
               IYF = (IY-1)*NSTEP + 1
               DO IX = 1, LC%NX

                  IXF = (IX-1)*NSTEP + 1
                  LC%IS_SOLID(IX,IY,IZ) = LF%IS_SOLID(IXF, IYF, IZF)

                  GC%NC_LOCAL(NM)  = GC%NC_LOCAL(NM) + 1
                  GC%CELL_NUMBER(IX,IY,IZ) = GC%NC_LOCAL(NM)

                  GC%ICX(GC%NC_LOCAL(NM)) = IX
                  GC%ICY(GC%NC_LOCAL(NM)) = IY
                  GC%ICZ(GC%NC_LOCAL(NM)) = IZ

               ENDDO
            ENDDO
         ENDDO

         GC%NC = GC%NC_LOCAL(NM)

      !
      ! Get cell numberings for coarser grid in case of unstructured discretization
      !
      CASE (NSCARC_GRID_UNSTRUCTURED)

         GC => LC%UNSTRUCTURED

         CALL SCARC_ALLOCATE_INT1(GC%NC_LOCAL , 1, NMESHES, NSCARC_INIT_ZERO, 'GC%NC_LOCAL')
         CALL SCARC_ALLOCATE_INT1(GC%NC_OFFSET, 1, NMESHES, NSCARC_INIT_ZERO, 'GC%NC_OFFSET')

         CALL SCARC_ALLOCATE_INT3(GC%CELL_NUMBER, 0, LC%NX+1, 0, LC%NY+1, 0, LC%NZ+1, &
                                  NSCARC_INIT_UNDEF, 'GC%CELL_NUMBER')

         NC = LC%NX * LC%NY * LC%NZ

         CALL SCARC_ALLOCATE_INT1(GC%ICX , 1, NC, NSCARC_INIT_UNDEF, 'GC%ICX')
         CALL SCARC_ALLOCATE_INT1(GC%ICY , 1, NC, NSCARC_INIT_UNDEF, 'GC%ICY')
         CALL SCARC_ALLOCATE_INT1(GC%ICZ , 1, NC, NSCARC_INIT_UNDEF, 'GC%ICZ')

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


! ------------------------------------------------------------------------------------------------
! Setup communication structure for data exchange along mesh interfaces
! and analyze neighborship structure of mesh subdivison 
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_NEIGHBORS
USE SCARC_POINTERS, ONLY: OS, OLF, OLC
INTEGER :: NM, NOM, NL

! Initialize communication counter for ScaRC, use same TAG for all communications
TAG   = 99
N_REQ =  0
N_EXCHANGES = 0

!
! Initialize level structures on neighboring meshes
!
LEVEL_MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   LEVEL_NEIGHBOR_LOOP: DO NOM = 1, NMESHES

      !
      ! On finest level point to exchange structures from surrounding FDS 
      !
      IF (.NOT. ARE_NEIGHBORS(NM, NOM)) CYCLE LEVEL_NEIGHBOR_LOOP

      N_EXCHANGES = N_EXCHANGES+1                                         ! count number of exchanges

      OS => SCARC(NM)%OSCARC(NOM)
      ALLOCATE (OS%LEVEL(NLEVEL_MIN:NLEVEL_MAX), STAT=IERROR)             ! allocate neighboring structures
      CALL CHKMEMERR ('SCARC_SETUP_NEIGHBORS', 'OS%LEVEL', IERROR)

      OLF => SCARC(NM)%OSCARC(NOM)%LEVEL(NLEVEL_MIN)                      ! point to neighbor on finest grid level

      OLF%NX = MESHES(NOM)%IBAR                                           ! number of cells in x-direction on neighbor
      OLF%NY = MESHES(NOM)%JBAR                                           ! number of cells in y-direction on neighbor
      OLF%NZ = MESHES(NOM)%KBAR                                           ! number of cells in z-direction on neighbor

      OLF%N_WALL_CELLS_EXT = MESHES(NOM)%N_EXTERNAL_WALL_CELLS            ! number of external wall cells on neighbor
      OLF%N_WALL_CELLS_INT = MESHES(NOM)%N_INTERNAL_WALL_CELLS            ! number of external wall cells on neighbor
      OLF%N_WALL_CELLS     = OLF%N_WALL_CELLS_EXT + OLF%N_WALL_CELLS_INT  ! number of walls cell on neighbor

      OLF%N_CELLS = OLF%NX*OLF%NY*OLF%NZ                                  ! number of cells on neighbor (structured)

      !
      ! In case of GMG with a predefined grid hierarchy define corresponding level-structures
      !
      IF (NLEVEL_MAX > NLEVEL_MIN) THEN                                   

         DO NL=NLEVEL_MIN+1,NLEVEL_MAX

            OLC => SCARC(NM)%OSCARC(NOM)%LEVEL(NL)                        ! OLF points to finer, OLC to coarser level

            OLC%NX = OLF%NX/2                                             ! use double grid width
            IF (TWO_D) THEN
               OLC%NY = 1
            ELSE
               OLC%NY = OLF%NY/2
            ENDIF
            OLC%NZ = OLF%NZ/2

            OLC%N_CELLS          = OLC%NX * OLC%NY * OLC%NZ               ! set new number of cells
            OLC%N_WALL_CELLS     = OLC%N_WALL_CELLS_EXT                   ! set new number of wall cells
            OLC%N_WALL_CELLS_EXT = 2 * (OLC%NX*OLC%NZ + OLC%NX*OLC%NY + OLC%NY*OLC%NZ)    ! TODO: CHECK!

         ENDDO
      ENDIF

   ENDDO LEVEL_NEIGHBOR_LOOP
ENDDO LEVEL_MESHES_LOOP

END SUBROUTINE SCARC_SETUP_NEIGHBORS


! ----------------------------------------------------------------------------------------------------
! Setup FACE related structures on finest grid level:
!   - get dimensions for each of the 6 faces of a mesh
!   - get grid width vector along face
!   - get information for adjacent neighbors
!   - allocate pointer arrays for data exchanges with neighbors
! ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_FACES
USE SCARC_POINTERS, ONLY: M, S, L, LC, F
INTEGER :: NL, NM, NOM
INTEGER :: IFACE, IOR0, JOR0, INBR, IWG, ICW
LOGICAL :: IS_KNOWN(-3:3)
INTEGER :: FACE_NEIGHBORS(-3:3, NSCARC_MAX_FACE_NEIGHBORS)
INTEGER :: MESH_NEIGHBORS(6*NSCARC_MAX_FACE_NEIGHBORS)
INTEGER :: N_FACE_NEIGHBORS(-3:3)
INTEGER :: N_MESH_NEIGHBORS

MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   !CALL SCARC_POINT_TO_LEVEL(NM, NLEVEL_MIN)             ! consider only finest grid level
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

   CALL SCARC_SETUP_FACE_BASICS(NM, NLEVEL_MIN)

   !
   ! Store first wall cell number for each face
   !
   ICW = 1
   FACE_ORDER_LOOP: DO IFACE = 1, 6
      IOR0 = FACE_ORIENTATION(IFACE)
      F => L%FACE(IOR0)
      F%NCW0 = ICW
      ICW = ICW + F%NCW
   ENDDO FACE_ORDER_LOOP

   !
   ! loop over external wall cells:
   ! store basic data and determine number of adajacent neighbors to each face
   !
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

   !
   ! Allocate array which stores numbers of all neighboring meshes
   !
   IF (N_MESH_NEIGHBORS /= 0) &
      CALL SCARC_ALLOCATE_INT1(S%NEIGHBORS, 1, N_MESH_NEIGHBORS, NSCARC_INIT_UNDEF, 'S%NEIGHBORS')
   S%N_NEIGHBORS = N_MESH_NEIGHBORS

   NEIGHBORS_OF_FACE_LOOP: DO IOR0 = -3, 3

      IF (IOR0 == 0) CYCLE NEIGHBORS_OF_FACE_LOOP

      ! if there are neighbors at face IOR0 store information about them
      F => L%FACE(IOR0)
      IF (N_FACE_NEIGHBORS(IOR0) /= 0) THEN

         ! allocate array for storing the numbers of the single neighbors
         F%N_NEIGHBORS = N_FACE_NEIGHBORS(IOR0)
         CALL SCARC_ALLOCATE_INT1(F%NEIGHBORS, 1, N_FACE_NEIGHBORS(IOR0), NSCARC_INIT_NONE, 'F%NEIGHBORS')

         ! store every neighbor and allocate corresponding administration arrays on finest level
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
! Setup subdivision information 
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_SUBDIVISION
USE SCARC_POINTERS, ONLY: SUB
INTEGER, ALLOCATABLE, DIMENSION(:) :: BUFFER_INT
INTEGER, ALLOCATABLE, DIMENSION(:) :: COUNTS_NBR   
INTEGER, ALLOCATABLE, DIMENSION(:) :: DISPLS_NBR    
INTEGER :: N, NM, INBR, IP, MAX_NBR

! Determine number of neighbors for each mesh and make them available in a global array
SUB => SUBDIVISION

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) '================= DETERMINE SUB%NEIGHBORS'
#endif

CALL SCARC_ALLOCATE_INT1 (SUB%N_NEIGHBORS, 1, NMESHES, NSCARC_INIT_ZERO, 'SUB%N_NEIGHBORS')
DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   SUB%N_NEIGHBORS(NM) = SCARC(NM)%N_NEIGHBORS
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'NM, SUB%N_NEIGHBORS=', NM, SUB%N_NEIGHBORS(NM)
#endif
ENDDO

IF (N_MPI_PROCESSES > 1) &
   CALL MPI_ALLGATHERV(MPI_IN_PLACE,1,MPI_INTEGER,SUB%N_NEIGHBORS,COUNTS,DISPLS, MPI_INTEGER,MPI_COMM_WORLD,IERROR)
SUB%N_NEIGHBORS_TOTAL = SUM(SUB%N_NEIGHBORS)

CALL SCARC_ALLOCATE_INT1 (COUNTS_NBR, 0, N_MPI_PROCESSES-1, NSCARC_INIT_ZERO, 'COUNTS_NBR')
CALL SCARC_ALLOCATE_INT1 (DISPLS_NBR, 0, N_MPI_PROCESSES-1, NSCARC_INIT_ZERO, 'DISPLS_NBR')

DO N = 0, N_MPI_PROCESSES - 1
   DO NM = 1, NMESHES
      IF (PROCESS(NM) == N) COUNTS_NBR(N) = COUNTS_NBR(N) + SUB%N_NEIGHBORS(NM)
   ENDDO
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'N, COUNTS_NBR=', N, COUNTS_NBR(N)
#endif
ENDDO
DO N = 1, N_MPI_PROCESSES -1
   DISPLS_NBR(N) = COUNTS_NBR(N-1) + DISPLS_NBR(N-1)
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'N, DISPLS_NBR=', N, DISPLS_NBR(N)
#endif
ENDDO

CALL SCARC_ALLOCATE_INT1 (BUFFER_INT, 1, SUB%N_NEIGHBORS_TOTAL, NSCARC_INIT_ZERO, 'BUFFER_INT')
IP = 1
DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   DO INBR = 1, SCARC(NM)%N_NEIGHBORS
      BUFFER_INT(DISPLS_NBR(PROCESS(NM)) + IP) = SCARC(NM)%NEIGHBORS(INBR)
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'NM, INBR, BUFFER_INT=', NM, INBR, SCARC(NM)%NEIGHBORS(INBR)
#endif
      IP = IP + 1
   ENDDO
ENDDO

IF (N_MPI_PROCESSES > 1) &
   CALL MPI_ALLGATHERV(MPI_IN_PLACE,1,MPI_INTEGER,BUFFER_INT,COUNTS_NBR,DISPLS_NBR, MPI_INTEGER,MPI_COMM_WORLD,IERROR)
MAX_NBR = MAXVAL(SUB%N_NEIGHBORS)

CALL SCARC_ALLOCATE_INT2 (SUB%NEIGHBORS, 1, MAX_NBR, 1, NMESHES,  NSCARC_INIT_ZERO, 'SUB%NEIGHBORS')
DO NM = 1, NMESHES
   DO INBR = 1, SUB%N_NEIGHBORS(NM)
      SUB%NEIGHBORS(INBR, NM) = BUFFER_INT(DISPLS_NBR(PROCESS(NM)) + INBR)
   ENDDO
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'MESH ', NM,': SUB%NEIGHBORS:', SUB%NEIGHBORS(1:SUB%N_NEIGHBORS(NM),NM)
#endif
ENDDO

DEALLOCATE(BUFFER_INT)
DEALLOCATE(COUNTS_NBR)
DEALLOCATE(DISPLS_NBR)

END SUBROUTINE SCARC_SETUP_SUBDIVISION

! ----------------------------------------------------------------------------------------------------
! Determine basic data for single faces (orientation, dimensions, numbers)
! ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_FACE_BASICS(NM, NL)
USE SCARC_POINTERS, ONLY: L, F
INTEGER, INTENT(IN) :: NM, NL
INTEGER:: IOR0

L => SCARC(NM)%LEVEL(NL)

FACES_OF_MESH_LOOP: DO IOR0 = -3, 3

   IF (IOR0 == 0) CYCLE FACES_OF_MESH_LOOP

   
   F => L%FACE(IOR0)
   
   SELECT CASE (ABS(IOR0))

      ! 
      ! ---------- Faces in x-direction
      ! 
      CASE (1)

         F%NOP =  L%NX                           ! number of cells between opposite mesh faces

         F%NX  =  1                              ! number of cells in x-direction
         F%NY  =  L%NY                           ! number of cells in y-direction
         F%NZ  =  L%NZ                           ! number of cells in z-direction

         F%NCW =  L%NY*L%NZ                      ! number of wall cells at that face
         F%DH  => L%DXL                          ! step size vector between opposite mesh faces

         F%SCAL_BOUNDARY  = L%DXI2               ! contribution due to boundary condition 
         F%SCAL_DIRICHLET = -2.0_EB * L%DXI2

         F%NOFFY =  0
         F%NOFFZ =  0
         IF (IOR0 > 0) THEN
            F%SCAL_NEUMANN = L%DXI
            F%SCAL_FACE    = 2.0_EB/(F%DH(0)*(F%DH(0)+F%DH(1)))
            F%SCAL_INSIDE  = 2.0_EB/(F%DH(1)*(F%DH(1)+F%DH(2)))
            F%NOFFX =  1                           ! offset to next internal cell in that direction
         ELSE
            F%SCAL_NEUMANN = -L%DXI
            F%SCAL_FACE    = 2.0_EB/(F%DH(F%NOP)  *(F%DH(F%NOP-1)+F%DH(F%NOP)))
            F%SCAL_INSIDE  = 2.0_EB/(F%DH(F%NOP-1)*(F%DH(F%NOP-2)+F%DH(F%NOP-1)))
            F%NOFFX = -1
         ENDIF     

      ! 
      ! ---------- Faces in y-direction
      ! 
      CASE (2)

         F%NOP =  L%NY                   ! dito

         F%NX  =  L%NX
         F%NY  =  1
         F%NZ  =  L%NZ

         F%NCW =  L%NX*L%NZ
         F%DH  => L%DYL

         F%SCAL_BOUNDARY  = L%DYI2
         F%SCAL_DIRICHLET = -2.0_EB * L%DYI2

         F%NOFFX =  0                           ! offset to next internal cell in that direction
         F%NOFFY =  0
         F%NOFFZ =  0
         IF (IOR0>0) THEN
            F%SCAL_NEUMANN = L%DYI
            IF (.NOT.TWO_D) THEN
               F%SCAL_FACE   = 2.0_EB/(F%DH(0)*(F%DH(0)+F%DH(1)))
               F%SCAL_INSIDE = 2.0_EB/(F%DH(1)*(F%DH(1)+F%DH(2)))
               F%NOFFY = 1
            ENDIF
         ELSE
            F%SCAL_NEUMANN = -L%DYI
            IF (.NOT.TWO_D) THEN
               F%SCAL_FACE   =  2.0_EB/(F%DH(F%NOP)  *(F%DH(F%NOP-1)+F%DH(F%NOP)))
               F%SCAL_INSIDE =  2.0_EB/(F%DH(F%NOP-1)*(F%DH(F%NOP-2)+F%DH(F%NOP-1)))
               F%NOFFY = -1
            ENDIF
         ENDIF

      ! 
      ! ---------- Faces in z-direction
      ! 
      CASE (3)

         F%NOP =  L%NZ                   ! dito

         F%NX  =  L%NX
         F%NY  =  L%NY
         F%NZ  =  1

         F%NCW =  L%NX*L%NY
         F%DH  => L%DZL

         F%NX  = L%NX
         F%SCAL_BOUNDARY  = L%DZI2
         F%SCAL_DIRICHLET = -2.0_EB * L%DZI2

         F%NOFFX =  0
         F%NOFFY =  0
         IF (IOR0>0) THEN
            F%SCAL_NEUMANN = L%DZI
            F%SCAL_FACE    = 2.0_EB/(F%DH(0)*(F%DH(0)+F%DH(1)))
            F%SCAL_INSIDE  = 2.0_EB/(F%DH(1)*(F%DH(1)+F%DH(2)))
            F%NOFFZ =  1
         ELSE
            F%SCAL_NEUMANN = -L%DZI
            F%SCAL_FACE    =  2.0_EB/(F%DH(F%NOP)  *(F%DH(F%NOP-1)+F%DH(F%NOP)))
            F%SCAL_INSIDE  =  2.0_EB/(F%DH(F%NOP-1)*(F%DH(F%NOP-2)+F%DH(F%NOP-1)))
            F%NOFFZ = -1
         ENDIF
   END SELECT

ENDDO FACES_OF_MESH_LOOP

END SUBROUTINE SCARC_SETUP_FACE_BASICS


! ----------------------------------------------------------------------------------------------------
! Setup WALL related structures and boundary conditions
! ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_WALLS
USE SCARC_POINTERS, ONLY: M, L, LF, LC, FF, FC, OL, OLF, OLC, G, GC, GF, WC, WF, &
                          OGC, OGF, GWC, MWC, EWC
INTEGER :: NL, NM, NOM
INTEGER :: IREFINE, IFACE, IOR0, JOR0, INBR, IWG, IWC, ICW
LOGICAL :: IS_KNOWN(-3:3), IS_DIRIC, IS_OPEN

!
! -------- Get dimensionings for wall cells
!
MESHES_LOOP1: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID (NM, NLEVEL_MIN)

   !
   ! First loop over external wall cells:
   ! Determine number of adajacent neighbors to each face with corresponding number of IW's
   ! Store neighbors, orientation and number of couplings for a single wall cell
   !
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

         G%NCE  = G%NCE  + 1                                          ! increase number of extended grid cells
         IF (IS_AMG .OR. IS_CG_SAMG) G%NCE2 = G%NCE2 + 2              ! increase number of extended grid cells type2
         CALL SCARC_POINT_TO_OTHER_GRID(NM, NOM, NLEVEL_MIN)

         IF (ANY(IS_KNOWN)) OG%NCG = OG%NCG + 1                       ! increase counter for local ghost cells
         IF (OL%GHOST_FIRSTW(IOR0) == 0) OL%GHOST_FIRSTW(IOR0) = OG%NCG               ! save first ghost cell for -IOR0
         IF (OL%GHOST_FIRSTE(IOR0) == 0) OL%GHOST_FIRSTE(IOR0) = OG%NCG               ! save first ghost cell for -IOR0
         OL%GHOST_LASTW(IOR0) = OG%NCG                                     
         OL%GHOST_LASTE(IOR0) = OG%NCG                                     

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'NM, NOM, IWG, IOR0, NCG, NCE2:', NM, NOM, IWG, IOR0, OG%NCG, G%NCE2
#endif

      ENDIF

   ENDDO EXTERNAL_WALL_CELLS_LOOP1
   IF (IS_AMG .OR. IS_CG_SAMG) G%ICE2 = G%NCE                           ! initialize counter for second layer ghost cells

   !
   ! Then process internal wall cells
   !
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

   !
   ! Allocate corresponding pointer arrays for data exchanges with neighbors
   !
   IF (G%NCE > G%NC) THEN
      G%N_FINE = G%NCE
      CALL SCARC_ALLOCATE_INT1 (G%ICE_TO_IWG, G%NC+1, G%NCE, NSCARC_INIT_ZERO, 'G%ICE_TO_IWG')
      CALL SCARC_ALLOCATE_INT1 (G%ICE_TO_ICN, G%NC+1, G%NCE, NSCARC_INIT_ZERO, 'G%ICE_TO_IWG')
   ENDIF

   FACE_NEIGHBORS_LOOP: DO IOR0 = -3, 3

      IF (IOR0 == 0) CYCLE FACE_NEIGHBORS_LOOP
      DO INBR = 1, L%FACE(IOR0)%N_NEIGHBORS

         NOM = L%FACE(IOR0)%NEIGHBORS(INBR)
         CALL SCARC_POINT_TO_OTHER_GRID(NM, NOM, NLEVEL_MIN)

         CALL SCARC_ALLOCATE_INT1 (OG%ICG_TO_IWG, 1, OG%NCG, NSCARC_INIT_ZERO, 'OG%ICG_TO_IWG')
         IF (IS_AMG .OR. IS_CG_SAMG) THEN
#ifdef WITH_SCARC_DEBU
WRITE(MSG%LU_DEBUG,*) 'ALLOCATING ICG_TO_ICE FOR A: ', NOM,': in length ', OG%NCG, 2
#endif
            CALL SCARC_ALLOCATE_INT2 (OG%ICG_TO_ICW, 1, OG%NCG, 1, 2, NSCARC_INIT_ZERO, 'OG%ICG_TO_ICW')
            CALL SCARC_ALLOCATE_INT2 (OG%ICG_TO_ICE, 1, OG%NCG, 1, 2, NSCARC_INIT_ZERO, 'OG%ICG_TO_ICE')
         ELSE
#ifdef WITH_SCARC_DEBU
WRITE(MSG%LU_DEBUG,*) 'ALLOCATING ICG_TO_ICE FOR B: ', NOM,': in length ', OG%NCG, 1
#endif
            CALL SCARC_ALLOCATE_INT2 (OG%ICG_TO_ICW, 1, OG%NCG, 1, 1, NSCARC_INIT_ZERO, 'OG%ICG_TO_ICW')
            CALL SCARC_ALLOCATE_INT2 (OG%ICG_TO_ICE, 1, OG%NCG, 1, 1, NSCARC_INIT_ZERO, 'OG%ICG_TO_ICE')
         ENDIF

      ENDDO

   ENDDO FACE_NEIGHBORS_LOOP

   !
   ! Second loop over external wall cells:
   ! Store detailed coordinate and cell data and get type of boundary condition
   !
   G%ICE = G%NC
   WALL_CELLS_LOOP2: DO IWG = 1, L%N_WALL_CELLS_EXT

      NOM  = G%WALL(IWG)%NOM
      IOR0 = G%WALL(IWG)%IOR

      MWC => M%WALL(IWG)
      EWC => M%EXTERNAL_WALL(IWG)

      !
      ! Preset ScaRC's boundary type indicator BTYPE
      ! INTERNAL  : the global Poisson problem is solved, so no BC's along mesh interfaces are needed
      ! DIRICHLET : - in the structured case face-wise BC-settings are used ccording to original FFT-solver
      !               (this also allows to use FFT as local preconditioner)
      !             - in the unstructured case Dirichlet BC's are only used for open boundary cells
      ! NEUMANN   : is used for the rest
      !
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
      IF (NOM /= 0) THEN
         CALL SCARC_SETUP_WALL_NEIGHBOR(EWC%IIO_MIN, EWC%IIO_MAX, &
                                        EWC%JJO_MIN, EWC%JJO_MAX, &
                                        EWC%KKO_MIN, EWC%KKO_MAX, &
                                        IWG, NM, NOM, NLEVEL_MIN)
      ENDIF

   ENDDO WALL_CELLS_LOOP2

ENDDO MESHES_LOOP1

!
! Set dimensions on finest level for requested type(s) of discretization
! and mapping from local to global cell numbering
!
CALL SCARC_SETUP_DIMENSIONS(NLEVEL_MIN)

!
! -------- For multi-level variants get discretization information and dimensions on coarser levels
!
DO NL = NLEVEL_MIN+1, NLEVEL_MAX
   CALL SCARC_SETUP_GRID_LEVEL(NL)
   CALL SCARC_SETUP_DIMENSIONS(NL)
ENDDO

!
! -------- Check whether there are no Dirichlet BC's available - TODO: Check !!!
!
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


!
! -------- Only for multi-level variants 
! (twolevel-CG or GMG-method as main solver or preconditioner):
! Determine WALL, FACE and OSCARC types for coarser levels
!
MULTI_LEVEL_IF: IF (HAS_MULTI_LEVELS .AND. .NOT.IS_AMG .AND. .NOT.IS_CG_SAMG) THEN

   MESHES_LOOP3: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

      CALL SCARC_POINT_TO_MESH(NM)

      IREFINE=1
      LEVEL_GMG_LEVEL_LOOP: DO NL = NLEVEL_MIN+1, NLEVEL_MAX

         LF => SCARC(NM)%LEVEL(NL-1)                 ! LF points to finer level
         LC => SCARC(NM)%LEVEL(NL)                   ! LC points to coarser level

         SELECT CASE(TYPE_GRID)
            CASE (NSCARC_GRID_STRUCTURED)
               GF => LF%STRUCTURED
               GC => LC%STRUCTURED
            CASE (NSCARC_GRID_UNSTRUCTURED)
               GF => LF%UNSTRUCTURED
               GC => LC%UNSTRUCTURED
         END SELECT

         IREFINE=IREFINE*2
         IF (IS_GMG) CALL SCARC_CHECK_DIVISIBILITY(GF%NCE-GF%NC, 'GF%NCE')

         !
         ! Initialize counts for overlapping and wall cells
         !
         GC%NCE = GC%NC + (GF%NCE-GF%NC)/2
         GC%ICE = GC%NC

         LC%N_WALL_CELLS_EXT = SCARC_COUNT_EXTERNAL_WALL_CELLS(NM, NL)
         LC%N_WALL_CELLS_INT = SCARC_COUNT_INTERNAL_WALL_CELLS(NM, NL)
         LC%N_WALL_CELLS     = LC%N_WALL_CELLS_EXT + LC%N_WALL_CELLS_INT

         SELECT CASE(TYPE_GRID)
            CASE (NSCARC_GRID_STRUCTURED)
               GC%NW = LC%N_WALL_CELLS_EXT
            CASE (NSCARC_GRID_UNSTRUCTURED)
               GC%NW = LC%N_WALL_CELLS_EXT + LC%N_WALL_CELLS_INT
         END SELECT

         ALLOCATE(GC%WALL(LC%N_WALL_CELLS), STAT=IERROR)
         CALL ChkMemErr('SCARC_SETUP_NEIGHBORS','WALL',IERROR)

         !
         ! First allocate administrative mapping arrays for own mesh
         !
         IF (GC%NCE > GC%NC) THEN
            CALL SCARC_ALLOCATE_INT1 (GC%ICE_TO_IWG , GC%NC+1, GC%NCE, NSCARC_INIT_ZERO, 'GC%ICE_TO_IWG')
            CALL SCARC_ALLOCATE_INT1 (GC%ICE_TO_ICN , GC%NC+1, GC%NCE, NSCARC_INIT_ZERO, 'GC%ICE_TO_ICN')
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'ALLOCATING ICE_TO_ICN ON COARSE LEVEL IN LENGTH ', GC%NC+1, GC%NCE, 1, 1
#endif
         ENDIF

         !
         ! set basic and wall information for all faces of coarser level
         !
         CALL SCARC_SETUP_FACE_BASICS(NM, NL)

         WC => GC%WALL
         WF => GF%WALL

         IWC = 1
         IWG = 1
         ICW = 1
         DO IFACE = 1, 6

            IOR0 = FACE_ORIENTATION(IFACE)

            FF => LF%FACE(IOR0)
            FC => LC%FACE(IOR0)

            ! initialize FACE type for coarser mesh
            FC%NCW0 = ICW
            FC%N_NEIGHBORS = FF%N_NEIGHBORS

            IF (FC%N_NEIGHBORS /= 0) &
               CALL SCARC_ALLOCATE_INT1(FC%NEIGHBORS, 1, FC%N_NEIGHBORS, NSCARC_INIT_NONE, 'FC%FACE_NEIGHBORS')
            DO INBR= 1, FC%N_NEIGHBORS
               FC%NEIGHBORS(INBR) = FF%NEIGHBORS(INBR)
            ENDDO
            
            FC%NCW = FC%NX * FC%NY * FC%NZ                                ! get number of wall cells for that face
            ICW = ICW + FC%NCW                                            ! increase global wall cell counter

            ! get related data and pointer structures for every mesh neighbor
            IF (LF%FACE(IOR0)%N_NEIGHBORS /= 0) THEN
               DO INBR = 1, LF%FACE(IOR0)%N_NEIGHBORS

                  NOM = LF%FACE(IOR0)%NEIGHBORS(INBR)

                  OLF => SCARC(NM)%OSCARC(NOM)%LEVEL(NL-1)
                  OLC => SCARC(NM)%OSCARC(NOM)%LEVEL(NL)

                  SELECT CASE(TYPE_GRID)
                     CASE (NSCARC_GRID_STRUCTURED)
                        OGF => OLF%STRUCTURED
                        OGC => OLC%STRUCTURED
                     CASE (NSCARC_GRID_UNSTRUCTURED)
                        OGF => OLF%UNSTRUCTURED
                        OGC => OLC%UNSTRUCTURED
                  END SELECT

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

                  CALL SCARC_SETUP_EXCHANGE_DIMENSIONS(IREFINE, NM, NOM, NL)

                  CALL SCARC_ALLOCATE_INT1(OGC%ICG_TO_IWG, 1, OGC%NCG, NSCARC_INIT_ZERO, 'OGC%ICG_TO_IWG')
                  CALL SCARC_ALLOCATE_INT2(OGC%ICG_TO_ICW, 1, OGC%NCG, 1, 1, NSCARC_INIT_ZERO, 'OGC%ICG_TO_ICW')
                  CALL SCARC_ALLOCATE_INT2(OGC%ICG_TO_ICE, 1, OGC%NCG, 1, 1, NSCARC_INIT_ZERO, 'OGC%ICG_TO_ICE')

               ENDDO
            ENDIF

            ! setup complete face information for coarser mesh
            CALL SCARC_SETUP_WALL_LEVEL(IOR0, IWC, IREFINE, NM, NL)

         ENDDO

         CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G
         CALL SCARC_SETUP_CELL_INDEX (NM, NL)
         CALL SCARC_SETUP_WALL_COORDS(NM, NL)
         CALL SCARC_SETUP_WALL_INDEX (NM, NL)
         CALL SCARC_SETUP_WALL_GHOSTS(NM, NL)

      ENDDO LEVEL_GMG_LEVEL_LOOP
   ENDDO MESHES_LOOP3
ENDIF MULTI_LEVEL_IF


! -------------------------------------------------------------------------------------------
! Correct boundary types for cells adjacent to obstructions on ghost cells
! -------------------------------------------------------------------------------------------
DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   IF (IS_UNSTRUCTURED) THEN
      CALL SCARC_IDENTIFY_INTERNAL_NEUMANNS(NM, NLEVEL_MIN)
      IF (.NOT.IS_AMG .AND. .NOT.IS_CG_AMG .AND. .NOT.IS_CG_SAMG) THEN
         DO NL = NLEVEL_MIN+1, NLEVEL_MAX
            CALL SCARC_IDENTIFY_INTERNAL_NEUMANNS(NM, NL)
         ENDDO
      ENDIF
   ENDIF
ENDDO

! -------------------------------------------------------------------------------------------
! Debug FACE, WALL and DISCRET structures - only if directive SCARC_DEBUG is set
! -------------------------------------------------------------------------------------------
#ifdef WITH_SCARC_DEBUG
CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_STACK, NLEVEL_MIN, 'STACK')
CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_FACE , NLEVEL_MIN, 'FACE')
CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_WALL , NLEVEL_MIN, 'WALL')
CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_GRID , NLEVEL_MIN, 'DISCRET')
#endif

END SUBROUTINE SCARC_SETUP_WALLS


! ------------------------------------------------------------------------------------------------
! Allocate workspace for data exchanges of different data types and sizes
! Perform first basic data exchanges to fix some dimensionings
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_EXCHANGES
USE SCARC_POINTERS, ONLY:  OS, OG
INTEGER :: NL, NM, NOM
INTEGER :: INBR

!
! Allocate request array for data exchanges
! Exchange basic information about wall sizes (needed for the dimensioning of the exchange buffers)
!
IF (N_MPI_PROCESSES>1) THEN
   ALLOCATE (REQ(N_EXCHANGES*40), STAT=IERROR)
   CALL CHKMEMERR ('SCARC_SETUP_EXCHANGES', 'REQ', IERROR)
   REQ = MPI_REQUEST_NULL
ENDIF
CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_BASIC_SIZES, NSCARC_NONE, NLEVEL_MIN)

!
! Allocate send and receive buffers (real and integer) in correct lengths
! These are allocated with sizes according to the requirements of the finest grid level
! In case of a multi-level method, they are also used for the coarser levels (with shorter exchange sizes)
!
DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   DO INBR = 1, SCARC(NM)%N_NEIGHBORS

      CALL SCARC_POINT_TO_GRID(NM, NLEVEL_MIN)                          

      NOM = S%NEIGHBORS(INBR)
      CALL SCARC_POINT_TO_OTHER_GRID(NM, NOM, NLEVEL_MIN)

      ! define maximum lengths for the single overlapping buffers
      OS%NLEN_MAX_BUFFER_LAYER1  = OG%NCG
      OS%NLEN_MAX_BUFFER_LAYER2  = OG%NCG * 2
      OS%NLEN_MAX_BUFFER_LAYER4  = OG%NCG * 4
      OS%NLEN_MAX_BUFFER_STENCIL = OG%NCG * NSCARC_MAX_STENCIL 
      OS%NLEN_MAX_BUFFER_FULL    = OG%NCG * NSCARC_MAX_STENCIL * 4

      ! Allocate buffers in maximum length for finest grid level (same buffers are used in shorter length on coarse levels, too)
      OG%NLEN_BUFFER_LAYER1  = OS%NLEN_MAX_BUFFER_LAYER1
      OG%NLEN_BUFFER_LAYER2  = OS%NLEN_MAX_BUFFER_LAYER2
      OG%NLEN_BUFFER_LAYER4  = OS%NLEN_MAX_BUFFER_LAYER4
      OG%NLEN_BUFFER_STENCIL = OS%NLEN_MAX_BUFFER_STENCIL
      OG%NLEN_BUFFER_FULL    = OS%NLEN_MAX_BUFFER_FULL
     
      CALL SCARC_ALLOCATE_INT1 (OS%SEND_BUFFER_INT , 1, 2*OS%NLEN_MAX_BUFFER_FULL, NSCARC_INIT_UNDEF, 'OS%SEND_BUFFER_INT')
      CALL SCARC_ALLOCATE_INT1 (OS%RECV_BUFFER_INT , 1, 2*OS%NLEN_MAX_BUFFER_FULL, NSCARC_INIT_UNDEF, 'OS%RECV_BUFFER_INT')
      CALL SCARC_ALLOCATE_REAL1(OS%SEND_BUFFER_REAL, 1, 2*OS%NLEN_MAX_BUFFER_FULL, NSCARC_INIT_UNDEF, 'OS%SEND_BUFFER_REAL')
      CALL SCARC_ALLOCATE_REAL1(OS%RECV_BUFFER_REAL, 1, 2*OS%NLEN_MAX_BUFFER_FULL, NSCARC_INIT_UNDEF, 'OS%RECV_BUFFER_REAL')

      ! neighboring wall structures for common wall cells
      ALLOCATE (OG%WALL(OG%NCG), STAT=IERROR)
      CALL CHKMEMERR ('SCARC_SETUP_NEIGHBORS', 'OG%WALL', IERROR)

      ! In case of GMG with a predefined grid hierarchy allocate corresponding level-structures
      IF (NLEVEL_MAX > NLEVEL_MIN) THEN
         DO NL=NLEVEL_MIN+1,NLEVEL_MAX
            CALL SCARC_POINT_TO_OTHER_GRID(NM, NOM, NL)
            ALLOCATE (OG%WALL(OG%NCG), STAT=IERROR)
            CALL CHKMEMERR ('SCARC_SETUP_NEIGHBORS', 'OG%WALL', IERROR)
         ENDDO
      ENDIF

   ENDDO
ENDDO

!
! If there is more than 1 mesh, initialize communication structures on finest level 
! and setup mapping from local to global numbering
!
IF (NMESHES > 1) THEN
   CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_CELL_NUMBERS, NSCARC_NONE, NLEVEL_MIN)
   CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_CELL_SIZES,   NSCARC_NONE, NLEVEL_MIN)

   IF (HAS_MULTI_LEVELS .AND. .NOT.IS_AMG .AND. .NOT.IS_CG_SAMG) THEN
      DO NL = NLEVEL_MIN+1, NLEVEL_MAX
         CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_CELL_NUMBERS, NSCARC_NONE, NL)
         CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_CELL_SIZES,   NSCARC_NONE, NL)
      ENDDO
   ENDIF
ENDIF

END SUBROUTINE SCARC_SETUP_EXCHANGES


! -----------------------------------------------------------------------------------------
! --- Store neighbors of mesh
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
! --- Setup CELL_INDEX array on coarser grid levels in case of MG-method
! -----------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_CELL_INDEX(NM, NL)
USE SCARC_POINTERS, ONLY: M, L, OB
INTEGER, INTENT(IN) :: NM, NL
INTEGER :: I, J, K, NOBST

CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G

! if finest level, the corresponding CELL_INDEX array is already available by surrounding routines
IF (NL == NLEVEL_MIN) THEN
   L%CELL_INDEX_PTR => M%CELL_INDEX

   ! on coarser levels, it must still be computed
ELSE

   CALL SCARC_ALLOCATE_INT3(L%CELL_INDEX, 0, L%NX+1, 0, L%NY+1, 0, L%NZ+1, NSCARC_INIT_ZERO, 'L%CELL_INDEX')
   L%CELL_INDEX_PTR => L%CELL_INDEX
   L%N_CELL_INDEX = 0

   !
   ! Preset it for all grid cells
   !
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

   !
   ! Consider cells in obstructions
   !
   DO NOBST=1,L%N_OBST
      OB => SCARC(NM)%LEVEL(NL)%OBST(NOBST)
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
! In case of an MG-method:
! Setup WALL_INDEX array on coarser grid levels
! -----------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_WALL_INDEX(NM, NL)
USE SCARC_POINTERS, ONLY: M, L, G
INTEGER, INTENT(IN) :: NM, NL
INTEGER :: I, J, K, ICG, IW, IOR0


CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G

! if on finest level, the array WALL_INDEX is already available by surrounding routines
IF (NL == NLEVEL_MIN) THEN

   L%WALL_INDEX_PTR => M%WALL_INDEX

   ! if on coarser levels, it must still be computed
ELSE

   CALL SCARC_ALLOCATE_INT2(L%WALL_INDEX, 1, L%N_CELL_INDEX, -3, 3, NSCARC_INIT_ZERO, 'L%WALL_INDEX')
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

! -----------------------------------------------------------------------------------------
! In case of an MG-method:
! Setup ghost cells order 
! -----------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_WALL_GHOSTS(NM, NL)
USE SCARC_POINTERS, ONLY: OL, OG, GWC
INTEGER, INTENT(IN) :: NM, NL
INTEGER :: IW, IOR0, NOM

CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G
DO IW = 1, L%N_WALL_CELLS

  GWC => G%WALL(IW)
  NOM  = GWC%NOM
  IF (NOM == 0) CYCLE

  IOR0 = GWC%IOR
  CALL SCARC_POINT_TO_OTHER_GRID(NM, NOM, NL)

  OG%ICG2 = OG%ICG2 + 1
  IF (OL%GHOST_FIRSTW(IOR0) == 0) OL%GHOST_FIRSTW(IOR0) = OG%ICG2
  IF (OL%GHOST_FIRSTE(IOR0) == 0) OL%GHOST_FIRSTE(IOR0) = OG%ICG2
  OL%GHOST_LASTE(IOR0) = OG%ICG2

ENDDO

END SUBROUTINE SCARC_SETUP_WALL_GHOSTS

! -------------------------------------------------------------------------------------------------
! Setup all necessary information for a wall cell with neighbor in case of MG-method
! Number of obstructions on coarse level is the same as on fine level
! TODO: Only works for special cases which run for GMG, must still be extended!!
! -------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_WALL_COORDS(NM, NL)
USE SCARC_POINTERS, ONLY: L, G, OB
INTEGER, INTENT(IN) :: NM, NL
INTEGER :: IC, IO, IWC
INTEGER :: I, J, K

CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G
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
! Correct BTYPE related to internal obstructions on ghost cells
! TODO: to check again !!!
! -------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_IDENTIFY_INTERNAL_NEUMANNS(NM, NL)
USE SCARC_POINTERS, ONLY: L, G, GWC
INTEGER, INTENT(IN) :: NM, NL
INTEGER :: IWG
INTEGER :: IX, IY, IZ, IOR0, BTYPE0

CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G

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
! Setup all necessary information for a wall cell with neighbor
! -------------------------------------------------------------------------------------------------
INTEGER FUNCTION SCARC_COUNT_EXTERNAL_WALL_CELLS(NM, NL)
USE SCARC_POINTERS, ONLY: LF, LC
INTEGER, INTENT(IN) :: NM, NL
INTEGER :: IXC, IYC, IZC
INTEGER :: IXF, IYF, IZF
INTEGER :: NLF, NLC
INTEGER :: IWC, ICF(4)=0, IWF(4)=0, IOR0

NLF = NL-1
NLC = NL

LF => SCARC(NM)%LEVEL(NLF)
LC => SCARC(NM)%LEVEL(NLC)

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
      IF (IS_EXTERNAL_WALLCELL(IOR0, ICF, 2, NM, NLF)) IWC = IWC + 1
   ENDDO

   ! IOR = -1
   IOR0 = -1
   DO IZC = 1, LC%NZ
      IZF = 2*IZC - 1
      ICF(1) = LF%CELL_INDEX_PTR(LF%NX, IYF  , IZF  )
      ICF(2) = LF%CELL_INDEX_PTR(LF%NX, IYF  , IZF+1)
      IF (IS_EXTERNAL_WALLCELL(IOR0, ICF, 2, NM, NLF)) IWC = IWC + 1
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
         IF (IS_EXTERNAL_WALLCELL(IOR0, ICF, 4, NM, NLF)) IWC = IWC + 1
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
         IF (IS_EXTERNAL_WALLCELL(IOR0, ICF, 4, NM, NLF)) IWC = IWC + 1
      ENDDO
   ENDDO

   ! IOR = 3
   IOR0 = 3
   DO IXC = 1, LC%NX
      IXF = 2*IXC - 1
      ICF(1) = LF%CELL_INDEX_PTR(IXF    , IYF  , 1)
      ICF(2) = LF%CELL_INDEX_PTR(IXF+1  , IYF  , 1)
      IF (IS_EXTERNAL_WALLCELL(IOR0, ICF, 2, NM, NLF)) IWC = IWC + 1
   ENDDO

   ! IOR = -3
   IOR0 = -3
   DO IXC = 1, LC%NX
      IXF = 2*IXC - 1
      ICF(1) = LF%CELL_INDEX_PTR(IXF  , IYF  , LF%NZ)
      ICF(2) = LF%CELL_INDEX_PTR(IXF+1, IYF  , LF%NZ)
      IF (IS_EXTERNAL_WALLCELL(IOR0, ICF, 2, NM, NLF)) IWC = IWC + 1
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
         IF (IS_EXTERNAL_WALLCELL(IOR0, ICF, 4, NM, NLF)) IWC = IWC + 1
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
         IF (IS_EXTERNAL_WALLCELL(IOR0, ICF, 4, NM, NLF)) IWC = IWC + 1
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
         IF (IS_EXTERNAL_WALLCELL(IOR0, ICF, 4, NM, NLF)) IWC = IWC + 1
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
         IF (IS_EXTERNAL_WALLCELL(IOR0, ICF, 4, NM, NLF)) IWC = IWC + 1
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
         IF (IS_EXTERNAL_WALLCELL(IOR0, ICF, 4, NM, NLF)) IWC = IWC + 1
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
         IF (IS_EXTERNAL_WALLCELL(IOR0, ICF, 4, NM, NLF)) IWC = IWC + 1
      ENDDO
   ENDDO

ENDIF

SCARC_COUNT_EXTERNAL_WALL_CELLS = IWC
END FUNCTION SCARC_COUNT_EXTERNAL_WALL_CELLS


! -------------------------------------------------------------------------------------------------
! Setup all necessary information for a wall cell with neighbor
! -------------------------------------------------------------------------------------------------
INTEGER FUNCTION SCARC_COUNT_INTERNAL_WALL_CELLS(NM, NL)
USE SCARC_POINTERS, ONLY: LF, LC, GF, GC, OB
INTEGER, INTENT(IN) :: NM, NL
INTEGER :: IWC, NW_INT
INTEGER :: IC, IO
INTEGER :: I, J, K

LC => SCARC(NM)%LEVEL(NL)
LF => SCARC(NM)%LEVEL(NL-1)

SELECT CASE(TYPE_GRID)
   CASE (NSCARC_GRID_STRUCTURED)
      GC => LC%STRUCTURED
      GF => LF%STRUCTURED
   CASE (NSCARC_GRID_UNSTRUCTURED)
      GC => LC%UNSTRUCTURED
      GF => LF%UNSTRUCTURED
END SELECT

LC%N_OBST = LF%N_OBST                   ! Number of obstructions is the same on all levels

ALLOCATE(LC%OBST(LC%N_OBST), STAT=IERROR)
CALL ChkMemErr('SCARC_COUNT_INTERNAL_WALL_CELLS','OBST',IERROR)

NW_INT = 0
IWC = LC%N_WALL_CELLS_EXT + 1

DO IO = 1, LF%N_OBST

   OB => SCARC(NM)%LEVEL(NL)%OBST(IO)

   OB%I1 = (LF%OBST(IO)%I1+1)/2
   OB%I2 =  LF%OBST(IO)%I2/2

   IF (TWO_D) THEN
      OB%J1 = 0
      OB%J2 = 1
   ELSE
      OB%J1 = (LF%OBST(IO)%J1+1)/2
      OB%J2 =  LF%OBST(IO)%J2/2
   ENDIF

   OB%K1 = (LF%OBST(IO)%K1+1)/2
   OB%K2 =  LF%OBST(IO)%K2/2

   ! Analyze IOR = 1
   I = OB%I1
   DO K = OB%K1+1, OB%K2
      DO J = OB%J1+1, OB%J2
         IC = GC%CELL_NUMBER(I, J, K)
         IF (IC == NSCARC_UNDEF_INT) CYCLE
         NW_INT = NW_INT + 1
      ENDDO
   ENDDO

   ! Analyze IOR = -1
   I = OB%I2
   DO K = OB%K1+1, OB%K2
      DO J = OB%J1+1, OB%J2
         IC = GC%CELL_NUMBER(I+1, J, K)
         IF (IC == NSCARC_UNDEF_INT) CYCLE
         NW_INT = NW_INT + 1
      ENDDO
   ENDDO

   ! Analyze IOR = 2
   J = OB%J1
   DO K = OB%K1+1, OB%K2
      DO I = OB%I1+1, OB%I2
         IC = GC%CELL_NUMBER(I, J, K)
         IF (IC == NSCARC_UNDEF_INT) CYCLE
         NW_INT = NW_INT + 1
      ENDDO
   ENDDO

   ! Analyze IOR = -2
   J = OB%J2
   DO K = OB%K1+1, OB%K2
      DO I = OB%I1+1, OB%I2
         IC = GC%CELL_NUMBER(I, J+1, K)
         IF (IC == NSCARC_UNDEF_INT) CYCLE
         NW_INT = NW_INT + 1
      ENDDO
   ENDDO

   ! Analyze IOR = 3
   K = OB%K1
   DO J = OB%J1+1, OB%J2
      DO I = OB%I1+1, OB%I2
         IC = GC%CELL_NUMBER(I, J, K)
         IF (IC == NSCARC_UNDEF_INT) CYCLE
         NW_INT = NW_INT + 1
      ENDDO
   ENDDO

   ! Analyze IOR = -3
   K = OB%K2
   DO J = OB%J1+1, OB%J2
      DO I = OB%I1+1, OB%I2
         IC = GC%CELL_NUMBER(I, J, K+1)
         IF (IC == NSCARC_UNDEF_INT) CYCLE
         NW_INT = NW_INT + 1
      ENDDO
   ENDDO

ENDDO

SCARC_COUNT_INTERNAL_WALL_CELLS = NW_INT
END FUNCTION SCARC_COUNT_INTERNAL_WALL_CELLS


! -------------------------------------------------------------------------------------------------
! Count external wall cells on face IOR
! -------------------------------------------------------------------------------------------------
LOGICAL FUNCTION IS_EXTERNAL_WALLCELL(IOR0, ICF, NCNT, NM, NL)
USE SCARC_POINTERS, ONLY: L, G
INTEGER, INTENT(IN) :: IOR0, NCNT, NM, NL
INTEGER, DIMENSION(:), INTENT(IN) :: ICF
INTEGER :: I, IWF_LAST, IWF(4)=0
REAL(EB) :: BSUM

IS_EXTERNAL_WALLCELL = .FALSE.

CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G

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
! Setup all necessary information for a wall cell with neighbor
! -------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_WALL_NEIGHBOR(NX1, NX2, NY1, NY2, NZ1, NZ2, IWG, NM, NOM, NL)
USE SCARC_POINTERS, ONLY: G, OG
INTEGER, INTENT(IN) :: NX1, NX2, NY1, NY2, NZ1, NZ2
INTEGER, INTENT(IN) :: IWG, NM, NOM, NL
INTEGER :: NOMX, NOMY, NOMZ
INTEGER :: ICG, ICE, IX, IY, IZ, JL, IXW, IYW, IZW, IXG, IYG, IZG

CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G
CALL SCARC_POINT_TO_OTHER_GRID(NM, NOM, NL)

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
! Check divisibility by 2 of a given number of elements (in one grid direction)
! ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_CHECK_DIVISIBILITY(NN, CDIR)
INTEGER, INTENT(IN) :: NN
CHARACTER(*) , INTENT(IN) :: CDIR
IF (MOD(NN,2) /= 0) CALL SCARC_SHUTDOWN(NSCARC_ERROR_GRID_NUMBER, CDIR, NSCARC_NONE)
END SUBROUTINE SCARC_CHECK_DIVISIBILITY


! ----------------------------------------------------------------------------------------------------
! Set wall cell information on coarse level
! ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_WALL_LEVEL(IOR0, IWC, IREFINE, NM, NL)
USE SCARC_POINTERS, ONLY: LF, LC, FF, WF, WC, OLC, OGC
INTEGER, INTENT(INOUT) :: IWC
INTEGER, INTENT(IN) :: NM, NL
INTEGER, INTENT(IN) :: IOR0, IREFINE
INTEGER :: IWF(4) , IBCF(4), NOMF(4)
INTEGER :: IX,  IY,  IZ, I
INTEGER :: NX1, NY1, NZ1
INTEGER :: NX2, NY2, NZ2
INTEGER :: IX1, IY1, IZ1
INTEGER :: IX2, IY2, IZ2
INTEGER :: IDIFF, JDIFF, KDIFF


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

!
! Loop over all wall cells of face IOR0
!
DO IZ = NZ1, NZ2
   DO IY = NY1, NY2
      DO IX = NX1, NX2

         ! Set orientation of neiboring face, indices of ghost and adjacent cell for coarse IW
         WC(IWC)%IOR = IOR0

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
         ! 2D-version
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

            ! in case of an internal boundary set neighboring WALL cells
            IF (NOMF(1) > 0) THEN

               OLC => SCARC(NM)%OSCARC(NOMF(1))%LEVEL(NL)
               SELECT CASE(TYPE_GRID)
                  CASE (NSCARC_GRID_STRUCTURED)
                     OGC => OLC%STRUCTURED
                  CASE (NSCARC_GRID_UNSTRUCTURED)
                     OGC => OLC%UNSTRUCTURED
               END SELECT

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

               !
               ! Allocate and specify ICN and ICE arrays for OC
               !
               CALL SCARC_SETUP_WALL_NEIGHBOR(IX1, IX2, 1, 1, IZ1, IZ2, IWC, NM, NOMF(1), NL)

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

               OLC => SCARC(NM)%OSCARC(NOMF(1))%LEVEL(NL)

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

               CALL SCARC_SETUP_WALL_NEIGHBOR(IX1, IX2, IY1, IY2, IZ1, IZ2, IWC, NM, NOMF(1), NL)

            ENDIF
         ENDIF
         IWC = IWC + 1
      ENDDO
   ENDDO
ENDDO

END SUBROUTINE SCARC_SETUP_WALL_LEVEL


! -------------------------------------------------------------------------------------------------
! Set analogues to I_MIN, IMAX, J_MIN, J_MAX, K_MIN, K_MAX on coarse level (only GMG!)
! -------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_EXCHANGE_DIMENSIONS(IREFINE, NM, NOM, NL)
USE SCARC_POINTERS, ONLY: G, OG
INTEGER, INTENT(IN) :: IREFINE, NM, NOM, NL
INTEGER :: IMIN, IMAX, JMIN, JMAX, KMIN, KMAX, IW
LOGICAL :: FOUND

CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G
CALL SCARC_POINT_TO_OTHER_GRID(NM, NOM, NL)

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

   ! neighborship structure already known from finest level
   IF (G%WALL(IW)%NOM/=NOM) CYCLE SEARCH_LOOP
   FOUND = .TRUE.

   SELECT CASE (G%WALL(IW)%IOR)
      CASE ( 1)
         IMIN=MAX(IMIN,G%WALL(NM)%IXN(1)-1)
      CASE (-1)
         IMAX=MIN(IMAX,G%WALL(NM)%IXN(2)+1)
      CASE ( 2)
         JMIN=MAX(JMIN,G%WALL(NM)%IYN(1)-1)
      CASE (-2)
         JMAX=MIN(JMAX,G%WALL(NM)%IYN(2)+1)
      CASE ( 3)
         KMIN=MAX(KMIN,G%WALL(NM)%IZN(1)-1)
      CASE (-3)
         KMAX=MIN(KMAX,G%WALL(NM)%IZN(2)+1)
   END SELECT
ENDDO SEARCH_LOOP

N_EXCHANGES = N_EXCHANGES+1

END SUBROUTINE SCARC_SETUP_EXCHANGE_DIMENSIONS


! ----------------------------------------------------------------------------------------------------
! Allocate several global structures for data exchange
! --------------------------------------------------------------------------------------------------ss
SUBROUTINE SCARC_SETUP_GLOBALS
INTEGER :: NM, NP

IF (N_MPI_PROCESSES > 1) THEN

   ! Allocate and preset counter and displacement vector for global data exchanges
   CALL SCARC_ALLOCATE_INT1 (COUNTS, 0, N_MPI_PROCESSES-1, NSCARC_INIT_ZERO, 'COUNTS')
   CALL SCARC_ALLOCATE_INT1 (DISPLS, 0, N_MPI_PROCESSES-1, NSCARC_INIT_ZERO, 'DISPLS')

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

CALL SCARC_ALLOCATE_INT1 (MESH_INT , 1, NMESHES, NSCARC_INIT_ZERO, 'MESH_INT')
CALL SCARC_ALLOCATE_REAL1(MESH_REAL, 1, NMESHES, NSCARC_INIT_ZERO, 'MESH_REAL')

CALL SCARC_SETUP_DIMENSIONS(NLEVEL_MIN)

END SUBROUTINE SCARC_SETUP_GLOBALS


! -----------------------------------------------------------------------------
! Get information about global numbers of unknowns for unstructured case
! -----------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_DIMENSIONS(NL)
USE SCARC_POINTERS, ONLY: G
INTEGER, INTENT(IN) :: NL
INTEGER :: NM, NM2

! Preset communication array MESH_INT with local numbers of cells for all meshes depending on type of discretization
MESHES_LOOP1: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G

   DO NM2 = LOWER_MESH_INDEX, UPPER_MESH_INDEX
      MESH_INT(NM2) = G%NC_LOCAL(NM2)
   ENDDO

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
! Setup system of equation:
! Define matrix stencils and initialize matrices and boundary conditions on all needed levels
! ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_SYSTEMS
INTEGER :: NM, NOM, NL, INBR


! 
! ------ Setup sizes for system matrices
! 
! Get mapping array from local to global numbering
CALL SCARC_SETUP_GLOBAL_NUMBERING(NLEVEL_MIN)

SELECT_SCARC_METHOD_SIZES: SELECT CASE (TYPE_METHOD)

   ! -------- Global Krylov method
   CASE (NSCARC_METHOD_KRYLOV)
   
      CALL SCARC_ASSIGN_GRID_TYPE (TYPE_GRID)                      ! process specified discretization type
      CALL SCARC_SETUP_POISSON_SIZES (NLEVEL_MIN)                  ! setup sizes on finest level
   
      IF (HAS_TWO_LEVELS .AND. .NOT. IS_CG_SAMG) &                                                  
         CALL SCARC_SETUP_POISSON_SIZES (NLEVEL_MAX)               ! twolevel-precon: also setup size for coarse level
   
      IF (IS_CG_GMG) THEN                                                   
         DO NL=NLEVEL_MIN+1, NLEVEL_MAX
            CALL SCARC_SETUP_POISSON_SIZES (NL)                    ! GMG-precon: also setup size for all other levels
         ENDDO
      ENDIF
   
   ! -------- Global Multigrid method
   CASE (NSCARC_METHOD_MULTIGRID)
   
      CALL SCARC_ASSIGN_GRID_TYPE (TYPE_GRID)                      ! process specified discretization type
      SELECT CASE (TYPE_MULTIGRID)
         CASE (NSCARC_MULTIGRID_GEOMETRIC)                                   
            DO NL=NLEVEL_MIN, NLEVEL_MAX
               CALL SCARC_SETUP_POISSON_SIZES (NL)                 ! GMG: setup size for all levels
            ENDDO
         CASE (NSCARC_MULTIGRID_ALGEBRAIC)                                   
            CALL SCARC_SETUP_POISSON_SIZES (NLEVEL_MIN)            ! AMG: setup sizes only on finest level
      END SELECT
   
#ifdef WITH_SCARC_MGM
   ! -------- Global MGM method - currently just proof of concept
   CASE (NSCARC_METHOD_MGM)
   
      CALL SCARC_ASSIGN_GRID_TYPE (NSCARC_GRID_STRUCTURED)         ! First process structured discretization
      CALL SCARC_SETUP_POISSON_SIZES (NLEVEL_MIN)        
   
      CALL SCARC_ASSIGN_GRID_TYPE (NSCARC_GRID_UNSTRUCTURED)       ! Then process unstructured discretization
      CALL SCARC_SETUP_POISSON_SIZES (NLEVEL_MIN)         
   
#endif

END SELECT SELECT_SCARC_METHOD_SIZES


! 
! ------ Assemble system matrices on requested grid levels and set boundary conditions
! 
MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   SELECT_SCARC_METHOD: SELECT CASE (TYPE_METHOD)

      ! ---------- Krylov method (CG/BICG) as main solver, different preconditioners possible
      CASE (NSCARC_METHOD_KRYLOV)

         SELECT_KRYLOV_PRECON: SELECT CASE (TYPE_PRECON)

            ! in case of multigrid as preconditioner:
            CASE (NSCARC_RELAX_GMG)

               SELECT_KRYLOV_PRECON_MG: SELECT CASE (TYPE_MULTIGRID)

                  ! Geometric multigrid:
                  !    -  assemble standard n-point-matrix hierarchy on all levels (finest and all coarser)
                  CASE (NSCARC_MULTIGRID_GEOMETRIC)

                     DO NL = NLEVEL_MIN, NLEVEL_MAX
                        CALL SCARC_SETUP_POISSON (NM, NL)
                        CALL SCARC_SETUP_BOUNDARY(NM, NL)
#ifdef WITH_MKL
                        IF (TYPE_MKL(NL) == NSCARC_MKL_LOCAL .OR. TYPE_MKL(NL) == NSCARC_MKL_GLOBAL) &
                           CALL SCARC_SETUP_POISSON_MKL(NM, NL)
#endif

                     ENDDO

                  ! Algebraic multigrid:
                  !    -  use compact storage technique on all levels (no other choise possible!)
                  !    -  assemble standard n-point-matrix only on finest level
                  !    -  construct all coarser levels by requested coarsening strategy
                  CASE (NSCARC_MULTIGRID_ALGEBRAIC)
   
                     CALL SCARC_SETUP_POISSON (NM, NLEVEL_MIN)
                     CALL SCARC_SETUP_BOUNDARY(NM, NLEVEL_MIN)

               END SELECT SELECT_KRYLOV_PRECON_MG

#ifdef WITH_MKL
            ! in case of LU-decomposition as preconditioner
            CASE (NSCARC_RELAX_MKL)

               SELECT_KRYLOV_PRECON_MKL: SELECT CASE(TYPE_SCOPE(1))

                  ! Locally acting: PARDISO from MKL as preconditioners with possible coarse grid correction
                  CASE (NSCARC_SCOPE_LOCAL)
   
                     CALL SCARC_SETUP_POISSON (NM, NLEVEL_MIN)
                     CALL SCARC_SETUP_BOUNDARY(NM, NLEVEL_MIN)
                     CALL SCARC_SETUP_POISSON_MKL (NM, NLEVEL_MIN)
   
                     IF (HAS_TWO_LEVELS .AND. .NOT.IS_CG_AMG .AND. .NOT.IS_CG_SAMG) THEN
                        CALL SCARC_SETUP_POISSON (NM, NLEVEL_MAX)
                        CALL SCARC_SETUP_BOUNDARY(NM, NLEVEL_MAX)
                        IF (TYPE_COARSE == NSCARC_COARSE_DIRECT) CALL SCARC_SETUP_POISSON_MKL (NM, NLEVEL_MAX)
                     ENDIF

                  ! Globally acting: Cluster_Sparse_Solver from MKL as preconditioner
                  CASE (NSCARC_SCOPE_GLOBAL)
   
                     CALL SCARC_SETUP_POISSON (NM, NLEVEL_MIN)
                     CALL SCARC_SETUP_BOUNDARY(NM, NLEVEL_MIN)
                     CALL SCARC_SETUP_POISSON_MKL (NM, NLEVEL_MIN)
   
               END SELECT SELECT_KRYLOV_PRECON_MKL
#endif

            ! in case of one-level preconditioners (JACOBI/SSOR/FFT):
            ! assemble standard n-point-matrix on finest level with possible coarse grid correction
            CASE DEFAULT
   
               CALL SCARC_SETUP_POISSON (NM, NLEVEL_MIN)
               CALL SCARC_SETUP_BOUNDARY(NM, NLEVEL_MIN)

               IF (HAS_TWO_LEVELS .AND. .NOT.IS_CG_SAMG) THEN
                  CALL SCARC_SETUP_POISSON (NM, NLEVEL_MAX)
                  CALL SCARC_SETUP_BOUNDARY(NM, NLEVEL_MAX)
#ifdef WITH_MKL
                  IF (TYPE_COARSE == NSCARC_COARSE_DIRECT) CALL SCARC_SETUP_POISSON_MKL (NM, NL)
#endif
               ENDIF

         END SELECT SELECT_KRYLOV_PRECON


      ! ---------- Multigrid as main solver
      CASE (NSCARC_METHOD_MULTIGRID)
   
         SELECT_MULTIGRID_TYPE: SELECT CASE (TYPE_MULTIGRID)
   
            !  Geometric multigrid:
            !    -  assemble standard n-point-matrix hierarchy on all levels
            !    -  use MKL coarse grid solver if requested
            CASE (NSCARC_MULTIGRID_GEOMETRIC)

               DO NL = NLEVEL_MIN, NLEVEL_MAX
                  CALL SCARC_SETUP_POISSON (NM, NL)
                  CALL SCARC_SETUP_BOUNDARY(NM, NL)
               ENDDO
#ifdef WITH_MKL
               DO NL = NLEVEL_MIN, NLEVEL_MAX-1
                  CALL SCARC_SETUP_POISSON_MKL(NM, NL)
               ENDDO
               IF (TYPE_COARSE == NSCARC_COARSE_DIRECT) CALL SCARC_SETUP_POISSON_MKL(NM, NLEVEL_MAX)
#endif

            !  Algebraic multigrid:
            !    -  use compact storage technique (no other choice possible!)
            !    -  assemble standard n-point-matrix only on finest level
            !    -  construct all coarser levels later by requested coarsening strategy
            CASE (NSCARC_MULTIGRID_ALGEBRAIC)

               CALL SCARC_SETUP_POISSON (NM, NLEVEL_MIN)
               CALL SCARC_SETUP_BOUNDARY(NM, NLEVEL_MIN)

         END SELECT SELECT_MULTIGRID_TYPE

#ifdef WITH_MKL
      ! ---------- MKL-LU decomposition as main solver
      CASE (NSCARC_METHOD_LU)
  
         CALL SCARC_SETUP_POISSON (NM, NLEVEL_MIN)
         CALL SCARC_SETUP_BOUNDARY(NM, NLEVEL_MIN)
         CALL SCARC_SETUP_POISSON_MKL(NM, NL)
#endif

#ifdef WITH_SCARC_MGM
      ! ---------- McKenny-Greengard-Mayo method:
      ! Solving for the structured and unstructured Poisson matrix
      ! Assemble both, the structured and unstructured Poisson matrix
      ! temporarily they will be stored separately in matrices AC and ACU due to the different
      ! settings along internal boundary cells,
      ! in the medium term, a toggle mechanism will be implemented which only switches the corresponding
      ! entries while keeping the entries which are the same for both discretization types
      CASE (NSCARC_METHOD_MGM)
   
         ! First assemble unstructured matrix
         CALL SCARC_ASSIGN_GRID_TYPE (NSCARC_GRID_UNSTRUCTURED)
         CALL SCARC_SETUP_POISSON (NM, NLEVEL_MIN)
         CALL SCARC_SETUP_BOUNDARY(NM, NLEVEL_MIN)

         ! Then assemble structured matrix
         CALL SCARC_ASSIGN_GRID_TYPE (NSCARC_GRID_STRUCTURED)
         CALL SCARC_SETUP_POISSON (NM, NLEVEL_MIN)
         CALL SCARC_SETUP_BOUNDARY(NM, NLEVEL_MIN)
#endif

   END SELECT SELECT_SCARC_METHOD

ENDDO MESHES_LOOP

! Setup global numbering of matrix column relationships (Compact storage technique only)
CALL SCARC_SETUP_GLOBAL_MATRIX_COLUMNS(NSCARC_MATRIX_POISSON, NLEVEL_MIN)
 

! If there is more than one mesh, exchange matrix values in overlapping parts
! This must be done for all multilevel methods at least at the finest grid level
! Furthermore also at all higher levels except for the AMG method,
! in this case it will be done later in routine SETUP_SAMG
NMESHES_IF: IF (NMESHES > 1) THEN

   TYPE_MATRIX_IF: IF (TYPE_MATRIX == NSCARC_MATRIX_COMPACT) THEN

      CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_POISSON_COLS,  NSCARC_MATRIX_POISSON, NLEVEL_MIN)
      CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_POISSON_COLSG, NSCARC_MATRIX_POISSON, NLEVEL_MIN)
      CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_POISSON_VALS,  NSCARC_MATRIX_POISSON, NLEVEL_MIN)

      CALL SCARC_EXTRACT_MATRIX_OVERLAPS(NSCARC_MATRIX_POISSON, NLEVEL_MIN)

      MESHES_FINE_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

         CALL SCARC_POINT_TO_GRID (NM, NLEVEL_MIN)    
         A => SCARC_POINT_TO_CMATRIX (G, NSCARC_MATRIX_POISSON)
         CALL SCARC_REDUCE_CMATRIX(A, 'G%POISSONC')

         OMESHES_FINE_LOOP: DO INBR = 1, S%N_NEIGHBORS
            NOM = S%NEIGHBORS(INBR)
            CALL SCARC_POINT_TO_OTHER_GRID(NM, NOM, NLEVEL_MIN)
            OA => SCARC_POINT_TO_OTHER_CMATRIX (OG, NSCARC_MATRIX_POISSON)
            CALL SCARC_REDUCE_CMATRIX(OA, 'OG%POISSONC')
         ENDDO OMESHES_FINE_LOOP

#ifdef WITH_SCARC_DEBUG
         CALL SCARC_DEBUG_CMATRIX(A, 'G%POISSONC', 'SETUP_SYSTEMS: AFTER EXTRAXT_MATRIX_OVERLAPS')
#endif

      ENDDO MESHES_FINE_LOOP
    

      MULTI_LEVEL_IF: IF (HAS_MULTI_LEVELS .AND. .NOT.IS_AMG .AND. .NOT.IS_CG_SAMG) THEN

         DO NL = NLEVEL_MIN+1, NLEVEL_MAX

            CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_POISSON_COLS,  NSCARC_MATRIX_POISSON, NL)
            CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_POISSON_COLSG, NSCARC_MATRIX_POISSON, NL)
            CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_POISSON_VALS,  NSCARC_MATRIX_POISSON, NL)
   
            MESHES_COARSE_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

               CALL SCARC_POINT_TO_GRID (NM, NL)    
               A => SCARC_POINT_TO_CMATRIX (G, NSCARC_MATRIX_POISSON)
               CALL SCARC_REDUCE_CMATRIX(A, 'G%POISSONC')

               OMESHES_COARSE_LOOP: DO INBR = 1, S%N_NEIGHBORS
                  NOM = S%NEIGHBORS(INBR)
                  CALL SCARC_POINT_TO_OTHER_GRID(NM, NOM, NL)
                  OA => SCARC_POINT_TO_OTHER_CMATRIX (OG, NSCARC_MATRIX_POISSON)
                  CALL SCARC_REDUCE_CMATRIX(OA, 'OG%POISSONC')
               ENDDO OMESHES_COARSE_LOOP

#ifdef WITH_SCARC_DEBUG
               CALL SCARC_DEBUG_CMATRIX(A, 'G%POISSONC', 'SETUP_SYSTEMS: AFTER EXTRAXT_MATRIX_OVERLAPS')
#endif

            ENDDO MESHES_COARSE_LOOP
         ENDDO

      ENDIF MULTI_LEVEL_IF
   ENDIF TYPE_MATRIX_IF
ENDIF NMESHES_IF


! Debug matrix and wall structures - only if directive SCARC_DEBUG is set
#ifdef WITH_SCARC_DEBUG
CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_WALL  , NLEVEL_MIN, 'WALL')
CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_FACE  , NLEVEL_MIN, 'FACE_AFTER_SYSTEM')
CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_CMATRIX, NLEVEL_MIN, 'SYSTEM-MATRIX')
#endif

END SUBROUTINE SCARC_SETUP_SYSTEMS


! -------------------------------------------------------------------------------------------
! Setup mapping vector from local to global numbering
! -------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_GLOBAL_NUMBERING(NL)
USE SCARC_POINTERS, ONLY : G
INTEGER, INTENT(IN) :: NL
INTEGER :: NM, NOM, IC, IW, ICE, ICN

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID (NM, NL)             ! set grid pointer G

   IF (IS_AMG .OR. IS_CG_SAMG) THEN              ! AMG: consider two overlapping layers
      CALL SCARC_ALLOCATE_INT1 (G%LOCAL_TO_GLOBAL, 1, G%NCE2, NSCARC_INIT_ZERO, 'G%LOCAL_TO_GLOBAL')
   ELSE                                          ! else: only consider one overlapping layer
      CALL SCARC_ALLOCATE_INT1 (G%LOCAL_TO_GLOBAL, 1, G%NCE , NSCARC_INIT_ZERO, 'G%LOCAL_TO_GLOBAL')
   ENDIF

   DO IC = 1, G%NC
      G%LOCAL_TO_GLOBAL(IC) = IC + G%NC_OFFSET(NM)
   ENDDO

   DO IW = 1, G%NW
      NOM = G%WALL(IW)%NOM
      IF (NOM == 0) CYCLE
      ICE = G%WALL(IW)%ICE
      ICN = G%ICE_TO_ICN(ICE)
      G%LOCAL_TO_GLOBAL(ICE) = ICN + G%NC_OFFSET(NOM)
   ENDDO

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'SETUP_GLOBAL_NUMBERING: LOCAL_TO_GLOBAL'
WRITE(MSG%LU_DEBUG,'(8I6)') G%LOCAL_TO_GLOBAL
#endif

ENDDO

END SUBROUTINE SCARC_SETUP_GLOBAL_NUMBERING


! -------------------------------------------------------------------------------------------
! Extract overlapping matrix part and add to main matrix
! -------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_GLOBAL_MATRIX_COLUMNS(NMATRIX, NL)
USE SCARC_POINTERS, ONLY : G
INTEGER, INTENT(IN) :: NMATRIX, NL
INTEGER :: NM, IC, ICOL, JC

IF (NMESHES == 1) THEN
   CALL SCARC_POINT_TO_GRID (NMESHES, NL)                                   ! Sets grid pointer G
   A => SCARC_POINT_TO_CMATRIX (G, NMATRIX)
   A%COLG = A%COL
ELSE
   DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
      CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G
      A => SCARC_POINT_TO_CMATRIX (G, NMATRIX)
      DO IC = 1, G%NC
         DO ICOL = A%ROW(IC), A%ROW(IC+1) -1
            JC = A%COL(ICOL)
            !IF (JC <= G%NC) THEN
               A%COLG(ICOL) = G%LOCAL_TO_GLOBAL(JC)
            !ELSE
            !NDIF
         ENDDO
      ENDDO
   ENDDO
ENDIF

END SUBROUTINE SCARC_SETUP_GLOBAL_MATRIX_COLUMNS


! -------------------------------------------------------------------------------------------
! Extract overlapping matrix part and add to main matrix
! -------------------------------------------------------------------------------------------
SUBROUTINE SCARC_EXTRACT_MATRIX_OVERLAPS (NMATRIX, NL)
USE SCARC_POINTERS, ONLY : G, F, OL, OG, A
INTEGER, INTENT(IN) :: NL, NMATRIX
INTEGER :: NM, IFACE, NOM, IOR0, ICG, ICE, IP, ICOL, INBR, ICN, ICE1, ICE2

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G

   SELECT CASE (TYPE_GRID)

      ! ---------- Version for structured grids
      CASE (NSCARC_GRID_STRUCTURED)

         A => SCARC_POINT_TO_CMATRIX(G, NMATRIX)
 
         IP = A%ROW(G%NC+1)
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'EXTRACT_MATRIX_OVERLAPS: NM=', NM,': IP=',IP, G%POISSONC%N_ROW
#endif
         !==============================
         !==============================
         !==============================
         IF ( (NMATRIX /= NSCARC_MATRIX_POISSON) .OR. &
              (NMATRIX == NSCARC_MATRIX_POISSON .AND. .NOT.((IS_AMG.OR.IS_CG_SAMG) .AND. NL>NLEVEL_MIN))) THEN

         FACE_LOOP: DO IFACE = 1, 6                          ! check if this face has connection to last cell

            IOR0 = FACE_ORIENTATION(IFACE)
            F => SCARC(NM)%LEVEL(NLEVEL_MIN)%FACE(IOR0)
      
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'IFACE, IOR0:', IFACE, IOR0
#endif
            DO INBR = 1, F%N_NEIGHBORS
      
               NOM = F%NEIGHBORS(INBR)
               CALL SCARC_POINT_TO_OTHER_GRID(NM, NOM, NL)
               OA => SCARC_POINT_TO_OTHER_CMATRIX(OG, NMATRIX)
   
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,'(A,5I6)') 'IFACE, IOR0, INBR, NOM:', IFACE, IOR0, INBR, NOM, G%POISSONC%N_ROW
WRITE(MSG%LU_DEBUG,*) 'G%POISSONC%ROW:'
WRITE(MSG%LU_DEBUG,'( 8I6)') G%POISSONC%ROW(1:G%POISSONC%N_ROW)
WRITE(MSG%LU_DEBUG,*) 'OA%COL'
WRITE(MSG%LU_DEBUG,'( 8I6)') OA%COL
WRITE(MSG%LU_DEBUG,*) 'OA%VAL'
WRITE(MSG%LU_DEBUG,'(8E12.4)') OA%VAL
#endif
               ICOL = 1
               DO ICG = OL%GHOST_FIRSTE(IOR0), OL%GHOST_LASTE(IOR0)
   
                  ICE = OG%ICG_TO_ICE(ICG, 1)
                  A%ROW(ICE) = IP 

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) '====================== ICG, ICE, IP =', ICG, ICE, IP
#endif
                  IF (NMATRIX == NSCARC_MATRIX_POISSON .AND. .NOT.((IS_AMG.OR.IS_CG_SAMG) .AND. NL>NLEVEL_MIN)) THEN

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
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,'(A,7I6,E12.4)') '       ----> B: ICG, ICE, ICOL, IP, ROW, COL, VAL :', &
                                       ICG, ICE, ICOL, IP, A%ROW(ICE), A%COL(IP), A%COLG(IP), A%VAL(IP)
#endif
                        IP = IP + 1
                     ENDDO

                  ELSE

                     DO ICOL = OA%ROW(ICG), OA%ROW(ICG+1)-1
                        A%COL(IP) = -OA%COL(ICOL)   
                        A%COLG(IP) = ABS(OA%COLG(ICOL))      
                        A%VAL(IP) = OA%VAL(ICOL)
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,'(A,7I6,E12.4)') '       ----> B: ICG, ICE, ICOL, IP, ROW, COL, VAL :', &
                                       ICG, ICE, ICOL, IP, A%ROW(ICE), A%COL(IP), A%COLG(IP), A%VAL(IP)
#endif
                        IP = IP + 1
                     ENDDO
                  ENDIF
   
               ENDDO

               A%ROW(ICE+1) = IP 
               A%N_ROW = ICE + 1
               A%N_VAL = IP - 1
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'FINAL: ICE+1, ROW, N_ROW, N_VAL:', ICE+1, IP, ICE+1, IP -1
WRITE(MSG%LU_DEBUG,*) 'FINAL: A%ROW:'
WRITE(MSG%LU_DEBUG,'(9I6)') A%ROW
WRITE(MSG%LU_DEBUG,*) 'FINAL: A%COL:'
WRITE(MSG%LU_DEBUG,'(9I6)') A%COL
WRITE(MSG%LU_DEBUG,*) 'FINAL: A%COGL:'
WRITE(MSG%LU_DEBUG,'(9I6)') A%COLG
WRITE(MSG%LU_DEBUG,*) 'FINAL: A%VAL:'
WRITE(MSG%LU_DEBUG,'(6E12.4)') A%VAL
#endif
            ENDDO

         ENDDO FACE_LOOP


         !=============================
         !=============================
         !=============================
         ELSE IF (NMATRIX == NSCARC_MATRIX_POISSON .AND. (IS_AMG.OR.IS_CG_SAMG) .AND. NL > NLEVEL_MIN) THEN


         FACE_LOOP2: DO IFACE = 1, 6                          ! check if this face has connection to last cell

            IOR0 = FACE_ORIENTATION(IFACE)
            F => SCARC(NM)%LEVEL(NLEVEL_MIN)%FACE(IOR0)
      
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'HALLO HALLO: IFACE, IOR0:', IFACE, IOR0
#endif
            DO INBR = 1, F%N_NEIGHBORS
      
               NOM = F%NEIGHBORS(INBR)
               CALL SCARC_POINT_TO_OTHER_GRID(NM, NOM, NL)
               OA => SCARC_POINT_TO_OTHER_CMATRIX(OG, NMATRIX)
   
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,'(A,5I6)') 'IFACE, IOR0, INBR, NOM:', IFACE, IOR0, INBR, NOM, G%POISSONC%N_ROW
WRITE(MSG%LU_DEBUG,*) 'G%POISSONC%ROW:'
WRITE(MSG%LU_DEBUG,'( 8I6)') G%POISSONC%ROW(1:G%POISSONC%N_ROW)
WRITE(MSG%LU_DEBUG,*) 'OA%COL'
WRITE(MSG%LU_DEBUG,'( 8I6)') OA%COL
WRITE(MSG%LU_DEBUG,*) 'OA%VAL'
WRITE(MSG%LU_DEBUG,'(8E12.4)') OA%VAL
#endif
               ICOL = 1
               DO ICG = OL%GHOST_FIRSTE(IOR0), OL%GHOST_LASTE(IOR0)
   
                  ICE = OG%ICG_TO_ICE(ICG, 1)
                  IP = A%ROW(ICE) 


#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) '====================== ICG, ICE, IP =', ICG, ICE, IP
WRITE(MSG%LU_DEBUG,*) 'OA%ROW(ICG+1) - OA%ROW(ICG) =', OA%ROW(ICG+1) - OA%ROW(ICG)
WRITE(MSG%LU_DEBUG,*) ' A%ROW(ICE+1) -  A%ROW(ICE) =',  A%ROW(ICE+1) -  A%ROW(ICE)
#endif
                  DO ICOL = OA%ROW(ICG), OA%ROW(ICG+1)-1
                     A%COL(IP) = OA%COL(ICOL)   
                     A%COLG(IP) = ABS(OA%COLG(ICOL))      
                     A%VAL(IP) = OA%VAL(ICOL)
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,'(A,7I6,E12.4)') '       ----> B: ICG, ICE, ICOL, IP, ROW, COL, VAL :', &
                                       ICG, ICE, ICOL, IP, A%ROW(ICE), A%COL(IP), A%COLG(IP), A%VAL(IP)
#endif
                     IP = IP + 1
                  ENDDO
   
               ENDDO
               A%ROW(ICE+1) = IP 
               A%N_ROW = ICE + 1
               A%N_VAL = IP - 1

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'FINAL: ICE+1, ROW, N_ROW, N_VAL:', ICE+1, IP, ICE+1, IP -1
WRITE(MSG%LU_DEBUG,*) 'FINAL: A%N_VAL=', A%N_VAL
WRITE(MSG%LU_DEBUG,*) 'FINAL: A%N_ROW=', A%N_ROW
WRITE(MSG%LU_DEBUG,*) 'FINAL: A%ROW:'
WRITE(MSG%LU_DEBUG,'(9I6)') A%ROW
WRITE(MSG%LU_DEBUG,*) 'FINAL: A%COL:'
WRITE(MSG%LU_DEBUG,'(9I6)') A%COL
WRITE(MSG%LU_DEBUG,*) 'FINAL: A%COGL:'
WRITE(MSG%LU_DEBUG,'(9I6)') A%COLG
WRITE(MSG%LU_DEBUG,*) 'FINAL: A%VAL:'
WRITE(MSG%LU_DEBUG,'(6E12.4)') A%VAL
#endif
            ENDDO

         ENDDO FACE_LOOP2

         !===============================
         !===============================
         !===============================
         ENDIF

      ! ---------- Version for structured grids
      CASE (NSCARC_GRID_UNSTRUCTURED)
         WRITE(*,*) 'SCARC-ERROR: No unstructured version of EXTRACT_MATRIX_OVERLAPS available yet'

   END SELECT
   
#ifdef WITH_SCARC_DEBUG
CALL SCARC_DEBUG_CMATRIX (A, 'A', 'AFTER EXTRACT_MATRIX_OVERLAPS')
#endif
ENDDO
END SUBROUTINE SCARC_EXTRACT_MATRIX_OVERLAPS


! -------------------------------------------------------------------------------------------
! Extract overlapping galerkin part and add to main matrix
! -------------------------------------------------------------------------------------------
SUBROUTINE SCARC_EXTRACT_GALERKIN_OVERLAPS (NMATRIX, NL)
USE SCARC_POINTERS, ONLY : G, F, OL, OG, A
INTEGER, INTENT(IN) :: NL, NMATRIX
INTEGER :: NM, IFACE, NOM, IOR0, ICG, ICE, IP, ICOL, INBR

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G

   SELECT CASE (TYPE_GRID)

      ! ---------- Version for structured grids
      CASE (NSCARC_GRID_STRUCTURED)

         A => SCARC_POINT_TO_CMATRIX(G, NMATRIX)
 
         IP = A%ROW(G%NC+1)
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'EXTRACT_MATRIX_OVERLAPS: NM=', NM,': IP=',IP, G%POISSONC%N_ROW
#endif
      
         FACE_LOOP2: DO IFACE = 1, 6                          ! check if this face has connection to last cell

            IOR0 = FACE_ORIENTATION(IFACE)
            F => SCARC(NM)%LEVEL(NLEVEL_MIN)%FACE(IOR0)
      
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'HALLO HALLO: IFACE, IOR0:', IFACE, IOR0
#endif
            DO INBR = 1, F%N_NEIGHBORS
      
               NOM = F%NEIGHBORS(INBR)
               CALL SCARC_POINT_TO_OTHER_GRID(NM, NOM, NL)
               OA => SCARC_POINT_TO_OTHER_CMATRIX(OG, NMATRIX)
   
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,'(A,5I6)') 'IFACE, IOR0, INBR, NOM:', IFACE, IOR0, INBR, NOM, G%POISSONC%N_ROW
WRITE(MSG%LU_DEBUG,*) 'G%POISSONC%ROW:'
WRITE(MSG%LU_DEBUG,'( 8I6)') G%POISSONC%ROW(1:G%POISSONC%N_ROW)
WRITE(MSG%LU_DEBUG,*) 'OA%COL'
WRITE(MSG%LU_DEBUG,'( 8I6)') OA%COL
WRITE(MSG%LU_DEBUG,*) 'OA%VAL'
WRITE(MSG%LU_DEBUG,'(8E12.4)') OA%VAL
#endif
               ICOL = 1
               DO ICG = OL%GHOST_FIRSTE(IOR0), OL%GHOST_LASTE(IOR0)
   
                  ICE = OG%ICG_TO_ICE(ICG, 1)
                  IP = A%ROW(ICE) 


#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) '====================== ICG, ICE, IP =', ICG, ICE, IP
WRITE(MSG%LU_DEBUG,*) 'OA%ROW(ICG+1) - OA%ROW(ICG) =', OA%ROW(ICG+1) - OA%ROW(ICG)
WRITE(MSG%LU_DEBUG,*) ' A%ROW(ICE+1) -  A%ROW(ICE) =',  A%ROW(ICE+1) -  A%ROW(ICE)
#endif
                  DO ICOL = OA%ROW(ICG), OA%ROW(ICG+1)-1
                     A%COL(IP) = OA%COL(ICOL)   
                     A%COLG(IP) = ABS(OA%COLG(ICOL))      
                     A%VAL(IP) = OA%VAL(ICOL)
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,'(A,7I6,E12.4)') '       ----> B: ICG, ICE, ICOL, IP, ROW, COL, VAL :', &
                                       ICG, ICE, ICOL, IP, A%ROW(ICE), A%COL(IP), A%COLG(IP), A%VAL(IP)
#endif
                     IP = IP + 1
                  ENDDO
   
               ENDDO

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'FINAL: ICE+1, ROW, N_ROW, N_VAL:', ICE+1, IP, ICE+1, IP -1
WRITE(MSG%LU_DEBUG,*) 'FINAL: A%N_VAL=', A%N_VAL
WRITE(MSG%LU_DEBUG,*) 'FINAL: A%N_ROW=', A%N_ROW
WRITE(MSG%LU_DEBUG,*) 'FINAL: A%ROW:'
WRITE(MSG%LU_DEBUG,'(10I6)') A%ROW
WRITE(MSG%LU_DEBUG,*) 'FINAL: A%COL:'
WRITE(MSG%LU_DEBUG,'(10I6)') A%COL
WRITE(MSG%LU_DEBUG,*) 'FINAL: A%COGL:'
WRITE(MSG%LU_DEBUG,'(10I6)') A%COLG
WRITE(MSG%LU_DEBUG,*) 'FINAL: A%VAL:'
WRITE(MSG%LU_DEBUG,'(6E12.4)') A%VAL
#endif
            ENDDO

         ENDDO FACE_LOOP2


      ! ---------- Version for structured grids
      CASE (NSCARC_GRID_UNSTRUCTURED)
         WRITE(*,*) 'SCARC-ERROR: No unstructured version of EXTRACT_MATRIX_OVERLAPS available yet'

   END SELECT
   
#ifdef WITH_SCARC_DEBUG
CALL SCARC_DEBUG_CMATRIX (A, 'A', 'AFTER EXTRACT_MATRIX_OVERLAPS')
#endif
ENDDO
END SUBROUTINE SCARC_EXTRACT_GALERKIN_OVERLAPS


! ------------------------------------------------------------------------------------------------
! Check if cell is within a given mesh
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
! Allocate Poisson matrix for the usual 5-point-stencil (2D) or 7-point-stencil (3D)
! Compact storage technique (POISSONC)
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
USE SCARC_POINTERS, ONLY: L, G, OG, A, AB, OA, OAB
INTEGER, INTENT(IN) :: NM, NL
INTEGER :: IX, IY, IZ, IC, IP, INBR, NOM

!
! Compute single matrix entries and corresponding row and column pointers
! Along internal boundaries use placeholders for the neighboring matrix entries
! which will be communicated in a following step
!
SELECT_STORAGE_TYPE: SELECT CASE (SCARC_MATRIX_LEVEL(NL))

   !
   ! ---------- COMPACT Storage technique
   !
   CASE (NSCARC_MATRIX_COMPACT)
   
      ! Allocate main matrix on non-overlapping part of mesh
      CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G
      A => SCARC_POINT_TO_CMATRIX (G, NSCARC_MATRIX_POISSON)
      CALL SCARC_ALLOCATE_CMATRIX (A, NL, NSCARC_PRECISION_DOUBLE, NSCARC_MATRIX_FULL, 'G%POISSONC')

      ! For every neighbor allocate small matrix on overlapping part of mesh
      DO INBR = 1, SCARC(NM)%N_NEIGHBORS
         NOM = S%NEIGHBORS(INBR)
         CALL SCARC_POINT_TO_OTHER_GRID(NM, NOM, NL)
         OA => SCARC_POINT_TO_OTHER_CMATRIX (OG, NSCARC_MATRIX_POISSON)
         CALL SCARC_ALLOCATE_CMATRIX(OA, NL, NSCARC_PRECISION_DOUBLE, NSCARC_MATRIX_FULL, 'OG%POISSONC')
      ENDDO

      IP = 1
      DO IZ = 1, L%NZ
         DO IY = 1, L%NY
            DO IX = 1, L%NX
   
               IF (IS_UNSTRUCTURED .AND. L%IS_SOLID(IX, IY, IZ)) CYCLE
               IC = G%CELL_NUMBER(IX, IY, IZ)

               ! Main diagonal 
               CALL SCARC_SETUP_POISSONC_MAINDIAG (IC, IX, IY, IZ, IP)
   
               ! Lower subdiagonals
               IF (VALID_SUBDIAG(IX, IY, IZ,  3)) CALL SCARC_SETUP_POISSONC_SUBDIAG(IC, IX, IY, IZ, IX  , IY  , IZ-1, IP,  3)
               IF (VALID_SUBDIAG(IX, IY, IZ,  2)) CALL SCARC_SETUP_POISSONC_SUBDIAG(IC, IX, IY, IZ, IX  , IY-1, IZ  , IP,  2)
               IF (VALID_SUBDIAG(IX, IY, IZ,  1)) CALL SCARC_SETUP_POISSONC_SUBDIAG(IC, IX, IY, IZ, IX-1, IY  , IZ  , IP,  1)
   
               ! Upper subdiagonals
               IF (VALID_SUBDIAG(IX, IY, IZ, -1)) CALL SCARC_SETUP_POISSONC_SUBDIAG(IC, IX, IY, IZ, IX+1, IY  , IZ  , IP, -1)
               IF (VALID_SUBDIAG(IX, IY, IZ, -2)) CALL SCARC_SETUP_POISSONC_SUBDIAG(IC, IX, IY, IZ, IX  , IY+1, IZ  , IP, -2)
               IF (VALID_SUBDIAG(IX, IY, IZ, -3)) CALL SCARC_SETUP_POISSONC_SUBDIAG(IC, IX, IY, IZ, IX  , IY  , IZ+1, IP, -3)
   
            ENDDO
         ENDDO
      ENDDO
   
      A%ROW(G%NC+1) = IP
      A%N_VAL = IP
   
      CALL SCARC_GET_MATRIX_STENCIL_MAX(A, G%NC)

   
   !
   ! ---------- bandwise storage technique
   !
   CASE (NSCARC_MATRIX_BANDWISE)
   
      ! Allocate main matrix on non-overlapping part of mesh
      CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G
      AB => SCARC_POINT_TO_BMATRIX (G, NSCARC_MATRIX_POISSON)
      CALL SCARC_ALLOCATE_BMATRIX(AB, NL, 'G%POISSONB')
   
      ! For every neighbor allocate little matrix on overlapping part of mesh
      DO INBR = 1, SCARC(NM)%N_NEIGHBORS
         NOM = S%NEIGHBORS(INBR)
         CALL SCARC_POINT_TO_OTHER_GRID(NM, NOM, NL)
         OAB => SCARC_POINT_TO_BMATRIX (G, NSCARC_MATRIX_POISSON)
         CALL SCARC_ALLOCATE_BMATRIX(OAB, NL, 'OG%POISSONB')
      ENDDO
   
      IP  = 1
      DO IZ = 1, L%NZ
         DO IY = 1, L%NY
            DO IX = 1, L%NX
   
               IF (IS_UNSTRUCTURED .AND. L%IS_SOLID(IX, IY, IZ)) CYCLE
               IC = G%CELL_NUMBER(IX, IY, IZ)
   
               ! Lower subdiagonals
               IF (VALID_SUBDIAG(IX, IY, IZ,  3)) CALL SCARC_SETUP_POISSONB_SUBDIAG(IC, IX, IY, IZ, IX  , IY  , IZ-1,  3)
               IF (VALID_SUBDIAG(IX, IY, IZ,  2)) CALL SCARC_SETUP_POISSONB_SUBDIAG(IC, IX, IY, IZ, IX  , IY-1, IZ  ,  2)
               IF (VALID_SUBDIAG(IX, IY, IZ,  1)) CALL SCARC_SETUP_POISSONB_SUBDIAG(IC, IX, IY, IZ, IX-1, IY  , IZ  ,  1)
   
               ! Main diagonal
               CALL SCARC_SETUP_POISSONB_MAINDIAG (IC, IX, IY, IZ)

               ! Upper subdiagonals
               IF (VALID_SUBDIAG(IX, IY, IZ, -1)) CALL SCARC_SETUP_POISSONB_SUBDIAG(IC, IX, IY, IZ, IX+1, IY  , IZ  , -1)
               IF (VALID_SUBDIAG(IX, IY, IZ, -2)) CALL SCARC_SETUP_POISSONB_SUBDIAG(IC, IX, IY, IZ, IX  , IY+1, IZ  , -2)
               IF (VALID_SUBDIAG(IX, IY, IZ, -3)) CALL SCARC_SETUP_POISSONB_SUBDIAG(IC, IX, IY, IZ, IX  , IY  , IZ+1, -3)
   
            ENDDO
         ENDDO
      ENDDO
   
END SELECT SELECT_STORAGE_TYPE

END SUBROUTINE SCARC_SETUP_POISSON


! ------------------------------------------------------------------------------------------------
! Check if a subdiagonal entry must be computed in direction IOR0
! If a structured discretization is used, then subdiagonals are built in every direction
! Else check the type of the neighboring cell in direction IOR0
! ------------------------------------------------------------------------------------------------
LOGICAL FUNCTION VALID_SUBDIAG(IX, IY, IZ, IOR0)
USE SCARC_POINTERS, ONLY: L, G
INTEGER, INTENT(IN)  :: IX, IY, IZ, IOR0
INTEGER :: IC_INDEX, IW_INDEX

VALID_SUBDIAG = .FALSE.
IF (TWO_D .AND. ABS(IOR0) == 2) RETURN

SELECT CASE (TYPE_GRID)
   CASE (NSCARC_GRID_STRUCTURED)
      VALID_SUBDIAG = .TRUE.                                                ! always build subdiagonals
      RETURN
   CASE (NSCARC_GRID_UNSTRUCTURED)
      IC_INDEX = L%CELL_INDEX_PTR(IX, IY, IZ)                               ! cell index of corresponding cell
      IW_INDEX = 0
      IF (IC_INDEX /= 0) IW_INDEX  = L%WALL_INDEX_PTR(IC_INDEX, -IOR0)      ! check its wall index

      IF (IW_INDEX == 0) THEN                                               ! if zero, build subdiagonal 
         VALID_SUBDIAG = .TRUE.
         RETURN
      ELSE                                                                  ! if not, only build along interfaces
         IF (G%WALL(IW_INDEX)%BOUNDARY_TYPE== INTERPOLATED_BOUNDARY) THEN
            VALID_SUBDIAG = .TRUE.
            RETURN
         ENDIF
      ENDIF
END SELECT
RETURN

END FUNCTION VALID_SUBDIAG


! ------------------------------------------------------------------------------------------------
! Compact storage technique:
! Set main diagonal entry for matrix - full matrix of the global problem
! In case of an equidistant grid, we get the usual 5-point (2d) and 7-point (3d) stencil
! If two meshes with different step sizes meet, we get a weighted stencil along internal wall cells
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_POISSONC_MAINDIAG (IC, IX, IY, IZ, IP)
USE SCARC_POINTERS, ONLY: L, G, A
INTEGER, INTENT(IN) :: IC, IX, IY, IZ
INTEGER, INTENT(INOUT) :: IP

A => SCARC_POINT_TO_CMATRIX (G, NSCARC_MATRIX_POISSON)

A%VAL(IP) = - 2.0_EB/(L%DXL(IX-1)*L%DXL(IX))
IF (.NOT.TWO_D) A%VAL(IP) = A%VAL(IP) - 2.0_EB/(L%DYL(IY-1)*L%DYL(IY))
A%VAL(IP) = A%VAL(IP) - 2.0_EB/(L%DZL(IZ-1)*L%DZL(IZ))

A%ROW(IC) = IP
A%COL(IP) = IC

A%STENCIL(0) = A%VAL(IP)

IP = IP + 1
END SUBROUTINE SCARC_SETUP_POISSONC_MAINDIAG


! ------------------------------------------------------------------------------------------------
! Compact storage technique:
! Compute subdiagonal contribution in direction IOR0 
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_POISSONC_SUBDIAG (IC, IX1, IY1, IZ1, IX2, IY2, IZ2, IP, IOR0)
USE SCARC_POINTERS, ONLY: L, F, G, A
INTEGER, INTENT(IN) :: IC, IX1, IY1, IZ1, IX2, IY2, IZ2, IOR0
INTEGER, INTENT(INOUT) :: IP
INTEGER :: IW
LOGICAL :: IS_INTERNAL_CELL

A => SCARC_POINT_TO_CMATRIX (G, NSCARC_MATRIX_POISSON)

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

! if IC is an internal cell of the mesh, compute usual matrix contribution for corresponding subdiagonal
IF (IS_INTERNAL_CELL) THEN

   IF (IS_STRUCTURED .OR. .NOT.L%IS_SOLID(IX2, IY2, IZ2)) THEN
      A%VAL(IP) = A%VAL(IP) + F%SCAL_INSIDE
      A%COL(IP) = G%CELL_NUMBER(IX2, IY2, IZ2)
      A%STENCIL(-IOR0) = A%VAL(IP)
      IP = IP + 1
   ELSE
      CALL SCARC_SHUTDOWN(NSCARC_ERROR_MATRIX_SUBDIAG, SCARC_NONE, NSCARC_NONE)
   ENDIF

! if IC is a boundary cell of the mesh, compute matrix contribution only if there is a neighbor for that cell
ELSE IF (L%FACE(IOR0)%N_NEIGHBORS /= 0) THEN

   IW = SCARC_ASSIGN_SUBDIAG_TYPE (IC, IOR0)           ! get IW of a possibly suitable neighbor at face IOR0
   IF (IW > 0) then                                    ! if available, build corresponding subdiagonal entry
      A%VAL(IP) = A%VAL(IP) + F%SCAL_FACE
      A%COL(IP) = G%WALL(IW)%ICE                        ! store its extended number in matrix column pointers
      A%STENCIL(-IOR0) = A%VAL(IP)
      IP = IP + 1
   ENDIF

ENDIF

END SUBROUTINE SCARC_SETUP_POISSONC_SUBDIAG


! ------------------------------------------------------------------------------------------------
! Determine if cell IC has a neighbor and, if yes, return corresponding IW-value
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
! Get type of matrix storage scheme for level NL
! ------------------------------------------------------------------------------------------------
INTEGER FUNCTION SCARC_MATRIX_LEVEL(NL)
INTEGER, INTENT(IN) :: NL

IF (NL == NLEVEL_MAX .AND. TYPE_COARSE == NSCARC_COARSE_DIRECT) THEN
   SCARC_MATRIX_LEVEL = NSCARC_MATRIX_COMPACT
ELSE
   SCARC_MATRIX_LEVEL = TYPE_MATRIX
ENDIF

END FUNCTION SCARC_MATRIX_LEVEL


! ------------------------------------------------------------------------------------------------
! Bandwise storage technique:
! Set main diagonal entry for matrix - full matrix of the global problem
! In case of an equidistant grid, we get the usual 5-point (2d) and 7-point (3d) stencil
! If two meshes with different step sizes meet, we get a weighted stencil along internal wall cells
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_POISSONB_MAINDIAG (IC, IX, IY, IZ)
USE SCARC_POINTERS, ONLY: L, G, AB
INTEGER, INTENT(IN)  :: IC, IX, IY, IZ
INTEGER :: ID

AB => G%POISSONB
ID = AB%POS(0)               ! get column vector corresponding to matrix diagonal

AB%VAL(IC, ID) = AB%VAL(IC, ID) - 2.0_EB/(L%DXL(IX-1)*L%DXL(IX))
IF (.NOT.TWO_D)  AB%VAL(IC, ID) = AB%VAL(IC, ID) - 2.0_EB/(L%DYL(IY-1)*L%DYL(IY))
AB%VAL(IC, ID) = AB%VAL(IC, ID) - 2.0_EB/(L%DZL(IZ-1)*L%DZL(IZ))

AB%STENCIL(0) = AB%VAL(IC, ID)

END SUBROUTINE SCARC_SETUP_POISSONB_MAINDIAG


! ------------------------------------------------------------------------------------------------
! Bandwise storage technique:
! Compute subdiagonal contribution in direction IOR0 
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_POISSONB_SUBDIAG (IC, IX1, IY1, IZ1, IX2, IY2, IZ2, IOR0)
USE SCARC_POINTERS, ONLY: L, F, G
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

!
! if IC is an internal cell of the mesh, compute usual matrix contribution for corresponding subdiagonal
!
IF (IS_INTERNAL_CELL) THEN

   IF (IS_STRUCTURED .OR. .NOT.L%IS_SOLID(IX2, IY2, IZ2)) THEN
      AB%VAL(IC, ID)   = AB%VAL(IC, ID) + F%SCAL_INSIDE
      AB%STENCIL(IOR0) = AB%VAL(IC, ID)
   ELSE
      CALL SCARC_SHUTDOWN(NSCARC_ERROR_MATRIX_SUBDIAG, SCARC_NONE, NSCARC_NONE)
   ENDIF

!
! if IC is a boundary cell of the mesh, compute matrix contribution only if there is a neighbor for that cell
!
!ELSE IF (L%FACE(IOR0)%N_NEIGHBORS /= 0) THEN
ELSE IF (L%FACE(IOR0)%N_NEIGHBORS == 123456) THEN       ! CAUTION: TO FIX AGAIN, ONLY FOR TESTING, IMPOSSIBLE CONDITION

   IW = SCARC_ASSIGN_SUBDIAG_TYPE (IC, IOR0)            ! get IW of a possibly suitable neighbor at face IOR0
   IF (IW /= 0) THEN
      AB%VAL(IC, ID)   = AB%VAL(IC, ID) + F%SCAL_FACE
      AB%STENCIL(IOR0) = AB%VAL(IC, ID)
   ENDIF

ENDIF

END SUBROUTINE SCARC_SETUP_POISSONB_SUBDIAG

! ------------------------------------------------------------------------------------------------
! Get maximum stencil size in matrix NMATRIX (this is known to be 7 on finest level),
! but in SAMG-method this size results only in the course and can be much larger
! (required for dimensioning the coarse-level matrices)
! If NTYPE == 0, only internal matrix part is considered, if NTYPE == 1, also the overlap
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_GET_MATRIX_STENCIL_MAX (A, NLEN)
TYPE (SCARC_MATRIX_COMPACT_TYPE), INTENT(INOUT) :: A
INTEGER, INTENT(IN) :: NLEN
INTEGER :: IC

A%N_STENCIL_MAX = 0
DO IC = 1, NLEN
   A%N_STENCIL_MAX = MAX(A%N_STENCIL_MAX, A%ROW(IC+1)-A%ROW(IC))
ENDDO

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'GET_MATRIX_STENCIL_MAX:', A%N_STENCIL_MAX
#endif
END SUBROUTINE SCARC_GET_MATRIX_STENCIL_MAX


#ifdef WITH_MKL
! ------------------------------------------------------------------------------------------------
! Build system matrix for MKL solver in double precision
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_POISSON_MKL (NM, NL)
USE SCARC_POINTERS, ONLY: G, A, AS
INTEGER, INTENT(IN) :: NM, NL
INTEGER :: IC, JC, JC0, ICS, JCS, JCG
INTEGER :: ICOL, JCOL, IAS
INTEGER :: ISYM, JSYM, NSYM
REAL(EB) :: VAL = 0.0_EB, VALS = 0.0_EB, DIFF
LOGICAL  :: BSYM, BCHECK_SYMMETRY = .FALSE.
INTEGER, DIMENSION(:), ALLOCATABLE :: AUX_ICOL, AUX_IC

CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G
A  => SCARC_POINT_TO_CMATRIX (G, NSCARC_MATRIX_POISSON)
AS => SCARC_POINT_TO_CMATRIX (G, NSCARC_MATRIX_POISSON_SYM)

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'STARTING MATRIX_MKL FOR NM, NL:', NM, NL
WRITE(MSG%LU_DEBUG,*) 'A%N_ROW =', A%N_ROW
WRITE(MSG%LU_DEBUG,*) 'A%N_VAL =', A%N_VAL
#endif

#ifdef WITH_SCARC_DEBUG
CALL SCARC_DEBUG_CMATRIX (A, 'POISSON', 'SETUP_MATRIX_MKL_DOUBLE: BEGIN')
#endif

! 
! ---------- Store only symmetric parts of matrix (diagonal and upper part)
! 
IF (SCARC_MKL_MTYPE == 'SYMMETRIC') THEN

   IF (BCHECK_SYMMETRY) THEN
      ! First check whether symmetry of system matrix is guaranteed
      DO IC = 1, G%NC
         COLUMN_LOOP: DO ICOL = A%ROW(IC)+1, A%ROW(IC+1)-1
            ICS = A%COLG(ICOL)
            VAL = A%VAL(ICOL)
            IF (ICS > IC .AND. ICS <= G%NC) THEN
               BSYM = .FALSE.
               DO JCOL = A%ROW(ICS)+1, A%ROW(ICS+1)-1
                  JCS = A%COLG(JCOL)
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

   !
   ! Compute number of entries in symmetric matrix
   !
   AS%N_VAL = 0
   DO IC = 1, G%NC
      DO ICOL = A%ROW(IC), A%ROW(IC+1)-1
         IF (TYPE_MKL(NL) == NSCARC_MKL_LOCAL) THEN
            JC = A%COLG(ICOL)
            IF (JC >= IC .AND. JC <= G%NC) AS%N_VAL = AS%N_VAL+1
         ELSE IF (TYPE_MKL(NL) == NSCARC_MKL_GLOBAL) THEN
            IF (NL == NLEVEL_MIN) THEN
               JCG = G%LOCAL_TO_GLOBAL(A%COLG(ICOL))
            ELSE
               JCG = A%COLG(ICOL)
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

!
! allocate storage for symmetric matrix and its column and row pointers
! 
CALL SCARC_GET_MATRIX_STENCIL_MAX(A, G%NC)
AS%N_ROW = G%NC+1
AS%N_VAL = A%N_STENCIL_MAX * G%NC
CALL SCARC_ALLOCATE_CMATRIX (AS, NL, NSCARC_PRECISION_DOUBLE, NSCARC_MATRIX_FULL, 'G%AS')

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'A%N_STENCIL_MAX=',A%N_STENCIL_MAX
WRITE(MSG%LU_DEBUG,*) 'AS%N_ROW=',AS%N_ROW
WRITE(MSG%LU_DEBUG,*) 'AS%N_VAL=',AS%N_VAL
#endif

! if global MKL method is used, also allocate auxiliary space for computation of global numbering
IF (IS_MKL_LEVEL(NL)) THEN
   CALL SCARC_ALLOCATE_INT1(AUX_ICOL, 1, A%N_STENCIL_MAX, NSCARC_HUGE_INT, 'AUX_ICOL')
   CALL SCARC_ALLOCATE_INT1(AUX_IC  , 1, A%N_STENCIL_MAX, NSCARC_HUGE_INT, 'AUX_IC')
ENDIF

! 
! extract symmetric matrix part from usual system matrix
! 
IAS = 1
DO IC = 1, AS%N_ROW - 1
   AS%ROW(IC) = IAS

   TYPE_MKL_SELECT: SELECT CASE (TYPE_MKL(NL)) 

      ! blockwise use of local MKL solvers - no global numbering required
      CASE(NSCARC_MKL_LOCAL) 

         DO ICOL = A%ROW(IC), A%ROW(IC+1)-1
            JC = A%COLG(ICOL)
            IF (JC >= IC .AND. JC <= G%NC) THEN
               AS%COLG(IAS) = A%COLG(ICOL)
               SELECT CASE (TYPE_MKL_PRECISION)
                  CASE (NSCARC_PRECISION_DOUBLE)
                     AS%VAL(IAS) = A%VAL(ICOL)
                  CASE (NSCARC_PRECISION_SINGLE)
                     AS%VAL_FB(IAS) = REAL(A%VAL(ICOL),FB)
                  END SELECT
               IAS = IAS + 1
            ENDIF
         ENDDO

      ! global use of MKL solver - get global numbering of matrix elements
      CASE(NSCARC_MKL_GLOBAL) 

         ! store indices of all diagonal and upper-diagonal entries
         AUX_ICOL = 0
         AUX_IC   = NSCARC_HUGE_INT
         ISYM = 1
         JC0 = A%COLG(A%ROW(IC))
         DO ICOL = A%ROW(IC), A%ROW(IC+1)-1
             JC = A%COLG(ICOL)
#ifdef WITH_SCARC_DEBUG2
WRITE(MSG%LU_DEBUG,*) 'IC, ICOL, JC, JC0: ', IC, ICOL, JC, JC0
#endif
            IF (SCARC_MKL_MTYPE == 'SYMMETRIC') THEN
               IF (JC >= JC0) THEN
                  AUX_ICOL(ISYM) = ICOL
                  AUX_IC(ISYM) = JC
                  ISYM  = ISYM  + 1
               ENDIF
            ELSE
               AUX_ICOL(ISYM) = ICOL
               AUX_IC(ISYM) = JC
               ISYM  = ISYM  + 1
            ENDIF
         ENDDO

         NSYM = ISYM - 1
         JSYM = 1

#ifdef WITH_SCARC_DEBUG2
WRITE(MSG%LU_DEBUG,*) '======================================'
WRITE(MSG%LU_DEBUG,*) 'NSYM ', NSYM
WRITE(MSG%LU_DEBUG,*) 'AUX_ICOL:', AUX_ICOL(1:NSYM)
WRITE(MSG%LU_DEBUG,*) 'AUX_IC  :', AUX_IC(1:NSYM)
#endif

         ! sort them in increasing order (for the use of Cluster_Sparse_Solver and PARDISO functionality)
         SORT_LOOP: DO WHILE (JSYM <= NSYM)
            DO ISYM = 1, NSYM
               JC = AUX_IC(ISYM)
               IF (JC == NSCARC_HUGE_INT) CYCLE
               IF (JC <= MINVAL(ABS(AUX_IC(1:NSYM)))) THEN
                  ICOL = AUX_ICOL(ISYM)
                  SELECT CASE (TYPE_MKL_PRECISION)
                     CASE (NSCARC_PRECISION_DOUBLE)
                        AS%VAL(IAS) = A%VAL(ICOL)
                     CASE (NSCARC_PRECISION_SINGLE)
                        AS%VAL_FB(IAS) = REAL(A%VAL(ICOL), FB)
                  END SELECT
                  AS%COLG(IAS) = A%COLG(ICOL)
#ifdef WITH_SCARC_DEBUG2
WRITE(MSG%LU_DEBUG,'(A,7I6, E12.4)') '  IC, JSYM, ISYM, JC, ICOL, IAS, VAL, COLG: ', &
                                        IC, JSYM, ISYM, JC, ICOL, IAS, AS%COLG(IAS), AS%VAL(IAS)
#endif
                  AUX_IC(ISYM) = NSCARC_HUGE_INT            ! mark entry as already used
                  IAS  = IAS  + 1
               ENDIF
            ENDDO
            JSYM = JSYM + 1
         ENDDO SORT_LOOP

   END SELECT TYPE_MKL_SELECT
ENDDO

AS%ROW(AS%N_ROW) = IAS

IF (IS_MKL_LEVEL(NL)) THEN
   DEALLOCATE (AUX_ICOL, STAT=IERROR)
   DEALLOCATE (AUX_IC,   STAT=IERROR)
ENDIF

CALL SCARC_REDUCE_CMATRIX (AS, 'POISSON_SYM')

#ifdef WITH_SCARC_DEBUG
CALL SCARC_DEBUG_CMATRIX(AS, 'POISSON_SYM', 'SETUP_MATRIX_MKL: END')
#endif
END SUBROUTINE SCARC_SETUP_POISSON_MKL
#endif


! ------------------------------------------------------------------------------------------------
! Insert correct boundary conditions into system matrix
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
INTEGER, INTENT(IN) :: NM, NL
INTEGER :: I, J, K, IOR0, IW, IC, NOM, IP, ICO, ICOL

CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G

SELECT CASE (SCARC_MATRIX_LEVEL(NL))

   !
   ! ---------- Matrix in COMPACT storage technique
   !
   CASE (NSCARC_MATRIX_COMPACT)

      A => G%POISSONC

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
                  A%VAL(IP) = A%VAL(IP) - F%SCAL_BOUNDARY
               CASE (NEUMANN)
                  A%VAL(IP) = A%VAL(IP) + F%SCAL_BOUNDARY
            END SELECT

         ! purely Neumann matrix
         ELSE IF (GWC%BTYPE == NEUMANN) THEN
            IP = A%ROW(IC)
            A%VAL(IP) = A%VAL(IP) + F%SCAL_BOUNDARY
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


   !
   ! ---------- Matrix in bandwise storage technique
   !
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
                  AB%VAL(IC, AB%POS(0)) = AB%VAL(IC, AB%POS(0)) - F%SCAL_BOUNDARY
               CASE (NEUMANN)
                  AB%VAL(IC, AB%POS(0)) = AB%VAL(IC, AB%POS(0)) + F%SCAL_BOUNDARY
            END SELECT

         ! Purely Neumann matrix
         ELSE
            IF (GWC%BTYPE == NEUMANN) AB%VAL(IC, AB%POS(0)) = AB%VAL(IC, AB%POS(0)) + F%SCAL_BOUNDARY
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
! Compact matrix storage technique:
! Define switch entries for toggle between original and condensed values
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_CMATRIX_CONDENSED (NM)
USE SCARC_POINTERS, ONLY: L, G, A, ACO
INTEGER, INTENT(IN) :: NM
INTEGER :: ICO = 0, NC, NOM, IP, IC, JC, ICE, ICN, ICOL, IOR0, IW, I, J, K

A => G%POISSONC
LAST_CELL_IN_LAST_MESH_IF: IF (NM == NMESHES) THEN

   NC = G%NC_LOCAL(NMESHES)
   IP = A%ROW(NC)

   ! store column indices and values of diagonal and all off-diagonal entries in last row
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

   !
   ! within last mesh: check which other cells have a connection to the last cell;
   ! in each corresponding matrix row store the column index and value of just that matrix entry
   ! for each direction only one value has to be stored
   !
   JC = NC - 1
   DO IP = A%ROW(JC)+1, A%ROW(JC+1)-1
      IF (A%COL(IP) == NC) THEN
         ICO = ICO + 1
         ACO => A%CONDENSED(ICO)
         ACO%PTR(1)  = IP
         ACO%COL(1)  = JC
         ACO%VAL1(1) = A%VAL(IP)                  ! store original value of system matrix
         ACO%VAL2(1) = 0.0_EB                      ! store new value of condensed system matrix
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
         ACO%VAL1(1) = A%VAL(IP)                  ! store original value of system matrix
         ACO%VAL2(1) = 0.0_EB                      ! store new value of condensed system matrix
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
            ACO%VAL2(1) = 0.0_EB                      ! store new value of condensed system matrix
            ACO%N_COL   = 1
            EXIT
         ENDIF
      ENDDO
   ENDIF

ENDIF LAST_CELL_IN_LAST_MESH_IF

!
! cycle boundary cells to check if there is a periodic communication partner whose stencil is coupled
! with the last cell of last mesh;
! this can be a cell on the opposite side of the own mesh or on a different mesh
! if such a cell exists, store corresponding matrix entry
!
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
! BANDWISE matrix storage technique:
! Define switch entries for toggle between original and condensed values
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_BMATRIX_CONDENSED (NM)
USE SCARC_POINTERS, ONLY: L, G, AB, ABCO
INTEGER, INTENT(IN) :: NM
INTEGER :: ICO = 0, NC, NOM, IOR0, IC, JC, ICE, ICN, IW, I, J, K

AB => G%POISSONB
LAST_CELL_IN_LAST_MESH_BANDWISE_IF: IF (NM == NMESHES) THEN

   NC = G%NC_LOCAL(NMESHES)

   ! store column indices and values of diagonal and all off-diagonal entries in last row
   ! index '1' corresponds to main diagonal entry
   ICO = ICO + 1
   ABCO => AB%CONDENSED(ICO)

   ABCO%IOR0 = 0
   ABCO%ICO  = NC
   ABCO%VAL1(1:AB%N_STENCIL) = AB%VAL(NC, 1:AB%N_STENCIL)
   ABCO%VAL2(1:AB%N_STENCIL) = 0.0_EB
   ABCO%VAL2(AB%POS(0)) = 1.0_EB

   !
   ! within last mesh: check which other cells have a connection to the last cell;
   ! in each corresponding matrix row store the column index and value of just that matrix entry
   ! for each direction only one value has to be stored
   !
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

!
! cycle boundary cells to check if there is a periodic communication partner whose stencil is coupled
! with the last cell of last mesh;
! this can be a cell on the opposite side of the own mesh or a cell on a different mesh
! if such a cell exists, store corresponding matrix entry
!
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
! Setup condensed system in case of periodic or pure Neumann boundary conditions
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_SYSTEM_CONDENSED (NV, NL, ITYPE)
USE SCARC_POINTERS, ONLY: G, F, OL, VC, A, ACO, AB, ABCO
INTEGER, INTENT(IN) :: NV, NL, ITYPE
INTEGER :: NM, NOM, IFACE, ICN, ICE, ICW, JC, NC, ICO, IOR0, IP, ICG, INBR

IF (N_DIRIC_GLOBAL(NLEVEL_MIN) > 0 .OR. &
    TYPE_PRECON == NSCARC_RELAX_FFT .OR. TYPE_PRECON == NSCARC_RELAX_FFTO) RETURN

!
! In last mesh:
! Subtract B*RHS(end) for internal legs of stencil
!
MESH_REAL = 0.0_EB
IF (UPPER_MESH_INDEX == NMESHES) THEN

   CALL SCARC_POINT_TO_GRID (NMESHES, NL)

   NC =  G%NC_LOCAL(NMESHES)
   VC => SCARC_POINT_TO_VECTOR(NMESHES, NL, NV)

   !
   ! process last column entries of all rows except of last one
   ! for those rows only one matrix entry was stored, namely that one which connects to the last cell
   !
   SELECT CASE (SCARC_MATRIX_LEVEL(NL))

      CASE (NSCARC_MATRIX_COMPACT)
         A => G%POISSONC
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

!
! broadcast last RHS-value of last cell in last mesh to all meshes
!
IF (N_MPI_PROCESSES > 1) &
   CALL MPI_ALLGATHER(MPI_IN_PLACE, 1, MPI_DOUBLE_PRECISION, MESH_REAL, 1, MPI_DOUBLE_PRECISION,&
                      MPI_COMM_WORLD, IERROR)

DO NM = 1, NMESHES
   SCARC(NM)%RHS_END = MESH_REAL(NMESHES)
ENDDO

!
! Only in case of periodic BC's:
! Subtract B*RHS(end) for corresponding entries of all periodic communication partners
!
DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G

   SNODE = PROCESS(NM)
   RNODE = PROCESS(NMESHES)

   IF (.NOT. ARE_NEIGHBORS(NM, NMESHES)) CYCLE

   CALL SCARC_POINT_TO_OTHER_GRID(NM, NMESHES, NL)
   VC => SCARC_POINT_TO_VECTOR (NM, NL, NV)

   !
   ! subtract B*RHS(end) at corresponding positions
   !
   DO IFACE = 1, 6                                         ! check if this face has connection to last cell

      IOR0 = FACE_ORIENTATION(IFACE)
      F => L%FACE(IOR0)

      DO INBR = 1, F%N_NEIGHBORS

         NOM = F%NEIGHBORS(INBR)
         IF (NOM /= NMESHES) CYCLE                      ! only check for common matrix entries with last mesh
         CALL SCARC_POINT_TO_OTHER_GRID(NM, NOM, NL)

         DO ICG = OL%GHOST_FIRSTW(IOR0), OL%GHOST_LASTW(IOR0)

            ICW = OG%ICG_TO_ICW(ICG, 1)
            ICE = OG%ICG_TO_ICE(ICG, 1)
            ICN = G%ICE_TO_ICN(ICE)                        ! get column index of neighboring offdiagonal matrix entry

            IF (ICN /= SCARC(NMESHES)%NC) CYCLE            ! if no relation to last cell in last mesh, cycle

            VC(ICW) = VC(ICW) - F%SCAL_FACE * SCARC(NM)%RHS_END

         ENDDO

      ENDDO
   ENDDO
ENDDO

END SUBROUTINE SCARC_SETUP_SYSTEM_CONDENSED


! ----------------------------------------------------------------------------------------------------
! Initialize global solver methods (CG/BICG/GMG/AMG) and corresponding level structures
! ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_METHODS
INTEGER :: NSTACK

SELECT_METHOD: SELECT CASE(TYPE_METHOD)

   !
   ! ------------------ Global Krylov method -------------------------------------
   !
   CASE (NSCARC_METHOD_KRYLOV)

      !
      ! Setup basic CG solver
      !
      NSTACK = NSCARC_STACK_ROOT
      STACK(NSTACK)%SOLVER => MAIN_CG
      CALL SCARC_SETUP_KRYLOV(NSCARC_SOLVER_MAIN, NSCARC_SCOPE_GLOBAL, NSCARC_STAGE_ONE, NSTACK, NLEVEL_MIN, NLEVEL_MIN)

      !
      ! Setup preconditioner for Krylov solver
      !
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

         !
         ! Preconditioning by Geometric multigrid,
         ! either locally or Globally acting, depending on user specification stored in TYPE_SCOPE(1)
         !
         CASE (NSCARC_RELAX_GMG)

            STACK(NSTACK)%SOLVER => PRECON_GMG
            CALL SCARC_SETUP_PRECON(NSTACK, TYPE_SCOPE(1))
            CALL SCARC_SETUP_MULTIGRID(NSCARC_SOLVER_PRECON, TYPE_SCOPE(1), NSCARC_STAGE_TWO, NSTACK, &
                                       NLEVEL_MIN, NLEVEL_MAX)

            NSTACK = NSTACK + 1
            SELECT CASE (TYPE_SMOOTH)

               ! Jacobi-smoothing (acting locally by default)
               CASE (NSCARC_RELAX_JAC)
                  STACK(NSTACK)%SOLVER => SMOOTH_JACOBI
                  CALL SCARC_SETUP_SMOOTH(NSTACK, NSCARC_SCOPE_LOCAL)

               ! SSOR-smoothing (acting locally by default)
               CASE (NSCARC_RELAX_SSOR)
                  STACK(NSTACK)%SOLVER => SMOOTH_SSOR
                  CALL SCARC_SETUP_SMOOTH(NSTACK, NSCARC_SCOPE_LOCAL)

               ! JACOBI-preconditioning in matrix form (acting locally by default)
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

      !
      ! If two-level Krylov, allocate intermediate structures for interpolation and workspace for global coarse solver
      !
      IF (HAS_TWO_LEVELS) THEN

         IF (.NOT.IS_CG_SAMG) CALL SCARC_SETUP_INTERPOLATION(NSCARC_STAGE_ONE, NLEVEL_MIN+1, NLEVEL_MAX)

         NSTACK = NSTACK + 1
         CALL SCARC_SETUP_COARSE_SOLVER(NSCARC_STAGE_ONE, NSCARC_SCOPE_GLOBAL, NSTACK, NLEVEL_MAX, NLEVEL_MAX)

      ENDIF

    !
    ! ------------------ Global Multigrid method -------------------------------------
    !
    CASE (NSCARC_METHOD_MULTIGRID)

       NSTACK = NSCARC_STACK_ROOT
       STACK(NSTACK)%SOLVER => MAIN_GMG
       CALL SCARC_SETUP_MULTIGRID(NSCARC_SOLVER_MAIN, NSCARC_SCOPE_GLOBAL, NSCARC_STAGE_ONE, NSTACK, NLEVEL_MIN, NLEVEL_MAX)

       NSTACK = NSTACK + 1
       SELECT CASE(TYPE_SMOOTH)

          ! Jacobi-smoothing (acting locally by default)
          CASE (NSCARC_RELAX_JAC)
             STACK(NSTACK)%SOLVER => SMOOTH_JACOBI
             CALL SCARC_SETUP_SMOOTH(NSTACK, NSCARC_SCOPE_LOCAL)

          ! SSOR-smoothing (acting locally by default)
          CASE (NSCARC_RELAX_SSOR)
             STACK(NSTACK)%SOLVER => SMOOTH_SSOR
             CALL SCARC_SETUP_SMOOTH(NSTACK, NSCARC_SCOPE_LOCAL)

          ! JACOBI-preconditioning in matrix form (acting locally by default)
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
          ! smoothing by LU-decomposition
          CASE (NSCARC_RELAX_MKL)
             CALL SCARC_SETUP_SMOOTH(NSTACK, TYPE_SCOPE(2))

             SELECT CASE(TYPE_SCOPE(2))

                ! Globally acting - call global CLUSTER_SPARSE_SOLVER on MKL
                CASE (NSCARC_SCOPE_GLOBAL)
                   !WRITE(MSG%LU_VERBOSE,*) 'B: SELECT GLOBAL SMOOTHING'
                   STACK(NSTACK)%SOLVER => SMOOTH_MKL
                   CALL SCARC_SETUP_CLUSTER(NLEVEL_MIN, NLEVEL_MIN)

                ! locally acting - call local PARDISO solvers based on MKL
                CASE (NSCARC_SCOPE_LOCAL)
                   !WRITE(MSG%LU_VERBOSE,*) 'B: SELECT LOCAL SMOOTHING'
                   STACK(NSTACK)%SOLVER => SMOOTH_MKL
                   CALL SCARC_SETUP_PARDISO(NLEVEL_MIN, NLEVEL_MIN)

             END SELECT
#endif

       END SELECT

       ! Globally acting coarse grid solver
       NSTACK = NSTACK + 1
       CALL SCARC_SETUP_COARSE_SOLVER(NSCARC_STAGE_ONE, NSCARC_SCOPE_GLOBAL, NSTACK, NLEVEL_MAX, NLEVEL_MAX)

#ifdef WITH_SCARC_MGM
   !
   ! ---------------- McKenny-Greengard-Mayo method (MGM) --------------------
   !
   CASE (NSCARC_METHOD_MGM)

      ! Allocate velocity vectors along internal obstructions for the setting of internal BC's
      CALL SCARC_SETUP_MGM(NLEVEL_MIN, NLEVEL_MIN)

      ! ------- First part of method: Setup CG solver for inhomogeneous problem on structured discretization
      TYPE_GRID = NSCARC_GRID_STRUCTURED

      NSTACK = NSCARC_STACK_ROOT
      STACK(NSTACK)%SOLVER => MAIN_CG_STRUCTURED
      CALL SCARC_SETUP_KRYLOV(NSCARC_SOLVER_MAIN, NSCARC_SCOPE_GLOBAL, NSCARC_STAGE_ONE, NSTACK, NLEVEL_MIN, NLEVEL_MIN)

      ! use FFT-preconditioning
      NSTACK = NSTACK + 1
      STACK(NSTACK)%SOLVER => PRECON_FFT
      CALL SCARC_SETUP_PRECON(NSTACK, NSCARC_SCOPE_LOCAL)
      CALL SCARC_SETUP_FFT(NLEVEL_MIN, NLEVEL_MIN)

      ! ------- Second part of method: Setup CG solver for homogeneous problem on unstructured discretization
      TYPE_GRID = NSCARC_GRID_UNSTRUCTURED

      NSTACK = NSTACK + 1
      STACK(NSTACK)%SOLVER => MAIN_CG_UNSTRUCTURED
      CALL SCARC_SETUP_KRYLOV(NSCARC_SOLVER_MAIN, NSCARC_SCOPE_GLOBAL, NSCARC_STAGE_ONE, NSTACK, NLEVEL_MIN, NLEVEL_MIN)

      ! for a first proof of concept only use SSOR-preconditioning (may be extended later to other preconditioners)
      NSTACK = NSTACK + 1
      STACK(NSTACK)%SOLVER => PRECON_ILU
      CALL SCARC_SETUP_PRECON(NSTACK, NSCARC_SCOPE_LOCAL)

#endif

   !
   ! ------------------ MKL method -------------------------------------
   !
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
! Setup references to solution vectors related to used scope (main/relax)
! ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_REFERENCES(BX, BB, BD, BR, BV, BY, BZ, NSTACK)
USE SCARC_POINTERS, ONLY: SV
LOGICAL, INTENT(IN) :: BX, BB, BD, BR, BV, BY, BZ
INTEGER, INTENT(IN) :: NSTACK

SV  => STACK(NSTACK)%SOLVER

SELECT CASE (SV%TYPE_STAGE)

   !
   ! Solver from working stage ONE, e.g. main Krylov or MG solver
   !
   CASE (NSCARC_STAGE_ONE)
      IF (BX) SV%X = NSCARC_VECTOR_ONE_X
      IF (BB) SV%B = NSCARC_VECTOR_ONE_B
      IF (BD) SV%D = NSCARC_VECTOR_ONE_D
      IF (BR) SV%R = NSCARC_VECTOR_ONE_R
      IF (BV) SV%V = NSCARC_VECTOR_ONE_V
      IF (BY) SV%Y = NSCARC_VECTOR_ONE_Y
      IF (BZ) SV%Z = NSCARC_VECTOR_ONE_Z

#ifdef WITH_MKL
     IF (TYPE_MKL_PRECISION == NSCARC_PRECISION_SINGLE) THEN
        IF (BX) SV%X_FB = NSCARC_VECTOR_ONE_X
        IF (BB) SV%B_FB = NSCARC_VECTOR_ONE_B
        IF (BR) SV%R_FB = NSCARC_VECTOR_ONE_R
        IF (BV) SV%V_FB = NSCARC_VECTOR_ONE_V
     ENDIF
#endif

#ifdef WITH_SCARC_DEBUG
      SV%E = NSCARC_VECTOR_ONE_E
#endif

   !
   ! Solver from working stage TWO, e.g. MG solver as preconditioner
   !
   CASE (NSCARC_STAGE_TWO)
      IF (BX) SV%X = NSCARC_VECTOR_TWO_X
      IF (BB) SV%B = NSCARC_VECTOR_TWO_B
      IF (BD) SV%D = NSCARC_VECTOR_TWO_D
      IF (BR) SV%R = NSCARC_VECTOR_TWO_R
      IF (BV) SV%V = NSCARC_VECTOR_TWO_V
      IF (BY) SV%Y = NSCARC_VECTOR_TWO_Y
      IF (BZ) SV%Z = NSCARC_VECTOR_TWO_Z
#ifdef WITH_SCARC_DEBUG
      SV%E = NSCARC_VECTOR_TWO_E
#endif

#ifdef WITH_MKL
     IF (TYPE_MKL_PRECISION == NSCARC_PRECISION_SINGLE) THEN
        IF (BX) SV%X_FB = NSCARC_VECTOR_TWO_X
        IF (BB) SV%B_FB = NSCARC_VECTOR_TWO_B
        IF (BR) SV%R_FB = NSCARC_VECTOR_TWO_R
        IF (BV) SV%V_FB = NSCARC_VECTOR_TWO_V
     ENDIF
#endif

END SELECT

END SUBROUTINE SCARC_SETUP_REFERENCES


! ----------------------------------------------------------------------------------------------------
! Allocate and initialize vectors for Krylov method
! ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_VECTORS()
USE SCARC_POINTERS, ONLY: G, SV, ST
INTEGER :: NM, NSTACK, NL

DO NSTACK = 1, N_STACK_TOTAL

   SV  => STACK(NSTACK)%SOLVER

   DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
      DO NL = SV%TYPE_LEVEL(1), SV%TYPE_LEVEL(2)

         CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G
         ST => SCARC(NM)%LEVEL(NL)%STAGE(SV%TYPE_STAGE)

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'Allocating arrays for NSTACK, NM, NL, G%NCE:', NSTACK, NM, NL, G%NCE
#endif
         IF (SV%X /= NSCARC_UNDEF_INT) CALL SCARC_ALLOCATE_REAL1(ST%X, 1, G%NCE, NSCARC_INIT_ZERO, 'ST%X')
         IF (SV%B /= NSCARC_UNDEF_INT) CALL SCARC_ALLOCATE_REAL1(ST%B, 1, G%NCE, NSCARC_INIT_ZERO, 'ST%B')
         IF (SV%D /= NSCARC_UNDEF_INT) CALL SCARC_ALLOCATE_REAL1(ST%D, 1, G%NCE, NSCARC_INIT_ZERO, 'ST%D')
         IF (SV%R /= NSCARC_UNDEF_INT) CALL SCARC_ALLOCATE_REAL1(ST%R, 1, G%NCE, NSCARC_INIT_ZERO, 'ST%R')
         IF (SV%V /= NSCARC_UNDEF_INT) CALL SCARC_ALLOCATE_REAL1(ST%V, 1, G%NCE, NSCARC_INIT_ZERO, 'ST%V')
         IF (SV%Y /= NSCARC_UNDEF_INT) CALL SCARC_ALLOCATE_REAL1(ST%Y, 1, G%NCE, NSCARC_INIT_ZERO, 'ST%Y')
         IF (SV%Z /= NSCARC_UNDEF_INT) CALL SCARC_ALLOCATE_REAL1(ST%Z, 1, G%NCE, NSCARC_INIT_ZERO, 'ST%Z')

#ifdef WITH_SCARC_DEBUG
         IF (SV%E /= NSCARC_UNDEF_INT) CALL SCARC_ALLOCATE_REAL1(ST%E, 1, G%NCE, NSCARC_INIT_ZERO, 'ST%E')
#endif

         IF (TYPE_MKL_PRECISION == NSCARC_PRECISION_SINGLE) THEN
            IF (SV%X_FB /= NSCARC_UNDEF_INT) CALL SCARC_ALLOCATE_REAL1_FB(ST%X_FB, 1, G%NCE, NSCARC_INIT_ZERO, 'ST%X_FB')
            IF (SV%B_FB /= NSCARC_UNDEF_INT) CALL SCARC_ALLOCATE_REAL1_FB(ST%B_FB, 1, G%NCE, NSCARC_INIT_ZERO, 'ST%B_FB')
            IF (SV%R_FB /= NSCARC_UNDEF_INT) CALL SCARC_ALLOCATE_REAL1_FB(ST%R_FB, 1, G%NCE, NSCARC_INIT_ZERO, 'ST%R_FB')
            IF (SV%V_FB /= NSCARC_UNDEF_INT) CALL SCARC_ALLOCATE_REAL1_FB(ST%V_FB, 1, G%NCE, NSCARC_INIT_ZERO, 'ST%V_FB')
         ENDIF

      ENDDO
   ENDDO
ENDDO

END SUBROUTINE SCARC_SETUP_VECTORS


! ----------------------------------------------------------------------------------------------------
! Basic setup of environent on specified stack level
! ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_STACK(NSTACK)
USE SCARC_POINTERS, ONLY: SV
INTEGER, INTENT(IN):: NSTACK

SV => STACK(NSTACK)%SOLVER

SV%TYPE_ACCURACY      = TYPE_ACCURACY 
SV%TYPE_COARSE        = TYPE_COARSE
SV%TYPE_COARSENING    = TYPE_COARSENING
SV%TYPE_CYCLING       = TYPE_CYCLING
SV%TYPE_EXCHANGE      = TYPE_EXCHANGE
SV%TYPE_GRID          = TYPE_GRID
SV%TYPE_INTERPOL      = TYPE_INTERPOL
SV%TYPE_KRYLOV        = TYPE_KRYLOV
SV%TYPE_LEVEL         = TYPE_LEVEL 
SV%TYPE_MATRIX        = TYPE_MATRIX
SV%TYPE_METHOD        = TYPE_METHOD
SV%TYPE_MKL           = TYPE_MKL
SV%TYPE_MKL_PRECISION = TYPE_MKL_PRECISION
SV%TYPE_MULTIGRID     = TYPE_MULTIGRID
SV%TYPE_PARENT        = TYPE_PARENT
SV%TYPE_PRECON        = TYPE_PRECON
SV%TYPE_RELAX         = TYPE_RELAX
SV%TYPE_SCOPE         = TYPE_SCOPE
SV%TYPE_SMOOTH        = TYPE_SMOOTH
SV%TYPE_SOLVER        = TYPE_SOLVER
SV%TYPE_STAGE         = TYPE_STAGE
SV%TYPE_STENCIL       = TYPE_STENCIL
SV%TYPE_TWOLEVEL      = TYPE_TWOLEVEL
SV%TYPE_VECTOR        = TYPE_VECTOR

END SUBROUTINE SCARC_SETUP_STACK


! ----------------------------------------------------------------------------------------------------
! Basic setup of environent on specified stack level
! ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_RESTORE_STACK(NSTACK)
USE SCARC_POINTERS, ONLY: SV
INTEGER, INTENT(IN) :: NSTACK

SV => STACK(NSTACK)%SOLVER

TYPE_ACCURACY      = SV%TYPE_ACCURACY 
TYPE_COARSE        = SV%TYPE_COARSE
TYPE_COARSENING    = SV%TYPE_COARSENING
TYPE_CYCLING       = SV%TYPE_CYCLING
TYPE_EXCHANGE      = SV%TYPE_EXCHANGE
TYPE_GRID          = SV%TYPE_GRID
TYPE_INTERPOL      = SV%TYPE_INTERPOL
TYPE_KRYLOV        = SV%TYPE_KRYLOV
TYPE_LEVEL         = SV%TYPE_LEVEL 
TYPE_MATRIX        = SV%TYPE_MATRIX
TYPE_METHOD        = SV%TYPE_METHOD
TYPE_MKL           = SV%TYPE_MKL
TYPE_MKL_PRECISION = SV%TYPE_MKL_PRECISION
TYPE_MULTIGRID     = SV%TYPE_MULTIGRID
TYPE_PARENT        = SV%TYPE_PARENT
TYPE_PRECON        = SV%TYPE_PRECON
TYPE_RELAX         = SV%TYPE_RELAX
TYPE_SCOPE         = SV%TYPE_SCOPE
TYPE_SMOOTH        = SV%TYPE_SMOOTH
TYPE_SOLVER        = SV%TYPE_SOLVER
TYPE_STAGE         = SV%TYPE_STAGE
TYPE_STENCIL       = SV%TYPE_STENCIL
TYPE_TWOLEVEL      = SV%TYPE_TWOLEVEL
TYPE_VECTOR        = SV%TYPE_VECTOR

END SUBROUTINE SCARC_RESTORE_STACK


! ----------------------------------------------------------------------------------------------------
! Allocate and initialize vectors for Krylov method
! ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_KRYLOV(NSOLVER, NSCOPE, NSTAGE, NSTACK, NLMIN, NLMAX)
USE SCARC_POINTERS, ONLY: SV
INTEGER, INTENT(IN) :: NSOLVER, NSCOPE, NSTAGE, NSTACK, NLMIN, NLMAX

!
! Basic setup of stack information and types for Krylov method
!
CALL SCARC_SETUP_STACK(NSTACK)

SV  => STACK(NSTACK)%SOLVER
SV%TYPE_METHOD   = NSCARC_METHOD_KRYLOV
SV%TYPE_SOLVER   = NSOLVER
SV%TYPE_SCOPE(0) = NSCOPE
SV%TYPE_STAGE    = NSTAGE
SV%TYPE_LEVEL(1) = NLMIN
SV%TYPE_LEVEL(2) = NLMAX
SV%TYPE_MATRIX   = TYPE_MATRIX

!
! Preset iteration parameters for Krylov method
!
SELECT CASE(NSOLVER)

   ! -------------- Krylov method is used as main solver
   CASE (NSCARC_SOLVER_MAIN)
   
      SV%CNAME = 'SCARC_MAIN_KRYLOV'
   
      SV%EPS = SCARC_KRYLOV_ACCURACY
      SV%NIT = SCARC_KRYLOV_ITERATIONS
   
      SV%TYPE_RELAX    = TYPE_PRECON
      SV%TYPE_TWOLEVEL = TYPE_TWOLEVEL
   
   ! -------------- Krylov method is used as coarse grid solver solver
   CASE (NSCARC_SOLVER_COARSE)
   
      SV%CNAME = 'SCARC_COARSE_KRYLOV'
   
      SV%EPS = SCARC_COARSE_ACCURACY
      SV%NIT = SCARC_COARSE_ITERATIONS
   
      SV%TYPE_RELAX  = NSCARC_RELAX_SSOR               ! only use SSOR-preconditioning for coarse solver
      SV%TYPE_TWOLEVEL = NSCARC_TWOLEVEL_NONE          ! only use one level for coarse solver
   
   ! -------------- Otherwise: print error message
   CASE DEFAULT
   
      CALL SCARC_SHUTDOWN(NSCARC_ERROR_PARSE_INPUT, SCARC_NONE, NSOLVER)
   
END SELECT

!
! Point to solution vectors (in corresponding scope)
!
CALL SCARC_SETUP_REFERENCES(.TRUE.,.TRUE.,.TRUE.,.TRUE.,.TRUE.,.TRUE.,.TRUE., NSTACK)

END SUBROUTINE SCARC_SETUP_KRYLOV


! ----------------------------------------------------------------------------------------------------
! Allocate and initialize vectors for Geometric Multigrid method
! ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_MULTIGRID(NSOLVER, NSCOPE, NSTAGE, NSTACK, NLMIN, NLMAX)
USE SCARC_POINTERS, ONLY: SV
INTEGER, INTENT(IN) :: NSOLVER, NSCOPE, NSTAGE, NSTACK, NLMIN, NLMAX

!
! Basic setup of stack information and types for multigrid method
!
CALL SCARC_SETUP_STACK(NSTACK)

SV  => STACK(NSTACK)%SOLVER
SV%TYPE_METHOD        = NSCARC_METHOD_MULTIGRID
SV%TYPE_SOLVER        = NSOLVER
SV%TYPE_SCOPE(0)      = NSCOPE
SV%TYPE_STAGE         = NSTAGE
SV%TYPE_LEVEL(0)      = NLMAX-NLMIN+1
SV%TYPE_LEVEL(1)      = NLMIN
SV%TYPE_LEVEL(2)      = NLMAX
SV%TYPE_RELAX         = TYPE_SMOOTH
SV%TYPE_GRID          = TYPE_GRID
SV%TYPE_INTERPOL      = TYPE_INTERPOL
SV%TYPE_ACCURACY      = TYPE_ACCURACY
SV%TYPE_CYCLING       = TYPE_CYCLING
SV%TYPE_MATRIX        = TYPE_MATRIX
SV%TYPE_MKL_PRECISION = TYPE_MKL_PRECISION

SELECT CASE(NSOLVER)

   ! -------------- Multigrid method is used as main solver
   CASE (NSCARC_SOLVER_MAIN)
      SV%CNAME = 'SCARC_MAIN_GMG'

   ! -------------- Multigrid method is only used as preconditioner for global Krylov-method
   CASE (NSCARC_SOLVER_PRECON)
      SV%CNAME = 'SCARC_PRECON_GMG'

   ! -------------- Otherwise: print error message
   CASE DEFAULT
      CALL SCARC_SHUTDOWN(NSCARC_ERROR_PARSE_INPUT, SCARC_NONE, NSOLVER)

END SELECT

!
! Preset iteration parameters for Multigrid method
!
SV%EPS = SCARC_MULTIGRID_ACCURACY
SV%NIT = SCARC_MULTIGRID_ITERATIONS

!
! Point to solution vectors (in corresponding scope)
!
CALL SCARC_SETUP_REFERENCES(.TRUE.,.TRUE.,.FALSE.,.TRUE.,.TRUE.,.FALSE.,.TRUE., NSTACK)

END SUBROUTINE SCARC_SETUP_MULTIGRID


! ----------------------------------------------------------------------------------------------------
! Allocate and initialize vectors for MKL-methods
! ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_COARSE_SOLVER(NSTAGE, NSCOPE, NSTACK, NLMIN, NLMAX)
INTEGER, INTENT(IN)    :: NSCOPE, NSTAGE, NLMIN, NLMAX
INTEGER, INTENT(INOUT) :: NSTACK

SELECT_COARSE: SELECT CASE (TYPE_COARSE)

   ! -------------- CG-method is used as iterative coarse grid solver
   CASE (NSCARC_COARSE_ITERATIVE)

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'SETUP_COARSE_SOLVER: KRYLOV'
#endif
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
      IF (NSCOPE == NSCARC_SCOPE_GLOBAL .AND. NMESHES > 1) THEN
         STACK(NSTACK)%SOLVER => COARSE_CLUSTER
         CALL SCARC_SETUP_MKL(NSCARC_SOLVER_COARSE, NSCOPE, NSTAGE, NSTACK, NLMIN, NLMAX)
         CALL SCARC_SETUP_CLUSTER(NLMIN, NLMAX)

      ! Local scope:
      ! initialize current stack position as PARDISO solver
      ELSE
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


#ifdef WITH_MKL
! ----------------------------------------------------------------------------------------------------
! Allocate and initialize vectors for LU-solvers (based on MKL)
! ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_MKL(NSOLVER, NSCOPE, NSTAGE, NSTACK, NLMIN, NLMAX)
USE SCARC_POINTERS, ONLY: SV
INTEGER, INTENT(IN) :: NSOLVER, NSCOPE, NSTAGE, NSTACK, NLMIN, NLMAX

!
! Basic setup of stack information and types for MKL
!
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

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'SETUP_MKL: OMEGA=', SV%OMEGA
#endif

!
! Reset types for LU-decomposition method
!
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

!
! Point to solution vectors (in corresponding scope)
!
CALL SCARC_SETUP_REFERENCES(.TRUE.,.TRUE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.TRUE., NSTACK)

END SUBROUTINE SCARC_SETUP_MKL
#endif


! ----------------------------------------------------------------------------------------------------
! Allocate and initialize vectors for Krylov method
! ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_PRECON(NSTACK, NSCOPE)
USE SCARC_POINTERS, ONLY: SV, SVP
INTEGER, INTENT(IN) :: NSTACK, NSCOPE

!
! Point to stack entry of called preconditioner and types for preconditioner
!
SV => STACK(NSTACK)%SOLVER

SV%TYPE_SCOPE(0) = NSCOPE

SV%EPS   =  SCARC_PRECON_ACCURACY
SV%NIT   =  SCARC_PRECON_ITERATIONS
SV%OMEGA =  SCARC_PRECON_OMEGA

!
! Preset name of preconditioning method
!
SELECT CASE(TYPE_PRECON)
   CASE (NSCARC_RELAX_JAC)
      SV%CNAME = 'SCARC_PRECON_JAC'
   CASE (NSCARC_RELAX_SSOR)
      SV%CNAME = 'SCARC_PRECON_SSOR'
   CASE (NSCARC_RELAX_MJAC)
      SV%CNAME = 'SCARC_PRECON_MJAC'
   CASE (NSCARC_RELAX_MGS)
      SV%CNAME = 'SCARC_PRECON_MGS'
   CASE (NSCARC_RELAX_MSGS)
      SV%CNAME = 'SCARC_PRECON_MSGS'
   CASE (NSCARC_RELAX_MSOR)
      SV%CNAME = 'SCARC_PRECON_MSOR'
   CASE (NSCARC_RELAX_MSSOR)
      SV%CNAME = 'SCARC_PRECON_MSSOR'
   CASE (NSCARC_RELAX_ILU)
      SV%CNAME = 'SCARC_PRECON_ILU'
      SV%OMEGA = 1.0_EB
   CASE (NSCARC_RELAX_FFT)
      SV%CNAME = 'SCARC_PRECON_FFT'
      SV%OMEGA = 1.0_EB
   CASE (NSCARC_RELAX_FFTO)
      SV%CNAME = 'SCARC_PRECON_FFTO'
      SV%OMEGA = 1.0_EB
   CASE (NSCARC_RELAX_GMG)
      SV%CNAME = 'SCARC_PRECON_GMG'
#ifdef WITH_MKL
   CASE (NSCARC_RELAX_MKL)
      SV%OMEGA = 1.0_EB
      IF (NSCOPE == NSCARC_SCOPE_LOCAL) THEN
         SV%CNAME = 'SCARC_PRECON_PARDISO'
      ELSE
         SV%CNAME = 'SCARC_PRECON_CLUSTER'
         SV%NIT   = 1
      ENDIF
#endif
   CASE DEFAULT
      CALL SCARC_SHUTDOWN(NSCARC_ERROR_PARSE_INPUT, SCARC_NONE, TYPE_PRECON)
END SELECT

!
! Preset types for preconditioner (use same as for calling solver)
!
SVP => STACK(NSTACK-1)%SOLVER

SV%TYPE_ACCURACY      = SVP%TYPE_ACCURACY
SV%TYPE_GRID          = SVP%TYPE_GRID
SV%TYPE_LEVEL(0:2)    = SVP%TYPE_LEVEL(0:2)
SV%TYPE_INTERPOL      = SVP%TYPE_INTERPOL
SV%TYPE_MATRIX        = SVP%TYPE_MATRIX
SV%TYPE_MKL_PRECISION = SVP%TYPE_MKL_PRECISION
SV%TYPE_RELAX         = SVP%TYPE_RELAX
SV%TYPE_SOLVER        = SVP%TYPE_SOLVER
SV%TYPE_STAGE         = SVP%TYPE_STAGE

!
! Preset pointers for preconditioner (use same as for alling solver)
!
SV%X = SVP%X
SV%B = SVP%B
SV%D = SVP%D
SV%R = SVP%R
SV%V = SVP%V
SV%Y = SVP%Y
SV%Z = SVP%Z

IF (TYPE_MKL_PRECISION == NSCARC_PRECISION_SINGLE) THEN
   SV%X_FB = SVP%X_FB
   SV%B_FB = SVP%B_FB
   SV%R_FB = SVP%R_FB
   SV%V_FB = SVP%V_FB
ENDIF

#ifdef WITH_SCARC_DEBUG
SV%E = SVP%E
#endif

END SUBROUTINE SCARC_SETUP_PRECON


! ----------------------------------------------------------------------------------------------------
! Allocate and initialize vectors for Krylov method
! ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_SMOOTH(NSTACK, NSCOPE)
USE SCARC_POINTERS, ONLY: SV, SVP
INTEGER, INTENT(IN) :: NSTACK, NSCOPE

!
! Basic setup of stack information and types/names for smoother
!
CALL SCARC_SETUP_STACK(NSTACK)

SV => STACK(NSTACK)%SOLVER
SV%TYPE_SCOPE(0) = NSCOPE
SV%OMEGA =  SCARC_SMOOTH_OMEGA

SELECT CASE(TYPE_SMOOTH)
   CASE (NSCARC_RELAX_JAC)
      SV%CNAME = 'SCARC_SMOOTH_JACOBI'
   CASE (NSCARC_RELAX_SSOR)
      SV%CNAME = 'SCARC_SMOOTH_SSOR'
   CASE (NSCARC_RELAX_MJAC)
      SV%CNAME = 'SCARC_SMOOTH_MJAC'
   CASE (NSCARC_RELAX_MGS)
      SV%CNAME = 'SCARC_SMOOTH_MGS'
   CASE (NSCARC_RELAX_MSGS)
      SV%CNAME = 'SCARC_SMOOTH_MSGS'
   CASE (NSCARC_RELAX_MSOR)
      SV%CNAME = 'SCARC_SMOOTH_MSOR'
   CASE (NSCARC_RELAX_MSSOR)
      SV%CNAME = 'SCARC_SMOOTH_MSSOR'
   CASE (NSCARC_RELAX_ILU)
      SV%CNAME = 'SCARC_SMOOTH_ILU'
   CASE (NSCARC_RELAX_FFT)
      SV%CNAME = 'SCARC_SMOOTH_FFT'
      SV%OMEGA = 1.0_EB
   CASE (NSCARC_RELAX_FFTO)
      SV%CNAME = 'SCARC_SMOOTH_FFTO'
      SV%OMEGA = 1.0_EB
#ifdef WITH_MKL
   CASE (NSCARC_RELAX_MKL)
      SV%OMEGA = 1.0_EB
      IF (NSCOPE == NSCARC_SCOPE_GLOBAL) THEN
         SV%CNAME = 'SCARC_SMOOTH_CLUSTER'
      ELSE
         SV%CNAME = 'SCARC_SMOOTH_PARDISO'
      ENDIF
#endif
   CASE DEFAULT
      CALL SCARC_SHUTDOWN(NSCARC_ERROR_PARSE_INPUT, SCARC_NONE, TYPE_SMOOTH)
END SELECT

!
! Preset iteration parameters for Multigrid method
!
SV%EPS = SCARC_SMOOTH_ACCURACY
SV%NIT = SCARC_SMOOTH_ITERATIONS

!
! Preset types for preconditioner (use same descriptors as calling solver)
!
SVP => STACK(NSTACK-1)%SOLVER

SV%TYPE_SOLVER        = NSCARC_SOLVER_SMOOTH
SV%TYPE_STAGE         = SVP%TYPE_STAGE
SV%TYPE_LEVEL(0:2)    = SVP%TYPE_LEVEL(0:2)
SV%TYPE_RELAX         = SVP%TYPE_RELAX
SV%TYPE_GRID          = SVP%TYPE_GRID
SV%TYPE_INTERPOL      = SVP%TYPE_INTERPOL
SV%TYPE_ACCURACY      = SVP%TYPE_ACCURACY
SV%TYPE_MKL_PRECISION = SVP%TYPE_MKL_PRECISION
SV%TYPE_MATRIX        = SVP%TYPE_MATRIX

!
!
! Preset references for preconditioner (use same pointers as calling solver)
SV%X = SVP%X
SV%B = SVP%B
SV%D = SVP%D
SV%R = SVP%R
SV%V = SVP%V
SV%Y = SVP%Y
SV%Z = SVP%Z

IF (TYPE_MKL_PRECISION == NSCARC_PRECISION_SINGLE) THEN
   SV%X_FB = SVP%X_FB
   SV%B_FB = SVP%B_FB
   SV%R_FB = SVP%R_FB
   SV%V_FB = SVP%V_FB
ENDIF

#ifdef WITH_SCARC_DEBUG
SV%E = SVP%E
SV%R = SVP%R
#endif

END SUBROUTINE SCARC_SETUP_SMOOTH


! ----------------------------------------------------------------------------------------------------
! Allocate and initialize vectors for additive or multiplicative coarse grid
! (corresponding to Schwarz domain decomposition method)
! ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_INTERPOLATION(NSTAGE, NLMIN, NLMAX)
USE SCARC_POINTERS, ONLY: G, ST
INTEGER, INTENT(IN) :: NSTAGE, NLMIN, NLMAX
INTEGER :: NM, NL

MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   LEVEL_LOOP: DO NL = NLMIN, NLMAX

      CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G
      ST => SCARC(NM)%LEVEL(NL)%STAGE(NSTAGE)

      CALL SCARC_ALLOCATE_REAL1(ST%X, 1, G%NCE, NSCARC_INIT_ZERO, 'ST%X')
      CALL SCARC_ALLOCATE_REAL1(ST%B, 1, G%NCE, NSCARC_INIT_ZERO, 'ST%B')
      CALL SCARC_ALLOCATE_REAL1(ST%V, 1, G%NCE, NSCARC_INIT_ZERO, 'ST%Q')
      CALL SCARC_ALLOCATE_REAL1(ST%R, 1, G%NCE, NSCARC_INIT_ZERO, 'ST%W')
      CALL SCARC_ALLOCATE_REAL1(ST%Y, 1, G%NCE, NSCARC_INIT_ZERO, 'ST%Y')
      CALL SCARC_ALLOCATE_REAL1(ST%Z, 1, G%NCE, NSCARC_INIT_ZERO, 'ST%Z')

   ENDDO LEVEL_LOOP
ENDDO MESHES_LOOP

END SUBROUTINE SCARC_SETUP_INTERPOLATION


! ----------------------------------------------------------------------------------------------------
! Allocate and initialize vectors for blockwise FFT methods
! New here: Perform own initialization of FFT based on H2CZIS/H3CZIS and use own SAVE and WORK arrays
! ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_FFT(NLMIN, NLMAX)
USE SCARC_POINTERS, ONLY: M, S, L, FFT
USE POIS, ONLY: H2CZIS, H3CZIS
INTEGER, INTENT(IN) :: NLMIN, NLMAX
INTEGER :: NM, NL, IERR = 0

!
! Allocate working space for FFT routine
!
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

      CALL SCARC_ALLOCATE_REAL1(FFT%SAVE1, -3, FFT%LSAVE, NSCARC_INIT_ZERO, 'FFT%SAVE1')
      CALL SCARC_ALLOCATE_REAL1(FFT%WORK ,  1, FFT%LWORK, NSCARC_INIT_ZERO, 'FFT%WORK')

      ! Allocate stretching vector (set to 1)
      CALL SCARC_ALLOCATE_REAL1(FFT%HX, 0, FFT%ITRN, NSCARC_INIT_ONE, 'FFT%HX')

      ! Allocate RHS vector for FFT routine
      IF (L%NY == 1) THEN
         CALL SCARC_ALLOCATE_REAL3(FFT%PRHS, 1, FFT%ITRN, 1,        1, 1, FFT%KTRN, NSCARC_INIT_ZERO, 'FFT%PRHS')
      ELSE
         CALL SCARC_ALLOCATE_REAL3(FFT%PRHS, 1, FFT%ITRN, 1, FFT%JTRN, 1, FFT%KTRN, NSCARC_INIT_ZERO, 'FFT%PRHS')
      ENDIF

      ! Allocate boundary data vector for XS
      IF (L%NZ>1) THEN
         IF (L%NY  > 1) CALL SCARC_ALLOCATE_REAL2(FFT%BXS, 1, FFT%JTRN, 1, FFT%KTRN, NSCARC_INIT_ZERO, 'FFT%BXS')
         IF (L%NY == 1) CALL SCARC_ALLOCATE_REAL2(FFT%BXS, 1,        1, 1, FFT%KTRN, NSCARC_INIT_ZERO, 'FFT%BXS')
      ELSE
         CALL SCARC_ALLOCATE_REAL2(FFT%BXS, 1, FFT%JTRN, 1, 1, NSCARC_INIT_ZERO, 'FFT%BXS')
      ENDIF

      ! Allocate boundary data vector for XF
      IF (L%NZ>1) THEN
         IF (L%NY  > 1) CALL SCARC_ALLOCATE_REAL2(FFT%BXF, 1, FFT%JTRN, 1, FFT%KTRN, NSCARC_INIT_ZERO, 'FFT%BXF')
         IF (L%NY == 1) CALL SCARC_ALLOCATE_REAL2(FFT%BXF, 1,        1, 1, FFT%KTRN, NSCARC_INIT_ZERO, 'FFT%BXF')
      ELSE
         CALL SCARC_ALLOCATE_REAL2(FFT%BXF, 1, FFT%JTRN, 1, 1, NSCARC_INIT_ZERO, 'FFT%BXF')
      ENDIF

      ! Allocate boundary data vector for YS
      IF (L%NZ > 1) THEN
         CALL SCARC_ALLOCATE_REAL2(FFT%BYS, 1, FFT%ITRN,1, FFT%KTRN, NSCARC_INIT_ZERO, 'FFT%BYS')
      ELSE
         CALL SCARC_ALLOCATE_REAL2(FFT%BYS, 1, FFT%ITRN,1,        1, NSCARC_INIT_ZERO, 'FFT%BYS')
      ENDIF

      ! Allocate boundary data vector for YF
      IF (L%NZ > 1) THEN
         CALL SCARC_ALLOCATE_REAL2(FFT%BYF, 1, FFT%ITRN,1, FFT%KTRN, NSCARC_INIT_ZERO, 'FFT%BYF')
      ELSE
         CALL SCARC_ALLOCATE_REAL2(FFT%BYF, 1, FFT%ITRN,1,        1, NSCARC_INIT_ZERO, 'FFT%BYF')
      ENDIF

      ! Allocate boundary data vector for ZS
      IF (L%NY  > 1) CALL SCARC_ALLOCATE_REAL2(FFT%BZS, 1, FFT%ITRN, 1, FFT%JTRN, NSCARC_INIT_ZERO, 'FFT%BZS')
      IF (L%NY == 1) CALL SCARC_ALLOCATE_REAL2(FFT%BZS, 1, FFT%ITRN, 1,        1, NSCARC_INIT_ZERO, 'FFT%BZS')

      ! Allocate boundary data vector for ZF
      IF (L%NY  >1)  CALL SCARC_ALLOCATE_REAL2(FFT%BZF, 1, FFT%ITRN, 1, FFT%JTRN, NSCARC_INIT_ZERO, 'FFT%BZF')
      IF (L%NY == 1) CALL SCARC_ALLOCATE_REAL2(FFT%BZF, 1, FFT%ITRN, 1,        1, NSCARC_INIT_ZERO, 'FFT%BZF')

      IF (TWO_D) THEN
         CALL H2CZIS(FFT%XS,FFT%XF,FFT%IBAR,FFT%LBC,&
                     FFT%ZS,FFT%ZF,FFT%KBAR,FFT%NBC,&
                     FFT%HX,FFT%XLM,FFT%ITRN,IERR,FFT%SAVE1)
      ELSE
         CALL H3CZIS(FFT%XS,FFT%XF,FFT%IBAR,FFT%LBC,&
                     FFT%YS,FFT%YF,FFT%JBAR,FFT%MBC,&
                     FFT%ZS,FFT%ZF,FFT%KBAR,FFT%NBC,&
                     FFT%HX,FFT%XLM,FFT%ITRN,FFT%JTRN,IERR,FFT%SAVE1)
      ENDIF

   ENDDO LEVEL_LOOP
ENDDO MESHES_LOOP

END SUBROUTINE SCARC_SETUP_FFT

! ----------------------------------------------------------------------------------------------------
! Allocate and initialize vectors for blockwise FFT methods with overlap
! New here: Perform own initialization of FFT based on H2CZIS/H3CZIS and use own SAVE and WORK arrays
! ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_FFTO(NLMIN, NLMAX)
USE SCARC_POINTERS, ONLY: M, S, L, FFT
USE POIS, ONLY: H2CZIS, H3CZIS
INTEGER, INTENT(IN) :: NLMIN, NLMAX
INTEGER :: NM, NL, IERR = 0

!
! Allocate working space for FFT routine
!
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

      CALL SCARC_ALLOCATE_REAL1(FFT%SAVE1, -3, FFT%LSAVE, NSCARC_INIT_ZERO, 'FFT%SAVE1')
      CALL SCARC_ALLOCATE_REAL1(FFT%WORK ,  1, FFT%LWORK, NSCARC_INIT_ZERO, 'FFT%WORK')

      ! Allocate stretching vector (set to 1)
      CALL SCARC_ALLOCATE_REAL1(FFT%HX, 0, FFT%ITRN, NSCARC_INIT_ONE, 'FFT')

      ! Allocate RHS vector for FFT routine
      IF (L%NY == 1) THEN
         CALL SCARC_ALLOCATE_REAL3(FFT%PRHS, 1, FFT%ITRN, 1,        1, 1, FFT%KTRN, NSCARC_INIT_ZERO, 'FFT%PRHS')
      ELSE
         CALL SCARC_ALLOCATE_REAL3(FFT%PRHS, 1, FFT%ITRN, 1, FFT%JTRN, 1, FFT%KTRN, NSCARC_INIT_ZERO, 'FFT%PRHS')
      ENDIF

      ! Allocate boundary data vector for XS
      IF (L%NZ>1) THEN
         IF (L%NY  > 1) CALL SCARC_ALLOCATE_REAL2(FFT%BXS, 1, FFT%JTRN, 1, FFT%KTRN, NSCARC_INIT_ZERO, 'FFT%BXS')
         IF (L%NY == 1) CALL SCARC_ALLOCATE_REAL2(FFT%BXS, 1,        1, 1, FFT%KTRN, NSCARC_INIT_ZERO, 'FFT%BXS')
      ELSE
         CALL SCARC_ALLOCATE_REAL2(FFT%BXS, 1, FFT%JTRN, 1, 1, NSCARC_INIT_ZERO, 'FFT%BXS')
      ENDIF

      ! Allocate boundary data vector for XF
      IF (L%NZ>1) THEN
         IF (L%NY  > 1) CALL SCARC_ALLOCATE_REAL2(FFT%BXF, 1, FFT%JTRN, 1, FFT%KTRN, NSCARC_INIT_ZERO, 'FFT%BXF')
         IF (L%NY == 1) CALL SCARC_ALLOCATE_REAL2(FFT%BXF, 1,        1, 1, FFT%KTRN, NSCARC_INIT_ZERO, 'FFT%BXF')
      ELSE
         CALL SCARC_ALLOCATE_REAL2(FFT%BXF, 1, FFT%JTRN, 1, 1, NSCARC_INIT_ZERO, 'FFT%BXF')
      ENDIF

      ! Allocate boundary data vector for YS
      IF (L%NZ > 1) THEN
         CALL SCARC_ALLOCATE_REAL2(FFT%BYS, 1, FFT%ITRN, 1, FFT%KTRN, NSCARC_INIT_ZERO, 'FFT%BYS')
      ELSE
         CALL SCARC_ALLOCATE_REAL2(FFT%BYS, 1, FFT%ITRN, 1,        1, NSCARC_INIT_ZERO, 'FFT%BYS')
      ENDIF

      ! Allocate boundary data vector for YF
      IF (L%NZ > 1) THEN
         CALL SCARC_ALLOCATE_REAL2(FFT%BYF, 1, FFT%ITRN,1, FFT%KTRN, NSCARC_INIT_ZERO, 'FFT%BYF')
      ELSE
         CALL SCARC_ALLOCATE_REAL2(FFT%BYF, 1, FFT%ITRN,1,        1, NSCARC_INIT_ZERO, 'FFT%BYF')
      ENDIF

      ! Allocate boundary data vector for ZS
      IF (L%NY  > 1) CALL SCARC_ALLOCATE_REAL2(FFT%BZS, 1, FFT%ITRN, 1, FFT%JTRN, NSCARC_INIT_ZERO, 'FFT%BZS')
      IF (L%NY == 1) CALL SCARC_ALLOCATE_REAL2(FFT%BZS, 1, FFT%ITRN, 1,        1, NSCARC_INIT_ZERO, 'FFT%BZS')

      ! Allocate boundary data vector for ZF
      IF (L%NY  >1)  CALL SCARC_ALLOCATE_REAL2(FFT%BZF, 1, FFT%ITRN, 1, FFT%JTRN, NSCARC_INIT_ZERO, 'FFT%BZF')
      IF (L%NY == 1) CALL SCARC_ALLOCATE_REAL2(FFT%BZF, 1, FFT%ITRN, 1,        1, NSCARC_INIT_ZERO, 'FFT%BZF')

      IF (TWO_D) THEN
         CALL H2CZIS(FFT%XS,FFT%XF,FFT%IBAR,FFT%LBC,&
                     FFT%ZS,FFT%ZF,FFT%KBAR,FFT%NBC,&
                     FFT%HX,FFT%XLM,FFT%ITRN,IERR,FFT%SAVE1)
      ELSE
         CALL H3CZIS(FFT%XS,FFT%XF,FFT%IBAR,FFT%LBC,&
                     FFT%YS,FFT%YF,FFT%JBAR,FFT%MBC,&
                     FFT%ZS,FFT%ZF,FFT%KBAR,FFT%NBC,&
                     FFT%HX,FFT%XLM,FFT%ITRN,FFT%JTRN,IERR,FFT%SAVE1)
      ENDIF

   ENDDO LEVEL_LOOP
ENDDO MESHES_LOOP

END SUBROUTINE SCARC_SETUP_FFTO


! ----------------------------------------------------------------------------------------------------
! Store JACOBI preconditioner in matrix form
! Based on the following splitting of A = D - E - F
! where :   D is the diagonal part
! the MJAC-preconditioner in matrix form is defined
!           M_MJAC = D 
! ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_MJAC(NLMIN, NLMAX)
USE SCARC_POINTERS, ONLY: G, A, AB
INTEGER, INTENT(IN) :: NLMIN, NLMAX
INTEGER :: NM, NL, IC

MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   LEVEL_LOOP: DO NL = NLMIN, NLMAX

      CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G

      SELECT CASE(TYPE_MATRIX)
         CASE (NSCARC_MATRIX_COMPACT)
            A => G%POISSONC
            CALL SCARC_ALLOCATE_REAL1(A%RELAX, 1, G%NC, NSCARC_INIT_ZERO, 'G%POISSONC%RELAX')
            DO IC = 1, G%NC
               A%RELAX(IC) = 1.0_EB/A%VAL(A%ROW(IC))
            ENDDO 
         CASE (NSCARC_MATRIX_BANDWISE)
            AB => G%POISSONB
            CALL SCARC_ALLOCATE_REAL1(AB%RELAXD, 1, AB%N_DIAG,  NSCARC_INIT_ZERO, 'G%POISSONC%RELAXD')
            DO IC = 1, G%NC
               AB%RELAXD(IC) = 1.0_EB/AB%VAL(IC, AB%POS(0))
            ENDDO 
       END SELECT

   ENDDO LEVEL_LOOP
ENDDO MESHES_LOOP

END SUBROUTINE SCARC_SETUP_MJAC


! ----------------------------------------------------------------------------------------------------
! Store GS preconditioner in matrix form
! Based on the following splitting of A = D - E - F
! where :   D is the diagonal part
!          -E is the strictly lower part
!          -F is the strictly upper part
! the SGS-preconditioner in matrix form is defined
!           M_MGS = (D - E) = (I - E D^{-1}) D
! ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_MGS(NLMIN, NLMAX)
USE SCARC_POINTERS, ONLY: G, A, AB
INTEGER, INTENT(IN) :: NLMIN, NLMAX
INTEGER :: NM, NL, IC, JC, IPTR

MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   LEVEL_LOOP: DO NL = NLMIN, NLMAX

      CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G

      SELECT CASE (SCARC_MATRIX_LEVEL(NL))

         !
         ! matrix in COMPACT storage technique
         !
         CASE (NSCARC_MATRIX_COMPACT)

            A => G%POISSONC
            CALL SCARC_ALLOCATE_REAL1(A%RELAX, 1, A%N_VAL, NSCARC_INIT_ZERO, 'G%POISSONC%RELAX')

            DO IC = 1, G%NC
               DO IPTR = A%ROW(IC), A%ROW(IC+1)-1
                  JC = A%COL(IPTR)
                  IF (JC <  IC) A%RELAX(IPTR) = A%VAL(IPTR) / A%VAL(A%ROW(JC))
                  IF (JC == IC) A%RELAX(IPTR) = A%VAL(IPTR)
               ENDDO
            ENDDO 

         !
         ! matrix in bandwise storage technique
         !
         CASE (NSCARC_MATRIX_BANDWISE)

            AB => G%POISSONB
            CALL SCARC_ALLOCATE_REAL2(AB%RELAX, 1, AB%N_DIAG, 1, AB%N_STENCIL, NSCARC_INIT_ZERO, 'G%POISSONB%RELAX')
            WRITE(*,*) 'SCARC_SETUP_MGS: BANDWISE: NOT FINISHED YET'

      END SELECT

   ENDDO LEVEL_LOOP
ENDDO MESHES_LOOP

END SUBROUTINE SCARC_SETUP_MGS


! ----------------------------------------------------------------------------------------------------
! Store symmetric Gauss-Seidel preconditioner in matrix form
! Based on the following splitting of A = D - E - F
! where :   D is the diagonal part
!          -E is the strictly lower part
!          -F is the strictly upper part
! the SGS-preconditioner in matrix form is defined
!           M_MSGS = (D - E) D^{-1} (D - F)  =  (I - E D^{-1}) (D - F)
! ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_MSGS(NLMIN, NLMAX)
USE SCARC_POINTERS, ONLY: G, A, AB
INTEGER, INTENT(IN) :: NLMIN, NLMAX
INTEGER :: NM, NL, IC, JC, IPTR, I, IS, IL, IOR0

MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   LEVEL_LOOP: DO NL = NLMIN, NLMAX

      CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G

      SELECT CASE(TYPE_MATRIX)

         !
         ! matrix in COMPACT storage technique
         !
         CASE (NSCARC_MATRIX_COMPACT)

            A => G%POISSONC
            CALL SCARC_ALLOCATE_REAL1(A%RELAX, 1, A%N_VAL, NSCARC_INIT_ZERO, 'G%POISSONC%RELAX')
            A%RELAX = A%VAL

            DO IC = 1, G%NC
               !  l(i,j) = a(i,j)/a(j,j)
               DO IPTR = A%ROW(IC), A%ROW(IC+1)-1
                  JC = A%COL(IPTR)
                  IF (JC < IC)  A%RELAX(IPTR) = A%VAL(IPTR) / A%VAL(A%ROW(JC))
               ENDDO
            ENDDO 

         !
         ! matrix in bandwise storage technique
         !
         CASE (NSCARC_MATRIX_BANDWISE)

            AB => G%POISSONB
            CALL SCARC_ALLOCATE_REAL2(AB%RELAX, 1, AB%N_DIAG, 1, AB%N_STENCIL, NSCARC_INIT_ZERO, 'G%AB%RELAX')
            AB%RELAX = AB%VAL
            DO IOR0 = 3, 1, -1
               IS = AB%TARGET(IOR0)
               IL = AB%LENGTH(IOR0)
               DO I = 1, AB%LENGTH(IOR0)
                  AB%RELAX(IS+I-1, AB%POS(IOR0)) = AB%VAL(IS+I-1, AB%POS(IOR0))/AB%VAL(I,AB%POS(0))
               ENDDO
            ENDDO

      END SELECT

   ENDDO LEVEL_LOOP
ENDDO MESHES_LOOP
      
END SUBROUTINE SCARC_SETUP_MSGS


! ----------------------------------------------------------------------------------------------------
! Store SOR preconditioner in matrix form
! Based on the following splitting of A = D - E - F
! where :   D is the diagonal part
!          -E is the strictly lower part
!          -F is the strictly upper part
! the SOR-preconditioner in matrix form is defined
!           M_MSOR = (D−ωE) = (I−ωE D^{-1}) D
! ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_MSOR(NLMIN, NLMAX, NSTACK)
USE SCARC_POINTERS, ONLY: G, A, AB
INTEGER, INTENT(IN) :: NLMIN, NLMAX, NSTACK
REAL (EB) :: OMEGA
INTEGER :: NM, NL, IC, JC, IPTR

OMEGA = STACK(NSTACK)%SOLVER%OMEGA

MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   LEVEL_LOOP: DO NL = NLMIN, NLMAX

      CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G

      SELECT CASE(TYPE_MATRIX)

         !
         ! matrix in COMPACT storage technique
         !
         CASE (NSCARC_MATRIX_COMPACT)

            A => G%POISSONC
            CALL SCARC_ALLOCATE_REAL1(A%RELAX, 1, A%N_VAL, NSCARC_INIT_ZERO, 'G%POISSONC%RELAX')
      
            DO IC = 1, G%NC
      
               !DO IPTR = A%ROW(IC), A%ROW(IC+1)-1
               !   JC = A%COL(IPTR)
               !   IF (JC <  IC) A%RELAX(IPTR) = OMEGA * A%VAL(IPTR) / A%VAL(A%ROW(JC))
               !   IF (JC == IC) A%RELAX(IPTR) = A%VAL(IPTR)
               !ENDDO
               DO IPTR = A%ROW(IC), A%ROW(IC+1)-1
                  JC = A%COL(IPTR)
                  IF (JC < IC) THEN
                     A%RELAX(IPTR) = OMEGA * A%VAL(IPTR) / A%VAL(A%ROW(JC))
                  ELSE IF (JC == IC) THEN
                     A%RELAX(IPTR) = A%VAL(IPTR)
                  ELSE IF (JC > IC) THEN
                     A%RELAX(IPTR) = OMEGA * A%VAL(IPTR) 
                  ENDIF
               ENDDO
      
            ENDDO 

         !
         ! matrix in bandwise storage technique
         !
         CASE (NSCARC_MATRIX_BANDWISE)

            AB => G%POISSONB
            CALL SCARC_ALLOCATE_REAL2(AB%RELAX, 1, AB%N_DIAG, 1, AB%N_STENCIL, NSCARC_INIT_ZERO, 'G%POISSONB%RELAX')
            WRITE(*,*) 'SCARC_SETUP_MSOR: BANDWISE: NOT FINISHED YET'
      
      END SELECT

   ENDDO LEVEL_LOOP
ENDDO MESHES_LOOP

END SUBROUTINE SCARC_SETUP_MSOR


! ----------------------------------------------------------------------------------------------------
! Store SSOR preconditioner in matrix form
! Based on the following splitting of A = D - E - F
! where :   D is the diagonal part
!          -E is the strictly lower part
!          -F is the strictly upper part
! the SSOR-preconditioner in matrix form is defined
!           B_SSOR = 1/(omega * (2-omega)) * (D - omega * E) D^{-1} (D - omega * F)
!                  = (I - omega E D^{-1}) * [1/(omega * (2-omega)) * D -  1/(2-omega) * F]
! Defining the triangular matrices
!               L  = I - omega E D^{-1}
!               U  = 1/(omega * (2-omega)) * D -  1/(2-omega) * F
! the SSOR-preconditioning can be thought as the solution of two triangular systems
! Both matrices can be stored as a single matrix that occupies the same amount of storage as A
! where the same row and column pointers can be used as for A (identical pattern)
! Note that the diagonal elements of L are 1 (and are omitted, only the diagonal of U is stored there)
! ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_MSSOR(NLMIN, NLMAX, NSTACK)
USE SCARC_POINTERS, ONLY: G, A, AB
INTEGER, INTENT(IN) :: NLMIN, NLMAX, NSTACK
INTEGER :: NM, NL, IC, JC, IPTR, INCR, IS, IL, IOR0, I
REAL(EB) :: OMEGA, SCAL1, SCAL2

OMEGA = STACK(NSTACK)%SOLVER%OMEGA
SCAL1  = 1.0_EB / (OMEGA * (2.0_EB - OMEGA))
SCAL2  = 1.0_EB / (2.0_EB - OMEGA)

MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   LEVEL_LOOP: DO NL = NLMIN, NLMAX

      CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G

      SELECT CASE(TYPE_MATRIX)

         !
         ! matrix in COMPACT storage technique
         !
         CASE (NSCARC_MATRIX_COMPACT)

            A => G%POISSONC
            CALL SCARC_ALLOCATE_REAL1(A%RELAX, 1, A%N_VAL, NSCARC_INIT_ZERO, 'G%POISSONC%RELAX')
      
            DO IC = 1, G%NC
      
               DO IPTR = A%ROW(IC), A%ROW(IC+1)-1
                  JC = A%COL(IPTR)
                  IF (JC < IC) THEN
                     A%RELAX(IPTR) = OMEGA * A%VAL(IPTR) / A%VAL(A%ROW(JC))
                  ELSE IF (JC == IC) THEN
                     A%RELAX(IPTR) = SCAL1 * A%VAL(IPTR)
                  ELSE IF (JC > IC) THEN
                     A%RELAX(IPTR) = SCAL2 * A%VAL(IPTR)
                  ENDIF
               ENDDO

            ENDDO 

         !
         ! matrix in bandwise storage technique
         !
         CASE (NSCARC_MATRIX_BANDWISE)

            AB => G%POISSONB
            CALL SCARC_ALLOCATE_REAL2(AB%RELAX, 1, AB%N_DIAG, 1, AB%N_STENCIL, NSCARC_INIT_ZERO, 'G%POISSONB%RELAX')

            IF (TWO_D) THEN
               INCR = -2
            ELSE
               INCR = -1
            ENDIF
            
            DO IOR0 = 3, 1, INCR
               IS = AB%TARGET(IOR0)
               IL = AB%LENGTH(IOR0)
               DO I = 1, AB%LENGTH(IOR0)
                  AB%RELAX(IS+I-1, AB%POS(IOR0)) = OMEGA * AB%VAL(IS+I-1, AB%POS(IOR0))/AB%VAL(I,AB%POS(0))
               ENDDO
            ENDDO
            DO I = 1, AB%LENGTH(0)
               AB%RELAX(I, AB%POS(0)) = SCAL1 * AB%VAL(I, AB%POS(0))
            ENDDO
            DO IOR0 = -1, -3, INCR
               IS = AB%TARGET(IOR0)
               IL = AB%LENGTH(IOR0)
               DO I = IS, IS + AB%LENGTH(IOR0) -1
                  AB%RELAX(I, AB%POS(IOR0)) = SCAL2 * AB%VAL(I, AB%POS(IOR0))
               ENDDO
            ENDDO

!               L  = I - omega E D^{-1}
!               U  = 1/(omega * (2-omega)) * D -  1/(2-omega) * F

      END SELECT

   ENDDO LEVEL_LOOP
ENDDO MESHES_LOOP

!STACK(NSTACK)%SOLVER%OMEGA = 1.0_EB

END SUBROUTINE SCARC_SETUP_MSSOR


! ----------------------------------------------------------------------------------------------------
! Allocate and initialize ILU(0) decomposition of Poisson matrix
! L- and U-parts are stored in the same array, diagonal elements of L are supposed to be 1
! Based on Saad-algorithm 10.4 from 'Iterative Methods for Sparse Linear Systems':
!   for i = 2 , ... , n do
!      for k = 1 , ... , i-1 and for (i,k) in NZ(A) do
!         compute a_ik = a_ik / a_kk
!         for j = k+1 , ... , n and for (i,j) in NZ(A) do
!            compute a_ij = a_ij - a_ik a_kj
!         enddo
!      enddo
!   enddo
! ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_ILU(NLMIN, NLMAX)
USE SCARC_POINTERS, ONLY: G, A, AB
INTEGER, INTENT(IN) :: NLMIN, NLMAX
INTEGER :: NM, NL, IC, JC, KC, IPTR, JPTR, KPTR, KPTR0, IOR0, JOR0, KOR0
LOGICAL :: BFOUND


MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   LEVEL_LOOP: DO NL = NLMIN, NLMAX

      CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G

      SELECT CASE(TYPE_MATRIX)

         !
         ! matrix in COMPACT storage technique
         !
         CASE (NSCARC_MATRIX_COMPACT)

            A => G%POISSONC
            CALL SCARC_ALLOCATE_REAL1(A%RELAX, 1, A%N_VAL, NSCARC_INIT_ZERO, 'G%POISSONC%RELAX')
            A%RELAX = A%VAL
      
            CELL_LOOP: DO IC = 2, G%NC
      
               COLUMN_LOOP: DO IPTR = A%ROW(IC), A%ROW(IC+1)-1
      
                  KC = A%COL(IPTR)                        ! get number of neighboring cell
                  IF (KC >= IC) CYCLE                      ! only consider neighbors with lower cell numbers than IC
                  IF (A%RELAX(IPTR) == 0) CYCLE
      
                  KPTR = A%ROW(KC)                        ! get diagonal entry of neighbor
                  A%RELAX(IPTR) = A%RELAX(IPTR)/A%RELAX(KPTR)
      
                  DO JPTR = A%ROW(IC), A%ROW(IC+1)-1
      
                     JC = A%COL(JPTR)
                     IF (JC<=KC) CYCLE                     ! only consider neighbors with higher cell numbers than IC
                     IF (A%RELAX(JPTR) == 0) CYCLE
      
                     KPTR = -1
                     DO KPTR0 = A%ROW(KC), A%ROW(KC+1)-1
                        IF (A%COL(KPTR0) == JC) THEN
                          KPTR = KPTR0
                        ENDIF
                     ENDDO
                     IF (KPTR>0) A%RELAX(JPTR) = A%RELAX(JPTR) - A%RELAX(IPTR) * A%RELAX(KPTR)
      
                  ENDDO
      
               ENDDO COLUMN_LOOP
            ENDDO CELL_LOOP
      
         !
         ! matrix in bandwise storage technique
         !
         CASE (NSCARC_MATRIX_BANDWISE)

            AB => G%POISSONB
            CALL SCARC_ALLOCATE_REAL2(AB%RELAX, 1, AB%N_DIAG, 1, AB%N_STENCIL, NSCARC_INIT_ZERO, 'G%POISSONB%RELAX')
            AB%RELAX = AB%VAL
      
            CELL_BANDWISE_LOOP: DO IC = 2, G%NC
      
               COLUMN_BANDWISE_LOOP: DO IOR0 = 3, -3, -1 
      
                  IF (TWO_D .AND. ABS(IOR0) == 2) CYCLE

                  KC = IC + AB%OFFSET(IOR0)                           ! get number of neighboring cell
                  IF (KC<=0  .OR. KC >= IC) CYCLE                     ! only consider neighbors with lower cell numbers than IC
                  IF (AB%RELAX(IC, AB%POS(IOR0)) == 0) CYCLE
      
                  AB%RELAX(IC,AB%POS(IOR0)) = AB%RELAX(IC, AB%POS(IOR0))/AB%RELAX(KC,AB%POS(0))
      
                  DO JOR0 = 3, -3, -1

                     IF (TWO_D .AND. ABS(JOR0) == 2) CYCLE
      
                     JC = IC + AB%OFFSET(JOR0)                       ! get number of neighboring cell

                     IF (JC<=KC .OR. JC >G%NC) CYCLE                 ! only consider neighbors with higher cell numbers than IC
                     IF (AB%RELAX(IC, AB%POS(JOR0)) == 0) CYCLE
      
                     BFOUND = .FALSE.
                     DO KOR0 = 3, -3, -1
                        IF (TWO_D .AND. ABS(KOR0) == 2) CYCLE
                        IF (KC + AB%OFFSET(KOR0) == JC) THEN
                           BFOUND = .TRUE.
                           EXIT
                        ENDIF
                     ENDDO
                     IF (BFOUND) AB%RELAX(JC, AB%POS(JOR0)) = AB%RELAX(JC, AB%POS(JOR0)) - &
                                 AB%RELAX(IC, AB%POS(IOR0)) * AB%RELAX(KC, AB%POS(KOR0))
      
                  ENDDO
      
               ENDDO COLUMN_BANDWISE_LOOP
            ENDDO CELL_BANDWISE_LOOP
      
      END SELECT

   ENDDO LEVEL_LOOP
ENDDO MESHES_LOOP

END SUBROUTINE SCARC_SETUP_ILU


#ifdef WITH_MKL
! ------------------------------------------------------------------------------------------------
! Initialize Pardiso solver
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_CLUSTER(NLMIN, NLMAX)
USE SCARC_POINTERS, ONLY: L, G, AS, MKL
INTEGER, INTENT(IN) :: NLMIN, NLMAX
INTEGER :: NM, NL, I 
REAL (EB) :: TNOW
REAL (EB) :: DUMMY(1)=0.0_EB
REAL (FB) :: DUMMY_FB(1)=0.0_FB

TNOW = CURRENT_TIME()

MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   LEVEL_LOOP: DO NL = NLMIN, NLMAX

      CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G

      MKL => L%MKL
      AS  => G%POISSONC_SYM

      !
      ! Allocate workspace for parameters and pointers needed in MKL-routine
      !
      CALL SCARC_ALLOCATE_INT1(MKL%IPARM, 1, 64, NSCARC_INIT_ZERO, 'MKL%IPARM')

      IF (.NOT.ALLOCATED(MKL%CT)) THEN
         ALLOCATE(MKL%CT(64), STAT=IERROR)
         CALL CHKMEMERR ('SCARC', 'CT', IERROR)
         DO I=1,64
            MKL%CT(I)%DUMMY = 0
         ENDDO
      ENDIF

      !
      ! Define corresponding parameters
      ! Note: IPARM-vectory is allocate from 1:64, not from 0:63
      !
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

      ! 
      ! First perform only reordering and symbolic factorization
      ! Then perform only factorization
      ! 
      IF (TYPE_MKL_PRECISION == NSCARC_PRECISION_SINGLE) THEN

         MKL%IPARM(28) = 1         ! single precision
         MKL%PHASE = 11
         CALL CLUSTER_SPARSE_SOLVER_S(MKL%CT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, G%NC_GLOBAL, &
                                      AS%VAL_FB, AS%ROW, AS%COLG, MKL%PERM, MKL%NRHS, MKL%IPARM, &
                                      MKL%MSGLVL, DUMMY_FB, DUMMY_FB, MPI_COMM_WORLD, MKL%ERROR)
         MKL%PHASE = 22
         CALL CLUSTER_SPARSE_SOLVER_S(MKL%CT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, G%NC_GLOBAL, &
                                      AS%VAL_FB, AS%ROW, AS%COLG, MKL%PERM, MKL%NRHS, MKL%IPARM, &
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
                                      AS%VAL, AS%ROW, AS%COLG, MKL%PERM, MKL%NRHS, MKL%IPARM, &
                                      MKL%MSGLVL, DUMMY, DUMMY, MPI_COMM_WORLD, MKL%ERROR)
         MKL%PHASE = 22
         CALL CLUSTER_SPARSE_SOLVER_D(MKL%CT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, G%NC_GLOBAL, &
                                      AS%VAL, AS%ROW, AS%COLG, MKL%PERM, MKL%NRHS, MKL%IPARM, &
                                      MKL%MSGLVL, DUMMY, DUMMY, MPI_COMM_WORLD, MKL%ERROR)
         IF (MKL%ERROR /= 0) THEN
            WRITE(*,*) 'ERROR in MKL SETUP, MKL%ERROR=', MKL%ERROR
            CALL MPI_FINALIZE(IERROR)
            STOP
         ENDIF

      ENDIF

   ENDDO LEVEL_LOOP
ENDDO MESHES_LOOP

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) '========= AFTER SETUP_CLUSTER : NL=',NL
#endif

END SUBROUTINE SCARC_SETUP_CLUSTER


! ------------------------------------------------------------------------------------------------
! Initialize Pardiso solver
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_PARDISO(NLMIN, NLMAX)
USE SCARC_POINTERS, ONLY: L, G, AS, MKL
INTEGER, INTENT(IN) :: NLMIN, NLMAX
INTEGER :: NM, NL, I, IDUMMY(1)=0
REAL (EB) :: TNOW
REAL (EB) :: DUMMY(1)=0.0_EB
REAL (FB) :: DUMMY_FB(1)=0.0_FB

TNOW = CURRENT_TIME()

MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   LEVEL_LOOP: DO NL = NLMIN, NLMAX

      CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G
      MKL => L%MKL
      AS  => G%POISSONC_SYM

      !
      ! Allocate workspace for parameters nnd pointers eeded in MKL-routine
      !
      CALL SCARC_ALLOCATE_INT1(MKL%IPARM, 1, 64, NSCARC_INIT_ZERO, 'MKL%IPARM')

      IF (.NOT.ALLOCATED(MKL%PT)) THEN
         ALLOCATE(MKL%PT(64), STAT=IERROR)
         CALL CHKMEMERR ('SCARC', 'PT', IERROR)
         DO I=1,64
            MKL%PT(I)%DUMMY = 0
         ENDDO
      ENDIF

      !
      ! Define corresponding parameters
      ! Note: IPARM-vectory is allocate from 1:64, not from 0:63
      !
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

      ! 
      ! First perform only reordering and symbolic factorization
      ! Then perform only factorization
      ! 
      IF (TYPE_MKL_PRECISION == NSCARC_PRECISION_SINGLE) THEN

         MKL%IPARM(28) = 1         ! single precision
         MKL%PHASE = 11
         CALL PARDISO_S(MKL%PT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, G%NC, &
                        AS%VAL_FB, AS%ROW, AS%COLG, IDUMMY, &
                        MKL%NRHS, MKL%IPARM, MKL%MSGLVL, DUMMY_FB, DUMMY_FB, MKL%ERROR)
         MKL%PHASE = 22
         CALL PARDISO_S(MKL%PT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, G%NC, &
                        AS%VAL_FB, AS%ROW, AS%COLG, IDUMMY, &
                        MKL%NRHS, MKL%IPARM, MKL%MSGLVL, DUMMY_FB, DUMMY_FB, MKL%ERROR)

      ELSE

         MKL%IPARM(28) = 0         ! double precision
         MKL%PHASE = 11
         CALL PARDISO_D(MKL%PT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, G%NC, &
                        AS%VAL, AS%ROW, AS%COLG, IDUMMY, &
                        MKL%NRHS, MKL%IPARM, MKL%MSGLVL, DUMMY, DUMMY, MKL%ERROR)
         MKL%PHASE = 22
         CALL PARDISO_D(MKL%PT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, G%NC, &
                        AS%VAL, AS%ROW, AS%COLG, IDUMMY, &
                        MKL%NRHS, MKL%IPARM, MKL%MSGLVL, DUMMY, DUMMY, MKL%ERROR)
      ENDIF

   ENDDO LEVEL_LOOP
ENDDO MESHES_LOOP

END SUBROUTINE SCARC_SETUP_PARDISO
#endif

! ------------------------------------------------------------------------------------------------
! Define sizes for system matrix A (including extended regions related to overlaps)
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_POISSON_SIZES(NL)
USE SCARC_POINTERS, ONLY: L, G, OG, A, AB
INTEGER, INTENT(IN) :: NL
INTEGER :: NM, NOM, INBR

MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   
   CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G
   
   SELECT_MATRIX_TYPE: SELECT CASE (SCARC_MATRIX_LEVEL(NL))
   
      !
      ! -------- Matrix in compact storage technique
      !
      CASE (NSCARC_MATRIX_COMPACT)

         A => SCARC_POINT_TO_CMATRIX (G, NSCARC_MATRIX_POISSON)

         ! assign IOR settings to corresponding positions in stencil
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
WRITE(MSG%LU_DEBUG,*) 'INBR =',INBR
WRITE(MSG%LU_DEBUG,*) 'OA%N_STENCIL=',OA%N_STENCIL
WRITE(MSG%LU_DEBUG,*) 'OA%N_VAL=',OA%N_VAL
WRITE(MSG%LU_DEBUG,*) 'OA%N_ROW=',OA%N_ROW
#endif
         ENDDO

      !
      ! -------- Matrix in bandwise storage technique
      !
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
   
!
! -------- Exchange matrix sizes in case of a multi-mesh geometry
!
IF (NMESHES > 1) CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_POISSON_SIZES, NSCARC_MATRIX_POISSON, NL)

END SUBROUTINE SCARC_SETUP_POISSON_SIZES


! ------------------------------------------------------------------------------------------------------
! Usual ScaRC solver routine
! ------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SOLVER(DT_CURRENT)
USE SCARC_ITERATION_ENVIRONMENT
REAL (EB), INTENT(IN) :: DT_CURRENT
REAL (EB) :: TNOW

TNOW = CURRENT_TIME()

DT  = DT_CURRENT
DTI = 1.0_EB/DT_CURRENT

ITE_PRES = ITE_PRES + 1
ITE_GLOBAL = ICYC

SELECT_METHOD: SELECT CASE (TYPE_METHOD)

   !
   ! ---------------- Krylov method (CG/BICG) --------------------------------
   !
   CASE (NSCARC_METHOD_KRYLOV)
   
      SELECT_KRYLOV: SELECT CASE (TYPE_KRYLOV)
         CASE (NSCARC_KRYLOV_CG)
            CALL SCARC_METHOD_CG (NSCARC_STACK_ROOT, NSCARC_STACK_ZERO, NSCARC_RHS_INHOMOGENEOUS, NLEVEL_MIN)
         CASE (NSCARC_KRYLOV_BICG)
            CALL SCARC_SHUTDOWN(NSCARC_ERROR_BICG_DISABLED, SCARC_KRYLOV, NSCARC_NONE)
            !CALL SCARC_METHOD_BICG(NSCARC_STACK_ROOT, NSCARC_STACK_ZERO, NSCARC_RHS_INHOMOGENEOUS, NLEVEL_MIN)
      END SELECT SELECT_KRYLOV
   
   !
   ! ---------------- Multigrid method ---------------------------------------
   !
   CASE (NSCARC_METHOD_MULTIGRID)
   
      CALL SCARC_METHOD_MULTIGRID(NSCARC_STACK_ROOT, NSCARC_STACK_ZERO, NSCARC_RHS_INHOMOGENEOUS, NLEVEL_MIN)
   

#ifdef WITH_SCARC_MGM
   !
   ! ---------------- McKenny-Greengard-Mayo method (MGM) --------------------
   !
   CASE (NSCARC_METHOD_MGM)
   
      ! first solve inhomogeneous Poisson problem on structured grid with ScaRC (with Block-FFT)
      CALL SCARC_ASSIGN_GRID_TYPE(NSCARC_GRID_STRUCTURED)
      CALL SCARC_METHOD_CG (NSCARC_STACK_ROOT, NSCARC_STACK_ZERO, NSCARC_RHS_INHOMOGENEOUS, NLEVEL_MIN)
   
      CALL SCARC_STORE_MGM(NLEVEL_MIN, 1)
      CALL SCARC_MGM_INTERNAL_VELOCITY(NLEVEL_MIN)
      CALL SCARC_STORE_MGM(NLEVEL_MIN, 1)
   
      ! then solve homogeneous Poisson problem on unstructured grid with UScaRC (first with SSOR-preconditioning)
      ! later the preconditioning will be replaced by an individual LU-process
      CALL SCARC_ASSIGN_GRID_TYPE(NSCARC_GRID_UNSTRUCTURED)
      CALL SCARC_METHOD_CG (NSCARC_STACK_ROOT, NSCARC_STACK_ZERO, NSCARC_RHS_HOMOGENEOUS, NLEVEL_MIN)

      CALL SCARC_STORE_MGM(NLEVEL_MIN, 2)
      CALL SCARC_STORE_MGM(NLEVEL_MIN, 3)

#endif
   
   !
   ! ---------------- MKL method ---------------------------------------------
   !
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

IF (STOP_STATUS==SETUP_STOP) RETURN


T_USED(5)=T_USED(5)+CURRENT_TIME()-TNOW
CPU(MYID)%SOLVER =CPU(MYID)%SOLVER+CURRENT_TIME()-TNOW
CPU(MYID)%OVERALL=CPU(MYID)%OVERALL+CURRENT_TIME()-TNOW
END SUBROUTINE SCARC_SOLVER


! ------------------------------------------------------------------------------------------------------
! Perform preceding FFT method to improve start solution for ScaRC
! ------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_METHOD_FFT
USE MESH_POINTERS
USE POIS, ONLY: H2CZSS, H3CZSS
USE SCARC_ITERATION_ENVIRONMENT
REAL(EB), POINTER, DIMENSION(:,:,:) :: HP
INTEGER :: NM, I, J, K
LOGICAL :: WITH_BDRY = .FALSE.

MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL POINT_TO_MESH(NM)
   
   IF (PREDICTOR) THEN
      HP => H
   ELSE
      HP => HS
   ENDIF
   
   !
   ! Call the Poisson solver
   !
   IF (.NOT.TWO_D) CALL H3CZSS(BXS,BXF,BYS,BYF,BZS,BZF,ITRN,JTRN,PRHS,POIS_PTB,SAVE1,WORK,HX)
   IF (TWO_D .AND. .NOT.CYLINDRICAL) CALL H2CZSS(BXS,BXF,BZS,BZF,ITRN,PRHS,POIS_PTB,SAVE1,WORK,HX)
   
   DO K=1,KBAR
     DO J=1,JBAR
        DO I=1,IBAR
            HP(I,J,K) = PRHS(I,J,K)
        ENDDO
      ENDDO
   ENDDO
   
   !
   ! Apply boundary conditions to H
   !
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


! ------------------------------------------------------------------------------------------------
! Compute global matrix-vector product on grid level NL
!                     Y := A*X
! where NV1 is a reference to X and NV2 is a reference to Y
! including data exchange along internal boundaries
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_MATVEC_PRODUCT(NV1, NV2, NL)
USE SCARC_POINTERS, ONLY: L, G, F, OG, GWC, A, AB, V1, V2
INTEGER, INTENT(IN) :: NV1, NV2, NL           
REAL(EB) :: TNOW
INTEGER :: NM, NOM, IC, JC, IOR0, ICOL, INBR, ICE, ICW, ICG
INTEGER :: I, J, K, IW, IS=0, IT=0, IL=0, INUM1, INUM2
REAL(EB) :: TMP, VSAVE
#ifdef WITH_MKL
EXTERNAL :: DAXPBY
#endif

TNOW = CURRENT_TIME()

!
! If this call is related to a globally acting solver, exchange internal boundary values of
! vector1 such that the ghost values contain the corresponding overlapped values of adjacent neighbor
!
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'CALLING MATVEC_PRODUCT FOR ', NV1, NV2, NL
#endif
IF (TYPE_MATVEC == NSCARC_MATVEC_GLOBAL) CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_VECTOR_PLAIN, NV1, NL)

!
! Perform global matrix-vector product:
! Note: - matrix already contains subdiagonal values from neighbor along internal boundaries
!       - if vector1 contains neighboring values, then correct values of global matvec are achieved
!
SELECT_MATRIX_TYPE: SELECT CASE (SCARC_MATRIX_LEVEL(NL))

   !
   ! ------------- COMPACT storage technique
   !
   CASE (NSCARC_MATRIX_COMPACT)
   
      MESHES_COMPACT_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   
         CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G
         A => SCARC_POINT_TO_CMATRIX(G, NSCARC_MATRIX_POISSON)
         
         V1 => SCARC_POINT_TO_VECTOR (NM, NL, NV1)
         V2 => SCARC_POINT_TO_VECTOR (NM, NL, NV2)
   
         DO IC = 1, G%NC
   
            ICOL = A%ROW(IC)                                                ! diagonal entry
            JC   = A%COL(ICOL)
            V2(IC) = A%VAL(ICOL)* V1(JC)
            DO ICOL = A%ROW(IC)+1, A%ROW(IC+1)-1                            ! subdiagonal entries
               JC = A%COL(ICOL)
               V2(IC) =  V2(IC) + A%VAL(ICOL) * V1(JC)
            ENDDO
         ENDDO
   
      ENDDO MESHES_COMPACT_LOOP
   
   !
   ! ------------- bandwise storage technique
   ! matrix diagonals are supposed to be constant
   ! matrix-vector multiplication is based on daxpy-routines using the constant matrix stencil
   ! the 'wrong' entries due to boundary conditions and zero entries in subdiagonals are explicitly corrected 
   !
   CASE (NSCARC_MATRIX_BANDWISE)
   
      SELECT_STENCIL_TYPE: SELECT CASE (TYPE_STENCIL)

         !
         ! ------------- Variable entries with own implementation of daxpyv
         ! matrix-vector multiplication is based on variable matrix stencil
         !
         CASE (NSCARC_STENCIL_VARIABLE)
         
            MESHES_BANDWISE_VARIABLE_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
         
               CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G
               AB => SCARC_POINT_TO_BMATRIX(G, NSCARC_MATRIX_POISSON)
         
               V1 => SCARC_POINT_TO_VECTOR (NM, NL, NV1)               ! point to X-vector
               V2 => SCARC_POINT_TO_VECTOR (NM, NL, NV2)               ! point to Y-vector
               V2 = 0.0_EB
         
               !!$OMP PARALLEL default(none) PRIVATE(IOR0, IL, IS, IT) SHARED(AB, V1, V2) 
               DO IOR0 = 3, -3, -1

                  IF (AB%POS(IOR0) == 0) CYCLE

                  IL = AB%LENGTH(IOR0) 
                  IS = AB%SOURCE(IOR0)
                  IT = AB%TARGET(IOR0)

                  CALL SCARC_DAXPY_VARIABLE(IL, AB%VAL(IT:IT+IL-1,AB%POS(IOR0)), V1(IS:IS+IL-1), V2(IT:IT+IL-1))

               ENDDO
               !!$OMP END PARALLEL

               IF (NM==NMESHES) VSAVE = V2(G%NC)                        ! save value in case of pure Neumann bdry

               DO IOR0 = 3, -3, -1
                  F => L%FACE(IOR0)
                  DO INBR = 1, F%N_NEIGHBORS
                     NOM = F%NEIGHBORS(INBR)
                     CALL SCARC_POINT_TO_OTHER_GRID(NM, NOM, NL)
                     DO ICG = OL%GHOST_FIRSTW(IOR0), OL%GHOST_LASTW(IOR0)
                        ICW = OG%ICG_TO_ICW(ICG, 1)
                        ICE = OG%ICG_TO_ICE(ICG, 1)
                        V2(ICW) = V2(ICW) + F%SCAL_FACE * V1(ICE)
                     ENDDO
                  ENDDO
               ENDDO

               IF (IS_PURE_NEUMANN.AND.NM==NMESHES) V2(G%NC) = VSAVE   ! restore value in last cell

            ENDDO MESHES_BANDWISE_VARIABLE_LOOP

         !
         ! ------------- Storage of variable matrix entries with local multiplications
         ! matrix-vector-multiplication according to compact storage technique
         !
         CASE (NSCARC_STENCIL_VARIABLE_LOCAL)
         
            MESHES_BANDWISE_VARIABLE_LOCAL_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
         
               CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G
               AB => G%POISSONB
         
               V1 => SCARC_POINT_TO_VECTOR (NM, NL, NV1)               ! point to X-vector
               V2 => SCARC_POINT_TO_VECTOR (NM, NL, NV2)               ! point to Y-vector
         
               DO IC = 1, G%NC
                  V2(IC) = V1(IC)*AB%VAL(IC, AB%POS(0))                      ! main-diagonal contribution
                  DO IOR0 = 3, 1, -1                                         ! lower sub-diagonal contributions
                     IF (AB%POS(IOR0) == 0) CYCLE                            ! no contribution for y-direction in 2D
                     JC = IC + AB%OFFSET(IOR0)
                     IF (JC >= 1) V2(IC) = V2(IC) + V1(JC)*AB%VAL(IC, AB%POS(IOR0))
                  ENDDO
                  DO IOR0 = -1, -3, -1                                       ! upper sub-diagonal contributions
                     IF (AB%POS(IOR0) == 0) CYCLE                            ! no contribution for y-direction in 2D
                     JC = IC + AB%OFFSET(IOR0)
                     IF (JC <= G%NC) V2(IC) = V2(IC) + V1(JC)*AB%VAL(IC, AB%POS(IOR0))
                  ENDDO
               ENDDO
         
            ENDDO MESHES_BANDWISE_VARIABLE_LOCAL_LOOP
      
         !
         ! ------------- Storage of constant matrix entries - with corrections at subdiagonals and diagonal (BC's)
         !
         CASE (NSCARC_STENCIL_CONSTANT)

            MESHES_BANDWISE_CONSTANT_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
         
               CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G
               AB => G%POISSONB
         
               V1 => SCARC_POINT_TO_VECTOR (NM, NL, NV1)               ! point to X-vector
               V2 => SCARC_POINT_TO_VECTOR (NM, NL, NV2)               ! point to Y-vector
               V2 = 0.0_EB
         
               DO IOR0 = 3, -3, -1

                  IF (AB%POS(IOR0) == 0) CYCLE
                  AB%AUX(1:G%NC)=0.0_EB
                  AB%AUX(1:AB%LENGTH(IOR0)) = V1(AB%SOURCE(IOR0):AB%SOURCE(IOR0)+AB%LENGTH(IOR0)-1)

                  IF (ABS(IOR0) == 1) THEN
                     DO IC = 1, G%NC
                        IF (MOD(IC,L%NX)==0) AB%AUX(IC)=0.0_EB
                     ENDDO
                  ELSE IF (ABS(IOR0) == 2) THEN
                     INUM1 = L%NX*L%NY
                     INUM2 = L%NX*L%NY - L%NX
                     DO IC = 1, G%NC
                        IF (MOD(IC,INUM1) > INUM2) AB%AUX(IC)=0.0_EB
                     ENDDO
                  ENDIF

                  IS = AB%SOURCE(IOR0)
                  IT = AB%TARGET(IOR0)
                  IL = AB%LENGTH(IOR0)

#ifdef WITH_MKL
                  CALL DAXPY(IL, AB%STENCIL(IOR0), AB%AUX(1:IL), 1, V2(IT:IT+IL-1), 1)
#else
                  WRITE(*,*) 'TODO: MATVEC: CONSTANT: NO-MKL: CHECK HERE'
                  V2(IT:IT+IL-1) = AB%STENCIL(IOR0) * AB%AUX(1:IL)
#endif

               ENDDO
               WALL_CELLS_BANDWISE_LOOP: DO IW = 1, G%NW 

                  IOR0 = G%WALL(IW)%IOR
                  IF (TWO_D .AND. ABS(IOR0) == 2) CYCLE              ! cycle boundaries in y-direction for 2D-cases

                  GWC => G%WALL(IW)
                  F  => L%FACE(IOR0)
            
                  I = GWC%IXW
                  J = GWC%IYW
                  K = GWC%IZW
            
                  IF (IS_UNSTRUCTURED .AND. L%IS_SOLID(I, J, K)) CYCLE 
            
                  IC  = G%CELL_NUMBER(I, J, K) 
            
                  ! SPD-matrix with mixture of Dirichlet and Neumann BC's according to the SETTING of BTYPE
                  IF (N_DIRIC_GLOBAL(NLEVEL_MIN) > 0) THEN
                     TMP = V2(IC)
                     SELECT CASE (GWC%BTYPE)
                        CASE (DIRICHLET)
                           V2(IC) = V2(IC) - F%SCAL_BOUNDARY * V1(IC)
                        CASE (NEUMANN)
                           V2(IC) = V2(IC) + F%SCAL_BOUNDARY * V1(IC)
                     END SELECT
                  ENDIF 
            
               ENDDO WALL_CELLS_BANDWISE_LOOP
         
            ENDDO MESHES_BANDWISE_CONSTANT_LOOP
         !
         ! ------------- Constant entries with own implementation of daxpy
         ! 'Wrong' entries due to boundary conditions and zero entries in subdiagonals are locally corrected 
         !
         CASE (NSCARC_STENCIL_CONSTANT_LOCAL)
         
            MESHES_BANDWISE_CONSTANT_LOCAL_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
         
               CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G
               AB => G%POISSONB
         
               V1 => SCARC_POINT_TO_VECTOR (NM, NL, NV1)               ! point to X-vector
               V2 => SCARC_POINT_TO_VECTOR (NM, NL, NV2)               ! point to Y-vector
               V2 = 0.0_EB
         
               DO IOR0 = 3, -3, -1

                  IF (AB%POS(IOR0) == 0) CYCLE

                  IS = AB%SOURCE(IOR0)
                  IT = AB%TARGET(IOR0)
                  IL = AB%LENGTH(IOR0)

                  AB%AUX(1:G%NC)=0.0_EB
                  AB%AUX(1:IL) = V1(IS:IS+IL-1)

                  IF (ABS(IOR0) == 1) THEN
                     DO IC = 1, G%NC
                        IF (MOD(IC,L%NX)==0) AB%AUX(IC)=0.0_EB
                     ENDDO
                  ELSE IF (ABS(IOR0) == 2) THEN
                     INUM1 = L%NX*L%NY
                     INUM2 = L%NX*L%NY - L%NX
                     DO IC = 1, G%NC
                        IF (MOD(IC,INUM1) > INUM2) AB%AUX(IC)=0.0_EB
                     ENDDO
                  ENDIF

                  CALL SCARC_DAXPY_CONSTANT(IL, AB%STENCIL(IOR0), AB%AUX(1:IL), V2(IT:IT+IL-1))

               ENDDO
               WALL_CELLS_BANDWISE_OWN_LOOP: DO IW = 1, G%NW 

                  IOR0 = G%WALL(IW)%IOR
                  IF (TWO_D .AND. ABS(IOR0) == 2) CYCLE              ! cycle boundaries in y-direction for 2D-cases

                  GWC => G%WALL(IW)
                  F  => L%FACE(IOR0)
            
                  I = GWC%IXW
                  J = GWC%IYW
                  K = GWC%IZW
            
                  IF (IS_UNSTRUCTURED .AND. L%IS_SOLID(I, J, K)) CYCLE 
            
                  IC  = G%CELL_NUMBER(I, J, K) 
            
                  ! SPD-matrix with mixture of Dirichlet and Neumann BC's according to the SETTING of BTYPE
                  IF (N_DIRIC_GLOBAL(NLEVEL_MIN) > 0) THEN
                  TMP = V2(IC)
                     SELECT CASE (GWC%BTYPE)
                        CASE (DIRICHLET)
                           V2(IC) = V2(IC) - F%SCAL_BOUNDARY * V1(IC)
                        CASE (NEUMANN)
                           V2(IC) = V2(IC) + F%SCAL_BOUNDARY * V1(IC)
                     END SELECT
                  ENDIF 
            
               ENDDO WALL_CELLS_BANDWISE_OWN_LOOP
         
            ENDDO MESHES_BANDWISE_CONSTANT_LOCAL_LOOP

      END SELECT SELECT_STENCIL_TYPE
END SELECT SELECT_MATRIX_TYPE

CPU(MYID)%MATVEC_PRODUCT =CPU(MYID)%MATVEC_PRODUCT+CURRENT_TIME()-TNOW
END SUBROUTINE SCARC_MATVEC_PRODUCT


! ------------------------------------------------------------------------------------------------
! Compute global scalar-product (including global data exchange)
! ------------------------------------------------------------------------------------------------
REAL(EB) FUNCTION SCARC_SCALAR_PRODUCT(NV1, NV2, NL)
USE SCARC_POINTERS, ONLY: G, V1, V2
INTEGER, INTENT(IN) :: NV1, NV2, NL
REAL(EB) :: TNOW, RANK_REAL
INTEGER :: NM
#ifdef WITH_MKL
REAL(EB) :: DDOT
EXTERNAL :: DDOT
#else
INTEGER :: IC
#endif

TNOW = CURRENT_TIME()

RANK_REAL = 0.0_EB
MESH_REAL = 0.0_EB

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G

   V1 => SCARC_POINT_TO_VECTOR (NM, NL, NV1)
   V2 => SCARC_POINT_TO_VECTOR (NM, NL, NV2)

#ifdef WITH_MKL
   MESH_REAL(NM) = DDOT(G%NC, V1, 1, V2, 1)
#else
   MESH_REAL(NM) = 0.0_EB
   !$OMP PARALLEL DO PRIVATE(IC) SCHEDULE(STATIC)
   DO IC = 1, G%NC
      MESH_REAL(NM) = MESH_REAL(NM) + V1(IC) * V2(IC)
   ENDDO
   !$OMP END PARALLEL DO 
#endif

   RANK_REAL = RANK_REAL + MESH_REAL(NM)

ENDDO

!
! Compute global scalar product as sum of local scalar products
!
IF (N_MPI_PROCESSES>1) & 
   CALL MPI_ALLREDUCE(MPI_IN_PLACE, RANK_REAL, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, IERROR)

SCARC_SCALAR_PRODUCT = RANK_REAL

CPU(MYID)%SCALAR_PRODUCT = CPU(MYID)%SCALAR_PRODUCT + CURRENT_TIME()-TNOW
END FUNCTION SCARC_SCALAR_PRODUCT


! ------------------------------------------------------------------------------------------------
! Compute global L2-norm (including global data exchange)
! ------------------------------------------------------------------------------------------------
REAL(EB) FUNCTION SCARC_L2NORM(NV1, NL)
INTEGER, INTENT(IN) :: NV1, NL
REAL(EB) :: TNOW
TNOW = CURRENT_TIME()

GLOBAL_REAL = SCARC_SCALAR_PRODUCT(NV1, NV1, NL)
GLOBAL_REAL = SQRT (GLOBAL_REAL)

SCARC_L2NORM = GLOBAL_REAL

CPU(MYID)%L2NORM =CPU(MYID)%L2NORM+CURRENT_TIME()-TNOW
END FUNCTION SCARC_L2NORM


! ------------------------------------------------------------------------------------------------
! Compute linear combination of two vectors for bandwise storage technique
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_VECTOR_SUM(NV1, NV2, SCAL1, SCAL2, NL)
USE SCARC_POINTERS, ONLY: G, V1, V2
INTEGER, INTENT(IN) :: NV1, NV2, NL
REAL(EB), INTENT(IN) :: SCAL1, SCAL2
INTEGER :: NM
#ifdef WITH_MKL
EXTERNAL :: DAXPBY
#endif

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G

   V1 => SCARC_POINT_TO_VECTOR(NM, NL, NV1)
   V2 => SCARC_POINT_TO_VECTOR(NM, NL, NV2)

#ifdef WITH_MKL
   CALL DAXPBY(G%NCE, SCAL1, V1, 1, SCAL2, V2, 1)
#else
   CALL SCARC_DAXPY_CONSTANT_DOUBLE(G%NCE, SCAL1, V1, SCAL2, V2)
#endif

ENDDO

END SUBROUTINE SCARC_VECTOR_SUM


! ------------------------------------------------------------------------------------------------
! Define vector2 to be a scaled copy of vector 1
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_VECTOR_COPY(NV1, NV2, SCAL1, NL)
USE SCARC_POINTERS, ONLY: G, V1, V2
INTEGER, INTENT(IN) :: NV1, NV2, NL
REAL(EB), INTENT(IN) :: SCAL1
INTEGER :: NM
#ifdef WITH_MKL
EXTERNAL :: DCOPY, DSCAL
#endif

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G

   V1 => SCARC_POINT_TO_VECTOR(NM, NL, NV1)
   V2 => SCARC_POINT_TO_VECTOR(NM, NL, NV2)

#ifdef WITH_MKL
   CALL DCOPY(G%NCE, V1, 1, V2, 1)
   CALL DSCAL(G%NCE, SCAL1, V2, 1)
#else
   CALL SCARC_SCALING_CONSTANT(G%NCE, SCAL1, V1, V2)
#endif

ENDDO

END SUBROUTINE SCARC_VECTOR_COPY


! ------------------------------------------------------------------------------------------------
! Clear vector
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_VECTOR_CLEAR(NV, NL)
USE SCARC_POINTERS, ONLY: VC
INTEGER, INTENT(IN) :: NV, NL
INTEGER :: NM

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   VC => SCARC_POINT_TO_VECTOR(NM, NL, NV)
   VC =  0.0_EB
ENDDO

END SUBROUTINE SCARC_VECTOR_CLEAR


! ------------------------------------------------------------------------------------------------
! Preset vector with specified value
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_VECTOR_RANDOM_INIT (NV, NL)
USE SCARC_POINTERS, ONLY: L, G, VC
INTEGER, INTENT(IN) :: NV, NL
INTEGER :: IC, NM, I, J, K
REAL (EB) :: VAL

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G

   VC => SCARC_POINT_TO_VECTOR (NM, NL, NV)

   !$OMP PARALLEL DO PRIVATE(I, J, K, IC) SCHEDULE(STATIC)
   DO K = 1, L%NZ
      DO J = 1, L%NY
         DO I = 1, L%NX
            IF (IS_UNSTRUCTURED .AND. L%IS_SOLID(I, J, K)) CYCLE
            IC = G%CELL_NUMBER(I,J,K)
            CALL RANDOM_NUMBER(VAL)
            VC(IC) = VAL
         ENDDO
      ENDDO
   ENDDO
   !$OMP END PARALLEL DO 

ENDDO

END SUBROUTINE SCARC_VECTOR_RANDOM_INIT


! ------------------------------------------------------------------------------------------------
! Preset vector with specified value
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_VECTOR_INIT (NV, VAL, NL)
USE SCARC_POINTERS, ONLY: L, G, VC
INTEGER, INTENT(IN) :: NV, NL
REAL (EB), INTENT(IN) :: VAL
INTEGER :: IC, NM, I, J, K

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G

   VC => SCARC_POINT_TO_VECTOR (NM, NL, NV)

   !$OMP PARALLEL DO PRIVATE(I, J, K, IC) SCHEDULE(STATIC)
   DO K = 1, L%NZ
      DO J = 1, L%NY
         DO I = 1, L%NX
            IF (IS_UNSTRUCTURED .AND. L%IS_SOLID(I, J, K)) CYCLE
            IC = G%CELL_NUMBER(I,J,K)
            VC(IC) = VAL
         ENDDO
      ENDDO
   ENDDO
   !$OMP END PARALLEL DO 
   DO IC = G%NC+1, G%NCE
      VC(IC) = VAL
   ENDDO

ENDDO

END SUBROUTINE SCARC_VECTOR_INIT


! ------------------------------------------------------------------------------------------------
! Perform preconditioning
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_RELAXATION (NV1, NV2, NS, NP, NL)
USE SCARC_POINTERS, ONLY: L, G, A, AB, FFT, V1, V2
#ifdef WITH_MKL
USE SCARC_POINTERS, ONLY: AS, MKL, V1_FB, V2_FB
#endif
USE POIS, ONLY: H2CZSS, H3CZSS
REAL(EB) :: AUX, OMEGA_SSOR = 1.5_EB 
REAL (EB) :: TNOW
INTEGER, INTENT(IN) :: NV1, NV2, NS, NP, NL
INTEGER :: NM, IC, JC, ICOL, ITYPE, IDIAG, IPTR, INCR, IOR0, IC0, IY, IZ

TNOW = CURRENT_TIME()
ITYPE = STACK(NS-1)%SOLVER%TYPE_RELAX

SELECT CASE (ITYPE)

   ! 
   ! --------- Preconditioning by blockwise Jacobi
   ! 
   CASE (NSCARC_RELAX_JAC)

      JACOBI_MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

         CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G
         V2 => SCARC_POINT_TO_VECTOR(NM, NL, NV2)

         SELECT CASE (SCARC_MATRIX_LEVEL(NL))
            
            !
            ! matrix in COMPACT storage technique
            !
            CASE (NSCARC_MATRIX_COMPACT)

               A => SCARC_POINT_TO_CMATRIX(G, NSCARC_MATRIX_POISSON)
#ifdef WITH_SCARC_DEBUG
CALL SCARC_DEBUG_CMATRIX (A, 'A','JACOBI')
#endif
               !$OMP PARALLEL DO PRIVATE(IC) SCHEDULE(STATIC)
               DO IC = 1, G%NC
WRITE(MSG%LU_DEBUG,'(A,2I6, 2E12.4)') 'IC, V2, ROW, VAL:', IC, A%ROW(IC), V2(IC), A%VAL(A%ROW(IC))
                  V2(IC) = V2(IC) / A%VAL(A%ROW(IC))
               ENDDO
               !$OMP END PARALLEL DO

            !
            ! matrix in bandwise storage technique
            !
            CASE (NSCARC_MATRIX_BANDWISE)

               AB => SCARC_POINT_TO_BMATRIX(G, NSCARC_MATRIX_POISSON)
               !$OMP PARALLEL DO PRIVATE(IC) SCHEDULE(STATIC)
               DO IC = 1, G%NC
                  V2(IC) = V2(IC) / AB%VAL(IC, AB%POS(0))
               ENDDO
               !$OMP END PARALLEL DO

         END SELECT 

      ENDDO JACOBI_MESHES_LOOP

   ! 
   ! --------- Preconditioning by blockwise SSOR
   ! 
   CASE (NSCARC_RELAX_SSOR)

      SSOR_MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

         CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G
         V2 => SCARC_POINT_TO_VECTOR(NM, NL, NV2)

         SELECT CASE (SCARC_MATRIX_LEVEL(NL))


            !
            ! matrix in COMPACT storage technique
            !
            CASE (NSCARC_MATRIX_COMPACT)

               A => SCARC_POINT_TO_CMATRIX(G, NSCARC_MATRIX_POISSON)

               SSOR_FORWARD_COMPACT_LOOP: DO IC = 1, G%NC                                   ! forward SSOR step
                  AUX = 0.0_EB
                  DO ICOL = A%ROW(IC)+1, A%ROW(IC+1)-1
                     IF (A%COL(ICOL) >= IC) CYCLE                                          ! only process lower diags
                     IF (A%COL(ICOL) <= G%NC) THEN
                        AUX = AUX + A%VAL(ICOL) * V2(A%COL(ICOL))  ! ignore overlaps
                     ENDIF
                  ENDDO
                  V2(IC) = (V2(IC) - AUX * OMEGA_SSOR) / A%VAL(A%ROW(IC))
               ENDDO SSOR_FORWARD_COMPACT_LOOP

               SSOR_BACKWARD_COMPACT_LOOP: DO IC = G%NC-1, 1, -1                           ! backward SSOR step
                  AUX = 0.0_EB
                  DO ICOL = A%ROW(IC)+1, A%ROW(IC+1)-1
                     IF (A%COL(ICOL) <= IC) CYCLE                                         ! only process upper diags
                     IF (A%COL(ICOL) <= G%NC) THEN                                        ! ignore overlaps
                        AUX = AUX + A%VAL(ICOL) * V2(A%COL(ICOL))
                     ENDIF
                  ENDDO
                  V2(IC) = V2(IC) - AUX * OMEGA_SSOR / A%VAL(A%ROW(IC))
               ENDDO SSOR_BACKWARD_COMPACT_LOOP

         !
         ! matrix in bandwise storage technique
         !
         CASE (NSCARC_MATRIX_BANDWISE)

            AB => SCARC_POINT_TO_BMATRIX(G, NSCARC_MATRIX_POISSON)

            !
            ! 2D version
            !
            IF (TWO_D) THEN

               SSOR_FORWARD_BANDWISE_2D_LOOP: DO IC = 1, G%NC                 ! forward SSOR step
                  AUX = 0.0_EB
                  DO IOR0 = 1, 3, 2                                          ! only process lower x- and z-diag
                     JC = IC + AB%OFFSET(IOR0)
                     IF (JC >= 1 .AND. JC <= G%NC) AUX = AUX + AB%VAL(IC, AB%POS(IOR0)) * V2(JC)
                  ENDDO
                  V2(IC) = (V2(IC) - AUX * OMEGA_SSOR) / AB%VAL(IC, AB%POS(0))
               ENDDO SSOR_FORWARD_BANDWISE_2D_LOOP

               SSOR_BACKWARD_BANDWISE_2D_LOOP: DO IC = G%NC-1, 1, -1          ! backward SSOR step
                  AUX = 0.0_EB
                  DO IOR0 = -1, -3, -2                                       ! only process upper x- and z-diag
                     JC = IC + AB%OFFSET(IOR0)
                     IF (JC <= G%NC) THEN
                        AUX = AUX + AB%VAL(IC, AB%POS(IOR0)) * V2(JC)
                     ENDIF
                  ENDDO
                  V2(IC) = V2(IC) - AUX * OMEGA_SSOR / AB%VAL(IC, AB%POS(0))
               ENDDO SSOR_BACKWARD_BANDWISE_2D_LOOP

            !
            ! 3D version
            !
            ELSE

               SSOR_FORWARD_BANDWISE_3D_LOOP: DO IC = 1, G%NC                  ! forward SSOR step
                  AUX = 0.0_EB
                  DO IOR0 = 1, 3                                             ! only process lower diags
                     IF (AB%POS(IOR0) == 0) CYCLE                            ! no contribution for y-direction in 2D
                     JC = IC + AB%OFFSET(IOR0)
                     IF (JC >= 1 .AND. JC <= G%NC) AUX = AUX + AB%VAL(IC, AB%POS(IOR0)) * V2(JC)
                  ENDDO
                  V2(IC) = (V2(IC) - AUX * OMEGA_SSOR) / AB%VAL(IC, AB%POS(0))
               ENDDO SSOR_FORWARD_BANDWISE_3D_LOOP

               SSOR_BACKWARD_BANDWISE_3D_LOOP: DO IC = G%NC-1, 1, -1           ! backward SSOR step
                  AUX = 0.0_EB
                  DO IOR0 = -1, -3, -1                                       ! only process upper diags
                     IF (AB%POS(IOR0) == 0) CYCLE                            ! no contribution for y-direction in 2D
                     JC = IC + AB%OFFSET(IOR0)
                     IF (JC >= IC .AND. JC <= G%NC) AUX = AUX + AB%VAL(IC, AB%POS(IOR0)) * V2(JC)
                  ENDDO
                  V2(IC) = V2(IC) - AUX * OMEGA_SSOR / AB%VAL(IC, AB%POS(0))
               ENDDO SSOR_BACKWARD_BANDWISE_3D_LOOP

            ENDIF

         END SELECT 

      ENDDO SSOR_MESHES_LOOP

   ! 
   ! --------- Preconditioning by Jacobi in matrix form
   ! 
   CASE (NSCARC_RELAX_MJAC)

      MJAC_MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

         CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G

         V1 => SCARC_POINT_TO_VECTOR(NM, NL, NV1)
         V2 => SCARC_POINT_TO_VECTOR(NM, NL, NV2)

         SELECT CASE(TYPE_MATRIX)

            ! ------------ matrix in COMPACT storage technique
            CASE (NSCARC_MATRIX_COMPACT)
               A => SCARC_POINT_TO_CMATRIX(G, NSCARC_MATRIX_POISSON)
               CALL SCARC_SCALING_VARIABLE(G%NC, A%RELAX, V1, V2)

            ! ------------ matrix in bandwise storage technique
            CASE (NSCARC_MATRIX_BANDWISE)
               AB => SCARC_POINT_TO_BMATRIX(G, NSCARC_MATRIX_POISSON)
               CALL SCARC_SCALING_VARIABLE(G%NC, AB%RELAXD, V1, V2)

         END SELECT

      ENDDO MJAC_MESHES_LOOP

   ! 
   ! --------- Preconditioning by different matrix-form preconditioners
   ! in all cases the preconditioner is given as separate matrix which is based
   ! on the same storage technique as the matrix AC itself;
   ! two tridiagonal systems have to be solved
   ! V1 contains the RHS to be solved for, V2 will contain the solution
   ! 
   CASE (NSCARC_RELAX_MGS, NSCARC_RELAX_MSGS, NSCARC_RELAX_MSOR, NSCARC_RELAX_MSSOR, NSCARC_RELAX_ILU)

      LU_MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

         CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G

         V1 => SCARC_POINT_TO_VECTOR(NM, NL, NV1)
         V2 => SCARC_POINT_TO_VECTOR(NM, NL, NV2)
      
         SELECT CASE(TYPE_MATRIX)

            ! ------------ matrix in COMPACT storage technique
            CASE (NSCARC_MATRIX_COMPACT)

               A => SCARC_POINT_TO_CMATRIX(G, NSCARC_MATRIX_POISSON)
      
               ! Forward solve:   Solve V2 = L^-1 V1
               ! Compute sol(i) = rhs(i) - sum L(i,j) x sol(j)
               DO IC = 1, G%NC
                  V2(IC) = V1(IC)
                  DO IPTR = A%ROW(IC), A%ROW(IC+1)-1
                     JC = A%COL(IPTR)
                     IF (JC >= IC) CYCLE
                     V2(IC) = V2(IC) - A%RELAX(IPTR) * V2(JC)
                  ENDDO
               ENDDO
      
               ! If preconditioner is not symmetric, upper matrix U is zero and nothing has to be solved
               IF (ITYPE == NSCARC_RELAX_MGS .OR. ITYPE == NSCARC_RELAX_MSOR) CYCLE
      
               ! Backward solve
               ! Compute sol: inv(U) sol
               DO IC = G%NC, 1, -1
      
                  DO IPTR = A%ROW(IC), A%ROW(IC+1)-1
                     JC = A%COL(IPTR)
                     IF (JC <= IC) CYCLE
                     V2(IC) = V2(IC) - A%RELAX(IPTR) * V2(JC)
                  ENDDO
      
                  ! Compute sol(i) = sol(i)/U(i,i)
                  IDIAG = A%ROW(IC)
                  V2(IC) = V2(IC)/A%RELAX(IDIAG)
      
               ENDDO
      

            !
            ! matrix in bandwise storage technique
            !
            CASE (NSCARC_MATRIX_BANDWISE)

               AB => SCARC_POINT_TO_BMATRIX(G, NSCARC_MATRIX_POISSON)
      
               IF (TWO_D) THEN
                  INCR = -2
               ELSE 
                  INCR = -1
               ENDIF
               
               ! Forward solve:   Solve V2 = L^-1 V1
               ! Compute sol(i) = rhs(i) - sum L(i,j) x sol(j)
               !!$OMP PARALLEL DO PRIVATE(IC, JC, IOR0) SCHEDULE(STATIC)
               DO IC = 1, G%NC
                  V2(IC) = V1(IC)
                  DO IOR0 = 3, 1, INCR
                     JC = IC + AB%OFFSET(IOR0)
                     IF (JC <= 0) CYCLE
                     V2(IC) = V2(IC) - AB%RELAX(IC, AB%POS(IOR0)) * V2(JC)
                  ENDDO
               ENDDO
               !!$OMP END PARALLEL DO 
      
               ! If preconditioner is not symmetric, upper matrix U is zero and nothing has to be solved
               IF (ITYPE == NSCARC_RELAX_MGS .OR. ITYPE == NSCARC_RELAX_MSOR) CYCLE
      
               ! Backward solve
               ! Compute sol: inv(U) sol
               !!$OMP PARALLEL DO PRIVATE(IC, JC, IOR0) SCHEDULE(STATIC)
               DO IC = G%NC, 1, -1
      
                  DO IOR0 = -1, -3, INCR
                     JC = IC + AB%OFFSET(IOR0)
                     IF (JC > G%NC) CYCLE
                     V2(IC) = V2(IC) - AB%RELAX(IC, AB%POS(IOR0)) * V2(JC)
                  ENDDO
      
                  ! Compute sol(i) = sol(i)/U(i,i)
                  V2(IC) = V2(IC)/AB%RELAX(IC, AB%POS(0))
               ENDDO
               !!$OMP END PARALLEL DO 
      
         END SELECT
      
      ENDDO LU_MESHES_LOOP
      

   ! 
   ! --------- Preconditioning by blockwise Geometric Multigrid
   ! 
   CASE (NSCARC_RELAX_GMG)

      CALL SCARC_METHOD_MULTIGRID (NS, NP, NSCARC_RHS_DEFECT, NLEVEL_MIN)

   ! 
   ! --------- Preconditioning by blockwise FFT based on Crayfishpak
   ! 
   CASE (NSCARC_RELAX_FFT)

      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

         CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G
         FFT => L%FFT

         V1  => SCARC_POINT_TO_VECTOR(NM, NL, NV1)
         V2  => SCARC_POINT_TO_VECTOR(NM, NL, NV2)

         !$OMP PARALLEL DO PRIVATE(IC) SCHEDULE(STATIC)
         DO IC = 1, G%NC
            FFT%PRHS(G%ICX(IC), G%ICY(IC), G%ICZ(IC)) = V1(IC)
         ENDDO
         !$OMP END PARALLEL DO

         IF (TWO_D) THEN
            CALL H2CZSS (FFT%BXS,  FFT%BXF, FFT%BZS, FFT%BZF, FFT%ITRN, &
                         FFT%PRHS, FFT%POIS_PTB, FFT%SAVE1, FFT%WORK, FFT%HX)
         ELSE
            CALL H3CZSS (FFT%BXS,  FFT%BXF, FFT%BYS, FFT%BYF, FFT%BZS, FFT%BZF, FFT%ITRN, FFT%JTRN, &
                         FFT%PRHS, FFT%POIS_PTB, FFT%SAVE1, FFT%WORK, FFT%HX)
         ENDIF

         !$OMP PARALLEL DO PRIVATE(IC) SCHEDULE(STATIC)
         DO IC = 1, G%NC
            V2(IC) = FFT%PRHS(G%ICX(IC), G%ICY(IC), G%ICZ(IC)) 
         ENDDO
         !$OMP END PARALLEL DO 

      ENDDO

   ! 
   ! --------- Preconditioning by blockwise overlapping FFT based on Crayfishpak 
   !           still test-version for tunnel-shaped geometries of type Mx1
   ! 
   CASE (NSCARC_RELAX_FFTO)

      ! Exchange overlapping parts
      CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_VECTOR_PLAIN, NV1, NL)

      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

         CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G
         FFT => L%FFT

         V1  => SCARC_POINT_TO_VECTOR(NM, NL, NV1)
         V2  => SCARC_POINT_TO_VECTOR(NM, NL, NV2)

         IF (NM == 1) THEN

            !$OMP PARALLEL DO PRIVATE(IC) SCHEDULE(STATIC)
            DO IC = 1, G%NC
               FFT%PRHS(G%ICX(IC), G%ICY(IC), G%ICZ(IC)) = V1(IC)
            ENDDO
            !$OMP END PARALLEL DO

            IC0 = G%NC+1
            DO IZ = 1, L%NZ
               DO IY = 1, L%NY
                  FFT%PRHS(FFT%IBAR, IY, IZ) = V1(IC0)
               IC0 = IC0 + 1
               ENDDO
            ENDDO

         ELSE IF (NM == NMESHES) THEN

            !$OMP PARALLEL DO PRIVATE(IC) SCHEDULE(STATIC)
            DO IC = 1, G%NC
               FFT%PRHS(G%ICX(IC)+1, G%ICY(IC), G%ICZ(IC)) = V1(IC)
            ENDDO
            !$OMP END PARALLEL DO

            IC0 = G%NC+1
            DO IZ = 1, L%NZ
               DO IY = 1, L%NY
                  FFT%PRHS(1, IY, IZ) = V1(IC0)
                  IC0 = IC0 + 1
               ENDDO
            ENDDO

         ELSE 

            !$OMP PARALLEL DO PRIVATE(IC) SCHEDULE(STATIC)
            DO IC = 1, G%NC
               FFT%PRHS(G%ICX(IC)+1, G%ICY(IC), G%ICZ(IC)) = V1(IC)
            ENDDO
            !$OMP END PARALLEL DO

            IC0 = G%NC+1
            DO IZ = 1, L%NZ
               DO IY = 1, L%NY
                  FFT%PRHS(1, IY, IZ) = V1(IC0)
                  IC0 = IC0 + 1
               ENDDO
            ENDDO
            DO IZ = 1, L%NZ
               DO IY = 1, L%NY
                  FFT%PRHS(FFT%IBAR, IY, IZ) = V1(IC0)
                  IC0 = IC0 + 1
               ENDDO
            ENDDO

         ENDIF

         ! 
         ! Call corresponding FFT solver
         ! 
         IF (TWO_D) THEN
            CALL H2CZSS (FFT%BXS,  FFT%BXF, FFT%BZS, FFT%BZF, FFT%ITRN, &
                         FFT%PRHS, FFT%POIS_PTB, FFT%SAVE1, FFT%WORK, FFT%HX)
         ELSE
            CALL H3CZSS (FFT%BXS,  FFT%BXF, FFT%BYS, FFT%BYF, FFT%BZS, FFT%BZF, FFT%ITRN, FFT%JTRN, &
                         FFT%PRHS, FFT%POIS_PTB, FFT%SAVE1, FFT%WORK, FFT%HX)
         ENDIF

         ! 
         ! Extract computed data from FFT%PRHS
         ! 
         IF (NM == 1) THEN
            !$OMP PARALLEL DO PRIVATE(IC) SCHEDULE(STATIC)
            DO IC = 1, G%NC
               V2(IC) = FFT%PRHS(G%ICX(IC), G%ICY(IC), G%ICZ(IC)) 
            ENDDO
            !$OMP END PARALLEL DO 

            IC0 = G%NC+1
            DO IZ = 1, L%NZ
               DO IY = 1, L%NY
                  V2(IC0) = FFT%PRHS(FFT%IBAR, IY, IZ) 
                  IC0 = IC0 + 1
               ENDDO
            ENDDO

         ELSE IF (NM == NMESHES) THEN

            !$OMP PARALLEL DO PRIVATE(IC) SCHEDULE(STATIC)
            DO IC = 1, G%NC
               V2(IC) = FFT%PRHS(G%ICX(IC)+1, G%ICY(IC), G%ICZ(IC)) 
            ENDDO
            !$OMP END PARALLEL DO 

            IC0 = G%NC+1
            DO IZ = 1, L%NZ
               DO IY = 1, L%NY
                  V2(IC0) = FFT%PRHS(1, IY, IZ) 
                  IC0 = IC0 + 1
               ENDDO
            ENDDO

         ELSE

            !$OMP PARALLEL DO PRIVATE(IC) SCHEDULE(STATIC)
            DO IC = 1, G%NC
               V2(IC) = FFT%PRHS(G%ICX(IC)+1, G%ICY(IC), G%ICZ(IC)) 
            ENDDO
            !$OMP END PARALLEL DO 

            IC0 = G%NC+1
            DO IZ = 1, L%NZ
               DO IY = 1, L%NY
                  V2(IC0) = FFT%PRHS(1, IY, IZ) 
                  IC0 = IC0 + 1
               ENDDO
            ENDDO
            DO IZ = 1, L%NZ
               DO IY = 1, L%NY
                  V2(IC0) = FFT%PRHS(FFT%IBAR, IY, IZ) 
                  IC0 = IC0 + 1
               ENDDO
            ENDDO

         ENDIF

      ENDDO

      ! Exchange overlapping parts
      CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_VECTOR_MEAN, NV2, NL)


#ifdef WITH_MKL
   ! 
   ! --------- Preconditioning by LU-decomposition
   ! 
   CASE (NSCARC_RELAX_MKL)

      ! Preconditioning by Cluster Sparse Solver from MKL
      !
      MKL_SCOPE_IF: IF (STACK(NS)%SOLVER%TYPE_SCOPE(0) == NSCARC_SCOPE_GLOBAL) THEN

         MKL_SCOPE_GLOBAL_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

            CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G
            MKL => L%MKL
            AS => SCARC_POINT_TO_CMATRIX(G, NSCARC_MATRIX_POISSON_SYM)

            MKL%PHASE  = 33                            ! only solving

            V1 => SCARC_POINT_TO_VECTOR (NM, NL, NV1)
            V2 => SCARC_POINT_TO_VECTOR (NM, NL, NV2)

            IF (TYPE_MKL_PRECISION == NSCARC_PRECISION_SINGLE) THEN

               V1_FB => SCARC_POINT_TO_VECTOR_FB (NM, NL, NV1)
               V2_FB => SCARC_POINT_TO_VECTOR_FB (NM, NL, NV2)

               V1_FB(1:G%NC) = REAL(V1(1:G%NC), FB)
               V2_FB(1:G%NC) = 0.0_FB

               CALL CLUSTER_SPARSE_SOLVER_S(MKL%CT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, G%NC_GLOBAL, &
                                            AS%VAL_FB, AS%ROW, AS%COLG, MKL%PERM, MKL%NRHS, MKL%IPARM, &
                                            MKL%MSGLVL, V1_FB, V2_FB, MPI_COMM_WORLD, MKL%ERROR)

               V2(1:G%NC) = REAL(V2_FB(1:G%NC), EB)

            ELSE
               CALL CLUSTER_SPARSE_SOLVER_D(MKL%CT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, G%NC_GLOBAL, &
                                            AS%VAL, AS%ROW, AS%COLG, MKL%PERM, MKL%NRHS, MKL%IPARM, &
                                            MKL%MSGLVL, V1, V2, MPI_COMM_WORLD, MKL%ERROR)
            ENDIF
            IF (MKL%ERROR /= 0) CALL SCARC_SHUTDOWN(NSCARC_ERROR_MKL_INTERNAL, SCARC_NONE, MKL%ERROR)

         ENDDO MKL_SCOPE_GLOBAL_LOOP

      !
      ! Preconditioning by Pardiso Solver from MKL
      !
      ELSE MKL_SCOPE_IF

         MKL_SCOPE_LOCAL_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

            CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G
            MKL => L%MKL
            AS => SCARC_POINT_TO_CMATRIX(G, NSCARC_MATRIX_POISSON_SYM)

            MKL%PHASE  = 33                            ! only solving

            V1 => SCARC_POINT_TO_VECTOR (NM, NL, NV1)
            V2 => SCARC_POINT_TO_VECTOR (NM, NL, NV2)

            IF (TYPE_MKL_PRECISION == NSCARC_PRECISION_SINGLE) THEN

               V1_FB => SCARC_POINT_TO_VECTOR_FB (NM, NL, NV1)
               V2_FB => SCARC_POINT_TO_VECTOR_FB (NM, NL, NV2)

               V1_FB(1:G%NC) = REAL(V1(1:G%NC), FB)
               V2_FB(1:G%NC) = 0.0_FB
   
               CALL PARDISO_S(MKL%PT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, G%NC, &
                              AS%VAL_FB, AS%ROW, AS%COLG, &
                              MKL%PERM, MKL%NRHS, MKL%IPARM, MKL%MSGLVL, V1_FB, V2_FB, MKL%ERROR)

               V2(1:G%NC) = REAL(V2_FB(1:G%NC), EB)
            ELSE

               V1 => SCARC_POINT_TO_VECTOR (NM, NL, NV1)
               V2 => SCARC_POINT_TO_VECTOR (NM, NL, NV2)

               CALL PARDISO_D(MKL%PT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, G%NC, &
                              AS%VAL, AS%ROW, AS%COLG, &
                              MKL%PERM, MKL%NRHS, MKL%IPARM, MKL%MSGLVL, V1, V2, MKL%ERROR)

            ENDIF
            IF (MKL%ERROR /= 0) CALL SCARC_SHUTDOWN(NSCARC_ERROR_MKL_INTERNAL, SCARC_NONE, MKL%ERROR)

         ENDDO MKL_SCOPE_LOCAL_LOOP

      ENDIF MKL_SCOPE_IF

#endif

END SELECT

CPU(MYID)%RELAXATION =CPU(MYID)%RELAXATION+CURRENT_TIME()-TNOW
END SUBROUTINE SCARC_RELAXATION


#ifdef WITH_MKL
! ------------------------------------------------------------------------------------------------
! Perform global Pardiso-method based on MKL
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_METHOD_CLUSTER(NSTACK, NPARENT, NLEVEL)
USE SCARC_POINTERS, ONLY: L, G, MKL, V1, V2, AS, V1_FB, V2_FB
USE SCARC_ITERATION_ENVIRONMENT
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

!CALL SCARC_DEBUG_CMATRIX(AS, 'AS','IN CLUSTER')

   V1 => SCARC_POINT_TO_VECTOR (NM, NL, B)
   V2 => SCARC_POINT_TO_VECTOR (NM, NL, X)

   MKL => L%MKL
   MKL%PHASE  = 33                                ! only solving

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'CLUSTER, PRE, V1:'
WRITE(MSG%LU_DEBUG,'(6E12.4)') V1
WRITE(MSG%LU_DEBUG,*) 'CLUSTER, PRE, V2:'
WRITE(MSG%LU_DEBUG,'(6E12.4)') V2
#endif
   IF (TYPE_MKL_PRECISION == NSCARC_PRECISION_SINGLE) THEN

      V1_FB => SCARC_POINT_TO_VECTOR_FB (NM, NL, B)
      V2_FB => SCARC_POINT_TO_VECTOR_FB (NM, NL, X)

      V1_FB = REAL(V1, FB)
      V2_FB = 0.0_FB
      CALL CLUSTER_SPARSE_SOLVER_S(MKL%CT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, G%NC_GLOBAL, &
                                   AS%VAL_FB, AS%ROW, AS%COLG, MKL%PERM, MKL%NRHS, MKL%IPARM, &
                                   MKL%MSGLVL, V1_FB, V2_FB, MPI_COMM_WORLD, MKL%ERROR)
      V2 = REAL(V2_FB, EB)

   ELSE

      CALL CLUSTER_SPARSE_SOLVER_D(MKL%CT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, G%NC_GLOBAL, &
                                   AS%VAL, AS%ROW, AS%COLG, MKL%PERM, MKL%NRHS, MKL%IPARM, &
                                   MKL%MSGLVL, V1, V2, MPI_COMM_WORLD, MKL%ERROR)
   ENDIF

   IF (MKL%ERROR /= 0) CALL SCARC_SHUTDOWN(NSCARC_ERROR_MKL_INTERNAL, SCARC_NONE, MKL%ERROR)

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'CLUSTER, POST, V1:'
WRITE(MSG%LU_DEBUG,'(2E16.8)') V1
WRITE(MSG%LU_DEBUG,*) 'CLUSTER, POST, V2:'
WRITE(MSG%LU_DEBUG,'(2E16.8)') V2
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
! Perform global Pardiso-method based on MKL
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_METHOD_PARDISO(NSTACK, NPARENT, NLEVEL)
USE SCARC_POINTERS, ONLY: L, G, MKL, AS, V1, V2, V1_FB, V2_FB
USE SCARC_ITERATION_ENVIRONMENT
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

#ifdef WITH_SCARC_DEBUG
CALL SCARC_DEBUG_CMATRIX(AS, 'AS','IN PARDISO')
#endif

   V1 => SCARC_POINT_TO_VECTOR (NM, NL, B)
   V2 => SCARC_POINT_TO_VECTOR (NM, NL, X)

   MKL => L%MKL
   MKL%PHASE  = 33         ! only solving

   IF (TYPE_MKL_PRECISION == NSCARC_PRECISION_SINGLE) THEN

      V1_FB => SCARC_POINT_TO_VECTOR_FB (NM, NL, B)
      V2_FB => SCARC_POINT_TO_VECTOR_FB (NM, NL, X)

      V1_FB = REAL(V1, FB)
      V2_FB = 0.0_FB
      CALL PARDISO_S(MKL%PT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, G%NC, &
                     AS%VAL_FB, AS%ROW, AS%COLG, MKL%PERM, MKL%NRHS, MKL%IPARM, &
                     MKL%MSGLVL, V1_FB, V2_FB, MKL%ERROR)

      V2 = REAL(V2_FB, EB)

   ELSE

      V2 = 0.0_EB
      CALL PARDISO_D(MKL%PT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, G%NC, &
                     AS%VAL, AS%ROW, AS%COLG, MKL%PERM, MKL%NRHS, MKL%IPARM, &
                     MKL%MSGLVL, V1, V2, MKL%ERROR)
   ENDIF

   IF (MKL%ERROR /= 0) CALL SCARC_SHUTDOWN(NSCARC_ERROR_MKL_INTERNAL, SCARC_NONE, MKL%ERROR)

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'PARDISO, POST, V1:'
WRITE(MSG%LU_DEBUG,'(2E16.8)') V1
WRITE(MSG%LU_DEBUG,*) 'PARDISO, POST, V2:'
WRITE(MSG%LU_DEBUG,'(2E16.8)') V2
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
! Increase corresponding iteration count (just for visualization of convergence behavior)
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_INCREASE_ITERATION_COUNTS(ITE0)
USE SCARC_ITERATION_ENVIRONMENT
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


! ------------------------------------------------------------------------------------------------
! Perform global CG-method based on global Possion-matrix
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_METHOD_CG(NSTACK, NPARENT, NRHS, NLEVEL)
USE SCARC_ITERATION_ENVIRONMENT
INTEGER, INTENT(IN) :: NSTACK, NPARENT, NRHS, NLEVEL
INTEGER :: NSTATE, NS, NP, NL
REAL (EB) :: ALPHA, BETA, SIGMA, SIGMA0=0.0_EB
REAL (EB) :: TNOW, TNOWI

TNOW = CURRENT_TIME()
ITE_CG = 0

! get current and parent stack position, and current level
TYPE_MATVEC = NSCARC_MATVEC_GLOBAL
NS  = NSTACK
NP  = NPARENT
NL  = NLEVEL

#ifdef WITH_SCARC_STANDALONE
IF (ICYC == 1) THEN
   CALL SCARC_DUMP_SYSTEM(NS, NSCARC_DUMP_MESH)
   CALL SCARC_DUMP_SYSTEM(NS, NSCARC_DUMP_A)
ENDIF
#endif

! ------------------------------------------------------------------------------------------------
! Initialization:
!   - Get parameters for current scope (note: NL denotes the finest level)
!   - Get right hand side vector and clear solution vectors
! ------------------------------------------------------------------------------------------------
CALL SCARC_SETUP_SOLVER(NS, NP)
CALL SCARC_SETUP_WORKSPACE(NS, NL, NRHS)

! In case of pure Neumann boundary conditions setup condensed system
IF (N_DIRIC_GLOBAL(NLEVEL_MIN) == 0) THEN
   CALL SCARC_VECTOR_INIT (X, 0.0_EB, NL)                    
   CALL SCARC_FILTER_MEANVALUE(B, NL)                       
   CALL SCARC_SETUP_SYSTEM_CONDENSED (B, NL, 1)            
ENDIF

#ifdef WITH_SCARC_STANDALONE
CALL SCARC_DUMP_SYSTEM(NS, NSCARC_DUMP_B)
#endif

#ifdef WITH_SCARC_DEBUG
CALL SCARC_DEBUG_LEVEL (X, 'CG-METHOD: X INIT1 ', NL)
CALL SCARC_DEBUG_LEVEL (B, 'CG-METHOD: B INIT1 ', NL)
#endif

! Compute initial residual 
TYPE_MATVEC = NSCARC_MATVEC_GLOBAL
CALL SCARC_MATVEC_PRODUCT (X, R, NL)                         !  r^0 := A*x^0
CALL SCARC_VECTOR_SUM     (B, R, -1.0_EB, 1.0_EB, NL)        !  r^0 := r^0 - b     corresponds to  A*x^0 - b

RES    = SCARC_L2NORM (R, NL)                                !  res   := ||r^0||
RESIN  = RES                                                 !  resin := res
NSTATE = SCARC_CONVERGENCE_STATE (0, NS, NL)                 !  res < tolerance ?

#ifdef WITH_SCARC_DEBUG
CALL SCARC_DEBUG_LEVEL (X, 'CG-METHOD: X INIT1 ', NL)
CALL SCARC_DEBUG_LEVEL (B, 'CG-METHOD: B INIT1 ', NL)
#endif

! Perform initial preconditioning
IF (NSTATE /= NSCARC_STATE_CONV0) THEN                       !  if no convergence yet, start precon
   CALL SCARC_PRECONDITIONER(NS, NS, NL)                     !  v^0 := Precon(r^0)
   SIGMA0 = SCARC_SCALAR_PRODUCT(R, V, NL)                   !  SIGMA0 := (r^0,v^0)
   CALL SCARC_VECTOR_COPY (V, D, -1.0_EB, NL)                !  d^0 := -v^0
ENDIF

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'RESIN, RES, ITE, SIGMA0:', RESIN, RES, ITE, SIGMA0
CALL SCARC_DEBUG_LEVEL (D, 'CG-METHOD: D INIT1 ', NL)
#endif

! ------------------------------------------------------------------------------------------------
! Perform conjugate gradient looping
! ------------------------------------------------------------------------------------------------
CG_LOOP: DO ITE = 1, NIT

   TNOWI = CURRENT_TIME()
   CALL SCARC_INCREASE_ITERATION_COUNTS(ITE)

   !TYPE_MATVEC = NSCARC_MATVEC_LOCAL
   CALL SCARC_MATVEC_PRODUCT (D, Y, NL)                      !  y^k := A*d^k

   ALPHA = SCARC_SCALAR_PRODUCT (D, Y, NL)                   !  alpha := (d^k,y^k)     corresponds to   (d^k,A*d^k)
   ALPHA = SIGMA0/ALPHA                                      !  alpha := (r^k,v^k)/(d^k,A*d^k)

   CALL SCARC_VECTOR_SUM (D, X, ALPHA, 1.0_EB, NL)           !  x^{k+1} := x^k + alpha * d^k
   CALL SCARC_VECTOR_SUM (Y, R, ALPHA, 1.0_EB, NL)           !  r^{k+1} := r^k + alpha * y^k   ~  r^k + alpha * A * d^k

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'ITE, ITE_CG=', ITE, ITE_CG
CALL SCARC_DEBUG_LEVEL (X, 'CG-METHOD: X ITE ', NL)
CALL SCARC_DEBUG_LEVEL (Y, 'CG-METHOD: Y ITE ', NL)
CALL SCARC_DEBUG_LEVEL (R, 'CG-METHOD: R ITE ', NL)
#endif

   RES = SCARC_L2NORM (R, NL)                                !  res := ||r^{k+1}||
   NSTATE = SCARC_CONVERGENCE_STATE (0, NS, NL)              !  res < tolerance ??
   IF (NSTATE /= NSCARC_STATE_PROCEED) EXIT CG_LOOP

   CALL SCARC_PRECONDITIONER(NS, NS, NL)                     !  v^{k+1} := Precon(r^{k+1})

   SIGMA  = SCARC_SCALAR_PRODUCT (R, V, NL)                  !  sigma := (r^{k+1},v^{k+1})
   BETA   = SIGMA/SIGMA0                                     !  beta  := (r^{k+1},v^{k+1})/(r^k,v^k)
   SIGMA0 = SIGMA                                            !  save last sigma

#ifdef WITH_SCARC_DEBUG
CALL SCARC_DEBUG_LEVEL (V, 'CG-METHOD: V ITE ', NL)
CALL SCARC_DEBUG_LEVEL (D, 'CG-METHOD: D ITE ', NL)
#endif

   CALL SCARC_VECTOR_SUM (V, D, -1.0_EB, BETA, NL)           !  d^{k+1} := -v^{k+1} + beta * d^{k+1}

   CPU(MYID)%ITERATION=MAX(CPU(MYID)%ITERATION,CURRENT_TIME()-TNOWI)

ENDDO CG_LOOP

! ------------------------------------------------------------------------------------------------
! Determine convergence rate and print corresponding information
! In case of CG as main solver:
!   - Transfer ScaRC solution vector X to FDS pressure vector
!   - Set ghost cell values along external boundaries
!   - Exchange values along internal boundaries
! ------------------------------------------------------------------------------------------------
CALL SCARC_CONVERGENCE_RATE(NSTATE, NS, NL)

IF (N_DIRIC_GLOBAL(NLEVEL_MIN) == 0) THEN
   CALL SCARC_RESTORE_LAST_CELL(X, NL)
   CALL SCARC_FILTER_MEANVALUE(X, NL)
ENDIF

#ifdef WITH_SCARC_DEBUG
CALL SCARC_DEBUG_LEVEL (X, 'CG-METHOD: X FINAL', NL)
#endif

IF (TYPE_SOLVER == NSCARC_SOLVER_MAIN) THEN
   CALL SCARC_UPDATE_MAINCELLS(NLEVEL_MIN)
   CALL SCARC_UPDATE_GHOSTCELLS(NLEVEL_MIN)
#ifdef WITH_SCARC_STANDALONE
   CALL SCARC_PRESSURE_DIFFERENCE(NLEVEL_MIN)
   CALL SCARC_DUMP_SYSTEM(NS, NSCARC_DUMP_X)
#endif
ENDIF

CALL SCARC_RELEASE_SOLVER(NS, NP)

END SUBROUTINE SCARC_METHOD_CG


! -----------------------------------------------------------------------------------------------
! Preconditioning method which is based on the following input and output convention:
!  - the residual which has to be preconditioned is passed in via vector R
!  - the result of preconditioning is passed out via vector V
!  - for several variants Y and Z are used as auxiliary vectors
!  - in the comments: call is based on current grid level l (mostly the finest one)
!  -                  l=1 denotes the finest  grid level NLEVEL_MIN
!  -                  l=L denotes the coarset grid level NLEVEL_MAX
! -----------------------------------------------------------------------------------------------
SUBROUTINE SCARC_PRECONDITIONER(NS, NP, NL)
USE SCARC_ITERATION_ENVIRONMENT
INTEGER, INTENT(IN) :: NS, NP, NL     
INTEGER :: IL

SELECT_PRECON_TYPE: SELECT CASE (TYPE_TWOLEVEL)

   ! 
   ! ---------- Classical one-level preconditioning
   ! 
   CASE (NSCARC_TWOLEVEL_NONE)

      CALL SCARC_VECTOR_COPY (R, V, 1.0_EB, NL)                   !  v := r
      CALL SCARC_RELAXATION (R, V, NS+1, NP, NL)                  !  v := Relax(r)

   ! 
   ! ---------- Additive two-level preconditioning
   ! 
   CASE (NSCARC_TWOLEVEL_ADD, NSCARC_TWOLEVEL_SAMG)

      CALL SCARC_VECTOR_COPY (R, B, 1.0_EB, NL)                   !  Use r^l as right hand side for preconditioner
      DO IL = NL, NLEVEL_MAX-1                                    !  successively restrict to coarser levels up to coarsest
         CALL SCARC_RESTRICTION (B, B, IL, IL+1)                  !  b^{l+1} := Restriction(r^l)
      ENDDO
      CALL SCARC_METHOD_COARSE(NS+2, NS, NLEVEL_MAX)              !  solve A^L * x^L := b^L on coarsest level
#ifdef WITH_SCARC_DEBUG
CALL SCARC_DEBUG_LEVEL (X, 'PRECON-TWOLEVEL-SAMG: X AFTER COARSE ', NLEVEL_MAX)
#endif
      CALL SCARC_VECTOR_COPY (X, Z, 1.0_EB, NLEVEL_MAX)           !  z^L := x^L
#ifdef WITH_SCARC_DEBUG
CALL SCARC_DEBUG_LEVEL (Z, 'PRECON-TWOLEVEL-SAMG: Z AFTER COARSE ', NLEVEL_MAX)
#endif

      DO IL = NLEVEL_MAX-1, NL, -1                                !  successively interpolate to finer levels up to finest
         CALL SCARC_PROLONGATION(Z, Z, IL+1, IL)                  !  z^l := Prolongation(z^{l+1})
#ifdef WITH_SCARC_DEBUG
CALL SCARC_DEBUG_LEVEL (Z, 'PRECON-TWOLEVEL-SAMG: Z AFTER PROL ', IL)
#endif
      ENDDO
      CALL SCARC_VECTOR_COPY (R, V, 1.0_EB, NL)                   !  v^l := r^l
      CALL SCARC_RELAXATION (R, V, NS+1, NP, NL)                  !  v^l := Relax(r^l)
      CALL SCARC_VECTOR_SUM (Z, V, 1.0_EB, 1.0_EB, NL)            !  v^l := z^l + v^l

   ! 
   ! ---------- Multiplicative two-level preconditioning (coarse first, fine second)
   ! 
   CASE (NSCARC_TWOLEVEL_MUL)

      CALL SCARC_VECTOR_COPY (R, B, 1.0_EB, NL)                   !  Use r^l as right hand side for preconditioner

      DO IL = NL, NLEVEL_MAX-1
         CALL SCARC_RESTRICTION (B, B, IL, IL+1)                  !  b^{l+1} := Restriction(r^l)
      ENDDO

      CALL SCARC_METHOD_COARSE(NS+2, NS, NLEVEL_MAX)              !  solve A^L * x^L := b^L on coarsest level
      CALL SCARC_VECTOR_COPY (X, Y, 1.0_EB, NLEVEL_MAX)           !  y^L := x^L

      DO IL = NLEVEL_MAX-1, NL, -1
         CALL SCARC_PROLONGATION (Y, Y, NL+1, NL)                 !  y^l := Prolongation(y^{l+1})
      ENDDO
      CALL SCARC_MATVEC_PRODUCT (Y, Z, NL)                        !  z^l := A^l * y^l

      CALL SCARC_VECTOR_SUM (R, Z, 1.0_EB, -1.0_EB, NL)           !  z^l := r^l - z^l
      CALL SCARC_VECTOR_COPY (Z, V, 1.0_EB, NL)                   !  v^l := z^l
      CALL SCARC_RELAXATION (Z, V, NS+1, NP, NL)                  !  v^l := Relax(z^l)
      CALL SCARC_VECTOR_SUM (Y, V, 1.0_EB, 1.0_EB, NL)            !  v^l := y^l - z^l

   ! 
   ! ---------- Multiplicative two-level preconditioning (fine first, coarse second):
   ! coarse level is one level away from finest one (one coarsening step)
   ! 
   CASE (NSCARC_TWOLEVEL_MUL2)

      CALL SCARC_VECTOR_COPY (R, V, 1.0_EB, NL)                   !  v^l := r^l
      CALL SCARC_RELAXATION (R, V, NS+1, NP, NL)                  !  v^l := Relax(r^l)
      CALL SCARC_MATVEC_PRODUCT (V, Z, NL)                        !  z^l := A^{l} * v^l

      CALL SCARC_VECTOR_SUM (R, Z, 1.0_EB, -1.0_EB, NL)           !  z^l := r^l - z^l

      CALL SCARC_RESTRICTION (Z, B, NL, NL+1)                     !  b^{l+1} := rest(R^{l})
      CALL SCARC_METHOD_COARSE(NS+2, NS, NLEVEL_MAX)              !  x^{l+1} := A^{l+1}^{-1}(b^{l+1})
      CALL SCARC_PROLONGATION (X, Z, NL+1, NL)                    !  v^l := Prolongation(x^{l+1})
      CALL SCARC_VECTOR_SUM (Z, V, 1.0_EB, 1.0_EB, NL)            !  z^l := r^l - z^l

   ! 
   ! ---------- Only coarse grid preconditioner
   ! 
   CASE (NSCARC_TWOLEVEL_COARSE)

      CALL SCARC_VECTOR_COPY (R, B, 1.0_EB, NL)                   !  Use r^l as right hand side for preconditioner
      DO IL = NL, NLEVEL_MAX-1                                    !  successively restrict to coarser levels up to coarsest
         CALL SCARC_RESTRICTION (B, B, IL, IL+1)                  !  b^{l+1} := Restriction(b^l)
      ENDDO

      CALL SCARC_METHOD_COARSE(NS+2, NS, NLEVEL_MAX)              !  solve A^L * x^L := b^L on coarsest level
      CALL SCARC_VECTOR_COPY (X, Y, 1.0_EB, NLEVEL_MAX)           !  y^L := x^L

      DO IL = NLEVEL_MAX-1, NL, -1                                !  successively interpolate to finer levels up to finest
         CALL SCARC_PROLONGATION (Y, Y, NL+1, NL)                 !  y^l := Prolongation(y^{l+1})
      ENDDO
      CALL SCARC_VECTOR_COPY (Y, V, 1.0_EB, NL)                   !  v^l := y^l

END SELECT SELECT_PRECON_TYPE

END SUBROUTINE SCARC_PRECONDITIONER


! ------------------------------------------------------------------------------------------------
! Call requested coarse grid solver (iterative/direct)
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_METHOD_COARSE(NSTACK, NPARENT, NLEVEL)
INTEGER, INTENT(IN) :: NSTACK, NPARENT, NLEVEL

SELECT CASE (TYPE_COARSE)

   CASE (NSCARC_COARSE_ITERATIVE)
      CALL SCARC_METHOD_CG (NSTACK, NPARENT, NSCARC_RHS_DEFECT, NLEVEL)

   CASE (NSCARC_COARSE_DIRECT)
#ifdef WITH_MKL
      IF (STACK(NPARENT)%SOLVER%TYPE_SCOPE(0) == NSCARC_SCOPE_GLOBAL .AND. NMESHES > 1) THEN
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
! Perform geometric multigrid method based on global possion-matrix
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_METHOD_MULTIGRID(NSTACK, NPARENT, NRHS, NLEVEL)
USE SCARC_ITERATION_ENVIRONMENT
INTEGER, INTENT(IN) :: NSTACK, NPARENT, NRHS, NLEVEL
INTEGER :: NS, NP, NL
INTEGER :: NSTATE, ICYCLE
REAL (EB) :: TNOW, TNOW_COARSE

TNOW = CURRENT_TIME()
ITE_MG = 0

! store current and parent stack position and current level
TYPE_MATVEC = NSCARC_MATVEC_GLOBAL
NS = NSTACK
NP = NPARENT
NL = NLEVEL

!
! ---------- Initialization:
!   - Save SETTING (in case that subsequent solvers with different SETTING are called)
!   - Define parameters for current scope (note: NL denotes the finest level)
!   - Initialize solution, right hand side vector
! 
CALL SCARC_SETUP_SOLVER(NS, NP)
CALL SCARC_SETUP_WORKSPACE(NS, NL, NRHS)

! 
! ---------- Compute initial defect:  
!            RESIN := || B - A*X ||
!   - Initialize cycle counts for MG-iteration
!   - Perform initial matrix-vector product on finest level
!   - calculate norm of initial residual on finest level
! 
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

ICYCLE = SCARC_CYCLING_CONTROL(NSCARC_CYCLING_SETUP, NL)

! 
! ---------- Perform multigrid-looping (start each iteration on finest level)
! 
MULTIGRID_LOOP: DO ITE = 1, NIT

   CALL SCARC_INCREASE_ITERATION_COUNTS(ITE)

   NL = NLEVEL_MIN
   ICYCLE = SCARC_CYCLING_CONTROL(NSCARC_CYCLING_RESET, NL)

   CYCLE_LOOP: DO WHILE (ICYCLE /= NSCARC_CYCLING_EXIT)

      !
      ! Presmoothing  (smoothing/restriction till coarsest level is reached)
      ! initial and final residual are passed via vector V by default
      !
      PRESMOOTHING_LOOP: DO WHILE (NL < NLEVEL_MAX)
         CALL SCARC_SMOOTHER (NSCARC_CYCLING_PRESMOOTH, NS+1, NS, NL)         ! D_fine   := Smooth(defect)
         CALL SCARC_RESTRICTION (V, B, NL, NL+1)                              ! B_coarse := Rest(D_fine)
#ifdef WITH_SCARC_DEBUG
CALL SCARC_DEBUG_LEVEL (V, 'MG PRE: V', NL)
CALL SCARC_DEBUG_LEVEL (B, 'MG PRE: B', NL+1)
#endif
         CALL SCARC_VECTOR_CLEAR (X, NL+1)                                    ! use zero initial guess on coarse level
         NL = NL + 1                                                          ! set coarser level
      ENDDO PRESMOOTHING_LOOP

      !
      ! Coarse grid solver
      !
      TNOW_COARSE = CURRENT_TIME()
      CALL SCARC_METHOD_COARSE(NS+2, NS, NLEVEL_MAX)                          ! X_coarse := exact_sol(.)
      CPU(MYID)%COARSE =CPU(MYID)%COARSE+CURRENT_TIME()-TNOW_COARSE

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) '==================> AFTER SCARC_METHOD_COARSE'
CALL SCARC_DEBUG_LEVEL (X, 'MG COA: X', NLEVEL_MAX)
CALL SCARC_DEBUG_LEVEL (B, 'MG COA: B', NLEVEL_MAX)
#endif
      !
      ! Postsmoothing (smoothing/restriction till finest level is reached again)
      !
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

   ! 
   ! Compute norm of new residual on finest level and  leave loop correspondingly
   ! 
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

! 
! ---------- Determine convergence rate and print corresponding information:
! In case of MG as main solver:
!   - Transfer ScaRC solution vector X to FDS pressure vector
!   - Set ghost cell values along external boundaries
!   - Exchange values along internal boundaries (consistency!)
! 
CALL SCARC_CONVERGENCE_RATE(NSTATE, NS, NL)

#ifdef WITH_SCARC_DEBUG
CALL SCARC_DEBUG_LEVEL (X, 'MG-METHOD: FINAL X', NL)
#endif

SELECT CASE (TYPE_SOLVER)
   CASE (NSCARC_SOLVER_MAIN)
      CALL SCARC_UPDATE_MAINCELLS(NLEVEL_MIN)
      CALL SCARC_UPDATE_GHOSTCELLS(NLEVEL_MIN)
#ifdef WITH_SCARC_STANDALONE
      CALL SCARC_PRESSURE_DIFFERENCE(NLEVEL_MIN)
#endif
   CASE (NSCARC_SOLVER_PRECON)
      CALL SCARC_UPDATE_PRECONDITIONER(NLEVEL_MIN)
END SELECT

CALL SCARC_RELEASE_SOLVER(NS, NP)

END SUBROUTINE SCARC_METHOD_MULTIGRID


! ------------------------------------------------------------------------------------------------
! Control multigrid cycling (F/V/W)
! Note: NLEVEL_MIN corresponds to finest level, NLEVEL_MAX to coarsest level
! ------------------------------------------------------------------------------------------------
INTEGER FUNCTION SCARC_CYCLING_CONTROL(NTYPE, NL)
USE SCARC_POINTERS, ONLY: MG
INTEGER, INTENT(IN) :: NTYPE, NL
INTEGER :: NM, IL, ICYCLE

SELECT CASE (NTYPE)

   ! 
   ! initialize cycle counts at beginning of multigrid method
   ! 
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
   
   ! 
   ! reset cycle counts at beginning of each new multigrid iteration
   ! 
   CASE (NSCARC_CYCLING_RESET)
   
      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
         DO IL = NLEVEL_MIN, NLEVEL_MAX
            MG => SCARC(NM)%LEVEL(IL)%MG
            MG%CYCLING(1)=MG%CYCLING(2)
         ENDDO
      ENDDO
      ICYCLE = NSCARC_CYCLING_NEXT
   
   ! 
   ! determine where to proceed with cycling
   ! 
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


! ------------------------------------------------------------------------------------------------
! Perform smoothing
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SMOOTHER(NTYPE, NSTACK, NPARENT, NLEVEL)
USE SCARC_ITERATION_ENVIRONMENT
INTEGER, INTENT(IN) :: NTYPE, NSTACK, NPARENT, NLEVEL
INTEGER :: NSTATE=0, NS, NP, NL
REAL(EB) :: TNOW
LOGICAL :: BMATVEC, BL2NORM, BVERBOSE

!
! ---------- Initialization
!
TNOW = CURRENT_TIME()
TYPE_MATVEC = NSCARC_MATVEC_GLOBAL
NS = NSTACK
NP = NPARENT
NL = NLEVEL

CALL SCARC_SETUP_SOLVER(NS, NP)

!
! Calculate initial defect on l2-norm on level NL (only if BMATVEC and Bl2NORM are set to .TRUE.)
! Because initial vector in MG is set to zero, this defect corresponds to F
!
ITE = 0
BVERBOSE = .TRUE.
IF (BVERBOSE) THEN
   BL2NORM  = .TRUE.
   BMATVEC  = .TRUE.
ELSE
   BL2NORM  = .FALSE.
   BMATVEC  = .FALSE.
ENDIF

#ifdef WITH_SCARC_DEBUG
CALL SCARC_DEBUG_LEVEL (X, 'SMOOTH INIT X', NL)
CALL SCARC_DEBUG_LEVEL (B, 'SMOOTH INIT B', NL)
#endif

IF (BMATVEC) THEN
   CALL SCARC_MATVEC_PRODUCT (X, V, NL)                                  !  v := A*x
   CALL SCARC_VECTOR_SUM (B, V, 1.0_EB, -1.0_EB, NL)                     !  v := b - v    corresponds to   b - A*x
ENDIF

#ifdef WITH_SCARC_DEBUG
CALL SCARC_DEBUG_LEVEL (V, 'SMOOTH INIT V', NL)
#endif

IF (BL2NORM) THEN
   RESIN = SCARC_L2NORM (V, NL)                                          !  resin := ||v||
ELSE
   RESIN = SCARC_RESIDUAL
ENDIF
IF (BVERBOSE) NSTATE = SCARC_CONVERGENCE_STATE(NTYPE, NS, NL)

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'SCARC_SMOOTHER: OMEGA = ', OMEGA
WRITE(MSG%LU_DEBUG,*) 'SCARC_SMOOTHER: RESIN = ', RESIN
#endif

!
! ---------- Smoothing loop - only temporarily
!
!IF (NTYPE == NSCARC_CYCLING_PRESMOOTH) THEN
!   NIT = SCARC_MULTIGRID_PRESMOOTH
!ELSE
!   NIT = SCARC_MULTIGRID_POSTSMOOTH
!ENDIF

SMOOTH_LOOP: DO ITE=1, NIT

   CALL SCARC_INCREASE_ITERATION_COUNTS(ITE)

#ifdef WITH_SCARC_DEBUG
CALL SCARC_DEBUG_LEVEL (V, 'SMOOTH BEFORE RELAX: V', NL)
#endif

#ifdef WITH_MKL
   IF (TYPE_SMOOTH == NSCARC_RELAX_MKL) THEN
      CALL SCARC_VECTOR_COPY(V, Z, 1.0_EB, NL)                          !  use additional auxiliary vector Z
      CALL SCARC_RELAXATION (Z, V, NS, NP, NL)                          !  v := Relax(z)
   ELSE
      CALL SCARC_RELAXATION (V, V, NS, NP, NL)                          !  v := Relax(v)
   ENDIF
#else
   CALL SCARC_RELAXATION (V, V, NS, NP, NL)                             !  v := Relax(v)
#endif

#ifdef WITH_SCARC_DEBUG
CALL SCARC_DEBUG_LEVEL (V, 'SMOOTH AFTER RELAX: V', NL)
CALL SCARC_DEBUG_LEVEL (X, 'SMOOTH BEFORE ADD: X', NL)
#endif
   CALL SCARC_VECTOR_SUM      (V, X, OMEGA, 1.0_EB, NL)                 !  x := omega * v + x

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'OMEGA=',OMEGA
CALL SCARC_DEBUG_LEVEL (X, 'SMOOTH ITE X', NL)
#endif

   CALL SCARC_MATVEC_PRODUCT  (X, V, NL)                                !  v := A*x
#ifdef WITH_SCARC_DEBUG
CALL SCARC_DEBUG_LEVEL (V, 'SMOOTH AFTER MATVEC V', NL)
#endif
   CALL SCARC_VECTOR_SUM      (B, V, 1.0_EB, -1.0_EB, NL)               !  v := b - v

#ifdef WITH_SCARC_DEBUG
CALL SCARC_DEBUG_LEVEL (V, 'SMOOTH AFTER ADD V', NL)
#endif


   IF (BL2NORM) THEN
      RES = SCARC_L2NORM (V, NL)                                        !  res := ||v||
      IF (BVERBOSE) THEN
         NSTATE = SCARC_CONVERGENCE_STATE(NTYPE, NS, NL)
         IF (NSTATE /= NSCARC_STATE_PROCEED) EXIT SMOOTH_LOOP
      ENDIF
   ENDIF

ENDDO SMOOTH_LOOP

CALL SCARC_RELEASE_SOLVER(NS, NP)

CPU(MYID)%SMOOTHER = CPU(MYID)%SMOOTHER + CURRENT_TIME() - TNOW
END SUBROUTINE SCARC_SMOOTHER


! ------------------------------------------------------------------------------------------------
! Setup environement in every solver CALL (i.e. set pointers to used vectors) related to NSTACK
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_SOLVER(NS, NP)
USE SCARC_ITERATION_ENVIRONMENT
USE SCARC_POINTERS, ONLY: SV, SVP
INTEGER, INTENT(IN) :: NS, NP                          ! references to current stack and parent

!
! if not first solver in stack, store last iteration parameters of parent solver NP
!
IF (NP > 0) THEN
   SVP => STACK(NP)%SOLVER
   SVP%ITE   = ITE
   SVP%RES   = RES
   SVP%RESIN = RESIN
   SVP%ERR   = ERR
   SVP%CAPPA = CAPPA
ENDIF

!
! set new environment for solver on stack position NS
!
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
TYPE_KRYLOV         = SV%TYPE_KRYLOV
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
! Reset SETTING of calling CURRENT-routine
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_RELEASE_SOLVER(NS, NP)
USE SCARC_ITERATION_ENVIRONMENT
USE SCARC_POINTERS, ONLY: SV, SVP
INTEGER, INTENT(IN)  :: NS, NP                            ! references to current stack and parent

SV  => STACK(NS)%SOLVER

!
! store convergence information of preceding solver for FDS dump routine
!
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
SV%ERR   = ERR

!
! if not first solver in stack, reset environment of parent (calling) routine
!
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
   TYPE_KRYLOV         = SVP%TYPE_KRYLOV
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
! Set initial solution corresponding to boundary data in BXS, BXF, ...
! ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_WORKSPACE(NS, NL, NRHS)
USE SCARC_POINTERS, ONLY: M, L, F, G, SV, ST, STP, GWC, PRHS, HP
#ifdef WITH_SCARC_MGM
USE SCARC_POINTERS, ONLY: MGM
#endif
#ifdef WITH_SCARC_STANDALONE
USE SCARC_POINTERS, ONLY: P
#endif
INTEGER, INTENT(IN) :: NS, NL, NRHS
INTEGER :: NM, IW, IW1, IW2, IOR0, I, J, K, IC
REAL(EB) :: VAL

SV  => STACK(NS)%SOLVER

SELECT_SOLVER_TYPE: SELECT CASE (SV%TYPE_SOLVER)

   !
   ! --------------- If used as main solver use values from pressure-routine as initialization
   !
   CASE (NSCARC_SOLVER_MAIN)
   
      MAIN_MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   
         !CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G
         CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G
         ST => SCARC(NM)%LEVEL(NL)%STAGE(SV%TYPE_STAGE)
   
         PRHS => M%PRHS
         IF (PREDICTOR) THEN
            HP => M%H
         ELSE
            HP => M%HS
         ENDIF

#ifdef WITH_SCARC_STANDALONE
         IF (PREDICTOR) THEN
            P%H_OLD = P%H_NEW
         ELSE
            P%HS_OLD = P%HS_NEW
         ENDIF
         P%B_OLD = ST%B
#endif
   
         ! get right hand side (PRHS from pres.f90) and initial vector (H or HS from last time step)
         SELECT_RHS_TYPE: SELECT CASE (NRHS)
   
            ! Solve original problem with inhomegeneous boundary conditions
            CASE (NSCARC_RHS_INHOMOGENEOUS)
      
               !$OMP PARALLEL DO PRIVATE(IC) SCHEDULE(STATIC)
               DO IC = 1, G%NC
                  ST%X(IC) = HP(G%ICX(IC), G%ICY(IC), G%ICZ(IC))        ! use last iterate as initial solution
                  ST%B(IC) = PRHS(G%ICX(IC), G%ICY(IC), G%ICZ(IC))      ! get new RHS from surrounding code
               ENDDO                         
               !$OMP END PARALLEL DO
      
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
         
                     ! Dirichlet BC's:
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
         
                     ENDIF IF_DIRICHLET
         
                     ! Neumann BC's:
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
         
                     ENDIF IF_NEUMANN
      
                  ENDDO FACE_INHOMOGENEOUS_LOOP
                  !!$OMP END DO

               ENDDO MAIN_INHOMOGENEOUS_LOOP
               !!$OMP END PARALLEL 
   
#ifdef WITH_SCARC_MGM

            ! Solve problem with homegeneous boundary conditions (MGM only)
            CASE (NSCARC_RHS_HOMOGENEOUS)
   
               ST%B = 0.0_EB                                    ! set RHS to zero
               ST%X = 0.0_EB                                    ! use zero as initial vector
   
               MGM => L%MGM

               !$OMP PARALLEL DO PRIVATE(IW, GWC, I, J, K, IOR0, IC) SCHEDULE(STATIC)
               DO IW = L%N_WALL_CELLS_EXT+1, L%N_WALL_CELLS_EXT + L%N_WALL_CELLS_INT
   
                  GWC => G%WALL(IW)
   
                  I = GWC%IXW
                  J = GWC%IYW
                  K = GWC%IZW
   
                  IF (IS_UNSTRUCTURED .AND. L%IS_SOLID(I, J, K)) CYCLE
   
                  IOR0 = GWC%IOR
                  IC   = G%CELL_NUMBER(I,J,K)
   
                  SELECT CASE (ABS(IOR0))
                     CASE(1)
                        ST%B(IC) = L%DXI * DTI * MGM%US(IW)
                     CASE(2)
                        ST%B(IC) = L%DYI * DTI * MGM%VS(IW)
                     CASE(3)
                        ST%B(IC) = L%DZI * DTI * MGM%WS(IW)
                  END SELECT
   
               ENDDO
               !$OMP END PARALLEL DO
#endif
   
         END SELECT SELECT_RHS_TYPE
 
#ifdef WITH_SCARC_STANDALONE
         P%B_NEW = ST%B
#endif
   
      ENDDO MAIN_MESHES_LOOP
      
      ! In case of a Krylov method clear overlapping parts of auxiliary vectors
      IF (IS_CG.OR.HAS_TWO_LEVELS) THEN
         DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
            L  => SCARC(NM)%LEVEL(NL)
            ST => SCARC(NM)%LEVEL(NL)%STAGE(SV%TYPE_STAGE)
            ST%D = 0.0_EB
            ST%R = 0.0_EB
            ST%V = 0.0_EB
            ST%Y = 0.0_EB
            ST%Z = 0.0_EB
         ENDDO
      ENDIF
   
      ! In case of a multigrid method as main solver clear
      ! overlapping parts of auxiliary vectors and coarse grid solver vectors
      IF (IS_GMG) THEN
         DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
            L  => SCARC(NM)%LEVEL(NL)
            ST => SCARC(NM)%LEVEL(NL)%STAGE(SV%TYPE_STAGE)
            ST%V = 0.0_EB
            ST%Z = 0.0_EB
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
   
   !
   ! --------------- If MG is used as Krylov-preconditioner, vector G of main Krylov is the RHS for MG
   !
   CASE (NSCARC_SOLVER_PRECON)
   
      IF (IS_CG_GMG) THEN
         DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
            L   => SCARC(NM)%LEVEL(NL)
            ST  => SCARC(NM)%LEVEL(NL)%STAGE(SV%TYPE_STAGE)
            STP => SCARC(NM)%LEVEL(NL)%STAGE(NSCARC_STAGE_ONE)
            ST%X = 0.0_EB
            ST%B = STP%R
            ST%V = 0.0_EB
            ST%Z = 0.0_EB
         ENDDO
      ENDIF
   
   !
   ! --------------- If used as coarse grid solver start with zero initialization
   !
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

END SUBROUTINE SCARC_SETUP_WORKSPACE


! ------------------------------------------------------------------------------------------------
! Check if solver converges or diverges and print out residual information
! ------------------------------------------------------------------------------------------------
INTEGER FUNCTION SCARC_CONVERGENCE_STATE(ISM, NS, NL)
USE SCARC_ITERATION_ENVIRONMENT
INTEGER, INTENT(IN) :: NL, NS, ISM
INTEGER :: NSTATE

NSTATE = NSCARC_STATE_PROCEED

SELECT CASE (TYPE_ACCURACY)
   CASE (NSCARC_ACCURACY_RELATIVE)
      IF (RES <= RESIN*EPS .OR. RES <= NSCARC_THRESHOLD_CONVERGENCE) NSTATE = NSCARC_STATE_CONV
   CASE (NSCARC_ACCURACY_ABSOLUTE)
      IF (RES <= EPS .AND. RES <= RESIN) THEN
         IF (ITE == 0) THEN
            NSTATE = NSCARC_STATE_CONV0
         ELSE
            NSTATE = NSCARC_STATE_CONV
         ENDIF
         NIT = 0
      ENDIF
END SELECT
IF (RES > NSCARC_THRESHOLD_DIVGERGENCE) NSTATE = NSCARC_STATE_DIVG

SCARC_CONVERGENCE_STATE = NSTATE

IF (HAS_CSV_DUMP) CALL SCARC_DUMP_CSV(ISM, NS, NL)

#ifdef WITH_SCARC_VERBOSE
IF (TYPE_SOLVER == NSCARC_SOLVER_MAIN) &
   WRITE(MSG%LU_VERBOSE,1100) STACK(NS)%SOLVER%CNAME, NL, ITE, RES
1100 FORMAT (A30,': Level=',i4,': Iteration = ',i4,': Residual =',e12.4)
#endif

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG, 1000) STACK(NS)%SOLVER%CNAME, NL, ITE, RES
1000 FORMAT (A30,': Level=',i4,': Iteration = ',i4,': Residual =',e25.16)
#endif

END FUNCTION SCARC_CONVERGENCE_STATE


! ------------------------------------------------------------------------------------------------
! Compute convergence rate and print out residual information for final loop
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_CONVERGENCE_RATE(NSTATE, NS, NL)
USE SCARC_ITERATION_ENVIRONMENT
INTEGER, INTENT(IN) :: NSTATE, NS, NL

IF (NSTATE == NSCARC_STATE_DIVG) THEN
   ITE   = - 1
   CAPPA = 1.0_EB
ELSE
   IF (NSTATE == NSCARC_STATE_CONV0) THEN
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
         IF (NSTATE == NSCARC_STATE_CONV0) THEN
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

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,2000) STACK(NS)%SOLVER%CNAME, ITE, CAPPA
#endif

#ifdef WITH_SCARC_VERBOSE
IF (TYPE_SOLVER == NSCARC_SOLVER_MAIN) &
   WRITE(MSG%LU_VERBOSE,2000) STACK(NS)%SOLVER%CNAME, ITE, CAPPA
#endif

#if defined (WITH_SCARC_DEBUG) || defined (WITH_SCARC_VERBOSE)
2000 FORMAT (A30,': Iterations: ',i6,':   Convergence Rate =',e14.6,/)
#endif

END SUBROUTINE SCARC_CONVERGENCE_RATE


! ------------------------------------------------------------------------------------------------
! Perform restriction from finer to coarser grid in multigrid method
!    - 'VF' corresponds to vector on fine   grid
!    - 'VC' corresponds to vector on coarse grid
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_RESTRICTION (NVB, NVC, NLF, NLC)
USE SCARC_POINTERS, ONLY: LF, LC, GF, GC, VF, VC, R
INTEGER, INTENT(IN) :: NVB, NVC, NLF, NLC
REAL(EB) :: DSUM
INTEGER :: NM
INTEGER :: IXF, IYF, IZF, ICF(8)=0, ICFB(-2:2,-2:2)=0
INTEGER :: IXC, IYC, IZC, ICC, IC, ICOL

IF (IS_GMG) THEN
!
! ------------------ Twolevel-CG or Geometric multigrid (as main solver or preconditioner) --------------
!
IF (HAS_MULTI_LEVELS) THEN

   DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

      LF => SCARC(NM)%LEVEL(NLF)
      LC => SCARC(NM)%LEVEL(NLC)

      SELECT CASE(TYPE_GRID)
         CASE(NSCARC_GRID_STRUCTURED)
            GC => LC%STRUCTURED
            GF => LF%STRUCTURED
         CASE(NSCARC_GRID_UNSTRUCTURED)
            GC => LC%UNSTRUCTURED
            GF => LF%UNSTRUCTURED
      END SELECT

      VF => SCARC_POINT_TO_VECTOR(NM, NLF, NVB)
      VC => SCARC_POINT_TO_VECTOR(NM, NLC, NVC)

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

ELSE

   IF (NMESHES > 1) CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_VECTOR_PLAIN, NVB, NLF)

   DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

      LF => SCARC(NM)%LEVEL(NLF)
      SELECT CASE(TYPE_GRID)
         CASE(NSCARC_GRID_STRUCTURED)
            GF => LF%STRUCTURED
         CASE(NSCARC_GRID_UNSTRUCTURED)
            GF => LF%UNSTRUCTURED
      END SELECT

      VF => SCARC_POINT_TO_VECTOR(NM, NLF, NVB)
      VC => SCARC_POINT_TO_VECTOR(NM, NLC, NVC)
   
      R => SCARC_POINT_TO_CMATRIX(GF, NSCARC_MATRIX_RESTRICTION)

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'RESTRICTION_AMG: NM, NLF, NLC, N_COARSE:', NM, NLF, NLC, GF%N_COARSE
WRITE(MSG%LU_DEBUG,*) 'VF:'
WRITE(MSG%LU_DEBUG,'(6E12.4)') VF
WRITE(MSG%LU_DEBUG,*) 'VC:'
WRITE(MSG%LU_DEBUG,'(6E12.4)') VC
#endif

      DO IC = 1, GF%N_COARSE
         DSUM = 0.0_EB
         DO ICOL = R%ROW(IC), R%ROW(IC+1)-1                            
            DSUM =  DSUM + R%VAL(ICOL) * VF(R%COLG(ICOL))
#ifdef WITH_SCARC_DEBUG2
WRITE(MSG%LU_DEBUG,'(A,2I6,2E12.4)') 'SCARC_RESTRICTION: IC, ICOL, RVAL, DSUM:', IC, ICOL, R%VAL(ICOL), DSUM
#endif
         ENDDO
         VC(IC) = DSUM
#ifdef WITH_SCARC_DEBUG2
WRITE(MSG%LU_DEBUG,*) '---------------->: IC,  VC(IC):', IC, VC(IC)
#endif
      ENDDO

   ENDDO

ENDIF

END SUBROUTINE SCARC_RESTRICTION


! ------------------------------------------------------------------------------------------------
! Perform prolongation from coarser to finer grid
!    - 'VC' corresponds to coarser grid
!    - 'VF' corresponds to finer   grid
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_PROLONGATION (NVC, NVB, NLC, NLF)
USE SCARC_POINTERS, ONLY: LF, LC, GF, GC, VF, VC, P
INTEGER, INTENT(IN) :: NVC, NVB, NLC, NLF
REAL(EB) :: DSUM
INTEGER :: NM, I
INTEGER :: IXF, IYF, IZF, ICF(8)=0, ICFB(-1:1,-1:1)=0
INTEGER :: IXC, IYC, IZC, ICC, IC, ICOL

IF (IS_GMG) THEN
!
! ------------------ Twolevel CG or Geometric Multigrid 
!
IF (HAS_MULTI_LEVELS) THEN

   DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

      LC => SCARC(NM)%LEVEL(NLC)
      LF => SCARC(NM)%LEVEL(NLF)

      VC => SCARC_POINT_TO_VECTOR(NM, NLC, NVC)
      VF => SCARC_POINT_TO_VECTOR(NM, NLF, NVB)

      SELECT CASE(TYPE_GRID)
         CASE(NSCARC_GRID_STRUCTURED)
            GC => LC%STRUCTURED
            GF => LF%STRUCTURED
         CASE(NSCARC_GRID_UNSTRUCTURED)
            GC => LC%UNSTRUCTURED
            GF => LF%UNSTRUCTURED
      END SELECT

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

ELSE

   DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

      LF => SCARC(NM)%LEVEL(NLF)
      SELECT CASE(TYPE_GRID)
         CASE(NSCARC_GRID_STRUCTURED)
            GF => LF%STRUCTURED
         CASE(NSCARC_GRID_UNSTRUCTURED)
            GF => LF%UNSTRUCTURED
      END SELECT

      VC => SCARC_POINT_TO_VECTOR(NM, NLC, NVC)
      VF => SCARC_POINT_TO_VECTOR(NM, NLF, NVB)

      P  => SCARC_POINT_TO_CMATRIX(GF, NSCARC_MATRIX_PROLONGATION)

      !DO IC = 1, GF%N_FINE
      DO IC = 1, GF%NC
         DSUM = 0.0_EB
         DO ICOL = P%ROW(IC), P%ROW(IC+1)-1                            
            DSUM = DSUM + P%VAL(ICOL) * VC(P%COL(ICOL))
#ifdef WITH_SCARC_DEBUG2
WRITE(MSG%LU_DEBUG,*) 'SCARC_PROLONGATION: IC, ICOL, DSUM:', IC, ICOL, DSUM
#endif
         ENDDO
         VF(IC) = DSUM
#ifdef WITH_SCARC_DEBUG2
WRITE(MSG%LU_DEBUG,*) '----------------->: IC,  VF(IC):', IC, VF(IC)
#endif
      ENDDO

   ENDDO

   IF (NMESHES > 1) CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_VECTOR_PLAIN, NVB, NLF)

ENDIF

END SUBROUTINE SCARC_PROLONGATION


! ------------------------------------------------------------------------------------------------
! Copy final solution from GMG (as preconditioner) to corresponding vector of CG (as main solver)
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_UPDATE_PRECONDITIONER(NL)
INTEGER, INTENT(IN) :: NL
INTEGER :: NM
DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   SCARC(NM)%LEVEL(NL)%STAGE(NSCARC_STAGE_ONE)%V = SCARC(NM)%LEVEL(NL)%STAGE(NSCARC_STAGE_TWO)%X
ENDDO
END SUBROUTINE SCARC_UPDATE_PRECONDITIONER


! ------------------------------------------------------------------------------------------------
! Finalize data
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_UPDATE_MAINCELLS(NL)
USE SCARC_POINTERS, ONLY: M, G, ST, HP
INTEGER, INTENT(IN) :: NL
INTEGER :: NM, IC !, I, J, K

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G
   ST => SCARC(NM)%LEVEL(NL)%STAGE(NSCARC_STAGE_ONE)

   IF (PREDICTOR) THEN
      HP => M%H
   ELSE
      HP => M%HS
   ENDIF

   !DO K = 1, M%KBAR
   !   DO J = 1, M%JBAR
   !      DO I = 1, M%IBAR
   !         IF (IS_UNSTRUCTURED .AND. L%IS_SOLID(I, J, K)) CYCLE
   !         IC = G%CELL_NUMBER(I,J,K)
   !         HP(I, J, K) = ST%X(IC)
   !      ENDDO
   !   ENDDO
   !ENDDO

   !$OMP PARALLEL DO PRIVATE(IC) SCHEDULE(STATIC)
   DO IC = 1, G%NC
      HP (G%ICX(IC), G%ICY(IC), G%ICZ(IC)) = ST%X(IC)
   ENDDO
   !!$OMP END PARALLEL DO 

ENDDO

!IF (MYID == 0) THEN
!WRITE(*,*) 'UPDATE_MAIN:'
!WRITE(*,'(5E25.16)') (((HP(I, J, K), I=5,9), J=1,4), K=1,4)
!ENDIF

END SUBROUTINE SCARC_UPDATE_MAINCELLS


! ------------------------------------------------------------------------------------------------
! Set correct boundary values at external and internal boundaries
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_UPDATE_GHOSTCELLS(NL)
USE SCARC_POINTERS, ONLY: M, L, G, GWC, HP
INTEGER, INTENT(IN) :: NL
INTEGER :: NM, IW, IOR0, IXG, IYG, IZG, IXW, IYW, IZW !, I, J, K

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G

   IF (PREDICTOR) THEN
      HP => M%H
   ELSE
      HP => M%HS
   ENDIF

   !
   ! Compute ghost cell values
   !
   !$OMP PARALLEL DO SHARED(HP, M, L, G) PRIVATE(IW, IXG, IYG, IZG, IXW, IYW, IZW, IOR0, GWC) SCHEDULE(STATIC)
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

   ENDDO WALL_CELLS_LOOP
   !$OMP END PARALLEL DO

#ifdef WITH_SCARC_DEBUG2
   WRITE(MSG%LU_DEBUG,*) 'Updating ghostcells '
   WRITE(MSG%LU_DEBUG,*) 'HP'
   DO IZG=0,M%KBP1
      DO IYG=0,M%JBP1
         WRITE(MSG%LU_DEBUG,'(8E12.4)') (HP(IXG,IYG,IZG), IXG=0,M%IBP1)
      ENDDO
   ENDDO
#endif

ENDDO

! -----------------------------------------------------------------------------------------------
! Perform data exchange to achieve consistency of ghost values along internal boundaries
! Note: this is no longer necessary because MESH_EXCHANGE(5) is used after the call of ScaRC
! -----------------------------------------------------------------------------------------------
CALL SCARC_EXCHANGE(NSCARC_EXCHANGE_PRESSURE, NSCARC_NONE, NL)

!IF (MYID == 0) THEN
!WRITE(*,*) 'UPDATE_MAIN2:'
!WRITE(*,'(5E25.16)') (((HP(I, J, K), I=5,9), J=1,4), K=1,4)
!ENDIF

END SUBROUTINE SCARC_UPDATE_GHOSTCELLS


! ------------------------------------------------------------------------------------------------
! Perform data exchange corresponding to requested exchange type 
! 
! NSCARC_EXCHANGE_BASIC_SIZES     :  exchange initial information about interface sizes
! NSCARC_EXCHANGE_CELL_NEIGHBORS  :  exchange neighboring cell numbers on overlap
! NSCARC_EXCHANGE_CELL_NUMBERS    :  exchange neighboring grid widths on overlap
! NSCARC_EXCHANGE_CELL_SIZES      :  exchange neighboring grid widths on overlap
! NSCARC_EXCHANGE_POISSON_COLS    :  exchange columns of neighboring matrix on overlap
! NSCARC_EXCHANGE_POISSON_COLSG   :  exchange columns of neighboring matrix on overlap
! NSCARC_EXCHANGE_POISSON_DIAGS   :  exchange size of neighboring matrix 
! NSCARC_EXCHANGE_POISSON_SIZES   :  exchange size of neighboring matrix 
! NSCARC_EXCHANGE_POISSON_VALS    :  exchange values of neighboring matrix on overlap
! NSCARC_EXCHANGE_NULLSPACE       :  exchange sum of nullspace entries
! NSCARC_EXCHANGE_PRESSURE        :  exchange vector values along internal boundaries
! NSCARC_EXCHANGE_VECTOR_MEAN     :  exchange vector and build mean values with own data
! NSCARC_EXCHANGE_VECTOR_PLAIN    :  exchange plain vector (just use data from neighbor)
! NSCARC_EXCHANGE_ZONE_NEIGHBORS  :  exchange number of aggregation zones (AMG only)
! NSCARC_EXCHANGE_ZONE_TYPES      :  exchange aggregation zone types (AMG only)
! 
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_EXCHANGE (NTYPE, NPARAM, NL)
INTEGER, INTENT(IN) :: NTYPE, NPARAM, NL
REAL(EB) :: TNOW
INTEGER :: NM, NOM

N_REQ = 0
TYPE_EXCHANGE = NTYPE

!
! ----------------- Receive data from neighbors 
!
RECEIVE_MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID (NM, NL)             

   RECEIVE_OMESHES_LOOP: DO NOM = 1, NMESHES

      RNODE = PROCESS(NM)
      SNODE = PROCESS(NOM)

      IF (RNODE==SNODE .OR.  .NOT.ARE_NEIGHBORS(NM, NOM)) CYCLE RECEIVE_OMESHES_LOOP
      CALL SCARC_POINT_TO_OTHER_GRID(NM, NOM, NL)

      SELECT_EXCHANGE_TYPE: SELECT CASE (NTYPE)

         CASE (NSCARC_EXCHANGE_BASIC_SIZES)
            CALL SCARC_RECV_MESSAGE_INT (NM, NOM, NL, NSCARC_BUFFER_BASIC, 'BASIC SIZES')

         CASE (NSCARC_EXCHANGE_CELL_NEIGHBORS)
            CALL SCARC_RECV_MESSAGE_INT (NM, NOM, NL, NSCARC_BUFFER_LAYER4, 'CELL NEIGHBORS')

         CASE (NSCARC_EXCHANGE_CELL_NUMBERS)
            CALL SCARC_RECV_MESSAGE_INT (NM, NOM, NL, NSCARC_BUFFER_LAYER1, 'CELL NUMBERS')

         CASE (NSCARC_EXCHANGE_CELL_SIZES)
            CALL SCARC_RECV_MESSAGE_REAL (NM, NOM, NL, NSCARC_BUFFER_BASIC, 'CELL SIZES')

         CASE (NSCARC_EXCHANGE_AUX)
            CALL SCARC_RECV_MESSAGE_REAL (NM, NOM, NL, NSCARC_BUFFER_LAYER2, 'AUX')

         CASE (NSCARC_EXCHANGE_NULLSPACE)
            CALL SCARC_RECV_MESSAGE_REAL (NM, NOM, NL, NSCARC_BUFFER_LAYER2, 'NULLSPACE')

         CASE (NSCARC_EXCHANGE_POISSON_COLS)
            CALL SCARC_RECV_MESSAGE_INT (NM, NOM, NL, NSCARC_BUFFER_STENCIL, 'POISSON COLS')

         CASE (NSCARC_EXCHANGE_POISSON_COLSG)
            CALL SCARC_RECV_MESSAGE_INT (NM, NOM, NL, NSCARC_BUFFER_STENCIL, 'POISSON COLSG')

         CASE (NSCARC_EXCHANGE_POISSON_DIAGS)
            CALL SCARC_RECV_MESSAGE_REAL (NM, NOM, NL, NSCARC_BUFFER_LAYER1, 'POISSON DIAGS')

         CASE (NSCARC_EXCHANGE_POISSON_SIZES)
            CALL SCARC_RECV_MESSAGE_INT (NM, NOM, NL, NSCARC_BUFFER_BASIC, 'POISSON SIZES')

         CASE (NSCARC_EXCHANGE_POISSON_VALS)
            CALL SCARC_RECV_MESSAGE_REAL (NM, NOM, NL, NSCARC_BUFFER_STENCIL, 'POISSON VALS')

         CASE (NSCARC_EXCHANGE_PRESSURE)
            CALL SCARC_RECV_MESSAGE_REAL (NM, NOM, NL, NSCARC_BUFFER_LAYER1, 'PRESSURE')

         CASE (NSCARC_EXCHANGE_VECTOR_MEAN)
            CALL SCARC_RECV_MESSAGE_REAL (NM, NOM, NL, NSCARC_BUFFER_LAYER2, 'VECTOR MEAN')

         CASE (NSCARC_EXCHANGE_VECTOR_PLAIN)
            CALL SCARC_RECV_MESSAGE_REAL (NM, NOM, NL, NSCARC_BUFFER_LAYER1, 'VECTOR PLAIN')

         CASE (NSCARC_EXCHANGE_ZONE_NEIGHBORS)
            CALL SCARC_RECV_MESSAGE_INT (NM, NOM, NL, NSCARC_BUFFER_LAYER4, 'ZONE NEIGHBORS')

         CASE (NSCARC_EXCHANGE_ZONE_TYPES)
            CALL SCARC_RECV_MESSAGE_INT (NM, NOM, NL, NSCARC_BUFFER_LAYER2, 'ZONE TYPES')

         CASE DEFAULT
            CALL SCARC_SHUTDOWN(NSCARC_ERROR_EXCHANGE_RECV, SCARC_NONE, TYPE_EXCHANGE)

      END SELECT SELECT_EXCHANGE_TYPE

   ENDDO RECEIVE_OMESHES_LOOP
ENDDO RECEIVE_MESHES_LOOP


! 
! --- Pack data for requested exchange type in corresponding SEND-buffer
! 
TNOW = CURRENT_TIME()

SEND_PACK_MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G

   SEND_PACK_OMESHES_LOOP: DO NOM = 1, NMESHES

      IF (.NOT. ARE_NEIGHBORS(NM, NOM)) CYCLE SEND_PACK_OMESHES_LOOP
      CALL SCARC_POINT_TO_OTHER_GRID(NM, NOM, NL)

      SNODE = PROCESS(NOM)
      RNODE = PROCESS(NM)

      SEND_PACK_OMESHES_SELECT: SELECT CASE (NTYPE)

         CASE (NSCARC_EXCHANGE_BASIC_SIZES)
            CALL SCARC_PACK_BASIC_SIZES
            CALL SCARC_SEND_MESSAGE_INT (NM, NOM, NL, NSCARC_BUFFER_BASIC, 'BASIC SIZES')

         CASE (NSCARC_EXCHANGE_CELL_NEIGHBORS)
            CALL SCARC_PACK_CELL_NEIGHBORS
            CALL SCARC_SEND_MESSAGE_INT (NM, NOM, NL, NSCARC_BUFFER_LAYER4, 'CELL NEIGHBORS')

         CASE (NSCARC_EXCHANGE_CELL_NUMBERS)
            CALL SCARC_PACK_CELL_NUMBERS
            CALL SCARC_SEND_MESSAGE_INT (NM, NOM, NL, NSCARC_BUFFER_LAYER1, 'CELL NUMBERS')

         CASE (NSCARC_EXCHANGE_CELL_SIZES)
            CALL SCARC_PACK_CELL_SIZES
            CALL SCARC_SEND_MESSAGE_REAL (NM, NOM, NL, NSCARC_BUFFER_BASIC, 'CELL SIZES')

         CASE (NSCARC_EXCHANGE_AUX)
            CALL SCARC_PACK_AUX
            CALL SCARC_SEND_MESSAGE_REAL (NM, NOM, NL, NSCARC_BUFFER_LAYER2, 'AUX')

         CASE (NSCARC_EXCHANGE_NULLSPACE)
            CALL SCARC_PACK_NULLSPACE
            CALL SCARC_SEND_MESSAGE_REAL (NM, NOM, NL, NSCARC_BUFFER_LAYER2, 'NULLSPACE')

         CASE (NSCARC_EXCHANGE_POISSON_SIZES)
            CALL SCARC_PACK_POISSON_SIZES(NL)
            CALL SCARC_SEND_MESSAGE_INT (NM, NOM, NL, NSCARC_BUFFER_BASIC, 'POISSON SIZES')

         CASE (NSCARC_EXCHANGE_POISSON_DIAGS)
            CALL SCARC_PACK_POISSON_DIAGS(NSCARC_MATRIX_POISSON)
            CALL SCARC_SEND_MESSAGE_REAL (NM, NOM, NL, NSCARC_BUFFER_LAYER1, 'POISSON DIAGS')

         CASE (NSCARC_EXCHANGE_POISSON_COLS)
            CALL SCARC_PACK_POISSON_COLS(NPARAM)
            CALL SCARC_SEND_MESSAGE_INT (NM, NOM, NL, NSCARC_BUFFER_STENCIL, 'POISSON COLS')

         CASE (NSCARC_EXCHANGE_POISSON_COLSG)
            CALL SCARC_PACK_POISSON_COLSG(NPARAM)
            CALL SCARC_SEND_MESSAGE_INT (NM, NOM, NL, NSCARC_BUFFER_STENCIL, 'POISSON COLSG')

         CASE (NSCARC_EXCHANGE_POISSON_VALS)
            CALL SCARC_PACK_POISSON_VALS(NPARAM, NL)
            CALL SCARC_SEND_MESSAGE_REAL (NM, NOM, NL, NSCARC_BUFFER_STENCIL, 'POISSON VALS')

         CASE (NSCARC_EXCHANGE_PRESSURE)
            CALL SCARC_PACK_PRESSURE(NM)
            CALL SCARC_SEND_MESSAGE_REAL (NM, NOM, NL, NSCARC_BUFFER_LAYER1, 'PRESSURE')

         CASE (NSCARC_EXCHANGE_VECTOR_MEAN)
            CALL SCARC_PACK_VECTOR_MEAN(NM, NL, NPARAM)
            CALL SCARC_SEND_MESSAGE_REAL (NM, NOM, NL, NSCARC_BUFFER_LAYER2, 'VECTOR MEAN')

         CASE (NSCARC_EXCHANGE_VECTOR_PLAIN)
            CALL SCARC_PACK_VECTOR_PLAIN(NM, NL, NPARAM)
            CALL SCARC_SEND_MESSAGE_REAL (NM, NOM, NL, NSCARC_BUFFER_LAYER1, 'VECTOR PLAIN')

         CASE (NSCARC_EXCHANGE_ZONE_NEIGHBORS)
            CALL SCARC_PACK_ZONE_NEIGHBORS
            CALL SCARC_SEND_MESSAGE_INT (NM, NOM, NL, NSCARC_BUFFER_LAYER4, 'ZONE NEIGHBORS')

         CASE (NSCARC_EXCHANGE_ZONE_TYPES)
            CALL SCARC_PACK_ZONE_TYPES
            CALL SCARC_SEND_MESSAGE_INT (NM, NOM, NL, NSCARC_BUFFER_LAYER2, 'ZONE TYPES')

         CASE DEFAULT
            CALL SCARC_SHUTDOWN(NSCARC_ERROR_EXCHANGE_SEND, SCARC_NONE, TYPE_EXCHANGE)

      END SELECT SEND_PACK_OMESHES_SELECT
   ENDDO SEND_PACK_OMESHES_LOOP
ENDDO SEND_PACK_MESHES_LOOP

CPU(MYID)%BUFFER_PACKING = CPU(MYID)%BUFFER_PACKING   + CURRENT_TIME() - TNOW


!
! ----------- Wait for all meshes to have sent and received their data
!

IF (N_MPI_PROCESSES > 1 .AND. N_REQ /= 0) CALL MPI_WAITALL(N_REQ,REQ(1:N_REQ),MPI_STATUSES_IGNORE,IERROR)


!
! ----------- Unpack received data from corresponding RECEIVE-buffers
!
TNOW = CURRENT_TIME()
SEND_UNPACK_MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G

   SEND_UNPACK_OMESHES_LOOP: DO NOM = 1, NMESHES

      SNODE  = PROCESS(NM)
      RNODE  = PROCESS(NOM)

      IF (.NOT. ARE_NEIGHBORS(NM, NOM)) CYCLE SEND_UNPACK_OMESHES_LOOP
      CALL SCARC_POINT_TO_OTHER_GRID(NM, NOM, NL)

      SEND_UNPACK_OMESHES_SELECT: SELECT CASE (NTYPE)

         CASE (NSCARC_EXCHANGE_BASIC_SIZES)
            CALL SCARC_UNPACK_BASIC_SIZES (NM, NOM)

         CASE (NSCARC_EXCHANGE_CELL_NEIGHBORS)
            CALL SCARC_UNPACK_CELL_NEIGHBORS (NM, NOM)

         CASE (NSCARC_EXCHANGE_CELL_NUMBERS)
            CALL SCARC_UNPACK_CELL_NUMBERS (NM, NOM)

         CASE (NSCARC_EXCHANGE_CELL_SIZES)
            CALL SCARC_UNPACK_CELL_SIZES(NM, NOM)

         CASE (NSCARC_EXCHANGE_AUX)
            CALL SCARC_UNPACK_AUX (NM, NOM)

         CASE (NSCARC_EXCHANGE_NULLSPACE)
            CALL SCARC_UNPACK_NULLSPACE (NM, NOM)

         CASE (NSCARC_EXCHANGE_POISSON_SIZES)
            CALL SCARC_UNPACK_POISSON_SIZES (NM, NOM, NL)

         CASE (NSCARC_EXCHANGE_POISSON_DIAGS)
            CALL SCARC_UNPACK_POISSON_DIAGS (NM, NOM)

         CASE (NSCARC_EXCHANGE_POISSON_COLS)
            CALL SCARC_UNPACK_POISSON_COLS (NM, NOM, NPARAM)

         CASE (NSCARC_EXCHANGE_POISSON_COLSG)
            CALL SCARC_UNPACK_POISSON_COLSG (NM, NOM, NPARAM)

         CASE (NSCARC_EXCHANGE_POISSON_VALS)
            CALL SCARC_UNPACK_POISSON_VALS (NM, NOM, NL, NPARAM)

         CASE (NSCARC_EXCHANGE_PRESSURE)
            CALL SCARC_UNPACK_PRESSURE (NM, NOM)

         CASE (NSCARC_EXCHANGE_VECTOR_MEAN)
            CALL SCARC_UNPACK_VECTOR_MEAN (NM, NOM, NL, NPARAM)

         CASE (NSCARC_EXCHANGE_VECTOR_PLAIN)
            CALL SCARC_UNPACK_VECTOR_PLAIN (NM, NOM, NL, NPARAM)

         CASE (NSCARC_EXCHANGE_ZONE_NEIGHBORS)
            CALL SCARC_UNPACK_ZONE_NEIGHBORS (NM, NOM)

         CASE (NSCARC_EXCHANGE_ZONE_TYPES)
            CALL SCARC_UNPACK_ZONE_TYPES (NM, NOM)

         CASE DEFAULT
            CALL SCARC_SHUTDOWN(NSCARC_ERROR_EXCHANGE_SEND, SCARC_NONE, TYPE_EXCHANGE)

      END SELECT SEND_UNPACK_OMESHES_SELECT
   ENDDO SEND_UNPACK_OMESHES_LOOP
ENDDO SEND_UNPACK_MESHES_LOOP

CPU(MYID)%BUFFER_UNPACKING = CPU(MYID)%BUFFER_UNPACKING   + CURRENT_TIME() - TNOW

END SUBROUTINE SCARC_EXCHANGE


! ------------------------------------------------------------------------------------------------
! Receive data of type integer
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_RECV_MESSAGE_INT(NM, NOM, NL, NTYPE, CTEXT)
USE SCARC_POINTERS, ONLY: OS, OG
INTEGER, INTENT(IN) :: NM, NOM, NL, NTYPE
CHARACTER(*), INTENT(IN) :: CTEXT
!INTEGER,  DIMENSION(:), POINTER :: RECV_BUFFER_INT
INTEGER,  POINTER :: RECV_BUFFER_INT
INTEGER :: NLEN
#ifndef WITH_SCARC_VERBOSE2
INTEGER :: IDUMMY
CHARACTER(40) :: CDUMMY
#endif

IF (RNODE == SNODE) RETURN
N_REQ = N_REQ+1

SELECT CASE(NTYPE)
   CASE (NSCARC_BUFFER_BASIC)
      RECV_BUFFER_INT => OS%RECV_BUFFER_INT0(1);  NLEN = NSCARC_MAX_BUFFER0
   CASE (NSCARC_BUFFER_FULL)
      RECV_BUFFER_INT => OS%RECV_BUFFER_INT(1) ;  NLEN = OG%NLEN_BUFFER_FULL
   CASE (NSCARC_BUFFER_LAYER1)
      RECV_BUFFER_INT => OS%RECV_BUFFER_INT(1) ;  NLEN = OG%NLEN_BUFFER_LAYER1
   CASE (NSCARC_BUFFER_LAYER2)
      RECV_BUFFER_INT => OS%RECV_BUFFER_INT(1) ;  NLEN = OG%NLEN_BUFFER_LAYER2
   CASE (NSCARC_BUFFER_LAYER4)
      RECV_BUFFER_INT => OS%RECV_BUFFER_INT(1) ;  NLEN = OG%NLEN_BUFFER_LAYER4
   CASE (NSCARC_BUFFER_STENCIL)
      RECV_BUFFER_INT => OS%RECV_BUFFER_INT(1) ;  NLEN = OG%NLEN_BUFFER_STENCIL
END SELECT

#ifdef WITH_SCARC_VERBOSE2
WRITE(MSG%LU_VERBOSE,1000, ADVANCE = 'NO') CTEXT, NLEN, NM, NOM, NL
#else
IDUMMY = NM; IDUMMY = NOM; IDUMMY = NL                ! prevent compilation warning in case that VERBOSE flag is not set
CDUMMY = CTEXT
#endif

RECV_BUFFER_INT = NSCARC_HUGE_INT
CALL MPI_IRECV(RECV_BUFFER_INT, NLEN, MPI_INTEGER, SNODE, TAG, MPI_COMM_WORLD, REQ(N_REQ),IERROR)

#ifdef WITH_SCARC_VERBOSE2
WRITE(MSG%LU_VERBOSE,*) ' ...  done'
1000 FORMAT('SCARC_RECV_MESSAGE_INT  : Receiving ',A20, ' in length =', I8,' from ',I4, ' to ', I4, ' on level ', I4)
#endif
END SUBROUTINE SCARC_RECV_MESSAGE_INT

! ------------------------------------------------------------------------------------------------
! Receive data of type real
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_RECV_MESSAGE_REAL(NM, NOM, NL, NTYPE, CTEXT)
USE SCARC_POINTERS, ONLY: OS, OG
INTEGER, INTENT(IN) :: NM, NOM, NL, NTYPE
CHARACTER(*), INTENT(IN) :: CTEXT
!REAL(EB), DIMENSION(:), POINTER :: RECV_BUFFER_REAL
REAL(EB), POINTER :: RECV_BUFFER_REAL
INTEGER :: NLEN
#ifndef WITH_SCARC_VERBOSE2
INTEGER :: IDUMMY
CHARACTER(40) :: CDUMMY
#endif

IF (RNODE == SNODE) RETURN
N_REQ = N_REQ+1

SELECT CASE(NTYPE)
   CASE (NSCARC_BUFFER_BASIC)
      RECV_BUFFER_REAL => OS%RECV_BUFFER_REAL0(1);  NLEN = NSCARC_MAX_BUFFER0
   CASE (NSCARC_BUFFER_FULL)
      RECV_BUFFER_REAL => OS%RECV_BUFFER_REAL(1) ;  NLEN = OG%NLEN_BUFFER_FULL
   CASE (NSCARC_BUFFER_LAYER1)
      RECV_BUFFER_REAL => OS%RECV_BUFFER_REAL(1) ;  NLEN = OG%NLEN_BUFFER_LAYER1
   CASE (NSCARC_BUFFER_LAYER2)
      RECV_BUFFER_REAL => OS%RECV_BUFFER_REAL(1) ;  NLEN = OG%NLEN_BUFFER_LAYER2
   CASE (NSCARC_BUFFER_LAYER4)
      RECV_BUFFER_REAL => OS%RECV_BUFFER_REAL(1) ;  NLEN = OG%NLEN_BUFFER_LAYER4
   CASE (NSCARC_BUFFER_STENCIL)
      RECV_BUFFER_REAL => OS%RECV_BUFFER_REAL(1) ;  NLEN = OG%NLEN_BUFFER_STENCIL
END SELECT

#ifdef WITH_SCARC_VERBOSE2
WRITE(MSG%LU_VERBOSE,1000, ADVANCE = 'NO') CTEXT, NLEN, NM, NOM, NL
#else
IDUMMY = NM; IDUMMY = NOM; IDUMMY = NL                ! prevent compilation warning in case that VERBOSE flag is not set
CDUMMY = CTEXT
#endif

RECV_BUFFER_REAL = NSCARC_INIT_UNDEF
CALL MPI_IRECV(RECV_BUFFER_REAL, NLEN, MPI_DOUBLE_PRECISION, SNODE, TAG, MPI_COMM_WORLD, REQ(N_REQ), IERROR)

#ifdef WITH_SCARC_VERBOSE2
WRITE(MSG%LU_VERBOSE,*) ' ...  done'
1000 FORMAT('SCARC_RECV_MESSAGE_REAL : Receiving ',A20, ' in length =', I8,' from ',I4, ' to ', I4, ' on level ', I4)
#endif
END SUBROUTINE SCARC_RECV_MESSAGE_REAL


! ------------------------------------------------------------------------------------------------
! Send data of integer type
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SEND_MESSAGE_INT(NM, NOM, NL, NTYPE, CTEXT)
USE SCARC_POINTERS, ONLY: OS, OG
CHARACTER(*), INTENT(IN) :: CTEXT
INTEGER, INTENT(IN) :: NM, NOM, NL, NTYPE
!INTEGER,  DIMENSION(:), POINTER :: SEND_BUFFER_INT
INTEGER,  POINTER :: SEND_BUFFER_INT
INTEGER :: NLEN
#ifndef WITH_SCARC_VERBOSE2
INTEGER :: IDUMMY
CHARACTER(40) :: CDUMMY
#endif

IF (RNODE == SNODE) RETURN
N_REQ = N_REQ+1

SELECT CASE(NTYPE)
   CASE (NSCARC_BUFFER_BASIC)
      SEND_BUFFER_INT => OS%SEND_BUFFER_INT0(1);  NLEN = NSCARC_MAX_BUFFER0
   CASE (NSCARC_BUFFER_FULL)
      SEND_BUFFER_INT => OS%SEND_BUFFER_INT(1) ;  NLEN = OG%NLEN_BUFFER_FULL
   CASE (NSCARC_BUFFER_LAYER1)
      SEND_BUFFER_INT => OS%SEND_BUFFER_INT(1) ;  NLEN = OG%NLEN_BUFFER_LAYER1
   CASE (NSCARC_BUFFER_LAYER2)
      SEND_BUFFER_INT => OS%SEND_BUFFER_INT(1) ;  NLEN = OG%NLEN_BUFFER_LAYER2
   CASE (NSCARC_BUFFER_LAYER4)
      SEND_BUFFER_INT => OS%SEND_BUFFER_INT(1) ;  NLEN = OG%NLEN_BUFFER_LAYER4
   CASE (NSCARC_BUFFER_STENCIL)
      SEND_BUFFER_INT => OS%SEND_BUFFER_INT(1) ;  NLEN = OG%NLEN_BUFFER_STENCIL
END SELECT

#ifdef WITH_SCARC_VERBOSE2
WRITE(MSG%LU_VERBOSE,1000, ADVANCE = 'NO') CTEXT, NLEN, NM, NOM, NL
#else
IDUMMY = NM; IDUMMY = NOM; IDUMMY = NL                ! prevent compilation warning in case that VERBOSE flag is not set
CDUMMY = CTEXT
#endif

CALL MPI_ISEND(SEND_BUFFER_INT, NLEN, MPI_INTEGER, SNODE, TAG, MPI_COMM_WORLD, REQ(N_REQ),IERROR)
      
#ifdef WITH_SCARC_VERBOSE2
WRITE(MSG%LU_VERBOSE,*) ' ...  done'
1000 FORMAT('SCARC_SEND_MESSAGE_INT  : Sending   ',A20, ' in length =', I8,' from ',I4, ' to ', I4, ' on level ', I4)
#endif
END SUBROUTINE SCARC_SEND_MESSAGE_INT


! ------------------------------------------------------------------------------------------------
! Send data of real type
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SEND_MESSAGE_REAL(NM, NOM, NL, NTYPE, CTEXT)
USE SCARC_POINTERS, ONLY: OS, OG
CHARACTER(*), INTENT(IN) :: CTEXT
INTEGER, INTENT(IN) :: NM, NOM, NL, NTYPE
!REAL(EB), DIMENSION(:), POINTER :: SEND_BUFFER_REAL
REAL(EB), POINTER :: SEND_BUFFER_REAL
INTEGER :: NLEN
#ifndef WITH_SCARC_VERBOSE2
INTEGER :: IDUMMY
CHARACTER(40) :: CDUMMY
#endif

IF (RNODE == SNODE) RETURN
N_REQ = N_REQ+1

SELECT CASE(NTYPE)
   CASE (NSCARC_BUFFER_BASIC)
      SEND_BUFFER_REAL => OS%SEND_BUFFER_REAL0(1);  NLEN = NSCARC_MAX_BUFFER0
   CASE (NSCARC_BUFFER_FULL)
      SEND_BUFFER_REAL => OS%SEND_BUFFER_REAL(1) ;  NLEN = OG%NLEN_BUFFER_FULL
   CASE (NSCARC_BUFFER_LAYER1)
      SEND_BUFFER_REAL => OS%SEND_BUFFER_REAL(1) ;  NLEN = OG%NLEN_BUFFER_LAYER1
   CASE (NSCARC_BUFFER_LAYER2)
      SEND_BUFFER_REAL => OS%SEND_BUFFER_REAL(1) ;  NLEN = OG%NLEN_BUFFER_LAYER2
   CASE (NSCARC_BUFFER_LAYER4)
      SEND_BUFFER_REAL => OS%SEND_BUFFER_REAL(1) ;  NLEN = OG%NLEN_BUFFER_LAYER4
   CASE (NSCARC_BUFFER_STENCIL)
      SEND_BUFFER_REAL => OS%SEND_BUFFER_REAL(1) ;  NLEN = OG%NLEN_BUFFER_STENCIL
END SELECT

#ifdef WITH_SCARC_VERBOSE2
WRITE(MSG%LU_VERBOSE,1000, ADVANCE = 'NO') CTEXT, NLEN, NM, NOM, NL
#else
IDUMMY = NM; IDUMMY = NOM; IDUMMY = NL                ! prevent compilation warning in case that VERBOSE flag is not set
CDUMMY = CTEXT
#endif

CALL MPI_ISEND(SEND_BUFFER_REAL, NLEN, MPI_DOUBLE_PRECISION, SNODE, TAG, MPI_COMM_WORLD, REQ(N_REQ), IERROR)

#ifdef WITH_SCARC_VERBOSE2
WRITE(MSG%LU_VERBOSE,*) ' ...  done'
1000 FORMAT('SCARC_SEND_MESSAGE_REAL : Sending   ',A20, ' in length =', I8,' from ',I4, ' to ', I4, ' on level ', I4)
#endif
END SUBROUTINE SCARC_SEND_MESSAGE_REAL


! ------------------------------------------------------------------------------------------------
! Pack and unpack overlapping parts of cell numbers
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_PACK_CELL_NUMBERS
USE SCARC_POINTERS, ONLY: L, G, OL, OG
INTEGER :: IOR0, ICG, IWG, IXW, IYW, IZW
DO IOR0 = -3, 3
   IF (OL%GHOST_LASTW(IOR0) == 0) CYCLE
   DO ICG = OL%GHOST_FIRSTW(IOR0), OL%GHOST_LASTW(IOR0)
      IWG = OG%ICG_TO_IWG(ICG)
      IXW = G%WALL(IWG)%IXW
      IYW = G%WALL(IWG)%IYW
      IZW = G%WALL(IWG)%IZW
      IF (L%IS_SOLID(IXW, IYW, IZW)) THEN
         OS%SEND_BUFFER_INT(ICG) = -G%CELL_NUMBER(IXW, IYW, IZW)
      ELSE
         OS%SEND_BUFFER_INT(ICG) =  G%CELL_NUMBER(IXW, IYW, IZW)
      ENDIF
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,'(A,7I6)') 'PACK_CELL_NUMBERS: IOR0, ICG, IWG, IXW, IYW, IZW, SEND_BUFFER_INT', &
                               IOR0, ICG, IWG, IXW, IYW, IZW, OS%SEND_BUFFER_INT(ICG)
#endif
   ENDDO
ENDDO
END SUBROUTINE SCARC_PACK_CELL_NUMBERS

! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_UNPACK_CELL_NUMBERS (NM, NOM)
USE SCARC_POINTERS, ONLY: L, G, OL, OG, RECV_BUFFER_INT
INTEGER, INTENT(IN) :: NM, NOM
INTEGER :: LL, IOR0, ICG, ICE, IWG, IXG, IYG, IZG

RECV_BUFFER_INT => SCARC_POINT_TO_BUFFER_INT (NM, NOM, 1)
LL = 1
DO IOR0 = -3, 3
   IF (OL%GHOST_LASTW(IOR0) == 0) CYCLE
   DO ICG = OL%GHOST_FIRSTE(IOR0), OL%GHOST_LASTE(IOR0)
      ICE = OG%ICG_TO_ICE(ICG, 1)
      G%ICE_TO_ICN(ICE) = RECV_BUFFER_INT(LL)
      IWG = OG%ICG_TO_IWG(ICG)
      IXG = G%WALL(IWG)%IXG
      IYG = G%WALL(IWG)%IYG
      IZG = G%WALL(IWG)%IZG
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,'(A,8I6)') 'UNPACK_CELL_NUMBERS: IOR0, ICG, ICE, IWG, IXG, IYG, IZG, RECV_BUFFER_INT', &
                               IOR0, ICG, ICE, IWG, IXG, IYG, IZG, RECV_BUFFER_INT(LL)
#endif
      IF (RECV_BUFFER_INT(ICG) < 0) THEN
         L%IS_SOLID(IXG, IYG, IZG) = .TRUE.
      ELSE
         L%IS_SOLID(IXG, IYG, IZG) = .FALSE.
      ENDIF
      LL = LL + 1
   ENDDO
ENDDO
END SUBROUTINE SCARC_UNPACK_CELL_NUMBERS


! ------------------------------------------------------------------------------------------------
! Pack and unpack cell width information 
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_PACK_CELL_SIZES
USE SCARC_POINTERS, ONLY: L, OS

OS%SEND_BUFFER_REAL0(1) = L%DXL(0)
OS%SEND_BUFFER_REAL0(2) = L%DXL(L%NX)
OS%SEND_BUFFER_REAL0(3) = L%DYL(0)
OS%SEND_BUFFER_REAL0(4) = L%DYL(L%NY)
OS%SEND_BUFFER_REAL0(5) = L%DZL(0)
OS%SEND_BUFFER_REAL0(6) = L%DZL(L%NZ)

END SUBROUTINE SCARC_PACK_CELL_SIZES

! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_UNPACK_CELL_SIZES (NM, NOM)
USE SCARC_POINTERS, ONLY: L, OL, RECV_BUFFER_REAL           
INTEGER, INTENT(IN) :: NM, NOM

RECV_BUFFER_REAL => SCARC_POINT_TO_BUFFER_REAL (NM, NOM, 0)

IF (OL%GHOST_LASTW(-1) /= 0) L%DXL(0)    = 0.5_EB*(RECV_BUFFER_REAL(1) + L%DXL(0))
IF (OL%GHOST_LASTW( 1) /= 0) L%DXL(L%NX) = 0.5_EB*(RECV_BUFFER_REAL(2) + L%DXL(L%NX))
IF (OL%GHOST_LASTW(-2) /= 0) L%DYL(0)    = 0.5_EB*(RECV_BUFFER_REAL(3) + L%DYL(0))
IF (OL%GHOST_LASTW( 2) /= 0) L%DYL(L%NY) = 0.5_EB*(RECV_BUFFER_REAL(4) + L%DYL(L%NY))
IF (OL%GHOST_LASTW(-3) /= 0) L%DZL(0)    = 0.5_EB*(RECV_BUFFER_REAL(5) + L%DZL(0))
IF (OL%GHOST_LASTW( 3) /= 0) L%DZL(L%NZ) = 0.5_EB*(RECV_BUFFER_REAL(6) + L%DZL(L%NZ))

END SUBROUTINE SCARC_UNPACK_CELL_SIZES


! ------------------------------------------------------------------------------------------------
! Pack and unpack initial exchange sizes along interfaces
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_PACK_BASIC_SIZES
USE SCARC_POINTERS, ONLY: OS

OS%SEND_BUFFER_INT0(1)=OG%NCG
OS%SEND_BUFFER_INT0(2)=OG%NZG

END SUBROUTINE SCARC_PACK_BASIC_SIZES

! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_UNPACK_BASIC_SIZES (NM, NOM)
USE SCARC_POINTERS, ONLY: OG, RECV_BUFFER_INT
INTEGER, INTENT(IN) :: NM, NOM

RECV_BUFFER_INT => SCARC_POINT_TO_BUFFER_INT (NM, NOM, 0)

OG%NCG = RECV_BUFFER_INT(1)
OG%NZG = RECV_BUFFER_INT(2)

END SUBROUTINE SCARC_UNPACK_BASIC_SIZES


! ------------------------------------------------------------------------------------------------
! Pack and unpack overlapping parts of specified vector VC 
! Vector VC is numbered via I, J, K values
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_PACK_PRESSURE(NM)
USE SCARC_POINTERS, ONLY: OL, OG, OS
INTEGER, INTENT(IN) :: NM
REAL(EB), DIMENSION(:,:,:), POINTER :: VC
INTEGER :: IOR0, ICG, IWG

IF (PREDICTOR) THEN
   VC => POINT_TO_HVECTOR (NM, NSCARC_VECTOR_H)
ELSE
   VC => POINT_TO_HVECTOR (NM, NSCARC_VECTOR_HS)
ENDIF
DO IOR0 = -3, 3
   IF (OL%GHOST_LASTW(IOR0) == 0) CYCLE
   DO ICG = OL%GHOST_FIRSTW(IOR0), OL%GHOST_LASTW(IOR0)
      IWG = OG%ICG_TO_IWG(ICG)
      OS%SEND_BUFFER_REAL(ICG) = VC(G%WALL(IWG)%IXW, G%WALL(IWG)%IYW, G%WALL(IWG)%IZW)
   ENDDO
ENDDO
END SUBROUTINE SCARC_PACK_PRESSURE

! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_UNPACK_PRESSURE(NM, NOM)
USE SCARC_POINTERS, ONLY: OL, OG, RECV_BUFFER_REAL
INTEGER, INTENT(IN) :: NM, NOM
REAL(EB), DIMENSION(:,:,:), POINTER :: VC
INTEGER :: LL, IOR0, IWG, IXG, IYG, IZG, ICG

RECV_BUFFER_REAL => SCARC_POINT_TO_BUFFER_REAL (NM, NOM, 1)
IF (PREDICTOR) THEN
   VC => POINT_TO_HVECTOR (NM, NSCARC_VECTOR_H)
ELSE
   VC => POINT_TO_HVECTOR (NM, NSCARC_VECTOR_HS)
ENDIF
LL = 1
DO IOR0 = -3, 3
   IF (OL%GHOST_LASTW(IOR0) == 0) CYCLE
   UNPACK_PRESSURE: DO ICG = OL%GHOST_FIRSTE(IOR0), OL%GHOST_LASTE(IOR0)
      IWG = OG%ICG_TO_IWG(ICG)
      IXG=G%WALL(IWG)%IXG
      IYG=G%WALL(IWG)%IYG
      IZG=G%WALL(IWG)%IZG
      VC(IXG, IYG, IZG) = RECV_BUFFER_REAL(LL)
      LL = LL + 1
   ENDDO UNPACK_PRESSURE
ENDDO
END SUBROUTINE SCARC_UNPACK_PRESSURE

! ------------------------------------------------------------------------------------------------
! Pack and unpack overlapping auxiliary vector (currently contained in G%AUX1, to be revised)
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_PACK_AUX
USE SCARC_POINTERS, ONLY: OL, OG, OS, F
INTEGER :: IOR0, ICG, ICW1, ICW2, IWG, IXW, IYW, IZW, LL

LL = 1
OS%SEND_BUFFER_REAL = 0.0_EB
IF (IS_STRUCTURED) THEN

   DO IOR0 = -3, 3
      IF (OL%GHOST_LASTW(IOR0) == 0) CYCLE
      F => L%FACE(IOR0)
      DO ICG = OL%GHOST_FIRSTW(IOR0), OL%GHOST_LASTW(IOR0)

         IWG  = OG%ICG_TO_IWG(ICG)
         IXW  = G%WALL(IWG)%IXW + F%NOFFX
         IYW  = G%WALL(IWG)%IYW + F%NOFFY
         IZW  = G%WALL(IWG)%IZW + F%NOFFZ

         ICW1 = OG%ICG_TO_ICW(ICG, 1)
         ICW2 = G%CELL_NUMBER(IXW, IYW, IZW)

         OS%SEND_BUFFER_REAL(LL)   = G%AUX1(ICW1)
         OS%SEND_BUFFER_REAL(LL+1) = G%AUX1(ICW2)
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,'(A,7I6, 2E12.4)') 'PACK_AUX: IOR0, IWG, IXW, IYW, IZW, ICW1, ICW2: ', &
                                       IOR0, IWG, IXW, IYW, IZW, ICW1, ICW2, G%AUX1(ICW1), G%AUX1(ICW2)
#endif
         LL = LL + 2
      ENDDO
   ENDDO

ELSE

   DO IOR0 = -3, 3

      WRITE(*,*) 'TODO: CHECK!'
      IF (OL%GHOST_LASTW(IOR0) == 0) CYCLE
      F => L%FACE(IOR0)
      DO ICG = OL%GHOST_FIRSTW(IOR0), OL%GHOST_LASTW(IOR0)

         IWG  = OG%ICG_TO_IWG(ICG)
         IXW  = G%WALL(IWG)%IXW + F%NOFFX
         IYW  = G%WALL(IWG)%IYW + F%NOFFY
         IZW  = G%WALL(IWG)%IZW + F%NOFFZ

         ICW1 = OG%ICG_TO_ICW(ICG, 1)
         ICW2 = G%CELL_NUMBER(IXW, IYW, IZW)

         IF (ICW1 > 0) THEN
            OS%SEND_BUFFER_REAL(LL)   = G%AUX1(ICW1)
         ELSE
            OS%SEND_BUFFER_REAL(LL)   = NSCARC_HUGE_REAL
         ENDIF
         IF (ICW2 > 0) THEN
            OS%SEND_BUFFER_REAL(LL+1) = G%AUX1(ICW2)
         ELSE
            OS%SEND_BUFFER_REAL(LL+1) = NSCARC_HUGE_REAL
         ENDIF
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,'(A,7I6, 2E12.4)') 'PACK_AUX: IOR0, IWG, IXW, IYW, IZW, ICW1, ICW2: ', &
                                       IOR0, IWG, IXW, IYW, IZW, ICW1, ICW2, G%AUX1(ICW1), G%AUX1(ICW2)
#endif
         LL = LL + 2
      ENDDO
   ENDDO

ENDIF

END SUBROUTINE SCARC_PACK_AUX

! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_UNPACK_AUX (NM, NOM)
USE SCARC_POINTERS, ONLY: OL, OG, RECV_BUFFER_REAL
INTEGER, INTENT(IN) :: NM, NOM
INTEGER :: IOR0, ICG, ICE1, ICE2, LL

RECV_BUFFER_REAL => SCARC_POINT_TO_BUFFER_REAL (NM, NOM, 1)
LL = 1
DO IOR0 = -3, 3
   IF (OL%GHOST_LASTW(IOR0) == 0) CYCLE
   DO ICG = OL%GHOST_FIRSTE(IOR0), OL%GHOST_LASTE(IOR0)
      ICE1 = OG%ICG_TO_ICE(ICG,1)
      ICE2 = OG%ICG_TO_ICE(ICG,2)
      G%AUX1(ICE1) = RECV_BUFFER_REAL(LL)
      G%AUX1(ICE2) = RECV_BUFFER_REAL(LL+1)
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,'(A,3I6, 2E12.4)') 'UNPACK_AUX: IOR0, ICE1, ICE2, AUX1-1, AUX2-2:', &
                                       IOR0, ICE1, ICE2, G%AUX1(ICE1), G%AUX1(ICE2)
#endif
      LL = LL + 2
   ENDDO
ENDDO

END SUBROUTINE SCARC_UNPACK_AUX


! ------------------------------------------------------------------------------------------------
! Pack and unpack overlapping nullspace vector 
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_PACK_NULLSPACE
USE SCARC_POINTERS, ONLY: OL, OG, OS, F
INTEGER :: IOR0, ICG, ICW1, ICW2, IWG, IXW, IYW, IZW, LL

LL = 1
OS%SEND_BUFFER_REAL = 0.0_EB
IF (IS_STRUCTURED) THEN

   DO IOR0 = -3, 3
      IF (OL%GHOST_LASTW(IOR0) == 0) CYCLE
      F => L%FACE(IOR0)
      DO ICG = OL%GHOST_FIRSTW(IOR0), OL%GHOST_LASTW(IOR0)

         IWG  = OG%ICG_TO_IWG(ICG)
         IXW  = G%WALL(IWG)%IXW + F%NOFFX
         IYW  = G%WALL(IWG)%IYW + F%NOFFY
         IZW  = G%WALL(IWG)%IZW + F%NOFFZ

         ICW1 = OG%ICG_TO_ICW(ICG, 1)
         ICW2 = G%CELL_NUMBER(IXW, IYW, IZW)

         OS%SEND_BUFFER_REAL(LL)   = G%NULLSPACE(ICW1)
         OS%SEND_BUFFER_REAL(LL+1) = G%NULLSPACE(ICW2)
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,'(A,7I6, 2E12.4)') 'PACK_NULLSPACE: IOR0, IWG, IXW, IYW, IZW, ICW1, ICW2: ', &
                                       IOR0, IWG, IXW, IYW, IZW, ICW1, ICW2, G%NULLSPACE(ICW1), G%NULLSPACE(ICW2)
#endif
         LL = LL + 2
      ENDDO
   ENDDO

ELSE

   DO IOR0 = -3, 3

      WRITE(*,*) 'TODO: CHECK!'
      IF (OL%GHOST_LASTW(IOR0) == 0) CYCLE
      F => L%FACE(IOR0)
      DO ICG = OL%GHOST_FIRSTW(IOR0), OL%GHOST_LASTW(IOR0)

         IWG  = OG%ICG_TO_IWG(ICG)
         IXW  = G%WALL(IWG)%IXW + F%NOFFX
         IYW  = G%WALL(IWG)%IYW + F%NOFFY
         IZW  = G%WALL(IWG)%IZW + F%NOFFZ

         ICW1 = OG%ICG_TO_ICW(ICG, 1)
         ICW2 = G%CELL_NUMBER(IXW, IYW, IZW)

         IF (ICW1 > 0) THEN
            OS%SEND_BUFFER_REAL(LL)   = G%NULLSPACE(ICW1)
         ELSE
            OS%SEND_BUFFER_REAL(LL)   = NSCARC_HUGE_REAL
         ENDIF
         IF (ICW2 > 0) THEN
            OS%SEND_BUFFER_REAL(LL+1) = G%NULLSPACE(ICW2)
         ELSE
            OS%SEND_BUFFER_REAL(LL+1) = NSCARC_HUGE_REAL
         ENDIF
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,'(A,7I6, 2E12.4)') 'PACK_NULLSPACE: IOR0, IWG, IXW, IYW, IZW, ICW1, ICW2: ', &
                                       IOR0, IWG, IXW, IYW, IZW, ICW1, ICW2, G%NULLSPACE(ICW1), G%NULLSPACE(ICW2)
#endif
         LL = LL + 2
      ENDDO
   ENDDO

ENDIF

END SUBROUTINE SCARC_PACK_NULLSPACE

! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_UNPACK_NULLSPACE (NM, NOM)
USE SCARC_POINTERS, ONLY: OL, OG, RECV_BUFFER_REAL
INTEGER, INTENT(IN) :: NM, NOM
INTEGER :: IOR0, ICG, ICE1, ICE2, LL

RECV_BUFFER_REAL => SCARC_POINT_TO_BUFFER_REAL (NM, NOM, 1)
LL = 1
DO IOR0 = -3, 3
   IF (OL%GHOST_LASTW(IOR0) == 0) CYCLE
   DO ICG = OL%GHOST_FIRSTE(IOR0), OL%GHOST_LASTE(IOR0)
      ICE1 = OG%ICG_TO_ICE(ICG,1)
      ICE2 = OG%ICG_TO_ICE(ICG,2)
      G%NULLSPACE(ICE1) = RECV_BUFFER_REAL(LL)
      G%NULLSPACE(ICE2) = RECV_BUFFER_REAL(LL+1)
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,'(A,3I6, 2E12.4)') 'UNPACK_NULLSPACE: IOR0, ICE1, ICE2, NULLSPACE-1, AUX2-2:', &
                                       IOR0, ICE1, ICE2, G%NULLSPACE(ICE1), G%NULLSPACE(ICE2)
#endif
      LL = LL + 2
   ENDDO
ENDDO

END SUBROUTINE SCARC_UNPACK_NULLSPACE



! ------------------------------------------------------------------------------------------------
! Pack and unpack overlapping parts of specified vector VC (numbered via IC values)
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_PACK_VECTOR_PLAIN(NM, NL, NV)
USE SCARC_POINTERS, ONLY: OL, OG, OS
INTEGER, INTENT(IN) :: NM, NL, NV
REAL(EB), DIMENSION(:), POINTER :: VC
INTEGER :: IOR0, ICG, ICW

VC => SCARC_POINT_TO_VECTOR(NM, NL, NV)
#ifdef WITH_SCARC_DEBUG2
WRITE(MSG%LU_DEBUG,*) '====================================================='
WRITE(MSG%LU_DEBUG,*) ' PACK_VECTOR_PLAIN : MESH ', NM, ' LEVEL ', NL, ' LENGTH ', OG%NCG
WRITE(MSG%LU_DEBUG,'(6E14.6)') VC
#endif
IF (IS_STRUCTURED) THEN
   DO IOR0 = -3, 3
      IF (OL%GHOST_LASTW(IOR0) == 0) CYCLE
      DO ICG = OL%GHOST_FIRSTW(IOR0), OL%GHOST_LASTW(IOR0)
         ICW = OG%ICG_TO_ICW(ICG, 1)
         OS%SEND_BUFFER_REAL(ICG) = VC(ICW)
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,'(A,3I6,E14.6)') 'PACK_VECTOR_PLAIN: IOR0, ICG, ICW, VC:', IOR0, ICG, ICW, VC(ICW)
#endif
      ENDDO
   ENDDO
ELSE
   DO IOR0 = -3, 3
      IF (OL%GHOST_LASTW(IOR0) == 0) CYCLE
      DO ICG = OL%GHOST_FIRSTE(IOR0), OL%GHOST_LASTE(IOR0)
         ICW = OG%ICG_TO_ICW(ICG, 1)
         IF (ICW > 0) THEN
            OS%SEND_BUFFER_REAL(ICG) = VC(ICW)
         ELSE
            OS%SEND_BUFFER_REAL(ICG) = NSCARC_HUGE_REAL
         ENDIF
      ENDDO
   ENDDO
ENDIF

END SUBROUTINE SCARC_PACK_VECTOR_PLAIN

! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_UNPACK_VECTOR_PLAIN(NM, NOM, NL, NVECTOR)
USE SCARC_POINTERS, ONLY: OL, OG, RECV_BUFFER_REAL
INTEGER, INTENT(IN) :: NM, NOM, NL, NVECTOR
REAL(EB), DIMENSION(:), POINTER :: VC
INTEGER :: IOR0, ICG, LL


#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) '====================================================='
WRITE(MSG%LU_DEBUG,*) ' UNPACK_VECTOR_PLAIN : MESH ', NM, ' LEVEL ', NL, ' NEIGHBOR ', NOM, ' LENGTH ', OG%NCG
#endif

RECV_BUFFER_REAL => SCARC_POINT_TO_BUFFER_REAL (NM, NOM, 1)
VC => SCARC_POINT_TO_VECTOR(NM, NL, NVECTOR)

#ifdef WITH_SCARC_DEBUG2
WRITE(MSG%LU_DEBUG,*) 'UNPACK_VECTOR_PLAIN: BEFORE: NM, NOM, NVECTOR:', NM, NOM, NVECTOR
WRITE(MSG%LU_DEBUG,'(6E14.6)') VC
#endif

LL = 1
DO IOR0 = -3, 3
   IF (OL%GHOST_LASTW(IOR0) == 0) CYCLE
   DO ICG = OL%GHOST_FIRSTE(IOR0), OL%GHOST_LASTE(IOR0)
      VC(OG%ICG_TO_ICE(ICG, 1)) = RECV_BUFFER_REAL(LL)
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,'(A,3I6,E14.6)') 'UNPACK_VECTOR_PLAIN: IOR0, ICG, VC:', &
                                     IOR0, ICG, OG%ICG_TO_ICE(ICG,1), RECV_BUFFER_REAL(LL)
#endif
      LL = LL + 1
   ENDDO
ENDDO

#ifdef WITH_SCARC_DEBUG2
WRITE(MSG%LU_DEBUG,*) 'UNPACK_VECTOR_PLAIN: AFTER: NM, NOM, NVECTOR:', NM, NOM, NVECTOR
WRITE(MSG%LU_DEBUG,'(6E14.6)') VC
#endif

END SUBROUTINE SCARC_UNPACK_VECTOR_PLAIN


! ------------------------------------------------------------------------------------------------
! Pack and unpack overlapping and internal parts of specified vector
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_PACK_VECTOR_MEAN(NM, NL, NVECTOR)
USE SCARC_POINTERS, ONLY: OL, OG, OS
INTEGER, INTENT(IN) :: NM, NL, NVECTOR
REAL(EB), DIMENSION(:), POINTER :: VC
INTEGER :: IOR0, ICG, ICW, ICE, LL

VC => SCARC_POINT_TO_VECTOR(NM, NL, NVECTOR)
SELECT CASE (TYPE_GRID)
   CASE (NSCARC_GRID_STRUCTURED)
      LL = 1
      DO IOR0 = -3, 3
         IF (OL%GHOST_LASTW(IOR0) == 0) CYCLE
         DO ICG = OL%GHOST_FIRSTW(IOR0), OL%GHOST_LASTW(IOR0)
            ICW = OG%ICG_TO_ICW(ICG, 1)
            OS%SEND_BUFFER_REAL(LL) = VC(ICW)
            LL = LL + 1
         ENDDO
         DO ICG = OL%GHOST_FIRSTW(IOR0), OL%GHOST_LASTW(IOR0)
            ICE = OG%ICG_TO_ICE(ICG, 1)
            OS%SEND_BUFFER_REAL(LL) = VC(ICE)
            LL = LL + 1
         ENDDO
      ENDDO
   CASE (NSCARC_GRID_UNSTRUCTURED)
      DO IOR0 = -3, 3
         IF (OL%GHOST_LASTW(IOR0) == 0) CYCLE
         DO ICG = OL%GHOST_FIRSTW(IOR0), OL%GHOST_LASTW(IOR0)
            ICE = OG%ICG_TO_ICE(ICG, 1)
            IF (ICW > 0) THEN
               OS%SEND_BUFFER_REAL(ICG) = VC(ICE)
            ELSE
               OS%SEND_BUFFER_REAL(ICG) = NSCARC_HUGE_REAL
            ENDIF
         ENDDO
      ENDDO
END SELECT

END SUBROUTINE SCARC_PACK_VECTOR_MEAN

! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_UNPACK_VECTOR_MEAN(NM, NOM, NL, NVECTOR)
USE SCARC_POINTERS, ONLY: OL, OG, RECV_BUFFER_REAL
INTEGER, INTENT(IN) :: NM, NOM, NL, NVECTOR
REAL(EB), DIMENSION(:), POINTER :: VC
INTEGER :: IOR0, ICG, LL

RECV_BUFFER_REAL => SCARC_POINT_TO_BUFFER_REAL (NM, NOM, 1)
VC => SCARC_POINT_TO_VECTOR(NM, NL, NVECTOR)
LL = 1
DO IOR0 = -3, 3
   IF (OL%GHOST_LASTW(IOR0) == 0) CYCLE
   DO ICG = OL%GHOST_FIRSTE(IOR0), OL%GHOST_LASTE(IOR0)
      !VC(OG%ICG_TO_ICE(ICG, 1)) = 0.5_EB * (VC(OG%ICG_TO_ICE(ICG, 1)) + RECV_BUFFER_REAL(LL))
      VC(OG%ICG_TO_ICE(ICG, 1)) = 1.0_EB/3.0_EB * VC(OG%ICG_TO_ICE(ICG, 1)) + 2.0_EB/3.0_EB * RECV_BUFFER_REAL(LL)
      LL = LL + 1
   ENDDO
   DO ICG = OL%GHOST_FIRSTE(IOR0), OL%GHOST_LASTE(IOR0)
      !VC(OG%ICG_TO_ICW(ICG, 1)) = 0.5_EB * (VC(OG%ICG_TO_ICW(ICG, 1)) + RECV_BUFFER_REAL(LL))
      VC(OG%ICG_TO_ICW(ICG, 1)) = 2.0_EB/3.0_EB * VC(OG%ICG_TO_ICW(ICG, 1)) + 1.0_EB/3.0_EB * RECV_BUFFER_REAL(LL)
      LL = LL + 1
   ENDDO
ENDDO

END SUBROUTINE SCARC_UNPACK_VECTOR_MEAN


! ------------------------------------------------------------------------------------------------
! Pack and unpack overlapping information about matrix columns (compact storage technique only)
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_PACK_POISSON_COLS(NMATRIX)                
USE SCARC_POINTERS, ONLY: OS, OL, OG, A
INTEGER, INTENT(IN) :: NMATRIX
INTEGER :: IOR0, ICG, ICW, LL, ICOL

A => SCARC_POINT_TO_CMATRIX(G, NMATRIX)

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) '=============== EXCHANGE_COLS: TYPE_EXCHANGE=', TYPE_EXCHANGE, NMATRIX
WRITE(MSG%LU_DEBUG,*) 'SIZE(A%ROW)', SIZE(A%ROW)
WRITE(MSG%LU_DEBUG,*) 'SIZE(A%COL)', SIZE(A%COL)
WRITE(MSG%LU_DEBUG,*) 'SIZE(A%VAL)', SIZE(A%VAL)
#endif
         
LL = 1
OS%SEND_BUFFER_INT = NSCARC_INIT_UNDEF 

SELECT CASE (TYPE_GRID) 

   CASE (NSCARC_GRID_STRUCTURED)

      DO IOR0 = -3, 3
         IF (OL%GHOST_LASTW(IOR0) == 0) CYCLE
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) '============= PACK_POISSON_COLS: SEND: FIRST, LAST, IOR0:', &
                       OL%GHOST_FIRSTW(IOR0), OL%GHOST_LASTW(IOR0), IOR0
#endif
         DO ICG = OL%GHOST_FIRSTW(IOR0), OL%GHOST_LASTW(IOR0)
            ICW = OG%ICG_TO_ICW(ICG, 1)
            
            ICOL = A%ROW(ICW)
            OS%SEND_BUFFER_INT(LL) = -A%COL(ICOL)

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,'(A,5I5)') 'PACK_POISSON_COLS: SEND: IOR0, ICG, ICW, ICOL, VAL:', &
                               IOR0, ICG, ICW, A%ROW(ICW), OS%SEND_BUFFER_INT(LL)
#endif
            LL = LL + 1                                ! such that receiver can identify start of column
            DO ICOL = A%ROW(ICW)+1, A%ROW(ICW+1)-1
               OS%SEND_BUFFER_INT(LL) = A%COL(ICOL)
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,'(A,5I5)') 'PACK_POISSON_COLS: SEND: IOR0, ICG, ICW, ICOL, VAL:', &
                               IOR0, ICG, ICW, ICOL, OS%SEND_BUFFER_INT(LL)
#endif
               LL = LL + 1
            ENDDO
         ENDDO
      ENDDO

   CASE (NSCARC_GRID_UNSTRUCTURED)

      DO IOR0 = -3, 3
         IF (OL%GHOST_LASTW(IOR0) == 0) CYCLE

         DO ICG= OL%GHOST_FIRSTW(IOR0), OL%GHOST_LASTW(IOR0)
            ICW = OG%ICG_TO_ICW(ICG, 1)
            IF (ICW > 0) THEN
               ICOL = A%ROW(ICW)
               OS%SEND_BUFFER_INT(LL) = -A%COL(ICOL)           ! send first element with negative sign
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,'(A,5I5)') 'PACK_POISSON_COLS: SEND: IOR0, ICG, ICW, ICOL, VAL:', &
                                  IOR0, ICG, ICW, A%ROW(ICW), OS%SEND_BUFFER_INT(LL)
#endif
               LL = LL + 1                                ! such that receiver can identify start of column
               DO ICOL = A%ROW(ICW)+1, A%ROW(ICW+1)-1
                  OS%SEND_BUFFER_INT(LL) = A%COL(ICOL)
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,'(A,5I5)') 'PACK_POISSON_COLS: SEND: IOR0, ICG, ICW, ICOL, VAL:', &
                               IOR0, ICG, ICW, ICOL, OS%SEND_BUFFER_INT(LL)
#endif
                  LL = LL + 1
               ENDDO
            ELSE
               OS%SEND_BUFFER_INT(LL) = NSCARC_HUGE_INT
               LL = LL + 1
            ENDIF
         ENDDO
      ENDDO

END SELECT
#ifdef WITH_SCARC_DEBUG2
WRITE(MSG%LU_DEBUG,*) 'SEND_BUFFER:'
WRITE(MSG%LU_DEBUG,'(10I6)') OS%SEND_BUFFER_INT
#endif
END SUBROUTINE SCARC_PACK_POISSON_COLS

! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_UNPACK_POISSON_COLS(NM, NOM, NMATRIX)
USE SCARC_POINTERS, ONLY: OL, OG, OA, RECV_BUFFER_INT
INTEGER, INTENT(IN) :: NM, NOM, NMATRIX
INTEGER :: IOR0, ICG, LL, ICE, ICP

#ifdef WITH_SCARC_DEBUG2
WRITE(MSG%LU_DEBUG,*) 'UNPACK_POISSON_COLS: OA%N_ROW, OA%N_VAL=',OA%N_ROW, OA%N_VAL, NM, NOM, NMATRIX
WRITE(MSG%LU_DEBUG,*) 'UNPACK_POISSON_COLS: SIZE(OA%COL)=', SIZE(OA%COL)
WRITE(MSG%LU_DEBUG,*) 'UNPACK_POISSON_COLS: SIZE(OA%ROW)=', SIZE(OA%ROW)
WRITE(MSG%LU_DEBUG,*) 'RECV_BUFFER_INT:'
WRITE(MSG%LU_DEBUG,'(10I6)') OS%RECV_BUFFER_INT
#endif

RECV_BUFFER_INT => SCARC_POINT_TO_BUFFER_INT (NM, NOM, 1)
OA => SCARC_POINT_TO_OTHER_CMATRIX(OG, NMATRIX)

LL = 1                                 
ICP = 1
DO IOR0 = -3, 3
   IF (OL%GHOST_LASTW(IOR0) == 0) CYCLE
   OA%ROW(ICP) = 1
   DO ICG = OL%GHOST_FIRSTE(IOR0), OL%GHOST_LASTE(IOR0)
      ICE = OG%ICG_TO_ICE(ICG, 1)
      OA%COL(LL) = ABS(RECV_BUFFER_INT(LL))
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,'(A, 4I4)') 'UNPACK_POISSON_COLS: A: RECV: ICG, LL, RECV_BUFFER_INT, OA%ROW(ICG), OA%COL(LL):', &
                                ICG, LL, RECV_BUFFER_INT(LL), OA%COL(LL)
#endif
      DO WHILE (RECV_BUFFER_INT(LL+1) >= 0)
         LL = LL + 1
         OA%COL(LL) = ABS(RECV_BUFFER_INT(LL))
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,'(A, 4I4)') 'UNPACK_POISSON_COLS: B: RECV: ICG, LL, RECV_BUFFER_INT, OA%ROW(ICG), OA%COL(LL):', &
                                ICG, LL, RECV_BUFFER_INT(LL), OA%COL(LL)
#endif
      ENDDO
      LL = LL + 1
      ICP = ICP + 1
      OA%ROW(ICP) = LL
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) '----------> UNPACK_POISSON_COLS: OA%ROW(',ICP,')=',OA%ROW(ICP)
#endif
   ENDDO
   OA%N_ROW = ICP  
   OA%N_VAL = LL - 1
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,'(A, 4I4)') 'UNPACK_POISSON_COLS: OA%N_ROW=',OA%N_ROW
WRITE(MSG%LU_DEBUG,'(A, 4I4)') 'UNPACK_POISSON_COLS: OA%N_VAL=',OA%N_VAL
#endif
ENDDO

END SUBROUTINE SCARC_UNPACK_POISSON_COLS


! ------------------------------------------------------------------------------------------------
! Pack and unpack overlapping information about matrix columns (compact storage technique only)
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_PACK_POISSON_COLSG(NMATRIX)                
USE SCARC_POINTERS, ONLY: OS, OL, OG, A
INTEGER, INTENT(IN) :: NMATRIX
INTEGER :: IOR0, ICG, ICW, LL, ICOL

A => SCARC_POINT_TO_CMATRIX(G, NMATRIX)

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) '========= EXCHANGE_COLSG: TYPE_EXCHANGE=', TYPE_EXCHANGE
WRITE(MSG%LU_DEBUG,*) 'SIZE(A%ROW)', SIZE(A%ROW)
WRITE(MSG%LU_DEBUG,*) 'SIZE(A%COL)', SIZE(A%COL)
WRITE(MSG%LU_DEBUG,*) 'SIZE(A%VAL)', SIZE(A%VAL)
#endif
         
LL = 1
OS%SEND_BUFFER_INT = NSCARC_INIT_UNDEF

SELECT CASE (TYPE_GRID) 

   CASE (NSCARC_GRID_STRUCTURED)

      DO IOR0 = -3, 3
         IF (OL%GHOST_LASTW(IOR0) == 0) CYCLE
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) '======= PACK_POISSON_COLSG: SEND: GHOST_FIRSTW, NCG:', OL%GHOST_FIRSTW(IOR0), OL%GHOST_LASTW(IOR0)
#endif
         DO ICG = OL%GHOST_FIRSTW(IOR0), OL%GHOST_LASTW(IOR0)
            ICW = OG%ICG_TO_ICW(ICG, 1)
            
            ICOL = A%ROW(ICW)
            !OS%SEND_BUFFER_INT(LL) = -G%LOCAL_TO_GLOBAL(A%COL(ICOL))           ! send first element with negative sign
            OS%SEND_BUFFER_INT(LL) = -A%COLG(ICOL)          ! send first element with negative sign

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,'(A,5I5)') 'PACK_POISSON_COLSG: SEND: IOR0, ICG, ICW, ICOL, VAL:', &
                               IOR0, ICG, ICW, A%ROW(ICW), OS%SEND_BUFFER_INT(LL)
#endif
            LL = LL + 1                                ! such that receiver can identify start of column
            DO ICOL = A%ROW(ICW)+1, A%ROW(ICW+1)-1
               !OS%SEND_BUFFER_INT(LL) = G%LOCAL_TO_GLOBAL(A%COL(ICOL))           ! send first element with negative sign
               OS%SEND_BUFFER_INT(LL) = A%COLG(ICOL)           ! send first element with negative sign
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,'(A,5I5)') 'PACK_POISSON_COLSG: SEND: IOR0, ICG, ICW, ICOL, VAL:', &
                               IOR0, ICG, ICW, ICOL, OS%SEND_BUFFER_INT(LL)
#endif
               LL = LL + 1
            ENDDO
         ENDDO
      ENDDO

   CASE (NSCARC_GRID_UNSTRUCTURED)

      DO IOR0 = -3, 3
         IF (OL%GHOST_LASTW(IOR0) == 0) CYCLE

         DO ICG= OL%GHOST_FIRSTW(IOR0), OL%GHOST_LASTW(IOR0)
            ICW = OG%ICG_TO_ICW(ICG, 1)
            IF (ICW > 0) THEN
               ICOL = A%ROW(ICW)
               !OS%SEND_BUFFER_INT(LL) = -G%LOCAL_TO_GLOBAL(A%COL(ICOL))           ! send first element with negative sign
               OS%SEND_BUFFER_INT(LL) = -A%COLG(ICOL)           ! send first element with negative sign
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,'(A,5I5)') 'PACK_POISSON_COLSG: SEND: IOR0, ICG, ICW, ICOL, VAL:', &
                                  IOR0, ICG, ICW, A%ROW(ICW), OS%SEND_BUFFER_INT(LL)
#endif
               LL = LL + 1                                ! such that receiver can identify start of column
               DO ICOL = A%ROW(ICW)+1, A%ROW(ICW+1)-1
                  !OS%SEND_BUFFER_INT(LL) = G%LOCAL_TO_GLOBAL(A%COL(ICOL))           ! send first element with negative sign
                  OS%SEND_BUFFER_INT(LL) = A%COL(ICOL)           ! send first element with negative sign
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,'(A,5I5)') 'PACK_POISSON_COLSG: SEND: IOR0, ICG, ICW, ICOL, VAL:', &
                               IOR0, ICG, ICW, ICOL, OS%SEND_BUFFER_INT(LL)
#endif
                  LL = LL + 1
               ENDDO
            ELSE
               OS%SEND_BUFFER_INT(LL) = NSCARC_HUGE_INT
               LL = LL + 1
            ENDIF
         ENDDO
      ENDDO

END SELECT
END SUBROUTINE SCARC_PACK_POISSON_COLSG

! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_UNPACK_POISSON_COLSG(NM, NOM, NMATRIX)
USE SCARC_POINTERS, ONLY: OL, OG, OA, RECV_BUFFER_INT
INTEGER, INTENT(IN) :: NM, NOM, NMATRIX
INTEGER :: IOR0, ICG, LL, ICE, ICP

RECV_BUFFER_INT => SCARC_POINT_TO_BUFFER_INT (NM, NOM, 1)
OA => SCARC_POINT_TO_OTHER_CMATRIX(OG, NMATRIX)

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'UNPACK_POISSON_COLSG: OA%N_ROW, OA%N_VAL=',OA%N_ROW, OA%N_VAL
#endif
LL = 1                                 
ICP = 1
DO IOR0 = -3, 3
   IF (OL%GHOST_LASTW(IOR0) == 0) CYCLE
   OA%ROW(ICP) = 1
   DO ICG = OL%GHOST_FIRSTE(IOR0), OL%GHOST_LASTE(IOR0)
      ICE = OG%ICG_TO_ICE(ICG, 1)
      OA%COLG(LL) = ABS(RECV_BUFFER_INT(LL))
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,'(A, 4I4)') 'UNPACK_POISSON_COLSG: A: RECV: ICG, LL, RECV_BUFFER_INT, OA%ROW(ICG), OA%COLG(LL):', &
                                ICG, LL, RECV_BUFFER_INT(LL), OA%COLG(LL)
#endif
      DO WHILE (RECV_BUFFER_INT(LL+1) > 0)
         LL = LL + 1
         OA%COLG(LL) = ABS(RECV_BUFFER_INT(LL))
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,'(A, 4I4)') 'UNPACK_POISSON_COLSG: B: RECV: ICG, LL, RECV_BUFFER_INT, OA%ROW(ICG), OA%COLG(LL):', &
                                ICG, LL, RECV_BUFFER_INT(LL), OA%COLG(LL)
#endif
      ENDDO
      LL = LL + 1
      ICP = ICP + 1
      OA%ROW(ICP) = LL
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) '----------> UNPACK_POISSON_COLSG: OA%ROW(',ICP,')=',OA%ROW(ICP)
#endif
   ENDDO
   OA%N_ROW = ICP  
   OA%N_VAL = LL - 1
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,'(A, 4I4)') 'UNPACK_POISSON_COLSG: OA%N_ROW=',OA%N_ROW
WRITE(MSG%LU_DEBUG,'(A, 4I4)') 'UNPACK_POISSON_COLSG: OA%N_VAL=',OA%N_VAL
#endif
ENDDO

END SUBROUTINE SCARC_UNPACK_POISSON_COLSG


! ------------------------------------------------------------------------------------------------
! Pack and unpack overlapping information about matrix values (compact storage technique only)
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_PACK_POISSON_VALS(NMATRIX, NL)
USE SCARC_POINTERS, ONLY:  A, AB, OS, OL, OG
INTEGER, INTENT(IN) :: NMATRIX, NL
INTEGER :: IOR0, ICG, ICW, LL, ID, ICOL

A => SCARC_POINT_TO_CMATRIX(G, NMATRIX)

LL = 1
SELECT CASE (SCARC_MATRIX_LEVEL(NL))

   ! Bandwise matrix on level NL
   CASE (NSCARC_MATRIX_BANDWISE)                
      DO IOR0 = -3, 3
         IF (OL%GHOST_LASTW(IOR0) == 0) CYCLE
         DO ICG = OL%GHOST_FIRSTW(IOR0), OL%GHOST_LASTW(IOR0)
            ICW = OG%ICG_TO_ICW(ICG, 1)
            DO ID = 1, AB%N_STENCIL
               OS%SEND_BUFFER_REAL(LL) = AB%VAL(ICW, ID)
               LL = LL + 1
            ENDDO
         ENDDO
      ENDDO

   ! Compact matrix on level NL
   CASE (NSCARC_MATRIX_COMPACT)

      SELECT CASE (TYPE_GRID)

         CASE (NSCARC_GRID_STRUCTURED)

            DO IOR0 = -3, 3
               IF (OL%GHOST_LASTW(IOR0) == 0) CYCLE
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) '============= EXCHANGE_POISSON_VALS: GHOST_FIRSTW, NCG:', OL%GHOST_FIRSTW(IOR0), OL%GHOST_LASTW(IOR0)
#endif
               DO ICG = OL%GHOST_FIRSTW(IOR0), OL%GHOST_LASTW(IOR0)
                  ICW = OG%ICG_TO_ICW(ICG, 1)
                  DO ICOL = A%ROW(ICW), A%ROW(ICW+1)-1
                     OS%SEND_BUFFER_REAL(LL) = A%VAL(ICOL)
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,'(A, 4I5,E12.4)') 'EXCHANGE_POISSON_VALS: SEND: IOR0, ICG, ICW, ICOL, VAL:', IOR0, ICG, ICW, ICOL, A%VAL(ICOL)
#endif
                     LL = LL + 1
                  ENDDO
               ENDDO
            ENDDO

         CASE (NSCARC_GRID_UNSTRUCTURED)

            DO IOR0 = -3, 3
               IF (OL%GHOST_LASTW(IOR0) == 0) CYCLE
               DO ICG= OL%GHOST_FIRSTW(IOR0), OL%GHOST_LASTW(IOR0)
                  ICW = OG%ICG_TO_ICW(ICG, 1)
                  IF (ICW > 0) THEN
                     DO ICOL = A%ROW(ICW), A%ROW(ICW+1)-1
                        OS%SEND_BUFFER_REAL(LL) = A%VAL(A%COL(ICOL))
                        LL = LL + 1
                     ENDDO
                  ELSE
                     OS%SEND_BUFFER_REAL(LL) = NSCARC_HUGE_REAL
                     LL = LL + 1
                  ENDIF
               ENDDO
            ENDDO

         END SELECT
END SELECT

END SUBROUTINE SCARC_PACK_POISSON_VALS

! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_UNPACK_POISSON_VALS(NM, NOM, NL, NMATRIX)
USE SCARC_POINTERS, ONLY: OL, OA, OAB, RECV_BUFFER_REAL
INTEGER, INTENT(IN) :: NM, NOM, NL, NMATRIX
INTEGER :: IOR0, ICG, LL, ID

RECV_BUFFER_REAL => SCARC_POINT_TO_BUFFER_REAL (NM, NOM, 1)
OA => SCARC_POINT_TO_OTHER_CMATRIX(OG, NMATRIX)

SELECT CASE (SCARC_MATRIX_LEVEL(NL))

   ! Bandwise matrix on level NL
   CASE (NSCARC_MATRIX_BANDWISE)              

      LL = 1
      DO IOR0 = -3, 3
         IF (OL%GHOST_LASTW(IOR0) == 0) CYCLE
         DO ICG = OL%GHOST_FIRSTE(IOR0), OL%GHOST_LASTE(IOR0)
            DO ID = 1, OAB%N_STENCIL
               OAB%VAL(ICG, ID) = RECV_BUFFER_REAL(LL)
               LL = LL + 1
            ENDDO
         ENDDO
      ENDDO

   ! Compact matrix on level NL
   CASE (NSCARC_MATRIX_COMPACT)

      DO IOR0 = -3, 3
         IF (OL%GHOST_LASTW(IOR0) == 0) CYCLE
         DO ICG = OL%GHOST_FIRSTE(IOR0), OL%GHOST_LASTE(IOR0)
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'UNPACK_POISSON_VALS: ICG, OA%ROW(ICG)', ICG, OA%ROW(ICG)
#endif
            DO LL = OA%ROW(ICG), OA%ROW(ICG+1)-1
               OA%VAL(LL) = RECV_BUFFER_REAL(LL)
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,'(A,3I4,E12.4)') 'UNPACK_POISSON_VALS, RECV: IOR0, ICG, LL, OA%VAL', &
                                     IOR0, ICG, LL, OA%VAL(LL)
#endif
            ENDDO
         ENDDO
      ENDDO

END SELECT

END SUBROUTINE SCARC_UNPACK_POISSON_VALS


! ------------------------------------------------------------------------------------------------
! Pack and unpack information about matrix sizes into send vector
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_PACK_POISSON_SIZES(NL)
USE SCARC_POINTERS, ONLY: OS, OG
INTEGER, INTENT(IN) :: NL

SELECT CASE (SCARC_MATRIX_LEVEL(NL))
   CASE (NSCARC_MATRIX_BANDWISE)
      OS%SEND_BUFFER_INT0(1) = OG%POISSONB%N_VAL
      OS%SEND_BUFFER_INT0(2) = OG%POISSONB%N_DIAG
      OS%SEND_BUFFER_INT0(3) = OG%POISSONB%N_STENCIL
   CASE (NSCARC_MATRIX_COMPACT)
      OS%SEND_BUFFER_INT0(1) = OG%POISSONC%N_VAL
      OS%SEND_BUFFER_INT0(2) = OG%POISSONC%N_ROW
      OS%SEND_BUFFER_INT0(3) = OG%POISSONC%N_STENCIL
END SELECT

END SUBROUTINE SCARC_PACK_POISSON_SIZES
   
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_UNPACK_POISSON_SIZES(NM, NOM, NL)
USE SCARC_POINTERS, ONLY: OG, RECV_BUFFER_INT
INTEGER, INTENT(IN) :: NM, NOM, NL

RECV_BUFFER_INT => SCARC_POINT_TO_BUFFER_INT (NM, NOM, 0)
SELECT CASE (SCARC_MATRIX_LEVEL(NL))
   CASE (NSCARC_MATRIX_BANDWISE)
      OG%POISSONB%N_VAL     = RECV_BUFFER_INT(1)
      OG%POISSONB%N_DIAG    = RECV_BUFFER_INT(2)
      OG%POISSONB%N_STENCIL = RECV_BUFFER_INT(3)
   CASE (NSCARC_MATRIX_COMPACT)
      OG%POISSONC%N_VAL     = RECV_BUFFER_INT(1)
      OG%POISSONC%N_ROW     = RECV_BUFFER_INT(2)
      OG%POISSONC%N_STENCIL = RECV_BUFFER_INT(3)
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'OG%POISSONC%N_VAL=',OG%POISSONC%N_VAL
WRITE(MSG%LU_DEBUG,*) 'OG%POISSONC%N_ROW=',OG%POISSONC%N_ROW
WRITE(MSG%LU_DEBUG,*) 'OG%POISSONC%N_STENCIL=',OG%POISSONC%N_STENCIL
#endif
END SELECT

END SUBROUTINE SCARC_UNPACK_POISSON_SIZES


! ------------------------------------------------------------------------------------------------
! Pack and unpack overlapping information about matrix diagonals (compact storage technique only)
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_PACK_POISSON_DIAGS(NTYPE)
USE SCARC_POINTERS, ONLY:  A, OS, OL, OG
INTEGER, INTENT(IN) :: NTYPE
INTEGER :: IOR0, ICG, ICW, ICE, ICOL

A => SCARC_POINT_TO_CMATRIX(G, NTYPE)

IF (IS_STRUCTURED) THEN
   DO IOR0 = -3, 3
      IF (OL%GHOST_LASTW(IOR0) == 0) CYCLE
      DO ICG = OL%GHOST_FIRSTW(IOR0), OL%GHOST_LASTW(IOR0)
         ICW = OG%ICG_TO_ICW(ICG, 1)
         ICE = OG%ICG_TO_ICE(ICG, 1)
         ICOL = A%ROW(ICW)
         OS%SEND_BUFFER_REAL(ICG) = A%VAL(ICOL)
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,'(A, 5I4,E12.4)') 'PACK_POISSON_DIAGS: IOR0, ICG, ICW, ICE, ICOL, VAL:', &
                                      IOR0, ICG, ICW, ICE, ICOL, A%VAL(ICOL)
#endif
      ENDDO
   ENDDO
ELSE
   DO IOR0 = -3, 3
      IF (OL%GHOST_LASTW(IOR0) == 0) CYCLE
      DO ICG = OL%GHOST_FIRSTW(IOR0), OL%GHOST_LASTW(IOR0)
         ICW = OG%ICG_TO_ICW(ICG, 1)
         IF (ICW > 0) THEN
            ICE = OG%ICG_TO_ICE(ICG, 1)
            ICOL = A%ROW(ICW)
            OS%SEND_BUFFER_REAL(ICG) = A%VAL(ICOL)
         ELSE
            OS%SEND_BUFFER_REAL(ICG) = NSCARC_HUGE_REAL
         ENDIF
      ENDDO
   ENDDO
ENDIF
END SUBROUTINE SCARC_PACK_POISSON_DIAGS
      
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_UNPACK_POISSON_DIAGS(NM, NOM)
USE SCARC_POINTERS, ONLY: G, OL, OG, RECV_BUFFER_REAL
INTEGER, INTENT(IN) :: NM, NOM
INTEGER :: IOR0, ICG, ICE, LL

RECV_BUFFER_REAL => SCARC_POINT_TO_BUFFER_REAL (NM, NOM, 1)
LL = 1
DO IOR0 = -3, 3
   IF (OL%GHOST_LASTW(IOR0) == 0) CYCLE
   DO ICG = OL%GHOST_FIRSTE(IOR0), OL%GHOST_LASTE(IOR0)
      ICE = OG%ICG_TO_ICE(ICG, 1)
      G%DIAG(ICE) = RECV_BUFFER_REAL(LL)
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,'(A, 3I4,E12.4)') 'UNPACK_POISSON_DIAGS: NOM, IOR0, ICG, ICE, DIAG:', &
                                      IOR0, ICG, ICE, G%DIAG(ICE)
#endif
      LL = LL + 1
   ENDDO
ENDDO

END SUBROUTINE SCARC_UNPACK_POISSON_DIAGS


! ------------------------------------------------------------------------------------------------
! Pack and unpack zones numbers along interfaces
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_PACK_CELL_NEIGHBORS
USE SCARC_POINTERS, ONLY: L, G, OL, OG, A, F
INTEGER :: IOR0, ICG, ICW, ICWG, LL, IWG, IXW, IYW, IZW

A => SCARC_POINT_TO_CMATRIX(G, NSCARC_MATRIX_POISSON)
OS%SEND_BUFFER_INT = NSCARC_INIT_UNDEF
LL = 1
DO IOR0 = -3, 3
   IF (OL%GHOST_LASTW(IOR0) == 0) CYCLE
   DO ICG = OL%GHOST_FIRSTW(IOR0), OL%GHOST_LASTW(IOR0)
      ICW  = OG%ICG_TO_ICW(ICG, 1)                                 ! local cell number adjacent to interface
      ICWG = G%LOCAL_TO_GLOBAL(ICW)                            ! global cell number adjacent to interface
      OS%SEND_BUFFER_INT(LL)   = ICW
      OS%SEND_BUFFER_INT(LL+1) = ICWG
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,'(A,4I6)') 'PACK_CELL_NEIGHBORS:A: IOR0, ICG, ICW, ICWG:', IOR0, ICG, ICW, ICWG
#endif
      LL = LL + 2
   ENDDO
   F => L%FACE(IOR0)
   DO ICG = OL%GHOST_FIRSTW(IOR0), OL%GHOST_LASTW(IOR0)
      IWG  = OG%ICG_TO_IWG(ICG)
      IXW  = G%WALL(IWG)%IXW + F%NOFFX
      IYW  = G%WALL(IWG)%IYW + F%NOFFY
      IZW  = G%WALL(IWG)%IZW + F%NOFFZ
      ICW  = G%CELL_NUMBER(IXW, IYW, IZW)
      OG%ICG_TO_ICW(ICG, 2) = ICW                                 ! TODO: do it in separate routine?
      ICWG = G%LOCAL_TO_GLOBAL(ICW)
      OS%SEND_BUFFER_INT(LL)   = ICW
      OS%SEND_BUFFER_INT(LL+1) = ICWG
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,'(A,7I6)') 'PACK_CELL_NEIGHBORS:B: IOR0, IWG, IXW, IYW, IZW, ICW, ICWG: ', &
                                                      IOR0, IWG, IXW, IYW, IZW, ICW, ICWG
#endif
      LL = LL + 2
   ENDDO
ENDDO

END SUBROUTINE SCARC_PACK_CELL_NEIGHBORS


! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_UNPACK_CELL_NEIGHBORS(NM, NOM)
USE SCARC_POINTERS, ONLY: OL, OG, RECV_BUFFER_INT
INTEGER, INTENT(IN) :: NM, NOM
INTEGER :: ICG, IOR0, LL, IOFF

RECV_BUFFER_INT => SCARC_POINT_TO_BUFFER_INT (NM, NOM, 1)

#ifdef WITH_SCARC_DEBUG
   WRITE(MSG%LU_DEBUG,*) 'CELL_NEIGHBORS: RECV_BUFFER_INT: OG%NCG:', OG%NCG
   WRITE(MSG%LU_DEBUG,'(9I6)') RECV_BUFFER_INT
#endif
CALL SCARC_ALLOCATE_INT1(OG%ICG_TO_OCELL, 1, 2*OG%NCG, NSCARC_INIT_ZERO, 'G%ICG_TO_OCELL')
CALL SCARC_ALLOCATE_INT1(OG%ICG_TO_GCELL, 1, 2*OG%NCG, NSCARC_INIT_ZERO, 'G%ICG_TO_GCELL')

LL = 1
DO IOR0 = -3, 3
   IF (OL%GHOST_LASTW(IOR0) == 0) CYCLE
   DO ICG = OL%GHOST_FIRSTE(IOR0), OL%GHOST_LASTE(IOR0)
      OG%ICG_TO_OCELL(ICG) = RECV_BUFFER_INT(LL)
      OG%ICG_TO_GCELL(ICG) = RECV_BUFFER_INT(LL+1)
      LL = LL + 2
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'UNPACK_CELL_NEIGHBORS: IOR0, ICG_TO_OCELL(',ICG,'):',IOR0,OG%ICG_TO_OCELL(ICG)
WRITE(MSG%LU_DEBUG,*) 'UNPACK_CELL_NEIGHBORS: IOR0, ICG_TO_GCELL(',ICG,'):',IOR0,OG%ICG_TO_GCELL(ICG)
WRITE(MSG%LU_DEBUG,*) 'UNPACK_CELL_NEIGHBORS: IOR0, ICG_TO_ICE(',ICG,') OLD:',IOR0,OG%ICG_TO_ICE(ICG, 1)
#endif
   ENDDO
   IOFF = OL%GHOST_LASTW(IOR0)
#ifdef WITH_SCARC_DEBUG
   WRITE(MSG%LU_DEBUG,*) 'IOFF=',IOFF
#endif
   DO ICG = OL%GHOST_FIRSTE(IOR0), OL%GHOST_LASTE(IOR0)
      G%ICE2 = G%ICE2 + 1
      OG%ICG_TO_OCELL(ICG + IOFF) = RECV_BUFFER_INT(LL)
      OG%ICG_TO_GCELL(ICG + IOFF) = RECV_BUFFER_INT(LL+1)
      OG%ICG_TO_ICE(ICG, 2) = G%ICE2
      LL = LL + 2
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'UNPACK_CELL_NEIGHBORS: IOR0, ICG_TO_OCELL(',ICG+IOFF,'):',IOR0,OG%ICG_TO_OCELL(ICG+IOFF)
WRITE(MSG%LU_DEBUG,*) 'UNPACK_CELL_NEIGHBORS: IOR0, ICG_TO_GCELL(',ICG+IOFF,'):',IOR0,OG%ICG_TO_GCELL(ICG+IOFF)
WRITE(MSG%LU_DEBUG,*) 'UNPACK_CELL_NEIGHBORS: IOR0, ICG_TO_ICE2(',ICG+IOFF,') :',IOR0,OG%ICG_TO_ICE(ICG,2)
#endif
   ENDDO
#ifdef WITH_SCARC_DEBUG
   WRITE(MSG%LU_DEBUG,*) 'CELL_NEIGHBORS: IOR0=', IOR0
   WRITE(MSG%LU_DEBUG,*) '===============> LCELL:'
   WRITE(MSG%LU_DEBUG,'(6I6)') (OG%ICG_TO_OCELL(ICG), ICG = 1, 2*OG%NCG)
   WRITE(MSG%LU_DEBUG,*) '===============> GCELL:'
   WRITE(MSG%LU_DEBUG,'(6I6)') (OG%ICG_TO_GCELL(ICG), ICG = 1, 2*OG%NCG)
   WRITE(MSG%LU_DEBUG,*) '===============> ICE:'
   WRITE(MSG%LU_DEBUG,'(6I6)') (OG%ICG_TO_ICE(ICG, 1), ICG = 1, OG%NCG)
   WRITE(MSG%LU_DEBUG,'(6I6)') (OG%ICG_TO_ICE(ICG, 2), ICG = 1, OG%NCG)
#endif
ENDDO

END SUBROUTINE SCARC_UNPACK_CELL_NEIGHBORS


! ------------------------------------------------------------------------------------------------
! Pack and unpack zones numbers along interfaces
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_PACK_ZONE_NEIGHBORS
USE SCARC_POINTERS, ONLY: L, G, OL, OG, A, F
INTEGER :: IOR0, ICG, ICW, LL, IWG, IXW, IYW, IZW, IZWG

A => SCARC_POINT_TO_CMATRIX(G, NSCARC_MATRIX_POISSON)
OS%SEND_BUFFER_INT = NSCARC_INIT_UNDEF
LL = 1
DO IOR0 = -3, 3
   IF (OL%GHOST_LASTW(IOR0) == 0) CYCLE
   DO ICG = OL%GHOST_FIRSTW(IOR0), OL%GHOST_LASTW(IOR0)
      ICW  = OG%ICG_TO_ICW(ICG, 1)
      IZW  = G%ZONES_LOCAL(ICW)                                ! local zone number adjacent to interface
      IZWG = G%ZONES_GLOBAL(ICW)                               ! global zone number adjacent to interface
      OS%SEND_BUFFER_INT(LL  ) = IZW
      OS%SEND_BUFFER_INT(LL+1) = IZWG
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,'(A,5I6)') 'PACK_ZONE_NEIGHBORS:A: IOR0, ICG, ICW, IZL, IZG:', IOR0, ICG, ICW, IZW, IZWG
#endif
      LL = LL + 2
   ENDDO
   F => L%FACE(IOR0)
   DO ICG = OL%GHOST_FIRSTW(IOR0), OL%GHOST_LASTW(IOR0)
      IWG  = OG%ICG_TO_IWG(ICG)
      IXW  = G%WALL(IWG)%IXW + F%NOFFX
      IYW  = G%WALL(IWG)%IYW + F%NOFFY
      IZW  = G%WALL(IWG)%IZW + F%NOFFZ
      ICW  = G%CELL_NUMBER(IXW, IYW, IZW)
      IZW  = G%ZONES_LOCAL(ICW)                       
      IZWG = G%ZONES_GLOBAL(ICW)                      
      OS%SEND_BUFFER_INT(LL  ) = IZW
      OS%SEND_BUFFER_INT(LL+1) = IZWG
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,'(A,8I6)') 'PACK_CELL_NEIGHBORS:B: IOR0, IWG, IXW, IYW, IZW, ICW, IZW, IZWG: ', &
                                                      IOR0, IWG, IXW, IYW, IZW, ICW, IZW, IZWG
#endif
      LL = LL + 2
   ENDDO
ENDDO
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,'(8I6)') OS%SEND_BUFFER_INT(1:2*OG%NCG)
#endif

END SUBROUTINE SCARC_PACK_ZONE_NEIGHBORS

! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_UNPACK_ZONE_NEIGHBORS(NM, NOM)
USE SCARC_POINTERS, ONLY: OL, OG, RECV_BUFFER_INT
INTEGER, INTENT(IN) :: NM, NOM
INTEGER :: ICG, IOR0, LL, KK

RECV_BUFFER_INT => SCARC_POINT_TO_BUFFER_INT (NM, NOM, 1)

CALL SCARC_ALLOCATE_INT1(OG%ICG_TO_OZONE, 1, 2*OG%NCG, NSCARC_INIT_ZERO, 'G%ICG_TO_OZONE')
CALL SCARC_ALLOCATE_INT1(OG%ICG_TO_GZONE, 1, 2*OG%NCG, NSCARC_INIT_ZERO, 'G%ICG_TO_GZONE')

LL = 1
KK = 1
DO IOR0 = -3, 3
   IF (OL%GHOST_LASTW(IOR0) == 0) CYCLE
   DO ICG = OL%GHOST_FIRSTE(IOR0), OL%GHOST_LASTE(IOR0)
      OG%ICG_TO_OZONE(KK) = RECV_BUFFER_INT(LL)
      OG%ICG_TO_GZONE(KK) = RECV_BUFFER_INT(LL+1)
      KK = KK + 1
      LL = LL + 2
#ifdef WITH_SCARC_DEBUG2
WRITE(MSG%LU_DEBUG,'(A,6I6)') 'UNPACK_ZONE_NEIGHBORS: IOR0, ICG_TO_OZONE(.):',IOR0,OG%ICG_TO_OZONE(KK)
WRITE(MSG%LU_DEBUG,'(A,6I6)') 'UNPACK_ZONE_NEIGHBORS: IOR0, ICG_TO_GZONE(.):',IOR0,OG%ICG_TO_GZONE(KK)
#endif
   ENDDO
   DO ICG = OL%GHOST_FIRSTE(IOR0), OL%GHOST_LASTE(IOR0)
      OG%ICG_TO_OZONE(KK) = RECV_BUFFER_INT(LL)
      OG%ICG_TO_GZONE(KK) = RECV_BUFFER_INT(LL+1)
      KK = KK + 1
      LL = LL + 2
#ifdef WITH_SCARC_DEBUG2
WRITE(MSG%LU_DEBUG,'(A,6I6)') 'UNPACK_ZONE_NEIGHBORS: IOR0, ICG_TO_OZONE(.):',IOR0,OG%ICG_TO_OZONE(KK)
WRITE(MSG%LU_DEBUG,'(A,6I6)') 'UNPACK_ZONE_NEIGHBORS: IOR0, ICG_TO_GZONE(.):',IOR0,OG%ICG_TO_GZONE(KK)
#endif
   ENDDO
#ifdef WITH_SCARC_DEBUG
   WRITE(MSG%LU_DEBUG,*) 'ZONE_NEIGHBORS: IOR0=', IOR0
   WRITE(MSG%LU_DEBUG,*) '===============> OZONE:'
   WRITE(MSG%LU_DEBUG,'(10I6)') (OG%ICG_TO_OZONE(ICG), ICG = 1, 2*OG%NCG)
   WRITE(MSG%LU_DEBUG,*) '===============> GZONE:'
   WRITE(MSG%LU_DEBUG,'(10I6)') (OG%ICG_TO_GZONE(ICG), ICG = 1, 2*OG%NCG)
#endif
ENDDO
END SUBROUTINE SCARC_UNPACK_ZONE_NEIGHBORS


! ------------------------------------------------------------------------------------------------
! Pack ZONES information into send vector
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_PACK_ZONE_TYPES
USE SCARC_POINTERS, ONLY: G, OS, OL, OG
INTEGER :: IOR0, ICG, ICW, ICE, LL
LL = 1
OS%SEND_BUFFER_INT = NSCARC_INIT_UNDEF
DO IOR0 = -3, 3
   IF (OL%GHOST_LASTW(IOR0) == 0) CYCLE

   DO ICG = OL%GHOST_FIRSTW(IOR0), OL%GHOST_LASTW(IOR0)
      ICW = OG%ICG_TO_ICW(ICG, 1)
      ICE = OG%ICG_TO_ICE(ICG, 1)
      IF (G%ZONES_LOCAL(ICW) /= 0) THEN
!         IF (G%CPOINTS(G%ZONES_LOCAL(ICW)) == ICW) THEN
!            OS%SEND_BUFFER_INT(LL) = -G%ZONES_LOCAL(ICW)
!         ELSE
            OS%SEND_BUFFER_INT(LL) = G%ZONES_LOCAL(ICW)
!         ENDIF
         OS%SEND_BUFFER_INT(LL+1) =  G%ZONES_LOCAL(ICE)
      ELSE
         OS%SEND_BUFFER_INT(LL)   = 0
         OS%SEND_BUFFER_INT(LL+1) = 0
      ENDIF
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,'(A, 6I8)') 'PACK_ZONE_TYPES: IOR0, ICG, ICW, SEND(LL:LL+1):', &
                                IOR0, ICG, ICW, ICE, OS%SEND_BUFFER_INT(LL), OS%SEND_BUFFER_INT(LL+1)
#endif
      LL = LL + 2
   ENDDO
ENDDO
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'PACK_ZONE_TYPES: SENDING'
WRITE(MSG%LU_DEBUG,'(6I8)') OS%SEND_BUFFER_INT(1:2*OG%NCG)
#endif
END SUBROUTINE SCARC_PACK_ZONE_TYPES


! ------------------------------------------------------------------------------------------------
! Unpack aggregation zone information from receive vector
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_UNPACK_ZONE_TYPES(NM, NOM)
USE SCARC_POINTERS, ONLY: G, OL, OG, RECV_BUFFER_INT
INTEGER, INTENT(IN) :: NM, NOM
INTEGER :: IOR0, LL, ICG, ICW, ICE, IZW, IZE !, IFOUND 

RECV_BUFFER_INT => SCARC_POINT_TO_BUFFER_INT (NM, NOM, 1)
LL = 1
DO IOR0 = -3, 3
   IF (OL%GHOST_LASTW(IOR0) == 0) CYCLE
   DO ICG = OL%GHOST_FIRSTE(IOR0), OL%GHOST_LASTE(IOR0)
      ICE = OG%ICG_TO_ICE(ICG, 1)
      ICW = OG%ICG_TO_ICW(ICG, 1)
      IZE = RECV_BUFFER_INT(LL)
      IZW = RECV_BUFFER_INT(LL+1)
      IF (G%ZONES_LOCAL(ICW) == 0) G%ZONES_LOCAL(ICW) = IZW
      IF (G%ZONES_LOCAL(ICE) == 0) G%ZONES_LOCAL(ICE) = IZE
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,'(A, 7I8)') 'UNPACK_ZONE_TYPES: IOR0, ICG, ICW, ICE, IZW, IZE, N_ZONES:', &
              IOR0, ICG, ICW, ICE, IZW, IZE, G%N_ZONES
#endif
      LL = LL + 2
   ENDDO
ENDDO

END SUBROUTINE SCARC_UNPACK_ZONE_TYPES


! ------------------------------------------------------------------------------------------------
! Determine if cell should be considered during zone-numbers packaging
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
! Check if difference of two values is less than a given tolerance
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
! Filter out mean value
! --------------------------------------------------------------------------------------------------ss
SUBROUTINE SCARC_FILTER_MEANVALUE(NV, NL)
USE SCARC_POINTERS, ONLY: L, G, VC
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
! Restore last cell of last mesh
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_RESTORE_LAST_CELL (XX, NL)
USE SCARC_POINTERS, ONLY: S, VC
INTEGER, INTENT(IN) :: XX, NL

IF (UPPER_MESH_INDEX /= NMESHES .OR. TYPE_RELAX == NSCARC_RELAX_FFT) RETURN
S => SCARC(UPPER_MESH_INDEX)

VC => SCARC_POINT_TO_VECTOR (UPPER_MESH_INDEX, NL, XX)
VC(S%NC) = S%RHS_END

END SUBROUTINE SCARC_RESTORE_LAST_CELL



! ================================================================================================
! ================================================================================================
!  Bundle of different pointer routines returning pointers to specified structures
! ================================================================================================
! ================================================================================================

! -----------------------------------------------------------------------------
! Point to mesh
! -----------------------------------------------------------------------------
SUBROUTINE SCARC_POINT_TO_MESH(NM)
USE SCARC_POINTERS, ONLY: M, S
INTEGER, INTENT(IN) :: NM
M => MESHES(NM)
S => SCARC(NM)
END SUBROUTINE SCARC_POINT_TO_MESH

! -----------------------------------------------------------------------------
! This routine sets the requested combination of mesh and level 
! -----------------------------------------------------------------------------
SUBROUTINE SCARC_POINT_TO_LEVEL(NM, NL)
USE SCARC_POINTERS, ONLY: M, S, L
INTEGER, INTENT(IN) :: NM, NL
M => MESHES(NM)
S => SCARC(NM)
L => S%LEVEL(NL)
END SUBROUTINE SCARC_POINT_TO_LEVEL


! ----------------------------------------------------------------------------------------------------
! Unset all pointers - 
! mainly used to test the correctness of the pointer settings in the different routines
! ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_POINT_TO_NONE 
USE SCARC_POINTERS, ONLY: M, S, L, G, F, A, P, R, C, Z, W, &
                                LC, GC, FC, AC, PC, RC, CC, ZC, &
                                LF, GF, FF, AF, PF, RF, CF, ZF
M => NULL()                 
S => NULL()                 
W => NULL()                 
L => NULL();  LF => NULL();  LC => NULL()                 
G => NULL();  GF => NULL();  GC => NULL()                 
F => NULL();  FF => NULL();  FC => NULL()                 
A => NULL();  AF => NULL();  AC => NULL()                 
P => NULL();  PF => NULL();  PC => NULL()                 
R => NULL();  RF => NULL();  RC => NULL()                 
C => NULL();  CF => NULL();  CC => NULL()                 
Z => NULL();  ZF => NULL();  ZC => NULL()                 

END SUBROUTINE SCARC_POINT_TO_NONE

! ----------------------------------------------------------------------------------------------------
! This routine sets the requested combination of mesh, level and discretization
! ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_POINT_TO_GRID (NM, NL)                              
USE SCARC_POINTERS, ONLY: M, S, L, G, W
INTEGER, INTENT(IN) ::  NM, NL

CALL SCARC_POINT_TO_NONE

M => MESHES(NM)
S => SCARC(NM)
L => S%LEVEL(NL)
SELECT CASE(TYPE_GRID)
   CASE (NSCARC_GRID_STRUCTURED)
      G => L%STRUCTURED
      G%NW = L%N_WALL_CELLS_EXT 
   CASE (NSCARC_GRID_UNSTRUCTURED)
      G => L%UNSTRUCTURED
      G%NW = L%N_WALL_CELLS_EXT + L%N_WALL_CELLS_INT
END SELECT
W => G%WALL

END SUBROUTINE SCARC_POINT_TO_GRID

! ----------------------------------------------------------------------------------------------------
! This routine sets the requested combination of mesh, level and discretization
! ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_POINT_TO_MULTIGRID(NM, NL1, NL2)
USE SCARC_POINTERS, ONLY: M, S, LF, LC, GF, GC
INTEGER, INTENT(IN) ::  NM, NL1, NL2

CALL SCARC_POINT_TO_NONE

M => MESHES(NM)
S => SCARC(NM)
LF => SCARC(NM)%LEVEL(NL1)
LC => SCARC(NM)%LEVEL(NL2)
SELECT CASE(TYPE_GRID)
   CASE (NSCARC_GRID_STRUCTURED)
      GF => LF%STRUCTURED
      GC => LC%STRUCTURED
      GF%NW = LF%N_WALL_CELLS_EXT 
      GC%NW = LC%N_WALL_CELLS_EXT 
   CASE (NSCARC_GRID_UNSTRUCTURED)
      GF => LF%UNSTRUCTURED
      GC => LC%UNSTRUCTURED
      GF%NW = LF%N_WALL_CELLS_EXT + LF%N_WALL_CELLS_INT
      GC%NW = LC%N_WALL_CELLS_EXT + LC%N_WALL_CELLS_INT
END SELECT

END SUBROUTINE SCARC_POINT_TO_MULTIGRID

! -----------------------------------------------------------------------------
! This routine sets the requested combination of mesh and level 
! -----------------------------------------------------------------------------
SUBROUTINE SCARC_POINT_TO_OTHER_GRID(NM, NOM, NL)
USE SCARC_POINTERS, ONLY : OS, OL, OLF, OG, OGF
INTEGER, INTENT(IN) :: NM, NOM, NL

OS  => SCARC(NM)%OSCARC(NOM)
OL  => OS%LEVEL(NL)
OLF => OS%LEVEL(NL)

SELECT CASE(TYPE_GRID)
   CASE (NSCARC_GRID_STRUCTURED)
      OG  => OL%STRUCTURED
      OGF => OLF%STRUCTURED
   CASE (NSCARC_GRID_UNSTRUCTURED)
      OG  => OL%UNSTRUCTURED
      OGF => OLF%UNSTRUCTURED
END SELECT

END SUBROUTINE SCARC_POINT_TO_OTHER_GRID


! -----------------------------------------------------------------------------
! This routine sets the requested combination of mesh and level 
! -----------------------------------------------------------------------------
SUBROUTINE SCARC_POINT_TO_OTHER_MULTIGRID(NM, NOM, NL1, NL2)
USE SCARC_POINTERS, ONLY : OS, OLC, OLF, OGC, OGF
INTEGER, INTENT(IN) :: NM, NOM, NL1, NL2

OS  => SCARC(NM)%OSCARC(NOM)
OLF => OS%LEVEL(NL1)
OLC => OS%LEVEL(NL2)

SELECT CASE(TYPE_GRID)
   CASE (NSCARC_GRID_STRUCTURED)
      OGF => OLF%STRUCTURED
      OGC => OLC%STRUCTURED
   CASE (NSCARC_GRID_UNSTRUCTURED)
      OGF => OLF%UNSTRUCTURED
      OGC => OLC%UNSTRUCTURED
END SELECT

END SUBROUTINE SCARC_POINT_TO_OTHER_MULTIGRID


! ----------------------------------------------------------------------------------------------------
! This routine sets a pointer to the requested matrix of compact storage type
! ----------------------------------------------------------------------------------------------------
FUNCTION SCARC_POINT_TO_CMATRIX(G, NTYPE)
TYPE(SCARC_MATRIX_COMPACT_TYPE), POINTER :: SCARC_POINT_TO_CMATRIX
TYPE(SCARC_GRID_TYPE), POINTER, INTENT(IN) :: G
INTEGER, INTENT(IN) :: NTYPE

SELECT CASE(NTYPE)
   CASE (NSCARC_MATRIX_POISSON_PROL)
      SCARC_POINT_TO_CMATRIX => G%AP
   CASE (NSCARC_MATRIX_CONNECTION)
      SCARC_POINT_TO_CMATRIX => G%CONNECTION
   CASE (NSCARC_MATRIX_POISSON)
      SCARC_POINT_TO_CMATRIX => G%POISSONC
#ifdef WITH_MKL
   CASE (NSCARC_MATRIX_POISSON_SYM)
      SCARC_POINT_TO_CMATRIX => G%POISSONC_SYM
#endif
   CASE (NSCARC_MATRIX_PROLONGATION)
      SCARC_POINT_TO_CMATRIX => G%PROLONGATION
   CASE (NSCARC_MATRIX_RESTRICTION)
      SCARC_POINT_TO_CMATRIX => G%RESTRICTION
   CASE (NSCARC_MATRIX_ZONES)
      SCARC_POINT_TO_CMATRIX => G%ZONES
END SELECT

END FUNCTION SCARC_POINT_TO_CMATRIX


! ----------------------------------------------------------------------------------------------------
! This routine sets a pointer to the requested matrix of compact storage type
! ----------------------------------------------------------------------------------------------------
FUNCTION SCARC_POINT_TO_BMATRIX(G, NTYPE)
TYPE(SCARC_MATRIX_BANDWISE_TYPE), POINTER :: SCARC_POINT_TO_BMATRIX
TYPE(SCARC_GRID_TYPE), POINTER, INTENT(IN) :: G
INTEGER, INTENT(IN) :: NTYPE

SELECT CASE(NTYPE)
   CASE (NSCARC_MATRIX_POISSON)
      SCARC_POINT_TO_BMATRIX => G%POISSONB
   CASE DEFAULT
      WRITE(*,*) 'No other bandwise matrix available yet except of POISSONB'
END SELECT

END FUNCTION SCARC_POINT_TO_BMATRIX


! ----------------------------------------------------------------------------------------------------
! This routine sets a pointer to the requested matrix of compact storage type
! ----------------------------------------------------------------------------------------------------
FUNCTION SCARC_POINT_TO_OTHER_CMATRIX(OG, NTYPE)
TYPE(SCARC_MATRIX_COMPACT_TYPE), POINTER :: SCARC_POINT_TO_OTHER_CMATRIX
TYPE(SCARC_GRID_TYPE), POINTER, INTENT(IN) :: OG
INTEGER, INTENT(IN) :: NTYPE

SELECT CASE(NTYPE)
   CASE (NSCARC_MATRIX_POISSON_PROL)
      SCARC_POINT_TO_OTHER_CMATRIX => OG%AP
   CASE (NSCARC_MATRIX_CONNECTION)
      SCARC_POINT_TO_OTHER_CMATRIX => OG%CONNECTION
   CASE (NSCARC_MATRIX_POISSON)
      SCARC_POINT_TO_OTHER_CMATRIX => OG%POISSONC
#ifdef WITH_MKL
   CASE (NSCARC_MATRIX_POISSON_SYM)
      SCARC_POINT_TO_OTHER_CMATRIX => OG%POISSONC_SYM
#endif
   CASE (NSCARC_MATRIX_PROLONGATION)
      SCARC_POINT_TO_OTHER_CMATRIX => OG%PROLONGATION
   CASE (NSCARC_MATRIX_RESTRICTION)
      SCARC_POINT_TO_OTHER_CMATRIX => OG%RESTRICTION
   CASE (NSCARC_MATRIX_ZONES)
      SCARC_POINT_TO_OTHER_CMATRIX => OG%ZONES
END SELECT

END FUNCTION SCARC_POINT_TO_OTHER_CMATRIX


! ----------------------------------------------------------------------------------------------------
! This routine sets a pointer to the requested matrix of compact storage type
! ----------------------------------------------------------------------------------------------------
FUNCTION SCARC_POINT_TO_OTHER_BMATRIX(OG, NTYPE)
TYPE(SCARC_MATRIX_BANDWISE_TYPE), POINTER :: SCARC_POINT_TO_OTHER_BMATRIX
TYPE(SCARC_GRID_TYPE), POINTER, INTENT(IN) :: OG
INTEGER, INTENT(IN) :: NTYPE

SELECT CASE(NTYPE)
   CASE (NSCARC_MATRIX_POISSON)
      SCARC_POINT_TO_OTHER_BMATRIX => OG%POISSONB
   CASE DEFAULT
      WRITE(*,*) 'No other bandwise matrix available yet except of POISSONB'
END SELECT

END FUNCTION SCARC_POINT_TO_OTHER_BMATRIX


! -----------------------------------------------------------------------------
! Point to requested integer receive buffer
! -----------------------------------------------------------------------------
FUNCTION SCARC_POINT_TO_BUFFER_INT(NM, NOM, NTYPE)
INTEGER, DIMENSION(:), POINTER :: SCARC_POINT_TO_BUFFER_INT
INTEGER, INTENT(IN) ::  NM, NOM, NTYPE

SCARC_POINT_TO_BUFFER_INT => NULL()
SELECT CASE (NTYPE)
   CASE (0)
      IF (RNODE/=SNODE) THEN
         SCARC_POINT_TO_BUFFER_INT => SCARC(NM)%OSCARC(NOM)%RECV_BUFFER_INT0
      ELSE
         SCARC_POINT_TO_BUFFER_INT => SCARC(NOM)%OSCARC(NM)%SEND_BUFFER_INT0
      ENDIF
   CASE (1)
      IF (RNODE/=SNODE) THEN
         SCARC_POINT_TO_BUFFER_INT => SCARC(NM)%OSCARC(NOM)%RECV_BUFFER_INT
      ELSE
         SCARC_POINT_TO_BUFFER_INT => SCARC(NOM)%OSCARC(NM)%SEND_BUFFER_INT
      ENDIF
END SELECT

END FUNCTION SCARC_POINT_TO_BUFFER_INT


! -----------------------------------------------------------------------------
! Point to requested real receive buffer
! -----------------------------------------------------------------------------
FUNCTION SCARC_POINT_TO_BUFFER_REAL(NM, NOM, NTYPE)
REAL(EB), DIMENSION(:), POINTER :: SCARC_POINT_TO_BUFFER_REAL
INTEGER, INTENT(IN) ::  NM, NOM, NTYPE

SCARC_POINT_TO_BUFFER_REAL => NULL()
SELECT CASE (NTYPE)
   CASE (0)
      IF (RNODE/=SNODE) THEN
         SCARC_POINT_TO_BUFFER_REAL => SCARC(NM)%OSCARC(NOM)%RECV_BUFFER_REAL0
      ELSE
         SCARC_POINT_TO_BUFFER_REAL => SCARC(NOM)%OSCARC(NM)%SEND_BUFFER_REAL0
      ENDIF
   CASE (1)
      IF (RNODE/=SNODE) THEN
         SCARC_POINT_TO_BUFFER_REAL => SCARC(NM)%OSCARC(NOM)%RECV_BUFFER_REAL
      ELSE
         SCARC_POINT_TO_BUFFER_REAL => SCARC(NOM)%OSCARC(NM)%SEND_BUFFER_REAL
      ENDIF
END SELECT

END FUNCTION SCARC_POINT_TO_BUFFER_REAL


! -----------------------------------------------------------------------------
! Point to exchange structure
! -----------------------------------------------------------------------------
LOGICAL FUNCTION ARE_NEIGHBORS(NM, NOM)
USE SCARC_POINTERS, ONLY: OM
INTEGER, INTENT(IN) :: NM, NOM

ARE_NEIGHBORS = .TRUE.

OM => MESHES(NM)%OMESH(NOM)
IF (OM%NIC_R == 0 .AND. OM%NIC_S == 0) ARE_NEIGHBORS = .FALSE.

END FUNCTION ARE_NEIGHBORS

! ------------------------------------------------------------------------------------------------
! Set vector pointer on given combination of mesh and level
! ------------------------------------------------------------------------------------------------
FUNCTION SCARC_POINT_TO_VECTOR (NM, NL, NV)
USE SCARC_POINTERS, ONLY: L
REAL(EB), POINTER, DIMENSION(:) :: SCARC_POINT_TO_VECTOR
INTEGER, INTENT(IN) :: NM, NL, NV

L => SCARC(NM)%LEVEL(NL)
SELECT CASE (NV)

   !
   ! Stage one vectors (for methods on first hierarchical level)
   !
   CASE (NSCARC_VECTOR_ONE_X)
      SCARC_POINT_TO_VECTOR => L%STAGE(NSCARC_STAGE_ONE)%X
   CASE (NSCARC_VECTOR_ONE_B)
      SCARC_POINT_TO_VECTOR => L%STAGE(NSCARC_STAGE_ONE)%B
   CASE (NSCARC_VECTOR_ONE_D)
      SCARC_POINT_TO_VECTOR => L%STAGE(NSCARC_STAGE_ONE)%D
   CASE (NSCARC_VECTOR_ONE_R)
      SCARC_POINT_TO_VECTOR => L%STAGE(NSCARC_STAGE_ONE)%R
   CASE (NSCARC_VECTOR_ONE_V)
      SCARC_POINT_TO_VECTOR => L%STAGE(NSCARC_STAGE_ONE)%V
   CASE (NSCARC_VECTOR_ONE_Y)
      SCARC_POINT_TO_VECTOR => L%STAGE(NSCARC_STAGE_ONE)%Y
   CASE (NSCARC_VECTOR_ONE_Z)
      SCARC_POINT_TO_VECTOR => L%STAGE(NSCARC_STAGE_ONE)%Z
#ifdef WITH_SCARC_DEBUG
   CASE (NSCARC_VECTOR_ONE_E)
      SCARC_POINT_TO_VECTOR => L%STAGE(NSCARC_STAGE_ONE)%E
#endif

   !
   ! Stage two vectors (for methods on second hierarchical level)
   !
   CASE (NSCARC_VECTOR_TWO_X)
      SCARC_POINT_TO_VECTOR => L%STAGE(NSCARC_STAGE_TWO)%X
   CASE (NSCARC_VECTOR_TWO_B)
      SCARC_POINT_TO_VECTOR => L%STAGE(NSCARC_STAGE_TWO)%B
   CASE (NSCARC_VECTOR_TWO_D)
      SCARC_POINT_TO_VECTOR => L%STAGE(NSCARC_STAGE_TWO)%D
   CASE (NSCARC_VECTOR_TWO_R)
      SCARC_POINT_TO_VECTOR => L%STAGE(NSCARC_STAGE_TWO)%R
   CASE (NSCARC_VECTOR_TWO_V)
      SCARC_POINT_TO_VECTOR => L%STAGE(NSCARC_STAGE_TWO)%V
   CASE (NSCARC_VECTOR_TWO_Y)
      SCARC_POINT_TO_VECTOR => L%STAGE(NSCARC_STAGE_TWO)%Y
   CASE (NSCARC_VECTOR_TWO_Z)
      SCARC_POINT_TO_VECTOR => L%STAGE(NSCARC_STAGE_TWO)%Z
#ifdef WITH_SCARC_DEBUG
   CASE (NSCARC_VECTOR_TWO_E)
      SCARC_POINT_TO_VECTOR => L%STAGE(NSCARC_STAGE_TWO)%E
#endif
END SELECT

END FUNCTION SCARC_POINT_TO_VECTOR


! ------------------------------------------------------------------------------------------------
! Single precision: Set vector pointer on given combination of mesh and level
! ------------------------------------------------------------------------------------------------
FUNCTION SCARC_POINT_TO_VECTOR_FB(NM, NL, NV)
USE SCARC_POINTERS, ONLY: L
REAL(FB), POINTER, DIMENSION(:) :: SCARC_POINT_TO_VECTOR_FB
INTEGER, INTENT(IN) :: NM, NL, NV

L => SCARC(NM)%LEVEL(NL)
SELECT CASE (NV)
   CASE (NSCARC_VECTOR_ONE_X)
      SCARC_POINT_TO_VECTOR_FB => L%STAGE(NSCARC_STAGE_ONE)%X_FB
   CASE (NSCARC_VECTOR_ONE_B)
      SCARC_POINT_TO_VECTOR_FB => L%STAGE(NSCARC_STAGE_ONE)%B_FB
   CASE (NSCARC_VECTOR_ONE_R)
      SCARC_POINT_TO_VECTOR_FB => L%STAGE(NSCARC_STAGE_ONE)%R_FB
   CASE (NSCARC_VECTOR_ONE_V)
      SCARC_POINT_TO_VECTOR_FB => L%STAGE(NSCARC_STAGE_ONE)%V_FB
END SELECT

END FUNCTION SCARC_POINT_TO_VECTOR_FB


! ------------------------------------------------------------------------------------------------
! Set pointer to chosen vector for bandwise storage technique
! ------------------------------------------------------------------------------------------------
FUNCTION POINT_TO_HVECTOR(NM, NV)
REAL(EB), POINTER, DIMENSION(:,:,:) :: POINT_TO_HVECTOR
INTEGER, INTENT(IN) :: NM, NV
SELECT CASE (NV)
   CASE (NSCARC_VECTOR_H)
      POINT_TO_HVECTOR => MESHES(NM)%H
   CASE (NSCARC_VECTOR_HS)
      POINT_TO_HVECTOR => MESHES(NM)%HS
END SELECT
END FUNCTION POINT_TO_HVECTOR




! ================================================================================================
! ================================================================================================
! Bundle of routines for different vector-operations, also based on use of OpenMP
! ================================================================================================
! ================================================================================================

! ------------------------------------------------------------------------------------------------
! Vector multiplied with a constant scalar is added to another vector 
!     DY(I) = DA * DX(I) + DY(I) 
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_DAXPY_CONSTANT(N, DA, DX, DY)
REAL(EB), INTENT(IN):: DA
REAL(EB), INTENT(IN), DIMENSION(:):: DX
REAL(EB), INTENT(INOUT), DIMENSION(:):: DY
INTEGER, INTENT(IN)::  N
INTEGER::  I

!$OMP PARALLEL DO PRIVATE(I) SCHEDULE(STATIC)
DO I = 1, N
  DY(I) = DY(I) + DA * DX(I)
ENDDO
!$OMP END PARALLEL DO 

END SUBROUTINE SCARC_DAXPY_CONSTANT


! ------------------------------------------------------------------------------------------------
! Vector multiplied with a constant scalar is added to vector multiplied with another scalar
!     DY(I) = DA1 * DX(I) + DA2 * DY(I) 
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_DAXPY_CONSTANT_DOUBLE(N, DA1, DX, DA2, DY)
REAL(EB), INTENT(IN):: DA1, DA2
REAL(EB), INTENT(IN), DIMENSION(:):: DX
REAL(EB), INTENT(INOUT), DIMENSION(:):: DY
INTEGER, INTENT(IN)::  N
INTEGER::  I

!$OMP PARALLEL DO PRIVATE(I) SCHEDULE(STATIC)
DO I = 1, N
  DY(I) = DA1 * DX(I) + DA2 * DY(I)
ENDDO
!$OMP END PARALLEL DO 

END SUBROUTINE SCARC_DAXPY_CONSTANT_DOUBLE


! ------------------------------------------------------------------------------------------------
! Vector multiplied with variable scalars (componentwise) is added to another vector 
!     DY(I) = DA(I)*DX(I) + DY(I)
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_DAXPY_VARIABLE(N, DA, DX, DY)
REAL(EB), INTENT(IN), DIMENSION(:):: DA, DX
REAL(EB), INTENT(INOUT), DIMENSION(:):: DY
INTEGER, INTENT(IN)::  N
INTEGER::  I

!$OMP PARALLEL DO PRIVATE(I) SCHEDULE(STATIC)
DO I = 1, N
  DY(I) = DY(I) + DA(I) * DX(I)
#ifdef WITH_SCARC_DEBUG2
WRITE(MSG%LU_DEBUG,'(A, I6, 3E14.6)') 'I, DX, DA, DY:', I, DX(I), DA(I), DY(I)
#endif
ENDDO
!$OMP END PARALLEL DO 

END SUBROUTINE SCARC_DAXPY_VARIABLE


! ------------------------------------------------------------------------------------------------
! Vector is multiplied with a constant scalar 
!     DY(I) = DA(I)*DX(I) 
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SCALING_CONSTANT(N, SCAL, DX, DY)
REAL(EB), INTENT(IN), DIMENSION(:):: DX
REAL(EB), INTENT(INOUT), DIMENSION(:):: DY
REAL(EB), INTENT(IN) :: SCAL
INTEGER, INTENT(IN)::  N
INTEGER::  I

!$OMP PARALLEL DO PRIVATE(I) SCHEDULE(STATIC)
DO I = 1, N
  DY(I) =  SCAL * DX(I)
ENDDO
!$OMP END PARALLEL DO 

END SUBROUTINE SCARC_SCALING_CONSTANT


! ------------------------------------------------------------------------------------------------
! Vector is multiplied with variable scalars (componentwise)
!     DY(I) = DA(I)*DX(I) 
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SCALING_VARIABLE(N, DA, DX, DY)
REAL(EB), INTENT(IN), DIMENSION(:):: DA, DX
REAL(EB), INTENT(INOUT), DIMENSION(:):: DY
INTEGER, INTENT(IN)::  N
INTEGER::  I

!$OMP PARALLEL DO PRIVATE(I) SCHEDULE(STATIC)
DO I = 1, N
  DY(I) = DA(I) * DX(I)
ENDDO
!$OMP END PARALLEL DO

END SUBROUTINE SCARC_SCALING_VARIABLE




! ====================================================================================================
! ====================================================================================================
! Bundle of routines for the Smoothed Algebraic Multigrid Method
! ====================================================================================================
! ====================================================================================================

! ----------------------------------------------------------------------------------------------------
! Setup coarsening process (Algebraic Multigrid only)
! Allocate needed workspace for hierarchy of system matrices, prolongation, restriction, etc.
! Note: all used pointers end with either 'F' or 'C' where:
!     'F' corresponds to fine   level NL
!     'C' corresponds to coarse level NL+1
! Determine mesh hierarchy based on smoothed aggregation
! Compute QR-decomposition of nullspace vector in order to determine tentative prolongator 
! Set nullspace for next level and perform Jacobi relaxation to get the final prolongator
! If the maximum allowed level is not yet reached, set dimensions for next coarser level, 
! define its nullspace and perform relaxation to define the respective prolongation matrix
! Define Poisson matrix on coarser level by Galerkin approach: A_coarse = R * A_fine * P
! ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_SAMG
INTEGER :: NL
LOGICAL :: FURTHER_COARSENING_REQUIRED

IF (.NOT.IS_CG_SAMG .AND. .NOT.IS_AMG) RETURN
FURTHER_COARSENING_REQUIRED = .TRUE.

NL = NLEVEL_MIN
!CALL  SCARC_PYTHON_MATRIX(NL, 'A')

COARSENING_LOOP: DO WHILE (FURTHER_COARSENING_REQUIRED)

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) '========================================================'
WRITE(MSG%LU_DEBUG,*) ' SAMG : LEVEL ', NL
WRITE(MSG%LU_DEBUG,*) '========================================================'
#endif

   CALL SCARC_SETUP_AGGREGATION_ORDER               ! Determine the aggregation order among the meshes
   CALL SCARC_EXTRACT_MATRIX_DIAGONAL(NL)           ! Extract matrix diagonal from Poisson matrix
   CALL SCARC_SETUP_STRENGTH_OF_CONNECTION(NL)      ! Determine strength of connection matrix for Poisson matrix
   CALL SCARC_INVERT_MATRIX_DIAGONAL(NL)            ! Store inverted matrix diagonal for later processing
   CALL SCARC_SETUP_AGGREGATION_ZONES(NL)           ! Apply smoothed aggregation heuristic to specify aggregation zones
   CALL SCARC_RELAX_NULLSPACE(NL)                   ! Improve near null space by Jacobi relaxation step
   CALL SCARC_SETUP_ZONE_OPERATOR(NL)               ! Setup final aggregation zones operator
   CALL SCARC_SETUP_PROLONGATION(NL)                ! Setup prolongation matrix based on QR-decomposition
   CALL SCARC_SETUP_NULLSPACE_COARSE(NL)            ! Setup nullspace on coarser level
   CALL SCARC_SETUP_RESTRICTION(NL)                 ! Setup corresponding restriction matrix 
   CALL SCARC_SETUP_AP(NL)             
   CALL SCARC_SETUP_GALERKIN(NL)             
   CALL SCARC_CLEAN_WORKSPACE_SAMG(NL)
      
   NL = NL + 1
   IF (NL == NLEVEL_MAX) FURTHER_COARSENING_REQUIRED = .FALSE.

ENDDO COARSENING_LOOP

END SUBROUTINE SCARC_SETUP_SAMG


! ------------------------------------------------------------------------------------------------------
! Setup order in which aggregation is performed over mesh decomposition
! ------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_AGGREGATION_ORDER
USE SCARC_POINTERS, ONLY: SUB, S
LOGICAL, ALLOCATABLE, DIMENSION(:) :: NOT_AGGREGATED
INTEGER :: NM, NOM, INBR, ICYCLE

SUB => SUBDIVISION

CALL SCARC_ALLOCATE_INT2(SUB%ORDER, 1, NMESHES, 1, NMESHES, NSCARC_INIT_NONE, 'SUB%ORDER')
CALL SCARC_ALLOCATE_LOG1(NOT_AGGREGATED, 1, NMESHES, NSCARC_INIT_TRUE, 'NOT_AGGREGATED')

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

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'SUB%N_CYCLES :', SUB%N_CYCLES
DO NM = 1, NMESHES
   WRITE(MSG%LU_DEBUG,*) 'NM, SUB%N_NEIGHBORS =', SUB%N_NEIGHBORS(NM)
   WRITE(MSG%LU_DEBUG,*) 'NM, SUB%NEIGHBORS   =', SUB%NEIGHBORS(1:SUB%N_NEIGHBORS(NM), NM)
ENDDO
DO ICYCLE = 1, SUB%N_CYCLES
   WRITE(MSG%LU_DEBUG,*) SUB%ORDER(1:NMESHES, ICYCLE)
ENDDO
#endif

DEALLOCATE(NOT_AGGREGATED)

END SUBROUTINE SCARC_SETUP_AGGREGATION_ORDER


! ------------------------------------------------------------------------------------------------------
! Extract diagonal of A and store it in a separate vector for further use
! ------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_EXTRACT_MATRIX_DIAGONAL(NL)
USE SCARC_POINTERS, ONLY: G, A
INTEGER, INTENT(IN) :: NL
INTEGER :: NM, IC, JC, ICOL

MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G
   A => SCARC_POINT_TO_CMATRIX (G, NSCARC_MATRIX_POISSON)

   CALL SCARC_ALLOCATE_REAL1 (G%DIAG, 1, G%NCE, NSCARC_INIT_ZERO, 'G%DIAG')
   DO IC = 1, G%NC
      DO ICOL = A%ROW(IC), A%ROW(IC+1) - 1
         JC = A%COL(ICOL)
         IF (JC == IC) G%DIAG(IC) = A%VAL(ICOL)
      ENDDO
   ENDDO

#ifdef WITH_SCARC_DEBUG
   WRITE(MSG%LU_DEBUG,*) 'EXTRACT_MATRIX_DIAGONAL: DIAG:'
   WRITE(MSG%LU_DEBUG,'(8E12.4)') G%DIAG
#endif

ENDDO MESHES_LOOP

! If there are multiple meshes exchange diagonal matrix on overlapping parts
IF (NMESHES > 1) CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_POISSON_DIAGS, NSCARC_NONE, NL)

END SUBROUTINE SCARC_EXTRACT_MATRIX_DIAGONAL


! ------------------------------------------------------------------------------------------------------
! Invert matrix diagonal which is already stored in DIAG-vector (reuse workspace)
! Scale each matrix element with inverse of diagonal and approximate spectral radius (currently disabled) 
! ------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_INVERT_MATRIX_DIAGONAL(NL)
USE SCARC_POINTERS, ONLY: G
INTEGER, INTENT(IN) :: NL
INTEGER :: NM, IC

MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G
   DO IC = 1, G%NCE
      G%DIAG(IC) = 1.0_EB/G%DIAG(IC)
   ENDDO
#ifdef WITH_SCARC_DEBUG
   WRITE(MSG%LU_DEBUG,*) 'INVERT_MATRIX_DIAGONAL:'
   WRITE(MSG%LU_DEBUG,'(8E12.4)') G%DIAG
#endif
ENDDO MESHES_LOOP

END SUBROUTINE SCARC_INVERT_MATRIX_DIAGONAL


! ------------------------------------------------------------------------------------------------------
! Compute a strength of connection matrix based on symmetric smoothed aggregation heuristic. 
! A nonzero connection A[i,j] is considered strong if:
!
!     abs(A[i,j]) >= theta * sqrt( abs(A[i,i]) * abs(A[j,j]) )
!
! The strength matrix S corresponds to the set of nonzero entries of A that are strong connections
! based on a strength of connection tolerance theta
! ------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_STRENGTH_OF_CONNECTION(NL)
USE SCARC_POINTERS, ONLY: G, A, S, OG
INTEGER, INTENT(IN) :: NL
REAL(EB):: VAL, EPS, THETA = 0.05_EB, SCAL, CVAL_MAX
INTEGER :: NM, NOM, IC, JC, ICOL, IZONE, INBR

MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   
   CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G

   ! Allocate workspaces for coarse points, aggregation zones and strength matrix
   ! Since they all consist of a subsets of A, a conservative bound is to first allocate the same storage 
   ! used for A and reduce it later corresponding to their real sizes
   CALL SCARC_ALLOCATE_INT1 (G%CPOINTS,      1, G%NCE,  NSCARC_INIT_ZERO, 'G%CPOINTS')
   CALL SCARC_ALLOCATE_INT1 (G%ZONES_GLOBAL, 1, G%NCE2, NSCARC_INIT_ZERO, 'G%ZONES_GLOBAL')
   CALL SCARC_ALLOCATE_INT1 (G%ZONES_LOCAL,  1, G%NCE2, NSCARC_INIT_ZERO, 'G%ZONES_LOCAL')

   A => SCARC_POINT_TO_CMATRIX (G, NSCARC_MATRIX_POISSON)
   C => SCARC_POINT_TO_CMATRIX (G, NSCARC_MATRIX_CONNECTION)

   C%N_VAL = A%N_VAL           
   C%N_ROW = A%N_ROW
   CALL SCARC_ALLOCATE_CMATRIX(C, NL, NSCARC_PRECISION_DOUBLE, NSCARC_MATRIX_LIGHT, 'G%CONNECTION')

   DO INBR = 1, SCARC(NM)%N_NEIGHBORS

      NOM = S%NEIGHBORS(INBR)
      CALL SCARC_POINT_TO_OTHER_GRID(NM, NOM, NL)

      OA => SCARC_POINT_TO_OTHER_CMATRIX(OG, NSCARC_MATRIX_POISSON)
      OC => SCARC_POINT_TO_OTHER_CMATRIX(OG, NSCARC_MATRIX_CONNECTION)

      OC%N_VAL = 2*OA%N_VAL                         !TODO: CHECK LENGTH
      OC%N_ROW = OA%N_ROW           
      CALL SCARC_ALLOCATE_CMATRIX(OC, NL, NSCARC_PRECISION_DOUBLE, NSCARC_MATRIX_LIGHT, 'OG%CONNECTION')
      
   ENDDO 

#ifdef WITH_SCARC_DEBUG
   WRITE(MSG%LU_DEBUG,*) '================== STRENGTH0:'
   WRITE(MSG%LU_DEBUG,*) 'G%NC:', G%NC
   WRITE(MSG%LU_DEBUG,*) 'G%NCE:', G%NCE
   WRITE(MSG%LU_DEBUG,*) 'G%NCE2:', G%NCE2
   WRITE(MSG%LU_DEBUG,*) 'G%N_FINE:', G%N_FINE
#endif
   
   ! Check strength-of-connection criterion  
   IZONE = 1
   C%ROW(1) = 1
   DO IC = 1, G%NC
   
      EPS = THETA**2 * ABS(G%DIAG(IC))                  ! EPS = theta**2 * A_ii
      DO ICOL = A%ROW(IC), A%ROW(IC+1) - 1
         JC  = A%COL(ICOL)
         VAL = A%VAL(ICOL)                              ! VAL = A_ij
   
#ifdef WITH_SCARC_DEBUG2
   WRITE(MSG%LU_DEBUG,'(A,4I4,4E12.4)') 'IC, ICOL, JC, IZONE, VAL2, EPS, DIAG(JC), EPS*DIAG(JC):', &
                                         IC, ICOL, JC, IZONE, VAL**2, EPS, ABS(G%DIAG(JC)), EPS*ABS(G%DIAG(JC))
#endif
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
   CALL SCARC_DEBUG_CMATRIX(C, 'CONNECTION', 'STRENGTH OF CNNECTION - PASS 2')
#endif
   
ENDDO MESHES_LOOP
   
! If there are multiple meshes, exchange strength matrix on overlapping parts
IF (NMESHES > 1) THEN
   CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_POISSON_COLS,  NSCARC_MATRIX_CONNECTION, NL)
   CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_POISSON_VALS,  NSCARC_MATRIX_CONNECTION, NL)
ENDIF

END SUBROUTINE SCARC_SETUP_STRENGTH_OF_CONNECTION
 

! ------------------------------------------------------------------------------------------------------
! Typical SAMG aggregation prodecure based on strength of connection matrix
! ------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_AGGREGATION_SAMG(G, C)
TYPE (SCARC_MATRIX_COMPACT_TYPE), POINTER, INTENT(IN) :: C
TYPE (SCARC_GRID_TYPE), POINTER, INTENT(IN) :: G
INTEGER :: IC, ICOL, JC, IZONE, JZONE
LOGICAL :: HAS_NEIGHBORS, HAS_AGGREGATED_NEIGHBORS

! 
! Pass 1 of aggregation:  Setup aggregation zones on internal cells of active mesh
! 
PASS1_LOOP: DO IC = 1, G%NC

   IF (G%ZONES_LOCAL(IC) /= 0) CYCLE                   ! has cell already been aggregated?

   HAS_NEIGHBORS = .FALSE.
   HAS_AGGREGATED_NEIGHBORS = .FALSE.

   DO ICOL = C%ROW(IC), C%ROW(IC+1)-1                  ! are all neighbors free (not already aggregated)?
      JC = C%COL(ICOL)
      IF (IC /= JC .AND. JC <= G%NC) THEN              ! only consider internal cells here
         HAS_NEIGHBORS = .TRUE.
         IF (G%ZONES_LOCAL(JC) /= 0) THEN
            HAS_AGGREGATED_NEIGHBORS = .TRUE.
            EXIT
         ENDIF
      ENDIF
   ENDDO

   IF (.NOT. HAS_NEIGHBORS) THEN                       ! do not aggregate isolated cells
      G%ZONES_LOCAL(IC) = NSCARC_HUGE_INT
   ELSE IF (.NOT. HAS_AGGREGATED_NEIGHBORS) THEN       ! build aggregate of this cell and its neighbors
      G%N_ZONES = G%N_ZONES + 1
      G%ZONES_LOCAL(IC) = G%N_ZONES
      G%CPOINTS(G%N_ZONES) = IC                
      DO ICOL = C%ROW(IC), C%ROW(IC+1)-1 
         JC = C%COL(ICOL)
         IF (JC <= G%NC) G%ZONES_LOCAL(C%COL(ICOL)) = G%N_ZONES
      ENDDO
   ENDIF

#ifdef WITH_SCARC_DEBUG
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
         IF (JC >= G%NC .OR. G%CPOINTS(JZONE)>0) THEN
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
   G%CPOINTS(G%N_ZONES) = IC

   DO ICOL = C%ROW(IC), C%ROW(IC+1)-1
      JC = C%COL(ICOL)
      IF (JC <= G%NC .AND. G%ZONES_LOCAL(JC) == 0) G%ZONES_LOCAL(JC) = G%N_ZONES
   ENDDO
   G%N_ZONES = G%N_ZONES + 1

ENDDO PASS3_LOOP

IF (MINVAL(G%ZONES_LOCAL) < 0) THEN
   WRITE(*,*) 'CAUTION: CELL ',MINLOC(G%ZONES_LOCAL),' HAS NOT BEEN AGGREGATED DURING AGGREGATION'
   CALL MPI_FINALIZE(IERROR)
   STOP
ENDIF
      
#ifdef WITH_SCARC_DEBUG
CALL SCARC_DEBUG_ZONES(G, -1, 1, 'AFTER ACTIVE PASS3')
#endif

END SUBROUTINE SCARC_SETUP_AGGREGATION_SAMG


! ------------------------------------------------------------------------------------------------------
! GMG-like aggregation
! ------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_AGGREGATION_GMG(L, G)
TYPE (SCARC_LEVEL_TYPE), POINTER, INTENT(IN) :: L
TYPE (SCARC_GRID_TYPE), POINTER, INTENT(IN) :: G
INTEGER :: NXM, NYM, NZM, NXD, NYD, NZD, NXI, NYI, NZI
INTEGER :: IX, IY, IZ, IXZ, IYZ, IZZ, IXP, IYP, IZP, IX0, IY0, IZ0, IC
INTEGER, DIMENSION(:), ALLOCATABLE :: OFFX, OFFY, OFFZ
LOGICAL :: BFIRST

NXM = MOD(L%NX,2)
NXD = L%NX/2

IF (TWO_D) THEN
   NYM = 0
   NYD = 1
ELSE
   NYM = MOD(L%NY,2)
   NYD = L%NY/2
ENDIF

NZM = MOD(L%NZ,2)
NZD = L%NZ/2

! Temporarily - to prevent failure of following algorithm
IF ((L%NX < 4) .OR. (.NOT.TWO_D .AND. L%NY < 4) .OR. (L%NZ < 4)) THEN 
   WRITE(*,*) 'Grid dimensions too small für GMG-like aggregation'
   CALL MPI_FINALIZE(IERROR)
   STOP
ENDIF

CALL SCARC_ALLOCATE_INT1 (OFFX, 1, NXD, NSCARC_INIT_ZERO, 'OFFX')
CALL SCARC_ALLOCATE_INT1 (OFFY, 1, NXD, NSCARC_INIT_ZERO, 'OFFY')
CALL SCARC_ALLOCATE_INT1 (OFFZ, 1, NZD, NSCARC_INIT_ZERO, 'OFFZ')

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'SETUP_AGGREGATION_GMG:'
WRITE(MSG%LU_DEBUG,*) 'NXD, NXM =', NXD, NXM
WRITE(MSG%LU_DEBUG,*) 'NYD, NYM =', NYD, NYM
WRITE(MSG%LU_DEBUG,*) 'NZD, NZM =', NZD, NZM
#endif

! if even number of cells in x-direction, use patch length of 2 as in default GMG
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

! if even number of cells in y-direction, use patch length of 2 as in default GMG
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

! if even number of cells in x-direction, use patch length of 2 as in default GMG
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

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'OFFX:'
WRITE(MSG%LU_DEBUG,'(16I6)') OFFX
WRITE(MSG%LU_DEBUG,*) 'OFFY:'
WRITE(MSG%LU_DEBUG,'(16I6)') OFFY
WRITE(MSG%LU_DEBUG,*) 'OFFZ:'
WRITE(MSG%LU_DEBUG,'(16I6)') OFFZ
#endif


DIMENSION_IF: IF (TWO_D) THEN

   IZ0 = 1
   DO IZ = 1, NZD
#ifdef WITH_SCARC_DEBUG
      WRITE(MSG%LU_DEBUG,*) '--------- IZ, IZ0 =', IZ, IZ0
#endif
      IX0 = 1
      DO IX = 1, NXD

#ifdef WITH_SCARC_DEBUG
         WRITE(MSG%LU_DEBUG,*) '                 IX, IX0 =', IX, IX0
#endif
         BFIRST = .TRUE.
         DO IZZ = 0, OFFZ(IZ)-1
            DO IXZ = 0, OFFX(IX)-1
               IXP = IX0 + IXZ
               IZP = IZ0 + IZZ
               IF (IS_UNSTRUCTURED .AND. L%IS_SOLID(IXP, 1, IZP)) CYCLE
               IC = G%CELL_NUMBER(IXP, 1, IZP)
               IF (BFIRST) THEN
                  G%N_ZONES = G%N_ZONES + 1
                  BFIRST = .FALSE. 
                  G%CPOINTS(G%N_ZONES) = IC
               ENDIF
               G%ZONES_LOCAL(IC) = G%N_ZONES
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,'(A,9I6)') 'IX, 1, IZ, IXP, 1, IZP, IC, ZONES_LOCAL(IC), N_ZONES:', &
                               IX, 1, IZ, IXP, 1, IZP, IC, G%ZONES_LOCAL(IC), G%N_ZONES
#endif
            ENDDO
         ENDDO
         IX0 = IX0 + OFFX(IX)
      ENDDO
      IZ0 = IZ0 + OFFZ(IZ)
   ENDDO

ELSE

   IZ0 = 1
   DO IZ = 1, NZD
#ifdef WITH_SCARC_DEBUG
      WRITE(MSG%LU_DEBUG,*) '--------- IZ, IZ0 =', IZ, IZ0
#endif
      IY0 = 1
      DO IY = 1, NYD
#ifdef WITH_SCARC_DEBUG
         WRITE(MSG%LU_DEBUG,*) '             IY, IY0 =', IY, IY0
#endif
         IX0 = 1
         DO IX = 1, NXD

#ifdef WITH_SCARC_DEBUG
            WRITE(MSG%LU_DEBUG,*) '                 IX, IX0 =', IX, IX0
#endif
            BFIRST = .TRUE.
            DO IZZ = 0, OFFZ(IZ)-1
               DO IYZ = 0, OFFY(IY)-1
                  DO IXZ = 0, OFFX(IX)-1
                     IXP = IX0 + IXZ
                     IYP = IY0 + IYZ
                     IZP = IZ0 + IZZ
                     IF (IS_UNSTRUCTURED .AND. L%IS_SOLID(IXP, IYP, IZP)) CYCLE
                     IF (BFIRST) THEN
                        G%N_ZONES = G%N_ZONES + 1
                        BFIRST = .FALSE. 
                     ENDIF
                     IC = G%CELL_NUMBER(IXP, IYP, IZP)
                     G%ZONES_LOCAL(IC) = G%N_ZONES
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,'(A,9I6)') 'IX, IY, IZ, IXP, IYP, IZP, IC, ZONES_LOCAL(IC), N_ZONES:', &
                               IX, IY, IZ, IXP, IYP, IZP, IC, G%ZONES_LOCAL(IC), G%N_ZONES
#endif
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

DEALLOCATE(OFFX)
DEALLOCATE(OFFY)
DEALLOCATE(OFFZ)

#ifdef WITH_SCARC_DEBUG
   CALL SCARC_DEBUG_ZONES(G, IC, 1, 'AFTER GMG-AGG ')
#endif
END SUBROUTINE SCARC_SETUP_AGGREGATION_GMG



! ------------------------------------------------------------------------------------------------------
! Setup aggregation zones for Smoothed Aggregation Algebraic Multigrid Method
! ------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_AGGREGATION_ZONES(NL)
USE SCARC_POINTERS, ONLY: SUB, C, CF, G, GF, GC
INTEGER, INTENT(IN) :: NL
INTEGER :: NM, NM2, ICYCLE, IC, IZL

SUB => SUBDIVISION
MESH_INT = -1

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'TYPE_COARSENING:', TYPE_COARSENING
#endif

COARSENING_TYPE_SELECT: SELECT CASE (TYPE_COARSENING)

   !
   ! ---- Default aggregation procedure for SAMG
   !
   CASE (NSCARC_COARSENING_SAMG)

      CYCLES_LOOP1: DO ICYCLE = 1, SUB%N_CYCLES
   
#ifdef WITH_SCARC_DEBUG
         WRITE(MSG%LU_DEBUG,*) '================================= ICYCLE =', ICYCLE, SUB%N_CYCLES
         WRITE(MSG%LU_DEBUG,*) 'SUB%ORDER(.,ICYCLE):'
         WRITE(MSG%LU_DEBUG,'(20I6)') SUB%ORDER(1:NMESHES,ICYCLE)
#endif
   
         ! First aggregate on active meshes
         DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   
            IF (SUB%ORDER(NM, ICYCLE) == NSCARC_ORDER_ACTIVE) THEN
               CALL SCARC_POINT_TO_GRID(NM, NL)
               C => SCARC_POINT_TO_CMATRIX (G, NSCARC_MATRIX_CONNECTION)
               CALL SCARC_SETUP_AGGREGATION_SAMG(G, C)
               MESH_INT (NM) = G%N_ZONES
#ifdef WITH_SCARC_DEBUG
               WRITE(MSG%LU_DEBUG,*) 'ACTIVE: G%N_ZONES, MESH_INT:', G%N_ZONES, MESH_INT
#endif
            ENDIF
   
         ENDDO
      
         ! Exchange overlapping information of active meshes
         IF (NMESHES > 1)  CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_ZONE_TYPES, NSCARC_ZONES_SEND, NL)
   
         ! Then aggregate on passive meshes (taking into account overlapping aggregate information of active meshes)
         DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
      
            IF (SUB%ORDER (NM, ICYCLE) /= NSCARC_ORDER_ACTIVE) THEN
               CALL SCARC_POINT_TO_GRID(NM, NL)
               C => SCARC_POINT_TO_CMATRIX (G, NSCARC_MATRIX_CONNECTION)
               CALL SCARC_SETUP_AGGREGATION_SAMG(G, C)
               MESH_INT (NM) = G%N_ZONES
#ifdef WITH_SCARC_DEBUG
               WRITE(MSG%LU_DEBUG,*) 'PASSIVE: G%N_ZONES, MESH_INT:', G%N_ZONES, MESH_INT
#endif
            ENDIF
   
         ENDDO
   
         ! Exchange overlapping information of passive meshes
         IF (NMESHES > 1)  CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_ZONE_TYPES, NSCARC_ZONES_SEND, NL)
   
      ENDDO CYCLES_LOOP1
   
#ifdef WITH_SCARC_DEBUG
      WRITE(MSG%LU_DEBUG,*) 'MESH_INT=',MESH_INT
#endif
   
   !
   ! ---- GMG-like aggregation procedure 
   !      In case of even cell numbers this process corresponds to the usual GMG coarsening
   !      in case of uneven cell number in a coordinate direction, on patch with 3 cells is used, the rest with patches of 2
   !
   CASE (NSCARC_COARSENING_GMG)

      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
         CALL SCARC_POINT_TO_GRID(NM, NL)
         CALL SCARC_SETUP_AGGREGATION_GMG(L, G)
         MESH_INT (NM) = G%N_ZONES
      ENDDO
      
      ! Exchange overlapping information 
      IF (NMESHES > 1)  CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_ZONE_TYPES, NSCARC_ZONES_SEND, NL)
   
END SELECT COARSENING_TYPE_SELECT


! Broadcast number of zones of all meshes
IF (N_MPI_PROCESSES>1) &
   CALL MPI_ALLGATHERV(MPI_IN_PLACE, 1, MPI_INTEGER, MESH_INT, COUNTS, DISPLS, MPI_INTEGER, MPI_COMM_WORLD, IERROR)
      

! Only for debugging - Preset single-mesh case corresponding to some multi-mesh test cases
! such that the same aggregation takes place
#ifdef WITH_SCARC_DEBUG
IF (TRIM(CHID) == 'b11111111') THEN
   CALL SCARC_PRESET_B1_CASE
   MESH_INT(1) = 6
ENDIF
IF (TRIM(CHID) == 'b1411111111') THEN
   CALL SCARC_PRESET_B14_CASE
   MESH_INT(1) = 8
ENDIF
IF (TRIM(CHID) == 'b14big') THEN
   CALL SCARC_PRESET_b14big_CASE
   MESH_INT(1) = 24
ENDIF
#endif

!
! Prepare grid dimensions of coarse grid level
!
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
WRITE(MSG%LU_DEBUG,*) 'GC%NC_LOCAL(1:NMESHES) =', GC%NC_LOCAL(1:NMESHES)
WRITE(MSG%LU_DEBUG,*) 'GC%NC_GLOBAL ', GC%NC_GLOBAL
WRITE(MSG%LU_DEBUG,*) 'GC%NC ', GC%NC
WRITE(MSG%LU_DEBUG,*) 'GC%NCE', GC%NCE
WRITE(MSG%LU_DEBUG,*) 'GC%NC_OFFSET(1:NMESHES) ', GC%NC_OFFSET(1:NMESHES)
#endif
   ENDIF                   

! Again only for debugging - Preset single-mesh case corresponding to some multi-mesh test cases
! such that the same aggregation takes place
#ifdef WITH_SCARC_DEBUG
   IF (TRIM(CHID) == 'b11111111') THEN
      GF%NC = 24
      GF%NCE = 24
      GC%NC_LOCAL  = 6
      GC%NC_GLOBAL = 6
      GC%NC_OFFSET = 0
      GC%NC = 6
      GC%NCE = 6
   ENDIF
   IF (TRIM(CHID) == 'b1411111111') THEN
      GF%NC = 36
      GF%NCE = 36
      GC%NC_LOCAL  = 8
      GC%NC_GLOBAL = 8
      GC%NC_OFFSET = 0
      GC%NC = 8
      GC%NCE = 8
   ENDIF
   IF (TRIM(CHID) == 'b14big') THEN
      GF%NC = 108
      GF%NCE = 108
      GF%N_FINE = 108
      GF%N_COARSE = 24
      GC%NC_LOCAL  = 108
      GC%NC_GLOBAL = 108
      GC%NC_OFFSET = 0
      GC%N_FINE = 24
      GC%NC = 24
      GC%NCE = 24
   ENDIF
#endif


   ! Setup mapping from local zones to global zones
   ! TODO Reduce GC%LOCAL_TO_GLOBAL in length 
   CALL SCARC_ALLOCATE_INT1(GC%LOCAL_TO_GLOBAL, 1, GF%NCE, NSCARC_INIT_ZERO, 'G%LOCAL_TO_GLOBAL')
   DO IZL = 1, GC%NC
      GC%LOCAL_TO_GLOBAL(IZL) = IZL + GC%NC_OFFSET(NM)
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'GC%LOCAL_TO_GLOBAL(',IZL,')=',GC%LOCAL_TO_GLOBAL(IZL)
#endif
   ENDDO

   DO IC = 1, GF%NC
      GF%ZONES_GLOBAL(IC) = GF%ZONES_LOCAL(IC) + GC%NC_OFFSET(NM)
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'GF%ZONES_GLOBAL(',IC,')=',GF%ZONES_GLOBAL(IC)
#endif
   ENDDO

ENDDO

! Exchange zones information between meshes
IF (NMESHES > 1)  THEN

   CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_CELL_NEIGHBORS, NSCARC_NONE, NL)
   CALL SCARC_EXTRACT_CELL_OVERLAPS(NL)

   CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_ZONE_NEIGHBORS, NSCARC_NONE, NL)
   CALL SCARC_EXTRACT_ZONE_OVERLAPS(NL)
   CALL SCARC_EXTRACT_ZONE_POINTERS(NL)

ENDIF

! Determine final grid dimensions on coarser level and reduce zone arrays to correct length
DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_MULTIGRID(NM, NL, NL+1)         

   CF => SCARC_POINT_TO_CMATRIX (GF, NSCARC_MATRIX_CONNECTION)

   GC%NCE  = GF%N_ZONES
   GC%NCE2 = GF%N_ZONES
   GF%N_COARSE = GF%N_ZONES

   CALL SCARC_REDUCE_INT1(GF%ZONES_LOCAL, GF%NCE2, 'GC%LOCAL_TO_GLOBAL')
   CALL SCARC_REDUCE_INT1(GF%ZONES_GLOBAL, GF%NCE2, 'GC%LOCAL_TO_GLOBAL')
   CALL SCARC_REDUCE_INT1(GC%LOCAL_TO_GLOBAL, GC%NCE2, 'GC%LOCAL_TO_GLOBAL')

   CALL SCARC_DEALLOCATE_CMATRIX(CF, 'STRENGTH OF CONNECTION')

#ifdef WITH_SCARC_DEBUG
   WRITE(MSG%LU_DEBUG,*) '================== END OF SETUP AGGREGATION_ZONES '
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
   CALL SCARC_BLENDER_ZONES(NM, NL)
#endif

ENDDO

END SUBROUTINE SCARC_SETUP_AGGREGATION_ZONES



! ------------------------------------------------------------------------------------------------------
! Remove workspace on level 'NL' which will no longer be needed after matrix setup
! ------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_CLEAN_WORKSPACE_SYSTEM(NL)
INTEGER, INTENT(IN) :: NL
INTEGER:: NM

! TODO: deallocate arrays which are no longer used
DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   CALL SCARC_POINT_TO_GRID (NM, NL)
ENDDO

END SUBROUTINE SCARC_CLEAN_WORKSPACE_SYSTEM



! ------------------------------------------------------------------------------------------------------
! Remove workspace on level 'NL' which will no longer be needed in SAMG method
! ------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_CLEAN_WORKSPACE_SAMG(NL)
INTEGER, INTENT(IN) :: NL
INTEGER:: NM

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   CALL SCARC_POINT_TO_GRID (NM, NL)
   DEALLOCATE (G%CPOINTS, STAT = IERROR)
   DEALLOCATE (G%AUX1, STAT = IERROR)
   DEALLOCATE (G%AUX2, STAT = IERROR)
   DEALLOCATE (G%RR, STAT = IERROR)
   DEALLOCATE (G%QQ, STAT = IERROR)
ENDDO

END SUBROUTINE SCARC_CLEAN_WORKSPACE_SAMG



! ------------------------------------------------------------------------------------------------------
! Perform AMG Jacobi :.. x = x - omega D^{-1} (Ax-b)
! Near-null space vector is given in vector G%NULLSPACE --> corresponds to x
! vector b is zero
! ------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_RELAX_NULLSPACE(NL)
USE SCARC_POINTERS, ONLY: L, G, A, MG 
INTEGER, INTENT(IN) :: NL
INTEGER :: IC, JC, ICOL, IDIAG, NM

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G
   A => SCARC_POINT_TO_CMATRIX (G, NSCARC_MATRIX_POISSON)
   MG => L%MG

#ifdef WITH_SCARC_DEBUG
   CALL SCARC_DEBUG_CMATRIX (A, 'POISSON', 'START RELAX NULLSPACE')
#endif

   CALL SCARC_ALLOCATE_REAL1 (G%NULLSPACE, 1, G%NCE2,   NSCARC_INIT_ONE,  'G%NULLSPACE')
   CALL SCARC_ALLOCATE_REAL1 (G%AUX1,      1, G%NCE2+1, NSCARC_INIT_ZERO, 'G%AUX1')
   CALL SCARC_ALLOCATE_REAL1 (G%AUX2,      1, G%NCE2+1, NSCARC_INIT_ZERO, 'G%AUX2')

   MG%OMEGA = 1.0_EB                                   ! currently used default
   MG%OMEGA = MG%OMEGA/MG%APPROX_SPECTRAL_RADIUS
   G%AUX1 = G%NULLSPACE
   
   ! Compute defect to near-null-space vector: d = Ax - b, in this case
   !    'x' corresponds to nullspace vector consisting of only '1'-entries 
   !    'b' corresponds to zero vector
   DO IC = 1, G%NC
   
      IDIAG = A%ROW(IC)                                                
      JC = A%COL(IDIAG)
      G%AUX2(IC) = A%VAL(IDIAG)* G%AUX1(JC)
#ifdef WITH_SCARC_DEBUG2
   WRITE(MSG%LU_DEBUG,'(A, 3I4,3E12.4)') 'NULLSPACE: IC; IDIAG, JC, AUX1=', &
                                          IC, IDIAG, JC, A%VAL(IDIAG), G%AUX1(JC), G%AUX2(IC)
#endif
   
      DO ICOL = A%ROW(IC)+1, A%ROW(IC+1)-1                          
         JC = A%COL(ICOL)
         G%AUX2(IC) =  G%AUX2(IC) + A%VAL(ICOL) * G%AUX1(JC)
#ifdef WITH_SCARC_DEBUG2
   WRITE(MSG%LU_DEBUG,'(A, 3I4,3E12.4)') 'NULLSPACE: IC; IDIAG, JC, AUX1=', &
                                          IC, IDIAG, JC, A%VAL(IDIAG), G%AUX1(JC), G%AUX2(IC)
#endif
      ENDDO
   
   ENDDO
   
#ifdef WITH_SCARC_DEBUG
   WRITE(MSG%LU_DEBUG,*) 'RELAX_NULLSPACE:A: AUX2:'
   WRITE(MSG%LU_DEBUG,'(6E12.4)') G%AUX2(1:G%NCE2)
#endif
   
   ! scale it by omega and inverse of diagonal:    d = omega D^{-1} d
   DO IC = 1, G%NC
      G%AUX2(IC) = MG%OMEGA * G%DIAG(IC) * G%AUX2(IC) 
#ifdef WITH_SCARC_DEBUG
   WRITE(MSG%LU_DEBUG,'(A, I4,3E12.4)') 'DIAG: IC =', IC, MG%OMEGA, G%DIAG(IC), G%AUX2(IC)
#endif
   ENDDO
   
#ifdef WITH_SCARC_DEBUG
   WRITE(MSG%LU_DEBUG,*) 'RELAX_NULLSPACE:B: AUX2:'
   WRITE(MSG%LU_DEBUG,'(6E12.4)') G%AUX2
#endif
   
   ! get new iterate:   x = x - d
#ifdef WITH_MKL
   CALL DAXPBY(G%NC, -1.0_EB, G%AUX2, 1, 1.0_EB, G%NULLSPACE, 1)
#else
   CALL SCARC_DAXPY_CONSTANT_DOUBLE(G%NC, -1.0_EB, G%AUX2, 1.0_EB, G%NULLSPACE)
#endif
   
#ifdef WITH_SCARC_DEBUG
   WRITE(MSG%LU_DEBUG,*) '======================================================='
   WRITE(MSG%LU_DEBUG,*) 'RELAX_NULLSPACE: AUX1: '
   WRITE(MSG%LU_DEBUG,'(6E12.4)') G%AUX1(1: G%NC)
   WRITE(MSG%LU_DEBUG,*) '---------------------------------'
   WRITE(MSG%LU_DEBUG,'(6E12.4)') G%AUX1(G%NC+1:G%NCE2)
   WRITE(MSG%LU_DEBUG,*) '======================================================='
   WRITE(MSG%LU_DEBUG,*) 'RELAX_NULLSPACE: AUX2: '
   WRITE(MSG%LU_DEBUG,'(6E12.4)') G%AUX2(1: G%NC)
   WRITE(MSG%LU_DEBUG,*) '---------------------------------'
   WRITE(MSG%LU_DEBUG,'(6E12.4)') G%AUX2(G%NC+1:G%NCE2)
   WRITE(MSG%LU_DEBUG,*) '======================================================='
   WRITE(MSG%LU_DEBUG,*) 'RELAX_NULLSPACE: NULLSPACE: '
   WRITE(MSG%LU_DEBUG,'(6E12.4)') G%NULLSPACE(1: G%NC)
   WRITE(MSG%LU_DEBUG,*) '---------------------------------'
   WRITE(MSG%LU_DEBUG,'(6E12.4)') G%NULLSPACE(G%NC+1:G%NCE2)
#endif

ENDDO

IF (NMESHES > 1) CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_NULLSPACE, NSCARC_NONE, NL)

END SUBROUTINE SCARC_RELAX_NULLSPACE

! ------------------------------------------------------------------------------------------------------
! Setup basic structure of prolongation matrix
! This concerns the setting of the number of rows and the column pointers
! The values are still missing and are set in a later step
! ------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_ZONE_OPERATOR(NL)
USE SCARC_POINTERS, ONLY: GF, GC, AF, ZF
INTEGER, INTENT(IN) :: NL
INTEGER :: NM, ICC, ICCL, ICCG, IC, IP, N_ROW, N_VAL
LOGICAL :: IS_INCLUDED 

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_MULTIGRID (NM, NL, NL+1)                       ! Sets grid pointer GF and GC

   AF => SCARC_POINT_TO_CMATRIX (GF, NSCARC_MATRIX_POISSON)
   ZF => SCARC_POINT_TO_CMATRIX (GF, NSCARC_MATRIX_ZONES)

   ! first use very conservative bounds for the size of the zones operator matrix 
   ! reduce it later once the real size is known
   ZF%N_VAL = AF%N_VAL                  
   ZF%N_ROW = AF%N_VAL                  
   CALL SCARC_ALLOCATE_CMATRIX(ZF, NL, NSCARC_PRECISION_DOUBLE, NSCARC_MATRIX_MINIMAL, 'GF%ZONES')
  
   ! Again conservative upper bound for length - to be reduced later 
   CALL SCARC_ALLOCATE_INT1(GF%ZONES_LOCAL, 1, GF%NCE2, NSCARC_INIT_ZERO, 'GF%ZONES_LOCAL')


#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) '=================== SETUP_ZONE_OPERATOR:'
WRITE(MSG%LU_DEBUG,*) 'ZF%ROW:', ZF%N_ROW, ZF%N_VAL, SIZE(ZF%ROW), SIZE(ZF%VAL)
WRITE(MSG%LU_DEBUG,*) 'ZONES_LOCAL:'
WRITE(MSG%LU_DEBUG,'(6I6)') GF%ZONES_LOCAL(1:GF%NCE2)
WRITE(MSG%LU_DEBUG,*) 'ZONES_GLOBAL:'
WRITE(MSG%LU_DEBUG,'(6I6)') GF%ZONES_GLOBAL(1:GF%NCE2)
WRITE(MSG%LU_DEBUG,*) '=================== SETUP_ZONE_OPERATOR: FINE'
WRITE(MSG%LU_DEBUG,*) 'SETUP_ZONE_OPERATOR: GF%N_FINE:', GF%N_FINE
WRITE(MSG%LU_DEBUG,*) 'SETUP_ZONE_OPERATOR: GF%N_COARSE:', GF%N_COARSE
WRITE(MSG%LU_DEBUG,*) 'SETUP_ZONE_OPERATOR: GF%NC:', GF%NC
WRITE(MSG%LU_DEBUG,*) 'SETUP_ZONE_OPERATOR: GF%NCE:', GF%NCE
WRITE(MSG%LU_DEBUG,*) 'SETUP_ZONE_OPERATOR: GF%NCE2:', GF%NCE2
WRITE(MSG%LU_DEBUG,*) 'SETUP_ZONE_OPERATOR: GF%LOCAL_TO_GLOBAL:'
WRITE(MSG%LU_DEBUG,'(6I6)') GF%LOCAL_TO_GLOBAL(1:GF%NCE2)
WRITE(MSG%LU_DEBUG,*) '=================== SETUP_ZONE_OPERATOR: COARSE'
WRITE(MSG%LU_DEBUG,*) 'SETUP_ZONE_OPERATOR: GC%NC_LOCAL:', GC%NC_LOCAL
WRITE(MSG%LU_DEBUG,*) 'SETUP_ZONE_OPERATOR: GC%NC:', GC%NC
WRITE(MSG%LU_DEBUG,*) 'SETUP_ZONE_OPERATOR: GC%NCE:', GC%NCE
WRITE(MSG%LU_DEBUG,*) 'SETUP_ZONE_OPERATOR: GC%NCE2:', GC%NCE2
WRITE(MSG%LU_DEBUG,*) 'SETUP_ZONE_OPERATOR: GC%NC_OFFSET:', GC%NC_OFFSET
WRITE(MSG%LU_DEBUG,*) 'SETUP_ZONE_OPERATOR: GC%LOCAL_TO_GLOBAL:'
WRITE(MSG%LU_DEBUG,'(6I6)') GC%LOCAL_TO_GLOBAL(1:GC%NCE2)
#endif

   ! Based on global zone numbers determine local zone numbers within mesh
   IP = 1
   ICC = 1
   ZF%ROW(ICC) = 1
   DO ICCL = 1, GC%NCE2

      ICCG = GC%LOCAL_TO_GLOBAL(ICCL)
#ifdef WITH_SCARC_DEBUG2
WRITE(MSG%LU_DEBUG,'(A,3I4)') '=================>  ICCL, ICCG: ', ICCL, ICCG
#endif
      IS_INCLUDED = .FALSE.
      DO IC = 1, GF%NCE2
         IF (GF%ZONES_LOCAL(IC) /= ICCL) CYCLE
         IS_INCLUDED = .TRUE.
         ZF%COL(IP)  = IC
         !ZF%COLG(IP) = GF%LOCAL_TO_GLOBAL(IC)

#ifdef WITH_SCARC_DEBUG2
WRITE(MSG%LU_DEBUG,'(A,4I4)') '     bingo, IC, ICC, ICCL, ICCG:', IC, ICC, ICCL, ICCG
#endif
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

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) ' ICC =',ICC
WRITE(MSG%LU_DEBUG,*) ' IP =',IP
WRITE(MSG%LU_DEBUG,*) ' N_ROW =', N_ROW
WRITE(MSG%LU_DEBUG,*) ' N_VAL =', N_VAL
WRITE(MSG%LU_DEBUG,*) 'ZF%N_ROW =',ZF%N_ROW
WRITE(MSG%LU_DEBUG,*) 'ZF%N_VAL =',ZF%N_VAL
CALL SCARC_DEBUG_CMATRIX (ZF, 'ZONES','AFTER SETUP AGGREGATION ZONES 1 ')
#endif

   CALL SCARC_REDUCE_CMATRIX(ZF, 'ZF')

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'GF%NCE=',GF%NCE
WRITE(MSG%LU_DEBUG,*) 'GF%N_COARSE=',GF%N_COARSE
WRITE(MSG%LU_DEBUG,*) 'GC%NCE=',GC%NCE
WRITE(MSG%LU_DEBUG,*) 'ZF%N_ROW=',ZF%N_ROW
WRITE(MSG%LU_DEBUG,*) 'ZF%N_VAL=',ZF%N_VAL
CALL SCARC_DEBUG_CMATRIX (ZF, 'ZONES','AFTER SETUP AGGREGATION ZONES 2')
#endif

ENDDO

END SUBROUTINE SCARC_SETUP_ZONE_OPERATOR

! -------------------------------------------------------------------------------------------
! Extract overlapping cell information (including second layers)
! -------------------------------------------------------------------------------------------
SUBROUTINE SCARC_EXTRACT_CELL_OVERLAPS(NL)
USE SCARC_POINTERS, ONLY : G, OL
INTEGER, INTENT(IN) :: NL
INTEGER :: NM, INBR, IOR0, NOM, ICG, ICE, ICE2, ICW, ICW2, IOFF


DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G

   CALL SCARC_ALLOCATE_INT1 (G%ICE_TO_ECELL, G%NC+1, G%NCE2, NSCARC_INIT_ZERO, 'G%ICE_TO_ECELL')
   CALL SCARC_ALLOCATE_INT1 (G%ICE_TO_ICELL, G%NC+1, G%NCE2, NSCARC_INIT_ZERO, 'G%ICE_TO_ICELL')
   CALL SCARC_ALLOCATE_INT1 (G%ICE_TO_OCELL, G%NC+1, G%NCE2, NSCARC_INIT_ZERO, 'G%ICE_TO_OCELL')
   CALL SCARC_ALLOCATE_INT1 (G%ICE_TO_GCELL, G%NC+1, G%NCE2, NSCARC_INIT_ZERO, 'G%ICE_TO_GCELL')

#ifdef WITH_SCARC_DEBUG
   WRITE(MSG%LU_DEBUG,*) 'ALLOCATING ICE_TO_OCELL IN LENGTH ', G%NC+1, G%NCE, G%NCE2
#endif

   DO INBR = 1, SCARC(NM)%N_NEIGHBORS

      NOM = S%NEIGHBORS(INBR)
      CALL SCARC_POINT_TO_OTHER_GRID(NM, NOM, NL)
    
      DO IOR0 = -3, 3
         IF (OL%GHOST_LASTW(IOR0) == 0) CYCLE
         DO ICG = OL%GHOST_FIRSTE(IOR0), OL%GHOST_LASTE(IOR0)

            ICE  = OG%ICG_TO_ICE(ICG, 1)
            ICW  = OG%ICG_TO_ICW(ICG, 1)

            G%ICE_TO_ECELL(ICE) = ICE
            G%ICE_TO_ICELL(ICE) = ICW
            G%ICE_TO_OCELL(ICE) = OG%ICG_TO_OCELL(ICG)
            G%ICE_TO_GCELL(ICE) = OG%ICG_TO_GCELL(ICG)

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,'(A,8I6)') 'CELL_OVERLAPS: INBR, NOM, IOR0, ICG, ICE: ', INBR, NOM, IOR0, ICG, ICW, ICE, &
                               G%ICE_TO_OCELL(ICE), G%ICE_TO_GCELL(ICE)
#endif
            
         ENDDO
         IOFF = OL%GHOST_LASTW(IOR0)
         DO ICG = OL%GHOST_FIRSTE(IOR0), OL%GHOST_LASTE(IOR0)

            ICE2 = OG%ICG_TO_ICE(ICG, 2)
            ICW2 = OG%ICG_TO_ICW(ICG, 2)

            G%ICE_TO_ECELL(ICE2) = ICE2
            G%ICE_TO_ICELL(ICE2) = ICW2
            G%ICE_TO_OCELL(ICE2) = OG%ICG_TO_OCELL(ICG + IOFF)
            G%ICE_TO_GCELL(ICE2) = OG%ICG_TO_GCELL(ICG + IOFF)

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,'(A,8I6)') 'CELL_OVERLAPS: INBR, NOM, IOR0, ICG, ICE2: ', INBR, NOM, IOR0, ICG, ICW2, ICE2, &
                               G%ICE_TO_OCELL(ICE2), G%ICE_TO_GCELL(ICE2)
#endif
            
         ENDDO
      ENDDO
   ENDDO
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) '===================== EXTRACT_CELL_OVERLAPS:'
WRITE(MSG%LU_DEBUG,*) 'G%ICE_TO_ECELL:'
WRITE(MSG%LU_DEBUG,'(10I6)') G%ICE_TO_ECELL(G%NC+1:G%NCE2)
WRITE(MSG%LU_DEBUG,*) 'G%ICE_TO_ICELL:'
WRITE(MSG%LU_DEBUG,'(10I6)') G%ICE_TO_ICELL(G%NC+1:G%NCE2)
WRITE(MSG%LU_DEBUG,*) 'G%ICE_TO_GCELL:'
WRITE(MSG%LU_DEBUG,'(10I6)') G%ICE_TO_GCELL(G%NC+1:G%NCE2)
WRITE(MSG%LU_DEBUG,*) 'G%ICE_TO_OCELL:'
WRITE(MSG%LU_DEBUG,'(10I6)') G%ICE_TO_OCELL(G%NC+1:G%NCE2)
#endif

ENDDO

END SUBROUTINE SCARC_EXTRACT_CELL_OVERLAPS

! -------------------------------------------------------------------------------------------
! Extract overlapping zone information (including second layers)
! -------------------------------------------------------------------------------------------
SUBROUTINE SCARC_EXTRACT_ZONE_OVERLAPS(NL)
USE SCARC_POINTERS, ONLY : S, GC, GF, OLF, OGF
INTEGER, INTENT(IN) :: NL
INTEGER :: NM, INBR, IOR0, NOM, IZ, ICG, ICE, ICE2, ICW, ICW2, IFOUND, IZL_CURRENT, ZONE_PREV
INTEGER :: IZL1, IZL2, IZG1, IZG2

MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_MULTIGRID (NM, NL, NL+1)    

   CALL SCARC_ALLOCATE_INT1 (GF%ICE_TO_OZONE, GF%NC+1, GF%NCE2, NSCARC_INIT_ZERO, 'GF%ICE_TO_OZONE')
   CALL SCARC_ALLOCATE_INT1 (GF%ICE_TO_GZONE, GF%NC+1, GF%NCE2, NSCARC_INIT_ZERO, 'GF%ICE_TO_GZONE')
   CALL SCARC_ALLOCATE_INT1 (GF%ICE_TO_IZONE, GF%NC+1, GF%NCE2, NSCARC_INIT_ZERO, 'GF%ICE_TO_IZONE')
   CALL SCARC_ALLOCATE_INT1 (GF%ICE_TO_EZONE, GF%NC+1, GF%NCE2, NSCARC_INIT_ZERO, 'GF%ICE_TO_EZONE')

   ! Clear overlapping parts and fill them with the recently exchanged data
   GF%ZONES_LOCAL (GF%NC+1: GF%NCE2) = 0
   GF%ZONES_GLOBAL(GF%NC+1: GF%NCE2) = 0

   IZL_CURRENT = GF%N_ZONES + 1
   NEIGHBORS_LOOP: DO INBR = 1, SCARC(NM)%N_NEIGHBORS

      NOM = S%NEIGHBORS(INBR)
      CALL SCARC_POINT_TO_OTHER_MULTIGRID(NM, NOM, NL, NL+1)

      ZONE_PREV = -1
      DO IOR0 = -3, 3
         IF (OLF%GHOST_LASTW(IOR0) == 0) CYCLE

         IZ = 0
         DO ICG = OLF%GHOST_FIRSTE(IOR0), OLF%GHOST_LASTE(IOR0)

            ICE  = OGF%ICG_TO_ICE(ICG, 1)
            ICE2 = OGF%ICG_TO_ICE(ICG, 2)

            ICW  = OGF%ICG_TO_ICW(ICG, 1)
            ICW2 = OGF%ICG_TO_ICW(ICG, 2)

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,'(A,6I6)') 'ZONE_OVERLAPS: INBR, NOM, IOR0, ICG, ICE, ICE2:', &
                                              INBR, NOM, IOR0, ICG, ICE, ICE2
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
WRITE(MSG%LU_DEBUG,'(A,7I6)') 'A: ICG, IFOUND, N_ZONES IZL1, ICE, LOCAL(ICE), GLOBAL(ICE):', &
                       ICG, IFOUND, GF%N_ZONES, IZL1, ICE, GF%ZONES_LOCAL(ICE), GF%ZONES_GLOBAL(ICE)
#endif

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
WRITE(MSG%LU_DEBUG,'(A,7I6)') 'B: ICG, IFOUND, N_ZONES IZL2, ICE, LOCAL(ICE), GLOBAL(ICE):', &
                       ICG, IFOUND, GF%N_ZONES, IZL2, ICE2, GF%ZONES_LOCAL(ICE2), GF%ZONES_GLOBAL(ICE2)
#endif

            GF%ICE_TO_EZONE(ICE)  = GF%ZONES_LOCAL(ICE)
            GF%ICE_TO_EZONE(ICE2) = GF%ZONES_LOCAL(ICE2)

            GF%ICE_TO_IZONE(ICE)  = GF%ZONES_LOCAL(ICW)
            GF%ICE_TO_IZONE(ICE2) = GF%ZONES_LOCAL(ICW2)

            GF%ICE_TO_OZONE(ICE)  = OGF%ICG_TO_OZONE(ICG)
            GF%ICE_TO_OZONE(ICE2) = OGF%ICG_TO_OZONE(ICG + OGF%NCG)

            GF%ICE_TO_GZONE(ICE)  = OGF%ICG_TO_GZONE(ICG)
            GF%ICE_TO_GZONE(ICE2) = OGF%ICG_TO_GZONE(ICG + OGF%NCG)

         ENDDO
      ENDDO
   ENDDO NEIGHBORS_LOOP

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) '===================== EXTRACT_ZONE_OVERLAPS:'
WRITE(MSG%LU_DEBUG,*) 'GF%ICE_TO_IZONE:'
WRITE(MSG%LU_DEBUG,'(10I6)') GF%ICE_TO_IZONE(GF%NC+1:GF%NCE2)
WRITE(MSG%LU_DEBUG,*) 'GF%ICE_TO_EZONE:'
WRITE(MSG%LU_DEBUG,'(10I6)') GF%ICE_TO_EZONE(GF%NC+1:GF%NCE2)
WRITE(MSG%LU_DEBUG,*) 'GF%ICE_TO_GZONE:'
WRITE(MSG%LU_DEBUG,'(10I6)') GF%ICE_TO_GZONE(GF%NC+1:GF%NCE2)
WRITE(MSG%LU_DEBUG,*) 'GF%ICE_TO_OZONE:'
WRITE(MSG%LU_DEBUG,'(10I6)') GF%ICE_TO_OZONE(GF%NC+1:GF%NCE2)
#endif

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'FINAL: GC%NCE=',GC%NCE
   CALL SCARC_DEBUG_ZONES(GF, -1, 1, 'AFTER EXTRACT_ZONES')
   CALL SCARC_DEBUG_ZONES(GF, -1, 2, 'AFTER EXTRAXT_ZONES')
#endif

ENDDO MESHES_LOOP

END SUBROUTINE SCARC_EXTRACT_ZONE_OVERLAPS


! -------------------------------------------------------------------------------------------
! Setup pointers for overlapping zones between levels NL and NL + 1
! -------------------------------------------------------------------------------------------
SUBROUTINE SCARC_EXTRACT_ZONE_POINTERS(NL)
USE SCARC_POINTERS, ONLY : S, GF, OLF, OGF, OGC
INTEGER, INTENT(IN) :: NL
INTEGER :: NM, INBR, IOR0, NOM, IZ, ICW1, ICW2, ICE1, ICE2, IZL1, IZL2, ICG, IZW, IZE, NCGE_TOTAL = 0

MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_MULTIGRID (NM, NL, NL+1)    
   
   NEIGHBORS_LOOP: DO INBR = 1, SCARC(NM)%N_NEIGHBORS

      NOM = S%NEIGHBORS(INBR)
      CALL SCARC_POINT_TO_OTHER_MULTIGRID(NM, NOM, NL, NL+1)
    
      CALL SCARC_ALLOCATE_INT1(OGF%ICG_TO_IZONE, 1, 2*OGF%NCG, NSCARC_INIT_ZERO, 'OGF%ICG_TO_IZONE')
      CALL SCARC_ALLOCATE_INT1(OGF%ICG_TO_EZONE, 1, 2*OGF%NCG, NSCARC_INIT_ZERO, 'OGF%ICG_TO_EZONE')

      ORIENTATION_LOOP: DO IOR0 = -3, 3
         IF (OLF%GHOST_LASTW(IOR0) == 0) CYCLE

         IZ  = 0
         INTERNAL_ZONES_LOOP: DO ICG = OLF%GHOST_FIRSTE(IOR0), OLF%GHOST_LASTE(IOR0)

            ICW1 = OGF%ICG_TO_ICW(ICG, 1)
            IZL1 = GF%ZONES_LOCAL(ICW1)
            IF (FINDLOC(OGF%ICG_TO_IZONE, VALUE = IZL1, DIM = 1) == 0) THEN
               IZ = IZ + 1
               OGF%ICG_TO_IZONE(IZ) = IZL1
            ENDIF

#ifdef WITH_SCARC_DEBUG2
WRITE(MSG%LU_DEBUG,'(A, 7I6)') 'INBR, NOM, IOR0, ICG, ICW1, IZL1, IZ :', &
                                INBR, NOM, IOR0, ICG, ICW1, IZL1, IZ
#endif

            ICW2 = OGF%ICG_TO_ICW(ICG, 2)
            IZL2 = GF%ZONES_LOCAL(ICW2)
            IF (FINDLOC(OGF%ICG_TO_IZONE, VALUE = IZL2, DIM = 1) == 0) THEN
               IZ = IZ + 1
               OGF%ICG_TO_IZONE(IZ) = IZL2
            ENDIF

#ifdef WITH_SCARC_DEBUG2
WRITE(MSG%LU_DEBUG,'(A, 7I6)') 'INBR, NOM, IOR0, ICG, ICW2, IZL2, IZ :', &
                                INBR, NOM, IOR0, ICG, ICW2, IZL2, IZ
#endif

         ENDDO INTERNAL_ZONES_LOOP
         OGF%NCGI = IZ
         CALL SCARC_REDUCE_INT1(OGF%ICG_TO_IZONE, OGF%NCGI, 'OGF%ICG_TO_IZONE')

         ! First allocate in fine cell related length
         CALL SCARC_ALLOCATE_INT2(OGC%ICG_TO_ICW, 1, OGF%NCGI, 1, 1, NSCARC_INIT_ZERO, 'OGF%ICG_TO_ICW')
         
         IZ = 0
         DO ICG = 1, OGF%NCGI
            IZW = OGF%ICG_TO_IZONE(ICG)
            IF (FINDLOC(OGC%ICG_TO_ICW(1:OGF%NCGI,1), VALUE = IZW, DIM = 1) == 0) THEN
               IZ = IZ + 1
               OGC%ICG_TO_ICW(IZ, 1) = IZW
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'OGC%ICG_TO_ICW(',IZ,',1)=',IZW
#endif
            ENDIF
         ENDDO
         OGC%NCG  = IZ
         OGC%NCGI = IZ

         OLC%GHOST_FIRSTW(IOR0) = 1
         OLC%GHOST_LASTW(IOR0)  = IZ
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'SUSI: GHOST_FIRSTW(',IOR0,')=', OLC%GHOST_FIRSTW(IOR0)
WRITE(MSG%LU_DEBUG,*) 'SUSI: GHOST_LASTW(',IOR0,')=', OLC%GHOST_LASTW(IOR0)
#endif
         ! Then reduce to zone related length 
         CALL SCARC_REDUCE_INT2(OGC%ICG_TO_ICW, OGC%NCGI, 1, 'OGC%ICG_TO_ICW')


         IZ  = 0 
         EXTERNAL_ZONES_LOOP: DO ICG = OLF%GHOST_FIRSTE(IOR0), OLF%GHOST_LASTE(IOR0)

            ICE1 = OGF%ICG_TO_ICE(ICG, 1)
            IZL1 = GF%ZONES_LOCAL(ICE1)
            IF (FINDLOC(OGF%ICG_TO_EZONE, VALUE = IZL1, DIM = 1) == 0) THEN
               IZ = IZ + 1
               OGF%ICG_TO_EZONE(IZ) = IZL1
            ENDIF

#ifdef WITH_SCARC_DEBUG2
WRITE(MSG%LU_DEBUG,'(A, 7I6)') 'INBR, NOM, IOR0, ICG, ICE1, IZL1, IZ :', &
                                INBR, NOM, IOR0, ICG, ICE1, IZL1, IZ
#endif

            ICE2 = OGF%ICG_TO_ICE(ICG, 2)
            IZL2 = GF%ZONES_LOCAL(ICE2)
            IF (FINDLOC(OGF%ICG_TO_EZONE, VALUE = IZL2, DIM = 1) == 0) THEN
               IZ = IZ + 1
               OGF%ICG_TO_EZONE(IZ) = IZL2
            ENDIF

#ifdef WITH_SCARC_DEBUG2
WRITE(MSG%LU_DEBUG,'(A, 7I6)') 'INBR, NOM, IOR0, ICG, ICE2, IZL2, IZ :', &
                                INBR, NOM, IOR0, ICG, ICE2, IZL2, IZ
#endif

         ENDDO EXTERNAL_ZONES_LOOP
         OGF%NCGE = IZ
         CALL SCARC_REDUCE_INT1(OGF%ICG_TO_EZONE, OGF%NCGE, 'OGF%ICG_TO_IZONE')

         ! First allocate in fine cell related length
         CALL SCARC_ALLOCATE_INT2 (OGC%ICG_TO_ICE, 1, 2*OGF%NCGE, 1, 1, NSCARC_INIT_ZERO, 'OGC%ICG_TO_ICE')
         
         IZ = 0
         DO ICG = 1, OGF%NCGE
            IZE = OGF%ICG_TO_EZONE(ICG)
            IF (FINDLOC(OGC%ICG_TO_ICE(1:OGF%NCGE,1), VALUE = IZE, DIM = 1) == 0) THEN
               IZ = IZ + 1
               OGC%ICG_TO_ICE(IZ, 1) = IZE
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'OGC%ICG_TO_ICE(',IZ,',1)=',IZE
#endif
            ENDIF
         ENDDO
         OGC%NCGE = IZ
         NCGE_TOTAL = NCGE_TOTAL + OGC%NCGE

         OLC%GHOST_FIRSTE(IOR0) = 1
         OLC%GHOST_LASTE(IOR0)  = IZ
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'SUSI: GHOST_FIRSTE(',IOR0,')=', OLC%GHOST_FIRSTE(IOR0)
WRITE(MSG%LU_DEBUG,*) 'SUSI: GHOST_LASTE(',IOR0,')=', OLC%GHOST_LASTE(IOR0)
#endif
         ! Then reduce to zone related length 
         CALL SCARC_REDUCE_INT2(OGC%ICG_TO_ICE, OGC%NCGE, 1, 'OGC%ICG_TO_ICE')


         GC%N_STENCIL_MAX = 20                  ! TODO: ONLY TEMPORARILY
         OGC%NLEN_BUFFER_LAYER1  = MAX(OGC%NCGI, OGC%NCGE)
         OGC%NLEN_BUFFER_LAYER2  = OGC%NLEN_BUFFER_LAYER1 * 2
         OGC%NLEN_BUFFER_LAYER4  = OGC%NLEN_BUFFER_LAYER1 * 4
         OGC%NLEN_BUFFER_STENCIL = OGC%NLEN_BUFFER_LAYER1 * GC%N_STENCIL_MAX
         OGC%NLEN_BUFFER_FULL    = OGC%NLEN_BUFFER_LAYER1 * GC%N_STENCIL_MAX * 2
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'FULL LENGTH OF COARSE BUFFER :', NM, NOM, OGC%NLEN_BUFFER_FULL
#endif
#ifdef WITH_SCARC_VERBISE
WRITE(MSG%LU_VERBOSE,*) 'FULL LENGTH OF COARSE BUFFER :', NM, NOM, OGC%NLEN_BUFFER_FULL
#endif

      ENDDO ORIENTATION_LOOP

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) '===================== EXTRACT_ZONE_POINTERS:'
WRITE(MSG%LU_DEBUG,'(A,3I6)') 'ZONE_POINTERS INBR, NOM, IOR0:', INBR, NOM, IOR0
WRITE(MSG%LU_DEBUG,'(A,2I6)') 'EXCHANGE LENGTH WITH NEIGHBOR ',NOM, OGC%NLEN_BUFFER_LAYER1
WRITE(MSG%LU_DEBUG,*) '------------------------------------------------'
WRITE(MSG%LU_DEBUG,'(A,I6)') 'ICG_TO_IZONE(.): OGF%NCGI', OGF%NCGI
WRITE(MSG%LU_DEBUG,'(10I6)')  (OGF%ICG_TO_IZONE(IZ), IZ=1, OGF%NCGI)
WRITE(MSG%LU_DEBUG,*) '------------------------------------------------'
WRITE(MSG%LU_DEBUG,'(A,I6)') 'OGF:ICG_TO_EZONE(.):  OGF%NCGE', OGF%NCGE
WRITE(MSG%LU_DEBUG,'(10I6)') (OGF%ICG_TO_EZONE(IZ), IZ=1, OGF%NCGE)
WRITE(MSG%LU_DEBUG,*) '------------------------------------------------'
WRITE(MSG%LU_DEBUG,'(A,I6)') 'OGC%ICG_TO_ICW: OGC%NCGI', OGC%NCGI
WRITE(MSG%LU_DEBUG,'(10I6)')  OGC%ICG_TO_ICW(1:OGC%NCGI,1)
WRITE(MSG%LU_DEBUG,*) '------------------------------------------------'
WRITE(MSG%LU_DEBUG,'(A,I6)') 'OGC%ICG_TO_ICE: OGC%NCGE', OGC%NCGE
WRITE(MSG%LU_DEBUG,'(10I6)')  OGC%ICG_TO_ICE(1:OGC%NCGE,1)
WRITE(MSG%LU_DEBUG,*) '------------------------------------------------'
WRITE(MSG%LU_DEBUG,*) 'NCGE_TOTAL = ', NCGE_TOTAL
WRITE(MSG%LU_DEBUG,*) 'OGC%NCE  = ', OGC%NCE
WRITE(MSG%LU_DEBUG,*) 'OGC%NCE2 = ', OGC%NCE2
WRITE(MSG%LU_DEBUG,*) 'OGC%NLEN_BUFFER_LAYER1 = ', OGC%NLEN_BUFFER_LAYER1
WRITE(MSG%LU_DEBUG,*) 'OGC%NLEN_BUFFER_LAYER2 = ', OGC%NLEN_BUFFER_LAYER2
WRITE(MSG%LU_DEBUG,*) 'OGC%NLEN_BUFFER_LAYER4 = ', OGC%NLEN_BUFFER_LAYER4
WRITE(MSG%LU_DEBUG,*) 'OGC%NLEN_BUFFER_STENCIL= ', OGC%NLEN_BUFFER_STENCIL
WRITE(MSG%LU_DEBUG,*) 'OGC%NLEN_BUFFER_FULL   = ', OGC%NLEN_BUFFER_FULL
#endif

   ENDDO NEIGHBORS_LOOP
ENDDO MESHES_LOOP

END SUBROUTINE SCARC_EXTRACT_ZONE_POINTERS


! ------------------------------------------------------------------------------------------------------
! Determine tentative prolongator for current level by computing QR-decomposition of smoothed 
! nullspace vector and set nullspace for next level
! Compute the tentative prolongator, T, which is a tentative interpolation
! matrix from the coarse-grid to the fine-grid.  T exactly interpolates  B_fine = T B_coarse.
! ------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_PROLONGATION(NL)
USE SCARC_POINTERS, ONLY:  G, A, OA, P, OP, GF, PF
INTEGER, INTENT(IN) :: NL
REAL(EB):: DSUM, SCAL, PSAVE !, TOL = 1.0E-12_EB
INTEGER :: NM, NOM, IC, JC, ICC, ICOL, ICCOL, JCCOL, IP, JCC, IQ, INBR

!
! Allocate several workspaces (with conservative bounds which will be reduced later)
!    - prolongation matrix on internal part of mesh 
!    - for every neighbor small prolongation matrix for corresponding overlap
!    - vectors Q and R for QR-decomposition of aggregation zones operator
! Initialize QR-decomposition
!
DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G

   A => SCARC_POINT_TO_CMATRIX(G, NSCARC_MATRIX_POISSON)
   P => SCARC_POINT_TO_CMATRIX(G, NSCARC_MATRIX_PROLONGATION)
   Z => SCARC_POINT_TO_CMATRIX(G, NSCARC_MATRIX_ZONES)

   P%N_VAL = A%N_VAL + 1        
   P%N_ROW = G%NCE + 1                  
   CALL SCARC_ALLOCATE_CMATRIX(P, NL, NSCARC_PRECISION_DOUBLE, NSCARC_MATRIX_FULL, 'G%PROLONGATION')

   DO INBR = 1, SCARC(NM)%N_NEIGHBORS

      NOM = S%NEIGHBORS(INBR)
      CALL SCARC_POINT_TO_OTHER_GRID(NM, NOM, NL)

      OA => SCARC_POINT_TO_OTHER_CMATRIX(OG, NSCARC_MATRIX_POISSON)
      OP => SCARC_POINT_TO_OTHER_CMATRIX(OG, NSCARC_MATRIX_PROLONGATION)

      OP%N_VAL = OA%N_VAL + 1              ! TODO : CHECK : MUCH TOO BIG !!!
      OP%N_ROW = G%NCE + 1                 ! TODO : CHECK : MUCH TOO BIG !!!
      CALL SCARC_ALLOCATE_CMATRIX(OP, NL, NSCARC_PRECISION_DOUBLE, NSCARC_MATRIX_FULL, 'OG%PROLONGATION')

   ENDDO

   CALL SCARC_ALLOCATE_REAL1(G%QQ, 1, G%NCE2+1, NSCARC_INIT_ZERO, 'G%QQ')
   CALL SCARC_ALLOCATE_REAL1(G%RR, 1, G%NCE2+1, NSCARC_INIT_ZERO, 'G%RR')     ! TODO check length!


   !
   ! Copy blocks into Q according to aggregation zones and compute norms for single ZONES
   ! In each cell corresponding to a single zone, store square-sum of entries
   !
   IQ = 1
   G%AUX1 = 0.0_EB
   DO ICC = 1, Z%N_ROW-1
      DSUM = 0.0_EB
      DO ICOL = Z%ROW(ICC), Z%ROW(ICC+1)-1
         IC = Z%COL(ICOL)
         G%QQ(IQ) = G%NULLSPACE(IC)
         DSUM = DSUM + G%QQ(IQ)**2
#ifdef WITH_SCARC_DEBUG2
WRITE(MSG%LU_DEBUG,'(A,3I6,2E12.4)') 'QQ: ICC, IC, IQ, QQ(IQ), DSUM:', ICC, IC, IQ, G%QQ(IQ), DSUM
#endif
         IQ = IQ + 1
      ENDDO
      DO ICOL = Z%ROW(ICC), Z%ROW(ICC+1)-1
         G%AUX1(Z%COL(ICOL)) = DSUM
      ENDDO
      G%RR(ICC) = DSUM
   ENDDO

#ifdef WITH_SCARC_DEBUG
   WRITE(MSG%LU_DEBUG,*) '================================ PART 0 :'
   WRITE(MSG%LU_DEBUG,*) 'Z%N_ROW:', Z%N_ROW
   WRITE(MSG%LU_DEBUG,*) 'G%NULLSPACE:'
   WRITE(MSG%LU_DEBUG,'(6E12.4)') G%NULLSPACE(1:G%NCE2)
   WRITE(MSG%LU_DEBUG,*) 'G%AUX2:'
   WRITE(MSG%LU_DEBUG,'(6E12.4)') G%AUX2(1:G%NCE2)
   WRITE(MSG%LU_DEBUG,*) 'G%QQ:'
   WRITE(MSG%LU_DEBUG,'(6E12.4)') G%QQ(1:G%NCE2)
   WRITE(MSG%LU_DEBUG,*) 'G%AUX1:'
   WRITE(MSG%LU_DEBUG,'(6E12.4)') G%AUX1(1:G%NCE2)
   WRITE(MSG%LU_DEBUG,*) 'G%RR:'
   WRITE(MSG%LU_DEBUG,'(6E12.4)') G%RR
#endif

ENDDO

!
! Exchange sums of nullspace entries within single aggregation zones
!
IF (NMESHES > 1) CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_AUX, NSCARC_NONE, NL)
   
#ifdef WITH_SCARC_DEBUG
   WRITE(MSG%LU_DEBUG,*) '================================ PART 1 :'
   WRITE(MSG%LU_DEBUG,*) 'G%QQ:'
   WRITE(MSG%LU_DEBUG,'(6E12.4)') G%QQ(1:G%NCE2)
   WRITE(MSG%LU_DEBUG,*) 'G%AUX1:'
   WRITE(MSG%LU_DEBUG,'(6E12.4)') G%AUX1(1:G%NCE2)
#endif
!
! Build norms over single zones and scale Q-entries by norms
!
DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID (NM, NL)                          
   Z => SCARC_POINT_TO_CMATRIX(G, NSCARC_MATRIX_ZONES)

   DO ICC = 1, Z%N_ROW-1
      ICOL = Z%ROW(ICC)
      IC = Z%COL(ICOL)
      G%RR(ICC) = SQRT(G%AUX1(IC))
#ifdef WITH_SCARC_DEBUG2
WRITE(MSG%LU_DEBUG,*) 'G%RR(',ICC,')=', G%RR(ICC)
#endif
   ENDDO

   IQ = 1
   DO ICC = 1, Z%N_ROW-1
      DO ICOL = Z%ROW(ICC), Z%ROW(ICC+1)-1
         IC = Z%COL(ICOL)
         G%QQ(IQ) = G%QQ(IQ)/G%RR(ICC)
#ifdef WITH_SCARC_DEBUG2
WRITE(MSG%LU_DEBUG,*) 'G%QQ(',IQ,')=', G%QQ(IQ)
#endif
         IQ = IQ + 1
      ENDDO
   ENDDO

ENDDO
   
#ifdef WITH_SCARC_DEBUG
   WRITE(MSG%LU_DEBUG,*) '================================ PART 2 :'
   WRITE(MSG%LU_DEBUG,*) 'G%QQ:'
   WRITE(MSG%LU_DEBUG,'(6E12.4)') G%QQ(1:G%NCE2)
#endif

! ----------- Relax prolongator:
! Smooth the tentative prolongator, so that it's accuracy is greatly improved for algebraically smooth error.
! Compute:                P =  P - A_Dinv * P   
! with:                   A_Dinv = 4/3 * 1/rho * D^{-1} A   
!
! First step: Compute P_0: = A_Dinv * Q
!
DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID (NM, NL)                                ! Sets grid pointer G

   A => SCARC_POINT_TO_CMATRIX (G, NSCARC_MATRIX_POISSON)
   P => SCARC_POINT_TO_CMATRIX (G, NSCARC_MATRIX_PROLONGATION)
   Z => SCARC_POINT_TO_CMATRIX (G, NSCARC_MATRIX_ZONES)

   MG%OMEGA = 4.0_EB/3.0_EB                                         ! currently used default
   SCAL = MG%OMEGA/MG%APPROX_SPECTRAL_RADIUS                        ! for testing purposes rho is set to 2 currently
   
#ifdef WITH_SCARC_DEBUG
   WRITE(MSG%LU_DEBUG,*) '================================ PART 3 :'
   WRITE(MSG%LU_DEBUG,*) 'G%QQ:'
   WRITE(MSG%LU_DEBUG,'(6E12.4)') G%QQ
   CALL SCARC_DEBUG_CMATRIX(Z, 'ZONES','BEFORE FIRST RELAX STEP')
#endif
   
   IP = 1
   P%ROW(1) = IP
   DO IC = 1, G%NC
      DO ICC = 1, Z%N_ROW-1
   
         DSUM = 0.0_EB
         DO ICCOL = Z%ROW(ICC), Z%ROW(ICC+1)-1
            JC = Z%COL(ICCOL)
            ICOL = SCARC_MATCH_MATRIX_COLUMN(A, IC, JC)
            IF (ICOL /= -1) THEN
               DSUM = DSUM - SCAL * G%DIAG(IC) * A%VAL(ICOL) * G%QQ(ICCOL)
#ifdef WITH_SCARC_DEBUG2
WRITE(MSG%LU_DEBUG,'(A,5I6,4E12.4)') 'RELAX: IC, ICC, ICCOL, JC, ICOL, DSUM:', &
                              IC, ICC, ICCOL, JC, ICOL, DSUM, G%DIAG(IC), A%VAL(ICOL), G%QQ(ICCOL)
#endif
            ENDIF
         ENDDO
   
         IF (ABS(DSUM) /= 0.0_EB) THEN
            P%VAL(IP) = DSUM
            P%COL(IP) = ICC
#ifdef WITH_SCARC_DEBUG2
WRITE(MSG%LU_DEBUG,*) '-----------------> IP, VAL, COL,:', IP, DSUM, ICC
#endif
            IP = IP + 1
         ENDIF
   
      ENDDO
      P%ROW(IC+1) = IP
   ENDDO
   P%N_VAL = IP - 1

#ifdef WITH_SCARC_DEBUG
   CALL SCARC_DEBUG_CMATRIX(Z, 'ZONES','AFTER FIRST RELAX STEP')
   CALL SCARC_DEBUG_CMATRIX(P, 'PROLONGATION','AFTER FIRST RELAX STEP')
#endif

ENDDO
   
!
! Second step: Compute P: = P - P_0
!
DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G
   P => SCARC_POINT_TO_CMATRIX(G, NSCARC_MATRIX_PROLONGATION)
   Z => SCARC_POINT_TO_CMATRIX(G, NSCARC_MATRIX_ZONES)

#ifdef WITH_SCARC_DEBUG
CALL SCARC_DEBUG_CMATRIX(Z, 'ZONES','AFTER RESORT PROL ')
CALL SCARC_DEBUG_CMATRIX(P, 'PROLONGATION','AFTER RESORT PROL ')
#endif

   DO ICC = 1, Z%N_ROW-1
      DO ICCOL = Z%ROW(ICC), Z%ROW(ICC+1) - 1
         IC = Z%COL(ICCOL)

         IF (IC > G%NC) CYCLE
         DO JCCOL = P%ROW(IC), P%ROW(IC+1) - 1
            JCC = P%COL(JCCOL)
#ifdef WITH_SCARC_DEBUG2
   WRITE(MSG%LU_DEBUG,'(A, 5I5)') 'Subtracting: ICC, ICCOL, IC, JCCOL, JCC:', ICC, ICCOL, IC, JCCOL, JCC
#endif
            IF (JCC == ICC) THEN
               PSAVE = P%VAL(JCCOL)
               P%VAL(JCCOL) = P%VAL(JCCOL) + G%QQ(ICCOL)
#ifdef WITH_SCARC_DEBUG2
   WRITE(MSG%LU_DEBUG,'(A, 3E12.4)') '       PSAVE,Q,PF  ---->', PSAVE, G%QQ(ICCOL), P%VAL(JCCOL)
#endif
               CYCLE
            ENDIF
         ENDDO
      ENDDO
   ENDDO
   
#ifdef WITH_SCARC_DEBUG
   CALL SCARC_DEBUG_CMATRIX(P, 'G%PROLONGATION','SETUP_PROLONGATION: AFTER RELAX STEP, BEFORE EXCHANGE ')
#endif
ENDDO

!
! Determine global columns array for prolongation matrix
!
DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_MULTIGRID (NM, NL, NL+1)                
   PF => SCARC_POINT_TO_CMATRIX(GF, NSCARC_MATRIX_PROLONGATION)

   DO IC = 1, GF%NC
#ifdef WITH_SCARC_DEBUG2
WRITE(MSG%LU_DEBUG,*) 'PROL: ------------ IC =', IC
#endif
      DO JCCOL = PF%ROW(IC), PF%ROW(IC+1) - 1
         JCC = PF%COL(JCCOL)
         PF%COLG(JCCOL) = GC%LOCAL_TO_GLOBAL(JCC)
#ifdef WITH_SCARC_DEBUG2
WRITE(MSG%LU_DEBUG,*) 'PROL: LOCAL_TO_GLOBAL: IC, JCCOL, JCC, RES:', IC, JCCOL, JCC, PF%COLG(JCCOL)
#endif
      ENDDO
   ENDDO
   
#ifdef WITH_SCARC_DEBUG
   CALL SCARC_DEBUG_CMATRIX(PF, 'G%PROLONGATION','SETUP_PROLONGATION: AFTER LAST EXCHANGE')
#endif
ENDDO

!
! Exchange resulting columns and values of prolongation matrix and extract exchanged data from 
! overlapping parts with single neighbors and attach them to main matrix
!
IF (NMESHES > 1) THEN
   CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_POISSON_COLS,  NSCARC_MATRIX_PROLONGATION, NL)
   CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_POISSON_COLSG, NSCARC_MATRIX_PROLONGATION, NL)
   CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_POISSON_VALS,  NSCARC_MATRIX_PROLONGATION, NL)
   CALL SCARC_EXTRACT_MATRIX_OVERLAPS(NSCARC_MATRIX_PROLONGATION, NL)
ENDIF
   
DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G
   P => SCARC_POINT_TO_CMATRIX(G, NSCARC_MATRIX_PROLONGATION)
   CALL SCARC_REDUCE_CMATRIX(P, 'P%PROLONGATION')

   DO INBR = 1, S%N_NEIGHBORS
      NOM = S%NEIGHBORS(INBR)
      CALL SCARC_POINT_TO_OTHER_GRID(NM, NOM, NL)
      OP => SCARC_POINT_TO_OTHER_CMATRIX(OG, NSCARC_MATRIX_PROLONGATION)
      CALL SCARC_REDUCE_CMATRIX(OP, 'OP%PROLONGATION')
   ENDDO

#ifdef WITH_SCARC_DEBUG
   CALL SCARC_DEBUG_CMATRIX(P, 'G%PROLONGATION','SETUP_PROLONGATION: FINAL')
#endif

ENDDO

END SUBROUTINE SCARC_SETUP_PROLONGATION


! ------------------------------------------------------------------------------------------------------
! If coarsest level isn't reached yet, define nullspace for next coarser level
! ------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_NULLSPACE_COARSE(NL)
USE SCARC_POINTERS, ONLY:  GC, GF, ZF
INTEGER, INTENT(IN) :: NL
INTEGER :: NM

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_MULTIGRID(NM, NL, NL+1)                   ! Sets pointers GC and GF
   ZF => SCARC_POINT_TO_CMATRIX (GF, NSCARC_MATRIX_ZONES)

   IF (NL < NLEVEL_MAX) THEN
      !GC%N_FINE = GF%N_COARSE
      GC%N_FINE = GC%NC_LOCAL(NM)
      !GC%NC = GC%N_FINE
      !GC%NCE = ZF%N_ROW-1
      CALL SCARC_ALLOCATE_REAL1(GC%NULLSPACE, 1, GC%NCE, NSCARC_INIT_ZERO, 'GC%NULLSPACE')
      GC%NULLSPACE(1:GC%NCE) = GF%RR(1:GC%NCE)
      CALL SCARC_REDUCE_INT1(GC%LOCAL_TO_GLOBAL, GC%NCE, 'GC%LOCAL_TO_GLOBAL')
   ENDIF

#ifdef WITH_SCARC_DEBUG
   WRITE(MSG%LU_DEBUG,*) '============== SETUP_NULLSPACE_COARSE ================='
   WRITE(MSG%LU_DEBUG,*) 'GF%NULLSPACE:'
   WRITE(MSG%LU_DEBUG,'(6E12.4)') GF%NULLSPACE
   WRITE(MSG%LU_DEBUG,*) 'GF%RR'
   WRITE(MSG%LU_DEBUG,'(6E12.4)') GF%RR
   WRITE(MSG%LU_DEBUG,*) 'GF%QR:'
   WRITE(MSG%LU_DEBUG,'(6E12.4)') GF%QQ
   WRITE(MSG%LU_DEBUG,*) 'GC%LOCAL_TO_GLOBAL:'
   WRITE(MSG%LU_DEBUG,'(8I6)') GC%LOCAL_TO_GLOBAL
   WRITE(MSG%LU_DEBUG,*) 'GC%N_FINE:', GC%N_FINE
   WRITE(MSG%LU_DEBUG,*) 'GC%NC:', GC%NC
   WRITE(MSG%LU_DEBUG,*) 'GC%NCE:', GC%NCE
   WRITE(MSG%LU_DEBUG,*) '==============================================='
#endif

ENDDO
   
END SUBROUTINE SCARC_SETUP_NULLSPACE_COARSE


! ------------------------------------------------------------------------------------------------------
! Determine which columns of system matrix are involved in multiplication with tentative prolongator
! ------------------------------------------------------------------------------------------------------
INTEGER FUNCTION SCARC_MATCH_MATRIX_COLUMN(A, IC, JC)
TYPE (SCARC_MATRIX_COMPACT_TYPE), POINTER :: A
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
! Setup restriction matrix: Build transpose of prolongation matrix
! Compute the restriction matrix, R, which interpolates from the fine-grid to the coarse-grid.
! ------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_RESTRICTION(NL)
USE SCARC_POINTERS, ONLY: GC, GF, RF, PF
INTEGER, INTENT(IN) :: NL
INTEGER :: NM, IC, IRROW, IRCOL, IPCOL, IPC, ICCL, ICCG, IFOUND 
LOGICAL :: IS_INCLUDED

! Allocate restriction matrix R on internal mesh part
MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_MULTIGRID (NM, NL, NL+1)            ! Sets pointers GF and GC to fine and coarse level

   PF => SCARC_POINT_TO_CMATRIX (GF, NSCARC_MATRIX_PROLONGATION)
   RF => SCARC_POINT_TO_CMATRIX (GF, NSCARC_MATRIX_RESTRICTION)

   ! TODO: Resize R
   RF%N_VAL = PF%N_VAL + 100
   RF%N_ROW = PF%N_ROW + 100

#ifdef WITH_SCARC_DEBUG
   WRITE(MSG%LU_DEBUG,*) 'RESTRICTION:  N_VAL, N_ROW=', RF%N_VAL, RF%N_ROW 
   WRITE(MSG%LU_DEBUG,*) 'PROLONGATION: N_VAL, N_ROW=', PF%N_VAL, PF%N_ROW 
   WRITE(MSG%LU_DEBUG,*) 'GC%LOCAL_TO_GLOBAL:', GC%LOCAL_TO_GLOBAL
   WRITE(MSG%LU_DEBUG,*) 'GF%NC:', GF%NC
   WRITE(MSG%LU_DEBUG,*) 'GF%NCE:', GF%NCE
#endif
   CALL SCARC_ALLOCATE_CMATRIX(RF, NL, NSCARC_PRECISION_DOUBLE, NSCARC_MATRIX_FULL, 'GF%RESTRICTION')

   IRROW = 1                             ! counter of current row of restriction matrix - corresponds to coarse cells
   IRCOL = 1                             ! counter of current column of restriction matrix
   RF%ROW(IRROW) = IRCOL

   LOCAL_COARSE_CELLS_LOOP: DO ICCL = 1, GC%NCE

      ICCG = GC%LOCAL_TO_GLOBAL(ICCL)                     ! corresponding global coarse cell

      IFOUND = -1
      IFOUND = FINDLOC(PF%COLG, VALUE = ICCG, DIM=1)
      IF (IFOUND == -1) CYCLE

      IS_INCLUDED = .FALSE.

      FINE_CELLS_LOOP: DO IC = 1, GF%NCE                  ! counter of fine cell (including overlaps)

#ifdef WITH_SCARC_DEBUG2
   WRITE(MSG%LU_DEBUG,*) ' ---------->  ICCL, ICCG, IC = ', ICCL, ICCG, IC
#endif
         ROW_LOOP: DO IPCOL = PF%ROW(IC), PF%ROW(IC+1)-1
            IPC = PF%COLG(IPCOL)
            IF (IPC == ICCG) THEN
               IS_INCLUDED = .TRUE.
               RF%VAL(IRCOL) = PF%VAL(IPCOL)
               RF%COLG(IRCOL) = IC

#ifdef WITH_SCARC_DEBUG2
   WRITE(MSG%LU_DEBUG,'(A,4I6,E12.4)') '     IRCOL, IPCOL, IRROW, RF%COLG, RF%VAL:', &
                                             IRCOL, IPCOL, IRROW, RF%COLG(IRCOL), RF%VAL(IRCOL)
#endif
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

#ifdef WITH_SCARC_DEBUG
   CALL SCARC_DEBUG_CMATRIX (RF, 'GF%RESTRICTION','AFTER SETUP_RESTRICTION')
#endif
   CALL SCARC_REDUCE_CMATRIX (RF, 'GF%RESTRICTION')

ENDDO MESHES_LOOP

END SUBROUTINE SCARC_SETUP_RESTRICTION


! ------------------------------------------------------------------------------------------------------
! Find matching column index during matrix-matrix multiplication of compact matrices
! ------------------------------------------------------------------------------------------------------
INTEGER FUNCTION SCARC_FIND_MATCHING_COLUMN(A, JC, ICCG)
TYPE (SCARC_MATRIX_COMPACT_TYPE), POINTER :: A
INTEGER, INTENT(IN) :: JC, ICCG
INTEGER :: IPCOL

#ifdef WITH_SCARC_DEBUG2
WRITE(MSG%LU_DEBUG,*) 'FIND_MATCHING_COLUMN: ICCG, JC:', ICCG, JC
#endif
SCARC_FIND_MATCHING_COLUMN = -1
DO IPCOL = A%ROW(JC), A%ROW(JC+1)-1
   IF (A%COLG(IPCOL) == ICCG) THEN
      SCARC_FIND_MATCHING_COLUMN = IPCOL
#ifdef WITH_SCARC_DEBUG2
WRITE(MSG%LU_DEBUG,*) ' -----------> JC, IPCOL, ICCG, VAL:', JC, IPCOL, ICCG, A%VAL(IPCOL)
#endif
   ENDIF
ENDDO
RETURN 

END FUNCTION SCARC_FIND_MATCHING_COLUMN


! ------------------------------------------------------------------------------------------------------
! Determine coarse grid matrix ACC which is computed as Galerkin operator R*A*P
! ------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_AP(NL)
USE SCARC_POINTERS, ONLY: GC, GF, AF, PF, APF, OA, OAP
INTEGER, INTENT(IN) :: NL
INTEGER  :: NM, NOM, IC, JC, IACOL, IPCOL, IP, ICCG, INBR, IFOUND
REAL(EB) :: DSUM, TOL = 1E-12_EB

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_MULTIGRID (NM, NL, NL+1)

   AF  => SCARC_POINT_TO_CMATRIX (GF, NSCARC_MATRIX_POISSON)          
   PF  => SCARC_POINT_TO_CMATRIX (GF, NSCARC_MATRIX_PROLONGATION)     
   APF => SCARC_POINT_TO_CMATRIX (GF, NSCARC_MATRIX_POISSON_PROL)         

#ifdef WITH_SCARC_DEBUG
   CALL SCARC_DEBUG_CMATRIX (AF, 'POISSON-FINE', 'START OF SETUP_AP')
   CALL SCARC_DEBUG_CMATRIX (PF, 'PROLONGATION-FINE', 'START OF SETUP_AP')
#endif

   ! TODO  Reduce size of AP
   APF%N_ROW = AF%N_ROW
   APF%N_VAL = APF%N_ROW*20            ! TODO: only temporarily
   CALL SCARC_ALLOCATE_CMATRIX(APF, NL, NSCARC_PRECISION_DOUBLE, NSCARC_MATRIX_FULL, 'GF%AP')

   DO INBR = 1, SCARC(NM)%N_NEIGHBORS

      NOM = SCARC(NM)%NEIGHBORS(INBR)
      CALL SCARC_POINT_TO_OTHER_GRID(NM, NOM, NL)

      OA  => SCARC_POINT_TO_OTHER_CMATRIX(OGF, NSCARC_MATRIX_POISSON)
      OAP => SCARC_POINT_TO_OTHER_CMATRIX(OGF, NSCARC_MATRIX_POISSON_PROL)

      OAP%N_VAL = AF%N_VAL              ! TODO : CHECK : MUCH TOO BIG !!!
      OAP%N_ROW = GF%NCE + 1            ! TODO : CHECK : MUCH TOO BIG !!!
      CALL SCARC_ALLOCATE_CMATRIX(OAP, NL, NSCARC_PRECISION_DOUBLE, NSCARC_MATRIX_FULL, 'OGF%AP')

   ENDDO

   IP = 1
   APF%ROW(1) = IP

   DO IC = 1, GF%NC

#ifdef WITH_SCARC_DEBUG2
WRITE(MSG%LU_DEBUG,*) '=================== IC = ', IC ,' ======================='
#endif
      DO ICCG = 1, GC%NC_GLOBAL
         IFOUND = -1
         IFOUND = FINDLOC(PF%COLG, VALUE=ICCG, DIM=1)
#ifdef WITH_SCARC_DEBUG2
WRITE(MSG%LU_DEBUG,*) '------------------- ICCG = ', ICCG,' IFOUND =',IFOUND
#endif
         IF (IFOUND == -1) CYCLE
      

         DSUM = 0.0_EB
         DO IACOL = AF%ROW(IC), AF%ROW(IC+1)-1
            JC = AF%COL(IACOL)
            IF (JC < 0) CYCLE
            IPCOL = SCARC_FIND_MATCHING_COLUMN(PF, JC, ICCG) 
#ifdef WITH_SCARC_DEBUG2
WRITE(MSG%LU_DEBUG,*) '                JC, IPCOL =', JC, IPCOL
#endif
            IF (IPCOL /= -1) THEN
               DSUM = DSUM + AF%VAL(IACOL) * PF%VAL(IPCOL)
#ifdef WITH_SCARC_DEBUG2
WRITE(MSG%LU_DEBUG,'(A,2I6,3E12.4)') 'SETUP_AP:, JC, PF%COL(IPCOL), AF%VAL(IACOL), PF%VAL(IPCOL), DSUM:', &
                                     JC, PF%COL(IPCOL), AF%VAL(IACOL), PF%VAL(IPCOL), DSUM
#endif
            ENDIF
         ENDDO

         IF (ABS(DSUM) > TOL) THEN
            APF%COL(IP)  = ICCG
            APF%COLG(IP) = ICCG
            APF%VAL(IP)  = DSUM
#ifdef WITH_SCARC_DEBUG2
WRITE(MSG%LU_DEBUG,'(A,2I6,E12.4)') '      bingo: IP, COL, VAL :', IP, ICCG, DSUM
#endif
            IP = IP + 1
         ENDIF

      ENDDO
      APF%ROW(IC+1) = IP

   ENDDO

ENDDO

! Exchange overlapping parts of prolongation matrix and extract exchanged data
IF (NMESHES > 1) THEN
   CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_POISSON_COLS,  NSCARC_MATRIX_POISSON_PROL, NL)
   CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_POISSON_COLSG, NSCARC_MATRIX_POISSON_PROL, NL)
   CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_POISSON_VALS,  NSCARC_MATRIX_POISSON_PROL, NL)
   CALL SCARC_EXTRACT_MATRIX_OVERLAPS(NSCARC_MATRIX_POISSON_PROL, NL)
ENDIF

! Reduce workspace for prolongation matrix to really needed size
DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_MULTIGRID (NM, NL, NL+1)
   APF => SCARC_POINT_TO_CMATRIX (GF, NSCARC_MATRIX_POISSON_PROL)         

#ifdef WITH_SCARC_DEBUG
   CALL SCARC_DEBUG_CMATRIX (APF, 'AP-FINE','END SETUP_AP')
#endif

   CALL SCARC_REDUCE_CMATRIX (APF, 'POISSON-PROL')

ENDDO

END SUBROUTINE SCARC_SETUP_AP

! ------------------------------------------------------------------------------------------------------
! ------------------------------------------------------------------------------------------------------
INTEGER FUNCTION SCARC_MAP_GLOBAL_TO_LOCAL(G, ICG, NM, NL)
USE SCARC_POINTERS, ONLY: SUB
TYPE (SCARC_GRID_TYPE), POINTER, INTENT(IN) :: G
INTEGER, INTENT(IN) :: ICG, NM, NL
INTEGER :: NOM, IOR0

SUB => SUBDIVISION
SCARC_MAP_GLOBAL_TO_LOCAL = -1

MESHES_LOOP: DO NOM = 1, NMESHES
   IF (G%NC_OFFSET(NM) < ICG .AND. ICG <= G%NC_OFFSET(NM) + G%NC_LOCAL(NM)) THEN
      EXIT MESHES_LOOP
   ENDIF
ENDDO MESHES_LOOP

IF (ARE_NEIGHBORS(NM, NOM)) THEN
   DO IOR0 = -3, 3
      CALL SCARC_POINT_TO_OTHER_GRID(NM, NOM, NL)
      !DO ICG = OL%GHOST_FIRSTW(IOR0), OL%GHOST_LASTW(IOR0)
      !   ICE = OG%ICG_TO_ICE(ICG, 1)
      !   
      !ENDDO
   ENDDO
ENDIF


END FUNCTION SCARC_MAP_GLOBAL_TO_LOCAL

! ------------------------------------------------------------------------------------------------------
! Note: Matrix POISPROL corresponds to POISSON x PROLONGATION  ~ AP
! ------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_GALERKIN(NL)
USE SCARC_POINTERS, ONLY: GF, GC, AP, RF, AC, OAC
INTEGER, INTENT(IN) :: NL
INTEGER  :: NM, ICCL, ICCG, JCCL, JCCG, IRCOL, IAPCOL, JC, IP, INBR, NOM, IFOUND
REAL(EB) :: DSUM, TOL = 1E-12_EB

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_MULTIGRID (NM, NL, NL+1)

   RF => SCARC_POINT_TO_CMATRIX (GF, NSCARC_MATRIX_RESTRICTION)         
   AP => SCARC_POINT_TO_CMATRIX (GF, NSCARC_MATRIX_POISSON_PROL)         
   AC => SCARC_POINT_TO_CMATRIX (GC, NSCARC_MATRIX_POISSON)         

#ifdef WITH_SCARC_DEBUG
   CALL SCARC_DEBUG_CMATRIX (RF, 'RESTRICTION-FINE', 'START OF SETUP_GALERKIN')
#endif

   IF (.NOT.ALLOCATED (AC%VAL)) THEN
      AC%N_ROW = GC%NCE+1
      AC%N_VAL = AC%N_ROW**2             ! only temporarily TODO TOO BIG
      CALL SCARC_ALLOCATE_CMATRIX(AC, NL+1, NSCARC_PRECISION_DOUBLE, NSCARC_MATRIX_FULL, 'GC%POISSONC')
   ENDIF

   DO INBR = 1, SCARC(NM)%N_NEIGHBORS

      NOM = SCARC(NM)%NEIGHBORS(INBR)
      CALL SCARC_POINT_TO_OTHER_MULTIGRID(NM, NOM, NL, NL+1)

      OAC => SCARC_POINT_TO_OTHER_CMATRIX(OGC, NSCARC_MATRIX_POISSON)

      OAC%N_VAL = AC%N_VAL              ! TODO : CHECK : MUCH TOO BIG !!!
      OAC%N_ROW = GC%NCE2 + 1           ! TODO : CHECK : MUCH TOO BIG !!!
      CALL SCARC_ALLOCATE_CMATRIX(OAC, NL, NSCARC_PRECISION_DOUBLE, NSCARC_MATRIX_FULL, 'OAC%POISSONC')

   ENDDO

#ifdef WITH_SCARC_DEBUG
CALL SCARC_DEBUG_CMATRIX (SCARC(MYID+1)%LEVEL(1)%STRUCTURED%POISSONC, 'AF', 'A in SETUP_GALERKIN 1')
#endif

   IP = 1
   AC%ROW(1) = IP

   LOCAL_COARSE_CELLS_LOOP: DO ICCL = 1, GC%NC

      ICCG = GC%LOCAL_TO_GLOBAL(ICCL)               ! corresponding global coarse cell

       GLOBAL_COARSE_CELLS_LOOP: DO JCCG = 1, GC%NC_GLOBAL                    

          IF (NMESHES == 1) THEN
             JCCL = JCCG
          ELSE
             IFOUND = FINDLOC (GC%LOCAL_TO_GLOBAL, VALUE = JCCG, DIM = 1)
             IF (IFOUND /= 0) THEN
                JCCL = IFOUND
             ELSE
                JCCL = 0
             ENDIF
          ENDIF

#ifdef WITH_SCARC_DEBUG2
WRITE(MSG%LU_DEBUG,*) '--------------- ICCL, ICCG, JCCL, JCCG = ', ICCL, ICCG, JCCL, JCCG 
#endif
         DSUM = 0.0_EB
         DO IRCOL = RF%ROW(ICCL), RF%ROW(ICCL+1)-1
            JC = RF%COLG(IRCOL)
            IAPCOL = SCARC_FIND_MATCHING_COLUMN(AP, JC, JCCG) 
            IF (IAPCOL /= -1) THEN
               DSUM = DSUM + RF%VAL(IRCOL) * AP%VAL(IAPCOL)
#ifdef WITH_SCARC_DEBUG2
WRITE(MSG%LU_DEBUG,'(A,I6,2E12.4)') 'SETUP_GALERKIN:, JC, AP%VAL(IAPCOL), DSUM:', JC, AP%VAL(IAPCOL), DSUM
#endif
            ENDIF
         ENDDO

         IF (ABS(DSUM) > TOL) THEN
            AC%COL(IP)  = JCCL
            AC%COLG(IP) = JCCG
            AC%VAL(IP) = DSUM
#ifdef WITH_SCARC_DEBUG2
WRITE(MSG%LU_DEBUG,'(A,3I6,E12.4)') ' -------- bingo: IP, COL, VAL :', IP, JCCL, JCCG, DSUM
#endif
            IP = IP + 1
         ENDIF
      ENDDO GLOBAL_COARSE_CELLS_LOOP

      AC%ROW(ICCL+1) = IP

   ENDDO LOCAL_COARSE_CELLS_LOOP

#ifdef WITH_SCARC_DEBUG
   CALL SCARC_DEBUG_CMATRIX (AC, 'POISSON-COARSE','END OF RAP')
#endif

   CALL SCARC_REDUCE_CMATRIX (AC, 'POISSON-COARSE')

ENDDO

CALL SCARC_RESORT_MATRIX_ROWS(NL+1)             
 
MESH_INT = 0                            
RANK_INT = 0

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   CALL SCARC_POINT_TO_MULTIGRID (NM, NL, NL+1)
   AC => SCARC_POINT_TO_CMATRIX (GC, NSCARC_MATRIX_POISSON)         
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'SETUP_GALERKIN: AC%N_STENCIL_MAX =', AC%N_STENCIL_MAX
#endif
   CALL SCARC_GET_MATRIX_STENCIL_MAX(AC, GC%NC)
   MESH_INT(NM) = AC%N_STENCIL_MAX
   RANK_INT = MAX(RANK_INT, MESH_INT(NM))
ENDDO

IF (N_MPI_PROCESSES>1) &
   CALL MPI_ALLREDUCE(MPI_IN_PLACE, RANK_INT, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, IERROR)

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'N_STENCIL_MAX  OVERALL:', RANK_INT
WRITE(MSG%LU_DEBUG,*) 'TYPE_MKL(0:2):', TYPE_MKL(0:2)
#endif

#ifdef WITH_MKL
IF (TYPE_MKL(NL+1) == NSCARC_MKL_LOCAL .OR. TYPE_MKL(NL+1) == NSCARC_MKL_GLOBAL) THEN
   DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
      CALL SCARC_SETUP_POISSON_MKL(NM, NL+1)
   ENDDO
ENDIF
#endif

END SUBROUTINE SCARC_SETUP_GALERKIN


! ------------------------------------------------------------------------------------------------------
! Compute entry of A*P matrix at position (IC, ICC)
! This consists of a summation over the entries:    P(:,ICC)*A(IC,:) 
! Thus, it must be checked, if - for a given entry of A in row IC - the prolongation matrix
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
#ifdef WITH_SCARC_DEBUG2
WRITE(MSG%LU_DEBUG,'(A,5I4,3E12.4)') 'VALUE_AP: IC, JC, ICC, GLOB(ICC), JCC, AF%VAL(IA), PF%VAL(IP), DSUM:',&
                                  IC, JC, ICC, GC%LOCAL_TO_GLOBAL(ICC), JCC, AF%VAL(IA), PF%VAL(IP), DSUM
#endif
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
! Determine coarse grid matrix ACC which is computed as Galerkin operator R*A*P
! ------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_GALERKIN_OPERATOR(NL)
USE SCARC_POINTERS, ONLY: GF, GC, AF, AC, RF, PF
INTEGER, INTENT(IN) :: NL
INTEGER  :: NM, IPC, IRC, IC, IR, IA
REAL(EB) :: DRSUM, DAP, TOL = 1.0E-12_EB

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_MULTIGRID (NM, NL, NL+1)

   PF => SCARC_POINT_TO_CMATRIX (GF, NSCARC_MATRIX_PROLONGATION)     
   RF => SCARC_POINT_TO_CMATRIX (GF, NSCARC_MATRIX_RESTRICTION)
   
   AF => SCARC_POINT_TO_CMATRIX (GF, NSCARC_MATRIX_POISSON)          
   AC => SCARC_POINT_TO_CMATRIX (GC, NSCARC_MATRIX_POISSON)          

#ifdef WITH_SCARC_DEBUG
CALL SCARC_DEBUG_CMATRIX (AF, 'AF','START OF SETUP_GALERKIN_OPERATOR')
#endif

   IF (.NOT.ALLOCATED (AC%VAL)) THEN
      AC%N_ROW = GC%NCE+1
      AC%N_VAL = AC%N_ROW**2             ! only temporarily TODO TOO BIG
      CALL SCARC_ALLOCATE_CMATRIX(AC, NL+1, NSCARC_PRECISION_DOUBLE, NSCARC_MATRIX_FULL, 'AC%POISSONC')
   ENDIF

   IA = 1
   AC%ROW(1) = 1

   DO IRC = 1, GC%NCE                      ! rows of restriction matrix

#ifdef WITH_SCARC_DEBUG2
WRITE(MSG%LU_DEBUG,*) '======================== IRC = ', IRC,' ============================'
#endif
      ! First process diagonal entry
      IPC = IRC

#ifdef WITH_SCARC_DEBUG2
WRITE(MSG%LU_DEBUG,*) '------------------------ IPC:A: = ', IPC,' ----------------------------'
WRITE(MSG%LU_DEBUG,*) '    R%ROW(IRW)-R%ROW(IRC+1)-1: ', (RF%COLG(IR), IR=RF%ROW(IRC),RF%ROW(IRC+1)-1)
WRITE(MSG%LU_DEBUG,*) '--------------------------------------------------------------------'
#endif
      ! Compute entry AC(IRC,IPC) = SUM( R(IRC,IC)*A(IC,:)*P(:,IPC))
      DRSUM = 0.0_EB
      DO IR = RF%ROW(IRC), RF%ROW(IRC+1) - 1
         IC = RF%COLG(IR)
         DAP = SCARC_VALUE_RAP(IC, IPC)
         DRSUM = DRSUM + DAP*RF%VAL(IR)
#ifdef WITH_SCARC_DEBUG2
         WRITE(MSG%LU_DEBUG,'(A,3I4,E12.4)') 'GAL: IRC, IC, IPC, DRSUM:', IRC, IC, IPC, DRSUM
#endif
      ENDDO
      IF (ABS(DRSUM) > TOL) THEN
         AC%VAL(IA) = DRSUM
         AC%COL(IA) = IPC
#ifdef WITH_SCARC_DEBUG2
         WRITE(MSG%LU_DEBUG,'(A,5I5,E12.4)') 'GALERKIN: IRC, IC, IPC, IA, COL, VAL:',&
                                              IRC,IC,IPC, IA,AC%COL(IA),AC%VAL(IA)
#endif
         IA = IA + 1
      ENDIF

      ! Then process remaining entries in row
      DO IPC = 1, GC%NCE                   ! columns of prolongation matrix

         IF (IPC == IRC) CYCLE

#ifdef WITH_SCARC_DEBUG2
WRITE(MSG%LU_DEBUG,*) '------------------------ IPC:B: = ', IPC,' ----------------------------'
#endif
         ! Compute entry AC(ICC,JCC) = SUM( P(:,JCC)*A(IC,:)*R(ICC,IC) )
         DRSUM = 0.0_EB
         DO IR = RF%ROW(IRC), RF%ROW(IRC+1) - 1
            IC = RF%COLG(IR)
            DAP = SCARC_VALUE_RAP(IC, IPC)
            DRSUM = DRSUM + DAP*RF%VAL(IR)
#ifdef WITH_SCARC_DEBUG2
         WRITE(MSG%LU_DEBUG,'(A,3I4,E12.4)') 'GAL: IRC, IC, IPC, DRSUM:', IRC, IC, IPC, DRSUM
#endif
         ENDDO
         IF (ABS(DRSUM) > TOL) THEN
            AC%VAL(IA) = DRSUM
            AC%COL(IA) = IPC
#ifdef WITH_SCARC_DEBUG2
         WRITE(MSG%LU_DEBUG,'(A,5I5,E12.4)') 'GALERKIN: IRC, IC, IPC, IA, COL, VAL:',&
                                              IRC,IC,IPC, IA,AC%COL(IA),AC%VAL(IA)
#endif
            IA = IA + 1
         ENDIF

      ENDDO

      AC%ROW(IRC+1) = IA
#ifdef WITH_SCARC_DEBUG2
WRITE(MSG%LU_DEBUG,*) 'AC%ROW(',IRC+1,')=', AC%ROW(IRC+1)
#endif
   ENDDO

   ! TODO temporarily
   AC%N_STENCIL = AF%N_STENCIL

   GC%NC = AC%N_ROW - 1

#ifdef WITH_SCARC_DEBUG
CALL SCARC_DEBUG_CMATRIX(AC, 'POISSON', 'COARSE POISSON AT END OF SETUP_GALERKIN')
#endif

#ifdef WITH_MKL
   CALL SCARC_SETUP_POISSON_MKL(NM, NL+1)
#endif

ENDDO

END SUBROUTINE SCARC_SETUP_GALERKIN_OPERATOR


! ------------------------------------------------------------------------------------------------------
! Resort matrix entries such that diagonal entry comes first
! ------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_RESORT_MATRIX_ROWS(NL)
USE SCARC_POINTERS, ONLY: G, A
INTEGER, INTENT(IN) :: NL
INTEGER:: NM, NCOL, ICOL, JCOL, KCOL, IC
INTEGER:: AUX_COL(100) = 0, AUX_COLG(100)
REAL(EB) :: AUX_VAL(100) = 0.0_EB
LOGICAL :: COLG_IS_DEFINED = .FALSE.

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   CALL SCARC_POINT_TO_GRID (NM, NL)
   A => SCARC_POINT_TO_CMATRIX (G, NSCARC_MATRIX_POISSON)

   IF (ALLOCATED(A%COLG)) COLG_IS_DEFINED = .TRUE.

   !DO IC = 1, G%NCE2
   DO IC = 1, G%NC
      AUX_COL = 0
      JCOL = 1
      DO ICOL = A%ROW(IC), A%ROW(IC+1)-1
         AUX_COL(JCOL) = A%COL(ICOL)
         IF (COLG_IS_DEFINED)  AUX_COLG(JCOL) = A%COLG(ICOL)
         AUX_VAL(JCOL) = A%VAL(ICOL)
         JCOL = JCOL + 1
      ENDDO
      NCOL = JCOL - 1

#ifdef WITH_SCARC_DEBUG2
WRITE(MSG%LU_DEBUG,*) '================= IC =',IC
WRITE(MSG%LU_DEBUG,*) '     AUX_COL =', AUX_COL(1:NCOL)
WRITE(MSG%LU_DEBUG,*) '     AUX_COLG=', AUX_COLG(1:NCOL)
WRITE(MSG%LU_DEBUG,*) '     AUX_VAL =', AUX_VAL(1:NCOL)
WRITE(MSG%LU_DEBUG,*) '     NCOL    = ', NCOL
#endif
      ! Find column index of diagonal element
      JCOL = 0
      DO WHILE (JCOL <= NCOL)
        JCOL = JCOL + 1
        IF (AUX_COL(JCOL) == IC) EXIT
      ENDDO
#ifdef WITH_SCARC_DEBUG2
WRITE(MSG%LU_DEBUG,*) '    -----> diag entry:', JCOL
#endif

      ! Store corresponding index and value in first matrix element of that row
      ICOL = A%ROW(IC)
      A%COL(ICOL)  = AUX_COL(JCOL)
      A%COLG(ICOL) = AUX_COLG(JCOL)
      A%VAL(ICOL)  = AUX_VAL(JCOL)

      AUX_COL(JCOL) = 99999999
      AUX_COLG(JCOL) = 99999999

      IF (COLG_IS_DEFINED) THEN
         JCOL = MINLOC(AUX_COLG(1:NCOL), DIM=1)
      ELSE
         JCOL = MINLOC(AUX_COL(1:NCOL), DIM=1)
      ENDIF
      KCOL = 1
#ifdef WITH_SCARC_DEBUG2
WRITE(MSG%LU_DEBUG,*) '    ------------->A'
WRITE(MSG%LU_DEBUG,*) '    A%COL(',ICOL,')=',A%COL(ICOL),'  : A%VAL(',ICOL,')=', A%VAL(ICOL)
WRITE(MSG%LU_DEBUG,*) '    A%COLG(',ICOL,')=',A%COLG(ICOL),'  : A%VAL(',ICOL,')=', A%VAL(ICOL)
WRITE(MSG%LU_DEBUG,*) '    AUX_COL =', AUX_COL(1:NCOL)
WRITE(MSG%LU_DEBUG,*) '    JCOL =', JCOL
WRITE(MSG%LU_DEBUG,*) '    KCOL =', KCOL
#endif
      ICOL = ICOL + 1
      DO WHILE (KCOL < NCOL)
         A%COLG(ICOL) = AUX_COLG(JCOL)
         A%COL(ICOL) = AUX_COL(JCOL)
         A%VAL(ICOL) = AUX_VAL(JCOL)
         IF (COLG_IS_DEFINED) THEN
            AUX_COLG(JCOL) = 99999999
            JCOL = MINLOC(AUX_COLG(1:NCOL), DIM=1)
         ELSE
            AUX_COL(JCOL) = 99999999
            JCOL = MINLOC(AUX_COL(1:NCOL), DIM=1)
         ENDIF
         KCOL = KCOL + 1
#ifdef WITH_SCARC_DEBUG2
WRITE(MSG%LU_DEBUG,*) '    ------------->B'
WRITE(MSG%LU_DEBUG,*) '    A%COL(',ICOL,')=',A%COL(ICOL),'  : A%VAL(',ICOL,')=', A%VAL(ICOL)
WRITE(MSG%LU_DEBUG,*) '    A%COLG(',ICOL,')=',A%COLG(ICOL),'  : A%VAL(',ICOL,')=', A%VAL(ICOL)
WRITE(MSG%LU_DEBUG,*) '    AUX_COL =', AUX_COL(1:NCOL)
WRITE(MSG%LU_DEBUG,*) '    JCOL =', JCOL
WRITE(MSG%LU_DEBUG,*) '    KCOL =', KCOL
#endif
         ICOL = ICOL + 1
      ENDDO

      IF (ICOL /= A%ROW(IC+1)) WRITE(*,*) 'ERROR IN RESORT_MATRIX_ROWS'

   ENDDO
#ifdef WITH_SCARC_DEBUG
   CALL SCARC_DEBUG_CMATRIX (A, 'A','AFTER RESORT_MATRIX_ROWS')
#endif
ENDDO

END SUBROUTINE SCARC_RESORT_MATRIX_ROWS

! ================================================================================================
! ================================================================================================
! Bundle of routines for the administration of different working spaces
!  - This includes allocation, deallocation and resizing
!  - Verbose messages are displayed to chid.logX
! ================================================================================================
! ================================================================================================

! ------------------------------------------------------------------------------------------------
! Allocate and initialize integer array of dimension 1
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_ALLOCATE_INT1(WORKSPACE, NL1, NR1, NINIT, CTEXT)
USE MEMORY_FUNCTIONS, ONLY: CHKMEMERR
INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: WORKSPACE
INTEGER, INTENT(IN) :: NL1, NR1, NINIT
CHARACTER(*), INTENT(IN) :: CTEXT
IF (.NOT.ALLOCATED(WORKSPACE)) THEN
   ALLOCATE (WORKSPACE(NL1:NR1), STAT=IERROR)
   CALL CHKMEMERR ('SCARC_ALLOCATE_INT1', CTEXT, IERROR)
   SELECT CASE (NINIT)
      CASE (NSCARC_INIT_UNDEF)
         WORKSPACE = NSCARC_UNDEF_INT
      CASE (NSCARC_INIT_ZERO)
         WORKSPACE = NSCARC_ZERO_INT
      CASE (NSCARC_INIT_ONE)
         WORKSPACE = NSCARC_ONE_INT
   END SELECT
#ifdef WITH_SCARC_VERBOSE
   WRITE(MSG%LU_VERBOSE,1000) CTEXT, NL1, NR1
#endif
ELSE
   IF (SIZE(WORKSPACE,1) /= NR1-NL1+1) CALL SCARC_SHUTDOWN(NSCARC_ERROR_VECTOR_LENGTH, CTEXT, NSCARC_NONE)
ENDIF
#ifdef WITH_SCARC_VERBOSE
1000 FORMAT('Allocating INT1     array  ',A20,'(',I8,':',I8,')')
#endif
END SUBROUTINE SCARC_ALLOCATE_INT1


! ------------------------------------------------------------------------------------------------
! Allocate and initialize integer array of dimension 2
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_ALLOCATE_INT2(WORKSPACE, NL1, NR1, NL2, NR2, NINIT, CTEXT)
USE MEMORY_FUNCTIONS, ONLY: CHKMEMERR
INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: WORKSPACE
INTEGER, INTENT(IN) :: NL1, NR1, NL2, NR2, NINIT
CHARACTER(*), INTENT(IN) :: CTEXT
IF (.NOT.ALLOCATED(WORKSPACE)) THEN
   ALLOCATE (WORKSPACE(NL1:NR1, NL2:NR2), STAT=IERROR)
   CALL CHKMEMERR ('SCARC_ALLOCATE_INT2', CTEXT, IERROR)
   SELECT CASE (NINIT)
      CASE (NSCARC_INIT_UNDEF)
         WORKSPACE = NSCARC_UNDEF_INT
      CASE (NSCARC_INIT_ZERO)
         WORKSPACE = NSCARC_ZERO_INT
      CASE (NSCARC_INIT_ONE)
         WORKSPACE = NSCARC_ONE_INT
   END SELECT
#ifdef WITH_SCARC_VERBOSE
   WRITE(MSG%LU_VERBOSE,1000) CTEXT, NL1, NR1, NL2, NR2
#endif
ELSE
   IF (SIZE(WORKSPACE,1) /= NR1-NL1+1 .OR. &
      SIZE(WORKSPACE,2) /= NR2-NL2+1) THEN
      CALL SCARC_SHUTDOWN(NSCARC_ERROR_VECTOR_LENGTH, CTEXT, NSCARC_NONE)
   ENDIF
ENDIF
#ifdef WITH_SCARC_VERBOSE
1000 FORMAT('Allocating INT2     array  ',A20,'(',I8,':',I8,' , ',I8,':',I8,')')
#endif
END SUBROUTINE SCARC_ALLOCATE_INT2


! ------------------------------------------------------------------------------------------------
! Allocate and initialize integer array of dimension 3
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_ALLOCATE_INT3(WORKSPACE, NL1, NR1, NL2, NR2, NL3, NR3, NINIT, CTEXT)
USE MEMORY_FUNCTIONS, ONLY: CHKMEMERR
INTEGER, DIMENSION(:,:,:), ALLOCATABLE, INTENT(INOUT) :: WORKSPACE
INTEGER, INTENT(IN) :: NL1, NR1, NL2, NR2, NL3, NR3, NINIT
CHARACTER(*), INTENT(IN) :: CTEXT
IF (.NOT.ALLOCATED(WORKSPACE)) THEN
   ALLOCATE (WORKSPACE(NL1:NR1, NL2:NR2, NL3:NR3), STAT=IERROR)
   CALL CHKMEMERR ('SCARC_ALLOCATE_INT3', CTEXT, IERROR)
   SELECT CASE (NINIT)
      CASE (NSCARC_INIT_UNDEF)
         WORKSPACE = NSCARC_UNDEF_INT
      CASE (NSCARC_INIT_ZERO)
         WORKSPACE = NSCARC_ZERO_INT
      CASE (NSCARC_INIT_ONE)
         WORKSPACE = NSCARC_ONE_INT
   END SELECT
#ifdef WITH_SCARC_VERBOSE
   WRITE(MSG%LU_VERBOSE,1000) CTEXT, NL1, NR1, NL2, NR2, NL3, NR3
#endif
ELSE
   IF (SIZE(WORKSPACE,1) /= NR1-NL1+1 .OR. &
      SIZE(WORKSPACE,2) /= NR2-NL2+1 .OR. &
      SIZE(WORKSPACE,3) /= NR3-NL3+1) THEN
      CALL SCARC_SHUTDOWN(NSCARC_ERROR_VECTOR_LENGTH, CTEXT, NSCARC_NONE)
   ENDIF
ENDIF
#ifdef WITH_SCARC_VERBOSE
1000 FORMAT('Allocating INT3     array  ',A20,'(',I8,':',I8,' , ',I8,':',I8,' , ',I8,':',I8,')')
#endif
END SUBROUTINE SCARC_ALLOCATE_INT3


! ------------------------------------------------------------------------------------------------
! Allocate and initialize Logical array of dimension 1
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_ALLOCATE_LOG1(WORKSPACE, NL1, NR1, NINIT, CTEXT)
USE MEMORY_FUNCTIONS, ONLY: CHKMEMERR
LOGICAL, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: WORKSPACE
INTEGER, INTENT(IN) :: NL1, NR1, NINIT
CHARACTER(*), INTENT(IN) :: CTEXT
IF (.NOT.ALLOCATED(WORKSPACE)) THEN
   ALLOCATE (WORKSPACE(NL1:NR1), STAT=IERROR)
   CALL CHKMEMERR ('SCARC_ALLOCATE_INT1', CTEXT, IERROR)
   SELECT CASE (NINIT)
      CASE (NSCARC_INIT_TRUE)
         WORKSPACE = .TRUE.
      CASE (NSCARC_INIT_FALSE)
         WORKSPACE = .FALSE.
   END SELECT
#ifdef WITH_SCARC_VERBOSE
   WRITE(MSG%LU_VERBOSE,1000) CTEXT, NL1, NR1
#endif
ELSE
   IF (SIZE(WORKSPACE,1) /= NR1-NL1+1) THEN
      CALL SCARC_SHUTDOWN(NSCARC_ERROR_VECTOR_LENGTH, CTEXT, NSCARC_NONE)
   ENDIF
ENDIF
#ifdef WITH_SCARC_VERBOSE
1000 FORMAT('Allocating LOG3     array  ',A20,'(',I8,':',I8,')')
#endif
END SUBROUTINE SCARC_ALLOCATE_LOG1


! ------------------------------------------------------------------------------------------------
! Allocate and initialize Logical array of dimension 2
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_ALLOCATE_LOG2(WORKSPACE, NL1, NR1, NL2, NR2, NINIT, CTEXT)
USE MEMORY_FUNCTIONS, ONLY: CHKMEMERR
LOGICAL, DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: WORKSPACE
INTEGER, INTENT(IN) :: NL1, NR1, NL2, NR2, NINIT
CHARACTER(*), INTENT(IN) :: CTEXT
IF (.NOT.ALLOCATED(WORKSPACE)) THEN
   ALLOCATE (WORKSPACE(NL1:NR1, NL2:NR2), STAT=IERROR)
   CALL CHKMEMERR ('SCARC_ALLOCATE_INT3', CTEXT, IERROR)
   SELECT CASE (NINIT)
      CASE (NSCARC_INIT_TRUE)
         WORKSPACE = .TRUE.
      CASE (NSCARC_INIT_FALSE)
         WORKSPACE = .FALSE.
   END SELECT
#ifdef WITH_SCARC_VERBOSE
   WRITE(MSG%LU_VERBOSE,1000) CTEXT, NL1, NR1, NL2, NR2
#endif
ELSE
   IF (SIZE(WORKSPACE,1) /= NR1-NL1+1 .OR. &
      SIZE(WORKSPACE,2) /= NR2-NL2+1) THEN
      CALL SCARC_SHUTDOWN(NSCARC_ERROR_VECTOR_LENGTH, CTEXT, NSCARC_NONE)
   ENDIF
ENDIF
#ifdef WITH_SCARC_VERBOSE
1000 FORMAT('Allocating LOG3     array  ',A20,'(',I8,':',I8,' , ',I8,':',I8,')')
#endif
END SUBROUTINE SCARC_ALLOCATE_LOG2


! ------------------------------------------------------------------------------------------------
! Allocate and initialize Logical array of dimension 3
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_ALLOCATE_LOG3(WORKSPACE, NL1, NR1, NL2, NR2, NL3, NR3, NINIT, CTEXT)
USE MEMORY_FUNCTIONS, ONLY: CHKMEMERR
LOGICAL, DIMENSION(:,:,:), ALLOCATABLE, INTENT(INOUT) :: WORKSPACE
INTEGER, INTENT(IN) :: NL1, NR1, NL2, NR2, NL3, NR3, NINIT
CHARACTER(*), INTENT(IN) :: CTEXT
IF (.NOT.ALLOCATED(WORKSPACE)) THEN
   ALLOCATE (WORKSPACE(NL1:NR1, NL2:NR2, NL3:NR3), STAT=IERROR)
   CALL CHKMEMERR ('SCARC_ALLOCATE_INT3', CTEXT, IERROR)
   SELECT CASE (NINIT)
      CASE (NSCARC_INIT_TRUE)
         WORKSPACE = .TRUE.
      CASE (NSCARC_INIT_FALSE)
         WORKSPACE = .FALSE.
   END SELECT
#ifdef WITH_SCARC_VERBOSE
   WRITE(MSG%LU_VERBOSE,1000) CTEXT, NL1, NR1, NL2, NR2, NL3, NR3
#endif
ELSE
   IF (SIZE(WORKSPACE,1) /= NR1-NL1+1 .OR. &
      SIZE(WORKSPACE,2) /= NR2-NL2+1 .OR. &
      SIZE(WORKSPACE,3) /= NR3-NL3+1) THEN
      CALL SCARC_SHUTDOWN(NSCARC_ERROR_VECTOR_LENGTH, CTEXT, NSCARC_NONE)
   ENDIF
ENDIF
#ifdef WITH_SCARC_VERBOSE
1000 FORMAT('Allocating LOG3     array  ',A20,'(',I8,':',I8,' , ',I8,':',I8,' , ',I8,':',I8,')')
#endif
END SUBROUTINE SCARC_ALLOCATE_LOG3


! ------------------------------------------------------------------------------------------------
! Allocate and initialize real array of dimension 1
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_ALLOCATE_REAL1(WORKSPACE, NL1, NR1, NINIT, CTEXT)
USE PRECISION_PARAMETERS, ONLY: EB
USE MEMORY_FUNCTIONS, ONLY: CHKMEMERR
REAL(EB), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: WORKSPACE
INTEGER, INTENT(IN) :: NL1, NR1, NINIT
CHARACTER(*), INTENT(IN) :: CTEXT

IF (.NOT.ALLOCATED(WORKSPACE)) THEN

   ALLOCATE (WORKSPACE(NL1:NR1), STAT=IERROR)
   CALL CHKMEMERR ('SCARC_ALLOCATE_REAL1', CTEXT, IERROR)

   SELECT CASE (NINIT)
      CASE (NSCARC_INIT_UNDEF)
         WORKSPACE = NSCARC_UNDEF_REAL_EB
      CASE (NSCARC_INIT_ZERO)
         WORKSPACE = NSCARC_ZERO_REAL_EB
      CASE (NSCARC_INIT_ONE)
         WORKSPACE = NSCARC_ONE_REAL_EB
   END SELECT

#ifdef WITH_SCARC_VERBOSE
   WRITE(MSG%LU_VERBOSE,1000) CTEXT, NL1, NR1
#endif
ELSE
   IF (SIZE(WORKSPACE,1) /= NR1-NL1+1) &
      CALL SCARC_SHUTDOWN(NSCARC_ERROR_VECTOR_LENGTH, CTEXT, NSCARC_NONE)
ENDIF

#ifdef WITH_SCARC_VERBOSE
1000 FORMAT('Allocating REAL1    array  ',A20,'(',I8,':',I8,')')
#endif
END SUBROUTINE SCARC_ALLOCATE_REAL1


! ------------------------------------------------------------------------------------------------
! Allocate and initialize real array of dimension 1
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_ALLOCATE_REAL1_FB(WORKSPACE, NL1, NR1, NINIT, CTEXT)
USE PRECISION_PARAMETERS, ONLY: FB
USE MEMORY_FUNCTIONS, ONLY: CHKMEMERR
REAL(FB), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: WORKSPACE
INTEGER, INTENT(IN) :: NL1, NR1, NINIT
CHARACTER(*), INTENT(IN) :: CTEXT
IF (.NOT.ALLOCATED(WORKSPACE)) THEN
   ALLOCATE (WORKSPACE(NL1:NR1), STAT=IERROR)
   CALL CHKMEMERR ('SCARC_ALLOCATE_REAL1_FB', CTEXT, IERROR)
   SELECT CASE (NINIT)
      CASE (NSCARC_INIT_UNDEF)
         WORKSPACE = NSCARC_UNDEF_REAL_FB
      CASE (NSCARC_INIT_ZERO)
         WORKSPACE = NSCARC_ZERO_REAL_FB
      CASE (NSCARC_INIT_ONE)
         WORKSPACE = NSCARC_ONE_REAL_FB
   END SELECT
#ifdef WITH_SCARC_VERBOSE
   WRITE(MSG%LU_VERBOSE,1000) CTEXT, NL1, NR1
#endif
ELSE
   IF (SIZE(WORKSPACE,1) /= NR1-NL1+1) &
      CALL SCARC_SHUTDOWN(NSCARC_ERROR_VECTOR_LENGTH, CTEXT, NSCARC_NONE)
ENDIF
#ifdef WITH_SCARC_VERBOSE
1000 FORMAT('Allocating REAL1_FB array  ',A20,'(',I8,':',I8,')')
#endif
END SUBROUTINE SCARC_ALLOCATE_REAL1_FB


! ------------------------------------------------------------------------------------------------
! Allocate and initialize real array of dimension 2
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_ALLOCATE_REAL2(WORKSPACE, NL1, NR1, NL2, NR2, NINIT, CTEXT)
USE PRECISION_PARAMETERS, ONLY: EB
USE MEMORY_FUNCTIONS, ONLY: CHKMEMERR
REAL(EB), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: WORKSPACE
INTEGER, INTENT(IN) :: NL1, NR1, NL2, NR2, NINIT
CHARACTER(*), INTENT(IN) :: CTEXT
IF (.NOT.ALLOCATED(WORKSPACE)) THEN
   ALLOCATE (WORKSPACE(NL1:NR1, NL2:NR2), STAT=IERROR)
   CALL CHKMEMERR ('SCARC_ALLOCATE_REAL2', CTEXT, IERROR)
   SELECT CASE (NINIT)
      CASE (NSCARC_INIT_UNDEF)
         WORKSPACE = NSCARC_UNDEF_REAL_EB
      CASE (NSCARC_INIT_ZERO)
         WORKSPACE = NSCARC_ZERO_REAL_EB
      CASE (NSCARC_INIT_ONE)
         WORKSPACE = NSCARC_ONE_REAL_EB
   END SELECT
#ifdef WITH_SCARC_VERBOSE
   WRITE(MSG%LU_VERBOSE,1000) CTEXT, NL1, NR1, NL2, NR2
#endif
ELSE
   IF (SIZE(WORKSPACE,1) /= NR1-NL1+1 .OR. &
      SIZE(WORKSPACE,2) /= NR2-NL2+1) THEN
      CALL SCARC_SHUTDOWN(NSCARC_ERROR_VECTOR_LENGTH, CTEXT, NSCARC_NONE)
   ENDIF
ENDIF
#ifdef WITH_SCARC_VERBOSE
1000 FORMAT('Allocating REAL2    array  ',A20,'(',I8,':',I8,' , ',I8,':',I8,')')
#endif
END SUBROUTINE SCARC_ALLOCATE_REAL2


! ------------------------------------------------------------------------------------------------
! Allocate and initialize real array of dimension 3
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_ALLOCATE_REAL3(WORKSPACE, NL1, NR1, NL2, NR2, NL3, NR3, NINIT, CTEXT)
USE PRECISION_PARAMETERS, ONLY: EB
USE MEMORY_FUNCTIONS, ONLY: CHKMEMERR
REAL(EB), DIMENSION(:,:,:), ALLOCATABLE, INTENT(INOUT) :: WORKSPACE
INTEGER, INTENT(IN) :: NL1, NR1, NL2, NR2, NL3, NR3, NINIT
CHARACTER(*), INTENT(IN) :: CTEXT
IF (.NOT.ALLOCATED(WORKSPACE)) THEN
   ALLOCATE (WORKSPACE(NL1:NR1, NL2:NR2, NL3:NR3), STAT=IERROR)
   CALL CHKMEMERR ('SCARC_ALLOCATE_REAL3', CTEXT, IERROR)
   SELECT CASE (NINIT)
      CASE (NSCARC_INIT_UNDEF)
         WORKSPACE = NSCARC_UNDEF_REAL_EB
      CASE (NSCARC_INIT_ZERO)
         WORKSPACE = NSCARC_ZERO_REAL_EB
      CASE (NSCARC_INIT_ONE)
         WORKSPACE = NSCARC_ONE_REAL_EB
   END SELECT
#ifdef WITH_SCARC_VERBOSE
   WRITE(MSG%LU_VERBOSE,1000) CTEXT, NL1, NR1, NL2, NR2, NL3, NR3
#endif
ELSE
   IF (SIZE(WORKSPACE,1) /= NR1-NL1+1 .OR. &
      SIZE(WORKSPACE,2) /= NR2-NL2+1 .OR. &
      SIZE(WORKSPACE,3) /= NR3-NL3+1) THEN
      CALL SCARC_SHUTDOWN(NSCARC_ERROR_VECTOR_LENGTH, CTEXT, NSCARC_NONE)
   ENDIF
ENDIF
#ifdef WITH_SCARC_VERBOSE
1000 FORMAT('Allocating REAL3    array  ',A20,'(',I8,':',I8,' , ',I8,':',I8,' , ',I8,':',I8,')')
#endif
END SUBROUTINE SCARC_ALLOCATE_REAL3

! ------------------------------------------------------------------------------------------------
! Reduce size of integer vector
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_REDUCE_INT1(WORKSPACE, NSIZE, CNAME)
INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(INOUT) :: WORKSPACE
INTEGER, INTENT(IN) :: NSIZE
CHARACTER(*), INTENT(IN) :: CNAME
INTEGER, ALLOCATABLE, DIMENSION(:) :: AUX

IF (NSIZE == SIZE(WORKSPACE)) THEN
   RETURN                                  ! vector already has requested size
ELSE IF (NSIZE < SIZE(WORKSPACE)) THEN

   CALL SCARC_ALLOCATE_INT1(AUX, 1, NSIZE, NSCARC_INIT_NONE, TRIM(CNAME))
   AUX(1:NSIZE) = WORKSPACE(1:NSIZE)
   DEALLOCATE(WORKSPACE)

   CALL SCARC_ALLOCATE_INT1(WORKSPACE, 1, NSIZE, NSCARC_INIT_NONE, TRIM(CNAME))
   WORKSPACE(1:NSIZE) = AUX(1:NSIZE)
   DEALLOCATE(AUX)
ENDIF

#ifdef WITH_SCARC_VERBOSE
WRITE(MSG%LU_VERBOSE,1000) CNAME, NSIZE
1000 FORMAT('Resizing   integer vector ',A20,' to size ',I8)
#endif
END SUBROUTINE SCARC_REDUCE_INT1

! ------------------------------------------------------------------------------------------------
! Reduce size of integer array with dimension 2
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_REDUCE_INT2(WORKSPACE, NSIZE1, NSIZE2, CNAME)
INTEGER, ALLOCATABLE, DIMENSION(:,:), INTENT(INOUT) :: WORKSPACE
INTEGER, INTENT(IN) :: NSIZE1, NSIZE2
CHARACTER(*), INTENT(IN) :: CNAME
INTEGER :: NWORK1, NWORK2, I1, I2
INTEGER, ALLOCATABLE, DIMENSION(:,:) :: AUX

NWORK1 = SIZE(WORKSPACE, DIM = 1)
NWORK2 = SIZE(WORKSPACE, DIM = 2)
IF (NSIZE1 == NWORK1 .AND. NSIZE2 == NWORK2 ) THEN
   RETURN                                  ! vector already has requested size
ELSE IF (NSIZE1 <= NWORK1 .AND. NSIZE2 <= NWORK2) THEN

   CALL SCARC_ALLOCATE_INT2(AUX, 1, NWORK1, 1, NWORK2, NSCARC_INIT_NONE, TRIM(CNAME))
   DO I2 = 1, NWORK2
      DO I1 = 1, NWORK1
         AUX(I1, I2) = WORKSPACE(I1, I2)
      ENDDO
   ENDDO
   DEALLOCATE(WORKSPACE)

   CALL SCARC_ALLOCATE_INT2(WORKSPACE, 1, NSIZE1, 1, NSIZE2, NSCARC_INIT_NONE, TRIM(CNAME))
   DO I2 = 1, NSIZE2
      DO I1 = 1, NSIZE1
         WORKSPACE(I1, I2) = AUX(I1, I2)
      ENDDO
   ENDDO
   DEALLOCATE(AUX)
ELSE
   WRITE(*,'(A,2I6,A,2I6)') 'Error in SCARC_REDUCE_INT2, orig ', NWORK1, NWORK2,': new ',NSIZE1, NSIZE2
ENDIF

#ifdef WITH_SCARC_VERBOSE
WRITE(MSG%LU_VERBOSE,1000) CNAME, NSIZE1, NSIZE2
1000 FORMAT('Resizing  array ',A20,' to size ',I8, 'x', I8)
#endif
END SUBROUTINE SCARC_REDUCE_INT2


! ------------------------------------------------------------------------------------------------
! Reduce size of integer vector
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_EXPAND_INT1(WORKSPACE, WORKSPACE_ADD, NSIZE, NSIZE_ADD, CNAME)
INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(INOUT) :: WORKSPACE, WORKSPACE_ADD
INTEGER, INTENT(IN) :: NSIZE, NSIZE_ADD
CHARACTER(*), INTENT(IN) :: CNAME
INTEGER, ALLOCATABLE, DIMENSION(:) :: AUX

CALL SCARC_ALLOCATE_INT1(AUX, 1, NSIZE + NSIZE_ADD, NSCARC_INIT_NONE, TRIM(CNAME))
AUX(1:NSIZE) = WORKSPACE(1:NSIZE)
AUX(NSIZE+1:NSIZE+NSIZE_ADD) = WORKSPACE_ADD(NSIZE+1:NSIZE+NSIZE_ADD)
DEALLOCATE(WORKSPACE)

CALL SCARC_ALLOCATE_INT1(WORKSPACE, 1, NSIZE + NSIZE_ADD, NSCARC_INIT_NONE, TRIM(CNAME))
WORKSPACE(1:NSIZE+NSIZE_ADD) = AUX(1:NSIZE+NSIZE_ADD)
DEALLOCATE(AUX)

#ifdef WITH_SCARC_VERBOSE
WRITE(MSG%LU_VERBOSE,1000) CNAME, NSIZE + NSIZE_ADD
1000 FORMAT('Resizing   integer vector ',A20,' to size ',I8)
#endif
END SUBROUTINE SCARC_EXPAND_INT1


! ------------------------------------------------------------------------------------------------
! Reduce size of integer vector
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_REDUCE_REAL1(WORKSPACE, NSIZE, CNAME)
REAL(EB), ALLOCATABLE, DIMENSION(:), INTENT(INOUT) :: WORKSPACE
INTEGER, INTENT(IN) :: NSIZE
CHARACTER(*), INTENT(IN) :: CNAME
REAL(EB), ALLOCATABLE, DIMENSION(:) :: AUX

IF (NSIZE == SIZE(WORKSPACE)) THEN
   RETURN                                  ! vector already has requested size
ELSE IF (NSIZE < SIZE(WORKSPACE)) THEN

   CALL SCARC_ALLOCATE_REAL1(AUX, 1, NSIZE, NSCARC_INIT_NONE, TRIM(CNAME))
   AUX(1:NSIZE) = WORKSPACE(1:NSIZE)
   DEALLOCATE(WORKSPACE)

   CALL SCARC_ALLOCATE_REAL1(WORKSPACE, 1, NSIZE, NSCARC_INIT_NONE, TRIM(CNAME))
   WORKSPACE(1:NSIZE) = AUX(1:NSIZE)
   DEALLOCATE(AUX)
ENDIF

#ifdef WITH_SCARC_VERBOSE
WRITE(MSG%LU_VERBOSE,1000) CNAME, NSIZE
1000 FORMAT('Resizing   integer vector ',A20,' to size ',I8)
#endif
END SUBROUTINE SCARC_REDUCE_REAL1


! ------------------------------------------------------------------------------------------------
! COMPACT storage format:
! Allocate matrix with corresponding pointer and length structures
!    NTYPE == NSCARC_MATRIX_FULL    :  ALLOCATE VAL, COL and COLG
!    NTYPE == NSCARC_MATRIX_LIGHT   :  ALLOCATE VAL, COL 
!    NTYPE == NSCARC_MATRIX_MINIMAL :  ALLOCATE COL 
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_ALLOCATE_CMATRIX(A, NL, NPREC, NTYPE, CNAME)
TYPE (SCARC_MATRIX_COMPACT_TYPE), INTENT(INOUT) :: A
INTEGER, INTENT(IN) :: NPREC, NTYPE, NL
CHARACTER(*), INTENT(IN) :: CNAME

A%CNAME = TRIM(CNAME)
A%NTYPE = NTYPE
#ifdef WITH_SCARC_VERBOSE
WRITE(MSG%LU_VERBOSE,1000) CNAME, NL
#endif

CALL SCARC_ALLOCATE_INT1(A%ROW, 1, A%N_ROW, NSCARC_INIT_ZERO, 'A%ROW')
CALL SCARC_ALLOCATE_INT1(A%COL, 1, A%N_VAL, NSCARC_INIT_ZERO, 'A%COL')

IF (NTYPE == NSCARC_MATRIX_LIGHT .OR. NTYPE == NSCARC_MATRIX_FULL) &
   CALL SCARC_ALLOCATE_INT1(A%COLG, 1, A%N_VAL, NSCARC_INIT_ZERO, 'A%COLG')

IF (NTYPE /= NSCARC_MATRIX_MINIMAL) THEN
   IF (NPREC == NSCARC_PRECISION_SINGLE) THEN
      CALL SCARC_ALLOCATE_REAL1_FB(A%VAL_FB, 1, A%N_VAL, NSCARC_INIT_ZERO, 'A%VAL_FB')
   ELSE
      CALL SCARC_ALLOCATE_REAL1(A%VAL, 1, A%N_VAL, NSCARC_INIT_ZERO, 'A%VAL')
   ENDIF
ENDIF

#ifdef WITH_SCARC_VERBOSE
1000 FORMAT('Allocating COMPACT  matrix ',A20,' on level ',I6)
#endif
END SUBROUTINE SCARC_ALLOCATE_CMATRIX


! ------------------------------------------------------------------------------------------------
! COMPACT storage format:
! Deallocate matrix with corresponding pointer and length structures
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_DEALLOCATE_CMATRIX(A, CNAME)
TYPE (SCARC_MATRIX_COMPACT_TYPE), INTENT(INOUT) :: A
CHARACTER(*) :: CNAME

A%N_STENCIL   = 0
A%N_CONDENSED = 0
A%N_ROW       = 0
A%N_VAL       = 0
A%STENCIL     = 0
A%POS         = 0

IF (ALLOCATED(A%VAL))   DEALLOCATE(A%VAL)
IF (ALLOCATED(A%ROW))   DEALLOCATE(A%ROW)
IF (ALLOCATED(A%COL))   DEALLOCATE(A%COL)
IF (ALLOCATED(A%COLG))  DEALLOCATE(A%COLG)
IF (ALLOCATED(A%RELAX)) DEALLOCATE(A%RELAX)

#ifdef WITH_SCARC_VERBOSE
WRITE(MSG%LU_VERBOSE,1000) CNAME
1000 FORMAT('Removing   COMPACT matrix ', A20)
#else
CNAME = ''
#endif
END SUBROUTINE SCARC_DEALLOCATE_CMATRIX


! ------------------------------------------------------------------------------------------------
! COMPACT storage format:
! Reduce sizes of VAL, COL and COLG arrays of a given matrix to real size
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_REDUCE_CMATRIX(A, CNAME)
TYPE (SCARC_MATRIX_COMPACT_TYPE), INTENT(INOUT) :: A
CHARACTER(*), INTENT(IN) :: CNAME
REAL(EB), ALLOCATABLE, DIMENSION(:) :: VAL
INTEGER , ALLOCATABLE, DIMENSION(:) :: COL, COLG
INTEGER :: NVAL, NVAL_ALLOCATED

NVAL = A%N_VAL
NVAL_ALLOCATED = SIZE(A%COL)

#ifdef WITH_SCARC_VERBOSE
WRITE(MSG%LU_VERBOSE, 1000, ADVANCE='NO') CNAME, SIZE(A%COL), NVAL, A%N_ROW
#endif

! If the matrix already has the desired size or specified values are to small, return or shutdown
IF (SIZE(A%COL) == NVAL) THEN
#ifdef WITH_SCARC_VERBOSE
   WRITE(MSG%LU_VERBOSE,*) ' ...  nothing to do '
#endif
   RETURN
ELSE IF (SIZE(A%COL) < NVAL) THEN
#ifdef WITH_SCARC_VERBOSE
   WRITE(MSG%LU_VERBOSE,*) ' ...  failed '
#endif
   CALL SCARC_SHUTDOWN(NSCARC_ERROR_MATRIX_SIZE, CNAME, NSCARC_NONE)
ENDIF

! If the allocated size of the matrix values workspace is too large, reduce it to the real size
IF (NVAL < SIZE(A%COL)) THEN

#ifdef WITH_SCARC_VERBOSE
   WRITE(MSG%LU_VERBOSE,*) 
#endif

   CALL SCARC_ALLOCATE_INT1(COL,   1, NVAL, NSCARC_INIT_NONE, 'COL')
   COL(1:NVAL) = A%COL(1:NVAL)
   DEALLOCATE(A%COL)
   CALL SCARC_ALLOCATE_INT1(A%COL, 1, NVAL, NSCARC_INIT_NONE, 'A%COL')
   A%COL(1:NVAL) = COL(1:NVAL)
   DEALLOCATE(COL)

   IF (ALLOCATED(A%COLG)) THEN
      CALL SCARC_ALLOCATE_INT1(COLG,   1, NVAL, NSCARC_INIT_NONE, 'COLG')
      COLG(1:NVAL) = A%COLG(1:NVAL)
      DEALLOCATE(A%COLG)
      CALL SCARC_ALLOCATE_INT1(A%COLG, 1, NVAL, NSCARC_INIT_NONE, 'A%COLG')
      A%COLG(1:NVAL) = COLG(1:NVAL)
      DEALLOCATE(COLG)
   ENDIF

   IF (A%NTYPE /= NSCARC_MATRIX_MINIMAL) THEN
      CALL SCARC_ALLOCATE_REAL1(VAL,   1, NVAL, NSCARC_INIT_NONE, 'VAL')
      VAL(1:NVAL) = A%VAL(1:NVAL)
      DEALLOCATE(A%VAL)
      CALL SCARC_ALLOCATE_REAL1(A%VAL, 1, NVAL, NSCARC_INIT_NONE, 'A%VAL')
      A%VAL(1:NVAL) = VAL(1:NVAL)
      DEALLOCATE(VAL)
   ENDIF

ENDIF

#ifdef WITH_SCARC_VERBOSE
1000 FORMAT('Reducing compact matrix',A20,' from ', I8,' to ', I8, ': N_ROW = ', I8)
#endif
END SUBROUTINE SCARC_REDUCE_CMATRIX


! ------------------------------------------------------------------------------------------------
! BANDED storage format:
! Allocate matrix with corresponding pointer and length structures
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_ALLOCATE_BMATRIX(A, NL, CNAME)
TYPE (SCARC_MATRIX_BANDWISE_TYPE), INTENT(INOUT) :: A
CHARACTER(*), INTENT(IN) :: CNAME
INTEGER, INTENT(IN) :: NL
CHARACTER(40) :: CINFO

A%CNAME = CNAME

WRITE(CINFO,'(A,A,I2.2,A)') TRIM(CNAME),'_LEV',NL,'.AUX'
CALL SCARC_ALLOCATE_REAL1(A%AUX, 1, A%N_DIAG, NSCARC_INIT_ZERO, CINFO)
WRITE(CINFO,'(A,A,I2.2,A)') TRIM(CNAME),'_LEV',NL,'.VAL'
CALL SCARC_ALLOCATE_REAL2(A%VAL, 1, A%N_DIAG, 1, A%N_STENCIL, NSCARC_INIT_ZERO, CINFO)

#ifdef WITH_SCARC_VERBOSE
WRITE(MSG%LU_VERBOSE,1000) CNAME, A%N_VAL, A%N_DIAG
1000 FORMAT('Allocating BANDED  matrix ',A20,' with NA=',I8,' and NDIAG=',I8)
#endif
END SUBROUTINE SCARC_ALLOCATE_BMATRIX


! ------------------------------------------------------------------------------------------------
! BANDED storage format:
! Deallocate matrix with corresponding pointer and length structures
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_DEALLOCATE_BMATRIX(A)
TYPE (SCARC_MATRIX_BANDWISE_TYPE), INTENT(INOUT) :: A

A%N_STENCIL   = 0
A%N_CONDENSED = 0
A%N_VAL       = 0
A%N_DIAG      = 0

IF (ALLOCATED(A%AUX))    DEALLOCATE(A%AUX)
IF (ALLOCATED(A%VAL))    DEALLOCATE(A%VAL)
IF (ALLOCATED(A%RELAX))  DEALLOCATE(A%RELAX)
IF (ALLOCATED(A%RELAXD)) DEALLOCATE(A%RELAXD)
A%STENCIL = 0

A%OFFSET = 0
A%LENGTH = 0
A%SOURCE = 0
A%TARGET = 0

#ifdef WITH_SCARC_VERBOSE
WRITE(MSG%LU_VERBOSE,1000)
1000 FORMAT('Removing   COMPACT matrix ')
#endif
END SUBROUTINE SCARC_DEALLOCATE_BMATRIX



! ====================================================================================================
! ====================================================================================================
! Bundle of routines for the output of different debugging information
! ====================================================================================================
! ====================================================================================================

! ----------------------------------------------------------------------------------------------------
! Dump residual information
! ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_DUMP_CSV(ISM, NS, NL)
USE SCARC_ITERATION_ENVIRONMENT
INTEGER, INTENT(IN) :: ISM, NS, NL

IF (.NOT.HAS_CSV_DUMP .OR. MYID /= 0) RETURN
IF (ITE_TOTAL == 0 .AND. TYPE_SOLVER /= NSCARC_SOLVER_MAIN) RETURN
IF (TYPE_SOLVER == NSCARC_SOLVER_COARSE) RETURN
WRITE(MSG%LU_STAT,1000) ITE_PRES, NS, ITE_TOTAL, ITE_CG, ITE_MG, NL, ITE_SMOOTH, ISM, ITE_COARSE, ITE_LU, RES, CAPPA

1000 FORMAT(10(I8,','), E12.4,',',E12.4)
END SUBROUTINE SCARC_DUMP_CSV


! ----------------------------------------------------------------------------------------------------
! Dump CPU times
! ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_DUMP_TIMERS
INTEGER, PARAMETER :: LINE_LENGTH = 5 + 12*11
INTEGER :: N, STATUS(MPI_STATUS_SIZE)
CHARACTER(LEN=LINE_LENGTH) :: LINE
CHARACTER(LEN=LINE_LENGTH), DIMENSION(0:N_MPI_PROCESSES-1) :: LINE_ARRAY

! All MPI processes except root send their timings to the root process. The root process then writes them out to a file.
WRITE(LINE,'(I5,12(",",ES10.3))') MYID,                       &
                                  CPU(MYID)%OVERALL,          &
                                  CPU(MYID)%SETUP,            &
                                  CPU(MYID)%SOLVER,           &
                                  CPU(MYID)%ITERATION,        &
                                  CPU(MYID)%MATVEC_PRODUCT,   &
                                  CPU(MYID)%SCALAR_PRODUCT,   &
                                  CPU(MYID)%RELAXATION,       &
                                  CPU(MYID)%SMOOTHER,         &
                                  CPU(MYID)%COARSE,           &
                                  CPU(MYID)%EXCHANGE,         &
                                  CPU(MYID)%BUFFER_PACKING,   &
                                  CPU(MYID)%BUFFER_UNPACKING

IF (MYID>0) THEN
   CALL MPI_SEND(LINE,LINE_LENGTH,MPI_CHARACTER,0,MYID,MPI_COMM_WORLD,IERROR)
ELSE
   LINE_ARRAY(0) = LINE
   DO N=1,N_MPI_PROCESSES-1
      CALL MPI_RECV(LINE_ARRAY(N),LINE_LENGTH,MPI_CHARACTER,N,N,MPI_COMM_WORLD,STATUS,IERROR)
   ENDDO
   MSG%FILE_CPU = TRIM(CHID)//'_scarc_cpu.csv'
   OPEN (MSG%LU_CPU, FILE=MSG%FILE_CPU, STATUS='REPLACE',FORM='FORMATTED')
   WRITE(MSG%LU_CPU,'(A,A)') 'Rank,OVERALL,SETUP,SOLVER,ITERATION,MATVEC_PRODUCT,SCALAR_PRODUCT,',&
                             'RELAXATION,SMOOTHER,COARSE,EXCHANGE,PACKING,UNPACKING'
   DO N=0,N_MPI_PROCESSES-1
      WRITE(MSG%LU_CPU,'(A)') LINE_ARRAY(N)
   ENDDO
   CLOSE(MSG%LU_CPU)
ENDIF

END SUBROUTINE SCARC_DUMP_TIMERS



#ifdef WITH_SCARC_MGM
! ====================================================================================================
! Begin  SCARC_WITH_MGM - Part
! Collection of routines to experiment with McKeeney-Greengard-Mayo method 
! ====================================================================================================
! ----------------------------------------------------------------------------------------------------
! Store preliminary solution vector in MGM structure
! ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_STORE_MGM(NL, NTYPE)
USE SCARC_POINTERS, ONLY: ST, MGM
INTEGER, INTENT(IN) :: NL, NTYPE
INTEGER :: NM

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   ST=> SCARC(NM)%LEVEL(NL)%STAGE(NSCARC_STAGE_ONE)
   MGM => SCARC(NM)%LEVEL(NL)%MGM

   SELECT CASE(NTYPE)
      CASE (1)
         MGM%H1 = ST%X
      CASE (2)
         MGM%H2 = ST%X
      CASE (3)
         ST%X = MGM%H1 + MGM%H2
   END SELECT

ENDDO

END SUBROUTINE SCARC_STORE_MGM


! ------------------------------------------------------------------------------------------------
! Set internal boundary conditions for unstructured, homogeneous part of MGM method
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_MGM_INTERNAL_VELOCITY(NL)
USE SCARC_POINTERS, ONLY: M, L, MGM, UU, VV, WW, HP
INTEGER, INTENT(IN) :: NL
INTEGER :: NM, IW, I, J, K

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   M   => MESHES(NM)
   L   => SCARC(NM)%LEVEL(NL)
   MGM => L%MGM

   IF (PREDICTOR) THEN
      UU => M%U
      VV => M%V
      WW => M%W
      HP => M%H
   ELSE
      UU => M%US
      VV => M%VS
      WW => M%WS
      HP => M%HS
   ENDIF

   DO IW = L%N_WALL_CELLS_EXT+1, L%N_WALL_CELLS_EXT + L%N_WALL_CELLS_INT

      GWC => L%STRUCTURED%WALL(IW)

      I = GWC%IXW
      J = GWC%IYW
      K = GWC%IZW

      MGM%US(IW) = UU(I,J,K) - DT*( M%FVX(I,J,K) + M%RDXN(I)*(HP(I+1,J  ,K  ) - HP(I,J,K)) )
      MGM%VS(IW) = VV(I,J,K) - DT*( M%FVY(I,J,K) + M%RDYN(J)*(HP(I  ,J+1,K  ) - HP(I,J,K)) )
      MGM%WS(IW) = WW(I,J,K) - DT*( M%FVZ(I,J,K) + M%RDZN(K)*(HP(I  ,J  ,K+1) - HP(I,J,K)) )

#ifdef WITH_SCARC_DEBUG
      WRITE(MSG%LU_DEBUG,'(a,4i5,5F12.6)') 'MGM: I,J,K,IW,US,WS,DT,RDXN, RDZN:', &
         I, J, K, IW, MGM%US(IW), MGM%WS(IW), DT, M%RDXN(I), M%RDZN(I)
#endif

   ENDDO

ENDDO

END SUBROUTINE SCARC_MGM_INTERNAL_VELOCITY


! ----------------------------------------------------------------------------------------------------
! Allocate velocity vectors for the setting of internal BC's in MGM-method
! ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_MGM(NLMIN, NLMAX)
USE SCARC_POINTERS, ONLY: L, G, MGM, AC
INTEGER, INTENT(IN) :: NLMIN, NLMAX
INTEGER :: NM, NL, IC, IW, IP
LOGICAL, DIMENSION(:), ALLOCATABLE :: IS_PROCESSED

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   DO NL = NLMIN, NLMAX

      L => SCARC(NM)%LEVEL(NL)
      G => L%UNSTRUCTURED
      MGM => L%MGM

      MGM%NWE = L%N_WALL_CELLS_EXT
      MGM%NWI = L%N_WALL_CELLS_INT

      MGM%NW1 = L%N_WALL_CELLS_EXT + 1
      MGM%NW2 = L%N_WALL_CELLS_EXT + L%N_WALL_CELLS_INT

      CALL SCARC_ALLOCATE_REAL1(MGM%H1, 1, G%NC, NSCARC_INIT_ZERO, 'H1')
      CALL SCARC_ALLOCATE_REAL1(MGM%H2, 1, G%NC, NSCARC_INIT_ZERO, 'H2')

      CALL SCARC_ALLOCATE_REAL1(MGM%US, MGM%NW1, MGM%NW2, NSCARC_INIT_ZERO, 'UMGM')
      CALL SCARC_ALLOCATE_REAL1(MGM%VS, MGM%NW1, MGM%NW2, NSCARC_INIT_ZERO, 'VMGM')
      CALL SCARC_ALLOCATE_REAL1(MGM%WS, MGM%NW1, MGM%NW2, NSCARC_INIT_ZERO, 'WMGM')

      CALL SCARC_ALLOCATE_INT1 (MGM%PERM , 1, G%NC, NSCARC_INIT_ZERO, 'PERMUTATION')

      CALL SCARC_ALLOCATE_REAL2(MGM%A , 1, G%NC, 1, G%NC, NSCARC_INIT_ZERO, 'A')
      CALL SCARC_ALLOCATE_REAL2(MGM%L , 1, G%NC, 1, G%NC, NSCARC_INIT_ZERO, 'L')
      CALL SCARC_ALLOCATE_REAL2(MGM%U , 1, G%NC, 1, G%NC, NSCARC_INIT_ZERO, 'U')

      CALL SCARC_ALLOCATE_REAL2(MGM%IA, 1, G%NC, 1, G%NC, NSCARC_INIT_ZERO, 'IA')
      CALL SCARC_ALLOCATE_REAL2(MGM%IL, 1, G%NC, 1, G%NC, NSCARC_INIT_ZERO, 'IL')
      CALL SCARC_ALLOCATE_REAL2(MGM%IU, 1, G%NC, 1, G%NC, NSCARC_INIT_ZERO, 'IU')

      CALL SCARC_ALLOCATE_REAL1(MGM%B, 1, G%NC, NSCARC_INIT_ZERO, 'B')
      CALL SCARC_ALLOCATE_REAL1(MGM%Y, 1, G%NC, NSCARC_INIT_ZERO, 'Y')
      CALL SCARC_ALLOCATE_REAL1(MGM%X, 1, G%NC, NSCARC_INIT_ZERO, 'X')

      CALL SCARC_ALLOCATE_LOG1(IS_PROCESSED, 1, G%NC, NSCARC_INIT_FALSE, 'IS_PROCESSED')

      A => SCARC_POINT_TO_CMATRIX(G, NSCARC_MATRIX_POISSON)

      IP = G%NC - MGM%NWI + 1
      DO IW = MGM%NW1, MGM%NW2
         IC = G%WALL(IW)%ICW
         IS_PROCESSED(IC) = .TRUE.
         MGM%PERM(IP) = IC
         IP = IP + 1
      ENDDO

      IP = 1
      DO IC = 1, G%NC
         IF (IS_PROCESSED(IC)) CYCLE
         IS_PROCESSED(IC) = .TRUE.
         MGM%PERM(IP) = IC
         IP = IP + 1
      ENDDO

      DEALLOCATE (IS_PROCESSED)

      CALL SCARC_SETUP_MGM_LU(NM, NL)
      CALL SCARC_METHOD_MGM_LU(NM, NL)

      CALL SCARC_SETUP_MGM_ILU(NM, NL)
      CALL SCARC_METHOD_MGM_ILU(NM, NL)

   ENDDO
ENDDO

END SUBROUTINE SCARC_SETUP_MGM


! ---------------------------------------------------------------------------------------------
! Setup LU-decomposition for MGM-method
! ---------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_MGM_LU(NM, NL)
USE SCARC_POINTERS, ONLY: G, MGM, A 
INTEGER, INTENT(IN) :: NM, NL
INTEGER :: I, J, K, IC, JC, ICOL, IP, IVERSION
REAL (EB) :: SCAL

G   => SCARC(NM)%LEVEL(NL)%UNSTRUCTURED
MGM => SCARC(NM)%LEVEL(NL)%MGM
A   => G%POISSONC

!
! Temporarily extract full matrix from compact storage technique - just for proof of concept
! Consider permutation in MGM%PERM
!
DO JC = 1, G%NC
   IC = MGM%PERM(JC)
   IC = JC

#ifdef WITH_SCARC_DEBUG
   WRITE(MSG%LU_DEBUG,*) 'PERMUTATION  --- JC =', JC, '-----------> IC=',IC
#endif

   DO IP = A%ROW(IC), A%ROW(IC+1)-1
      ICOL = A%COL(IP)

#ifdef WITH_SCARC_DEBUG
      WRITE(MSG%LU_DEBUG,*) 'PERMUTATION ICOL ', ICOL, A%VAL(IP), ' TO ', JC, ICOL
#endif

      MGM%A(JC,ICOL) = A%VAL(IP)
   ENDDO
ENDDO

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) '------- MGM%A - Copy (1:24)'
DO IC = 1, G%NC
   WRITE(MSG%LU_DEBUG,'(24F11.4)') (MGM%A(IC, JC), JC=1, 24)
ENDDO
#endif

IVERSION = 1

! Version 1
!IF (IVERSION == 1) THEN
DO I = 1, G%NC
   MGM%U(I,I) = 1.0_EB
ENDDO

DO J = 1, G%NC
   DO I = J, G%NC
      SCAL = 0.0_EB
      DO K = 1, J-1
         SCAL = SCAL + MGM%L(I,K) * MGM%U(K,J)
      ENDDO
      MGM%L(I,J) = MGM%A(I,J) - SCAL
   ENDDO
   DO I = J, G%NC
      SCAL = 0.0_EB
      DO K = 1, J-1
         SCAL = SCAL + MGM%L(J,K) * MGM%U(K,I)
      ENDDO
      MGM%U(J,I) = (MGM%A(J,I) - SCAL) / MGM%L(J,J)
   ENDDO
ENDDO

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) '=============================== MGM%A (1:24) '
DO IC = 1, G%NC
   WRITE(MSG%LU_DEBUG,'(24F11.4)') (MGM%A(IC, JC), JC=1, 24)
ENDDO
WRITE(MSG%LU_DEBUG,*) 'MGM%L'
DO IC = 1, G%NC
   WRITE(MSG%LU_DEBUG,'(24F11.4)') (MGM%L(IC, JC), JC=1, 24)
ENDDO
WRITE(MSG%LU_DEBUG,*) 'MGM%U'
DO IC = 1, G%NC
   WRITE(MSG%LU_DEBUG,'(24F11.4)') (MGM%U(IC, JC), JC=1, 24)
ENDDO
#endif

! Version 2
!ELSE IF (IVERSION == 2) THEN

DO I = 1, G%NC

   MGM%L(I,I) = 1.0_EB

   DO J = I, G%NC

      SCAL = 0.0_EB
      DO K = 1, I-1
         SCAL = SCAL + MGM%L(I,K) * MGM%U(K,J)
      ENDDO
      MGM%U(I,J) = MGM%A(I,J) - SCAL

      SCAL = 0.0_EB
      DO K = 1, I-1
         SCAL = SCAL + MGM%L(J,K) * MGM%U(K,I)
      ENDDO
      MGM%L(J,I) = (MGM%A(J,I) - SCAL)/MGM%U(I,I)

   ENDDO
ENDDO

#ifdef WITH_SCARC_MGM
WRITE(MSG%LU_DEBUG,*) '=============================== MGM%A (1:24) '
DO IC = 1, G%NC
   WRITE(MSG%LU_DEBUG,'(24F11.4)') (MGM%A(IC, JC), JC=1, 24)
ENDDO
WRITE(MSG%LU_DEBUG,*) 'MGM%L'
DO IC = 1, G%NC
   WRITE(MSG%LU_DEBUG,'(24F11.4)') (MGM%L(IC, JC), JC=1, 24)
ENDDO
WRITE(MSG%LU_DEBUG,*) 'MGM%U'
DO IC = 1, G%NC
   WRITE(MSG%LU_DEBUG,'(24F11.4)') (MGM%U(IC, JC), JC=1, 24)
ENDDO
#endif

! Version 3
!ELSE IF (IVERSION == 3) THEN

!   DO I = 1, G%NC
!      MGM%L(I,1) = MGM%A(I,1)
!      MGM%U(I,I) = 1.0_EB
!   ENDDO
!
!   DO J = 2, G%NC
!      MGM%U(1,J) = MGM%A(1,J)/MGM%L(1,1)
!   ENDDO
!
!   DO I = 2, G%NC
!
!      DO J = 2, I
!         SCAL = 0.0_EB
!         DO K = 1, J-1
!            SCAL = SCAL + MGM%L(I,K)*MGM%U(K,J)
!         ENDDO
!         L (I,J) = MGM%A(I,J) - SCAL
!      ENDDO
!
!      DO J = I+1, G%NC
!         SCAL = 0.0_EB
!         DO K = 1, I-1
!            SCAL = SCAL + MGM%L(I,K)*MGM%U(K,J)/MGM%L(I,I)
!         ENDDO
!         MGM%U(I,J) = MGM%A(I,J) - SCAL
!      ENDDO
!   ENDDO

!ENDIF

DO I = 1, G%NC
   DO J = 1, I-1
      MGM%A(I,J) = MGM%L(I,J)
   ENDDO
   DO J = I, G%NC
      MGM%A(I,J) = MGM%U(I,J)
   ENDDO
ENDDO

END SUBROUTINE SCARC_SETUP_MGM_LU


! ---------------------------------------------------------------------------------------------
! Setup ILU-decomposition for MGM-method
! ---------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_MGM_ILU(NM, NL)
USE SCARC_POINTERS, ONLY: G, MGM, AC
INTEGER, INTENT(IN) :: NM, NL
INTEGER :: I, J, K, IC, JC, ICOL, IP

G   => SCARC(NM)%LEVEL(NL)%UNSTRUCTURED
MGM => SCARC(NM)%LEVEL(NL)%MGM
AC  => G%POISSONC

!
! Temporarily extract full matrix from compact storage technique - just for proof of concept
! Consider permutation in MGM%PERM
!
DO JC = 1, G%NC
   !   IC = MGM%PERM(JC)
   IC = JC
   DO IP = A%ROW(IC), A%ROW(IC+1)-1
      ICOL = A%COL(IP)
      MGM%IA(JC,ICOL) = A%VAL(IP)
   ENDDO
ENDDO

DO I = 2, G%NC
   DO K = 1, I-1
      IF (MGM%IA(I,K) /= 0 .AND. MGM%IA(K,K) /= 0) THEN
         MGM%IA(I,K) = MGM%IA(I,K)/MGM%IA(K,K)
         DO J = K+1, G%NC
            IF (MGM%IA(I,J) /= 0) THEN
               MGM%IA(I,J) = MGM%IA(I,J) - MGM%IA(I,K)*MGM%IA(K,J)
            ENDIF
         ENDDO
      ENDIF
   ENDDO
ENDDO

END SUBROUTINE SCARC_SETUP_MGM_ILU


! ---------------------------------------------------------------------------------------------
! Setup LU-decomposition for MGM-method
! ---------------------------------------------------------------------------------------------
SUBROUTINE SCARC_METHOD_MGM_LU(NM, NL)
USE SCARC_POINTERS, ONLY: G, MGM, AC
INTEGER, INTENT(IN) :: NM, NL
INTEGER :: I, J, K, IC, IW, N

G   => SCARC(NM)%LEVEL(NL)%UNSTRUCTURED
MGM => SCARC(NM)%LEVEL(NL)%MGM
AC  => G%POISSONC

N = G%NC

DO IW = MGM%NW1, MGM%NW2
   IC = G%WALL(IW)%ICW
   MGM%B(IC) = 1.0_EB
ENDDO

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) '=============================== A'
DO I = 1, N
   WRITE(MSG%LU_DEBUG,'(24F11.4)') (MGM%A(I, J), J=1, 24)
ENDDO
WRITE(MSG%LU_DEBUG,*) '=============================== L'
DO I = 1, N
   WRITE(MSG%LU_DEBUG,'(24F11.4)') (MGM%L(I, J), J=1, 24)
ENDDO
WRITE(MSG%LU_DEBUG,*) '=============================== U'
DO I = 1, N
   WRITE(MSG%LU_DEBUG,'(24F11.4)') (MGM%U(I, J), J=1, 24)
ENDDO
WRITE(MSG%LU_DEBUG,*) '=============================== B'
WRITE(MSG%LU_DEBUG,'(16F11.4)') (MGM%B(I), I=1, G%NC)
#endif

DO J = 1, N
   MGM%Y(J) = MGM%B(J)
   DO K = 1, J-1
      MGM%Y(J) = MGM%Y(J) - MGM%A(J,K)*MGM%Y(K)
   ENDDO
ENDDO

DO J = N, 1, -1
   MGM%X(J) = MGM%Y(J)
   DO K = J+1, N
      MGM%X(J) = MGM%X(J) - MGM%A(J,K)*MGM%X(K)
   ENDDO
   MGM%X(J) = MGM%X(J)/MGM%A(J,J)
ENDDO

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) '=============================== Y'
WRITE(MSG%LU_DEBUG,'(16F11.4)') (MGM%Y(I), I=1, G%NC)
WRITE(MSG%LU_DEBUG,*) '=============================== X'
WRITE(MSG%LU_DEBUG,'(16F11.4)') (MGM%X(I), I=1, G%NC)
#endif

END SUBROUTINE SCARC_METHOD_MGM_LU


! ---------------------------------------------------------------------------------------------
! Setup ILU-decomposition for MGM-method
! ---------------------------------------------------------------------------------------------
SUBROUTINE SCARC_METHOD_MGM_ILU(NM, NL)
USE SCARC_POINTERS, ONLY: G, MGM, AC
INTEGER, INTENT(IN) :: NM, NL
INTEGER :: I, J, K, IC, IW, N
REAL (EB), DIMENSION(:,:), POINTER :: IA, IL, IU
REAL (EB), DIMENSION(:),   POINTER :: B, X, Y

G   => SCARC(NM)%LEVEL(NL)%UNSTRUCTURED
MGM => SCARC(NM)%LEVEL(NL)%MGM
AC  => G%POISSONC

IA => MGM%IA
IL => MGM%IL
IU => MGM%IU

B => MGM%B
X => MGM%X
Y => MGM%Y

N = G%NC

DO IW = MGM%NW1, MGM%NW2
   IC = G%WALL(IW)%ICW
   MGM%B(IC) = 1.0_EB
ENDDO

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) '=============================== A'
DO I = 1, N
   WRITE(MSG%LU_DEBUG,'(24F11.4)') (MGM%IA(I, J), J=1, 24)
ENDDO
WRITE(MSG%LU_DEBUG,*) '=============================== B'
WRITE(MSG%LU_DEBUG,'(16F11.4)') (MGM%B(I), I=1, G%NC)
#endif

DO J = 1, N
   MGM%Y(J) = MGM%B(J)
   DO K = 1, J-1
      MGM%Y(J) = MGM%Y(J) - MGM%IA(J,K)*MGM%Y(K)
   ENDDO
ENDDO

DO J = N, 1, -1
   MGM%X(J) = MGM%Y(J)
   DO K = J+1, N
      MGM%X(J) = MGM%X(J) - MGM%IA(J,K)*MGM%X(K)
   ENDDO
   MGM%X(J) = MGM%X(J)/MGM%IA(J,J)
ENDDO

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) '=============================== Y'
WRITE(MSG%LU_DEBUG,'(16F11.4)') (MGM%Y(I), I=1, G%NC)
WRITE(MSG%LU_DEBUG,*) '=============================== X'
WRITE(MSG%LU_DEBUG,'(16F11.4)') (MGM%X(I), I=1, G%NC)
#endif

END SUBROUTINE SCARC_METHOD_MGM_ILU

! ====================================================================================================
! End SCARC_WITH_MGM - Part
! ====================================================================================================
#endif



#ifdef WITH_SCARC_DEBUG
! ------------------------------------------------------------------------------------------------------
! Temporarily preset B1 case with specified aggregation pattern
! ------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_PRESET_B1_CASE
USE SCARC_POINTERS, ONLY: GF, GC

CALL SCARC_POINT_TO_MULTIGRID(1, 1, 2)

GF%N_COARSE = 6
GF%N_FINE = 24

GF%ZONES_LOCAL = 0
GF%ZONES_LOCAL(1)  = 1
GF%ZONES_LOCAL(2)  = 1
GF%ZONES_LOCAL(3)  = 2
GF%ZONES_LOCAL(4)  = 4
GF%ZONES_LOCAL(5)  = 4
GF%ZONES_LOCAL(6)  = 5
GF%ZONES_LOCAL(7)  = 1
GF%ZONES_LOCAL(8)  = 2
GF%ZONES_LOCAL(9)  = 2
GF%ZONES_LOCAL(10) = 4
GF%ZONES_LOCAL(11) = 5
GF%ZONES_LOCAL(12) = 5
GF%ZONES_LOCAL(13) = 3
GF%ZONES_LOCAL(14) = 2
GF%ZONES_LOCAL(15) = 2
GF%ZONES_LOCAL(16) = 6
GF%ZONES_LOCAL(17) = 5
GF%ZONES_LOCAL(18) = 5
GF%ZONES_LOCAL(19) = 3
GF%ZONES_LOCAL(20) = 3
GF%ZONES_LOCAL(21) = 2
GF%ZONES_LOCAL(22) = 6
GF%ZONES_LOCAL(23) = 6
GF%ZONES_LOCAL(24) = 5
GF%ZONES_GLOBAL = GF%ZONES_LOCAL

GF%CPOINTS(1) =  1
GF%CPOINTS(2) =  9
GF%CPOINTS(3) = 19
GF%CPOINTS(4) =  4
GF%CPOINTS(5) = 12
GF%CPOINTS(6) = 22

GF%N_COARSE = 6
GF%N_ZONES = 6

GC%N_FINE = 6
GC%NC_LOCAL(1) = 6
GC%NC_GLOBAL = 6

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'PRESETTING B1: ZONES_GLOBAL:'
WRITE(MSG%LU_DEBUG,CFORM_INT) GF%ZONES_GLOBAL(1:GF%NC)
WRITE(MSG%LU_DEBUG,*) 'PRESETTING B1: CPOINTS:'
WRITE(MSG%LU_DEBUG,CFORM_INT) GF%CPOINTS
WRITE(MSG%LU_DEBUG,*) 'PRESETTING B1: GF%N_ZONES  =',GF%N_ZONES
WRITE(MSG%LU_DEBUG,*) 'PRESETTING B1: GF%N_FINE  =',GF%N_FINE
WRITE(MSG%LU_DEBUG,*) 'PRESETTING B1: GF%N_COARSE=',GF%N_COARSE
WRITE(MSG%LU_DEBUG,*) 'PRESETTING B1: GC%N_FINE=',GC%N_FINE
WRITE(MSG%LU_DEBUG,*) 'PRESETTING B1: GC%NC_LOCAL(1)=',GC%NC_LOCAL(1)
WRITE(MSG%LU_DEBUG,*) 'PRESETTING B1: GC%NC_GLOBAL=',GC%NC_GLOBAL
#endif

END SUBROUTINE SCARC_PRESET_B1_CASE
 

! ------------------------------------------------------------------------------------------------------
! Temporarily preset B1 case with specified aggregation pattern
! ------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_PRESET_B14_CASE
USE SCARC_POINTERS, ONLY: GF, GC

CALL SCARC_POINT_TO_MULTIGRID(1, 1, 2)

GF%N_COARSE = 8
GF%N_FINE = 36

GF%ZONES_LOCAL = 0
GF%ZONES_LOCAL(1)  = 1
GF%ZONES_LOCAL(2)  = 1
GF%ZONES_LOCAL(3)  = 2
GF%ZONES_LOCAL(4)  = 3
GF%ZONES_LOCAL(5)  = 3
GF%ZONES_LOCAL(6)  = 4
GF%ZONES_LOCAL(7)  = 1
GF%ZONES_LOCAL(8)  = 2
GF%ZONES_LOCAL(9)  = 2
GF%ZONES_LOCAL(10) = 3
GF%ZONES_LOCAL(11) = 4
GF%ZONES_LOCAL(12) = 4
GF%ZONES_LOCAL(13) = 1
GF%ZONES_LOCAL(14) = 2
GF%ZONES_LOCAL(15) = 2
GF%ZONES_LOCAL(16) = 3
GF%ZONES_LOCAL(17) = 4
GF%ZONES_LOCAL(18) = 4
GF%ZONES_LOCAL(19) = 5
GF%ZONES_LOCAL(20) = 5
GF%ZONES_LOCAL(21) = 6
GF%ZONES_LOCAL(22) = 7
GF%ZONES_LOCAL(23) = 7
GF%ZONES_LOCAL(24) = 8
GF%ZONES_LOCAL(25) = 5
GF%ZONES_LOCAL(26) = 6
GF%ZONES_LOCAL(27) = 6
GF%ZONES_LOCAL(28) = 7
GF%ZONES_LOCAL(29) = 8
GF%ZONES_LOCAL(30) = 8
GF%ZONES_LOCAL(31) = 5
GF%ZONES_LOCAL(32) = 6
GF%ZONES_LOCAL(33) = 6
GF%ZONES_LOCAL(34) = 7
GF%ZONES_LOCAL(35) = 8
GF%ZONES_LOCAL(36) = 8
GF%ZONES_GLOBAL = GF%ZONES_LOCAL

GF%CPOINTS(1) =  1
GF%CPOINTS(2) =  9
GF%CPOINTS(3) =  4
GF%CPOINTS(4) = 12
GF%CPOINTS(5) = 19
GF%CPOINTS(6) = 27
GF%CPOINTS(7) = 22
GF%CPOINTS(8) = 30

GF%N_COARSE = 8
GF%N_ZONES = 8

GC%N_FINE = 8
GC%NC_LOCAL(1) = 8
GC%NC_GLOBAL = 8

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'PRESETTING B14_CASE: ZONES_GLOBAL:'
WRITE(MSG%LU_DEBUG,CFORM_INT) GF%ZONES_GLOBAL(1:GF%NC)
WRITE(MSG%LU_DEBUG,*) 'PRESETTING B14_CASE: CPOINTS:'
WRITE(MSG%LU_DEBUG,CFORM_INT) GF%CPOINTS
WRITE(MSG%LU_DEBUG,*) 'PRESETTING B14_CASE: GF%N_ZONES  =',GF%N_ZONES
WRITE(MSG%LU_DEBUG,*) 'PRESETTING B14_CASE: GF%N_FINE  =',GF%N_FINE
WRITE(MSG%LU_DEBUG,*) 'PRESETTING B14_CASE: GF%N_COARSE=',GF%N_COARSE
WRITE(MSG%LU_DEBUG,*) 'PRESETTING B14_CASE: GC%N_FINE=',GC%N_FINE
WRITE(MSG%LU_DEBUG,*) 'PRESETTING B14_CASE: GC%NC_LOCAL(1)=',GC%NC_LOCAL(1)
WRITE(MSG%LU_DEBUG,*) 'PRESETTING B14_CASE: GC%NC_GLOBAL=',GC%NC_GLOBAL
#endif

END SUBROUTINE SCARC_PRESET_B14_CASE
 
! ------------------------------------------------------------------------------------------------------
! Temporarily preset B1 case with specified aggregation pattern
! ------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_PRESET_b14big_CASE
USE SCARC_POINTERS, ONLY: GF, GC

CALL SCARC_POINT_TO_MULTIGRID(1, 1, 2)

GF%N_COARSE = 24
GF%N_FINE = 108

GF%ZONES_LOCAL = 0
GF%ZONES_LOCAL(1)  = 1
GF%ZONES_LOCAL(2)  = 1
GF%ZONES_LOCAL(3)  = 2
GF%ZONES_LOCAL(4)  = 2
GF%ZONES_LOCAL(5)  = 2
GF%ZONES_LOCAL(6)  = 3
GF%ZONES_LOCAL(7)  = 3
GF%ZONES_LOCAL(8)  = 9
GF%ZONES_LOCAL(9)  = 9
GF%ZONES_LOCAL(10) = 10
GF%ZONES_LOCAL(11) = 10
GF%ZONES_LOCAL(12) = 10

GF%ZONES_LOCAL(13) = 1
GF%ZONES_LOCAL(14) = 4
GF%ZONES_LOCAL(15) = 2
GF%ZONES_LOCAL(16) = 2
GF%ZONES_LOCAL(17) = 5
GF%ZONES_LOCAL(18) = 3
GF%ZONES_LOCAL(19) = 3
GF%ZONES_LOCAL(20) = 9
GF%ZONES_LOCAL(21) = 11
GF%ZONES_LOCAL(22) = 10
GF%ZONES_LOCAL(23) = 10
GF%ZONES_LOCAL(24) = 12

GF%ZONES_LOCAL(25) = 4
GF%ZONES_LOCAL(26) = 4
GF%ZONES_LOCAL(27) = 4
GF%ZONES_LOCAL(28) = 5
GF%ZONES_LOCAL(29) = 5
GF%ZONES_LOCAL(30) = 5
GF%ZONES_LOCAL(31) = 6
GF%ZONES_LOCAL(32) = 11
GF%ZONES_LOCAL(33) = 11
GF%ZONES_LOCAL(34) = 11
GF%ZONES_LOCAL(35) = 12
GF%ZONES_LOCAL(36) = 12

GF%ZONES_LOCAL(37) = 7
GF%ZONES_LOCAL(38) = 4
GF%ZONES_LOCAL(39) = 4
GF%ZONES_LOCAL(40) = 8
GF%ZONES_LOCAL(41) = 5
GF%ZONES_LOCAL(42) = 6
GF%ZONES_LOCAL(43) = 6
GF%ZONES_LOCAL(44) = 13
GF%ZONES_LOCAL(45) = 11
GF%ZONES_LOCAL(46) = 11
GF%ZONES_LOCAL(47) = 14
GF%ZONES_LOCAL(48) = 12

GF%ZONES_LOCAL(49) = 7
GF%ZONES_LOCAL(50) = 7
GF%ZONES_LOCAL(51) = 8
GF%ZONES_LOCAL(52) = 8
GF%ZONES_LOCAL(53) = 8
GF%ZONES_LOCAL(54) = 6
GF%ZONES_LOCAL(55) = 6
GF%ZONES_LOCAL(56) = 13
GF%ZONES_LOCAL(57) = 13
GF%ZONES_LOCAL(58) = 14
GF%ZONES_LOCAL(59) = 14
GF%ZONES_LOCAL(60) = 14

GF%ZONES_LOCAL(61) = 15
GF%ZONES_LOCAL(62) = 15
GF%ZONES_LOCAL(63) = 16
GF%ZONES_LOCAL(64) = 16
GF%ZONES_LOCAL(65) = 16
GF%ZONES_LOCAL(66) = 17
GF%ZONES_LOCAL(67) = 17
GF%ZONES_LOCAL(68) = 21
GF%ZONES_LOCAL(69) = 21
GF%ZONES_LOCAL(70) = 22
GF%ZONES_LOCAL(71) = 22
GF%ZONES_LOCAL(72) = 22

GF%ZONES_LOCAL(73) = 15
GF%ZONES_LOCAL(74) = 18
GF%ZONES_LOCAL(75) = 16
GF%ZONES_LOCAL(76) = 16
GF%ZONES_LOCAL(77) = 19
GF%ZONES_LOCAL(78) = 17
GF%ZONES_LOCAL(79) = 17
GF%ZONES_LOCAL(80) = 21
GF%ZONES_LOCAL(81) = 23
GF%ZONES_LOCAL(82) = 22
GF%ZONES_LOCAL(83) = 22
GF%ZONES_LOCAL(84) = 24

GF%ZONES_LOCAL(85) = 18
GF%ZONES_LOCAL(86) = 18
GF%ZONES_LOCAL(87) = 18
GF%ZONES_LOCAL(88) = 19
GF%ZONES_LOCAL(89) = 19
GF%ZONES_LOCAL(90) = 19
GF%ZONES_LOCAL(91) = 20
GF%ZONES_LOCAL(92) = 23
GF%ZONES_LOCAL(93) = 23
GF%ZONES_LOCAL(94) = 23
GF%ZONES_LOCAL(95) = 24
GF%ZONES_LOCAL(96) = 24

GF%ZONES_LOCAL(97) = 18
GF%ZONES_LOCAL(98) = 18
GF%ZONES_LOCAL(99) = 18
GF%ZONES_LOCAL(100) = 19
GF%ZONES_LOCAL(101) = 19
GF%ZONES_LOCAL(102) = 20
GF%ZONES_LOCAL(103) = 20
GF%ZONES_LOCAL(104) = 23
GF%ZONES_LOCAL(105) = 23
GF%ZONES_LOCAL(106) = 23
GF%ZONES_LOCAL(107) = 24
GF%ZONES_LOCAL(108) = 24

GF%ZONES_GLOBAL = GF%ZONES_LOCAL

GF%CPOINTS(1) =  1
GF%CPOINTS(2) =  4
GF%CPOINTS(3) =  7
GF%CPOINTS(4) = 26
GF%CPOINTS(5) = 29
GF%CPOINTS(6) = 43
GF%CPOINTS(7) = 49
GF%CPOINTS(8) = 52
GF%CPOINTS(9) = 8
GF%CPOINTS(10) = 11
GF%CPOINTS(11) = 33
GF%CPOINTS(12) = 36
GF%CPOINTS(13) = 56
GF%CPOINTS(14) = 59
GF%CPOINTS(15) = 61
GF%CPOINTS(16) = 64
GF%CPOINTS(17) = 67
GF%CPOINTS(18) = 86
GF%CPOINTS(19) = 89
GF%CPOINTS(20) = 103
GF%CPOINTS(21) = 68
GF%CPOINTS(22) = 71
GF%CPOINTS(23) = 93
GF%CPOINTS(24) = 96

GF%N_COARSE = 24
GF%N_ZONES = 24

GC%N_FINE = 24
GC%NC_LOCAL(1) = 24
GC%NC_GLOBAL = 24

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'PRESETTING B14_CASE: ZONES_GLOBAL:'
WRITE(MSG%LU_DEBUG,CFORM_INT) GF%ZONES_GLOBAL(1:GF%NC)
WRITE(MSG%LU_DEBUG,*) 'PRESETTING B14_CASE: CPOINTS:'
WRITE(MSG%LU_DEBUG,CFORM_INT) GF%CPOINTS
WRITE(MSG%LU_DEBUG,*) 'PRESETTING B14_CASE: GF%N_ZONES  =',GF%N_ZONES
WRITE(MSG%LU_DEBUG,*) 'PRESETTING B14_CASE: GF%N_FINE  =',GF%N_FINE
WRITE(MSG%LU_DEBUG,*) 'PRESETTING B14_CASE: GF%N_COARSE=',GF%N_COARSE
WRITE(MSG%LU_DEBUG,*) 'PRESETTING B14_CASE: GC%N_FINE=',GC%N_FINE
WRITE(MSG%LU_DEBUG,*) 'PRESETTING B14_CASE: GC%NC_LOCAL(1)=',GC%NC_LOCAL(1)
WRITE(MSG%LU_DEBUG,*) 'PRESETTING B14_CASE: GC%NC_GLOBAL=',GC%NC_GLOBAL
#endif

END SUBROUTINE SCARC_PRESET_b14big_CASE
 
! ================================================================================================
! Start  WITH_SCARC_DEBUG  - Part
! Collection of routines which print out different quantities or allow to preset them
! ================================================================================================
! ------------------------------------------------------------------------------------------------
! Preset right hand side in such a way that exact solution is known
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_PRESET_EXACT (NE, NL)
USE SCARC_POINTERS, ONLY: M, L, G, VC, XMID, ZMID
USE SCARC_ITERATION_ENVIRONMENT
INTEGER, INTENT(IN) :: NE, NL
INTEGER :: IC, NM, I, K

IF (ITE_TOTAL == 0) WRITE(*,*) 'TODO: PRESET_EXACT is active !!!'

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G
   VC => SCARC_POINT_TO_VECTOR (NM, NL, NE)

   DO K = 1, L%NZ
      DO I = 1, L%NX
         IF (IS_UNSTRUCTURED .AND. L%IS_SOLID(I,1,K)) CYCLE
         IC = G%CELL_NUMBER(I,1,K)
         IF (NL == NLEVEL_MIN) THEN
            XMID => M%XC
            ZMID => M%ZC
         ELSE
            XMID => L%XMID
            ZMID => L%ZMID
         ENDIF
         VC(IC) = EXACT(XMID(I),ZMID(K))
         !WRITE(MSG%LU_DEBUG,'(A,i3,a,e10.2,a,e10.2,a,e12.4)') 'IC=',IC,':X=',XMID(i),':Z=',ZMID(k),': RHS=',VC(IC)
      ENDDO
   ENDDO
ENDDO

END SUBROUTINE SCARC_PRESET_EXACT


! ------------------------------------------------------------------------------------------------
! Preset vector
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_PRESET_VECTOR (NV, NL)
USE SCARC_POINTERS, ONLY: M, G, VC, XMID, ZMID
INTEGER, INTENT(IN) :: NV, NL
INTEGER :: IC, NM, I, K

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G
   VC => SCARC_POINT_TO_VECTOR (NM, NL, NV)

   DO K = 1, L%NZ
      DO I = 1, L%NX
         IF (IS_UNSTRUCTURED .AND. L%IS_SOLID(I,1,K)) CYCLE
         IC = G%CELL_NUMBER(I,1,K)
         IF (NL == NLEVEL_MIN) THEN
            XMID => M%XC
            ZMID => M%ZC
         ELSE
            XMID => L%XMID
            ZMID => L%ZMID
         ENDIF
         VC(IC) = XMID(I)
         !WRITE(MSG%LU_DEBUG,*) 'IC=',IC,':X=',XMID,':Z=',ZMID,': RHS=',VC(IC)
      ENDDO
   ENDDO
ENDDO

END SUBROUTINE SCARC_PRESET_VECTOR


! ------------------------------------------------------------------------------------------------
! Preset right hand side in such a way that exact solution is known
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_PRESET_RHS (NV, NL)
USE SCARC_POINTERS, ONLY: M, L, G, VC
INTEGER, INTENT(IN) :: NV, NL
INTEGER :: IC, NM, I, K
REAL (EB) :: X, Z

IF (NL > NLEVEL_MIN) WRITE(*,*) 'Wrong level for presetting RHS '

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   M => MESHES(NM)
   M%BXS = 0.0_EB
   M%BXF = 0.0_EB
   M%BZS = 0.0_EB
   M%BZF = 0.0_EB

   CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G
   VC => SCARC_POINT_TO_VECTOR (NM, NL, NV)

   DO K = 1, L%NZ
      DO I = 1, L%NX
         IF (IS_UNSTRUCTURED .AND. L%IS_SOLID(I,1,K)) CYCLE
         IC = G%CELL_NUMBER(I,1,K)
         X  = M%XC(I)
         Z  = M%ZC(K)
         !WRITE(MSG%LU_DEBUG,'(A,i3,a,e10.2,a,e10.2,a,e12.4)') 'IC=',IC,':X=',X,':Z=',Z,': RHS=',VC(IC)
         VC(IC) = RHS(X,Z)
      ENDDO
   ENDDO
ENDDO

END SUBROUTINE SCARC_PRESET_RHS


! ------------------------------------------------------------------------------------------------
! Set exact solution
! ------------------------------------------------------------------------------------------------
DOUBLE PRECISION FUNCTION EXACT(X,Z)
REAL (EB), INTENT(IN) :: X, Z
!EXACT = (X**2 - X**4) * (Z**4 - Z**2)                                    ! FUNCTION 1
!EXACT = (X**2 - 1) * (Z**2 - 1)                                         ! FUNCTION 2
!EXACT =  625.0_EB/16.0_EB * X * (0.8_EB - X) * Z * (0.8_EB - Z)        ! FUNCTION 3
EXACT = - X * (0.8_EB - X) * Z * (0.8_EB - Z)        ! FUNCTION 3
END FUNCTION EXACT


! ------------------------------------------------------------------------------------------------
! Set right hand side
! ------------------------------------------------------------------------------------------------
DOUBLE PRECISION FUNCTION RHS(X,Z)
REAL (EB), INTENT(IN) :: X, Z
!RHS = 2.0_EB*((1.0_EB - 6.0_EB*X**2)*Z**2*(1.0_EB-Z**2)+(1.0_EB-6.0_EB*Z**2)*X**2*(1.0_EB-X**2))
!RHS = -X**2 - Z**2 +2
!RHS = 625.0_EB/8.0_EB * (X * (0.8_EB - X) + Z * (0.8_EB - Z))
RHS = 2.0_EB * (X * (0.8_EB - X) + Z * (0.8_EB - Z))
END FUNCTION RHS


! ------------------------------------------------------------------------------------------------
! Save dump of vector in dump-directory
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_DUMP_QUANTITY (NV, CNAME, ISM, NS, NL)
USE SCARC_POINTERS, ONLY: G, VC
USE SCARC_ITERATION_ENVIRONMENT
INTEGER, INTENT(IN) :: NV, NS, NL, ISM
CHARACTER(*), INTENT(IN) :: CNAME
CHARACTER(80) :: FN_DUMP, CDIR
INTEGER :: LU_DUMP, IC, NM

IF (TRIM(CNAME) == 'RESIDUAL') THEN
   CDIR = 'res'
ELSE IF (TRIM(CNAME) == 'ERROR') THEN
   CDIR = 'err'
ELSE IF (TRIM(CNAME) == 'RHS') THEN
   CDIR = 'rhs'
ELSE IF (TRIM(CNAME) == 'EXACT') THEN
   CDIR = 'exa'
ELSE IF (TRIM(CNAME) == 'DISCRET') THEN
   CDIR = 'dis'
ENDIF

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   VC => SCARC_POINT_TO_VECTOR (NM, NL, NV)
   CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G
   if (ISM == 0) THEN
      WRITE (FN_DUMP, '(A,A3,A,i3.3,A,I3.3,A,I3.3,A,I3.3,A,I1,A,A)') &
         'dump/',CDIR,'_t',ITE_TOTAL,'_ALL',ITE_SMOOTH,'_mg',ITE_MG,'_cg',ITE_CG,'_level',NL,'_',&
         TRIM(STACK(NS)%SOLVER%CNAME)
   else if (ISM == NSCARC_CYCLING_PRESMOOTH) THEN
      WRITE (FN_DUMP, '(A,A3,A,i3.3,A,I3.3,A,I3.3,A,I3.3,A,I1,A,A)') &
         'dump/',CDIR,'_t',ITE_TOTAL,'_PRE',ITE_SMOOTH,'_mg',ITE_MG,'_cg',ITE_CG,'_level',NL,'_',&
         TRIM(STACK(NS)%SOLVER%CNAME)
   ELSE IF (ISM == NSCARC_CYCLING_POSTSMOOTH) THEN
      WRITE (FN_DUMP, '(A,A3,A,i3.3,A,I3.3,A,I3.3,A,I3.3,A,I1,A,A)') &
         'dump/',CDIR,'_t',ITE_TOTAL,'_POST',ITE_SMOOTH,'_mg',ITE_MG,'_cg',ITE_CG,'_level',NL,'_',&
         TRIM(STACK(NS)%SOLVER%CNAME)
   ENDIF
   LU_DUMP = GET_FILE_NUMBER()
   OPEN (LU_DUMP, FILE=FN_DUMP)
   DO IC = 1, G%NC
      !WRITE(LU_DUMP,'(F25.16)') VC(IC)
      WRITE(LU_DUMP,*) VC(IC)
   ENDDO
   CLOSE(LU_DUMP)
ENDDO

FN_DUMP = CNAME
END SUBROUTINE SCARC_DUMP_QUANTITY

! ------------------------------------------------------------------------------------------------------
! Debug INT1 array
! ------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_DEBUG_INT1(ARR, I1, I2, CNAME, CTEXT)
INTEGER, DIMENSION(:), INTENT(IN) :: ARR
INTEGER, INTENT(IN) :: I1, I2
CHARACTER(*), INTENT(IN) :: CNAME, CTEXT
INTEGER :: IC
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) '============ DEBUGGING INT1 ARRAY ', CNAME, ' AT ', TRIM(CTEXT)
WRITE(MSG%LU_DEBUG,'(8I6)') (ARR(IC), IC=I1, I2)
#endif
END SUBROUTINE SCARC_DEBUG_INT1


! ------------------------------------------------------------------------------------------------------
! Debug REAL1 array
! ------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_DEBUG_REAL1(ARR, I1, I2, CNAME, CTEXT)
REAL(EB), DIMENSION(:), INTENT(IN) :: ARR
CHARACTER(*), INTENT(IN) :: CNAME, CTEXT
INTEGER, INTENT(IN) :: I1, I2
INTEGER :: IC
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) '============ DEBUGGING REAL1 ARRAY ', CNAME, ' AT ', TRIM(CTEXT)
WRITE(MSG%LU_DEBUG,'(8E12.4)') (ARR(IC), IC=I1, I2)
#endif
END SUBROUTINE SCARC_DEBUG_REAL1

! ------------------------------------------------------------------------------------------------------
! Debug zones information
! ------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_DEBUG_ZONES(G, IC, ITYPE, CTEXT)
TYPE (SCARC_GRID_TYPE), POINTER, INTENT(IN) :: G
INTEGER, INTENT(IN) :: IC, ITYPE
CHARACTER(*), INTENT(IN) :: CTEXT
WRITE(MSG%LU_DEBUG,*) '================= DEBUG_ZONES at ', CTEXT
IF (IC /= -1) WRITE(MSG%LU_DEBUG,*) ' IC = ', IC
IF (ITYPE == 1) THEN
   WRITE(MSG%LU_DEBUG,*) '-------------- ZONES_LOCAL: internal:'
   WRITE(MSG%LU_DEBUG,'(16I7)') G%ZONES_LOCAL(1:G%NC)
   WRITE(MSG%LU_DEBUG,*) '-------------- ZONES_LOCAL: overlap:'
   WRITE(MSG%LU_DEBUG,'(16I7)') G%ZONES_LOCAL(G%NC+1: G%NCE2)
ELSE
   WRITE(MSG%LU_DEBUG,*) '-------------- ZONES_GLOBAL: internal:'
   WRITE(MSG%LU_DEBUG,'(16I7)') G%ZONES_GLOBAL(1:G%NC)
   WRITE(MSG%LU_DEBUG,*) '-------------- ZONES_GLOBAL: overlap:'
   WRITE(MSG%LU_DEBUG,'(16I7)') G%ZONES_GLOBAL(G%NC+1: G%NCE2)
ENDIF
WRITE(MSG%LU_DEBUG,*) '-------------- CPOINTS'
WRITE(MSG%LU_DEBUG,'(16I4)') G%CPOINTS
END SUBROUTINE SCARC_DEBUG_ZONES

! ------------------------------------------------------------------------------------------------------
! Debug compact matrix
! ------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_DEBUG_CMATRIX(A, CNAME, CTEXT)
CHARACTER(*), INTENT(IN) :: CNAME, CTEXT
TYPE (SCARC_MATRIX_COMPACT_TYPE), INTENT(INOUT) :: A      
INTEGER :: IC, ICOL
CHARACTER(40) :: CFORM

WRITE(MSG%LU_DEBUG,*)
WRITE(MSG%LU_DEBUG,*) '============ START DEBUGGING MATRIX ', CNAME, ' AT ', TRIM(CTEXT)
WRITE(MSG%LU_DEBUG,*) 'INTERNAL NAME OF MATRIX :', A%CNAME
WRITE(MSG%LU_DEBUG,*) 'REQUESTED SIZES N_ROW, N_VAL:', A%N_ROW, A%N_VAL
WRITE(MSG%LU_DEBUG,*) 'ALLOCATED SIZES N_ROW, N_VAL:', SIZE(A%ROW), SIZE(A%VAL)

#ifdef WITH_SCARC_VERBOSE
IF (A%NTYPE /= NSCARC_MATRIX_MINIMAL) THEN
   IF (A%N_ROW /= SIZE(A%ROW)) WRITE(MSG%LU_VERBOSE,*) &
         '!!! CAUTION : different sizes for ',TRIM(CNAME),' N_ROW=',A%N_ROW,' SIZE(A%ROW)=',SIZE(A%ROW), ' AT', CTEXT
   IF (A%N_VAL /= SIZE(A%VAL)) WRITE(MSG%LU_VERBOSE,*) &
         '!!! CAUTION : different sizes for ',TRIM(CNAME),' N_VAL=',A%N_VAL,' SIZE(A%VAL)=',SIZE(A%VAL), ' AT', CTEXT
ENDIF
#endif

IF (A%NTYPE /= NSCARC_MATRIX_MINIMAL) THEN
   IF (A%N_ROW /= SIZE(A%ROW)) WRITE(MSG%LU_DEBUG,*) &
         '!!! CAUTION : different sizes for ',TRIM(CNAME),' N_ROW=',A%N_ROW,' SIZE(A%ROW)=',SIZE(A%ROW), ' AT', CTEXT
   IF (A%N_VAL /= SIZE(A%VAL)) WRITE(MSG%LU_DEBUG,*) &
         '!!! CAUTION : different sizes for ',TRIM(CNAME),' N_VAL=',A%N_VAL,' SIZE(A%VAL)=',SIZE(A%VAL), ' AT', CTEXT
ENDIF

WRITE(MSG%LU_DEBUG,*)
WRITE(MSG%LU_DEBUG,*) "------------->", TRIM(CNAME),'%ROW:'
WRITE(MSG%LU_DEBUG,'(8I6)') (A%ROW(IC), IC=1, A%N_ROW)
WRITE(MSG%LU_DEBUG,*) "------------->", TRIM(CNAME),'%COL:'
DO IC = 1, A%N_ROW-1
   IF (A%ROW(IC) == 0) CYCLE
   IF (A%ROW(IC+1)-A%ROW(IC) < 10) THEN
      CFORM = "(I4,A,10I6)"
   ELSE IF (A%ROW(IC+1)-A%ROW(IC) < 20) THEN
      CFORM = "(I4,A,20I6)"
   ELSE IF (A%ROW(IC+1)-A%ROW(IC)  < 30) THEN
      CFORM = "(I4,A,30I6)"
   ELSE
      CFORM = "(I4,A,40I6)"
   ENDIF
   WRITE(MSG%LU_DEBUG,CFORM) IC,':', (A%COL(ICOL), ICOL=A%ROW(IC), A%ROW(IC+1)-1)
ENDDO
IF (ALLOCATED(A%COLG)) THEN
   WRITE(MSG%LU_DEBUG,*) "------------->", TRIM(CNAME),'%COLG:'
   DO IC = 1, A%N_ROW-1
      IF (A%ROW(IC) == 0) CYCLE
      WRITE (CFORM, '(A,I3,A)' ) "(I4,A,", A%ROW(IC+1)-A%ROW(IC),"I5)"
      IF (A%ROW(IC+1)-A%ROW(IC) < 10) THEN
         CFORM = "(I4,A,10I6)"
      ELSE IF (A%ROW(IC+1)-A%ROW(IC) < 20) THEN
         CFORM = "(I4,A,20I6)"
      ELSE IF (A%ROW(IC+1)-A%ROW(IC) < 30) THEN
         CFORM = "(I4,A,30I6)"
      ELSE
         CFORM = "(I4,A,40I6)"
      ENDIF
      WRITE(MSG%LU_DEBUG,CFORM) IC,':', (A%COLG(ICOL), ICOL=A%ROW(IC), A%ROW(IC+1)-1)
   ENDDO
ENDIF
IF (ALLOCATED(A%VAL)) THEN
WRITE(MSG%LU_DEBUG,*) "------------->", TRIM(CNAME),'%VAL:'
DO IC = 1, A%N_ROW-1
   IF (A%ROW(IC) == 0) CYCLE
      IF (A%ROW(IC+1)-A%ROW(IC) < 10) THEN
         CFORM = "(I4,A,10E12.3)"
      ELSE IF (A%ROW(IC+1)-A%ROW(IC) < 20) THEN
         CFORM = "(I4,A,20E12.3)"
      ELSE IF (A%ROW(IC+1)-A%ROW(IC) < 30) THEN
         CFORM = "(I4,A,30E12.3)"
      ELSE
         CFORM = "(I4,A,40E12.3)"
      ENDIF
   WRITE(MSG%LU_DEBUG,CFORM) IC,':', (A%VAL(ICOL), ICOL=A%ROW(IC), A%ROW(IC+1)-1)
ENDDO
WRITE(MSG%LU_DEBUG,*) '============ END DEBUGGING MATRIX ', CNAME, ' AT ', TRIM(CTEXT)
ENDIF

END SUBROUTINE SCARC_DEBUG_CMATRIX


! ------------------------------------------------------------------------------------------------
! Print out vector information on level NL
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_DEBUG_LEVEL (NV, CVEC, NL)
USE SCARC_POINTERS, ONLY: L, G, VC
INTEGER, INTENT(IN) :: NV, NL
REAL (EB) :: VALUES(0:100)
INTEGER :: NM, II, JJ, KK, IC, NNX, NNY, NNZ
CHARACTER (*), INTENT(IN) :: CVEC

!IF (TYPE_SOLVER /= NSCARC_SOLVER_MAIN) RETURN
DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G
   VC => SCARC_POINT_TO_VECTOR (NM, NL, NV)

   NNX=MIN(10,L%NX)
   NNY=MIN(10,L%NY)
   NNZ=MIN(10,L%NZ)

   WRITE(MSG%LU_DEBUG,*) '=========================================================='
   WRITE(MSG%LU_DEBUG,2001) CVEC, NM, NL
   WRITE(MSG%LU_DEBUG,2002) G%NC, NNX, NNY, NNZ, NV, SIZE(VC)
   WRITE(MSG%LU_DEBUG,*) '=========================================================='
   IF ((IS_AMG.OR.IS_CG_SAMG) .AND. NL > NLEVEL_MIN) THEN

      WRITE(MSG%LU_DEBUG, '(6E14.6)') VC

   ELSE
   !IF (NL == NLEVEL_MIN) THEN
   DO KK = NNZ, 1, - 1
      DO JJ = NNY, 1, - 1
         DO II=1, NNX
            IF (IS_UNSTRUCTURED.AND.L%IS_SOLID(II,JJ,KK)) THEN
               VALUES(II)=-999999.0_EB
            ELSE
               IC=G%CELL_NUMBER(II,JJ,KK)
               IF (ABS(VC(IC))<1.0E-14_EB) THEN
                  VALUES(II)=0.0_EB
               ELSE
                  VALUES(II)=VC(IC)
               ENDIF
            ENDIF
         ENDDO
         !WRITE(MSG%LU_DEBUG, '(5E23.14)') (VALUES(II), II=4, NNX)
         WRITE(MSG%LU_DEBUG, '(12E14.6)') (VALUES(II), II=1, NNX)
         !WRITE(MSG%LU_DEBUG, '(12E11.3)') (VALUES(II), II=1, NNX
      ENDDO
      IF (.NOT. TWO_D) WRITE(MSG%LU_DEBUG, *) '----------------'
   ENDDO
   !ENDIF
   WRITE(MSG%LU_DEBUG, *) '---------------- Overlap ----------------'
   WRITE(MSG%LU_DEBUG, '(4E14.6)') (VC(IC), IC = G%NC+1, G%NCE)
   ENDIF
ENDDO

!CALL SCARC_MATLAB_VECTOR(NV, CVEC, NL)

!2000 FORMAT('=== ',A,' : ',A,' on mesh ',I8,' on level ',I8, ': NX, NY, NZ=',3I8,': NV=',I8)
2001 FORMAT('=== ',A,' on mesh ',I8,' on level ',I8)
2002 FORMAT('=== NC = ',I6, ': NX, NY, NZ=',3I6,': NV=',I6,': Size=',I8)
END SUBROUTINE SCARC_DEBUG_LEVEL


! ------------------------------------------------------------------------------------------------
! Print out vector information on level NL
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_DEBUG_LEVEL2 (NV, CVEC, NL)
USE SCARC_POINTERS, ONLY:  G, VC
INTEGER, INTENT(IN) :: NV, NL
INTEGER :: NM, IC
CHARACTER (*), INTENT(IN) :: CVEC

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G
   VC => SCARC_POINT_TO_VECTOR (NM, NL, NV)

   WRITE(MSG%LU_DEBUG,*) '=========================================================='
   WRITE(MSG%LU_DEBUG,*) CVEC, NM, NL
   WRITE(MSG%LU_DEBUG,*) G%NC, NV, SIZE(VC)
   WRITE(MSG%LU_DEBUG,*) '=========================================================='
   DO IC=1, G%NC
      WRITE(MSG%LU_DEBUG, '(I4,E25.16)') IC, VC(IC)
   ENDDO
ENDDO

END SUBROUTINE SCARC_DEBUG_LEVEL2


! ------------------------------------------------------------------------------------------------
! Debug specified quantity
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_DEBUG_QUANTITY(NTYPE, NL, CQUANTITY)
USE SCARC_POINTERS, ONLY: M, L, OL, G, OG, SV, A, AB
#ifdef WITH_MKL
USE SCARC_POINTERS, ONLY: AS
#endif
INTEGER, INTENT(IN) :: NTYPE, NL
INTEGER :: NM, NOM, IP, IC, ID, IW, I, J, IOR0, INBR, III, JJJ, KKK, IWG, NNX, NNY, NNZ
CHARACTER (*), INTENT(IN) :: CQUANTITY

SELECT CASE (NTYPE)

   ! ------------------------------------------------------------------------------------------------
   ! Debug system matrix A (corresponding to system type)
   ! ------------------------------------------------------------------------------------------------
   CASE (NSCARC_DEBUG_CMATRIX)

      SELECT CASE (SCARC_MATRIX_LEVEL(NL))

         CASE (NSCARC_MATRIX_COMPACT)

            DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
               CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G
               A => G%POISSONC
               WRITE(MSG%LU_DEBUG,1000) CQUANTITY, NM, NL
               WRITE(MSG%LU_DEBUG,*) '----------- SHOWING FULL COMPACT MATRIX ENTRIES'
               WRITE(MSG%LU_DEBUG,*) 'G%NC =',G%NC
               WRITE(MSG%LU_DEBUG,*) '=========================================='
               WRITE(MSG%LU_DEBUG,*) 'A%N_VAL =',A%N_VAL
               WRITE(MSG%LU_DEBUG,*) 'A%N_ROW =',A%N_ROW
               WRITE(MSG%LU_DEBUG,*) 'SIZE(A%VAL) =',SIZE(A%VAL)
               WRITE(MSG%LU_DEBUG,*) 'SIZE(A%COL) =',SIZE(A%COL)
               WRITE(MSG%LU_DEBUG,*) 'SIZE(A%ROW) =',SIZE(A%ROW)
               WRITE(MSG%LU_DEBUG,*) '---------------------- A%POS:'
               WRITE(MSG%LU_DEBUG,'(7i10)') (A%POS(IC), IC=-3,3)
               WRITE(MSG%LU_DEBUG,*) '---------------------- A%STENCIL:'
               WRITE(MSG%LU_DEBUG,'(7F10.3)') (A%STENCIL(IC), IC=-3,3)
               WRITE(MSG%LU_DEBUG,*) '---------------------- A%ROW:'
               WRITE(MSG%LU_DEBUG,'(7i10)') (A%ROW(IC), IC=1,A%N_ROW)
               WRITE(MSG%LU_DEBUG,*) '---------------------- A%COL:'
               DO IC = 1, A%N_ROW-1
                  WRITE(MSG%LU_DEBUG,'(i5,a,20i9)') IC,':',(A%COL(IP),IP=A%ROW(IC),A%ROW(IC+1)-1)
               ENDDO
               WRITE(MSG%LU_DEBUG,*) '---------------------- AC(',NM,',',NL,'):'
               DO IC = 1, A%N_ROW-1
                  WRITE(MSG%LU_DEBUG,'(i5,a,20f8.1)') IC,':',(A%VAL(IP),IP=A%ROW(IC),A%ROW(IC+1)-1)
               ENDDO
               !IF (ALLOCATED(A%RELAX)) THEN
               !   WRITE(MSG%LU_DEBUG,*) '---------------------- RELAX:'
               !   DO IC = 1, A%N_ROW-1
               !      WRITE(MSG%LU_DEBUG,'(i5,a,20f15.6)') IC,':',(A%RELAX(IP),IP=A%ROW(IC),A%ROW(IC+1)-1)
               !   ENDDO
               !ENDIF
#ifdef WITH_MKL
               !IF ((TYPE_METHOD == NSCARC_METHOD_LU) .AND. (TYPE_MKL(0) == NSCARC_MKL_GLOBAL)) THEN
               !   WRITE(MSG%LU_DEBUG,*) '---------------------- AG_COL:'
               !   DO IC = 1, G%NC
               !      WRITE(MSG%LU_DEBUG,'(i5,a,20i9)') IC,':',(A%COLG(IP),IP=A%ROW(IC),A%ROW(IC+1)-1)
               !      WRITE(MSG%LU_DEBUG,*)  IC,':',(A%COLG(IP),IP=A%ROW(IC),A%ROW(IC+1)-1)
               !   ENDDO
               !ENDIF
               !DO IC=1,A%N_ROW-1
               !   DO IP=A%ROW(IC),A%ROW(IC+1)-1
               !       WRITE(MSG%LU_DEBUG,'(2i10,F24.12)') IC,A%COL(IP),A%VAL(IP)
               !       WRITE(MSG%LU_DEBUG,*) IC,A%COL(IP),A%VAL(IP)
               !   ENDDO
               !   WRITE(MSG%LU_DEBUG,*)
               !ENDDO
               !DO IC=1,A%N_ROW-1
               !   DO IP=A%ROW(IC),A%ROW(IC+1)-1
               !       IF (IC == A%COL(IP)) WRITE(MSG%LU_DEBUG,'(2i10,F24.12)') IC,A%COL(IP),A%VAL(IP)
               !   ENDDO
               !ENDDO
#endif
               !IF (IS_STRUCTURED) THEN
               !   CALL SCARC_MATLAB_MATRIX(A%VAL, A%ROW, A%COL, G%NC, G%NC, NM, NL, 'A_struct')
               !ELSE
               !   CALL SCARC_MATLAB_MATRIX(A%VAL, A%ROW, A%COL, G%NC, G%NC, NM, NL, 'A_unstruct')
               !ENDIF

            ENDDO

         CASE (NSCARC_MATRIX_BANDWISE)

            DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
               CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G
               AB => G%POISSONB
               WRITE(MSG%LU_DEBUG,1000) CQUANTITY, NM, NL
               WRITE(MSG%LU_DEBUG,*) '----------- SHOWING FULL BANDWISE MATRIX ENTRIES'
               WRITE(MSG%LU_DEBUG,*) 'NC =',G%NC
               WRITE(MSG%LU_DEBUG,*) 'N_VAL  =',AB%N_VAL
               WRITE(MSG%LU_DEBUG,*) 'N_DIAG =',AB%N_DIAG
               WRITE(MSG%LU_DEBUG,*) 'N_STENCIL =',AB%N_STENCIL
               WRITE(MSG%LU_DEBUG,*) 'SIZE(AB%VAL) =',SIZE(AB%VAL)
               WRITE(MSG%LU_DEBUG,*) '---------------------- AB%POS:'
               WRITE(MSG%LU_DEBUG,'(7i10)') (AB%POS(IC), IC=-3,3)
               WRITE(MSG%LU_DEBUG,*) '---------------------- AB%STENCIL:'
               WRITE(MSG%LU_DEBUG,'(7F10.3)') (AB%STENCIL(IC), IC=-3,3)
               WRITE(MSG%LU_DEBUG,*) '---------------------- AB%OFFSET:'
               WRITE(MSG%LU_DEBUG,'(7i10)') (AB%OFFSET(IC), IC=-3,3)
               WRITE(MSG%LU_DEBUG,*) '---------------------- AB%LENGTH:'
               WRITE(MSG%LU_DEBUG,'(7i10)') (AB%LENGTH(IC), IC=-3,3)
               WRITE(MSG%LU_DEBUG,*) '---------------------- AB%SOURCE:'
               WRITE(MSG%LU_DEBUG,'(7i10)') (AB%SOURCE(IC), IC=-3,3)
               WRITE(MSG%LU_DEBUG,*) '---------------------- AB%TARGET:'
               WRITE(MSG%LU_DEBUG,'(7i10)') (AB%TARGET(IC), IC=-3,3)
               WRITE(MSG%LU_DEBUG,*) '---------------------- AB%VAL:'
               DO IC = 1, AB%N_DIAG
                 WRITE(MSG%LU_DEBUG,'(I5,A,7(F12.3))') IC,':',(AB%VAL(IC, ID),ID=1,AB%N_STENCIL)
               ENDDO
            ENDDO

      END SELECT

   ! ------------------------------------------------------------------------------------------------
   ! Debug symmetric system matrix AS
   ! ------------------------------------------------------------------------------------------------
   CASE (NSCARC_DEBUG_MATRIX_SYM)

      IF (NL /= NLEVEL_MAX) RETURN

#ifdef WITH_MKL
      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
         WRITE(MSG%LU_DEBUG,1000) CQUANTITY, NM, NL
         CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G
         AS => G%POISSONC_SYM
         WRITE(MSG%LU_DEBUG,*) '----------- SHOWING SYMMETRIC MATRIX ENTRIES'
         WRITE(MSG%LU_DEBUG,*) 'G%NC=', G%NC
         WRITE(MSG%LU_DEBUG,*) 'AS%ROW(G%NC)=', AS%ROW(G%NC)
         WRITE(MSG%LU_DEBUG,*) 'AS%ROW(G%NC+1)=', AS%ROW(G%NC+1)
         WRITE(MSG%LU_DEBUG,*) 'SIZE(ASC) =',    SIZE(AS%VAL)
         WRITE(MSG%LU_DEBUG,*) 'SIZE(AS%COLG) =',SIZE(AS%COLG)
         WRITE(MSG%LU_DEBUG,*) 'SIZE(AS%ROW) =',SIZE(AS%ROW)
         WRITE(MSG%LU_DEBUG,*) '---------------------- AS%ROW:', G%NC
         WRITE(MSG%LU_DEBUG,*) (AS%ROW(IC), IC=1,G%NC+1)
         WRITE(MSG%LU_DEBUG,*) '---------------------- AS%COLG:'
         DO IC = 1, G%NC
            WRITE(MSG%LU_DEBUG,*) IC,':',(AS%COLG(IP),IP=AS%ROW(IC),AS%ROW(IC+1)-1)
         ENDDO
         DO IC = 1, G%NC
            WRITE(MSG%LU_DEBUG,*) IC,':',(AS%VAL(IP),IP=AS%ROW(IC),AS%ROW(IC+1)-1)
         ENDDO
         DO IC=1,G%NC
            DO IP=AS%ROW(IC),AS%ROW(IC+1)-1
               WRITE(MSG%LU_DEBUG,*) IC,AS%COLG(IP),AS%VAL(IP)
            ENDDO
            WRITE(MSG%LU_DEBUG,*)
         ENDDO

      ENDDO
#endif

   ! ------------------------------------------------------------------------------------------------
   ! Debug strength matrix 
   ! ------------------------------------------------------------------------------------------------
   CASE (NSCARC_DEBUG_CONNECTION)

      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
         WRITE(MSG%LU_DEBUG,1000) CQUANTITY, NM, NL
         CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G
         WRITE(MSG%LU_DEBUG,*) '----------- SHOWING STRENGTH MATRIX ENTRIES'
         WRITE(MSG%LU_DEBUG,*) 'G%NC=', G%NC
         WRITE(MSG%LU_DEBUG,*) 'C%ROW(G%NC)=', C%ROW(G%NC)
         WRITE(MSG%LU_DEBUG,*) 'C%ROW(G%NC+1)=', C%ROW(G%NC+1)
         WRITE(MSG%LU_DEBUG,*) 'SIZE(SVAL) =',SIZE(C%VAL)
         WRITE(MSG%LU_DEBUG,*) 'SIZE(SCOL) =',SIZE(C%COL)
         WRITE(MSG%LU_DEBUG,*) 'SIZE(SROW) =',SIZE(C%ROW)
         WRITE(MSG%LU_DEBUG,*) '---------------------- SROW:', G%NC
         WRITE(MSG%LU_DEBUG,*) (C%ROW(IC), IC=1,G%NC+1)
         WRITE(MSG%LU_DEBUG,*) '---------------------- SCOL:'
         DO IC = 1, G%NC
            WRITE(MSG%LU_DEBUG,*) IC,':',(C%COL(IP),IP=C%ROW(IC),C%ROW(IC+1)-1)
         ENDDO
         DO IC = 1, G%NC
            WRITE(MSG%LU_DEBUG,*) IC,':',(C%VAL(IP),IP=C%ROW(IC),C%ROW(IC+1)-1)
         ENDDO

      ENDDO

   ! ------------------------------------------------------------------------------------------------
   ! Debug prolongation matrix 
   ! ------------------------------------------------------------------------------------------------
   CASE (NSCARC_DEBUG_PROLONGATION)

      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
         WRITE(MSG%LU_DEBUG,1000) CQUANTITY, NM, NL
         CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G
         WRITE(MSG%LU_DEBUG,*) '----------- SHOWING PROLONGATION MATRIX ENTRIES'
         WRITE(MSG%LU_DEBUG,*) 'G%NC=', G%NC
         WRITE(MSG%LU_DEBUG,*) 'P%ROW(G%NC)=', P%ROW(G%NC)
         WRITE(MSG%LU_DEBUG,*) 'P%ROW(G%NC+1)=', P%ROW(G%NC+1)
         WRITE(MSG%LU_DEBUG,*) 'SIZE(PVAL) =',SIZE(P%VAL)
         WRITE(MSG%LU_DEBUG,*) 'SIZE(PCOL) =',SIZE(P%COL)
         WRITE(MSG%LU_DEBUG,*) 'SIZE(PROW) =',SIZE(P%ROW)
         WRITE(MSG%LU_DEBUG,*) '---------------------- PROW:', G%NC
         WRITE(MSG%LU_DEBUG,*) (P%ROW(IC), IC=1,G%NC+1)
         WRITE(MSG%LU_DEBUG,*) '---------------------- PCOL:'
         DO IC = 1, G%NC
            WRITE(MSG%LU_DEBUG,*) IC,':',(P%COL(IP),IP=P%ROW(IC),P%ROW(IC+1)-1)
         ENDDO
         DO IC = 1, G%NC
            WRITE(MSG%LU_DEBUG,*) IC,':',(P%VAL(IP),IP=P%ROW(IC),P%ROW(IC+1)-1)
         ENDDO

      ENDDO

   ! ------------------------------------------------------------------------------------------------
   ! Debug strength matrix SC
   ! ------------------------------------------------------------------------------------------------
   CASE (NSCARC_DEBUG_RESTRICTION)

      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
         WRITE(MSG%LU_DEBUG,1000) CQUANTITY, NM, NL
         CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G
         WRITE(MSG%LU_DEBUG,*) '----------- SHOWING RESTRICTION MATRIX ENTRIES'
         WRITE(MSG%LU_DEBUG,*) 'G%NC=', G%NC
         WRITE(MSG%LU_DEBUG,*) 'R%ROW(G%NC)=', R%ROW(G%NC)
         WRITE(MSG%LU_DEBUG,*) 'R%ROW(G%NC+1)=', R%ROW(G%NC+1)
         WRITE(MSG%LU_DEBUG,*) 'SIZE(R%VAL) =',SIZE(R%VAL)
         WRITE(MSG%LU_DEBUG,*) 'SIZE(R%COL) =',SIZE(R%COL)
         WRITE(MSG%LU_DEBUG,*) 'SIZE(R%ROW) =',SIZE(R%ROW)
         WRITE(MSG%LU_DEBUG,*) '---------------------- R%ROW:', G%NC
         WRITE(MSG%LU_DEBUG,*) (R%ROW(IC), IC=1,G%NC+1)
         WRITE(MSG%LU_DEBUG,*) '---------------------- R%COL:'
         DO IC = 1, G%NC
            WRITE(MSG%LU_DEBUG,*) IC,':',(R%COLG(IP),IP=R%ROW(IC),R%ROW(IC+1)-1)
         ENDDO
         DO IC = 1, G%NC
            WRITE(MSG%LU_DEBUG,*) IC,':',(R%VAL(IP),IP=R%ROW(IC),R%ROW(IC+1)-1)
         ENDDO

      ENDDO

   ! ------------------------------------------------------------------------------------------------
   ! Debug FACEINFO
   ! ------------------------------------------------------------------------------------------------
   CASE (NSCARC_DEBUG_FACE)

      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
         L => SCARC(NM)%LEVEL(NL)
         WRITE(MSG%LU_DEBUG,1000) CQUANTITY, NM, NL
         DO IOR0 = -3, 3
            IF (IOR0 == 0) CYCLE
            WRITE(MSG%LU_DEBUG,*) '========================================='
            WRITE(MSG%LU_DEBUG,*) '============= DEBUGGING FACE(',IOR0,'): FOR LEVEL ', NL
            WRITE(MSG%LU_DEBUG,*) '========================================='
            WRITE(MSG%LU_DEBUG,*) 'FACE(.)%DH: '
            WRITE(MSG%LU_DEBUG,'(12F8.2)')  L%FACE(IOR0)%DH
            WRITE(MSG%LU_DEBUG,*) 'FACE(.)%NOP:', L%FACE(IOR0)%NOP
            WRITE(MSG%LU_DEBUG,*) 'FACE(.)%NCW:', L%FACE(IOR0)%NCW
            WRITE(MSG%LU_DEBUG,*) 'FACE(.)%NX: ', L%FACE(IOR0)%NX
            WRITE(MSG%LU_DEBUG,*) 'FACE(.)%NY: ', L%FACE(IOR0)%NY
            WRITE(MSG%LU_DEBUG,*) 'FACE(.)%NZ: ', L%FACE(IOR0)%NZ
            WRITE(MSG%LU_DEBUG,*) 'FACE(.)%NW0:', L%FACE(IOR0)%NCW0
            WRITE(MSG%LU_DEBUG,*) 'FACE(.)%N_NEIGHBORS: ', L%FACE(IOR0)%N_NEIGHBORS
            WRITE(MSG%LU_DEBUG,*) 'FACE(.)%SCAL_BOUNDARY: ', L%FACE(IOR0)%SCAL_BOUNDARY
            WRITE(MSG%LU_DEBUG,*) 'FACE(.)%SCAL_DIRICHLET: ', L%FACE(IOR0)%SCAL_DIRICHLET
            WRITE(MSG%LU_DEBUG,*) 'FACE(.)%SCAL_NEUMANN: ', L%FACE(IOR0)%SCAL_NEUMANN
            WRITE(MSG%LU_DEBUG,*) 'FACE(.)%SCAL_FACE: ', L%FACE(IOR0)%SCAL_FACE
            WRITE(MSG%LU_DEBUG,*) '----------------------------------------------'
            WRITE(MSG%LU_DEBUG,*)
            IF (L%FACE(IOR0)%N_NEIGHBORS /= 0) THEN
               DO INBR=1,L%FACE(IOR0)%N_NEIGHBORS
                  NOM = L%FACE(IOR0)%NEIGHBORS(INBR)
                  OL => SCARC(NM)%OSCARC(NOM)%LEVEL(NL)
                  SELECT CASE(TYPE_GRID)
                  CASE (NSCARC_GRID_STRUCTURED)
                     OG => OL%STRUCTURED
                  CASE (NSCARC_GRID_UNSTRUCTURED)
                     OG => OL%UNSTRUCTURED
                  END SELECT
                  WRITE(MSG%LU_DEBUG,*) 'N_NEIGHBORS:', L%FACE(IOR0)%N_NEIGHBORS
                  WRITE(MSG%LU_DEBUG,*) 'NOM:', NOM
                  WRITE(MSG%LU_DEBUG,*) 'SIZE(OG%ICG_TO_IWG)=',SIZE(OG%ICG_TO_IWG)
                  WRITE(MSG%LU_DEBUG,*) 'SIZE(OG%ICG_TO_ICE)=',SIZE(OG%ICG_TO_ICE)
                  WRITE(MSG%LU_DEBUG,'(a,i8,a,2i8)') '---OG(',NOM,')%GHOST_LASTW(.):',OG%NCG
                  WRITE(MSG%LU_DEBUG,*)
                  WRITE(MSG%LU_DEBUG,'(a,i8,a)') '------OG(',NOM,')%ICG_TO_IWG:'
                  DO IW = 1, OG%NCG
                     WRITE(MSG%LU_DEBUG,'(16i8)') OG%ICG_TO_IWG(IW)
                  ENDDO
                  WRITE(MSG%LU_DEBUG,*)
                  WRITE(MSG%LU_DEBUG,'(a,i8,a)') '------OG(',NOM,')%ICG_TO_ICE:'
                  DO IW = 1, OG%NCG
                     WRITE(MSG%LU_DEBUG,'(16i8)') OG%ICG_TO_ICE(IW, 1)
                  ENDDO
                  WRITE(MSG%LU_DEBUG,*)
                  WRITE(MSG%LU_DEBUG,'(a,i8,a)') '------OG(',NOM,')%ICG_TO_ICW:'
                  DO IW = 1, OG%NCG
                     WRITE(MSG%LU_DEBUG,'(16i8)') OG%ICG_TO_ICW(IW, 1)
                  ENDDO
               ENDDO
            ENDIF
         ENDDO
      ENDDO

   ! ------------------------------------------------------------------------------------------------
   ! Debug Pressure information 
   ! ------------------------------------------------------------------------------------------------
   CASE (NSCARC_DEBUG_PRESSURE)

      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

         CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G

         NNX=MIN(8,L%NX)
         NNY=MIN(8,L%NY)
         NNZ=MIN(8,L%NZ)

         WRITE(MSG%LU_DEBUG,*) '========= PRESSURE ========= ', NM, NL, CQUANTITY
         IF (PREDICTOR) THEN
            WRITE(MSG%LU_DEBUG,*) 'RHO'
            DO KKK = NNZ+1, 0, -1
               WRITE(MSG%LU_DEBUG,'(I4,10E14.6)') KKK, (M%RHO(III, 1, KKK), III=0,NNX+1)
            ENDDO
            WRITE(MSG%LU_DEBUG,*) 'H'
            DO KKK = NNZ+1, 0, -1
               WRITE(MSG%LU_DEBUG,'(I4,10E14.6)') KKK, (M%H(III, 1, KKK), III=0,NNX+1)
            ENDDO
         ELSE
            WRITE(MSG%LU_DEBUG,*) 'RHOS'
            DO KKK = NNZ+1, 0, -1
               WRITE(MSG%LU_DEBUG,'(I4,10E14.6)') KKK, (M%RHOS(III, 1, KKK), III=0,NNX+1)
            ENDDO
            WRITE(MSG%LU_DEBUG,*) 'HS'
            DO KKK = NNZ+1, 0, -1
               WRITE(MSG%LU_DEBUG,'(I4,10E14.6)') KKK, (M%HS(III, 1, KKK), III=0,NNX+1)
            ENDDO
         ENDIF
         WRITE(MSG%LU_DEBUG,*) 'FVX'
         DO KKK = NNZ+1, 0, -1
            WRITE(MSG%LU_DEBUG,'(I4,10E14.6)') KKK, (M%FVX(III, 1, KKK), III=0,NNX+1)
         ENDDO
         WRITE(MSG%LU_DEBUG,*) 'FVY'
         DO KKK = NNZ+1, 0, -1
            WRITE(MSG%LU_DEBUG,'(I4,10E14.6)') KKK, (M%FVY(III, 1, KKK), III=0,NNX+1)
         ENDDO
         WRITE(MSG%LU_DEBUG,*) 'FVZ'
         DO KKK = NNZ+1, 0, -1
            WRITE(MSG%LU_DEBUG,'(I4,10E14.6)') KKK, (M%FVZ(III, 1, KKK), III=0,NNX+1)
         ENDDO
         WRITE(MSG%LU_DEBUG,*) 'KRES'
         DO KKK = NNZ+1, 0, -1
            WRITE(MSG%LU_DEBUG,'(I4,10E14.6)') KKK, (M%KRES(III, 1, KKK), III=0,NNX+1)
         ENDDO
         WRITE(MSG%LU_DEBUG,*) 'M%PRHS'
         DO KKK = NNZ+1, 1, -1
            WRITE(MSG%LU_DEBUG,'(I4,9E14.6)') KKK, (M%PRHS(III, 1, KKK), III=1,NNX+1)
         ENDDO
         !WRITE(MSG%LU_DEBUG,*) 'P%PRHS'
         !DO KKK = NNZ+1, 1, -1
         !   WRITE(MSG%LU_DEBUG,'(I4,9E14.6)') KKK, (P%PRHS(III, 1, KKK), III=1,NNX+1)
         !ENDDO

      ENDDO

   ! ------------------------------------------------------------------------------------------------
   ! Debug stack information
   ! ------------------------------------------------------------------------------------------------
   CASE (NSCARC_DEBUG_STACK)

      WRITE(MSG%LU_DEBUG,*) 'N_STACK_TOTAL=',N_STACK_TOTAL
      DO I = 1, N_STACK_TOTAL
         SV  => STACK(I)%SOLVER
         WRITE(MSG%LU_DEBUG,*) '===================== STACK ', I,' ======================'
         WRITE(MSG%LU_DEBUG,*) '-------------------- SOLVER:'
         WRITE(MSG%LU_DEBUG,*) 'NAME=',SV%CNAME
         WRITE(MSG%LU_DEBUG,*) '-- SETTING:'
         WRITE(MSG%LU_DEBUG,*) 'EPS   = ',SV%EPS
         WRITE(MSG%LU_DEBUG,*) 'RES   = ',SV%RES
         WRITE(MSG%LU_DEBUG,*) 'RESIN = ',SV%RESIN
         WRITE(MSG%LU_DEBUG,*) 'OMEGA = ',SV%OMEGA
         WRITE(MSG%LU_DEBUG,*) 'ITE   = ',SV%ITE
         WRITE(MSG%LU_DEBUG,*) 'NIT   = ',SV%NIT
         WRITE(MSG%LU_DEBUG,*) '-- TYPES:'
         WRITE(MSG%LU_DEBUG,*) 'TYPE_PARENT   = ',SV%TYPE_PARENT
         WRITE(MSG%LU_DEBUG,*) 'TYPE_SOLVER   = ',SV%TYPE_SOLVER
         WRITE(MSG%LU_DEBUG,*) 'TYPE_STAGE    = ',SV%TYPE_STAGE
         WRITE(MSG%LU_DEBUG,*) 'TYPE_SCOPE(0) = ',SV%TYPE_SCOPE(0)
         WRITE(MSG%LU_DEBUG,*) 'TYPE_PRECON   = ',SV%TYPE_RELAX
         WRITE(MSG%LU_DEBUG,*) 'TYPE_ACCURACY = ',SV%TYPE_ACCURACY
         WRITE(MSG%LU_DEBUG,*) 'TYPE_INTERPOL = ',SV%TYPE_INTERPOL
         WRITE(MSG%LU_DEBUG,*) 'TYPE_CYCLING  = ',SV%TYPE_CYCLING
         WRITE(MSG%LU_DEBUG,*) 'TYPE_TWOLEVEL = ',SV%TYPE_TWOLEVEL
         WRITE(MSG%LU_DEBUG,*) '-- POINTERS:'
         WRITE(MSG%LU_DEBUG,*) 'X   = ',SV%X
         WRITE(MSG%LU_DEBUG,*) 'B   = ',SV%B
         WRITE(MSG%LU_DEBUG,*) 'D   = ',SV%D
         WRITE(MSG%LU_DEBUG,*) 'E   = ',SV%E
         WRITE(MSG%LU_DEBUG,*) 'R   = ',SV%R
         WRITE(MSG%LU_DEBUG,*) 'V   = ',SV%V
         WRITE(MSG%LU_DEBUG,*) 'Y   = ',SV%Y
         WRITE(MSG%LU_DEBUG,*) 'Z   = ',SV%Z
      ENDDO

   ! ------------------------------------------------------------------------------------------------
   ! Debug WALLINFO
   ! ------------------------------------------------------------------------------------------------
   CASE (NSCARC_DEBUG_WALL)

      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
         CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G
         WRITE(MSG%LU_DEBUG,1000) CQUANTITY, NM, NL
         WRITE(MSG%LU_DEBUG,*) 'SIZE(G%ICE_TO_IWG)=',SIZE(G%ICE_TO_IWG)
         WRITE(MSG%LU_DEBUG,*) 'NM  =',NM
         WRITE(MSG%LU_DEBUG,*) 'NL  =',NL
         WRITE(MSG%LU_DEBUG,*) 'NC  =',G%NC
         WRITE(MSG%LU_DEBUG,*) 'NCE =',G%NCE
         WRITE(MSG%LU_DEBUG,*) 'N_WALL_CELLS  =',L%N_WALL_CELLS
         WRITE(MSG%LU_DEBUG,*)
         WRITE(MSG%LU_DEBUG,'(a,I8,a)') '------G%ICE_TO_IWG:'
         WRITE(MSG%LU_DEBUG,'(8I8)') (G%ICE_TO_IWG(IW), IW = G%NC+1, G%NCE)
         WRITE(MSG%LU_DEBUG,*)
         IF (N_DIRIC_GLOBAL(NLEVEL_MIN) == 0) THEN
            WRITE(MSG%LU_DEBUG,'(a,I8,a)') '------G%ICE_TO_ICN:'
            WRITE(MSG%LU_DEBUG,'(8I8)') (G%ICE_TO_ICN(IW), IW = G%NC+1, G%NCE)
         ENDIF
         WRITE(MSG%LU_DEBUG,*)
         IF (NL == 1) THEN
            WRITE(MSG%LU_DEBUG,*) '============= WALL(.)%IXG:', NM
            WRITE(MSG%LU_DEBUG,'(16I6)') (G%WALL(IW)%IXG, IW=1,L%N_WALL_CELLS)
            WRITE(MSG%LU_DEBUG,*) '============= WALL(.)%IYG:', NM
            WRITE(MSG%LU_DEBUG,'(16I6)') (G%WALL(IW)%IYG, IW=1,L%N_WALL_CELLS)
            WRITE(MSG%LU_DEBUG,*) '============= WALL(.)%IZG:', NM
            WRITE(MSG%LU_DEBUG,'(16I6)') (G%WALL(IW)%IZG, IW=1,L%N_WALL_CELLS)
            WRITE(MSG%LU_DEBUG,*)
            WRITE(MSG%LU_DEBUG,*) '============= WALL(.)%IXW:', NM
            WRITE(MSG%LU_DEBUG,'(16I6)') (G%WALL(IW)%IXW, IW=1,L%N_WALL_CELLS)
            WRITE(MSG%LU_DEBUG,*) '============= WALL(.)%IYW:', NM
            WRITE(MSG%LU_DEBUG,'(16I6)') (G%WALL(IW)%IYW, IW=1,L%N_WALL_CELLS)
            WRITE(MSG%LU_DEBUG,*) '============= WALL(.)%IZW:', NM
            WRITE(MSG%LU_DEBUG,'(16I6)') (G%WALL(IW)%IZW, IW=1,L%N_WALL_CELLS)
            WRITE(MSG%LU_DEBUG,*)
            WRITE(MSG%LU_DEBUG,*) '============= WALL(.)%IXN(1):', NM
            WRITE(MSG%LU_DEBUG,'(16I6)') (G%WALL(IW)%IXN(1), IW=1,L%N_WALL_CELLS)
            WRITE(MSG%LU_DEBUG,*) '============= WALL(.)%IYN(1):', NM
            WRITE(MSG%LU_DEBUG,'(16I6)') (G%WALL(IW)%IYN(1), IW=1,L%N_WALL_CELLS)
            WRITE(MSG%LU_DEBUG,*) '============= WALL(.)%IZN(1):', NM
            WRITE(MSG%LU_DEBUG,'(16I6)') (G%WALL(IW)%IZN(1), IW=1,L%N_WALL_CELLS)
            WRITE(MSG%LU_DEBUG,*)
            WRITE(MSG%LU_DEBUG,*) '============= WALL(.)%IXN(2):', NM
            WRITE(MSG%LU_DEBUG,'(16I6)') (G%WALL(IW)%IXN(2), IW=1,L%N_WALL_CELLS)
            WRITE(MSG%LU_DEBUG,*) '============= WALL(.)%IYN(2):', NM
            WRITE(MSG%LU_DEBUG,'(16I6)') (G%WALL(IW)%IYN(2), IW=1,L%N_WALL_CELLS)
            WRITE(MSG%LU_DEBUG,*) '============= WALL(.)%IZN(2):', NM
            WRITE(MSG%LU_DEBUG,'(16I6)') (G%WALL(IW)%IZN(2), IW=1,L%N_WALL_CELLS)
         ENDIF
         WRITE(MSG%LU_DEBUG,*)
         WRITE(MSG%LU_DEBUG,*) '============= WALL(.)%BTYPE:', NM
         WRITE(MSG%LU_DEBUG,'(16I6)') (G%WALL(IW)%BTYPE, IW=1,L%N_WALL_CELLS)
         WRITE(MSG%LU_DEBUG,*)
         WRITE(MSG%LU_DEBUG,*) '============= WALL(.)%IOR:', NM
         WRITE(MSG%LU_DEBUG,'(16I6)') (G%WALL(IW)%IOR, IW=1,L%N_WALL_CELLS)
         WRITE(MSG%LU_DEBUG,*)
         WRITE(MSG%LU_DEBUG,*) '============= WALL(.)%NOM:', NM
         WRITE(MSG%LU_DEBUG,'(16I6)') (G%WALL(IW)%NOM, IW=1,L%N_WALL_CELLS)
         WRITE(MSG%LU_DEBUG,*)
         WRITE(MSG%LU_DEBUG,*)
         WRITE(MSG%LU_DEBUG,*) '============= WALL(.)%ICW:', NM
         DO IW=1,L%N_WALL_CELLS
            WRITE(MSG%LU_DEBUG,'(a,I8, a,I8)') 'IW=',IW,':',G%WALL(IW)%ICW
         ENDDO
         WRITE(MSG%LU_DEBUG,*)
         WRITE(MSG%LU_DEBUG,*) '============= WALL(.)%ICE:', NM
         DO IW=1,L%N_WALL_CELLS
            IF (G%WALL(IW)%NOM /=0) &
               WRITE(MSG%LU_DEBUG,'(a,I8, a,I8)') 'IW=',IW,':',G%WALL(IW)%ICE
         ENDDO
         WRITE(MSG%LU_DEBUG,*)
         WRITE(MSG%LU_DEBUG,*)
         WRITE(MSG%LU_DEBUG,*) '============= WALL(.)%ICG:', NM
         DO IW=1,L%N_WALL_CELLS
            IF (G%WALL(IW)%NOM /=0) &
               WRITE(MSG%LU_DEBUG,'(a,I8, a,I8)') 'IW=',IW,':',G%WALL(IW)%ICG
         ENDDO
         WRITE(MSG%LU_DEBUG,*) '====================================================='
         WRITE(MSG%LU_DEBUG,*) ' Plotting out M%WALL-structure'
         WRITE(MSG%LU_DEBUG,*) '====================================================='
         WRITE(MSG%LU_DEBUG,*) ' M%WALL(.)%WALL_INDEX'
         WRITE(MSG%LU_DEBUG,'(16I6)') (M%WALL(IW)%WALL_INDEX, IW=1,L%N_WALL_CELLS)
         WRITE(MSG%LU_DEBUG,*) ' M%WALL(.)%SURF_INDEX'
         WRITE(MSG%LU_DEBUG,'(16I6)') (M%WALL(IW)%SURF_INDEX, IW=1,L%N_WALL_CELLS)
         WRITE(MSG%LU_DEBUG,*) ' M%WALL(.)%BACK_INDEX'
         WRITE(MSG%LU_DEBUG,'(16I6)') (M%WALL(IW)%BACK_INDEX, IW=1,L%N_WALL_CELLS)
         WRITE(MSG%LU_DEBUG,*) ' M%WALL(.)%BOUNDARY_TYPE'
         WRITE(MSG%LU_DEBUG,'(16I6)') (M%WALL(IW)%BOUNDARY_TYPE, IW=1,L%N_WALL_CELLS)
         WRITE(MSG%LU_DEBUG,*) ' M%WALL(.)%OBST_INDEX'
         WRITE(MSG%LU_DEBUG,'(16I6)') (M%WALL(IW)%OBST_INDEX, IW=1,L%N_WALL_CELLS)
         WRITE(MSG%LU_DEBUG,*) ' M%WALL(.)%PRESSURE_BC_INDEX'
         WRITE(MSG%LU_DEBUG,'(16I6)') (M%WALL(IW)%PRESSURE_BC_INDEX, IW=1,L%N_WALL_CELLS)
         WRITE(MSG%LU_DEBUG,*) ' M%WALL(.)%ONE_D%PRESSURE_ZONE'
         WRITE(MSG%LU_DEBUG,'(16I6)') (M%WALL(IW)%ONE_D%PRESSURE_ZONE, IW=1,L%N_WALL_CELLS)
         WRITE(MSG%LU_DEBUG,*) ' M%WALL(.)%VENT_INDEX'
         WRITE(MSG%LU_DEBUG,'(16I6)') (M%WALL(IW)%VENT_INDEX, IW=1,L%N_WALL_CELLS)
         WRITE(MSG%LU_DEBUG,*) ' M%WALL(.)%NOM'
         WRITE(MSG%LU_DEBUG,'(16I6)') (M%EXTERNAL_WALL(IW)%NOM, IW=1,L%N_WALL_CELLS_EXT)
         WRITE(MSG%LU_DEBUG,*) ' M%WALL(.)%ONE_D%II'
         WRITE(MSG%LU_DEBUG,'(16I6)') (M%WALL(IW)%ONE_D%II, IW=1,L%N_WALL_CELLS)
         WRITE(MSG%LU_DEBUG,*) ' M%WALL(.)%ONE_D%JJ'
         WRITE(MSG%LU_DEBUG,'(16I6)') (M%WALL(IW)%ONE_D%JJ, IW=1,L%N_WALL_CELLS)
         WRITE(MSG%LU_DEBUG,*) ' M%WALL(.)%ONE_D%KK'
         WRITE(MSG%LU_DEBUG,'(16I6)') (M%WALL(IW)%ONE_D%KK, IW=1,L%N_WALL_CELLS)
         WRITE(MSG%LU_DEBUG,*) ' M%WALL(.)%ONE_D%IIG'
         WRITE(MSG%LU_DEBUG,'(16I6)') (M%WALL(IW)%ONE_D%IIG, IW=1,L%N_WALL_CELLS)
         WRITE(MSG%LU_DEBUG,*) ' M%WALL(.)%ONE_D%JJG'
         WRITE(MSG%LU_DEBUG,'(16I6)') (M%WALL(IW)%ONE_D%JJG, IW=1,L%N_WALL_CELLS)
         WRITE(MSG%LU_DEBUG,*) ' M%WALL(.)%ONE_D%KKG'
         WRITE(MSG%LU_DEBUG,'(16I6)') (M%WALL(IW)%ONE_D%KKG, IW=1,L%N_WALL_CELLS)
         WRITE(MSG%LU_DEBUG,*) ' M%WALL(.)%ONE_D%N_LAYER_CELLS'
         WRITE(MSG%LU_DEBUG,'(16I6)') (M%WALL(IW)%ONE_D%IOR, IW=1,L%N_WALL_CELLS)
         WRITE(MSG%LU_DEBUG,*) ' M%WALL(.)%ONE_D%IOR'
         !WRITE(MSG%LU_DEBUG,'(16I6)') (M%WALL(IW)%ONE_D%N_LAYER_CELLS, IW=1,L%N_WALL_CELLS)

      ENDDO
      !ENDIF

   ! ------------------------------------------------------------------------------------------------
   ! Debug complete grid information
   ! ------------------------------------------------------------------------------------------------
   CASE (NSCARC_DEBUG_GRID)
      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
         CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G
         WRITE(MSG%LU_DEBUG,*) 'M%N_OBST=',M%N_OBST
         WRITE(MSG%LU_DEBUG,*) 'M%OBST ... I1, I2, J1, J2, K1, K2'
         DO IWG = 1, M%N_OBST
            WRITE(MSG%LU_DEBUG,'(6I8)') M%OBSTRUCTION(IWG)%I1,M%OBSTRUCTION(IWG)%I2,&
               M%OBSTRUCTION(IWG)%J1,M%OBSTRUCTION(IWG)%J2,&
               M%OBSTRUCTION(IWG)%K1,M%OBSTRUCTION(IWG)%K2
         ENDDO
         WRITE(MSG%LU_DEBUG,*) 'M%N_OBST=',M%N_OBST
         WRITE(MSG%LU_DEBUG,*) 'M%N_WALL_CELLS=',M%N_WALL_CELLS
         WRITE(MSG%LU_DEBUG,*) 'M%CELL_INDEX:'
         DO JJJ = L%NY+1,0,-1
            WRITE(MSG%LU_DEBUG,*) ' ------------- JJJ = ', JJJ
            DO KKK = L%NZ+1,0,-1
               WRITE(MSG%LU_DEBUG,*) (M%CELL_INDEX(III,JJJ,KKK), III=0,L%NX+1)
            ENDDO
         ENDDO
         WRITE(MSG%LU_DEBUG,*) 'M%WALL(.)%BOUNDARY_TYPE:'
         WRITE(MSG%LU_DEBUG,'(16i6)') (M%WALL(IWG)%BOUNDARY_TYPE, IWG=1, L%N_WALL_CELLS)

         WRITE(MSG%LU_DEBUG,*) 'CELL_COUNT:', CELL_COUNT(NM)

         WRITE(MSG%LU_DEBUG,*) 'M%WALL_INDEX:'
         WRITE(MSG%LU_DEBUG,'(i8,a,6i8)') ( IWG, ' : ', &
            M%WALL_INDEX(IWG, 1),  &
            M%WALL_INDEX(IWG,-1),  &
            M%WALL_INDEX(IWG, 2),  &
            M%WALL_INDEX(IWG,-2),  &
            M%WALL_INDEX(IWG, 3),  &
            M%WALL_INDEX(IWG,-3),  &
            IWG=1, CELL_COUNT(NM))

         WRITE(MSG%LU_DEBUG,*) 'M%WALL(.)%ONE_D% IOR,II,JJ,KK, BOUNDARY_TYPE, BTYPE, PRESSURE_BC_INDEX:'
         DO IWG = 1, L%N_WALL_CELLS
            WRITE(MSG%LU_DEBUG,'(9I8)') &
               IWG,M%WALL(IWG)%ONE_D%IOR,M%WALL(IWG)%ONE_D%II,M%WALL(IWG)%ONE_D%JJ,M%WALL(IWG)%ONE_D%KK,&
               M%WALL(IWG)%BOUNDARY_TYPE, G%WALL(IWG)%BTYPE, M%WALL(IWG)%PRESSURE_BC_INDEX
         ENDDO

         WRITE(MSG%LU_DEBUG,*) 'GRIG%CELL_NUMBER(...)'
         DO JJJ=L%NY+1,0,-1
            WRITE(MSG%LU_DEBUG,*) ' ------------- JJJ = ', JJJ
            DO KKK = L%NZ+1,0,-1
               WRITE(MSG%LU_DEBUG,*) (G%CELL_NUMBER(III,JJJ,KKK),III=0,L%NX+1)
            ENDDO
         ENDDO
         WRITE(MSG%LU_DEBUG,*) 'GRIL%IS_SOLID(...)'
         DO JJJ=L%NY+1,0,-1
            WRITE(MSG%LU_DEBUG,*) ' ------------- JJJ = ', JJJ
            DO KKK = L%NZ+1,0,-1
               WRITE(MSG%LU_DEBUG,*) (L%IS_SOLID(III,JJJ,KKK),III=0,L%NX+1)
            ENDDO
         ENDDO
         WRITE(MSG%LU_DEBUG,*) 'M%CELL_INDEX:'
         DO JJJ=L%NY+1,0,-1
            WRITE(MSG%LU_DEBUG,*) ' ------------- JJJ = ', JJJ
            DO KKK = L%NZ+1,0,-1
               WRITE(MSG%LU_DEBUG,*) (M%CELL_INDEX(III,JJJ,KKK),III=0,L%NX+1)
            ENDDO
         ENDDO
      ENDDO

   ! ------------------------------------------------------------------------------------------------
   ! Debug PRESSURE_BC_INDEX
   ! ------------------------------------------------------------------------------------------------
   CASE (NSCARC_DEBUG_BDRY)

      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
         L => SCARC(NM)%LEVEL(NL)
         WRITE(MSG%LU_DEBUG,1000) CQUANTITY, NM, NL
         WRITE(MSG%LU_DEBUG, '(16I6)') (G%WALL(J)%BTYPE, J=1,L%N_WALL_CELLS)
      ENDDO

END SELECT

1000 FORMAT('======================================================================================',/, &
   '=== ', A30,' for mesh ',i3,' on level ', i3, /, &
   '======================================================================================')

END SUBROUTINE SCARC_DEBUG_QUANTITY


! ------------------------------------------------------------------------------------------------
! Print out vector information on level NL for matlab
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_MATLAB_VECTOR (NV, CVEC, NL)
USE SCARC_POINTERS, ONLY: G, VC
INTEGER, INTENT(IN) :: NV, NL
INTEGER :: NM 
CHARACTER (*), INTENT(IN) :: CVEC
INTEGER :: JC, MVEC
CHARACTER(60) :: CNAME, CFORM

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G

   VC => SCARC_POINT_TO_VECTOR (NM, NL, NV)
   WRITE (CNAME, '(A,A1,A,i2.2,A,i2.2,A)') 'matlab/',CVEC,'_mesh',NM,'_level',NL,'_vec.txt'
   WRITE (CFORM, '(I3, A)' ) G%NC-1, "(F7.2,;),F7.2"
   MVEC=GET_FILE_NUMBER()

   OPEN(MVEC,FILE=CNAME)
   WRITE(MVEC, *) CVEC, ' = ['
   WRITE(MVEC,'(8F12.2)') (VC(JC),JC=1,G%NC)
   WRITE(MVEC, *) ' ]'
   CLOSE(MVEC)

ENDDO

END SUBROUTINE SCARC_MATLAB_VECTOR


! ------------------------------------------------------------------------------------------------
! Print out matrix information on level NL for matlab
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_MATLAB_MATRIX(VAL, ROW, COL, NC1, NC2, NM, NL, CNAME)
REAL(EB), DIMENSION(:), INTENT(IN) :: VAL
INTEGER, DIMENSION(:), INTENT(IN) :: ROW
INTEGER, DIMENSION(:), INTENT(IN) :: COL
INTEGER, INTENT(IN) :: NM, NL, NC1, NC2
CHARACTER(*), INTENT(IN) :: CNAME
INTEGER :: IC, JC, ICOL, MMATRIX
CHARACTER(60) :: CFILE, CFORM
REAL(EB) :: MATRIX_LINE(1000)

RETURN
WRITE (CFILE, '(A,A,A,i2.2,A,i2.2,A)') 'matlab/',TRIM(CNAME),'_mesh',NM,'_level',NL,'_mat.txt'
!WRITE (CFORM, '(A,I3, 2A)' ) "(", NC2-1, "(F9.3,','),F9.3,';')"
!WRITE (CFORM, '(A,I3, 2A)' ) "(", NC2-1, "(F9.3,' '),F9.3,' ')"
WRITE (CFORM, '(A,I3, 2A)' ) "(", NC2-1, "(F9.3,' '),F9.3,' ')"
MMATRIX=GET_FILE_NUMBER()
OPEN(MMATRIX,FILE=CFILE)
!WRITE(MMATRIX, *) CNAME, ' = ['
DO IC = 1, NC1
   MATRIX_LINE=0.0_EB
   DO JC = 1, NC2
      DO  ICOL= ROW(IC), ROW(IC+1)-1
         !IF (COL(ICOL)==JC) MATRIX_LINE(JC)=VAL(ICOL)/25.0_EB
         IF (COL(ICOL)==JC) MATRIX_LINE(JC)=VAL(ICOL)
      ENDDO
   ENDDO
   WRITE(MMATRIX, CFORM) (MATRIX_LINE(JC),JC=1,NC2)
ENDDO
!WRITE(MMATRIX, *) ' ]'
CLOSE(MMATRIX)

END SUBROUTINE SCARC_MATLAB_MATRIX

! ------------------------------------------------------------------------------------------------
! Print out matrix information on level NL for matlab
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_PYTHON_MATRIX(NL, CNAME)
USE SCARC_POINTERS, ONLY : G
INTEGER, INTENT(IN) :: NL
CHARACTER(*), INTENT(IN) :: CNAME
INTEGER :: NM, IC, JC, ICOL, MVAL, MCOL, MROW, I, J, NLEN
CHARACTER(60) :: CVAL, CROW, CCOL
INTEGER :: STENCIL(7) = 99999999, COLUMNS(7) = 99999999

if (NL /= NLEVEL_MIN) RETURN

#ifdef WITH_SCARC_DEBUG
MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G
   A => SCARC_POINT_TO_CMATRIX(G, NSCARC_MATRIX_POISSON)

   WRITE (CVAL, '(5A,i2.2,A,i2.2,A)') 'python/',TRIM(CHID),'_',TRIM(CNAME),'_m',NM,'_l',NL,'.val'
   WRITE (CCOL, '(5A,i2.2,A,i2.2,A)') 'python/',TRIM(CHID),'_',TRIM(CNAME),'_m',NM,'_l',NL,'.col'
   WRITE (CROW, '(5A,i2.2,A,i2.2,A)') 'python/',TRIM(CHID),'_',TRIM(CNAME),'_m',NM,'_l',NL,'.row'

   MVAL=GET_FILE_NUMBER()
   MCOL=GET_FILE_NUMBER()
   MROW=GET_FILE_NUMBER()

   OPEN(MVAL,FILE=CVAL)
   OPEN(MCOL,FILE=CCOL)
   OPEN(MROW,FILE=CROW)
   
   !WRITE(MVAL, *) '['
   !WRITE(MCOL, *) '['
   !WRITE(MROW, *) '['
   
   DO IC = 1, G%NC
      I = 1
      COLUMNS = 0
      STENCIL = NSCARC_HUGE_INT
      DO ICOL= A%ROW(IC), A%ROW(IC+1)-1
         COLUMNS(I) = ICOL
         STENCIL(I) = A%COL(ICOL)
         I = I + 1
      ENDDO
      NLEN = I - 1

#ifdef WITH_SCARC_DEBUG2
WRITE(MSG%LU_DEBUG,*) 'IC, NLEN, COLUMNS1, STENCIL1:', IC, NLEN, COLUMNS(1:NLEN), STENCIL(1:NLEN)
#endif

      DO I = NLEN, 2, -1
         DO J = 1, I-1
            IF (STENCIL(J) > STENCIL(J+1)) THEN
               JC = STENCIL(J+1)
               STENCIL(J+1) = STENCIL(J)
               STENCIL(J) = JC
               JC = COLUMNS(J+1)
               COLUMNS(J+1) = COLUMNS(J)
               COLUMNS(J) = JC
            ENDIF
         ENDDO
      ENDDO

#ifdef WITH_SCARC_DEBUG2
WRITE(MSG%LU_DEBUG,*) '          COLUMNS2, STENCIL2:', IC, NLEN, COLUMNS(1:NLEN), STENCIL(1:NLEN)
#endif

      DO I = 1, NLEN
         ICOL = COLUMNS(I)
         WRITE(MCOL,1001, ADVANCE="NO") A%COL(ICOL)-1
         WRITE(MVAL,1002, ADVANCE="NO") A%VAL(ICOL)
      ENDDO
      WRITE(MROW,1001, ADVANCE="NO") A%ROW(IC)-1

   ENDDO
   WRITE(MROW,1001, ADVANCE="NO") A%ROW(IC)-1

   !WRITE(MVAL, *) ' ]'
   !WRITE(MCOL, *) ' ]'
   !WRITE(MROW, *) ' ]'

   CLOSE(MVAL)
   CLOSE(MCOL)
   CLOSE(MROW)

ENDDO MESHES_LOOP
#endif

1001 FORMAT(I4,',')
1002 FORMAT(E10.2,',')
END SUBROUTINE SCARC_PYTHON_MATRIX


! ------------------------------------------------------------------------------------------------
! Print out matrix information on level NL for matlab
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_PYTHON_ZONES(NM, NL, CNAME)
USE SCARC_POINTERS, ONLY: G
INTEGER, INTENT(IN) :: NM, NL
CHARACTER(*), INTENT(IN) :: CNAME
INTEGER :: IC, MAGG
CHARACTER(60) :: CAGG

CALL SCARC_POINT_TO_GRID(NM, NL)

WRITE (CAGG, '(5A,i2.2,A,i2.2,A)') 'python/',TRIM(CHID),'_',TRIM(CNAME),'_m',NM,'_l',NL,'.val'
MAGG=GET_FILE_NUMBER()
OPEN(MAGG,FILE=CAGG)
!WRITE(MAGG, *) '['
DO IC = 1, G%N_FINE-1
   WRITE(MAGG,1001, ADVANCE="NO") G%ZONES_GLOBAL(IC)
ENDDO
WRITE(MAGG,1002, ADVANCE="NO") G%ZONES_GLOBAL(G%N_FINE)
!WRITE(MAGG, *) ' ]'
CLOSE(MAGG)

1001 FORMAT(I4,',')
1002 FORMAT(I4)
END SUBROUTINE SCARC_PYTHON_ZONES

#endif


! ------------------------------------------------------------------------------------------------
! Print out matrix information on level NL for matlab
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_BLENDER_ZONES(NM, NL)
USE SCARC_POINTERS, ONLY: M, L, OL, G, OG, F, XCOR, YCOR, ZCOR
INTEGER, INTENT(IN) :: NM, NL
REAL(EB) :: INCX, INCY, INCZ, XCOR0, YCOR0, ZCOR0
INTEGER :: IC, IZL, IZG, ICPT, MAGG, II, JJ, KK, IFACE, IOR0, NOM, ICG, INBR, ICW, ICE, ITYPE = 0
CHARACTER(60) :: CAGG

CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G

IF (NL == NLEVEL_MIN) THEN
   XCOR => M%X; YCOR => M%Y; ZCOR => M%Z
ELSE
   XCOR => L%XCOR; YCOR => L%YCOR; ZCOR => L%ZCOR
ENDIF

WRITE (CAGG, '(3A,i3.3,A,i3.3,A)') 'python/',TRIM(CHID),'/m',NM,'_l',NL,'.val'
MAGG=GET_FILE_NUMBER()
OPEN(MAGG,FILE=CAGG)

!CALL SCARC_ALLOCATE_INT1(VALUES, 1, L%NX, NSCARC_INIT_ZERO, 'VALUES')

IF (ITYPE == 0) THEN
   WRITE(MAGG,1001) G%NC, G%NC
ELSE
   WRITE(MAGG,1001) G%NC, G%NCE
ENDIF
WRITE(MAGG,1002) L%NX, L%NY, L%NZ
WRITE(MAGG,1003) L%DX, L%DY, L%DZ

IF (NL == NLEVEL_MIN) THEN
   DO KK = 1, L%NZ
      DO JJ = 1, L%NY
         DO II=1, L%NX
            IF (IS_UNSTRUCTURED.AND.L%IS_SOLID(II,JJ,KK)) THEN
               CYCLE
            ELSE
               IC=G%CELL_NUMBER(II,JJ,KK)
               IZL = ABS(G%ZONES_LOCAL(IC))
               IZG = ABS(G%ZONES_GLOBAL(IC))
               ICPT = G%CPOINTS(IZL)      
    
               XCOR0 = XCOR(II-1)
               ZCOR0 = ZCOR(KK-1)
   
               IF (TWO_D) THEN
                  YCOR0 = 0.0_EB
               ELSE
                  YCOR0 = YCOR(JJ-1)
               ENDIF
   
               IF (IC == ICPT) THEN
                  WRITE(MAGG,1000) IC, -IZG, XCOR0, YCOR0, ZCOR0
               ELSE
                  WRITE(MAGG,1000) IC,  IZG, XCOR0, YCOR0, ZCOR0
               ENDIF
   
            ENDIF
         ENDDO
      ENDDO
   ENDDO
ENDIF

! If ITYPE == 1 also plot overlapping information
IF (ITYPE == 1) THEN
   DO IFACE = 1, 6                                  
   
      IOR0 = FACE_ORIENTATION(IFACE)
      F => L%FACE(IOR0)
   
      INCX = 0
      INCY = 0
      INCZ = 0
   
      SELECT CASE(IOR0)
         CASE ( 1)
            INCX = -L%DX
         CASE (-1)
            INCX =  L%DX
         CASE ( 2)
            INCY = -L%DY
         CASE (-2)
            INCY =  L%DY
         CASE ( 3)
            INCZ = -L%DZ
         CASE (-3)
            INCZ =  L%DZ
      END SELECT
   
      DO INBR = 1, F%N_NEIGHBORS
   
         NOM = F%NEIGHBORS(INBR)
         CALL SCARC_POINT_TO_OTHER_GRID(NM, NOM, NL)
   
         DO ICG = OL%GHOST_FIRSTW(IOR0), OL%GHOST_LASTW(IOR0)
            ICW = OG%ICG_TO_ICW(ICG, 1)
            ICE = OG%ICG_TO_ICE(ICG, 1)
            II = G%ICX(ICW) 
            JJ = G%ICY(ICW) 
            KK = G%ICZ(ICW) 
            XCOR0 = XCOR(II-1) + INCX
            ZCOR0 = ZCOR(KK-1) + INCZ
            IF (TWO_D) THEN
               YCOR0 = 0.0_EB
            ELSE
               YCOR0 = YCOR(JJ-1) + INCY
            ENDIF
   WRITE(*,*) 'TODO: FIX BLENDER OUTPUT'
            IZG = ABS(G%ZONES_GLOBAL(ICE))
            WRITE(MAGG,1000) ICE, IZG, XCOR0, YCOR0, ZCOR0
         ENDDO
   
      ENDDO
   ENDDO
ENDIF

!DEALLOCATE(VALUES)
CLOSE(MAGG)

1001 FORMAT(I8,',', I8)
1002 FORMAT(I8,',', I8,',', I8)
1003 FORMAT(E12.4,',',  E12.4,',', E12.4)
1000 FORMAT(I8,',', I8,',', E12.4,',',  E12.4,',', E12.4)
END SUBROUTINE SCARC_BLENDER_ZONES

! ================================================================================================
! End  WITH_SCARC_DEBUG  - Part
! ================================================================================================

#ifdef WITH_SCARC_STANDALONE
! ================================================================================================
! Begin  WITH_SCARC_STANDALONE  - PART
! ================================================================================================
! ------------------------------------------------------------------------------------------------
! Dump matrix and vectors belonging to pressure system 
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_DUMP_SYSTEM (NSTACK, ITYPE)
USE SCARC_POINTERS, ONLY: SV, ST, G
INTEGER, INTENT(IN) :: NSTACK, ITYPE
INTEGER  :: NM, IC, JC, JCG, IP, IW, IOR0, N
INTEGER  :: COLUMNSL(7), COLUMNSG(7)
REAL(EB) :: VALUES(7), VAL

SV  => STACK(NSTACK)%SOLVER

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   ST => SCARC(NM)%LEVEL(NLEVEL_MIN)%STAGE(SV%TYPE_STAGE)
   CALL SCARC_POINT_TO_GRID (NM, NLEVEL_MIN)

   SELECT CASE(ITYPE)

      CASE(NSCARC_DUMP_MESH)          

         MSG%LU_CNN = GET_FILE_NUMBER()
         WRITE (MSG%FILE_CNN, '(A,I3.3,A,I5.5)') 'CNN/M',NM,'/Mesh.dat'
         OPEN (MSG%LU_CNN, FILE=MSG%FILE_CNN, FORM = 'FORMATTED')
         WRITE(MSG%LU_CNN,*) 'Total number of cells over all meshes'
         WRITE(MSG%LU_CNN,*) NC_GLOBAL(NLEVEL_MIN)
         WRITE(MSG%LU_CNN,*) 'Local number of cells in current mesh'
         WRITE(MSG%LU_CNN,*) G%NC_LOCAL(NM)
         WRITE(MSG%LU_CNN,*) 'Local number of cells in x, y and z in current mesh'
         WRITE(MSG%LU_CNN,*) L%NX, L%NY, L%NZ
         WRITE(MSG%LU_CNN,*) 'Section of local mesh in overall geometry with respect to global numbering'
         WRITE(MSG%LU_CNN,*) G%NC_OFFSET(NM)+1, G%NC_OFFSET(NM)+G%NC_LOCAL(NM)
         WRITE(MSG%LU_CNN,*) 'Coordinates of mesh x1, x2, y1, y2, z1, z2'
         WRITE(MSG%LU_CNN,*) S%XS, S%XF, S%YS, S%YF, S%ZS, S%ZF
         WRITE(MSG%LU_CNN,*) 'Local grid sizes in x, y and z in current mesh'
         WRITE(MSG%LU_CNN,*) L%DX, L%DY, L%DZ
         CLOSE(MSG%LU_CNN)

      CASE(NSCARC_DUMP_A)          

         A => G%POISSONC
         MSG%LU_CNN1 = GET_FILE_NUMBER() 
         MSG%LU_CNN2 = GET_FILE_NUMBER() 
         MSG%LU_CNN3 = GET_FILE_NUMBER() 
         WRITE (MSG%FILE_CNN1, '(A,I3.3,A)') 'CNN/M',NM,'/A/values.dat'
         WRITE (MSG%FILE_CNN2, '(A,I3.3,A)') 'CNN/M',NM,'/A/stencilLocal.dat'
         WRITE (MSG%FILE_CNN3, '(A,I3.3,A)') 'CNN/M',NM,'/A/stencilGlobal.dat'
         OPEN (NEWUNIT=MSG%LU_CNN1, FILE=MSG%FILE_CNN1, ACTION='READWRITE', FORM = 'FORMATTED')
         OPEN (NEWUNIT=MSG%LU_CNN2, FILE=MSG%FILE_CNN2, ACTION='READWRITE', FORM = 'FORMATTED')
         OPEN (NEWUNIT=MSG%LU_CNN3, FILE=MSG%FILE_CNN3, ACTION='READWRITE', FORM = 'FORMATTED')
         DO IC = 1, A%N_ROW - 1
            COLUMNSL = 0
            COLUMNSG = 0
            VALUES = 0.0_EB
            DO IP = A%ROW(IC), A%ROW(IC+1)-1
               JC  = A%COL(IP)
               JCG = G%LOCAL_TO_GLOBAL(JC)
               VAL = A%VAL(IP)
               N = -1
               IW = -1
               IOR0 = 0
               IF (TWO_D) THEN
                  IF (IC - JC > 1 .AND. JC <= G%NC)  THEN
                     N = 1
                  ELSE IF (IC - JC == 1)  THEN
                     N = 2
                  ELSE IF (IC - JC == 0)  THEN
                     N = 3
                  ELSE IF (IC - JC == -1)  THEN
                     N = 4
                  ELSE IF ( JC <= G%NC .AND. IC - JC < -1)  THEN
                     N = 5
                  ELSE IF (JC > G%NC) THEN
                     IW = G%ICE_TO_IWG(JC)
                     IOR0 = G%WALL(IW)%IOR
                     SELECT CASE (IOR0)
                        CASE (-3)
                           N = 5
                        CASE (-1)
                           N = 4
                        CASE ( 1)
                           N = 2
                        CASE ( 3)
                           N = 1
                     END SELECT
                  ENDIF
               ELSE
                  IF (IC - JC == 0)  THEN
                     N = 4
                  ELSE IF (IC - JC == 1)  THEN
                     N = 3
                  ELSE IF (IC - JC == -1)  THEN
                     N = 5
                  ELSE IF (1 < IC - JC .AND. IC - JC <= L%NX)  THEN
                     N = 2
                  ELSE IF (-L%NX <= IC - JC  .AND. IC - JC < -1 )  THEN
                     N = 6
                  ELSE IF (L%NX < IC - JC .AND. JC <= G%NC)  THEN
                     N = 1
                  ELSE IF (JC <= G%NC .AND. IC - JC <= -L%NX)  THEN
                     N = 7
                  ELSE IF (JC > G%NC) THEN
                     IW = G%ICE_TO_IWG(JC)
                     IOR0 = G%WALL(IW)%IOR
                     SELECT CASE (IOR0)
                        CASE (-3)
                           N = 7
                        CASE (-2)
                           N = 6
                        CASE (-1)
                           N = 5
                        CASE ( 1)
                           N = 3
                        CASE ( 2)
                           N = 2
                        CASE ( 3)
                           N = 1
                     END SELECT
                  ENDIF
               ENDIF
               IF (N > 0) THEN
                  VALUES(N) = VAL
                  COLUMNSL(N) = JC
                  COLUMNSG(N) = JCG
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'IP, IC, JC, IW, IOR0, VAL, N:', IP, IC, JC, JCG, IW, IOR0, VAL, N
#endif
               ENDIF
            ENDDO
            IF (TWO_D) THEN
               WRITE(MSG%LU_CNN1,*) VALUES(1:5)
               WRITE(MSG%LU_CNN2,*) COLUMNSL(1:5)
               WRITE(MSG%LU_CNN3,*) COLUMNSG(1:5)
            ELSE
               WRITE(MSG%LU_CNN1,*) VALUES(1:7)
               WRITE(MSG%LU_CNN2,*) COLUMNSL(1:7)
               WRITE(MSG%LU_CNN3,*) COLUMNSG(1:7)
#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,'(A,7I84,7E12.4)') 'COLUMNS, VALUES:', COLUMNSL(1:7), VALUES(1:7)
#endif
            ENDIF
         ENDDO
         CLOSE(MSG%LU_CNN1)
         CLOSE(MSG%LU_CNN2)
         CLOSE(MSG%LU_CNN3)

      CASE(NSCARC_DUMP_B)

         MSG%LU_CNN = GET_FILE_NUMBER() 
         WRITE (MSG%FILE_CNN, '(A,I3.3,A,I6.6,A)') 'CNN/M',NM,'/B/values_t',ITE_GLOBAL,'.dat'
         OPEN (MSG%LU_CNN, FILE=MSG%FILE_CNN, FORM = 'FORMATTED')
         DO IC = 1, G%NC
            WRITE(MSG%LU_CNN,*) ST%B(IC)
         ENDDO
         CLOSE(MSG%LU_CNN)

      CASE(NSCARC_DUMP_X)

         MSG%LU_CNN = GET_FILE_NUMBER() 
         WRITE (MSG%FILE_CNN, '(A,I3.3,A,I6.6,A)') 'CNN/M',NM,'/X/values_t',ITE_GLOBAL,'.dat'
         OPEN (MSG%LU_CNN, FILE=MSG%FILE_CNN, FORM = 'FORMATTED')
         DO IC = 1, G%NC
            WRITE(MSG%LU_CNN,*) ST%X(IC)
         ENDDO
         CLOSE(MSG%LU_CNN)

   END SELECT

ENDDO

END SUBROUTINE SCARC_DUMP_SYSTEM

! ------------------------------------------------------------------------------------------------
! Dump several arrays and structures needed for the standalone version of ScaRC 
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_DUMP_ENVIRONMENT
USE SCARC_POINTERS, ONLY: M
INTEGER :: NM, NOM, I, J, K, IW

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   M => MESHES(NM)

   ! Open SCARC-file for mesh NM
   WRITE (MSG%FILE_SCARC, '(A,A,A,i3.3)') 'dump/',TRIM(CHID),'_elk.dump',NM
   MSG%LU_SCARC = GET_FILE_NUMBER()
   OPEN (MSG%LU_SCARC, FILE=MSG%FILE_SCARC, ACTION = 'readwrite')

   ! Current process number
   WRITE(MSG%LU_SCARC,*) 'PROCESS(NM)'
   WRITE(MSG%LU_SCARC,*) PROCESS(NM)

   ! Dump information about all necessary MESHES-quantities
   WRITE(MSG%LU_SCARC,*) 'M%IBAR, M%JBAR, M%KBAR'
   WRITE(MSG%LU_SCARC,*) M%IBAR, M%JBAR, M%KBAR
   WRITE(MSG%LU_SCARC,*) 'M%IBP1, M%JBP1, M%KBP1'
   WRITE(MSG%LU_SCARC,*) M%IBP1, M%JBP1, M%KBP1
   WRITE(MSG%LU_SCARC,*) 'M%ITRN, M%JTRN, M%KTRN'
   WRITE(MSG%LU_SCARC,*) M%ITRN, M%JTRN, M%KTRN
   WRITE(MSG%LU_SCARC,*) 'M%LBC, M%MBC, M%NBC'
   WRITE(MSG%LU_SCARC,*) M%LBC, M%MBC, M%NBC
   WRITE(MSG%LU_SCARC,*) 'M%N_WALL_CELLS'
   WRITE(MSG%LU_SCARC,*) M%N_WALL_CELLS
   WRITE(MSG%LU_SCARC,*) 'M%N_EXTERNAL_WALL_CELLS'
   WRITE(MSG%LU_SCARC,*) M%N_EXTERNAL_WALL_CELLS
   WRITE(MSG%LU_SCARC,*) 'M%N_INTERNAL_WALL_CELLS'
   WRITE(MSG%LU_SCARC,*) M%N_INTERNAL_WALL_CELLS
   WRITE(MSG%LU_SCARC,*) 'M%IPS, M%LSAVE, M%LWORK'
   WRITE(MSG%LU_SCARC,*) M%IPS, M%LSAVE, M%LWORK
   WRITE(MSG%LU_SCARC,*) 'M%POIS_PTB'
   WRITE(MSG%LU_SCARC,*) M%POIS_PTB
   WRITE(MSG%LU_SCARC,*) 'M%XS, M%XF, M%YS, M%YF, M%ZS, M%ZF'
   WRITE(MSG%LU_SCARC,*) M%XS, M%XF, M%YS, M%YF, M%ZS, M%ZF
   WRITE(MSG%LU_SCARC,*) 'M%DXI, M%DETA, M%DZETA'
   WRITE(MSG%LU_SCARC,*) M%DXI, M%DETA, M%DZETA
   WRITE(MSG%LU_SCARC,*) 'M%RDXI, M%RDETA, M%RDZETA'
   WRITE(MSG%LU_SCARC,*) M%RDXI, M%RDETA, M%RDZETA
   WRITE(MSG%LU_SCARC,*) 'M%X'
   WRITE(MSG%LU_SCARC,*) (M%X(I), I=0, M%IBAR)
   WRITE(MSG%LU_SCARC,*) 'M%Y'
   WRITE(MSG%LU_SCARC,*) (M%Y(I), I=0, M%JBAR)
   WRITE(MSG%LU_SCARC,*) 'M%Z'
   WRITE(MSG%LU_SCARC,*) (M%Z(I), I=0, M%KBAR)
   WRITE(MSG%LU_SCARC,*) 'M%XC'
   WRITE(MSG%LU_SCARC,*) (M%XC(I), I=0, M%IBP1)
   WRITE(MSG%LU_SCARC,*) 'M%YC'
   WRITE(MSG%LU_SCARC,*) (M%YC(I), I=0, M%JBP1)
   WRITE(MSG%LU_SCARC,*) 'M%ZC'
   WRITE(MSG%LU_SCARC,*) (M%ZC(I), I=0, M%KBP1)
   WRITE(MSG%LU_SCARC,*) 'M%HX'
   WRITE(MSG%LU_SCARC,*) (M%HX(I), I=0, M%IBP1)
   WRITE(MSG%LU_SCARC,*) 'M%HY'
   WRITE(MSG%LU_SCARC,*) (M%HY(I), I=0, M%JBP1)
   WRITE(MSG%LU_SCARC,*) 'M%HZ'
   WRITE(MSG%LU_SCARC,*) (M%HZ(I), I=0, M%KBP1)

   DO NOM = 1, NMESHES
      WRITE(MSG%LU_SCARC,*) 'M%OMESH(',NOM,')%NIC_R, NICS'
      WRITE(MSG%LU_SCARC,*) M%OMESH(NOM)%NIC_R, M%OMESH(NOM)%NIC_S
   ENDDO

   IF (PREDICTOR) THEN
      WRITE(MSG%LU_SCARC,*) 'M%H'
      WRITE(MSG%LU_SCARC,*) (((M%H(I,J,K), I=0, M%IBP1), J=0, M%JBP1), K=0, M%KBP1)
   ELSE
      WRITE(MSG%LU_SCARC,*) 'M%HS'
      WRITE(MSG%LU_SCARC,*) (((M%HS(I,J,K), I=0, M%IBP1), J=0, M%JBP1), K=0, M%KBP1)
   ENDIF

   WRITE(MSG%LU_SCARC,*) 'M%PRHS'
   IF (TWO_D) THEN
      WRITE(MSG%LU_SCARC,*) ((M%PRHS(I,1,K), I=1, M%IBP1), K=1, M%KBP1)
   ELSE
      WRITE(MSG%LU_SCARC,*) (((M%PRHS(I,J,K), I=1, M%IBP1), J=1, M%JBP1), K=1, M%KBP1)
   ENDIF

   IF (TWO_D) THEN
      WRITE(MSG%LU_SCARC,*) 'M%BXS'
      WRITE(MSG%LU_SCARC,*) (M%BXS(1,K), K=1, M%KBP1)
      WRITE(MSG%LU_SCARC,*) 'M%BXF'
      WRITE(MSG%LU_SCARC,*) (M%BXF(1,K), K=1, M%KBP1)
   ELSE
      WRITE(MSG%LU_SCARC,*) 'M%BXS'
      WRITE(MSG%LU_SCARC,*) ((M%BXS(J,K), J=1, M%JBP1), K=1, M%KBP1)
      WRITE(MSG%LU_SCARC,*) 'M%BXF'
      WRITE(MSG%LU_SCARC,*) ((M%BXF(J,K), J=1, M%JBP1), K=1, M%KBP1)
   ENDIF

   WRITE(MSG%LU_SCARC,*) 'M%BYS'
   WRITE(MSG%LU_SCARC,*) ((M%BYS(I,K), I=1, M%IBP1), K=1, M%KBP1)
   WRITE(MSG%LU_SCARC,*) 'M%BYF'
   WRITE(MSG%LU_SCARC,*) ((M%BYF(I,K), I=1, M%IBP1), K=1, M%KBP1)

   IF (TWO_D) THEN
      WRITE(MSG%LU_SCARC,*) 'M%BZS'
      WRITE(MSG%LU_SCARC,*) (M%BZS(I,1), I=1, M%IBP1)
      WRITE(MSG%LU_SCARC,*) 'M%BZF'
      WRITE(MSG%LU_SCARC,*) (M%BZF(I,1), I=1, M%IBP1)
   ELSE
      WRITE(MSG%LU_SCARC,*) 'M%BZS'
      WRITE(MSG%LU_SCARC,*) ((M%BZS(I,J), I=1, M%IBP1), J=1, M%JBP1)
      WRITE(MSG%LU_SCARC,*) 'M%BZF'
      WRITE(MSG%LU_SCARC,*) ((M%BZF(I,J), I=1, M%IBP1), J=1, M%JBP1)
   ENDIF

   WRITE(MSG%LU_SCARC,*) 'M%WALL(.)%PRESSURE_BC_INDEX'
   WRITE(MSG%LU_SCARC,*) (M%WALL(IW)%PRESSURE_BC_INDEX, IW=1, M%N_WALL_CELLS)

   WRITE(MSG%LU_SCARC,*) 'M%WALL(.)%BOUNDARY_TYPE'
   WRITE(MSG%LU_SCARC,*) (M%WALL(IW)%BOUNDARY_TYPE, IW=1, M%N_WALL_CELLS)

   WRITE(MSG%LU_SCARC,*) 'M%WALL(.)%WALL_INDEX'
   WRITE(MSG%LU_SCARC,*) (M%WALL(IW)%WALL_INDEX, IW=1, M%N_WALL_CELLS)

   WRITE(MSG%LU_SCARC,*) 'M%WALL(.)%SURF_INDEX'
   WRITE(MSG%LU_SCARC,*) (M%WALL(IW)%SURF_INDEX, IW=1, M%N_WALL_CELLS)

   WRITE(MSG%LU_SCARC,*) 'M%WALL(.)%BACK_INDEX'
   WRITE(MSG%LU_SCARC,*) (M%WALL(IW)%BACK_INDEX, IW=1, M%N_WALL_CELLS)

   WRITE(MSG%LU_SCARC,*) 'M%WALL(.)%OBST_INDEX'
   WRITE(MSG%LU_SCARC,*) (M%WALL(IW)%OBST_INDEX, IW=1, M%N_WALL_CELLS)

   WRITE(MSG%LU_SCARC,*) 'M%WALL(.)%VENT_INDEX'
   WRITE(MSG%LU_SCARC,*) (M%WALL(IW)%VENT_INDEX, IW=1, M%N_WALL_CELLS)

   WRITE(MSG%LU_SCARC,*) 'M%WALL(.)%ONE_D%IOR'
   WRITE(MSG%LU_SCARC,*) (M%WALL(IW)%ONE_D%IOR, IW=1, M%N_EXTERNAL_WALL_CELLS)

   WRITE(MSG%LU_SCARC,*) 'M%WALL(.)%ONE_D%II'
   WRITE(MSG%LU_SCARC,*) (M%WALL(IW)%ONE_D%II, IW=1, M%N_EXTERNAL_WALL_CELLS)

   WRITE(MSG%LU_SCARC,*) 'M%WALL(.)%ONE_D%JJ'
   WRITE(MSG%LU_SCARC,*) (M%WALL(IW)%ONE_D%JJ, IW=1, M%N_EXTERNAL_WALL_CELLS)

   WRITE(MSG%LU_SCARC,*) 'M%WALL(.)%ONE_D%KK'
   WRITE(MSG%LU_SCARC,*) (M%WALL(IW)%ONE_D%KK, IW=1, M%N_EXTERNAL_WALL_CELLS)

   WRITE(MSG%LU_SCARC,*) 'M%WALL(.)%ONE_D%IIG'
   WRITE(MSG%LU_SCARC,*) (M%WALL(IW)%ONE_D%IIG, IW=1, M%N_EXTERNAL_WALL_CELLS)

   WRITE(MSG%LU_SCARC,*) 'M%WALL(.)%ONE_D%JJG'
   WRITE(MSG%LU_SCARC,*) (M%WALL(IW)%ONE_D%JJG, IW=1, M%N_EXTERNAL_WALL_CELLS)

   WRITE(MSG%LU_SCARC,*) 'M%WALL(.)%ONE_D%KKG'
   WRITE(MSG%LU_SCARC,*) (M%WALL(IW)%ONE_D%KKG, IW=1, M%N_EXTERNAL_WALL_CELLS)

   WRITE(MSG%LU_SCARC,*) 'M%EXTERNAL_WALL(.)%NOM'
   WRITE(MSG%LU_SCARC,*) (M%EXTERNAL_WALL(IW)%NOM, IW=1, M%N_EXTERNAL_WALL_CELLS)

   WRITE(MSG%LU_SCARC,*) 'M%EXTERNAL_WALL(.)%IIO_MIN'
   WRITE(MSG%LU_SCARC,*) (M%EXTERNAL_WALL(IW)%IIO_MIN, IW=1, M%N_EXTERNAL_WALL_CELLS)

   WRITE(MSG%LU_SCARC,*) 'M%EXTERNAL_WALL(.)%IIO_MAX'
   WRITE(MSG%LU_SCARC,*) (M%EXTERNAL_WALL(IW)%IIO_MAX, IW=1, M%N_EXTERNAL_WALL_CELLS)

   WRITE(MSG%LU_SCARC,*) 'M%EXTERNAL_WALL(.)%JJO_MIN'
   WRITE(MSG%LU_SCARC,*) (M%EXTERNAL_WALL(IW)%JJO_MIN, IW=1, M%N_EXTERNAL_WALL_CELLS)

   WRITE(MSG%LU_SCARC,*) 'M%EXTERNAL_WALL(.)%JJO_MAX'
   WRITE(MSG%LU_SCARC,*) (M%EXTERNAL_WALL(IW)%JJO_MAX, IW=1, M%N_EXTERNAL_WALL_CELLS)

   WRITE(MSG%LU_SCARC,*) 'M%EXTERNAL_WALL(.)%KKO_MIN'
   WRITE(MSG%LU_SCARC,*) (M%EXTERNAL_WALL(IW)%KKO_MIN, IW=1, M%N_EXTERNAL_WALL_CELLS)

   WRITE(MSG%LU_SCARC,*) 'M%EXTERNAL_WALL(.)%KKO_MAX'
   WRITE(MSG%LU_SCARC,*) (M%EXTERNAL_WALL(IW)%KKO_MAX, IW=1, M%N_EXTERNAL_WALL_CELLS)

   WRITE(MSG%LU_SCARC,*) 'M%N_OBST'
   WRITE(MSG%LU_SCARC,*) M%N_OBST
   IF (M%N_OBST > 0) THEN
      WRITE(MSG%LU_SCARC,*) 'M%OBSTRUCTION(.)%I1, I2, J1, J2, K1, K2'
      DO I = 1, M%N_OBST
         WRITE(MSG%LU_SCARC,*) M%OBSTRUCTION(I)%I1, M%OBSTRUCTION(I)%I2, &
                               M%OBSTRUCTION(I)%J1, M%OBSTRUCTION(I)%J2, &
                               M%OBSTRUCTION(I)%K1, M%OBSTRUCTION(I)%K2
      ENDDO
   ENDIF

   WRITE(MSG%LU_SCARC,*) 'CELL_COUNT(',NM,')'
   WRITE(MSG%LU_SCARC,*) CELL_COUNT(NM)

   WRITE(MSG%LU_SCARC,*) 'M%SOLID'
   WRITE(MSG%LU_SCARC,*) (M%SOLID(I), I=1, CELL_COUNT(NM))

   WRITE(MSG%LU_SCARC,*) 'M%CELL_INDEX'
   WRITE(MSG%LU_SCARC,*) (((M%CELL_INDEX(I,J,K), I=0, M%IBAR+1), J=0, M%JBAR+1), K=0, M%KBAR+1)

   WRITE(MSG%LU_SCARC,*) 'M%WALL_INDEX'
   DO I = 0, CELL_COUNT(NM)
      WRITE(MSG%LU_SCARC,'(7I8)') (M%WALL_INDEX(I,J), J=-3,3)
   ENDDO

   WRITE(MSG%LU_SCARC,*) 'M%SAVE1'
   WRITE(MSG%LU_SCARC,*) (M%SAVE1(I), I=-3, M%LSAVE)

   WRITE(MSG%LU_SCARC,*) 'M%WORK'
   WRITE(MSG%LU_SCARC,*) (M%WORK(I), I=1, M%LWORK)

   CLOSE (MSG%LU_SCARC)

ENDDO

SCARC_DUMP = .FALSE.
END SUBROUTINE SCARC_DUMP_ENVIRONMENT


! ------------------------------------------------------------------------------------------------
! Dump several arrays and structures needed for ScaRC - standalone version
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_RESTORE_ENVIRONMENT
USE SCARC_POINTERS, ONLY: M
INTEGER :: NM, NOM, I, J, K, IW, LU, IERR=0

!CPATH = 'D:\GIT\HHP\A_ScaRC\VisualStudio\Cases\'
!CPATH = '../VisualStudio/Cases/'

ALLOCATE(CELL_COUNT(NMESHES), STAT= IERR)
CALL CHKMEMERR ('SCARC_RESTORE_ENVIRONMENT', 'CELL_COUNT', IERR)

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   M => MESHES(NM)

   MSG%LU_SCARC = 99
   WRITE (MSG%FILE_SCARC, '(A,A,A,i3.3)') 'dump/',TRIM(CHID),'.dump',NM
   write(*,*) MSG%FILE_SCARC
   MSG%FILE_SCARC = TRIM(MSG%FILE_SCARC)

   LU = MSG%LU_SCARC
   OPEN (LU, FILE=MSG%FILE_SCARC)

   ! Current process number
   READ(LU,*)
   READ(LU,*) PROCESS(NM)

   ! Dump information about all necessary MESHES-quantities
   READ(LU,*)
   READ(LU,*) M%IBAR, M%JBAR, M%KBAR
   READ(LU,*) 
   READ(LU,*) M%IBP1, M%JBP1, M%KBP1
   READ(LU,*) 
   READ(LU,*) M%ITRN, M%JTRN, M%KTRN
   READ(LU,*) 
   READ(LU,*) M%LBC, M%MBC, M%NBC
   READ(LU,*) 
   READ(LU,*) M%N_WALL_CELLS
   READ(LU,*)
   READ(LU,*) M%N_EXTERNAL_WALL_CELLS
   READ(LU,*)
   READ(LU,*) M%N_INTERNAL_WALL_CELLS
   READ(LU,*)
   READ(LU,*) M%IPS, M%LSAVE, M%LWORK
   READ(LU,*) 
   READ(LU,*) M%POIS_PTB
   READ(LU,*)
   READ(LU,*) M%XS, M%XF, M%YS, M%YF, M%ZS, M%ZF
   READ(LU,*) 
   READ(LU,*) M%DXI, M%DETA, M%DZETA
   READ(LU,*) 
   READ(LU,*) M%RDXI, M%RDETA, M%RDZETA

   ALLOCATE(M%X(0:M%IBAR), STAT= IERR)
   CALL CHKMEMERR ('SCARC_RESTORE_ENVIRONMENT', 'M%X', IERR)

   ALLOCATE(M%Y(0:M%JBAR), STAT= IERR)
   CALL CHKMEMERR ('SCARC_RESTORE_ENVIRONMENT', 'M%Y', IERR)

   ALLOCATE(M%Z(0:M%KBAR), STAT= IERR)
   CALL CHKMEMERR ('SCARC_RESTORE_ENVIRONMENT', 'M%Z', IERR)

   ALLOCATE(M%XC(0:M%IBP1), STAT= IERR)
   CALL CHKMEMERR ('SCARC_RESTORE_ENVIRONMENT', 'M%XC', IERR)

   ALLOCATE(M%YC(0:M%JBP1), STAT= IERR)
   CALL CHKMEMERR ('SCARC_RESTORE_ENVIRONMENT', 'M%YC', IERR)

   ALLOCATE(M%ZC(0:M%KBP1), STAT= IERR)
   CALL CHKMEMERR ('SCARC_RESTORE_ENVIRONMENT', 'M%ZC', IERR)

   ALLOCATE(M%HX(0:M%IBP1), STAT= IERR)
   CALL CHKMEMERR ('SCARC_RESTORE_ENVIRONMENT', 'M%HX', IERR)

   ALLOCATE(M%HY(0:M%JBP1), STAT= IERR)
   CALL CHKMEMERR ('SCARC_RESTORE_ENVIRONMENT', 'M%HY', IERR)

   ALLOCATE(M%HZ(0:M%KBP1), STAT= IERR)
   CALL CHKMEMERR ('SCARC_RESTORE_ENVIRONMENT', 'M%HZ', IERR)

   READ(LU,*)
   READ(LU,*) (M%X(I), I=0, M%IBAR)
   READ(LU,*)
   READ(LU,*) (M%Y(I), I=0, M%JBAR)
   READ(LU,*)
   READ(LU,*) (M%Z(I), I=0, M%KBAR)
   READ(LU,*)
   READ(LU,*) (M%XC(I), I=0, M%IBP1)
   READ(LU,*)
   READ(LU,*) (M%YC(I), I=0, M%JBP1)
   READ(LU,*)
   READ(LU,*) (M%ZC(I), I=0, M%KBP1)

   READ(LU,*) 
   READ(LU,*) (M%HX(I), I=0, M%IBP1)
   READ(LU,*) 
   READ(LU,*) (M%HY(I), I=0, M%JBP1)
   READ(LU,*) 
   READ(LU,*) (M%HZ(I), I=0, M%KBP1)

   ALLOCATE(M%OMESH(NMESHES), STAT= IERR)
   CALL CHKMEMERR ('SCARC_RESTORE_ENVIRONMENT', 'M%OMESH', IERR)

   DO NOM = 1, NMESHES
      READ(LU,*)
      READ(LU,*) M%OMESH(NOM)%NIC_R, M%OMESH(NOM)%NIC_S
   ENDDO

   IF (PREDICTOR) THEN

      ALLOCATE(M%H(0:M%IBP1, 0:M%JBP1, 0:M%KBP1), STAT= IERR)
      CALL CHKMEMERR ('SCARC_RESTORE_ENVIRONMENT', 'M%H', IERR)
      M%H = 0.0_EB

      READ(LU,*)
      READ(LU,*) (((M%H(I,J,K), I=0, M%IBP1), J=0, M%JBP1), K=0, M%KBP1)

   ELSE

      ALLOCATE(M%HS(0:M%IBP1, 0:M%JBP1, 0:M%KBP1), STAT= IERR)
      CALL CHKMEMERR ('SCARC_RESTORE_ENVIRONMENT', 'M%H', IERR)
      M%HS = 0.0_EB

      READ(LU,*)
      READ(LU,*) (((M%HS(I,J,K), I=0, M%IBP1), J=0, M%JBP1), K=0, M%KBP1)

   ENDIF
   IF (TWO_D) THEN

      ALLOCATE(M%PRHS(M%IBP1, 1, M%KBP1), STAT= IERR)
      CALL CHKMEMERR ('SCARC_RESTORE_ENVIRONMENT', 'M%H', IERR)

      READ(LU,*)
      READ(LU,*) ((M%PRHS(I,1,K), I=1, M%IBP1), K=1, M%KBP1)

   ELSE

      ALLOCATE(M%PRHS(M%IBP1, M%JBP1, M%KBP1), STAT= IERR)
      CALL CHKMEMERR ('SCARC_RESTORE_ENVIRONMENT', 'M%H', IERR)

      READ(LU,*)
      READ(LU,*) (((M%PRHS(I,J,K), I=1, M%IBP1), J=1, M%JBP1), K=1, M%KBP1)

   ENDIF

   IF (TWO_D) THEN

      ALLOCATE(M%BXS(1,M%KBP1), STAT= IERR)
      CALL CHKMEMERR ('READ_SCARC', 'M%BXS', IERR)

      ALLOCATE(M%BXF(1,M%KBP1), STAT= IERR)
      CALL CHKMEMERR ('READ_SCARC', 'M%BXF', IERR)

      ALLOCATE(M%BZS(M%IBP1,1), STAT= IERR)
      CALL CHKMEMERR ('READ_SCARC', 'M%BZS', IERR)

      ALLOCATE(M%BZF(M%IBP1,1), STAT= IERR)
      CALL CHKMEMERR ('READ_SCARC', 'M%BZF', IERR)

   ELSE

      ALLOCATE(M%BXS(M%JBP1,M%KBP1), STAT= IERR)
      CALL CHKMEMERR ('READ_SCARC', 'M%BXS', IERR)

      ALLOCATE(M%BXF(M%JBP1,M%KBP1), STAT= IERR)
      CALL CHKMEMERR ('READ_SCARC', 'M%BXF', IERR)

      ALLOCATE(M%BZS(M%IBP1,M%JBP1), STAT= IERR)
      CALL CHKMEMERR ('READ_SCARC', 'M%BZS', IERR)

      ALLOCATE(M%BZF(M%IBP1,M%JBP1), STAT= IERR)
      CALL CHKMEMERR ('READ_SCARC', 'M%BZF', IERR)

   ENDIF
   
   ALLOCATE(M%BYS(M%IBP1,M%KBP1), STAT= IERR)
   CALL CHKMEMERR ('READ_SCARC', 'M%BYS', IERR)

   ALLOCATE(M%BYF(M%IBP1,M%KBP1), STAT= IERR)
   CALL CHKMEMERR ('READ_SCARC', 'M%BYF', IERR)

   IF (TWO_D) THEN
      READ(LU,*)
      READ(LU,*) (M%BXS(1,K), K=1, M%KBP1)
      READ(LU,*)
      READ(LU,*) (M%BXF(1,K), K=1, M%KBP1)
   ELSE
      READ(LU,*)
      READ(LU,*) ((M%BXS(J,K), J=1, M%JBP1), K=1, M%KBP1)
      READ(LU,*)
      READ(LU,*) ((M%BXF(J,K), J=1, M%JBP1), K=1, M%KBP1)
   ENDIF
   
   READ(LU,*)
   READ(LU,*) ((M%BYS(I,K), I=1, M%IBP1), K=1, M%KBP1)
   READ(LU,*)
   READ(LU,*) ((M%BYF(I,K), I=1, M%IBP1), K=1, M%KBP1)

   IF (TWO_D) THEN
      READ(LU,*)
      READ(LU,*) (M%BZS(I,1), I=1, M%IBP1)
      READ(LU,*)
      READ(LU,*) (M%BZF(I,1), I=1, M%IBP1)
   ELSE
      READ(LU,*)
      READ(LU,*) ((M%BZS(I,J), I=1, M%IBP1), J=1, M%JBP1)
      READ(LU,*)
      READ(LU,*) ((M%BZF(I,J), I=1, M%IBP1), J=1, M%JBP1)
   ENDIF

   ALLOCATE(M%WALL(M%N_WALL_CELLS), STAT= IERR)
   CALL CHKMEMERR ('READ_SCARC', 'M%WALL', IERR)

   ALLOCATE(M%EXTERNAL_WALL(M%N_EXTERNAL_WALL_CELLS), STAT= IERR)
   CALL CHKMEMERR ('READ_SCARC', 'M%EXTERNAL_WALL', IERR)

   READ(LU,*)
   READ(LU,*) (M%WALL(IW)%PRESSURE_BC_INDEX, IW=1, M%N_WALL_CELLS)

   READ(LU,*)
   READ(LU,*) (M%WALL(IW)%BOUNDARY_TYPE, IW=1, M%N_WALL_CELLS)

   READ(LU,*)
   READ(LU,*) (M%WALL(IW)%WALL_INDEX, IW=1, M%N_WALL_CELLS)

   READ(LU,*)
   READ(LU,*) (M%WALL(IW)%SURF_INDEX, IW=1, M%N_WALL_CELLS)

   READ(LU,*)
   READ(LU,*) (M%WALL(IW)%BACK_INDEX, IW=1, M%N_WALL_CELLS)

   READ(LU,*)
   READ(LU,*) (M%WALL(IW)%OBST_INDEX, IW=1, M%N_WALL_CELLS)

   READ(LU,*)
   READ(LU,*) (M%WALL(IW)%VENT_INDEX, IW=1, M%N_WALL_CELLS)

   READ(LU,*)
   READ(LU,*) (M%WALL(IW)%ONE_D%IOR, IW=1, M%N_EXTERNAL_WALL_CELLS)

   READ(LU,*)
   READ(LU,*) (M%WALL(IW)%ONE_D%II, IW=1, M%N_EXTERNAL_WALL_CELLS)

   READ(LU,*)
   READ(LU,*) (M%WALL(IW)%ONE_D%JJ, IW=1, M%N_EXTERNAL_WALL_CELLS)

   READ(LU,*)
   READ(LU,*) (M%WALL(IW)%ONE_D%KK, IW=1, M%N_EXTERNAL_WALL_CELLS)

   READ(LU,*)
   READ(LU,*) (M%WALL(IW)%ONE_D%IIG, IW=1, M%N_EXTERNAL_WALL_CELLS)

   READ(LU,*)
   READ(LU,*) (M%WALL(IW)%ONE_D%JJG, IW=1, M%N_EXTERNAL_WALL_CELLS)

   READ(LU,*)
   READ(LU,*) (M%WALL(IW)%ONE_D%KKG, IW=1, M%N_EXTERNAL_WALL_CELLS)

   READ(LU,*) 
   READ(LU,*) (M%EXTERNAL_WALL(IW)%NOM, IW=1, M%N_EXTERNAL_WALL_CELLS)

   READ(LU,*)
   READ(LU,*) (M%EXTERNAL_WALL(IW)%IIO_MIN, IW=1, M%N_EXTERNAL_WALL_CELLS)

   READ(LU,*)
   READ(LU,*) (M%EXTERNAL_WALL(IW)%IIO_MAX, IW=1, M%N_EXTERNAL_WALL_CELLS)

   READ(LU,*) 
   READ(LU,*) (M%EXTERNAL_WALL(IW)%JJO_MIN, IW=1, M%N_EXTERNAL_WALL_CELLS)

   READ(LU,*) 
   READ(LU,*) (M%EXTERNAL_WALL(IW)%JJO_MAX, IW=1, M%N_EXTERNAL_WALL_CELLS)

   READ(LU,*) 
   READ(LU,*) (M%EXTERNAL_WALL(IW)%KKO_MIN, IW=1, M%N_EXTERNAL_WALL_CELLS)

   READ(LU,*)
   READ(LU,*) (M%EXTERNAL_WALL(IW)%KKO_MAX, IW=1, M%N_EXTERNAL_WALL_CELLS)

   READ(LU,*)
   READ(LU,*) M%N_OBST

   IF (M%N_OBST > 0) THEN

      ALLOCATE(M%OBSTRUCTION(M%N_OBST), STAT= IERR)
      CALL CHKMEMERR ('READ_SCARC', 'M%OBSTRUCTION', IERR)

      READ(LU,*)
      DO I = 1, M%N_OBST
         READ(LU,*) M%OBSTRUCTION(I)%I1, M%OBSTRUCTION(I)%I2, &
                    M%OBSTRUCTION(I)%J1, M%OBSTRUCTION(I)%J2, &
                    M%OBSTRUCTION(I)%K1, M%OBSTRUCTION(I)%K2
      ENDDO
   ENDIF

   READ(LU,*)
   READ(LU,*) CELL_COUNT(NM)

   ALLOCATE(M%SOLID(0:CELL_COUNT(NM)), STAT= IERR)
   CALL CHKMEMERR ('READ_SCARC', 'M%SOLID', IERR)
   M%SOLID=.FALSE.

   ALLOCATE(M%CELL_INDEX(0:M%IBP1,0:M%JBP1,0:M%KBP1), STAT= IERR)
   CALL CHKMEMERR ('READ_SCARC', 'M%CELL_INDEX', IERR)
   M%CELL_INDEX = 0

   ALLOCATE(M%WALL_INDEX(0:CELL_COUNT(NM),-3:3), STAT= IERR)
   CALL CHKMEMERR ('READ_SCARC', 'M%WALL_INDEX', IERR)

   READ(LU,*)
   READ(LU,*) (M%SOLID(I), I=1, CELL_COUNT(NM))

   READ(LU,*)
   READ(LU,*) (((M%CELL_INDEX(I,J,K), I=0, M%IBP1), J=0, M%JBP1), K=0, M%KBP1)

   READ(LU,*)
   DO I = 0, CELL_COUNT(NM)
      READ(LU,'(7I8)') (M%WALL_INDEX(I,J), J=-3,3)
   ENDDO

   ALLOCATE(M%SAVE1(-3:M%LSAVE), STAT= IERR)
   CALL CHKMEMERR ('SCARC_RESTORE_ENVIRONMENT', 'M%SAVE1', IERR)

   ALLOCATE(M%WORK(M%LWORK), STAT= IERR)
   CALL CHKMEMERR ('SCARC_RESTORE_ENVIRONMENT', 'M%WORK', IERR)

   READ(LU,*) 
   READ(LU,*) (M%SAVE1(I), I=-3, M%LSAVE)

   READ(LU,*) 
   READ(LU,*) (M%WORK(I), I=1, M%LWORK)

   CLOSE (LU)

ENDDO

END SUBROUTINE SCARC_RESTORE_ENVIRONMENT


! ----------------------------------------------------------------------------------------------------
! Allocate and initialize vectors pressure diagnostics - only for developping purposes
! ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_PRESSURE()
USE SCARC_POINTERS, ONLY: L, G, PRES
INTEGER :: NM

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID (NM, NLEVEL_MIN)

   CALL SCARC_ALLOCATE_REAL1(PRES%B_OLD, 1, G%NC, NSCARC_INIT_ZERO, 'PRES%B_OLD')
   CALL SCARC_ALLOCATE_REAL1(PRES%B_NEW, 1, G%NC, NSCARC_INIT_ZERO, 'PRES%B_NEW')

   CALL SCARC_ALLOCATE_REAL3(PRES%H_OLD, 0, L%NX+1, 0, L%NY+1, 0, L%NZ+1, NSCARC_INIT_ZERO, 'PRES%H_OLD')
   CALL SCARC_ALLOCATE_REAL3(PRES%H_NEW, 0, L%NX+1, 0, L%NY+1, 0, L%NZ+1, NSCARC_INIT_ZERO, 'PRES%H_NEW')

   CALL SCARC_ALLOCATE_REAL3(PRES%HS_OLD, 0, L%NX+1, 0, L%NY+1, 0, L%NZ+1, NSCARC_INIT_ZERO, 'PRES%HS_OLD')
   CALL SCARC_ALLOCATE_REAL3(PRES%HS_NEW, 0, L%NX+1, 0, L%NY+1, 0, L%NZ+1, NSCARC_INIT_ZERO, 'PRES%HS_NEW')

   !CALL SCARC_ALLOCATE_REAL3(PRES%PRHS, 1, L%NX+1, 1, L%NY+1, 1, L%NZ+1, NSCARC_INIT_ZERO, 'PRES%PRHS')

ENDDO

END SUBROUTINE SCARC_SETUP_PRESSURE


! ------------------------------------------------------------------------------------------------
! Compute Differences between old and new pressure solutions - only for developping purposes
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_PRESSURE_DIFFERENCE(NL)
USE SCARC_POINTERS, ONLY: L, G, PRES
INTEGER, INTENT(IN) :: NL
INTEGER :: NM, IX, IY, IZ

! -----------------------------------------------------------------------------------------------
! Store new pressure vector for comparison in next time step
! -----------------------------------------------------------------------------------------------
DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_GRID (NM, NL)                                   ! Sets grid pointer G

   IF (PREDICTOR) THEN
   
      PRES%H_NEW  = M%H
      PRES%DIFF_H = 0.0_EB
      DO IZ = 1, L%NZ
         DO IY = 1, L%NY
            DO IX = 1, L%NX
               PRES%DIFF_H = PRES%DIFF_H + (PRES%H_NEW(IX, IY, IZ) - PRES%H_OLD(IX, IY, IZ))**2         
            ENDDO
         ENDDO
      ENDDO
      PRES%DIFF_H = PRES%DIFF_H / REAL(L%N_CELLS, EB)

   ELSE

      PRES%HS_NEW = M%HS
      PRES%DIFF_HS = 0.0_EB
      DO IZ = 1, L%NZ
         DO IY = 1, L%NY
            DO IX = 1, L%NX
               PRES%DIFF_HS = PRES%DIFF_HS + (PRES%HS_NEW(IX, IY, IZ) - PRES%HS_OLD(IX, IY, IZ))**2         
            ENDDO
         ENDDO
      ENDDO
      PRES%DIFF_HS = PRES%DIFF_HS / REAL(L%N_CELLS, EB)

   ENDIF

ENDDO

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'Differences of pressure vectors on mesh ', NM,' : ', PRES%DIFF_H, PRES%DIFF_HS
#endif

END SUBROUTINE SCARC_PRESSURE_DIFFERENCE

! ================================================================================================
! End  WITH_SCARC_STANDALONE  - PART
! ================================================================================================
#endif

END MODULE SCRC




