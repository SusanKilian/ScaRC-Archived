!> --------------------------------------------------------------------------------------------------------
!> Use of different directives possible
!>     - WITH_SCARC_VERBOSE : print more detailed information about ScaRC iterations
!>     - WITH_SCARC_DEBUG   : print detaild debugging info (only for code development)
!>     - WITH_SCARC_CGBARO  : test modified CG-routine which considers baroclinic effect 
!>     - WITH_MKL           : use MKL routines PARDISO, CLUSTER_SPARSE_SOLVER, DDOT, DAXPBY, DCOPY, DSCAL
!>     - WITH_MKL_FB        : perform MKL routines for LU-decomposition only in single precision
!> --------------------------------------------------------------------------------------------------------
#define WITH_SCARC_VERBOSE
#define WITH_SCARC_DEBUG
#undef WITH_SCARC_CGBARO
#undef WITH_MKL_FB

MODULE SCRC

USE PRECISION_PARAMETERS
USE GLOBAL_CONSTANTS
USE MESH_VARIABLES
USE MEMORY_FUNCTIONS, ONLY: CHKMEMERR
USE COMP_FUNCTIONS, ONLY: CURRENT_TIME, GET_FILE_NUMBER, SHUTDOWN
USE MPI

#ifdef WITH_MKL
  USE MKL_PARDISO
  USE MKL_CLUSTER_SPARSE_SOLVER
#endif

IMPLICIT NONE

!> ------------------------------------------------------------------------------------------------
!> Public subroutines (initialization, solver, time measurement and revisioning)
!> ------------------------------------------------------------------------------------------------
PUBLIC SCARC_SETUP                             !> Setup routine for ScaRC, needed in main.f90
PUBLIC SCARC_SOLVER                            !> Call of basic ScaRC solver
PUBLIC SCARC_TIMINGS                           !> Call of time measurements for ScaRC

!> ------------------------------------------------------------------------------------------------
!> Public variables   (explanations in declaration part below)
!> Note: For input parameters in character format corresponding INTEGER type-parameters will
!> Be introduced later to simplify inquiries
!> ------------------------------------------------------------------------------------------------
PUBLIC SCARC_METHOD                            !> ScaRC method
PUBLIC SCARC_DISCRETIZATION                    !> ScaRC discretization type
PUBLIC SCARC_MATRIX_FORMAT                     !> Storage format for matrices
PUBLIC SCARC_RESIDUAL                          !> Residual of iterative solver
PUBLIC SCARC_ITERATIONS                        !> Number of iterations
PUBLIC SCARC_CAPPA                             !> Convergence rate
PUBLIC SCARC_ACCURACY                          !> Chosen accuracy type (relative/absolute)
PUBLIC SCARC_PRECISION                         !> Single/double precision for preconditioning or LU-decomposition
PUBLIC SCARC_TWOLEVEL                          !> Type of Twolevel-CG-method
PUBLIC SCARC_ERROR_FILE                        !> Print ScaRC error file

PUBLIC SCARC_KRYLOV                            !> Type of Krylov method
PUBLIC SCARC_KRYLOV_ITERATIONS                 !> Maximum number of iterations for Krylov method
PUBLIC SCARC_KRYLOV_INTERPOL                   !> Interpolation type for twolevel methods
PUBLIC SCARC_KRYLOV_ACCURACY                   !> Requested accuracy for Krylov method

PUBLIC SCARC_PRECON                            !> Preconditioner of Krylov method
PUBLIC SCARC_PRECON_SCOPE                      !> Scope of preconditioner (local/global)
PUBLIC SCARC_PRECON_ITERATIONS                 !> Maximum number of iterations for perconditioning method
PUBLIC SCARC_PRECON_ACCURACY                   !> Requested accuracy for preconditioning method
PUBLIC SCARC_PRECON_OMEGA                      !> Relaxation parameter for preconditioning method

PUBLIC SCARC_MULTIGRID                         !> Type of multigrid method
PUBLIC SCARC_MULTIGRID_LEVEL                   !> Multigrid level
PUBLIC SCARC_MULTIGRID_ITERATIONS              !> Maximum number of iterations for multigrid method
PUBLIC SCARC_MULTIGRID_INTERPOL                !> Interpolation method for multigrid (AMG only)
PUBLIC SCARC_MULTIGRID_ACCURACY                !> Requested accuracy for multigrid method
PUBLIC SCARC_MULTIGRID_CYCLE                   !> Type of multigrid cycle
PUBLIC SCARC_MULTIGRID_COARSENING              !> Coarsening method for multigrid (AMG only)

PUBLIC SCARC_SMOOTH                            !> Smoother for multigrid method
PUBLIC SCARC_SMOOTH_SCOPE                      !> Scope of smoother (local/global)
PUBLIC SCARC_SMOOTH_ITERATIONS                 !> Maximum number of iterations for smoothing method
PUBLIC SCARC_SMOOTH_ACCURACY                   !> Requested accuracy for smoothing method
PUBLIC SCARC_SMOOTH_OMEGA                      !> Damping parameter for smoothing method

PUBLIC SCARC_COARSE                            !> Coarse grid solver for multigrid method
PUBLIC SCARC_COARSE_ITERATIONS                 !> Maximum number of iterations for coarse grid solver
PUBLIC SCARC_COARSE_ACCURACY                   !> Requested accuracy for coarse grid solver
PUBLIC SCARC_COARSE_OMEGA                      !> Relaxation parameter for coarse grid solver
PUBLIC SCARC_COARSE_LEVEL                      !> Coarse grid level

!> ------------------------------------------------------------------------------------------------
!> Miscellaneous declarations
!> ------------------------------------------------------------------------------------------------
!> General definitions
CHARACTER(40) :: SCARC_METHOD               = 'NONE'        !> Requested solver method (KRYLOV/MULTIGRID)
CHARACTER(40) :: SCARC_TWOLEVEL             = 'NONE'        !> Type of two-level method (NONE/ADDITIVE/MULTIPLICATIVE)
CHARACTER(40) :: SCARC_MATRIX_FORMAT        = 'COMPACT'     !> Type of matrix storage (COMPACT/BANDED)
CHARACTER(40) :: SCARC_DISCRETIZATION       = 'STRUCTURED'  !> Type of discretization (STRUCTURED/UNSTRUCTURED)

!> General iteration parameters
INTEGER       :: SCARC_ITERATIONS           =  0            !> Number of iterations of selected ScaRC solver
REAL (EB)     :: SCARC_RESIDUAL             =  0.0_EB       !> Residual of global selected solver
REAL (EB)     :: SCARC_CAPPA                =  0.0_EB       !> Convergence rate of selected ScarC solver
CHARACTER(40) :: SCARC_ACCURACY             = 'ABSOLUTE'    !> Accuracy type (ABSOLUTE/RELATIVE)
CHARACTER(6)  :: SCARC_PRECISION            = 'DOUBLE'      !> Single/double precision for preconditioning or LU-decomposition

!> Parameters for multigrid-type methods
CHARACTER(40) :: SCARC_MULTIGRID            = 'GEOMETRIC'   !> Type of MG-method (GEOMETRIC/ALGEBRAIC)
CHARACTER(40) :: SCARC_MULTIGRID_COARSENING = 'FALGOUT'     !> Coarsening strategy  (Falgout/RS3/A1/A2/...)
CHARACTER(1)  :: SCARC_MULTIGRID_CYCLE      = 'V'           !> Cycling type  (F/V/W)
INTEGER       :: SCARC_MULTIGRID_LEVEL      = -1            !> User defined number of MG-levels (optionally, otherwise maximum)
INTEGER       :: SCARC_MULTIGRID_ITERATIONS = 100           !> Max number of iterations
REAL (EB)     :: SCARC_MULTIGRID_ACCURACY   = 1.E-8_EB      !> Requested accuracy for convergence
CHARACTER(40) :: SCARC_MULTIGRID_INTERPOL   = 'CONSTANT'    !> Interpolation strategy (CONSTANT/BILINEAR)

!> Parameters for Krylov-type methods
CHARACTER(40) :: SCARC_KRYLOV            = 'CG'             !> Type of Krylov-method (CG/BICG)
INTEGER       :: SCARC_KRYLOV_ITERATIONS = 1000             !> Max number of iterations
REAL (EB)     :: SCARC_KRYLOV_ACCURACY   = 1.E-8_EB         !> Requested accuracy for convergence
CHARACTER(40) :: SCARC_KRYLOV_INTERPOL   = 'CONSTANT'       !> twolevel-interpolation (CONSTANT/BILINEAR)

!> Parameters for smoothing method (used in multigrids-methods)
CHARACTER(40) :: SCARC_SMOOTH            = 'SSOR'           !> Smoother for MG (JACOBI/SSOR)
CHARACTER(40) :: SCARC_SMOOTH_SCOPE      = 'GLOBAL'         !> Scope of action (LOCAL/GLOBAL)
INTEGER       :: SCARC_SMOOTH_ITERATIONS = 5                !> Max number of iterations
REAL (EB)     :: SCARC_SMOOTH_ACCURACY   = 1.E-8_EB         !> Requested accuracy for convergence
REAL (EB)     :: SCARC_SMOOTH_OMEGA      = 0.80E+0_EB       !> Relaxation parameter

!> Parameters for preconditioning method (used in Krylov-methods)
CHARACTER(40) :: SCARC_PRECON            = 'NONE'           !> Preconditioner for CG/BICG (JACOBI/SSOR/FFT/PARDISO/MG)
CHARACTER(40) :: SCARC_PRECON_SCOPE      = 'LOCAL'          !> Scope of action (LOCAL/GLOBAL)
INTEGER       :: SCARC_PRECON_ITERATIONS = 100              !> Max number of iterations
REAL (EB)     :: SCARC_PRECON_ACCURACY   = 1.E-10_EB        !> Requested accuracy for convergence
REAL (EB)     :: SCARC_PRECON_OMEGA      = 0.80E+0_EB       !> Relaxation parameter

!> Parameters for coarse grid method
CHARACTER(40) :: SCARC_COARSE            = 'ITERATIVE'      !> Type of coarse grid solver (ITERATIVE/DIRECT)
INTEGER       :: SCARC_COARSE_ITERATIONS = 100              !> Max number of iterations for iterative solver
REAL (EB)     :: SCARC_COARSE_ACCURACY   = 1.E-14_EB        !> Requested accuracy for iterative solver
REAL (EB)     :: SCARC_COARSE_OMEGA      = 0.80E+0_EB       !> Relaxation parameter
INTEGER       :: SCARC_COARSE_LEVEL      =  1               !> Coarse grid level for twolevel-Krylov method (default minimum level)

#ifdef WITH_MKL
!> Parameter for MKL solver
CHARACTER(40) :: SCARC_MKL       = 'GLOBAL'                 !> Type of MKL solver (LOCAL->Pardiso/GLOBAL->Cluster_Sparse_solver)
CHARACTER(40) :: SCARC_MKL_MTYPE = 'SYMMETRIC'              !> Type of MKL matrix (SYMMETRIC/UNSYMMETRIC)
#endif

!> Debugging parameters
CHARACTER(40) :: SCARC_ERROR_FILE = '.FALSE.'               !> Print ScaRC statistics into chid_scarc.csv (TRUE/FALSE)

INTEGER :: IERROR = 0

PRIVATE
!> ------------------------------------------------------------------------------------------------
!> Global constants
!> ------------------------------------------------------------------------------------------------
INTEGER, PARAMETER :: NSCARC_DISCRET_STRUCTURED      =  1, &    !> structured discretization 
                      NSCARC_DISCRET_UNSTRUCTURED    =  2       !> unstructured discretization

INTEGER, PARAMETER :: NSCARC_CELL_GASPHASE           =  1, &    !> gasphase cell
                      NSCARC_CELL_SOLID              =  2       !> solid cell

INTEGER, PARAMETER :: NSCARC_RHS_HOMOGENEOUS        =  1, &     !> homogeneous boundary conditions for global problem
                      NSCARC_RHS_INHOMOGENEOUS      =  2, &     !> inhomogeneous boundary conditions for global problem
                      NSCARC_RHS_DEFECT             =  3        !> RHS of precond. is set to defect of main iteration

INTEGER, PARAMETER :: NSCARC_STAGE_ONE               =  1, &    !> primary scope for solution vectors
                      NSCARC_STAGE_TWO               =  2       !> secondary scope for solution vectors

INTEGER, PARAMETER :: NSCARC_METHOD_KRYLOV           =  1, &    !> Krylov   -method as global solver
                      NSCARC_METHOD_MULTIGRID        =  2, &    !> multigrid-method as global solver
                      NSCARC_METHOD_LU               =  3, &    !> multigrid-method as global solver
                      NSCARC_METHOD_MGM              =  4       !> LU-decomposition based on MKL solver

INTEGER, PARAMETER :: NSCARC_KRYLOV_CG               =  1, &    !> CG   as Krylov solver
                      NSCARC_KRYLOV_CGBARO           =  2, &    !> CGBARO (Krylov with baroclinic effect), only testwise
                      NSCARC_KRYLOV_BICG             =  3, &    !> BICG as Krylov solver
                      NSCARC_KRYLOV_MAIN             =  1, &    !> Krylov solver as main solver
                      NSCARC_KRYLOV_COARSE           =  2       !> Krylov solver as coarse grid solver

INTEGER, PARAMETER :: NSCARC_MULTIGRID_GEOMETRIC     =  1, &    !> geometric multigrid
                      NSCARC_MULTIGRID_ALGEBRAIC     =  2, &    !> algebraic multigrid
                      NSCARC_MULTIGRID_MAIN          =  1, &    !> multigrid solver as main solver
                      NSCARC_MULTIGRID_PRECON        =  2       !> multigrid solver as preconditioner

INTEGER, PARAMETER :: NSCARC_MKL_LOCAL               =  1, &    !> local LU-decompositions (PARDISO)
                      NSCARC_MKL_GLOBAL              =  2, &    !> global LU-decomposition(CLUSTER_SPARSE_SOLVER) 
                      NSCARC_MKL_COARSE              =  3       !> LU-decomposition on coarse grid level

INTEGER, PARAMETER :: NSCARC_EXCHANGE_BASIC          =  1, &    !> initialize wall information
                      NSCARC_EXCHANGE_VECTOR         =  2, &    !> matrix-vector communication
                      NSCARC_EXCHANGE_PRESSURE       =  3, &    !> vector values along internal boundaries
                      NSCARC_EXCHANGE_DISCRET        =  4, &    !> exchange discretization information
                      NSCARC_EXCHANGE_WALL_INFO      =  5, &    !> initialize wall information
                      NSCARC_EXCHANGE_CELL_INDEX     =  6, &    !> internal transfer weights
                      NSCARC_EXCHANGE_CELL_WIDTH     =  7, &    !> neighboring grid resolution
                      NSCARC_EXCHANGE_CELL_NUMBER    =  8, &    !> neighboring mesh information
                      NSCARC_EXCHANGE_MATRIX_SIZE    =  9, &    !> neighboring matrix size
                      NSCARC_EXCHANGE_MATRIX_VALUE   = 10       !> neighboring matrix size

INTEGER, PARAMETER :: NSCARC_BROADCAST_SUM           =  1, &     !> broadcast local value and deliver sum of all
                      NSCARC_BROADCAST_PRODUCT       =  2, &     !> broadcast local value and deliver product of all
                      NSCARC_BROADCAST_MEAN          =  3, &     !> broadcast local value and deliver mean of all
                      NSCARC_BROADCAST_FIRST         =  4, &     !> broadcast local value and deliver first
                      NSCARC_BROADCAST_LAST          =  5, &     !> broadcast local value and deliver last
                      NSCARC_BROADCAST_MIN           =  6, &     !> broadcast local value and deliver minimum
                      NSCARC_BROADCAST_MAX           =  7, &     !> broadcast local value and deliver maximum
                      NSCARC_BROADCAST_MINLOC        =  8, &     !> broadcast local value and deliver minimum
                      NSCARC_BROADCAST_MAXLOC        =  9        !> broadcast local value and deliver maximum

INTEGER, PARAMETER :: NSCARC_RELAX_JACOBI            =  1, &     !> preconditioning by local JACOBI-methods
                      NSCARC_RELAX_SSOR              =  2, &     !> preconditioning by local SSOR-methods
                      NSCARC_RELAX_GSM               =  3, &     !> preconditioning by local GS-methods (matrix form)
                      NSCARC_RELAX_SGSM              =  4, &     !> preconditioning by local SGS-methods (matrix form)
                      NSCARC_RELAX_SORM              =  5, &     !> preconditioning by local SOR-methods (matrix form)
                      NSCARC_RELAX_SSORM             =  6, &     !> preconditioning by local SSOR-methods (matrix form)
                      NSCARC_RELAX_FFT               =  7, &     !> preconditioning by local FFT-methods
                      NSCARC_RELAX_GMG               =  8, &     !> preconditioning by local GMG-methods
                      NSCARC_RELAX_ILU               =  9, &     !> preconditioning by local ILU-decompositions (own)
                      NSCARC_RELAX_MKL               = 10        !> preconditioning by local LU-decompositions (MKL)
 
INTEGER, PARAMETER :: NSCARC_SCOPE_GLOBAL            =  1, &     !> scope of defect correction is global
                      NSCARC_SCOPE_LOCAL             =  2        !> scope of defect correction is local

INTEGER, PARAMETER :: NSCARC_TWOLEVEL_NONE           =  0, &     !> no two levels, only one level
                      NSCARC_TWOLEVEL_ADD            =  1, &     !> additive 2-level method
                      NSCARC_TWOLEVEL_MUL            =  2, &     !> multiplicative 2-level method
                      NSCARC_TWOLEVEL_MUL2           =  3, &     !> multiplicative 2-level method
                      NSCARC_TWOLEVEL_COARSE         =  4        !> only coarse grid

INTEGER, PARAMETER :: NSCARC_CYCLING_F               =  0, &     !> F-cycle for mg-method
                      NSCARC_CYCLING_V               =  1, &     !> V-cycle for mg-method
                      NSCARC_CYCLING_W               =  2, &     !> W-cycle for mg-method
                      NSCARC_CYCLING_SETUP           =  3, &     !> initialize cycle counts
                      NSCARC_CYCLING_RESET           =  4, &     !> reset cycle counts
                      NSCARC_CYCLING_PROCEED         =  5, &     !> proceed cycle counts
                      NSCARC_CYCLING_NEXT            =  6, &     !> perform next cycling loop
                      NSCARC_CYCLING_EXIT            =  7, &     !> exit cycling loop
                      NSCARC_CYCLING_PRESMOOTH       = -1, &     !> presmoothing cycle
                      NSCARC_CYCLING_POSTSMOOTH      =  1        !> postsmoothing cycle

INTEGER, PARAMETER :: NSCARC_STATE_PROCEED           =  0, &     !> proceed loop
                      NSCARC_STATE_CONV0             =  1, &     !> check convergence already for initial residual
                      NSCARC_STATE_CONV              =  2, &     !> check convergence for residual
                      NSCARC_STATE_DIVG              =  3        !> check divergence for residual

INTEGER, PARAMETER :: NSCARC_DEBUG_STACK             =  1, &     !> show matrix
                      NSCARC_DEBUG_DISCRET           =  2, &     !> show discretization information
                      NSCARC_DEBUG_MATRIX            =  3, &     !> show matrix
                      NSCARC_DEBUG_MATRIX_SYM        =  4, &     !> show symmetric matrix
                      NSCARC_DEBUG_WALL              =  5, &     !> show wall information
                      NSCARC_DEBUG_FACE              =  6, &     !> show face information
                      NSCARC_DEBUG_BDRY              =  7, &     !> show pressure_bc_index
                      NSCARC_DEBUG_SUBDIVISION       =  8        !> show subdivision

INTEGER, PARAMETER :: NSCARC_COARSENING_BASIC        =  1, &     !> basic coarsening
                      NSCARC_COARSENING_FALGOUT      =  2, &     !> parallel Falgout
                      NSCARC_COARSENING_RS3          =  3, &     !> parallel RS3
                      NSCARC_COARSENING_A1           =  4, &     !> aggressive 1 (path =1, length =2)
                      NSCARC_COARSENING_A2           =  5, &     !> aggressive 2 (path =2, length =2)
                      NSCARC_COARSENING_BDRY         =  6        !> FDSA2  : FDS variant similar to A2
 
INTEGER, PARAMETER :: NSCARC_COARSE_ITERATIVE        =  1, &     !> iterative solution of coarse grid problem
                      NSCARC_COARSE_DIRECT           =  2, &     !> direct solution of coarse grid problem
                      NSCARC_COARSE_KRYLOV           =  3, &     !> direct solution of coarse grid problem
                      NSCARC_COARSE_PARDISO          =  4, &     !> direct solution of coarse grid problem
                      NSCARC_COARSE_CLUSTER          =  5        !> direct solution of coarse grid problem

INTEGER, PARAMETER :: NSCARC_SIZE_MATRIX             =  1, &     !> size of system matrix for compact system
                      NSCARC_SIZE_TRANSFER           =  2        !> size of transfer matrices for compact system

INTEGER, PARAMETER :: NSCARC_SOLVER_MAIN             =  1, &     !> Krylov solver as main solver
                      NSCARC_SOLVER_PRECON           =  2, &     !> Multigrid solver as main solver
                      NSCARC_SOLVER_SMOOTH           =  3, &     !> Cluster sparse solver as main solver
                      NSCARC_SOLVER_COARSE           =  4        !> Cluster sparse solver as main solver

INTEGER, PARAMETER :: NSCARC_VECTOR_ONE_X            =  1, &     !> selection parameter for 1-stage vector X
                      NSCARC_VECTOR_ONE_B            =  2, &     !> selection parameter for 1-stage vector F
                      NSCARC_VECTOR_ONE_D            =  3, &     !> selection parameter for 1-stage vector G
                      NSCARC_VECTOR_ONE_R            =  4, &     !> selection parameter for 1-stage vector D
                      NSCARC_VECTOR_ONE_V            =  5, &     !> selection parameter for 1-stage vector R
                      NSCARC_VECTOR_ONE_Y            =  6, &     !> selection parameter for 1-stage vector Y
                      NSCARC_VECTOR_ONE_Z            =  7, &     !> selection parameter for 1-stage vector Z
                      NSCARC_VECTOR_ONE_E            =  8, &     !> selection parameter for 1-stage vector Z
                      NSCARC_VECTOR_TWO_X            =  9, &     !> selection parameter for 2-stage vector X
                      NSCARC_VECTOR_TWO_B            = 10, &     !> selection parameter for 2-stage vector F
                      NSCARC_VECTOR_TWO_D            = 11, &     !> selection parameter for 2-stage vector G
                      NSCARC_VECTOR_TWO_R            = 12, &     !> selection parameter for 2-stage vector D
                      NSCARC_VECTOR_TWO_V            = 13, &     !> selection parameter for 2-stage vector R
                      NSCARC_VECTOR_TWO_Y            = 14, &     !> selection parameter for 2-stage vector Y
                      NSCARC_VECTOR_TWO_Z            = 15, &     !> selection parameter for 2-stage vector Z
                      NSCARC_VECTOR_TWO_E            = 16, &     !> selection parameter for 2-stage vector Z
                      NSCARC_VECTOR_H                = 17, &     !> selection parameter for vector H
                      NSCARC_VECTOR_HS               = 18        !> selection parameter for vector HS

INTEGER, PARAMETER :: NSCARC_MATRIX_COMPACT          =  1, &     !> matrix in BANDED storage format
                      NSCARC_MATRIX_BANDED           =  2, &     !> matrix in diagonal storage format
                      NSCARC_MATRIX_FREE             =  3, &     !> matrix in free storage format
                      NSCARC_MATRIX_CONDENSED        =  4, &     !> matrix condensing for Neumann problems
                      NSCARC_MATRIX_GENERAL          =  5, &     !> general matrix (not necessarily symmetric)
                      NSCARC_MATRIX_SYMMETRIC        =  6        !> symmetric matrix

INTEGER, PARAMETER :: NSCARC_ACCURACY_ABSOLUTE       =  1, &     !> absolute accuracy must be reached
                      NSCARC_ACCURACY_RELATIVE       =  2        !> relative accuracy must be reached

INTEGER, PARAMETER :: NSCARC_ERROR_PARSE_INPUT       =  1, &     !> wrong input parameter
                      NSCARC_ERROR_MKL_INTERNAL      =  2, &     !> pardiso solver not available
                      NSCARC_ERROR_MKL_PARDISO       =  3, &     !> pardiso solver not available
                      NSCARC_ERROR_MKL_CLUSTER       =  4, &     !> cluster_sparse_solver not available
                      NSCARC_ERROR_MKL_STORAGE       =  5, &     !> wrong storage scheme for MKL-solvers
                      NSCARC_ERROR_BOUNDARY_SUM      =  6, &     !> wrong sum of elements along boundary
                      NSCARC_ERROR_BOUNDARY_TYPE     =  7, &     !> wrong boundary type
                      NSCARC_ERROR_GRID_NUMBER       =  8, &     !> cell number not divisable by two
                      NSCARC_ERROR_GRID_NUMBERX      =  9, &     !> cell number not divisable by two
                      NSCARC_ERROR_GRID_NUMBERY      = 10, &     !> cell number not divisable by two
                      NSCARC_ERROR_GRID_NUMBERZ      = 11, &     !> cell number not divisable by two
                      NSCARC_ERROR_GRID_RESOLUTION   = 12, &     !> error with grid resolution
                      NSCARC_ERROR_GRID_INDEX        = 13, &     !> error with grid index
                      NSCARC_ERROR_NEIGHBOR_TYPE     = 14, &     !> wrong neighbor type
                      NSCARC_ERROR_NEIGHBOR_NUMBER   = 15, &     !> wrong neighbor number
                      NSCARC_ERROR_MATRIX_SETUP      = 16, &     !> error in matrix setup
                      NSCARC_ERROR_MATRIX_SIZE       = 17, &     !> error in matrix size
                      NSCARC_ERROR_MATRIX_SYMMETRY   = 18, &     !> matrix not symmetric
                      NSCARC_ERROR_MATRIX_SUBDIAG    = 19, &     !> subdiagonal missing
                      NSCARC_ERROR_MATRIX_COPY       = 20, &     !> subdiagonal missing
                      NSCARC_ERROR_STACK_SOLVER      = 21, &     !> error in solver stack
                      NSCARC_ERROR_STACK_MESSAGE     = 22, &     !> error with stack message
                      NSCARC_ERROR_FFT_DISCRET       = 23, &     !> wrong unstructured discretization for FFT
                      NSCARC_ERROR_VECTOR_LENGTH     = 24, &     !> error in vector length
                      NSCARC_ERROR_MULTIGRID_LEVEL   = 25, &     !> wrong multigrid level
                      NSCARC_ERROR_AMG_MISSING       = 26, &     !> AMG-method is currently missing
                      NSCARC_ERROR_RHS_SETUP         = 27, &     !> error with rhs setup
                      NSCARC_ERROR_POC_STOP          = 28, &     !> temporary stop, only proof of concept
                      NSCARC_ERROR_BICG_DISABLED     = 29        !> BICG temporarily disabled
  
INTEGER, PARAMETER :: NSCARC_PRECISION_FB            =  1, &     !> single precision for preconditioning or LU-decomposition
                      NSCARC_PRECISION_EB            =  2        !> double precision for preconditioning or LU-decomposition

INTEGER, PARAMETER :: NSCARC_INTERPOL_STANDARD       =  1, &     !> standard interpolation
                      NSCARC_INTERPOL_CONSTANT       =  2, &     !> standard interpolation
                      NSCARC_INTERPOL_BILINEAR       =  3, &     !> standard interpolation
                      NSCARC_INTERPOL_CLASSICAL      =  4, &     !> classical interpolation
                      NSCARC_INTERPOL_CLASSICAL2     =  5, &     !> classical interpolation
                      NSCARC_INTERPOL_DIRECT         =  6, &     !> direct interpolation
                      NSCARC_INTERPOL_DIRECT_BDRY    =  7        !> direct interpolation with special boundary

INTEGER, PARAMETER :: NSCARC_LEVEL_MIN               =  0, &     !> minimum multigrid level
                      NSCARC_LEVEL_MAX               = 15, &     !> maximum multigrid level
                      NSCARC_LEVEL_SINGLE            =  1, &     !> only one grid level needed 
                      NSCARC_LEVEL_MULTI             =  2, &     !> multiple grid levels needed for GMG meethod
                      NSCARC_LEVEL_AMG               =  3        !> number of grid levels specified by AMG method

INTEGER, PARAMETER :: NSCARC_DUMP_RHS                =  1, &     !> dump rhs
                      NSCARC_DUMP_PRES               =  2        !> dump pressure

INTEGER, PARAMETER :: NSCARC_DIAG_MAIN               =  1, &     !> standard 5- or 7-point stencil
                      NSCARC_DIAG_LOWER              =  2, &     !> standard 5- or 7-point stencil
                      NSCARC_DIAG_UPPER              =  3        !> standard 5- or 7-point stencil
  
INTEGER, PARAMETER :: NSCARC_COUPLING_MAX            = 10        !> maximum of possible couplings in stencil

INTEGER, PARAMETER :: NSCARC_MAX_FACE_NEIGHBORS      = 10        !> max number neighbors per mesh face
INTEGER, PARAMETER :: NSCARC_MAX_STENCIL             =  7

INTEGER, PARAMETER :: NSCARC_UNDEFINED_INT           = -1, &     !> undefined integer value
                      NSCARC_ZERO_INT                =  0, &     !> zero integer value
                      NSCARC_ONE_INT                 =  1        !> one integer value

REAL(EB), PARAMETER:: NSCARC_UNDEFINED_REAL_EB       = -1.0_EB, &    !> undefined real value
                      NSCARC_ZERO_REAL_EB            =  0.0_EB, &    !> zero real value
                      NSCARC_ONE_REAL_EB             =  1.0_EB       !> one real value

REAL(EB), PARAMETER:: NSCARC_UNDEFINED_REAL_FB       = -1.0_FB, &    !> undefined real value
                      NSCARC_ZERO_REAL_FB            =  0.0_FB, &    !> zero real value
                      NSCARC_ONE_REAL_FB             =  1.0_FB       !> one real value

INTEGER, PARAMETER :: NSCARC_INIT_UNDEF          = -99, &        !> initialize allocated array as undefined
                      NSCARC_INIT_NONE               =  -1, &        !> do not initialize allocated arrays
                      NSCARC_INIT_ZERO               =   0, &        !> initialize allocated array with zero
                      NSCARC_INIT_ONE                =   1           !> initialize allocated array with one

INTEGER, PARAMETER :: NSCARC_HUGE_INT                = -999999999    !> undefined integer value
REAL(EB), PARAMETER:: NSCARC_HUGE_REAL               = -999999999_EB !> undefined integer value

INTEGER, PARAMETER :: NSCARC_STACK_ZERO              =   0, &        !> zero stack 
                      NSCARC_STACK_ROOT              =   1, &        !> root stack
                      NSCARC_STACK_MAX               =  10, &        !> maximum stack
                      NSCARC_STACK_NOPARENT          = -99           !> no stack information

REAL(EB), PARAMETER:: NSCARC_THRESHOLD_CONVERGENCE   = 1.0E-15_EB, & !> threshold for convergence
                      NSCARC_THRESHOLD_DIVGERGENCE   = 1.0E+15_EB    !> threshold for divergence

INTEGER, PARAMETER :: NSCARC_NONE                    = -123456789    !> dummy integer value for requests
CHARACTER(40), PARAMETER :: SCARC_NONE = 'NONE'                      !> dummy character value for requests


!> --------------------------------------------------------------------------------------------
!> Logical indicators for different methods and mechanisms
!> --------------------------------------------------------------------------------------------
LOGICAL :: IS_STRUCTURED   = .FALSE.                      !> is structured discretization?
LOGICAL :: IS_UNSTRUCTURED = .FALSE.                      !> is unstructured discretization?
LOGICAL :: IS_CG           = .FALSE.                      !> is Krylov-method?
LOGICAL :: IS_CG_GMG       = .FALSE.                      !> is Krylov-method with GMG-preconditioning?
LOGICAL :: IS_CG_ADD       = .FALSE.                      !> is additive twolevel-Krylov-method?
LOGICAL :: IS_CG_MUL       = .FALSE.                      !> is multiplicative Twolevel-Krylov-method ?
LOGICAL :: IS_CG_COARSE    = .FALSE.                      !> is only coarse grid preconditiner?
LOGICAL :: IS_MG           = .FALSE.                      !> is Multigrid-method?
LOGICAL :: IS_GMG          = .FALSE.                      !> is Geometric Multigrid-method?
LOGICAL :: IS_FFT          = .FALSE.                      !> is FFT-method?
LOGICAL :: IS_MKL          = .FALSE.                      !> is MKL-method?
LOGICAL :: IS_MKL_LEVEL(NSCARC_LEVEL_MAX)  = .FALSE.      !> is level-dependent MKL method?

LOGICAL :: HAS_CSV           = .FALSE.                    !> has CSV-file to be dumped out
LOGICAL :: HAS_TWO_LEVELS     = .FALSE.                    !> has two grid levels?
LOGICAL :: HAS_MULTIPLE_LEVELS   = .FALSE.                    !> has multiple grid levels?
LOGICAL :: HAS_MULTIPLE_DISCRETIZATIONS = .FALSE.                    !> has multiple discretization types?


!> --------------------------------------------------------------------------------------------
!> Globally used types for description of different solvers
!> --------------------------------------------------------------------------------------------
INTEGER :: TYPE_DISCRET      = NSCARC_DISCRET_STRUCTURED    !> default type of discretization (structured/unstructured)
INTEGER :: TYPE_METHOD       = NSCARC_METHOD_KRYLOV         !> default type of ScaRC method
INTEGER :: TYPE_SOLVER       = NSCARC_SOLVER_MAIN           !> default type of surrounding solver stage
INTEGER :: TYPE_STAGE        = NSCARC_STAGE_ONE             !> default type of surrounding solver stage
INTEGER :: TYPE_TWOLEVEL     = NSCARC_TWOLEVEL_NONE         !> default type of two-level method
INTEGER :: TYPE_INTERPOL     = NSCARC_INTERPOL_CONSTANT     !> default type of interpolation method
INTEGER :: TYPE_KRYLOV       = NSCARC_KRYLOV_CG             !> default type of Krylov method (CG/BICG)
INTEGER :: TYPE_MULTIGRID    = NSCARC_MULTIGRID_GEOMETRIC   !> default type of multigrid method (GMG/AMG)
INTEGER :: TYPE_ACCURACY     = NSCARC_ACCURACY_ABSOLUTE     !> default type of requested accuracy
INTEGER :: TYPE_CYCLING      = NSCARC_CYCLING_V             !> default type of cycling for multigrid method
INTEGER :: TYPE_COARSE       = NSCARC_COARSE_ITERATIVE      !> default type of coarse grid solver for multigrid method
INTEGER :: TYPE_MATRIX       = NSCARC_MATRIX_COMPACT        !> default type of storage for matrix
INTEGER :: TYPE_COARSENING   = NSCARC_UNDEFINED_INT         !> no default type of coarsening algorithm for AMG
INTEGER :: TYPE_RELAX        = NSCARC_UNDEFINED_INT         !> no default type of preconditioner for iterative solver
INTEGER :: TYPE_PRECON       = NSCARC_UNDEFINED_INT         !> no default type of preconditioner for iterative solver
INTEGER :: TYPE_PRECON_SCOPE = NSCARC_SCOPE_LOCAL           !> default type of preconditioning is local
INTEGER :: TYPE_PRECISION    = NSCARC_UNDEFINED_INT         !> no default type of preconditioner for iterative solver
INTEGER :: TYPE_SMOOTH       = NSCARC_UNDEFINED_INT         !> no default type of smoother for multigrid method
INTEGER :: TYPE_SMOOTH_SCOPE = NSCARC_SCOPE_LOCAL           !> default type of smoothing is local
INTEGER :: TYPE_EXCHANGE     = NSCARC_UNDEFINED_INT         !> no default type of data exchange
INTEGER :: TYPE_VECTOR       = NSCARC_UNDEFINED_INT         !> no default type of vector to point to
INTEGER :: TYPE_PARENT       = NSCARC_UNDEFINED_INT         !> no default type of parent (calling) solver
INTEGER :: TYPE_MKL          = NSCARC_UNDEFINED_INT         !> no default type of MKL method (PARDISO/CLUSTER_SPARSE_SOLVER)
INTEGER :: TYPE_MKL_LEVEL(NSCARC_LEVEL_MAX) = NSCARC_UNDEFINED_INT  !> no default type of MKL for single levels

!> ------------------------------------------------------------------------------------------------
!> Globally used parameters
!> ------------------------------------------------------------------------------------------------
!> discretization information
INTEGER :: NLEVEL_MIN, NLEVEL_MAX                         !> minimum and maximum number of multigrid levels
INTEGER :: N_CELLS_GLOBAL(NSCARC_LEVEL_MAX) = 0           !> number of global cells
INTEGER :: N_DIRIC_GLOBAL(NSCARC_LEVEL_MAX) = 0           !> global number of Dirichlet BCs

!> stack information
INTEGER :: N_STACK_TOTAL                                  !> maximum number of used solvers in stack

!> communication parameters
INTEGER :: N_REQ, N_EXCHANGE, TAG                         !> Variables for data exchange
INTEGER :: SNODE, RNODE                                   !> Process identifier for data exchange

INTEGER,  ALLOCATABLE, DIMENSION (:)  :: REQ              !> Request array for data exchange
INTEGER,  ALLOCATABLE, DIMENSION (:)  :: COUNTS           !> Counter array for data exchange
INTEGER,  ALLOCATABLE, DIMENSION (:)  :: DISPLS           !> Displacement array for data exchange
INTEGER,  ALLOCATABLE, DIMENSION (:)  :: LOCAL_INT        !> Local integer data array for data exchange
REAL(EB), ALLOCATABLE, DIMENSION (:)  :: LOCAL_REAL       !> Local real data array for data exchange

INTEGER :: GLOBAL_INT
REAL(EB) :: GLOBAL_REAL

!> Face information
INTEGER :: FACE_ORIENTATION(6) = (/1,-1,2,-2,3,-3/)
INTEGER :: FACE_ORDER_XYZ(6) = (/1,-1,2,-2,3,-3/)           !> Coordinate direction related order of mesh faces

INTEGER :: FACE_NEIGHBORS(-3:3, NSCARC_MAX_FACE_NEIGHBORS)
INTEGER :: NUM_FACE_NEIGHBORS(-3:3)

!> Neighbor information
INTEGER :: MESH_NEIGHBORS(6*NSCARC_MAX_FACE_NEIGHBORS)
INTEGER :: NUM_MESH_NEIGHBORS

!> --------------------------------------------------------------------------------------------
!> Some auxiliary parameters used for transfer operators
!> --------------------------------------------------------------------------------------------
REAL(EB), PARAMETER :: SCALR  = 0.015625_EB
REAL(EB), PARAMETER :: SCALP  = 0.0625_EB
REAL(EB), PARAMETER :: W1     =  1.0_EB
REAL(EB), PARAMETER :: W3     =  3.0_EB
REAL(EB), PARAMETER :: W4     =  4.0_EB
REAL(EB), PARAMETER :: W9     =  9.0_EB
REAL(EB), PARAMETER :: W12    = 12.0_EB
REAL(EB), PARAMETER :: W16    = 16.0_EB

!> -----------------------------------------------------------------------------------------------
!> Time measurement
!> ------------------------------------------------------------------------------------------------
TYPE SCARC_TIME_TYPE
REAL(EB) :: OVERALL   = 0.0_EB               !> complete time
REAL(EB) :: SOLVER    = 0.0_EB               !> time for solver (general version)
REAL(EB) :: CLUSTER   = 0.0_EB               !> time for cluster solver
REAL(EB) :: PARDISO   = 0.0_EB               !> time for pardiso solver
REAL(EB) :: KRYLOV    = 0.0_EB               !> time for krylov solver
REAL(EB) :: MULTIGRID = 0.0_EB               !> time for multigrid solver
REAL(EB) :: MATVEC    = 0.0_EB               !> time for matrix vector multiplication
REAL(EB) :: SCALPROD  = 0.0_EB               !> time for scalar product
REAL(EB) :: L2NORM    = 0.0_EB               !> time for l2-norm
REAL(EB) :: PRECON    = 0.0_EB               !> time for preconditioner
REAL(EB) :: SMOOTH    = 0.0_EB               !> time for smoother
REAL(EB) :: COARSE    = 0.0_EB               !> time for coarse grid solver
REAL(EB) :: EXCHANGE  = 0.0_EB               !> time for data exchange
END TYPE SCARC_TIME_TYPE

!> -----------------------------------------------------------------------------------------------
!> Messaging and debugging mechanisms
!> ------------------------------------------------------------------------------------------------
TYPE SCARC_MESSAGE_TYPE
CHARACTER(100) :: TEXT
CHARACTER(60) :: FILE_DEBUG, FILE_STAT, FILE_DUMP, FILE_VERBOSE
INTEGER :: LU_DEBUG = 0, LU_STAT = 0, LU_DUMP = 0, LU_VERBOSE = 0
END TYPE SCARC_MESSAGE_TYPE

!> --------------------------------------------------------------------------------------------
!> Face information related to wall cells and neighbors
!> --------------------------------------------------------------------------------------------
TYPE SCARC_FACE_TYPE
REAL(EB), POINTER, DIMENSION(:) :: DH              !> adjacent grid sizes
INTEGER, ALLOCATABLE, DIMENSION(:) :: NEIGHBORS    !> adjacent neighbors
INTEGER :: N_NEIGHBORS = 0                         !> number of adjacent neighbors
INTEGER :: NFC, NFW                                !> number of cells and wall cells along face
INTEGER :: NFX, NFY, NFZ                           !> local face dimensions
INTEGER :: NCPL                                    !> number of adjacent couplings
INTEGER :: IWG_PTR                                 !> first (global) IW number to that face
INTEGER :: IOFFSET_WALL   = 0                      !> counter for wall cells over all faces
END TYPE SCARC_FACE_TYPE

!> --------------------------------------------------------------------------------------------
!> Wall information related to neighbors and BC's
!> --------------------------------------------------------------------------------------------
TYPE SCARC_WALL_TYPE

!> state/boundary/orientation/neighboring information of wall cells
INTEGER :: STATE                                  !> Degree of freedom and state of related cell (gasphase/solid)
INTEGER :: BTYPE                                  !> type of wall cell (Dirichlet/Neumann/Internal)
INTEGER :: BOUNDARY_TYPE = 0                      !> state of wall cell (Solid/Interpolated/Open))
INTEGER :: IOR = 0                                !> orientation of wall cell
INTEGER :: NCPL = 1                               !> number of couplings at wall cell (depending on resolution of neighbor)
INTEGER :: NOM = 0                                !> adjacent neighbor at wall cell

!> counters for different elements
INTEGER :: ICW = NSCARC_UNDEFINED_INT             !> internal wall cell for IW
INTEGER :: ICO = NSCARC_UNDEFINED_INT             !> overlapping cell for IW
INTEGER :: IWL = NSCARC_UNDEFINED_INT             !> corresponding local wall cell number for neighbor NOM

!> coordinates in x-,y- and z-direction and extended/ghost cell information
INTEGER :: IXG, IYG, IZG                          !> x-, y- and z-indices of ghost cells
INTEGER :: IXW, IYW, IZW                          !> x-, y- and z-indices of (internal) wall cells
INTEGER :: IXN(2), IYN(2), IZN(2)                 !> x-, y- and z-indices of neighboring cells
INTEGER, ALLOCATABLE, DIMENSION(:) :: ICE         !> extended cell for IW
INTEGER, ALLOCATABLE, DIMENSION(:) :: ICG         !> ghost cell for IW
END TYPE SCARC_WALL_TYPE

!> --------------------------------------------------------------------------------------------
!> Mappings between different discretization description arrays
!> --------------------------------------------------------------------------------------------
TYPE SCARC_POINTER_TYPE
INTEGER :: ICG = 0                                         !> ghost cell pointer
INTEGER :: ICO = 0                                         !> overlapping cell pointer
INTEGER :: ICE = 0                                         !> extended cell pointer
INTEGER :: IWL = 0                                         !> local wall cell pointer
INTEGER, ALLOCATABLE, DIMENSION (:)  :: ICE_TO_IWG        !> mapping from ICE to IWG
INTEGER, ALLOCATABLE, DIMENSION (:)  :: ICE_TO_IWL        !> mapping from ICE to IWL
INTEGER, ALLOCATABLE, DIMENSION (:)  :: ICE_TO_ICG        !> mapping from ICE to ICG
INTEGER, ALLOCATABLE, DIMENSION (:)  :: ICE_TO_ICN        !> mapping from ICE to ICN
INTEGER, ALLOCATABLE, DIMENSION (:)  :: ICG_TO_IWG        !> mapping from ICG to IWG
INTEGER, ALLOCATABLE, DIMENSION (:)  :: ICG_TO_ICE        !> mapping from ICG to ICE
INTEGER, ALLOCATABLE, DIMENSION (:)  :: ICG_TO_ICO        !> mapping from ICG to ICE
INTEGER, ALLOCATABLE, DIMENSION (:)  :: IWL_TO_IWG        !> mapping from IWL to IWG
INTEGER, ALLOCATABLE, DIMENSION (:)  :: IWL_TO_ICW        !> mapping from IWL to ICW
INTEGER, ALLOCATABLE, DIMENSION (:)  :: IWL_TO_ICO        !> mapping from IWL to ICO
INTEGER, ALLOCATABLE, DIMENSION (:,:) :: IWL_TO_ICG        !> mapping from IWL to ICG
REAL(EB), ALLOCATABLE, DIMENSION (:)  :: ICE_TO_VAL        !> mapping from ICE to VAL
!INTEGER, ALLOCATABLE, DIMENSION (:)  :: ICN_TO_ICE        !> mapping from ICN to ICE    ! CAUTION: HUGE
END TYPE SCARC_POINTER_TYPE

!> --------------------------------------------------------------------------------------------
!> Obstruction information
!> --------------------------------------------------------------------------------------------
TYPE SCARC_OBST_TYPE
INTEGER :: I1, I2, J1, J2, K1, K2                          !> cell indices of obstructions
END TYPE SCARC_OBST_TYPE

!> --------------------------------------------------------------------------------------------
!> Matrix entries which will be stored and exchanged during generation of condensed system
!> --------------------------------------------------------------------------------------------
TYPE SCARC_CONDENSED_TYPE
INTEGER :: NROW, NCOL
#ifdef WITH_MKL_FB
REAL(FB) :: VAL_ORIG(NSCARC_MAX_STENCIL) = 0.0_EB          !> original values (double precision)
REAL(FB) :: VAL_COND(NSCARC_MAX_STENCIL) = 0.0_EB          !> condensed values (double precision)
#else
REAL(EB) :: VAL_ORIG(NSCARC_MAX_STENCIL) = 0.0_EB          !> original values (single precision)
REAL(EB) :: VAL_COND(NSCARC_MAX_STENCIL) = 0.0_EB          !> condensed values (single precision)
#endif
INTEGER :: COL(NSCARC_MAX_STENCIL) = 0                    !> column pointers
INTEGER :: PTR(NSCARC_MAX_STENCIL) = 0                    !> storage pointer
END TYPE SCARC_CONDENSED_TYPE

!> --------------------------------------------------------------------------------------------
!> Compact sparse row (COMPACT) storage type for matrices
!> Is based on three arrays:
!>    - non-zero matrix values
!>    - corresponding columns pointers
!>    - row pointers
!> --------------------------------------------------------------------------------------------
TYPE SCARC_MATRIX_COMPACT_TYPE
INTEGER :: POS(-3:3)                                      !> Position of IOR's in STENCIL
INTEGER :: NSTENCIL                                       !> number of points in matrix stencil
INTEGER :: NA                                             !> number of matrix values 
INTEGER :: NR                                             !> number of matrix rows
INTEGER :: NSTORE = 0
#ifdef WITH_MKL_FB
REAL(FB), ALLOCATABLE, DIMENSION (:) :: VAL_FB             !> values of matrix (single precision)
REAL(FB), DIMENSION (-3:3)           :: STENCIL_FB         !> store basic stencil information in single precision
#else
REAL(EB), ALLOCATABLE, DIMENSION (:) :: VAL                !> values of matrix (real precision)
REAL(EB), ALLOCATABLE, DIMENSION (:) :: ILU                !> ILU-decomposition
REAL(EB), ALLOCATABLE, DIMENSION (:) :: GSM                !> GS-decomposition in matrix form
REAL(EB), ALLOCATABLE, DIMENSION (:) :: SGSM               !> symmetric GS-decomposition in matrix form
REAL(EB), ALLOCATABLE, DIMENSION (:) :: SORM               !> SOR-decomposition in matrix form
REAL(EB), ALLOCATABLE, DIMENSION (:) :: SSORM              !> SSOR-decomposition in matrix form
REAL(EB), DIMENSION (-3:3)           :: STENCIL            !> store basic stencil information in single precision
#endif
INTEGER,  ALLOCATABLE, DIMENSION (:) :: ROW                !> row pointer
INTEGER,  ALLOCATABLE, DIMENSION (:) :: COL                !> column pointers
INTEGER,  ALLOCATABLE, DIMENSION (:) :: COL_GLOBAL         !> column pointer for global numbering
TYPE (SCARC_CONDENSED_TYPE) :: STORE(NSCARC_MAX_STENCIL)
END TYPE SCARC_MATRIX_COMPACT_TYPE

!> --------------------------------------------------------------------------------------------
!> Banded storage type for matrices
!> The entries are stored one diagonal after the other
!> Missing entries of subdiagonals are filled with zero
!> Is based on two arrays:
!>      - non-zero matrix entries diagonal-wise
!>      - the offsets from the main diagonal
!> --------------------------------------------------------------------------------------------
TYPE SCARC_MATRIX_BANDED_TYPE
INTEGER :: NSTENCIL                                    !> number of points in matrix stencil
INTEGER :: POS(-3:3)                                   !> position of IOR's in STENCIL and in matrix storage array
INTEGER :: NA, NAS                                     !> number of matrix values in general and symmetric cass
INTEGER :: NDIAG, NLEN                                 !> length of main diagonal
INTEGER :: NLO, NUP                                    !> number of lower and upper diagonals
#ifdef WITH_MKL_FB
REAL(FB), ALLOCATABLE, DIMENSION (:,:) :: VAL_FB        !> values of matrix (double precision)
REAL(FB), DIMENSION (-3:3)             :: STENCIL_FB    !> store basic stencil information in single precision
#else
REAL(EB), ALLOCATABLE, DIMENSION (:,:) :: VAL           !> values of matrix (single precision)
REAL(EB), DIMENSION (-3:3)             :: STENCIL       !> store basic stencil information in double precision
#endif
INTEGER,  DIMENSION (-3:3)             :: OFFSET        !> offset pointers
END TYPE SCARC_MATRIX_BANDED_TYPE

!> --------------------------------------------------------------------------------------------
!> Workspace for FFT preconditioners
!> --------------------------------------------------------------------------------------------
TYPE SCARC_FFT_TYPE
INTEGER :: LSAVE, LWORK
INTEGER :: LBC, MBC, NBC
INTEGER :: ITRN, JTRN, KTRN
INTEGER :: IBAR0, JBAR0, KBAR0
INTEGER :: ITRN0, JTRN0, KTRN0
INTEGER :: IS0=0, IF0=0, JS0=0, JF0=0, KS0=0, KF0=0
REAL(EB) :: XS0, XF0, YS0, YF0, ZS0, ZF0
REAL(EB) :: POIS_PTB = 0.0_EB, XLM = 0.0_EB
REAL(EB), ALLOCATABLE, DIMENSION (:)       :: SAVE1, WORK, HX
REAL(EB), ALLOCATABLE, DIMENSION (:, :)    :: BXS, BXF, BYS, BYF, BZS, BZF
REAL(EB), ALLOCATABLE, DIMENSION (:, :, :) :: PRHS
END TYPE SCARC_FFT_TYPE

#ifdef WITH_MKL
!> --------------------------------------------------------------------------------------------
!> MKL information
!> --------------------------------------------------------------------------------------------
TYPE SCARC_MKL_TYPE
INTEGER, ALLOCATABLE :: IPARM(:)
INTEGER :: MAXFCT, MNUM, MTYPE, PHASE, NRHS, ERROR, MSGLVL
INTEGER :: PERM(1)
TYPE (MKL_PARDISO_HANDLE),               ALLOCATABLE :: PT_H(:), PT(:)
TYPE (MKL_CLUSTER_SPARSE_SOLVER_HANDLE), ALLOCATABLE :: CT_H(:), CT(:)
END TYPE SCARC_MKL_TYPE
#endif

!> --------------------------------------------------------------------------------------------
!> Different scopes for solution, rhs and auxiliary vectors of different solvers
!> --------------------------------------------------------------------------------------------
TYPE SCARC_STAGE_TYPE
REAL (EB), ALLOCATABLE, DIMENSION (:) :: X          !> solution vector double precision
REAL (EB), ALLOCATABLE, DIMENSION (:) :: B          !> right hand side vector double precision
REAL (EB), ALLOCATABLE, DIMENSION (:) :: R          !> residual vector double precision
REAL (EB), ALLOCATABLE, DIMENSION (:) :: D          !> auxiliary vector double precision
REAL (EB), ALLOCATABLE, DIMENSION (:) :: V          !> auxiliary vector double precision
REAL (EB), ALLOCATABLE, DIMENSION (:) :: Y          !> auxiliary vector double precision
REAL (EB), ALLOCATABLE, DIMENSION (:) :: Z          !> auxiliary vector double precision
#ifdef WITH_SCARC_DEBUG
REAL (EB), ALLOCATABLE, DIMENSION (:) :: E          !> error vector double precision
#endif
#ifdef WITH_MKL_FB
REAL (FB), ALLOCATABLE, DIMENSION (:) :: X_FB       !> solution vector vector single precision
REAL (FB), ALLOCATABLE, DIMENSION (:) :: B_FB       !> right hand side vector single precision
REAL (FB), ALLOCATABLE, DIMENSION (:) :: R_FB       !> residual vector single precision
REAL (FB), ALLOCATABLE, DIMENSION (:) :: D_FB       !> auxiliary vector single precision
REAL (FB), ALLOCATABLE, DIMENSION (:) :: V_FB       !> auxiliary vector single precision
REAL (FB), ALLOCATABLE, DIMENSION (:) :: Y_FB       !> auxiliary vector single precision
REAL (FB), ALLOCATABLE, DIMENSION (:) :: Z_FB       !> auxiliary vector single precision
#endif
END TYPE SCARC_STAGE_TYPE

!> --------------------------------------------------------------------------------------------
!> Information related to discretization type
!> --------------------------------------------------------------------------------------------
TYPE SCARC_DISCRET_TYPE

INTEGER :: NC                                               !> number of cells

!> Some basic dimensions
INTEGER :: N_CELLS = 0                                      !> number of cells needed for matrix 
INTEGER :: N_DIRIC = 0                                      !> number of Dirichlet BCs
INTEGER :: N_NEUMANN = 0                                    !> number of Neumann BCs

!> Global and local numbers of cells and offsets between meshes
INTEGER :: NC_GLOBAL = 0                                    !> global number of cells in all meshes
INTEGER, ALLOCATABLE, DIMENSION (:) :: NC_LOCAL             !> number of cells in local meshes
INTEGER, ALLOCATABLE, DIMENSION (:) :: NC_OFFSET            !> offset in cell numbering between meshes

!> number and state of single cells
INTEGER, ALLOCATABLE, DIMENSION (:,:,:) :: CELL_NUMBER      !> numbering of single cells

!> Matrices in different storage techniques
TYPE (SCARC_MATRIX_BANDED_TYPE)  :: AB                      !> Poisson matrix in banded storage technique
TYPE (SCARC_MATRIX_COMPACT_TYPE) :: AC                      !> Poisson matrix in compact storage technique 
#ifdef WITH_MKL
TYPE (SCARC_MATRIX_COMPACT_TYPE) :: AC_SYM                  !> symmetric part of Poisson matrix (only for MKL)
#endif

TYPE (SCARC_WALL_TYPE), ALLOCATABLE, DIMENSION(:) :: WALL   !> wall information

!> Different mapping arrays for communication of overlapping and wall cells
TYPE (SCARC_POINTER_TYPE) :: PTR                            !> pointers between different meshes

END TYPE SCARC_DISCRET_TYPE

!> --------------------------------------------------------------------------------------------
!> Multigrid type - to be extended for algebraic multigrid
!> --------------------------------------------------------------------------------------------
TYPE SCARC_MULTIGRID_TYPE
INTEGER :: CYCLING(2) = 0                                  !> Counter for multigrid cycling
END TYPE SCARC_MULTIGRID_TYPE

!> --------------------------------------------------------------------------------------------
!> Multigrid type - to be extended for algebraic multigrid
!> --------------------------------------------------------------------------------------------
TYPE SCARC_MGM_TYPE
REAL(EB), DIMENSION(:), ALLOCATABLE :: US, VS, WS          !> velocity components in predictor
REAL(EB), DIMENSION(:), ALLOCATABLE :: U , V , W           !> velocity components in corrector
END TYPE SCARC_MGM_TYPE

!> --------------------------------------------------------------------------------------------
!> Collection of grid level related information
!> --------------------------------------------------------------------------------------------
TYPE SCARC_LEVEL_TYPE

INTEGER :: IOR
INTEGER :: SUBDIVISION(3,-3:3)=0                            !> basic information related to single faces

!> Set stage for used working vectors
TYPE (SCARC_STAGE_TYPE), DIMENSION(2) :: STAGE              !> stage information (hierarchy of solvers)

!> Discretization information
TYPE (SCARC_DISCRET_TYPE) :: SD, UD                   !> Structured or unstructured discretization

!> Administrative structures for obstructions and faces 
TYPE (SCARC_FACE_TYPE), ALLOCATABLE, DIMENSION(:) :: FACE   !> face information
TYPE (SCARC_OBST_TYPE), ALLOCATABLE, DIMENSION(:) :: OBST   !> obstruction information

!> Different dimensions and lengths (for both structured and unstructured discretizations equally required)
INTEGER :: N_DISCRET                                        !> number of discretization types used
INTEGER :: N_OBST                                           !> number of obstructions
INTEGER :: N_WALL_CELLS                                     !> number of wall cells
INTEGER :: N_WALL_CELLS_EXT                                 !> number of external cells
INTEGER :: N_WALL_CELLS_INT                                 !> number of internal cells
INTEGER :: N_CELL_INDEX                                     !> number of entries in CELL_INDEX array

!> Cell and wall index information:
!>    - on the finest level, the original arrays from FDS are used
!>    - separate arrays will only be allocated for coarser levels
!>    - to address them on all levels, corresponding pointers are used 
INTEGER, POINTER, DIMENSION (:,:,:)     :: CELL_INDEX_PTR   !> Pointer to CELL_INDEX 
INTEGER, ALLOCATABLE, DIMENSION (:,:,:) :: CELL_INDEX       !> cell index list (only allocated for coarser levels)
INTEGER, ALLOCATABLE, DIMENSION (:,:,:) :: CELL_STATE       !> state of single cells (gasphase/solid)
INTEGER, POINTER, DIMENSION (:,:)       :: WALL_INDEX_PTR   !> Pointer to WALL_INDEX
INTEGER, ALLOCATABLE, DIMENSION (:,:)   :: WALL_INDEX       !> wall index list (only allocated for coarser levels)

!> Different counters for some grid dimensions and lengths
INTEGER :: NX, NY, NZ                                       !> number of grid cells in x-, y- and z-direction
INTEGER :: NW, NWL                                          !> number of global and local wall cells
INTEGER :: NC, NCG, NCE                                     !> number of cells, ghost cells and extended cells
INTEGER :: NCPL=1, NCPL_MAX=-NSCARC_UNDEFINED_INT           !> number of couplings
INTEGER :: NCPLS, NCPLR                                     !> number of couplings to send and read

!> Different coordinate information
REAL(EB) :: DX , DY , DZ                                    !> step sizes in x-, y- and z-direction
REAL(EB) :: DXI, DYI, DZI                                   !> inversed of step sizes in x-, y- and z-direction
REAL(EB) :: DXI2, DYI2, DZI2                                !> squared and inversed step sizes in x-, y- and z-direction
REAL(EB) :: DH = 0.0_EB                                     !> local step sizes
REAL(EB), ALLOCATABLE, DIMENSION (:) :: XCOR, YCOR, ZCOR    !> coordinate vectors in x-, y- and z-direction
REAL(EB), ALLOCATABLE, DIMENSION (:) :: XMID, YMID, ZMID    !> midpoint vectors in x-, y- and z-direction
REAL(EB), ALLOCATABLE, DIMENSION (:) :: DXL, DYL, DZL       !> step size vectors in x-, y- and z-direction

!> Additonal structures for different solvers/preconditioners
TYPE (SCARC_MULTIGRID_TYPE) :: MG                            !> Multigrid type
TYPE (SCARC_MGM_TYPE) :: MGM                                 !> McKenney-Greengar-Mayo method
TYPE (SCARC_FFT_TYPE) :: FFT                                 !> FFT preconditioner based on CRAYFISHPAK
#ifdef WITH_MKL
TYPE (SCARC_MKL_TYPE) :: MKL                                 !> MKL preconditioner based on Intel MKL 
#endif

END TYPE SCARC_LEVEL_TYPE

!> --------------------------------------------------------------------------------------------
!> Sample sequence of used solvers in stack
!> --------------------------------------------------------------------------------------------
TYPE SCARC_SOLVER_TYPE

CHARACTER(30) :: CNAME = 'NONE'                      !> name of current solver

!> Types of different solver components
INTEGER :: TYPE_METHOD    = NSCARC_UNDEFINED_INT    !> current global method
INTEGER :: TYPE_SOLVER    = NSCARC_UNDEFINED_INT    !> current solver
INTEGER :: TYPE_DISCRET   = NSCARC_UNDEFINED_INT    !> current discretization
INTEGER :: TYPE_PARENT    = NSCARC_UNDEFINED_INT    !> parent (calling) solver
INTEGER :: TYPE_STAGE     = NSCARC_UNDEFINED_INT    !> stage of working vectors 
INTEGER :: TYPE_SCOPE     = NSCARC_UNDEFINED_INT    !> scope of acting 
INTEGER :: TYPE_NLMIN     = NSCARC_UNDEFINED_INT    !> minimum level for that solver
INTEGER :: TYPE_NLMAX     = NSCARC_UNDEFINED_INT    !> maximum level for that solver
INTEGER :: TYPE_RELAX     = NSCARC_UNDEFINED_INT    !> relaxation method
INTEGER :: TYPE_TWOLEVEL  = NSCARC_UNDEFINED_INT    !> two-level method used?
INTEGER :: TYPE_INTERPOL  = NSCARC_UNDEFINED_INT    !> interpolation type
INTEGER :: TYPE_ACCURACY  = NSCARC_UNDEFINED_INT    !> accuracy requirements
INTEGER :: TYPE_PRECISION = NSCARC_UNDEFINED_INT    !> precision type 
INTEGER :: TYPE_CYCLING   = NSCARC_UNDEFINED_INT    !> multigrid cycle

!> References to different vectors which are needed for the current solver
INTEGER :: X = NSCARC_UNDEFINED_INT                !> reference to local X-vector, double precision
INTEGER :: B = NSCARC_UNDEFINED_INT                !> reference to local B-vector, double precision
INTEGER :: D = NSCARC_UNDEFINED_INT                !> reference to local D-vector, double precision
INTEGER :: R = NSCARC_UNDEFINED_INT                !> reference to local R-vector, double precision
INTEGER :: V = NSCARC_UNDEFINED_INT                !> reference to local V-vector, double precision
INTEGER :: Y = NSCARC_UNDEFINED_INT                !> reference to local Y-vector, double precision
INTEGER :: Z = NSCARC_UNDEFINED_INT                !> reference to local Z-vector, double precision
#ifdef WITH_SCARC_DEBUG
INTEGER :: E = NSCARC_UNDEFINED_INT                !> reference to local E-vector, double precision
#endif
#ifdef WITH_MKL_FB
INTEGER :: X_FB = NSCARC_UNDEFINED_INT             !> reference to local X-vector, single precision
INTEGER :: B_FB = NSCARC_UNDEFINED_INT             !> reference to local B-vector, single precision
INTEGER :: D_FB = NSCARC_UNDEFINED_INT             !> reference to local D-vector, single precision
INTEGER :: R_FB = NSCARC_UNDEFINED_INT             !> reference to local R-vector, single precision
INTEGER :: V_FB = NSCARC_UNDEFINED_INT             !> reference to local V-vector, single precision
INTEGER :: Y_FB = NSCARC_UNDEFINED_INT             !> reference to local Y-vector, single precision
INTEGER :: W_FB = NSCARC_UNDEFINED_INT             !> reference to local W-vector, single precision
#endif

!> Converegence requirements for current solver
INTEGER :: NIT   = NSCARC_UNDEFINED_INT            !> maximum iteration number
INTEGER :: ITE   = NSCARC_UNDEFINED_INT            !> current iteration number
REAL(EB) :: EPS   = NSCARC_UNDEFINED_REAL_EB        !> required accuracy
REAL(EB) :: RES   = NSCARC_UNDEFINED_REAL_EB        !> current residual
REAL(EB) :: RESIN = NSCARC_UNDEFINED_REAL_EB        !> initial residual
REAL(EB) :: ERR   = NSCARC_UNDEFINED_REAL_EB        !> initial residual
REAL(EB) :: OMEGA = NSCARC_UNDEFINED_REAL_EB        !> relaxation parameter
REAL(EB) :: CAPPA = NSCARC_UNDEFINED_REAL_EB        !> convergence rate

END TYPE SCARC_SOLVER_TYPE

!> --------------------------------------------------------------------------------------------
!> Stack type
!> --------------------------------------------------------------------------------------------
TYPE SCARC_STACK_TYPE
TYPE (SCARC_SOLVER_TYPE), POINTER :: SOLVER                    !> type of current solver
END TYPE SCARC_STACK_TYPE

!> --------------------------------------------------------------------------------------------
!> Administration other mesh data needed for the coupling of adjacent neighbors
!> --------------------------------------------------------------------------------------------
TYPE SCARC_NEIGHBOR_TYPE

!> Sizes and vectors for data exchange between meshes
INTEGER :: NICMAX_R=0, NICMAX_S=0, NIC_R=0, NIC_S=0
REAL(EB) :: SEND_REAL_BASIC(50), RECV_REAL_BASIC(50)             !> initial real send and receive buffers
INTEGER :: SEND_INT_BASIC(50) , RECV_INT_BASIC(50)              !> initial integer send and receive buffers
REAL (EB), ALLOCATABLE, DIMENSION (:) :: SEND_REAL, RECV_REAL    !> main real send and receive buffers
INTEGER  , ALLOCATABLE, DIMENSION (:) :: SEND_INT, RECV_INT      !> main integer send and receive buffers

!> level dependent information of mesh neighbor
TYPE (SCARC_LEVEL_TYPE),  ALLOCATABLE, DIMENSION(:) :: LEVEL     !> Description of level related information

END TYPE SCARC_NEIGHBOR_TYPE

!> --------------------------------------------------------------------------------------------
!> Basic administration type for ScaRC-method
!> --------------------------------------------------------------------------------------------
TYPE SCARC_TYPE
INTEGER :: N_NEIGHBORS = 0                                           !> number of adjacent neighbors of whole mesh
INTEGER :: N_CELLS = 0                                               !> total number of cells on that mesh
INTEGER :: IBAR, JBAR, KBAR                                          !> number of cells (corresponding to main prg)
REAL(EB) :: XS, XF, YS, YF, ZS, ZF                                    !> x-, y- and z-bounds (corresponding to main prg)
REAL(EB) :: RHS_END = 0.0_EB                                          !> very last RHS entry, needed for matrix condensing
INTEGER, ALLOCATABLE, DIMENSION(:) :: NEIGHBORS                       !> List of adjacent neighbors of whole mesh
TYPE (SCARC_NEIGHBOR_TYPE) , ALLOCATABLE, DIMENSION(:)   :: OSCARC    !> ScaRC type on other mesh
TYPE (SCARC_LEVEL_TYPE)    , ALLOCATABLE, DIMENSION(:)   :: LEVEL     !> of level related information
END TYPE SCARC_TYPE

!> -----------------------------------------------------------------------------------------------
!> Iteration parameters vectors and iteration parameters
!> -----------------------------------------------------------------------------------------------
REAL (EB) :: DT, RDT
REAL (EB) :: EPS, RES, RESIN, OMEGA, CAPPA, ERR
INTEGER :: X, B, D, R, V, Y, Z
#ifdef WITH_SCARC_DEBUG
INTEGER :: E
#endif
INTEGER :: NIT, ITE
INTEGER :: ITE_PRES=0, ITE_TOTAL=0, ITE_CG=0, ITE_MG=0, ITE_LU=0, ITE_SMOOTH=0, ITE_COARSE=0
CHARACTER(60) :: CNAME

!> -----------------------------------------------------------------------------------------------
!> Preliminary definitions for proof of CONPept for CG_Baroclinic routine
!> -----------------------------------------------------------------------------------------------
REAL(EB) :: T
INTEGER :: N_REQ7
INTEGER, ALLOCATABLE, DIMENSION(:) :: REQ7

!> --------------------------------------------------------------------------------------------
!> Basic ScaRC type, different solver types and stack of used solvers
!> --------------------------------------------------------------------------------------------
TYPE (SCARC_TYPE), SAVE, DIMENSION(:), ALLOCATABLE, TARGET :: SCARC

TYPE (SCARC_SOLVER_TYPE) , SAVE, TARGET :: MAIN_CG, MAIN_CG_STRUCTURED, MAIN_CG_UNSTRUCTURED
TYPE (SCARC_SOLVER_TYPE) , SAVE, TARGET :: MAIN_GMG
TYPE (SCARC_SOLVER_TYPE) , SAVE, TARGET :: COARSE_KRYLOV
#ifdef WITH_MKL
TYPE (SCARC_SOLVER_TYPE) , SAVE, TARGET :: MAIN_LU      
TYPE (SCARC_SOLVER_TYPE) , SAVE, TARGET :: COARSE_CLUSTER, COARSE_PARDISO
#endif
TYPE (SCARC_SOLVER_TYPE) , SAVE, TARGET :: PRECON_JACOBI, PRECON_SSOR, PRECON_ILU, PRECON_FFT, PRECON_GMG
TYPE (SCARC_SOLVER_TYPE) , SAVE, TARGET :: PRECON_GSM, PRECON_SGSM, PRECON_SORM, PRECON_SSORM
TYPE (SCARC_SOLVER_TYPE) , SAVE, TARGET :: SMOOTH_JACOBI, SMOOTH_SSOR, SMOOTH_FFT
#ifdef WITH_MKL
TYPE (SCARC_SOLVER_TYPE) , SAVE, TARGET :: PRECON_MKL, SMOOTH_MKL
#endif
TYPE (SCARC_STACK_TYPE)  , SAVE, DIMENSION(:), ALLOCATABLE, TARGET :: STACK
TYPE (SCARC_TIME_TYPE)   , SAVE, DIMENSION(:), ALLOCATABLE :: TSETUP, TSUM, TSTEP
TYPE (SCARC_MESSAGE_TYPE), SAVE :: MSG

CONTAINS


!> ------------------------------------------------------------------------------------------------
!> Initialize ScaRC structures based on SCARC-input parameters from &PRES namelist
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP
REAL(EB) :: TNOW

TNOW = CURRENT_TIME()

!> Setup messaging/debugging mechanisms and time measurements
CALL SCARC_SETUP_MESSAGING                                                 
CALL SCARC_SETUP_TIMING                                                    

!> Parse all ScaRC parameters which have been read in read.f90
CALL SCARC_PARSE_INPUT                   ; IF (STOP_STATUS==SETUP_STOP) RETURN     

!> Setup information on grid levels, ScaRC-types, discretization, neighbors, global variables and faces
CALL SCARC_SETUP_LEVELS                  ; IF (STOP_STATUS==SETUP_STOP) RETURN     
CALL SCARC_SETUP_TYPES                   ; IF (STOP_STATUS==SETUP_STOP) RETURN     
CALL SCARC_SETUP_DISCRETIZATION          ; IF (STOP_STATUS==SETUP_STOP) RETURN     
CALL SCARC_SETUP_NEIGHBORS               ; IF (STOP_STATUS==SETUP_STOP) RETURN     
CALL SCARC_SETUP_GLOBALS                 ; IF (STOP_STATUS==SETUP_STOP) RETURN       
CALL SCARC_SETUP_FACES                   ; IF (STOP_STATUS==SETUP_STOP) RETURN              

!> Setup wall information according to specified discretization type/method
IF (HAS_MULTIPLE_DISCRETIZATIONS) THEN
   TYPE_DISCRET = NSCARC_DISCRET_UNSTRUCTURED                             
   CALL SCARC_SETUP_WALLS                ; IF (STOP_STATUS==SETUP_STOP) RETURN       
   TYPE_DISCRET = NSCARC_DISCRET_STRUCTURED                               
   CALL SCARC_SETUP_WALLS                ; IF (STOP_STATUS==SETUP_STOP) RETURN       
ELSE                                                                      
   CALL SCARC_SETUP_WALLS                ; IF (STOP_STATUS==SETUP_STOP) RETURN       
ENDIF

!> Setup information for data exchange, matrix systems, used solvers and vectors
CALL SCARC_SETUP_EXCHANGE                ; IF (STOP_STATUS==SETUP_STOP) RETURN       
CALL SCARC_SETUP_SYSTEM                  ; IF (STOP_STATUS==SETUP_STOP) RETURN       
CALL SCARC_SETUP_METHODS                 ; IF (STOP_STATUS==SETUP_STOP) RETURN       
CALL SCARC_SETUP_VECTORS                 ; IF (STOP_STATUS==SETUP_STOP) RETURN       

TSETUP(MYID+1)%OVERALL = TSETUP(MYID+1)%OVERALL + CURRENT_TIME() - TNOW
END SUBROUTINE SCARC_SETUP


!> ------------------------------------------------------------------------------------------------
!> Setup time measurements
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_TIMING

ALLOCATE (TSETUP(LOWER_MESH_INDEX:UPPER_MESH_INDEX), STAT = IERROR)
CALL CHKMEMERR ('SCARC_SETUP_TIMING', 'TSETUP', IERROR)

ALLOCATE (TSTEP(LOWER_MESH_INDEX:UPPER_MESH_INDEX), STAT = IERROR)
CALL CHKMEMERR ('SCARC_SETUP_TIMING', 'TSTEP', IERROR)

ALLOCATE (TSUM(LOWER_MESH_INDEX:UPPER_MESH_INDEX), STAT = IERROR)
CALL CHKMEMERR ('SCARC_SETUP_TIMING', 'TSUM', IERROR)

END SUBROUTINE SCARC_SETUP_TIMING


!> ------------------------------------------------------------------------------------------------
!> Setup debug file if requested
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_MESSAGING
#if defined(WITH_SCARC_VERBOSE) || defined(WITH_SCARC_DEBUG)
INTEGER :: NM, LASTID
#endif

IF (TRIM(SCARC_ERROR_FILE) == '.TRUE.') HAS_CSV = .TRUE.

!> If requested, open file for CSV-information about convergence of different solvers
IF (HAS_CSV) THEN
   IF (MYID == 0) THEN
      WRITE (MSG%FILE_STAT, '(A,A)') TRIM(CHID),'_scarc.csv'
      MSG%LU_STAT = GET_FILE_NUMBER()
      OPEN (MSG%LU_STAT, FILE=MSG%FILE_STAT)
      WRITE(MSG%LU_STAT,*) '  #Pres,   Stack,  #ScaRC,     #CG,     #MG,   Level, #Smooth, SmoType, #Coarse,', &
                           '     #LU,    Residual,   Cappa'
   ENDIF
ENDIF

#ifdef WITH_SCARC_VERBOSE
!> If requested, open file for log-information
LASTID = -99999
DO NM=LOWER_MESH_INDEX, UPPER_MESH_INDEX
   IF (MYID == LASTID) CYCLE
   WRITE (MSG%FILE_VERBOSE, '(A,A,i3.3)') TRIM(CHID),'.log',MYID+1
   MSG%LU_VERBOSE = GET_FILE_NUMBER()
   OPEN (MSG%LU_VERBOSE, FILE=MSG%FILE_VERBOSE, ACTION = 'readwrite')
   LASTID = MYID
ENDDO
#endif 

#ifdef WITH_SCARC_DEBUG
!> If requested, open file for debug messages
LASTID = -99999
DO NM=LOWER_MESH_INDEX, UPPER_MESH_INDEX
   IF (MYID == LASTID) CYCLE
   WRITE (MSG%FILE_DEBUG, '(A,A,i3.3)') TRIM(CHID),'.debug',MYID+1
   MSG%LU_DEBUG = GET_FILE_NUMBER()
   OPEN (MSG%LU_DEBUG, FILE=MSG%FILE_DEBUG, ACTION = 'readwrite')
   LASTID = MYID
ENDDO
#endif 

END SUBROUTINE SCARC_SETUP_MESSAGING


!> ------------------------------------------------------------------------------------------------
!> Shutdown ScaRC with error message
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SHUTDOWN(NERROR, CPARAM, NPARAM)
CHARACTER(*), INTENT(IN) :: CPARAM
INTEGER, INTENT(IN) :: NERROR, NPARAM
CHARACTER(80) :: CERROR

!>
!> Assign according error message
!>
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
   CASE (NSCARC_ERROR_BOUNDARY_TYPE)
      CERROR = 'Wrong boundary type'
   CASE (NSCARC_ERROR_NEIGHBOR_TYPE)
      CERROR = 'Wrong neighbor'
   CASE (NSCARC_ERROR_NEIGHBOR_NUMBER)
      CERROR = 'More than 20 neighbors along one face not allowed'
   CASE (NSCARC_ERROR_MATRIX_SYMMETRY)
      CERROR = 'Matrix not symmetric for mesh'
   CASE (NSCARC_ERROR_MATRIX_SUBDIAG)
      CERROR = 'Subdiagonal missing for system matrix'
   CASE (NSCARC_ERROR_MATRIX_SETUP)
      CERROR = 'Matrix setup failed for level type'
   CASE (NSCARC_ERROR_MATRIX_SIZE)
      CERROR = 'Matrix resized failed due to too big new length'
   CASE (NSCARC_ERROR_MATRIX_COPY)
      CERROR = 'Matrix copy failed due to too already existing array'
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
   CASE (NSCARC_ERROR_AMG_MISSING)
      CERROR = 'Algebraic multigrid is currently missing'
   CASE (NSCARC_ERROR_VECTOR_LENGTH)
      CERROR = 'Inconsistent length for vector allocation'
   CASE (NSCARC_ERROR_POC_STOP)
      CERROR = 'Only one call of solver due to proof of concept'
   CASE (NSCARC_ERROR_BICG_DISABLED)
      CERROR = 'Krylov solver BICG temporarily disabled'
END SELECT

IF (CPARAM /= SCARC_NONE) THEN
   IF (MYID == 0) WRITE(LU_ERR,1000)  CERROR, CPARAM, TRIM(CHID)
ELSE IF (NPARAM /= NSCARC_NONE) THEN
   IF (MYID == 0) WRITE(LU_ERR,2000)  CERROR, NPARAM, TRIM(CHID)
ELSE
   IF (MYID == 0) WRITE(LU_ERR,3000)  CERROR, TRIM(CHID)
ENDIF

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


!> ----------------------------------------------------------------------------------------------------
!> Determine types of input parameters
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_PARSE_INPUT

ITERATE_PRESSURE = .TRUE.  ! Although there is no need to do pressure iterations to drive down velocity error
                           ! leave it .TRUE. to write out velocity error diagnostics.

!>
!> ------------- Set type of discretization
!>
SELECT CASE (TRIM(SCARC_DISCRETIZATION))
   CASE ('STRUCTURED')
      TYPE_DISCRET = NSCARC_DISCRET_STRUCTURED
   CASE ('UNSTRUCTURED')
      TYPE_DISCRET = NSCARC_DISCRET_UNSTRUCTURED
   CASE DEFAULT
      CALL SCARC_SHUTDOWN(NSCARC_ERROR_PARSE_INPUT, SCARC_DISCRETIZATION, NSCARC_NONE)
END SELECT

!>
!> ------------ Set type of matrix storage (COMPACT/BANDED/FREE)
!>
SELECT CASE (TRIM(SCARC_MATRIX_FORMAT))
   CASE ('COMPACT')
      TYPE_MATRIX = NSCARC_MATRIX_COMPACT
   CASE ('BANDED')
      TYPE_MATRIX = NSCARC_MATRIX_BANDED
   CASE ('FREE')
      TYPE_MATRIX = NSCARC_MATRIX_FREE
   CASE DEFAULT
      CALL SCARC_SHUTDOWN(NSCARC_ERROR_PARSE_INPUT, SCARC_MATRIX_FORMAT, NSCARC_NONE)
END SELECT

!>
!> ------------ Set type of global solver
!>
SELECT CASE (TRIM(SCARC_METHOD))

   !> ------------------------- McKenny-Greengard-Mayo solver -------------------------
   CASE ('MGM')

      !> Just preset some values for proof of concept
      TYPE_METHOD       = NSCARC_METHOD_MGM
      TYPE_KRYLOV       = NSCARC_KRYLOV_CG
      TYPE_TWOLEVEL     = NSCARC_TWOLEVEL_NONE
      TYPE_PRECON       = NSCARC_RELAX_FFT
      TYPE_PRECON_SCOPE = NSCARC_SCOPE_LOCAL
     

   !> ------------------------- Global Krylov solver ----------------------------------
   CASE ('KRYLOV')

      TYPE_METHOD = NSCARC_METHOD_KRYLOV

      !> Set type of Krylov-method (CG/BICG)
      SELECT CASE (TRIM(SCARC_KRYLOV))
         CASE ('CG')
            TYPE_KRYLOV = NSCARC_KRYLOV_CG
#ifdef WITH_SCARC_CGBARO
         CASE ('CGBARO')
            TYPE_KRYLOV = NSCARC_KRYLOV_CGBARO
#endif
         CASE ('BICG')
            TYPE_KRYLOV = NSCARC_KRYLOV_BICG
            CALL SCARC_SHUTDOWN(NSCARC_ERROR_BICG_DISABLED, SCARC_KRYLOV, NSCARC_NONE)
         CASE DEFAULT
            CALL SCARC_SHUTDOWN(NSCARC_ERROR_PARSE_INPUT, SCARC_KRYLOV, NSCARC_NONE)
      END SELECT

      !> Set type of two-level method
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
         CASE DEFAULT
            CALL SCARC_SHUTDOWN(NSCARC_ERROR_PARSE_INPUT, SCARC_TWOLEVEL, NSCARC_NONE)
      END SELECT

      !> Set type of interpolation for two-level Krylov method
      SELECT CASE (TRIM(SCARC_KRYLOV_INTERPOL))
         CASE ('NONE')
            TYPE_INTERPOL = NSCARC_UNDEFINED_INT
         CASE ('CONSTANT')
            TYPE_INTERPOL = NSCARC_INTERPOL_CONSTANT
         CASE ('BILINEAR')
            TYPE_INTERPOL = NSCARC_INTERPOL_BILINEAR
         CASE DEFAULT
            CALL SCARC_SHUTDOWN(NSCARC_ERROR_PARSE_INPUT, SCARC_KRYLOV_INTERPOL, NSCARC_NONE)
      END SELECT

      !> Set type of preconditioner (JACOBI/SSOR/GSM/SGSM/SORM/SSORM/ILU/FFT/GMG/PARDISO/CLUSTER)
      SELECT CASE (TRIM(SCARC_PRECON))
         CASE ('JACOBI')                                    !> Jacobi preconditioner
            TYPE_PRECON = NSCARC_RELAX_JACOBI
         CASE ('SSOR')                                      !> Symmetric SOR preconditioner
            TYPE_PRECON = NSCARC_RELAX_SSOR
         CASE ('GSM')                                       !> Gauss-Seidel preconditioner in matrix form
            TYPE_PRECON = NSCARC_RELAX_GSM
         CASE ('SGSM')                                      !> Symmetric Gauss-Seidel preconditioner in matrix form
            TYPE_PRECON = NSCARC_RELAX_SGSM
         CASE ('SORM')                                      !> SOR preconditioner in matrix form
            TYPE_PRECON = NSCARC_RELAX_SORM
         CASE ('SSORM')                                     !> Symmetric SOR preconditioner in matrix form
            TYPE_PRECON = NSCARC_RELAX_SSORM
         CASE ('ILU')                                       !> ILU preconditioner
            TYPE_PRECON = NSCARC_RELAX_ILU
         CASE ('GMG')                                       !> Geometric multigrid preconditioner
            TYPE_PRECON = NSCARC_RELAX_GMG
            SELECT CASE (TRIM(SCARC_SMOOTH))
               CASE ('JACOBI')
                  TYPE_SMOOTH = NSCARC_RELAX_JACOBI
               CASE ('SSOR')
                  TYPE_SMOOTH = NSCARC_RELAX_SSOR
               CASE ('ILU')
                  TYPE_SMOOTH = NSCARC_RELAX_ILU
               CASE ('FFT')
                  IF (TYPE_DISCRET == NSCARC_DISCRET_UNSTRUCTURED) &
                     CALL SCARC_SHUTDOWN(NSCARC_ERROR_FFT_DISCRET, SCARC_NONE, NSCARC_NONE)
                  TYPE_PRECON = NSCARC_RELAX_FFT
               CASE ('PARDISO')
#ifdef WITH_MKL
                  IF (TYPE_MATRIX == NSCARC_MATRIX_COMPACT) THEN
                     TYPE_SMOOTH = NSCARC_RELAX_MKL
                  ELSE
                     CALL SCARC_SHUTDOWN(NSCARC_ERROR_MKL_PARDISO, SCARC_NONE, NSCARC_NONE)
                  ENDIF
#else
                  CALL SCARC_SHUTDOWN(NSCARC_ERROR_MKL_PARDISO, SCARC_NONE, NSCARC_NONE)
#endif

               CASE ('CLUSTER')
#ifdef WITH_MKL
                  IF (TYPE_MATRIX == NSCARC_MATRIX_COMPACT) THEN
                     TYPE_SMOOTH = NSCARC_RELAX_MKL
                  ELSE
                     CALL SCARC_SHUTDOWN(NSCARC_ERROR_MKL_STORAGE, SCARC_NONE, NSCARC_NONE)
                  ENDIF
#else
                  CALL SCARC_SHUTDOWN(NSCARC_ERROR_MKL_STORAGE, SCARC_NONE, NSCARC_NONE)
#endif
            END SELECT
         CASE ('FFT')                                                !> FFT preconditioner
            IF (TYPE_DISCRET == NSCARC_DISCRET_UNSTRUCTURED) &
               CALL SCARC_SHUTDOWN(NSCARC_ERROR_FFT_DISCRET, SCARC_NONE, NSCARC_NONE)
            TYPE_PRECON = NSCARC_RELAX_FFT
         CASE ('PARDISO')                                            !> LU preconditioner based on MKL-PARDISO
#ifdef WITH_MKL
            IF (TYPE_MATRIX == NSCARC_MATRIX_COMPACT) THEN
               TYPE_PRECON   = NSCARC_RELAX_MKL
               TYPE_MKL      = NSCARC_MKL_LOCAL
            ELSE
               CALL SCARC_SHUTDOWN(NSCARC_ERROR_MKL_STORAGE, SCARC_NONE, NSCARC_NONE)
            ENDIF                                    
#else
            CALL SCARC_SHUTDOWN(NSCARC_ERROR_MKL_PARDISO, SCARC_NONE, NSCARC_NONE)
#endif
         CASE ('CLUSTER')                            !>  LU-preconditioner based on MKL Cluster_Sparse_Solver
#ifdef WITH_MKL
            IF (TYPE_MATRIX == NSCARC_MATRIX_COMPACT) THEN
               TYPE_PRECON   = NSCARC_RELAX_MKL
               TYPE_MKL      = NSCARC_MKL_GLOBAL
            ELSE
               CALL SCARC_SHUTDOWN(NSCARC_ERROR_MKL_STORAGE, SCARC_NONE, NSCARC_NONE)
            ENDIF
#else
            CALL SCARC_SHUTDOWN(NSCARC_ERROR_MKL_CLUSTER, SCARC_NONE, NSCARC_NONE)
#endif
         CASE DEFAULT
            CALL SCARC_SHUTDOWN(NSCARC_ERROR_PARSE_INPUT, SCARC_PRECON, NSCARC_NONE)
      END SELECT

      !> set type scope for preconditioner (GLOBAL/LOCAL)
      SELECT CASE (TRIM(SCARC_PRECON_SCOPE))
         CASE ('GLOBAL')
            TYPE_PRECON_SCOPE = NSCARC_SCOPE_GLOBAL
         CASE ('LOCAL')
            TYPE_PRECON_SCOPE = NSCARC_SCOPE_LOCAL
         CASE DEFAULT
            CALL SCARC_SHUTDOWN(NSCARC_ERROR_PARSE_INPUT, SCARC_PRECON_SCOPE, NSCARC_NONE)
      END SELECT

   !> ------------------------- Global geometric multigrid solver -------------------------------
   CASE ('MULTIGRID')

      TYPE_METHOD = NSCARC_METHOD_MULTIGRID

      !> Set type of multigrid method (GEOMETRIC/ALGEBRAIC)
      SELECT CASE (TRIM(SCARC_MULTIGRID))
         CASE ('GEOMETRIC')
            TYPE_MULTIGRID = NSCARC_MULTIGRID_GEOMETRIC
         CASE ('ALGEBRAIC')
            TYPE_MULTIGRID = NSCARC_MULTIGRID_ALGEBRAIC
         CASE DEFAULT
            CALL SCARC_SHUTDOWN(NSCARC_ERROR_PARSE_INPUT, SCARC_MULTIGRID, NSCARC_NONE)
      END SELECT

      !> Set type of smoother (JACOBI/SGS/SSOR/SSORM/ILU/PARDISO/CLUSTER)
      SELECT CASE (TRIM(SCARC_SMOOTH))                        !> use same parameters as for preconditioner
         CASE ('JACOBI')
            TYPE_SMOOTH = NSCARC_RELAX_JACOBI
         CASE ('SSOR')
            TYPE_SMOOTH = NSCARC_RELAX_SSOR
         CASE ('ILU')
            TYPE_SMOOTH = NSCARC_RELAX_ILU
         CASE ('FFT')
            IF (TYPE_DISCRET == NSCARC_DISCRET_UNSTRUCTURED) &
               CALL SCARC_SHUTDOWN(NSCARC_ERROR_FFT_DISCRET, SCARC_NONE, NSCARC_NONE)
            TYPE_SMOOTH = NSCARC_RELAX_FFT
         CASE ('PARDISO')
#ifdef WITH_MKL
            IF (TYPE_MATRIX == NSCARC_MATRIX_COMPACT) THEN
               TYPE_SMOOTH = NSCARC_RELAX_MKL
            ELSE
               CALL SCARC_SHUTDOWN(NSCARC_ERROR_MKL_STORAGE, SCARC_NONE, NSCARC_NONE)
            ENDIF
#else
            CALL SCARC_SHUTDOWN(NSCARC_ERROR_MKL_PARDISO, SCARC_NONE, NSCARC_NONE)
#endif
         CASE ('CLUSTER')
#ifdef WITH_MKL
            IF (TYPE_MATRIX == NSCARC_MATRIX_COMPACT) THEN
               TYPE_SMOOTH = NSCARC_RELAX_MKL
            ELSE
               CALL SCARC_SHUTDOWN(NSCARC_ERROR_MKL_STORAGE, SCARC_NONE, NSCARC_NONE)
            ENDIF
#else
            CALL SCARC_SHUTDOWN(NSCARC_ERROR_MKL_CLUSTER, SCARC_NONE, NSCARC_NONE)
#endif
         CASE DEFAULT
            CALL SCARC_SHUTDOWN(NSCARC_ERROR_PARSE_INPUT, SCARC_SMOOTH, NSCARC_NONE)
      END SELECT

      !> set type scope for smoother (GLOBAL/LOCAL)
      SELECT CASE (TRIM(SCARC_SMOOTH_SCOPE))
         CASE ('GLOBAL')
            TYPE_SMOOTH_SCOPE = NSCARC_SCOPE_GLOBAL
         CASE ('LOCAL')
            TYPE_SMOOTH_SCOPE = NSCARC_SCOPE_LOCAL
         CASE DEFAULT
            CALL SCARC_SHUTDOWN(NSCARC_ERROR_PARSE_INPUT, SCARC_SMOOTH_SCOPE, NSCARC_NONE)
      END SELECT
WRITE(MSG%LU_DEBUG,*) 'C: TYPE_SMOOTH_SCOPE =',TYPE_SMOOTH_SCOPE

   !> ------------------------- Global LU-decomposition solver -------------------------------
#ifdef WITH_MKL
   CASE ('MKL')

      TYPE_METHOD  = NSCARC_METHOD_LU      

      !> Set type of MKL method (global/local)
      SELECT CASE (TRIM(SCARC_MKL))                      !Achtung, hier noch nacharbeiten!
         CASE ('GLOBAL')
#ifdef WITH_MKL
            IF (TYPE_MATRIX == NSCARC_MATRIX_COMPACT) THEN
               TYPE_MKL     = NSCARC_MKL_GLOBAL
            ELSE
               CALL SCARC_SHUTDOWN(NSCARC_ERROR_MKL_STORAGE, SCARC_NONE, NSCARC_NONE)
            ENDIF
#else
            CALL SCARC_SHUTDOWN(NSCARC_ERROR_MKL_CLUSTER, SCARC_NONE, NSCARC_NONE)
#endif
         CASE ('LOCAL')
#ifdef WITH_MKL
            IF (TYPE_MATRIX == NSCARC_MATRIX_COMPACT) THEN
               TYPE_MKL      = NSCARC_MKL_LOCAL
            ELSE
               CALL SCARC_SHUTDOWN(NSCARC_ERROR_MKL_STORAGE, SCARC_NONE, NSCARC_NONE)
            ENDIF
#else
            CALL SCARC_SHUTDOWN(NSCARC_ERROR_MKL_PARDISO, SCARC_NONE, NSCARC_NONE)
#endif
         CASE DEFAULT
            CALL SCARC_SHUTDOWN(NSCARC_ERROR_PARSE_INPUT, SCARC_MKL, NSCARC_NONE)
      END SELECT
#endif

   CASE DEFAULT
      CALL SCARC_SHUTDOWN(NSCARC_ERROR_PARSE_INPUT, SCARC_METHOD, NSCARC_NONE)

END SELECT

!>
!> If a multigrid solver is used (either as main solver or as preconditioner)
!> set types for multigrid, coarse grid solver and cycling pattern
!>
IF (TYPE_METHOD == NSCARC_METHOD_MULTIGRID .OR. TYPE_PRECON == NSCARC_RELAX_GMG) THEN

   !> Set type of multigrid (GEOMETRIC/ALGEBRAIC with corresponding coarsening strategy)
   SELECT CASE (TRIM(SCARC_MULTIGRID))

      CASE ('GEOMETRIC')
         TYPE_MULTIGRID = NSCARC_MULTIGRID_GEOMETRIC

      CASE ('ALGEBRAIC')
         TYPE_MULTIGRID = NSCARC_MULTIGRID_ALGEBRAIC

         !> Set type of coarsening strategy (STANDARD/AGGRESSIVE)
         SELECT CASE (TRIM(SCARC_MULTIGRID_COARSENING))
            CASE ('BASIC')
               TYPE_COARSENING = NSCARC_COARSENING_BASIC
            CASE ('FALGOUT')
               TYPE_COARSENING = NSCARC_COARSENING_FALGOUT
            CASE ('RS3')
               TYPE_COARSENING = NSCARC_COARSENING_RS3
            CASE ('A1')
               TYPE_COARSENING = NSCARC_COARSENING_A1
            CASE ('A2')
               TYPE_COARSENING = NSCARC_COARSENING_A2
            CASE DEFAULT
            CALL SCARC_SHUTDOWN(NSCARC_ERROR_PARSE_INPUT, SCARC_MULTIGRID_COARSENING, NSCARC_NONE)
         END SELECT

      CASE DEFAULT
         CALL SCARC_SHUTDOWN(NSCARC_ERROR_PARSE_INPUT, SCARC_MULTIGRID, NSCARC_NONE)
   END SELECT

   !> Set type of cycling pattern (F/V/W)
   SELECT CASE (TRIM(SCARC_MULTIGRID_CYCLE))
      CASE ('F')
         TYPE_CYCLING = NSCARC_CYCLING_F
      CASE ('V')
         TYPE_CYCLING = NSCARC_CYCLING_V
      CASE ('W')
         TYPE_CYCLING = NSCARC_CYCLING_W
      CASE DEFAULT
         CALL SCARC_SHUTDOWN(NSCARC_ERROR_PARSE_INPUT, SCARC_MULTIGRID_CYCLE, NSCARC_NONE)
   END SELECT

   !> Set type of interpolation (STANDARD/DIRECT/MULTIPASS)
   SELECT CASE (TRIM(SCARC_MULTIGRID_INTERPOL))
      CASE ('STANDARD')
         TYPE_INTERPOL = NSCARC_INTERPOL_STANDARD
      CASE ('CONSTANT')
         TYPE_INTERPOL = NSCARC_INTERPOL_CONSTANT
      CASE ('BILINEAR')
         TYPE_INTERPOL = NSCARC_INTERPOL_BILINEAR
      CASE ('CLASSICAL')
         TYPE_INTERPOL = NSCARC_INTERPOL_CLASSICAL
      CASE ('CLASSICAL2')
         TYPE_INTERPOL = NSCARC_INTERPOL_CLASSICAL2
      CASE ('DIRECT')
         TYPE_INTERPOL = NSCARC_INTERPOL_DIRECT
      CASE ('DIRECT_BDRY')
         TYPE_INTERPOL = NSCARC_INTERPOL_DIRECT_BDRY
      CASE DEFAULT
         CALL SCARC_SHUTDOWN(NSCARC_ERROR_PARSE_INPUT, SCARC_MULTIGRID_INTERPOL, NSCARC_NONE)
   END SELECT

ENDIF

!> Set type of coarse grid solver
SELECT CASE (TRIM(SCARC_COARSE))
   CASE ('ITERATIVE')
      TYPE_COARSE = NSCARC_COARSE_ITERATIVE
      TYPE_KRYLOV = NSCARC_KRYLOV_CG
   CASE ('DIRECT')
#ifdef WITH_MKL
      IF (TYPE_MATRIX == NSCARC_MATRIX_COMPACT) THEN
         TYPE_COARSE  = NSCARC_COARSE_DIRECT
         TYPE_MKL     = NSCARC_MKL_COARSE
      ELSE
         CALL SCARC_SHUTDOWN(NSCARC_ERROR_MKL_STORAGE, SCARC_NONE, NSCARC_NONE)
      ENDIF
#else
      CALL SCARC_SHUTDOWN(NSCARC_ERROR_MKL_CLUSTER, SCARC_NONE, NSCARC_NONE)
#endif
   CASE DEFAULT
      CALL SCARC_SHUTDOWN(NSCARC_ERROR_PARSE_INPUT, SCARC_COARSE, NSCARC_NONE)
END SELECT


!>
!> Set type of accuracy (ABSOLUTE/RELATIVE)
!>
SELECT CASE (TRIM(SCARC_ACCURACY))
   CASE ('ABSOLUTE')
      TYPE_ACCURACY = NSCARC_ACCURACY_ABSOLUTE
   CASE ('RELATIVE')
      TYPE_ACCURACY = NSCARC_ACCURACY_RELATIVE
   CASE DEFAULT
      CALL SCARC_SHUTDOWN(NSCARC_ERROR_PARSE_INPUT, SCARC_ACCURACY, NSCARC_NONE)
END SELECT

!>
!> Set type of precision for preconditioner (SINGLE/DOUBLE)
!>
SELECT CASE (TRIM(SCARC_PRECISION))
   CASE ('SINGLE')
      TYPE_PRECISION = NSCARC_PRECISION_FB
   CASE ('DOUBLE')
      TYPE_PRECISION = NSCARC_PRECISION_EB
   CASE DEFAULT
      CALL SCARC_SHUTDOWN(NSCARC_ERROR_PARSE_INPUT, SCARC_PRECISION, NSCARC_NONE)
END SELECT

!>
!> define some logical variables - just for notational convenience
!>
IS_STRUCTURED   = TYPE_DISCRET == NSCARC_DISCRET_STRUCTURED
IS_UNSTRUCTURED = TYPE_DISCRET == NSCARC_DISCRET_UNSTRUCTURED

IS_CG      = TYPE_METHOD == NSCARC_METHOD_KRYLOV
IS_CG_GMG  = IS_CG .AND. TYPE_PRECON == NSCARC_RELAX_GMG .AND. TYPE_MULTIGRID == NSCARC_MULTIGRID_GEOMETRIC

IS_MG  = TYPE_METHOD == NSCARC_METHOD_MULTIGRID
IS_GMG = IS_MG .AND. TYPE_MULTIGRID == NSCARC_MULTIGRID_GEOMETRIC

IS_CG_ADD    = HAS_TWO_LEVELS .AND. TYPE_TWOLEVEL == NSCARC_TWOLEVEL_ADD
IS_CG_MUL    = HAS_TWO_LEVELS .AND. TYPE_TWOLEVEL == NSCARC_TWOLEVEL_MUL
IS_CG_COARSE = HAS_TWO_LEVELS .AND. TYPE_TWOLEVEL == NSCARC_TWOLEVEL_COARSE

HAS_TWO_LEVELS   = IS_CG .AND. TYPE_PRECON /= NSCARC_RELAX_GMG .AND.  TYPE_TWOLEVEL > NSCARC_TWOLEVEL_NONE
HAS_MULTIPLE_LEVELS = IS_GMG .OR. IS_CG_GMG .OR. HAS_TWO_LEVELS

IS_FFT =  TYPE_PRECON == NSCARC_RELAX_FFT .OR.  TYPE_SMOOTH == NSCARC_RELAX_FFT
IS_MKL = (TYPE_PRECON >= NSCARC_RELAX_MKL) .OR. (TYPE_SMOOTH >= NSCARC_RELAX_MKL) .OR. &
         (HAS_MULTIPLE_LEVELS .AND. TYPE_COARSE == NSCARC_COARSE_DIRECT)

END SUBROUTINE SCARC_PARSE_INPUT


!> ------------------------------------------------------------------------------------------------
!> Determine number of grid levels  (1 for CG/BICG-method, NLEVEL for MG-method)
!> Note: NLEVEL_MIN corresponds to finest grid resolution, NLEVEL_MAX to coarsest resolution
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_LEVELS
#ifdef WITH_MKL
INTEGER :: NL
#endif

SELECT_METHOD: SELECT CASE (TYPE_METHOD)

   !>
   !> ----------------------- Global McKenney-Greengard-Mayo method - only finest level --------------------
   !>
   CASE (NSCARC_METHOD_MGM)

      CALL SCARC_GET_NUMBER_OF_LEVELS(NSCARC_LEVEL_SINGLE) 

   !>
   !> ----------------------- Global data-parallel Krylov method --------------------------------------------
   !>
   CASE (NSCARC_METHOD_KRYLOV)

      SELECT_KRYLOV_PRECON: SELECT CASE (TYPE_PRECON)

#ifdef WITH_MKL
         !>
         !> Preconditioning by defect correction based on LU-decomposition
         !> if two-level method, also use coarse grid level, otherwise only use single (finest) grid level
         !>
         CASE (NSCARC_RELAX_MKL)                                

            IF (HAS_TWO_LEVELS) THEN                                      
               CALL SCARC_GET_NUMBER_OF_LEVELS(NSCARC_LEVEL_MULTI)  
               IF (TYPE_COARSE == NSCARC_COARSE_DIRECT) TYPE_MKL_LEVEL(NLEVEL_MAX) = NSCARC_MKL_GLOBAL
            ELSE                                                   
               CALL SCARC_GET_NUMBER_OF_LEVELS(NSCARC_LEVEL_SINGLE)
            ENDIF

            !> Either using a global CLUSTER_SPARSE_SOLVER or local PARDISO solvers from MKL
            SELECT CASE (TYPE_PRECON_SCOPE)
               CASE(NSCARC_SCOPE_GLOBAL) 
                  TYPE_MKL_LEVEL(NLEVEL_MIN) = NSCARC_MKL_GLOBAL
               CASE(NSCARC_SCOPE_LOCAL) 
                  TYPE_MKL_LEVEL(NLEVEL_MIN) = NSCARC_MKL_LOCAL
            END SELECT

#endif

         !>
         !> Preconditioning by defect correction based on geometric multigrid method,
         !> use specified hierarchy of grid levels
         !>
         CASE (NSCARC_RELAX_GMG)                             
            CALL SCARC_GET_NUMBER_OF_LEVELS(NSCARC_LEVEL_MULTI)     
#ifdef WITH_MKL
            IF (TYPE_COARSE == NSCARC_COARSE_DIRECT) TYPE_MKL_LEVEL(NLEVEL_MAX) = NSCARC_MKL_GLOBAL
#endif

         !>
         !> Preconditioning by defect correction based on local basic iterations (JACOBI/SSOR),
         !> if two-level method, also use coarse grid, otherwise only use single (finest) grid level
         !>
         CASE DEFAULT
            IF (HAS_TWO_LEVELS) THEN                                     
               CALL SCARC_GET_NUMBER_OF_LEVELS(NSCARC_LEVEL_MULTI)
#ifdef WITH_MKL
               IF (TYPE_COARSE == NSCARC_COARSE_DIRECT) TYPE_MKL_LEVEL(NLEVEL_MAX) = NSCARC_MKL_GLOBAL
#endif
            ELSE                                       
               CALL SCARC_GET_NUMBER_OF_LEVELS(NSCARC_LEVEL_SINGLE) 
            ENDIF

      END SELECT SELECT_KRYLOV_PRECON

   !>
   !> ----------------------- Global data-parallel Multigrid method -------------------------------------
   !>
   CASE (NSCARC_METHOD_MULTIGRID)

      SELECT_MULTIGRID_TYPE: SELECT CASE (TYPE_MULTIGRID)

         !> 
         !> Use of geometric multigrid method
         !> If not specified by user, determine number of possible grid levels
         !> 
         CASE (NSCARC_MULTIGRID_GEOMETRIC)
            CALL SCARC_GET_NUMBER_OF_LEVELS(NSCARC_LEVEL_MULTI)  

#ifdef WITH_MKL
            !> in case of smoothing by different MKL solvers, mark levels for the use of MKL,
            !> either by locally acting PARDISO solvers or globally acting CLUSTER_SPARSE_SOLVER
            IF (TYPE_SMOOTH == NSCARC_RELAX_MKL) THEN

               LU_SCOPE_SELECT: SELECT CASE (TYPE_SMOOTH_SCOPE) 
                  CASE (NSCARC_SCOPE_LOCAL)
WRITE(MSG%LU_DEBUG,*) 'SELECT LOCAL SMOOTHING'
                     DO NL = NLEVEL_MIN, NLEVEL_MAX-1
                        TYPE_MKL_LEVEL(NL) = NSCARC_MKL_LOCAL
                     ENDDO
                  CASE (NSCARC_SCOPE_GLOBAL)
WRITE(MSG%LU_DEBUG,*) 'SELECT GLOBAL SMOOTHING'
                     DO NL = NLEVEL_MIN, NLEVEL_MAX-1
                        TYPE_MKL_LEVEL(NL) = NSCARC_MKL_GLOBAL
                     ENDDO
               END SELECT LU_SCOPE_SELECT

            ENDIF

            IF (TYPE_MKL == NSCARC_MKL_COARSE) TYPE_MKL_LEVEL(NLEVEL_MAX) = NSCARC_MKL_GLOBAL
#endif

         !> 
         !> Use of algebraic multigrid method - currently disabled, still in progres
         !> first, only finest level is set, further levels are defined during coarsening process
         !> 
         CASE (NSCARC_MULTIGRID_ALGEBRAIC)
            CALL SCARC_SHUTDOWN(NSCARC_ERROR_AMG_MISSING, SCARC_NONE, NSCARC_NONE)
            !CALL SCARC_GET_NUMBER_OF_LEVELS(NSCARC_LEVEL_AMG)

      END SELECT SELECT_MULTIGRID_TYPE

   !>
   !> ----------------------- Global LU-decomposition -----------------------------------------
   !>
   CASE (NSCARC_METHOD_LU)

      !> Only use single (finest) grid level
      CALL SCARC_GET_NUMBER_OF_LEVELS(NSCARC_LEVEL_SINGLE)

      !> and mark this level for the use of MKL methods
      SELECT_MKL_SCOPE: SELECT CASE (TYPE_MKL)
         CASE (NSCARC_MKL_LOCAL)
            TYPE_MKL_LEVEL(NLEVEL_MIN) = NSCARC_MKL_LOCAL
         CASE (NSCARC_MKL_GLOBAL)
            TYPE_MKL_LEVEL(NLEVEL_MIN) = NSCARC_MKL_GLOBAL
      END SELECT SELECT_MKL_SCOPE

END SELECT SELECT_METHOD

!>
!> Define remaining logical short names based on computed number of levels
!>
#ifdef WITH_MKL
DO NL = NLEVEL_MIN, NLEVEL_MAX
   IS_MKL_LEVEL(NL) = (TYPE_MKL == NSCARC_MKL_GLOBAL .AND. NL == NLEVEL_MIN) .OR. &
                    (TYPE_MKL == NSCARC_MKL_COARSE .AND. NL == NLEVEL_MAX) .OR. &
                    (TYPE_MKL_LEVEL(NL) == NSCARC_MKL_GLOBAL) 
ENDDO
#endif

END SUBROUTINE SCARC_SETUP_LEVELS


!> ------------------------------------------------------------------------------------------------
!> Setup single level in case of default Krylov method
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_GET_NUMBER_OF_LEVELS(NTYPE)
INTEGER, INTENT(IN) :: NTYPE
INTEGER :: KLEVEL(3), KLEVEL_MIN, NM, NLEVEL
TYPE (MESH_TYPE), POINTER :: M=>NULL()

SELECT_LEVEL_TYPE: SELECT CASE (NTYPE)

   !> only use the finest grid level
   CASE(NSCARC_LEVEL_SINGLE)

      NLEVEL     = 1
      NLEVEL_MIN = 1
      NLEVEL_MAX = 1

   !> determine maximum number of possible levels based on number of grid cells
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
      IF (IS_GMG.OR.IS_CG_GMG) THEN
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

   !> use user specified number of grid levels - currently disabled, still in progress
   CASE(NSCARC_LEVEL_AMG)
       CALL SCARC_SHUTDOWN(NSCARC_ERROR_AMG_MISSING, SCARC_NONE, NSCARC_NONE)
   !   NLEVEL_MIN = 1
   !   IF (SCARC_MULTIGRID_LEVEL /= -1) THEN
   !      NLEVEL_MAX  = SCARC_MULTIGRID_LEVEL
   !   ELSE
   !      NLEVEL_MAX  = NSCARC_LEVEL_MAX
   !   ENDIF
   !   NLEVEL = NLEVEL_MAX

   END SELECT SELECT_LEVEL_TYPE

END SUBROUTINE SCARC_GET_NUMBER_OF_LEVELS


!> ------------------------------------------------------------------------------------------------
!> Determine maximum number of possible levels on direction IOR0 of mesh NM
!> In case of GMG- or 2-Level-method, NC must be divisable by 2 at least one time
!> ------------------------------------------------------------------------------------------------
INTEGER FUNCTION SCARC_GET_MAX_LEVEL(NC, IOR0)
INTEGER, INTENT(IN) :: NC, IOR0
INTEGER :: NC0, NL

IF (HAS_MULTIPLE_LEVELS .AND.  MOD(NC,2)/=0) THEN
   SELECT CASE (ABS(IOR0))
      CASE (1)
         CALL SCARC_SHUTDOWN(NSCARC_ERROR_GRID_NUMBERX, SCARC_NONE, NC)
      CASE (2)
         CALL SCARC_SHUTDOWN(NSCARC_ERROR_GRID_NUMBERY, SCARC_NONE, NC)
      CASE (3)
         CALL SCARC_SHUTDOWN(NSCARC_ERROR_GRID_NUMBERZ, SCARC_NONE, NC)
   END SELECT

ENDIF

!> Divide by 2 as often as possible or until user defined max-level is reached
NC0=NC
DO NL=1,NSCARC_LEVEL_MAX
   NC0=NC0/2
   IF (MOD(NC0,2)/=0) EXIT                !> if no longer divisable by two, leave loop ...
   IF (NL==SCARC_MULTIGRID_LEVEL) EXIT    !> if max possible number of levels reached, leave loop ...
   IF (NC0==1) EXIT                       !> if corresponding power of two has been found, leave loop ...
ENDDO

SCARC_GET_MAX_LEVEL=NL

RETURN
END FUNCTION SCARC_GET_MAX_LEVEL


!> ------------------------------------------------------------------------------------------------
!> Allocate basic ScaRC-structures for all needed levels
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_TYPES
INTEGER :: NM
TYPE (SCARC_TYPE), POINTER :: S=>NULL()

!> Basic information for all requested grid levels
ALLOCATE (SCARC(NMESHES), STAT=IERROR)
CALL CHKMEMERR ('SCARC_SETUP', 'SCARC', IERROR)

!> Basic solver stack
ALLOCATE (STACK(NSCARC_STACK_MAX), STAT=IERROR)
CALL CHKMEMERR ('SCARC_SETUP', 'STACK', IERROR)

MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   S => SCARC(NM)

   !> Needed information about other meshes
   ALLOCATE (S%OSCARC(NMESHES), STAT=IERROR)
   CALL CHKMEMERR ('SCARC_SETUP_TYPES', 'OSCARC', IERROR)

   !> Information for single grid levels
   ALLOCATE (S%LEVEL(NLEVEL_MIN:NLEVEL_MAX), STAT=IERROR)
   CALL CHKMEMERR ('SCARC_SETUP_TYPES', 'LEVEL', IERROR)

ENDDO MESHES_LOOP
END SUBROUTINE SCARC_SETUP_TYPES


!> -----------------------------------------------------------------------------
!> Setup discretization information 
!> -----------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_DISCRETIZATION
INTEGER :: NL, NM, IX, IY, IZ, IO, I, J, K
INTEGER :: NX, NY, NZ
REAL(EB), DIMENSION(:), POINTER :: XCOR, YCOR, ZCOR
REAL(EB), DIMENSION(:), POINTER :: XMID, YMID, ZMID
TYPE (MESH_TYPE), POINTER :: M=>NULL()
TYPE (SCARC_TYPE), POINTER :: S=>NULL()
TYPE (SCARC_LEVEL_TYPE), POINTER :: L=>NULL()
TYPE (SCARC_DISCRET_TYPE), POINTER :: D=>NULL()

!>
!> -------------------- On all grid levels --------------------------------------
!> Specify general mesh related geometry information
!>
MESHES_LOOP1: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   M => MESHES(NM)
   S => SCARC(NM)

   !> store bounds of mesh in SCARC-structure
   S%XS = M%XS
   S%XF = M%XF
   S%YS = M%YS
   S%YF = M%YF
   S%ZS = M%ZS
   S%ZF = M%ZF

   S%IBAR = M%IBAR
   S%JBAR = M%JBAR
   S%KBAR = M%KBAR

   NX = M%IBAR
   NY = M%JBAR
   NZ = M%KBAR

   LEVEL_LOOP1: DO NL = NLEVEL_MIN, NLEVEL_MAX

      L => SCARC(NM)%LEVEL(NL)

      L%NX = NX
      L%NY = NY
      L%NZ = NZ

      L%NC  = L%NX * L%NY * L%NZ
      L%NCE = L%NC

      L%N_WALL_CELLS_EXT = M%N_EXTERNAL_WALL_CELLS
      L%N_WALL_CELLS_INT = M%N_INTERNAL_WALL_CELLS
      L%N_WALL_CELLS     = L%N_WALL_CELLS_EXT + L%N_WALL_CELLS_INT

      !> get coordination information
      L%DX = (S%XF-S%XS)/REAL(L%NX,EB)
      L%DY = (S%YF-S%YS)/REAL(L%NY,EB)
      L%DZ = (S%ZF-S%ZS)/REAL(L%NZ,EB)

      L%DXI = 1.0_EB/L%DX
      L%DYI = 1.0_EB/L%DY
      L%DZI = 1.0_EB/L%DZ

      L%DXI2 = L%DXI**2
      L%DYI2 = L%DYI**2
      L%DZI2 = L%DZI**2

      !> needed in case of GMG with multiple grid levels
      NX=NX/2
      IF (.NOT.TWO_D) NY=NY/2
      NZ=NZ/2

      IF (NL == NLEVEL_MIN) THEN

         !> On finest level store information about obstructions
         L%N_OBST = M%N_OBST
         ALLOCATE(L%OBST(L%N_OBST), STAT=IERROR)
         CALL ChkMemErr('SCARC_SETUP_DISCRETIZATION','OBST',IERROR)

         DO IO = 1, L%N_OBST
            L%OBST(IO)%I1  = M%OBSTRUCTION(IO)%I1
            L%OBST(IO)%I2  = M%OBSTRUCTION(IO)%I2
            L%OBST(IO)%J1  = M%OBSTRUCTION(IO)%J1
            L%OBST(IO)%J2  = M%OBSTRUCTION(IO)%J2
            L%OBST(IO)%K1  = M%OBSTRUCTION(IO)%K1
            L%OBST(IO)%K2  = M%OBSTRUCTION(IO)%K2
         ENDDO

         !> Point to already existing arrays from main FDS program
         XCOR => M%X
         YCOR => M%Y
         ZCOR => M%Z

         XMID => M%XC
         YMID => M%YC
         ZMID => M%ZC

      ELSE

         !> Allocate and compute coordinate information for coarser levels
         CALL SCARC_ALLOCATE_REAL1(L%XCOR, 0, L%NX, NSCARC_INIT_NONE, 'XCOR')
         CALL SCARC_ALLOCATE_REAL1(L%YCOR, 0, L%NY, NSCARC_INIT_NONE, 'YCOR')
         CALL SCARC_ALLOCATE_REAL1(L%ZCOR, 0, L%NZ, NSCARC_INIT_NONE, 'ZCOR')

         DO IX = 0, L%NX
            L%XCOR(IX) = S%XS + IX*L%DX
         ENDDO
         DO IY = 0, L%NY
            L%YCOR(IY) = S%YS + IY*L%DY
         ENDDO
         DO IZ = 0, L%NZ
            L%ZCOR(IZ) = S%ZS + IZ*L%DZ
         ENDDO

         !> Allocate and compute midpoint information for coarser levels
         CALL SCARC_ALLOCATE_REAL1(L%XMID, 0, L%NX+1, NSCARC_INIT_NONE, 'XMID')
         CALL SCARC_ALLOCATE_REAL1(L%YMID, 0, L%NY+1, NSCARC_INIT_NONE, 'YMID')
         CALL SCARC_ALLOCATE_REAL1(L%ZMID, 0, L%NZ+1, NSCARC_INIT_NONE, 'ZMID')

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

         XCOR => L%XCOR
         YCOR => L%YCOR
         ZCOR => L%ZCOR

         XMID => L%XMID
         YMID => L%YMID
         ZMID => L%ZMID

      ENDIF

      !> Allocate vectors for step sizes in different directions
      CALL SCARC_ALLOCATE_REAL1(L%DXL, 0, L%NX, NSCARC_INIT_ZERO, 'DXL')
      CALL SCARC_ALLOCATE_REAL1(L%DYL, 0, L%NY, NSCARC_INIT_ZERO, 'DYL')
      CALL SCARC_ALLOCATE_REAL1(L%DZL, 0, L%NZ, NSCARC_INIT_ZERO, 'DZL')

      !> set step sizes between cell midpoints, use interior step sizes for ghost cells as initial values
      !> correct sizes for ghost cells are exchanged later
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

!>
!> ---------------------- On finest grid level -------------------------------------------------
!> Allocate several arrays for the administration of discretization related data
!>
MESHES_LOOP2: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   M => MESHES(NM)
   S => SCARC(NM)
   L => SCARC(NM)%LEVEL(NLEVEL_MIN)

   !> Set pointers to already existing CELL_INDEX WALL_INDEX arrays from main program (on finest level)
   CALL SCARC_SETUP_CELL_INDEX(NM, NLEVEL_MIN)
   CALL SCARC_SETUP_WALL_INDEX(NM, NLEVEL_MIN)

   !> Allocate and initialize CELL_STATE array which indicates the state of a cell (gasphase/solid)
   CALL SCARC_ALLOCATE_INT3(L%CELL_STATE, 0, L%NX+1, 0, L%NY+1, 0, L%NZ+1, NSCARC_INIT_NONE, 'CELL_STATE')
   L%CELL_STATE  = NSCARC_UNDEFINED_INT
   L%CELL_STATE(1:L%NX, 1:L%NY, 1:L%NZ) = NSCARC_CELL_GASPHASE

   !> Identify and mark solid obstruction cells in CELL_STATE-part of the discretization
   DO K = 1, L%NZ
      DO J = 1, L%NY
         DO I = 1, L%NX
            IF (M%SOLID(M%CELL_INDEX(I, J, K))) L%CELL_STATE(I, J, K) = NSCARC_CELL_SOLID
         ENDDO
      ENDDO
   ENDDO

   !>
   !> If both discretization types (structured/unstructured) must be administrated (MGM method only):
   !> Allocate all arrays which are related to a specific discretization type
   !>
   IF (HAS_MULTIPLE_DISCRETIZATIONS) THEN

      !> ---------------- First process structured discretization 
      D => L%SD

      !> Number of local cells per mesh
      CALL SCARC_ALLOCATE_INT1(D%NC_LOCAL , 1, NMESHES, NSCARC_INIT_ZERO, 'NC_LOCAL')
      CALL SCARC_ALLOCATE_INT1(D%NC_OFFSET, 1, NMESHES, NSCARC_INIT_ZERO, 'NC_OFFSET')

      !> Allocate wall information array 
      ALLOCATE(D%WALL(L%N_WALL_CELLS), STAT=IERROR)
      CALL ChkMemErr('SCARC_SETUP_DISCRETIZATION','WALL',IERROR)

      !> Allocate and preset cell numbers array 
      CALL SCARC_ALLOCATE_INT3(D%CELL_NUMBER, 0, L%NX+1, 0, L%NY+1, 0, L%NZ+1, NSCARC_INIT_NONE, 'CELL_NUMBER')
      D%CELL_NUMBER = NSCARC_UNDEFINED_INT

      !> Define local cell numbers for Poisson equation
      DO K=1,L%NZ
         DO J=1,L%NY
            DO I=1,L%NX
               D%NC_LOCAL(NM) = D%NC_LOCAL(NM) + 1
               D%CELL_NUMBER(I,J,K) = D%NC_LOCAL(NM)          
            ENDDO
         ENDDO
      ENDDO
      D%NC = D%NC_LOCAL(NM)


      !> ---------------- Then process unstructured discretization
      D => L%UD

      !> Also allocate and preset cell numbers and state arrays for unstructured discretization
      CALL SCARC_ALLOCATE_INT3(D%CELL_NUMBER, 0, L%NX+1, 0, L%NY+1, 0, L%NZ+1, NSCARC_INIT_NONE, 'CELL_NUMBER')
      D%CELL_NUMBER = NSCARC_UNDEFINED_INT

      !> Number of local cells per mesh
      CALL SCARC_ALLOCATE_INT1(D%NC_LOCAL , 1, NMESHES, NSCARC_INIT_ZERO, 'NC_LOCAL')
      CALL SCARC_ALLOCATE_INT1(D%NC_OFFSET, 1, NMESHES, NSCARC_INIT_ZERO, 'NC_OFFSET')

      !> Allocate wall information array 
      ALLOCATE(D%WALL(L%N_WALL_CELLS), STAT=IERROR)
      CALL ChkMemErr('SCARC_SETUP_DISCRETIZATION','WALL',IERROR)

      !> Define local cell numbers for Poisson equation
      DO K=1,L%NZ
         DO J=1,L%NY
            DO I=1,L%NX
               IF (L%CELL_STATE(I,J,K) == NSCARC_CELL_GASPHASE ) THEN
                  D%NC_LOCAL(NM) = D%NC_LOCAL(NM) + 1
                  D%CELL_NUMBER(I,J,K) = D%NC_LOCAL(NM)      
               ENDIF
            ENDDO
         ENDDO
      ENDDO
      D%NC = D%NC_LOCAL(NM)

   !>
   !> If only one specified type of discretization must be admistrated:
   !> Allocate and preset cell numbers and state arrays for requested type of discretization
   !>
   ELSE

      !> ---------------- Only process specified type of discretization
      SELECT CASE(TYPE_DISCRET)
         CASE (NSCARC_DISCRET_STRUCTURED)
            D => L%SD
         CASE (NSCARC_DISCRET_UNSTRUCTURED)
            D => L%UD
      END SELECT

      !> Also allocate and preset cell numbers and state arrays for unstructured discretization
      CALL SCARC_ALLOCATE_INT3(D%CELL_NUMBER, 0, L%NX+1, 0, L%NY+1, 0, L%NZ+1, NSCARC_INIT_NONE, 'CELL_NUMBER')
      D%CELL_NUMBER = NSCARC_UNDEFINED_INT

      !> Number of local cells per mesh
      CALL SCARC_ALLOCATE_INT1(D%NC_LOCAL , 1, NMESHES, NSCARC_INIT_ZERO, 'NC_LOCAL')
      CALL SCARC_ALLOCATE_INT1(D%NC_OFFSET, 1, NMESHES, NSCARC_INIT_ZERO, 'NC_OFFSET')

      !> Allocate wall information array 
      ALLOCATE(D%WALL(L%N_WALL_CELLS), STAT=IERROR)
      CALL ChkMemErr('SCARC_SETUP_DISCRETIZATION','WALL',IERROR)

      !> Define local cell numbers for Poisson equation
      DO K=1,L%NZ
         DO J=1,L%NY
            DO I=1,L%NX
               IF (L%CELL_STATE(I,J,K) == NSCARC_CELL_GASPHASE ) THEN
                  D%NC_LOCAL(NM) = D%NC_LOCAL(NM) + 1
                  D%CELL_NUMBER(I,J,K) = D%NC_LOCAL(NM)      
               ENDIF
            ENDDO
         ENDDO
      ENDDO
      D%NC = D%NC_LOCAL(NM)

   ENDIF

ENDDO MESHES_LOOP2

!> -------------------------------------------------------------------------------------------
!> Debug DISCRET structures - only if directive SCARC_DEBUG is set
!> -------------------------------------------------------------------------------------------
#ifdef WITH_SCARC_DEBUG
CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_DISCRET, NLEVEL_MIN, 'DISCRET')
#endif
END SUBROUTINE SCARC_SETUP_DISCRETIZATION


!> -----------------------------------------------------------------------------
!> Setup discretization information on coarser levels
!> -----------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_DISCRETIZATION_LEVEL(NL)

INTEGER, INTENT(IN) :: NL
INTEGER :: NM, IXF, IYF, IZF, IX, IY, IZ, NSTEP
TYPE (SCARC_LEVEL_TYPE), POINTER :: LF=>NULL(), LC=>NULL()
TYPE (SCARC_DISCRET_TYPE), POINTER :: DC=>NULL()

MESHES_LOOP1: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   LF => SCARC(NM)%LEVEL(NLEVEL_MIN)
   LC => SCARC(NM)%LEVEL(NL)

   CALL SCARC_ALLOCATE_INT3(LC%CELL_STATE, 0, LC%NX+1, 0, LC%NY+1, 0, LC%NZ+1, NSCARC_INIT_NONE, 'CELL_STATE')
   LC%CELL_STATE  = NSCARC_UNDEFINED_INT
   LC%CELL_STATE(1:LC%NX, 1:LC%NY, 1:LC%NZ) = NSCARC_CELL_GASPHASE

   NSTEP = 2**(NL - NLEVEL_MIN)

   SELECT CASE(TYPE_DISCRET)

      !> Get cell numberings for coarser grid in case of structured discretization
      CASE (NSCARC_DISCRET_STRUCTURED)

         DC => LC%SD

         CALL SCARC_ALLOCATE_INT3(DC%CELL_NUMBER, 0, LC%NX+1, 0, LC%NY+1, 0, LC%NZ+1, NSCARC_INIT_NONE, 'CELL_NUMBER')
         DC%CELL_NUMBER = NSCARC_UNDEFINED_INT

         DO IZ = 1, LC%NZ
            IZF = (IZ-1)*NSTEP + 1
            DO IY = 1, LC%NY
               IYF = (IY-1)*NSTEP + 1
               DO IX = 1, LC%NX
                  IXF = (IX-1)*NSTEP + 1
                  LC%CELL_STATE(IX,IY,IZ) = LF%CELL_STATE(IXF, IYF, IZF)
                  DC%NC_LOCAL(NM)  = DC%NC_LOCAL(NM) + 1
                  DC%CELL_NUMBER(IX,IY,IZ) = DC%NC_LOCAL(NM)
               ENDDO
            ENDDO
         ENDDO

         DC%NC = DC%NC_LOCAL(NM)

      !> Get cell numberings for coarser grid in case of unstructured discretization
      CASE (NSCARC_DISCRET_UNSTRUCTURED)

         DC => LC%UD

         CALL SCARC_ALLOCATE_INT3(DC%CELL_NUMBER, 0, LC%NX+1, 0, LC%NY+1, 0, LC%NZ+1, NSCARC_INIT_NONE, 'CELL_NUMBER')
         DC%CELL_NUMBER = NSCARC_UNDEFINED_INT

         DO IZ = 1, LC%NZ
            IZF = (IZ-1)*NSTEP + 1
            DO IY = 1, LC%NY
               IYF = (IY-1)*NSTEP + 1
               DO IX = 1, LC%NX
                  IXF = (IX-1)*NSTEP + 1
                  LC%CELL_STATE(IX,IY,IZ) = LF%CELL_STATE(IXF, IYF, IZF)
                  IF (LF%CELL_STATE(IXF, IYF, IZF) == NSCARC_CELL_GASPHASE) THEN
                     DC%NC_LOCAL(NM) = DC%NC_LOCAL(NM) + 1
                     DC%CELL_NUMBER(IX,IY,IZ) = DC%NC_LOCAL(NM)
                  ENDIF
               ENDDO
            ENDDO
         ENDDO

         DC%NC = DC%NC_LOCAL(NM)

   END SELECT

ENDDO MESHES_LOOP1

END SUBROUTINE SCARC_SETUP_DISCRETIZATION_LEVEL

!> ------------------------------------------------------------------------------------------------
!> Setup communication structure for data exchange along mesh interfaces
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_NEIGHBORS
INTEGER :: NM, NOM, NL
TYPE (MESH_TYPE), POINTER :: M=>NULL(), OM=>NULL()
TYPE (SCARC_NEIGHBOR_TYPE), POINTER :: OS=>NULL()
TYPE (SCARC_LEVEL_TYPE), POINTER :: OLF=>NULL(), OLC=>NULL()
TYPE (SCARC_DISCRET_TYPE), POINTER :: ODF=>NULL(), ODC=>NULL()

!> Initialize communication counter for ScaRC, use same TAG for all communications
TAG   = 99
N_REQ =  0
N_EXCHANGE =  0

!> Allocate basic WALL and FACE types on mesh NM for all requested grid levels
MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   !> Get communication lengths for other meshes from main program
   M => MESHES(NM)
   OTHER_MESHES_LOOP: DO NOM = 1, NMESHES

      !IF (NOM == NM) CYCLE OTHER_MESHES_LOOP
      OS => SCARC(NM)%OSCARC(NOM)

      OS%NIC_R    = M%OMESH(NOM)%NIC_R
      OS%NIC_S    = M%OMESH(NOM)%NIC_S

      OS%NICMAX_R = M%OMESH(NOM)%NIC_R
      OS%NICMAX_S = M%OMESH(NOM)%NIC_S

      IF (OS%NICMAX_R==0 .AND. OS%NICMAX_S==0) CYCLE OTHER_MESHES_LOOP
      N_EXCHANGE  = N_EXCHANGE+1

   ENDDO OTHER_MESHES_LOOP

ENDDO MESHES_LOOP

!> Initialize level structures on neighboring meshes
LEVEL_MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   LEVEL_OTHER_MESHES_LOOP: DO NOM = 1, NMESHES

      OS  => SCARC(NM)%OSCARC(NOM)
      IF (OS%NICMAX_S==0 .AND. OS%NICMAX_R==0)  CYCLE LEVEL_OTHER_MESHES_LOOP

      ALLOCATE (OS%LEVEL(NLEVEL_MIN:NLEVEL_MAX), STAT=IERROR)
      CALL CHKMEMERR ('SCARC_SETUP_NEIGHBORS', 'OS%LEVEL', IERROR)

      OM  => MESHES(NOM)
      OLF => SCARC(NM)%OSCARC(NOM)%LEVEL(NLEVEL_MIN)

      OLF%NX = OM%IBAR                                              !> number of cells in x-direction on other mesh
      OLF%NY = OM%JBAR                                              !> number of cells in y-direction on other mesh
      OLF%NZ = OM%KBAR                                              !> number of cells in z-direction on other mesh

      OLF%N_WALL_CELLS_EXT = OM%N_EXTERNAL_WALL_CELLS               !> number of external wall cells on other mesh
      OLF%N_WALL_CELLS_INT = OM%N_INTERNAL_WALL_CELLS               !> number of external wall cells on other mesh

      OLF%NC  = OLF%NX*OLF%NY*OLF%NZ                                !> number of cells on other mesh
      OLF%NW  = OLF%N_WALL_CELLS_EXT + OLF%N_WALL_CELLS_INT         !> number of walls cell on other mesh
      OLF%NCG = 0                                                   !> number of ghost cells on other mesh

      SELECT CASE (TYPE_MATRIX)
         CASE (NSCARC_MATRIX_COMPACT)
            ODF%AC%NA = 0                                        !> number of matrix values
            ODF%AC%NR = 0                                        !> number of row pointers
         CASE (NSCARC_MATRIX_BANDED)
            ODF%AB%NA   = 0                                      !> number of matrix values
            ODF%AB%NDIAG  = 0                                      !> number of bands to store
      END SELECT

      IF (OS%NICMAX_S == 0 .AND. OS%NICMAX_R == 0) CYCLE LEVEL_OTHER_MESHES_LOOP

      !> In case of GMG with a predefined grid hierarchy allocate corresponding level-structures
      IF (NLEVEL_MAX > NLEVEL_MIN) THEN

         DO NL=NLEVEL_MIN+1,NLEVEL_MAX

            OLC => SCARC(NM)%OSCARC(NOM)%LEVEL(NL)         !> Note: OLF points to finer, OLC to coarser level
            
            OLC%NX = OLF%NX/2
            IF (TWO_D) THEN
               OLC%NY = 1
            ELSE
               OLC%NY = OLF%NY/2
            ENDIF
            OLC%NZ = OLF%NZ/2

            OLC%N_WALL_CELLS_EXT = 4 * (OLC%NX * OLC%NZ + OLC%NX * OLC%NY + OLC%NY * OLC%NZ)   !> TOFIX : why 4?

            OLC%NC  = OLC%NX * OLC%NY * OLC%NZ                     !> see above
            OLC%NW  = OLC%N_WALL_CELLS_EXT
            OLC%NCG = 0

            SELECT CASE (TYPE_MATRIX)
               CASE (NSCARC_MATRIX_COMPACT)
                  ODC%AC%NA = 0                                            !> number of matrix values
                  ODC%AC%NR = 0                                            !> number of row pointers
               CASE (NSCARC_MATRIX_BANDED)
                  ODC%AB%NA = 0                                            !> number of matrix values
                  ODC%AB%NDIAG = 0                                         !> number of bands to store
            END SELECT

         ENDDO
      ENDIF

   ENDDO LEVEL_OTHER_MESHES_LOOP
ENDDO LEVEL_MESHES_LOOP

END SUBROUTINE SCARC_SETUP_NEIGHBORS

!> ----------------------------------------------------------------------------------------------------
!> Setup FACE related structures 
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_FACES
INTEGER :: NL, NM, NOM
INTEGER :: IFACE, IOR0, JOR0, INBR, IWG, IWL
LOGICAL :: BKNOWN(-3:3)
TYPE (MESH_TYPE), POINTER :: M
TYPE (SCARC_TYPE), POINTER :: S=>NULL()
TYPE (SCARC_LEVEL_TYPE), POINTER :: L=>NULL(),  LC=>NULL()

MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   M => MESHES(NM)
   S => SCARC(NM)
   L => SCARC(NM)%LEVEL(NLEVEL_MIN)

   !> Allocate FACE arrays on different grid levels
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

   NUM_FACE_NEIGHBORS = 0
   NUM_MESH_NEIGHBORS = 0

   !> Determine basic data for single faces (orientation, dimensions, numbers)
   FACES_OF_MESH_INDEX1: DO IOR0 = -3, 3

      IF (IOR0 == 0) CYCLE FACES_OF_MESH_INDEX1

      !> information about face orientation and local dimensions
      SELECT CASE (ABS(IOR0))
         CASE (1)
            L%FACE(IOR0)%NFC  =  L%NX                   !> number of cells between opposite mesh faces
            L%FACE(IOR0)%NFX  =  1                      !> number of cells in x-direction
            L%FACE(IOR0)%NFY  =  L%NY                   !> number of cells in y-direction
            L%FACE(IOR0)%NFZ  =  L%NZ                   !> number of cells in z-direction
            L%FACE(IOR0)%NFW  =  L%NY*L%NZ              !> number of wall cells at that face
            L%FACE(IOR0)%DH   => L%DXL                  !> step size vector between opposite mesh faces
         CASE (2)
            L%FACE(IOR0)%NFC  =  L%NY                   !> dito
            L%FACE(IOR0)%NFX  =  L%NX
            L%FACE(IOR0)%NFY  =  1
            L%FACE(IOR0)%NFZ  =  L%NZ
            L%FACE(IOR0)%NFW  =  L%NX*L%NZ
            L%FACE(IOR0)%DH   => L%DYL
         CASE (3)
            L%FACE(IOR0)%NFC  =  L%NZ                   !> dito
            L%FACE(IOR0)%NFX  =  L%NX
            L%FACE(IOR0)%NFY  =  L%NY
            L%FACE(IOR0)%NFZ  =  1
            L%FACE(IOR0)%NFW  =  L%NX*L%NY
            L%FACE(IOR0)%DH   => L%DZL
      END SELECT

   ENDDO FACES_OF_MESH_INDEX1

   !> Store first IWG-number for each face
   IWG = 1
   FACE_ORDER_LOOP: DO IFACE = 1, 6
      IOR0 = FACE_ORDER_XYZ(IFACE)
      L%FACE(IOR0)%IWG_PTR = IWG
      IWG = IWG + L%FACE(IOR0)%NFW
   ENDDO FACE_ORDER_LOOP

   !> loop over external IW's:
   !> store basic data and determine number of adajacent neighbors to each face
   IWL = 0
   EXTERNAL_WALL_CELLS_LOOP: DO IWG = 1, L%N_WALL_CELLS_EXT

      NOM  = M%EXTERNAL_WALL(IWG)%NOM
      IOR0 = M%WALL(IWG)%ONE_D%IOR

      IF (NOM /= 0) THEN
         BKNOWN = .FALSE.
         DO JOR0 = -3, 3                                                !> neighbor already known?
            IF (JOR0 == 0) CYCLE
            DO INBR = 1, NUM_FACE_NEIGHBORS(JOR0)
               IF (FACE_NEIGHBORS(JOR0, INBR) == NOM) THEN
                  BKNOWN(JOR0) = .TRUE.
                  EXIT
               ENDIF
            ENDDO
         ENDDO
         IF (.NOT.BKNOWN(IOR0)) THEN
            NUM_FACE_NEIGHBORS(IOR0) = NUM_FACE_NEIGHBORS(IOR0) + 1     !> increase neighbor counter for face
            FACE_NEIGHBORS(IOR0, NUM_FACE_NEIGHBORS(IOR0)) = NOM        !> store number of neighbor for face
         ENDIF
         IF (.NOT.ANY(BKNOWN)) THEN
            NUM_MESH_NEIGHBORS = NUM_MESH_NEIGHBORS + 1                 !> increase neighbor counter for mesh
            MESH_NEIGHBORS(NUM_FACE_NEIGHBORS(IOR0)) = NOM              !> store number of neighbor for mesh
         ENDIF
      ENDIF

   ENDDO EXTERNAL_WALL_CELLS_LOOP

   !> Allocate array which stores numbers of all neighboring meshes
   IF (NUM_MESH_NEIGHBORS /= 0) &
      CALL SCARC_ALLOCATE_INT1(S%NEIGHBORS, 1, NUM_MESH_NEIGHBORS, NSCARC_INIT_UNDEF, 'NEIGHBORS')
   S%N_NEIGHBORS = NUM_MESH_NEIGHBORS

   !> Store information about adjacent neighbors on different faces
   !> Allocate corresponding index arrays in OSCARC-structures
   !> First allocate administrative mapping arrays for own mesh
   CALL SCARC_SETUP_POINTERS(NM, NLEVEL_MIN)

   NEIGHBORS_OF_FACE_LOOP: DO IOR0 = -3, 3

      IF (IOR0 == 0) CYCLE NEIGHBORS_OF_FACE_LOOP

      !> if there are neighbors at face IOR0 store information about them
      IF (NUM_FACE_NEIGHBORS(IOR0) /= 0) THEN

         !> allocate array for storing the numbers of the single neighbors
         L%FACE(IOR0)%N_NEIGHBORS = NUM_FACE_NEIGHBORS(IOR0)       
         CALL SCARC_ALLOCATE_INT1(L%FACE(IOR0)%NEIGHBORS, 1, NUM_FACE_NEIGHBORS(IOR0), NSCARC_INIT_NONE, 'NEIGHBORS')

         !> store every neighbor and allocate corresponding administration arrays
         DO INBR = 1, NUM_FACE_NEIGHBORS(IOR0)
            NOM = FACE_NEIGHBORS(IOR0, INBR)
            L%FACE(IOR0)%NEIGHBORS(INBR) = NOM               !> store NOM as a neighbor of that face and if
            CALL SCARC_UPDATE_NEIGHBORS(NM, NOM)             !> not already done also as mesh neighbor itself
            CALL SCARC_SETUP_OPOINTERS(NM, NLEVEL_MIN, NOM)  !> allocate pointers for neighboring meshes
         ENDDO

      ENDIF
   ENDDO NEIGHBORS_OF_FACE_LOOP

ENDDO MESHES_LOOP

!> -------------------------------------------------------------------------------------------
!> Debug FACE structures - only if directive SCARC_DEBUG is set
!> -------------------------------------------------------------------------------------------
#ifdef WITH_SCARC_DEBUG
CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_FACE, NLEVEL_MIN, 'FACE')
#endif

END SUBROUTINE SCARC_SETUP_FACES


!> ----------------------------------------------------------------------------------------------------
!> Setup FACE and WALL related structures and boundary conditions
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_WALLS
USE GEOMETRY_FUNCTIONS, ONLY: SEARCH_OTHER_MESHES
INTEGER :: NL, NM, NM2, NOM, NCPL
INTEGER :: IREFINE, IFACE, IOR0, INBR, IWG, IWL, IWC
LOGICAL :: BKNOWN(-3:3), IS_BC_DIRICHLET, IS_OPEN_BOUNDARY
TYPE (MESH_TYPE), POINTER :: M
TYPE (SCARC_TYPE), POINTER :: S=>NULL()
TYPE (SCARC_LEVEL_TYPE), POINTER :: L=>NULL(),  LF=>NULL(),  LC=>NULL()
TYPE (SCARC_LEVEL_TYPE), POINTER :: OL=>NULL(), OLF=>NULL(), OLC=>NULL()
TYPE (SCARC_DISCRET_TYPE), POINTER :: D=>NULL(), DC=>NULL()
TYPE (SCARC_POINTER_TYPE), POINTER :: PTR=>NULL()
TYPE (SCARC_WALL_TYPE), POINTER :: DWC
TYPE (WALL_TYPE), POINTER :: MWC
TYPE (EXTERNAL_WALL_TYPE), POINTER :: EWC

MESHES_LOOP1: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   M => MESHES(NM)
   S => SCARC(NM)
   L => SCARC(NM)%LEVEL(NLEVEL_MIN)

   SELECT CASE(TYPE_DISCRET)
      CASE (NSCARC_DISCRET_STRUCTURED)
         D => L%SD
      CASE (NSCARC_DISCRET_UNSTRUCTURED)
         D => L%UD
   END SELECT
   PTR => D%PTR

   !> First loop over external wall cells:
   !> Determine number of adajacent neighbors to each face with corresponding number of IW's
   IWL = 0
   EXTERNAL_WALL_CELLS_LOOP1: DO IWG = 1, L%N_WALL_CELLS_EXT

      MWC  => M%WALL(IWG)
      EWC => M%EXTERNAL_WALL(IWG)

      !> Determine and store neighbors, orientation and number of couplings for a single wall cell
      NOM  =  EWC%NOM
      IOR0 =  MWC%ONE_D%IOR
      NCPL = (EWC%IIO_MAX - EWC%IIO_MIN + 1) * &
             (EWC%JJO_MAX - EWC%JJO_MIN + 1) * &
             (EWC%KKO_MAX - EWC%KKO_MIN + 1)

      D%WALL(IWG)%NOM  = NOM                            !> store number of neighbor in wall cell
      D%WALL(IWG)%IOR  = IOR0                           !> store orientation of that cell

      IWL = IWL + 1                                     !> count local wall cells for that face

      IF (NOM /= 0) THEN
         D%WALL(IWG)%NCPL = NCPL                        !> store number of couplings for that cell
      ELSE
         D%WALL(IWG)%NCPL = 0                           !> no couplings
      ENDIF

      IF (NOM /= 0) THEN

         L%NCE = L%NCE + NCPL                           !> increase number of extended grid cells

         OL => SCARC(NM)%OSCARC(NOM)%LEVEL(NLEVEL_MIN)
         OL%NCPL = NCPL                                 !> initialize own counter for local wall cells
         OL%IOR  = -IOR0                                !> initialize own orientation variable

         IF (ANY(BKNOWN)) THEN
            OL%NWL = OL%NWL + 1                         !> increase own counter for local wall cells
            OL%NCG = OL%NCG + NCPL                      !> increase counter for local ghost cells
         ELSE
            OL%NWL = 1                                  !> initialize own counter for local wall cells
            OL%NCG = NCPL                               !> initialize counter for local ghost cells
            L%NCPL_MAX  = MAX(L%NCPL_MAX, NCPL)         !> get max NCPL ever used on this mesh
         ENDIF
      ENDIF

   ENDDO EXTERNAL_WALL_CELLS_LOOP1

   !> Then process internal wall cells
   INTERNAL_WALL_CELLS_LOOP1: DO IWG = L%N_WALL_CELLS_EXT+1, L%N_WALL_CELLS_EXT+L%N_WALL_CELLS_INT

      MWC => M%WALL(IWG)
      DWC => D%WALL(IWG)

      DWC%IOR  = MWC%ONE_D%IOR
      DWC%NOM  = 0

      DWC%BTYPE = NEUMANN
      DWC%BOUNDARY_TYPE = M%WALL(IWG)%BOUNDARY_TYPE

      DWC%IXG =  MWC%ONE_D%II                        !> ghost cell indices
      DWC%IYG =  MWC%ONE_D%JJ
      DWC%IZG =  MWC%ONE_D%KK

      DWC%IXW =  MWC%ONE_D%IIG                       !> (internal) wall cell indices
      DWC%IYW =  MWC%ONE_D%JJG
      DWC%IZW =  MWC%ONE_D%KKG

   ENDDO INTERNAL_WALL_CELLS_LOOP1


   !> Second loop over external wall cells:
   !> Store detailed coordinate and cell data and get type of boundary condition
   PTR%ICE = L%NC
   PTR%ICO = L%NC

   WALL_CELLS_LOOP2: DO IWG = 1, L%N_WALL_CELLS_EXT

      NOM  = D%WALL(IWG)%NOM
      IOR0 = D%WALL(IWG)%IOR
      NCPL = D%WALL(IWG)%NCPL

      MWC => M%WALL(IWG)
      EWC => M%EXTERNAL_WALL(IWG)

      !>
      !> Preset ScaRC's boundary type indicator BTYPE
      !> INTERNAL  : the global Poisson problem is solved, no need to impose BC's along mesh interfaces
      !> DIRICHLET : in the structured case face-wise BC-SETTING are used ccording to FFT-SETTING
      !>             (this also allows to use FFT as local preconditioner)
      !>             in the unstructured case Dirichlet BCs are only used for open boundary cells
      !> NEUMANN   : is used for the rest
      !>
      IS_BC_DIRICHLET  = MWC%PRESSURE_BC_INDEX == DIRICHLET
      IS_OPEN_BOUNDARY = MWC%BOUNDARY_TYPE     == OPEN_BOUNDARY

      DWC => D%WALL(IWG)

      IF (EWC%NOM /= 0) THEN
         DWC%BTYPE = INTERNAL
      ELSE IF ((TYPE_DISCRET == NSCARC_DISCRET_STRUCTURED   .AND. IS_BC_DIRICHLET) .OR. &
               (TYPE_DISCRET == NSCARC_DISCRET_UNSTRUCTURED .AND. IS_OPEN_BOUNDARY)) THEN
         DWC%BTYPE = DIRICHLET
         D%N_DIRIC = D%N_DIRIC + 1
      ELSE
         DWC%BTYPE = NEUMANN
         D%N_NEUMANN = D%N_NEUMANN + 1
      ENDIF

      DWC%BOUNDARY_TYPE  = MWC%BOUNDARY_TYPE

      DWC%IXG =  MWC%ONE_D%II                                 !> ghost cell indices
      DWC%IYG =  MWC%ONE_D%JJ
      DWC%IZG =  MWC%ONE_D%KK

      DWC%IXW =  MWC%ONE_D%IIG                                !> (internal) wall cell indices
      DWC%IYW =  MWC%ONE_D%JJG
      DWC%IZW =  MWC%ONE_D%KKG

      !> If there exists a neighbor for that wall cell, setup corresponding neighborship information
      IF (NOM /= 0) THEN
         CALL SCARC_SETUP_WALLCELL_NEIGHBOR(EWC%IIO_MIN, EWC%IIO_MAX, &
                                            EWC%JJO_MIN, EWC%JJO_MAX, &
                                            EWC%KKO_MIN, EWC%KKO_MAX, &
                                            IWG, IOR0, NM, NOM, NLEVEL_MIN)
      ENDIF

   ENDDO WALL_CELLS_LOOP2

ENDDO MESHES_LOOP1

CALL SCARC_SETUP_SUBDIVISION(NLEVEL_MIN)

!>
!> Set dimensions on finest level for requested type(s) of discretization 
!>
CALL SCARC_SETUP_DIMENSIONS(NLEVEL_MIN)

!> If necessary (in case of MG-method) also get discretization information and dimensions on coarser levels
DO NL = NLEVEL_MIN+1, NLEVEL_MAX
   CALL SCARC_SETUP_DISCRETIZATION_LEVEL(NL)
   CALL SCARC_SETUP_DIMENSIONS(NL)
ENDDO

!> Count number of Dirichlet BCs on finest level 
LOCAL_INT = 0
DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   DO NM2 = LOWER_MESH_INDEX, UPPER_MESH_INDEX
      L => SCARC(NM2)%LEVEL(NLEVEL_MIN)
      LOCAL_INT(NM2) = D%N_DIRIC
   ENDDO
ENDDO
N_DIRIC_GLOBAL(NLEVEL_MIN) = SCARC_BROADCAST_INT(NSCARC_BROADCAST_SUM)

!> In case that there are no Dirichlet boundary cells, 
!> allocate mapping information for overlapped cells for later definition of condensed matrix system
IF (N_DIRIC_GLOBAL(NLEVEL_MIN) == 0) THEN
   DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
      L => SCARC(NM)%LEVEL(NLEVEL_MIN)
      SELECT CASE(TYPE_DISCRET)
         CASE (NSCARC_DISCRET_STRUCTURED)
            D => L%SD
         CASE (NSCARC_DISCRET_UNSTRUCTURED)
            D => L%UD
      END SELECT
      CALL SCARC_ALLOCATE_INT1 (D%PTR%ICE_TO_ICN, L%NC+1, L%NCE, NSCARC_INIT_ZERO, 'ICE_TO_ICN')
      CALL SCARC_ALLOCATE_REAL1(D%PTR%ICE_TO_VAL, L%NC+1, L%NCE, NSCARC_INIT_ZERO, 'ICE_TO_VAL')
   ENDDO
ENDIF


!>
!> Only in case of Twolevel-CG- or GMG-method (as main solver or preconditioner):
!> Determine WALL, FACE and OSCARC types for coarser levels
!>
MULTI_LEVEL_IF: IF (HAS_MULTIPLE_LEVELS) THEN

   MESHES_LOOP2: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

      S => SCARC(NM)

      IREFINE=1
      LEVEL_GMG_LEVEL_LOOP: DO NL = NLEVEL_MIN+1, NLEVEL_MAX

         LF => SCARC(NM)%LEVEL(NL-1)                 !> LF points to finer level
         LC => SCARC(NM)%LEVEL(NL)                   !> LC points to coarser level

         IREFINE=IREFINE*2
         CALL SCARC_CHECK_DIVISIBILITY(LF%NCE-LF%NC, 'LF%NCE')

         LC%NCE = LC%NC + (LF%NCE-LF%NC)/2

         SELECT CASE(TYPE_DISCRET)
            CASE (NSCARC_DISCRET_STRUCTURED)
               DC => LC%SD
            CASE (NSCARC_DISCRET_UNSTRUCTURED)
               DC => LC%UD
         END SELECT

         DC%PTR%ICE = LC%NC
         DC%PTR%ICO = LC%NC

         LC%N_WALL_CELLS_EXT = SCARC_COUNT_EXTERNAL_WALL_CELLS(NM, NL)
         LC%N_WALL_CELLS_INT = SCARC_COUNT_INTERNAL_WALL_CELLS(NM, NL)

         LC%NW = LC%N_WALL_CELLS_EXT + LC%N_WALL_CELLS_INT
         ALLOCATE(DC%WALL(LC%NW), STAT=IERROR)
         CALL ChkMemErr('SCARC_SETUP_NEIGHBORS','WALL',IERROR)

         CALL SCARC_SETUP_POINTERS(NM, NL)

         !> compute FACE and WALL information for all faces of coarser level
         IWC = 1
         IWG = 1
         DO IFACE = 1, 6

            IOR0 = FACE_ORDER_XYZ(IFACE)

            !> compute mesh dimensions of coarser mesh level
            CALL SCARC_SETUP_FACE_DIMENSIONS(IOR0, IWG, NM, NL)

            !> for every neighbor do:
            IF (LF%FACE(IOR0)%N_NEIGHBORS /= 0) THEN
               DO INBR = 1, LF%FACE(IOR0)%N_NEIGHBORS

                  NOM = LF%FACE(IOR0)%NEIGHBORS(INBR)
                  
                  OLF => SCARC(NM)%OSCARC(NOM)%LEVEL(NL-1)
                  OLC => SCARC(NM)%OSCARC(NOM)%LEVEL(NL)

                  OLC%NCPL = OLF%NCPL

                  CALL SCARC_CHECK_DIVISIBILITY(OLF%NWL, 'OLF%NWL')
                  CALL SCARC_CHECK_DIVISIBILITY(OLF%NCG, 'OLF%NCG')

                  IF (.NOT.TWO_D) THEN
                     OLC%NWL = OLF%NWL/4
                     OLC%NCG = OLF%NCG/4
                  ELSE
                     OLC%NWL = OLF%NWL/2
                     OLC%NCG = OLF%NCG/2
                  ENDIF

                  CALL SCARC_SETUP_EXCHANGE_DIMENSIONS(IREFINE, NM, NOM, NL)
                  CALL SCARC_SETUP_OPOINTERS(NM, NL, NOM)

               ENDDO
            ENDIF

            !> setup complete face information for coarser mesh
            CALL SCARC_SETUP_FACE(IOR0, IWC, IREFINE, NM, NL)

         ENDDO
         CALL SCARC_SETUP_CELL_INDEX(NM, NL)
         CALL SCARC_SETUP_INTERNAL_WALL_COORDS(NM, NL)
         CALL SCARC_SETUP_WALL_INDEX(NM, NL)

      ENDDO LEVEL_GMG_LEVEL_LOOP
   ENDDO MESHES_LOOP2
ENDIF MULTI_LEVEL_IF

!> -------------------------------------------------------------------------------------------
!> Debug WALL structures - only if directive SCARC_DEBUG is set
!> -------------------------------------------------------------------------------------------
#ifdef WITH_SCARC_DEBUG
CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_WALL, NLEVEL_MIN, 'WALL')
#endif

END SUBROUTINE SCARC_SETUP_WALLS


!> ------------------------------------------------------------------------------------------------
!> Initialize arrays for data exchange
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_EXCHANGE
INTEGER :: NL, NM, NOM, NLEN_SEND, NLEN_RECV
INTEGER :: NLMIN, NLMAX
INTEGER :: INBR
TYPE (SCARC_TYPE), POINTER :: S=>NULL()
TYPE (SCARC_NEIGHBOR_TYPE), POINTER :: OS=>NULL()
TYPE (SCARC_LEVEL_TYPE), POINTER :: L=>NULL(), OL=>NULL()
TYPE (SCARC_DISCRET_TYPE), POINTER :: OD=>NULL()

!>
!> Allocate request array for data exchanges
!> Exchange basic information about wall sizes (needed for the dimensioning of the exchange buffers)
!>
IF (N_MPI_PROCESSES>1) THEN
   ALLOCATE (REQ(N_EXCHANGE*40), STAT=IERROR)
   CALL CHKMEMERR ('SCARC_SETUP_EXCHANGE', 'REQ', IERROR)
   REQ = MPI_REQUEST_NULL
ENDIF

CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_BASIC, NLEVEL_MIN)

!>
!> Allocate send and receive buffers (real and integer) in correct lengths
!>
DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   S => SCARC(NM)
   DO INBR = 1, S%N_NEIGHBORS

      NOM = S%NEIGHBORS(INBR)

      L  => SCARC(NM)%LEVEL(NLEVEL_MIN)

      OS => SCARC(NM)%OSCARC(NOM)
      OL => SCARC(NM)%OSCARC(NOM)%LEVEL(NLEVEL_MIN)
      SELECT CASE(TYPE_DISCRET)
         CASE (NSCARC_DISCRET_STRUCTURED)
            OD => OL%SD
         CASE (NSCARC_DISCRET_UNSTRUCTURED)
            OD => OL%UD
      END SELECT

      !> allocate real buffers in maximum needed length
      NLEN_SEND = NSCARC_MAX_STENCIL*MAX(OL%NWL, OL%NCG)+10
      NLEN_RECV = NLEN_SEND

      CALL SCARC_ALLOCATE_REAL1(OS%SEND_REAL, 1, NLEN_SEND, NSCARC_INIT_ZERO, 'SEND_REAL')
      CALL SCARC_ALLOCATE_REAL1(OS%RECV_REAL, 1, NLEN_RECV, NSCARC_INIT_ZERO, 'RECV_REAL')

      !> allocate integer buffers in maximum needed length
      OL%NCPLS = OL%NCPL
      IF (OL%NCG == OL%NWL) THEN
         NLEN_SEND = 15*OL%NWL
         NLEN_RECV = 15*OL%NCG
         OL%NCPLR = OL%NCPL
      ELSE IF (OL%NCG == 2*OL%NWL) THEN
         NLEN_SEND = 13*OL%NWL + 2*OL%NCG
         NLEN_RECV = 15*OL%NCG
         OL%NCPLR =  1
      ELSE IF (OL%NWL == 2*OL%NCG) THEN
         NLEN_SEND = 15*OL%NWL
         NLEN_RECV = 13*OL%NCG + 2*OL%NWL
         OL%NCPLR =  2
      ENDIF

      CALL SCARC_ALLOCATE_INT1(OS%SEND_INT, 1, NLEN_SEND, NSCARC_INIT_ZERO, 'SEND_INT')
      CALL SCARC_ALLOCATE_INT1(OS%RECV_INT, 1, NLEN_RECV, NSCARC_INIT_ZERO, 'RECV_INT')

      !> neighboring wall structures for common wall cells
      ALLOCATE (OD%WALL(OL%NCG), STAT=IERROR)
      CALL CHKMEMERR ('SCARC_SETUP_NEIGHBORS', 'OD%WALL', IERROR)

      !> In case of GMG with a predefined grid hierarchy allocate corresponding level-structures
      IF (NLEVEL_MAX > NLEVEL_MIN) THEN
         DO NL=NLEVEL_MIN+1,NLEVEL_MAX
            OL => SCARC(NM)%OSCARC(NOM)%LEVEL(NL)
            ALLOCATE (OD%WALL(OL%NCG), STAT=IERROR)
            CALL CHKMEMERR ('SCARC_SETUP_NEIGHBORS', 'OD%WALL', IERROR)
         ENDDO
      ENDIF

   ENDDO
ENDDO

!>
!> Initialize communication structures on finest level (if there is more than 1 mesh)
!>
IF (N_MPI_PROCESSES > 1) THEN
   NLMIN = NLEVEL_MIN
   NLMAX = NLEVEL_MAX
   DO NL = NLMIN, NLMAX
      CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_DISCRET, NL)
      CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_CELL_WIDTH, NL)
!      CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_WALL_INFO, NL)
   ENDDO
ENDIF

!>
!> Correct boundary types for cells adjacent to obstructions on ghost cells
!>
DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   IF (TYPE_DISCRET == NSCARC_DISCRET_UNSTRUCTURED) CALL SCARC_IDENTIFY_INTERNAL_NEUMANNS(NM, NLEVEL_MIN)
   DO NL = NLEVEL_MIN+1, NLEVEL_MAX
      CALL SCARC_IDENTIFY_INTERNAL_NEUMANNS(NM, NL)
   ENDDO
ENDDO

END SUBROUTINE SCARC_SETUP_EXCHANGE

!> -----------------------------------------------------------------------------------------
!> --- Store neighbors of mesh
!> -----------------------------------------------------------------------------------------
SUBROUTINE SCARC_UPDATE_NEIGHBORS(NM, NOM)
INTEGER, INTENT(IN) :: NM, NOM
INTEGER :: INBR
DO INBR = 1, SCARC(NM)%N_NEIGHBORS
   IF (SCARC(NM)%NEIGHBORS(INBR) == NSCARC_UNDEFINED_INT) EXIT      !> not found, to be stored
   IF (SCARC(NM)%NEIGHBORS(INBR) == NOM) RETURN                     !> nothing to do, already stored
ENDDO
SCARC(NM)%NEIGHBORS(INBR) = NOM
RETURN
END SUBROUTINE SCARC_UPDATE_NEIGHBORS


!> -----------------------------------------------------------------------------------------
!> --- Setup CELL_INDEX array on coarser grid levels in case of MG-method
!> -----------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_CELL_INDEX(NM, NL)
INTEGER, INTENT(IN) :: NM, NL
INTEGER :: I, J, K, NOBST
TYPE (MESH_TYPE), POINTER :: M=>NULL()
TYPE (SCARC_LEVEL_TYPE), POINTER :: L=>NULL()
TYPE (SCARC_OBST_TYPE), POINTER :: OB=>NULL()

M => MESHES(NM)
L => SCARC(NM)%LEVEL(NL)

!> if finest level, the corresponding CELL_INDEX array is already available by surrounding routines
IF (NL == NLEVEL_MIN) THEN
   L%CELL_INDEX_PTR => M%CELL_INDEX              

!> on coarser levels, it must still be computed
ELSE

   CALL SCARC_ALLOCATE_INT3(L%CELL_INDEX, 0, L%NX+1, 0, L%NY+1, 0, L%NZ+1, NSCARC_INIT_ZERO, 'CELL_INDEX')
   L%CELL_INDEX_PTR => L%CELL_INDEX              

   L%N_CELL_INDEX = 0
   
   !>
   !> Preset it for all grid cells
   !>
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
   
   !>
   !> Consider cells in obstructions
   !>
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


!> -----------------------------------------------------------------------------------------
!> In case of an MG-method: 
!> Setup WALL_INDEX array on coarser grid levels 
!> -----------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_WALL_INDEX(NM, NL)
INTEGER, INTENT(IN) :: NM, NL
INTEGER :: I, J, K, ICG, IW, IOR0
TYPE (MESH_TYPE), POINTER :: M=>NULL()
TYPE (SCARC_LEVEL_TYPE), POINTER :: L=>NULL()
TYPE (SCARC_DISCRET_TYPE), POINTER :: D=>NULL()

M => MESHES(NM)
L => SCARC(NM)%LEVEL(NL)

!> if on finest level, the array WALL_INDEX is already available by surrounding routines
IF (NL == NLEVEL_MIN) THEN

   L%WALL_INDEX_PTR => M%WALL_INDEX

!> if on coarser levels, it must still be computed
ELSE

   CALL SCARC_ALLOCATE_INT2(L%WALL_INDEX, 1, L%N_CELL_INDEX, -3, 3, NSCARC_INIT_ZERO, 'WALL_INDEX')
   L%WALL_INDEX_PTR => L%WALL_INDEX
   
   SELECT CASE(TYPE_DISCRET)
      CASE (NSCARC_DISCRET_STRUCTURED)
         D => L%SD
      CASE (NSCARC_DISCRET_UNSTRUCTURED)
         D => L%UD
   END SELECT

   DO IW = 1, L%N_WALL_CELLS
   
      I = D%WALL(IW)%IXW
      J = D%WALL(IW)%IYW
      K = D%WALL(IW)%IZW
   
      IOR0 = D%WALL(IW)%IOR
      ICG  = L%CELL_INDEX(I,J,K)
   
      L%WALL_INDEX(ICG,-IOR0) = IW
   
   ENDDO

ENDIF

END SUBROUTINE SCARC_SETUP_WALL_INDEX


!> -------------------------------------------------------------------------------------------------
!> Setup all necessary information for a wall cell with neighbor in case of MG-method
!> Number of obstructions on coarse level is the same as on fine level
!> ACHTUNG: FUNKTIONIERT NUR FUER SPEZIALFAELLE, DIE AUCH FUER GMG LAUFEN !>!
!> ACHTUNG: MUSS DRINGEND NOCH ERWEITERT WERDEN
!> -------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_INTERNAL_WALL_COORDS(NM, NL)
INTEGER, INTENT(IN) :: NM, NL
INTEGER :: IC, IO, IWC
INTEGER :: I, J, K
TYPE (SCARC_LEVEL_TYPE), POINTER :: LC=>NULL()
TYPE (SCARC_OBST_TYPE), POINTER :: OB=>NULL()
TYPE (SCARC_DISCRET_TYPE), POINTER :: DC=>NULL()

LC => SCARC(NM)%LEVEL(NL)
SELECT CASE(TYPE_DISCRET)
   CASE (NSCARC_DISCRET_STRUCTURED)
      DC => LC%SD
   CASE (NSCARC_DISCRET_UNSTRUCTURED)
      DC => LC%UD
END SELECT

IWC = LC%N_WALL_CELLS_EXT + 1

DO IO = 1, LC%N_OBST

   OB => SCARC(NM)%LEVEL(NL)%OBST(IO)

   !> Analyze IOR = 1
   I = OB%I1
   DO K = OB%K1+1, OB%K2
      DO J = OB%J1+1, OB%J2
         IC = DC%CELL_NUMBER(I, J, K)
         IF (IC == NSCARC_UNDEFINED_INT) CYCLE
         DC%WALL(IWC)%IXW = I+1
         DC%WALL(IWC)%IYW = J
         DC%WALL(IWC)%IZW = K
         DC%WALL(IWC)%IXG = I
         DC%WALL(IWC)%IYG = J
         DC%WALL(IWC)%IZG = K
         DC%WALL(IWC)%IOR = 1
         DC%WALL(IWC)%BTYPE = NEUMANN
         DC%WALL(IWC)%BOUNDARY_TYPE = SOLID_BOUNDARY
         IWC = IWC + 1
      ENDDO
   ENDDO

   !> Analyze IOR = -1
   I = OB%I2
   DO K = OB%K1+1, OB%K2
      DO J = OB%J1+1, OB%J2
         IC = DC%CELL_NUMBER(I+1, J, K)
         IF (IC == NSCARC_UNDEFINED_INT) CYCLE
         DC%WALL(IWC)%IXW = I
         DC%WALL(IWC)%IYW = J
         DC%WALL(IWC)%IZW = K
         DC%WALL(IWC)%IXG = I+1
         DC%WALL(IWC)%IYG = J
         DC%WALL(IWC)%IZG = K
         DC%WALL(IWC)%IOR =-1
         DC%WALL(IWC)%BTYPE = NEUMANN
         DC%WALL(IWC)%BOUNDARY_TYPE = SOLID_BOUNDARY
         IWC = IWC + 1
      ENDDO
   ENDDO

   !> Analyze IOR = 2
   J = OB%J1
   DO K = OB%K1+1, OB%K2
      DO I = OB%I1+1, OB%I2
         IC = DC%CELL_NUMBER(I, J, K)
         IF (IC == NSCARC_UNDEFINED_INT) CYCLE
         DC%WALL(IWC)%IXW = I
         DC%WALL(IWC)%IYW = J+1
         DC%WALL(IWC)%IZW = K
         DC%WALL(IWC)%IXG = I
         DC%WALL(IWC)%IYG = J
         DC%WALL(IWC)%IZG = K
         DC%WALL(IWC)%IOR = 2
         DC%WALL(IWC)%BTYPE = NEUMANN
         DC%WALL(IWC)%BOUNDARY_TYPE = SOLID_BOUNDARY
         IWC = IWC + 1
      ENDDO
   ENDDO

   !> Analyze IOR = -2
   J = OB%J2
   DO K = OB%K1+1, OB%K2
      DO I = OB%I1+1, OB%I2
         IC = DC%CELL_NUMBER(I, J+1, K)
         IF (IC == NSCARC_UNDEFINED_INT) CYCLE
         DC%WALL(IWC)%IXW = I
         DC%WALL(IWC)%IYW = J
         DC%WALL(IWC)%IZW = K
         DC%WALL(IWC)%IXG = I
         DC%WALL(IWC)%IYG = J+1
         DC%WALL(IWC)%IZG = K
         DC%WALL(IWC)%IOR =-2
         DC%WALL(IWC)%BTYPE = NEUMANN
         DC%WALL(IWC)%BOUNDARY_TYPE = SOLID_BOUNDARY
         IWC = IWC + 1
      ENDDO
   ENDDO

   !> Analyze IOR = 3
   K = OB%K1
   DO J = OB%J1+1, OB%J2
      DO I = OB%I1+1, OB%I2
         IC = DC%CELL_NUMBER(I, J, K)
         IF (IC == NSCARC_UNDEFINED_INT) CYCLE
         DC%WALL(IWC)%IXW = I
         DC%WALL(IWC)%IYW = J
         DC%WALL(IWC)%IZW = K+1
         DC%WALL(IWC)%IXG = I
         DC%WALL(IWC)%IYG = J
         DC%WALL(IWC)%IZG = K
         DC%WALL(IWC)%IOR = 3
         DC%WALL(IWC)%BTYPE = NEUMANN
         DC%WALL(IWC)%BOUNDARY_TYPE = SOLID_BOUNDARY
         IWC = IWC + 1
      ENDDO
   ENDDO

   !> Analyze IOR = -3
   K = OB%K2
   DO J = OB%J1+1, OB%J2
      DO I = OB%I1+1, OB%I2
         IC = DC%CELL_NUMBER(I, J, K+1)
         IF (IC == NSCARC_UNDEFINED_INT) CYCLE
         DC%WALL(IWC)%IXW = I
         DC%WALL(IWC)%IYW = J
         DC%WALL(IWC)%IZW = K
         DC%WALL(IWC)%IXG = I
         DC%WALL(IWC)%IYG = J
         DC%WALL(IWC)%IZG = K+1
         DC%WALL(IWC)%IOR =-3
         DC%WALL(IWC)%BTYPE = NEUMANN
         DC%WALL(IWC)%BOUNDARY_TYPE = SOLID_BOUNDARY
         IWC = IWC + 1
      ENDDO
   ENDDO

ENDDO

END SUBROUTINE SCARC_SETUP_INTERNAL_WALL_COORDS


!> -------------------------------------------------------------------------------------------------
!> Correct BTYPE related to internal obstructions on ghost cells
!> ACHTUNG: NOCHMAL UEBERARBEITEN !
!> -------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_IDENTIFY_INTERNAL_NEUMANNS(NM, NL)
INTEGER, INTENT(IN) :: NM, NL
INTEGER :: IWG
INTEGER :: IX, IY, IZ, IOR0, BTYPE0
TYPE (SCARC_LEVEL_TYPE), POINTER :: L=>NULL()
TYPE (SCARC_DISCRET_TYPE), POINTER :: D=>NULL()
TYPE (SCARC_WALL_TYPE), POINTER :: WC=>NULL()

L => SCARC(NM)%LEVEL(NL)
SELECT CASE(TYPE_DISCRET)
   CASE (NSCARC_DISCRET_STRUCTURED)
      D => L%SD
   CASE (NSCARC_DISCRET_UNSTRUCTURED)
      D => L%UD
END SELECT

DO IWG = 1, L%N_WALL_CELLS_EXT

   WC => D%WALL(IWG)

   IF (WC%NOM == 0) CYCLE                    !> TOFIX: GLEICH ODER UNGLEICH NULL????????

   IX = WC%IXW
   IY = WC%IYW
   IZ = WC%IZW

   IOR0   = WC%IOR
   BTYPE0 = WC%BTYPE

  ! ICG = L%CELL_INDEX_PTR(IX, IY, IZ)
  ! IWG = L%WALL_INDEX_PTR(ICG, IOR0)
  ! IF (WC%BOUNDARY_TYPE /= INTERPOLATED_BOUNDARY) WC%BTYPE=NEUMANN
  ! IF (L%CELL_STATE(IX, IY, IZ) == NSCARC_CELL_SOLID) WC%BTYPE=NEUMANN

   IF (WC%BOUNDARY_TYPE == SOLID_BOUNDARY) WC%BTYPE=NEUMANN

ENDDO

END SUBROUTINE SCARC_IDENTIFY_INTERNAL_NEUMANNS


!> -------------------------------------------------------------------------------------------------
!> Setup all necessary information for a wall cell with neighbor
!> -------------------------------------------------------------------------------------------------
INTEGER FUNCTION SCARC_COUNT_EXTERNAL_WALL_CELLS(NM, NL)
INTEGER, INTENT(IN) :: NM, NL
INTEGER :: IXC, IYC, IZC
INTEGER :: IXF, IYF, IZF
INTEGER :: NLF, NLC
INTEGER :: IWC, ICF(4)=0, IWF(4)=0, IOR0
TYPE (SCARC_LEVEL_TYPE), POINTER :: LF=>NULL(), LC=>NULL()

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

   !> IOR = 1
   IOR0 = 1
   DO IZC = 1, LC%NZ
      IZF = 2*IZC - 1
      ICF(1) = LF%CELL_INDEX_PTR(1  , IYF  , IZF  )
      ICF(2) = LF%CELL_INDEX_PTR(1  , IYF  , IZF+1)
      IF (IS_EXTERNAL_WALLCELL(IOR0, ICF, 2, NM, NLF)) IWC = IWC + 1
   ENDDO

   !> IOR = -1
   IOR0 = -1
   DO IZC = 1, LC%NZ
      IZF = 2*IZC - 1
      ICF(1) = LF%CELL_INDEX_PTR(LF%NX, IYF  , IZF  )
      ICF(2) = LF%CELL_INDEX_PTR(LF%NX, IYF  , IZF+1)
      IF (IS_EXTERNAL_WALLCELL(IOR0, ICF, 2, NM, NLF)) IWC = IWC + 1
   ENDDO

   !> IOR = 2
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

   !> IOR = -2
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

   !> IOR = 3
   IOR0 = 3
   DO IXC = 1, LC%NX
      IXF = 2*IXC - 1
      ICF(1) = LF%CELL_INDEX_PTR(IXF    , IYF  , 1)
      ICF(2) = LF%CELL_INDEX_PTR(IXF+1  , IYF  , 1)
      IF (IS_EXTERNAL_WALLCELL(IOR0, ICF, 2, NM, NLF)) IWC = IWC + 1
   ENDDO

   !> IOR = -3
   IOR0 = -3
   DO IXC = 1, LC%NX
      IXF = 2*IXC - 1
      ICF(1) = LF%CELL_INDEX_PTR(IXF  , IYF  , LF%NZ)
      ICF(2) = LF%CELL_INDEX_PTR(IXF+1, IYF  , LF%NZ)
      IF (IS_EXTERNAL_WALLCELL(IOR0, ICF, 2, NM, NLF)) IWC = IWC + 1
   ENDDO

ELSE

   !> IOR = 1
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

   !> IOR = -1
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

   !> IOR = 2
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

   !> IOR = -2
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

   !> IOR = 3
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

   !> IOR = -3
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

!> -------------------------------------------------------------------------------------------------
!> Setup all necessary information for a wall cell with neighbor
!> -------------------------------------------------------------------------------------------------
INTEGER FUNCTION SCARC_COUNT_INTERNAL_WALL_CELLS(NM, NL)
INTEGER, INTENT(IN) :: NM, NL
INTEGER :: IWC, NW_INT
INTEGER :: IC, IO
INTEGER :: I, J, K
TYPE (SCARC_LEVEL_TYPE), POINTER :: LF=>NULL(), LC=>NULL()
TYPE (SCARC_OBST_TYPE), POINTER :: OB=>NULL()
TYPE (SCARC_DISCRET_TYPE), POINTER :: DF=>NULL(), DC=>NULL()

LC => SCARC(NM)%LEVEL(NL)
LF => SCARC(NM)%LEVEL(NL-1)

SELECT CASE(TYPE_DISCRET)
   CASE (NSCARC_DISCRET_STRUCTURED)
      DC => LC%SD
      DF => LF%SD
   CASE (NSCARC_DISCRET_UNSTRUCTURED)
      DC => LC%UD
      DF => LF%UD
END SELECT

LC%N_OBST = LF%N_OBST                   !> Number of obstructions is the same on all levels

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

   !> Analyze IOR = 1
   I = OB%I1
   DO K = OB%K1+1, OB%K2
      DO J = OB%J1+1, OB%J2
         IC = DC%CELL_NUMBER(I, J, K)
         IF (IC == NSCARC_UNDEFINED_INT) CYCLE
         NW_INT = NW_INT + 1
      ENDDO
   ENDDO

   !> Analyze IOR = -1
   I = OB%I2
   DO K = OB%K1+1, OB%K2
      DO J = OB%J1+1, OB%J2
         IC = DC%CELL_NUMBER(I+1, J, K)
         IF (IC == NSCARC_UNDEFINED_INT) CYCLE
         NW_INT = NW_INT + 1
      ENDDO
   ENDDO

   !> Analyze IOR = 2
   J = OB%J1
   DO K = OB%K1+1, OB%K2
      DO I = OB%I1+1, OB%I2
         IC = DC%CELL_NUMBER(I, J, K)
         IF (IC == NSCARC_UNDEFINED_INT) CYCLE
         NW_INT = NW_INT + 1
      ENDDO
   ENDDO

   !> Analyze IOR = -2
   J = OB%J2
   DO K = OB%K1+1, OB%K2
      DO I = OB%I1+1, OB%I2
         IC = DC%CELL_NUMBER(I, J+1, K)
         IF (IC == NSCARC_UNDEFINED_INT) CYCLE
         NW_INT = NW_INT + 1
      ENDDO
   ENDDO

   !> Analyze IOR = 3
   K = OB%K1
   DO J = OB%J1+1, OB%J2
      DO I = OB%I1+1, OB%I2
         IC = DC%CELL_NUMBER(I, J, K)
         IF (IC == NSCARC_UNDEFINED_INT) CYCLE
         NW_INT = NW_INT + 1
      ENDDO
   ENDDO

   !> Analyze IOR = -3
   K = OB%K2
   DO J = OB%J1+1, OB%J2
      DO I = OB%I1+1, OB%I2
         IC = DC%CELL_NUMBER(I, J, K+1)
         IF (IC == NSCARC_UNDEFINED_INT) CYCLE
         NW_INT = NW_INT + 1
      ENDDO
   ENDDO

ENDDO

SCARC_COUNT_INTERNAL_WALL_CELLS = NW_INT
END FUNCTION SCARC_COUNT_INTERNAL_WALL_CELLS



!> -------------------------------------------------------------------------------------------------
!> Count external wall cells on face IOR
!> -------------------------------------------------------------------------------------------------
LOGICAL FUNCTION IS_EXTERNAL_WALLCELL(IOR0, ICF, NCNT, NM, NL)
INTEGER, INTENT(IN) :: IOR0, NCNT, NM, NL
INTEGER, DIMENSION(:), INTENT(IN) :: ICF
INTEGER :: I, IWF_LAST, IWF(4)=0
REAL(EB) :: BSUM
TYPE (SCARC_LEVEL_TYPE), POINTER :: L=>NULL()
TYPE (SCARC_DISCRET_TYPE), POINTER :: D=>NULL()

IS_EXTERNAL_WALLCELL = .FALSE.

L => SCARC(NM)%LEVEL(NL)
SELECT CASE(TYPE_DISCRET)
   CASE (NSCARC_DISCRET_STRUCTURED)
      D => L%SD
   CASE (NSCARC_DISCRET_UNSTRUCTURED)
      D => L%UD
END SELECT


DO I = 1, NCNT
   IWF(I) = L%WALL_INDEX_PTR(ICF(I), -IOR0)
ENDDO

BSUM = 0.0_EB
IWF_LAST = 0

DO I = 1, NCNT
   IF (IWF(I)>0) THEN
      BSUM = BSUM + REAL(D%WALL(IWF(I))%BTYPE,EB)
      IWF_LAST = IWF(I)
   ENDIF
ENDDO

IF (IWF_LAST == 0) RETURN
IF (ABS(BSUM/REAL(NCNT,EB) - REAL(D%WALL(IWF_LAST)%BTYPE,EB)) < 1E-12) THEN
   IS_EXTERNAL_WALLCELL = .TRUE.
   RETURN
ELSE
   CALL SCARC_SHUTDOWN(NSCARC_ERROR_BOUNDARY_SUM, SCARC_NONE, IOR0)
ENDIF

END FUNCTION IS_EXTERNAL_WALLCELL


!> -------------------------------------------------------------------------------------------------
!> Setup all necessary information for a wall cell with neighbor
!> -------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_WALLCELL_NEIGHBOR(NX1, NX2, NY1, NY2, NZ1, NZ2, IWG, IOR0, NM, NOM, NL)
INTEGER, INTENT(IN) :: NX1, NX2, NY1, NY2, NZ1, NZ2
INTEGER, INTENT(IN) :: IWG, IOR0, NM, NOM, NL
INTEGER :: NOMX, NOMY, NOMZ
INTEGER :: ICG, ICO, ICE, IWL, ICPL, IX, IY, IZ, JL
TYPE (SCARC_NEIGHBOR_TYPE), POINTER :: OS=>NULL()
TYPE (SCARC_LEVEL_TYPE), POINTER :: L=>NULL(), OL=>NULL()
TYPE (SCARC_DISCRET_TYPE), POINTER :: D=>NULL(), OD=>NULL()
TYPE (SCARC_POINTER_TYPE), POINTER :: PTR=>NULL(), OPTR=>NULL()
TYPE (SCARC_WALL_TYPE), POINTER :: WC=>NULL()

L  => SCARC(NM)%LEVEL(NL)
OS => SCARC(NM)%OSCARC(NOM)
OL => OS%LEVEL(NL)

SELECT CASE(TYPE_DISCRET)
   CASE (NSCARC_DISCRET_STRUCTURED)
      D  => L%SD
      OD => OL%SD
   CASE (NSCARC_DISCRET_UNSTRUCTURED)
      D =>  L%UD
      OD => OL%UD
END SELECT

PTR  => D%PTR
OPTR => OD%PTR

!> store basic information about neighbor and orientation
OL%IOR  = IOR0

ICO  = PTR%ICO
ICE  = PTR%ICE

ICG  = OPTR%ICG
IWL  = OPTR%IWL

OL%NCPLS = OL%NCPL
IF (OL%NCG == OL%NWL) THEN
   OL%NCPLR = OL%NCPL
ELSE IF (OL%NCG == 2*OL%NWL) THEN
   OL%NCPLR =  1
ELSE IF (OL%NWL == 2*OL%NCG) THEN
   OL%NCPLR =  2
ENDIF

!> set neighboring coordinates
WC => D%WALL(IWG)

WC%IXN(1) = NX1
WC%IXN(2) = NX2
WC%IYN(1) = NY1
WC%IYN(2) = NY2
WC%IZN(1) = NZ1
WC%IZN(2) = NZ2

!> allocate pointer arrays for extended, ghost and neighboring cells
CALL SCARC_ALLOCATE_INT1(WC%ICE, 1, OL%NCPL, NSCARC_INIT_UNDEF, 'ICE')
CALL SCARC_ALLOCATE_INT1(WC%ICG, 1, OL%NCPL, NSCARC_INIT_UNDEF, 'ICG')

IWL = IWL + 1
ICO = ICO + 1

IF (OL%SUBDIVISION(1, IOR0) == NSCARC_ZERO_INT) OL%SUBDIVISION(1, IOR0) =  IWL
OL%SUBDIVISION(2, IOR0) = OL%SUBDIVISION(2, IOR0) + 1

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

!> store information about overlapped cells and set mapping arrays
ICPL = 0
DO IZ = NZ1, NZ2
   DO IY = NY1, NY2
      DO IX = NX1, NX2

         ICPL = ICPL + 1
         ICG  = ICG  + 1
         ICE  = ICE  + 1

         WC%ICE(ICPL)  = ICE                              !> number of extended grid cell
         WC%ICG(ICPL)  = ICG                              !> number of ghost grid cell

         PTR%ICE_TO_IWG(ICE)  = IWG                       !> map extended cell to global wall cell
         PTR%ICE_TO_IWL(ICE)  = IWL                       !> map extended cell to local wall cell
         PTR%ICE_TO_ICG(ICE)  = ICG                       !> map extended cell to ghost cell
         PTR%ICE_TO_ICG(ICE)  = ICG                       !> map extended cell to ghost cell

         OPTR%ICG_TO_IWG(ICG) = IWG                       !> map ghost cell to global wall cell
         OPTR%ICG_TO_ICO(ICG) = ICO                       !> map ghost cell to extended grid cell
         OPTR%ICG_TO_ICE(ICG) = D%WALL(IWG)%ICE(ICPL)     !> map ghost cell to extended grid cell

      ENDDO
   ENDDO
ENDDO

WC%ICO = ICO                                                !> number of overlapping cell
WC%IWL = IWL                                                !> number of local wall cell

PTR%ICO  = ICO                                            !> store overlapping cell counter
PTR%ICE  = ICE                                            !> store extended cell counter

OPTR%IWL = IWL                                            !> store local wall cell pointer
OPTR%ICG = ICG                                            !> store ghost cell counter

OPTR%IWL_TO_IWG(IWL) = IWG                                !> map local wall cell to global wall cell
OPTR%IWL_TO_ICO(IWL) = ICO                                !> map local wall cell to internal grid cell
!OPTR%IWL_TO_ICW(IWL) = D%WALL(IWG)%ICW                   !> map local wall cell to internal grid cell (AMG only)

END SUBROUTINE SCARC_SETUP_WALLCELL_NEIGHBOR


!> -------------------------------------------------------------------------------------------------
!> Allocate needed administration arrays for SCARC(NM)%LEVEL(NL)
!> -------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_POINTERS(NM, NL)
INTEGER, INTENT(IN) :: NM, NL
TYPE (SCARC_LEVEL_TYPE), POINTER :: L=>NULL()
TYPE (SCARC_DISCRET_TYPE), POINTER :: D=>NULL()
TYPE (SCARC_POINTER_TYPE), POINTER :: PTR=>NULL()

L   => SCARC(NM)%LEVEL(NL)
SELECT CASE(TYPE_DISCRET)
   CASE (NSCARC_DISCRET_STRUCTURED)
      D => L%SD
   CASE (NSCARC_DISCRET_UNSTRUCTURED)
      D => L%UD
END SELECT
PTR => D%PTR

IF (L%NCE > L%NC) THEN
   CALL SCARC_ALLOCATE_INT1(PTR%ICE_TO_IWG, L%NC+1, L%NCE, NSCARC_INIT_ZERO, 'ICE_TO_IWG')
   CALL SCARC_ALLOCATE_INT1(PTR%ICE_TO_IWL, L%NC+1, L%NCE, NSCARC_INIT_ZERO, 'ICL_TO_IWG')
   CALL SCARC_ALLOCATE_INT1(PTR%ICE_TO_ICG, L%NC+1, L%NCE, NSCARC_INIT_ZERO, 'ICE_TO_ICG')
ENDIF

END SUBROUTINE SCARC_SETUP_POINTERS


!> -------------------------------------------------------------------------------------------------
!> Allocate needed administration arrays for SCARC(NM)%OSCARC(NOM)%LEVEL(NL)
!> -------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_OPOINTERS(NM, NL, NOM)
INTEGER, INTENT(IN) :: NM, NL, NOM
TYPE (SCARC_LEVEL_TYPE), POINTER :: OL=>NULL()
TYPE (SCARC_POINTER_TYPE), POINTER :: OPTR=>NULL()
TYPE (SCARC_DISCRET_TYPE), POINTER :: OD=>NULL()

OL   => SCARC(NM)%OSCARC(NOM)%LEVEL(NL)
SELECT CASE(TYPE_DISCRET)
   CASE (NSCARC_DISCRET_STRUCTURED)
      OD => OL%SD
   CASE (NSCARC_DISCRET_UNSTRUCTURED)
      OD => OL%UD
END SELECT
OPTR => OD%PTR

CALL SCARC_ALLOCATE_INT1(OPTR%IWL_TO_IWG, 1, OL%NWL, NSCARC_INIT_ZERO, 'IWL_TO_IWG')
CALL SCARC_ALLOCATE_INT1(OPTR%IWL_TO_ICO, 1, OL%NWL, NSCARC_INIT_ZERO, 'IWL_TO_ICO')

CALL SCARC_ALLOCATE_INT1(OPTR%ICG_TO_IWG, 1, OL%NCG, NSCARC_INIT_ZERO, 'ICG_TO_IWG')
CALL SCARC_ALLOCATE_INT1(OPTR%ICG_TO_ICO, 1, OL%NCG, NSCARC_INIT_ZERO, 'ICG_TO_ICO')
CALL SCARC_ALLOCATE_INT1(OPTR%ICG_TO_ICE, 1, OL%NCG, NSCARC_INIT_ZERO, 'ICG_TO_ICE')

!CALL SCARC_ALLOCATE_INT1(OL%MAP%IWL_TO_ICW, 1, OL%NWL, NSCARC_INIT_ZERO, 'IWL_TO_ICW')   !> TOFIX: still needed?
!IF (TYPE_DISCRET == NSCARC_DISCRET_STRUCTURED.AND.OL%NC>0) &                             !> TOFIX: huge!
!   CALL SCARC_ALLOCATE_INT1(OPTR%ICN_TO_ICE, 1, OL%NC, NSCARC_INIT_ZERO,'ICN_TO_ICE')

END SUBROUTINE SCARC_SETUP_OPOINTERS


!> ----------------------------------------------------------------------------------------------------
!> Check divisibility by 2 of a given number of elements (in one grid direction)
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_CHECK_DIVISIBILITY(NN, CDIR)
INTEGER, INTENT(IN) :: NN
CHARACTER(*) , INTENT(IN) :: CDIR
IF (MOD(NN,2) /= 0) CALL SCARC_SHUTDOWN(NSCARC_ERROR_GRID_NUMBER, CDIR, NSCARC_NONE)
END SUBROUTINE SCARC_CHECK_DIVISIBILITY


!> ----------------------------------------------------------------------------------------------------
!> Set wall cell information on coarse level
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_FACE_DIMENSIONS(IOR0, IWG, NM, NL)
INTEGER, INTENT(IN)    :: IOR0, NM, NL
INTEGER, INTENT(INOUT) :: IWG
INTEGER :: INBR
TYPE (SCARC_LEVEL_TYPE), POINTER :: LF=>NULL(), LC=>NULL()
TYPE (SCARC_FACE_TYPE), POINTER :: FF=>NULL(), FC=>NULL()

!> reference coarse and fine LEVEL type
LF => SCARC(NM)%LEVEL(NL-1)
LC => SCARC(NM)%LEVEL(NL)

FF => SCARC(NM)%LEVEL(NL-1)%FACE(IOR0)
FC => SCARC(NM)%LEVEL(NL)%FACE(IOR0)

!> initialize FACE type for coarser mesh
FC%IWG_PTR = IWG

FC%N_NEIGHBORS = FF%N_NEIGHBORS
IF (FC%N_NEIGHBORS /= 0) CALL SCARC_ALLOCATE_INT1(FC%NEIGHBORS, 1, FC%N_NEIGHBORS, NSCARC_INIT_NONE, 'FACE_NEIGHBORS')
DO INBR= 1, FC%N_NEIGHBORS
   FC%NEIGHBORS(INBR) = FF%NEIGHBORS(INBR)
ENDDO

SELECT CASE (ABS(IOR0))

   CASE (1)
      FC%DH => LC%DXL
      FC%NFX = 1                                              !> no extension in x-direction
      IF (.NOT.TWO_D) THEN                                    !> only subdivide y-direction in 3D-case
         CALL SCARC_CHECK_DIVISIBILITY(FF%NFY, 'Y')
         FC%NFY = FF%NFY/2
      ELSE
         FC%NFY = FF%NFY
      ENDIF
      CALL SCARC_CHECK_DIVISIBILITY(FF%NFZ, 'Z')              !> number of z-cells divisible by 2?
      FC%NFZ = FF%NFZ/2
      FC%NFC = LC%NX                                          !> number of cells between opposite faces
   CASE (2)
      FC%DH => LC%DYL
      CALL SCARC_CHECK_DIVISIBILITY(FF%NFX, 'X')              !> number of x-cells divisible by 2?
      FC%NFX = FF%NFX/2
      FC%NFY = 1                                              !> no extension in y-direction
      CALL SCARC_CHECK_DIVISIBILITY(FF%NFZ, 'Z')              !> number of z-cells divisible by 2?
      FC%NFZ = FF%NFZ/2
      FC%NFC = LC%NY                                          !> number of cells between opposite faces
   CASE (3)
      FC%DH => LC%DZL
      CALL SCARC_CHECK_DIVISIBILITY(FF%NFX, 'X')              !> number of x-cells divisible by 2?
      FC%NFX = FF%NFX/2
      IF (.NOT.TWO_D) THEN                                    !> only subdivide y-direction in 3D-case
         CALL SCARC_CHECK_DIVISIBILITY(FF%NFY, 'Y')
         FC%NFY = FF%NFY/2
      ELSE
         FC%NFY = FF%NFY
      ENDIF
      FC%NFZ = 1                                              !> no extension in y-direction
      FC%NFC = LC%NZ                                          !> number of cells between opposite faces
END SELECT

FC%NFW = FC%NFX * FC%NFY * FC%NFZ                             !> get number of wall cells for that face
IWG = IWG + FC%NFW                                            !> increase global wall cell counter

END SUBROUTINE SCARC_SETUP_FACE_DIMENSIONS


!> ----------------------------------------------------------------------------------------------------
!> Set wall cell information on coarse level
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_FACE(IOR0, IWC, IREFINE, NM, NL)
INTEGER, INTENT(INOUT) :: IWC
INTEGER, INTENT(IN) :: NM, NL
INTEGER, INTENT(IN) :: IOR0, IREFINE
INTEGER :: IWF(4) , IBCF(4), NOMF(4)
INTEGER :: NCPL
INTEGER :: IX,  IY,  IZ, I
INTEGER :: NX1, NY1, NZ1
INTEGER :: NX2, NY2, NZ2
INTEGER :: IX1, IY1, IZ1
INTEGER :: IX2, IY2, IZ2
INTEGER :: IDIFF, JDIFF, KDIFF
TYPE (SCARC_LEVEL_TYPE), POINTER :: LF=>NULL(), LC=>NULL(), OLC=>NULL()
TYPE (SCARC_FACE_TYPE), POINTER :: FF=>NULL(), FC=>NULL()
TYPE (SCARC_DISCRET_TYPE), POINTER :: DF=>NULL(), DC=>NULL()
TYPE (SCARC_WALL_TYPE), POINTER, DIMENSION(:) :: WC=>NULL(), WF=>NULL()

!> reference coarse and fine LEVEL type
LC => SCARC(NM)%LEVEL(NL)
LF => SCARC(NM)%LEVEL(NL-1)

SELECT CASE(TYPE_DISCRET)
   CASE (NSCARC_DISCRET_STRUCTURED)
      DC => LC%SD
      DF => LF%SD
   CASE (NSCARC_DISCRET_UNSTRUCTURED)
      DC => LC%UD
      DF => LF%UD
END SELECT

!> reference coarse and fine WALL type
WC => DC%WALL
WF => DF%WALL

!> reference coarse and fine FACE type
FC => LC%FACE(IOR0)
FF => LF%FACE(IOR0)

!> set coordinate dimensions for correspoding face
SELECT CASE (ABS(IOR0))
   CASE (1)
      IF (IOR0 > 0) THEN                                           !> set dimensions for wall cell counting
         NX1 = 0
         NX2 = 0
      ELSE
         NX1 = LC%NX+1
         NX2 = LC%NX+1
      ENDIF
      NY1 = 1
      NY2 = LC%NY
      NZ1 = 1
      NZ2 = LC%NZ
   CASE (2)
      NX1 = 1                                                      !> set dimensions for wall cell counting
      NX2 = LC%NX
      IF (IOR0 > 0) THEN
         NY1 = 0
         NY2 = 0
      ELSE
         NY1 = LC%NY+1
         NY2 = LC%NY+1
      ENDIF
      NZ1 = 1
      NZ2 = LC%NZ
   CASE (3)
      NX1 = 1                                                      !> set dimensions for wall cell counting
      NX2 = LC%NX
      NY1 = 1
      NY2 = LC%NY
      IF (IOR0 > 0) THEN
         NZ1 = 0
         NZ2 = 0
      ELSE
         NZ1 =LC%NZ+1
         NZ2 =LC%NZ+1
      ENDIF
END SELECT

!>
!> Loop over all wall cells of face IOR0
!>
DO IZ = NZ1, NZ2
   DO IY = NY1, NY2
      DO IX = NX1, NX2

         !> Set orientation of neiboring face, indices of ghost and adjacent cell for coarse IW
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

         !> ------------------------------------------------------------
         !> 2D-version
         !> ------------------------------------------------------------
         IF (TWO_D) THEN

            !> determine fine IW's, which must be merged to one coarse IW
            SELECT CASE (ABS(IOR0))
               CASE ( 1)
                  IWF(1) = FF%IWG_PTR + 2*(IZ-1)
               CASE ( 2)
                  IWF(1) = FF%IWG_PTR + 2*(IZ-1)*LF%NX + 2*(IX - 1)
               CASE ( 3)
                  IWF(1) = FF%IWG_PTR + 2*(IX-1)
            END SELECT
            IWF(2) = IWF(1)+1

            !> set fine cell neighbors (they must be the same for all fine IW's)
            NOMF(1) = WF(IWF(1))%NOM
            NOMF(2) = WF(IWF(2))%NOM
            IF (NOMF(1) /= NOMF(2)) CALL SCARC_SHUTDOWN(NSCARC_ERROR_NEIGHBOR_TYPE, SCARC_NONE, NOMF(1))

            WC(IWC)%NOM = NOMF(1)

            !> set corresponding pressure_bc_index on coarser level
            IBCF(1) = WF(IWF(1))%BTYPE
            IBCF(2) = WF(IWF(2))%BTYPE
            IF (IBCF(1) == INTERNAL .OR. IBCF(2) == INTERNAL) THEN
               WC(IWC)%BTYPE = INTERNAL
            ELSE IF (IBCF(1) == DIRICHLET .OR. IBCF(2) == DIRICHLET) THEN
               WC(IWC)%BTYPE = DIRICHLET
            ELSE
               WC(IWC)%BTYPE = NEUMANN
            ENDIF

            !> set corresponding pressure_bc_index on coarser level
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

            !> in case of an internal boundary set neighboring WALL cells
            IF (NOMF(1) > 0) THEN

               OLC => SCARC(NM)%OSCARC(NOMF(1))%LEVEL(NL)

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

               !>!
               !> Allocate and specify ICN and ICE arrays for OC
               !>!
               NCPL = (IZ2-IZ1+1)*(IX2-IX1+1)
               OLC%NCPL = NCPL

               CALL SCARC_SETUP_WALLCELL_NEIGHBOR(IX1, IX2, 1, 1, IZ1, IZ2, IWC, IOR0, NM, NOMF(1), NL)

            ENDIF


         !> ------------------------------------------------------------
         !> 3D-version
         !> ------------------------------------------------------------
         ELSE

            !> determine fine IW's, which must be merged to one coarse IW
            SELECT CASE (ABS(IOR0))
               CASE (1)
                  IWF(1) = FF%IWG_PTR + (2*IZ-2)*LF%NY + 2*IY - 2
                  IWF(3) = FF%IWG_PTR + (2*IZ-1)*LF%NY + 2*IY - 2
               CASE (2)
                  IWF(1) = FF%IWG_PTR + (2*IZ-2)*LF%NX + 2*IX - 2
                  IWF(3) = FF%IWG_PTR + (2*IZ-1)*LF%NX + 2*IX - 2
               CASE (3)
                  IWF(1) = FF%IWG_PTR + (2*IY-2)*LF%NX + 2*IX - 2
                  IWF(3) = FF%IWG_PTR + (2*IY-1)*LF%NX + 2*IX - 2
            END SELECT
            IWF(2) = IWF(1)+1
            IWF(4) = IWF(3)+1

            !> set fine cell neighbors (they must be the same for all fine IW's)
            DO I=1,4
               NOMF(I) = WF(IWF(I))%NOM
            ENDDO

            IF (NOMF(1) /= NOMF(2) .OR. NOMF(1) /= NOMF(3) .OR. NOMF(1) /= NOMF(4)) &
               CALL SCARC_SHUTDOWN(NSCARC_ERROR_NEIGHBOR_TYPE, SCARC_NONE, IOR0)
            WC(IWC)%NOM = NOMF(1)

            !> set corresponding pressure_bc_index on coarser level
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

            !> set corresponding pressure_bc_index on coarser level
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

            !> in case of an internal boundary set WALL(10:15,IWC)
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

               !> Allocate and specify ICN and ICE arrays for OLC
               NCPL = (IZ2-IZ1+1)*(IY2-IY1+1)*(IX2-IX1+1)
               OLC%NCPL = NCPL

               CALL SCARC_SETUP_WALLCELL_NEIGHBOR(IX1, IX2, IY1, IY2, IZ1, IZ2, IWC, IOR0, NM, NOMF(1), NL)

            ENDIF
         ENDIF
         IWC = IWC + 1
      ENDDO
   ENDDO
ENDDO

FF%IOFFSET_WALL = FF%IOFFSET_WALL + FF%NFW


END SUBROUTINE SCARC_SETUP_FACE


!> -------------------------------------------------------------------------------------------------
!> Set analogues to I_MIN, IMAX, J_MIN, J_MAX, K_MIN, K_MAX on coarse level (only GMG!)
!> -------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_EXCHANGE_DIMENSIONS(IREFINE, NM, NOM, NL)
INTEGER, INTENT(IN) :: IREFINE, NM, NOM, NL
INTEGER :: IMIN, IMAX, JMIN, JMAX, KMIN, KMAX, IW
LOGICAL :: FOUND
TYPE (SCARC_LEVEL_TYPE), POINTER :: L=>NULL(), OL=>NULL()
TYPE (SCARC_NEIGHBOR_TYPE), POINTER :: OS=>NULL()
TYPE (SCARC_DISCRET_TYPE), POINTER :: D=>NULL()

L  => SCARC(NM)%LEVEL(NL)
SELECT CASE(TYPE_DISCRET)
   CASE (NSCARC_DISCRET_STRUCTURED)
      D => L%SD
   CASE (NSCARC_DISCRET_UNSTRUCTURED)
      D => L%UD
END SELECT

OS => SCARC(NM)%OSCARC(NOM)
OL => SCARC(NM)%OSCARC(NOM)%LEVEL(NL)

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

OS%NIC_R = 0
FOUND = .FALSE.

SEARCH_LOOP: DO IW=1, OL%NCG

   !> neighborship structure already known from finest level
   IF (D%WALL(IW)%NOM/=NOM) CYCLE SEARCH_LOOP
   OS%NIC_R = OS%NIC_R + 1
   FOUND = .TRUE.

   SELECT CASE (D%WALL(IW)%IOR)
      CASE ( 1)
         IMIN=MAX(IMIN,D%WALL(NM)%IXN(1)-1)
      CASE (-1)
         IMAX=MIN(IMAX,D%WALL(NM)%IXN(2)+1)
      CASE ( 2)
         JMIN=MAX(JMIN,D%WALL(NM)%IYN(1)-1)
      CASE (-2)
         JMAX=MIN(JMAX,D%WALL(NM)%IYN(2)+1)
      CASE ( 3)
         KMIN=MAX(KMIN,D%WALL(NM)%IZN(1)-1)
      CASE (-3)
         KMAX=MIN(KMAX,D%WALL(NM)%IZN(2)+1)
   END SELECT
ENDDO SEARCH_LOOP

N_EXCHANGE = N_EXCHANGE+1

END SUBROUTINE SCARC_SETUP_EXCHANGE_DIMENSIONS

!> ----------------------------------------------------------------------------------------------------
!> Allocate several global structures for data exchange
!> --------------------------------------------------------------------------------------------------ss
SUBROUTINE SCARC_SETUP_GLOBALS
INTEGER :: NM, NP

!> Allocate and preset counter and displacement vector for global data exchanges
CALL SCARC_ALLOCATE_INT1 (COUNTS, 0, N_MPI_PROCESSES-1, NSCARC_INIT_ZERO, 'COUNTS')
CALL SCARC_ALLOCATE_INT1 (DISPLS, 0, N_MPI_PROCESSES-1, NSCARC_INIT_ZERO, 'DISPLS')

CALL SCARC_ALLOCATE_INT1 (LOCAL_INT , 1, NMESHES, NSCARC_INIT_ZERO, 'LOCAL_INT')
CALL SCARC_ALLOCATE_REAL1(LOCAL_REAL, 1, NMESHES, NSCARC_INIT_ZERO, 'LOCAL_REAL')

!> Get number of data to send per process
DO NP = 0, N_MPI_PROCESSES-1
   DO NM = 1, NMESHES
      IF (PROCESS(NM)==NP) COUNTS(NP) = COUNTS(NP) + 1
   ENDDO
ENDDO

!> Get displacements on communication vector for all meshes
DO NP = 1, N_MPI_PROCESSES-1
   DISPLS(NP) = COUNTS(NP-1) + DISPLS(NP-1)
ENDDO

CALL SCARC_SETUP_DIMENSIONS(NLEVEL_MIN)

END SUBROUTINE SCARC_SETUP_GLOBALS


!> --------------------------------------------------------------------------------------------------
!> Broadcast local integer information to all meshes 
!> Process received data due to specified NTYPE (SUM/PRODUCT/FIRST/LAST/MAX/MIN)
!> Corresponding local data is passed in array LOCAL_INT by default
!> --------------------------------------------------------------------------------------------------
INTEGER FUNCTION SCARC_BROADCAST_INT(NTYPE)
INTEGER, INTENT(IN) :: NTYPE

IF (N_MPI_PROCESSES > 1) &
   CALL MPI_ALLGATHERV(MPI_IN_PLACE,1,MPI_INTEGER,LOCAL_INT,COUNTS,DISPLS, MPI_INTEGER,MPI_COMM_WORLD,IERROR)

GLOBAL_INT = 0
SELECT CASE (NTYPE)
   CASE (NSCARC_BROADCAST_SUM)
      GLOBAL_INT = SUM(LOCAL_INT(1:NMESHES))
   CASE (NSCARC_BROADCAST_PRODUCT)
      GLOBAL_INT = PRODUCT(LOCAL_INT(1:NMESHES))
   CASE (NSCARC_BROADCAST_MIN)
      GLOBAL_INT = MINVAL(LOCAL_INT(1:NMESHES))
   CASE (NSCARC_BROADCAST_MAX)
      GLOBAL_INT = MAXVAL(LOCAL_INT(1:NMESHES))
   CASE (NSCARC_BROADCAST_FIRST)
      GLOBAL_INT = LOCAL_INT(1)
   CASE (NSCARC_BROADCAST_LAST)
      GLOBAL_INT = LOCAL_INT(NMESHES)
END SELECT

SCARC_BROADCAST_INT = GLOBAL_INT
RETURN
END FUNCTION SCARC_BROADCAST_INT


!> --------------------------------------------------------------------------------------------------
!> Broadcast local double precision information to all meshes 
!> Process received data due to specified NTYPE (SUM/PRODUCT/FIRST/LAST/MAX/MIN)
!> Corresponding local data is passed in array LOCAL_INT by default
!> --------------------------------------------------------------------------------------------------
DOUBLE PRECISION FUNCTION SCARC_BROADCAST_REAL(NTYPE)
INTEGER, INTENT(IN) :: NTYPE

IF (N_MPI_PROCESSES > 1) THEN
   CALL MPI_ALLGATHERV(MPI_IN_PLACE,1,MPI_DOUBLE_PRECISION,LOCAL_REAL,COUNTS,DISPLS,&
                       MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,IERROR)
ENDIF

GLOBAL_REAL = 0.0_EB
SELECT CASE (NTYPE)
   CASE (NSCARC_BROADCAST_SUM)
      GLOBAL_REAL = SUM(LOCAL_REAL(1:NMESHES))
   CASE (NSCARC_BROADCAST_PRODUCT)
      GLOBAL_REAL = PRODUCT(LOCAL_REAL(1:NMESHES))
   CASE (NSCARC_BROADCAST_MIN)
      GLOBAL_REAL = MINVAL(LOCAL_REAL(1:NMESHES))
   CASE (NSCARC_BROADCAST_MAX)
      GLOBAL_REAL = MAXVAL(LOCAL_REAL(1:NMESHES))
   CASE (NSCARC_BROADCAST_FIRST)
      GLOBAL_REAL = LOCAL_REAL(1)
   CASE (NSCARC_BROADCAST_LAST)
      GLOBAL_REAL = LOCAL_REAL(NMESHES)
END SELECT

SCARC_BROADCAST_REAL = GLOBAL_REAL
RETURN
END FUNCTION SCARC_BROADCAST_REAL


!> -----------------------------------------------------------------------------
!> Get information about global numbers of unknowns for unstructured case
!> -----------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_DIMENSIONS(NL)
INTEGER, INTENT(IN) :: NL
INTEGER :: NM, NM2
TYPE (SCARC_LEVEL_TYPE), POINTER :: L=>NULL()
TYPE (SCARC_DISCRET_TYPE), POINTER :: D=>NULL()

!> Preset communication array LOCAL_INT with local numbers of cells for all meshes depending on type of discretization
MESHES_LOOP1: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   L => SCARC(NM)%LEVEL(NL)
   SELECT CASE(TYPE_DISCRET)
      CASE (NSCARC_DISCRET_STRUCTURED)
         D => L%SD                                 
      CASE (NSCARC_DISCRET_UNSTRUCTURED)
         D => L%UD                                
   END SELECT

   DO NM2 = LOWER_MESH_INDEX, UPPER_MESH_INDEX
      LOCAL_INT(NM2) = D%NC_LOCAL(NM2)
   ENDDO

ENDDO MESHES_LOOP1

!> Broadcast number of local mesh cells on level NL to all and build global sum
N_CELLS_GLOBAL(NL) = SCARC_BROADCAST_INT(NSCARC_BROADCAST_SUM)

!> Store information on local and global cells numbers on data structure of corresponding discretization type
MESHES_LOOP2: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   L => SCARC(NM)%LEVEL(NL)
   SELECT CASE(TYPE_DISCRET)
      CASE (NSCARC_DISCRET_STRUCTURED)
         D => L%SD                               
      CASE (NSCARC_DISCRET_UNSTRUCTURED)
         D => L%UD                              
   END SELECT

   D%NC_LOCAL(1:NMESHES) = LOCAL_INT(1:NMESHES)
   D%NC_GLOBAL = SUM(LOCAL_INT(1:NMESHES))

   !> compute offset between local grid numberings
   IF (NMESHES > 1) THEN
      DO NM2=2,NMESHES
         D%NC_OFFSET(NM2) = D%NC_OFFSET(NM2-1) + D%NC_LOCAL(NM2-1)
      ENDDO
   ENDIF

ENDDO MESHES_LOOP2

IF (NL == NLEVEL_MIN) THEN
   DO NM = 1, NMESHES
      SCARC(NM)%N_CELLS = LOCAL_INT(NM)
   ENDDO
ENDIF

END SUBROUTINE SCARC_SETUP_DIMENSIONS

!> ----------------------------------------------------------------------------------------------------
!> Setup system of equation:
!> Define matrix stencils and initialize matrices and boundary conditions on all needed levels
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_SYSTEM
INTEGER :: NM, NL

!> ---------------------------------------------------------------------------------------------
!> Setup sizes for system matrices
!> ---------------------------------------------------------------------------------------------
SELECT_METHOD: SELECT CASE (TYPE_METHOD)

   !> Global MGM method - proof of concept
   CASE (NSCARC_METHOD_MGM)

      TYPE_DISCRET = NSCARC_DISCRET_STRUCTURED                                   !> use structured discretization
      CALL SCARC_SETUP_MATRIX_SIZES (NSCARC_SIZE_MATRIX, NLEVEL_MIN)             !> setup size for fine grid

      TYPE_DISCRET = NSCARC_DISCRET_UNSTRUCTURED                                 !> use unstructured discretization
      CALL SCARC_SETUP_MATRIX_SIZES (NSCARC_SIZE_MATRIX, NLEVEL_MIN)             !> setup size for fine grid

   !> Global Krylov method
   CASE (NSCARC_METHOD_KRYLOV)

      CALL SCARC_SETUP_MATRIX_SIZES (NSCARC_SIZE_MATRIX, NLEVEL_MIN)             !> use specified discretization

      IF (HAS_TWO_LEVELS) &                                                       !> two grid levels?
         CALL SCARC_SETUP_MATRIX_SIZES (NSCARC_SIZE_MATRIX, NLEVEL_MAX)          !> setup size for coarse grid

      IF (IS_CG_GMG) THEN                                                        !> GMG-preconditioning?
         DO NL=NLEVEL_MIN+1, NLEVEL_MAX
            CALL SCARC_SETUP_MATRIX_SIZES (NSCARC_SIZE_MATRIX , NL)              !> setup size for all levels
         ENDDO
      ENDIF

   !> Global Multigrid  method
   CASE (NSCARC_METHOD_MULTIGRID)

      SELECT CASE (TYPE_MULTIGRID)                                               !> use specified discretization
         CASE (NSCARC_MULTIGRID_GEOMETRIC)                                       !> Geometric MG?
            DO NL=NLEVEL_MIN, NLEVEL_MAX
               CALL SCARC_SETUP_MATRIX_SIZES (NSCARC_SIZE_MATRIX, NL)            !> setup size for all levels
            ENDDO
         CASE (NSCARC_MULTIGRID_ALGEBRAIC)                                       !> Algebraic MG?
            CALL SCARC_SETUP_MATRIX_SIZES (NSCARC_SIZE_MATRIX, NLEVEL_MIN)       !> setup sizes for AMG
      END SELECT

END SELECT SELECT_METHOD


MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   SELECT_SOLVER: SELECT CASE (TYPE_METHOD)

      !> ---------------------------------------------------------------------------------------------
      !> McKenny-Greengard-Mayo method, solving for the structured and unstructured Poisson matrix
      !> Assemble both, the structured and unstructured Poisson matrix
      !> temporarily they will be stored separately in matrices AC and ACU due to the different
      !> settings along internal boundary cells,
      !> in the medium term, a toggle mechanism will be implemented which only switches the corresponding
      !> entries while keeping the entries which are the same for both discretization types
      !> ---------------------------------------------------------------------------------------------
      CASE (NSCARC_METHOD_MGM)

         !> First assemble unstructured matrix 
         TYPE_DISCRET = NSCARC_DISCRET_UNSTRUCTURED
         CALL SCARC_SETUP_MATRIX  (NM, NLEVEL_MIN)
         CALL SCARC_SETUP_BOUNDARY(NM, NLEVEL_MIN)

         !> Then assemble structured matrix 
         TYPE_DISCRET = NSCARC_DISCRET_STRUCTURED
         CALL SCARC_SETUP_MATRIX  (NM, NLEVEL_MIN)
         CALL SCARC_SETUP_BOUNDARY(NM, NLEVEL_MIN)


      !> ---------------------------------------------------------------------------------------------
      !> Krylov method (CG/BICG) as main solver, different preconditioners possible
      !> ---------------------------------------------------------------------------------------------
      CASE (NSCARC_METHOD_KRYLOV)

         SELECT_PRECON: SELECT CASE (TYPE_PRECON)

            !> ---------------------------------------------------------------------------------------
            !> in case of multigrid as preconditioner:
            !> ---------------------------------------------------------------------------------------
            CASE (NSCARC_RELAX_GMG)

               SELECT_PRECON_MG: SELECT CASE (TYPE_MULTIGRID)

                  !> Geometric multigrid:
                  !>    -  assemble standard n-point-matrix hierarchy on all levels (finest and all coarser)
                  CASE (NSCARC_MULTIGRID_GEOMETRIC)

                     DO NL = NLEVEL_MIN, NLEVEL_MAX
                        CALL SCARC_SETUP_MATRIX  (NM, NL)
                        CALL SCARC_SETUP_BOUNDARY(NM, NL)
#ifdef WITH_MKL
                        IF (TYPE_MKL_LEVEL(NL) == NSCARC_MKL_LOCAL) CALL SCARC_SETUP_MATRIX_MKL(NM, NL)
                        IF (TYPE_MKL_LEVEL(NL) == NSCARC_MKL_GLOBAL) CALL SCARC_SETUP_MATRIX_MKL(NM, NL)
#endif
                     ENDDO


                  !> Algebraic multigrid:
                  !>    -  use compact storage technique on all levels (no other choise possible!)
                  !>    -  assemble standard n-point-matrix only on finest level
                  !>    -  construct all coarser levels by requested coarsening strategy
                  CASE (NSCARC_MULTIGRID_ALGEBRAIC)

                     CALL SCARC_SETUP_MATRIX  (NM, NLEVEL_MIN)
                     CALL SCARC_SETUP_BOUNDARY(NM, NLEVEL_MIN)

               END SELECT SELECT_PRECON_MG

#ifdef WITH_MKL
            !> ---------------------------------------------------------------------------------------
            !> in case of LU-decomposition as preconditioner
            !> ---------------------------------------------------------------------------------------
            CASE (NSCARC_RELAX_MKL)

               SELECT CASE(TYPE_PRECON_SCOPE)

                  !> locally acting: PARDISO from MKL as preconditioners with possible coarse grid correction
                  CASE (NSCARC_SCOPE_LOCAL) 

                     CALL SCARC_SETUP_MATRIX  (NM, NLEVEL_MIN)
                     CALL SCARC_SETUP_BOUNDARY(NM, NLEVEL_MIN)
                     CALL SCARC_SETUP_MATRIX_MKL(NM, NLEVEL_MIN)
                     IF (HAS_TWO_LEVELS) THEN
                        CALL SCARC_SETUP_MATRIX  (NM, NLEVEL_MAX)
                        CALL SCARC_SETUP_BOUNDARY(NM, NLEVEL_MAX)
                        IF (TYPE_COARSE == NSCARC_COARSE_DIRECT) CALL SCARC_SETUP_MATRIX_MKL(NM, NLEVEL_MAX)
                     ENDIF

                  !> globally acting: Cluster_Sparse_Solver from MKL as preconditioner
                  CASE (NSCARC_SCOPE_GLOBAL)

                     CALL SCARC_SETUP_MATRIX  (NM, NLEVEL_MIN)
                     CALL SCARC_SETUP_BOUNDARY(NM, NLEVEL_MIN)
                     CALL SCARC_SETUP_MATRIX_MKL(NM, NLEVEL_MIN)

               END SELECT
#endif

            !> ---------------------------------------------------------------------------------------
            !> in case of one-level preconditioners (JACOBI/SSOR/FFT):
            !> assemble standard n-point-matrix on finest level with possible coarse grid correction
            !> ---------------------------------------------------------------------------------------
            CASE DEFAULT

               CALL SCARC_SETUP_MATRIX  (NM, NLEVEL_MIN)
               CALL SCARC_SETUP_BOUNDARY(NM, NLEVEL_MIN)

               IF (HAS_TWO_LEVELS) THEN
                  CALL SCARC_SETUP_MATRIX  (NM, NLEVEL_MAX)
                  CALL SCARC_SETUP_BOUNDARY(NM, NLEVEL_MAX)
#ifdef WITH_MKL
                  IF (TYPE_COARSE == NSCARC_COARSE_DIRECT) CALL SCARC_SETUP_MATRIX_MKL(NM, NLEVEL_MAX)
#endif
               ENDIF


         END SELECT SELECT_PRECON

      !> ---------------------------------------------------------------------------------------------
      !> Multigrid as main solver
      !> ---------------------------------------------------------------------------------------------
      CASE (NSCARC_METHOD_MULTIGRID)

         SELECT_MG: SELECT CASE (TYPE_MULTIGRID)

            !> ---------------------------------------------------------------------------------------
            !> Geometric multigrid:
            !>    -  assemble standard n-point-matrix hierarchy on all levels
            !>    -  use MKL coarse grid solver if requested
            !> ---------------------------------------------------------------------------------------
            CASE (NSCARC_MULTIGRID_GEOMETRIC)

               DO NL = NLEVEL_MIN, NLEVEL_MAX
                  CALL SCARC_SETUP_MATRIX  (NM, NL)
                  CALL SCARC_SETUP_BOUNDARY(NM, NL)
               ENDDO

#ifdef WITH_MKL
               DO NL = NLEVEL_MIN, NLEVEL_MAX-1
                  IF (TYPE_MKL_LEVEL(NL) == NSCARC_MKL_LOCAL   .OR. &
                      TYPE_MKL_LEVEL(NL) == NSCARC_MKL_GLOBAL) THEN
                     CALL SCARC_SETUP_MATRIX_MKL(NM, NL)
                  ENDIF
               ENDDO

               IF (TYPE_COARSE == NSCARC_COARSE_DIRECT) CALL SCARC_SETUP_MATRIX_MKL(NM, NLEVEL_MAX)
#endif

            !> ---------------------------------------------------------------------------------------
            !> Algebraic multigrid:
            !>    -  use compact storage technique (no other choice possible!)
            !>    -  assemble standard n-point-matrix only on finest level
            !>    -  construct all coarser levels later by requested coarsening strategy
            !> ---------------------------------------------------------------------------------------
            CASE (NSCARC_MULTIGRID_ALGEBRAIC)

               CALL SCARC_SETUP_MATRIX  (NM, NLEVEL_MIN)
               CALL SCARC_SETUP_BOUNDARY(NM, NLEVEL_MIN)

         END SELECT SELECT_MG



#ifdef WITH_MKL
      !> ---------------------------------------------------------------------------------------------
      !> MKL-Pardiso method
      !> ---------------------------------------------------------------------------------------------
      CASE (NSCARC_METHOD_LU)

         CALL SCARC_SETUP_MATRIX(NM, NLEVEL_MIN)
         CALL SCARC_SETUP_BOUNDARY(NM, NLEVEL_MIN)
         CALL SCARC_SETUP_MATRIX_MKL(NM, NLEVEL_MIN)
#endif

   END SELECT SELECT_SOLVER

ENDDO MESHES_LOOP

!> -------------------------------------------------------------------------------------------
!> Debug matrix and wall structures - only if directive SCARC_DEBUG is set
!> -------------------------------------------------------------------------------------------
#ifdef WITH_SCARC_DEBUG
IF (TYPE_METHOD == NSCARC_METHOD_MGM) THEN
   TYPE_DISCRET = NSCARC_DISCRET_UNSTRUCTURED
   CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_MATRIX , NL, 'SYSTEM-MATRIX-UNSTRUCTURED')
   TYPE_DISCRET = NSCARC_DISCRET_STRUCTURED
   CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_MATRIX , NL, 'SYSTEM-MATRIX-STRUCTURED')
ELSE
   DO NL = NLEVEL_MIN, NLEVEL_MAX
      CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_MATRIX , NL, 'SYSTEM-MATRIX')
   ENDDO
ENDIF
IF (TYPE_PRECON == NSCARC_RELAX_MKL) &
CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_MATRIX_SYM , NLEVEL_MIN, 'SYSTEM-MATRIX-SYMM')
CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_WALL, NLEVEL_MIN, 'SYSTEM-WALL')
#endif

END SUBROUTINE SCARC_SETUP_SYSTEM

!> ------------------------------------------------------------------------------------------------
!> Allocate matrix for the usual 5-point-stencil (2D) or 7-point-stencil (3D)
!> Use compact storage technique:
!>
!> Compression technique to store sparse matrices, non-zero entries are stored
!> in a 1D-vector B(.), row after row,
!> Each row starts with its diagonal entry followed by the other non-zero entries
!> In order to identify each element, pointer arrays ROW and COL are needed,
!> ROW points to the several diagonal entries in vector B(.),
!> COL points to the columns which non-zero entries in the matrix stencil
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_MATRIX (NM, NL)
INTEGER, INTENT(IN) :: NM, NL
INTEGER :: IX, IY, IZ, IC, IP
TYPE (SCARC_LEVEL_TYPE), POINTER :: L=>NULL()
TYPE (SCARC_DISCRET_TYPE), POINTER :: D=>NULL()
TYPE (SCARC_MATRIX_COMPACT_TYPE), POINTER :: AC=>NULL()
TYPE (SCARC_MATRIX_BANDED_TYPE),  POINTER :: AB=>NULL()

!> Initialize data structures
L => SCARC(NM)%LEVEL(NL)

!> Compute single matrix entries and corresponding row and column pointers
!> Along internal boundaries use placeholders for the neighboring matrix entries
!> which will be communicated in a following step
SELECT_STORAGE_TYPE: SELECT CASE (TYPE_MATRIX)

   !>
   !> ---------- COMPACT Storage technique
   !>
   CASE (NSCARC_MATRIX_COMPACT)

      SELECT CASE(TYPE_DISCRET)
         CASE (NSCARC_DISCRET_STRUCTURED)
            D => L%SD
         CASE (NSCARC_DISCRET_UNSTRUCTURED)
            D => L%UD
      END SELECT
      AC => D%AC
      CALL SCARC_ALLOCATE_MATRIX_COMPACT (AC, 'AC', NL, NSCARC_INIT_ZERO)

WRITE(MSG%LU_DEBUG,*) 'TYPE_DISCRET =',TYPE_DISCRET
      IP  = 1
      DO IZ = 1, L%NZ
         DO IY = 1, L%NY
            DO IX = 1, L%NX
   
               IF (TYPE_DISCRET == NSCARC_DISCRET_UNSTRUCTURED .AND. &
                   L%CELL_STATE(IX, IY, IZ) /= NSCARC_CELL_GASPHASE) CYCLE

               !> Main diagonal entry
               IC = D%CELL_NUMBER(IX, IY, IZ)
               CALL SCARC_SETUP_MATRIX_MAINDIAG_COMPACT (IC, IX, IY, IZ, IP, NM, NL)

               !> Subdiagonals in x-direction
               IF (VALID_SUBDIAG(IX, IY, IZ,  1, NM, NL)) &
                  CALL SCARC_SETUP_MATRIX_SUBDIAG_COMPACT(IC, IX, IY, IZ, IX-1, IY, IZ, IP,  1, NM, NL)
               IF (VALID_SUBDIAG(IX, IY, IZ, -1, NM, NL)) &
                  CALL SCARC_SETUP_MATRIX_SUBDIAG_COMPACT(IC, IX, IY, IZ, IX+1, IY, IZ, IP, -1, NM, NL)
          
               !> Subdiagonals in y-direction
               IF (VALID_SUBDIAG(IX, IY, IZ,  2, NM, NL)) &
                  CALL SCARC_SETUP_MATRIX_SUBDIAG_COMPACT(IC, IX, IY, IZ, IX, IY-1, IZ, IP,  2, NM, NL)
               IF (VALID_SUBDIAG(IX, IY, IZ, -2, NM, NL)) &
                  CALL SCARC_SETUP_MATRIX_SUBDIAG_COMPACT(IC, IX, IY, IZ, IX, IY+1, IZ, IP, -2, NM, NL)

               !> Subdiagonals in z-direction
               IF (VALID_SUBDIAG(IX, IY, IZ,  3, NM, NL)) &
                  CALL SCARC_SETUP_MATRIX_SUBDIAG_COMPACT(IC, IX, IY, IZ, IX, IY, IZ-1, IP,  3, NM, NL)
               IF (VALID_SUBDIAG(IX, IY, IZ, -3, NM, NL)) &
                  CALL SCARC_SETUP_MATRIX_SUBDIAG_COMPACT(IC, IX, IY, IZ, IX, IY, IZ+1, IP, -3, NM, NL)

            ENDDO
         ENDDO
      ENDDO
      
      AC%ROW(AC%NR) = IP
      AC%NA         = IP-1                         !> set correct number of matrix entries

      WRITE(MSG%LU_DEBUG,*) 'SYSTEM: AC%NA=',AC%NA
      WRITE(MSG%LU_DEBUG,*) 'SYSTEM: AC%NR=',AC%NR

      CALL SCARC_RESIZE_MATRIX_COMPACT(AC, AC%NA, 'Resized System-Matrix')


   !> 
   !> ---------- BANDED Storage technique
   !>
   CASE (NSCARC_MATRIX_BANDED)

      SELECT CASE(TYPE_DISCRET)
         CASE (NSCARC_DISCRET_STRUCTURED)
            D => L%SD
         CASE (NSCARC_DISCRET_UNSTRUCTURED)
            D => L%UD
      END SELECT
      AB => D%AB
      CALL SCARC_ALLOCATE_MATRIX_BANDED(AB, 'AB', NL, NSCARC_INIT_ZERO)

      IP  = 1
      DO IZ = 1, L%NZ
         DO IY = 1, L%NY
            DO IX = 1, L%NX
   
               IF (TYPE_DISCRET == NSCARC_DISCRET_UNSTRUCTURED .AND. &
                   L%CELL_STATE(IX, IY, IZ) /= NSCARC_CELL_GASPHASE) CYCLE

               IC = D%CELL_NUMBER(IX, IY, IZ)

               !> Lower diagonal entries
               IF (VALID_SUBDIAG(IX, IY, IZ,  3, NM, NL)) &
                  CALL SCARC_SETUP_MATRIX_SUBDIAG_BANDED(IC, IX, IY, IZ, IX, IY, IZ-1, 3, NM, NL)
               IF (VALID_SUBDIAG(IX, IY, IZ,  2, NM, NL)) &
                  CALL SCARC_SETUP_MATRIX_SUBDIAG_BANDED(IC, IX, IY, IZ, IX, IY-1, IZ, 2, NM, NL)
               IF (VALID_SUBDIAG(IX, IY, IZ,  1, NM, NL)) &
                  CALL SCARC_SETUP_MATRIX_SUBDIAG_BANDED(IC, IX, IY, IZ, IX-1, IY, IZ, 1, NM, NL)

               !> Main diagonal entry
               CALL SCARC_SETUP_MATRIX_MAINDIAG_BANDED (IC, IX, IY, IZ, NM, NL)

               !> Upper diagonal entries
               IF (VALID_SUBDIAG(IX, IY, IZ, -1, NM, NL)) &
                  CALL SCARC_SETUP_MATRIX_SUBDIAG_BANDED(IC, IX, IY, IZ, IX+1, IY, IZ, -1, NM, NL)
               IF (VALID_SUBDIAG(IX, IY, IZ, -2, NM, NL)) &
                  CALL SCARC_SETUP_MATRIX_SUBDIAG_BANDED(IC, IX, IY, IZ, IX, IY+1, IZ, -2, NM, NL)
               IF (VALID_SUBDIAG(IX, IY, IZ, -3, NM, NL)) &
                  CALL SCARC_SETUP_MATRIX_SUBDIAG_BANDED(IC, IX, IY, IZ, IX, IY, IZ+1, -3, NM, NL)

            ENDDO
         ENDDO
      ENDDO
      
END SELECT SELECT_STORAGE_TYPE

END SUBROUTINE SCARC_SETUP_MATRIX

!> ------------------------------------------------------------------------------------------------
!> Check if a subdiagonal entry must be computed in direction IOR0
!> ------------------------------------------------------------------------------------------------
LOGICAL FUNCTION VALID_SUBDIAG(IX, IY, IZ, IOR0, NM, NL)
INTEGER, INTENT(IN)  :: IX, IY, IZ, IOR0, NM, NL
INTEGER :: IC_INDEX, IW_INDEX
TYPE (SCARC_LEVEL_TYPE), POINTER :: L=>NULL()
TYPE (SCARC_DISCRET_TYPE), POINTER :: D=>NULL()

VALID_SUBDIAG = .FALSE.
IF (TWO_D .AND. ABS(IOR0) == 2) RETURN

L => SCARC(NM)%LEVEL(NL)
SELECT CASE(TYPE_DISCRET)
   CASE (NSCARC_DISCRET_STRUCTURED)
      D => L%SD
   CASE (NSCARC_DISCRET_UNSTRUCTURED)
      D => L%UD
END SELECT

SELECT CASE (TYPE_DISCRET)

   !> If a structured discretization is used, then subdiagonals are built in every direction
   CASE (NSCARC_DISCRET_STRUCTURED)

      VALID_SUBDIAG = .TRUE.
      RETURN

   !> Else check the type of the neighboring cell in direction IOR0
   CASE (NSCARC_DISCRET_UNSTRUCTURED)

      !> get cell index of corresponding cell and check its wall index 
      IC_INDEX = L%CELL_INDEX_PTR(IX, IY, IZ)
      IW_INDEX = 0
      IF (IC_INDEX /= 0) IW_INDEX  = L%WALL_INDEX_PTR(IC_INDEX, -IOR0)

      !> If this wall index is zero, build a subdiagonal in this direction
      IF (IW_INDEX == 0) THEN
         VALID_SUBDIAG = .TRUE.
         RETURN

      !> Else build the subdiagonal entry only if the wall cell is of interpolated type
      ELSE
         IF (D%WALL(IW_INDEX)%BOUNDARY_TYPE== INTERPOLATED_BOUNDARY) THEN
            VALID_SUBDIAG = .TRUE.
            RETURN
         ENDIF
      ENDIF

END SELECT
RETURN
   
END FUNCTION VALID_SUBDIAG

!> ------------------------------------------------------------------------------------------------
!> Set main diagonal entry for matrix - full matrix of the global problem
!> In case of an equidistant grid, we get the usual 5-point (2d) and 7-point (3d) stencil
!> If two meshes with different step sizes meet, we get a weighted stencil along internal wall cells
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_MATRIX_MAINDIAG_COMPACT (IC, IX, IY, IZ, IP, NM, NL)
INTEGER, INTENT(IN) :: IC, IX, IY, IZ, NM, NL
INTEGER, INTENT(INOUT) :: IP
TYPE (SCARC_LEVEL_TYPE), POINTER :: L=>NULL()
TYPE (SCARC_DISCRET_TYPE), POINTER :: D=>NULL()
TYPE (SCARC_MATRIX_COMPACT_TYPE), POINTER :: AC=>NULL()

L  => SCARC(NM)%LEVEL(NL)
SELECT CASE(TYPE_DISCRET)
   CASE (NSCARC_DISCRET_STRUCTURED)
      AC => L%SD%AC
   CASE (NSCARC_DISCRET_UNSTRUCTURED)
      AC => L%UD%AC
END SELECT

WRITE(MSG%LU_DEBUG,*) 'MATRIX_MAINDIAG_COMPACT: IC, IX, IY, IZ:', IC, IX, IY, IZ

AC%VAL(IP) = - 2.0_EB/(L%DXL(IX-1)*L%DXL(IX))
IF (.NOT.TWO_D) AC%VAL(IP) = AC%VAL(IP) - 2.0_EB/(L%DYL(IY-1)*L%DYL(IY))
AC%VAL(IP) = AC%VAL(IP) - 2.0_EB/(L%DZL(IZ-1)*L%DZL(IZ))

AC%ROW(IC) = IP
AC%COL(IP) = IC

#ifdef WITH_MKL
IF (IS_MKL_LEVEL(NL)) AC%COL_GLOBAL(IP) = AC%COL(IP) + D%NC_OFFSET(NM)
#endif

AC%STENCIL(0) = AC%VAL(IP)

IP = IP + 1
END SUBROUTINE SCARC_SETUP_MATRIX_MAINDIAG_COMPACT


!> ------------------------------------------------------------------------------------------------
!> Compute subdiagonal contribution in direction IOR0 in case of a COMPACT matrix
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_MATRIX_SUBDIAG_COMPACT (IC, IX1, IY1, IZ1, IX2, IY2, IZ2, IP, IOR0, NM, NL)
INTEGER, INTENT(IN) :: IC, IX1, IY1, IZ1, IX2, IY2, IZ2, IOR0, NM, NL
INTEGER, INTENT(INOUT) :: IP
INTEGER :: IW, IPTR
#ifdef WITH_MKL
INTEGER :: IX, IY, IZ
#endif
REAL(EB) :: DSCAL, DH1, DH2
LOGICAL  :: BINTERNAL
TYPE (SCARC_LEVEL_TYPE), POINTER :: L=>NULL()
TYPE (SCARC_DISCRET_TYPE), POINTER :: D=>NULL()
TYPE (SCARC_MATRIX_COMPACT_TYPE), POINTER :: AC=>NULL()

L => SCARC(NM)%LEVEL(NL)
SELECT CASE(TYPE_DISCRET)
   CASE (NSCARC_DISCRET_STRUCTURED)
      D => L%SD
   CASE (NSCARC_DISCRET_UNSTRUCTURED)
      D => L%UD
END SELECT
AC => D%AC

SELECT CASE (ABS(IOR0))
   CASE ( 1)
      IPTR = IX1
   CASE ( 2)
      IPTR = IY1
   casE ( 3)
      IPTR = IZ1
END SELECT

DH1 = L%FACE(IOR0)%DH(IPTR-1)
DH2 = L%FACE(IOR0)%DH(IPTR)

!> set bounds and step sizes depending on orientation of face (lower or upper subdiagonal)
IF (IOR0 > 0) THEN
   BINTERNAL = IPTR > 1
   DSCAL = 2.0_EB/(DH1*(DH1+DH2))
ELSE
   BINTERNAL = IPTR < L%FACE(IOR0)%NFC
   DSCAL = 2.0_EB/(DH2*(DH1+DH2))
ENDIF

!> if IC is an internal cell of the mesh, compute usual matrix contribution for corresponding subdiagonal
IF (BINTERNAL) THEN

   IF (TYPE_DISCRET == NSCARC_DISCRET_STRUCTURED .OR. L%CELL_STATE(IX2, IY2, IZ2) == NSCARC_CELL_GASPHASE) THEN
WRITE(MSG%LU_DEBUG,*) 'MATRIX_SUBDIAG_COMPACT: IC, IP, IPTR, IOR0:', IC, IP, IPTR, IOR0
      AC%VAL(IP) = AC%VAL(IP) + DSCAL
      AC%COL(IP) = D%CELL_NUMBER(IX2, IY2, IZ2)

      AC%STENCIL(-IOR0) = AC%VAL(IP)

#ifdef WITH_MKL
      IF (IS_MKL_LEVEL(NL)) AC%COL_GLOBAL(IP)= AC%COL(IP) + D%NC_OFFSET(NM)
#endif
      IP = IP + 1
   ELSE
     CALL SCARC_SHUTDOWN(NSCARC_ERROR_MATRIX_SUBDIAG, SCARC_NONE, NSCARC_NONE)
   ENDIF

!> if IC is a boundary cell of the mesh, compute matrix contribution only if there is a neighbor for that cell
ELSE IF (L%FACE(IOR0)%N_NEIGHBORS /= 0) THEN

   IW = ASSIGN_SUBDIAG_TYPE (IC, IOR0, NM, NL)         !> get IW of a possibly suitable neighbor at face IOR0
   IF (IW /= 0) then                                   !> if available, build corresponding subdiagonal entry

      AC%VAL(IP) = AC%VAL(IP) + DSCAL
      AC%COL(IP) = D%WALL(IW)%ICE(1)                   !> store its extended number in matrix column pointers
WRITE(MSG%LU_DEBUG,*) 'MATRIX_SUBDIAG: IP, IW:', IP, IW, AC%COL(IP)

      AC%STENCIL(-IOR0) = AC%VAL(IP)

#ifdef WITH_MKL
      IF (IS_MKL_LEVEL(NL)) THEN                         !> if MKL method used, also store its global number
         IX = D%WALL(IW)%IXG
         IY = D%WALL(IW)%IYG
         IZ = D%WALL(IW)%IZG
         AC%COL_GLOBAL(IP) = D%CELL_NUMBER(IX, IY, IZ) + D%NC_OFFSET(D%WALL(IW)%NOM)
      ENDIF
#endif

      IP = IP + 1
   ENDIF

ENDIF

END SUBROUTINE SCARC_SETUP_MATRIX_SUBDIAG_COMPACT

!> ------------------------------------------------------------------------------------------------
!> Determine if cell IC has a neighbor and, if yes, return corresponding IW-value
!> ------------------------------------------------------------------------------------------------
INTEGER FUNCTION ASSIGN_SUBDIAG_TYPE (IC, IOR0, NM, NL) 
INTEGER, INTENT(IN) :: IC, IOR0, NM, NL
INTEGER :: IXW, IYW, IZW
INTEGER :: IXG, IYG, IZG
INTEGER :: IW
TYPE (SCARC_LEVEL_TYPE), POINTER :: L=>NULL()
TYPE (SCARC_DISCRET_TYPE), POINTER :: D=>NULL()

L => SCARC(NM)%LEVEL(NL)
SELECT CASE(TYPE_DISCRET)
   CASE (NSCARC_DISCRET_STRUCTURED)
      D => L%SD
   CASE (NSCARC_DISCRET_UNSTRUCTURED)
      D => L%UD
END SELECT

ASSIGN_SUBDIAG_TYPE = -1
SEARCH_WALL_CELLS_LOOP: DO IW = L%FACE(IOR0)%IWG_PTR, L%FACE(IOR0)%IWG_PTR+L%FACE(IOR0)%NFW

  IF (D%WALL(IW)%NOM == 0) CYCLE

  IXW = D%WALL(IW)%IXW
  IYW = D%WALL(IW)%IYW
  IZW = D%WALL(IW)%IZW

  IF (D%CELL_NUMBER(IXW, IYW, IZW) /= IC) CYCLE

  IXG = D%WALL(IW)%IXG
  IYG = D%WALL(IW)%IYG
  IZG = D%WALL(IW)%IZG

  IF (L%CELL_STATE(IXG, IYG, IZG) /= NSCARC_CELL_SOLID) THEN
     ASSIGN_SUBDIAG_TYPE = IW
     RETURN
  ENDIF

ENDDO SEARCH_WALL_CELLS_LOOP

END FUNCTION ASSIGN_SUBDIAG_TYPE

!> ------------------------------------------------------------------------------------------------
!> Set main diagonal entry for matrix - full matrix of the global problem
!> In case of an equidistant grid, we get the usual 5-point (2d) and 7-point (3d) stencil
!> If two meshes with different step sizes meet, we get a weighted stencil along internal wall cells
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_MATRIX_MAINDIAG_BANDED (IC, IX, IY, IZ, NM, NL)
INTEGER, INTENT(IN)  :: IC, IX, IY, IZ, NM, NL
INTEGER :: ID
TYPE (SCARC_LEVEL_TYPE), POINTER :: L=>NULL()
TYPE (SCARC_DISCRET_TYPE), POINTER :: D=>NULL()
TYPE (SCARC_MATRIX_BANDED_TYPE), POINTER :: AB=>NULL()

L  => SCARC(NM)%LEVEL(NL)
SELECT CASE(TYPE_DISCRET)
   CASE (NSCARC_DISCRET_STRUCTURED)
      D => L%SD
   CASE (NSCARC_DISCRET_UNSTRUCTURED)
      D => L%UD
END SELECT
AB => D%AB

ID = AB%POS(0)               !> get column vector corresponding to matrix diagonal

AB%VAL(ID, IC) = AB%VAL(ID, IC) - 2.0_EB/(L%DXL(IX-1)*L%DXL(IX))
IF (.NOT.TWO_D)  AB%VAL(ID, IC) = AB%VAL(ID, IC) - 2.0_EB/(L%DYL(IY-1)*L%DYL(IY))
AB%VAL(ID, IC) = AB%VAL(ID, IC) - 2.0_EB/(L%DZL(IZ-1)*L%DZL(IZ))

AB%STENCIL(0) = AB%VAL(ID, IC)

END SUBROUTINE SCARC_SETUP_MATRIX_MAINDIAG_BANDED

!> ------------------------------------------------------------------------------------------------
!> Compute subdiagonal contribution in direction IOR0 in case of a COMPACT matrix
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_MATRIX_SUBDIAG_BANDED (IC, IX1, IY1, IZ1, IX2, IY2, IZ2, IOR0, NM, NL)
INTEGER, INTENT(IN) :: IC, IX1, IY1, IZ1, IX2, IY2, IZ2, IOR0, NM, NL
INTEGER :: IW, IPTR, ID
REAL(EB) :: DSCAL, DH1, DH2
LOGICAL  :: BINTERNAL
TYPE (SCARC_LEVEL_TYPE), POINTER :: L=>NULL()
TYPE (SCARC_DISCRET_TYPE), POINTER :: D=>NULL()
TYPE (SCARC_MATRIX_BANDED_TYPE), POINTER :: AB=>NULL()

L  => SCARC(NM)%LEVEL(NL)
SELECT CASE(TYPE_DISCRET)
   CASE (NSCARC_DISCRET_STRUCTURED)
      D => L%SD
   CASE (NSCARC_DISCRET_UNSTRUCTURED)
      D => L%UD
END SELECT
AB => D%AB

SELECT CASE (ABS(IOR0))
   CASE ( 1)
      IPTR = IX1
   CASE ( 2)
      IPTR = IY1
   CASE ( 3)
      IPTR = IZ1
END SELECT

DH1 = L%FACE(IOR0)%DH(IPTR-1)
DH2 = L%FACE(IOR0)%DH(IPTR)

!> set bounds and step sizes depending on orientation of face (lower or upper subdiagonal)
IF (IOR0 > 0) THEN
   BINTERNAL = IPTR > 1
   DSCAL = 2.0_EB/(DH1*(DH1+DH2))
ELSE
   BINTERNAL = IPTR < L%FACE(IOR0)%NFC
   DSCAL = 2.0_EB/(DH2*(DH1+DH2))
ENDIF

ID  = AB%POS(IOR0)                                !> get column vector corresponding to diagonal of IOR0

!> if IC is an internal cell of the mesh, compute usual matrix contribution for corresponding subdiagonal
IF (BINTERNAL) THEN

   IF (TYPE_DISCRET == NSCARC_DISCRET_STRUCTURED .OR. L%CELL_STATE(IX2, IY2, IZ2) == NSCARC_CELL_GASPHASE) THEN
      AB%VAL(ID, IC)   = AB%VAL(ID, IC) + DSCAL
      AB%STENCIL(IOR0) = AB%VAL(ID, IC)
   ELSE
     CALL SCARC_SHUTDOWN(NSCARC_ERROR_MATRIX_SUBDIAG, SCARC_NONE, NSCARC_NONE)
   ENDIF

!> if IC is a boundary cell of the mesh, compute matrix contribution only if there is a neighbor for that cell
ELSE IF (L%FACE(IOR0)%N_NEIGHBORS /= 0) THEN

   IW = ASSIGN_SUBDIAG_TYPE (IC, IOR0, NM, NL)         !> get IW of a possibly suitable neighbor at face IOR0
   IF (IW /= 0) THEN
      AB%VAL(ID, IC)   = AB%VAL(ID, IC) + DSCAL
      AB%STENCIL(IOR0) = AB%VAL(ID, IC)
   ENDIF
   
ENDIF

END SUBROUTINE SCARC_SETUP_MATRIX_SUBDIAG_BANDED

#ifdef WITH_MKL
!> ------------------------------------------------------------------------------------------------
!> Build system matrix for MKL solver
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_MATRIX_MKL (NM, NL)
INTEGER, INTENT(IN) :: NM, NL
INTEGER :: IC, JC, JC0
INTEGER :: ICS, JCS
INTEGER :: ICOL, JCOL, IAS
INTEGER :: ISYM, JSYM, NSYM
REAL(EB) :: VAL = 0.0_EB, VALS = 0.0_EB, DIFF
LOGICAL  :: BSYM
INTEGER, DIMENSION(:), ALLOCATABLE :: KCOL_AUX, KC_AUX
TYPE (SCARC_LEVEL_TYPE), POINTER :: L=>NULL()
TYPE (SCARC_DISCRET_TYPE), POINTER :: D=>NULL()
TYPE (SCARC_MATRIX_COMPACT_TYPE), POINTER :: AC=>NULL(), AC_SYM=>NULL()

L   => SCARC(NM)%LEVEL(NL)
SELECT CASE(TYPE_DISCRET)
   CASE (NSCARC_DISCRET_STRUCTURED)
      D => L%SD
   CASE (NSCARC_DISCRET_UNSTRUCTURED)
      D => L%UD
END SELECT

AC     => D%AC
AC_SYM => D%AC_SYM

!> ------------------------------------------------------------------------------------------------
!> Store only symmetric parts of matrix (diagonal and upper part)
!> ------------------------------------------------------------------------------------------------
IF (SCARC_MKL_MTYPE == 'SYMMETRIC') THEN

   !> First check whether symmetry of system matrix is guaranteed
   DO IC = 1, D%NC

      COLUMN_LOOP: DO ICOL = AC%ROW(IC)+1, AC%ROW(IC+1)-1
         ICS = AC%COL(ICOL)
         VAL = AC%VAL(ICOL)
         IF (ICS > IC .AND. ICS <= D%NC) THEN
            BSYM = .FALSE.
            DO JCOL = AC%ROW(ICS)+1, AC%ROW(ICS+1)-1
               JCS = AC%COL(JCOL)
               IF (JCS == IC) THEN
                  VALS = AC%VAL(JCOL)
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

   !> Compute number of entries in symmetric matrix
   AC_SYM%NA = 0
   DO IC = 1, D%NC
      DO ICOL = AC%ROW(IC), AC%ROW(IC+1)-1
         IF (TYPE_MKL_LEVEL(NL) == NSCARC_MKL_LOCAL) THEN
            JC = AC%COL(ICOL)
            IF (JC >= IC .AND. JC <= D%NC) AC_SYM%NA = AC_SYM%NA+1  
         ELSE IF (TYPE_MKL_LEVEL(NL) == NSCARC_MKL_GLOBAL) THEN
            JC = AC%COL_GLOBAL(ICOL)
            IF (JC >= IC + D%NC_OFFSET(NM)) AC_SYM%NA = AC_SYM%NA+1
         ELSE
            CALL SCARC_SHUTDOWN(NSCARC_ERROR_MATRIX_SETUP, SCARC_NONE, TYPE_MKL_LEVEL(NL))
         ENDIF
      ENDDO
   ENDDO

ELSE
   AC_SYM%NA = AC%NA
ENDIF

!> ------------------------------------------------------------------------------------------------
!> allocate storage for symmetric matrix and its column and row pointers
!> ------------------------------------------------------------------------------------------------
AC_SYM%NR = AC%NR
CALL SCARC_ALLOCATE_MATRIX_COMPACT(AC_SYM, 'AC_SYM', NL, NSCARC_INIT_ZERO)

!> if global MKL method is used, also allocate auxiliary space for computation of global numbering
IF (IS_MKL_LEVEL(NL)) THEN
   CALL SCARC_ALLOCATE_INT1(KCOL_AUX, 1, AC%NSTENCIL, NSCARC_INIT_NONE, 'KCOL_AUX')
   CALL SCARC_ALLOCATE_INT1(KC_AUX  , 1, AC%NSTENCIL, NSCARC_INIT_NONE, 'KC_AUX')
ENDIF

!> ------------------------------------------------------------------------------------------------
!> extract symmetric matrix part from usual system matrix
!> ------------------------------------------------------------------------------------------------
IAS = 1
DO IC = 1, D%NC
   AC_SYM%ROW(IC) = IAS

   !> blockwise use of local MKL solvers - no global numbering required
   IF (TYPE_MKL_LEVEL(NL) == NSCARC_MKL_LOCAL) THEN    

      DO ICOL = AC%ROW(IC), AC%ROW(IC+1)-1
         JC = AC%COL(ICOL)

            IF (JC >= IC .AND. JC <= D%NC) THEN
               AC_SYM%COL(IAS) = AC%COL(ICOL)
#ifdef WITH_MKL_FB
               AC_SYM%VAL(IAS) = REAL(AC%VAL(ICOL),FB)
#else
               AC_SYM%VAL(IAS) = AC%VAL(ICOL)
#endif
               IAS = IAS + 1
            ENDIF
      ENDDO

   !> global use of MKL solver - get global numbering of matrix elements
   ELSE IF (TYPE_MKL_LEVEL(NL) == NSCARC_MKL_GLOBAL) THEN

      !> store indices of all diagonal and upper-diagonal entries
      KCOL_AUX = 0
      KC_AUX   = 99999999
      ISYM = 1
      JC0 = AC%COL_GLOBAL(AC%ROW(IC))
      DO ICOL = AC%ROW(IC), AC%ROW(IC+1)-1
         !JC = AC%COL(ICOL)
         JC = AC%COL_GLOBAL(ICOL)
         IF (SCARC_MKL_MTYPE == 'SYMMETRIC') THEN
            IF (JC >= JC0) THEN
               KCOL_AUX(ISYM) = ICOL
               KC_AUX(ISYM) = JC
               ISYM  = ISYM  + 1
            ENDIF
         ELSE
            KCOL_AUX(ISYM) = ICOL
            KC_AUX(ISYM) = JC
            ISYM  = ISYM  + 1
         ENDIF
      ENDDO
      NSYM = ISYM - 1

      !> sort them in increasing order (for the use of DSS-PARDISO functionality)
      JSYM = 1
      SORT_LOOP: DO WHILE (JSYM <= NSYM)
         DO ISYM = 1, NSYM
            JC = KC_AUX(ISYM)
            IF (JC == 99999999) CYCLE
            IF (JC <= MINVAL(KC_AUX)) THEN
               ICOL = KCOL_AUX(ISYM)
#ifdef WITH_MKL_FB
               AC_SYM%VAL(IAS) = REAL(AC%VAL(ICOL),FB)
#else
               AC_SYM%VAL(IAS) = AC%VAL(ICOL)
#endif
               !AC_SYM%COL(IAS) = AC%COL(ICOL)
               AC_SYM%COL(IAS) = AC%COL_GLOBAL(ICOL)
               KC_AUX(ISYM) = 99999999            ! mark entry as already used
               IAS  = IAS  + 1
            ENDIF
         ENDDO
         JSYM = JSYM + 1
      ENDDO SORT_LOOP
   ENDIF
ENDDO

AC_SYM%ROW(D%NC+1) = IAS

IF (IS_MKL_LEVEL(NL)) THEN
   DEALLOCATE (KCOL_AUX, STAT=IERROR)
   DEALLOCATE (KC_AUX,   STAT=IERROR)
   DEALLOCATE (AC%COL_GLOBAL, STAT=IERROR)
ENDIF

END SUBROUTINE SCARC_SETUP_MATRIX_MKL
#endif


!> ------------------------------------------------------------------------------------------------
!> Set pointer for different structures on level NL
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_BOUNDARY (NM, NL)
INTEGER, INTENT(IN) :: NM, NL
INTEGER :: I, J, K, IOR0, IW, IC, NOM, IP, NW, NC, ICN, ICE, JC, ICOL, IS=0
REAL(EB) :: DBC
TYPE (SCARC_LEVEL_TYPE), POINTER :: L=>NULL()
TYPE (SCARC_DISCRET_TYPE), POINTER :: D=>NULL()
TYPE (SCARC_MATRIX_COMPACT_TYPE), POINTER :: AC=>NULL()
TYPE (SCARC_MATRIX_BANDED_TYPE),  POINTER :: AB=>NULL()

L  => SCARC(NM)%LEVEL(NL)

SELECT CASE(TYPE_DISCRET)
   CASE (NSCARC_DISCRET_STRUCTURED)
      D => L%SD
      NW = L%N_WALL_CELLS_EXT
   CASE (NSCARC_DISCRET_UNSTRUCTURED)
      D => L%UD
      NW = L%N_WALL_CELLS_EXT + L%N_WALL_CELLS_INT
END SELECT

MATRIX_TYPE_SELECT: SELECT CASE (TYPE_MATRIX)
   
   !>
   !> ---------- Matrix in COMPACT storage technique
   !>
   CASE (NSCARC_MATRIX_COMPACT)

      !>
      !> If A is a pure Neumann matrix, get neighboring cell indices of communicated stencil legs for condensed system
      !> also save values and column indices of last matrix row of last mesh
      !>
      AC => D%AC
      NO_DIRIC_COMPACT_IF: IF (N_DIRIC_GLOBAL(NLEVEL_MIN) == 0) THEN
      
         CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_CELL_INDEX  , NLEVEL_MIN)
         CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_MATRIX_VALUE, NLEVEL_MIN)
      
         LAST_CELL_IN_LAST_MESH_COMPACT_IF: IF (NM == NMESHES) THEN
      
            NC = D%NC_LOCAL(NMESHES)
            IP = AC%ROW(NC)
      
            !> store column indices and values of diagonal and all off-diagonal entries in last row
            !> index '1' correspondings to main diagonal entry
            IS = IS + 1
            ICOL = 1
            AC%STORE(1)%PTR(ICOL)      = IP
            AC%STORE(1)%COL(ICOL)      = AC%COL(IP)
            AC%STORE(1)%VAL_ORIG(ICOL) = AC%VAL(IP)
            AC%STORE(1)%VAL_COND(ICOL) = 1.0_EB
      
            DO IP = AC%ROW(NC)+1, AC%ROW(NC+1)-1
               ICOL = ICOL + 1
               AC%STORE(1)%PTR(ICOL)      = IP
               AC%STORE(1)%COL(ICOL)      = AC%COL(IP)
               AC%STORE(1)%VAL_ORIG(ICOL) = AC%VAL(IP)
               AC%STORE(1)%VAL_COND(ICOL) = 0.0_EB
            ENDDO
      
            AC%STORE(1)%NROW = NC              !> row index of last row
            AC%STORE(1)%NCOL = ICOL            !> number of stored columns
      
            !> within last mesh: check which other cells have a connection to the last cell;
            !> in each corresponding matrix row store the column index and value of just that matrix entry
            !> for each direction only one value has to be stored
            JC = NC - 1
            DO IP = AC%ROW(JC)+1, AC%ROW(JC+1)-1
               IF (AC%COL(IP) == NC) THEN
                  IS = IS + 1
                  AC%STORE(IS)%PTR(1)      = IP
                  AC%STORE(IS)%COL(1)      = JC
                  AC%STORE(IS)%VAL_ORIG(1) = AC%VAL(IP)
                  AC%STORE(IS)%VAL_COND(1) = 0.0_EB
                  AC%STORE(IS)%NROW        = JC
                  AC%STORE(IS)%NCOL        = 1
                  EXIT
               ENDIF
            ENDDO
      
            JC = NC - L%NX
            DO IP = AC%ROW(JC)+1, AC%ROW(JC+1)-1
               IF (AC%COL(IP) == NC) THEN
                  IS = IS + 1
                  AC%STORE(IS)%PTR(1)      = IP
                  AC%STORE(IS)%COL(1)      = JC
                  AC%STORE(IS)%VAL_ORIG(1) = AC%VAL(IP)
                  AC%STORE(IS)%VAL_COND(1) = 0.0_EB
                  AC%STORE(IS)%NROW        = JC
                  AC%STORE(IS)%NCOL        = 1
                  EXIT
               ENDIF
            ENDDO
      
            IF (.NOT.TWO_D) THEN
               JC = NC - L%NX * L%NY
               DO IP = AC%ROW(JC)+1, AC%ROW(JC+1)-1
                  IF (AC%COL(IP) == NC) THEN
                     IS = IS + 1
                     AC%STORE(IS)%PTR(1)      = IP
                     AC%STORE(IS)%COL(1)      = JC
                     AC%STORE(IS)%VAL_ORIG(1) = AC%VAL(IP)
                     AC%STORE(IS)%VAL_COND(1) = 0.0_EB
                     AC%STORE(IS)%NROW        = JC
                     AC%STORE(IS)%NCOL        = 1
                     EXIT
                  ENDIF
               ENDDO
            ENDIF
      
         ENDIF LAST_CELL_IN_LAST_MESH_COMPACT_IF
      
         !> cycle boundary cells to check if there is a periodic communication partner whose stencil is coupled
         !> with the last cell of last mesh;
         !> this can be a cell on the opposite side of the own mesh or on a different mesh
         !> if such a cell exists, store corresponding matrix entry
         WALL_CELLS_COMPACT_LOOP1: DO IW = 1, NW
      
            IOR0 = D%WALL(IW)%IOR
            IF (TWO_D .AND. ABS(IOR0) == 2) CYCLE
      
            I    = D%WALL(IW)%IXW
            J    = D%WALL(IW)%IYW
            K    = D%WALL(IW)%IZW
      
            IF (TYPE_DISCRET == NSCARC_DISCRET_UNSTRUCTURED .AND. L%CELL_STATE(I, J, K) /= NSCARC_CELL_GASPHASE) CYCLE
      
            NOM  = D%WALL(IW)%NOM
            D%WALL(IW)%ICW =D%CELL_NUMBER(I, J, K)
            IC =D%CELL_NUMBER(I, J, K)
      
            PERIODIC_NBR_COMPACT_IF1: IF (NOM == NMESHES) THEN
               ICE = D%WALL(IW)%ICE(1)                            !> adjacent ghost cell number
               ICN = D%PTR%ICE_TO_ICN(ICE)                        !> get column index of neighboring offdiagonal matrix entry
               IF (ICN /= SCARC(NMESHES)%N_CELLS) CYCLE           !> if no relation to last cell in last mesh, cycle
               DO IP = AC%ROW(ICN)+1, AC%ROW(ICN+1)-1
                  IF (AC%COL(IP) == ICE) THEN
                     IS = IS + 1
                     AC%STORE(IS)%PTR(1)      = IP
                     AC%STORE(IS)%COL(1)      = ICN
                     AC%STORE(IS)%VAL_ORIG(1) = AC%VAL(IP)
                     AC%STORE(IS)%VAL_COND(1) = 0.0_EB
                     AC%STORE(IS)%NROW        = JC
                     AC%STORE(IS)%NCOL        = 1
                     EXIT
                  ENDIF
               ENDDO
            ENDIF PERIODIC_NBR_COMPACT_IF1
         ENDDO WALL_CELLS_COMPACT_LOOP1
         AC%NSTORE = IS
      
      ENDIF NO_DIRIC_COMPACT_IF
      
      !>
      !> Set correct boundary conditions for system matrix
      !> Take care of whether the structured or unstructured discretization is used
      !>
      WALL_CELLS_COMPACT_LOOP2: DO IW = 1, NW
      
         IOR0 = D%WALL(IW)%IOR
         IF (TWO_D .AND. ABS(IOR0) == 2) CYCLE              !> cycle boundaries in y-direction for 2D-cases
      
         I    = D%WALL(IW)%IXW
         J    = D%WALL(IW)%IYW
         K    = D%WALL(IW)%IZW
      
         IF (TYPE_DISCRET == NSCARC_DISCRET_UNSTRUCTURED .AND. L%CELL_STATE(I, J, K) /= NSCARC_CELL_GASPHASE) CYCLE
      
         NOM  = D%WALL(IW)%NOM
         D%WALL(IW)%ICW =D%CELL_NUMBER(I, J, K)
         IC =D%CELL_NUMBER(I, J, K)
      
         SELECT CASE (ABS(IOR0))
            CASE (1)
               DBC= L%DXI2           ! Achtung: Wirklich richtig oder Mittelwert?
            CASE (2)
               DBC= L%DYI2
            CASE (3)
               DBC= L%DZI2
         END SELECT
      
         !> SPD-matrix with mixture of Dirichlet and Neumann BC's according to the SETTING of BTYPE
         IF (N_DIRIC_GLOBAL(NLEVEL_MIN) > 0) THEN
      
            IP = AC%ROW(IC)
            SELECT CASE (D%WALL(IW)%BTYPE)
               CASE (DIRICHLET)
                  AC%VAL(IP) = AC%VAL(IP) - DBC
               CASE (NEUMANN)
                  AC%VAL(IP) = AC%VAL(IP) + DBC
            END SELECT
      
         !> purely Neumann matrix
         ELSE
            IP = AC%ROW(IC)
            IF (D%WALL(IW)%BTYPE == NEUMANN) AC%VAL(IP) = AC%VAL(IP) + DBC
         ENDIF
      
      ENDDO WALL_CELLS_COMPACT_LOOP2
      
      !> If there are no Dirichlet BC's transform sytem into condensed one by replacing the
      !> matrix entries in last column and row by the stored ones (zeros and one at diaonal position)
      SETUP_CONDENSED_COMPACT_IF: IF (N_DIRIC_GLOBAL(NLEVEL_MIN) == 0 .AND. TYPE_PRECON /= NSCARC_RELAX_FFT) THEN
         DO IS = 1, AC%NSTORE
            DO ICOL = 1, AC%STORE(IS)%NCOL
               IP = AC%STORE(IS)%PTR(ICOL)
               AC%VAL(IP) = AC%STORE(IS)%VAL_COND(ICOL)
            ENDDO
         ENDDO
      ENDIF SETUP_CONDENSED_COMPACT_IF
      
   !>
   !> ---------- Matrix in Banded storage technique
   !>
   CASE (NSCARC_MATRIX_BANDED)

      AB => D%AB
      
      !>
      !> If A is a pure Neumann matrix, get neighboring cell indices of communicated stencil legs for condensed system
      !> also save values and column indices of last matrix row of last mesh
      !>
      NO_DIRIC_BANDED_IF: IF (N_DIRIC_GLOBAL(NLEVEL_MIN) == 0) THEN
         WRITE(*,*) 'TOFIX: NO_DIRIC_BANDED_IF'
      ENDIF NO_DIRIC_BANDED_IF

      !>
      !> Set correct boundary conditions for system matrix
      !> Take care of whether the structured or unstructured discretization is used
      !>
      WALL_CELLS_LOOP2: DO IW = 1, NW
      
         IOR0 = D%WALL(IW)%IOR
         IF (TWO_D .AND. ABS(IOR0) == 2) CYCLE              !> cycle boundaries in y-direction for 2D-cases
      
         I    = D%WALL(IW)%IXW
         J    = D%WALL(IW)%IYW
         K    = D%WALL(IW)%IZW
      
         IF (TYPE_DISCRET == NSCARC_DISCRET_UNSTRUCTURED .AND. L%CELL_STATE(I, J, K) /= NSCARC_CELL_GASPHASE) CYCLE
      
         NOM  = D%WALL(IW)%NOM
         D%WALL(IW)%ICW =D%CELL_NUMBER(I, J, K)
         IC =D%CELL_NUMBER(I, J, K)
      
         SELECT CASE (ABS(IOR0))
            CASE (1)
               DBC= L%DXI2           ! Achtung: Wirklich richtig oder Mittelwert?
            CASE (2)
               DBC= L%DYI2
            CASE (3)
               DBC= L%DZI2
         END SELECT
      
         !> SPD-matrix with mixture of Dirichlet and Neumann BC's according to the SETTING of BTYPE
         IF (N_DIRIC_GLOBAL(NLEVEL_MIN) > 0) THEN
      
            SELECT CASE (D%WALL(IW)%BTYPE)
               CASE (DIRICHLET)
                  AB%VAL(AB%POS(0), IC) = AB%VAL(AB%POS(0), IC) - DBC
               CASE (NEUMANN)
                  AB%VAL(AB%POS(0), IC) = AB%VAL(AB%POS(0), IC) + DBC
            END SELECT
      
         !> purely Neumann matrix
         ELSE
            IF (D%WALL(IW)%BTYPE == NEUMANN) AB%VAL(AB%POS(0),IC) = AB%VAL(AB%POS(0),IC) + DBC
         ENDIF
      
      ENDDO WALL_CELLS_LOOP2
      
      !> If there are no Dirichlet BC's transform sytem into condensed one by replacing the
      !> matrix entries in last column and row by the stored ones (zeros and one at diaonal position)
      SETUP_CONDENSED_BANDED_IF: IF (N_DIRIC_GLOBAL(NLEVEL_MIN) == 0 .AND. TYPE_PRECON /= NSCARC_RELAX_FFT) THEN
         WRITE(*,*) 'TOFIX: NO_DIRIC_BANDED_IF2'
      ENDIF SETUP_CONDENSED_BANDED_IF

END SELECT MATRIX_TYPE_SELECT

END SUBROUTINE SCARC_SETUP_BOUNDARY


!> ------------------------------------------------------------------------------------------------
!> Setup condensed system in case of periodic or pure Neumann boundary conditions
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_CONDENSING (NV, NL, ITYPE)
INTEGER, INTENT(IN) :: NV, NL, ITYPE
INTEGER :: NM, NOM, IW, IWG, IWL, NWL, IFACE, ICN, ICE, ICW, JC, NC, IS
REAL(EB) :: RHS_END
REAL(EB), POINTER, DIMENSION(:) :: VC
TYPE (SCARC_LEVEL_TYPE), POINTER :: L=>NULL(), OL=>NULL()
TYPE (SCARC_DISCRET_TYPE), POINTER :: D=>NULL(), OD=>NULL()
TYPE (SCARC_NEIGHBOR_TYPE), POINTER :: OS=>NULL()
TYPE (SCARC_POINTER_TYPE), POINTER :: PTR=>NULL(), OPTR=>NULL()
TYPE (SCARC_MATRIX_COMPACT_TYPE), POINTER :: AC=>NULL()
TYPE (SCARC_WALL_TYPE), POINTER :: WC=>NULL()

IF (N_DIRIC_GLOBAL(NLEVEL_MIN) > 0 .OR. TYPE_PRECON == NSCARC_RELAX_FFT) RETURN

!> In last mesh:
!> Subtract B*RHS(end) for internal legs of stencil
LOCAL_REAL = 0.0_EB
IF (UPPER_MESH_INDEX == NMESHES) THEN

   L  => SCARC(NMESHES)%LEVEL(NL)
   SELECT CASE(TYPE_DISCRET)
      CASE (NSCARC_DISCRET_STRUCTURED)
         D => L%SD
      CASE (NSCARC_DISCRET_UNSTRUCTURED)
         D => L%UD
   END SELECT

   AC => D%AC
   NC =  D%NC_LOCAL(NMESHES)

   VC => POINT_TO_VECTOR(NMESHES, NL, NV)

   !> process last column entries of all rows except of last one
   !> for those rows only one matrix entry was stored, namely that one which connects to the last cell
   DO IS = 2, AC%NSTORE
      JC = AC%STORE(IS)%COL(1)
      IF (JC < NC) THEN
         VC(JC) = VC(JC) - AC%STORE(IS)%VAL_ORIG(1)*VC(NC)
      ENDIF
   ENDDO

   LOCAL_REAL(NMESHES) = VC(NC)    !> store last entry of RHS
   VC(NC) = 0.0_EB                 !> set last entry of last mesh to zero

ENDIF

IF (ITYPE == 0) RETURN

!> broadcast last RHS-value of last cell in last mesh to all meshes
RHS_END = SCARC_BROADCAST_REAL(NSCARC_BROADCAST_LAST)
DO NM = 1, NMESHES
   SCARC(NM)%RHS_END = RHS_END
ENDDO

!>
!> Only in case of periodic BC's:
!> Subtract B*RHS(end) for corresponding entries of all periodic communication partners
!>
DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   L => SCARC(NM)%LEVEL(NL)
   SELECT CASE(TYPE_DISCRET)
      CASE (NSCARC_DISCRET_STRUCTURED)
         D => L%SD
      CASE (NSCARC_DISCRET_UNSTRUCTURED)
         D => L%UD
   END SELECT
   PTR => D%PTR

   SNODE = PROCESS(NM)
   RNODE = PROCESS(NMESHES)

   OS => SCARC(NM)%OSCARC(NMESHES)
   IF (OS%NICMAX_S==0 .OR. OS%NICMAX_R==0) CYCLE             !> if there is no communication with last mesh, cycle

   OL => SCARC(NM)%OSCARC(NMESHES)%LEVEL(NL)
   SELECT CASE(TYPE_DISCRET)
      CASE (NSCARC_DISCRET_STRUCTURED)
         OD => OL%SD
      CASE (NSCARC_DISCRET_UNSTRUCTURED)
         OD => OL%UD
   END SELECT
   OPTR => OD%PTR
   VC => POINT_TO_VECTOR (NM, NL, NV)

   !> subtract B*RHS(end) at corresponding positions
   DO IFACE = 1, 6                                          !> check if this face has connection to last cell

      IWL = OL%SUBDIVISION(1, FACE_ORIENTATION(IFACE))      !> first wall cell on this face
      NWL = OL%SUBDIVISION(2, FACE_ORIENTATION(IFACE))      !> number of wall cells on this face

      DO IW = IWL, IWL + NWL - 1

         WC => D%WALL(IWG)

         IWG = OPTR%IWL_TO_IWG(IW)                          !> corresponding global wall cell number
         ICE = WC%ICE(1)                           !> adjacent ghost cell number
         ICW = WC%ICW                              !> adjacent internal cell number
         NOM = WC%NOM                              !> neighbor for that wall cell

         IF (NOM > 0) THEN
            IF (NOM /= NMESHES) CYCLE                       !> only check for common matrix entries with last mesh
            ICN = PTR%ICE_TO_ICN(ICE)                       !> get column index of neighboring offdiagonal matrix entry
            IF (ICN /= SCARC(NMESHES)%N_CELLS) CYCLE        !> if no relation to last cell in last mesh, cycle
            VC(ICW) = VC(ICW) - PTR%ICE_TO_VAL(ICE)*SCARC(NM)%RHS_END
         ENDIF

      ENDDO
   ENDDO
ENDDO

END SUBROUTINE SCARC_SETUP_CONDENSING


!> ----------------------------------------------------------------------------------------------------
!> Initialize global solver methods (CG/BICG/GMG/AMG) and corresponding level structures
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_METHODS
INTEGER :: NSTACK

SELECT_METHOD: SELECT CASE(TYPE_METHOD)

   !> ---------------- McKenny-Greengard-Mayo method (MGM) --------------------
   CASE (NSCARC_METHOD_MGM)

      !> Allocate velocity vectors along internal obstructions for the setting of internal BC's
      CALL SCARC_SETUP_MGM(NLEVEL_MIN, NLEVEL_MIN)

      !> ------- First part of method: Setup CG solver for inhomogeneous problem on structured discretization
      TYPE_DISCRET = NSCARC_DISCRET_STRUCTURED

      NSTACK = NSCARC_STACK_ROOT
      STACK(NSTACK)%SOLVER => MAIN_CG_STRUCTURED
      CALL SCARC_SETUP_KRYLOV(NSCARC_SOLVER_MAIN, NSCARC_SCOPE_GLOBAL, NSCARC_STAGE_ONE, NSTACK, NLEVEL_MIN, NLEVEL_MIN)

      !> use FFT-preconditioning
      NSTACK = NSTACK + 1
      STACK(NSTACK)%SOLVER => PRECON_FFT
      CALL SCARC_SETUP_PRECON(NSTACK, NSCARC_SCOPE_LOCAL)

      !> ------- Second part of method: Setup CG solver for homogeneous problem on unstructured discretization
      TYPE_DISCRET = NSCARC_DISCRET_UNSTRUCTURED

      NSTACK = NSTACK + 1
      STACK(NSTACK)%SOLVER => MAIN_CG_UNSTRUCTURED
      CALL SCARC_SETUP_KRYLOV(NSCARC_SOLVER_MAIN, NSCARC_SCOPE_GLOBAL, NSCARC_STAGE_ONE, NSTACK, NLEVEL_MIN, NLEVEL_MIN)

      !> for a first proof of concept only use SSOR-preconditioning (may be extended later to other preconditioners)
      NSTACK = NSTACK + 1
      STACK(NSTACK)%SOLVER => PRECON_SSOR
      CALL SCARC_SETUP_PRECON(NSTACK, NSCARC_SCOPE_LOCAL)

   !>
   !> ------------------ Global Krylov method -------------------------------------
   !>
   CASE (NSCARC_METHOD_KRYLOV)

      !> Setup basic CG solver
      NSTACK = NSCARC_STACK_ROOT
      STACK(NSTACK)%SOLVER => MAIN_CG
      CALL SCARC_SETUP_KRYLOV(NSCARC_SOLVER_MAIN, NSCARC_SCOPE_GLOBAL, NSCARC_STAGE_ONE, NSTACK, NLEVEL_MIN, NLEVEL_MIN)

      !> Setup preconditioner for Krylov solver
      NSTACK = NSTACK + 1
      SELECT_KRYLOV_PRECON: SELECT CASE (TYPE_PRECON)

         !> Jacobi-preconditioning (acting locally by default)
         CASE (NSCARC_RELAX_JACOBI)                                
            STACK(NSTACK)%SOLVER => PRECON_JACOBI
            CALL SCARC_SETUP_PRECON(NSTACK, NSCARC_SCOPE_LOCAL)

         !> SSOR-preconditioning (acting locally by default)
         CASE (NSCARC_RELAX_SSOR)                                
            STACK(NSTACK)%SOLVER => PRECON_SSOR
            CALL SCARC_SETUP_PRECON(NSTACK, NSCARC_SCOPE_LOCAL)

         !> GS-preconditioning in matrix form (acting locally by default)
         CASE (NSCARC_RELAX_GSM)                                
            STACK(NSTACK)%SOLVER => PRECON_GSM
            CALL SCARC_SETUP_PRECON(NSTACK, NSCARC_SCOPE_LOCAL)
            CALL SCARC_SETUP_GSM(NLEVEL_MIN, NLEVEL_MAX)

         !> SGS-preconditioning in matrix form (acting locally by default)
         CASE (NSCARC_RELAX_SGSM)                                
            STACK(NSTACK)%SOLVER => PRECON_SGSM
            CALL SCARC_SETUP_PRECON(NSTACK, NSCARC_SCOPE_LOCAL)
            CALL SCARC_SETUP_SGSM(NLEVEL_MIN, NLEVEL_MAX)

         !> SOR-preconditioning in matrix form (acting locally by default)
         CASE (NSCARC_RELAX_SORM)                                
            STACK(NSTACK)%SOLVER => PRECON_SORM
            CALL SCARC_SETUP_PRECON(NSTACK, NSCARC_SCOPE_LOCAL)
            CALL SCARC_SETUP_SORM(NLEVEL_MIN, NLEVEL_MAX, NSTACK)

         !> SSOR-preconditioning in matrix form (acting locally by default)
         CASE (NSCARC_RELAX_SSORM)                                
            STACK(NSTACK)%SOLVER => PRECON_SSORM
            CALL SCARC_SETUP_PRECON(NSTACK, NSCARC_SCOPE_LOCAL)
            CALL SCARC_SETUP_SSORM(NLEVEL_MIN, NLEVEL_MAX, NSTACK)

         !> ILU(0)-preconditioning in matrix form (acting locally by default)
         CASE (NSCARC_RELAX_ILU)                                
            STACK(NSTACK)%SOLVER => PRECON_ILU
            CALL SCARC_SETUP_PRECON(NSTACK, NSCARC_SCOPE_LOCAL)
            CALL SCARC_SETUP_ILU(NLEVEL_MIN, NLEVEL_MAX)

         !> FFT-preconditioning (acting locally by default)
         CASE (NSCARC_RELAX_FFT)                                
            STACK(NSTACK)%SOLVER => PRECON_FFT
            CALL SCARC_SETUP_PRECON(NSTACK, NSCARC_SCOPE_LOCAL)
            CALL SCARC_SETUP_FFT(NLEVEL_MIN, NLEVEL_MIN)

#ifdef WITH_MKL
         !> LU-preconditioning based on MKL (either locally or globally acting depending on user specification)
         CASE (NSCARC_RELAX_MKL)                                 
            STACK(NSTACK)%SOLVER => PRECON_MKL

            SELECT CASE(TYPE_PRECON_SCOPE)

               !> globally acting - call global CLUSTER_SPARSE_SOLVER from MKL
               CASE (NSCARC_SCOPE_GLOBAL)
                  CALL SCARC_SETUP_PRECON(NSTACK, NSCARC_SCOPE_GLOBAL)
                  CALL SCARC_SETUP_CLUSTER(NLEVEL_MIN, NLEVEL_MIN)

               !> locally acting - call global PARDISO solver from MKL
               CASE (NSCARC_SCOPE_LOCAL)
                  CALL SCARC_SETUP_PRECON(NSTACK, NSCARC_SCOPE_LOCAL)
                  CALL SCARC_SETUP_PARDISO(NLEVEL_MIN, NLEVEL_MIN)

            END SELECT
#endif

         !> Preconditioning by Geometric multigrid,
         !> either locally or globally acting, depending on user specification stored in TYPE_PRECON_SCOPE
         CASE (NSCARC_RELAX_GMG)                                

            STACK(NSTACK)%SOLVER => PRECON_GMG
            CALL SCARC_SETUP_PRECON(NSTACK, TYPE_PRECON_SCOPE)
            CALL SCARC_SETUP_MULTIGRID(NSCARC_SOLVER_PRECON, TYPE_PRECON_SCOPE, NSCARC_STAGE_TWO, NSTACK, &
                                       NLEVEL_MIN, NLEVEL_MAX)

            NSTACK = NSTACK + 1
            SELECT CASE (TYPE_SMOOTH)

               !> Jacobi-smoothing (acting locally by default)
               CASE (NSCARC_RELAX_JACOBI)                              
                  STACK(NSTACK)%SOLVER => SMOOTH_JACOBI
                  CALL SCARC_SETUP_SMOOTH(NSTACK, NSCARC_SCOPE_LOCAL)

               !> SSOR-smoothing (acting locally by default)
               CASE (NSCARC_RELAX_SSOR)                                
                  STACK(NSTACK)%SOLVER => SMOOTH_SSOR
                  CALL SCARC_SETUP_SMOOTH(NSTACK, NSCARC_SCOPE_LOCAL)

               !> FFT-smoothing (acting locally by default)
               CASE (NSCARC_RELAX_FFT)                            
                  STACK(NSTACK)%SOLVER => SMOOTH_FFT
                  CALL SCARC_SETUP_SMOOTH(NSTACK, NSCARC_SCOPE_LOCAL)
                  CALL SCARC_SETUP_FFT(NLEVEL_MIN, NLEVEL_MAX)
#ifdef WITH_MKL
               !> LU-smoothing (acting locally by default)
               CASE (NSCARC_RELAX_MKL)                                
                  STACK(NSTACK)%SOLVER => SMOOTH_MKL
                  CALL SCARC_SETUP_SMOOTH(NSTACK, NSCARC_SCOPE_LOCAL)
                  CALL SCARC_SETUP_PARDISO(NLEVEL_MIN, NLEVEL_MIN)
#endif
            END SELECT

            !> Coarse grid solver (same scope of action as calling GMG)
            NSTACK = NSTACK + 1                                        
            CALL SCARC_SETUP_COARSE_SOLVER(NSCARC_STAGE_TWO, TYPE_PRECON_SCOPE, NSTACK, NLEVEL_MAX, NLEVEL_MAX)

      END SELECT SELECT_KRYLOV_PRECON

      !> If two-level Krylov, allocate intermediate structures for interpolation and workspace for global coarse solver
      IF (HAS_TWO_LEVELS) THEN

         CALL SCARC_SETUP_INTERPOLATION(NSCARC_STAGE_ONE, NLEVEL_MIN+1, NLEVEL_MAX)

         NSTACK = NSTACK + 1
         CALL SCARC_SETUP_COARSE_SOLVER(NSCARC_STAGE_ONE, NSCARC_SCOPE_GLOBAL, NSTACK, NLEVEL_MAX, NLEVEL_MAX)

      ENDIF

   !>
   !> ------------------ Global Multigrid method -------------------------------------
   !>
   CASE (NSCARC_METHOD_MULTIGRID)

      NSTACK = NSCARC_STACK_ROOT
      STACK(NSTACK)%SOLVER => MAIN_GMG
      CALL SCARC_SETUP_MULTIGRID(NSCARC_SOLVER_MAIN, NSCARC_SCOPE_GLOBAL, NSCARC_STAGE_ONE, NSTACK, NLEVEL_MIN, NLEVEL_MAX)

      NSTACK = NSTACK + 1
      SELECT CASE(TYPE_SMOOTH)

         !> Jacobi-smoothing (acting locally by default)
         CASE (NSCARC_RELAX_JACOBI)                         
            STACK(NSTACK)%SOLVER => SMOOTH_JACOBI
            CALL SCARC_SETUP_SMOOTH(NSTACK, NSCARC_SCOPE_LOCAL)

         !> SSOR-smoothing (acting locally by default)
         CASE (NSCARC_RELAX_SSOR)                                     
            STACK(NSTACK)%SOLVER => SMOOTH_SSOR
            CALL SCARC_SETUP_SMOOTH(NSTACK, NSCARC_SCOPE_LOCAL)

         !> FFT-smoothing (acting locally by default)
         CASE (NSCARC_RELAX_FFT)                                     
            STACK(NSTACK)%SOLVER => SMOOTH_FFT
            CALL SCARC_SETUP_SMOOTH(NSTACK, NSCARC_SCOPE_LOCAL)
            CALL SCARC_SETUP_FFT(NLEVEL_MIN, NLEVEL_MAX)

#ifdef WITH_MKL
         !> smoothing by LU-decomposition
         CASE (NSCARC_RELAX_MKL)                   
            CALL SCARC_SETUP_SMOOTH(NSTACK, TYPE_SMOOTH_SCOPE)
            
            SELECT CASE(TYPE_SMOOTH_SCOPE)

               !> globally acting - call global CLUSTER_SPARSE_SOLVER on MKL
               CASE (NSCARC_SCOPE_GLOBAL)
WRITE(MSG%LU_VERBOSE,*) 'B: SELECT GLOBAL SMOOTHING'
                  STACK(NSTACK)%SOLVER => SMOOTH_MKL
                  CALL SCARC_SETUP_CLUSTER(NLEVEL_MIN, NLEVEL_MIN)

               !> locally acting - call local PARDISO solvers based on MKL
               CASE (NSCARC_SCOPE_LOCAL)
WRITE(MSG%LU_VERBOSE,*) 'B: SELECT LOCAL SMOOTHING'
                  STACK(NSTACK)%SOLVER => SMOOTH_MKL
                  CALL SCARC_SETUP_PARDISO(NLEVEL_MIN, NLEVEL_MIN)
   
            END SELECT
#endif

      END SELECT

      !> Globally acting coarse grid solver
      NSTACK = NSTACK + 1                                              
      CALL SCARC_SETUP_COARSE_SOLVER(NSCARC_STAGE_ONE, NSCARC_SCOPE_GLOBAL, NSTACK, NLEVEL_MAX, NLEVEL_MAX)


#ifdef WITH_MKL
   !> ------------------ MKL method -------------------------------------
   CASE (NSCARC_METHOD_LU)

      NSTACK = NSCARC_STACK_ROOT
      STACK(NSTACK)%SOLVER => MAIN_LU

      CALL SCARC_SETUP_MKL(NSCARC_SOLVER_MAIN, NSCARC_SCOPE_GLOBAL, NSCARC_STAGE_ONE, NSTACK, NLEVEL_MIN, NLEVEL_MIN)

      !> In the multi-mesh case use CLUSTER_SPARSE_SOLVER, else PARDISO solver (only on finest grid level)
      IF (NMESHES > 1) THEN
         CALL SCARC_SETUP_CLUSTER(NLEVEL_MIN, NLEVEL_MIN)
      ELSE
         CALL SCARC_SETUP_PARDISO(NLEVEL_MIN, NLEVEL_MIN)
      ENDIF

#endif

END SELECT SELECT_METHOD

!> Store total number of stack entries (used solvers)
N_STACK_TOTAL = NSTACK

END SUBROUTINE SCARC_SETUP_METHODS

!> ----------------------------------------------------------------------------------------------------
!> Setup references to solution vectors related to used scope (main/relax)
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_REFERENCES(BX, BB, BD, BR, BV, BY, BZ, NSTACK)
LOGICAL, INTENT(IN) :: BX, BB, BD, BR, BV, BY, BZ
INTEGER, INTENT(IN) :: NSTACK
TYPE (SCARC_SOLVER_TYPE), POINTER :: SV=>NULL()

SV  => STACK(NSTACK)%SOLVER

SELECT CASE (SV%TYPE_STAGE)

   !> Solver from working stage ONE, e.g. main Krylov or MG solver
   CASE (NSCARC_STAGE_ONE)
      IF (BX) SV%X = NSCARC_VECTOR_ONE_X
      IF (BB) SV%B = NSCARC_VECTOR_ONE_B
      IF (BD) SV%D = NSCARC_VECTOR_ONE_D
      IF (BR) SV%R = NSCARC_VECTOR_ONE_R
      IF (BV) SV%V = NSCARC_VECTOR_ONE_V
      IF (BY) SV%Y = NSCARC_VECTOR_ONE_Y
      IF (BZ) SV%Z = NSCARC_VECTOR_ONE_Z
#ifdef WITH_SCARC_DEBUG
      SV%E = NSCARC_VECTOR_ONE_E
#endif

   !> Solver from working stage TWO, e.g. MG solver as preconditioner
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
END SELECT

END SUBROUTINE SCARC_SETUP_REFERENCES

!> ----------------------------------------------------------------------------------------------------
!> Allocate and initialize vectors for Krylov method
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_VECTORS()
INTEGER :: NM, NSTACK, NL
TYPE (SCARC_LEVEL_TYPE) , POINTER :: L=>NULL()
TYPE (SCARC_SOLVER_TYPE), POINTER :: SV=>NULL()
TYPE (SCARC_STAGE_TYPE) , POINTER :: ST=>NULL()

DO NSTACK = 1, N_STACK_TOTAL

   SV  => STACK(NSTACK)%SOLVER

   DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
      DO NL = SV%TYPE_NLMIN, SV%TYPE_NLMAX

         L  => SCARC(NM)%LEVEL(NL)
         ST => SCARC(NM)%LEVEL(NL)%STAGE(SV%TYPE_STAGE)

         IF (SV%X /= NSCARC_UNDEFINED_INT) CALL SCARC_ALLOCATE_REAL1(ST%X, 1, L%NCE, NSCARC_INIT_ZERO, 'X')
         IF (SV%B /= NSCARC_UNDEFINED_INT) CALL SCARC_ALLOCATE_REAL1(ST%B, 1, L%NCE, NSCARC_INIT_ZERO, 'B')
         IF (SV%D /= NSCARC_UNDEFINED_INT) CALL SCARC_ALLOCATE_REAL1(ST%D, 1, L%NCE, NSCARC_INIT_ZERO, 'D')
         IF (SV%R /= NSCARC_UNDEFINED_INT) CALL SCARC_ALLOCATE_REAL1(ST%R, 1, L%NCE, NSCARC_INIT_ZERO, 'R')
         IF (SV%V /= NSCARC_UNDEFINED_INT) CALL SCARC_ALLOCATE_REAL1(ST%V, 1, L%NCE, NSCARC_INIT_ZERO, 'V')
         IF (SV%Y /= NSCARC_UNDEFINED_INT) CALL SCARC_ALLOCATE_REAL1(ST%Y, 1, L%NCE, NSCARC_INIT_ZERO, 'Y')
         IF (SV%Z /= NSCARC_UNDEFINED_INT) CALL SCARC_ALLOCATE_REAL1(ST%Z, 1, L%NCE, NSCARC_INIT_ZERO, 'Z')

#ifdef WITH_SCARC_DEBUG
         IF (SV%E /= NSCARC_UNDEFINED_INT) CALL SCARC_ALLOCATE_REAL1(ST%E, 1, L%NCE, NSCARC_INIT_ZERO, 'E')
#endif

#ifdef WITH_MKL_FB
         IF (SV%D_FB /= NSCARC_UNDEFINED_INT) CALL SCARC_ALLOCATE_REAL1_FB(ST%D_FB, 1, L%NCE, NSCARC_INIT_ZERO, 'D_FB')
         IF (SV%R_FB /= NSCARC_UNDEFINED_INT) CALL SCARC_ALLOCATE_REAL1_FB(ST%R_FB, 1, L%NCE, NSCARC_INIT_ZERO, 'R_FB')
#endif

      ENDDO
   ENDDO
ENDDO

END SUBROUTINE SCARC_SETUP_VECTORS

!> ----------------------------------------------------------------------------------------------------
!> Allocate velocity vectors for the setting of internal BC's in MGM-method
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_MGM(NLMIN, NLMAX)
INTEGER, INTENT(IN) :: NLMIN, NLMAX
INTEGER :: NM, NL, NW1, NW2
TYPE (SCARC_LEVEL_TYPE), POINTER :: L=>NULL()
TYPE (SCARC_MGM_TYPE), POINTER :: MGM=>NULL()

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   DO NL = NLMIN, NLMAX

      L   => SCARC(NM)%LEVEL(NL)
      MGM => L%MGM
      
      NW1 = L%N_WALL_CELLS_EXT + 1
      NW2 = L%N_WALL_CELLS_EXT + L%N_WALL_CELLS_INT

      CALL SCARC_ALLOCATE_REAL1(MGM%US, NW1, NW2, NSCARC_INIT_ZERO, 'UM') 
      CALL SCARC_ALLOCATE_REAL1(MGM%VS, NW1, NW2, NSCARC_INIT_ZERO, 'VM') 
      CALL SCARC_ALLOCATE_REAL1(MGM%WS, NW1, NW2, NSCARC_INIT_ZERO, 'WM') 

      CALL SCARC_ALLOCATE_REAL1(MGM%U , NW1, NW2, NSCARC_INIT_ZERO, 'UM') 
      CALL SCARC_ALLOCATE_REAL1(MGM%V , NW1, NW2, NSCARC_INIT_ZERO, 'VM') 
      CALL SCARC_ALLOCATE_REAL1(MGM%W , NW1, NW2, NSCARC_INIT_ZERO, 'WM') 

   ENDDO
ENDDO

END SUBROUTINE SCARC_SETUP_MGM

!> ----------------------------------------------------------------------------------------------------
!> Allocate and initialize vectors for Krylov method
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_KRYLOV(NSOLVER, NSCOPE, NSTAGE, NSTACK, NLMIN, NLMAX)
INTEGER, INTENT(IN) :: NSOLVER, NSCOPE, NSTAGE, NSTACK, NLMIN, NLMAX
TYPE (SCARC_SOLVER_TYPE), POINTER :: SV=>NULL()

SV  => STACK(NSTACK)%SOLVER

!> Preset types for Krylov method
SV%TYPE_METHOD    = NSCARC_METHOD_KRYLOV
SV%TYPE_SOLVER    = NSOLVER
SV%TYPE_SCOPE     = NSCOPE
SV%TYPE_STAGE     = NSTAGE
SV%TYPE_NLMIN     = NLMIN
SV%TYPE_NLMAX     = NLMAX
SV%TYPE_DISCRET   = TYPE_DISCRET
SV%TYPE_INTERPOL  = TYPE_INTERPOL
SV%TYPE_ACCURACY  = TYPE_ACCURACY
SV%TYPE_PRECISION = TYPE_PRECISION

!> Preset iteration parameters for Krylov method
SELECT CASE(NSOLVER)

   !> -------------- Krylov method is used as main solver
   CASE (NSCARC_SOLVER_MAIN)                            

      SV%CNAME = 'SCARC_MAIN_KRYLOV'

      SV%EPS = SCARC_KRYLOV_ACCURACY
      SV%NIT = SCARC_KRYLOV_ITERATIONS

      SV%TYPE_RELAX   = TYPE_PRECON                    
      SV%TYPE_TWOLEVEL = TYPE_TWOLEVEL                  

   !> -------------- Krylov method is used as coarse grid solver solver
   CASE (NSCARC_SOLVER_COARSE)                 

      SV%CNAME = 'SCARC_COARSE_KRYLOV'

      SV%EPS = SCARC_COARSE_ACCURACY
      SV%NIT = SCARC_COARSE_ITERATIONS

      SV%TYPE_RELAX   = NSCARC_RELAX_SSOR               !> only use SSOR-preconditioning for coarse solver
      SV%TYPE_TWOLEVEL = NSCARC_TWOLEVEL_NONE             !> only use one level for coarse solver

   CASE DEFAULT                                          

      CALL SCARC_SHUTDOWN(NSCARC_ERROR_PARSE_INPUT, SCARC_NONE, NSOLVER)

END SELECT

!> Point to solution vectors (in corresponding scope)
CALL SCARC_SETUP_REFERENCES(.TRUE.,.TRUE.,.TRUE.,.TRUE.,.TRUE.,.TRUE.,.TRUE., NSTACK)
END SUBROUTINE SCARC_SETUP_KRYLOV

!> ----------------------------------------------------------------------------------------------------
!> Allocate and initialize vectors for Geometric Multigrid method
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_MULTIGRID(NSOLVER, NSCOPE, NSTAGE, NSTACK, NLMIN, NLMAX)
INTEGER, INTENT(IN) :: NSOLVER, NSCOPE, NSTAGE, NSTACK, NLMIN, NLMAX
TYPE (SCARC_SOLVER_TYPE), POINTER :: SV=>NULL()

SV  => STACK(NSTACK)%SOLVER

!> Preset types for Multigrid method
SV%TYPE_METHOD    = NSCARC_METHOD_MULTIGRID
SV%TYPE_SOLVER    = NSOLVER
SV%TYPE_SCOPE     = NSCOPE
SV%TYPE_STAGE     = NSTAGE
SV%TYPE_NLMIN     = NLMIN
SV%TYPE_NLMAX     = NLMAX
SV%TYPE_RELAX     = TYPE_SMOOTH
SV%TYPE_DISCRET   = TYPE_DISCRET
SV%TYPE_INTERPOL  = TYPE_INTERPOL
SV%TYPE_ACCURACY  = TYPE_ACCURACY
SV%TYPE_CYCLING   = TYPE_CYCLING
SV%TYPE_PRECISION = TYPE_PRECISION

SELECT CASE(NSOLVER)

   !> -------------- Multigrid method is used as main solver
   CASE (NSCARC_SOLVER_MAIN)                 
      SV%CNAME = 'SCARC_MAIN_GMG'

   !> -------------- Multigrid method is only used as preconditioner for global Krylov-method
   CASE (NSCARC_SOLVER_PRECON)                              
      SV%CNAME = 'SCARC_PRECON_GMG'

   CASE DEFAULT                                            
      CALL SCARC_SHUTDOWN(NSCARC_ERROR_PARSE_INPUT, SCARC_NONE, NSOLVER)

END SELECT

!> Preset iteration parameters for Multigrid method
SV%EPS = SCARC_MULTIGRID_ACCURACY
SV%NIT = SCARC_MULTIGRID_ITERATIONS

!> Point to solution vectors (in corresponding scope)
CALL SCARC_SETUP_REFERENCES(.TRUE.,.TRUE.,.FALSE.,.TRUE.,.TRUE.,.FALSE.,.TRUE., NSTACK)
END SUBROUTINE SCARC_SETUP_MULTIGRID

!> ----------------------------------------------------------------------------------------------------
!> Allocate and initialize vectors for MKL-methods
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_COARSE_SOLVER(NSTAGE, NSCOPE, NSTACK, NLMIN, NLMAX)
INTEGER, INTENT(IN)    :: NSCOPE, NSTAGE, NLMIN, NLMAX
INTEGER, INTENT(INOUT) :: NSTACK

SELECT_COARSE: SELECT CASE (TYPE_COARSE)

   !> -------------- CG-method is used as iterative coarse grid solver
   CASE (NSCARC_COARSE_ITERATIVE)

      !> initialize current stack position as CG-method
      STACK(NSTACK)%SOLVER => COARSE_KRYLOV           
      CALL SCARC_SETUP_KRYLOV(NSCARC_SOLVER_COARSE, NSCOPE, NSTAGE, NSTACK, NLMIN, NLMAX)

      !> and next stack position as its SSOR-preconditioner
      NSTACK = NSTACK + 1
      TYPE_PRECON = NSCARC_RELAX_SSOR
      STACK(NSTACK)%SOLVER => PRECON_SSOR                     
      CALL SCARC_SETUP_PRECON(NSTACK, NSCOPE)

   !> -------------- LU-decomposition (from MKL) is used as direct coarse grid solver
#ifdef WITH_MKL
   CASE (NSCARC_COARSE_DIRECT)

      !> Global scope in the multi-mesh case:
      !> initialize current stack position as global CLUSTER_SPARSE_SOLVER
      IF (NSCOPE == NSCARC_SCOPE_GLOBAL .AND. N_MPI_PROCESSES > 1) THEN
         STACK(NSTACK)%SOLVER => COARSE_CLUSTER               
         CALL SCARC_SETUP_MKL(NSCARC_SOLVER_COARSE, NSCOPE, NSTAGE, NSTACK, NLMIN, NLMAX)
         CALL SCARC_SETUP_CLUSTER(NLMIN, NLMAX)

      !> Local scope:
      !> initialize current stack position as PARDISO solver
      ELSE
         STACK(NSTACK)%SOLVER => COARSE_PARDISO
         CALL SCARC_SETUP_MKL(NSCARC_SOLVER_COARSE, NSCOPE, NSTAGE, NSTACK, NLMIN, NLMAX)
         CALL SCARC_SETUP_PARDISO(NLMIN, NLMAX)
      ENDIF
#endif
   CASE DEFAULT
      CALL SCARC_SHUTDOWN(NSCARC_ERROR_PARSE_INPUT, SCARC_NONE, TYPE_COARSE)
END SELECT SELECT_COARSE
END SUBROUTINE SCARC_SETUP_COARSE_SOLVER


#ifdef WITH_MKL
!> ----------------------------------------------------------------------------------------------------
!> Allocate and initialize vectors for LU-solvers (based on MKL)
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_MKL(NSOLVER, NSCOPE, NSTAGE, NSTACK, NLMIN, NLMAX)
INTEGER, INTENT(IN) :: NSOLVER, NSCOPE, NSTAGE, NSTACK, NLMIN, NLMAX
TYPE (SCARC_SOLVER_TYPE), POINTER :: SV=>NULL()

SV  => STACK(NSTACK)%SOLVER

SELECT CASE (NSOLVER)
   CASE (NSCARC_SOLVER_MAIN)
      SV%CNAME = 'SCARC_MAIN_LU'
   CASE (NSCARC_SOLVER_COARSE)
      SV%CNAME = 'SCARC_COARSE_LU'
   CASE DEFAULT
      CALL SCARC_SHUTDOWN(NSCARC_ERROR_PARSE_INPUT, SCARC_NONE, NSOLVER)
END SELECT

!> Preset types for LU-decomposition method
SV%TYPE_METHOD    = NSCARC_METHOD_LU
SV%TYPE_SOLVER    = NSOLVER
SV%TYPE_SCOPE     = NSCOPE
SV%TYPE_STAGE     = NSTAGE
SV%TYPE_NLMIN     = NLMIN
SV%TYPE_NLMAX     = NLMAX
SV%TYPE_PRECISION = TYPE_PRECISION
SV%TYPE_DISCRET   = TYPE_DISCRET

!> Point to solution vectors (in corresponding scope)
CALL SCARC_SETUP_REFERENCES(.TRUE.,.TRUE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE., NSTACK)

END SUBROUTINE SCARC_SETUP_MKL
#endif

!> ----------------------------------------------------------------------------------------------------
!> Allocate and initialize vectors for Krylov method
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_PRECON(NSTACK, NSCOPE)
INTEGER, INTENT(IN) :: NSTACK, NSCOPE
TYPE (SCARC_SOLVER_TYPE), POINTER :: SV=>NULL(), SVP=>NULL()

!> Point to stack entry of called preconditioner
SV => STACK(NSTACK)%SOLVER
SV%TYPE_SCOPE = NSCOPE
SV%EPS   =  SCARC_PRECON_ACCURACY
SV%NIT   =  SCARC_PRECON_ITERATIONS
SV%OMEGA =  SCARC_PRECON_OMEGA

!> Preset name of preconditioning method
SELECT CASE(TYPE_PRECON)
   CASE (NSCARC_RELAX_JACOBI)
      SV%CNAME = 'SCARC_PRECON_JACOBI'
   CASE (NSCARC_RELAX_SSOR)
      SV%CNAME = 'SCARC_PRECON_SSOR'
   CASE (NSCARC_RELAX_GSM)
      SV%CNAME = 'SCARC_PRECON_GSM'
   CASE (NSCARC_RELAX_SGSM)
      SV%CNAME = 'SCARC_PRECON_SGSM'
   CASE (NSCARC_RELAX_SORM)
      SV%CNAME = 'SCARC_PRECON_SORM'
   CASE (NSCARC_RELAX_SSORM)
      SV%CNAME = 'SCARC_PRECON_SSORM'
   CASE (NSCARC_RELAX_ILU)
      SV%CNAME = 'SCARC_PRECON_ILU'
      SV%OMEGA = 1.0_EB
   CASE (NSCARC_RELAX_FFT)
      SV%CNAME = 'SCARC_PRECON_FFT'
      SV%OMEGA = 1.0_EB
   CASE (NSCARC_RELAX_GMG)
      SV%CNAME = 'SCARC_PRECON_GMG'
#ifdef WITH_MKL
   CASE (NSCARC_RELAX_MKL)
      IF (NSCOPE == NSCARC_SCOPE_LOCAL) THEN
         SV%CNAME = 'SCARC_PRECON_PARDISO'
         SV%OMEGA = 1.0_EB
      ELSE
         SV%CNAME = 'SCARC_PRECON_CLUSTER'
         SV%NIT   = 1
         SV%OMEGA = 1.0_EB
      ENDIF
#endif
   CASE DEFAULT
      CALL SCARC_SHUTDOWN(NSCARC_ERROR_PARSE_INPUT, SCARC_NONE, TYPE_PRECON)
END SELECT

!> Preset types for preconditioner (use same as for calling solver)
SVP => STACK(NSTACK-1)%SOLVER

SV%TYPE_SOLVER    = SVP%TYPE_SOLVER
SV%TYPE_STAGE     = SVP%TYPE_STAGE
SV%TYPE_NLMIN     = SVP%TYPE_NLMIN
SV%TYPE_NLMAX     = SVP%TYPE_NLMAX
SV%TYPE_RELAX     = SVP%TYPE_RELAX
SV%TYPE_DISCRET   = SVP%TYPE_DISCRET
SV%TYPE_INTERPOL  = SVP%TYPE_INTERPOL
SV%TYPE_ACCURACY  = SVP%TYPE_ACCURACY
SV%TYPE_PRECISION = SVP%TYPE_PRECISION

!> Preset pointers for preconditioner (use same as for alling solver)
SV%X = SVP%X                                    
SV%B = SVP%B                             
SV%D = SVP%D
SV%R = SVP%R
SV%V = SVP%V
SV%Y = SVP%Y
SV%Z = SVP%Z

#ifdef WITH_SCARC_DEBUG
SV%E = SVP%E
#endif

#ifdef WITH_MKL_FB
SV%X_FB = SVP%X_FB
SV%B_FB = SVP%B_FB
SV%D_FB = SVP%D_FB
SV%R_FB = SVP%R_FB
SV%V_FB = SVP%V_FB
SV%Y_FB = SVP%Y_FB                              
SV%Z_FB = SVP%Z_FB                              
#endif

END SUBROUTINE SCARC_SETUP_PRECON

!> ----------------------------------------------------------------------------------------------------
!> Allocate and initialize vectors for Krylov method
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_SMOOTH(NSTACK, NSCOPE)
INTEGER, INTENT(IN) :: NSTACK, NSCOPE
TYPE (SCARC_SOLVER_TYPE), POINTER :: SV=>NULL(), SVP=>NULL()

!> Point to stack entry of called smoother
SV => STACK(NSTACK)%SOLVER
SV%TYPE_SCOPE = NSCOPE

!> Preset name of smoothing method
SELECT CASE(TYPE_SMOOTH)
   CASE (NSCARC_RELAX_JACOBI)
      SV%CNAME = 'SCARC_SMOOTH_JACOBI'
   CASE (NSCARC_RELAX_SSOR)
      SV%CNAME = 'SCARC_SMOOTH_SSOR'
   CASE (NSCARC_RELAX_ILU)
      SV%CNAME = 'SCARC_SMOOTH_ILU'
   CASE (NSCARC_RELAX_FFT)
      SV%CNAME = 'SCARC_SMOOTH_FFT'
#ifdef WITH_MKL
   CASE (NSCARC_RELAX_MKL)
      IF (NSCOPE == NSCARC_SCOPE_GLOBAL) THEN
         SV%CNAME = 'SCARC_SMOOTH_CLUSTER'
      ELSE
         SV%CNAME = 'SCARC_SMOOTH_PARDISO'
      ENDIF
#endif
   CASE DEFAULT
      CALL SCARC_SHUTDOWN(NSCARC_ERROR_PARSE_INPUT, SCARC_NONE, TYPE_SMOOTH)
END SELECT

!> Preset iteration parameters for Multigrid method
SV%EPS   = SCARC_SMOOTH_ACCURACY                 
SV%NIT   = SCARC_SMOOTH_ITERATIONS
SV%OMEGA = SCARC_SMOOTH_OMEGA

!> Preset types for preconditioner (use same descriptors as calling solver)
SVP => STACK(NSTACK-1)%SOLVER

SV%TYPE_SOLVER    = NSCARC_SOLVER_SMOOTH
SV%TYPE_STAGE     = SVP%TYPE_STAGE
SV%TYPE_NLMIN     = SVP%TYPE_NLMIN
SV%TYPE_NLMAX     = SVP%TYPE_NLMAX
SV%TYPE_RELAX     = SVP%TYPE_RELAX
SV%TYPE_DISCRET   = SVP%TYPE_DISCRET
SV%TYPE_INTERPOL  = SVP%TYPE_INTERPOL
SV%TYPE_ACCURACY  = SVP%TYPE_ACCURACY
SV%TYPE_PRECISION = SVP%TYPE_PRECISION

!> Preset references for preconditioner (use same pointers as calling solver)
SV%X = SVP%X                          
SV%B = SVP%B
SV%D = SVP%D
SV%R = SVP%R
SV%V = SVP%V
SV%Y = SVP%Y
SV%Z = SVP%Z

#ifdef WITH_SCARC_DEBUG
SV%E = SVP%E
SV%R = SVP%R
#endif

#ifdef WITH_MKL_FB
SV%X_FB = SVP%X_FB                        
SV%B_FB = SVP%B_FB
SV%D_FB = SVP%D_FB
SV%R_FB = SVP%R_FB
SV%V_FB = SVP%V_FB
SV%Y_FB = SVP%Y_FB
SV%Z_FB = SVP%Z_FB
#endif

END SUBROUTINE SCARC_SETUP_SMOOTH

!> ----------------------------------------------------------------------------------------------------
!> Allocate and initialize vectors for additive or multiplicative coarse grid
!> (corresponding to Schwarz domain decomposition method)
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_INTERPOLATION(NSTAGE, NLMIN, NLMAX)
INTEGER, INTENT(IN) :: NSTAGE, NLMIN, NLMAX
INTEGER :: NM, NL
TYPE (SCARC_LEVEL_TYPE), POINTER :: L=>NULL()
TYPE (SCARC_STAGE_TYPE), POINTER :: ST=>NULL()

MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   LEVEL_LOOP: DO NL = NLMIN, NLMAX

      L  => SCARC(NM)%LEVEL(NL)
      ST => SCARC(NM)%LEVEL(NL)%STAGE(NSTAGE)

      CALL SCARC_ALLOCATE_REAL1(ST%X, 1, L%NCE, NSCARC_INIT_ZERO, 'X')
      CALL SCARC_ALLOCATE_REAL1(ST%B, 1, L%NCE, NSCARC_INIT_ZERO, 'B')
      CALL SCARC_ALLOCATE_REAL1(ST%V, 1, L%NCE, NSCARC_INIT_ZERO, 'Q')
      CALL SCARC_ALLOCATE_REAL1(ST%R, 1, L%NCE, NSCARC_INIT_ZERO, 'W')
      CALL SCARC_ALLOCATE_REAL1(ST%Y, 1, L%NCE, NSCARC_INIT_ZERO, 'Y')
      CALL SCARC_ALLOCATE_REAL1(ST%Z, 1, L%NCE, NSCARC_INIT_ZERO, 'Z')

   ENDDO LEVEL_LOOP
ENDDO MESHES_LOOP

END SUBROUTINE SCARC_SETUP_INTERPOLATION

!> ----------------------------------------------------------------------------------------------------
!> Allocate and initialize vectors for blockwise FFT methods
!> New here: Perform own initialization of FFT based on H2CZIS/H3CZIS and use own SAVE and WORK arrays
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_FFT(NLMIN, NLMAX)
USE POIS, ONLY: H2CZIS, H3CZIS
INTEGER, INTENT(IN) :: NLMIN, NLMAX
INTEGER :: NM, NL, IERR = 0
TYPE (MESH_TYPE), POINTER :: M=>NULL()
TYPE (SCARC_TYPE), POINTER :: S=>NULL()
TYPE (SCARC_LEVEL_TYPE), POINTER :: L=>NULL()
TYPE (SCARC_FFT_TYPE), POINTER :: FFT=>NULL()

MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   LEVEL_LOOP: DO NL = NLMIN, NLMAX
                                 
      M   => MESHES(NM)
      S   => SCARC(NM)
      L   => SCARC(NM)%LEVEL(NL)
      FFT => SCARC(NM)%LEVEL(NL)%FFT

      !> Allocate working space for FFT routine
      FFT%LBC = M%LBC
      FFT%MBC = M%MBC
      FFT%NBC = M%NBC

WRITE(MSG%LU_DEBUG,*) 'FFT: --------> NM, NL, LBC, MBC, NBC =',NM, NL, M%LBC, M%MBC, M%NBC

      FFT%ITRN = L%NX+1
      IF (TWO_D) THEN
         FFT%JTRN = 1
      ELSE
         FFT%JTRN = L%NY+1
      ENDIF
      FFT%KTRN = L%NZ+1

      FFT%LSAVE = (FFT%ITRN+1)*FFT%JTRN*FFT%KTRN+7*FFT%ITRN+5*FFT%JTRN+6*FFT%KTRN+56
      FFT%LWORK = (FFT%ITRN+1)*FFT%JTRN*FFT%KTRN

      CALL SCARC_ALLOCATE_REAL1(FFT%SAVE1, -3, FFT%LSAVE, NSCARC_INIT_ZERO, 'FFT')
      CALL SCARC_ALLOCATE_REAL1(FFT%WORK ,  1, FFT%LWORK, NSCARC_INIT_ZERO, 'FFT')

      !> Allocate stretching vector (set to 1)
      CALL SCARC_ALLOCATE_REAL1(FFT%HX, 1, L%NX+1, NSCARC_INIT_ONE, 'FFT')

      !> Allocate RHS vector for FFT routine
      IF (L%NY == 1) THEN
         CALL SCARC_ALLOCATE_REAL3(FFT%PRHS, 1, L%NX+1, 1, 1,      1, L%NZ+1, NSCARC_INIT_ZERO, 'FFT')
      ELSE
         CALL SCARC_ALLOCATE_REAL3(FFT%PRHS, 1, L%NX+1, 1, L%NY+1, 1, L%NZ+1, NSCARC_INIT_ZERO, 'FFT')
      ENDIF

      !> Allocate boundary data vector for XS
      IF (L%NZ>1) THEN
         IF (L%NY  > 1) CALL SCARC_ALLOCATE_REAL2(FFT%BXS, 1, L%NY+1, 1, L%NZ+1, NSCARC_INIT_ZERO, 'BXS')
         IF (L%NY == 1) CALL SCARC_ALLOCATE_REAL2(FFT%BXS, 1,      1, 1, L%NZ+1, NSCARC_INIT_ZERO, 'BXS')
      ELSE
         CALL SCARC_ALLOCATE_REAL2(FFT%BXS, 1, L%NY+1, 1, 1, NSCARC_INIT_ZERO, 'BXS')
      ENDIF

      !> Allocate boundary data vector for XF
      IF (L%NZ>1) THEN
         IF (L%NY  > 1) CALL SCARC_ALLOCATE_REAL2(FFT%BXF, 1, L%NY+1, 1, L%NZ+1, NSCARC_INIT_ZERO, 'BXF')
         IF (L%NY == 1) CALL SCARC_ALLOCATE_REAL2(FFT%BXF, 1,      1, 1, L%NZ+1, NSCARC_INIT_ZERO, 'BXF')
      ELSE
         CALL SCARC_ALLOCATE_REAL2(FFT%BXF, 1, L%NY+1, 1, 1, NSCARC_INIT_ZERO, 'BXF')
      ENDIF

      !> Allocate boundary data vector for YS
      IF (L%NZ > 1) THEN
         CALL SCARC_ALLOCATE_REAL2(FFT%BYS, 1, L%NX+1,1, L%NZ+1, NSCARC_INIT_ZERO, 'BYS')
      ELSE
         CALL SCARC_ALLOCATE_REAL2(FFT%BYS, 1, L%NX+1,1,      1, NSCARC_INIT_ZERO, 'BYS')
      ENDIF

      !> Allocate boundary data vector for YF
      IF (L%NZ > 1) THEN
         CALL SCARC_ALLOCATE_REAL2(FFT%BYF, 1, L%NX+1,1, L%NZ+1, NSCARC_INIT_ZERO, 'BYF')
      ELSE
         CALL SCARC_ALLOCATE_REAL2(FFT%BYF, 1, L%NX+1,1,      1, NSCARC_INIT_ZERO, 'BYF')
      ENDIF

      !> Allocate boundary data vector for ZS
      IF (L%NY  > 1) CALL SCARC_ALLOCATE_REAL2(FFT%BZS, 1, L%NX+1, 1, L%NY+1, NSCARC_INIT_ZERO, 'BZS')
      IF (L%NY == 1) CALL SCARC_ALLOCATE_REAL2(FFT%BZS, 1, L%NX+1, 1,      1, NSCARC_INIT_ZERO, 'BZS')

      !> Allocate boundary data vector for ZF
      IF (L%NY  >1)  CALL SCARC_ALLOCATE_REAL2(FFT%BZF, 1, L%NX+1, 1, L%NY+1, NSCARC_INIT_ZERO, 'BZF')
      IF (L%NY == 1) CALL SCARC_ALLOCATE_REAL2(FFT%BZF, 1, L%NX+1, 1,      1, NSCARC_INIT_ZERO, 'BZF')

      IF (TWO_D) THEN
         CALL H2CZIS(S%XS,S%XF,S%IBAR,FFT%LBC,S%ZS,S%ZF,S%KBAR,FFT%NBC,FFT%HX,FFT%XLM,FFT%ITRN,IERR,FFT%SAVE1)
      ELSE
         CALL H3CZIS(S%XS,S%XF,S%IBAR,FFT%LBC,S%YS,S%YF,S%JBAR,FFT%MBC,S%ZS,S%ZF,S%KBAR,FFT%NBC,&
                     FFT%HX,FFT%XLM,FFT%ITRN,FFT%JTRN,IERR,FFT%SAVE1)
      ENDIF

   ENDDO LEVEL_LOOP
ENDDO MESHES_LOOP

END SUBROUTINE SCARC_SETUP_FFT

!> ----------------------------------------------------------------------------------------------------
!> Store SGS preconditioner in matrix form
!> Based on the following splitting of A = D - E - F
!> where :   D is the diagonal part
!>          -E is the strictly lower part
!>          -F is the strictly upper part
!> the SSOR-preconditioner in matrix form is defined
!>           B_GSM = (D - E) 
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_GSM(NLMIN, NLMAX)
INTEGER, INTENT(IN) :: NLMIN, NLMAX
INTEGER :: NM, NL, IC, JC, IPTR
TYPE (SCARC_LEVEL_TYPE), POINTER :: L=>NULL()
TYPE (SCARC_MATRIX_COMPACT_TYPE), POINTER :: AC=>NULL()
TYPE (SCARC_DISCRET_TYPE), POINTER :: D=>NULL()

MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   LEVEL_LOOP: DO NL = NLMIN, NLMAX

      L  => SCARC(NM)%LEVEL(NL)
      SELECT CASE(TYPE_DISCRET)
         CASE (NSCARC_DISCRET_STRUCTURED)
            D => L%SD
         CASE (NSCARC_DISCRET_UNSTRUCTURED)
            D => L%UD
      END SELECT
      AC => D%AC

      !> Allocate workspace for GSM decomposition of Poisson matrix 
      CALL SCARC_ALLOCATE_REAL1(AC%GSM, 1, AC%NA, NSCARC_INIT_ZERO, 'ILU')

      CELL_LOOP: DO IC = 1, L%NC
        
         DO IPTR = AC%ROW(IC), AC%ROW(IC+1)-1
            JC = AC%COL(IPTR)
            IF (JC > IC) CYCLE
            AC%GSM(IPTR) = AC%VAL(IPTR) 
         ENDDO
            
      ENDDO CELL_LOOP

   ENDDO LEVEL_LOOP
ENDDO MESHES_LOOP

#ifdef WITH_SCARC_DEBUG
DO NL = NLEVEL_MIN, NLEVEL_MAX
   CALL SCARC_DEBUG_QUANTITY(NSCARC_DEBUG_MATRIX, NL, 'GSM-Decomposition')
ENDDO
#endif

END SUBROUTINE SCARC_SETUP_GSM

!> ----------------------------------------------------------------------------------------------------
!> Store symmetric Gauss-Seidel preconditioner in matrix form
!> Based on the following splitting of A = D - E - F
!> where :   D is the diagonal part
!>          -E is the strictly lower part
!>          -F is the strictly upper part
!> the SGS-preconditioner in matrix form is defined
!>           B_SGSM = (D - E) D^{-1} (D - F)  =  (I - E D^{-1}) (D - F)
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_SGSM(NLMIN, NLMAX)
INTEGER, INTENT(IN) :: NLMIN, NLMAX
INTEGER :: NM, NL, IC, JC, IPTR
TYPE (SCARC_LEVEL_TYPE), POINTER :: L=>NULL()
TYPE (SCARC_DISCRET_TYPE), POINTER :: D=>NULL()
TYPE (SCARC_MATRIX_COMPACT_TYPE), POINTER :: AC=>NULL()

MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   LEVEL_LOOP: DO NL = NLMIN, NLMAX

      L  => SCARC(NM)%LEVEL(NL)
      SELECT CASE(TYPE_DISCRET)
         CASE (NSCARC_DISCRET_STRUCTURED)
            D => L%SD
         CASE (NSCARC_DISCRET_UNSTRUCTURED)
            D => L%UD
      END SELECT
      AC => D%AC

      !> Allocate workspace for SGSM decomposition of Poisson matrix 
      CALL SCARC_ALLOCATE_REAL1(AC%SGSM, 1, AC%NA, NSCARC_INIT_ZERO, 'ILU')

      CELL_LOOP: DO IC = 1, L%NC
        
         DO IPTR = AC%ROW(IC), AC%ROW(IC+1)-1
            JC = AC%COL(IPTR)

            !> Definition of lower elements
            !>  l(i,j) = omega * a(i,j)/a(j,j)
            IF (JC < IC)  AC%SGSM(IPTR) = AC%VAL(IPTR) / AC%VAL(AC%ROW(JC))

            !> Definition of diagonal element
            !> u(i,i) = a(i,i)
            IF (JC == IC) AC%SGSM(IPTR) = AC%VAL(IPTR)

            !> Definition of upper elements
            !> u(i,j) = omega * a(i,j)
            IF (JC > IC)  AC%SGSM(IPTR) = AC%VAL(IPTR) 
         ENDDO
            
      ENDDO CELL_LOOP

   ENDDO LEVEL_LOOP
ENDDO MESHES_LOOP

#ifdef WITH_SCARC_DEBUG
DO NL = NLEVEL_MIN, NLEVEL_MAX
   CALL SCARC_DEBUG_QUANTITY(NSCARC_DEBUG_MATRIX, NL, 'SGSM-Decomposition')
ENDDO
#endif

END SUBROUTINE SCARC_SETUP_SGSM

!> ----------------------------------------------------------------------------------------------------
!> Store SOR preconditioner in matrix form
!> Based on the following splitting of A = D - E - F
!> where :   D is the diagonal part
!>          -E is the strictly lower part
!>          -F is the strictly upper part
!> the SOR-preconditioner in matrix form is defined
!>           B_SORM = 1/omega * (D - omega * E) 
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_SORM(NLMIN, NLMAX, NSTACK)
INTEGER, INTENT(IN) :: NLMIN, NLMAX, NSTACK
REAL (EB) :: OMEGA
INTEGER :: NM, NL, IC, JC, IPTR
TYPE (SCARC_LEVEL_TYPE), POINTER :: L=>NULL()
TYPE (SCARC_DISCRET_TYPE), POINTER :: D=>NULL()
TYPE (SCARC_MATRIX_COMPACT_TYPE), POINTER :: AC=>NULL()

OMEGA = STACK(NSTACK)%SOLVER%OMEGA

MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   LEVEL_LOOP: DO NL = NLMIN, NLMAX

      L  => SCARC(NM)%LEVEL(NL)
      SELECT CASE(TYPE_DISCRET)
         CASE (NSCARC_DISCRET_STRUCTURED)
            D => L%SD
         CASE (NSCARC_DISCRET_UNSTRUCTURED)
            D => L%UD
      END SELECT
      AC => D%AC

      !> Allocate workspace for SORM decomposition of Poisson matrix 
      CALL SCARC_ALLOCATE_REAL1(AC%SORM, 1, AC%NA, NSCARC_INIT_ZERO, 'ILU')

      CELL_LOOP: DO IC = 1, L%NC
        
         DO IPTR = AC%ROW(IC), AC%ROW(IC+1)-1
            JC = AC%COL(IPTR)

            !> Definition of lower elements
            !>  l(i,j) = omega * a(i,j)/a(j,j)
            IF (JC < IC)  AC%SORM(IPTR) = AC%VAL(IPTR) 

            !> Definition of diagonal element
            !> u(i,i) = a(i,i)
            IF (JC == IC) AC%SORM(IPTR) = AC%VAL(IPTR)/OMEGA

         ENDDO
            
      ENDDO CELL_LOOP

   ENDDO LEVEL_LOOP
ENDDO MESHES_LOOP

#ifdef WITH_SCARC_DEBUG
DO NL = NLEVEL_MIN, NLEVEL_MAX
   CALL SCARC_DEBUG_QUANTITY(NSCARC_DEBUG_MATRIX, NL, 'SORM-Decomposition')
ENDDO
#endif

END SUBROUTINE SCARC_SETUP_SORM

!> ----------------------------------------------------------------------------------------------------
!> Store SSOR preconditioner in matrix form
!> Based on the following splitting of A = D - E - F
!> where :   D is the diagonal part
!>          -E is the strictly lower part
!>          -F is the strictly upper part
!> the SSOR-preconditioner in matrix form is defined
!>           B_SSOR = 1/(omega * (2-omega)) * (D - omega * E) D^{-1} (D - omega * F) 
!>                  = (I - omega E D^{-1}) * [1/(omega * (2-omega)) * D -  1/(2-omega) * F]
!> Defining the triangular matrices
!>               L  = I - omega E D^{-1} 
!>               U  = 1/(omega * (2-omega)) * D -  1/(2-omega) * F
!> the SSOR-preconditioning can be thought as the solution of two triangular systems
!> Both matrices can be stored as a single matrix that occupies the same amount of storage as A
!> where the same row and column pointers can be used as for A (identical pattern)
!> Note that the diagonal elements of L are 1 (and are omitted, only the diagonal of U is stored there)
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_SSORM(NLMIN, NLMAX, NSTACK)
INTEGER, INTENT(IN) :: NLMIN, NLMAX, NSTACK
INTEGER :: NM, NL, IC, JC, IPTR
REAL(EB) :: OMEGA, SCAL1, SCAL2
TYPE (SCARC_LEVEL_TYPE), POINTER :: L=>NULL()
TYPE (SCARC_DISCRET_TYPE), POINTER :: D=>NULL()
TYPE (SCARC_MATRIX_COMPACT_TYPE), POINTER :: AC=>NULL()

OMEGA = STACK(NSTACK)%SOLVER%OMEGA
SCAL1  = 1.0_EB / (OMEGA * (2.0_EB - OMEGA))
SCAL2  = 1.0_EB / (2.0_EB - OMEGA)

MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   LEVEL_LOOP: DO NL = NLMIN, NLMAX

      L  => SCARC(NM)%LEVEL(NL)
      SELECT CASE(TYPE_DISCRET)
         CASE (NSCARC_DISCRET_STRUCTURED)
            D => L%SD
         CASE (NSCARC_DISCRET_UNSTRUCTURED)
            D => L%UD
      END SELECT
      AC => D%AC

      !> Allocate workspace for SSORM decomposition of Poisson matrix 
      CALL SCARC_ALLOCATE_REAL1(AC%SSORM, 1, AC%NA, NSCARC_INIT_ZERO, 'ILU')

      CELL_LOOP: DO IC = 1, L%NC
        
         DO IPTR = AC%ROW(IC), AC%ROW(IC+1)-1

            JC = AC%COL(IPTR)

            !> Definition of lower elements:  l(i,j) = omega * a(i,j)/a(j,j)
            IF (JC < IC) AC%SSORM(IPTR) = OMEGA * AC%VAL(IPTR) / AC%VAL(AC%ROW(JC))

            !> Definition of diagonal element:  u(i,i) = a(i,i)
            IF (JC == IC) AC%SSORM(IPTR) = SCAL1 * AC%VAL(IPTR)

            !> Definition of upper elements:  u(i,j) = omega * a(i,j)
            IF (JC > IC) AC%SSORM(IPTR) = SCAL2 * AC%VAL(IPTR) 

         ENDDO
            
      ENDDO CELL_LOOP

   ENDDO LEVEL_LOOP
ENDDO MESHES_LOOP

!STACK(NSTACK)%SOLVER%OMEGA = 1.0_EB

#ifdef WITH_SCARC_DEBUG
DO NL = NLEVEL_MIN, NLEVEL_MAX
   CALL SCARC_DEBUG_QUANTITY(NSCARC_DEBUG_MATRIX, NL, 'SSORM-Decomposition')
ENDDO
#endif

END SUBROUTINE SCARC_SETUP_SSORM


!> ----------------------------------------------------------------------------------------------------
!> Allocate and initialize ILU(0) decomposition of Poisson matrix 
!> L- and U-parts are stored in the same array, diagonal elements of L are supposed to be 1
!> Based on Saad-algorithm 10.4 from 'Iterative Methods for Sparse Linear Systems':
!>   for i = 2 , ... , n do
!>      for k = 1 , ... , i-1 and for (i,k) in NZ(A) do
!>         compute a_ik = a_ik / a_kk
!>         for j = k+1 , ... , n and for (i,j) in NZ(A) do
!>            compute a_ij = a_ij - a_ik a_kj
!>         enddo
!>      enddo
!>   enddo
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_ILU(NLMIN, NLMAX)
INTEGER, INTENT(IN) :: NLMIN, NLMAX
INTEGER :: NM, NL, IC, JC, KC, IPTR, JPTR, KPTR, KPTR0
TYPE (SCARC_LEVEL_TYPE), POINTER :: L=>NULL()
TYPE (SCARC_DISCRET_TYPE), POINTER :: D=>NULL()
TYPE (SCARC_MATRIX_COMPACT_TYPE), POINTER :: AC=>NULL()

MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   LEVEL_LOOP: DO NL = NLMIN, NLMAX

      L  => SCARC(NM)%LEVEL(NL)
      SELECT CASE(TYPE_DISCRET)
         CASE (NSCARC_DISCRET_STRUCTURED)
            D => L%SD
         CASE (NSCARC_DISCRET_UNSTRUCTURED)
            D => L%UD
      END SELECT
      AC => D%AC

      !> Allocate ILU-part of Poisson matrix and preset it with Poisson matrix itself
      CALL SCARC_ALLOCATE_REAL1(AC%ILU, 1, AC%NA, NSCARC_INIT_ZERO, 'ILU')
      AC%ILU = AC%VAL

      CELL_LOOP: DO IC = 2, L%NC

         COLUMN_LOOP: DO IPTR = AC%ROW(IC), AC%ROW(IC+1)-1

            KC = AC%COL(IPTR)                        !> get number of neighboring cell
            IF (KC >= IC) CYCLE                      !> only consider neighbors with lower cell numbers than IC
            IF (AC%ILU(IPTR) == 0) CYCLE

            KPTR = AC%ROW(KC)                        !> get diagonal entry of neighbor
            AC%ILU(IPTR) = AC%ILU(IPTR)/AC%ILU(KPTR)           

            DO JPTR = AC%ROW(IC), AC%ROW(IC+1)-1

               JC = AC%COL(JPTR)
               IF (JC<=KC) CYCLE                     !> only consider neighbors with higher cell numbers than IC
               IF (AC%ILU(JPTR) == 0) CYCLE

               KPTR = -1
               DO KPTR0 = AC%ROW(KC), AC%ROW(KC+1)-1
                  IF (AC%COL(KPTR0) == JC) KPTR = KPTR0
               ENDDO
               IF (KPTR>0) AC%ILU(JPTR) = AC%ILU(JPTR) - AC%ILU(IPTR) * AC%ILU(KPTR)

            ENDDO

         ENDDO COLUMN_LOOP
      ENDDO CELL_LOOP

   ENDDO LEVEL_LOOP
ENDDO MESHES_LOOP

#ifdef WITH_SCARC_DEBUG
DO NL = NLEVEL_MIN, NLEVEL_MAX
   CALL SCARC_DEBUG_QUANTITY(NSCARC_DEBUG_MATRIX, NL, 'LU-Decomposition')
ENDDO
#endif

END SUBROUTINE SCARC_SETUP_ILU



#ifdef WITH_MKL
!> ------------------------------------------------------------------------------------------------
!> Initialize Pardiso solver
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_CLUSTER(NLMIN, NLMAX)
INTEGER, INTENT(IN) :: NLMIN, NLMAX
INTEGER :: NM, NL, I !, IC, IP
REAL (EB) :: TNOW
#ifdef WITH_MKL_FB
REAL (FB) :: DUMMY(1)=0.0_FB
#else
REAL (EB) :: DUMMY(1)=0.0_EB
#endif
TYPE (SCARC_LEVEL_TYPE), POINTER :: L=>NULL()
TYPE (SCARC_DISCRET_TYPE), POINTER :: D=>NULL()
TYPE (SCARC_MKL_TYPE), POINTER :: MKL=>NULL()
TYPE (SCARC_MATRIX_COMPACT_TYPE), POINTER :: AC_SYM=>NULL()

TNOW = CURRENT_TIME()

MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   LEVEL_LOOP: DO NL = NLMIN, NLMAX

      L   => SCARC(NM)%LEVEL(NL)
      MKL => L%MKL

      SELECT CASE(TYPE_DISCRET)
         CASE (NSCARC_DISCRET_STRUCTURED)
            D => L%SD
         CASE (NSCARC_DISCRET_UNSTRUCTURED)
            D => L%UD
      END SELECT
      AC_SYM => D%AC_SYM

      !> Allocate workspace for parameters needed in MKL-routine
      CALL SCARC_ALLOCATE_INT1(MKL%IPARM, 1, 64, NSCARC_INIT_ZERO, 'MKL_IPARM')

      !> Allocate workspace for pointers needed in MKL-routine
      IF (.NOT.ALLOCATED(MKL%CT)) THEN
         ALLOCATE(MKL%CT(64), STAT=IERROR)
         CALL CHKMEMERR ('SCARC', 'CT', IERROR)
         DO I=1,64
            MKL%CT(I)%DUMMY = 0
         ENDDO
      ENDIF

      !> Define corresponding parameters
      !> Note: IPARM-vectory is allocate from 1:64, not from 0:63
      MKL%NRHS   =  1         ! one right hand side
      MKL%MAXFCT =  1         ! one matrix
      MKL%MNUM   =  1         ! number of matrix to be factorized
      MKL%ERROR  =  0         ! initialize error flag
      MKL%MSGLVL =  0         ! do not print statistical information
      IF (SCARC_MKL_MTYPE == 'SYMMETRIC') THEN
         MKL%MTYPE  = -2         ! Matrix type real and symmetric indefinite
      ELSE
         MKL%MTYPE  = 11         ! Matrix type real and non-symmetric
      ENDIF
      MKL%IPARM(1)  =  1      ! supply own parameters
      MKL%IPARM(2)  =  3      ! supply own parameters
      MKL%IPARM(4)  =  0      ! supply own parameters
      MKL%IPARM(5)  =  0      ! supply own parameters
      MKL%IPARM(6)  =  0      ! write solution to x
      MKL%IPARM(8)  =  2      ! automatic setting of iterative refinement steps
      MKL%IPARM(10) = 13      ! pivoting perturbation
      MKL%IPARM(11) =  1      ! pivoting perturbation
      MKL%IPARM(13) =  1      ! pivoting perturbation
      MKL%IPARM(14) =  0      ! pivoting perturbation
      MKL%IPARM(18) =  0      ! pivoting perturbation
      MKL%IPARM(19) =  0      ! pivoting perturbation
      MKL%IPARM(20) =  0      ! pivoting perturbation
      MKL%IPARM(21) =  1      ! Bunch-Kaufman pivoting which is default in case of IPARM(0)=0
      MKL%IPARM(24) =  0      ! Bunch-Kaufman pivoting which is default in case of IPARM(0)=0
      MKL%IPARM(27) =  1      ! use matrix checker
      MKL%IPARM(40) = 2                                             ! provide matrix in distributed format
      MKL%IPARM(41) = D%NC_OFFSET(NM) + 1                      ! first global cell number for mesh NM
      MKL%IPARM(42) = D%NC_OFFSET(NM) + D%NC_LOCAL(NM)    ! last global cell number for mesh NM
      !MKL%IPARM(39) = 2                                            ! provide matrix in distributed format
      !MKL%IPARM(40) = D%NC_OFFSET(NM)+1                       ! first global cell number for mesh NM
      !MKL%IPARM(41) = D%NC_OFFSET(NM)+D%NC_LOCAL(NM)     ! last global cell number for mesh NM

#ifdef WITH_MKL_FB
      MKL%IPARM(28)=1         ! single precision

      ! perform only reordering and symbolic factorization
      WRITE(*,*) 'CHECK AC_SYM%VAL !!'
      MKL%PHASE = 11
      CALL CLUSTER_SPARSE_SOLVER_S(MKL%CT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, D%NC_GLOBAL, &
                                   AC_SYM%VAL_FB, AC_SYM%ROW, AC_SYM%COL, MKL%PERM, MKL%NRHS, MKL%IPARM, &
                                   MKL%MSGLVL, DUMMY, DUMMY, MPI_COMM_WORLD, MKL%ERROR)

      ! perform only factorization
      MKL%PHASE = 22
      CALL CLUSTER_SPARSE_SOLVER_S(MKL%CT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, D%NC_GLOBAL, &
                                   AC_SYM%VAL_FB, AC_SYM%ROW, AC_SYM%COL, MKL%PERM, MKL%NRHS, MKL%IPARM, &
                                   MKL%MSGLVL, DUMMY, DUMMY, MPI_COMM_WORLD, MKL%ERROR)

#else
      MKL%IPARM(28)=0         ! double precision

      ! perform only reordering and symbolic factorization
      MKL%PHASE = 11
      CALL CLUSTER_SPARSE_SOLVER_D(MKL%CT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, D%NC_GLOBAL, &
                                   AC_SYM%VAL, AC_SYM%ROW, AC_SYM%COL, MKL%PERM, MKL%NRHS, MKL%IPARM, &
                                   MKL%MSGLVL, DUMMY, DUMMY, MPI_COMM_WORLD, MKL%ERROR)

      ! perform only factorization
      MKL%PHASE = 22
      CALL CLUSTER_SPARSE_SOLVER_D(MKL%CT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, D%NC_GLOBAL, &
                                   AC_SYM%VAL, AC_SYM%ROW, AC_SYM%COL, MKL%PERM, MKL%NRHS, MKL%IPARM, &
                                   MKL%MSGLVL, DUMMY, DUMMY, MPI_COMM_WORLD, MKL%ERROR)

#endif

   ENDDO LEVEL_LOOP
ENDDO MESHES_LOOP

TSETUP(MYID+1)%CLUSTER=TSETUP(MYID+1)%CLUSTER+CURRENT_TIME()-TNOW
END SUBROUTINE SCARC_SETUP_CLUSTER

!> ------------------------------------------------------------------------------------------------
!> Initialize Pardiso solver
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_PARDISO(NLMIN, NLMAX)
INTEGER, INTENT(IN) :: NLMIN, NLMAX
INTEGER :: NM, NL, I, IDUMMY(1)=0
REAL (EB) :: TNOW
#ifdef WITH_MKL_FB
REAL (FB) :: DUMMY(1)=0.0_FB
#else
REAL (EB) :: DUMMY(1)=0.0_EB
#endif
TYPE (SCARC_LEVEL_TYPE), POINTER :: L=>NULL()
TYPE (SCARC_DISCRET_TYPE), POINTER :: D=>NULL()
TYPE (SCARC_MKL_TYPE), POINTER :: MKL=>NULL()
TYPE (SCARC_MATRIX_COMPACT_TYPE), POINTER :: AC_SYM=>NULL()

TNOW = CURRENT_TIME()

MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   LEVEL_LOOP: DO NL = NLMIN, NLMAX

      L   => SCARC(NM)%LEVEL(NL)
      MKL => L%MKL

      SELECT CASE(TYPE_DISCRET)
         CASE (NSCARC_DISCRET_STRUCTURED)
            D => L%SD
         CASE (NSCARC_DISCRET_UNSTRUCTURED)
            D => L%UD
      END SELECT
      AC_SYM => D%AC_SYM

      !> Allocate workspace for parameters needed in MKL-routine
      CALL SCARC_ALLOCATE_INT1(MKL%IPARM, 1, 64, NSCARC_INIT_ZERO, 'MKL_IPARM')

      !> Allocate workspace for pointers needed in MKL-routine
      IF (.NOT.ALLOCATED(MKL%PT)) THEN
         ALLOCATE(MKL%PT(64), STAT=IERROR)
         CALL CHKMEMERR ('SCARC', 'PT', IERROR)
         DO I=1,64
            MKL%PT(I)%DUMMY = 0
         ENDDO
      ENDIF

      !> Define corresponding parameters
      !> Note: IPARM-vectory is allocate from 1:64, not from 0:63
      MKL%NRHS   = 1
      MKL%MAXFCT = 1
      MKL%MNUM   = 1

      MKL%IPARM(1)  =  1      ! no solver default
      MKL%IPARM(2)  =  3      ! parallel (OpenMP) version of the nested dissection algorithm
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

      MKL%ERROR  =  0         ! initialize error flag
      MKL%MSGLVL =  0         ! do not print statistical information
      MKL%MTYPE  = -2         ! Matrix type real non-symmetric

#ifdef WITH_MKL_FB
      MKL%IPARM(28)=1         ! single precision

      ! perform only reordering and symbolic factorization
      MKL%PHASE = 11
      CALL PARDISO_S(MKL%PT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, L%NC, &
                     AC_SYM%VAL_FB, AC_SYM%ROW, AC_SYM%COL, IDUMMY, &
                     MKL%NRHS, MKL%IPARM, MKL%MSGLVL, DUMMY, DUMMY, MKL%ERROR)

      ! perform only Factorization
      MKL%PHASE = 22
      CALL PARDISO_S(MKL%PT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, L%NC, &
                     AC_SYM%VAL_FB, AC_SYM%ROW, AC_SYM%COL, IDUMMY, &
                     MKL%NRHS, MKL%IPARM, MKL%MSGLVL, DUMMY, DUMMY, MKL%ERROR)

#else
      MKL%IPARM(28)=0         ! double precision

      ! perform only reordering and symbolic factorization
      MKL%PHASE = 11
      CALL PARDISO_D(MKL%PT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, L%NC, &
                     AC_SYM%VAL, AC_SYM%ROW, AC_SYM%COL, IDUMMY, &
                     MKL%NRHS, MKL%IPARM, MKL%MSGLVL, DUMMY, DUMMY, MKL%ERROR)

      ! perform only Factorization
      MKL%PHASE = 22
      CALL PARDISO_D(MKL%PT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, L%NC, &
                     AC_SYM%VAL, AC_SYM%ROW, AC_SYM%COL, IDUMMY, &
                     MKL%NRHS, MKL%IPARM, MKL%MSGLVL, DUMMY, DUMMY, MKL%ERROR)

#endif

   ENDDO LEVEL_LOOP
ENDDO MESHES_LOOP

TSETUP(MYID+1)%PARDISO=TSETUP(MYID+1)%PARDISO+CURRENT_TIME()-TNOW
END SUBROUTINE SCARC_SETUP_PARDISO
#endif


!> ------------------------------------------------------------------------------------------------
!> Set sizes for transfer matrices
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_MATRIX_SIZES(NTYPE, NL)
INTEGER, INTENT(IN) :: NTYPE, NL
INTEGER :: IW
INTEGER :: NM, NOM
TYPE (SCARC_LEVEL_TYPE), POINTER :: L=>NULL(), OL=>NULL()
TYPE (SCARC_DISCRET_TYPE), POINTER :: D=>NULL(), OD=>NULL()
TYPE (SCARC_MATRIX_COMPACT_TYPE), POINTER :: AC=>NULL()
TYPE (SCARC_MATRIX_BANDED_TYPE), POINTER :: AB=>NULL()

SELECT CASE (NTYPE)

   !> --------------------------------------------------------------------------------------------------
   !> Define sizes for system matrix A (including extended regions related to overlaps)
   !> --------------------------------------------------------------------------------------------------
   CASE (NSCARC_SIZE_MATRIX)

      LEVEL_SYSTEM_MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

         L => SCARC(NM)%LEVEL(NL)
         SELECT CASE(TYPE_DISCRET)
            CASE (NSCARC_DISCRET_STRUCTURED)
               D => L%SD
            CASE (NSCARC_DISCRET_UNSTRUCTURED)
               D => L%UD
         END SELECT

         SELECT CASE (TYPE_MATRIX)

            !> Matrix in compact storage technique
            CASE (NSCARC_MATRIX_COMPACT)

               AC => D%AC
               IF (TWO_D) THEN
                  AC%NSTENCIL = 5
                  AC%POS(-3:3) = (/1,0,2,3,4,0,5/)     !> assignment of IOR settings to position in stencil
               ELSE
                  AC%NSTENCIL = 7
                  AC%POS(-3:3) = (/1,2,3,4,5,6,7/)
               ENDIF

               AC%NA  = D%NC * AC%NSTENCIL
               AC%NR  = D%NC + 1
WRITE(MSG%LU_DEBUG,*) 'MATRIX_SIZED: AC%NA=',AC%NA,': AC%NR=',AC%NR

               !> Determine sizes of overlapped parts for later communication with corresponding neighbors
               DO IW = 1, L%N_WALL_CELLS_EXT
                  NOM = D%WALL(IW)%NOM
                  IF (NOM /= 0) THEN
                     OL  => SCARC(NM)%OSCARC(NOM)%LEVEL(NL)
                     SELECT CASE(TYPE_DISCRET)
                        CASE (NSCARC_DISCRET_STRUCTURED)
                           OD => L%SD
                        CASE (NSCARC_DISCRET_UNSTRUCTURED)
                           OD => L%UD
                     END SELECT
                     OD%AC%NA = OD%AC%NA + AC%NSTENCIL
                  ENDIF
               ENDDO

            !> Matrix in banded storage technique
            CASE (NSCARC_MATRIX_BANDED)

               AB => D%AB

               IF (TWO_D) THEN
                  AB%NSTENCIL   = 5                      !> 5-point Laplacian
                  AB%NLO        = 2                      !> number of lower diagonals
                  AB%NUP        = 2                      !> number of upper diagonals
                  AB%POS(-3:3)  = (/5,0,4,3,2,0,1/)      !> assignment of IOR settings to columns in matrix array
                  AB%OFFSET( 3) = -L%NX                  !> lower z
                  AB%OFFSET( 1) = -1                     !> lower x
                  AB%OFFSET( 0) =  0                     !> diag
                  AB%OFFSET(-1) =  1                     !> upper x
                  AB%OFFSET(-3) =  L%NX                  !> upper z
               ELSE
                  AB%NSTENCIL   = 7                      !> 7-point Laplacian
                  AB%NLO        = 3                      !> number of lower diagonals
                  AB%NUP        = 3                      !> number of upper diagonals
                  AB%POS(-3:3)  = (/7,6,5,4,3,2,1/)      !> assignment of IOR settings to columns in matrix array
                  AB%OFFSET( 3) = -L%NX*L%NY             !> lower z
                  AB%OFFSET( 2) = -L%NX                  !> lower y
                  AB%OFFSET( 1) = -1                     !> lower x
                  AB%OFFSET( 0) =  0                     !> diag
                  AB%OFFSET(-1) =  1                     !> upper x
                  AB%OFFSET(-2) =  L%NX                  !> upper y
                  AB%OFFSET(-3) =  L%NX*L%NY             !> upper z
               ENDIF

               AB%NA    = L%NC * AB%NSTENCIL
               AB%NDIAG = L%NC

               !> Determine sizes of overlapped parts for later communication with corresponding neighbors
               DO IW = 1, L%N_WALL_CELLS_EXT
                  NOM = D%WALL(IW)%NOM
                  IF (NOM /= 0) THEN
                     OL  => SCARC(NM)%OSCARC(NOM)%LEVEL(NL)
                     SELECT CASE(TYPE_DISCRET)
                        CASE (NSCARC_DISCRET_STRUCTURED)
                           OD => L%SD
                        CASE (NSCARC_DISCRET_UNSTRUCTURED)
                           OD => L%UD
                     END SELECT
                     OD%AB%NA = OD%AB%NA + AB%NSTENCIL
                  ENDIF
               ENDDO

         END SELECT

      ENDDO LEVEL_SYSTEM_MESHES_LOOP

      !> Exchange sizes for matrices and vectors
      IF (NMESHES > 1) CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_MATRIX_SIZE, NL)

END SELECT

END SUBROUTINE SCARC_SETUP_MATRIX_SIZES


!> ----------------------------------------------------------------------------------------------------
!> Store subdivision information
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_SUBDIVISION(NL)
INTEGER, INTENT(IN) :: NL
INTEGER :: NM, IW, IOR0, IOR_LAST, NOM, INBR
INTEGER :: NEIGHBORS(20,-3:3)
TYPE (SCARC_LEVEL_TYPE), POINTER :: L=>NULL()
TYPE (SCARC_DISCRET_TYPE), POINTER :: D=>NULL()

IOR_LAST    = 0
NEIGHBORS   = 0

MESHES_LOOP1: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   L => SCARC(NM)%LEVEL(NL)
   SELECT CASE(TYPE_DISCRET)
      CASE (NSCARC_DISCRET_STRUCTURED)
         D => L%SD
      CASE (NSCARC_DISCRET_UNSTRUCTURED)
         D => L%UD
   END SELECT

   L%SUBDIVISION = 0

   WALL_CELLS_LOOP: DO IW = 1, L%N_WALL_CELLS_EXT

      IOR0 = D%WALL(IW)%IOR

      IF (IOR_LAST /= IOR0)   L%SUBDIVISION(1,IOR0) = IW
      L%SUBDIVISION(2,IOR0) = L%SUBDIVISION(2,IOR0) + 1

      NOM= D%WALL(IW)%NOM

      IF (NOM /= 0) THEN
         NEIGHBOR_LOOP: DO INBR = 1, 20
            IF (NOM == NEIGHBORS(INBR, IOR0)) THEN
               EXIT NEIGHBOR_LOOP
            ELSE IF (NEIGHBORS(INBR, IOR0) /= 0) THEN
               CYCLE NEIGHBOR_LOOP
            ELSE IF (NEIGHBORS(INBR, IOR0) == 0) THEN
               NEIGHBORS(INBR, IOR0) = NOM
               L%SUBDIVISION(3,IOR0) = L%SUBDIVISION(3,IOR0) + 1
               EXIT NEIGHBOR_LOOP
            ELSE
               CALL SCARC_SHUTDOWN(NSCARC_ERROR_NEIGHBOR_NUMBER, SCARC_NONE, NSCARC_NONE)
            ENDIF
         ENDDO NEIGHBOR_LOOP
      ENDIF

      IOR_LAST = IOR0

   ENDDO WALL_CELLS_LOOP
ENDDO MESHES_LOOP1

END SUBROUTINE SCARC_SETUP_SUBDIVISION


!> ------------------------------------------------------------------------------------------------
!> Get cell number of internal neighbor for ghost cell IW
!> ------------------------------------------------------------------------------------------------
INTEGER FUNCTION GET_WALL_CELL(WALL, IW, NX, NY)
TYPE (SCARC_WALL_TYPE), DIMENSION(:), INTENT(IN) :: WALL
INTEGER, INTENT(IN) :: IW, NX, NY
INTEGER :: IC

IF (TWO_D) THEN
   IF (ABS(WALL(IW)%NOM) == 2) THEN
      IC = -1
   ELSE
      IC = (WALL(IW)%IZW-1)*NX + WALL(IW)%IXW
   ENDIF
ELSE
   IC = (WALL(IW)%IZW-1)*NX*NY + (WALL(IW)%IYW-1)*NX + WALL(IW)%IXW
ENDIF

GET_WALL_CELL = IC
RETURN

END FUNCTION GET_WALL_CELL


!> -----------------------------------------------------------------------------------------------
!> Interface for the call of different ScaRC-solvers
!> -----------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SOLVER(DT_CURRENT)
REAL (EB), INTENT(IN) :: DT_CURRENT
REAL (EB) :: TNOW

TNOW = CURRENT_TIME()

DT  = DT_CURRENT
RDT = 1.0_EB/DT_CURRENT


ITE_PRES = ITE_PRES + 1

SELECT_METHOD: SELECT CASE (TYPE_METHOD)

   !> ---------------- McKenny-Greengard-Mayo method (MGM) --------------------
   CASE (NSCARC_METHOD_MGM)

      !> first solve inhomogeneous Poisson problem on structured grid with ScaRC (with Block-FFT)
      TYPE_DISCRET = NSCARC_DISCRET_STRUCTURED
      CALL SCARC_METHOD_CG (NSCARC_STACK_ROOT, NSCARC_STACK_ZERO, NSCARC_RHS_INHOMOGENEOUS, NLEVEL_MIN)

      CALL SCARC_VELOCITY_PREDICTOR(NLEVEL_MIN)

      !> then solve homogeneous Poisson problem on unstructured grid with UScaRC (first with SSOR-preconditioning)
      !> later the preconditioning will be replaced by an individual LU-process
      TYPE_DISCRET = NSCARC_DISCRET_UNSTRUCTURED
      CALL SCARC_METHOD_CG (NSCARC_STACK_ROOT, NSCARC_STACK_ZERO, NSCARC_RHS_HOMOGENEOUS, NLEVEL_MIN)
         

   !> ---------------- Krylov method (CG/BICG) --------------------------------
   CASE (NSCARC_METHOD_KRYLOV)

      SELECT_KRYLOV: SELECT CASE (TYPE_KRYLOV)
         CASE (NSCARC_KRYLOV_CG)
            CALL SCARC_METHOD_CG (NSCARC_STACK_ROOT, NSCARC_STACK_ZERO, NSCARC_RHS_INHOMOGENEOUS, NLEVEL_MIN)
#ifdef WITH_SCARC_CGBARO
         CASE (NSCARC_KRYLOV_CGBARO)
            CALL SCARC_METHOD_CGBARO (NSCARC_STACK_ROOT, NSCARC_STACK_ZERO, NSCARC_RHS_INHOMOGENEOUS, NLEVEL_MIN)
#endif
         CASE (NSCARC_KRYLOV_BICG)
            CALL SCARC_SHUTDOWN(NSCARC_ERROR_BICG_DISABLED, SCARC_KRYLOV, NSCARC_NONE)
!            CALL SCARC_METHOD_BICG(NSCARC_STACK_ROOT, NSCARC_STACK_ZERO, NSCARC_RHS_INHOMOGENEOUS, NLEVEL_MIN)
      END SELECT SELECT_KRYLOV


   !> ---------------- Multigrid method ---------------------------------------
   CASE (NSCARC_METHOD_MULTIGRID)

      CALL SCARC_METHOD_MULTIGRID(NSCARC_STACK_ROOT, NSCARC_STACK_ZERO, NSCARC_RHS_INHOMOGENEOUS, NLEVEL_MIN)


   !> ---------------- MKL method ---------------------------------------------
#ifdef WITH_MKL
   CASE (NSCARC_METHOD_LU)

      SELECT_MKL: SELECT CASE (TYPE_MKL)
         CASE (NSCARC_MKL_GLOBAL)
            CALL SCARC_METHOD_CLUSTER(NSCARC_STACK_ROOT, NSCARC_STACK_ZERO, NLEVEL_MIN)
         CASE (NSCARC_MKL_LOCAL)
            CALL SCARC_METHOD_PARDISO(NSCARC_STACK_ROOT, NSCARC_STACK_ZERO, NLEVEL_MIN)
      END SELECT SELECT_MKL
#endif

END SELECT SELECT_METHOD

IF (STOP_STATUS==SETUP_STOP) RETURN

#ifdef WITH_SCARC_DEBUG
CALL SCARC_SHUTDOWN(NSCARC_ERROR_POC_STOP, SCARC_NONE, NSCARC_NONE)
#endif

T_USED(5)=T_USED(5)+CURRENT_TIME()-TNOW
TSTEP(MYID+1)%SOLVER=MAX(TSTEP(MYID+1)%SOLVER,CURRENT_TIME()-TNOW)
TSUM(MYID+1)%SOLVER =TSUM(MYID+1)%SOLVER+CURRENT_TIME()-TNOW
END SUBROUTINE SCARC_SOLVER

                                     
!> ------------------------------------------------------------------------------------------------
!> Compute global matrix-vector product on grid level NL  
!>                     Y := A*X         
!> where NV1 is a reference to X and NV2 is a reference to Y
!> including data exchange along internal boundaries
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_MATVEC_PRODUCT(NV1, NV2, NS, NL)
INTEGER, INTENT(IN) :: NV1, NV2, NS, NL           !> references to vector1, vector2, current stack position and level
REAL(EB) :: TNOW
INTEGER :: NM, IC, JC, IOR0, ICOL
REAL(EB), POINTER, DIMENSION(:) :: V1, V2 
TYPE (SCARC_LEVEL_TYPE), POINTER :: L=>NULL()
TYPE (SCARC_MATRIX_COMPACT_TYPE), POINTER :: AC=>NULL()
TYPE (SCARC_MATRIX_BANDED_TYPE), POINTER :: AB=>NULL()
TYPE (SCARC_DISCRET_TYPE), POINTER :: D=>NULL()

TNOW = CURRENT_TIME()

#ifdef WITH_SCARC_DEBUG
!CALL SCARC_DEBUG_LEVEL (NV1, 'A: MATVEC 1 ', NL)
!CALL SCARC_DEBUG_LEVEL (NV2, 'A: MATVEC 2 ', NL)
#endif

!>
!> If this call is related to a globally acting solver, exchange internal boundary values of 
!> vector1 such that the ghost values contain the corresponding overlapped values of adjacent neighbor
!>
IF (STACK(NS)%SOLVER%TYPE_SCOPE == NSCARC_SCOPE_GLOBAL) THEN
   TYPE_VECTOR = NV1
   CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_VECTOR, NL)
ENDIF

#ifdef WITH_SCARC_DEBUG
!CALL SCARC_DEBUG_LEVEL (NV1, 'B: MATVEC 1 ', NL)
!CALL SCARC_DEBUG_LEVEL (NV2, 'B: MATVEC 2 ', NL)
#endif

!>
!> Perform global matrix-vector product:
!> Note: - matrix already contains subdiagonal values from neighbor along internal boundaries
!>       - if vector1 contains neighboring values, then correct values of global matvec are achieved
!>
SELECT CASE (TYPE_MATRIX)

   !>
   !> ------------- COMPACT storage technique
   !>
   CASE (NSCARC_MATRIX_COMPACT)

      MESHES_COMPACT_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

         L => SCARC(NM)%LEVEL(NL)                  

         SELECT CASE(TYPE_DISCRET)
            CASE (NSCARC_DISCRET_STRUCTURED)
               D => L%SD
            CASE (NSCARC_DISCRET_UNSTRUCTURED)
               D => L%UD
         END SELECT
         AC => D%AC

         V1 => POINT_TO_VECTOR (NM, NL, NV1)
         V2 => POINT_TO_VECTOR (NM, NL, NV2)

         DO IC = 1, D%NC

            ICOL = AC%ROW(IC)                                                !> diagonal entry
            JC   = AC%COL(ICOL)
            V2(IC) = AC%VAL(ICOL)* V1(JC)

            DO ICOL = AC%ROW(IC)+1, AC%ROW(IC+1)-1                            !> subdiagonal entries
               JC = AC%COL(ICOL)
               V2(IC) =  V2(IC) + AC%VAL(ICOL)* V1(JC)
!               WRITE(MSG%LU_DEBUG,'(A,3i4,3e14.6)') 'IC, ICOL, JC, V1(IC), V2(IC), A(ICOL):', &
!                                IC, ICOL, JC, V1(IC), V2(IC), AC%VAL(ICOL)

            ENDDO
         ENDDO

      ENDDO MESHES_COMPACT_LOOP

   !>
   !> ------------- BANDED storage technique
   !>
   CASE (NSCARC_MATRIX_BANDED)

      MESHES_BANDED_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

         L  => SCARC(NM)%LEVEL(NL)                  
         SELECT CASE(TYPE_DISCRET)
            CASE (NSCARC_DISCRET_STRUCTURED)
               D => L%SD
            CASE (NSCARC_DISCRET_UNSTRUCTURED)
               D => L%UD
         END SELECT
         AB => D%AB

         V1 => POINT_TO_VECTOR (NM, NL, NV1)               !> point to X-vector
         V2 => POINT_TO_VECTOR (NM, NL, NV2)               !> point to Y-vector

         DO IC = 1, D%NC
            V2(IC) = V1(IC)*AB%VAL(AB%POS(0),IC)                       !> main-diagonal contribution
            DO IOR0 = 3, 1, -1                                         !> lower sub-diagonal contributions
               IF (AB%POS(IOR0) == 0) CYCLE                            !> no contribution for y-direction in 2D
               JC = IC + AB%OFFSET(IOR0)         
               IF (JC >= 1) V2(IC) = V2(IC) + V1(JC)*AB%VAL(AB%POS(IOR0),IC)           
            ENDDO
            DO IOR0 = -1, -3, -1                                       !> upper sub-diagonal contributions    
               IF (AB%POS(IOR0) == 0) CYCLE                            !> no contribution for y-direction in 2D
               JC = IC + AB%OFFSET(IOR0)         
               IF (JC <= L%NC) V2(IC) = V2(IC) + V1(JC)*AB%VAL(AB%POS(IOR0),IC)           
            ENDDO
         ENDDO

      ENDDO MESHES_BANDED_LOOP

END SELECT

#ifdef WITH_SCARC_DEBUG
!CALL SCARC_DEBUG_LEVEL (NV1, 'C: MATVEC 1 ', NL)
!CALL SCARC_DEBUG_LEVEL (NV2, 'C: MATVEC 2 ', NL)
#endif

TSTEP(MYID+1)%MATVEC=MAX(TSTEP(MYID+1)%MATVEC,CURRENT_TIME()-TNOW)
TSUM(MYID+1)%MATVEC =TSUM(MYID+1)%MATVEC+CURRENT_TIME()-TNOW
END SUBROUTINE SCARC_MATVEC_PRODUCT


!> ------------------------------------------------------------------------------------------------
!> Compute global scalarproductt (including global data exchange)
!> ------------------------------------------------------------------------------------------------
REAL(EB) FUNCTION SCARC_SCALAR_PRODUCT(NV1, NV2, NL)
INTEGER, INTENT(IN) :: NV1, NV2, NL
REAL(EB) :: TNOW
INTEGER :: NM
#ifdef WITH_MKL
REAL(EB) :: DDOT
EXTERNAL :: DDOT
#else
INTEGER :: IC
#endif
REAL(EB), DIMENSION(:), POINTER ::  V1, V2
TYPE (SCARC_LEVEL_TYPE), POINTER :: L=>NULL()

TNOW = CURRENT_TIME()
LOCAL_REAL = 0.0_EB

!> Compute local scalar products on single meshes und process group
DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   L => SCARC(NM)%LEVEL(NL)

   V1 => POINT_TO_VECTOR (NM, NL, NV1)
   V2 => POINT_TO_VECTOR (NM, NL, NV2)

#ifdef WITH_MKL
   LOCAL_REAL(NM) = DDOT(L%NC, V1, 1, V2, 1)
#else
   LOCAL_REAL(NM) = 0.0_EB
   DO IC = 1, L%NC
      LOCAL_REAL(NM) = LOCAL_REAL(NM) + V1(IC) * V2(IC)
   ENDDO
#endif

ENDDO

!> Compute global scalar product as sum of local scalar products
GLOBAL_REAL = 0.0_EB

IF (N_MPI_PROCESSES>1) THEN
   CALL MPI_ALLGATHERV(MPI_IN_PLACE,1,MPI_DOUBLE_PRECISION,LOCAL_REAL,COUNTS,DISPLS,&
                       MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,IERROR)
ENDIF
GLOBAL_REAL = SUM(LOCAL_REAL(1:NMESHES))

SCARC_SCALAR_PRODUCT = GLOBAL_REAL

TSTEP(MYID+1)%SCALPROD=MAX(TSTEP(MYID+1)%SCALPROD,CURRENT_TIME()-TNOW)
TSUM(MYID+1)%SCALPROD =MAX(TSUM(MYID+1)%SCALPROD,CURRENT_TIME()-TNOW)

RETURN
END FUNCTION SCARC_SCALAR_PRODUCT


!> ------------------------------------------------------------------------------------------------
!> Compute global L2-norm (including global data exchange)
!> ------------------------------------------------------------------------------------------------
REAL(EB) FUNCTION SCARC_L2NORM(NV1, NL)
INTEGER, INTENT(IN) :: NV1, NL
REAL(EB) :: TNOW
TNOW = CURRENT_TIME()

GLOBAL_REAL = SCARC_SCALAR_PRODUCT(NV1, NV1, NL)
GLOBAL_REAL = SQRT (GLOBAL_REAL)

SCARC_L2NORM = GLOBAL_REAL

TSTEP(MYID+1)%L2NORM=MAX(TSTEP(MYID+1)%L2NORM,CURRENT_TIME()-TNOW)
TSUM(MYID+1)%L2NORM =TSUM(MYID+1)%L2NORM+CURRENT_TIME()-TNOW

RETURN
END FUNCTION SCARC_L2NORM


!> ------------------------------------------------------------------------------------------------
!> Compute linear combination of two vectors for banded storage technique
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_VECTOR_SUM(NV1, NV2, SCAL1, SCAL2, NL)
INTEGER, INTENT(IN) :: NV1, NV2, NL
REAL(EB), INTENT(IN) :: SCAL1, SCAL2
INTEGER :: NM
REAL(EB), DIMENSION(:), POINTER ::  V1, V2
#ifdef WITH_MKL
EXTERNAL :: DAXPBY
TYPE (SCARC_LEVEL_TYPE), POINTER :: L=>NULL()
#endif

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   V1 => POINT_TO_VECTOR(NM, NL, NV1)
   V2 => POINT_TO_VECTOR(NM, NL, NV2)
#ifdef WITH_MKL
   L => SCARC(NM)%LEVEL(NL)
   CALL DAXPBY(L%NC, SCAL1, V1, 1, SCAL2, V2, 1)
#else
   V2 = SCAL1 * V1 + SCAL2 * V2
#endif
ENDDO

END SUBROUTINE SCARC_VECTOR_SUM


!> ------------------------------------------------------------------------------------------------
!> Define vector2 to be a scaled copy of vector 1
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_VECTOR_COPY(NV1, NV2, SCAL1, NL)
INTEGER, INTENT(IN) :: NV1, NV2, NL
REAL(EB), INTENT(IN) :: SCAL1
INTEGER :: NM
REAL(EB), DIMENSION(:), POINTER ::  V1, V2
#ifdef WITH_MKL
EXTERNAL :: DCOPY, DSCAL
TYPE (SCARC_LEVEL_TYPE), POINTER :: L=>NULL()
#endif

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   V1 => POINT_TO_VECTOR(NM, NL, NV1)
   V2 => POINT_TO_VECTOR(NM, NL, NV2)

#ifdef WITH_MKL
   L => SCARC(NM)%LEVEL(NL)
   CALL DCOPY(L%NC, V1, 1, V2, 1)
   CALL DSCAL(L%NC, SCAL1, V2, 1)
#else
   V2 = SCAL1 * V1
#endif

ENDDO

END SUBROUTINE SCARC_VECTOR_COPY


!> ------------------------------------------------------------------------------------------------
!> Clear vector
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_VECTOR_CLEAR(NV, NL)
INTEGER, INTENT(IN) :: NV, NL
INTEGER :: NM
REAL(EB), DIMENSION(:), POINTER ::  VC

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   VC => POINT_TO_VECTOR(NM, NL, NV)
   VC =  0.0_EB
ENDDO

END SUBROUTINE SCARC_VECTOR_CLEAR


!> ------------------------------------------------------------------------------------------------
!> Preset vector with specified value
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_VECTOR_INIT (NV, VAL, NL)
INTEGER, INTENT(IN) :: NV, NL
REAL (EB), INTENT(IN) :: VAL
INTEGER :: IC, NM, I, J, K
REAL (EB), POINTER, DIMENSION(:) :: VC
TYPE (SCARC_LEVEL_TYPE), POINTER :: L=>NULL()
TYPE (SCARC_DISCRET_TYPE), POINTER :: D=>NULL()

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   L => SCARC(NM)%LEVEL(NL)
   SELECT CASE(TYPE_DISCRET)
      CASE (NSCARC_DISCRET_STRUCTURED)
         D => L%SD
      CASE (NSCARC_DISCRET_UNSTRUCTURED)
         D => L%UD
   END SELECT
   VC => POINT_TO_VECTOR (NM, NL, NV)

   DO K = 1, L%NZ
      DO J = 1, L%NY
         DO I = 1, L%NX
            IF (TYPE_DISCRET == NSCARC_DISCRET_UNSTRUCTURED .AND. L%CELL_STATE(I, J, K) /= NSCARC_CELL_GASPHASE) CYCLE
            IC = D%CELL_NUMBER(I,J,K)
            VC(IC) = VAL
         ENDDO
      ENDDO
   ENDDO

ENDDO

END SUBROUTINE SCARC_VECTOR_INIT

!> ------------------------------------------------------------------------------------------------
!> Perform preconditioning
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_RELAXATION (NV1, NV2, NS, NP, NL)
USE POIS, ONLY: H2CZSS, H3CZSS
INTEGER, INTENT(IN) :: NV1, NV2, NS, NP, NL
INTEGER :: NM, IC, JC, IW, I, J, K, ICOL, ITYPE, IDIAG
INTEGER :: IXW, IYW, IZW, ICW, IOR0
REAL(EB) :: AUX, OMEGA_SSOR = 1.5_EB, VAL
REAL(EB), DIMENSION(:), POINTER ::  V1, V2, LU
#ifdef WITH_MKL_FB
REAL(FB), DIMENSION(:), POINTER ::  V1_FB, V2_FB
#endif
TYPE (SCARC_LEVEL_TYPE), POINTER :: L=>NULL()
TYPE (SCARC_DISCRET_TYPE), POINTER :: D=>NULL()
TYPE (SCARC_FFT_TYPE), POINTER :: FFT=>NULL()
#ifdef WITH_MKL
TYPE (SCARC_MKL_TYPE), POINTER :: MKL=>NULL()
TYPE (SCARC_MATRIX_COMPACT_TYPE), POINTER :: AC_SYM=>NULL()
#endif
TYPE (SCARC_MATRIX_COMPACT_TYPE), POINTER :: AC=>NULL()
TYPE (SCARC_MATRIX_BANDED_TYPE), POINTER :: AB=>NULL()

REAL (EB) :: TNOW

TNOW = CURRENT_TIME()

#ifdef WITH_SCARC_DEBUG
CALL SCARC_DEBUG_LEVEL (NV1, 'A: BLOCK_SOLVER 1', NL)
CALL SCARC_DEBUG_LEVEL (NV2, 'A: BLOCK_SOLVER 2', NL)
#endif

ITYPE = STACK(NS-1)%SOLVER%TYPE_RELAX
SELECT CASE (ITYPE)

   !> ----------------------------------------------------------------------------------------
   !> Preconditioning by blockwise Jacobi
   !> ----------------------------------------------------------------------------------------
   CASE (NSCARC_RELAX_JACOBI)

      JACOBI_MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

         L  => SCARC(NM)%LEVEL(NL)
         SELECT CASE(TYPE_DISCRET)
            CASE (NSCARC_DISCRET_STRUCTURED)
               D => L%SD
            CASE (NSCARC_DISCRET_UNSTRUCTURED)
               D => L%UD
         END SELECT
         V2 => POINT_TO_VECTOR(NM, NL, NV2)

         JACOBI_MATRIX_CASE: SELECT CASE (TYPE_MATRIX)

            !> ------------ matrix in COMPACT storage technique
            CASE (NSCARC_MATRIX_COMPACT)

               AC => D%AC

               DO IC = 1, L%NC
                  V2(IC) = V2(IC) / AC%VAL(AC%ROW(IC))
               ENDDO

            !> ------------ matrix in BANDED storage technique
            CASE (NSCARC_MATRIX_BANDED)

               AB => D%AB
               DO IC = 1, L%NC
                  V2(IC) = V2(IC) / AB%VAL(AB%POS(0),IC)
               ENDDO

         END SELECT JACOBI_MATRIX_CASE

      ENDDO JACOBI_MESHES_LOOP

   !> ----------------------------------------------------------------------------------------
   !> Preconditioning by blockwise SSOR
   !> ----------------------------------------------------------------------------------------
   CASE (NSCARC_RELAX_SSOR)

!WRITE(MSG%LU_DEBUG,*) 'RELAX_SSOR: OMEGA = ', OMEGA_SSOR

      SSOR_MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

         L  => SCARC(NM)%LEVEL(NL)
         SELECT CASE(TYPE_DISCRET)
            CASE (NSCARC_DISCRET_STRUCTURED)
               D => L%SD
            CASE (NSCARC_DISCRET_UNSTRUCTURED)
               D => L%UD
         END SELECT
         V2 => POINT_TO_VECTOR(NM, NL, NV2)

         SSOR_MATRIX_CASE: SELECT CASE (TYPE_MATRIX)

            !> ------------ matrix in COMPACT storage technique
            CASE (NSCARC_MATRIX_COMPACT)

               AC => D%AC

               SSOR_FORWARD_COMPACT_LOOP: DO IC = 1, L%NC                        !> forward SSOR step
                  AUX = 0.0_EB
                  DO ICOL = AC%ROW(IC)+1, AC%ROW(IC+1)-1                  
                     IF (AC%COL(ICOL) >= IC) CYCLE                            !> only process lower diags
                     IF (AC%COL(ICOL) <= L%NC) AUX = AUX + AC%VAL(ICOL) * V2(AC%COL(ICOL))
                  ENDDO 
                  V2(IC) = (V2(IC) - AUX * OMEGA_SSOR) / AC%VAL(AC%ROW(IC))
!WRITE(MSG%LU_DEBUG,*) 'SSOR1: ',IC, AUX, AC%VAL(AC%ROW(IC)), V2(IC)
               ENDDO SSOR_FORWARD_COMPACT_LOOP
      
               SSOR_BACKWARD_COMPACT_LOOP: DO IC = L%NC-1, 1, -1                 !> backward SSOR step
                  AUX = 0.0_EB
                  DO ICOL = AC%ROW(IC)+1, AC%ROW(IC+1)-1
                     IF (AC%COL(ICOL) <= IC) CYCLE                            !> only process upper diags
                     IF (AC%COL(ICOL) <= L%NC) THEN
                       AUX = AUX + AC%VAL(ICOL) * V2(AC%COL(ICOL))
!WRITE(MSG%LU_DEBUG,*) 'SSOR2: ',IC, AUX, AC%VAL(AC%ROW(IC)), V2(IC), AC%COL(ICOL)
                     ENDIF
                  ENDDO 
                  V2(IC) = V2(IC) - AUX * OMEGA_SSOR / AC%VAL(AC%ROW(IC))
!WRITE(MSG%LU_DEBUG,*) 'SSOR3: ',IC, AUX, AC%VAL(AC%ROW(IC)), V2(IC)
               ENDDO SSOR_BACKWARD_COMPACT_LOOP


            !> ------------ matrix in BANDED storage technique
            CASE (NSCARC_MATRIX_BANDED)

               AB => D%AB

               !> ---------- 2D version
               IF (TWO_D) THEN

                  SSOR_FORWARD_BANDED_2D_LOOP: DO IC = 1, L%NC                 !> forward SSOR step
                     AUX = 0.0_EB
                     DO IOR0 = 1, 3, 2                                          !> only process lower x- and z-diag
                        JC = IC + AB%OFFSET(IOR0)         
                        IF (JC >= 1 .AND. JC <= L%NC) AUX = AUX + AB%VAL(AB%POS(IOR0),IC) * V2(JC)
                     ENDDO 
                     V2(IC) = (V2(IC) - AUX * OMEGA_SSOR) / AB%VAL(AB%POS(0),IC)
!WRITE(MSG%LU_DEBUG,*) 'SSOR1: ',IC, AUX, AB%VAL(AB%POS(0),IC), V2(IC)
                  ENDDO SSOR_FORWARD_BANDED_2D_LOOP
      
                  SSOR_BACKWARD_BANDED_2D_LOOP: DO IC = L%NC-1, 1, -1          !> backward SSOR step
                     AUX = 0.0_EB
                     DO IOR0 = -1, -3, -2                                       !> only process upper x- and z-diag
                       JC = IC + AB%OFFSET(IOR0)         
                        IF (JC <= L%NC) THEN
                           AUX = AUX + AB%VAL(AB%POS(IOR0),IC) * V2(JC)
!WRITE(MSG%LU_DEBUG,*) 'SSOR2: ',IC, AUX, AB%VAL(AB%POS(0),IC), V2(IC), JC
                        ENDIF
                     ENDDO 
                     V2(IC) = V2(IC) - AUX * OMEGA_SSOR / AB%VAL(AB%POS(0),IC)
!WRITE(MSG%LU_DEBUG,*) 'SSOR3: ',IC, AUX, AB%VAL(AB%POS(0),IC), V2(IC)
                  ENDDO SSOR_BACKWARD_BANDED_2D_LOOP

               !> ---------- 3D version
               ELSE

                  SSOR_FORWARD_BANDED_3D_LOOP: DO IC = 1, L%NC                    !> forward SSOR step
                     AUX = 0.0_EB
                     DO IOR0 = 1, 3                                             !> only process lower diags
                        IF (AB%POS(IOR0) == 0) CYCLE                            !> no contribution for y-direction in 2D
                        JC = IC + AB%OFFSET(IOR0)         
                        IF (JC >= 1 .AND. JC <= L%NC) AUX = AUX + AB%VAL(AB%POS(IOR0),IC) * V2(JC)
                     ENDDO 
                     V2(IC) = (V2(IC) - AUX * OMEGA_SSOR) / AB%VAL(AB%POS(0),IC)
!WRITE(MSG%LU_DEBUG,*) 'SSOR1: ',IC, AUX, AB%VAL(AB%POS(0),IC), V2(IC)
                  ENDDO SSOR_FORWARD_BANDED_3D_LOOP
      
                  SSOR_BACKWARD_BANDED_3D_LOOP: DO IC = L%NC-1, 1, -1             !> backward SSOR step
                     AUX = 0.0_EB
                     DO IOR0 = -1, -3, -1                                       !> only process upper diags
                        IF (AB%POS(IOR0) == 0) CYCLE                            !> no contribution for y-direction in 2D
                        JC = IC + AB%OFFSET(IOR0)         
!WRITE(MSG%LU_DEBUG,*) 'SSOR2: Testing ',IC, JC
!WRITE(MSG%LU_DEBUG,*) 'SSOR2: yes'
                        IF (JC >= IC .AND. JC <= L%NC) AUX = AUX + AB%VAL(AB%POS(IOR0),IC) * V2(JC)
!WRITE(MSG%LU_DEBUG,*) 'SSOR2: ',IC, AUX, AB%VAL(AB%POS(0),IC), V2(IC), JC
                     ENDDO 
                     V2(IC) = V2(IC) - AUX * OMEGA_SSOR / AB%VAL(AB%POS(0),IC)
!WRITE(MSG%LU_DEBUG,*) 'SSOR3: ',IC, AUX, AB%VAL(AB%POS(0),IC), V2(IC)
                  ENDDO SSOR_BACKWARD_BANDED_3D_LOOP

               ENDIF

         END SELECT SSOR_MATRIX_CASE

      ENDDO SSOR_MESHES_LOOP

   !> ----------------------------------------------------------------------------------------
   !> Preconditioning by different matrix-form preconditioners
   !> in all cases the preconditioner is given as separate matrix which is based
   !> on the same storage technique as the matrix AC itself;
   !> two tridiagonal systems have to be solved
   !> V1 contains the RHS to be solved for, V2 will contain the solution
   !> ----------------------------------------------------------------------------------------
   CASE (NSCARC_RELAX_GSM, NSCARC_RELAX_SGSM, NSCARC_RELAX_SORM, NSCARC_RELAX_SSORM, NSCARC_RELAX_ILU)

      LU_MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

         L  => SCARC(NM)%LEVEL(NL)
         SELECT CASE(TYPE_DISCRET)
            CASE (NSCARC_DISCRET_STRUCTURED)
               D => L%SD
            CASE (NSCARC_DISCRET_UNSTRUCTURED)
               D => L%UD
         END SELECT
         AC => D%AC

         SELECT CASE(ITYPE)
            CASE(NSCARC_RELAX_GSM)
               LU => AC%GSM
            CASE(NSCARC_RELAX_SGSM)
               LU => AC%SGSM
            CASE(NSCARC_RELAX_SORM)
               LU => AC%SORM
            CASE(NSCARC_RELAX_SSORM)
               LU => AC%SSORM
            CASE(NSCARC_RELAX_ILU)
               LU => AC%ILU
         END SELECT

         V1 => POINT_TO_VECTOR(NM, NL, NV1)
         V2 => POINT_TO_VECTOR(NM, NL, NV2)

WRITE(MSG%LU_DEBUG,*) 'LU-SOLVE, IN: V1:'
WRITE(MSG%LU_DEBUG,'(8E14.6)') V1

         !> Forward solve: 
         !> Compute sol(i) = rhs(i) - sum L(i,j) x sol(j)
         DO IC = 1, L%NC
            V2(IC) = V1(IC)
            DO ICOL = AC%ROW(IC), AC%ROW(IC+1)-1
               JC = AC%COL(ICOL)
               IF (JC >= IC) CYCLE
               V2(IC) = V2(IC) - LU(ICOL) * V2(JC)
WRITE(MSG%LU_DEBUG,*) 'A: V2(',IC,')=',V2(IC), JC
            ENDDO
         ENDDO

         !> If preconditioner is not symmetric, upper matrix U is zero and nothing has to be solved
         IF (ITYPE == NSCARC_RELAX_GSM .OR. ITYPE == NSCARC_RELAX_SORM) CYCLE

         !> Backward solve
         !> Compute sol: inv(U) sol
         DO IC = L%NC, 1, -1

            DO ICOL = AC%ROW(IC), AC%ROW(IC+1)-1
               JC = AC%COL(ICOL)
               IF (JC <= IC) CYCLE
               V2(IC) = V2(IC) - LU(ICOL) * V2(JC)
WRITE(MSG%LU_DEBUG,*) 'B: V2(',IC,')=',V2(IC), JC
            ENDDO

            !> Compute sol(i) = sol(i)/U(i,i)
            IDIAG = AC%ROW(IC)
            V2(IC) = V2(IC)/LU(IDIAG) 

WRITE(MSG%LU_DEBUG,*) 'C: V2(',IC,')=',V2(IC), IDIAG
         ENDDO

WRITE(MSG%LU_DEBUG,*) 'LU-SOLVE: OUT: V2:'
WRITE(MSG%LU_DEBUG,'(8E14.6)') V2

      ENDDO LU_MESHES_LOOP


   !> ----------------------------------------------------------------------------------------
   !> Preconditioning by blockwise Geometric Multigrid
   !> ----------------------------------------------------------------------------------------
   CASE (NSCARC_RELAX_GMG)

      CALL SCARC_METHOD_MULTIGRID (NS, NP, NSCARC_RHS_DEFECT, NLEVEL_MIN)


   !> ----------------------------------------------------------------------------------------
   !> Preconditioning by blockwise FFT based on Crayfishpak
   !> ----------------------------------------------------------------------------------------
   CASE (NSCARC_RELAX_FFT)

      TYPE_VECTOR = NV1
      CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_VECTOR, NL)

#ifdef WITH_SCARC_DEBUG
      CALL SCARC_DEBUG_LEVEL (NV1, 'B: FFT-PRECON 1', NL)
      CALL SCARC_DEBUG_LEVEL (NV2, 'B: FFT-PRECON 2', NL)
#endif

      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

         L   => SCARC(NM)%LEVEL(NL)
         FFT => L%FFT

         V1  => POINT_TO_VECTOR(NM, NL, NV1)
         V2  => POINT_TO_VECTOR(NM, NL, NV2)

         DO K = 1, L%NZ
            DO J = 1, L%NY
               DO I = 1, L%NX
                  IC = (K-1) * L%NX * L%NY + (J-1) * L%NX + I
                  FFT%PRHS(I, J, K) = V1(IC)
               ENDDO
            ENDDO
         ENDDO

         DO IW = 1, L%N_WALL_CELLS_EXT

            IXW = D%WALL(IW)%IXW
            IYW = D%WALL(IW)%IYW
            IZW = D%WALL(IW)%IZW

            ICW = L%SD%CELL_NUMBER(IXW,IYW,IZW)       !> only for structured discretization available

            VAL = 0.0_EB                              !> Use zero BC's 
            SELECT CASE(D%WALL(IW)%IOR)
               CASE( 1)
                  FFT%BXS(IYW,IZW) = VAL
               CASE(-1)
                  FFT%BXF(IYW,IZW) = VAL
               CASE( 2)
                  FFT%BYS(IXW,IZW) = VAL
               CASE(-2)
                  FFT%BYF(IXW,IZW) = VAL
               CASE( 3)
                  FFT%BZS(IXW,IYW) = VAL
               CASE(-3)
                  FFT%BZF(IXW,IYW) = VAL
            END SELECT

         ENDDO

WRITE(MSG%LU_DEBUG,*) 'FFT: IN: L%PRHS:'
WRITE(MSG%LU_DEBUG,'(8E14.6)') FFT%PRHS

         IF (TWO_D) THEN
            CALL H2CZSS (FFT%BXS,  FFT%BXF, FFT%BZS, FFT%BZF, L%NX+1, &
                         FFT%PRHS, FFT%POIS_PTB, FFT%SAVE1, FFT%WORK, FFT%HX)
         ELSE
            CALL H3CZSS (FFT%BXS,  FFT%BXF, FFT%BYS, FFT%BYF, FFT%BZS, FFT%BZF, L%NX+1, L%NY+1, &
                         FFT%PRHS, FFT%POIS_PTB, FFT%SAVE1, FFT%WORK, FFT%HX)
         ENDIF

         DO K = 1, L%NZ
            DO J = 1, L%NY
               DO I = 1, L%NX
                  IC = (K-1) * L%NX * L%NY + (J-1) * L%NX + I
                  V2(IC) = FFT%PRHS(I, J, K)
               ENDDO
            ENDDO
         ENDDO

WRITE(MSG%LU_DEBUG,*) 'FFT: OUT: V2:'
WRITE(MSG%LU_DEBUG,'(8E14.6)') V2

      ENDDO

#ifdef WITH_SCARC_DEBUG
CALL SCARC_DEBUG_LEVEL (NV1, 'C: FFT-PRECON 1', NL)
CALL SCARC_DEBUG_LEVEL (NV2, 'C: FFT-PRECON 2', NL)
#endif


#ifdef WITH_MKL
   !> ----------------------------------------------------------------------------------------
   !> Preconditioning by LU-decomposition
   !> ----------------------------------------------------------------------------------------
   CASE (NSCARC_RELAX_MKL)

      !> 
      !> ------------- Preconditioning by Cluster Sparse Solver from MKL
      !> 
      LU_SCOPE_IF: IF (STACK(NS)%SOLVER%TYPE_SCOPE == NSCARC_SCOPE_GLOBAL) THEN
      
         LU_SCOPE_GLOBAL_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   
            L   => SCARC(NM)%LEVEL(NL)
            MKL => L%MKL

            SELECT CASE(TYPE_DISCRET)
               CASE (NSCARC_DISCRET_STRUCTURED)
                  D => L%SD
               CASE (NSCARC_DISCRET_UNSTRUCTURED)
                  D => L%UD
            END SELECT
            AC_SYM => D%AC_SYM
   
            MKL%PHASE  = 33                            ! only solving
   
            V1 => POINT_TO_VECTOR (NM, NL, NV1)
            V2 => POINT_TO_VECTOR (NM, NL, NV2)
   
#ifdef WITH_MKL_FB
            V1_FB => POINT_TO_VECTOR_FB (NM, NL, NV1)
            V2_FB => POINT_TO_VECTOR_FB (NM, NL, NV2)
   
            V1_FB(1:L%NC) = REAL(V1(1:L%NC), FB)
            V2_FB(1:L%NC) = 0.0_FB
   
            CALL CLUSTER_SPARSE_SOLVER_S(MKL%CT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, D%NC_GLOBAL, &
                                         AC_SYM%VAL_FB, AC_SYM%ROW, AC_SYM%COL, MKL%PERM, MKL%NRHS, MKL%IPARM, &
                                         MKL%MSGLVL, V1_FB, V2_FB, MPI_COMM_WORLD, MKL%ERROR)
   
            V2(1:L%NC) = REAL(V2_FB(1:L%NC), EB)
   
#else
            CALL CLUSTER_SPARSE_SOLVER_D(MKL%CT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, D%NC_GLOBAL, &
                                         AC_SYM%VAL, AC_SYM%ROW, AC_SYM%COL, MKL%PERM, MKL%NRHS, MKL%IPARM, &
                                         MKL%MSGLVL, V1, V2, MPI_COMM_WORLD, MKL%ERROR)
#endif
            IF (MKL%ERROR /= 0) CALL SCARC_SHUTDOWN(NSCARC_ERROR_MKL_INTERNAL, SCARC_NONE, MKL%ERROR)
   
         ENDDO LU_SCOPE_GLOBAL_LOOP
   

      !> 
      !> -------------- Preconditioning by Pardiso Solver from MKL
      !> 
      ELSE LU_SCOPE_IF
   
         LU_SCOPE_LOCAL_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   
            L   => SCARC(NM)%LEVEL(NL)
            MKL => L%MKL

            SELECT CASE(TYPE_DISCRET)
               CASE (NSCARC_DISCRET_STRUCTURED)
                  D => L%SD
               CASE (NSCARC_DISCRET_UNSTRUCTURED)
                  D => L%UD
            END SELECT
            AC_SYM => D%AC_SYM
   
            MKL%PHASE  = 33                            ! only solving
   
            V1 => POINT_TO_VECTOR (NM, NL, NV1)
            V2 => POINT_TO_VECTOR (NM, NL, NV2)
   
#ifdef WITH_MKL_FB
   
            V1_FB => POINT_TO_VECTOR_FB (NM, NL, NV1)
            V2_FB => POINT_TO_VECTOR_FB (NM, NL, NV2)
   
            V1_FB(1:L%NC) = REAL(V1(1:L%NC), FB)
            V2_FB(1:L%NC) = 0.0_FB
   
            CALL PARDISO_S(MKL%PT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, L%NC, AC_SYM%VAL_FB, AC_SYM%ROW, AC_SYM%COL, &
                           MKL%PERM, MKL%NRHS, MKL%IPARM, MKL%MSGLVL, V1_FB, V2_FB, MKL%ERROR)
   
            V2(1:L%NC) = REAL(V2_FB(1:L%NC), EB)
#else
   
            V1 => POINT_TO_VECTOR (NM, NL, NV1)
            V2 => POINT_TO_VECTOR (NM, NL, NV2)
   
            CALL PARDISO_D(MKL%PT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, L%NC, AC_SYM%VAL, AC_SYM%ROW, AC_SYM%COL, &
                           MKL%PERM, MKL%NRHS, MKL%IPARM, MKL%MSGLVL, V1, V2, MKL%ERROR)
   
#endif
            IF (MKL%ERROR /= 0) CALL SCARC_SHUTDOWN(NSCARC_ERROR_MKL_INTERNAL, SCARC_NONE, MKL%ERROR)
   
         ENDDO LU_SCOPE_LOCAL_LOOP

      ENDIF LU_SCOPE_IF
   
#endif

END SELECT

#ifdef WITH_SCARC_DEBUG
CALL SCARC_DEBUG_LEVEL (NV1, 'B: BLOCK_SOLVER 1', NL)
CALL SCARC_DEBUG_LEVEL (NV2, 'B: BLOCK_SOLVER 2', NL)
#endif

TSTEP(MYID+1)%PRECON=MAX(TSTEP(MYID+1)%PRECON,CURRENT_TIME()-TNOW)
TSUM(MYID+1)%PRECON =TSUM(MYID+1)%PRECON+CURRENT_TIME()-TNOW

END SUBROUTINE SCARC_RELAXATION


#ifdef WITH_MKL
!> ------------------------------------------------------------------------------------------------
!> Perform global Pardiso-method based on MKL
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_METHOD_CLUSTER(NSTACK, NPARENT, NLEVEL)
INTEGER, INTENT(IN) :: NSTACK, NPARENT, NLEVEL
INTEGER ::  NM, NS, NP, NL
REAL (EB) :: TNOW
REAL(EB), POINTER, DIMENSION(:) :: V1, V2
#ifdef WITH_MKL_FB
REAL(FB), POINTER, DIMENSION(:) :: V2_FB, V1_FB
#endif
TYPE (SCARC_LEVEL_TYPE), POINTER :: L=>NULL()
TYPE (SCARC_DISCRET_TYPE), POINTER :: D=>NULL()
TYPE (SCARC_MKL_TYPE), POINTER :: MKL=>NULL()
TYPE (SCARC_MATRIX_COMPACT_TYPE), POINTER :: AC_SYM=>NULL()

NS = NSTACK
NP = NPARENT
NL = NLEVEL

TNOW = CURRENT_TIME()

CALL SCARC_SETUP_SOLVER(NS, NP)
CALL SCARC_SETUP_WORKSPACE(NS, NL, NSCARC_RHS_INHOMOGENEOUS)

MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   L   => SCARC(NM)%LEVEL(NL)
   SELECT CASE(TYPE_DISCRET)
      CASE (NSCARC_DISCRET_STRUCTURED)
         D => L%SD
      CASE (NSCARC_DISCRET_UNSTRUCTURED)
         D => L%UD
   END SELECT
   AC_SYM => D%AC_SYM

   V1  => POINT_TO_VECTOR (NM, NL, B)
   V2  => POINT_TO_VECTOR (NM, NL, X)

   MKL => L%MKL
   MKL%PHASE  = 33                                !> only solving

#ifdef WITH_MKL_FB

   V1_FB => POINT_TO_VECTOR_FB (NM, NL, B)
   V2_FB => POINT_TO_VECTOR_FB (NM, NL, X)

   V1_FB = REAL(V1, FB)
   V2_FB = 0.0_FB

   CALL CLUSTER_SPARSE_SOLVER_S(MKL%CT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, D%NC_GLOBAL, &
                                AC_SYM%VAL_FB, AC_SYM%ROW, AC_SYM%COL, MKL%PERM, MKL%NRHS, MKL%IPARM, &
                                MKL%MSGLVL, V1_FB, V2_FB, MPI_COMM_WORLD, MKL%ERROR)
   V2 = REAL(V2_FB, EB)

#else

   V1 => POINT_TO_VECTOR (NM, NL, B)

   CALL CLUSTER_SPARSE_SOLVER_D(MKL%CT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, D%NC_GLOBAL, &
                                AC_SYM%VAL, AC_SYM%ROW, AC_SYM%COL, MKL%PERM, MKL%NRHS, MKL%IPARM, &
                                MKL%MSGLVL, V1, V2, MPI_COMM_WORLD, MKL%ERROR)

#endif

   IF (MKL%ERROR /= 0) CALL SCARC_SHUTDOWN(NSCARC_ERROR_MKL_INTERNAL, SCARC_NONE, MKL%ERROR)

ENDDO MESHES_LOOP

TYPE_VECTOR = X
CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_VECTOR, NL)

IF (TYPE_SOLVER == NSCARC_SOLVER_MAIN) THEN
   CALL SCARC_UPDATE_PRESSURE_MAINCELLS (NLEVEL_MIN)
   CALL SCARC_UPDATE_PRESSURE_GHOSTCELLS(NLEVEL_MIN)
ENDIF

CALL SCARC_RELEASE_SOLVER(NS, NP)

TSTEP(MYID+1)%CLUSTER=MAX(TSTEP(MYID+1)%CLUSTER,CURRENT_TIME()-TNOW)
TSUM(MYID+1)%CLUSTER =TSUM(MYID+1)%CLUSTER+CURRENT_TIME()-TNOW
RETURN

END SUBROUTINE SCARC_METHOD_CLUSTER
#endif


#ifdef WITH_MKL
!> ------------------------------------------------------------------------------------------------
!> Perform global Pardiso-method based on MKL
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_METHOD_PARDISO(NSTACK, NPARENT, NLEVEL)
INTEGER, INTENT(IN) :: NSTACK, NPARENT, NLEVEL
INTEGER ::  NM, NS, NP, NL
REAL (EB) :: TNOW
REAL(EB), POINTER, DIMENSION(:) :: V1, V2
#ifdef WITH_MKL_FB
REAL(FB), POINTER, DIMENSION(:) :: V2_FB, V1_FB
#endif
TYPE (SCARC_LEVEL_TYPE), POINTER :: L=>NULL()
TYPE (SCARC_DISCRET_TYPE), POINTER :: D=>NULL()
TYPE (SCARC_MKL_TYPE), POINTER :: MKL=>NULL()
TYPE (SCARC_MATRIX_COMPACT_TYPE), POINTER :: AC_SYM=>NULL()

TNOW = CURRENT_TIME()

NS = NSTACK
NP = NPARENT
NL = NLEVEL

CALL SCARC_SETUP_SOLVER(NS, NP)
CALL SCARC_SETUP_WORKSPACE(NS, NL, NSCARC_RHS_INHOMOGENEOUS)

MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   L => SCARC(NM)%LEVEL(NL)

   SELECT CASE(TYPE_DISCRET)
      CASE (NSCARC_DISCRET_STRUCTURED)
         D => L%SD
      CASE (NSCARC_DISCRET_UNSTRUCTURED)
         D => L%UD
   END SELECT
   AC_SYM => D%AC_SYM

   V1  => POINT_TO_VECTOR (NM, NL, B)
   V2  => POINT_TO_VECTOR (NM, NL, X)

   MKL => L%MKL
   MKL%PHASE  = 33         ! only solving

#ifdef WITH_MKL_FB

   V1_FB => POINT_TO_VECTOR_FB (NM, NL, B)
   V2_FB => POINT_TO_VECTOR_FB (NM, NL, X)

   V1_FB = REAL(V1, FB)
   V2_FB = 0.0_FB

   CALL PARDISO_S(MKL%PT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, L%NC, &
                  AC_SYM%VAL_FB, AC_SYM%ROW, AC_SYM%COL, MKL%PERM, MKL%NRHS, MKL%IPARM, &
                  MKL%MSGLVL, V1_FB, V2_FB, MKL%ERROR)

   V2 = REAL(V2_FB, EB)

#else

   V1 => POINT_TO_VECTOR (NM, NL, B)
   V2 => POINT_TO_VECTOR (NM, NL, X)

   V2 = 0.0_EB

   CALL PARDISO_D(MKL%PT, MKL%MAXFCT, MKL%MNUM, MKL%MTYPE, MKL%PHASE, L%NC, &
                  AC_SYM%VAL, AC_SYM%ROW, AC_SYM%COL, MKL%PERM, MKL%NRHS, MKL%IPARM, &
                  MKL%MSGLVL, V1, V2, MKL%ERROR)

#endif

   IF (MKL%ERROR /= 0) CALL SCARC_SHUTDOWN(NSCARC_ERROR_MKL_INTERNAL, SCARC_NONE, MKL%ERROR)

ENDDO MESHES_LOOP

IF (TYPE_SOLVER == NSCARC_SOLVER_MAIN) THEN
   CALL SCARC_UPDATE_PRESSURE_MAINCELLS (NLEVEL_MIN)
   CALL SCARC_UPDATE_PRESSURE_GHOSTCELLS(NLEVEL_MIN)
ENDIF

CALL SCARC_RELEASE_SOLVER(NSTACK, NPARENT)

TSTEP(MYID+1)%PARDISO=MAX(TSTEP(MYID+1)%PARDISO,CURRENT_TIME()-TNOW)
TSUM(MYID+1)%PARDISO =TSUM(MYID+1)%PARDISO+CURRENT_TIME()-TNOW
RETURN

END SUBROUTINE SCARC_METHOD_PARDISO
#endif

!> ------------------------------------------------------------------------------------------------
!> Increase corresponding iteration count (just for visualization of convergence behavior)
!> ------------------------------------------------------------------------------------------------
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


!> ------------------------------------------------------------------------------------------------
!> Perform global CG-method based on global Possion-matrix
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_METHOD_CG(NSTACK, NPARENT, NRHS, NLEVEL)
INTEGER, INTENT(IN) :: NSTACK, NPARENT, NRHS, NLEVEL
INTEGER :: NSTATE, NS, NP, NL
REAL (EB) :: ALPHA, BETA, SIGMA, SIGMA_LAST
REAL (EB) :: TNOW                                                                                       !>

TNOW = CURRENT_TIME()
ITE_CG = 0

!> get current and parent stack position, and current level
NS = NSTACK
NP = NPARENT
NL = NLEVEL


#ifdef WITH_SCARC_VERBOSE
WRITE(MSG%LU_VERBOSE,*) '==========================================================================='
WRITE(MSG%LU_VERBOSE,*) '==================== Entering CG method on level ', NLEVEL
WRITE(MSG%LU_VERBOSE,*) '==========================================================================='
#endif
!> ------------------------------------------------------------------------------------------------
!> Initialization:
!>   - Get parameters for current scope (note: NL denotes the finest level)
!>   - Get right hand side vector and clear solution vectors
!> ------------------------------------------------------------------------------------------------
CALL SCARC_SETUP_SOLVER(NS, NP)
CALL SCARC_SETUP_WORKSPACE(NS, NL, NRHS)

#ifdef WITH_SCARC_DEBUG
CALL SCARC_DEBUG_LEVEL (X, 'X INIT0', NL)
CALL SCARC_DEBUG_LEVEL (B, 'B INIT0', NL)
#endif

IF (N_DIRIC_GLOBAL(NLEVEL_MIN) == 0) THEN
   CALL SCARC_VECTOR_INIT (X, 0.0_EB, NL)                    ! set x to zero
   CALL SCARC_FILTER_MEANVALUE(B, NL)                        ! filter out mean value of B
   CALL SCARC_SETUP_CONDENSING (B, NL, 1)                    ! setup condensed system
ENDIF

#ifdef WITH_SCARC_DEBUG
CALL SCARC_DEBUG_LEVEL (X, 'X INIT1', NL)
CALL SCARC_DEBUG_LEVEL (B, 'B INIT1', NL)
#endif

!> ------------------------------------------------------------------------------------------------
!> Compute initial residual and perform initial preconditioning
!> ------------------------------------------------------------------------------------------------
CALL SCARC_MATVEC_PRODUCT (X, R, NS, NL)                     !>  r^0 := A*x^0
CALL SCARC_VECTOR_SUM     (B, R, -1.0_EB, 1.0_EB, NL)        !>  r^0 := r^0 - b     corresponds to  A*x^0 - b

#ifdef WITH_SCARC_DEBUG
CALL SCARC_DEBUG_LEVEL (X, 'X INIT2', NL)
CALL SCARC_DEBUG_LEVEL (B, 'B INIT2', NL)
CALL SCARC_DEBUG_LEVEL (R, 'R INIT2', NL)
#endif

RES   = SCARC_L2NORM (R, NL)                                 !>  res   := ||r^0||
RESIN = RES                                                  !>  resin := res

NSTATE = SCARC_CONVERGENCE_STATE (0, NS, NL)                 !>  res < tolerance ?

IF (NSTATE /= NSCARC_STATE_CONV0) THEN                       !>  if no convergence yet, start precon
   CALL SCARC_PRECONDITIONER(NS, NS, NL)                     !>  v^0 := Precon(r^0)
   SIGMA_LAST = SCARC_SCALAR_PRODUCT(R, V, NL)               !>  sigma_last := (r^0,v^0)
   CALL SCARC_VECTOR_COPY (V, D, -1.0_EB, NL)                !>  d^0 := -v^0
ENDIF

!> ------------------------------------------------------------------------------------------------
!> Perform conjugate gradient looping
!> ------------------------------------------------------------------------------------------------
CG_LOOP: DO ITE = 1, NIT

   CALL SCARC_INCREASE_ITERATION_COUNTS(ITE)
WRITE(*,*) 'CG_METHOD , ITE_CG =',ITE_CG, TYPE_SOLVER

   CALL SCARC_MATVEC_PRODUCT (D, Y, NS, NL)                  !>  y^k := A*d^k

#ifdef WITH_SCARC_DEBUG
   CALL SCARC_DEBUG_LEVEL (D, 'V ITE1', NL)
   CALL SCARC_DEBUG_LEVEL (Y, 'Y ITE1', NL)
#endif

   ALPHA = SCARC_SCALAR_PRODUCT (D, Y, NL)                   !>  alpha := (d^k,y^k)     corresponds to   (d^k,A*d^k)
   ALPHA = SIGMA_LAST/ALPHA                                  !>  alpha := (r^k,v^k)/(d^k,A*d^k)

   CALL SCARC_VECTOR_SUM (D, X, ALPHA, 1.0_EB, NL)           !>  x^{k+1} := x^k + alpha * d^k 
   CALL SCARC_VECTOR_SUM (Y, R, ALPHA, 1.0_EB, NL)           !>  r^{k+1} := r^k + alpha * y^k   ~  r^k + alpha * A * d^k

#ifdef WITH_SCARC_DEBUG
   CALL SCARC_DEBUG_LEVEL (X, 'X ITE2', NL)
   CALL SCARC_DEBUG_LEVEL (R, 'Y ITE2', NL)
   WRITE(MSG%LU_DEBUG,*) 'SIGMA_LAST=',SIGMA_LAST
   WRITE(MSG%LU_DEBUG,*) 'ALPHA=',ALPHA
#endif

   RES = SCARC_L2NORM (R, NL)                               !>  res := ||r^{k+1}||
   NSTATE = SCARC_CONVERGENCE_STATE (0, NS, NL)             !>  res < tolerance ??
   IF (NSTATE /= NSCARC_STATE_PROCEED) EXIT CG_LOOP

   CALL SCARC_PRECONDITIONER(NS, NS, NL)                    !>  v^{k+1} := Precon(r^{k+1})

   SIGMA = SCARC_SCALAR_PRODUCT (R, V, NL)                  !>  sigma := (r^{k+1},v^{k+1})
   BETA  = SIGMA/SIGMA_LAST                                 !>  beta  := (r^{k+1},v^{k+1})/(r^k,v^k)
   SIGMA_LAST = SIGMA                                       !>  save last sigma

   CALL SCARC_VECTOR_SUM (V, D, -1.0_EB, BETA, NL)          !>  d^{k+1} := -v^{k+1} + beta * d^{k+1}

#ifdef WITH_SCARC_DEBUG
   CALL SCARC_DEBUG_LEVEL (V, 'V ITE2', NL)
   CALL SCARC_DEBUG_LEVEL (D, 'D ITE2', NL)
   WRITE(MSG%LU_DEBUG,*) 'SIGMA=',SIGMA
   WRITE(MSG%LU_DEBUG,*) 'BETA=',BETA
#endif

ENDDO CG_LOOP

!> ------------------------------------------------------------------------------------------------
!> Determine convergence rate and print corresponding information
!> In case of CG as main solver:
!>   - Transfer ScaRC solution vector X to FDS pressure vector
!>   - Set ghost cell values along external boundaries
!>   - Exchange values along internal boundaries
!> ------------------------------------------------------------------------------------------------
CALL SCARC_CONVERGENCE_RATE(NSTATE, NS, NL)

#ifdef WITH_SCARC_DEBUG
CALL SCARC_DEBUG_LEVEL (X, 'X FINAL', NL)
#endif

IF (N_DIRIC_GLOBAL(NLEVEL_MIN) == 0) THEN
   CALL SCARC_RESTORE_LAST_CELL(X, NL)
   CALL SCARC_FILTER_MEANVALUE(X, NL)
ENDIF

IF (TYPE_STAGE == NSCARC_SOLVER_MAIN) THEN
   CALL SCARC_UPDATE_PRESSURE_MAINCELLS  (NLEVEL_MIN)
   CALL SCARC_UPDATE_PRESSURE_GHOSTCELLS (NLEVEL_MIN)
ENDIF

CALL SCARC_RELEASE_SOLVER(NS, NP)

#ifdef WITH_SCARC_VERBOSE
WRITE(MSG%LU_VERBOSE,*) '==================== Leaving CG method on level ', NLEVEL
#endif

TSTEP(MYID+1)%KRYLOV=MAX(TSTEP(MYID+1)%KRYLOV,CURRENT_TIME()-TNOW)
TSUM(MYID+1)%KRYLOV =TSUM(MYID+1)%KRYLOV+CURRENT_TIME()-TNOW
1000 FORMAT(A30,' : ', A10,' = ', I10,' : ', A15,' = ', E14.5)
END SUBROUTINE SCARC_METHOD_CG


!> -----------------------------------------------------------------------------------------------
!> Preconditioning method which is based on the following input and output convention:
!>  - the residual which has to be preconditioned is passed in via vector R 
!>  - the result of preconditioning is passed out via vector V 
!>  - for several variants Y and Z are used as auxiliary vectors
!>  - in the comments: call is based on current grid level l (mostly the finest one)
!>  -                  l=1 denotes the finest  grid level NLEVEL_MIN
!>  -                  l=L denotes the coarset grid level NLEVEL_MAX
!> -----------------------------------------------------------------------------------------------
SUBROUTINE SCARC_PRECONDITIONER(NS, NP, NL)
INTEGER, INTENT(IN) :: NS, NP, NL                                 !> references to current stack, parent and level
INTEGER :: IL

SELECT_PRECON_TYPE: SELECT CASE (TYPE_TWOLEVEL)

   !> --------------------------------------------------------------------
   !> Classical one-level preconditioning
   !> --------------------------------------------------------------------
   CASE (NSCARC_TWOLEVEL_NONE)

      CALL SCARC_VECTOR_COPY (R, V, 1.0_EB, NL)                   !>  v := r
      CALL SCARC_RELAXATION (R, V, NS+1, NP, NL)                  !>  v := Relax(r)

   !> --------------------------------------------------------------------
   !> Additive two-level preconditioning
   !> --------------------------------------------------------------------
   CASE (NSCARC_TWOLEVEL_ADD)

      CALL SCARC_VECTOR_COPY (R, B, 1.0_EB, NL)                   !>  Use r^l as right hand side for preconditioner
      DO IL = NL, NLEVEL_MAX-1                                    !>  successively restrict it to coarser levels up to coarsest
         CALL SCARC_RESTRICTION (B, B, IL, IL+1)                  !>  b^{l+1} := Restriction(r^l)
      ENDDO

      CALL SCARC_METHOD_COARSE(NS+2, NS, NLEVEL_MAX)              !>  solve A^L * x^L := b^L on coarsest level
      CALL SCARC_VECTOR_COPY (X, Z, 1.0_EB, NLEVEL_MAX)           !>  z^L := x^L

      DO IL = NLEVEL_MAX-1, NL, -1                                !>  successively interpolate to finer levels up to finest
         CALL SCARC_PROLONGATION(Z, Z, IL+1, IL)                  !>  z^l := Prolongation(z^{l+1})
      ENDDO

      CALL SCARC_VECTOR_COPY (R, V, 1.0_EB, NL)                   !>  v^l := r^l
      CALL SCARC_RELAXATION (R, V, NS+1, NP, NL)                  !>  v^l := Relax(r^l)
      CALL SCARC_VECTOR_SUM (Z, V, 1.0_EB, 1.0_EB, NL)            !>  v^l := z^l + v^l

   !> --------------------------------------------------------------------
   !> Multiplicative two-level preconditioning (coarse first, fine second)
   !> --------------------------------------------------------------------
   CASE (NSCARC_TWOLEVEL_MUL)

      CALL SCARC_VECTOR_COPY (R, B, 1.0_EB, NL)                   !>  Use r^l as right hand side for preconditioner

      DO IL = NL, NLEVEL_MAX-1
         CALL SCARC_RESTRICTION (B, B, IL, IL+1)                  !>  b^{l+1} := Restriction(r^l)
      ENDDO

      CALL SCARC_METHOD_COARSE(NS+2, NS, NLEVEL_MAX)              !>  solve A^L * x^L := b^L on coarsest level
      CALL SCARC_VECTOR_COPY (X, Y, 1.0_EB, NLEVEL_MAX)           !>  y^L := x^L

      DO IL = NLEVEL_MAX-1, NL, -1
         CALL SCARC_PROLONGATION (Y, Y, NL+1, NL)                 !>  y^l := Prolongation(y^{l+1})
      ENDDO
      CALL SCARC_MATVEC_PRODUCT (Y, Z, NS, NL)                    !>  z^l := A^l * y^l

      CALL SCARC_VECTOR_SUM (R, Z, 1.0_EB, -1.0_EB, NL)           !>  z^l := r^l - z^l
      CALL SCARC_VECTOR_COPY (Z, V, 1.0_EB, NL)                   !>  v^l := z^l
      CALL SCARC_RELAXATION (Z, V, NS+1, NP, NL)                  !>  v^l := Relax(z^l)
      CALL SCARC_VECTOR_SUM (Y, V, 1.0_EB, 1.0_EB, NL)            !>  v^l := y^l - z^l

   !> --------------------------------------------------------------------
   !> Multiplicative two-level preconditioning (fine first, coarse second),
   !> coarse level is one level away from finest one (one coarsening step)
   !> --------------------------------------------------------------------
   CASE (NSCARC_TWOLEVEL_MUL2)

      WRITE(*,*) 'TOFIX: where is restriction?'

      CALL SCARC_VECTOR_COPY (R, V, 1.0_EB, NL)                   !>  v^l := r^l
      CALL SCARC_RELAXATION (R, V, NS+1, NP, NL)                  !>  v^l := Relax(r^l)
      CALL SCARC_MATVEC_PRODUCT (V, Z, NS, NL)                    !>  z^l := A^{l} * v^l

      CALL SCARC_VECTOR_SUM (R, Z, 1.0_EB, -1.0_EB, NL)           !>  z^l := r^l - z^l

      CALL SCARC_RESTRICTION (Z, B, NL, NL+1)                     !>  b^{l+1} := rest(R^{l})
      CALL SCARC_METHOD_COARSE(NS+2, NS, NLEVEL_MAX)              !>  x^{l+1} := A^{l+1}^{-1}(b^{l+1})
      CALL SCARC_PROLONGATION (X, Z, NL+1, NL)                    !>  v^l := Prolongation(x^{l+1})
      CALL SCARC_VECTOR_SUM (Z, V, 1.0_EB, 1.0_EB, NL)            !>  z^l := r^l - z^l


   !> --------------------------------------------------------------------
   !> Only coarse grid preconditioner
   !> --------------------------------------------------------------------
   CASE (NSCARC_TWOLEVEL_COARSE)

      CALL SCARC_VECTOR_COPY (R, B, 1.0_EB, NL)                   !>  Use r^l as right hand side for preconditioner
      DO IL = NL, NLEVEL_MAX-1                                    !>  successively restrict it to coarser levels up to coarsest
         CALL SCARC_RESTRICTION (B, B, IL, IL+1)                  !>  b^{l+1} := Restriction(b^l)
      ENDDO

      CALL SCARC_METHOD_COARSE(NS+2, NS, NLEVEL_MAX)              !>  solve A^L * x^L := b^L on coarsest level
      CALL SCARC_VECTOR_COPY (X, Y, 1.0_EB, NLEVEL_MAX)           !>  y^L := x^L

      DO IL = NLEVEL_MAX-1, NL, -1                                !>  successively interpolate to finer levels up to finest
         CALL SCARC_PROLONGATION (Y, Y, NL+1, NL)                 !>  y^l := Prolongation(y^{l+1})
      ENDDO
      CALL SCARC_VECTOR_COPY (Y, V, 1.0_EB, NL)                   !>  v^l := y^l

END SELECT SELECT_PRECON_TYPE

END SUBROUTINE SCARC_PRECONDITIONER


!> ------------------------------------------------------------------------------------------------
!> Call requested coarse grid solver (iterative/direct)
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_METHOD_COARSE(NSTACK, NPARENT, NLEVEL)
INTEGER, INTENT(IN) :: NSTACK, NPARENT, NLEVEL

#ifdef WITH_SCARC_VERBOSE
WRITE(MSG%LU_VERBOSE,*) '==================== Solving coarse grid on level ', NLEVEL_MAX
#endif

SELECT CASE (TYPE_COARSE)
   CASE (NSCARC_COARSE_ITERATIVE)
WRITE(MSG%LU_VERBOSE,*) ' --------> Solving coarse grid iteratively'
      CALL SCARC_METHOD_CG (NSTACK, NPARENT, NSCARC_RHS_DEFECT, NLEVEL)
   CASE (NSCARC_COARSE_DIRECT)
#ifdef WITH_MKL
      IF (STACK(NPARENT)%SOLVER%TYPE_SCOPE == NSCARC_SCOPE_GLOBAL .AND. N_MPI_PROCESSES > 1) THEN
WRITE(MSG%LU_VERBOSE,*) ' --------> Solving coarse grid with cluster_sparse_solver'
         CALL SCARC_METHOD_CLUSTER (NSTACK, NPARENT, NLEVEL)
      ELSE
WRITE(MSG%LU_VERBOSE,*) ' --------> Solving coarse grid with pardiso solver'
         CALL SCARC_METHOD_PARDISO (NSTACK, NPARENT, NLEVEL)
      ENDIF
#else
      WRITE(*,*) 'SCARC_METHOD_DIRECT not working yet '
#endif
END SELECT

END SUBROUTINE SCARC_METHOD_COARSE


!> ------------------------------------------------------------------------------------------------
!> Perform geometric multigrid method based on global possion-matrix
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_METHOD_MULTIGRID(NSTACK, NPARENT, NRHS, NLEVEL)
INTEGER, INTENT(IN) :: NSTACK, NPARENT, NRHS, NLEVEL
INTEGER :: NS, NP, NL
INTEGER :: NSTATE, ICYCLE
REAL (EB) :: TNOW, TNOW_COARSE

TNOW = CURRENT_TIME()
ITE_MG = 0

!> store current and parent stack position and current level
NS = NSTACK
NP = NPARENT
NL = NLEVEL

#ifdef WITH_SCARC_VERBOSE
WRITE(MSG%LU_VERBOSE,*) '------------------------------------------------------------------------------'
WRITE(MSG%LU_VERBOSE,*) '--------------------------  Entering Multigrid method on level ', NLEVEL
WRITE(MSG%LU_VERBOSE,*) '------------------------------------------------------------------------------'
#endif
!> ------------------------------------------------------------------------------------------------
!> Initialization:
!>   - Save SETTING (in case that subsequent solvers with different SETTING are called)
!>   - Define parameters for current scope (note: NL denotes the finest level)
!>   - Initialize solution, right hand side vector
!> ------------------------------------------------------------------------------------------------
CALL SCARC_SETUP_SOLVER(NS, NP)
CALL SCARC_SETUP_WORKSPACE(NS, NL, NRHS)

#ifdef WITH_SCARC_DEBUG

!> BEGIN: CAUTION CAUTION CAUTION - ONLY TEMPORARILY
!CALL SCARC_PRESET_RHS(F, NL)
!CALL SCARC_DUMP_QUANTITY(F, 'RHS', 0, NS, NL)
!CALL SCARC_PRESET_EXACT(E, NL)
!CALL SCARC_DUMP_QUANTITY(E, 'EXACT', 0, NS, NL)
!CALL SCARC_PRESET_EXACT(X, NL)
!> END: CAUTION CAUTION CAUTION - ONLY TEMPORARILY

CALL SCARC_DEBUG_LEVEL (X, 'X INIT0', NL)
CALL SCARC_DEBUG_LEVEL (B, 'B INIT0', NL)
#endif

!> ------------------------------------------------------------------------------------------------
!> Compute initial defect:  RESIN := || B - A*X ||
!>   - Initialize cycle counts for MG-iteration
!>   - Perform initial matrix-vector product on finest level
!>   - calculate norm of initial residual on finest level
!> ------------------------------------------------------------------------------------------------
!IF (TYPE_SOLVER == NSCARC_SOLVER_MAIN) THEN
   CALL SCARC_MATVEC_PRODUCT (X, V, NS, NL)                              !>  V := A*X
   CALL SCARC_VECTOR_SUM (B, V, 1.0_EB, -1.0_EB, NL)                     !>  V := B - V
   RESIN  = SCARC_L2NORM (V, NL)                                         !>  RESIN := ||V||
   NSTATE = SCARC_CONVERGENCE_STATE (0, NS, NL)                          !>  RES < TOL already ??
!ELSE
!   CALL SCARC_VECTOR_COPY(B, V, 1.0_EB, NL)                              !>  B already contains matvec
!ENDIF

#ifdef WITH_SCARC_DEBUG
CALL SCARC_DEBUG_LEVEL (V, 'V INIT0', NL)
#endif
#ifdef WITH_SCARC_VERBOSE
WRITE(MSG%LU_VERBOSE,*) '-----------> SCARC_MULTIGRID:  NIT=',NIT
#endif

ICYCLE = SCARC_CYCLING_CONTROL(NSCARC_CYCLING_SETUP, NL)

!> ------------------------------------------------------------------------------------------------
!> Perform multigrid-looping (start each iteration on finest level)
!> ------------------------------------------------------------------------------------------------
MULTIGRID_LOOP: DO ITE = 1, NIT

   CALL SCARC_INCREASE_ITERATION_COUNTS(ITE)
WRITE(*,*) 'MG-METHOD, ITE_MG =',ITE_MG, TYPE_SOLVER

   NL = NLEVEL_MIN
   ICYCLE = SCARC_CYCLING_CONTROL(NSCARC_CYCLING_RESET, NL)

   CYCLE_LOOP: DO WHILE (ICYCLE /= NSCARC_CYCLING_EXIT)

      !>
      !> Presmoothing  (smoothing/restriction till coarsest level is reached)
      !> initial and final residual are passed via vector V by default
      !>
      PRESMOOTHING_LOOP: DO WHILE (NL < NLEVEL_MAX)
         CALL SCARC_SMOOTHER (NSCARC_CYCLING_PRESMOOTH, NS+1, NS, NL)         !> D_fine   := Smooth(defect)

#ifdef WITH_SCARC_DEBUG
!CALL SCARC_PRESET_VECTOR(D, NL)
#endif
         CALL SCARC_RESTRICTION (V, B, NL, NL+1)                              !> B_coarse := Rest(D_fine)

#ifdef WITH_SCARC_DEBUG
CALL SCARC_DEBUG_LEVEL (V, 'V before restriction', NL)
CALL SCARC_DEBUG_LEVEL (B, 'B after restriction', NL+1)
#endif
         CALL SCARC_VECTOR_CLEAR (X, NL+1)                                    !> use zero initial guess on coarse level
         NL = NL + 1                                                          !> set coarser level
      ENDDO PRESMOOTHING_LOOP

      !>
      !> Coarse grid solver
      !>
      TNOW_COARSE = CURRENT_TIME()
      CALL SCARC_METHOD_COARSE(NS+2, NS, NLEVEL_MAX)                          !> X_coarse := exact_sol(.)
      TSTEP(MYID+1)%COARSE=MAX(TSTEP(MYID+1)%COARSE,CURRENT_TIME()-TNOW_COARSE)
      TSUM(MYID+1)%COARSE =TSUM(MYID+1)%COARSE+CURRENT_TIME()-TNOW_COARSE

#ifdef WITH_SCARC_DEBUG
!CALL SCARC_VECTOR_CLEAR (X, NLEVEL_MAX)                                    !> X_coarse := 0.0, only for testing
!CALL SCARC_PRESET_VECTOR(X, NLEVEL_MAX)
CALL SCARC_DEBUG_LEVEL (X, 'X Coarse', NLEVEL_MAX)
#endif

      !>
      !> Postsmoothing (smoothing/restriction till finest level is reached again)
      !>
      POSTSMOOTHING_LOOP: DO WHILE (NL > NLEVEL_MIN)
         NL=NL-1

         CALL SCARC_PROLONGATION (X, V, NL+1, NL)                             !> V_fine := Prol(X_coarse)
#ifdef WITH_SCARC_DEBUG
CALL SCARC_DEBUG_LEVEL (X, 'before prolongation', NL+1)
CALL SCARC_DEBUG_LEVEL (V, 'after prolongation', NL)
#endif
         CALL SCARC_VECTOR_SUM (V, X, 0.8_EB, 1.0_EB, NL)                     !> X_fine := V_fine + X_fine

#ifdef WITH_SCARC_DEBUG
CALL SCARC_DEBUG_LEVEL (X, 'X New Iterate', NL)
#endif
         CALL SCARC_SMOOTHER (NSCARC_CYCLING_POSTSMOOTH, NS+1, NS, NL)        !> V_fine := Smooth(defect)
         ICYCLE = SCARC_CYCLING_CONTROL(NSCARC_CYCLING_PROCEED, NL)           !> perform requested cycle
         IF (ICYCLE /= NSCARC_CYCLING_POSTSMOOTH) CYCLE CYCLE_LOOP
      ENDDO POSTSMOOTHING_LOOP

   ENDDO CYCLE_LOOP

   IF (NL /= NLEVEL_MIN) CALL SCARC_SHUTDOWN(NSCARC_ERROR_MULTIGRID_LEVEL, SCARC_NONE, NL)


 !> ---------------------------------------------------------------------------------------------
 !> Compute norm of new residual on finest level and  leave loop correspondingly
 !> ---------------------------------------------------------------------------------------------
   CALL SCARC_MATVEC_PRODUCT (X, V, NS, NL)                                   !> V := A*X
   CALL SCARC_VECTOR_SUM (B, V, 1.0_EB, -1.0_EB, NL)                          !> V := F - V

#ifdef WITH_SCARC_VERBOSE
WRITE(MSG%LU_VERBOSE,*) '--------------> SCARC_MULTIGRID: AFTER CYCLING: ', ITE
#endif

!   IF (TYPE_SOLVER == NSCARC_SOLVER_MAIN) THEN
      RES = SCARC_L2NORM (V, NL)                                              !> RES := ||V||
      NSTATE = SCARC_CONVERGENCE_STATE(0, NS, NL)                             !> convergence ?
      IF (NSTATE /= NSCARC_STATE_PROCEED) EXIT MULTIGRID_LOOP
!   ENDIF

#ifdef WITH_SCARC_DEBUG
CALL SCARC_DEBUG_LEVEL (X, 'X ITE', NL)
CALL SCARC_DEBUG_LEVEL (V, 'V ITE', NL)
WRITE(MSG%LU_DEBUG,1000) 'SCARC_MG-Iteration','ITE',ITE,'Residual=',RES
#endif

ENDDO MULTIGRID_LOOP

!> ------------------------------------------------------------------------------------------------
!> Determine convergence rate and print corresponding information
!> In case of MG as main solver:
!>   - Transfer ScaRC solution vector X to FDS pressure vector
!>   - Set ghost cell values along external boundaries
!>   - Exchange values along internal boundaries (consistency!)
!> ------------------------------------------------------------------------------------------------
!IF (TYPE_SOLVER == NSCARC_SOLVER_MAIN) CALL SCARC_CONVERGENCE_RATE(NSTATE, NS, NL)
CALL SCARC_CONVERGENCE_RATE(NSTATE, NS, NL)

SELECT CASE (TYPE_STAGE)
   CASE (NSCARC_SOLVER_MAIN)
      CALL SCARC_UPDATE_PRESSURE_MAINCELLS(NLEVEL_MIN)
      CALL SCARC_UPDATE_PRESSURE_GHOSTCELLS(NLEVEL_MIN)
   CASE (NSCARC_SOLVER_PRECON)
      CALL SCARC_UPDATE_PRECONDITIONER(NLEVEL_MIN)
END SELECT

CALL SCARC_RELEASE_SOLVER(NS, NP)

#ifdef WITH_SCARC_DEBUG
!IF (TYPE_SOLVER == NSCARC_SOLVER_MAIN) STOP
#endif

#ifdef WITH_SCARC_VERBOSE
WRITE(MSG%LU_VERBOSE,*) '==================== Leaving Multigrid method on level ', NLEVEL
#endif

TSTEP(MYID+1)%MULTIGRID=MAX(TSTEP(MYID+1)%MULTIGRID,CURRENT_TIME()-TNOW)
TSUM(MYID+1)%MULTIGRID =TSUM(MYID+1)%MULTIGRID+CURRENT_TIME()-TNOW
1000 FORMAT(A30,' : ', A10,' = ', I10,' : ', A15,' = ', E14.5)
END SUBROUTINE SCARC_METHOD_MULTIGRID


!> ------------------------------------------------------------------------------------------------
!> Control multigrid cycling (F/V/W)
!> Note: NLEVEL_MIN corresponds to finest level, NLEVEL_MAX to coarsest level
!> ------------------------------------------------------------------------------------------------
INTEGER FUNCTION SCARC_CYCLING_CONTROL(NTYPE, NL)
INTEGER, INTENT(IN) :: NTYPE, NL
INTEGER :: NM, IL, ICYCLE
TYPE (SCARC_MULTIGRID_TYPE), POINTER :: MG=>NULL()

SELECT CASE (NTYPE)

 !> ---------------------------------------------------------------------------------------------
 !> initialize cycle counts at beginning of multigrid method
 !> ---------------------------------------------------------------------------------------------
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

 !> ---------------------------------------------------------------------------------------------
 !> reset cycle counts at beginning of each new multigrid iteration
 !> ---------------------------------------------------------------------------------------------
   CASE (NSCARC_CYCLING_RESET)

      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
         DO IL = NLEVEL_MIN, NLEVEL_MAX
            MG => SCARC(NM)%LEVEL(IL)%MG
            MG%CYCLING(1)=MG%CYCLING(2)
         ENDDO
      ENDDO
      ICYCLE = NSCARC_CYCLING_NEXT

 !> ---------------------------------------------------------------------------------------------
 !> determine where to proceed with cycling
 !> ---------------------------------------------------------------------------------------------
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


!> ------------------------------------------------------------------------------------------------
!> Perform smoothing
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SMOOTHER(NTYPE, NSTACK, NPARENT, NLEVEL)
INTEGER, INTENT(IN) :: NTYPE, NSTACK, NPARENT, NLEVEL
INTEGER :: NSTATE=0, NS, NP, NL
REAL(EB) :: TNOW
LOGICAL :: BMATVEC, BL2NORM

TNOW = CURRENT_TIME()

!> get current and parent stack position and current level
NS = NSTACK
NP = NPARENT
NL = NLEVEL

CALL SCARC_SETUP_SOLVER(NS, NP)

!> 
!> Calculate initial defect on l2-norm on level NL (only if BMATVEC and Bl2NORM are set to .TRUE.)
!> Because initial vector in MG is set to zero, this defect corresponds to F
!> 
ITE = 0
BL2NORM  = .TRUE.
BMATVEC  = .TRUE.

IF (BMATVEC) THEN
   CALL SCARC_MATVEC_PRODUCT (X, V, NS, NL)                              !>  v := A*x
   CALL SCARC_VECTOR_SUM (B, V, 1.0_EB, -1.0_EB, NL)                     !>  v := b - v    corresponds to   b - A*x 
ENDIF

IF (BL2NORM) THEN
   RESIN = SCARC_L2NORM (V, NL)                                          !>  resin := ||v||
ELSE
   RESIN = SCARC_RESIDUAL
ENDIF
!NSTATE = SCARC_CONVERGENCE_STATE(NTYPE, NS, NL)

!>
!> Smoothing loop
!>
SMOOTH_LOOP: DO ITE=1, NIT

   CALL SCARC_INCREASE_ITERATION_COUNTS(ITE)
WRITE(*,*) 'SM-METHOD, ITE_SMOOTH =',ITE_SMOOTH, TYPE_SOLVER

#ifdef WITH_MKL
   IF (TYPE_SMOOTH == NSCARC_RELAX_MKL) THEN
      CALL SCARC_VECTOR_COPY(V, Z, 1.0_EB, NL)                          !>  use additional auxiliary vector Z
      CALL SCARC_RELAXATION (Z, V, NS, NP, NL)                          !>  v := Relax(z)
   ELSE
      CALL SCARC_RELAXATION (V, V, NS, NP, NL)                          !>  v := Relax(v)
   ENDIF
#else
   CALL SCARC_RELAXATION (V, V, NS, NP, NL)                             !>  v := Relax(v)
#endif

   CALL SCARC_VECTOR_SUM      (V, X, OMEGA, 1.0_EB, NL)                 !>  x := omega * v + x

   CALL SCARC_MATVEC_PRODUCT  (X, V, NS, NL)                            !>  v := A*x
   CALL SCARC_VECTOR_SUM      (B, V, 1.0_EB, -1.0_EB, NL)               !>  v := b - v

   !IF (BL2NORM.OR.ITE==NIT) THEN
   IF (BL2NORM) THEN
      RES = SCARC_L2NORM (V, NL)                                        !>  res := ||v||
      NSTATE = SCARC_CONVERGENCE_STATE(NTYPE, NS, NL)
!      IF (NSTATE /= NSCARC_STATE_PROCEED) EXIT SMOOTH_LOOP             !>  res < tol ?
   ENDIF

ENDDO SMOOTH_LOOP

CALL SCARC_RELEASE_SOLVER(NS, NP)

TSTEP(MYID+1)%SMOOTH=MAX(TSTEP(MYID+1)%SMOOTH,CURRENT_TIME()-TNOW)
TSUM(MYID+1)%SMOOTH =TSUM(MYID+1)%SMOOTH+CURRENT_TIME()-TNOW
END SUBROUTINE SCARC_SMOOTHER


!> ------------------------------------------------------------------------------------------------
!> Setup environement in every solver CALL (i.e. set pointers to used vectors) related to NSTACK
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_SOLVER(NS, NP)
INTEGER, INTENT(IN) :: NS, NP                          !> references to current stack and parent
TYPE (SCARC_SOLVER_TYPE), POINTER :: SV=>NULL(), SVP=>NULL()

SV => STACK(NS)%SOLVER

CNAME = SV%CNAME
NIT   = SV%NIT
EPS   = SV%EPS
OMEGA = SV%OMEGA
CAPPA = -1.0

TYPE_PARENT   = NP
TYPE_METHOD   = SV%TYPE_METHOD
TYPE_SOLVER   = SV%TYPE_SOLVER
TYPE_STAGE    = SV%TYPE_STAGE
TYPE_RELAX   = SV%TYPE_RELAX
TYPE_INTERPOL = SV%TYPE_INTERPOL
TYPE_TWOLEVEL = SV%TYPE_TWOLEVEL
TYPE_CYCLING  = SV%TYPE_CYCLING
TYPE_ACCURACY = SV%TYPE_ACCURACY

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

!> if not first solver in stack, restore last iteration parameters of predecessor
IF (NP > 0) THEN
   SVP => STACK(NS-1)%SOLVER
   SVP%ITE   = ITE
   SVP%RES   = RES
   SVP%RESIN = RESIN
   SVP%ERR   = ERR
   SVP%CAPPA = CAPPA
ENDIF

#ifdef SCARC_VERBOSE
WRITE(MSG%LU_VERBOSE,*)
#endif
END SUBROUTINE SCARC_SETUP_SOLVER


!> ------------------------------------------------------------------------------------------------
!> Reset SETTING of calling CURRENT-routine
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_RELEASE_SOLVER(NS, NP)
INTEGER, INTENT(IN)  :: NS, NP                            !> references to current stack and parent
TYPE (SCARC_SOLVER_TYPE), POINTER :: SV=>NULL(), SVP=>NULL()

SV  => STACK(NS)%SOLVER

!> store convergence information of preceding solver for FDS dump routine
IF (TYPE_SOLVER == NSCARC_SOLVER_MAIN) THEN
   SCARC_CAPPA      = CAPPA
   SCARC_RESIDUAL   = RES
   SCARC_ITERATIONS = ITE
   !CALL SCARC_CLOSE_CSV_FILE()
ENDIF

SV%RESIN = RESIN
SV%RES   = RES
SV%ITE   = ITE
SV%ERR   = ERR

!> if not first solver in stack, reset environment of parent (calling) routine
IF (NP > 0) THEN

   SVP => STACK(NP)%SOLVER

   ITE   = SVP%ITE
   NIT   = SVP%NIT
   EPS   = SVP%EPS
   RESIN = SVP%RESIN
   RES   = SVP%RES
   OMEGA = SVP%OMEGA
   CAPPA = SVP%CAPPA

   TYPE_PARENT   = SVP%TYPE_PARENT
   TYPE_METHOD   = SVP%TYPE_METHOD
   TYPE_SOLVER   = SVP%TYPE_SOLVER
   TYPE_STAGE    = SVP%TYPE_STAGE
   TYPE_RELAX    = SVP%TYPE_RELAX
   TYPE_INTERPOL = SVP%TYPE_INTERPOL
   TYPE_TWOLEVEL = SVP%TYPE_TWOLEVEL
   TYPE_CYCLING  = SVP%TYPE_CYCLING
   TYPE_ACCURACY = SVP%TYPE_ACCURACY

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

!> ----------------------------------------------------------------------------------------------------
!> Set initial solution corresponding to boundary data in BXS, BXF, ...
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_WORKSPACE(NS, NL, NRHS)
INTEGER, INTENT(IN) :: NS, NL, NRHS
INTEGER :: NM, IW, IOR0, I, J, K, IC
REAL(EB) :: VAL
REAL(EB), POINTER, DIMENSION(:,:,:) :: PRHS, HP
TYPE (MESH_TYPE), POINTER :: M=>NULL()
TYPE (SCARC_LEVEL_TYPE), POINTER :: L=>NULL()
TYPE (SCARC_DISCRET_TYPE), POINTER :: D=>NULL()
TYPE (SCARC_SOLVER_TYPE), POINTER :: SV=>NULL()
TYPE (SCARC_STAGE_TYPE), POINTER :: ST=>NULL(), STP=>NULL()
TYPE (SCARC_WALL_TYPE), POINTER :: WC=>NULL()
TYPE (SCARC_MGM_TYPE), POINTER :: MGM=>NULL()

SV  => STACK(NS)%SOLVER

SELECT CASE (SV%TYPE_SOLVER)

   !> --------------- If used as main solver use values from pressure-routine as initialization
   CASE (NSCARC_SOLVER_MAIN)

      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

         M  => MESHES(NM)
         L  => SCARC(NM)%LEVEL(NL)
         SELECT CASE(TYPE_DISCRET)
            CASE (NSCARC_DISCRET_STRUCTURED)
               D => L%SD
            CASE (NSCARC_DISCRET_UNSTRUCTURED)
               D => L%UD
         END SELECT
         ST => SCARC(NM)%LEVEL(NL)%STAGE(SV%TYPE_STAGE)

         PRHS => M%PRHS
         IF (PREDICTOR) THEN
            HP => M%H
         ELSE
            HP => M%HS
         ENDIF

         !> get right hand side (PRHS from pres.f90) and initial vector (H or HS from last time step)
         SELECT CASE (NRHS)

            !> Solve original problem with inhomegeneous boundary conditions
            CASE (NSCARC_RHS_INHOMOGENEOUS)

               DO K = 1, M%KBAR
                  DO J = 1, M%JBAR
                     DO I = 1, M%IBAR
                        IF (TYPE_DISCRET == NSCARC_DISCRET_UNSTRUCTURED .AND. &
                            L%CELL_STATE(I, J, K) /= NSCARC_CELL_GASPHASE) CYCLE
                        IC = D%CELL_NUMBER(I,J,K)
                        ST%B(IC) = PRHS(I, J, K)                 !> use right hand side from pres-routine
                        ST%X(IC) = HP(I, J, K)                   !> use last iterate as initial solution
                     ENDDO
                  ENDDO
               ENDDO

            !> Solve problem with homegeneous boundary conditions (MGM only)
            CASE (NSCARC_RHS_HOMOGENEOUS)

               ST%B = 0.0_EB                                    !> set RHS to zero
               ST%X = 0.0_EB                                    !> use zero as initial vector

               RDT = 1._EB/DT

               IF (TYPE_METHOD == NSCARC_METHOD_MGM) THEN
                  MGM => L%MGM
                  DO IW = L%N_WALL_CELLS_EXT+1, L%N_WALL_CELLS_EXT + L%N_WALL_CELLS_INT
                     WC => L%SD%WALL(IW)
                     IOR0 = WC%IOR
                     SELECT CASE (ABS(IOR0))
                        CASE(1)
                           ST%B(IC) = RDT * MGM%US(IW)
                        CASE(2)
                           ST%B(IC) = RDT * MGM%VS(IW)
                        CASE(3)
                           ST%B(IC) = RDT * MGM%WS(IW)
                     END SELECT
                  ENDDO
               ENDIF

         END SELECT

         LEVEL_WALL_CELLS_LOOP: DO IW = 1, M%N_EXTERNAL_WALL_CELLS

            I = D%WALL(IW)%IXW
            J = D%WALL(IW)%IYW
            K = D%WALL(IW)%IZW

            IF (TWO_D .AND. J /= 1) CALL SCARC_SHUTDOWN(NSCARC_ERROR_GRID_INDEX, SCARC_NONE, J)

            IF (TYPE_DISCRET == NSCARC_DISCRET_UNSTRUCTURED .AND. L%CELL_STATE(I, J, K) /= NSCARC_CELL_GASPHASE) CYCLE

            IOR0 = D%WALL(IW)%IOR
            IC   = D%CELL_NUMBER(I,J,K)

            !> Dirichlet BC's:
            !> these are based on the SETTING in BTYPE
            !> in the structured case this corresponds to the face-wise SETTING according to the FFT
            !> (this allows to use local FFT's as preconditioners)
            !> in the unstructured case only open boundary cells lead to Dirichlet BC's
            IF_DIRICHLET: IF (D%WALL(IW)%BTYPE == DIRICHLET) THEN

               SELECT CASE (IOR0)
                  CASE (1)
                     VAL = - 2.0_EB * L%DXI2 * M%BXS(J,K)
                  CASE (-1)
                     VAL = - 2.0_EB * L%DXI2 * M%BXF(J,K)
                  CASE (2)
                     VAL = - 2.0_EB * L%DYI2 * M%BYS(I,K)
                  CASE (-2)
                     VAL = - 2.0_EB * L%DYI2 * M%BYF(I,K)
                  CASE (3)
                     VAL = - 2.0_EB * L%DZI2 * M%BZS(I,J)
                  CASE (-3)
                     VAL = - 2.0_EB * L%DZI2 * M%BZF(I,J)
               END SELECT

               ST%B(IC) = ST%B(IC) + VAL

            ENDIF IF_DIRICHLET

            !> Neumann BC's:
            !> Note for the unstructured case only:
            !> Here, the matrix also contains Neumann BC's for those cells which have a
            !> PRESSURE_BC_INDEX == DIRICHLET but are NOT open; these cells must be excluded below,
            !> because BXS, BXF, ... contain the Dirichlet information from pres.f90 there;
            !> excluding them corresponds to a homogeneous Neumann condition for these cells
            IF_NEUMANN: IF (D%WALL(IW)%BTYPE == NEUMANN) THEN

               IF (TYPE_DISCRET == NSCARC_DISCRET_UNSTRUCTURED .AND. M%WALL(IW)%PRESSURE_BC_INDEX /= NEUMANN) CYCLE

               SELECT CASE (IOR0)
                  CASE (1)
                     VAL =   L%DXI * M%BXS(J,K)
                  CASE (-1)
                     VAL = - L%DXI * M%BXF(J,K)
                  CASE (2)
                     VAL =   L%DYI * M%BYS(I,K)
                  CASE (-2)
                     VAL = - L%DYI * M%BYF(I,K)
                  CASE (3)
                     VAL =   L%DZI * M%BZS(I,J)
                  CASE (-3)
                     VAL = - L%DZI * M%BZF(I,J)
               END SELECT

               ST%B(IC) = ST%B(IC) + VAL

            ENDIF IF_NEUMANN

         ENDDO LEVEL_WALL_CELLS_LOOP
      ENDDO

      !> In case of a Krylov method clear overlapping parts of auxiliary vectors
      IF (IS_CG.OR.HAS_TWO_LEVELS) THEN
         DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
            L  => SCARC(NM)%LEVEL(NL)
            ST => SCARC(NM)%LEVEL(NL)%STAGE(SV%TYPE_STAGE)
            ST%D(L%NC+1:L%NCE) = 0.0_EB
            ST%R(L%NC+1:L%NCE) = 0.0_EB
            ST%V(L%NC+1:L%NCE) = 0.0_EB
            ST%Y(L%NC+1:L%NCE) = 0.0_EB
            ST%Z(L%NC+1:L%NCE) = 0.0_EB
         ENDDO
      ENDIF

      !> In case of a multigrid method as main solver clear
      !> overlapping parts of auxiliary vectors and coarse grid solver vectors
      IF (IS_GMG) THEN
         DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
            L  => SCARC(NM)%LEVEL(NL)
            ST => SCARC(NM)%LEVEL(NL)%STAGE(SV%TYPE_STAGE)
            ST%V(L%NC+1:L%NCE) = 0.0_EB
            ST%Z(L%NC+1:L%NCE) = 0.0_EB
         ENDDO
      ENDIF

      !> In case of pure Neumann or periodic BCs, broadcast RHS(end) from last mesh
      !> to all and store it on all meshes
      IF (N_DIRIC_GLOBAL(NLEVEL_MIN) == 0) THEN
         IF (UPPER_MESH_INDEX == NMESHES) THEN
            L  => SCARC(NM)%LEVEL(NL)
            ST => SCARC(NM)%LEVEL(NL)%STAGE(SV%TYPE_STAGE)
            LOCAL_REAL = ST%B(L%NC)
         ELSE
            LOCAL_REAL = 0.0_EB
         ENDIF
         GLOBAL_REAL = SCARC_BROADCAST_REAL(NSCARC_BROADCAST_LAST)
         DO NM = 1, NMESHES
            SCARC(NM)%RHS_END = GLOBAL_REAL
         ENDDO
      ENDIF


   !> --------------- If MG is used as Krylov-preconditioner, vector G of main Krylov is the RHS for MG
   CASE (NSCARC_SOLVER_PRECON)

      IF (IS_CG_GMG) THEN
         DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
            L   => SCARC(NM)%LEVEL(NL)
            ST  => SCARC(NM)%LEVEL(NL)%STAGE(SV%TYPE_STAGE)
            STP => SCARC(NM)%LEVEL(NL)%STAGE(NSCARC_STAGE_ONE)
            ST%X = 0.0_EB
            ST%B = STP%R
            ST%V(L%NC+1:L%NCE) = 0.0_EB
            ST%Z(L%NC+1:L%NCE) = 0.0_EB
         ENDDO
      ENDIF

   !> --------------- If used as coarse grid solver start with zero initialization
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

END SELECT

END SUBROUTINE SCARC_SETUP_WORKSPACE

!> ------------------------------------------------------------------------------------------------
!> Check if solver converges or diverges and print out residual information
!> ------------------------------------------------------------------------------------------------
INTEGER FUNCTION SCARC_CONVERGENCE_STATE(ISM, NS, NL)
INTEGER, INTENT(IN) :: NL, NS, ISM
INTEGER :: NSTATE

NSTATE = NSCARC_STATE_PROCEED

!>
!> ONLY TEMPORARILY - START
!> Dump residual, exact and approximate solution
!>
#ifdef WITH_SCARC_DEBUG
!CALL SCARC_PRESET_EXACT(E, NL)

!call SCARC_DEBUG_LEVEL (R, 'RESIDUAL', NL)
!call SCARC_DEBUG_LEVEL (E, 'EXACT'   , NL)
!call SCARC_DEBUG_LEVEL (X, 'DISCRET' , NL)
!
!CALL SCARC_DUMP_QUANTITY(R, 'RESIDUAL', ISM, NS, NL)
!CALL SCARC_DUMP_QUANTITY(X, 'DISCRET' , ISM, NS, NL)
!
!CALL SCARC_VECTOR_SUM(X, E, 1.0_EB, -1.0_EB, NL)
!ERR = SCARC_L2NORM (E, NL)
!
!call SCARC_DEBUG_LEVEL (E, 'ERROR', NL)
!CALL SCARC_DUMP_QUANTITY(E, 'ERROR', ISM, NS, NL)

#endif
!>
!> ONLY TEMPORARILY - END
!>

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
      ENDIF
END SELECT
IF (RES > NSCARC_THRESHOLD_DIVGERGENCE) NSTATE = NSCARC_STATE_DIVG

SCARC_CONVERGENCE_STATE = NSTATE

IF (HAS_CSV) CALL SCARC_DUMP_CSV(ISM, NS, NL)

#ifdef WITH_SCARC_VERBOSE
!IF (MYID == 0) WRITE(LU_OUTPUT,1000) TRIM(CNAME), NL, ITE, RES
WRITE(MSG%LU_VERBOSE,1000) STACK(NS)%SOLVER%CNAME, NL, ITE, RES
#endif

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG, 1000) STACK(NS)%SOLVER%CNAME, NL, ITE, RES
#endif

RETURN

1000 FORMAT (A30,': Level=',i4,': Iteration = ',i4,': Residual =',e14.6)
END FUNCTION SCARC_CONVERGENCE_STATE


!> ------------------------------------------------------------------------------------------------
!> Compute convergence rate and print out residual information for final loop
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_CONVERGENCE_RATE(NSTATE, NS, NL)
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

#ifdef WITH_SCARC_VERBOSE
IF (TRIM(CNAME) /= 'SCARC_COARSE_CG') THEN
   WRITE(MSG%LU_VERBOSE,2000) STACK(NS)%SOLVER%CNAME, ITE, CAPPA
ENDIF
#endif

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,2000) STACK(NS)%SOLVER%CNAME, ITE, CAPPA
#endif

2000 FORMAT (A30,': Iterations: ',i6,':   Convergence Rate =',e14.6,/)
END SUBROUTINE SCARC_CONVERGENCE_RATE


!> ------------------------------------------------------------------------------------------------
!> Perform restriction from finer to coarser grid in multigrid method
!>    - 'VF' corresponds to vector on fine   grid
!>    - 'VC' corresponds to vector on coarse grid
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_RESTRICTION (NVB, NVC, NLF, NLC)
INTEGER, INTENT(IN) :: NVB, NVC, NLF, NLC
INTEGER :: NM
INTEGER :: IXF, IYF, IZF, ICF(8)=0, ICFB(-2:2,-2:2)=0
INTEGER :: IXC, IYC, IZC, ICC
REAL(EB), POINTER, DIMENSION(:) :: VC, VF
TYPE (SCARC_LEVEL_TYPE), POINTER :: LF=>NULL(), LC=>NULL()
TYPE (SCARC_DISCRET_TYPE), POINTER :: DF=>NULL(), DC=>NULL()

!>
!> ------------------ Twolevel-CG or Geometric multigrid (as main solver or preconditioner) --------------
!>
IF (HAS_MULTIPLE_LEVELS) THEN

   DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

      LF => SCARC(NM)%LEVEL(NLF)
      LC => SCARC(NM)%LEVEL(NLC)

      SELECT CASE(TYPE_DISCRET)
         CASE(NSCARC_DISCRET_STRUCTURED)
            DC => LC%SD 
            DF => LF%SD 
         CASE(NSCARC_DISCRET_UNSTRUCTURED)
            DC => LC%UD 
            DF => LF%UD 
      END SELECT

      VF => POINT_TO_VECTOR(NM, NLF, NVB)
      VC => POINT_TO_VECTOR(NM, NLC, NVC)

      IF (TWO_D) THEN

         SELECT_INTERPOL: SELECT CASE (TYPE_INTERPOL)

            !> ---------- Constant Interpolation
            CASE (NSCARC_INTERPOL_CONSTANT)

               DO IZC = 1, LC%NZ
                  DO IXC = 1, LC%NX

                     IF (TYPE_DISCRET == NSCARC_DISCRET_UNSTRUCTURED .AND. &
                         LC%CELL_STATE(IXC, 1, IZC) /= NSCARC_CELL_GASPHASE) CYCLE

                     IXF = 2*IXC
                     IZF = 2*IZC

                     ICC = DC%CELL_NUMBER(IXC, 1, IZC)

                     ICF(1) = DF%CELL_NUMBER(IXF-1, 1, IZF-1)
                     ICF(2) = DF%CELL_NUMBER(IXF-1, 1, IZF  )
                     ICF(3) = DF%CELL_NUMBER(IXF  , 1, IZF-1)
                     ICF(4) = DF%CELL_NUMBER(IXF  , 1, IZF  )

                     VC(ICC) = 0.25_EB * (  VF(ICF(1)) &
                                          + VF(ICF(2)) &
                                          + VF(ICF(3)) &
                                          + VF(ICF(4)) )
                  ENDDO
               ENDDO

            !> ---------- Bilinear Interpolation
            CASE (NSCARC_INTERPOL_BILINEAR)

               VC=0.0_EB

               DO IZC = 1, LC%NZ
                  DO IXC = 1, LC%NX

                     IF (TYPE_DISCRET == NSCARC_DISCRET_UNSTRUCTURED .AND. &
                         LC%CELL_STATE(IXC, 1, IZC) /= NSCARC_CELL_GASPHASE) CYCLE

                     IXF = 2*IXC
                     IZF = 2*IZC

                     ICC = DC%CELL_NUMBER(IXC, 1, IZC)

                     ICFB(-2,-2) = DF%CELL_NUMBER(IXF-2, 1, IZF-2)
                     ICFB(-1,-2) = DF%CELL_NUMBER(IXF-1, 1, IZF-2)
                     ICFB( 1,-2) = DF%CELL_NUMBER(IXF  , 1, IZF-2)
                     ICFB( 2,-2) = DF%CELL_NUMBER(IXF+1, 1, IZF-2)

                     ICFB(-2,-1) = DF%CELL_NUMBER(IXF-2, 1, IZF-1)
                     ICFB(-1,-1) = DF%CELL_NUMBER(IXF-1, 1, IZF-1)
                     ICFB( 1,-1) = DF%CELL_NUMBER(IXF  , 1, IZF-1)
                     ICFB( 2,-1) = DF%CELL_NUMBER(IXF+1, 1, IZF-1)

                     ICFB(-2, 1) = DF%CELL_NUMBER(IXF-2, 1, IZF)
                     ICFB(-1, 1) = DF%CELL_NUMBER(IXF-1, 1, IZF)
                     ICFB( 1, 1) = DF%CELL_NUMBER(IXF  , 1, IZF)
                     ICFB( 2, 1) = DF%CELL_NUMBER(IXF+1, 1, IZF)

                     ICFB(-2, 2) = DF%CELL_NUMBER(IXF-2, 1, IZF+1)
                     ICFB(-1, 2) = DF%CELL_NUMBER(IXF-1, 1, IZF+1)
                     ICFB( 1, 2) = DF%CELL_NUMBER(IXF  , 1, IZF+1)
                     ICFB( 2, 2) = DF%CELL_NUMBER(IXF+1, 1, IZF+1)

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

         END SELECT SELECT_INTERPOL

      ELSE

         !> ---------- Constant Interpolation (Note: 3D-bilinear case is still missing)
         DO IZC = 1, LC%NZ
            DO IYC = 1, LC%NY
               DO IXC = 1, LC%NX

                  IF (TYPE_DISCRET == NSCARC_DISCRET_UNSTRUCTURED .AND. &
                      LC%CELL_STATE(IXC, IYC, IZC) /= NSCARC_CELL_GASPHASE) CYCLE

                  IXF = 2*IXC
                  IYF = 2*IYC
                  IZF = 2*IZC

                  ICC = DC%CELL_NUMBER(IXC, IYC, IZC)

                  ICF(1) = DF%CELL_NUMBER(IXF-1, IYF-1, IZF-1)
                  ICF(2) = DF%CELL_NUMBER(IXF-1, IYF-1, IZF  )
                  ICF(3) = DF%CELL_NUMBER(IXF-1, IYF  , IZF-1)
                  ICF(4) = DF%CELL_NUMBER(IXF-1, IYF  , IZF  )
                  ICF(5) = DF%CELL_NUMBER(IXF  , IYF-1, IZF-1)
                  ICF(6) = DF%CELL_NUMBER(IXF  , IYF-1, IZF  )
                  ICF(7) = DF%CELL_NUMBER(IXF  , IYF  , IZF-1)
                  ICF(8) = DF%CELL_NUMBER(IXF  , IYF  , IZF  )

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

      ENDIF

   ENDDO

ENDIF

END SUBROUTINE SCARC_RESTRICTION


!> ------------------------------------------------------------------------------------------------
!> Perform prolongation from coarser to finer grid
!>    - 'VC' corresponds to coarser grid
!>    - 'VF' corresponds to finer   grid
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_PROLONGATION (NVC, NVB, NLC, NLF)
INTEGER, INTENT(IN) :: NVC, NVB, NLC, NLF
INTEGER :: NM, I
INTEGER :: IXF, IYF, IZF, ICF(8)=0, ICFB(-1:1,-1:1)=0
INTEGER :: IXC, IYC, IZC, ICC
REAL(EB), POINTER, DIMENSION(:) :: VC, VF
TYPE (SCARC_LEVEL_TYPE), POINTER :: LF=>NULL(), LC=>NULL()
TYPE (SCARC_DISCRET_TYPE), POINTER :: DF=>NULL(), DC=>NULL()

!>
!> ------------------ Twolevel CG or Geometric Multigrid -------------------------------------------
!>
IF (HAS_MULTIPLE_LEVELS) THEN

   DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

      LC => SCARC(NM)%LEVEL(NLC)
      LF => SCARC(NM)%LEVEL(NLF)

      VC => POINT_TO_VECTOR(NM, NLC, NVC)
      VF => POINT_TO_VECTOR(NM, NLF, NVB)

      SELECT CASE(TYPE_DISCRET)
         CASE(NSCARC_DISCRET_STRUCTURED)
            DC => LC%SD
            DF => LF%SD
         CASE(NSCARC_DISCRET_UNSTRUCTURED)
            DC => LC%UD
            DF => LF%UD
      END SELECT

      IF (TWO_D) THEN

         SELECT_INTERPOL: SELECT CASE (TYPE_INTERPOL)

            CASE (NSCARC_INTERPOL_CONSTANT)

               DO IZC = 1, LC%NZ
                  DO IXC = 1, LC%NX

                     IF (TYPE_DISCRET == NSCARC_DISCRET_UNSTRUCTURED .AND. &
                         LC%CELL_STATE(IXC, 1, IZC) /= NSCARC_CELL_GASPHASE) CYCLE

                     IXF = 2*IXC
                     IZF = 2*IZC

                     ICC = DC%CELL_NUMBER(IXC, 1, IZC)

                     ICF(1) = DF%CELL_NUMBER(IXF-1, 1, IZF-1)
                     ICF(2) = DF%CELL_NUMBER(IXF-1, 1, IZF  )
                     ICF(3) = DF%CELL_NUMBER(IXF  , 1, IZF-1)
                     ICF(4) = DF%CELL_NUMBER(IXF  , 1, IZF  )

                     DO I = 1, 4
                        VF(ICF(I)) = VC(ICC)
                     ENDDO
                  ENDDO
               ENDDO

            CASE (NSCARC_INTERPOL_BILINEAR)

               DO IZC = 1, LC%NZ
                  DO IXC = 1, LC%NX

                     IF (TYPE_DISCRET == NSCARC_DISCRET_UNSTRUCTURED .AND. &
                         LC%CELL_STATE(IXC, 1, IZC) /= NSCARC_CELL_GASPHASE) CYCLE

                     IXF = 2*IXC
                     IZF = 2*IZC

                     ICC = DC%CELL_NUMBER(IXC, 1, IZC)

                     ICFB(-1,-1) = DF%CELL_NUMBER(IXF-1, 1, IZF-1)
                     ICFB(-1, 1) = DF%CELL_NUMBER(IXF-1, 1, IZF  )
                     ICFB( 1,-1) = DF%CELL_NUMBER(IXF  , 1, IZF-1)
                     ICFB( 1, 1) = DF%CELL_NUMBER(IXF  , 1, IZF  )

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

         END SELECT SELECT_INTERPOL

      ELSE

         !> ---------- Constant Interpolation (Note: 3D-bilinear case is still missing)
         DO IZC = 1, LC%NZ
            DO IYC = 1, LC%NY
               DO IXC = 1, LC%NX

                  IF (TYPE_DISCRET == NSCARC_DISCRET_UNSTRUCTURED .AND. &
                      LC%CELL_STATE(IXC, IYC, IZC) /= NSCARC_CELL_GASPHASE) CYCLE

                  IXF = 2*IXC
                  IYF = 2*IYC
                  IZF = 2*IZC

                  ICC = DC%CELL_NUMBER(IXC, IYC, IZC)

                  ICF(1) = DF%CELL_NUMBER(IXF-1, IYF-1, IZF-1)
                  ICF(2) = DF%CELL_NUMBER(IXF-1, IYF-1, IZF  )
                  ICF(3) = DF%CELL_NUMBER(IXF-1, IYF  , IZF-1)
                  ICF(4) = DF%CELL_NUMBER(IXF-1, IYF  , IZF  )
                  ICF(5) = DF%CELL_NUMBER(IXF  , IYF-1, IZF-1)
                  ICF(6) = DF%CELL_NUMBER(IXF  , IYF-1, IZF  )
                  ICF(7) = DF%CELL_NUMBER(IXF  , IYF  , IZF-1)
                  ICF(8) = DF%CELL_NUMBER(IXF  , IYF  , IZF  )

                  DO I = 1, 8
                     VF(ICF(I)) = VC(ICC)
                  ENDDO

               ENDDO
            ENDDO
         ENDDO
      ENDIF
   ENDDO
ENDIF

END SUBROUTINE SCARC_PROLONGATION

!> ------------------------------------------------------------------------------------------------
!> Copy final solution from GMG (as preconditioner) to corresponding vector of CG (as main solver)
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_UPDATE_PRECONDITIONER(NL)
INTEGER, INTENT(IN) :: NL
INTEGER :: NM
DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   SCARC(NM)%LEVEL(NL)%STAGE(NSCARC_STAGE_ONE)%V = SCARC(NM)%LEVEL(NL)%STAGE(NSCARC_STAGE_TWO)%X
ENDDO
END SUBROUTINE SCARC_UPDATE_PRECONDITIONER

!> ------------------------------------------------------------------------------------------------
!> Finalize data
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_UPDATE_PRESSURE_MAINCELLS(NL)
INTEGER, INTENT(IN) :: NL
INTEGER :: NM, I, J, K, IC
REAL(EB), POINTER, DIMENSION(:,:,:) :: HP
TYPE (MESH_TYPE), POINTER :: M=>NULL()
TYPE (SCARC_LEVEL_TYPE), POINTER :: L=>NULL()
TYPE (SCARC_DISCRET_TYPE), POINTER :: D=>NULL()
TYPE (SCARC_STAGE_TYPE), POINTER :: ST=>NULL()

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   M  => MESHES(NM)
   L  => SCARC(NM)%LEVEL(NL)
   SELECT CASE(TYPE_DISCRET)
      CASE (NSCARC_DISCRET_STRUCTURED)
         D => L%SD
      CASE (NSCARC_DISCRET_UNSTRUCTURED)
         D => L%UD
   END SELECT
   ST => SCARC(NM)%LEVEL(NL)%STAGE(NSCARC_STAGE_ONE)

   IF (PREDICTOR) THEN
      HP => M%H
   ELSE
      HP => M%HS
   ENDIF

   DO K = 1, M%KBAR
      DO J = 1, M%JBAR
         DO I = 1, M%IBAR
            IF (TYPE_DISCRET == NSCARC_DISCRET_UNSTRUCTURED .AND. L%CELL_STATE(I, J, K) /= NSCARC_CELL_GASPHASE) CYCLE
            IC = D%CELL_NUMBER(I,J,K)
            HP(I, J, K) = ST%X(IC)
         ENDDO
      ENDDO
   ENDDO

ENDDO

END SUBROUTINE SCARC_UPDATE_PRESSURE_MAINCELLS


!> ------------------------------------------------------------------------------------------------
!> Set correct boundary values at external and internal boundaries
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_UPDATE_PRESSURE_GHOSTCELLS(NL)
INTEGER, INTENT(IN) :: NL
INTEGER :: NM, IW, IOR0, IXG, IYG, IZG, IXW, IYW, IZW
REAL(EB), POINTER, DIMENSION(:,:,:) :: HP=>NULL()
TYPE (MESH_TYPE), POINTER :: M=>NULL()
TYPE (SCARC_LEVEL_TYPE), POINTER :: L=>NULL()
TYPE (SCARC_DISCRET_TYPE), POINTER :: D=>NULL()
TYPE (SCARC_WALL_TYPE), POINTER :: WC=>NULL()

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   M => MESHES(NM)
   L => SCARC(NM)%LEVEL(NL)
   SELECT CASE(TYPE_DISCRET)
      CASE (NSCARC_DISCRET_STRUCTURED)
         D => L%SD
      CASE (NSCARC_DISCRET_UNSTRUCTURED)
         D => L%UD
   END SELECT

   IF (PREDICTOR) THEN
      HP => M%H
   ELSE
      HP => M%HS
   ENDIF

   !> compute ghost cell values
   WALL_CELLS_LOOP: DO IW = 1, L%N_WALL_CELLS_EXT

      WC => D%WALL(IW)

      IXG = WC%IXG
      IYG = WC%IYG
      IZG = WC%IZG

      IXW = WC%IXW
      IYW = WC%IYW
      IZW = WC%IZW

      IOR0 = WC%IOR

      SELECT CASE (IOR0)
         CASE ( 1)
            IF (WC%BTYPE==DIRICHLET) THEN
               HP(IXG,IYW,IZW) = -HP(IXW,IYW,IZW) + 2.0_EB * M%BXS(IYW,IZW)
            ELSE IF (WC%BTYPE==NEUMANN) THEN
               HP(IXG,IYW,IZW) =  HP(IXW,IYW,IZW) - L%DX *M%BXS(IYW,IZW)
            ENDIF
         CASE (-1)
            IF (WC%BTYPE==DIRICHLET) THEN
               HP(IXG,IYW,IZW) = -HP(IXW,IYW,IZW) + 2.0_EB * M%BXF(IYW,IZW)
            ELSE IF (WC%BTYPE==NEUMANN) THEN
               HP(IXG,IYW,IZW) =  HP(IXW,IYW,IZW) + L%DX *M%BXF(IYW,IZW)
            ENDIF
         CASE ( 2)
            IF (WC%BTYPE==DIRICHLET) THEN
               HP(IXW,IYG,IZW) = -HP(IXW,IYW,IZW) + 2.0_EB * M%BYS(IXW,IZW)
            ELSE IF (WC%BTYPE==NEUMANN) THEN
               HP(IXW,IYG,IZW) =  HP(IXW,IYW,IZW) - L%DY *M%BYS(IXW,IZW)
            ENDIF
         CASE (-2)
            IF (WC%BTYPE==DIRICHLET) THEN
               HP(IXW,IYG,IZW) = -HP(IXW,IYW,IZW) + 2.0_EB * M%BYF(IXW,IZW)
            ELSE IF (WC%BTYPE==NEUMANN) THEN
               HP(IXW,IYG,IZW) =  HP(IXW,IYW,IZW) + L%DY *M%BYF(IXW,IZW)
            ENDIF
         CASE ( 3)
            IF (WC%BTYPE==DIRICHLET) THEN
               HP(IXW,IYW,IZG) = -HP(IXW,IYW,IZW) + 2.0_EB * M%BZS(IXW,IYW)
            ELSE IF (WC%BTYPE==NEUMANN) THEN
               HP(IXW,IYW,IZG) =  HP(IXW,IYW,IZW) - L%DZ *M%BZS(IXW,IYW)
            ENDIF
         CASE (-3)
         IF (WC%BTYPE==DIRICHLET) THEN
               HP(IXW,IYW,IZG) = -HP(IXW,IYW,IZW) + 2.0_EB * M%BZF(IXW,IYW)
            ELSE IF (WC%BTYPE==NEUMANN) THEN
               HP(IXW,IYW,IZG) =  HP(IXW,IYW,IZW) + L%DZ *M%BZF(IXW,IYW)
            ENDIF
      END SELECT
   ENDDO WALL_CELLS_LOOP

ENDDO

!> -----------------------------------------------------------------------------------------------
!> Perform data exchange to achieve consistency of ghost values along internal boundaries
!> -----------------------------------------------------------------------------------------------
CALL SCARC_EXCHANGE(NSCARC_EXCHANGE_PRESSURE, NL)

END SUBROUTINE SCARC_UPDATE_PRESSURE_GHOSTCELLS


!> ------------------------------------------------------------------------------------------------
!>  Perform data exchange corresponding to requested exchange type (CALL receive and send-routines)
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_EXCHANGE (NTYPE, NL)
INTEGER, INTENT(IN) :: NTYPE, NL

N_REQ = 0
TYPE_EXCHANGE = NTYPE

CALL SCARC_EXCHANGE_RECEIVE(NL)
CALL SCARC_EXCHANGE_SEND(NL)

END SUBROUTINE SCARC_EXCHANGE


!> ------------------------------------------------------------------------------------------------
!>  Receive data from neighbors (corresponds to POST_RECEIVES)
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_EXCHANGE_RECEIVE (NL)
INTEGER, INTENT(IN) :: NL
INTEGER :: NM, NOM
TYPE (SCARC_LEVEL_TYPE), POINTER :: L=>NULL()
TYPE (SCARC_NEIGHBOR_TYPE), POINTER :: OS=>NULL()

RECEIVE_MESH_INDEX: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   L => SCARC(NM)%LEVEL(NL)
   RECEIVE_OMESH_INDEX: DO NOM = 1, NMESHES

      SNODE = PROCESS(NOM)
      IF (PROCESS(NM)==SNODE) CYCLE RECEIVE_OMESH_INDEX

      OS => SCARC(NM)%OSCARC(NOM)
      IF (OS%NICMAX_S==0 .AND. OS%NICMAX_R==0) CYCLE RECEIVE_OMESH_INDEX

      SELECT_EXCHANGE_TYPE: SELECT CASE (TYPE_EXCHANGE)


         !> ---------------------------------------------------------------------------------------
         !> Exchange information about neighboring step size along internal boundary
         !> ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_BASIC)

            N_REQ = N_REQ+1
            CALL MPI_IRECV(OS%RECV_INT_BASIC(1),1,MPI_INTEGER,SNODE, &
                           TAG,MPI_COMM_WORLD,REQ(N_REQ),IERROR)

         !> ---------------------------------------------------------------------------------------
         !> Exchange information about neighboring wall data
         !> ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_WALL_INFO)

            N_REQ = N_REQ+1
            CALL MPI_IRECV(OS%RECV_INT(1),SIZE(OS%RECV_INT),MPI_INTEGER,SNODE, &
                           TAG,MPI_COMM_WORLD,REQ(N_REQ),IERROR)

         !> ---------------------------------------------------------------------------------------
         !> Exchange information about neighboring wall data
         !> ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_DISCRET)

            N_REQ = N_REQ+1
            OS%RECV_INT = 0
            CALL MPI_IRECV(OS%RECV_INT(1),SIZE(OS%RECV_INT),MPI_INTEGER,SNODE, &
                           TAG,MPI_COMM_WORLD,REQ(N_REQ),IERROR)

         !> ---------------------------------------------------------------------------------------
         !> Exchange information about neighboring grid dimensions
         !> ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_CELL_NUMBER)

            N_REQ = N_REQ+1
            CALL MPI_IRECV(OS%RECV_INT_BASIC(1),8,MPI_INTEGER,SNODE, &
                           TAG,MPI_COMM_WORLD,REQ(N_REQ),IERROR)

         !> ---------------------------------------------------------------------------------------
         !> Exchange information about neighboring step size along internal boundary
         !> ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_CELL_WIDTH)

            N_REQ = N_REQ+1
            CALL MPI_IRECV(OS%RECV_REAL_BASIC(1),1,MPI_DOUBLE_PRECISION,SNODE, &
                           TAG,MPI_COMM_WORLD,REQ(N_REQ),IERROR)

         !> ---------------------------------------------------------------------------------------
         !> Exchange neighboring grid type or CELL_INDEX
         !> ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_CELL_INDEX)

            N_REQ = N_REQ+1
            CALL MPI_IRECV(OS%RECV_INT(1),SIZE(OS%RECV_INT),MPI_INTEGER,&
                           SNODE,TAG,MPI_COMM_WORLD,REQ(N_REQ),IERROR)

         !> ---------------------------------------------------------------------------------------
         !> Perform exchanges for
         !>    - internal values for matrix-vector multiplication
         !>    - internal boundary values
         !>    - internal subdiagonal matrix values
         !>    - internal subdiagonal or ghost matrix values
         !> ---------------------------------------------------------------------------------------
         CASE DEFAULT

            N_REQ = N_REQ+1
            CALL MPI_IRECV(OS%RECV_REAL(1),SIZE(OS%RECV_REAL),MPI_DOUBLE_PRECISION,&
                           SNODE,TAG,MPI_COMM_WORLD,REQ(N_REQ),IERROR)

      END SELECT SELECT_EXCHANGE_TYPE
   ENDDO RECEIVE_OMESH_INDEX
ENDDO RECEIVE_MESH_INDEX

END SUBROUTINE SCARC_EXCHANGE_RECEIVE


!> ------------------------------------------------------------------------------------------------
!> Send data to neighbors (corresponds to MESH_EXCHANGE)
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_EXCHANGE_SEND (NL)
INTEGER, INTENT(IN) :: NL
INTEGER :: NM, NOM, NWL
INTEGER :: IW, IWL, IWG, ICOL, IPTR, ICPL
INTEGER :: IOR_NBR, IOR_OWN, IC, ICE, ICG, IFACE
INTEGER :: IX, IY, IZ
INTEGER :: IXW, IYW, IZW
INTEGER :: IXG, IYG, IZG
!INTEGER :: IXN, IYN, IZN
INTEGER :: I, J, K, LL
REAL(EB) :: ZSUM
INTEGER, POINTER, DIMENSION(:)     ::  RECV_INT, RECV_INT_BASIC
REAL(EB), POINTER, DIMENSION(:)     ::  RECV_REAL, RECV_REAL_BASIC
REAL(EB), POINTER, DIMENSION(:)     ::  VECTOR
REAL(EB), POINTER, DIMENSION(:,:,:) :: HVECTOR
TYPE (SCARC_TYPE), POINTER :: S=>NULL()
TYPE (SCARC_NEIGHBOR_TYPE), POINTER :: OS=>NULL()
TYPE (SCARC_POINTER_TYPE), POINTER :: PTR=>NULL(), OPTR=>NULL()
TYPE (SCARC_LEVEL_TYPE), POINTER :: L=>NULL(), OL=>NULL()
TYPE (SCARC_WALL_TYPE), POINTER :: WC=>NULL()
TYPE (SCARC_DISCRET_TYPE), POINTER :: D=>NULL(), OD=>NULL()
TYPE (SCARC_MATRIX_COMPACT_TYPE), POINTER :: AC=>NULL()

!> ------------------------------------------------------------------------------------------------
!> Collect data for sending corresponding to requested exchange type
!> ------------------------------------------------------------------------------------------------
MESH_PACK_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   L   => SCARC(NM)%LEVEL(NL)
   SELECT CASE(TYPE_DISCRET)
      CASE (NSCARC_DISCRET_STRUCTURED)
         D => L%SD
      CASE (NSCARC_DISCRET_UNSTRUCTURED)
         D => L%UD
   END SELECT
   PTR => D%PTR

   OMESH_PACK_LOOP: DO NOM = 1, NMESHES

      OS => SCARC(NM)%OSCARC(NOM)
      IF (OS%NICMAX_S == 0 .AND. OS%NICMAX_R == 0) CYCLE OMESH_PACK_LOOP

      OL   => SCARC(NM)%OSCARC(NOM)%LEVEL(NL)
      SELECT CASE(TYPE_DISCRET)
         CASE (NSCARC_DISCRET_STRUCTURED)
            OD => OL%SD
         CASE (NSCARC_DISCRET_UNSTRUCTURED)
            OD => OL%UD
      END SELECT
      OPTR => OD%PTR

      SNODE = PROCESS(NOM)
      RNODE = PROCESS(NM)

      OMESH_PACK_SELECT: SELECT CASE (TYPE_EXCHANGE)

         !> ---------------------------------------------------------------------------------------
         !> send needed size of exchange vector
         !> ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_BASIC)

            OS%SEND_INT_BASIC(1)=OL%NWL

            IF (RNODE /= SNODE) THEN
               N_REQ = N_REQ+1
               CALL MPI_ISEND(OS%SEND_INT_BASIC(1),1,MPI_INTEGER,SNODE, &
                              TAG,MPI_COMM_WORLD,REQ(N_REQ),IERROR)
            ENDIF

         !> ---------------------------------------------------------------------------------------
         !> send overlapped parts of cell widths
         !> ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_CELL_WIDTH)

            SELECT CASE(OL%IOR)
               CASE (1)
                  OL%DH = L%DXL(0)
               CASE (-1)
                  OL%DH = L%DXL(L%NX)
                CASE (2)
                  OL%DH = L%DYL(0)
               CASE (-2)
                  OL%DH = L%DYL(L%NY)
               CASE (3)
                  OL%DH = L%DZL(0)
               CASE (-3)
                  OL%DH = L%DZL(L%NZ)
            END SELECT
            OS%SEND_REAL_BASIC(1) = OL%DH

            IF (RNODE /= SNODE) THEN
               N_REQ = N_REQ+1
               CALL MPI_ISEND(OS%SEND_REAL_BASIC(1),1,MPI_DOUBLE_PRECISION,SNODE, &
                              TAG,MPI_COMM_WORLD,REQ(N_REQ),IERROR)
            ENDIF

         !> ---------------------------------------------------------------------------------------
         !> send overlapped parts of cell numbers
         !> ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_CELL_NUMBER)

            OS%SEND_INT_BASIC(1)=L%NX
            OS%SEND_INT_BASIC(2)=L%NY
            OS%SEND_INT_BASIC(3)=L%NZ
            OS%SEND_INT_BASIC(4)=L%NC
            OS%SEND_INT_BASIC(5)=L%NC
            OS%SEND_INT_BASIC(6)=L%N_WALL_CELLS
            OS%SEND_INT_BASIC(7)=L%N_WALL_CELLS_EXT
            OS%SEND_INT_BASIC(8)=L%N_WALL_CELLS_INT

            IF (RNODE /= SNODE) THEN
               N_REQ = N_REQ+1
               CALL MPI_ISEND(OS%SEND_INT_BASIC(1),8,MPI_INTEGER,SNODE, &
                              TAG,MPI_COMM_WORLD,REQ(N_REQ),IERROR)
            ENDIF


         !> ---------------------------------------------------------------------------------------
         !> send overlapped parts of wall information
         !> ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_WALL_INFO)

            IPTR=1
            DO IWL = 1, OL%NWL
               IWG = OPTR%IWL_TO_IWG(IWL)
               WC => D%WALL(IWG)
               OS%SEND_INT(IPTR   ) = WC%IXG
               OS%SEND_INT(IPTR+ 1) = WC%IYG
               OS%SEND_INT(IPTR+ 2) = WC%IZG
               OS%SEND_INT(IPTR+ 3) = WC%IXW
               OS%SEND_INT(IPTR+ 4) = WC%IYW
               OS%SEND_INT(IPTR+ 5) = WC%IZW
               OS%SEND_INT(IPTR+ 6) = WC%IXN(1)
               OS%SEND_INT(IPTR+ 7) = WC%IXN(2)
               OS%SEND_INT(IPTR+ 8) = WC%IYN(1)
               OS%SEND_INT(IPTR+ 9) = WC%IYN(2)
               OS%SEND_INT(IPTR+10) = WC%IZN(1)
               OS%SEND_INT(IPTR+11) = WC%IZN(2)
               OS%SEND_INT(IPTR+12) = WC%NOM
               IPTR = IPTR + 13
               DO ICPL=1,OL%NCPLS
                  OS%SEND_INT(IPTR)=WC%ICE(ICPL)
                  IPTR = IPTR + 1
               ENDDO
            ENDDO

            IF (RNODE /= SNODE) THEN
               N_REQ = N_REQ+1
               CALL MPI_ISEND(OS%SEND_INT(1),SIZE(OS%SEND_INT),MPI_INTEGER,SNODE, &
                              TAG,MPI_COMM_WORLD,REQ(N_REQ),IERROR)
            ENDIF

         !> ---------------------------------------------------------------------------------------
         !> send overlapped parts of grid information
         !> ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_DISCRET)

            IPTR=1
            OS%SEND_INT=0
            DO IWL = 1, OL%NWL
               IWG = OPTR%IWL_TO_IWG(IWL)
               WC => D%WALL(IWG)
               OS%SEND_INT(IPTR  ) = L%CELL_STATE(WC%IXW, WC%IYW, WC%IZW)
               OS%SEND_INT(IPTR+1) = D%CELL_NUMBER(WC%IXW, WC%IYW, WC%IZW)
               IPTR = IPTR + 2
            ENDDO

            IF (RNODE /= SNODE) THEN
               N_REQ = N_REQ+1
               CALL MPI_ISEND(OS%SEND_INT(1),SIZE(OS%SEND_INT),MPI_INTEGER,SNODE, &
                              TAG,MPI_COMM_WORLD,REQ(N_REQ),IERROR)
            ENDIF



         !> ---------------------------------------------------------------------------------------
         !> send overapped parts of pressure vector (H or HS) from predictor/corrector
         !> ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_PRESSURE)

            IF (PREDICTOR) THEN
               HVECTOR => POINT_TO_HVECTOR (NM, NSCARC_VECTOR_H)
            ELSE
               HVECTOR => POINT_TO_HVECTOR (NM, NSCARC_VECTOR_HS)
            ENDIF

            LL  = 1
            DO ICG=1, OL%NCG
               IWG = OPTR%ICG_TO_IWG(ICG)
               WC => D%WALL(IWG)
               IX  = WC%IXW
               IY  = WC%IYW
               IZ  = WC%IZW
               OS%SEND_REAL(LL) = HVECTOR(IX,IY,IZ)
               LL = LL+1
            ENDDO

            IF (RNODE/=SNODE) THEN
               N_REQ=N_REQ+1
               CALL MPI_ISEND(OS%SEND_REAL(1), SIZE(OS%SEND_REAL), MPI_DOUBLE_PRECISION, SNODE, &
                              TAG, MPI_COMM_WORLD, REQ(N_REQ), IERROR)
            ENDIF


         !> ---------------------------------------------------------------------------------------
         !> send overapped parts of specified vector
         !> ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_VECTOR)

            VECTOR => POINT_TO_VECTOR(NM, NL, TYPE_VECTOR)

            LL = 1
            DO ICG= 1, OL%NCG
               !ZSUM = 0.0_EB
               IWG = OPTR%ICG_TO_IWG(ICG)
               WC => D%WALL(IWG)
               IX  = WC%IXW
               IY  = WC%IYW
               IZ  = WC%IZW

               !ZSUM = VECTOR(D%CELL_NUMBER(IX, IY, IZ))
               !OS%SEND_REAL(LL) = ZSUM/REAL(OL%NCPLR,EB)

               IF (TYPE_DISCRET == NSCARC_DISCRET_UNSTRUCTURED .AND. L%CELL_STATE(IX, IY, IZ) /= NSCARC_CELL_GASPHASE) THEN
                  OS%SEND_REAL(LL) = NSCARC_HUGE_REAL
               ELSE
                  OS%SEND_REAL(LL) = VECTOR(D%CELL_NUMBER(IX, IY, IZ))
               ENDIF

               LL = LL + 1
            ENDDO

            IF (RNODE/=SNODE) THEN
               N_REQ=N_REQ+1
               CALL MPI_ISEND(OS%SEND_REAL(1), SIZE(OS%SEND_REAL), MPI_DOUBLE_PRECISION, SNODE, &
                              TAG, MPI_COMM_WORLD, REQ(N_REQ), IERROR)
            ENDIF

         !> ---------------------------------------------------------------------------------------
         !> send overapped cell indices
         !> ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_CELL_INDEX)

            OS%SEND_INT = 0
            LL = 1
            DO ICG= 1, OL%NCG
               IWG = OPTR%ICG_TO_IWG(ICG)
               WC => D%WALL(IWG)
               IX  = WC%IXW
               IY  = WC%IYW
               IZ  = WC%IZW

               IF (TYPE_DISCRET == NSCARC_DISCRET_UNSTRUCTURED .AND. L%CELL_STATE(IX, IY, IZ) /= NSCARC_CELL_GASPHASE) THEN
                  OS%SEND_INT(LL) = NSCARC_HUGE_INT
               ELSE
                  OS%SEND_INT(LL) = D%CELL_NUMBER(IX, IY, IZ)
               ENDIF

               LL = LL + 1
            ENDDO

            IF (RNODE/=SNODE) THEN
               N_REQ=N_REQ+1
               CALL MPI_ISEND(OS%SEND_INT(1), SIZE(OS%SEND_INT), MPI_INTEGER, SNODE, &
                              TAG, MPI_COMM_WORLD, REQ(N_REQ), IERROR)
            ENDIF

         !> ---------------------------------------------------------------------------------------
         !> send overlapped values of neighboring matrix 
         !> ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_MATRIX_VALUE)

            LL = 1
            DO ICG= 1, OL%NCG
               IWG = OPTR%ICG_TO_IWG(ICG)
               WC => D%WALL(IWG)
               IX  = WC%IXW
               IY  = WC%IYW
               IZ  = WC%IZW

               IF (TYPE_DISCRET == NSCARC_DISCRET_UNSTRUCTURED .AND. L%CELL_STATE(IX, IY, IZ) /= NSCARC_CELL_GASPHASE) THEN
                  OS%SEND_REAL(LL) = NSCARC_HUGE_REAL
               ELSE
                  IC = D%CELL_NUMBER(IX, IY, IZ)
                  ICE = OPTR%ICG_TO_ICE(ICG)
                  DO ICOL = AC%ROW(IC)+1, AC%ROW(IC+1)-1
                     IF (AC%COL(ICOL) == ICE) THEN
                       OS%SEND_REAL(LL) = AC%VAL(ICOL)
                       EXIT
                     ENDIF
                  ENDDO
               ENDIF

               LL = LL + 1
            ENDDO

            IF (RNODE/=SNODE) THEN
               N_REQ=N_REQ+1
               CALL MPI_ISEND(OS%SEND_REAL(1), SIZE(OS%SEND_REAL), MPI_DOUBLE_PRECISION, SNODE, &
                              TAG, MPI_COMM_WORLD, REQ(N_REQ), IERROR)
            ENDIF


         !> ---------------------------------------------------------------------------------------
         !> send size of neighboring matrix
         !> ---------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_MATRIX_SIZE)

            OS%SEND_INT_BASIC(1) = OD%AC%NA

            IF (RNODE /= SNODE) THEN
               N_REQ = N_REQ+1
               CALL MPI_ISEND(OS%SEND_INT_BASIC(1),1,MPI_INTEGER,SNODE,TAG,MPI_COMM_WORLD,REQ(N_REQ),IERROR)
            ENDIF

      END SELECT OMESH_PACK_SELECT
   ENDDO OMESH_PACK_LOOP
ENDDO MESH_PACK_LOOP


!> ------------------------------------------------------------------------------------------------
!> Information from Mesh NM is received by Mesh NOM  (NOM receiver, NM sender)
!> ------------------------------------------------------------------------------------------------
IF (N_MPI_PROCESSES>1.AND.N_REQ/=0) &
   CALL MPI_WAITALL(N_REQ,REQ(1:N_REQ),MPI_STATUSES_IGNORE,IERROR)

!> ------------------------------------------------------------------------------------------------
!> Extract communication data from corresponding RECEIVE-buffers
!> ------------------------------------------------------------------------------------------------
MESH_UNPACK_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   S   => SCARC(NM)
   L   => SCARC(NM)%LEVEL(NL)
   SELECT CASE(TYPE_DISCRET)
      CASE (NSCARC_DISCRET_STRUCTURED)
         D => L%SD
      CASE (NSCARC_DISCRET_UNSTRUCTURED)
         D => L%UD
   END SELECT
   PTR => D%PTR

   OMESH_UNPACK_LOOP: DO NOM=1,NMESHES

      SNODE  = PROCESS(NM)
      RNODE  = PROCESS(NOM)

      OS => SCARC(NM)%OSCARC(NOM)

      OMESH_UNPACK_IF: IF (OS%NICMAX_S/=0 .AND. OS%NICMAX_R/=0) THEN

         IF (RNODE/=SNODE) THEN
            RECV_INT        => SCARC(NM)%OSCARC(NOM)%RECV_INT
            RECV_INT_BASIC  => SCARC(NM)%OSCARC(NOM)%RECV_INT_BASIC
            RECV_REAL       => SCARC(NM)%OSCARC(NOM)%RECV_REAL
            RECV_REAL_BASIC => SCARC(NM)%OSCARC(NOM)%RECV_REAL_BASIC
         ELSE
            RECV_INT        => SCARC(NOM)%OSCARC(NM)%SEND_INT
            RECV_INT_BASIC  => SCARC(NOM)%OSCARC(NM)%SEND_INT_BASIC
            RECV_REAL       => SCARC(NOM)%OSCARC(NM)%SEND_REAL
            RECV_REAL_BASIC => SCARC(NOM)%OSCARC(NM)%SEND_REAL_BASIC
         ENDIF

         OL   => SCARC(NM)%OSCARC(NOM)%LEVEL(NL)
         SELECT CASE(TYPE_DISCRET)
            CASE (NSCARC_DISCRET_STRUCTURED)
               OD => OL%SD
            CASE (NSCARC_DISCRET_UNSTRUCTURED)
               OD => OL%UD
         END SELECT
         OPTR => OD%PTR

         OMESH_UNPACK_SELECT: SELECT CASE (TYPE_EXCHANGE)

            !> ------------------------------------------------------------------------------------
            !> unpack information about neighboring exchange size
            !> ------------------------------------------------------------------------------------
            CASE (NSCARC_EXCHANGE_BASIC)

               OL%NCG  = RECV_INT_BASIC(1)

            !> ------------------------------------------------------------------------------------
            !> unpack information about neighboring cell widths
            !> ------------------------------------------------------------------------------------
            CASE (NSCARC_EXCHANGE_CELL_WIDTH)

               OL%DH = RECV_REAL_BASIC(1)

               SELECT CASE (OL%IOR)
                  CASE ( 1)
                     L%DXL(0)    = 0.5_EB*(OL%DH + L%DXL(0))
                  CASE (-1)
                     L%DXL(L%NX) = 0.5_EB*(OL%DH + L%DXL(L%NX))
                  CASE ( 2)
                     L%DYL(0)    = 0.5_EB*(OL%DH + L%DYL(0))
                  CASE (-2)
                     L%DYL(L%NY) = 0.5_EB*(OL%DH + L%DYL(L%NY))
                  CASE ( 3)
                     L%DZL(0)    = 0.5_EB*(OL%DH + L%DZL(0))
                  CASE (-3)
                     L%DZL(L%NZ) = 0.5_EB*(OL%DH + L%DZL(L%NZ))
               END SELECT

            !> ------------------------------------------------------------------------------------
            !> unpack information about neighboring walls
            !> ------------------------------------------------------------------------------------
            CASE (NSCARC_EXCHANGE_WALL_INFO)
               IPTR=1
               DO ICG = 1, OL%NCG
                  WC => OD%WALL(ICG)
                  WC%IXG    = RECV_INT(IPTR    )
                  WC%IYG    = RECV_INT(IPTR + 1)
                  WC%IZG    = RECV_INT(IPTR + 2)
                  WC%IXW    = RECV_INT(IPTR + 3)
                  WC%IYW    = RECV_INT(IPTR + 4)
                  WC%IZW    = RECV_INT(IPTR + 5)
                  WC%IXN(1) = RECV_INT(IPTR + 6)
                  WC%IXN(2) = RECV_INT(IPTR + 7)
                  WC%IYN(1) = RECV_INT(IPTR + 8)
                  WC%IYN(2) = RECV_INT(IPTR + 9)
                  WC%IZN(1) = RECV_INT(IPTR +10)
                  WC%IZN(2) = RECV_INT(IPTR +11)
                  WC%NOM    = RECV_INT(IPTR +12)
                  IPTR = IPTR + 13
                  ALLOCATE (WC%ICE(OL%NCPLR))
                  DO ICPL=1,OL%NCPLR
                     WC%ICE(ICPL) = RECV_INT(IPTR)
                     IPTR = IPTR + 1
                  ENDDO
               ENDDO

            !> ------------------------------------------------------------------------------------
            !> unpack information about neighboring grids
            !> ------------------------------------------------------------------------------------
            CASE (NSCARC_EXCHANGE_DISCRET)

               DO ICG = 1, OL%NCG

                  OD%WALL(ICG)%STATE = RECV_INT(IPTR  )
                  IPTR = IPTR + 1

                  IWG = OPTR%ICG_TO_IWG(ICG)
                  WC => D%WALL(IWG)

                  IXG = WC%IXG
                  IYG = WC%IYG
                  IZG = WC%IZG

                  IXW = WC%IXW
                  IYW = WC%IYW
                  IZW = WC%IZW

                  !IF (OD%WALL(ICG)%STATE == NSCARC_CELL_GASPHASE.AND.L%CELL_STATE(IXW,IYW,IZW) /= NSCARC_CELL_SOLID) THEN
                  IF (OD%WALL(ICG)%STATE == NSCARC_CELL_GASPHASE) L%CELL_STATE(IXG, IYG, IZG) = OD%WALL(ICG)%STATE

               ENDDO


            !> ------------------------------------------------------------------------------------
            !> unpack information about neighboring cell numbers
            !> ------------------------------------------------------------------------------------
            CASE (NSCARC_EXCHANGE_CELL_NUMBER)

               OL%NX = RECV_INT_BASIC(1)
               OL%NY = RECV_INT_BASIC(2)
               OL%NZ = RECV_INT_BASIC(3)
               OL%NC = RECV_INT_BASIC(4)
               OL%NC = RECV_INT_BASIC(5)
               OL%N_WALL_CELLS     = RECV_INT_BASIC(6)
               OL%N_WALL_CELLS_EXT = RECV_INT_BASIC(7)
               OL%N_WALL_CELLS_INT = RECV_INT_BASIC(8)

            !> ------------------------------------------------------------------------------------
            !> unpack information about neighboring matrix sizes
            !> ------------------------------------------------------------------------------------
            CASE (NSCARC_EXCHANGE_MATRIX_SIZE)

               SELECT CASE (TYPE_MATRIX)
                  CASE (NSCARC_MATRIX_COMPACT)
                     OD%AC%NA = RECV_INT_BASIC(1)
                  CASE (NSCARC_MATRIX_BANDED)
                     OD%AB%NA = RECV_INT_BASIC(1)
               END SELECT

            !> ------------------------------------------------------------------------------------
            !> unpack overlapping parts of pressure vectors H or HS from predictor/corrector
            !> ------------------------------------------------------------------------------------
            CASE (NSCARC_EXCHANGE_PRESSURE)

               IF (PREDICTOR) THEN
                  HVECTOR => POINT_TO_HVECTOR (NM, NSCARC_VECTOR_H)
               ELSE
                  HVECTOR => POINT_TO_HVECTOR (NM, NSCARC_VECTOR_HS)
               ENDIF

               LL = 1
               DO IFACE = 1, 6
                  IOR_NBR = FACE_ORIENTATION(IFACE)
                  IWL = OL%SUBDIVISION(1, -IOR_NBR)
                  NWL = OL%SUBDIVISION(2, -IOR_NBR)

                  UNPACK_PRESSURE: DO IW = IWL, IWL + NWL - 1

                     ZSUM=0.0_EB
                     IWG = OPTR%IWL_TO_IWG(IW)
                     DO ICPL = 1, OL%NCPL
                        ZSUM=ZSUM+RECV_REAL(LL)
                        LL = LL+1
                     ENDDO

                     WC => D%WALL(IWG)
                     I=WC%IXG
                     J=WC%IYG
                     K=WC%IZG

                     HVECTOR(I, J, K) = ZSUM/REAL(OL%NCPL,EB)

                  ENDDO UNPACK_PRESSURE
               ENDDO

            !> ------------------------------------------------------------------------------------
            !> unpack overlapping parts of a specified vector
            !> ------------------------------------------------------------------------------------
            CASE (NSCARC_EXCHANGE_VECTOR)

               VECTOR => POINT_TO_VECTOR (NM, NL, TYPE_VECTOR)

               LL = 1
               DO IFACE = 1, 6
                  IOR_NBR = FACE_ORIENTATION(IFACE)
                  IWL = OL%SUBDIVISION(1, -IOR_NBR)
                  NWL = OL%SUBDIVISION(2, -IOR_NBR)

                  DO IW = IWL, IWL + NWL - 1
                     IWG = OPTR%IWL_TO_IWG(IW)

                     ICE     = D%WALL(IWG)%ICE(1)
                     IOR_OWN = D%WALL(IWG)%IOR

                     VECTOR(ICE) = RECV_REAL(LL)

                     LL = LL + 1
                  ENDDO
               ENDDO

            !> ------------------------------------------------------------------------------------
            !> unpack information about neighboring cell indices
            !> ------------------------------------------------------------------------------------
            CASE (NSCARC_EXCHANGE_CELL_INDEX)

               LL = 1
               DO IFACE = 1, 6
                  IOR_NBR = FACE_ORIENTATION(IFACE)
                  IWL = OL%SUBDIVISION(1, -IOR_NBR)
                  NWL = OL%SUBDIVISION(2, -IOR_NBR)

                  DO IW = IWL, IWL + NWL - 1
                     IWG = OPTR%IWL_TO_IWG(IW)
                     ICE     = D%WALL(IWG)%ICE(1)
                     IOR_OWN = D%WALL(IWG)%IOR

                     PTR%ICE_TO_ICN(ICE) = RECV_INT(LL)
                     LL = LL + 1
                  ENDDO
               ENDDO

            !> ------------------------------------------------------------------------------------
            !> unpack information about neighboring matrix values
            !> ------------------------------------------------------------------------------------
            CASE (NSCARC_EXCHANGE_MATRIX_VALUE)

               LL = 1
               DO IFACE = 1, 6
                  IOR_NBR = FACE_ORIENTATION(IFACE)
                  IWL = OL%SUBDIVISION(1, -IOR_NBR)
                  NWL = OL%SUBDIVISION(2, -IOR_NBR)

                  DO IW = IWL, IWL + NWL - 1
                     IWG = OPTR%IWL_TO_IWG(IW)
                     ICE     = D%WALL(IWG)%ICE(1)
                     IOR_OWN = D%WALL(IWG)%IOR

                     PTR%ICE_TO_VAL(ICE) = RECV_REAL(LL)

                     LL = LL + 1
                  ENDDO
               ENDDO

          END SELECT OMESH_UNPACK_SELECT
      ENDIF OMESH_UNPACK_IF
   ENDDO OMESH_UNPACK_LOOP
ENDDO MESH_UNPACK_LOOP

END SUBROUTINE SCARC_EXCHANGE_SEND


!> ------------------------------------------------------------------------------------------------
!> Check if difference of two values is less than a given tolerance
!> ------------------------------------------------------------------------------------------------
LOGICAL FUNCTION MATCH (VAL1, VAL2)
REAL (EB), INTENT(IN) :: VAL1, VAL2
REAL (EB) :: TOL
TOL = 1.0E-10_EB
MATCH = .FALSE.
IF (Abs(VAL1-VAL2) <= TOL) MATCH = .TRUE.
RETURN
END FUNCTION MATCH

!> ------------------------------------------------------------------------------------------------
!> Print out timings for ScaRC - not updated at the moment
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_TIMINGS

#ifdef WITH_SCARC_VERBOSE
IF (MYID == 0) THEN
   IF (N_MPI_PROCESSES == 1) THEN
      WRITE(LU_OUTPUT,1002)
      WRITE(LU_OUTPUT,1001)
      WRITE(LU_OUTPUT,1002)
   ELSE
      WRITE(LU_OUTPUT,2002)
      WRITE(LU_OUTPUT,2001)
      WRITE(LU_OUTPUT,2002)
   ENDIF
ENDIF

IF (MYID == 0) THEN
   IF (N_MPI_PROCESSES == 1) THEN
      WRITE(LU_OUTPUT,1002)
   ELSE
      WRITE(LU_OUTPUT,2002)
   ENDIF
ENDIF
#endif

1001 FORMAT('| Scarc routine       |    Time (s)    |')
1002 FORMAT('|---------------------|----------------|')
2001 FORMAT('| Scarc routine       |  Time_min (s)  |  Time_max (s)  | Time_mean (s)  |')
2002 FORMAT('|---------------------|----------------|----------------|----------------|')
END SUBROUTINE SCARC_TIMINGS



!> ----------------------------------------------------------------------------------------------------
!> Filter out mean value
!> --------------------------------------------------------------------------------------------------ss
SUBROUTINE SCARC_FILTER_MEANVALUE(NV, NL)
INTEGER, INTENT(IN) :: NV, NL
INTEGER :: NM, IC, I, J, K
REAL(EB), DIMENSION(:) , POINTER :: VC
TYPE (SCARC_LEVEL_TYPE), POINTER :: L=>NULL()
TYPE (SCARC_DISCRET_TYPE), POINTER :: D=>NULL()

LOCAL_REAL = 0.0_EB
DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   L => SCARC(NM)%LEVEL(NL)
   SELECT CASE(TYPE_DISCRET)
      CASE (NSCARC_DISCRET_STRUCTURED)
         D => L%SD
      CASE (NSCARC_DISCRET_UNSTRUCTURED)
         D => L%UD
   END SELECT

   VC => POINT_TO_VECTOR(NM, NL, NV)
   DO K = 1, L%NZ
      DO J = 1, L%NY
         DO I = 1, L%NX
            IF (TYPE_DISCRET == NSCARC_DISCRET_UNSTRUCTURED .AND. L%CELL_STATE(I, J, K) /= NSCARC_CELL_GASPHASE) CYCLE
            IC = D%CELL_NUMBER(I,J,K)
            LOCAL_REAL(NM) = LOCAL_REAL(NM) + VC(IC)
         ENDDO
      ENDDO
   ENDDO
ENDDO

IF (N_MPI_PROCESSES > 1) &
   CALL MPI_ALLGATHERV(MPI_IN_PLACE,1,MPI_INTEGER,LOCAL_REAL,COUNTS,DISPLS,&
                       MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,IERROR)

GLOBAL_REAL = SUM(LOCAL_REAL(1:NMESHES))/REAL(N_CELLS_GLOBAL(NL))

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   L   => SCARC(NM)%LEVEL(NL)
   SELECT CASE(TYPE_DISCRET)
      CASE (NSCARC_DISCRET_STRUCTURED)
         D => L%SD
      CASE (NSCARC_DISCRET_UNSTRUCTURED)
         D => L%UD
   END SELECT
   VC  => POINT_TO_VECTOR(NM, NL, NV)

   DO K = 1, L%NZ
      DO J = 1, L%NY
         DO I = 1, L%NX
            IF (TYPE_DISCRET == NSCARC_DISCRET_UNSTRUCTURED .AND. L%CELL_STATE(I, J, K) /= NSCARC_CELL_GASPHASE) CYCLE
            IC = D%CELL_NUMBER(I,J,K)
            VC(IC) = VC(IC) - GLOBAL_REAL
         ENDDO
      ENDDO
   ENDDO

ENDDO

END SUBROUTINE SCARC_FILTER_MEANVALUE

!> ------------------------------------------------------------------------------------------------
!> Velocity predictor for MGM-method
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_VELOCITY_PREDICTOR(NL)
INTEGER, INTENT(IN) :: NL
INTEGER :: NM, IW, I, J, K, IC
TYPE (MESH_TYPE), POINTER :: M=>NULL()
TYPE (SCARC_LEVEL_TYPE), POINTER :: L=>NULL()
TYPE (SCARC_WALL_TYPE), POINTER :: WC=>NULL()
TYPE (SCARC_MGM_TYPE), POINTER :: MGM=>NULL()

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   M   => MESHES(NM)
   L   => SCARC(NM)%LEVEL(NL)
   MGM => L%MGM

   DO IW = L%N_WALL_CELLS_EXT+1, L%N_WALL_CELLS_EXT + L%N_WALL_CELLS_INT
   
      WC => L%SD%WALL(IW)                        !> based on structured discretization

      I = WC%IXW 
      J = WC%IYW 
      K = WC%IZW 

      IC = WC%ICW

      MGM%US(IC) = M%U(I,J,K) - DT*( M%FVX(I,J,K) + M%RDXN(I)*(M%H(I+1,J  ,K  ) - M%H(I,J,K)) )
      MGM%VS(IC) = M%V(I,J,K) - DT*( M%FVY(I,J,K) + M%RDYN(J)*(M%H(I  ,J+1,K  ) - M%H(I,J,K)) )
      MGM%WS(IC) = M%W(I,J,K) - DT*( M%FVZ(I,J,K) + M%RDZN(K)*(M%H(I  ,J  ,K+1) - M%H(I,J,K)) )
 
   ENDDO

ENDDO

END SUBROUTINE SCARC_VELOCITY_PREDICTOR


!> ------------------------------------------------------------------------------------------------
!> Restore last cell of last mesh
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_RESTORE_LAST_CELL (XX, NL)
INTEGER, INTENT(IN) :: XX, NL
REAL (EB), POINTER, DIMENSION(:) :: VX
TYPE (SCARC_TYPE), POINTER :: S=>NULL()

IF (UPPER_MESH_INDEX /= NMESHES .OR. TYPE_RELAX == NSCARC_RELAX_FFT) RETURN
S => SCARC(UPPER_MESH_INDEX)

VX => POINT_TO_VECTOR (UPPER_MESH_INDEX, NL, XX)
VX(S%N_CELLS) = S%RHS_END

END SUBROUTINE SCARC_RESTORE_LAST_CELL


!> ================================================================================================
!>  Set of different pointer functions which return pointers to specified structures
!> ================================================================================================

!> ------------------------------------------------------------------------------------------------
!> Set pointer to chosen vector for compact storage technique
!> ------------------------------------------------------------------------------------------------
FUNCTION POINT_TO_VECTOR(NM, NL, NV)
REAL(EB), POINTER, DIMENSION(:) :: POINT_TO_VECTOR
INTEGER, INTENT(IN) :: NM, NL, NV
TYPE (SCARC_LEVEL_TYPE), POINTER :: L=>NULL()

L => SCARC(NM)%LEVEL(NL)
SELECT CASE (NV)
   !>
   !> Stage one vectors
   !>
   CASE (NSCARC_VECTOR_ONE_X)
      POINT_TO_VECTOR => L%STAGE(NSCARC_STAGE_ONE)%X
   CASE (NSCARC_VECTOR_ONE_B)
      POINT_TO_VECTOR => L%STAGE(NSCARC_STAGE_ONE)%B
   CASE (NSCARC_VECTOR_ONE_D)
      POINT_TO_VECTOR => L%STAGE(NSCARC_STAGE_ONE)%D
   CASE (NSCARC_VECTOR_ONE_R)
      POINT_TO_VECTOR => L%STAGE(NSCARC_STAGE_ONE)%R
   CASE (NSCARC_VECTOR_ONE_V)
      POINT_TO_VECTOR => L%STAGE(NSCARC_STAGE_ONE)%V
   CASE (NSCARC_VECTOR_ONE_Y)
      POINT_TO_VECTOR => L%STAGE(NSCARC_STAGE_ONE)%Y
   CASE (NSCARC_VECTOR_ONE_Z)
      POINT_TO_VECTOR => L%STAGE(NSCARC_STAGE_ONE)%Z
#ifdef WITH_SCARC_DEBUG
   CASE (NSCARC_VECTOR_ONE_E)
      POINT_TO_VECTOR => L%STAGE(NSCARC_STAGE_ONE)%E
#endif
   !>
   !> Stage two vectors
   !>
   CASE (NSCARC_VECTOR_TWO_X)
      POINT_TO_VECTOR => L%STAGE(NSCARC_STAGE_TWO)%X
   CASE (NSCARC_VECTOR_TWO_B)
      POINT_TO_VECTOR => L%STAGE(NSCARC_STAGE_TWO)%B
   CASE (NSCARC_VECTOR_TWO_D)
      POINT_TO_VECTOR => L%STAGE(NSCARC_STAGE_TWO)%D
   CASE (NSCARC_VECTOR_TWO_R)
      POINT_TO_VECTOR => L%STAGE(NSCARC_STAGE_TWO)%R
   CASE (NSCARC_VECTOR_TWO_V)
      POINT_TO_VECTOR => L%STAGE(NSCARC_STAGE_TWO)%V
   CASE (NSCARC_VECTOR_TWO_Y)
      POINT_TO_VECTOR => L%STAGE(NSCARC_STAGE_TWO)%Y
   CASE (NSCARC_VECTOR_TWO_Z)
      POINT_TO_VECTOR => L%STAGE(NSCARC_STAGE_TWO)%Z
#ifdef WITH_SCARC_DEBUG
   CASE (NSCARC_VECTOR_TWO_E)
      POINT_TO_VECTOR => L%STAGE(NSCARC_STAGE_TWO)%E
#endif
END SELECT

END FUNCTION POINT_TO_VECTOR

#ifdef WITH_MKL_FB
!> ------------------------------------------------------------------------------------------------
!> Set pointer to chosen vector for compact storage technique
!> ------------------------------------------------------------------------------------------------
FUNCTION POINT_TO_VECTOR_FB(NM, NL, NV)
REAL(FB), POINTER, DIMENSION(:) :: POINT_TO_VECTOR_FB
INTEGER, INTENT(IN) :: NM, NL, NV
TYPE (SCARC_LEVEL_TYPE), POINTER :: L=>NULL()

WRITE(*,*) 'CHECK NAMES OF FB-VECTORS !!!!'

L => SCARC(NM)%LEVEL(NL)
SELECT CASE (NV)
   CASE (NSCARC_VECTOR_ONE_X)
      POINT_TO_VECTOR_FB => L%STAGE(NSCARC_STAGE_ONE)%X_FB
   CASE (NSCARC_VECTOR_ONE_B)
      POINT_TO_VECTOR_FB => L%STAGE(NSCARC_STAGE_ONE)%B_FB
   CASE (NSCARC_VECTOR_ONE_D)
      POINT_TO_VECTOR_FB => L%STAGE(NSCARC_STAGE_ONE)%D_FB
   CASE (NSCARC_VECTOR_ONE_R)
      POINT_TO_VECTOR_FB => L%STAGE(NSCARC_STAGE_ONE)%R_FB
   CASE (NSCARC_VECTOR_ONE_V)
      POINT_TO_VECTOR_FB => L%STAGE(NSCARC_STAGE_ONE)%V_FB
END SELECT

END FUNCTION POINT_TO_VECTOR_FB
#endif

!> ------------------------------------------------------------------------------------------------
!> Set pointer to chosen vector for banded storage technique
!> ------------------------------------------------------------------------------------------------
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


!> ================================================================================================
!>> Routines for the allocation of different work spaces
!> ================================================================================================
!> ------------------------------------------------------------------------------------------------
!> Allocate and initialize integer array of dimension 1
!> ------------------------------------------------------------------------------------------------
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
         WORKSPACE = NSCARC_UNDEFINED_INT
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
1000 FORMAT('Allocating INT1  array ',A,' in length (',I8,':',I8,')')
END SUBROUTINE SCARC_ALLOCATE_INT1

!> ------------------------------------------------------------------------------------------------
!> Allocate and initialize integer array of dimension 2
!> ------------------------------------------------------------------------------------------------
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
         WORKSPACE = NSCARC_UNDEFINED_INT
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
1000 FORMAT('Allocating INT2  array ',A,' in length (',I8,':',I8,' , ',I8,':',I8,')')
END SUBROUTINE SCARC_ALLOCATE_INT2

!> ------------------------------------------------------------------------------------------------
!> Allocate and initialize integer array of dimension 3
!> ------------------------------------------------------------------------------------------------
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
         WORKSPACE = NSCARC_UNDEFINED_INT
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
1000 FORMAT('Allocating INT3  array ',A,' in length (',I8,':',I8,' , ',I8,':',I8,' , ',I8,':',I8,')')
END SUBROUTINE SCARC_ALLOCATE_INT3

!> ------------------------------------------------------------------------------------------------
!> Allocate and initialize real array of dimension 1
!> ------------------------------------------------------------------------------------------------
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
         WORKSPACE = NSCARC_UNDEFINED_REAL_EB
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
1000 FORMAT('Allocating REAL1 array ',A,' in length (',I8,':',I8,')')
END SUBROUTINE SCARC_ALLOCATE_REAL1


#ifdef WITH_MKL_FB
!> ------------------------------------------------------------------------------------------------
!> Allocate and initialize real array of dimension 1
!> ------------------------------------------------------------------------------------------------
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
         WORKSPACE = NSCARC_UNDEFINED_REAL_FB
      CASE (NSCARC_INIT_ZERO)
         WORKSPACE = NSCARC_ZERO_REAL_FB
      CASE (NSCARC_INIT_ONE)
         WORKSPACE = NSCARC_ONE_REAL_FB
   END SELECT
ELSE
   IF (SIZE(WORKSPACE,1) /= NR1-NL1+1) &
      CALL SCARC_SHUTDOWN(NSCARC_ERROR_VECTOR_LENGTH, CTEXT, NSCARC_NONE)
ENDIF
END SUBROUTINE SCARC_ALLOCATE_REAL1_FB
#endif

!> ------------------------------------------------------------------------------------------------
!> Allocate and initialize real array of dimension 2
!> ------------------------------------------------------------------------------------------------
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
         WORKSPACE = NSCARC_UNDEFINED_REAL_EB
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
1000 FORMAT('Allocating REAL2 array ',A,' in length (',I8,':',I8,' , ',I8,':',I8,')')
END SUBROUTINE SCARC_ALLOCATE_REAL2

!> ------------------------------------------------------------------------------------------------
!> Allocate and initialize real array of dimension 3
!> ------------------------------------------------------------------------------------------------
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
         WORKSPACE = NSCARC_UNDEFINED_REAL_EB
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
1000 FORMAT('Allocating REAL3 array ',A,' in length (',I8,':',I8,' , ',I8,':',I8,' , ',I8,':',I8,')')
END SUBROUTINE SCARC_ALLOCATE_REAL3

!> ------------------------------------------------------------------------------------------------
!> Allocate matrix with corresponding pointer and length structures
!> COMPACT storage format
!>    - allocate arrays for matrix values, column and row pointers
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_ALLOCATE_MATRIX_COMPACT(AC, CMATRIX, NL, NINIT)
TYPE (SCARC_MATRIX_COMPACT_TYPE), INTENT(INOUT) :: AC
INTEGER, INTENT(IN) :: NINIT, NL
CHARACTER(*), INTENT(IN) :: CMATRIX
CHARACTER(40) :: CINFO

#ifdef WITH_MKL_FB
WRITE(CINFO,'(A,A,I2.2,A)') TRIM(CMATRIX),'_LEV',NL,'.VAL_FB'
CALL SCARC_ALLOCATE_REAL1_FB(AC%VAL_FB, 1, AC%NA, NINIT, CINFO)
#else
WRITE(CINFO,'(A,A,I2.2,A)') TRIM(CMATRIX),'_LEV',NL,'.VAL'
CALL SCARC_ALLOCATE_REAL1(AC%VAL, 1, AC%NA, NINIT, CINFO)
#endif

WRITE(CINFO,'(A,A,I2.2,A)') TRIM(CMATRIX),'_LEV',NL,'.ROW'
CALL SCARC_ALLOCATE_INT1(AC%ROW, 1, AC%NR, NINIT, CINFO)

WRITE(CINFO,'(A,A,I2.2,A)') TRIM(CMATRIX),'_LEV',NL,'.COL'
CALL SCARC_ALLOCATE_INT1(AC%COL, 1, AC%NA, NINIT, CINFO)

#ifdef WITH_MKL
IF (IS_MKL_LEVEL(NL)) THEN
   WRITE(CINFO,'(A,A)') TRIM(CMATRIX),'.COL_GLOBAL'
   CALL SCARC_ALLOCATE_INT1 (AC%COL_GLOBAL, 1, AC%NA, NINIT, CINFO)
ENDIF
#endif

#ifdef WITH_SCARC_VERBOSE
WRITE(MSG%LU_VERBOSE,1000) CMATRIX, AC%NA, AC%NR
#endif
1000 FORMAT('Allocating COMPACT matrix ',A,' with NA=',I8,' and NR=',I8)
END SUBROUTINE SCARC_ALLOCATE_MATRIX_COMPACT

!> ------------------------------------------------------------------------------------------------
!> Get a full copy of matrix structure AC1 in AC2
!> If requested, the original matrix structure AC1 is completely removed
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_COPY_MATRIX_COMPACT(AC1, AC2, CMATRIX1, CMATRIX2, BREMOVE, NL)
TYPE (SCARC_MATRIX_COMPACT_TYPE), INTENT(INOUT) :: AC1, AC2
INTEGER, INTENT(IN) :: NL
LOGICAL :: BREMOVE
CHARACTER(*), INTENT(IN) :: CMATRIX1, CMATRIX2
CHARACTER(40) :: CINFO

#ifdef WITH_SCARC_VERBOSE
WRITE(MSG%LU_VERBOSE,*) 'Trying to copy matrix ', CMATRIX1, ' into ',CMATRIX2
#endif

AC2%NA       = AC1%NA
AC2%NR       = AC1%NR 
AC2%NSTORE   = AC1%NSTORE 
AC2%POS      = AC1%POS

AC2%NSTENCIL = AC1%NSTENCIL
AC2%STENCIL  = AC1%STENCIL

IF (.NOT.ALLOCATED(AC2%VAL)) THEN
   WRITE(CINFO,'(A,A,I2.2,A)') TRIM(CMATRIX2),'_LEV',NL,'.VAL'
   CALL SCARC_ALLOCATE_REAL1(AC2%VAL, 1, AC2%NA, NSCARC_INIT_NONE, CINFO)
ELSE IF (SIZE(AC2%VAL) /= AC1%NA) THEN
   CALL SCARC_SHUTDOWN(NSCARC_ERROR_MATRIX_COPY, CMATRIX2, NSCARC_NONE)
ENDIF
AC2%VAL = AC1%VAL
IF (BREMOVE) DEALLOCATE(AC1%VAL)

IF (.NOT.ALLOCATED(AC2%ROW)) THEN
   WRITE(CINFO,'(A,A,I2.2,A)') TRIM(CMATRIX2),'_LEV',NL,'.ROW'
   CALL SCARC_ALLOCATE_INT1(AC2%ROW, 1, AC2%NR, NSCARC_INIT_NONE, CINFO)
ELSE IF (SIZE(AC2%ROW) /= AC1%NR+1) THEN
   CALL SCARC_SHUTDOWN(NSCARC_ERROR_MATRIX_COPY, CMATRIX2, NSCARC_NONE)
ENDIF
AC2%ROW = AC1%ROW
IF (BREMOVE) DEALLOCATE(AC1%ROW)

IF (.NOT.ALLOCATED(AC2%VAL)) THEN
   WRITE(CINFO,'(A,A,I2.2,A)') TRIM(CMATRIX2),'_LEV',NL,'.COL'
   CALL SCARC_ALLOCATE_INT1(AC2%COL, 1, AC2%NA, NSCARC_INIT_NONE, CINFO)
ELSE IF (SIZE(AC2%COL) /= AC1%NA) THEN
   CALL SCARC_SHUTDOWN(NSCARC_ERROR_MATRIX_COPY, CMATRIX2, NSCARC_NONE)
ENDIF
AC2%COL = AC1%COL
IF (BREMOVE) DEALLOCATE(AC1%COL)

#ifdef WITH_MKL
IF (.NOT.ALLOCATED(AC2%COL_GLOBAL) .AND. IS_MKL_LEVEL(NL)) THEN
   WRITE(CINFO,'(A,A)') TRIM(CMATRIX2),'.COL_GLOBAL'
   CALL SCARC_ALLOCATE_INT1 (AC2%COL_GLOBAL, 1, AC2%NA, NSCARC_INIT_NONE, CINFO)
ELSE IF (SIZE(AC2%COL) /= AC1%NA) THEN
   CALL SCARC_SHUTDOWN(NSCARC_ERROR_MATRIX_COPY, CMATRIX2, NSCARC_NONE)
ENDIF
AC2%COL_GLOBAL = AC1%COL_GLOBAL
IF (BREMOVE) DEALLOCATE(AC1%COL_GLOBAL)
#endif

END SUBROUTINE SCARC_COPY_MATRIX_COMPACT


!> ------------------------------------------------------------------------------------------------
!> Allocate matrix with corresponding pointer and length structures
!> Diagoal storage format
!>    - allocate arrays for matrix values and offset vectors (in ascending order)
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_ALLOCATE_MATRIX_BANDED(AB, CMATRIX, NL, NINIT)
TYPE (SCARC_MATRIX_BANDED_TYPE), INTENT(INOUT) :: AB
INTEGER, INTENT(IN) :: NINIT, NL
CHARACTER(*), INTENT(IN) :: CMATRIX
CHARACTER(40) :: CINFO

#ifdef WITH_MKL_FB
WRITE(CINFO,'(A,A,I2.2,A)') TRIM(CMATRIX),'_LEV',NL,'.VAL_FB'
CALL SCARC_ALLOCATE_REAL2_FB(AB%VAL_FB, 1, AB%NSTENCIL , 1, AB%NDIAG, NINIT, CINFO)
#else
WRITE(CINFO,'(A,A,I2.2,A)') TRIM(CMATRIX),'_LEV',NL,'.VAL'
CALL SCARC_ALLOCATE_REAL2(AB%VAL, 1, AB%NSTENCIL , 1, AB%NDIAG, NINIT, CINFO)
#endif

END SUBROUTINE SCARC_ALLOCATE_MATRIX_BANDED

!> ------------------------------------------------------------------------------------------------
!> Allocate matrix with corresponding pointer and length structures
!> COMPACT storage format
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_DEALLOCATE_MATRIX_COMPACT(AC)
TYPE (SCARC_MATRIX_COMPACT_TYPE), INTENT(INOUT) :: AC

AC%NSTENCIL = 0
AC%STENCIL = 0

#ifdef WITH_MKL_FB
DEALLOCATE(AC%VAL_FB)
#else
DEALLOCATE(AC%VAL)
#endif
DEALLOCATE(AC%COL)
DEALLOCATE(AC%ROW)

AC%NA  = 0
AC%NR  = 0

END SUBROUTINE SCARC_DEALLOCATE_MATRIX_COMPACT


!> ------------------------------------------------------------------------------------------------
!> Reduce size of matrix to specified size
!> COMPACT storage format
!>    - reduce values and column pointer arrays
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_RESIZE_MATRIX_COMPACT(AC, NSIZE, CMATRIX)
TYPE (SCARC_MATRIX_COMPACT_TYPE), INTENT(INOUT) :: AC
INTEGER, INTENT(IN) :: NSIZE
CHARACTER(*), INTENT(IN) :: CMATRIX
CHARACTER(40) :: CINFO
#ifdef WITH_MKL_FB
REAL(FB), ALLOCATABLE, DIMENSION(:) :: VAL_FB
#else
REAL(EB), ALLOCATABLE, DIMENSION(:) :: VAL
#endif
INTEGER, ALLOCATABLE, DIMENSION(:) :: COL

IF (AC%NA == NSIZE) THEN
   RETURN                                  !> matrix has already desired size
ELSE IF (AC%NA > NSIZE) THEN

   AC%NA = NSIZE

   WRITE(CINFO,'(A,A)') TRIM(CMATRIX),'.VAL'

#ifdef WITH_MKL_FB

   CALL SCARC_ALLOCATE_REAL1_FB(VAL_FB, 1, NSIZE, NSCARC_INIT_NONE, CINFO)
   VAL_FB(1:NSIZE) = AC%VAL_FB(1:NSIZE)
   DEALLOCATE(AC%VAL_FB)

   CALL SCARC_ALLOCATE_REAL1_FB(AC%VAL_FB, 1, NSIZE, NSCARC_INIT_NONE, CINFO)
   AC%VAL_FB(1:NSIZE) = VAL_FB(1:NSIZE)
   DEALLOCATE(VAL_FB)

#else

   CALL SCARC_ALLOCATE_REAL1(VAL, 1, NSIZE, NSCARC_INIT_NONE, CINFO)
   VAL(1:NSIZE) = AC%VAL(1:NSIZE)
   DEALLOCATE(AC%VAL)

   CALL SCARC_ALLOCATE_REAL1(AC%VAL, 1, NSIZE, NSCARC_INIT_NONE, CINFO)
   AC%VAL(1:NSIZE) = VAL(1:NSIZE)
   DEALLOCATE(VAL)

#endif

   WRITE(CINFO,'(A,A)') TRIM(CMATRIX),'.COL'
   
   CALL SCARC_ALLOCATE_INT1(COL, 1, NSIZE, NSCARC_INIT_NONE, CINFO)
   COL(1:NSIZE) = AC%COL(1:NSIZE)
   DEALLOCATE(AC%COL)

   CALL SCARC_ALLOCATE_INT1(AC%COL, 1, NSIZE, NSCARC_INIT_NONE, CINFO)
   AC%COL(1:NSIZE) = COL(1:NSIZE)
   DEALLOCATE(COL)

ELSE
   CALL SCARC_SHUTDOWN(NSCARC_ERROR_MATRIX_SIZE, SCARC_NONE, NSIZE)
ENDIF

END SUBROUTINE SCARC_RESIZE_MATRIX_COMPACT

!> ------------------------------------------------------------------------------------------------
!> Send timeout to MPI
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_TIMEOUT(RNAME,NR,RR)

REAL(EB) :: START_TIME,WAIT_TIME
INTEGER, INTENT(IN) :: NR
INTEGER, DIMENSION(:) :: RR,STATUSS(MPI_STATUS_SIZE)
LOGICAL :: FLAG,FLAG2,FLAG3
CHARACTER(*) :: RNAME
INTEGER :: NNN, IERR = 0

IF (.NOT.PROFILING) THEN

   ! Normally, PROFILING=F and this branch continually tests the communication and cancels the requests if too much time elapses.

   START_TIME = MPI_WTIME()
   FLAG = .FALSE.
   DO WHILE(.NOT.FLAG)
      CALL MPI_TESTALL(NR,RR(1:NR),FLAG,MPI_STATUSES_IGNORE,IERR)
      WAIT_TIME = MPI_WTIME() - START_TIME
      IF (WAIT_TIME>MPI_TIMEOUT) THEN
         WRITE(LU_ERR,'(A,A,I6)') TRIM(RNAME),' timed out for MPI process ',MYID
         FLAG = .TRUE.
         DO NNN=1,NR
            CALL MPI_CANCEL(RR(NNN),IERR)
            CALL MPI_TEST(RR(NNN),FLAG2,MPI_STATUS_IGNORE,IERR)
            CALL MPI_TEST_CANCELLED(STATUSS,FLAG3,IERR)
         ENDDO
      ENDIF
   ENDDO

ELSE

   ! If PROFILING=T, do not do MPI_TESTALL because too many calls to this routine swamps the tracing and profiling.

   CALL MPI_WAITALL(NR,RR(1:NR),MPI_STATUSES_IGNORE,IERR)

ENDIF

END SUBROUTINE SCARC_TIMEOUT


!> ----------------------------------------------------------------------------------------------------
!> Dump residual information
!> ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_DUMP_CSV(ISM, NS, NL)
INTEGER, INTENT(IN) :: ISM, NS, NL

IF (.NOT.HAS_CSV .OR. MYID /= 0) RETURN
IF (ITE_TOTAL == 0 .AND. TYPE_SOLVER /= NSCARC_SOLVER_MAIN) RETURN
WRITE(MSG%LU_STAT,1000) ITE_PRES, NS, ITE_TOTAL, ITE_CG, ITE_MG, NL, ITE_SMOOTH, ISM, ITE_COARSE, ITE_LU, RES, CAPPA

1000 FORMAT(10(I8,','), E12.4,',',E12.4)
END SUBROUTINE SCARC_DUMP_CSV



!> ================================================================================================
!> START DEBUGGING PART  
!> Collection of routines which print out different quantities or allow to preset them
!> Only enabled if directive SCARC_DEBUG is set
!> ================================================================================================
#ifdef WITH_SCARC_DEBUG

!> ------------------------------------------------------------------------------------------------
!> Preset right hand side in such a way that exact solution is known
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_PRESET_EXACT (NE, NL)
INTEGER, INTENT(IN) :: NE, NL
REAL (EB), POINTER, DIMENSION(:) :: VE
INTEGER :: IC, NM, I, K
REAL(EB), DIMENSION(:), POINTER :: XMID, ZMID
TYPE (MESH_TYPE), POINTER :: M=>NULL()
TYPE (SCARC_LEVEL_TYPE), POINTER :: L=>NULL()
TYPE (SCARC_DISCRET_TYPE), POINTER :: D=>NULL()

IF (ITE_TOTAL == 0) WRITE(*,*) 'CAUTION: PRESET_EXACT is active !!!'

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   M => MESHES(NM)
   L => SCARC(NM)%LEVEL(NL)
   SELECT CASE(TYPE_DISCRET)
      CASE (NSCARC_DISCRET_STRUCTURED)
         D => L%SD
      CASE (NSCARC_DISCRET_UNSTRUCTURED)
         D => L%UD
   END SELECT
   VE => POINT_TO_VECTOR (NM, NL, NE)

   DO K = 1, L%NZ
      DO I = 1, L%NX
         IF (TYPE_DISCRET == NSCARC_DISCRET_UNSTRUCTURED.AND.L%CELL_STATE(I,1,K) /= NSCARC_CELL_GASPHASE) CYCLE
         IC = D%CELL_NUMBER(I,1,K)
         IF (NL == NLEVEL_MIN) THEN
            XMID => M%XC
            ZMID => M%ZC
         ELSE
            XMID => L%XMID
            ZMID => L%ZMID
         ENDIF
         VE(IC) = EXACT(XMID(I),ZMID(K))
         WRITE(MSG%LU_DEBUG,'(A,i3,a,e10.2,a,e10.2,a,e12.4)') 'IC=',IC,':X=',XMID(i),':Z=',ZMID(k),': RHS=',VE(IC)
      ENDDO
   ENDDO
ENDDO

END SUBROUTINE SCARC_PRESET_EXACT

!> ------------------------------------------------------------------------------------------------
!> Preset vector 
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_PRESET_VECTOR (NV, NL)
INTEGER, INTENT(IN) :: NV, NL
REAL (EB), POINTER, DIMENSION(:) :: VEC
INTEGER :: IC, NM, I, K
REAL(EB), DIMENSION(:), POINTER :: XMID, ZMID
TYPE (MESH_TYPE), POINTER :: M=>NULL()
TYPE (SCARC_LEVEL_TYPE), POINTER :: L=>NULL()
TYPE (SCARC_DISCRET_TYPE), POINTER :: D=>NULL()

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   M => MESHES(NM)
   L => SCARC(NM)%LEVEL(NL)
   SELECT CASE(TYPE_DISCRET)
      CASE (NSCARC_DISCRET_STRUCTURED)
         D => L%SD
      CASE (NSCARC_DISCRET_UNSTRUCTURED)
         D => L%UD
   END SELECT
   VEC => POINT_TO_VECTOR (NM, NL, NV)

   DO K = 1, L%NZ
      DO I = 1, L%NX
         IF (TYPE_DISCRET == NSCARC_DISCRET_UNSTRUCTURED.AND.L%CELL_STATE(I,1,K) /= NSCARC_CELL_GASPHASE) CYCLE
         IC = D%CELL_NUMBER(I,1,K)
         IF (NL == NLEVEL_MIN) THEN
            XMID => M%XC
            ZMID => M%ZC
         ELSE
            XMID => L%XMID
            ZMID => L%ZMID
         ENDIF
         VEC(IC) = XMID(I)
         WRITE(MSG%LU_DEBUG,'(A,i3,a,e10.2,a,e10.2,a,e12.4)') 'IC=',IC,':X=',XMID,':Z=',ZMID,': RHS=',VEC(IC)
      ENDDO
   ENDDO
ENDDO

END SUBROUTINE SCARC_PRESET_VECTOR


!> ------------------------------------------------------------------------------------------------
!> Preset right hand side in such a way that exact solution is known
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_PRESET_RHS (NV, NL)
INTEGER, INTENT(IN) :: NV, NL
REAL (EB), POINTER, DIMENSION(:) :: VC
INTEGER :: IC, NM, I, K
REAL (EB) :: X, Z
TYPE (MESH_TYPE), POINTER :: M=>NULL()
TYPE (SCARC_LEVEL_TYPE), POINTER :: L=>NULL()
TYPE (SCARC_DISCRET_TYPE), POINTER :: D=>NULL()

IF (NL > NLEVEL_MIN) WRITE(*,*) 'Wrong level for presetting RHS '

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   M => MESHES(NM)
   M%BXS = 0.0_EB
   M%BXF = 0.0_EB
   M%BZS = 0.0_EB
   M%BZF = 0.0_EB

   L => SCARC(NM)%LEVEL(NL)
   SELECT CASE(TYPE_DISCRET)
      CASE (NSCARC_DISCRET_STRUCTURED)
         D => L%SD
      CASE (NSCARC_DISCRET_UNSTRUCTURED)
         D => L%UD
   END SELECT
   VC => POINT_TO_VECTOR (NM, NL, NV)

   DO K = 1, L%NZ
      DO I = 1, L%NX
         IF (TYPE_DISCRET == NSCARC_DISCRET_UNSTRUCTURED.AND.L%CELL_STATE(I,1,K) /= NSCARC_CELL_GASPHASE) CYCLE
         IC = D%CELL_NUMBER(I,1,K)
         X  = M%XC(I)
         Z  = M%ZC(K)
         WRITE(MSG%LU_DEBUG,'(A,i3,a,e10.2,a,e10.2,a,e12.4)') 'IC=',IC,':X=',X,':Z=',Z,': RHS=',VC(IC)
         VC(IC) = RHS(X,Z)
      ENDDO
   ENDDO
ENDDO

END SUBROUTINE SCARC_PRESET_RHS

!> ------------------------------------------------------------------------------------------------
!> Set exact solution 
!> ------------------------------------------------------------------------------------------------
DOUBLE PRECISION FUNCTION EXACT(X,Z)
REAL (EB), INTENT(IN) :: X, Z
!EXACT = (X**2 - X**4) * (Z**4 - Z**2)                                    !> FUNCTION 1
!EXACT = (X**2 - 1) * (Z**2 - 1)                                         !> FUNCTION 2
!EXACT =  625.0_EB/16.0_EB * X * (0.8_EB - X) * Z * (0.8_EB - Z)        !> FUNCTION 3
EXACT = - X * (0.8_EB - X) * Z * (0.8_EB - Z)        !> FUNCTION 3
END FUNCTION EXACT


!> ------------------------------------------------------------------------------------------------
!> Set right hand side
!> ------------------------------------------------------------------------------------------------
DOUBLE PRECISION FUNCTION RHS(X,Z)
REAL (EB), INTENT(IN) :: X, Z
!RHS = 2.0_EB*((1.0_EB - 6.0_EB*X**2)*Z**2*(1.0_EB-Z**2)+(1.0_EB-6.0_EB*Z**2)*X**2*(1.0_EB-X**2))
!RHS = -X**2 - Z**2 +2
!RHS = 625.0_EB/8.0_EB * (X * (0.8_EB - X) + Z * (0.8_EB - Z))
RHS = 2.0_EB * (X * (0.8_EB - X) + Z * (0.8_EB - Z))
END FUNCTION RHS

!> ------------------------------------------------------------------------------------------------
!> Save dump of vector in dump-directory
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_DUMP_QUANTITY (NV, CNAME, ISM, NS, NL)
INTEGER, INTENT(IN) :: NV, NS, NL, ISM
CHARACTER(*), INTENT(IN) :: CNAME
REAL (EB), POINTER, DIMENSION(:)     :: VEC
CHARACTER(80) :: FN_DUMP, CDIR
INTEGER :: LU_DUMP, IC, NM
TYPE (SCARC_LEVEL_TYPE), POINTER :: L=>NULL()
TYPE (SCARC_DISCRET_TYPE), POINTER :: D=>NULL()

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
   VEC => POINT_TO_VECTOR (NM, NL, NV)
   L => SCARC(NM)%LEVEL(NL)
   SELECT CASE(TYPE_DISCRET)
      CASE (NSCARC_DISCRET_STRUCTURED)
         D => L%SD
      CASE (NSCARC_DISCRET_UNSTRUCTURED)
         D => L%UD
   END SELECT
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
   DO IC = 1, D%NC
      !WRITE(LU_DUMP,'(F25.16)') VEC(IC)
      WRITE(LU_DUMP,*) VEC(IC)
   ENDDO
   CLOSE(LU_DUMP)
ENDDO

FN_DUMP = CNAME
END SUBROUTINE SCARC_DUMP_QUANTITY


!> ------------------------------------------------------------------------------------------------
!> Print out vector information on level NL
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_DEBUG_LEVEL (NV, CVEC, NL)
INTEGER, INTENT(IN) :: NV, NL
REAL (EB) :: VALUES(0:100)
INTEGER :: NM, II, JJ, KK, IC, NX8, NY8, NZ8
CHARACTER (*), INTENT(IN) :: CVEC
REAL (EB), POINTER, DIMENSION(:)     :: VC
TYPE (SCARC_LEVEL_TYPE), POINTER :: L=>NULL()
TYPE (SCARC_DISCRET_TYPE), POINTER :: D=>NULL()

!IF (TYPE_SOLVER /= NSCARC_SOLVER_MAIN) RETURN
DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   L  => SCARC(NM)%LEVEL(NL)
   SELECT CASE(TYPE_DISCRET)
      CASE (NSCARC_DISCRET_STRUCTURED)
         D => L%SD
      CASE (NSCARC_DISCRET_UNSTRUCTURED)
         D => L%UD
   END SELECT
   VC => POINT_TO_VECTOR (NM, NL, NV)

   NX8=MIN(12,L%NX)
   NY8=MIN(12,L%NY)
   NZ8=MIN(12,L%NZ)

   WRITE(MSG%LU_DEBUG,*) '=========================================================='
   WRITE(MSG%LU_DEBUG,2001) CVEC, NM, NL
   WRITE(MSG%LU_DEBUG,2002) D%NC, NX8, NY8, NZ8, NV, SIZE(VC)
   WRITE(MSG%LU_DEBUG,*) '=========================================================='
   !IF (NL == NLEVEL_MIN) THEN
         DO KK = NZ8, 1, - 1
            DO JJ = NY8, 1, - 1
               DO II=1,NX8
                  IF (L%CELL_STATE(II,JJ,KK) == NSCARC_CELL_GASPHASE) THEN
                     IC=D%CELL_NUMBER(II,JJ,KK)
                     IF (ABS(VC(IC))<1.0E-14_EB) THEN
                        VALUES(II)=0.0_EB
                     ELSE
                        VALUES(II)=VC(IC)
                     ENDIF
                  ELSE
                     VALUES(II)=-999999.0_EB
                  ENDIF
               ENDDO
               WRITE(MSG%LU_DEBUG, '(12E14.6)') (VALUES(II), II=1, NX8)
            ENDDO
         ENDDO
   !ENDIF
ENDDO

!CALL SCARC_MATLAB_VECTOR(NV, CVEC, NL)

!2000 FORMAT('=== ',A,' : ',A,' on mesh ',I8,' on level ',I8, ': NX, NY, NZ=',3I8,': NV=',I8)
2001 FORMAT('=== ',A,' on mesh ',I8,' on level ',I8)
2002 FORMAT('=== NC = ',I4, ': NX, NY, NZ=',3I4,': NV=',I3,': Size=',I8)
END SUBROUTINE SCARC_DEBUG_LEVEL

!> ------------------------------------------------------------------------------------------------
!> Debug specified quantity
!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_DEBUG_QUANTITY(NTYPE, NL, CQUANTITY)
INTEGER, INTENT(IN) :: NTYPE, NL
INTEGER :: NM, NOM, IP, IC, ID, IW, I, J, IOR0, INBR, III, JJJ, KKK, IWG
CHARACTER (*), INTENT(IN) :: CQUANTITY
TYPE (MESH_TYPE), POINTER :: M=>NULL()
TYPE (SCARC_LEVEL_TYPE), POINTER :: L=>NULL(), OL=>NULL()
TYPE (SCARC_POINTER_TYPE), POINTER :: PTR=>NULL(), OPTR=>NULL()
TYPE (SCARC_DISCRET_TYPE), POINTER :: D=>NULL(), OD=>NULL()
TYPE (SCARC_SOLVER_TYPE), POINTER :: SV=>NULL()
TYPE (SCARC_MATRIX_COMPACT_TYPE), POINTER :: AC=>NULL()
TYPE (SCARC_MATRIX_BANDED_TYPE), POINTER :: AB=>NULL()
#ifdef WITH_MKL
TYPE (SCARC_MATRIX_COMPACT_TYPE), POINTER :: AC_SYM=>NULL()
#endif

SELECT CASE (NTYPE)

   !> ------------------------------------------------------------------------------------------------
   !> Debug system matrix A (corresponding to system type)
   !> ------------------------------------------------------------------------------------------------
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
         WRITE(MSG%LU_DEBUG,*) 'TYPE_SCOPE    = ',SV%TYPE_SCOPE 
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


   !> ------------------------------------------------------------------------------------------------
   !> Debug system matrix A (corresponding to system type)
   !> ------------------------------------------------------------------------------------------------
   CASE (NSCARC_DEBUG_MATRIX)

      SELECT CASE (TYPE_MATRIX)

         CASE (NSCARC_MATRIX_COMPACT)

            DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
               L  => SCARC(NM)%LEVEL(NL)
               SELECT CASE(TYPE_DISCRET)
                  CASE (NSCARC_DISCRET_STRUCTURED)
                     D => L%SD
                  CASE (NSCARC_DISCRET_UNSTRUCTURED)
                     D => L%UD
               END SELECT
               AC => D%AC
               WRITE(MSG%LU_DEBUG,1000) CQUANTITY, NM, NL
               WRITE(MSG%LU_DEBUG,*) '----------- SHOWING FULL COMPACT MATRIX ENTRIES'
               WRITE(MSG%LU_DEBUG,*) 'NC =',L%NC
               WRITE(MSG%LU_DEBUG,*) 'NV =',AC%NA
               WRITE(MSG%LU_DEBUG,*) 'NR =',AC%NR
               WRITE(MSG%LU_DEBUG,*) 'SIZE(AC%VAL) =',SIZE(AC%VAL)
               WRITE(MSG%LU_DEBUG,*) 'SIZE(AC%COL) =',SIZE(AC%COL)
               WRITE(MSG%LU_DEBUG,*) 'SIZE(AC%ROW) =',SIZE(AC%ROW)
               WRITE(MSG%LU_DEBUG,*) '---------------------- AC%POS:'
               WRITE(MSG%LU_DEBUG,'(7i8)') (AC%POS(IC), IC=-3,3)
               WRITE(MSG%LU_DEBUG,*) '---------------------- AC%STENCIL:'
               WRITE(MSG%LU_DEBUG,'(7F8.1)') (AC%STENCIL(IC), IC=-3,3)
               WRITE(MSG%LU_DEBUG,*) '---------------------- AC%ROW:'
               WRITE(MSG%LU_DEBUG,'(7i8)') (AC%ROW(IC), IC=1,AC%NR)
               WRITE(MSG%LU_DEBUG,*) '---------------------- AC%COL:'
               DO IC = 1, AC%NR-1
                  WRITE(MSG%LU_DEBUG,'(i5,a,20i9)') IC,':',(AC%COL(IP),IP=AC%ROW(IC),AC%ROW(IC+1)-1)
               ENDDO
               WRITE(MSG%LU_DEBUG,*) '---------------------- A:'
               DO IC = 1, AC%NR-1
                  WRITE(MSG%LU_DEBUG,'(i5,a,20f8.1)') IC,':',(AC%VAL(IP),IP=AC%ROW(IC),AC%ROW(IC+1)-1)
               ENDDO
               IF (ALLOCATED(AC%GSM)) THEN
                  WRITE(MSG%LU_DEBUG,*) '---------------------- GSM:'
                  DO IC = 1, AC%NR-1
                     WRITE(MSG%LU_DEBUG,'(i5,a,20f15.6)') IC,':',(AC%GSM(IP),IP=AC%ROW(IC),AC%ROW(IC+1)-1)
                  ENDDO
               ENDIF
               IF (ALLOCATED(AC%SGSM)) THEN
                  WRITE(MSG%LU_DEBUG,*) '---------------------- SGSM:'
                  DO IC = 1, AC%NR-1
                     WRITE(MSG%LU_DEBUG,'(i5,a,20f15.6)') IC,':',(AC%SGSM(IP),IP=AC%ROW(IC),AC%ROW(IC+1)-1)
                  ENDDO
               ENDIF
               IF (ALLOCATED(AC%SORM)) THEN
                  WRITE(MSG%LU_DEBUG,*) '---------------------- SORM:'
                  DO IC = 1, AC%NR-1
                     WRITE(MSG%LU_DEBUG,'(i5,a,20f15.6)') IC,':',(AC%SORM(IP),IP=AC%ROW(IC),AC%ROW(IC+1)-1)
                  ENDDO
               ENDIF
               IF (ALLOCATED(AC%SSORM)) THEN
                  WRITE(MSG%LU_DEBUG,*) '---------------------- SSORM:'
                  DO IC = 1, AC%NR-1
                     WRITE(MSG%LU_DEBUG,'(i5,a,20f15.6)') IC,':',(AC%SSORM(IP),IP=AC%ROW(IC),AC%ROW(IC+1)-1)
                  ENDDO
               ENDIF
               IF (ALLOCATED(AC%ILU)) THEN
                  WRITE(MSG%LU_DEBUG,*) '---------------------- ILU:'
                  DO IC = 1, AC%NR-1
                     WRITE(MSG%LU_DEBUG,'(i5,a,20f15.6)') IC,':',(AC%ILU(IP),IP=AC%ROW(IC),AC%ROW(IC+1)-1)
                  ENDDO
               ENDIF
#ifdef WITH_MKL
               !IF ((TYPE_METHOD == NSCARC_METHOD_LU) .AND. (TYPE_MKL == NSCARC_MKL_GLOBAL)) THEN
               !   WRITE(MSG%LU_DEBUG,*) '---------------------- AG_COL:'
               !   DO IC = 1, L%NC
               !      WRITE(MSG%LU_DEBUG,'(i5,a,20i9)') IC,':',(AC%COLG(IP),IP=AC%ROW(IC),AC%ROW(IC+1)-1)
               !      WRITE(MSG%LU_DEBUG,*)  IC,':',(AC%COLG(IP),IP=AC%ROW(IC),AC%ROW(IC+1)-1)
               !   ENDDO
               !ENDIF
               !DO IC=1,AC%NR-1
               !   DO IP=AC%ROW(IC),AC%ROW(IC+1)-1
               !       WRITE(MSG%LU_DEBUG,'(2I8,F24.12)') IC,AC%COL(IP),AC%VAL(IP)
               !       WRITE(MSG%LU_DEBUG,*) IC,AC%COL(IP),AC%VAL(IP)
               !   ENDDO
               !   WRITE(MSG%LU_DEBUG,*)
               !ENDDO
               !DO IC=1,AC%NR-1
               !   DO IP=AC%ROW(IC),AC%ROW(IC+1)-1
               !       IF (IC == AC%COL(IP)) WRITE(MSG%LU_DEBUG,'(2I8,F24.12)') IC,AC%COL(IP),AC%VAL(IP)
               !   ENDDO
               !ENDDO
#endif
               CALL SCARC_MATLAB_MATRIX(AC%VAL, AC%ROW, AC%COL, D%NC, D%NC, NM, NL, 'A')

            ENDDO

         CASE (NSCARC_MATRIX_BANDED)

            DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
               L => SCARC(NM)%LEVEL(NL)
               SELECT CASE(TYPE_DISCRET)
                  CASE (NSCARC_DISCRET_STRUCTURED)
                     AB => L%SD%AB
                  CASE (NSCARC_DISCRET_UNSTRUCTURED)
                     AB => L%UD%AB
               END SELECT
               WRITE(MSG%LU_DEBUG,1000) CQUANTITY, NM, NL
               WRITE(MSG%LU_DEBUG,*) '----------- SHOWING FULL BANDED MATRIX ENTRIES'
               WRITE(MSG%LU_DEBUG,*) 'NC =',L%NC
               WRITE(MSG%LU_DEBUG,*) 'NA  =',AB%NA
               WRITE(MSG%LU_DEBUG,*) 'NDIAG =',AB%NDIAG
               WRITE(MSG%LU_DEBUG,*) 'SIZE(AB%VAL) =',SIZE(AB%VAL)
               WRITE(MSG%LU_DEBUG,*) '---------------------- AB%POS:'
               WRITE(MSG%LU_DEBUG,'(7i8)') (AB%POS(IC), IC=-3,3)
               WRITE(MSG%LU_DEBUG,*) '---------------------- AB%STENCIL:'
               WRITE(MSG%LU_DEBUG,'(7F8.1)') (AB%STENCIL(IC), IC=-3,3)
               WRITE(MSG%LU_DEBUG,*) '---------------------- AB%OFFSET:'
               WRITE(MSG%LU_DEBUG,'(7i8)') (AB%OFFSET(IC), IC=-3,3)
               WRITE(MSG%LU_DEBUG,*) '---------------------- AB%VAL:'
               DO ID = 1, AB%NSTENCIL
                  WRITE(MSG%LU_DEBUG,*) ID,':',(AB%VAL(ID, IC),IC=1,AB%NDIAG)
               ENDDO
            ENDDO

      END SELECT


   !> ------------------------------------------------------------------------------------------------
   !> Debug symmetric system matrix AS
   !> ------------------------------------------------------------------------------------------------
   CASE (NSCARC_DEBUG_MATRIX_SYM)

      IF (NL > NLEVEL_MIN) RETURN

#ifdef WITH_MKL
      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
         WRITE(MSG%LU_DEBUG,1000) CQUANTITY, NM, NL
         L   => SCARC(NM)%LEVEL(NL)
         SELECT CASE(TYPE_DISCRET)
            CASE (NSCARC_DISCRET_STRUCTURED)
               AC_SYM => L%SD%AC_SYM
            CASE (NSCARC_DISCRET_UNSTRUCTURED)
               AC_SYM => L%UD%AC_SYM
         END SELECT
         WRITE(MSG%LU_DEBUG,*) '----------- SHOWING SYMMETRIC MATRIX ENTRIES'
         WRITE(MSG%LU_DEBUG,*) 'L%NC=', L%NC
         WRITE(MSG%LU_DEBUG,*) 'AC_SYM%ROW(L%NC)=', AC_SYM%ROW(L%NC)
         WRITE(MSG%LU_DEBUG,*) 'AC_SYM%ROW(L%NC+1)=', AC_SYM%ROW(L%NC+1)
         WRITE(MSG%LU_DEBUG,*) 'SIZE(ASC) =',    SIZE(AC_SYM%VAL)
         WRITE(MSG%LU_DEBUG,*) 'SIZE(AC_SYM%COL) =',SIZE(AC_SYM%COL)
         WRITE(MSG%LU_DEBUG,*) 'SIZE(AC_SYM%ROW) =',SIZE(AC_SYM%ROW)
         WRITE(MSG%LU_DEBUG,*) '---------------------- AC_SYM%ROW:', L%NC
         WRITE(MSG%LU_DEBUG,*) (AC_SYM%ROW(IC), IC=1,L%NC+1)
         WRITE(MSG%LU_DEBUG,*) '---------------------- AC_SYM%COL:'
         DO IC = 1, L%NC
            WRITE(MSG%LU_DEBUG,*) IC,':',(AC_SYM%COL(IP),IP=AC_SYM%ROW(IC),AC_SYM%ROW(IC+1)-1)
         ENDDO
         DO IC = 1, L%NC
            WRITE(MSG%LU_DEBUG,*) IC,':',(AC_SYM%VAL(IP),IP=AC_SYM%ROW(IC),AC_SYM%ROW(IC+1)-1)
         ENDDO
         DO IC=1,L%NC
            DO IP=AC_SYM%ROW(IC),AC_SYM%ROW(IC+1)-1
                WRITE(MSG%LU_DEBUG,*) IC,AC_SYM%COL(IP),AC_SYM%VAL(IP)
            ENDDO
            WRITE(MSG%LU_DEBUG,*)
         ENDDO
            WRITE(MSG%LU_DEBUG,*)
            WRITE(MSG%LU_DEBUG,*)
         DO IC=1,L%NC
            DO IP=AC_SYM%ROW(IC),AC_SYM%ROW(IC+1)-1
                IF (IC == AC_SYM%COL(IP)) WRITE(MSG%LU_DEBUG,'(2I8,F24.12)') IC,AC_SYM%COL(IP),AC_SYM%VAL(IP)
            ENDDO
            WRITE(MSG%LU_DEBUG,*)
         ENDDO

      ENDDO
#endif

 !> ------------------------------------------------------------------------------------------------
 !> Debug FACEINFO
 !> ------------------------------------------------------------------------------------------------
   CASE (NSCARC_DEBUG_FACE)

      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
         L => SCARC(NM)%LEVEL(NL)
         WRITE(MSG%LU_DEBUG,1000) CQUANTITY, NM, NL
         DO IOR0 = -3, 3
            IF (IOR0 == 0) CYCLE
            WRITE(MSG%LU_DEBUG,*) '========================================='
            WRITE(MSG%LU_DEBUG,*) '============= DEBUGGING FACE(',IOR0,'): FOR LEVEL ', NL
            WRITE(MSG%LU_DEBUG,*) '========================================='
            WRITE(MSG%LU_DEBUG,*) 'FACE(.)%NFX:', L%FACE(IOR0)%NFX
            WRITE(MSG%LU_DEBUG,*) 'FACE(.)%NFY:', L%FACE(IOR0)%NFY
            WRITE(MSG%LU_DEBUG,*) 'FACE(.)%NFZ:', L%FACE(IOR0)%NFZ
            WRITE(MSG%LU_DEBUG,*) 'FACE(.)%NWL:', L%FACE(IOR0)%NFW
            WRITE(MSG%LU_DEBUG,*) 'FACE(.)%IWG_PTR:', L%FACE(IOR0)%IWG_PTR
            WRITE(MSG%LU_DEBUG,*) '----------------------------------------------'
            WRITE(MSG%LU_DEBUG,*)
            !WRITE(MSG%LU_DEBUG,*) 'FACE(.)%DH :', L%FACE(IOR0)%DH
            WRITE(MSG%LU_DEBUG,*)
            IF (L%FACE(IOR0)%N_NEIGHBORS /= 0) THEN
               DO INBR=1,L%FACE(IOR0)%N_NEIGHBORS
                  NOM = L%FACE(IOR0)%NEIGHBORS(INBR)
                  OL => SCARC(NM)%OSCARC(NOM)%LEVEL(NL)
                  SELECT CASE(TYPE_DISCRET)
                     CASE (NSCARC_DISCRET_STRUCTURED)
                        OD => OL%SD
                     CASE (NSCARC_DISCRET_UNSTRUCTURED)
                        OD => OL%UD
                  END SELECT
                  OPTR => OD%PTR
                  WRITE(MSG%LU_DEBUG,*) 'N_NEIGHBORS:', L%FACE(IOR0)%N_NEIGHBORS
                  WRITE(MSG%LU_DEBUG,*) 'NOM:', NOM
                  WRITE(MSG%LU_DEBUG,*) 'SIZE(OPTR%ICG_TO_IWG)=',SIZE(OPTR%ICG_TO_IWG)
                  WRITE(MSG%LU_DEBUG,*) 'SIZE(OPTR%ICG_TO_ICE)=',SIZE(OPTR%ICG_TO_ICE)
                  WRITE(MSG%LU_DEBUG,*) 'SIZE(OPTR%ICG_TO_ICO)=',SIZE(OPTR%ICG_TO_ICO)
                  !WRITE(MSG%LU_DEBUG,*) 'SIZE(OPTR%IWL_TO_ICW)=',SIZE(OPTR%IWL_TO_ICW)     !> AMG only
                  WRITE(MSG%LU_DEBUG,*) 'SIZE(OPTR%IWL_TO_ICO)=',SIZE(OPTR%IWL_TO_ICO)
                  !IF (BAMG) WRITE(MSG%LU_DEBUG,*) 'SIZE(OPTR%IWL_TO_ICG)=',SIZE(OPTR%IWL_TO_ICG)
                  WRITE(MSG%LU_DEBUG,'(a,i8,a,2f12.6)') '---OC(',NOM,')%DH :',OL%DH
                  WRITE(MSG%LU_DEBUG,'(a,i8,a,2i8)') '---OL(',NOM,')%IOR   :',OL%IOR
                  WRITE(MSG%LU_DEBUG,'(a,i8,a,2i8)') '---OL(',NOM,')%NWL(.):',OL%NWL
                  WRITE(MSG%LU_DEBUG,'(a,i8,a,2i8)') '---OL(',NOM,')%NCG(.):',OL%NCG
                  WRITE(MSG%LU_DEBUG,*)
                  WRITE(MSG%LU_DEBUG,'(a,i8,a)') '------OL(',NOM,')%ICG_TO_IWG:'
                  DO IW = 1, OL%NCG
                     WRITE(MSG%LU_DEBUG,'(16i8)') OPTR%ICG_TO_IWG(IW)
                  ENDDO
                  WRITE(MSG%LU_DEBUG,*)
                  WRITE(MSG%LU_DEBUG,'(a,i8,a)') '------OL(',NOM,')%ICG_TO_ICO:'
                  DO IW = 1, OL%NCG
                     WRITE(MSG%LU_DEBUG,'(16i8)') OPTR%ICG_TO_ICO(IW)
                  ENDDO
                  WRITE(MSG%LU_DEBUG,*)
                  WRITE(MSG%LU_DEBUG,'(a,i8,a)') '------OL(',NOM,')%ICG_TO_ICE:'
                  DO IW = 1, OL%NCG
                     WRITE(MSG%LU_DEBUG,'(16i8)') OPTR%ICG_TO_ICE(IW)
                  ENDDO
                  WRITE(MSG%LU_DEBUG,*)
                  WRITE(MSG%LU_DEBUG,*)
                  WRITE(MSG%LU_DEBUG,'(a,i8,a)') '------OL(',NOM,')%IWL_TO_IWG:'
                  DO IW = 1, OL%NWL
                     WRITE(MSG%LU_DEBUG,'(16i8)') OPTR%IWL_TO_IWG(IW)
                  ENDDO
                  !WRITE(MSG%LU_DEBUG,*)
                  !WRITE(MSG%LU_DEBUG,'(a,i8,a)') '------OL(',NOM,')%IWL_TO_ICW:'
                  !DO IW = 1, OL%NWL
                  !   WRITE(MSG%LU_DEBUG,'(16i8)') OPTR%IWL_TO_ICW(IW)
                  !ENDDO
                  WRITE(MSG%LU_DEBUG,*)
                  WRITE(MSG%LU_DEBUG,'(a,i8,a)') '------OL(',NOM,')%IWL_TO_ICO:'
                  DO IW = 1, OL%NWL
                     WRITE(MSG%LU_DEBUG,'(16i8)') OPTR%IWL_TO_ICO(IW)
                  ENDDO
                  !IF (BAMG) THEN
                  !   WRITE(MSG%LU_DEBUG,*)
                  !   WRITE(MSG%LU_DEBUG,'(a,i4,a)') '------OL(',NOM,')%IWL_TO_ICG:'
                  !   DO IW = 1, OPTR%NWL
                  !      WRITE(MSG%LU_DEBUG,'(32i8)') OPTR%IWL_TO_ICG(IW, 1:OPTR%NCPL)
                  !   ENDDO
                  !ENDIF
               ENDDO
            ENDIF
         ENDDO
      ENDDO

 !> ------------------------------------------------------------------------------------------------
 !> Debug WALLINFO
 !> ------------------------------------------------------------------------------------------------
   CASE (NSCARC_DEBUG_WALL)

      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
         M => MESHES(NM)
         L => SCARC(NM)%LEVEL(NL)
         SELECT CASE(TYPE_DISCRET)
            CASE (NSCARC_DISCRET_STRUCTURED)
               D => L%SD
            CASE (NSCARC_DISCRET_UNSTRUCTURED)
               D => L%UD
         END SELECT
         PTR => D%PTR
         WRITE(MSG%LU_DEBUG,1000) CQUANTITY, NM, NL
         WRITE(MSG%LU_DEBUG,*) 'SIZE(PTR%ICE_TO_IWG)=',SIZE(PTR%ICE_TO_IWG)
         WRITE(MSG%LU_DEBUG,*) 'SIZE(PTR%ICE_TO_IWL)=',SIZE(PTR%ICE_TO_IWL)
         WRITE(MSG%LU_DEBUG,*) 'SIZE(PTR%ICE_TO_ICG)=',SIZE(PTR%ICE_TO_ICG)
         WRITE(MSG%LU_DEBUG,*) 'NM  =',NM
         WRITE(MSG%LU_DEBUG,*) 'NL  =',NL
         WRITE(MSG%LU_DEBUG,*) 'NC  =',L%NC
         WRITE(MSG%LU_DEBUG,*) 'NCE =',L%NCE
         WRITE(MSG%LU_DEBUG,*) 'NW  =',L%N_WALL_CELLS
         WRITE(MSG%LU_DEBUG,*)
         WRITE(MSG%LU_DEBUG,'(a,I8,a)') '------PTR%ICE_TO_IWG:'
         WRITE(MSG%LU_DEBUG,'(8I8)') (PTR%ICE_TO_IWG(IW), IW = L%NC+1, L%NCE)
         WRITE(MSG%LU_DEBUG,'(a,I8,a)') '------PTR%ICE_TO_IWL:'
         WRITE(MSG%LU_DEBUG,'(8I8)') (PTR%ICE_TO_IWL(IW), IW = L%NC+1, L%NCE)
         WRITE(MSG%LU_DEBUG,*)
         WRITE(MSG%LU_DEBUG,'(a,I8,a)') '------PTR%ICE_TO_ICG:'
         WRITE(MSG%LU_DEBUG,'(8I8)') (PTR%ICE_TO_ICG(IW), IW = L%NC+1, L%NCE)
         IF (N_DIRIC_GLOBAL(NLEVEL_MIN) == 0) THEN
            WRITE(MSG%LU_DEBUG,'(a,I8,a)') '------PTR%ICE_TO_ICN:'
            WRITE(MSG%LU_DEBUG,'(8I8)') (PTR%ICE_TO_ICN(IW), IW = L%NC+1, L%NCE)
            WRITE(MSG%LU_DEBUG,'(a,I8,a)') '------PTR%ICE_TO_VAL:'
            WRITE(MSG%LU_DEBUG,'(8e12.4)') (PTR%ICE_TO_VAL(IW), IW = L%NC+1, L%NCE)
         ENDIF
         WRITE(MSG%LU_DEBUG,*)
         IF (NL == 1) THEN
         WRITE(MSG%LU_DEBUG,*) '============= WALL(.)%IXG:', NM
         WRITE(MSG%LU_DEBUG,'(16I6)') (D%WALL(IW)%IXG, IW=1,L%N_WALL_CELLS)
         WRITE(MSG%LU_DEBUG,*) '============= WALL(.)%IYG:', NM
         WRITE(MSG%LU_DEBUG,'(16I6)') (D%WALL(IW)%IYG, IW=1,L%N_WALL_CELLS)
         WRITE(MSG%LU_DEBUG,*) '============= WALL(.)%IZG:', NM
         WRITE(MSG%LU_DEBUG,'(16I6)') (D%WALL(IW)%IZG, IW=1,L%N_WALL_CELLS)
         WRITE(MSG%LU_DEBUG,*)
         WRITE(MSG%LU_DEBUG,*) '============= WALL(.)%IXW:', NM
         WRITE(MSG%LU_DEBUG,'(16I6)') (D%WALL(IW)%IXW, IW=1,L%N_WALL_CELLS)
         WRITE(MSG%LU_DEBUG,*) '============= WALL(.)%IYW:', NM
         WRITE(MSG%LU_DEBUG,'(16I6)') (D%WALL(IW)%IYW, IW=1,L%N_WALL_CELLS)
         WRITE(MSG%LU_DEBUG,*) '============= WALL(.)%IZW:', NM
         WRITE(MSG%LU_DEBUG,'(16I6)') (D%WALL(IW)%IZW, IW=1,L%N_WALL_CELLS)
         WRITE(MSG%LU_DEBUG,*)
         WRITE(MSG%LU_DEBUG,*) '============= WALL(.)%IXN(1):', NM
         WRITE(MSG%LU_DEBUG,'(16I6)') (D%WALL(IW)%IXN(1), IW=1,L%N_WALL_CELLS)
         WRITE(MSG%LU_DEBUG,*) '============= WALL(.)%IYN(1):', NM
         WRITE(MSG%LU_DEBUG,'(16I6)') (D%WALL(IW)%IYN(1), IW=1,L%N_WALL_CELLS)
         WRITE(MSG%LU_DEBUG,*) '============= WALL(.)%IZN(1):', NM
         WRITE(MSG%LU_DEBUG,'(16I6)') (D%WALL(IW)%IZN(1), IW=1,L%N_WALL_CELLS)
         WRITE(MSG%LU_DEBUG,*)
         WRITE(MSG%LU_DEBUG,*) '============= WALL(.)%IXN(2):', NM
         WRITE(MSG%LU_DEBUG,'(16I6)') (D%WALL(IW)%IXN(2), IW=1,L%N_WALL_CELLS)
         WRITE(MSG%LU_DEBUG,*) '============= WALL(.)%IYN(2):', NM
         WRITE(MSG%LU_DEBUG,'(16I6)') (D%WALL(IW)%IYN(2), IW=1,L%N_WALL_CELLS)
         WRITE(MSG%LU_DEBUG,*) '============= WALL(.)%IZN(2):', NM
         WRITE(MSG%LU_DEBUG,'(16I6)') (D%WALL(IW)%IZN(2), IW=1,L%N_WALL_CELLS)
         ENDIF
         WRITE(MSG%LU_DEBUG,*)
         WRITE(MSG%LU_DEBUG,*) '============= WALL(.)%BTYPE:', NM
         WRITE(MSG%LU_DEBUG,'(16I6)') (D%WALL(IW)%BTYPE, IW=1,L%N_WALL_CELLS)
         WRITE(MSG%LU_DEBUG,*)
         WRITE(MSG%LU_DEBUG,*) '============= WALL(.)%IOR:', NM
         WRITE(MSG%LU_DEBUG,'(16I6)') (D%WALL(IW)%IOR, IW=1,L%N_WALL_CELLS)
         WRITE(MSG%LU_DEBUG,*)
         WRITE(MSG%LU_DEBUG,*) '============= WALL(.)%NOM:', NM
         WRITE(MSG%LU_DEBUG,'(16I6)') (D%WALL(IW)%NOM, IW=1,L%N_WALL_CELLS)
         WRITE(MSG%LU_DEBUG,*)
         !IF (NL == 1) THEN
         !WRITE(MSG%LU_DEBUG,*) '============= WALL(.)%NCPL:', NM
         !WRITE(MSG%LU_DEBUG,'(16i10)') (D%WALL(IW)%NCPL, IW=1,L%N_WALL_CELLS)
         !WRITE(MSG%LU_DEBUG,*) (D%WALL(IW)%NCPL, IW=1,L%N_WALL_CELLS)
         !ENDIF
         WRITE(MSG%LU_DEBUG,*)
         WRITE(MSG%LU_DEBUG,*) '============= WALL(.)%ICW:', NM
         DO IW=1,L%N_WALL_CELLS
            WRITE(MSG%LU_DEBUG,'(a,I8, a,I8)') 'IW=',IW,':',D%WALL(IW)%ICW
         ENDDO
         WRITE(MSG%LU_DEBUG,*)
         WRITE(MSG%LU_DEBUG,*) '============= WALL(.)%ICO:', NM
         DO IW=1,L%N_WALL_CELLS
            IF (D%WALL(IW)%NOM /=0) &
               WRITE(MSG%LU_DEBUG,'(a,I8, a,I8)') 'IW=',IW,':',D%WALL(IW)%ICO
         ENDDO
         WRITE(MSG%LU_DEBUG,*)
         WRITE(MSG%LU_DEBUG,*) '============= WALL(.)%ICE:', NM
         DO IW=1,L%N_WALL_CELLS
            IF (D%WALL(IW)%NOM /=0) &
               WRITE(MSG%LU_DEBUG,'(a,I8, a,I8)') 'IW=',IW,':',D%WALL(IW)%ICE
         ENDDO
         WRITE(MSG%LU_DEBUG,*)
         WRITE(MSG%LU_DEBUG,*)
         WRITE(MSG%LU_DEBUG,*) '============= WALL(.)%ICG:', NM
         DO IW=1,L%N_WALL_CELLS
            IF (D%WALL(IW)%NOM /=0) &
               WRITE(MSG%LU_DEBUG,'(a,I8, a,I8)') 'IW=',IW,':',D%WALL(IW)%ICG
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

 !> ------------------------------------------------------------------------------------------------
 !> Debug complete grid information
 !> ------------------------------------------------------------------------------------------------
   CASE (NSCARC_DEBUG_DISCRET)
      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
         M => MESHES(NM)
         L => SCARC(NM)%LEVEL(NL)
         SELECT CASE(TYPE_DISCRET)
            CASE (NSCARC_DISCRET_STRUCTURED)
               D => L%SD
            CASE (NSCARC_DISCRET_UNSTRUCTURED)
               D => L%UD
         END SELECT
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
  
         WRITE(MSG%LU_DEBUG,*) 'M%WALL_INDEX:', CELL_COUNT(NM)
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
               M%WALL(IWG)%BOUNDARY_TYPE, D%WALL(IWG)%BTYPE, M%WALL(IWG)%PRESSURE_BC_INDEX
         ENDDO
         WRITE(MSG%LU_DEBUG,*) 'GRID%CELL_NUMBER(...)'
         DO JJJ=L%NY+1,0,-1
            WRITE(MSG%LU_DEBUG,*) ' ------------- JJJ = ', JJJ
            DO KKK = L%NZ+1,0,-1
               WRITE(MSG%LU_DEBUG,*) (D%CELL_NUMBER(III,JJJ,KKK),III=0,L%NX+1)
            ENDDO
         ENDDO
         WRITE(MSG%LU_DEBUG,*) 'GRIL%CELL_STATE(...)'
         DO JJJ=L%NY+1,0,-1
            WRITE(MSG%LU_DEBUG,*) ' ------------- JJJ = ', JJJ
            DO KKK = L%NZ+1,0,-1
               WRITE(MSG%LU_DEBUG,*) (L%CELL_STATE(III,JJJ,KKK),III=0,L%NX+1)
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
    
 !> ------------------------------------------------------------------------------------------------
 !> Debug PRESSURE_BC_INDEX
 !> ------------------------------------------------------------------------------------------------
   CASE (NSCARC_DEBUG_BDRY)

      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
         L => SCARC(NM)%LEVEL(NL)
         WRITE(MSG%LU_DEBUG,1000) CQUANTITY, NM, NL
         WRITE(MSG%LU_DEBUG, '(16I6)') (D%WALL(J)%BTYPE, J=1,L%N_WALL_CELLS)
      ENDDO

 !> ------------------------------------------------------------------------------------------------
 !> Debug SUBDIVISION
 !> ------------------------------------------------------------------------------------------------
   CASE (NSCARC_DEBUG_SUBDIVISION)

      DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
         L => SCARC(NM)%LEVEL(NL)
         WRITE(MSG%LU_DEBUG,1000) CQUANTITY, NM, NL
         WRITE(MSG%LU_DEBUG,*) 'SUBDIVISION IOR= 1 '
         WRITE(MSG%LU_DEBUG,'(3i8)') (L%SUBDIVISION(I, 1), I=1,3)
         WRITE(MSG%LU_DEBUG,*) 'SUBDIVISION IOR=-1 '
         WRITE(MSG%LU_DEBUG,'(3i8)') (L%SUBDIVISION(I,-1), I=1,3)
         WRITE(MSG%LU_DEBUG,*) 'SUBDIVISION IOR= 2 '
         WRITE(MSG%LU_DEBUG,'(3i8)') (L%SUBDIVISION(I, 2), I=1,3)
         WRITE(MSG%LU_DEBUG,*) 'SUBDIVISION IOR=-2 '
         WRITE(MSG%LU_DEBUG,'(3i8)') (L%SUBDIVISION(I,-2), I=1,3)
         WRITE(MSG%LU_DEBUG,*) 'SUBDIVISION IOR= 3 '
         WRITE(MSG%LU_DEBUG,'(3i8)') (L%SUBDIVISION(I, 3), I=1,3)
         WRITE(MSG%LU_DEBUG,*) 'SUBDIVISION IOR=-3 '
         WRITE(MSG%LU_DEBUG,'(3i8)') (L%SUBDIVISION(I,-3), I=1,3)
      ENDDO

END SELECT

1000 FORMAT('======================================================================================',/, &
            '=== ', A30,' for mesh ',i3,' on level ', i3, /, &
            '======================================================================================')

END SUBROUTINE SCARC_DEBUG_QUANTITY


!!> ------------------------------------------------------------------------------------------------
!!> Print out vector information on level NL for matlab
!!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_MATLAB_VECTOR (NV, CVEC, NL)
INTEGER, INTENT(IN) :: NV, NL
INTEGER :: NM
CHARACTER (*), INTENT(IN) :: CVEC
INTEGER :: JC, MVEC
CHARACTER(60) :: CNAME, CFORM
REAL (EB), POINTER, DIMENSION(:) :: VC
TYPE (SCARC_LEVEL_TYPE), POINTER :: L=>NULL()

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   L  => SCARC(NM)%LEVEL(NL)
   VC => POINT_TO_VECTOR (NM, NL, NV)
   WRITE (CNAME, '(A,A1,A,i2.2,A,i2.2,A)') 'matlab/',CVEC,'_mesh',NM,'_level',NL,'_vec.txt'
   WRITE (CFORM, '(I3, A)' ) L%NC-1, "(F24.16,;),F24.16" 
   MVEC=GET_FILE_NUMBER()
   OPEN(MVEC,FILE=CNAME)
   WRITE(MVEC, *) CVEC, ' = ['
   WRITE(MVEC,CFORM) (VC(JC),JC=1,L%NC)
   WRITE(MVEC, *) ' ]'
   CLOSE(MVEC)
ENDDO

END SUBROUTINE SCARC_MATLAB_VECTOR


!!> ------------------------------------------------------------------------------------------------
!!> Print out matrix information on level NL for matlab
!!> ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_MATLAB_MATRIX(VAL, ROW, COL, NC1, NC2, NM, NL, CMATRIX)
REAL(EB), DIMENSION(:), INTENT(IN) :: VAL
INTEGER, DIMENSION(:), INTENT(IN) :: ROW
INTEGER, DIMENSION(:), INTENT(IN) :: COL
INTEGER, INTENT(IN) :: NM, NL, NC1, NC2
CHARACTER(*), INTENT(IN) :: CMATRIX
INTEGER :: IC, JC, ICOL, MMATRIX
CHARACTER(60) :: CNAME, CFORM
REAL(EB) :: MATRIX_LINE(1000)

WRITE (CNAME, '(A,A1,A,i2.2,A,i2.2,A)') 'matlab/',CMATRIX,'_mesh',NM,'_level',NL,'_mat.txt'
!WRITE (CFORM, '(A,I3, 2A)' ) "(", NC2-1, "(F7.1,','),F7.1,';')" 
WRITE (CFORM, '(A,I3, 2A)' ) "(", NC2-1, "(F7.1,' '),F7.1,' ')" 
MMATRIX=GET_FILE_NUMBER()
OPEN(MMATRIX,FILE=CNAME)
!WRITE(MMATRIX, *) CMATRIX, ' = ['
DO IC = 1, NC1
   MATRIX_LINE=0.0_EB
   DO JC = 1, NC2
      DO  ICOL= ROW(IC), ROW(IC+1)-1
         IF (COL(ICOL)==JC) MATRIX_LINE(JC)=VAL(ICOL)
      ENDDO
   ENDDO
   WRITE(MMATRIX, CFORM) (MATRIX_LINE(JC),JC=1,NC2)
ENDDO
!WRITE(MMATRIX, *) ' ]'
CLOSE(MMATRIX)
 
END SUBROUTINE SCARC_MATLAB_MATRIX

#endif
!!> ================================================================================================
!!> END DEBUGGING PART  -  only enabled if directive is set
!!> ================================================================================================


END MODULE SCRC

!WRITE(*,*) 'TEST:', __LINE__,__FILE__
