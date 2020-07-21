MODULE SCRC
            
USE PRECISION_PARAMETERS
USE GLOBAL_CONSTANTS
USE MESH_POINTERS
USE MEMORY_FUNCTIONS, ONLY: CHKMEMERR
USE COMP_FUNCTIONS, ONLY: SECOND, GET_FILE_NUMBER
USE MPI
 
IMPLICIT NONE

CHARACTER(255), PARAMETER :: scrchid='$Id: scrc.f90 9718 2011-12-30 17:49:06Z drjfloyd $'
CHARACTER(255), PARAMETER :: scrcrev='$Revision: 9718 $'
CHARACTER(255), PARAMETER :: scrcdate='$Date: 2011-12-30 18:49:06 +0100 (Fr, 30 Dez 2011) $'

PRIVATE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Public structures (needed in main, read, divg, dump)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!----------------------------------------------------------------------------------------------------
!!! Public subroutines (initialization, solver, time measurement and revisioning)
!!!----------------------------------------------------------------------------------------------------
PUBLIC SCARC_SETUP                   
PUBLIC SCARC_SOLVER                 
PUBLIC SCARC_TIMINGS               
PUBLIC GET_REV_SCRC               

!!!----------------------------------------------------------------------------------------------------
!!! Public variables   (explanations in declaration part below)
!!! Note: For input parameters in character format corresponding INTEGER type-parameters will
!!! be introduced later to simplify inquiries
!!!----------------------------------------------------------------------------------------------------
PUBLIC SCARC_METHOD   
PUBLIC SCARC_SYSTEM  
PUBLIC SCARC_INITIAL  

PUBLIC SCARC_RESIDUAL  
PUBLIC SCARC_ITERATIONS       
PUBLIC SCARC_CAPPA     
PUBLIC SCARC_ACCURACY        
PUBLIC SCARC_ACCURACY_DIVERGENCE 
PUBLIC SCARC_ACCURACY_RELATIVE 

PUBLIC SCARC_MULTIGRID      
PUBLIC SCARC_MULTIGRID_LEVEL
PUBLIC SCARC_MULTIGRID_CYCLE
PUBLIC SCARC_MULTIGRID_COARSENING 
PUBLIC SCARC_MULTIGRID_INTERPOL 
PUBLIC SCARC_MULTIGRID_ITERATIONS       
PUBLIC SCARC_MULTIGRID_ACCURACY      

PUBLIC SCARC_KRYLOV      
PUBLIC SCARC_KRYLOV_ITERATIONS
PUBLIC SCARC_KRYLOV_ACCURACY

PUBLIC SCARC_SMOOTH     
PUBLIC SCARC_SMOOTH_ITERATIONS 
PUBLIC SCARC_SMOOTH_ACCURACY 
PUBLIC SCARC_SMOOTH_OMEGA 

PUBLIC SCARC_PRECON      
PUBLIC SCARC_PRECON_ITERATIONS 
PUBLIC SCARC_PRECON_ACCURACY
PUBLIC SCARC_PRECON_OMEGA 

PUBLIC SCARC_COARSE      
PUBLIC SCARC_COARSE_ITERATIONS  
PUBLIC SCARC_COARSE_ACCURACY  
PUBLIC SCARC_COARSE_OMEGA
PUBLIC SCARC_COARSE_PRECON 
 
PUBLIC SCARC_DEBUG 
PUBLIC SCARC_MKL        !!! MKL code

!!!----------------------------------------------------------------------------------------------------
!!! corresponding declarations (with default settings)
!!!----------------------------------------------------------------------------------------------------

!!! General definitions
CHARACTER(20) :: SCARC_METHOD    = 'null'                     ! requested solver method (KRYLOV/MULTIGRID)
CHARACTER(20) :: SCARC_SYSTEM    = 'null'                     ! matrix storage technique (BANDED/COMPACT)
CHARACTER(20) :: SCARC_INITIAL   = 'null'                     ! initial solution (currently only default is used)

!!! General iteration parameters
REAL (EB)     :: SCARC_RESIDUAL             = -1.0_EB         ! residual of global selected solver
INTEGER       :: SCARC_ITERATIONS           =  0              ! number of iterations of selected ScaRC solver
REAL (EB)     :: SCARC_CAPPA                =  1.0_EB         ! convergence rate of selected ScarC solver
REAL (EB)     :: SCARC_ACCURACY_DIVERGENCE  = 1.E+6_EB        ! divergence epsilon for all solvers
REAL (EB)     :: SCARC_ACCURACY_RELATIVE    = 1.E-2_EB        ! minimum relative accuracy for all solvers
CHARACTER(20) :: SCARC_ACCURACY             = 'ABSOLUTE'      ! accuracy type (ABSOLUTE/RELATIVE)

!!! Parameters for multigrid-type methods
CHARACTER(20) :: SCARC_MULTIGRID               = 'GEOMETRIC'  ! type of MG-method (GEOMETRIC/ALGEBRAIC)
INTEGER       :: SCARC_MULTIGRID_LEVEL         = -1           ! User defined number of MG-levels (optionally)
CHARACTER(1)  :: SCARC_MULTIGRID_CYCLE         = 'V'          ! Cycling type  (F/V/W)
CHARACTER(20) :: SCARC_MULTIGRID_COARSENING    = 'RS3'        ! Coarsening strategy  (RS3/A1/A2/PMIS/FDS...)
CHARACTER(20) :: SCARC_MULTIGRID_INTERPOL      = 'DIRECT'     ! Interpolation strategy (DIRECT/RS/STANDARD)
INTEGER       :: SCARC_MULTIGRID_ITERATIONS    = 1000         ! max number of iterations
REAL (EB)     :: SCARC_MULTIGRID_ACCURACY      = 1.E-15_EB    ! requested accuracy for convergence

!!! Parameters for Krylov-type methods
CHARACTER(20) :: SCARC_KRYLOV            = 'CG'               ! type of Krylov-method (CG/BICG)
INTEGER       :: SCARC_KRYLOV_ITERATIONS = 1000               ! max number of iterations 
REAL (EB)     :: SCARC_KRYLOV_ACCURACY   = 1.E-15_EB          ! requested accuracy for convergence

!!! Parameters for smoothing method (used in multigrids-methods)
CHARACTER(20) :: SCARC_SMOOTH            = 'SSOR'             ! smoother for MG (JACOBI/SSOR/GSTRIX)
INTEGER       :: SCARC_SMOOTH_ITERATIONS = 1000               ! max number of iterations 
REAL (EB)     :: SCARC_SMOOTH_ACCURACY   = 1.E-15_EB          ! requested accuracy for convergence
REAL (EB)     :: SCARC_SMOOTH_OMEGA      = 0.90E+0_EB         ! relaxation parameter 

!!! Parameters for preconditioning method (used in Krylov-methods)
CHARACTER(20) :: SCARC_PRECON            = 'SSOR'             ! preconditioner for CG/BICG (JACOBI/SSOR/GSTRIX/MG)
INTEGER       :: SCARC_PRECON_ITERATIONS = 1000               ! max number of iterations 
REAL (EB)     :: SCARC_PRECON_ACCURACY   = 1.E-15_EB          ! requested accuracy for convergence
REAL (EB)     :: SCARC_PRECON_OMEGA      = 0.9E+0_EB         ! relaxation parameter 

!!! Parameters for coarse grid method
CHARACTER(20) :: SCARC_COARSE            = 'ITERATIVE'        ! coarse grid solver (iterative/direct)
INTEGER       :: SCARC_COARSE_ITERATIONS = 100                ! max number of iterations for iterative variant
REAL (EB)     :: SCARC_COARSE_ACCURACY   = 1.E-12_EB          ! requested accuracy for convergencefor iterative variant
REAL (EB)     :: SCARC_COARSE_OMEGA      = 1.5E+0_EB          ! relaxation parameter for iterative variant
CHARACTER(20) :: SCARC_COARSE_PRECON     = 'SSOR'             ! preconditioner for iterative variant
 
!!! debugging parameters
CHARACTER(20) :: SCARC_DEBUG = 'NONE'                         ! debugging level (NONE/LESS/MEDIUM/MUCH)
CHARACTER(40) :: SCARC_FN                                     ! file name for ScaRC debug messages
INTEGER       :: SCARC_LU                                     ! unit number for ScaRC debug file

!!! MKL library
LOGICAL       :: SCARC_MKL = .FALSE.


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Private structures
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!----------------------------------------------------------------------------------------------------
!!! Global constants
!!!----------------------------------------------------------------------------------------------------
INTEGER, PARAMETER :: NSCARC_DIMENSION_NONE         = -1, &
                      NSCARC_DIMENSION_TWO          =  2, &      ! two-dimensional problem
                      NSCARC_DIMENSION_THREE        =  3         ! three-dimensional problem

INTEGER, PARAMETER :: NSCARC_SCOPE_NONE             = -1, &         
                      NSCARC_SCOPE_MAIN             =  1, &      ! method used as main solver
                      NSCARC_SCOPE_SMOOTH           =  2, &      ! method used as smoother
                      NSCARC_SCOPE_PRECON           =  3, &      ! method used as preconditiner
                      NSCARC_SCOPE_COARSE           =  4         ! method used as coarse grid solver

INTEGER, PARAMETER :: NSCARC_METHOD_NONE            = -1, &
                      NSCARC_METHOD_KRYLOV          =  1, &      ! Krylov   -method as global solver
                      NSCARC_METHOD_MULTIGRID       =  2         ! multigrid-method as global solver

INTEGER, PARAMETER :: NSCARC_KRYLOV_NONE            = -1, &
                      NSCARC_KRYLOV_CG              =  1, &      ! CG   as Krylov solver
                      NSCARC_KRYLOV_BICG            =  2         ! BICG as Krylov solver

INTEGER, PARAMETER :: NSCARC_MULTIGRID_NONE         = -1, &
                      NSCARC_MULTIGRID_GEOMETRIC    =  1, &      ! geometric multigrid
                      NSCARC_MULTIGRID_ALGEBRAIC    =  2         ! algebraic multigrid

INTEGER, PARAMETER :: NSCARC_SYSTEM_NONE            = -1, &
                      NSCARC_SYSTEM_BANDED          =  1, &      ! system in banded storage technique
                      NSCARC_SYSTEM_COMPACT         =  2         ! system in compact  storage technique

INTEGER, PARAMETER :: NSCARC_EXCHANGE_NONE          = -1, &
                      NSCARC_EXCHANGE_IJKW          =  1, &      ! initialize communication
                      NSCARC_EXCHANGE_VECTOR        =  2, &      ! matrix-vector communication 
                      NSCARC_EXCHANGE_BDRY          =  3, &      ! vector values along internal boundaries
                      NSCARC_EXCHANGE_MATRIX        =  4, &      ! internal subdiagonal matrix values 
                      NSCARC_EXCHANGE_MEASURE       =  5, &      ! measure values along internal boundaries
                      NSCARC_EXCHANGE_CELLTYPE      =  6, &      ! cell types along internal boundaries
                      NSCARC_EXCHANGE_CELLTYPE2     =  7, &      ! cell types along internal boundaries, coarse
                      NSCARC_EXCHANGE_WEIGHTS       =  8, &      ! internal transfer weights
                      NSCARC_EXCHANGE_AMG           =  9, &      ! AMG-information on neighboring meshes
                      NSCARC_EXCHANGE_GRID          = 10         ! neighboring grid dimensions

INTEGER, PARAMETER :: NSCARC_SMOOTH_NONE            = -1, &
                      NSCARC_SMOOTH_JACOBI          =  1, &      ! smoothing by JACOBI-method
                      NSCARC_SMOOTH_SSOR            =  2, &      ! smoothing by SSOR-method
                      NSCARC_SMOOTH_GSTRIX          =  3, &      ! smoothing by GSTRIX-method
                      NSCARC_SMOOTH_FFT             =  4         ! smoothing by FFT-method

INTEGER, PARAMETER :: NSCARC_PRECON_NONE            = -1, &
                      NSCARC_PRECON_JACOBI          =  1, &      ! preconditioning by JACOBI-method
                      NSCARC_PRECON_SSOR            =  2, &      ! preconditioning by SSOR-method
                      NSCARC_PRECON_GSTRIX          =  3, &      ! preconditioning by GSTRIX-method
                      NSCARC_PRECON_FFT             =  4, &      ! preconditioning by FFT-method
                      NSCARC_PRECON_MULTIGRID       =  5         ! preconditioning by MG-method

INTEGER, PARAMETER :: NSCARC_CYCLE_NONE             = -1, &
                      NSCARC_CYCLE_F                =  0, &      ! F-cycle for mg-method
                      NSCARC_CYCLE_V                =  1, &      ! V-cycle for mg-method
                      NSCARC_CYCLE_W                =  2, &      ! W-cycle for mg-method
                      NSCARC_CYCLE_SETUP            =  3, &      ! initialize cycle counts
                      NSCARC_CYCLE_RESET            =  4, &      ! reset cycle counts
                      NSCARC_CYCLE_PROCEED          =  5, &      ! proceed cycle counts
                      NSCARC_CYCLE_PRESMOOTH        =  6, &      ! presmoothing cycle
                      NSCARC_CYCLE_POSTSMOOTH       =  7, &      ! postsmoothing cycle
                      NSCARC_CYCLE_NEXT             =  8, &      ! perform next cycling loop
                      NSCARC_CYCLE_EXIT             =  9         ! exit cycling loop

INTEGER, PARAMETER :: NSCARC_STATE_PROCEED           =  0, &      ! proceed loop
                      NSCARC_STATE_CONV              =  1, &      ! convergence
                      NSCARC_STATE_DIVG              =  2         ! divergence

INTEGER, PARAMETER :: NSCARC_DEBUG_NONE             = -1, &      ! no debugging requested
                      NSCARC_DEBUG_LESS             =  1, &      ! low    level of debugging requested
                      NSCARC_DEBUG_MEDIUM           =  2, &      ! medium level of debugging requested
                      NSCARC_DEBUG_MUCH             =  3, &      ! strong level of debugging requested
                      NSCARC_DEBUG_MATRIX           =  4, &      ! show matrix
                      NSCARC_DEBUG_MATRIXE          =  5, &      ! show matrix
                      NSCARC_DEBUG_IJKW             =  6, &      ! show IJKW
                      NSCARC_DEBUG_BCINDEX          =  7, &      ! show PRESSURE_BC_INDEX
                      NSCARC_DEBUG_ACELL            =  8, &      ! show ADJACENT_CELL
                      NSCARC_DEBUG_GCELL            =  9, &      ! show GHOST_CELL
                      NSCARC_DEBUG_SUBDIVISION      = 10, &      ! show SUBDIVISION
                      NSCARC_DEBUG_MEASURE          = 11, &      ! show MEASURE
                      NSCARC_DEBUG_CELLTYPE         = 12, &      ! show CELLTYPE
                      NSCARC_DEBUG_COARSE           = 13, &      ! show coarse grid
                      NSCARC_DEBUG_PROLONGATION     = 14, &      ! show prolongation matrix
                      NSCARC_DEBUG_RESTRICTION      = 15         ! show restriction matrix 

INTEGER, PARAMETER :: NSCARC_COARSENING_NONE        = -1, &
                      NSCARC_COARSENING_RS3         =  1, &      ! parallel Ruge-StÃ¼ben 
                      NSCARC_COARSENING_A1          =  2, &      ! aggressive 1 (path=1, length=2)
                      NSCARC_COARSENING_A2          =  3, &      ! aggressive 2 (path=2, length=2)
                      NSCARC_COARSENING_PMIS        =  4, &      ! PMIS 
                      NSCARC_COARSENING_FDSRS3      =  5, &      ! FDSRS3 : ScaRC-FDS variant similar to RS3
                      NSCARC_COARSENING_FDSA1       =  6, &      ! FDSA1  : ScaRC-FDS variant similar to A1
                      NSCARC_COARSENING_FDSA2       =  7, &      ! FDSA2  : ScaRC-FDS variant similar to A2
                      NSCARC_COARSENING_FDSPMIS     =  8, &      ! FDSPMIS: ScaRC-FDS variant similar to PMIS
                      NSCARC_COARSENING_PMIS2       =  8, &      ! PMIS2 (modified PMIS)
                      NSCARC_COARSENING_BDRY        = 10, &      ! own boundary test
                      NSCARC_COARSENING_TEST        = 11         ! own boundary test

INTEGER, PARAMETER :: NSCARC_COARSE_NONE            = -1, &
                      NSCARC_COARSE_ITERATIVE       =  1, &      ! iterative solution of coarse grid problem
                      NSCARC_COARSE_DIRECT          =  2         ! direct solution of coarse grid problem

INTEGER, PARAMETER :: NSCARC_DIRECT_NONE            = -1, &
                      NSCARC_DIRECT_GE              =  1, &      ! direct solution by Gaussian elmination
                      NSCARC_DIRECT_LU              =  2         ! direct solution by LU-decomposition 

INTEGER, PARAMETER :: NSCARC_VECTOR_NONE            = -1, &
                      NSCARC_VECTOR_X               =  1, &      ! selection parameter for vector X
                      NSCARC_VECTOR_F               =  2, &      ! selection parameter for vector F
                      NSCARC_VECTOR_Y               =  3, &      ! selection parameter for vector Y
                      NSCARC_VECTOR_G               =  4, &      ! selection parameter for vector G
                      NSCARC_VECTOR_W               =  5, &      ! selection parameter for vector R
                      NSCARC_VECTOR_D               =  6, &      ! selection parameter for vector D
                      NSCARC_VECTOR_Z               =  7, &      ! selection parameter for vector Z
                      NSCARC_VECTOR_X2              =  8, &      ! selection parameter for vector X2
                      NSCARC_VECTOR_D2              =  9, &      ! selection parameter for vector D2
                      NSCARC_VECTOR_W2              = 10, &      ! selection parameter for vector R2
                      NSCARC_VECTOR_Y2              = 11, &      ! selection parameter for vector Y2
                      NSCARC_VECTOR_H               = 12, &      ! selection parameter for vector Y2
                      NSCARC_VECTOR_HS              = 13, &      ! selection parameter for vector Y2
                      NSCARC_VECTOR_MEASURE         = 14, &      ! selection parameter for vector MEASURE
                      NSCARC_VECTOR_CELLTYPE        = 15         ! selection parameter for vector CELLTYPE

INTEGER, PARAMETER :: NSCARC_MATRIX_NONE            = -1, &
                      NSCARC_MATRIX_SUBDIAG         =  1, &      ! exchange subdiagonal matrix entries
                      NSCARC_MATRIX_SYSTEM1         =  2, &      ! exchange system matrix entries on ghost cells
                      NSCARC_MATRIX_SYSTEM2         =  3, &      ! exchange system matrix entries on ghost cells
                      NSCARC_MATRIX_PROLONGATION    =  4, &      ! exchange prolongation matrix on ghost cells
                      NSCARC_MATRIX_RESTRICTION     =  5         ! exchange restriction matrix on ghost cells

INTEGER, PARAMETER :: NSCARC_ACCURACY_NONE          = -1, &
                      NSCARC_ACCURACY_ABSOLUTE      =  1, &      ! absolute accuracy must be reached
                      NSCARC_ACCURACY_RELATIVE      =  2         ! relative accuracy must be reached

REAL(EB), PARAMETER:: NSCARC_MEASURE_NONE           =  0.0_EB, &
                      NSCARC_MEASURE_COARSE         =  6.0_EB, & ! coarse-grid cell
                      NSCARC_MEASURE_FINE           =  5.0_EB, & ! fine-grid cell
                      NSCARC_MEASURE_SFINE          =  5.0_EB, & ! strongly coupled fine-grid cell
                      NSCARC_MEASURE_WFINE          =  4.0_EB    ! weakly   coupled fine-grid cell

INTEGER, PARAMETER :: NSCARC_CELLTYPE_NONE          =  0, &
                      NSCARC_CELLTYPE_COARSE        =  1, &      ! coarse-grid cell
                      NSCARC_CELLTYPE_FINE          = -2, &      ! fine-grid cell
                      NSCARC_CELLTYPE_SFINE         = -2, &      ! strongly coupled fine-grid cell
                      NSCARC_CELLTYPE_WFINE         = -3         ! weakly   coupled fine-grid cell

INTEGER, PARAMETER :: NSCARC_INTERPOL_NONE          =  0, &
                      NSCARC_INTERPOL_STANDARD      =  1, &      ! standard interpolation
                      NSCARC_INTERPOL_CLASSICAL     =  2, &      ! classical interpolation
                      NSCARC_INTERPOL_DIRECT        =  3, &      ! direct interpolation
                      NSCARC_INTERPOL_MULTIPASS     =  4         ! multipass interpolation

INTEGER, PARAMETER :: NSCARC_LATEX_NONE             = -1, &      ! no latex information requested
                      NSCARC_LATEX_STAGGERED        =  1, &      ! show staggered latex information
                      NSCARC_LATEX_EQUAL            =  2, &      ! show equal latex information
                      NSCARC_LATEX_NUMBER           =  3         ! show number latex information

INTEGER, PARAMETER :: NSCARC_TIME_NONE              = -1, &
                      NSCARC_TIME_TOTAL             =  1, &      ! time for complete ScaRC part of FDS
                      NSCARC_TIME_SETUP             =  2, &      ! time for setup phase
                      NSCARC_TIME_SOLVER            =  3, &      ! time for ScaRC solver 
                      NSCARC_TIME_KRYLOV            =  4, &      ! time for Krylov solver
                      NSCARC_TIME_MULTIGRID         =  5, &      ! time for multigrid solver
                      NSCARC_TIME_PRECON            =  6, &      ! time for preconditioner
                      NSCARC_TIME_SMOOTH            =  7, &      ! time for smoother
                      NSCARC_TIME_COARSE            =  8, &      ! time for coarse grid solver
                      NSCARC_TIME_MATVEC            =  9, &      ! time for matrix-vector product
                      NSCARC_TIME_SCALPROD          = 10, &      ! time for scalar product
                      NSCARC_TIME_L2NORM            = 11, &      ! time for l2norm
                      NSCARC_TIME_EXCH_INIT         = 12, &      ! time for exchange initialization
                      NSCARC_TIME_EXCH_VECTOR       = 13, &      ! time for exchange of internal boundary
                      NSCARC_TIME_EXCH_MATRIX       = 14, &      ! time for exchange of internal matrix
                      NSCARC_TIME_EXCH_MEASURE      = 15         ! time for exchange of internal measure

INTEGER, PARAMETER :: NSCARC_LEVEL_NONE             = -1, &      ! no predefined initial solution used
                      NSCARC_LEVEL_MIN              =  0, &      ! minimum multigrid level 
                      NSCARC_LEVEL_MAX              =  15        ! maximum multigrid level 

INTEGER, PARAMETER :: NSCARC_INITIAL_NONE           = -1         ! another initial function ?
 
INTEGER, PARAMETER :: NSCARC_DUMMY                  = -1         ! dummy variable (needed at several places)


!!!----------------------------------------------------------------------------------------------------
!!! Global variables 
!!!----------------------------------------------------------------------------------------------------
!!! use integer types for the user defined input data (based on SCARC_TYPE_... variables)
INTEGER :: TYPE_DIMENSION  = NSCARC_DIMENSION_NONE
INTEGER :: TYPE_SCOPE      = NSCARC_SCOPE_NONE
INTEGER :: TYPE_METHOD     = NSCARC_METHOD_NONE
INTEGER :: TYPE_KRYLOV     = NSCARC_KRYLOV_NONE
INTEGER :: TYPE_MULTIGRID  = NSCARC_MULTIGRID_NONE
INTEGER :: TYPE_SYSTEM     = NSCARC_SYSTEM_NONE
INTEGER :: TYPE_MATRIX     = NSCARC_MATRIX_NONE
INTEGER :: TYPE_ACCURACY   = NSCARC_ACCURACY_NONE
INTEGER :: TYPE_SMOOTH     = NSCARC_SMOOTH_NONE
INTEGER :: TYPE_PRECON     = NSCARC_PRECON_NONE
INTEGER :: TYPE_CYCLE      = NSCARC_CYCLE_NONE
INTEGER :: TYPE_COARSENING = NSCARC_COARSENING_NONE
INTEGER :: TYPE_INTERPOL   = NSCARC_INTERPOL_NONE
INTEGER :: TYPE_COARSE     = NSCARC_COARSE_NONE
INTEGER :: TYPE_DEBUG      = NSCARC_DEBUG_NONE   
INTEGER :: TYPE_INITIAL    = NSCARC_INITIAL_NONE
INTEGER :: TYPE_EXCHANGE   = NSCARC_EXCHANGE_NONE
INTEGER :: TYPE_VECTOR     = NSCARC_VECTOR_NONE
INTEGER :: TYPE_DIRECT     = NSCARC_DIRECT_NONE

!!! range of meshes which must be processed for MYID
INTEGER :: NMESHES_MIN, NMESHES_MAX                 
 
!!! total, minimum and maximum number of multigrid levels
INTEGER :: NLEVEL, NLEVEL_MAX, NLEVEL_MIN                              

INTEGER :: NMASTER, NC_COARSE0

!!! additional arrays for data exchange
INTEGER :: NREQ_SCARC, N_EXCHANGES, TAG_SCARC, SNODE, RNODE, STATUS2_SCARC(MPI_STATUS_SIZE)
INTEGER, ALLOCATABLE, DIMENSION (:)    :: REQ_SCARC

!!! time measurements with ScaRC
INTEGER, PARAMETER :: N_TIMERS_SCARC=18         
REAL(EB), ALLOCATABLE, DIMENSION(:,:) :: TUSED_SCARC

!!! auxiliary variables for global and local scalproducts, number of global cells and global dimensions
REAL(EB), ALLOCATABLE, DIMENSION (:) :: SP_LOCAL
INTEGER,  ALLOCATABLE, DIMENSION (:) :: NC_GLOBAL, NC_LOCAL
INTEGER,  ALLOCATABLE, DIMENSION (:) :: NA_COARSE, NC_COARSE, NX_COARSE, NY_COARSE, NZ_COARSE, NW_COARSE
REAL(EB), ALLOCATABLE, DIMENSION (:) :: DXI2_COARSE, DYI2_COARSE, DZI2_COARSE

!!! number of couplings in given matrix stencil and pointer to indices in matrix stencil on given level 
INTEGER :: ID, ILX, ILY, ILZ, IUX, IUY, IUZ

!!! Index numbers of different vector types (used in different ScaRC-solvers)
INTEGER :: VEC_NONE, VEC_X, VEC_F, VEC_Y, VEC_G, VEC_W, VEC_D, VEC_Z, VEC_X2, VEC_D2, VEC_W2, VEC_Y2

!!! counter and displacement arrays for global data exchanges
INTEGER, ALLOCATABLE, DIMENSION(:) :: COUNTS_SCARC, DISPLS_SCARC

 
!!! Private type declarations

!!!----------------------------------------------------------------------------------------------------
!!! OSCARC type for banded system information on other mesh
!!!----------------------------------------------------------------------------------------------------
TYPE OSCARC_BANDED_TYPE
INTEGER :: NX, NY, NZ, NC, NW
INTEGER, POINTER, DIMENSION (:, :) :: IJKW
END TYPE OSCARC_BANDED_TYPE


!!!----------------------------------------------------------------------------------------------------
!!! OSCARC type for mesh information on other mesh
!!! local numbers of cells (also per direction), wall cells, matrix entries
!!! neighborship structures 
!!!----------------------------------------------------------------------------------------------------
TYPE OSCARC_COMPACT_TYPE
INTEGER :: NX, NY, NZ, NC, NCE, NW
INTEGER, POINTER, DIMENSION (:, :) :: IJKW
INTEGER, POINTER, DIMENSION (:)    :: ADJACENT_CELL, GHOST_CELL, GHOST_PTR
END TYPE OSCARC_COMPACT_TYPE


!!!----------------------------------------------------------------------------------------------------
!!! SCARC type for matrix stencil information
!!! offset in stencil, measure of cells, coarse and fine cells, couplings
!!!----------------------------------------------------------------------------------------------------
TYPE SCARC_PRECON_TYPE
REAL (EB), POINTER, DIMENSION (:)       :: MDX, LDX, UDX, MDY, MDZ, DWORK, PERIOD
REAL (EB), POINTER, DIMENSION (:, :, :) :: FFT
END TYPE SCARC_PRECON_TYPE


!!!----------------------------------------------------------------------------------------------------
!!! SCARC types for different matrix storage technique 
!!!   - NX, NY, NZ : correspond to IBAR, JBAR, KBAR
!!!   - NC, NCE    : number of cells and cells_extended (including ghost cells)
!!!   - NA, NP, NR : number of matrix entries in A, P and R
!!!   - NCPL       : number of couplings in matrix stencil on finest level (2D: 5, 3D: 7)
!!!   - NCC, NCF   : number of coarse and fine cells (only used in AMG)
!!!----------------------------------------------------------------------------------------------------
!!! Banded storage technique 
TYPE SCARC_BANDED_TYPE

INTEGER :: SUBDIVISION(3,-3:3)
INTEGER  , POINTER, DIMENSION (:, :)    :: IJKW
INTEGER  , POINTER, DIMENSION (:)       :: BC_INDEX, ADJACENT_CELL, GHOST_CELL
INTEGER  , POINTER, DIMENSION (:, :, :) :: CELLTYPE
REAL (EB), POINTER, DIMENSION (:, :)    :: A
REAL (EB), POINTER, DIMENSION (:, :, :) :: X , F , D , Y , G , W, Z
REAL (EB), POINTER, DIMENSION (:, :, :) :: X2, F2, D2, Y2, G2, W2
INTEGER   :: NX, NY, NZ, NW
INTEGER   :: NC, NCE, NCG, NCEP
INTEGER   :: NA, NCPL
REAL(EB)  :: DXI, DYI, DZI, DXI2, DYI2, DZI2, DI2(3)

END TYPE SCARC_BANDED_TYPE


!!! Compact storage technique 
TYPE SCARC_COMPACT_TYPE

INTEGER :: SUBDIVISION(3,-3:3)
INTEGER,   POINTER, DIMENSION (:, :) :: IJKW
INTEGER,   POINTER, DIMENSION (:, :) :: WC_INDICES, P_PTR, W_PTR
INTEGER,   POINTER, DIMENSION (:)    :: CELLTYPE, BC_INDEX, ADJACENT_CELL, GHOST_CELL, GHOST_PTR
REAL(EB),  POINTER, DIMENSION (:)    :: MEASURE
INTEGER,   POINTER, DIMENSION (:)    :: A_ROW, A_COL      ! row and column pointers for system matrix A
INTEGER,   POINTER, DIMENSION (:)    :: P_ROW, P_COL      ! row and column pointers for prolongation matrix P
INTEGER,   POINTER, DIMENSION (:)    :: R_ROW, R_COL      ! row and column pointers for restriction matrix A
REAL (EB), POINTER, DIMENSION (:)    :: A , P , R
REAL (EB), POINTER, DIMENSION (:)    :: X , F , D , Y , G , W, Z
REAL (EB), POINTER, DIMENSION (:)    :: X2, F2, D2, Y2, G2, W2
INTEGER   :: NX, NY, NZ, NW
INTEGER   :: NC, NCE, NCG, NCEP
INTEGER   :: NA, NCPL
INTEGER   :: NCC, NCCE, NCCP, NCF, NCW
INTEGER   :: NP, NR
REAL(EB)  :: DXI, DYI, DZI, DXI2, DYI2, DZI2, DI2(3)

END TYPE SCARC_COMPACT_TYPE


!!!----------------------------------------------------------------------------------------------------
!!! General ScaRC type with pointers to different structures on all grid levels
!!!----------------------------------------------------------------------------------------------------
TYPE SCARC_TYPE

INTEGER :: CYCLE_COUNT(2, NSCARC_LEVEL_MAX) = 0
REAL (EB), POINTER, DIMENSION (:,:) :: A_COARSE
REAL (EB), POINTER, DIMENSION (:)   :: X_COARSE, X_BUF
INTEGER  , POINTER, DIMENSION (:,:) :: PIVOT
INTEGER  , POINTER, DIMENSION (:)   :: OFFSET, COUNTS1, COUNTS2, DISPLS1, DISPLS2

TYPE (SCARC_PRECON_TYPE), POINTER, DIMENSION(:) :: PRECON

TYPE (SCARC_BANDED_TYPE) , POINTER, DIMENSION(:) :: BANDED
TYPE (SCARC_COMPACT_TYPE), POINTER, DIMENSION(:) :: COMPACT

TYPE (OSCARC_TYPE), POINTER, DIMENSION(:) :: OSCARC

END TYPE SCARC_TYPE


!!!----------------------------------------------------------------------------------------------------
!!! General OSCARC type on other mesh with mesh and exchange structures
!!!----------------------------------------------------------------------------------------------------
TYPE OSCARC_TYPE

REAL (EB), POINTER, DIMENSION (:) :: SEND_BUF, RECV_BUF
INTEGER :: NICMAX_R=0, NICMAX_S=0
INTEGER :: I_MIN_R=-10,I_MAX_R=-10,J_MIN_R=-10,J_MAX_R=-10,K_MIN_R=-10,K_MAX_R=-10,NIC_R=0, &
           I_MIN_S=-10,I_MAX_S=-10,J_MIN_S=-10,J_MAX_S=-10,K_MIN_S=-10,K_MAX_S=-10,NIC_S=0

TYPE (OSCARC_BANDED_TYPE) , POINTER, DIMENSION(:) :: BANDED
TYPE (OSCARC_COMPACT_TYPE), POINTER, DIMENSION(:) :: COMPACT

END TYPE OSCARC_TYPE


!!! globally used types
TYPE ( SCARC_TYPE), SAVE, DIMENSION(:), ALLOCATABLE, TARGET ::  SCARC
TYPE (OSCARC_TYPE), SAVE, DIMENSION(:), ALLOCATABLE, TARGET :: OSCARC


CONTAINS
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! SCARC_SETUP : Initialize ScaRC structures 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SETUP
INTEGER :: IERR
REAL(EB):: TNOW_SETUP

!!!----------------------------------------------------------------------------------------------------
!!! Initialize time measurement 
!!!----------------------------------------------------------------------------------------------------
IERR = 0
TNOW_SETUP = SECOND()

ALLOCATE(TUSED_SCARC(0:N_TIMERS_SCARC,NMESHES),STAT=IERR)
CALL ChkMemErr('SCARC_SETUP','TUSED_SCARC',IERR)

TUSED_SCARC = 0._EB
TUSED_SCARC(NSCARC_TIME_TOTAL,:) = SECOND()


!!!----------------------------------------------------------------------------------------------------
!!! Parse input parameters
!!!----------------------------------------------------------------------------------------------------
CALL SCARC_INPUT_PARSER               ! parse input parameters (corresponding to FDS-file)

WRITE(*,*) 'SCARC_MKL=',SCARC_MKL


!!!----------------------------------------------------------------------------------------------------
!!! Setup different components of ScaRC 
!!!----------------------------------------------------------------------------------------------------
CALL SCARC_SETUP_DIMENSION            ! define dimension of underlying problem
CALL SCARC_SETUP_DEBUGGING            ! open debug file if requested
CALL SCARC_SETUP_PROCESSES            ! determine set of meshes which must be processed on MYID
CALL SCARC_SETUP_LEVELS               ! define number of necessary grid levels 
CALL SCARC_SETUP_TYPES                ! allocate requested ScaRC-types for all necessary grid levels
CALL SCARC_SETUP_MESHES               ! set mesh information
CALL SCARC_SETUP_WALLCELLS            ! set wall cell information
CALL SCARC_SETUP_EXCHANGE             ! set information for data exchange
CALL SCARC_SETUP_SYSTEM               ! assemble system matrix with boundary conditions and solver vectors
CALL SCARC_SETUP_COARSENING           ! perform coarsening on different grid levels if requested (AMG only)
CALL SCARC_SETUP_VECTORS              ! allocate solution and auxiliary vectors on all needed grid levels
CALL SCARC_SETUP_GLOBALS              ! define some global variables
CALL SCARC_SETUP_NEIGHBORS            ! compute information about abutting neighbors on coarser levels


!!!----------------------------------------------------------------------------------------------------
!!! Measure time for setup routine
!!!----------------------------------------------------------------------------------------------------
TUSED_SCARC(NSCARC_TIME_SETUP,:)=TUSED_SCARC(NSCARC_TIME_SETUP,:)+SECOND()-TNOW_SETUP
TUSED_SCARC(NSCARC_TIME_TOTAL,:)=TUSED_SCARC(NSCARC_TIME_TOTAL,:)+SECOND()-TNOW_SETUP

END SUBROUTINE SCARC_SETUP


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Define dimension of problem
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SETUP_DIMENSION
IF (TWO_D) THEN
   TYPE_DIMENSION = NSCARC_DIMENSION_TWO
ELSE
   TYPE_DIMENSION = NSCARC_DIMENSION_THREE
ENDIF
END SUBROUTINE SCARC_SETUP_DIMENSION


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Determine types of input parameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_INPUT_PARSER
CHARACTER(300):: CMESSAGE

!!!----------------------------------------------------------------------------------------------------
!!! set type of global solver
!!!----------------------------------------------------------------------------------------------------
SELECT CASE (TRIM(SCARC_METHOD))

   CASE ('KRYLOV')

      TYPE_METHOD = NSCARC_METHOD_KRYLOV

      !!! set type of Krylov-method (CG/BICG)
      SELECT CASE (TRIM(SCARC_KRYLOV))
         CASE ('CG')
            TYPE_KRYLOV = NSCARC_KRYLOV_CG
         CASE ('BICG')
            TYPE_KRYLOV = NSCARC_KRYLOV_BICG
         CASE DEFAULT
            WRITE(CMESSAGE,1002) 'Krylov-solver',TRIM(SCARC_KRYLOV),'ScaRC','CG','BICG'
            CALL SCARC_SHUTDOWN(CMESSAGE)
      END SELECT 

      !!! set type of preconditioner (JACOBI/SSOR/GSTRIX/MG)
      SELECT CASE (TRIM(SCARC_PRECON))
         CASE ('JACOBI')
            TYPE_PRECON = NSCARC_PRECON_JACOBI
         CASE ('SSOR')
            TYPE_PRECON = NSCARC_PRECON_SSOR
         CASE ('GSTRIX')
            TYPE_PRECON = NSCARC_PRECON_GSTRIX
         CASE ('MULTIGRID')
            TYPE_PRECON = NSCARC_PRECON_MULTIGRID
         CASE ('FFT')
            TYPE_PRECON = NSCARC_PRECON_FFT
         CASE DEFAULT
            WRITE(CMESSAGE,1005) 'preconditioner',TRIM(SCARC_PRECON),&
                                 'Krylov-method','JACOBI','SSOR','GSTRIX','MG','FFT'
            CALL SCARC_SHUTDOWN(CMESSAGE)
      END SELECT

   CASE ('MULTIGRID')

      TYPE_METHOD = NSCARC_METHOD_MULTIGRID

      !!! set type of multigrid method (GEOMETRIC/ALGEBRAIC)
      SELECT CASE (TRIM(SCARC_MULTIGRID))
         CASE ('GEOMETRIC')
            TYPE_MULTIGRID = NSCARC_MULTIGRID_GEOMETRIC
         CASE ('ALGEBRAIC')
            TYPE_MULTIGRID = NSCARC_MULTIGRID_ALGEBRAIC
         CASE DEFAULT
            WRITE(CMESSAGE,1002) 'multigrid',TRIM(SCARC_MULTIGRID),&
                                 'ScaRC','GEOMETRIC','ALGEBRAIC'
            CALL SCARC_SHUTDOWN(CMESSAGE)
      END SELECT 

      !!! set type of multigrid method (GEOMETRIC/ALGEBRAIC)
      SELECT CASE (TRIM(SCARC_COARSE))
         CASE ('ITERATIVE')
            TYPE_COARSE = NSCARC_COARSE_ITERATIVE
            TYPE_KRYLOV = NSCARC_KRYLOV_CG
         CASE ('DIRECT')
            TYPE_COARSE = NSCARC_COARSE_DIRECT
            TYPE_DIRECT = NSCARC_DIRECT_GE
      END SELECT 

      !!! set type of smoother (JACOBI/SSOR/GSTRIX)
      SELECT CASE (TRIM(SCARC_SMOOTH))                          ! use same parameters as for preconditioner
         CASE ('JACOBI')
            TYPE_PRECON = NSCARC_PRECON_JACOBI
         CASE ('SSOR')
            TYPE_PRECON = NSCARC_PRECON_SSOR
         CASE ('GSTRIX')
            TYPE_PRECON = NSCARC_PRECON_GSTRIX
         CASE DEFAULT
            WRITE(CMESSAGE,1003) 'smoother',TRIM(SCARC_SMOOTH),'multigrid','JACOBI','SSOR','GSTRIX'
            CALL SCARC_SHUTDOWN(CMESSAGE)
      END SELECT

   CASE DEFAULT
      WRITE(CMESSAGE,1002) 'solver',TRIM(SCARC_METHOD),'ScaRC','KRYLOV','MULTIGRID'
      CALL SCARC_SHUTDOWN(CMESSAGE)

END SELECT 

!!!----------------------------------------------------------------------------------------------------
!!! if a multigrid solver is used (either as main solver or as preconditioner)
!!! set types for multigrid, coarse grid solver and cycling pattern
!!!----------------------------------------------------------------------------------------------------
IF (TYPE_METHOD == NSCARC_METHOD_MULTIGRID .OR. TYPE_PRECON == NSCARC_PRECON_MULTIGRID) THEN

   !!! set type of multigrid (GEOMETRIC/ALGEBRAIC with corresponding coarsening strategy)
   SELECT CASE (TRIM(SCARC_MULTIGRID))

      CASE ('GEOMETRIC')
         TYPE_MULTIGRID = NSCARC_MULTIGRID_GEOMETRIC

      CASE ('ALGEBRAIC')
         TYPE_MULTIGRID = NSCARC_MULTIGRID_ALGEBRAIC

         !!! set type of coarsening strategy (STANDARD/AGGRESSIVE)
         SELECT CASE (TRIM(SCARC_MULTIGRID_COARSENING))
            CASE ('RS3')
               TYPE_COARSENING = NSCARC_COARSENING_RS3
            CASE ('A1')
               TYPE_COARSENING = NSCARC_COARSENING_A1
            CASE ('A2')
               TYPE_COARSENING = NSCARC_COARSENING_A2
            CASE ('PMIS')
               TYPE_COARSENING = NSCARC_COARSENING_PMIS
            CASE ('FDSRS3')
               TYPE_COARSENING = NSCARC_COARSENING_FDSRS3
            CASE ('FDSA1')
               TYPE_COARSENING = NSCARC_COARSENING_FDSA1
            CASE ('FDSA2')
               TYPE_COARSENING = NSCARC_COARSENING_FDSA2
            CASE ('FDSPMIS')
               TYPE_COARSENING = NSCARC_COARSENING_FDSPMIS
            CASE ('PMIS2')
               TYPE_COARSENING = NSCARC_COARSENING_PMIS2
            CASE ('BDRY')
               TYPE_COARSENING = NSCARC_COARSENING_BDRY
            CASE ('TEST')
               TYPE_COARSENING = NSCARC_COARSENING_TEST
            CASE DEFAULT
               WRITE(CMESSAGE,1005) 'coarsening',TRIM(SCARC_MULTIGRID_COARSENING),&
                                    'algebraic multigrid','RS3'   ,'A1'   ,'A2'   ,'PMIS',&
                                                          'FDSRS3','FDSA1','FDSA1','FDSPMIS'
               CALL SCARC_SHUTDOWN(CMESSAGE)
         END SELECT

      CASE DEFAULT
         WRITE(CMESSAGE,1002) 'multigrid',TRIM(SCARC_MULTIGRID),&
                              'Krylov-preconditioner','GEOMETRIC','ALGEBRAIC'
         CALL SCARC_SHUTDOWN(CMESSAGE)
   END SELECT

   !!! set type of coarse grid solver (CG/GE)
   SELECT CASE (TRIM(SCARC_COARSE))
      CASE ('ITERATIVE')
         TYPE_COARSE = NSCARC_COARSE_ITERATIVE
         TYPE_KRYLOV = NSCARC_KRYLOV_CG
      CASE ('DIRECT')
         TYPE_COARSE = NSCARC_COARSE_DIRECT
         TYPE_DIRECT = NSCARC_DIRECT_GE
      CASE DEFAULT
         WRITE(CMESSAGE,1002) 'coarse grid solver',TRIM(SCARC_COARSE),'multigrid','ITERATIVE','DIRECT'
         CALL SCARC_SHUTDOWN(CMESSAGE)
   END SELECT

   !!! set type of cycling pattern (F/V/W)
   SELECT CASE (TRIM(SCARC_MULTIGRID_CYCLE))
      CASE ('F')
         TYPE_CYCLE = NSCARC_CYCLE_F
      CASE ('V')
         TYPE_CYCLE = NSCARC_CYCLE_V
      CASE ('W')
         TYPE_CYCLE = NSCARC_CYCLE_W
      CASE DEFAULT
         WRITE(CMESSAGE,1003) 'cycling ',TRIM(SCARC_MULTIGRID_CYCLE),'multigrid','F','V','W'
         CALL SCARC_SHUTDOWN(CMESSAGE)
   END SELECT

   !!! set type of interpolation (STANDARD/DIRECT/MULTIPASS)
   SELECT CASE (TRIM(SCARC_MULTIGRID_INTERPOL))
      CASE ('STANDARD')
         TYPE_INTERPOL = NSCARC_INTERPOL_STANDARD
      CASE ('CLASSICAL')
         TYPE_INTERPOL = NSCARC_INTERPOL_CLASSICAL
      CASE ('DIRECT')
         TYPE_INTERPOL = NSCARC_INTERPOL_DIRECT
      CASE ('MULTIPASS')
         TYPE_INTERPOL = NSCARC_INTERPOL_MULTIPASS
      CASE DEFAULT
         WRITE(CMESSAGE,1003) 'cycling ',TRIM(SCARC_MULTIGRID_INTERPOL), &
                              'multigrid','STANDARD','DIRECT','MULTIPASS'
         CALL SCARC_SHUTDOWN(CMESSAGE)
   END SELECT

ENDIF

!!!----------------------------------------------------------------------------------------------------
!!! set storage type (BANDED/COMPACT)
!!!----------------------------------------------------------------------------------------------------
SELECT CASE (TRIM(SCARC_SYSTEM))
   CASE ('BANDED')
      IF (TYPE_MULTIGRID == NSCARC_MULTIGRID_ALGEBRAIC) THEN
         WRITE(CMESSAGE,2002) 'system',TRIM(SCARC_SYSTEM),'ScaRC','TYPE_MULTIGRID','ALGEBRAIC',&
                              'TYPE_SYSTEM','COMPACT'
         TYPE_SYSTEM = NSCARC_SYSTEM_COMPACT
      ELSE
         TYPE_SYSTEM = NSCARC_SYSTEM_BANDED
      ENDIF
   CASE ('COMPACT')
      TYPE_SYSTEM = NSCARC_SYSTEM_COMPACT
   CASE ('null')

      SELECT_METHOD: SELECT CASE (TYPE_METHOD)

         CASE (NSCARC_METHOD_KRYLOV)

            SELECT_PRECON: SELECT CASE (TYPE_PRECON)
               CASE (NSCARC_PRECON_MULTIGRID) 
                  IF (TYPE_MULTIGRID == NSCARC_MULTIGRID_GEOMETRIC) THEN
                     TYPE_SYSTEM = NSCARC_SYSTEM_BANDED
                  ELSE
                     TYPE_SYSTEM = NSCARC_SYSTEM_COMPACT
                  ENDIF
               CASE DEFAULT
                  TYPE_SYSTEM = NSCARC_SYSTEM_BANDED
            END SELECT SELECT_PRECON

         CASE (NSCARC_METHOD_MULTIGRID)

            SELECT_MULTIGRID: SELECT CASE (TYPE_MULTIGRID)
               CASE (NSCARC_MULTIGRID_GEOMETRIC) 
                  TYPE_SYSTEM = NSCARC_SYSTEM_BANDED
               CASE (NSCARC_MULTIGRID_ALGEBRAIC) 
                  TYPE_SYSTEM = NSCARC_SYSTEM_COMPACT
            END SELECT SELECT_MULTIGRID

      END SELECT SELECT_METHOD

END SELECT


!!!----------------------------------------------------------------------------------------------------
!!! set type of accuracy (ABSOLUTE/RELATIVE)
!!!----------------------------------------------------------------------------------------------------
SELECT CASE (TRIM(SCARC_ACCURACY))
   CASE ('ABSOLUTE')
      TYPE_ACCURACY = NSCARC_ACCURACY_ABSOLUTE
   CASE ('RELATIVE')
      TYPE_ACCURACY = NSCARC_ACCURACY_RELATIVE
   CASE DEFAULT
      WRITE(CMESSAGE,1002) 'accuracy',TRIM(SCARC_ACCURACY),'ScaRC','ABSOLUTE','RELATIVE'
      CALL SCARC_SHUTDOWN(CMESSAGE)
END SELECT

!!!----------------------------------------------------------------------------------------------------
!!! set level of debugging (NONE/LESS/MEDIUM/MUCH)
!!!----------------------------------------------------------------------------------------------------
SELECT CASE (TRIM(SCARC_DEBUG))
   CASE ('NONE')
      TYPE_DEBUG = NSCARC_DEBUG_NONE
   CASE ('LESS')
      TYPE_DEBUG = NSCARC_DEBUG_LESS
   CASE ('MEDIUM')
      TYPE_DEBUG = NSCARC_DEBUG_MEDIUM
   CASE ('MUCH')
      TYPE_DEBUG = NSCARC_DEBUG_MUCH
   CASE DEFAULT
      WRITE(CMESSAGE,1004) 'debugging',TRIM(SCARC_DEBUG),'ScaRC','NONE','LESS','MEDIUM','MUCH'
      CALL SCARC_SHUTDOWN(CMESSAGE)
END SELECT

!!!----------------------------------------------------------------------------------------------------
!!! set type of initial solution (not yet used, may be used to define own initial vector)
!!!----------------------------------------------------------------------------------------------------
SELECT CASE (TRIM(SCARC_INITIAL))
   CASE ('null')
      TYPE_INITIAL = NSCARC_INITIAL_NONE
   CASE DEFAULT
      TYPE_INITIAL = NSCARC_INITIAL_NONE
END SELECT

1002 FORMAT ('ScaRC: Wrong ',A,' type ',A,' for ',A,/, &
             'ScaRC: Possible choices are ',A,'\/',A,/,&
             'ScaRC: Aborting program ...')
1003 FORMAT ('ScaRC: Wrong ',A,' type ',A,' for ',A,/, &
             'ScaRC: Possible choices are ',A,'\/',A,'\/',A,/,&
             'ScaRC: Aborting program ...')
1004 FORMAT ('ScaRC: Wrong ',A,' type ',A,' for ',A,/, &
             'ScaRC: Possible choices are ',A,'\/',A,'\/',A,'\/',A,/,&
             'ScaRC: Aborting program ...')
1005 FORMAT ('ScaRC: Wrong ',A,' type ',A,' for ',A,/, &
             'ScaRC: Possible choices are ',A,'\/',A,'\/',A,'\/',A,'\/',A,/,&
             'ScaRC: Aborting program ...')
!1006 FORMAT ('ScaRC: Wrong ',A,' type ',A,' for ',A,/, &
!             'ScaRC: Possible choices are ',A,'\/',A,'\/',A,'\/',A,'\/',A,'\/',A,/,&
!             'ScaRC: Aborting program ...')
!1007 FORMAT ('ScaRC: Wrong ',A,' type ',A,' for ',A,/, &
!             'ScaRC: Possible choices are ',A,'\/',A,'\/',A,'\/',A,'\/',A,'\/',A,'\/',A,/,&
!             'ScaRC: Aborting program ...')
!1008 FORMAT ('ScaRC: Wrong ',A,' type ',A,' for ',A,/, &
!             'ScaRC: Possible choices are ',A,'\/',A,'\/',A,'\/',A,'\/',A,'\/',A,'\/',A,'\/',A,/,&
!             'ScaRC: Aborting program ...')
2002 FORMAT ('ScaRC: Wrong ',A,' type ',A,' for ',A,/, ' in case of ',A,'=',A,/, &
             'ScaRC: Redefining ',A,' to ',A,' !',/)
END SUBROUTINE SCARC_INPUT_PARSER


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Setup debug file if requested
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SETUP_DEBUGGING
INTEGER:: NM

IF (TYPE_DEBUG /= NSCARC_DEBUG_NONE) THEN
   IF (USE_MPI) THEN
      DO NM=1,NMESHES
         IF (PROCESS(NM)/=MYID) CYCLE
         WRITE (SCARC_FN, '(A,A,i3.3)') TRIM(CHID),'.scarc',MYID+1
         SCARC_LU = GET_FILE_NUMBER()
         OPEN (SCARC_LU, FILE=SCARC_FN)
      ENDDO
   ELSE
      WRITE (SCARC_FN, '(A,A,i3.3)') TRIM(CHID),'.scarc',MYID+1
      SCARC_LU = GET_FILE_NUMBER()
      OPEN (SCARC_LU, FILE=SCARC_FN)
   ENDIF
ENDIF

END SUBROUTINE SCARC_SETUP_DEBUGGING


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Serial or parallel version ? Decide which meshes must be processed
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SETUP_PROCESSES
IF (USE_MPI) THEN
   NMESHES_MIN = MYID+1
   NMESHES_MAX = MYID+1
ELSE
   NMESHES_MIN = 1
   NMESHES_MAX = NMESHES
ENDIF
END SUBROUTINE SCARC_SETUP_PROCESSES

 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Determine number of grid levels  (1 for CG/BICG-method, NLEVEL for MG-method)
!!! Note: NLEVEL_MIN corresponds to finest grid resolution, NLEVEL_MAX to coarsest resolution
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SETUP_LEVELS
INTEGER :: KLEVEL(3), KLEVEL_MIN, NM

SELECT CASE (TYPE_MULTIGRID)

   !!!------------------------------------------------------------------------------------------------------
   !!! predefined hierarchy of levels in case of geometric multigrid-method
   !!!------------------------------------------------------------------------------------------------------
   CASE (NSCARC_MULTIGRID_GEOMETRIC)

      NLEVEL = NSCARC_LEVEL_MAX
      KLEVEL = NSCARC_LEVEL_MAX

      DO NM=1,NMESHES

         KLEVEL(1)=SCARC_GET_MAXLEVEL(MESHES(NM)%IBAR,1)
         IF (TYPE_DIMENSION == NSCARC_DIMENSION_THREE) KLEVEL(2)=SCARC_GET_MAXLEVEL(MESHES(NM)%JBAR,2)
         KLEVEL(3)=SCARC_GET_MAXLEVEL(MESHES(NM)%KBAR,3)

         KLEVEL_MIN = MINVAL(KLEVEL)
         IF (KLEVEL_MIN<NLEVEL) NLEVEL=KLEVEL_MIN

      ENDDO 
      NLEVEL_MIN  = 1
      IF (TYPE_MULTIGRID == NSCARC_MULTIGRID_GEOMETRIC .AND. TYPE_COARSE == NSCARC_COARSE_DIRECT) THEN
         NLEVEL_MAX  = NLEVEL + 1
      ELSE IF (SCARC_MULTIGRID_LEVEL /= -1) THEN
         NLEVEL_MAX  = NLEVEL_MIN + SCARC_MULTIGRID_LEVEL - 1
      ELSE
         NLEVEL_MAX  = NLEVEL
      ENDIF
      NLEVEL = NLEVEL_MAX

   !!!------------------------------------------------------------------------------------------------------
   !!! first, only finest level is set, further levels are defined during coarsening process
   !!!------------------------------------------------------------------------------------------------------
   CASE (NSCARC_MULTIGRID_ALGEBRAIC)
 
      NLEVEL_MIN = 1
      IF (SCARC_MULTIGRID_LEVEL /= -1) THEN
         NLEVEL_MAX  = SCARC_MULTIGRID_LEVEL 
      ELSE
         NLEVEL_MAX  = NSCARC_LEVEL_MAX
      ENDIF
      NLEVEL = SCARC_MULTIGRID_LEVEL

      IF (TYPE_DEBUG > NSCARC_LEVEL_NONE) THEN
         WRITE(SCARC_LU,*) 'NLEVEL_MIN=',NLEVEL_MIN
         WRITE(SCARC_LU,*) 'NLEVEL_MAX=',NLEVEL_MAX
         WRITE(SCARC_LU,*) 'NLEVEL=',NLEVEL
      ENDIF

   !!!------------------------------------------------------------------------------------------------------
   !!! no multigrid-hierachy needed in case of a pure Krylov-method: use only one level
   !!!------------------------------------------------------------------------------------------------------
   CASE DEFAULT

      NLEVEL     = 1
      NLEVEL_MIN = 1
      NLEVEL_MAX = 1

END SELECT

END SUBROUTINE SCARC_SETUP_LEVELS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Determine maximum number of possible levels on direction IOR0 of mesh NM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
INTEGER FUNCTION SCARC_GET_MAXLEVEL(NC, IOR0)
INTEGER, INTENT(IN) :: NC, IOR0
INTEGER :: NC0, NL

!!! In case of the GMG-method, NC must be divisable by 2 at least one time 
IF (MOD(NC,2)/=0 .AND. TYPE_MULTIGRID == NSCARC_MULTIGRID_GEOMETRIC) THEN
   SELECT CASE (IOR0)
      CASE (1)
         WRITE(*,1000) 'IBAR', NC
      CASE (2)
         WRITE(*,1000) 'JBAR', NC
      CASE (3)
         WRITE(*,1000) 'KBAR', NC
   END SELECT
   STOP
ENDIF

!!! divide by 2 as often as possible or till user defined max-level is reached
NC0=NC
DO NL=1,NSCARC_LEVEL_MAX
   NC0=NC0/2
   IF (MOD(NC0,2)/=0) EXIT                  ! NC no longer divisable by two
   IF (NL==SCARC_MULTIGRID_LEVEL) EXIT      ! max number of levels defined by user
   IF (NC0==1) EXIT                         ! NC is power of two, minimum has been reached
ENDDO

SCARC_GET_MAXLEVEL=NL
RETURN
1000 FORMAT(A,'=',I3,' must be divisable by two for ScaRC-Multigrid!')
END FUNCTION SCARC_GET_MAXLEVEL


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Allocate ScaRC-structures for all needed levels
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SETUP_TYPES
INTEGER :: IERR, NM
TYPE (SCARC_TYPE), POINTER :: S

!!! --------------------------------------------------------------------------------------------------------
!!! Allocate global ScaRC-structure
!!! --------------------------------------------------------------------------------------------------------
ALLOCATE (SCARC(NMESHES), STAT=IERR)
CALL CHKMEMERR ('SCARC_SETUP', 'SCARC', IERR)

!!! --------------------------------------------------------------------------------------------------------
!!! Allocate local OSCARC, MESHES, SYSTEM (COMPACT/BANDED) and PRECON structures for single meshes
!!! --------------------------------------------------------------------------------------------------------
MESHES_LOOP: DO NM = NMESHES_MIN, NMESHES_MAX

   S => SCARC(NM)

   ALLOCATE (S%OSCARC(NMESHES), STAT=IERR)
   CALL CHKMEMERR ('SCARC_SETUP_TYPES', 'OSCARC', IERR)
      
   SELECT CASE (TYPE_SYSTEM)
      CASE (NSCARC_SYSTEM_BANDED)
         ALLOCATE (SCARC(NM)%BANDED(NLEVEL_MIN:NLEVEL_MAX), STAT=IERR)
         CALL CHKMEMERR ('SCARC_SETUP_TYPES', 'BANDED', IERR)
      CASE (NSCARC_SYSTEM_COMPACT)
         ALLOCATE (SCARC(NM)%COMPACT(NLEVEL_MIN:NLEVEL_MAX), STAT=IERR)
         CALL CHKMEMERR ('SCARC_SETUP_TYPES', 'COMPACT', IERR)
   END SELECT

   IF (TYPE_PRECON == NSCARC_PRECON_FFT .OR. TYPE_PRECON == NSCARC_PRECON_GSTRIX) THEN
      ALLOCATE (SCARC(NM)%PRECON(NLEVEL_MIN:NLEVEL_MAX), STAT=IERR)
      CALL CHKMEMERR ('SCARC_SETUP_TYPES', 'PRECON', IERR)
   ENDIF

ENDDO MESHES_LOOP

END SUBROUTINE SCARC_SETUP_TYPES


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Setup geometry information for mesh NM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SETUP_MESHES
INTEGER  :: IERR, NL, NM
INTEGER  :: NX0, NY0, NZ0
TYPE (MESH_TYPE)         , POINTER :: M
TYPE (SCARC_TYPE)        , POINTER :: S
TYPE (SCARC_BANDED_TYPE) , POINTER :: SB
TYPE (SCARC_COMPACT_TYPE), POINTER :: SC

IERR=0

SELECT_SYSTEM: SELECT CASE (TYPE_SYSTEM)
 
   !!!-------------------------------------------------------------------------------------------------
   !!! Banded system
   !!!-------------------------------------------------------------------------------------------------
   CASE (NSCARC_SYSTEM_BANDED)
      
      BANDED_MESHES_LOOP: DO NM = NMESHES_MIN, NMESHES_MAX
      
         M => MESHES(NM)
         S => SCARC(NM)
         
         !!!
         !!!  define hierarchy of meshes depending on the chosen method
         !!!
         NX0=MESHES(NM)%IBAR
         NY0=MESHES(NM)%JBAR
         NZ0=MESHES(NM)%KBAR
         
         BANDED_LEVEL_LOOP: DO NL = NLEVEL_MIN, NLEVEL_MAX
         
            IF (NL > NLEVEL_MIN .AND. TYPE_MULTIGRID /= NSCARC_MULTIGRID_GEOMETRIC) EXIT BANDED_LEVEL_LOOP

            !!! let SM point to SCARC(NM)%MESHES(NL)
            SB => S%BANDED(NL)
            
            !!! numbers of cells in x-, y- and z-direction for level 'NL'
            SB%NX = NX0
            SB%NY = NY0
            SB%NZ = NZ0
            
            !!! step widths in x-, y- and z-direction for level 'NL'
            SB%DXI = REAL(SB%NX,EB)/(M%XF-M%XS)
            SB%DYI = REAL(SB%NY,EB)/(M%YF-M%YS)
            SB%DZI = REAL(SB%NZ,EB)/(M%ZF-M%ZS)
            
            SB%DXI2 = SB%DXI**2
            SB%DYI2 = SB%DYI**2
            SB%DZI2 = SB%DZI**2
          
            SB%DI2(1) = SB%DXI2
            SB%DI2(2) = SB%DYI2
            SB%DI2(3) = SB%DZI2

            !!! Get global number of grid cells (internal and including ghost cells)
            SB%NC  = SB%NX * SB%NY * SB%NZ 
            SB%NCE = SB%NC
            SELECT CASE (TYPE_DIMENSION)
               CASE (NSCARC_DIMENSION_TWO)
                  SB%NCG = (SB%NX+2) * (SB%NZ+2) 
               CASE (NSCARC_DIMENSION_THREE)
                  SB%NCG = (SB%NX+2) * (SB%NY+2) * (SB%NZ+2) 
            END SELECT
             
            NX0=NX0/2
            IF (TYPE_DIMENSION == NSCARC_DIMENSION_THREE) NY0=NY0/2
            NZ0=NZ0/2
         
         ENDDO BANDED_LEVEL_LOOP

      ENDDO BANDED_MESHES_LOOP
      
   !!!-------------------------------------------------------------------------------------------------
   !!! Compact system
   !!!-------------------------------------------------------------------------------------------------
   CASE (NSCARC_SYSTEM_COMPACT)

      COMPACT_MESHES_LOOP: DO NM = NMESHES_MIN, NMESHES_MAX
      
         M => MESHES(NM)
         S => SCARC(NM)
         
         !!!
         !!!  define hierarchy of meshes depending on the chosen method
         !!!
         NX0=MESHES(NM)%IBAR
         NY0=MESHES(NM)%JBAR
         NZ0=MESHES(NM)%KBAR
         
         COMPACT_LEVEL_LOOP: DO NL = NLEVEL_MIN, NLEVEL_MAX
         
            IF (NL > NLEVEL_MIN .AND. TYPE_MULTIGRID /= NSCARC_MULTIGRID_GEOMETRIC) EXIT COMPACT_LEVEL_LOOP

            !!! let SM point to SCARC(NM)%MESHES(NL)
            SC => S%COMPACT(NL)
            
            !!! numbers of cells in x-, y- and z-direction for level 'NL'
            SC%NX = NX0
            SC%NY = NY0
            SC%NZ = NZ0
            
            !!! step widths in x-, y- and z-direction for level 'NL'
            SC%DXI = REAL(SC%NX,EB)/(M%XF-M%XS)
            SC%DYI = REAL(SC%NY,EB)/(M%YF-M%YS)
            SC%DZI = REAL(SC%NZ,EB)/(M%ZF-M%ZS)
            
            SC%DXI2 = SC%DXI**2
            SC%DYI2 = SC%DYI**2
            SC%DZI2 = SC%DZI**2
          
            SC%DI2(1) = SC%DXI2
            SC%DI2(2) = SC%DYI2
            SC%DI2(3) = SC%DZI2

            !!! Get global number of grid cells (internal and including ghost cells)
            SC%NC  = SC%NX * SC%NY * SC%NZ 
            SC%NCE = SC%NC
            SELECT CASE (TYPE_DIMENSION)
               CASE (NSCARC_DIMENSION_TWO)
                  SC%NCG = (SC%NX+2) * (SC%NZ+2) 
               CASE (NSCARC_DIMENSION_THREE)
                  SC%NCG = (SC%NX+2) * (SC%NY+2) * (SC%NZ+2) 
            END SELECT
             
            NX0=NX0/2
            IF (TYPE_DIMENSION == NSCARC_DIMENSION_THREE) NY0=NY0/2
            NZ0=NZ0/2
         
         ENDDO COMPACT_LEVEL_LOOP

      ENDDO COMPACT_MESHES_LOOP

END SELECT SELECT_SYSTEM

END SUBROUTINE SCARC_SETUP_MESHES
 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Setup communication structure for data exchange 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SETUP_EXCHANGE 
INTEGER :: NOM, NM, NL, NLMIN, NLMAX, N
INTEGER :: IERR
TYPE (MESH_TYPE)          , POINTER :: M
TYPE (SCARC_TYPE)         , POINTER :: S
TYPE (OMESH_TYPE)         , POINTER :: OM
TYPE (OSCARC_TYPE)        , POINTER :: OS
TYPE (OSCARC_BANDED_TYPE) , POINTER :: OSBF, OSBC
TYPE (OSCARC_COMPACT_TYPE), POINTER :: OSCF, OSCC

IERR = 0

!!! Initialize communication counter for ScaRC, use same TAG for all communications
N_EXCHANGES = 0
TAG_SCARC  = 99

!!!----------------------------------------------------------------------------------------------------
!!! Store communication counters from FDS-code
!!!----------------------------------------------------------------------------------------------------
MESHES_LOOP: DO NM = NMESHES_MIN, NMESHES_MAX
   
   M => MESHES(NM)
   S => SCARC(NM)
   
   !!! Initialize level structures on neighboring meshes
   !!! Note that OSMF%IJKW corresponds to MESHES(NM)%OMESH(NOM)%IJKW
   OTHER_MESHES_LOOP: DO NOM = 1, NMESHES
      
      IF (NOM == NM) CYCLE OTHER_MESHES_LOOP
   
      !!! Search for this neighbor in IJKW to determine the IOR-value of the boundary
  !    FOUND = .FALSE.
  !    SEARCH_NEIGHBOR_LOOP: DO IW = 1, SMF%N_EXTERNAL_WALL_CELLS
  !       IF (SMF%IJKW(9,IW) == NOM) THEN
  !          FOUND = .TRUE.
  !          EXIT SEARCH_NEIGHBOR_LOOP
  !       ENDIF
  !    ENDDO SEARCH_NEIGHBOR_LOOP
  !
  !    IF (.NOT.FOUND) CYCLE OTHER_MESHES_LOOP

      OM => M%OMESH(NOM)
      OS => S%OSCARC(NOM)
   
      OS%I_MIN_S  = OM%I_MIN_S
      OS%I_MAX_S  = OM%I_MAX_S
      OS%J_MIN_S  = OM%J_MIN_S
      OS%J_MAX_S  = OM%J_MAX_S
      OS%K_MIN_S  = OM%K_MIN_S
      OS%K_MAX_S  = OM%K_MAX_S
      OS%NIC_S    = OM%NIC_S
      OS%NICMAX_S = OM%NIC_S
       
      OS%I_MIN_R  = OM%I_MIN_R
      OS%I_MAX_R  = OM%I_MAX_R
      OS%J_MIN_R  = OM%J_MIN_R
      OS%J_MAX_R  = OM%J_MAX_R
      OS%K_MIN_R  = OM%K_MIN_R
      OS%K_MAX_R  = OM%K_MAX_R
      OS%NIC_R    = OM%NIC_R
      OS%NICMAX_R = OM%NIC_R
   
      IF (OS%NICMAX_S==0 .AND. OS%NICMAX_R==0)  CYCLE OTHER_MESHES_LOOP

      N_EXCHANGES  = N_EXCHANGES+1
   
   ENDDO OTHER_MESHES_LOOP
ENDDO MESHES_LOOP

SELECT_SYSTEM: SELECT CASE (TYPE_SYSTEM)

   !!!-------------------------------------------------------------------------------------------------
   !!! Banded system
   !!!-------------------------------------------------------------------------------------------------
   CASE (NSCARC_SYSTEM_BANDED)

      BANDED_MESHES_LOOP: DO NM = NMESHES_MIN, NMESHES_MAX
         
         M => MESHES(NM)
         S => SCARC(NM)
         
         !!! Initialize level structures on neighboring meshes
         !!! Note that OSMF%IJKW corresponds to MESHES(NM)%OMESH(NOM)%IJKW
         BANDED_OTHER_MESHES_LOOP: DO NOM = 1, NMESHES
            
            IF (NOM == NM) CYCLE BANDED_OTHER_MESHES_LOOP
         
            !!! Allocate OSCARC grid structure for mesh NM
            ALLOCATE (S%OSCARC(NOM)%BANDED(NLEVEL_MIN:NLEVEL_MAX), STAT=IERR)
            CALL CHKMEMERR ('SCARC_SETUP_EXCHANGE', 'OS%BANDED', IERR)
         
            !!! point to grid structure of OSCARC(NM) on finest level
            OSBF => S%OSCARC(NOM)%BANDED(NLEVEL_MIN)

            OSBF%IJKW => MESHES(NOM)%OMESH(NOM)%IJKW  

            OSBF%NX =  MESHES(NOM)%IBAR 
            OSBF%NY =  MESHES(NOM)%JBAR 
            OSBF%NZ =  MESHES(NOM)%KBAR 

            OSBF%NW = 2*OSBF%NX*OSBF%NY + 2*OSBF%NX*OSBF%NZ + 2*OSBF%NY*OSBF%NZ  
         
            !!! In case of GMG with a predefined grid hierarchy allocate corresponding level-structures
            IF (NLEVEL_MAX > NLEVEL_MIN .AND. TYPE_MULTIGRID==NSCARC_MULTIGRID_GEOMETRIC) THEN                   
               DO NL=NLEVEL_MIN+1,NLEVEL_MAX
         
                  OSBC => S%OSCARC(NOM)%BANDED(NL)                            ! pointer to coarser level
                  OSBF => S%OSCARC(NOM)%BANDED(NL-1)                          ! pointer to finer level
         
                  ! get number of internal cells and external wall cells on neighbor NOM for level NL
                  OSBC%NX=OSBF%NX/2
                  SELECT CASE (TYPE_DIMENSION)
                     CASE (NSCARC_DIMENSION_TWO)
                        OSBC%NY=1
                     CASE (NSCARC_DIMENSION_THREE)
                        OSBC%NY=OSBF%NY/2
                  END SELECT
                  OSBC%NZ=OSBF%NZ/2
         
                  OSBC%NC = OSBC%NX * OSBC%NY * OSBC%NZ
                  OSBC%NW = 2*OSBC%NX*OSBC%NY + 2*OSBC%NX*OSBC%NZ + 2*OSBC%NY*OSBC%NZ  
         
                 ALLOCATE(OSBC%IJKW(15,OSBC%NW), STAT=IERR)
                  CALL ChkMemErr('SCARC_SETUP_EXCHANGE','IJKW',IERR)
                  OSBC%IJKW = 0
                         
               ENDDO
            ENDIF
               
         ENDDO BANDED_OTHER_MESHES_LOOP
      ENDDO BANDED_MESHES_LOOP
   
   !!!-------------------------------------------------------------------------------------------------
   !!! Compact system
   !!!-------------------------------------------------------------------------------------------------
   CASE (NSCARC_SYSTEM_COMPACT)

      COMPACT_MESHES_LOOP: DO NM = NMESHES_MIN, NMESHES_MAX
         
         M => MESHES(NM)
         S => SCARC(NM)
         
         !!! Initialize level structures on neighboring meshes
         !!! Note that OSMF%IJKW corresponds to MESHES(NM)%OMESH(NOM)%IJKW
         COMPACT_OTHER_MESHES_LOOP: DO NOM = 1, NMESHES
            
            IF (NOM == NM) CYCLE COMPACT_OTHER_MESHES_LOOP
         
            !!! Allocate OSCARC grid structure for mesh NM
            ALLOCATE (S%OSCARC(NOM)%COMPACT(NLEVEL_MIN:NLEVEL_MAX), STAT=IERR)
            CALL CHKMEMERR ('SCARC_SETUP_EXCHANGE', 'OS%COMPACT', IERR)
         
            !!! point to grid structure of OSCARC(NM) on finest level
            OSCF => S%OSCARC(NOM)%COMPACT(NLEVEL_MIN)

            OSCF%IJKW => MESHES(NOM)%OMESH(NOM)%IJKW   !****CHECK

            OSCF%NX = MESHES(NOM)%IBAR 
            OSCF%NY = MESHES(NOM)%JBAR 
            OSCF%NZ = MESHES(NOM)%KBAR 

            OSCF%NW = 2*OSCF%NX*OSCF%NY + 2*OSCF%NX*OSCF%NZ + 2*OSCF%NY*OSCF%NZ  
         
            ALLOCATE(OSCF%IJKW(15,OSCF%NW), STAT=IERR)
            CALL ChkMemErr('SCARC_SETUP_EXCHANGE','IJKW',IERR)
            OSCF%IJKW = 0
         
            !!! In case of GMG with a predefined grid hierarchy allocate corresponding level-structures
            IF (NLEVEL_MAX > NLEVEL_MIN .AND. TYPE_MULTIGRID==NSCARC_MULTIGRID_GEOMETRIC) THEN                   
               DO NL=NLEVEL_MIN+1,NLEVEL_MAX
         
                  OSCC => S%OSCARC(NOM)%COMPACT(NL)                            ! pointer to coarser level
                  OSCF => S%OSCARC(NOM)%COMPACT(NL-1)                          ! pointer to finer level
         
                  ! get number of internal cells and external wall cells on neighbor NOM for level NL
                  OSCC%NX=OSCF%NX/2
                  SELECT CASE (TYPE_DIMENSION)
                     CASE (NSCARC_DIMENSION_TWO)
                        OSCC%NY=1
                     CASE (NSCARC_DIMENSION_THREE)
                        OSCC%NY=OSCF%NY/2
                  END SELECT
                  OSCC%NZ=OSCF%NZ/2
         
                  OSCC%NC=OSCC%NX * OSCC%NY * OSCC%NZ
                  OSCC%NW= 2*OSCC%NX*OSCC%NY + 2*OSCC%NX*OSCC%NZ + 2*OSCC%NY*OSCC%NZ  
         
                  ALLOCATE(OSCC%IJKW(15,OSCC%NW), STAT=IERR)
                  CALL ChkMemErr('SCARC_SETUP_EXCHANGE','IJKW',IERR)
                  OSCC%IJKW = 0
               
               ENDDO
            ENDIF
               
         ENDDO COMPACT_OTHER_MESHES_LOOP
      ENDDO COMPACT_MESHES_LOOP
   
END SELECT SELECT_SYSTEM

!!!
!!! Allocate request array for data exchanges
!!!
IERR = 0
IF (NMESHES>1) THEN
   ALLOCATE (REQ_SCARC(N_EXCHANGES*40))
   CALL CHKMEMERR ('SCARC_SETUP_GLOBAL', 'REQ_SCARC', IERR)
   REQ_SCARC = MPI_REQUEST_NULL
ENDIF

!!!
!!! Allocate counter and displacement vector for global data exchanges
!!!
ALLOCATE(COUNTS_SCARC(0:NUMPROCS-1))
ALLOCATE(DISPLS_SCARC(0:NUMPROCS-1))

COUNTS_SCARC = 0
DO N=0,NUMPROCS-1
   DO NM=1,NMESHES
      IF (PROCESS(NM)==N) COUNTS_SCARC(N) = COUNTS_SCARC(N) + 1
   ENDDO
ENDDO
DISPLS_SCARC(0) = 0
DO N=1,NUMPROCS-1
   DISPLS_SCARC(N) = COUNTS_SCARC(N-1) + DISPLS_SCARC(N-1)
ENDDO


!!!-------------------------------------------------------------------------------------------------------
!!! Initialize communication structures on finest level (if there is more than 1 mesh) 
!!!-------------------------------------------------------------------------------------------------------
IF (NMESHES>1) THEN

   SELECT CASE (TYPE_MULTIGRID)
      CASE (NSCARC_MULTIGRID_GEOMETRIC)
         NLMIN = NLEVEL_MIN
         NLMAX = NLEVEL_MAX
      CASE DEFAULT
         NLMIN = NLEVEL_MIN
         NLMAX = NLEVEL_MIN
   END SELECT


   DO NL = NLMIN, NLMAX
      CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_IJKW, NL)
   ENDDO

ENDIF

END SUBROUTINE SCARC_SETUP_EXCHANGE

 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Setup neighborship structures and boundary conditions on finest level
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SETUP_WALLCELLS
USE GEOMETRY_FUNCTIONS, ONLY: SEARCH_OTHER_MESHES
INTEGER :: NM, NOM, NL
INTEGER :: IOFFSET, IXC, IYC, IZC, IWF, IWC, IREFINE, IERR, IOR0, IOR_LAST, INBR, INBR0, INBR1
INTEGER :: NEIGHBORS(20,-3:3)
LOGICAL :: FOUND
TYPE (MESH_TYPE) , POINTER :: M
TYPE (SCARC_TYPE), POINTER :: S
TYPE (SCARC_BANDED_TYPE) , POINTER :: SBF, SBC
TYPE (SCARC_COMPACT_TYPE), POINTER :: SCF, SCC

IERR=0

SELECT_SYSTEM: SELECT CASE (TYPE_SYSTEM)
 
   !!!----------------------------------------------------------------------------------------------------
   !!! Banded system
   !!!----------------------------------------------------------------------------------------------------
   CASE (NSCARC_SYSTEM_BANDED)

      MESHES_LOOP_BANDED: DO NM = NMESHES_MIN, NMESHES_MAX
      
         M => MESHES(NM)
         S => SCARC(NM)
            
         !!! For all solver: 
         !!! Determine array IJKW and PRESSURE_BC_INDEX on finest level
         SBF => S%BANDED(NLEVEL_MIN)
         
         SBF%IJKW     => M%OMESH(NM)%IJKW   !****CHECK
         SBF%BC_INDEX => M%PRESSURE_BC_INDEX

         SBF%NW =  M%N_EXTERNAL_WALL_CELLS

         ALLOCATE(SBF%ADJACENT_CELL(SBF%NW), STAT=IERR)
         CALL ChkMemErr('SCARC_SETUP_WALLCELLS','ADJACENT_CELL',IERR)
         SBF%ADJACENT_CELL = 0

         ALLOCATE(SBF%GHOST_CELL(SBF%NW), STAT=IERR)
         CALL ChkMemErr('SCARC_SETUP_WALLCELLS','GHOST_CELL',IERR)
         SBF%GHOST_CELL = 0

         DO IWF = 1, SBF%NW
      
            !!! Determine boundary type for IWF
            IF (M%WALL(IWF)%BOUNDARY_TYPE == OPEN_BOUNDARY) THEN
               SBF%BC_INDEX(IWF) = DIRICHLET
            ELSE IF (M%OMESH(NM)%IJKW(9,IWF) /= 0) THEN !****CHECK
               SBF%BC_INDEX(IWF) = INTERNAL
            !ELSE IF (M%BOUNDARY_TYPE(IWF) == NULL_BOUNDARY) THEN
            !   SBF%BC_INDEX(IWF) = DIRICHLET
            ELSE
               SBF%BC_INDEX(IWF) = NEUMANN
            ENDIF

         ENDDO

         CALL SCARC_SETUP_SUBDIVISION(SBF%IJKW, SBF%SUBDIVISION, SBF%NW, NM)
         
         
         !!! Only in case of MG-method:
         !!! Determine arrays IJKW for coarser levels
         BANDED_GMG_IF: IF (TYPE_MULTIGRID == NSCARC_MULTIGRID_GEOMETRIC) THEN
         
            IREFINE=1
            BANDED_GMG_LEVEL_LOOP: DO NL = NLEVEL_MIN+1, NLEVEL_MAX
            
               !!! point to SCARC grid structures on coarser and finer level
               SBC => S%BANDED(NL)
               SBF => S%BANDED(NL-1)
         
               IREFINE=IREFINE*2

               SBC%NW = 2*SBC%NX*SBC%NY + 2*SBC%NX*SBC%NZ + 2*SBC%NY*SBC%NZ

               ALLOCATE(SBC%BC_INDEX(SBC%NW), STAT=IERR)
               CALL ChkMemErr('SCARC_SETUP_WALLCELLS','BC_INDEX',IERR)
               SBC%BC_INDEX = 0
      
               ALLOCATE(SBC%IJKW(15,SBC%NW), STAT=IERR)
               CALL ChkMemErr('SCARC_SETUP_WALLCELLS','IJKW',IERR)
               SBC%IJKW = 0
            
               ALLOCATE(SBC%ADJACENT_CELL(SBC%NW), STAT=IERR)
               CALL ChkMemErr('SCARC_SETUP_WALLCELLS','ADJACENT_CELL',IERR)
               SBC%ADJACENT_CELL = 0
            
               ALLOCATE(SBC%GHOST_CELL(SBC%NW), STAT=IERR)
               CALL ChkMemErr('SCARC_SETUP_WALLCELLS','GHOST_CELL',IERR)
               SBC%GHOST_CELL = 0
            
               !!! 
               !!! set wall cells for coarser grid and define corresponding IJKW
               !!! 
               IWC=1
         
               !!! wall cells along IOR=1
               IOFFSET = 0
               DO IZC=1,SBC%NZ
                  DO IYC=1,SBC%NY
                     CALL SCARC_SETUP_FACE(IWC,  1, IOFFSET, IREFINE, 0, IYC, IZC, SBF%NY, NM, NL)
                  ENDDO
               ENDDO
         
               !!! wall cells along IOR=-1
              IOFFSET = SBF%NY*SBF%NZ
               DO IZC=1,SBC%NZ
                  DO IYC=1,SBC%NY
                     CALL SCARC_SETUP_FACE(IWC, -1, IOFFSET, IREFINE, SBC%NX+1, IYC, IZC, SBF%NY, NM, NL)
                  ENDDO
               ENDDO
            
               !!! wall cells along IOR=2
               IOFFSET = 2*SBF%NY*SBF%NZ
               DO IZC=1,SBC%NZ
                  DO IXC=1,SBC%NX
                     CALL SCARC_SETUP_FACE(IWC,  2, IOFFSET, IREFINE, IXC, 0, IZC, SBF%NX, NM, NL)
                  ENDDO
               ENDDO
         
               !!! wall cells along IOR=-2
               IOFFSET = 2*SBF%NY*SBF%NZ + SBF%NX*SBF%NZ
               DO IZC=1,SBC%NZ
                  DO IXC=1,SBC%NX
                     CALL SCARC_SETUP_FACE(IWC, -2, IOFFSET, IREFINE, IXC, SBC%NY+1, IZC, SBF%NX, NM, NL)
                  ENDDO
               ENDDO
         
               !!! wall cells along IOR=3
               IOFFSET = 2*SBF%NY*SBF%NZ + 2*SBF%NX*SBF%NZ
               DO IYC=1,SBC%NY
                  DO IXC=1,SBC%NX
                     CALL SCARC_SETUP_FACE(IWC,  3, IOFFSET, IREFINE, IXC, IYC, 0, SBF%NX, NM, NL)
                  ENDDO
               ENDDO
         
               !!! wall cells along IOR=-3
               IOFFSET = 2*SBF%NY*SBF%NZ + 2*SBF%NX*SBF%NZ + SBF%NX*SBF%NY
               DO IYC=1,SBC%NY
                  DO IXC=1,SBC%NX
                     CALL SCARC_SETUP_FACE(IWC, -3, IOFFSET, IREFINE, IXC, IYC, SBC%NZ+1, SBF%NX, NM, NL)
                  ENDDO
               ENDDO
         
               !!!
               !!! compute analogues to NIC, I_MIN, I_MAX, K_MIN, K_MAX on level 'NL'
               !!!
      
               !CALL SCARC_SETUP_COARSE_DIMENSIONS(IREFINE, NM, NL)
         
               !!! Store subdivision information on coarser level
               CALL SCARC_SETUP_SUBDIVISION(SBC%IJKW, SBC%SUBDIVISION, SBC%NW, NM)
         
            ENDDO BANDED_GMG_LEVEL_LOOP
         ENDIF BANDED_GMG_IF
      
      ENDDO MESHES_LOOP_BANDED

   !!!----------------------------------------------------------------------------------------------------
   !!! Compact system
   !!!----------------------------------------------------------------------------------------------------
   CASE (NSCARC_SYSTEM_COMPACT)

      MESHES_LOOP_COMPACT: DO NM = NMESHES_MIN, NMESHES_MAX

         M => MESHES(NM)
         S => SCARC(NM)
            
         !!! For all solver: 
         !!! Determine array IJKW and PRESSURE_BC_INDEX on finest level
         SCF => S%COMPACT(NLEVEL_MIN)
         
         SCF%IJKW     => M%OMESH(NM)%IJKW!****CHECK
         SCF%BC_INDEX => M%PRESSURE_BC_INDEX

         SCF%NW =  M%N_EXTERNAL_WALL_CELLS

         ALLOCATE(SCF%ADJACENT_CELL(SCF%NW), STAT=IERR)
         CALL ChkMemErr('SCARC_SETUP_WALLCELLS','ADJACENT_CELL',IERR)
         SCF%ADJACENT_CELL = 0

         ALLOCATE(SCF%GHOST_CELL(SCF%NW), STAT=IERR)
         CALL ChkMemErr('SCARC_SETUP_WALLCELLS','GHOST_CELL',IERR)
         SCF%GHOST_CELL = 0

         DO IWF = 1, SCF%NW
      
            !!! Determine boundary type for IWF
            IF (M%WALL(IWF)%BOUNDARY_TYPE == OPEN_BOUNDARY) THEN
               SCF%BC_INDEX(IWF) = DIRICHLET
            ELSE IF (M%OMESH(NM)%IJKW(9,IWF) /= 0) THEN !****CHECK
               SCF%BC_INDEX(IWF) = INTERNAL
            !ELSE IF (M%BOUNDARY_TYPE(IWF) == NULL_BOUNDARY) THEN
            !   SCF%BC_INDEX(IWF) = DIRICHLET
            ELSE
               SCF%BC_INDEX(IWF) = NEUMANN
            ENDIF

         ENDDO

         CALL SCARC_SETUP_SUBDIVISION(SCF%IJKW, SCF%SUBDIVISION, SCF%NW, NM)


         !!! Only in case of MG-method:
         !!! Determine arrays IJKW for coarser levels
         COMPACT_GMG_IF: IF (TYPE_MULTIGRID == NSCARC_MULTIGRID_GEOMETRIC) THEN
         
            IREFINE=1
            COMPACT_GMG_LEVEL_LOOP: DO NL = NLEVEL_MIN+1, NLEVEL_MAX
            
               !!! point to SCARC grid structures on coarser and finer level
               SCC => S%COMPACT(NL)
               SCF => S%COMPACT(NL-1)
         
               IREFINE=IREFINE*2

               SCC%NW = 2*SCC%NX*SCC%NY + 2*SCC%NX*SCC%NZ + 2*SCC%NY*SCC%NZ

               ALLOCATE(SCC%BC_INDEX(SCC%NW), STAT=IERR)
               CALL ChkMemErr('SCARC_SETUP_WALLCELLS','BC_INDEX',IERR)
               SCC%BC_INDEX = 0
      
               ALLOCATE(SCC%IJKW(15,SCC%NW), STAT=IERR)
               CALL ChkMemErr('SCARC_SETUP_WALLCELLS','IJKW',IERR)
               SCC%IJKW = 0
            
               ALLOCATE(SCC%ADJACENT_CELL(SCC%NW), STAT=IERR)
               CALL ChkMemErr('SCARC_SETUP_WALLCELLS','ADJACENT_CELL',IERR)
               SCC%ADJACENT_CELL = 0
            
               ALLOCATE(SCC%GHOST_CELL(SCC%NW), STAT=IERR)
               CALL ChkMemErr('SCARC_SETUP_WALLCELLS','GHOST_CELL',IERR)
               SCC%GHOST_CELL = 0
            
               !!! 
               !!! set wall cells for coarser grid and define corresponding IJKW
               !!! 
               IWC=1
         
               !!! wall cells along IOR=1
               IOFFSET = 0
               DO IZC=1,SCC%NZ
                  DO IYC=1,SCC%NY
                     CALL SCARC_SETUP_FACE(IWC,  1, IOFFSET, IREFINE, 0, IYC, IZC, SCF%NY, NM, NL)
                  ENDDO
               ENDDO
         
               !!! wall cells along IOR=-1
              IOFFSET = SCF%NY*SCF%NZ
               DO IZC=1,SCC%NZ
                  DO IYC=1,SCC%NY
                     CALL SCARC_SETUP_FACE(IWC, -1, IOFFSET, IREFINE, SCC%NX+1, IYC, IZC, SCF%NY, NM, NL)
                  ENDDO
               ENDDO
            
               !!! wall cells along IOR=2
               IOFFSET = 2*SCF%NY*SCF%NZ
               DO IZC=1,SCC%NZ
                  DO IXC=1,SCC%NX
                     CALL SCARC_SETUP_FACE(IWC,  2, IOFFSET, IREFINE, IXC, 0, IZC, SCF%NX, NM, NL)
                  ENDDO
               ENDDO
         
               !!! wall cells along IOR=-2
               IOFFSET = 2*SCF%NY*SCF%NZ + SCF%NX*SCF%NZ
               DO IZC=1,SCC%NZ
                  DO IXC=1,SCC%NX
                     CALL SCARC_SETUP_FACE(IWC, -2, IOFFSET, IREFINE, IXC, SCC%NY+1, IZC, SCF%NX, NM, NL)
                  ENDDO
               ENDDO
         
               !!! wall cells along IOR=3
               IOFFSET = 2*SCF%NY*SCF%NZ + 2*SCF%NX*SCF%NZ
               DO IYC=1,SCC%NY
                  DO IXC=1,SCC%NX
                     CALL SCARC_SETUP_FACE(IWC,  3, IOFFSET, IREFINE, IXC, IYC, 0, SCF%NX, NM, NL)
                  ENDDO
               ENDDO
         
               !!! wall cells along IOR=-3
               IOFFSET = 2*SCF%NY*SCF%NZ + 2*SCF%NX*SCF%NZ + SCF%NX*SCF%NY
               DO IYC=1,SCC%NY
                  DO IXC=1,SCC%NX
                     CALL SCARC_SETUP_FACE(IWC, -3, IOFFSET, IREFINE, IXC, IYC, SCC%NZ+1, SCF%NX, NM, NL)
                  ENDDO
               ENDDO
         
               !!!
               !!! compute analogues to NIC, I_MIN, I_MAX, K_MIN, K_MAX on level 'NL'
               !!!
      
               !CALL SCARC_SETUP_COARSE_DIMENSIONS(IREFINE, NM, NL)
         
               !!! Store subdivision information on coarser level
               CALL SCARC_SETUP_SUBDIVISION(SCC%IJKW, SCC%SUBDIVISION, SCC%NW, NM)
         
            ENDDO COMPACT_GMG_LEVEL_LOOP
         ENDIF COMPACT_GMG_IF
      
         
      ENDDO MESHES_LOOP_COMPACT

END SELECT SELECT_SYSTEM

END SUBROUTINE SCARC_SETUP_WALLCELLS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Store subdivision information
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
SUBROUTINE SCARC_SETUP_SUBDIVISION(IJKW, SUBDIVISION, NW, NM)
INTEGER, DIMENSION(:,:), INTENT(IN)    :: IJKW
INTEGER, DIMENSION(3,-3:3), INTENT(OUT) :: SUBDIVISION
INTEGER, INTENT(IN) :: NW, NM
INTEGER :: IW, IOR0, IOR_LAST, NOM, INBR
INTEGER :: NEIGHBORS(20,-3:3)

IOR_LAST    = 0
NEIGHBORS   = 0
SUBDIVISION = 0
         
WALLCELL_LOOP: DO IW = 1, NW

   IOR0 = IJKW(4,IW)
   
   IF (IOR_LAST /= IOR0) SUBDIVISION(1,IOR0) = IW
!   IF (IOR_LAST /= IOR0) WRITE(SCARC_LU,'(a,i3,a,i3,a,i3,a,2i3)') &
!         '============NM=',NM,': NOM=',NOM,': SUBDIVISION(1,',IOR0,')=',SUBDIVISION(1,IOR0), IOR0
   SUBDIVISION(2,IOR0) = SUBDIVISION(2,IOR0) + 1
   
   NOM= IJKW(9,IW)
!WRITE(SCARC_LU,'(a,i3,a,i3,a,i3,a,3i3)') '    NM=',NM,': NOM=',NOM,': SUBDIVISION(2,',IOR0,')=',SUBDIVISION(2,IOR0), IOR0, SUBDIVISION(2,-2)

   IF (NOM /= 0) THEN
      NEIGHBOR_LOOP: DO INBR = 1, 20
!WRITE(SCARC_LU,'(10i3)') NEIGHBORS(1:10,IOR0)
         IF (NOM == NEIGHBORS(INBR, IOR0)) THEN
            EXIT NEIGHBOR_LOOP
         ELSE IF (NEIGHBORS(INBR, IOR0) /= 0) THEN
            CYCLE NEIGHBOR_LOOP
         ELSE IF (NEIGHBORS(INBR, IOR0) == 0) THEN
            NEIGHBORS(INBR, IOR0) = NOM
            SUBDIVISION(3,IOR0) = SUBDIVISION(3,IOR0) + 1
!WRITE(SCARC_LU,'(a,i3,a,i3,a,i3,a,3i3)') '--->NM=',NM,': NOM=',NOM,': SUBDIVISION(3,',IOR0,')=',SUBDIVISION(3,IOR0),IOR0,SUBDIVISION(2,-2)
            EXIT NEIGHBOR_LOOP
         ELSE
            WRITE(*,*) 'More than 20 neighbors at one face not allowed yet!'
            STOP
         ENDIF
      ENDDO NEIGHBOR_LOOP
   ENDIF
   
   IOR_LAST = IOR0

ENDDO WALLCELL_LOOP
         
!WRITE(SCARC_LU,'(a,7i5)') 'SUBDIVISION(1,:)',SUBDIVISION(1,-3:3)
!WRITE(SCARC_LU,'(a,7i5)') 'SUBDIVISION(2,:)',SUBDIVISION(2,-3:3)
!WRITE(SCARC_LU,'(a,7i5)') 'SUBDIVISION(3,:)',SUBDIVISION(3,-3:3)
END SUBROUTINE SCARC_SETUP_SUBDIVISION


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Set wall cell information on coarse level
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SETUP_FACE(IW_CO, IOR0, IOFFSET, IREFINE, IP, JP, KP, ILEN_FI, NM, NL)
INTEGER, INTENT(INOUT) :: IW_CO
INTEGER, INTENT(IN) :: IOR0, IOFFSET, IREFINE, NM, NL
INTEGER, INTENT(IN) :: IP, JP, KP, ILEN_FI
INTEGER :: IW_FI(4) , IBC_FI(4), NOM_FI(4)
INTEGER :: I, I1, I2, J1, J2, K1, K2
INTEGER :: IDIFF, JDIFF, KDIFF
INTEGER, DIMENSION(:)  , POINTER:: BC_CO, BC_FI
INTEGER, DIMENSION(:,:), POINTER:: IJKW_CO, IJKW_FI

SELECT CASE (TYPE_SYSTEM)
   CASE (NSCARC_SYSTEM_BANDED)
      IJKW_CO => SCARC(NM)%BANDED(NL)%IJKW
      IJKW_FI => SCARC(NM)%BANDED(NL-1)%IJKW
      BC_CO   => SCARC(NM)%BANDED(NL)%BC_INDEX
      BC_FI   => SCARC(NM)%BANDED(NL-1)%BC_INDEX
   CASE (NSCARC_SYSTEM_COMPACT)
      IJKW_CO => SCARC(NM)%COMPACT(NL)%IJKW
      IJKW_FI => SCARC(NM)%COMPACT(NL-1)%IJKW
      BC_CO   => SCARC(NM)%COMPACT(NL)%BC_INDEX
      BC_FI   => SCARC(NM)%COMPACT(NL-1)%BC_INDEX
END SELECT

!!! Set orientation of neiboring face, indices of ghost and adjacent cell for coarse IW
!!!
IJKW_CO (4, IW_CO) = IOR0

IJKW_CO (1, IW_CO) = IP
IJKW_CO (2, IW_CO) = JP
IJKW_CO (3, IW_CO) = KP

SELECT CASE (IOR0)
   CASE (1)
      IJKW_CO (6, IW_CO) = IP+1
      IJKW_CO (7, IW_CO) = JP
      IJKW_CO (8, IW_CO) = KP
   CASE (-1)
      IJKW_CO (6, IW_CO) = IP-1
      IJKW_CO (7, IW_CO) = JP
      IJKW_CO (8, IW_CO) = KP
   CASE (2)
      IJKW_CO (6, IW_CO) = IP
      IJKW_CO (7, IW_CO) = JP+1
      IJKW_CO (8, IW_CO) = KP
   CASE (-2)
      IJKW_CO (6, IW_CO) = IP
      IJKW_CO (7, IW_CO) = JP-1
      IJKW_CO (8, IW_CO) = KP
   CASE (3)
      IJKW_CO (6, IW_CO) = IP
      IJKW_CO (7, IW_CO) = JP
      IJKW_CO (8, IW_CO) = KP+1
   CASE (-3)
      IJKW_CO (6, IW_CO) = IP
      IJKW_CO (7, IW_CO) = JP
      IJKW_CO (8, IW_CO) = KP-1
END SELECT


SELECT CASE (TYPE_DIMENSION)

   !!!----------------------------------------------------------------------------------------------------
   !!! 2D-version
   !!!----------------------------------------------------------------------------------------------------
   CASE (NSCARC_DIMENSION_TWO)

      !!! set IJKW(1:8,IW_CO) for coarser grid IW_CO
      SELECT CASE (ABS(IOR0))
         CASE ( 1)
            IW_FI(1) = IOFFSET +  2*KP-1
         CASE ( 2)
            IW_FI(1) = IOFFSET + (2*KP-1)*ILEN_FI + 2*IP - 1
         CASE ( 3)
            IW_FI(1) = IOFFSET +  2*IP-1
      END SELECT
      IW_FI(2) = IW_FI(1)+1
   
      !!! set neighbors IJKW(9,IW_CO) for coarser grid IW_CO
      NOM_FI(1) = IJKW_FI (9, IW_FI(1))
      NOM_FI(2) = IJKW_FI (9, IW_FI(2))
      IF (NOM_FI(1) /= NOM_FI(2)) THEN
         WRITE(*,*) 'SCARC_SETUP_FACE: Inconsistent neighbors on IOR=', IOR0,' not allowed!'
         STOP
      ENDIF
   
      IJKW_CO (9, IW_CO)=NOM_FI(1) 
   
      !!! set corresponding pressure_bc_index on coarser level
      IBC_FI(1) = BC_FI(IW_FI(1))
      IBC_FI(2) = BC_FI(IW_FI(2))
      IF (IBC_FI(1) == INTERNAL .OR. IBC_FI(2) == INTERNAL) THEN
         BC_CO(IW_CO) = INTERNAL
      ELSE IF (IBC_FI(1) == DIRICHLET .OR. IBC_FI(2) == DIRICHLET) THEN
         BC_CO(IW_CO) = DIRICHLET
      ELSE
         BC_CO(IW_CO) = NEUMANN
      ENDIF
   
      !!! in case of an internal boundary set IJKW(10:15,IW_CO)
      IF (NOM_FI(1) > 0) THEN   
   
         J1 = 1
         J2 = 1
         SELECT CASE (ABS(IOR0))
            CASE (1)
               KDIFF = IJKW_FI (12, IW_FI(2)) - IJKW_FI (12, IW_FI(1))
               IF (KDIFF == 1) THEN
                  K1 = IJKW_FI (15, IW_FI(2))/2
                  K2 = K1
               ELSE IF (KDIFF == 2) THEN
                  K1 = IJKW_FI (15, IW_FI(1))/2     
                  K2 = IJKW_FI (15, IW_FI(2))/2    
               ELSE IF (KDIFF == 0) THEN
                  K1 = (IJKW_FI (15, IW_FI(1))+1)/2
                  K2 = K1
               ELSE
                  WRITE(*,*) 'WRONG resolutions in SCARC_SETUP_FACE, IOR0=',IOR0
                  STOP
               ENDIF
            CASE (3)
               IDIFF = IJKW_FI (10, IW_FI(2)) - IJKW_FI (10, IW_FI(1))
               IF (IDIFF == 1) THEN
                  I1 = IJKW_FI (13, IW_FI(2))/2
                  I2 = I1
               ELSE IF (IDIFF == 2) THEN
                  I1 = IJKW_FI (13, IW_FI(1))/2
                  I1 = IJKW_FI (13, IW_FI(2))/2
               ELSE IF (IDIFF == 0) THEN
                  I1 = (IJKW_FI (13, IW_FI(1))+1)/2
                  I2 = I1
               ELSE
                  WRITE(*,*) 'WRONG resolutions in SCARC_SETUP_FACE, IOR0=',IOR0
                  STOP
               ENDIF
         END SELECT
   
         SELECT CASE (IOR0)
            CASE (1)
               I1 = MESHES(NOM_FI(1))%IBAR/IREFINE
               I2 = I1
            CASE (-1)
               I1 = 1
               I2 = 1
            CASE (3)
               K1 = MESHES(NOM_FI(1))%KBAR/IREFINE
               K2 = K1
            CASE (-3)
               K1 = 1
               K2 = 1
         END SELECT
   
         !!!
         !!! Set ranges for data exchange
         !!!
         IJKW_CO (10, IW_CO) = I1
         IJKW_CO (11, IW_CO) = J1
         IJKW_CO (12, IW_CO) = K1
         IJKW_CO (13, IW_CO) = I2
         IJKW_CO (14, IW_CO) = J2
         IJKW_CO (15, IW_CO) = K2
   
      ENDIF
         

!!!----------------------------------------------------------------------------------------------------
!!! 3D-version
!!!----------------------------------------------------------------------------------------------------
   CASE (NSCARC_DIMENSION_THREE)

      !!! set IJKW(1:8,IW_CO) for coarser grid IW_CO
      SELECT CASE (ABS(IOR0))
         CASE (1)
            IW_FI(1) = IOFFSET + (2*KP-2)*ILEN_FI + 2*JP - 1
            IW_FI(3) = IOFFSET + (2*KP-1)*ILEN_FI + 2*JP - 1
         CASE (2)
            IW_FI(1) = IOFFSET + (2*KP-2)*ILEN_FI + 2*IP - 1
            IW_FI(3) = IOFFSET + (2*KP-1)*ILEN_FI + 2*IP - 1
         CASE (3)
            IW_FI(1) = IOFFSET + (2*JP-2)*ILEN_FI + 2*IP - 1
            IW_FI(3) = IOFFSET + (2*JP-1)*ILEN_FI + 2*IP - 1
      END SELECT
      IW_FI(2) = IW_FI(1)+1
      IW_FI(4) = IW_FI(3)+1
   
      !!! set neighbors IJKW(9,IW_CO) for coarser grid IW_CO
      DO I=1,4
         NOM_FI(I) = IJKW_FI(9,IW_FI(I))
      ENDDO
         
      IF (NOM_FI(1)/=NOM_FI(2) .OR. NOM_FI(1)/=NOM_FI(3) .OR. NOM_FI(1)/=NOM_FI(4)) THEN
         WRITE(*,*) 'SCARC_SETUP_FACE: Inconsistent neighbors on IOR=', IOR0,' not allowed!'
         STOP
      ENDIF
      IJKW_CO (9, IW_CO)=NOM_FI(1) 
   
      !!! set corresponding pressure_bc_index on coarser level
      DO I=1,4
         IBC_FI(I) = BC_FI(IW_FI(I))
      ENDDO
      IF (IBC_FI(1)==INTERNAL.OR.IBC_FI(2)==INTERNAL.OR.&
          IBC_FI(3)==INTERNAL.OR.IBC_FI(4)==INTERNAL) THEN
         BC_CO(IW_CO)=INTERNAL
      ELSE IF (IBC_FI(1)==DIRICHLET.OR.IBC_FI(2)==DIRICHLET.OR.&
               IBC_FI(3)==DIRICHLET.OR.IBC_FI(4)==DIRICHLET) THEN
         BC_CO(IW_CO)=DIRICHLET
      ELSE
         BC_CO(IW_CO)=NEUMANN
      ENDIF
   
      !!! in case of an internal boundary set IJKW(10:15,IW_CO)
      IF (NOM_FI(1) > 0) THEN   
   
         SELECT CASE (ABS(IOR0))
            CASE (1)
               JDIFF = IJKW_FI (11, IW_FI(2)) - IJKW_FI (11, IW_FI(1))
               KDIFF = IJKW_FI (12, IW_FI(3)) - IJKW_FI (12, IW_FI(1))
               IF (JDIFF==1 .AND. KDIFF==1) THEN
                  J1 = IJKW_FI (14, IW_FI(2))/2
                  J2 = J1
                  K1 = IJKW_FI (15, IW_FI(3))/2
                  K2 = K1
               ELSE IF (JDIFF==2 .AND. KDIFF==2) THEN
                  J1 = IJKW_FI (14, IW_FI(1))/2
                  J2 = IJKW_FI (14, IW_FI(2))/2
                  K1 = IJKW_FI (15, IW_FI(1))/2
                  K1 = IJKW_FI (15, IW_FI(3))/2
               ELSE IF (JDIFF==0 .AND. KDIFF==0) THEN
                  J1 = IJKW_FI (14, IW_FI(1))/2
                  J2 = J1
                  K1 = IJKW_FI (15, IW_FI(1))/2
                  K2 = K1
               ELSE
                  WRITE(*,*) 'WRONG resolutions in SCARC_SETUP_FACE, IOR0=',IOR0
                  STOP
               ENDIF
            CASE (2)
               IDIFF = IJKW_FI (10, IW_FI(2)) - IJKW_FI (10, IW_FI(1))
               KDIFF = IJKW_FI (12, IW_FI(3)) - IJKW_FI (12, IW_FI(1))
               IF (IDIFF==1 .AND. KDIFF==1) THEN
                  I1 = IJKW_FI (13, IW_FI(2))/2
                  I2 = I1
                  K1 = IJKW_FI (15, IW_FI(3))/2
                  K2 = K1
               ELSE IF (IDIFF==2 .AND. KDIFF==2) THEN
                  I1 = IJKW_FI (13, IW_FI(1))/2
                  I2 = IJKW_FI (13, IW_FI(2))/2
                  K1 = IJKW_FI (15, IW_FI(1))/2
                  K1 = IJKW_FI (15, IW_FI(3))/2
               ELSE IF (IDIFF==0 .AND. KDIFF==0) THEN
                  I1 = IJKW_FI (13, IW_FI(1))/2
                  I2 = I1
                  K1 = IJKW_FI (15, IW_FI(1))/2
                  K2 = K1
               ELSE
                  WRITE(*,*) 'WRONG resolutions in SCARC_SETUP_FACE, IOR0=',IOR0
                  STOP
               ENDIF
            CASE (3)
               IDIFF = IJKW_FI (10, IW_FI(2)) - IJKW_FI (10, IW_FI(1))
               JDIFF = IJKW_FI (11, IW_FI(3)) - IJKW_FI (11, IW_FI(1))
               IF (IDIFF==1 .AND. JDIFF==1) THEN
                  I1 = IJKW_FI (13, IW_FI(2))/2
                  I2 = I1
                  J1 = IJKW_FI (14, IW_FI(3))/2
                  J2 = J1
               ELSE IF (IDIFF==2 .AND. JDIFF==2) THEN
                  I1 = IJKW_FI (13, IW_FI(1))/2
                  I2 = IJKW_FI (13, IW_FI(2))/2
                  J1 = IJKW_FI (14, IW_FI(1))/2
                  J1 = IJKW_FI (14, IW_FI(3))/2
               ELSE IF (IDIFF==0 .AND. JDIFF==0) THEN
                  I1 = IJKW_FI (13, IW_FI(2))/2
                  I2 = I1
                  J1 = IJKW_FI (14, IW_FI(3))/2
                  J2 = J1
               ELSE
                  WRITE(*,*) 'WRONG resolutions in SCARC_SETUP_FACE, IOR0=',IOR0
                  STOP
               ENDIF
         END SELECT
   
         SELECT CASE (IOR0)
            CASE (1)
               I1 = MESHES(NOM_FI(1))%IBAR/IREFINE
               I2 = I1
            CASE (-1)
               I1 = 1
               I2 = I1
            CASE (2)
               J1 = MESHES(NOM_FI(1))%JBAR/IREFINE
               J2 = J1
            CASE (-2)
               J1 = 1
               J2 = J1
            CASE (3)
               K1 = MESHES(NOM_FI(1))%KBAR/IREFINE
               K2 = K1
            CASE (-3)
               K1 = 1
               K2 = K1
         END SELECT
   
         !!!
         !!! Set IJKW_CO(10:15, IW_CO)
         !!!
         IJKW_CO (10, IW_CO) = I1
         IJKW_CO (11, IW_CO) = J1
         IJKW_CO (12, IW_CO) = K1
         IJKW_CO (13, IW_CO) = I2
         IJKW_CO (14, IW_CO) = J2
         IJKW_CO (15, IW_CO) = K2
   
      ENDIF
   
END SELECT

IW_CO = IW_CO + 1

END SUBROUTINE SCARC_SETUP_FACE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Set analogues to I_MIN, IMAX, J_MIN, J_MAX, K_MIN, K_MAX on coarse level (only GMG!)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SETUP_COARSE_DIMENSIONS(IREFINE, NM, NL)
INTEGER, INTENT(IN) :: IREFINE, NM, NL
INTEGER :: NOM, IMIN, IMAX, JMIN, JMAX, KMIN, KMAX, IERR, IW, NW
LOGICAL :: FOUND
TYPE (SCARC_TYPE) , POINTER :: S
TYPE (OSCARC_TYPE), POINTER :: OS

IERR = 0

S => SCARC(NM)

OTHER_MESH_LOOP: DO NOM = 1, NMESHES

   IF (NOM == NM) CYCLE OTHER_MESH_LOOP

   !!! let OS point to SCARC structure on neighboring mesh
   OS => S%OSCARC(NOM)
    
   IMIN=0
   IMAX=MESHES(NOM)%IBAR/IREFINE+1

   SELECT CASE (TYPE_DIMENSION)
      CASE (NSCARC_DIMENSION_TWO)
         JMIN=0
         JMAX=2
      CASE (NSCARC_DIMENSION_THREE)
         JMIN=0
         JMAX=MESHES(NOM)%JBAR/IREFINE+1
   END SELECT

   KMIN=0
   KMAX=MESHES(NOM)%KBAR/IREFINE+1
   
   OS%NIC_S = 0
   FOUND = .FALSE.

   SELECT CASE (TYPE_SYSTEM)
      CASE (NSCARC_SYSTEM_BANDED)
         NW = S%BANDED(NL)%NW
      CASE (NSCARC_SYSTEM_COMPACT)
         NW = S%COMPACT(NL)%NW
   END SELECT

   SEARCH_LOOP: DO IW=1,NW
   
      ! neighborship structure already known from finest level
      IF (MESHES(NM)%OMESH(NM)%IJKW(9,IW)/=NOM) CYCLE SEARCH_LOOP !****CHECK
      OS%NIC_S = OS%NIC_S + 1
      FOUND = .TRUE.
   
      SELECT CASE (MESHES(NM)%OMESH(NM)%IJKW(4,IW))!****CHECK
         CASE ( 1)
            IMIN=MAX(IMIN,MESHES(NM)%OMESH(NM)%IJKW(10,IW)-1)!****CHECK
         CASE (-1) 
            IMAX=MIN(IMAX,MESHES(NM)%OMESH(NM)%IJKW(13,IW))!****CHECK
         CASE ( 2) 
            JMIN=MAX(JMIN,MESHES(NM)%OMESH(NM)%IJKW(11,IW)-1)!****CHECK
         CASE (-2) 
            JMAX=MIN(JMAX,MESHES(NM)%OMESH(NM)%IJKW(14,IW))!****CHECK
         CASE ( 3) 
            KMIN=MAX(KMIN,MESHES(NM)%OMESH(NM)%IJKW(12,IW)-1)!****CHECK
         CASE (-3)
            KMAX=MIN(KMAX,MESHES(NM)%OMESH(NM)%IJKW(15,IW))!****CHECK
      END SELECT
   ENDDO SEARCH_LOOP
   
   IF (.NOT.FOUND) CYCLE OTHER_MESH_LOOP

   N_EXCHANGES = N_EXCHANGES+1

   OS%I_MIN_R = IMIN
   OS%I_MAX_R = IMAX
   OS%J_MIN_R = JMIN
   OS%J_MAX_R = JMAX
   OS%K_MIN_R = KMIN
   OS%K_MAX_R = KMAX
   
ENDDO OTHER_MESH_LOOP

END SUBROUTINE SCARC_SETUP_COARSE_DIMENSIONS

 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Allocate several global structures for data exchange 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SETUP_GLOBALS
INTEGER :: NM, NL, IERR
INTEGER :: IREFINE

IERR = 0

!!!-------------------------------------------------------------------------------------------------------
!!! Allocate arrays which are used for (global) communciations
!!!-------------------------------------------------------------------------------------------------------
ALLOCATE(NC_GLOBAL(NLEVEL_MIN:NLEVEL_MAX), STAT=IERR)
CALL CHKMEMERR ('SCARC_SETUP', 'NC_GLOBAL', IERR)
NC_GLOBAL = 0

ALLOCATE(NC_LOCAL(NMESHES), STAT=IERR)
CALL CHKMEMERR ('SCARC_SETUP', 'NC_LOCAL', IERR)
NC_LOCAL = 0

ALLOCATE(SP_LOCAL(NMESHES), STAT=IERR)
CALL CHKMEMERR ('SCARC_SETUP_GLOBAL', 'SP_LOCAL', IERR)
SP_LOCAL = 0.0_EB

!!!-------------------------------------------------------------------------------------------------------
!!! Compute global number of cells for all levels
!!!-------------------------------------------------------------------------------------------------------
IREFINE=0
LEVEL_LOOP: DO NL = NLEVEL_MIN, NLEVEL_MAX

   MESHES_LOOP: DO NM = NMESHES_MIN, NMESHES_MAX
      SELECT CASE (TYPE_SYSTEM)
         CASE (NSCARC_SYSTEM_BANDED)
            NC_LOCAL(NM)=SCARC(NM)%BANDED(NL)%NC
         CASE (NSCARC_SYSTEM_COMPACT)
            NC_LOCAL(NM)=SCARC(NM)%COMPACT(NL)%NC
      END SELECT
   ENDDO MESHES_LOOP

   !!! Determine global number of cells for all levels 
   IF (NMESHES>1) THEN
      IF (USE_MPI) THEN
         CALL MPI_ALLREDUCE(NC_LOCAL(MYID+1),NC_GLOBAL(NL),1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,IERR)
      ELSE
         NC_GLOBAL(NL)=0
         DO NM=1,NMESHES
            NC_GLOBAL(NL) = NC_GLOBAL(NL) + NC_LOCAL(NM)
         ENDDO
      ENDIF
   ELSE
      NC_GLOBAL(NL) = NC_LOCAL(1)
   ENDIF

ENDDO LEVEL_LOOP

END SUBROUTINE SCARC_SETUP_GLOBALS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Allocate several global structures for data exchange 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SETUP_NEIGHBORS
INTEGER :: NM, NOM, NL, IERR
TYPE (OSCARC_TYPE), POINTER :: OS, OSO


IF (NMESHES < 0) THEN

   LEVEL_LOOP: DO NL = NLEVEL_MIN+1, NLEVEL_MAX
   
      MESHES_LOOP: DO NM=NMESHES_MIN,NMESHES_MAX
         OMESHES_LOOP: DO NOM=1,NMESHES
   
            OS  => SCARC(NM)%OSCARC(NOM)
            OSO => SCARC(NOM)%OSCARC(NM)
   
            IF (USE_MPI) THEN
               CALL MPI_SEND(OS%I_MIN_R,1,MPI_INTEGER,PROCESS(NOM),1,MPI_COMM_WORLD,IERR)
               CALL MPI_SEND(OS%I_MAX_R,1,MPI_INTEGER,PROCESS(NOM),1,MPI_COMM_WORLD,IERR)
               CALL MPI_SEND(OS%J_MIN_R,1,MPI_INTEGER,PROCESS(NOM),1,MPI_COMM_WORLD,IERR)
               CALL MPI_SEND(OS%J_MAX_R,1,MPI_INTEGER,PROCESS(NOM),1,MPI_COMM_WORLD,IERR)
               CALL MPI_SEND(OS%K_MIN_R,1,MPI_INTEGER,PROCESS(NOM),1,MPI_COMM_WORLD,IERR)
               CALL MPI_SEND(OS%K_MAX_R,1,MPI_INTEGER,PROCESS(NOM),1,MPI_COMM_WORLD,IERR)
               CALL MPI_SEND(OS%NIC_S,  1,MPI_INTEGER,PROCESS(NOM),1,MPI_COMM_WORLD,IERR)
            ELSE
               OSO%I_MIN_S = OS%I_MIN_R
               OSO%I_MAX_S = OS%I_MAX_R
               OSO%J_MIN_S = OS%J_MIN_R
               OSO%J_MAX_S = OS%J_MAX_R
               OSO%K_MIN_S = OS%K_MIN_R
               OSO%K_MAX_S = OS%K_MAX_R
               OSO%NIC_R   = OS%NIC_S
            ENDIF

         ENDDO OMESHES_LOOP
      ENDDO MESHES_LOOP
   
      MESHES_LOOP2: DO NM=1,NMESHES
         OMESHES_LOOP2: DO NOM=NMESHES_MIN,NMESHES_MAX
   
            OSO => SCARC(NOM)%OSCARC(NM)
   
            IF (USE_MPI) THEN
               CALL MPI_RECV(OSO%I_MIN_S,1,MPI_INTEGER,PROCESS(NM),1,MPI_COMM_WORLD,STATUS2_SCARC,IERR)
               CALL MPI_RECV(OSO%I_MAX_S,1,MPI_INTEGER,PROCESS(NM),1,MPI_COMM_WORLD,STATUS2_SCARC,IERR)
               CALL MPI_RECV(OSO%J_MIN_S,1,MPI_INTEGER,PROCESS(NM),1,MPI_COMM_WORLD,STATUS2_SCARC,IERR)
               CALL MPI_RECV(OSO%J_MAX_S,1,MPI_INTEGER,PROCESS(NM),1,MPI_COMM_WORLD,STATUS2_SCARC,IERR)
               CALL MPI_RECV(OSO%K_MIN_S,1,MPI_INTEGER,PROCESS(NM),1,MPI_COMM_WORLD,STATUS2_SCARC,IERR)
               CALL MPI_RECV(OSO%K_MAX_S,1,MPI_INTEGER,PROCESS(NM),1,MPI_COMM_WORLD,STATUS2_SCARC,IERR)
               CALL MPI_RECV(OSO%NIC_R,  1,MPI_INTEGER,PROCESS(NM),1,MPI_COMM_WORLD,STATUS2_SCARC,IERR)
            ENDIF
   
         ENDDO OMESHES_LOOP2
      ENDDO MESHES_LOOP2
   
   ENDDO LEVEL_LOOP

ENDIF

END SUBROUTINE SCARC_SETUP_NEIGHBORS

 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Setup system of equation:
!!! Define matrix stencils and initialize matrices and boundary conditions on all needed levels
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SETUP_SYSTEM
INTEGER :: NM, NL, IERR
TYPE (SCARC_TYPE), POINTER :: S

IERR = 0
 
MESHES_LOOP: DO NM = NMESHES_MIN, NMESHES_MAX

   S  => SCARC(NM)

   SELECT_SOLVER: SELECT CASE (TYPE_METHOD)
   
      !!!-------------------------------------------------------------------------------------------------
      !!! Krylov method (CG/BICG) as main solver, different preconditioners possible
      !!!-------------------------------------------------------------------------------------------------
      CASE (NSCARC_METHOD_KRYLOV)
   
         SELECT_PRECON: SELECT CASE (TYPE_PRECON)
   
            !!!-------------------------------------------------------------------------------------------
            !!! in case of multigrid as preconditioner:
            !!!-------------------------------------------------------------------------------------------
            CASE (NSCARC_PRECON_MULTIGRID)
   
               SELECT_PRECON_MG: SELECT CASE (TYPE_MULTIGRID)
   
                  !!! geometric multigrid:
                  !!!    -  use banded storage technique on all levels unless otherwise specified
                  !!!    -  assemble standard n-point-matrix hierarchy on all levels (finest and all coarser)
                  CASE (NSCARC_MULTIGRID_GEOMETRIC)
                     
                     DO NL = NLEVEL_MIN, NLEVEL_MAX
                        CALL SCARC_SETUP_MATRIX  (NM, NL)
                        CALL SCARC_SETUP_BOUNDARY(NM, NL)
                     ENDDO 
   
                  !!! algebraic multigrid:
                  !!!    -  use compact storage technique on all levels (no other choise possible!)
                  !!!    -  assemble standard n-point-matrix only on finest level 
                  !!!    -  construct all coarser levels by requested coarsening strategy
                  CASE (NSCARC_MULTIGRID_ALGEBRAIC)
   
                     CALL SCARC_SETUP_MATRIX  (NM, NLEVEL_MIN)
                     CALL SCARC_SETUP_BOUNDARY(NM, NLEVEL_MIN)
   
               END SELECT SELECT_PRECON_MG
   
            !!!-------------------------------------------------------------------------------------------
            !!! in case of one-level preconditioners (JACOBI/SSOR/GSTRIX/FFT)
            !!!    -  use banded storage technique on finest level unless otherwise specified 
            !!!    -  assemble standard n-point-matrix on finest level 
            !!!-------------------------------------------------------------------------------------------
            CASE DEFAULT
   
               CALL SCARC_SETUP_MATRIX  (NM, NLEVEL_MIN)
               CALL SCARC_SETUP_BOUNDARY(NM, NLEVEL_MIN)
   
         END SELECT SELECT_PRECON
   
      !!!-------------------------------------------------------------------------------------------------
      !!! Multigrid as main solver
      !!!-------------------------------------------------------------------------------------------------
      CASE (NSCARC_METHOD_MULTIGRID)
   
         SELECT_MG: SELECT CASE (TYPE_MULTIGRID)
   
            !!!-------------------------------------------------------------------------------------------
            !!! geometric multigrid:
            !!!    -  use banded storage technique on all levels unless otherwise specified 
            !!!    -  assemble standard n-point-matrix hierarchy on all levels 
            !!!-------------------------------------------------------------------------------------------
            CASE (NSCARC_MULTIGRID_GEOMETRIC)
   
               DO NL = NLEVEL_MIN, NLEVEL_MAX
                  CALL SCARC_SETUP_MATRIX  (NM, NL)
                  CALL SCARC_SETUP_BOUNDARY(NM, NL)
               ENDDO 

   
            !!!-------------------------------------------------------------------------------------------
            !!! algebraic multigrid:
            !!!    -  use compact storage technique (no other choice possible!)
            !!!    -  assemble standard n-point-matrix only on finest level
            !!!    -  construct all coarser levels later by requested coarsening strategy
            !!!-------------------------------------------------------------------------------------------
            CASE (NSCARC_MULTIGRID_ALGEBRAIC)
   
               CALL SCARC_SETUP_MATRIX  (NM, NLEVEL_MIN)
               CALL SCARC_SETUP_BOUNDARY(NM, NLEVEL_MIN)

         END SELECT SELECT_MG
   
   END SELECT SELECT_SOLVER
   
ENDDO MESHES_LOOP

DO NL=NLEVEL_MIN, NLEVEL_MIN
   CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_MATRIX , NL, 'SETUP_SYSTEM0', 'MATRIX')
   CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_IJKW   , NL, 'SETUP_SYSTEM0', 'IJKW')
ENDDO

!!!-------------------------------------------------------------------------------------------------------
!!! Exchange matrix entries along internal boundaries:
!!!  - in case of GMG as solver or preconditioner: matrices of all levels must be exchanged
!!!  - in all other cases: only matrix of finest level must be exchanged
!!!-------------------------------------------------------------------------------------------------------
IF (NMESHES>1) THEN

   SELECT CASE (TYPE_MULTIGRID)
      CASE (NSCARC_MULTIGRID_GEOMETRIC)
         TYPE_MATRIX = NSCARC_MATRIX_SUBDIAG
         DO NL = NLEVEL_MIN, NLEVEL_MAX
            CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_MATRIX, NL)
         ENDDO
      CASE DEFAULT
         TYPE_MATRIX = NSCARC_MATRIX_SUBDIAG
         CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_MATRIX, NLEVEL_MIN)
         !TYPE_MATRIX = NSCARC_MATRIX_SYSTEM
         !CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_MATRIX, NLEVEL_MIN)
   END SELECT

ENDIF

DO NL=NLEVEL_MIN, NLEVEL_MIN
   CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_MATRIX , NL, 'SETUP_SYSTEM', 'MATRIX')
   CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_IJKW   , NL, 'SETUP_SYSTEM', 'IJKW')
ENDDO


!!!-------------------------------------------------------------------------------------------------------
!!! If multigrid is used, setup global coarse grid matrix in case of a direct coarse grid solver
!!!-------------------------------------------------------------------------------------------------------
IF ((TYPE_METHOD == NSCARC_METHOD_MULTIGRID .OR. TYPE_PRECON == NSCARC_PRECON_MULTIGRID) .AND. &
     TYPE_COARSE == NSCARC_COARSE_DIRECT) THEN
   CALL SCARC_SETUP_COARSE_MATRIX(NLEVEL_MAX)
ENDIF

DO NL=NLEVEL_MIN, NLEVEL_MIN
   CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_MATRIX , NL, 'SETUP_SYSTEM', 'MATRIX000')
ENDDO


CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_IJKW   , NLEVEL_MIN, 'SETUP_SYSTEM', 'IJKW')
CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_MATRIX , NLEVEL_MIN, 'SETUP_SYSTEM', 'MATRIX')
CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_MATRIXE, NLEVEL_MIN, 'SETUP_SYSTEM', 'MATRIXE')
CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_BCINDEX, NLEVEL_MIN, 'SETUP_SYSTEM', 'PRESSURE_BC_INDEX')
CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_ACELL  , NLEVEL_MIN, 'SETUP_SYSTEM', 'ADJACENT_CELL')
CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_GCELL  , NLEVEL_MIN, 'SETUP_SYSTEM', 'GHOST_CELL')
END SUBROUTINE SCARC_SETUP_SYSTEM


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Allocate matrix for the usual 5-point-stencil (2D) or 7-point-stencil (3D)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SETUP_COARSE_MATRIX (NL)
INTEGER, INTENT(IN) :: NL
TYPE (SCARC_TYPE) , POINTER :: SM
TYPE (SCARC_BANDED_TYPE) , POINTER :: SB, SBO
TYPE (SCARC_COMPACT_TYPE), POINTER :: SC, SCO
INTEGER , ALLOCATABLE, DIMENSION(:) :: BUF_INT, AC_COL, AC_ROW, IJKWC
REAL(EB), ALLOCATABLE, DIMENSION(:) :: BUF_REAL, AC
REAL(EB), ALLOCATABLE, DIMENSION(:,:,:) :: ACG
INTEGER , ALLOCATABLE, DIMENSION(:,:,:) :: IJKWG
!INTEGER, ALLOCATABLE, DIMENSION(:)  :: COUNTS1, DISPLS1
!INTEGER, ALLOCATABLE, DIMENSION(:)  :: COUNTS2, DISPLS2
INTEGER :: NM, NOM, NA_SUM, NC_SUM, NW_SUM, DISP, CNT, NCPL, NX_OTHER, NY_OTHER
INTEGER :: NM0, NX0, NY0, NZ0, NW0, NC0
INTEGER :: IW, IC, IC0, IC1, IC2, ICO, IOFFSET, IERR, IOR0, ICOL, IA, ICPL, I
INTEGER :: II0, JJ0, KK0, II1, JJ1, KK1, II2, JJ2, KK2
REAL(EB) :: VAL

IERR=0

!!! For the moment allocate coarse grid matrix on mesh 1
!!! Later a strategy for chosing the mesh with minimal cell load will be regarded ...
NMASTER = 1
SM => SCARC(NMASTER)

ALLOCATE(NA_COARSE(NMESHES), STAT=IERR)
CALL ChkMemErr('SCARC_SETUP_COARSE_MATRIX','NA_COARSE',IERR)
ALLOCATE(NC_COARSE(NMESHES), STAT=IERR)
CALL ChkMemErr('SCARC_SETUP_COARSE_MATRIX','NC_COARSE',IERR)
ALLOCATE(NW_COARSE(NMESHES), STAT=IERR)
CALL ChkMemErr('SCARC_SETUP_COARSE_MATRIX','NW_COARSE',IERR)

SELECT CASE(TYPE_SYSTEM)

   CASE (NSCARC_SYSTEM_BANDED) 

      ALLOCATE(NX_COARSE(NMESHES), STAT=IERR)
      CALL ChkMemErr('SCARC_SETUP_COARSE_MATRIX','NX_COARSE',IERR)
      ALLOCATE(NY_COARSE(NMESHES), STAT=IERR)
      CALL ChkMemErr('SCARC_SETUP_COARSE_MATRIX','NY_COARSE',IERR)
      ALLOCATE(NZ_COARSE(NMESHES), STAT=IERR)
      CALL ChkMemErr('SCARC_SETUP_COARSE_MATRIX','NZ_COARSE',IERR)
      ALLOCATE(DXI2_COARSE(NMESHES), STAT=IERR)
      CALL ChkMemErr('SCARC_SETUP_COARSE_MATRIX','DXI2_COARSE',IERR)
      ALLOCATE(DYI2_COARSE(NMESHES), STAT=IERR)
      CALL ChkMemErr('SCARC_SETUP_COARSE_MATRIX','DYI2_COARSE',IERR)
      ALLOCATE(DZI2_COARSE(NMESHES), STAT=IERR)
      CALL ChkMemErr('SCARC_SETUP_COARSE_MATRIX','DZI2_COARSE',IERR)

      DO NM = NMESHES_MIN, NMESHES_MAX
         SB => SCARC(NM)%BANDED(NL)
         NA_COARSE(NM)   = SB%NA
         NC_COARSE(NM)   = SB%NC
         NW_COARSE(NM)   = SB%NW
         NX_COARSE(NM)   = SB%NX
         NY_COARSE(NM)   = SB%NY
         NZ_COARSE(NM)   = SB%NZ
         DXI2_COARSE(NM) = SB%DXI2
         DYI2_COARSE(NM) = SB%DYI2
         DZI2_COARSE(NM) = SB%DZI2
      ENDDO

   CASE (NSCARC_SYSTEM_COMPACT) 

      DO NM = NMESHES_MIN, NMESHES_MAX
         SC => SCARC(NM)%COMPACT(NL)
         NA_COARSE(NM) = SC%NA
         NC_COARSE(NM) = SC%NC
         NW_COARSE(NM) = SC%NW
      ENDDO

END SELECT

!CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_GRID, NL)

IF (USE_MPI) THEN

   ALLOCATE(BUF_INT(NMESHES), STAT=IERR)
   CALL ChkMemErr('SCARC_SETUP_COARSE_MATRIX','BUF_INT',IERR)

   ALLOCATE(BUF_REAL(NMESHES), STAT=IERR)
   CALL ChkMemErr('SCARC_SETUP_COARSE_MATRIX','BUF_REAL',IERR)

   SELECT CASE(TYPE_SYSTEM)

      CASE (NSCARC_SYSTEM_BANDED) 

         BUF_INT = NA_COARSE
         CALL MPI_ALLGATHERV(BUF_INT(DISPLS_SCARC(MYID)+1),COUNTS_SCARC(MYID), MPI_INTEGER, &
                             NA_COARSE, COUNTS_SCARC, DISPLS_SCARC, MPI_INTEGER, MPI_COMM_WORLD, IERR)
         BUF_INT = NC_COARSE
         CALL MPI_ALLGATHERV(BUF_INT(DISPLS_SCARC(MYID)+1),COUNTS_SCARC(MYID), MPI_INTEGER, &
                             NC_COARSE, COUNTS_SCARC, DISPLS_SCARC, MPI_INTEGER, MPI_COMM_WORLD, IERR)
         BUF_INT = NW_COARSE
         CALL MPI_ALLGATHERV(BUF_INT(DISPLS_SCARC(MYID)+1),COUNTS_SCARC(MYID), MPI_INTEGER, &
                             NW_COARSE, COUNTS_SCARC, DISPLS_SCARC, MPI_INTEGER, MPI_COMM_WORLD, IERR)
         BUF_INT = NX_COARSE
         CALL MPI_ALLGATHERV(BUF_INT(DISPLS_SCARC(MYID)+1),COUNTS_SCARC(MYID), MPI_INTEGER, &
                             NX_COARSE, COUNTS_SCARC, DISPLS_SCARC, MPI_INTEGER, MPI_COMM_WORLD, IERR)
         BUF_INT = NY_COARSE
         CALL MPI_ALLGATHERV(BUF_INT(DISPLS_SCARC(MYID)+1),COUNTS_SCARC(MYID), MPI_INTEGER, &
                             NY_COARSE, COUNTS_SCARC, DISPLS_SCARC, MPI_INTEGER, MPI_COMM_WORLD, IERR)
         BUF_INT = NZ_COARSE
         CALL MPI_ALLGATHERV(BUF_INT(DISPLS_SCARC(MYID)+1),COUNTS_SCARC(MYID), MPI_INTEGER, &
                             NZ_COARSE, COUNTS_SCARC, DISPLS_SCARC, MPI_INTEGER, MPI_COMM_WORLD, IERR)
         BUF_REAL = DXI2_COARSE
         CALL MPI_ALLGATHERV(BUF_REAL(DISPLS_SCARC(MYID)+1),COUNTS_SCARC(MYID), MPI_DOUBLE_PRECISION, &
                             DXI2_COARSE, COUNTS_SCARC, DISPLS_SCARC, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, IERR)
         BUF_REAL = DYI2_COARSE
         CALL MPI_ALLGATHERV(BUF_REAL(DISPLS_SCARC(MYID)+1),COUNTS_SCARC(MYID), MPI_DOUBLE_PRECISION, &
                             DYI2_COARSE, COUNTS_SCARC, DISPLS_SCARC, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, IERR)
         BUF_REAL = DZI2_COARSE
         CALL MPI_ALLGATHERV(BUF_REAL(DISPLS_SCARC(MYID)+1),COUNTS_SCARC(MYID), MPI_DOUBLE_PRECISION, &
                             DZI2_COARSE, COUNTS_SCARC, DISPLS_SCARC, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, IERR)
         DEALLOCATE(BUF_REAL)

      CASE (NSCARC_SYSTEM_COMPACT) 
         BUF_INT = NA_COARSE
         CALL MPI_ALLGATHERV(BUF_INT(DISPLS_SCARC(MYID)+1),COUNTS_SCARC(MYID), MPI_INTEGER, &
                             NA_COARSE, COUNTS_SCARC, DISPLS_SCARC, MPI_INTEGER, MPI_COMM_WORLD, IERR)
         BUF_INT = NC_COARSE
         CALL MPI_ALLGATHERV(BUF_INT(DISPLS_SCARC(MYID)+1),COUNTS_SCARC(MYID), MPI_INTEGER, &
                             NC_COARSE, COUNTS_SCARC, DISPLS_SCARC, MPI_INTEGER, MPI_COMM_WORLD, IERR)
         BUF_INT = NW_COARSE
         CALL MPI_ALLGATHERV(BUF_INT(DISPLS_SCARC(MYID)+1),COUNTS_SCARC(MYID), MPI_INTEGER, &
                             NW_COARSE, COUNTS_SCARC, DISPLS_SCARC, MPI_INTEGER, MPI_COMM_WORLD, IERR)
   END SELECT

   DEALLOCATE(BUF_INT)

ENDIF

IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) WRITE(SCARC_LU,*) 'NA_COARSE:',(NA_COARSE(1:NMESHES))
IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) WRITE(SCARC_LU,*) 'NC_COARSE:',(NC_COARSE(1:NMESHES))
IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) WRITE(SCARC_LU,*) 'NW_COARSE:',(NW_COARSE(1:NMESHES))

IF (USE_MPI) THEN

   ALLOCATE(SM%COUNTS1(0:NUMPROCS-1), STAT=IERR)
   CALL ChkMemErr('SCARC_SETUP_COARSE_MATRIX','COUNTS1',IERR)
   ALLOCATE(SM%COUNTS2(0:NUMPROCS-1), STAT=IERR)
   CALL ChkMemErr('SCARC_SETUP_COARSE_MATRIX','COUNTS2',IERR)
   ALLOCATE(SM%DISPLS1(0:NUMPROCS-1), STAT=IERR)
   CALL ChkMemErr('SCARC_SETUP_COARSE_MATRIX','DISPLS1',IERR)
   ALLOCATE(SM%DISPLS2(0:NUMPROCS-1), STAT=IERR)
   CALL ChkMemErr('SCARC_SETUP_COARSE_MATRIX','DISPLS2',IERR)

   SELECT CASE (TYPE_SYSTEM)

      CASE (NSCARC_SYSTEM_BANDED) 

         NA_SUM = 0
         NC_SUM = 0
         NW_SUM = 0
         DO NM = 1, NMESHES

            SM%COUNTS1(NM-1) = NA_COARSE(NM)
            SM%DISPLS1(NM-1) = NA_SUM
            SM%COUNTS2(NM-1) = NW_COARSE(NM)*15
            SM%DISPLS2(NM-1) = NW_SUM*15
   
            NA_SUM = NA_SUM + NA_COARSE(NM)
            NC_SUM = NC_SUM + NC_COARSE(NM) 
            NW_SUM = NW_SUM + NW_COARSE(NM)

         ENDDO

IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) THEN
   WRITE(SCARC_LU,*) 'NA_SUM ', NA_SUM
   WRITE(SCARC_LU,*) 'NC_SUM ', NC_SUM
   WRITE(SCARC_LU,*) 'NW_SUM ', NW_SUM
   WRITE(SCARC_LU,*) 'COUNTS1 ', (SM%COUNTS1(0:NUMPROCS-1))
   WRITE(SCARC_LU,*) 'DISPLS1 ', (SM%DISPLS1(0:NUMPROCS-1))
   WRITE(SCARC_LU,*) 'COUNTS2 ', (SM%COUNTS2(0:NUMPROCS-1))
   WRITE(SCARC_LU,*) 'DISPLS2 ', (SM%DISPLS2(0:NUMPROCS-1))
ENDIF

         NM = MYID+1 

         !!! -----------------------------------------------------------------------------
         !!! Sample system matrix
         !!! -----------------------------------------------------------------------------
         ALLOCATE(AC(NA_SUM), STAT=IERR)
         CALL ChkMemErr('SCARC_SETUP_COARSE_MATRIX','AC',IERR)
         AC = 0.0_EB

         ALLOCATE(BUF_REAL(NA_SUM), STAT=IERR)
         CALL ChkMemErr('SCARC_SETUP_COARSE_MATRIX','BUF_REAL',IERR)
         BUF_REAL = 0.0_EB

         IF (TYPE_DIMENSION == NSCARC_DIMENSION_TWO) THEN
            NCPL = 5
         ELSE
            NCPL = 7
         ENDIF

         ALLOCATE(ACG(MAXVAL(NC_COARSE), NCPL, NMESHES), STAT=IERR)
         CALL ChkMemErr('SCARC_SETUP_COARSE_MATRIX','ACG',IERR)
         ACG = 0.0_EB

         DISP = SM%DISPLS1(MYID)+1
         DO IC = 1, NC_COARSE(NM)
            DO ICPL = 1, NCPL  
               BUF_REAL(DISP) = SCARC(NM)%BANDED(NL)%A(IC, ICPL)
               DISP = DISP + 1
            ENDDO
         ENDDO
         DISP = SM%DISPLS1(MYID)+1
         CALL MPI_ALLGATHERV(BUF_REAL(DISP),SM%COUNTS1(MYID),MPI_DOUBLE_PRECISION, &
                             AC,SM%COUNTS1,SM%DISPLS1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,IERR)

IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) WRITE(SCARC_LU,*) 'BUF_REAL (A): after'
IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) WRITE(SCARC_LU,'(7f12.5)') (AC(1:NA_SUM))


         !!! -----------------------------------------------------------------------------
         !!! Sample ijkw
         !!! -----------------------------------------------------------------------------
         ALLOCATE(IJKWC(NW_SUM*15), STAT=IERR)
         CALL ChkMemErr('SCARC_SETUP_COARSE_MATRIX','IJKWC',IERR)
         IJKWC = 0

         ALLOCATE(BUF_INT(NW_SUM*15), STAT=IERR)
         CALL ChkMemErr('SCARC_SETUP_COARSE_MATRIX','BUF_INT',IERR)
         BUF_INT = 0

         ALLOCATE(IJKWG(15,MAXVAL(NW_COARSE), NMESHES), STAT=IERR)
         CALL ChkMemErr('SCARC_SETUP_COARSE_MATRIX','IJKWG',IERR)
         IJKWG = 0

         DISP = SM%DISPLS2(MYID)+1
         DO IW = 1, NW_COARSE(NM)
            DO I = 1, 15  
               BUF_INT(DISP) = SCARC(NM)%BANDED(NL)%IJKW(I, IW)
               DISP = DISP + 1
            ENDDO
         ENDDO
         DISP = SM%DISPLS2(MYID)+1
         CALL MPI_ALLGATHERV(BUF_INT(DISP),SM%COUNTS2(MYID),MPI_INTEGER, &
                             IJKWC,SM%COUNTS2,SM%DISPLS2,MPI_INTEGER,MPI_COMM_WORLD,IERR)

IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) WRITE(SCARC_LU,*) 'BUF_INT (IJKW): after'
IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) WRITE(SCARC_LU,'(15i4)') (IJKWC(1:NW_SUM*15))

         DISP = 1
         DO NM = 1, NMESHES
            DO IC = 1, NC_COARSE(NM)
               DO ICPL = 1, NCPL  
                  ACG(IC, ICPL, NM) = AC(DISP)  
                  DISP = DISP + 1
               ENDDO
!               WRITE(SCARC_LU,'(a,i3,a,i3,a,i3,a,7f9.2)') 'ACG(',IC,',',ICPL,',',NM,')=',(ACG(IC,ICPL,NM),ICPL=1,NCPL)
            ENDDO
         ENDDO

         DISP = 1
         DO NM = 1, NMESHES
            DO IW = 1, NW_COARSE(NM)
               DO I = 1, 15
                  IJKWG(I, IW, NM) = IJKWC(DISP)  
                  DISP = DISP + 1
               ENDDO
!               WRITE(SCARC_LU,'(a,i3,a,i3,a,i3,a,15i5)') 'IJKWG(',I,',',IW,',',NM,')=',(IJKWG(I,IW,NM),I=1,15)
            ENDDO
         ENDDO

      CASE (NSCARC_SYSTEM_COMPACT) 

         NA_SUM = 0
         NC_SUM = 0
         DO NM = 1, NMESHES

            SM%COUNTS1(NM-1) = NA_COARSE(NM)+1
            SM%DISPLS1(NM-1) = NA_SUM
            SM%COUNTS2(NM-1) = NC_COARSE(NM)+1
            SM%DISPLS2(NM-1) = NC_SUM
   
            NA_SUM = NA_SUM + NA_COARSE(NM) + 1
            NC_SUM = NC_SUM + NC_COARSE(NM) + 1

         ENDDO

IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) THEN
   WRITE(SCARC_LU,*) 'NA_SUM ', NA_SUM
   WRITE(SCARC_LU,*) 'NC_SUM ', NC_SUM
   WRITE(SCARC_LU,*) 'COUNTS1 ', (SM%COUNTS1(0:NUMPROCS-1))
   WRITE(SCARC_LU,*) 'DISPLS1 ', (SM%COUNTS1(0:NUMPROCS-1))
   WRITE(SCARC_LU,*) 'COUNTS2 ', (SM%COUNTS1(0:NUMPROCS-1))
   WRITE(SCARC_LU,*) 'DISPLS2 ', (SM%COUNTS1(0:NUMPROCS-1))
ENDIF

         ALLOCATE(AC(NA_SUM), STAT=IERR)
         CALL ChkMemErr('SCARC_SETUP_COARSE_MATRIX','AC',IERR)
         AC = 0.0_EB

         ALLOCATE(AC_COL(NA_SUM), STAT=IERR)
         CALL ChkMemErr('SCARC_SETUP_COARSE_MATRIX','AC_COL',IERR)
         AC_COL = 0

         ALLOCATE(AC_ROW(NC_SUM), STAT=IERR)
         CALL ChkMemErr('SCARC_SETUP_COARSE_MATRIX','AC_ROW',IERR)
         AC_ROW = 0

         ALLOCATE(BUF_REAL(NA_SUM), STAT=IERR)
         CALL ChkMemErr('SCARC_SETUP_COARSE_MATRIX','BUF_REAL',IERR)
         BUF_REAL = 0

         ALLOCATE(BUF_INT(NA_SUM), STAT=IERR)
         CALL ChkMemErr('SCARC_SETUP_COARSE_MATRIX','BUF_INT',IERR)
         BUF_INT = 0

         DISP = SM%DISPLS2(MYID)+1
         CNT  = SM%COUNTS2(MYID)+1

         WRITE(*,*) 'NOCH NICHT FERTIG !!!'
         !BUF_INT(DISP:DISP+CNT) = SCARC(MYID+1)%BANDED(NL)%A_ROW(1:NA_COARSE(MYID+1)+1)          
         CALL MPI_GATHERV(BUF_INT(DISP),SM%COUNTS2(MYID),MPI_INTEGER, &
                          AC_ROW,SM%COUNTS2,SM%DISPLS2,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)

IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) WRITE(SCARC_LU,*) 'BUF_INT (A_ROW):',(BUF_INT(1:NC_SUM))

         DISP = SM%DISPLS1(MYID)+1
         CNT  = SM%COUNTS1(MYID)+1

         !BUF_INT(DISP:DISP+CNT) = SCARC(MYID+1)%BANDED(NL)%A_COL(1:NA_COARSE(MYID+1)+1)          
         CALL MPI_GATHERV(BUF_INT(DISP),SM%COUNTS1(MYID),MPI_INTEGER, &
                          AC_ROW,SM%COUNTS1,SM%DISPLS1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)

IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) WRITE(SCARC_LU,*) 'BUF_INT (A_COL):',(BUF_INT(1:NA_SUM))

         BUF_REAL(DISP:DISP+CNT) = MYID+1
         !BUF_REAL(DISP:DISP+CNT) = SCARC(MYID+1)%BANDED(NL)%A(1:NA_COARSE(MYID+1)+1)          
         CALL MPI_GATHERV(BUF_REAL(DISP),SM%COUNTS1(MYID),MPI_DOUBLE_PRECISION, &
                          AC,SM%COUNTS1,SM%DISPLS1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)

IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) WRITE(SCARC_LU,*) 'BUF_REAL (A):'
IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) WRITE(SCARC_LU,'(5f12.5)') (BUF_INT(1:NA_SUM))
         WRITE(*,*) 'SETUP_COARSE_MATRIX: COMPACT NOT YET IMPLEMENTED'
         STOP

   END SELECT
 
ENDIF

ALLOCATE(SM%OFFSET(NMESHES), STAT=IERR)
CALL ChkMemErr('SCARC_SETUP_COARSE_MATRIX','OFFSET',IERR)

NC_COARSE0 = 0
IOFFSET    = 0
DO NM = 1, NMESHES

   SM%OFFSET(NM) = IOFFSET

   NC_COARSE0 = NC_COARSE0 + NC_COARSE(NM)
   IOFFSET    = IOFFSET + NC_COARSE(NM) 

ENDDO

IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) WRITE(SCARC_LU,*) 'ALLOCATING A_COARSE(',NC_COARSE0,',',NC_COARSE0,')'

ALLOCATE(SM%A_COARSE(NC_COARSE0, NC_COARSE0), STAT=IERR)
CALL ChkMemErr('SCARC_SETUP_COARSE_MATRIX','A_COARSE',IERR)
SM%A_COARSE = 0.0_EB

ALLOCATE(SM%PIVOT(NC_COARSE0, NC_COARSE0), STAT=IERR)
CALL ChkMemErr('SCARC_SETUP_COARSE_MATRIX','PIVOT',IERR)
SM%PIVOT = 0

ALLOCATE(SM%X_COARSE(NC_COARSE0), STAT=IERR)
CALL ChkMemErr('SCARC_SETUP_COARSE_MATRIX','X_COARSE',IERR)
SM%X_COARSE = 0.0_EB

ALLOCATE(SM%X_BUF(NC_COARSE0), STAT=IERR)
CALL ChkMemErr('SCARC_SETUP_COARSE_MATRIX','X_BUF',IERR)
SM%X_BUF = 0.0_EB

IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) WRITE(SCARC_LU,*) 'Ready', USE_MPI
 
!!! Sample local coarse matrices and put them together to a global coarse matrix

!!! ---------------------------------------------------------------------------------------------
!!! Parallel version
!!! ---------------------------------------------------------------------------------------------
IF (USE_MPI) THEN
   
      NM0 = MYID+1
!      IF (NM0 == NMASTER) THEN
   
      SELECT CASE(TYPE_SYSTEM)
         CASE (NSCARC_SYSTEM_BANDED)
   
            !!! chose dimension
            SELECT CASE (TYPE_DIMENSION)
   
               CASE (NSCARC_DIMENSION_TWO)
   
                  SB => SCARC(NM0)%BANDED(NL)

                  DO NM = 1, NMESHES
   
                     NX0 = NX_COARSE(NM)
                     NY0 = NY_COARSE(NM)
                     NZ0 = NZ_COARSE(NM)
                     NC0 = NC_COARSE(NM)
                     NW0 = NW_COARSE(NM)
   
                     INTERNALCELL_BLOOP2D_1: DO IC = 1, NC_COARSE(NM)
   
                        IC0 = IC + SM%OFFSET(NM)
   
                        SM%A_COARSE(IC0, IC0) = ACG(IC,ID,NM)
   
                        IF (IC > 1        .AND. ACG(IC,ILX,NM)>0) SM%A_COARSE(IC0, IC0-1)     = ACG(IC,ILX,NM)
                        IF (IC < NC0      .AND. ACG(IC,IUX,NM)>0) SM%A_COARSE(IC0, IC0+1)     = ACG(IC,IUX,NM)
                        IF (IC > NX0      .AND. ACG(IC,ILZ,NM)>0) SM%A_COARSE(IC0, IC0-SB%NX) = ACG(IC,ILZ,NM)
                        IF (IC < NC0-NX0+1.AND. ACG(IC,IUZ,NM)>0) SM%A_COARSE(IC0, IC0+SB%NX) = ACG(IC,IUZ,NM)
   
                     ENDDO INTERNALCELL_BLOOP2D_1
   
                     WALLCELL_BLOOP2D_1: DO IW = 1, NW0
   
                        NOM = IJKWG(9,IW,NM)
                        IF (NOM == 0) CYCLE WALLCELL_BLOOP2D_1
   
                        IOR0 = IJKWG(4,IW,NM)
   
                        II0 = IJKWG(6 ,IW,NM) 
                        II1 = IJKWG(10,IW,NM) 
                        II2 = IJKWG(13,IW,NM) 
   
                        KK0 = IJKWG(8 ,IW,NM) 
                        KK1 = IJKWG(12,IW,NM) 
                        KK2 = IJKWG(15,IW,NM) 
   
                        IF (TYPE_DEBUG >= NSCARC_DEBUG_MEDIUM) WRITE(SCARC_LU,*) II0, II1, II2, KK0, KK1, KK2
                       
                        !!! For the moment: write error message in case of different resolutions
                        IF (II2-II1.NE.0.OR.KK2-KK1.NE.0) THEN
                           WRITE(*,*) 'Error in SCARC_SETUP_COARSE_MATRIX: Different resolutions on meshes ',NM, NOM
                           STOP
                        ENDIF
   
                        IC0 = (KK0-1)*NC_COARSE(NM)  + II0  
                        IC  = IC0 + SM%OFFSET(NM)      
   
                        IC1 = (KK1-1)*NX_COARSE(NOM) + II1  
                        ICO = IC1 + SM%OFFSET(NOM)        
                       
                        SELECT CASE (ABS(IOR0))
                           CASE (1)
                              VAL = 0.5_EB*(DXI2_COARSE(NM)+DXI2_COARSE(NOM))
                           CASE (3)
                              VAL = 0.5_EB*(DZI2_COARSE(NM)+DZI2_COARSE(NOM))
                        END SELECT
                        SM%A_COARSE(IC, ICO) = VAL
   
                     ENDDO WALLCELL_BLOOP2D_1
   
                     IOFFSET = IOFFSET + SB%NC
   
                  ENDDO
   
               CASE (NSCARC_DIMENSION_THREE)
   
                  SB => SCARC(NM0)%BANDED(NL)

                  DO NM = 1, NMESHES
   
                     NX0 = NX_COARSE(NM)
                     NY0 = NY_COARSE(NM)
                     NZ0 = NZ_COARSE(NM)
                     NC0 = NC_COARSE(NM)
                     NW0 = NW_COARSE(NM)
   
   
                     INTERNALCELL_BLOOP3D_1: DO IC = 1, NC_COARSE(NM)
   
                        IC0 = IC + SM%OFFSET(NM)
   
                        SM%A_COARSE(IC0, IC0) = ACG(IC,ID,NM)
   
                        IF (IC > 1             .AND. ACG(IC,ILX,NM)>0) SM%A_COARSE(IC0, IC0-1)       = ACG(IC, ILX,NM)
                        IF (IC < NC0           .AND. ACG(IC,IUX,NM)>0) SM%A_COARSE(IC0, IC0+1)       = ACG(IC, IUX,NM)
                        IF (IC > NX0           .AND. ACG(IC,ILY,NM)>0) SM%A_COARSE(IC0, IC0-NX0)     = ACG(IC, ILY,NM)
                        IF (IC < NC0-NX0+1     .AND. ACG(IC,IUY,NM)>0) SM%A_COARSE(IC0, IC0+NX0)     = ACG(IC, IUY,NM)
                        IF (IC > NX0*NY0       .AND. ACG(IC,ILZ,NM)>0) SM%A_COARSE(IC0, IC0-NX0*NY0) = ACG(IC, ILZ,NM)
                        IF (IC < NC0-NX0*NY0+1 .AND. ACG(IC,IUZ,NM)>0) SM%A_COARSE(IC0, IC0+NX0*NY0) = ACG(IC, IUZ,NM)
   
                     ENDDO INTERNALCELL_BLOOP3D_1
   
                     WALLCELL_BLOOP3D_1: DO IW = 1, NW0
   
                        NOM = IJKWG( 9, IW, NM)
                        IF (NOM == 0) CYCLE WALLCELL_BLOOP3D_1
   
                        IOR0 = IJKWG( 4, IW, NM)
   
                        II0 = IJKWG( 6 , IW, NM)
                        II1 = IJKWG( 10, IW, NM)
                        II2 = IJKWG( 13, IW, NM)
    
                        JJ0 = IJKWG( 7 , IW, NM)
                        JJ1 = IJKWG( 11, IW, NM)
                        JJ2 = IJKWG( 14, IW, NM)
   
                        KK0 = IJKWG( 8 , IW, NM)
                        KK1 = IJKWG( 12, IW, NM)
                        KK2 = IJKWG( 15, IW, NM)
   
                        !!! For the moment: write error message in case of different resolutions
                        IF (II2-II1.NE.0.OR.KK2-KK1.NE.0) THEN
                           WRITE(*,*) 'Error in SCARC_SETUP_COARSE_MATRIX: Different resolutions on meshes ',NM, NOM
                           STOP
                        ENDIF
   
                        IC0 = (KK0-1)*NX_COARSE(NM)*NY_COARSE(NM)   + (JJ0-1)*NX_COARSE(NM)  + II0
                        IC  = IC0 + SM%OFFSET(NM)
   
                        IC1 = (KK1-1)*NX_COARSE(NOM)*NY_COARSE(NOM) + (JJ1-1)*NX_COARSE(NOM) + II1
                        ICO = IC1 + SM%OFFSET(NOM)
   
                        SELECT CASE (ABS(IOR0))
                           CASE (1)
                              VAL = 0.5_EB*(DXI2_COARSE(NM)+DXI2_COARSE(NOM))
                           CASE (2)
                              VAL = 0.5_EB*(DYI2_COARSE(NM)+DYI2_COARSE(NOM))
                           CASE (3)
                              VAL = 0.5_EB*(DZI2_COARSE(NM)+DZI2_COARSE(NOM))
                        END SELECT
                        SM%A_COARSE(IC, ICO) = VAL
   
                     ENDDO WALLCELL_BLOOP3D_1
   
                     IOFFSET = IOFFSET + SB%NC

                  ENDDO
      
               END SELECT 
      
         CASE (NSCARC_SYSTEM_COMPACT)
            WRITE(*,*) 'NOT YET IMPLEMENTED 2'
            STOP

      END SELECT 
   
   !ENDIF
   
!!! ---------------------------------------------------------------------------------------------
!!! Serial version
!!! ---------------------------------------------------------------------------------------------
ELSE

IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) WRITE(SCARC_LU,*) 'Choose system', TYPE_SYSTEM, TYPE_DIMENSION
   !!! chose system type
   SELECT_SYSTEM: SELECT CASE(TYPE_SYSTEM)
      CASE (NSCARC_SYSTEM_BANDED)

         !!! chose dimension
         SELECT_DIMENSION: SELECT CASE (TYPE_DIMENSION)

            CASE (NSCARC_DIMENSION_TWO)

IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) WRITE(SCARC_LU,*) 'DIM=', 2
               DO NM = 1, NMESHES

                  SB => SCARC(NM)%BANDED(NL)

IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) WRITE(SCARC_LU,*) 'NM=', NM

                  INTERNALCELL_BLOOP2D: DO IC = 1, SB%NC

                     IC0 = IC + SM%OFFSET(NM)

                     SM%A_COARSE(IC0, IC0) = SB%A(IC, ID)

                     IF (IC > 1             .AND. SB%A(IC,ILX)>0) SM%A_COARSE(IC0, IC0-1)     = SB%A(IC, ILX)
                     IF (IC < SB%NC         .AND. SB%A(IC,IUX)>0) SM%A_COARSE(IC0, IC0+1)     = SB%A(IC, IUX)
                     IF (IC > SB%NX         .AND. SB%A(IC,ILZ)>0) SM%A_COARSE(IC0, IC0-SB%NX) = SB%A(IC, ILZ)
                     IF (IC < SB%NC-SB%NX+1 .AND. SB%A(IC,IUZ)>0) SM%A_COARSE(IC0, IC0+SB%NX) = SB%A(IC, IUZ)

                  ENDDO INTERNALCELL_BLOOP2D

                  WALLCELL_BLOOP2D: DO IW = 1, SB%NW

                     NOM = SB%IJKW(9,IW)
                     IF (NOM == 0) CYCLE WALLCELL_BLOOP2D

                     IOR0 = SB%IJKW(4,IW)

                     II0 = SB%IJKW(6 ,IW) 
                     II1 = SB%IJKW(10,IW) 
                     II2 = SB%IJKW(13,IW) 

                     KK0 = SB%IJKW(8 ,IW) 
                     KK1 = SB%IJKW(12,IW) 
                     KK2 = SB%IJKW(15,IW) 

                     IF (TYPE_DEBUG >= NSCARC_DEBUG_MEDIUM) WRITE(SCARC_LU,*) II0, II1, II2, KK0, KK1, KK2
                    
                     !!! For the moment: write error message in case of different resolutions
                     IF (II2-II1.NE.0.OR.KK2-KK1.NE.0) THEN
                        WRITE(*,*) 'Error in SCARC_SETUP_COARSE_MATRIX: Different resolutions on meshes ',NM, NOM
                        STOP
                     ENDIF

                     IC0 = (KK0-1)*NC_COARSE(NM)  + II0  
                     IC  = IC0 + SM%OFFSET(NM)      

                     IC1 = (KK1-1)*NX_COARSE(NOM) + II1  
                     ICO = IC1 + SM%OFFSET(NOM)        
                    
                     SELECT CASE (ABS(IOR0))
                        CASE (1)
                           VAL = 0.5_EB*(DXI2_COARSE(NM)+DXI2_COARSE(NOM))
                        CASE (3)
                           VAL = 0.5_EB*(DZI2_COARSE(NM)+DZI2_COARSE(NOM))
                     END SELECT
                     SM%A_COARSE(IC, ICO) = VAL

                  ENDDO WALLCELL_BLOOP2D

                  IOFFSET = IOFFSET + SB%NC

               ENDDO

            CASE (NSCARC_DIMENSION_THREE)

               DO NM = 1, NMESHES

                  SB => SCARC(NM)%BANDED(NL)

                  INTERNALCELL_BLOOP3D: DO IC = 1, SB%NC

                     IC0 = IC + SM%OFFSET(NM)

                     SM%A_COARSE(IC0, IC0) = SB%A(IC, ID)

                     IF (IC > 1                   .AND. SB%A(IC,ILX)>0) SM%A_COARSE(IC0, IC0-1)     =       SB%A(IC, ILX)
                     IF (IC < SB%NC               .AND. SB%A(IC,IUX)>0) SM%A_COARSE(IC0, IC0+1)           = SB%A(IC, IUX)
                     IF (IC > SB%NX               .AND. SB%A(IC,ILY)>0) SM%A_COARSE(IC0, IC0-SB%NX)       = SB%A(IC, ILY)
                     IF (IC < SB%NC-SB%NX+1       .AND. SB%A(IC,IUY)>0) SM%A_COARSE(IC0, IC0+SB%NX)       = SB%A(IC, IUY)
                     IF (IC > SB%NX*SB%NY         .AND. SB%A(IC,ILZ)>0) SM%A_COARSE(IC0, IC0-SB%NX*SB%NY) = SB%A(IC, ILZ)
                     IF (IC < SB%NC-SB%NX*SB%NY+1 .AND. SB%A(IC,IUZ)>0) SM%A_COARSE(IC0, IC0+SB%NX*SB%NY) = SB%A(IC, IUZ)

                  ENDDO INTERNALCELL_BLOOP3D

                  WALLCELL_BLOOP3D: DO IW = 1, SB%NW

                     NOM = SB%IJKW(9,IW)
                     IF (NOM == 0) CYCLE WALLCELL_BLOOP3D

                     IOR0 = SB%IJKW(4,IW)

                     II0 = SB%IJKW(6 ,IW) 
                     II1 = SB%IJKW(10,IW) 
                     II2 = SB%IJKW(13,IW) 

                     JJ0 = SB%IJKW(7 ,IW) 
                     JJ1 = SB%IJKW(11,IW) 
                     JJ2 = SB%IJKW(14,IW) 
                    
                     KK0 = SB%IJKW(8 ,IW) 
                     KK1 = SB%IJKW(12,IW) 
                     KK2 = SB%IJKW(15,IW) 
                    
                     !!! For the moment: write error message in case of different resolutions
                     IF (II2-II1/=.0 .OR. JJ2-JJ1/=0 .OR. KK2-KK1/=.0) THEN
                        WRITE(*,*) 'Error in SCARC_SETUP_COARSE_MATRIX: Different resolutions on meshes ',NM, NOM
                        STOP
                     ENDIF

                     IC0 = (KK0-1)*NX_COARSE(NM)*NY_COARSE(NM)   + (JJ0-1)*NX_COARSE(NM)  + II0  
                     IC  = IC0 + SM%OFFSET(NM)      

                     IC1 = (KK1-1)*NX_COARSE(NOM)*NY_COARSE(NOM) + (JJ1-1)*NX_COARSE(NOM) + II1  
                     ICO = IC1 + SM%OFFSET(NOM)        
                    
                     IF (TYPE_DEBUG >= NSCARC_DEBUG_MEDIUM) &
                        WRITE(SCARC_LU,'(a,i3,a,2i3,a,2i3,a,2i3,a,2i3)') 'IW=',IW,': ', II1,II2, &
                                                                         ':', JJ1, JJ2, ':', KK1, KK2,':', IC, ICO
                     SELECT CASE (ABS(IOR0))
                        CASE (1)
                           VAL = 0.5_EB*(DXI2_COARSE(NM)+DXI2_COARSE(NOM))
                        CASE (2)
                           VAL = 0.5_EB*(DYI2_COARSE(NM)+DYI2_COARSE(NOM))
                        CASE (3)
                           VAL = 0.5_EB*(DZI2_COARSE(NM)+DZI2_COARSE(NOM))
                     END SELECT
                     SM%A_COARSE(IC, ICO) = VAL

                  ENDDO WALLCELL_BLOOP3D

                  IOFFSET = IOFFSET + SB%NC

               ENDDO

         END SELECT SELECT_DIMENSION

      CASE (NSCARC_SYSTEM_COMPACT)

         DO NM = 1, NMESHES

            SC => SCARC(NM)%COMPACT(NL)

            INTERNALCELL_CLOOP2D: DO IC = 1, SC%NC

               IC0 = IC + SM%OFFSET(NM)

               DO ICOL = SC%A_ROW(IC), SC%A_ROW(IC+1)-1
                  IC2 = SC%A_COL(ICOL)
                  SM%A_COARSE(IC, IC2) = SC%A(ICOL)
               ENDDO

            ENDDO INTERNALCELL_CLOOP2D

            WALLCELL_CLOOP2D: DO IW = 1, SC%NW

               NOM = SC%IJKW(9,IW)
               IF (NOM == 0) CYCLE WALLCELL_CLOOP2D

               SCO => SCARC(NOM)%COMPACT(NL)

               IOR0 = SC%IJKW(4,IW)

               II0 = SC%IJKW(6 ,IW) 
               II1 = SC%IJKW(10,IW) 
               II2 = SC%IJKW(13,IW) 

               KK0 = SC%IJKW(8 ,IW) 
               KK1 = SC%IJKW(12,IW) 
               KK2 = SC%IJKW(15,IW) 
              
               !!! For the moment: write error message in case of different resolutions
               IF (II2-II1.NE.0.OR.KK2-KK1.NE.0) THEN
                  WRITE(*,*) 'Error in SCARC_SETUP_COARSE_MATRIX: Different resolutions on meshes ',NM, NOM
                  STOP
               ENDIF

               IC0 = (KK0-1)*SC%NX  + II0  
               IC  = IC0 + SM%OFFSET(NM)      

               IC1 = (KK1-1)*SCO%NX + II1  
               ICO = IC1 + SM%OFFSET(NOM)        
              
               SELECT CASE (IOR0)
                  CASE (1)
                     !SM%A_COARSE(IC , ICO) = SCO%A(IC1,ILX)    !! noch falsch !!!
                  CASE (-1)
                     !SM%A_COARSE(IC , ICO) = SCO%A(IC1,IUX)    !! noch falsch !!!
                  CASE (3)
                     !SM%A_COARSE(IC , ICO) = SCO%A(IC1,ILZ)    !! noch falsch !!!
                  CASE (-3)
                     !SM%A_COARSE(IC , ICO) = SCO%A(IC1,IUZ)    !! noch falsch !!!
               END SELECT

            ENDDO WALLCELL_CLOOP2D

            IOFFSET = IOFFSET + SC%NC

         ENDDO

   END SELECT SELECT_SYSTEM

ENDIF
      
IF (TYPE_DEBUG >= NSCARC_DEBUG_MEDIUM) THEN
   WRITE(SCARC_LU,*) 'NC_COARSE0:', NC_COARSE0
   WRITE(SCARC_LU,*) 'NA_COARSE:'
   WRITE(SCARC_LU,'(16i4)') (NA_COARSE(NM),NM=1,NMESHES)
   WRITE(SCARC_LU,*) 'NC_COARSE:'
   WRITE(SCARC_LU,'(16i4)') (NC_COARSE(NM),NM=1,NMESHES)
   IF (TYPE_SYSTEM == NSCARC_SYSTEM_BANDED) THEN
      WRITE(SCARC_LU,*) 'NX_COARSE:'
      WRITE(SCARC_LU,'(16i4)') (NX_COARSE(NM),NM=1,NMESHES)
      WRITE(SCARC_LU,*) 'NY_COARSE:'
      WRITE(SCARC_LU,'(16i4)') (NY_COARSE(NM),NM=1,NMESHES)
      WRITE(SCARC_LU,*) 'NZ_COARSE:'
      WRITE(SCARC_LU,'(16i4)') (NZ_COARSE(NM),NM=1,NMESHES)
   ENDIF
   WRITE(SCARC_LU,*) 'A_COARSE before LU'
   DO IC = 1, NC_COARSE0
      WRITE(SCARC_LU,'(25f7.2)') (SM%A_COARSE(IC,IC0), IC0=1, NC_COARSE0)
   ENDDO
ENDIF

!!! ---------------------------------------------------------------------------------------------------
!!! Get LU-decomposition of coarse matrix
!!! ---------------------------------------------------------------------------------------------------
CALL DGETRF(NC_COARSE0, NC_COARSE0, SM%A_COARSE, NC_COARSE0, SM%PIVOT, IERR)
IF (IERR/=0) THEN
   WRITE(*,*) 'Error ', IERR,' in DGETRF while LU-factorization of coarse grid matrix, aborting program  ...'
   STOP
ENDIF


IF (TYPE_DEBUG >= NSCARC_DEBUG_MEDIUM) THEN
   WRITE(SCARC_LU,*) 'A_COARSE after LU:'
   DO IC = 1, NC_COARSE0
      WRITE(SCARC_LU,'(25f7.2)') (SM%A_COARSE(IC,IC0), IC0=1, NC_COARSE0)
   ENDDO
   WRITE(SCARC_LU,*) 'PIVOT:'
   DO IC = 1, NC_COARSE0
      WRITE(SCARC_LU,'(25i7)') (SM%PIVOT(IC,IC0), IC0=1, NC_COARSE0)
   ENDDO
ENDIF

END SUBROUTINE SCARC_SETUP_COARSE_MATRIX



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Allocate matrix for the usual 5-point-stencil (2D) or 7-point-stencil (3D)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SETUP_MATRIX (NM, NL)
INTEGER, INTENT(IN) :: NM, NL
INTEGER :: I, J, K, IC, IP, IW=0, IW0(-3:3), IL0(-3:3), IERR
TYPE (SCARC_BANDED_TYPE) , POINTER :: SB
TYPE (SCARC_COMPACT_TYPE), POINTER :: SC
TYPE (MESH_TYPE), POINTER :: M

M  => MESHES(NM)


SELECT CASE (TYPE_SYSTEM)

!!!----------------------------------------------------------------------------------------------------
!!! Bandwise storage technique:
!!!
!!! The matrix is stored "banded", namely diagonal for diagonal.
!!! Although the sub- and superdiagonals are shorter than the main diagonal,
!!! all bands are stored in the same length and are simply filled with zeros
!!! at redundant positions. The 'wasting' of this storage space is justified
!!! by the possibility to use a much more efficient matrix-vector-multiplication
!!! (which must not use an expensive referencing logic)
!!!
!!! A(.,ILZ)  : lower subdiagonal corresponding to z-coordinate
!!! A(.,ILY)  : lower subdiagonal corresponding to y-coordinate  (only 3D)
!!! A(.,ILX)  : lower subdiagonal corresponding to x-coordinate
!!! A(.,ID )  : main diagonal
!!! A(.,IUX)  : upper subdiagonal corresponding to x-coordinate
!!! A(.,IUY)  : upper subdiagonal corresponding to y-coordinate  (only 3D)
!!! A(.,IUZ)  : upper subdiagonal corresponding to z-coordinate
!!!----------------------------------------------------------------------------------------------------
   CASE (NSCARC_SYSTEM_BANDED)

      SB => SCARC(NM)%BANDED(NL)

      SELECT CASE (TYPE_DIMENSION)

         !!! --------------------- 2D ---------------------
         CASE (NSCARC_DIMENSION_TWO)
      
            SB%NCPL = 5                                              ! number of couplings in matrix stencil
         
            ALLOCATE (SB%A(SB%NC, SB%NCPL), STAT=IERR)
            CALL CHKMEMERR ('SCARC_SETUP_MATRIX', 'SB%A', IERR)
            SB%A = 0.0_EB
      
            ID  = 1                                                  ! pointer to the single subdiagonals
            ILZ = 2
            ILX = 3
            IUX = 4
            IUZ = 5
         
            J = 1
            DO K = 1, SB%NZ
               DO I = 1, SB%NX
              
                  IC = (K-1) * SB%NX + I
              
                  !!! main diagonal
                  SB%A(IC,ID) = - 2.0_EB*SB%DXI2 - 2.0_EB*SB%DZI2   
          
                  ! lower subdiagonal in z-direction
                  IF (K > 1) THEN
                     SB%A(IC,ILZ) = SB%DZI2                       
                  ELSE IF (SB%SUBDIVISION(3,3) > 0) THEN
                     IW0(3) = SB%SUBDIVISION(1,3)
                     IL0(3) = SB%SUBDIVISION(2,3)
                     IF (CELL_WITH_NEIGHBOR(I,J,K,IW,IW0(3),IL0(3),SB%IJKW)) SB%A(IC,ILZ) = -IW
                  ENDIF

                  !!! lower subdiagonal in x-direction
                  IF (I > 1) THEN
                     SB%A(IC,ILX) = SB%DXI2                      
                  ELSE IF (SB%SUBDIVISION(3,1) > 0) THEN
                     IW0(1) = SB%SUBDIVISION(1,1)
                     IL0(1) = SB%SUBDIVISION(2,1)
                     IF (CELL_WITH_NEIGHBOR(I,J,K,IW,IW0(1),IL0(1),SB%IJKW)) SB%A(IC,ILX) = -IW
                  ENDIF

                  ! upper subdiagonal in x-direction
                  IF (I < SB%NX) THEN
                     SB%A(IC,IUX) = SB%DXI2                     
                  ELSE IF (SB%SUBDIVISION(3,-1) > 0) THEN
                     IW0(-1) = SB%SUBDIVISION(1,-1)
                     IL0(-1) = SB%SUBDIVISION(2,-1)
                     IF (CELL_WITH_NEIGHBOR(I,J,K,IW,IW0(-1),IL0(-1),SB%IJKW)) SB%A(IC,IUX) = -IW
                  ENDIF

                  ! upper subdiagonal in z-direction
                  IF (K < SB%NZ) THEN
                     SB%A(IC,IUZ) = SB%DZI2                    
                  ELSE  IF (SB%SUBDIVISION(3,-3) > 0) THEN
                     IW0(-3) = SB%SUBDIVISION(1,-3)
                     IL0(-3) = SB%SUBDIVISION(2,-3)
                     IF (CELL_WITH_NEIGHBOR(I,J,K,IW,IW0(-3),IL0(-3),SB%IJKW)) SB%A(IC,IUZ) = -IW
                  ENDIF
              
               ENDDO
            ENDDO
      
         !!! --------------------- 3D ---------------------
         CASE (NSCARC_DIMENSION_THREE)
    
            SB%NCPL = 7                                              ! number of couplings in matrix stencil

            ALLOCATE (SB%A(SB%NC, SB%NCPL), STAT=IERR)
            CALL CHKMEMERR ('SCARC_SETUP_MATRIX', 'SB%A', IERR)
            SB%A = 0.0_EB
      
            ID  = 1                                                  ! pointer to the single subdiagonals
            ILZ = 2
            ILY = 3
            ILX = 4
            IUX = 5
            IUY = 6
            IUZ = 7

            DO K = 1, SB%NZ
               DO J = 1, SB%NY
                  DO I = 1, SB%NX
             
                     IC = (K-1) * SB%NX * SB%NY + (J-1) * SB%NX + I
             
                     !!! main diagonal
                     SB%A(IC,ID) = - 2.0_EB*SB%DXI2 - 2.0_EB*SB%DYI2 - 2.0_EB*SB%DZI2   
             
                     ! lower subdiagonal in z-direction
                     IF (K > 1) THEN
                        SB%A(IC,ILZ) = SB%DZI2                       
                     ELSE IF (SB%SUBDIVISION(3,3) > 0) THEN
                        IW0(3) = SB%SUBDIVISION(1,3)
                        IL0(3) = SB%SUBDIVISION(2,3)
                        IF (CELL_WITH_NEIGHBOR(I,J,K,IW,IW0(3),IL0(3),SB%IJKW)) SB%A(IC,ILZ) = -IW
IF (TYPE_DEBUG >= NSCARC_DEBUG_MEDIUM) WRITE(SCARC_LU,*) ' 3: ',IW0(3), IL0(3), SB%A(IC,ILZ), IC
                     ENDIF
   
                     ! lower subdiagonal in y-direction
                     IF (J > 1) THEN
                        SB%A(IC,ILY) = SB%DYI2                       
                     ELSE IF (SB%SUBDIVISION(3,2) > 0) THEN
                        IW0(2) = SB%SUBDIVISION(1,2)
                        IL0(2) = SB%SUBDIVISION(2,2)
                        IF (CELL_WITH_NEIGHBOR(I,J,K,IW,IW0(2),IL0(2),SB%IJKW)) SB%A(IC,ILY) = -IW
IF (TYPE_DEBUG >= NSCARC_DEBUG_MEDIUM) WRITE(SCARC_LU,*) ' 2: ',IW0(2), IL0(2), SB%A(IC,ILY), IC
                     ENDIF
   
                     !!! lower subdiagonal in x-direction
                     IF (I > 1) THEN
                        SB%A(IC,ILX) = SB%DXI2                      
                     ELSE IF (SB%SUBDIVISION(3,1) > 0) THEN
                        IW0(1) = SB%SUBDIVISION(1,1)
                        IL0(1) = SB%SUBDIVISION(2,1)
                        IF (CELL_WITH_NEIGHBOR(I,J,K,IW,IW0(1),IL0(1),SB%IJKW)) SB%A(IC,ILX) = -IW
IF (TYPE_DEBUG >= NSCARC_DEBUG_MEDIUM) WRITE(SCARC_LU,*) ' 1: ',IW0(1), IL0(1), SB%A(IC,ILX), IC
                     ENDIF

                     ! upper subdiagonal in x-direction
                     IF (I < SB%NX) THEN
                        SB%A(IC,IUX) = SB%DXI2                     
                     ELSE IF (SB%SUBDIVISION(3,-1) > 0) THEN
                        IW0(-1) = SB%SUBDIVISION(1,-1)
                        IL0(-1) = SB%SUBDIVISION(2,-1)
                        IF (CELL_WITH_NEIGHBOR(I,J,K,IW,IW0(-1),IL0(-1),SB%IJKW)) SB%A(IC,IUX) = -IW
IF (TYPE_DEBUG >= NSCARC_DEBUG_MEDIUM) WRITE(SCARC_LU,*) '-1: ',IW0(-1), IL0(-1), SB%A(IC,IUX), IC
                     ENDIF
   
                     ! upper subdiagonal in y-direction
                     IF (J < SB%NY) THEN
                        SB%A(IC,IUY) = SB%DYI2                     
                     ELSE IF (SB%SUBDIVISION(3,-2) > 0) THEN
                        IW0(-2) = SB%SUBDIVISION(1,-2)
                        IL0(-2) = SB%SUBDIVISION(2,-2)
                        IF (CELL_WITH_NEIGHBOR(I,J,K,IW,IW0(-2),IL0(-2),SB%IJKW)) SB%A(IC,IUY) = -IW
IF (TYPE_DEBUG >= NSCARC_DEBUG_MEDIUM) WRITE(SCARC_LU,*) '-2: ',IW0(-2), IL0(-2), SB%A(IC,IUY), IC
                     ENDIF
   
                     ! upper subdiagonal in z-direction
                     IF (K < SB%NZ) THEN
                        SB%A(IC,IUZ) = SB%DZI2                    
                     ELSE IF (SB%SUBDIVISION(3,-3) > 0) THEN
                        IW0(-3) = SB%SUBDIVISION(1,-3)
                        IL0(-3) = SB%SUBDIVISION(2,-3)
                        IF (CELL_WITH_NEIGHBOR(I,J,K,IW,IW0(-3),IL0(-3),SB%IJKW)) SB%A(IC,IUZ) = -IW
IF (TYPE_DEBUG >= NSCARC_DEBUG_MEDIUM) WRITE(SCARC_LU,*) '-3: ',IW0(-3), IL0(-3), SB%A(IC,IUZ), IC
                     ENDIF
             
                  ENDDO
               ENDDO
            ENDDO
             
      END SELECT
             
      SB%NA   = SB%NC * SB%NCPL                                      ! total number of matrix entries


!!!----------------------------------------------------------------------------------------------------
!!! Compact storage technique:
!!!
!!! Compression technique to store sparse matrices, non-zero entries are stored
!!! in a 1D-vector B(.), row after row, 
!!! each row starts with its diagonal entry followed by the other non-zero entries
!!! In order to identify each element, pointer arrays ROW and COL are needed,
!!! ROW points to the several diagonal entries in vector B(.), 
!!! COL points to the columns which non-zero entries in the matrix stencil
!!!----------------------------------------------------------------------------------------------------
   CASE (NSCARC_SYSTEM_COMPACT)

      SC => SCARC(NM)%COMPACT(NL)

!! ONLY TEMPORARILY !!!
!SC%DXI2=1.0_EB
!SC%DYI2=1.0_EB
!SC%DZI2=1.0_EB

      ALLOCATE (SC%A_ROW(SC%NCG+1000), STAT=IERR)
      CALL CHKMEMERR ('SCARC_SETUP_MATRIX', 'SC%A_ROW', IERR)
      SC%A_ROW = 0
      
      SELECT CASE (TYPE_DIMENSION)
      
         !!! --------------------- 2D ---------------------
         CASE (NSCARC_DIMENSION_TWO)

            SC%NCPL = 5                                        ! number of couplings in matrix stencil
            SC%NA   = SC%NCG * SC%NCPL                         ! first guess for number of matrix entries

            ALLOCATE (SC%A(SC%NA+1000), STAT=IERR)
            CALL CHKMEMERR ('SCARC_SETUP_MATRIX', 'SC%A', IERR)
            SC%A = 0.0_EB
            
            ALLOCATE (SC%A_COL(SC%NA+1000), STAT=IERR)
            CALL CHKMEMERR ('SCARC_SETUP_MATRIX', 'SC%A_COL', IERR)
            SC%A_COL = 0.0_EB
      
            J   = 1
            IP  = 1
            IW0 = 0
            DO K = 1, SC%NZ
               DO I = 1, SC%NX
             
                  IC = (K-1) * SC%NX + I
  
                  !!! main diagonal
                  SC%A(IP)     = - 2.0_EB*SC%DXI2 - 2.0_EB*SC%DZI2      
                  SC%A_ROW(IC) = IP                                     
                  SC%A_COL(IP) = IC                                     
                  IP = IP + 1
         
                  ! lower subdiagonal in z-direction
                  IF (K > 1) THEN
                     SC%A(IP)    = SC%DZI2           
                     SC%A_COL(IP) = IC - SC%NX
                     IP = IP + 1
                  ELSE IF (SC%SUBDIVISION(3,3) > 0) THEN      !!! neighbor in lower k direction ?
                     IW0(3) = SC%SUBDIVISION(1,3)
                     IL0(3) = SC%SUBDIVISION(2,3)
                     IF (CELL_WITH_NEIGHBOR(I,J,K,IW,IW0(3),IL0(3),SC%IJKW)) THEN
                        SC%A_COL(IP) = -IW
                        IP = IP + 1
                     ENDIF
                  ENDIF
   
                  !!! lower subdiagonal in x-direction
                  IF (I > 1) THEN
                     SC%A(IP)    = SC%DXI2                            
                     SC%A_COL(IP) = IC - 1
                     IP = IP + 1
                  ELSE IF (SC%SUBDIVISION(3,1) > 0) THEN      !!! neighbor in lower x direction ?
                     IW0(1) = SC%SUBDIVISION(1,1)
                     IL0(1) = SC%SUBDIVISION(2,1)
                     IF (CELL_WITH_NEIGHBOR(I,J,K,IW,IW0(1),IL0(1),SC%IJKW)) THEN
                        SC%A_COL(IP) = -IW
                        IP = IP + 1
                     ENDIF
                  ENDIF
   
                  ! upper subdiagonal in x-direction
                  IF (I < SC%NX) THEN
                     SC%A(IP)    = SC%DXI2           
                     SC%A_COL(IP) = IC + 1
                     IP = IP + 1
                  ELSE IF (SC%SUBDIVISION(3,-1) > 0) THEN      !!! neighbor in upper x direction ?
                     IW0(-1) = SC%SUBDIVISION(1,-1)
                     IL0(-1) = SC%SUBDIVISION(2,-1)
                     IF (CELL_WITH_NEIGHBOR(I,J,K,IW,IW0(-1),IL0(-1),SC%IJKW)) THEN
                        SC%A_COL(IP) = -IW
                        IP = IP + 1
                     ENDIF
                  ENDIF
   
                  ! upper subdiagonal in z-direction      
                  IF (K < SC%NZ) THEN
                     SC%A(IP)    = SC%DZI2           
                     SC%A_COL(IP) = IC + SC%NX
                     IP = IP + 1
                  ELSE IF (SC%SUBDIVISION(3,-3) > 0) THEN      !!! neighbor in upper k direction ?
                     IW0(-3) = SC%SUBDIVISION(1,-3)
                     IL0(-3) = SC%SUBDIVISION(2,-3)
                     IF (CELL_WITH_NEIGHBOR(I,J,K,IW,IW0(-3),IL0(-3),SC%IJKW)) THEN
                        SC%A_COL(IP) = -IW
                        IP = IP + 1
                     ENDIF
                  ENDIF
   
               ENDDO
            ENDDO
         
         !!! --------------------- 3D ---------------------
         CASE (NSCARC_DIMENSION_THREE)
         
            SC%NCPL = 7                                        ! number of couplings in matrix stencil
            SC%NA   = SC%NCG * SC%NCPL                         ! first guess for number of matrix entries

            ALLOCATE (SC%A(SC%NA+1000), STAT=IERR)
            CALL CHKMEMERR ('SCARC_SETUP_MATRIX', 'SC%A', IERR)
            SC%A = 0.0_EB
            
            ALLOCATE (SC%A_COL(SC%NA+1000), STAT=IERR)
            CALL CHKMEMERR ('SCARC_SETUP_MATRIX', 'SC%A_COL', IERR)
            SC%A_COL = 0.0_EB
      

            !!! Compute single matrix entries and corresponding row and column pointers
            !!! Along internal boundaries use placeholders for the neighboring matrix entries
            !!! which will be communicated in a following step
            IP  = 1
            IW0 = 0
            DO K = 1, SC%NZ
               DO J = 1, SC%NY
                  DO I = 1, SC%NX
             
                     IC = (K-1) * SC%NX * SC%NY + (J-1) * SC%NX + I
             
                     !!! main diagonal
                     SC%A(IP)     = - 2.0_EB*SC%DXI2 - 2.0*SC%DYI2 - 2.0_EB*SC%DZI2  
                     SC%A_ROW(IC) = IP                                     
                     SC%A_COL(IP) = IC 
                     IP = IP + 1
            
                     ! lower subdiagonal in z-direction      !!! neighbor in lower k direction ?
                     IF (K > 1) THEN
                        SC%A(IP) = SC%DZI2           
                        SC%A_COL(IP) = IC - SC%NX * SC%NY
                        IP = IP + 1
                     ELSE IF (SC%SUBDIVISION(3,3) > 0) THEN
                        IW0(3) = SC%SUBDIVISION(1,3)
                        IL0(3) = SC%SUBDIVISION(2,3)
                        IF (CELL_WITH_NEIGHBOR(I,J,K,IW,IW0(3),IL0(3),SC%IJKW)) THEN
                           SC%A_COL(IP) = -IW
                           IP = IP + 1
                        ENDIF
                     ENDIF
     
                     !!! lower subdiagonal in y-direction
                     IF (J > 1) THEN
                        SC%A(IP) =  SC%DYI2                            
                        SC%A_COL(IP) = IC - SC%NX
                        IP = IP + 1
                     ELSE IF (SC%SUBDIVISION(3,2) > 0) THEN      !!! neighbor in lower y direction ?
                        IW0(2) = SC%SUBDIVISION(1,2)
                        IL0(2) = SC%SUBDIVISION(2,2)
                        IF (CELL_WITH_NEIGHBOR(I,J,K,IW,IW0(2),IL0(2),SC%IJKW)) THEN
                           SC%A_COL(IP) = -IW
                           IP = IP + 1
                        ENDIF
                     ENDIF
     
                     !!! lower subdiagonal in x-direction
                     IF (I > 1) THEN
                        SC%A(IP) = SC%DXI2                            
                        SC%A_COL(IP) = IC - 1
                        IP = IP + 1
                     ELSE IF (SC%SUBDIVISION(3,1) > 0) THEN      !!! neighbor in lower x direction ?
                        IW0(1) = SC%SUBDIVISION(1,1)
                        IL0(1) = SC%SUBDIVISION(2,1)
                        IF (CELL_WITH_NEIGHBOR(I,J,K,IW,IW0(1),IL0(1),SC%IJKW)) THEN
                           SC%A_COL(IP) = -IW
                           IP = IP + 1
                        ENDIF
                     ENDIF
     
                     ! upper subdiagonal in x-direction
                     IF (I < SC%NX) THEN
                        SC%A(IP) = SC%DXI2           
                        SC%A_COL(IP) = IC + 1
                        IP = IP + 1
                     ELSE IF (SC%SUBDIVISION(3,-1) > 0) THEN      !!! neighbor in upper x direction ?
                        IW0(-1) = SC%SUBDIVISION(1,-1)
                        IL0(-1) = SC%SUBDIVISION(2,-1)
                        IF (CELL_WITH_NEIGHBOR(I,J,K,IW,IW0(-1),IL0(-1),SC%IJKW)) THEN
                           SC%A_COL(IP) = -IW
                           IP = IP + 1
                        ENDIF
                     ENDIF
     
                     ! upper subdiagonal in y-direction
                     IF (J < SC%NY) THEN
                        SC%A(IP) = SC%DYI2           
                        SC%A_COL(IP) = IC + SC%NX
                        IP = IP + 1
                     ELSE IF (SC%SUBDIVISION(3,-2) > 0) THEN      !!! neighbor in upper y direction ?
                        IW0(-2) = SC%SUBDIVISION(1,-2)
                        IL0(-2) = SC%SUBDIVISION(2,-2)
                        IF (CELL_WITH_NEIGHBOR(I,J,K,IW,IW0(-2),IL0(-2),SC%IJKW)) THEN
                           SC%A_COL(IP) = -IW
                           IP = IP + 1
                        ENDIF
                     ENDIF
     
                     ! upper subdiagonal in z-direction      
                     IF (K < SC%NZ) THEN
                        SC%A(IP) = SC%DZI2           
                        SC%A_COL(IP) = IC + SC%NX * SC%NY
                        IP = IP + 1
                     ELSE IF (SC%SUBDIVISION(3,-3) > 0) THEN      !!! neighbor in upper z direction ?
                        IW0(-3) = SC%SUBDIVISION(1,-3)
                        IL0(-3) = SC%SUBDIVISION(2,-3)
                        IF (CELL_WITH_NEIGHBOR(I,J,K,IW,IW0(-3),IL0(-3),SC%IJKW)) THEN
                           SC%A_COL(IP) = -IW
                           IP = IP + 1
                        ENDIF
                     ENDIF
   
                  ENDDO
               ENDDO
            ENDDO
          
      END SELECT
      SC%A_ROW(SC%NC+1) = IP
      SC%NA            = IP -1                                     ! set correct number of matrix entries

IF (TYPE_DEBUG >= NSCARC_DEBUG_MEDIUM) THEN
   WRITE(SCARC_LU,*) '======================== NM=',NM
   WRITE(SCARC_LU,*) 'SUBDIVISION(1,:):'
   WRITE(SCARC_LU,'(7i8)') (SC%SUBDIVISION(1,I), I=-3,3)
   WRITE(SCARC_LU,*) 'SUBDIVISION(2,:):'
   WRITE(SCARC_LU,'(7i8)') (SC%SUBDIVISION(2,I), I=-3,3)
   WRITE(SCARC_LU,*) 'SUBDIVISION(3,:):'
   WRITE(SCARC_LU,'(7i8)') (SC%SUBDIVISION(3,I), I=-3,3)
ENDIF

END SELECT

END SUBROUTINE SCARC_SETUP_MATRIX


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Determine if cell (I,J,K) has a neighbor and, if yes, save corresponding IW-value
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
LOGICAL FUNCTION CELL_WITH_NEIGHBOR(I, J, K, IW, IW0, IL0, IJKW0)
INTEGER, INTENT(IN)    :: I, J, K, IW0, IL0
INTEGER, INTENT(INOUT) :: IW
INTEGER, DIMENSION(:,:), INTENT(IN) :: IJKW0

CELL_WITH_NEIGHBOR = .FALSE.

SEARCH_WALLCELL_LOOP: DO IW = IW0, IW0+IL0-1
  IF (I == IJKW0(6,IW) .AND. J == IJKW0(7,IW) .AND. K == IJKW0(8,IW) .AND. IJKW0(9,IW) /= 0) THEN
     CELL_WITH_NEIGHBOR = .TRUE.
     EXIT SEARCH_WALLCELL_LOOP
  ENDIF
ENDDO SEARCH_WALLCELL_LOOP

RETURN
END FUNCTION CELL_WITH_NEIGHBOR


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Set pointer for different structures on level NL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SETUP_BOUNDARY (NM, NL)
INTEGER, INTENT(IN) :: NM, NL
INTEGER :: I, J, K, IOR0, IW, IC, NOM, IP
REAL(EB) :: DBC
TYPE (SCARC_BANDED_TYPE) , POINTER :: SB
TYPE (SCARC_COMPACT_TYPE), POINTER :: SC


SELECT_SYSTEM: SELECT CASE (TYPE_SYSTEM)

   !!!----------------------------------------------------------------------------------------------------
   !!! Banded system
   !!!----------------------------------------------------------------------------------------------------
   CASE (NSCARC_SYSTEM_BANDED)

      SB  => SCARC(NM)%BANDED(NL)

      SELECT_DIMENSION_BANDED: SELECT CASE (TYPE_DIMENSION)
      
         !!!------------------ 2D ---------------------------
         CASE (NSCARC_DIMENSION_TWO)
      
            WALLCELL_LOOP2D_BANDED: DO IW = 1, SB%NW
             
               IOR0 = SB%IJKW(4,IW)
               IF (ABS(IOR0) == 2) CYCLE            ! 2D: in case of y-boundary cycle
         
               I   = SB%IJKW(6, IW)
               K   = SB%IJKW(8, IW)

               NOM = SB%IJKW(9,IW)
         
               SB%ADJACENT_CELL(IW) = (K-1)*SB%NX + I
      
               SELECT CASE (IOR0)
                  CASE (1)
                     IC = (K-1) * SB%NX + I
                     DBC= SB%DXI2
                  CASE (-1)
                     IC = K * SB%NX
                     DBC= SB%DXI2
                  CASE (3)
                     IC = I
                     DBC= SB%DZI2
                  CASE (-3)
                     IC = (SB%NZ-1) * SB%NX + I
                     DBC= SB%DZI2
               END SELECT
         
               SELECT CASE (SB%BC_INDEX(IW))
                  CASE (DIRICHLET)                        ! set Dirichlet BC's along open boundary cells
                    SB%A(IC,ID) = SB%A(IC,ID) - DBC
                  !CASE (INTERNAL)                        ! do nothing along internal boundaries (only debugging)
                  CASE (NEUMANN)                          ! set Neumann BC's at all other nodes
                    SB%A(IC,ID) = SB%A(IC,ID) + DBC
               END SELECT
         
            ENDDO WALLCELL_LOOP2D_BANDED

        
         !!!------------------ 3D ---------------------------
         CASE (NSCARC_DIMENSION_THREE)
      
            WALLCELL_LOOP3D_BANDED: DO IW = 1, SB%NW
         
               IOR0 = SB%IJKW(4,IW)

               I    = SB%IJKW(6, IW)
               J    = SB%IJKW(7, IW)
               K    = SB%IJKW(8, IW)

               NOM  = SB%IJKW(9,IW)
         
               SB%ADJACENT_CELL(IW) = (K-1)*SB%NX*SB%NY + (J-1)*SB%NX + I
      
               SELECT CASE (IOR0)
                  CASE (1)
                     IC = (K-1) * SB%NX * SB%NY + (J-1) * SB%NX + I
                     DBC= SB%DXI2
                  CASE (-1)
                     IC = (K-1) * SB%NX * SB%NY + J * SB%NX 
                     DBC= SB%DXI2
                  CASE (2)
                     IC = (K-1) * SB%NX * SB%NY + I
                     DBC= SB%DYI2
                  CASE (-2)
                     IC = (K-1) * SB%NX * SB%NY + (SB%NY-1) * SB%NX + I
                     DBC= SB%DYI2
                  CASE (3)
                     IC = (J-1) * SB%NX + I
                     DBC= SB%DZI2
                  CASE (-3)
                     IC = (SB%NZ-1) * SB%NX * SB%NY + (J-1) * SB%NX + I
                     DBC= SB%DZI2
               END SELECT
         
               SELECT CASE (SB%BC_INDEX(IW))
                  CASE (DIRICHLET)                     ! set Dirichlet BC's at open and null boundary cells
                     SB%A(IC,ID) = SB%A(IC,ID) - DBC
                  !CASE (INTERNAL)                     ! do nothing along internal boundaries (only debugging)
                  CASE (NEUMANN)                       ! set Neumann BC's at all other cells
                     SB%A(IC,ID) = SB%A(IC,ID) + DBC
               END SELECT
         
            ENDDO WALLCELL_LOOP3D_BANDED
           
      END SELECT SELECT_DIMENSION_BANDED


   !!!----------------------------------------------------------------------------------------------------
   !!! Compact system
   !!!----------------------------------------------------------------------------------------------------
   CASE (NSCARC_SYSTEM_COMPACT)

      SC => SCARC(NM)%COMPACT(NL)

      SELECT_DIMENSION_COMPACT: SELECT CASE (TYPE_DIMENSION)
      
         !!!------------------ 2D ---------------------------
         CASE (NSCARC_DIMENSION_TWO)
      
            WALLCELL_LOOP2D_COMPACT: DO IW = 1, SC%NW
             
               IOR0 = SC%IJKW(4,IW)
               IF (ABS(IOR0) == 2) CYCLE            ! 2D: in case of y-boundary cycle
         
               I    = SC%IJKW(6, IW)
               K    = SC%IJKW(8, IW)

               NOM  = SC%IJKW(9,IW)
         
               SC%ADJACENT_CELL(IW) = (K-1)*SC%NX + I
      
               SELECT CASE (IOR0)
                  CASE (1)
                     IC = (K-1) * SC%NX + I
                     DBC= SC%DXI2
                  CASE (-1)
                     IC = K * SC%NX
                     DBC= SC%DXI2
                  CASE (3)
                     IC = I
                     DBC= SC%DZI2
                  CASE (-3)
                     IC = (SC%NZ-1) * SC%NX + I
                     DBC= SC%DZI2
               END SELECT
         
               IP = SC%A_ROW(IC)
               SELECT CASE (SC%BC_INDEX(IW))
                  CASE (DIRICHLET)                        ! set Dirichlet BC's along open boundary cells
                     SC%A(IP) = SC%A(IP) - DBC
                  !CASE (INTERNAL)                        ! do nothing along internal boundaries (only debugging)
                  CASE (NEUMANN)                          ! set Neumann BC's at all other nodes
                     SC%A(IP) = SC%A(IP) + DBC
               END SELECT
         
            ENDDO WALLCELL_LOOP2D_COMPACT
         
      
         !!!------------------ 2D ---------------------------
         CASE (NSCARC_DIMENSION_THREE)
      
            WALLCELL_LOOP3D_COMPACT: DO IW = 1, SC%NW
         
               IOR0 = SC%IJKW(4,IW)

               I    = SC%IJKW(6, IW)
               J    = SC%IJKW(7, IW)
               K    = SC%IJKW(8, IW)

               NOM  = SC%IJKW(9,IW)
         
               SC%ADJACENT_CELL(IW) = (K-1)*SC%NX*SC%NY + (J-1)*SC%NX + I
      
               SELECT CASE (IOR0)
                  CASE (1)
                     IC = (K-1) * SC%NX * SC%NY + (J-1) * SC%NX + I
                     DBC= SC%DXI2
                  CASE (-1)
                     IC = (K-1) * SC%NX * SC%NY + J * SC%NX 
                     DBC= SC%DXI2
                  CASE (2)
                     IC = (K-1) * SC%NX * SC%NY + I
                     DBC= SC%DYI2
                  CASE (-2)
                     IC = (K-1) * SC%NX * SC%NY + (SC%NY-1) * SC%NX + I
                     DBC= SC%DYI2
                  CASE (3)
                     IC = (J-1) * SC%NX + I
                     DBC= SC%DZI2
                  CASE (-3)
                     IC = (SC%NZ-1) * SC%NX * SC%NY + (J-1) * SC%NX + I
                     DBC= SC%DZI2
               END SELECT
         
               IP = SC%A_ROW(IC)
               SELECT CASE (SC%BC_INDEX(IW))
                  CASE (DIRICHLET)                        ! set Dirichlet BC's at open and null boundary cells
                     SC%A(IP) = SC%A(IP) - DBC
                  !CASE (INTERNAL)                        ! do nothing along internal boundaries (only debugging)
                  CASE (NEUMANN)                          ! set Neumann BC's at all other cells
                     SC%A(IP) = SC%A(IP) + DBC
               END SELECT
         
            ENDDO WALLCELL_LOOP3D_COMPACT
           
      END SELECT SELECT_DIMENSION_COMPACT

END SELECT SELECT_SYSTEM

END SUBROUTINE SCARC_SETUP_BOUNDARY


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Initialize global solver methods (CG/BICG/GMG/AMG) and corresponding level structures
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SETUP_VECTORS
INTEGER :: IERR, IBP1, JBP1, KBP1
TYPE (SCARC_BANDED_TYPE) , POINTER :: SB
TYPE (SCARC_COMPACT_TYPE), POINTER :: SC
TYPE (SCARC_PRECON_TYPE) , POINTER :: SP
INTEGER :: NM, NL

IERR = 0

MESHES_LOOP: DO NM = NMESHES_MIN, NMESHES_MAX
   LEVEL_LOOP: DO NL = NLEVEL_MIN, NLEVEL_MAX

      SELECT_STORAGE: SELECT CASE (TYPE_SYSTEM) 
      
         !!!-------------------------------------------------------------------------------------------------
         !!! Bandwise storage technique
         !!!-------------------------------------------------------------------------------------------------
         CASE (NSCARC_SYSTEM_BANDED)
            
            !!! let SB point to SYSTEM-BANDED structure
            SB => SCARC(NM)%BANDED(NL)
         
            IBP1 = SB%NX + 1
            JBP1 = SB%NY + 1
            KBP1 = SB%NZ + 1
      
            SELECT_METHOD_BANDED: SELECT CASE (TYPE_METHOD)
         
               !!! working and auxiliary vectors for global CG/BICG-method
               CASE (NSCARC_METHOD_KRYLOV)
         
                  ALLOCATE (SB%X(0:IBP1, 0:JBP1, 0:KBP1), STAT=IERR)
                  CALL CHKMEMERR ('SCARC', 'X', IERR)
                  SB%X = 0.0_EB
             
                  ALLOCATE (SB%F(0:IBP1, 0:JBP1, 0:KBP1), STAT=IERR)
                  CALL CHKMEMERR ('SCARC', 'B', IERR)
                  SB%F = 0.0_EB
             
                  ALLOCATE (SB%D(0:IBP1, 0:JBP1, 0:KBP1), STAT=IERR)
                  CALL CHKMEMERR ('SCARC', 'D', IERR)
                  SB%D = 0.0_EB
         
                  ALLOCATE (SB%W(0:IBP1, 0:JBP1, 0:KBP1), STAT=IERR)
                  CALL CHKMEMERR ('SCARC', 'W', IERR)
                  SB%W = 0.0_EB
             
                  ALLOCATE (SB%G(0:IBP1, 0:JBP1, 0:KBP1), STAT=IERR)
                  CALL CHKMEMERR ('SCARC', 'G', IERR)
                  SB%G = 0.0_EB
             
                  ALLOCATE (SB%Y(0:IBP1, 0:JBP1, 0:KBP1), STAT=IERR)
                  CALL CHKMEMERR ('SCARC', 'Y', IERR)
                  SB%Y = 0.0_EB
         
                  ALLOCATE (SB%Z(0:IBP1, 0:JBP1, 0:KBP1), STAT=IERR)
                  CALL CHKMEMERR ('SCARC', 'Z', IERR)
                  SB%Z = 0.0_EB
         
         
                  IF (TYPE_PRECON == NSCARC_PRECON_MULTIGRID) THEN
               
                     ALLOCATE (SB%X2(0:IBP1, 0:JBP1, 0:KBP1), STAT=IERR)
                     CALL CHKMEMERR ('SCARC', 'X2', IERR)
                     SB%X2 = 0.0_EB
               
                     ALLOCATE (SB%F2(0:IBP1, 0:JBP1, 0:KBP1), STAT=IERR)
                     CALL CHKMEMERR ('SCARC', 'F2', IERR)
                     SB%F2 = 0.0_EB
               
                     ALLOCATE (SB%D2(0:IBP1, 0:JBP1, 0:KBP1), STAT=IERR)
                     CALL CHKMEMERR ('SCARC', 'D2', IERR)
                     SB%D2 = 0.0_EB
         
                     IF (NL==NLEVEL_MAX) THEN
            
                        ALLOCATE (SB%W2(0:IBP1, 0:JBP1, 0:KBP1), STAT=IERR)
                        CALL CHKMEMERR ('SCARC', 'W2', IERR)
                        SB%W2 = 0.0_EB
                
                        ALLOCATE (SB%G2(0:IBP1, 0:JBP1, 0:KBP1), STAT=IERR)
                        CALL CHKMEMERR ('SCARC', 'G2', IERR)
                        SB%G2 = 0.0_EB
                
                        ALLOCATE (SB%Y2(0:IBP1, 0:JBP1, 0:KBP1), STAT=IERR)
                        CALL CHKMEMERR ('SCARC', 'Y2', IERR)
                        SB%Y2 = 0.0_EB
            
                     ENDIF
                   
                  ENDIF
         
                  IF (TYPE_PRECON == NSCARC_PRECON_FFT) THEN 
                     SP => SCARC(NM)%PRECON(NL)
                     ALLOCATE (SP%FFT(1:IBP1, 1:JBP1, 1:KBP1), STAT=IERR)
                     CALL CHKMEMERR ('SCARC', 'FFT', IERR)
                     SP%FFT = 0.0_EB
                  ENDIF
         
               !!! working and auxiliary vectors for global GMG/AMG-method
               CASE (NSCARC_METHOD_MULTIGRID)
         
                  ALLOCATE (SB%X(0:IBP1, 0:JBP1, 0:KBP1), STAT=IERR)
                  CALL CHKMEMERR ('SCARC', 'X', IERR)
                  SB%X = 0.0_EB
             
                  ALLOCATE (SB%F(0:IBP1, 0:JBP1, 0:KBP1), STAT=IERR)
                  CALL CHKMEMERR ('SCARC', 'B', IERR)
                  SB%F = 0.0_EB
             
                  ALLOCATE (SB%D(0:IBP1, 0:JBP1, 0:KBP1), STAT=IERR)
                  CALL CHKMEMERR ('SCARC', 'D', IERR)
                  SB%D = 0.0_EB
         
                  IF (NL==NLEVEL_MAX) THEN
         
                     ALLOCATE (SB%W(0:IBP1, 0:JBP1, 0:KBP1), STAT=IERR)
                     CALL CHKMEMERR ('SCARC', 'W', IERR)
                     SB%W = 0.0_EB
             
                     ALLOCATE (SB%G(0:IBP1, 0:JBP1, 0:KBP1), STAT=IERR)
                     CALL CHKMEMERR ('SCARC', 'G', IERR)
                     SB%G = 0.0_EB
             
                     ALLOCATE (SB%Y(0:IBP1, 0:JBP1, 0:KBP1), STAT=IERR)
                     CALL CHKMEMERR ('SCARC', 'Y', IERR)
                     SB%Y = 0.0_EB
         
                  ENDIF
         
                  IF (TYPE_PRECON == NSCARC_PRECON_FFT) THEN
                     SP => SCARC(NM)%PRECON(NL)
                     ALLOCATE (SP%FFT(1:IBP1, 1:JBP1, 1:KBP1), STAT=IERR)
                     CALL CHKMEMERR ('SCARC', 'FFT', IERR)
                     SP%FFT = 0.0_EB
                  ENDIF
         
            END SELECT SELECT_METHOD_BANDED
            
      
         !!!-------------------------------------------------------------------------------------------------
         !!! Compact storage technique
         !!!-------------------------------------------------------------------------------------------------
         CASE (NSCARC_SYSTEM_COMPACT)
      
            SC => SCARC(NM)%COMPACT(NL)
         
            SELECT_METHOD_COMPACT: SELECT CASE (TRIM(SCARC_METHOD))
         
               !!! working and auxiliary vectors for global CG/BICG-method
               CASE ('KRYLOV')
         
                  ALLOCATE (SC%X(SC%NCE), STAT=IERR)
                  CALL CHKMEMERR ('SCARC', 'X', IERR)
                  SC%X = 0.0_EB
             
                  ALLOCATE (SC%F(SC%NCE), STAT=IERR)
                  CALL CHKMEMERR ('SCARC', 'B', IERR)
                  SC%F = 0.0_EB
             
                  ALLOCATE (SC%D(SC%NCE), STAT=IERR)
                  CALL CHKMEMERR ('SCARC', 'D', IERR)
                  SC%D = 0.0_EB
         
                  ALLOCATE (SC%W(SC%NCE), STAT=IERR)
                  CALL CHKMEMERR ('SCARC', 'W', IERR)
                  SC%W = 0.0_EB
             
                  ALLOCATE (SC%G(SC%NCE), STAT=IERR)
                  CALL CHKMEMERR ('SCARC', 'G', IERR)
                  SC%G = 0.0_EB
             
                  ALLOCATE (SC%Y(SC%NCE), STAT=IERR)
                  CALL CHKMEMERR ('SCARC', 'Y', IERR)
                  SC%Y = 0.0_EB
         
                  ALLOCATE (SC%Z(SC%NCE), STAT=IERR)
                  CALL CHKMEMERR ('SCARC', 'Z', IERR)
                  SC%Z = 0.0_EB
      
                  IF (TYPE_PRECON == NSCARC_PRECON_MULTIGRID) THEN
      
                     ALLOCATE (SC%X2(SC%NCE), STAT=IERR)
                     CALL CHKMEMERR ('SCARC', 'X2', IERR)
                     SC%X2 = 0.0_EB
               
                     ALLOCATE (SC%F2(SC%NCE), STAT=IERR)
                     CALL CHKMEMERR ('SCARC', 'F2', IERR)
                     SC%F2 = 0.0_EB
               
                     ALLOCATE (SC%D2(SC%NCE), STAT=IERR)
                     CALL CHKMEMERR ('SCARC', 'D2', IERR)
                     SC%D2 = 0.0_EB
         
                     IF (NL==NLEVEL_MAX) THEN
            
                        ALLOCATE (SC%W2(SC%NCE), STAT=IERR)
                        CALL CHKMEMERR ('SCARC', 'W2', IERR)
                        SC%W2 = 0.0_EB
                
                        ALLOCATE (SC%G2(SC%NCE), STAT=IERR)
                        CALL CHKMEMERR ('SCARC', 'G2', IERR)
                        SC%G2 = 0.0_EB
                
                        ALLOCATE (SC%Y2(SC%NCE), STAT=IERR)
                        CALL CHKMEMERR ('SCARC', 'Y2', IERR)
                        SC%Y2 = 0.0_EB
            
                     ENDIF
                   
                  ENDIF
         
                  IF (TYPE_PRECON == NSCARC_PRECON_FFT) THEN
                     SP => SCARC(NM)%PRECON(NL)
                     ALLOCATE (SP%FFT(1:SC%NX+1, 1:SC%NY+1, 1:SC%NZ+1), STAT=IERR)
                     CALL CHKMEMERR ('SCARC', 'FFT', IERR)
                     SP%FFT = 0.0_EB
                  ENDIF
         
               !!! working and auxiliary vectors for global GMG/AMG-method
               CASE ('MULTIGRID')
         
                  ALLOCATE (SC%X(SC%NCE), STAT=IERR)
                  CALL CHKMEMERR ('SCARC', 'X', IERR)
                  SC%X = 0.0_EB
             
                  ALLOCATE (SC%F(SC%NCE), STAT=IERR)
                  CALL CHKMEMERR ('SCARC', 'F', IERR)
                  SC%F = 0.0_EB
             
                  ALLOCATE (SC%D(SC%NCE), STAT=IERR)
                  CALL CHKMEMERR ('SCARC', 'D', IERR)
                  SC%D = 0.0_EB
         
                  IF (NL==NLEVEL_MAX) THEN
         
                     ALLOCATE (SC%W(SC%NCE), STAT=IERR)
                     CALL CHKMEMERR ('SCARC', 'W', IERR)
                     SC%W = 0.0_EB
             
                     ALLOCATE (SC%G(SC%NCE), STAT=IERR)
                     CALL CHKMEMERR ('SCARC', 'G', IERR)
                     SC%G = 0.0_EB
             
                     ALLOCATE (SC%Y(SC%NCE), STAT=IERR)
                     CALL CHKMEMERR ('SCARC', 'Y', IERR)
                     SC%Y = 0.0_EB
         
                  ENDIF
         
                  IF (TYPE_PRECON == NSCARC_PRECON_FFT) THEN
                     SP => SCARC(NM)%PRECON(NL)
                     ALLOCATE (SP%FFT(1:SC%NX+1, 1:SC%NY+1, 1:SC%NZ+1), STAT=IERR)
                     CALL CHKMEMERR ('SCARC', 'FFT', IERR)
                     SP%FFT = 0.0_EB
                  ENDIF
         
            END SELECT SELECT_METHOD_COMPACT
      
      END SELECT SELECT_STORAGE

   ENDDO LEVEL_LOOP
ENDDO  MESHES_LOOP

END SUBROUTINE SCARC_SETUP_VECTORS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Initialize global 3D-solver methods (cg/mg)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SETUP_COARSENING 
INTEGER :: IERR, NL, NM, NOM, IC, ICC, NLEN(3), ILEN
TYPE (SCARC_TYPE), POINTER :: S
TYPE (SCARC_COMPACT_TYPE), POINTER :: SCF, SCC
TYPE (OSCARC_COMPACT_TYPE), POINTER :: OS, OSCF, OSCC

IERR = 0

IF (TYPE_MULTIGRID /= NSCARC_MULTIGRID_ALGEBRAIC) RETURN

CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_ACELL , NLEVEL_MIN,'SETUP_COARSENING','ADJACENT_CELL')
CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_MATRIX, NLEVEL_MIN,'SETUP_COARSENING', 'MATRIX')


!!! Determine number of multigrid levels
LEVEL_LOOP: DO NL = NLEVEL_MIN, NLEVEL_MAX-1
   
IF (TYPE_DEBUG >= NSCARC_DEBUG_MEDIUM) WRITE(SCARC_LU,*) '=============== PREPARING LEVEL ', NL

   !!!-------------------------------------------------------------------------------------------------
   !!! Determine coarser meshes corresponding to requested coarsening strategy
   !!!  --- allocate necessary arrays
   !!!  --- setup measures of single cells
   !!!  --- setup celltypes of single cells
   !!!  --- setup sizes of transformation matrices (prolongation/restriction)
   !!!-------------------------------------------------------------------------------------------------
   DO NM = NMESHES_MIN, NMESHES_MAX

      SCF => SCARC(NM)%COMPACT(NL)                          ! system compact on fine level

      ALLOCATE (SCF%MEASURE(1:SCF%NCE), STAT=IERR)
      CALL CHKMEMERR ('SCARC_SETUP_COARSENING', 'SCF%MEASURE', IERR)
      SCF%MEASURE = 0.0_EB

      ALLOCATE (SCF%CELLTYPE(1:SCF%NCE), STAT=IERR)
      CALL CHKMEMERR ('SCARC_SETUP_COARSENING', 'SCF%CELLTYPE', IERR)
      SCF%CELLTYPE = NSCARC_CELLTYPE_NONE

   ENDDO
    
   CALL SCARC_SETUP_MEASURES (NL)
   CALL SCARC_SETUP_CELLTYPES (NL)

   CALL SCARC_LATEX_INFORMATION(NSCARC_LATEX_STAGGERED,NL)
   CALL SCARC_LATEX_INFORMATION(NSCARC_LATEX_EQUAL,NL)
   !CALL SCARC_LATEX_INFORMATION(NSCARC_LATEX_NUMBER,NL)

   CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_CELLTYPE,  NL)
IF (TYPE_DEBUG >= NSCARC_DEBUG_MEDIUM) WRITE(SCARC_LU,*) 'CELLTYPE BEFORE PROLONG '
   CALL SCARC_DEBUG_QUANTITY(NSCARC_DEBUG_CELLTYPE, NL, 'SETUP_CELLTYPES', 'CELLTYPEE')

   !!!-------------------------------------------------------------------------------------------------
   !!! Allocate and define grid transfer matrices
   !!!-------------------------------------------------------------------------------------------------
   !CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_WEIGHTS,  NL)

   DO NM = NMESHES_MIN, NMESHES_MAX

WRITE(SCARC_LU,*) '==================== SUSISUSI: NM=',NM
      SCF => SCARC(NM)%COMPACT(NL)            

      SCF%NCCE = SCF%NCC
      SCF%NCW  = SCF%NCE-SCF%NC
WRITE(SCARC_LU,*) 'HALLOSUSI: ',SCF%NCW, SCF%NCE, SCF%NC, SCF%NCCE
   DO IC = SCF%NC+1, SCF%NCE
WRITE(SCARC_LU,*) 'HALLO: IC=',IC,SCF%CELLTYPE(IC), SCF%NCW, SCF%NCE, SCF%NC
      SCF%NCCE = SCF%NCCE + 1
      IF (SCF%CELLTYPE(IC)>0) THEN
         SCF%NP = SCF%NP + 1
      ELSE
         SCF%NP = SCF%NP + 8                      !!! ONLY TEMPORARILY
      ENDIF
      SCF%CELLTYPE(IC)=IC
WRITE(SCARC_LU,*) 'HALLO: CELLTYPE(',IC,')=',SCF%CELLTYPE(IC), SCF%NCCE, SCF%NP
   ENDDO
   SCF%NR = SCF%NP

WRITE(SCARC_LU,*) 'HALLO1: NCCE=',SCF%NCCE, NM

      !!! allocate prolongation matrix including row and column pointers
      ALLOCATE (SCF%P(SCF%NP), STAT=IERR)
      CALL CHKMEMERR ('SCARC_SETUP_COARSENING', 'P', IERR)
      SCF%P = 0.0_EB
   
WRITE(*,*) 'HALLO: DIMENSIONIERUNG CHECKEN !!!'
      ALLOCATE (SCF%P_ROW(2*SCF%NCE-SCF%NC), STAT=IERR)
      CALL CHKMEMERR ('SCARC_SETUP_COARSENING', 'P_ROW', IERR)
      SCF%P_ROW = 0.0_EB
   
      ALLOCATE (SCF%P_COL(SCF%NP), STAT=IERR)
      CALL CHKMEMERR ('SCARC_SETUP_COARSENING', 'P_COL', IERR)
      SCF%P_COL = 0.0_EB
   
WRITE(*,*) 'HALLO2: AUCH HIER !!!'
      ALLOCATE (SCF%P_PTR(2*SCF%NCE-SCF%NC,2), STAT=IERR)
      CALL CHKMEMERR ('SCARC_SETUP_COARSENING', 'P_PTR', IERR)
      SCF%P_PTR = 0.0_EB
   
      !!! allocate restriction matrix including row and column pointers
      ALLOCATE (SCF%R(SCF%NR), STAT=IERR)
      CALL CHKMEMERR ('SCARC_SETUP_COARSENING', 'R', IERR)
      SCF%R = 0.0_EB
   
      ALLOCATE (SCF%R_ROW(2*SCF%NCCE-SCF%NCC), STAT=IERR)
      CALL CHKMEMERR ('SCARC_SETUP_COARSENING', 'R_ROW', IERR)
      SCF%R_ROW = 0.0_EB
   
      ALLOCATE (SCF%R_COL(SCF%NR), STAT=IERR)
      CALL CHKMEMERR ('SCARC_SETUP_COARSENING', 'R_COL', IERR)
      SCF%R_COL = 0.0_EB
   
      NLEN(1) = SCF%NX
      NLEN(2) = SCF%NY
      NLEN(3) = SCF%NZ
      ILEN = MAXVAL(NLEN)
      ALLOCATE (SCF%WC_INDICES(-3:3, ILEN), STAT=IERR)
      CALL CHKMEMERR ('SCARC_SETUP_COARSENING', 'WC_INDICES', IERR)
      SCF%WC_INDICES = 0

WRITE(SCARC_LU,*) 'HALLO2: NCCE=',SCF%NCCE, NM
   ENDDO


   CALL SCARC_SETUP_PROLONGATION(NL)

   !CALL SCARC_SETUP_PROLONGATION2(NL)

   CALL SCARC_SETUP_RESTRICTION(NL)

   IF (NL>=1) THEN
   !CALL PRINT_MATRIX(SCF%P, SCF%P_ROW, SCF%P_COL, SCF%NC, SCF%NCC, NM, NL, 'P')
   !CALL PRINT_MATRIX(SCF%R, SCF%R_ROW, SCF%R_COL, SCF%NCC, SCF%NC, NM, NL, 'R')
   CALL PRINT_MATRIX2(SCF%P, SCF%P_ROW, SCF%P_COL, SCF%NC, SCF%NCC, NM, NL, 'P')
   CALL PRINT_MATRIX2(SCF%R, SCF%R_ROW, SCF%R_COL, SCF%NCC, SCF%NC, NM, NL, 'R')
   ENDIF

   !!!-------------------------------------------------------------------------------------------------
   !!! Allocate coarse grid matrix including pointer arrays
   !!! Note: number of cells on coarse level corresponds to number of c-points on fine level
   !!! Compute coarse grid matrix by multiplication with restriction and prolongation matrix:
   !!!  A_coarse := R * A_fine * P
   !!!-------------------------------------------------------------------------------------------------
   DO NM = NMESHES_MIN, NMESHES_MAX

      SCF => SCARC(NM)%COMPACT(NL)                 ! Pointer to fine level
      SCC => SCARC(NM)%COMPACT(NL+1)               ! Pointer to coarse level

      SCC%NC  = SCF%NCC
      SCC%NCE = SCF%NCCE
      SCC%NW  = SCF%NCW
   
WRITE(SCARC_LU,*) 'NC on level ', NL+1,' = ', SCC%NC
WRITE(SCARC_LU,*) 'NW on level ', NL+1,' = ', SCC%NW
WRITE(SCARC_LU,*) 'HALLO4: NCCE=',SCF%NCCE, NM

   ENDDO

   !!! get neighboring NC's and NW's
   CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_AMG, NL+1)

   DO NM = NMESHES_MIN, NMESHES_MAX

      SCF => SCARC(NM)%COMPACT(NL)                 ! Pointer to fine level
      SCC => SCARC(NM)%COMPACT(NL+1)               ! Pointer to coarse level

WRITE(SCARC_LU,*) 'HALLO5: NCCE=',SCF%NCCE, NM, SCC%NW
      !!! Allocate coarse matrix structures
      ALLOCATE (SCC%A(SCF%NCCE*10), STAT=IERR)
      CALL CHKMEMERR ('SCARC_SETUP_COARSENING', 'SCC%A', IERR)
      SCC%A = 0.0_EB
   
      ALLOCATE (SCC%A_ROW(SCF%NCCE+1), STAT=IERR)
      CALL CHKMEMERR ('SCARC_SETUP_COARSENING', 'SCC%A_ROW', IERR)
      SCC%A_ROW = 0
   
      ALLOCATE (SCC%A_COL(SCF%NCCE*10), STAT=IERR)
      CALL CHKMEMERR ('SCARC_SETUP_COARSENING', 'SCC%A_COL', IERR)
      SCC%A_COL = 0

      !!! Allocate IJKW-information for coarser grids (own and neighboring grids)
      ALLOCATE(SCC%IJKW(15,SCC%NW), STAT=IERR)
      CALL ChkMemErr('SCARC_SETUP_COARSENING','IJKW',IERR)
      SCC%IJKW = 0

      ALLOCATE(SCC%ADJACENT_CELL(SCC%NW), STAT=IERR)
      CALL ChkMemErr('SCARC_SETUP_COARSENING','ADJACENT_CELL',IERR)
      SCC%ADJACENT_CELL = 0

      ALLOCATE(SCC%GHOST_CELL(SCC%NW), STAT=IERR)
      CALL ChkMemErr('SCARC_SETUP_COARSENING','GHOST_CELL',IERR)
      SCC%GHOST_CELL = 0

      ALLOCATE(SCC%GHOST_PTR(SCC%NC+1:SCC%NCE), STAT=IERR)
      CALL ChkMemErr('SCARC_SETUP_COARSENING','GHOST_PTR',IERR)
      SCC%GHOST_PTR = 0
  
WRITE(SCARC_LU,*) 'SCF%NCE=',SCF%NCE, SCF%NC
WRITE(SCARC_LU,*) 'HALLO6: NCCE=',SCF%NCCE, NM
      DO IC = SCF%NC+1,SCF%NCE
WRITE(SCARC_LU,*) 'IC=',IC
         ICC = SCF%CELLTYPE(IC)
WRITE(SCARC_LU,*) 'IC=',IC,': ICC=',ICC
!         IF (ICC >= NSCARC_CELLTYPE_COARSE) SCC%GHOST_PTR(ICC)=IC
!         IF (ICC >= NSCARC_CELLTYPE_COARSE) WRITE(SCARC_LU,*) 'SCC%GHOST_PTR(',ICC,')=',SCC%GHOST_PTR(ICC)
      ENDDO
WRITE(SCARC_LU,*) 'READY'

      OTHER_MESHES_LOOP: DO NOM = 1, NMESHES

         IF (NOM == NM) CYCLE OTHER_MESHES_LOOP

         OSCC => SCARC(NM)%OSCARC(NOM)%COMPACT(NL+1)                     

WRITE(SCARC_LU,*) 'OSCC%NC=',OSCC%NC
WRITE(SCARC_LU,*) 'OSCC%NW=',OSCC%NW

         ALLOCATE(OSCC%IJKW(15,OSCC%NW), STAT=IERR)
         CALL ChkMemErr('SCARC_SETUP_EXCHANGE','IJKW',IERR)
         OSCC%IJKW = 0

         ALLOCATE(OSCC%ADJACENT_CELL(OSCC%NW), STAT=IERR)
         CALL ChkMemErr('SCARC_SETUP_EXCHANGE','ADJACENT_CELL',IERR)
         OSCC%ADJACENT_CELL = 0

         ALLOCATE(OSCC%GHOST_CELL(OSCC%NW), STAT=IERR)
         CALL ChkMemErr('SCARC_SETUP_EXCHANGE','GHOST_CELL',IERR)
         OSCC%GHOST_CELL = 0

         ALLOCATE(OSCC%GHOST_CELL(OSCC%NC:OSCC%NCE), STAT=IERR)
         CALL ChkMemErr('SCARC_SETUP_EXCHANGE','GHOST_CELL',IERR)
         OSCC%GHOST_CELL = 0

      ENDDO OTHER_MESHES_LOOP

   ENDDO

   !CALL SCARC_SETUP_SUBDIVISION_AMG (NL)
   IF (NMESHES > 1) CALL SCARC_SETUP_WALLCELLS2_AMG (NL)
   CALL SCARC_TRANSFER_MATRIX (NL)

IF (TYPE_DEBUG >= NSCARC_DEBUG_MEDIUM) WRITE(SCARC_LU,*) 'CELLTYPE FINAL AFTER EXCHANGE'
   CALL SCARC_DEBUG_QUANTITY(NSCARC_DEBUG_CELLTYPE, NL, 'SETUP_CELLTYPES', 'CELLTYPE AFTER EXCHANGE')

IF (TYPE_DEBUG >= NSCARC_DEBUG_MEDIUM) WRITE(SCARC_LU,*) 'SIZE MATRIX on level ', NL+1,'= ',SCC%NC

   !IF (SCC%NC<5) EXIT LEVEL_LOOP

IF (TYPE_DEBUG >= NSCARC_DEBUG_MEDIUM) WRITE(SCARC_LU,*) '=============== FINISHING LEVEL ', NL
ENDDO LEVEL_LOOP

DO NL = NLEVEL_MIN, NLEVEL_MAX
   DO NM = NMESHES_MIN, NMESHES_MAX
      SCF => SCARC(NM)%COMPACT(NL)  
      !CALL PRINT_MATRIX(SCF%A, SCF%A_ROW, SCF%A_COL, SCF%NC, SCF%NC, NM, NL, 'A')
      CALL PRINT_MATRIX2(SCF%A, SCF%A_ROW, SCF%A_COL, SCF%NC, SCF%NC, NM, NL, 'A')
      !DEALLOCATE(SCF%MEASURE)
      !DEALLOCATE(SCF%CELLTYPE)
   ENDDO
ENDDO

END SUBROUTINE SCARC_SETUP_COARSENING



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Store subdivision information
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SETUP_SUBDIVISION_AMG(NL)
INTEGER, INTENT(IN) :: NL
INTEGER :: NM, IW, IOR0, IOR_LAST, NOM, INBR
INTEGER :: NEIGHBORS(20,-3:3)
TYPE (SCARC_COMPACT_TYPE), POINTER :: SCF, SCC

IOR_LAST    = 0
NEIGHBORS   = 0
         
MESHES_LOOP1: DO NM = NMESHES_MIN, NMESHES_MAX
   
   SCF => SCARC(NM)%COMPACT(NL)
   SCC => SCARC(NM)%COMPACT(NL+1)

   SCC%SUBDIVISION = 0

   WALLCELL_LOOP: DO IW = 1, SCC%NW
   
      IOR0 = SCF%IJKW(4,IW)
      
      IF (IOR_LAST /= IOR0) SCC%SUBDIVISION(1,IOR0) = IW
      SCC%SUBDIVISION(2,IOR0) = SCC%SUBDIVISION(2,IOR0) + 1
      
      NOM= SCF%IJKW(9,IW)
   !WRITE(SCARC_LU,'(a,i3,a,i3,a,i3,a,3i3)') '    NM=',NM,': NOM=',NOM,': SUBDIVISION(2,',IOR0,')=',SUBDIVISION(2,IOR0), IOR0, SUBDIVISION(2,-2)
   
      IF (NOM /= 0) THEN
         NEIGHBOR_LOOP: DO INBR = 1, 20
   !WRITE(SCARC_LU,'(10i3)') NEIGHBORS(1:10,IOR0)
            IF (NOM == NEIGHBORS(INBR, IOR0)) THEN
               EXIT NEIGHBOR_LOOP
            ELSE IF (NEIGHBORS(INBR, IOR0) /= 0) THEN
               CYCLE NEIGHBOR_LOOP
            ELSE IF (NEIGHBORS(INBR, IOR0) == 0) THEN
               NEIGHBORS(INBR, IOR0) = NOM
               SCC%SUBDIVISION(3,IOR0) = SCC%SUBDIVISION(3,IOR0) + 1
   !WRITE(SCARC_LU,'(a,i3,a,i3,a,i3,a,3i3)') '--->NM=',NM,': NOM=',NOM,': SUBDIVISION(3,',IOR0,')=',SUBDIVISION(3,IOR0),IOR0,SUBDIVISION(2,-2)
               EXIT NEIGHBOR_LOOP
            ELSE
               WRITE(*,*) 'More than 20 neighbors at one face not allowed yet!'
               STOP
            ENDIF
         ENDDO NEIGHBOR_LOOP
      ENDIF
      
      IOR_LAST = IOR0
   
   ENDDO WALLCELL_LOOP

ENDDO MESHES_LOOP1
         
!WRITE(SCARC_LU,'(a,7i5)') 'SUBDIVISION(1,:)',SUBDIVISION(1,-3:3)
!WRITE(SCARC_LU,'(a,7i5)') 'SUBDIVISION(2,:)',SUBDIVISION(2,-3:3)
!WRITE(SCARC_LU,'(a,7i5)') 'SUBDIVISION(3,:)',SUBDIVISION(3,-3:3)
END SUBROUTINE SCARC_SETUP_SUBDIVISION_AMG



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Determine if a cell IC is strongly coupled to another cell JC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SETUP_WALLCELLS2_AMG(NL)
INTEGER , INTENT(IN) :: NL
INTEGER :: NM, NOM, ICPL, ICO, JCC, ICOL, ICOL0, IOR0
INTEGER :: II, JJ, KK, INBR, NNUM, IPTR, ITYPE, ICF, ICA, ICC, ICS, IC2, ICG, IWF, IWC, ICN
LOGICAL :: COARSE_BDRY_CELL
TYPE (SCARC_COMPACT_TYPE), POINTER :: SCF, SCC
TYPE (OSCARC_COMPACT_TYPE), POINTER :: OSCF, OSCC


IF (TYPE_DEBUG >= NSCARC_DEBUG_LESS) WRITE(SCARC_LU,*) '=============== STARTING WALLCELLS2', NL
CALL SCARC_DEBUG_QUANTITY(NSCARC_DEBUG_CELLTYPE, NL, 'SETUP_CELLTYPES', 'CELLTYPE AFTER EXCHANGE')

SELECT CASE (NL)

   !!! -----------------------------------------------------------------------------------------
   !!! Take original IJKW from finest level
   !!! -----------------------------------------------------------------------------------------
   CASE (1)

      WRITE(SCARC_LU,*) '================== NL=',NL, NMESHES_MIN, NMESHES_MAX
      MESHES_LOOP1: DO NM = NMESHES_MIN, NMESHES_MAX
   
        WRITE(SCARC_LU,*) '================== NM=',NM

         SCF => SCARC(NM)%COMPACT(NL)
         SCC => SCARC(NM)%COMPACT(NL+1)

WRITE(SCARC_LU,*) 'ADJACENT_CELL:'
WRITE(SCARC_LU,*) (SCF%ADJACENT_CELL(IWF), IWF=1, SCF%NW)
         IWC = 0
         !!! -----------------------------------------------------------------------------
         !!! Loop over all boundary cells IW of fine grid
         !!! get corresponding adjacent and ghost cell
         !!! Store the following information in IJKW for coarser mesh
         !!! IJKW(1 ,.) : new coarse number of cell
         !!! IJKW(2 ,.) : corresponding IOR
         !!! IJKW(2 ,.) : neighboring NOM
         !!! IJKW(4 ,.) : number of connections in neighboring mesh
         !!! IJKW(5+,.) : coarse cell numbers of connected cells in neighboring mesh
         !!! -----------------------------------------------------------------------------
         DO IWF = 1, SCF%NW
            IF (SCF%IJKW(9,IWF)==0) CYCLE
            ICA = SCF%ADJACENT_CELL(IWF)
            ICG = SCF%GHOST_CELL(IWF)
WRITE(SCARC_LU,*) '===================================================================='
WRITE(SCARC_LU,*) 'IWF=',IWF,': ICA=',ICA, ICG
WRITE(SCARC_LU,*) '       ---->', SCF%P_ROW(ICG), SCF%P_ROW(ICG+1),SCF%P_ROW(ICG+1)-SCF%P_ROW(ICG)

            !!! 
            !!! If adjacent cell is a coarse cell, then get it's coarse neighbors in neighboring mesh
            !!! directly from relations in prolongation matrix
            !!! 
            IF (SCF%CELLTYPE(ICA) >= NSCARC_CELLTYPE_COARSE) THEN
               IWC = IWC + 1
               SCC%IJKW(1,IWC) = SCF%CELLTYPE(ICA)
               SCC%IJKW(2,IWC) = SCF%IJKW(4,IWF)
               SCC%IJKW(3,IWC) = SCF%IJKW(9,IWF)
               IF (ICG>0) THEN
                  IPTR=0
                  DO ICOL = SCF%P_ROW(ICG), SCF%P_ROW(ICG+1)-1
                    IF (SCF%P_COL(ICOL)<0) THEN
                       SCC%IJKW(5+IPTR,IWC) = ABS(SCF%P_COL(ICOL))
                       IPTR = IPTR + 1
WRITE(SCARC_LU,*) 'SUSI: IJKW(',5+IPTR-1,',',IWC,')=',ABS(SCF%P_COL(ICOL)), IPTR
                    ENDIF
                  ENDDO
                  SCC%IJKW(4,IWC) = IPTR
                  IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) &
                     WRITE(SCARC_LU,'(a,i3,a,i3,a,i3,a,i3,15i4)') &
                        'A:IWF=',IWF,' IWC=',IWC,' ICA=',ICA,' ICG=',ICG, (SCC%IJKW(II,IWC), II=1,5+IPTR-1)
               ENDIF
            ELSE
               COARSE_BDRY_CELL = .FALSE.
               PROLONGATION_LOOP1: DO ICOL = SCF%P_ROW(ICA), SCF%P_ROW(ICA+1)-1
                  IF (SCF%P_COL(ICOL)<0) THEN
                     IWC = IWC + 1
                     COARSE_BDRY_CELL = .TRUE.
                     EXIT PROLONGATION_LOOP1
                  ENDIF
               ENDDO PROLONGATION_LOOP1
               ICOL0 = -1
               IF (COARSE_BDRY_CELL) THEN
                  SCC%IJKW(2,IWC) = SCF%IJKW(4,IWF)
                  SCC%IJKW(3,IWC) = SCF%IJKW(9,IWF)
                  DO ICOL = SCF%A_ROW(ICA)+1, SCF%A_ROW(ICA+1)-1
                     IF (SCF%A_COL(ICOL) == SECOND_LINE_NBR1(IWF,SCF%IJKW,SCF%NX,SCF%NY)) ICOL0 = ICOL
                  ENDDO
                  IPTR=0
                  DO ICOL = SCF%P_ROW(ICA), SCF%P_ROW(ICA+1)-1
                     ICN = SCF%P_COL(ICOL)
                     IF (ICN < 0) THEN
                        SCC%IJKW(5+IPTR,IWC) = ABS(ICN)
                        IPTR = IPTR + 1
                     ENDIF
    WRITE(SCARC_LU,'(a,i3,a,i3,a,i3,a,i3,a,i3,a,i3,a,i3)') &
                    'B: ICOL=',ICOL,': ICOL0=',ICOL0, 'ICN=',ICN,': IWF=',IWF,': ICA=',ICA,': ICOL0=',ICOL0
                  ENDDO
                  IF (ICOL0 > 0) THEN
                     SCC%IJKW(1,IWC) = SCF%CELLTYPE(SCF%A_COL(ICOL0))
                  ELSE
                     WRITE(*,*) 'Error during initialization of neighboring prolongation cells, stopping!'
                     STOP
                  ENDIF
                  SCC%IJKW(4,IWC) = IPTR
                  IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) &
                     WRITE(SCARC_LU,'(a,i3,a,i3,a,i3,a,i3,15i4)') &
                        'C: IWF=',IWF,' IWC=',IWC,' ICA=',ICA,' ICG=',ICG, (SCC%IJKW(II,IWC), II=1,5+IPTR-1)
               ENDIF
            ENDIF
         ENDDO


      WRITE(SCARC_LU,*) '================== WALLCELLS-FINAL:'
      DO IPTR=1, IWC
      WRITE(SCARC_LU,'(15i3)') (SCC%IJKW(IWF, IPTR), IWF=1,15)
      ENDDO

      ENDDO MESHES_LOOP1

END SELECT

END SUBROUTINE SCARC_SETUP_WALLCELLS2_AMG


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Determine if a cell IC is strongly coupled to another cell JC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SETUP_WALLCELLS_AMG(NL)
INTEGER , INTENT(IN) :: NL
INTEGER :: NM, NOM, ICPL, IWC, IWF, ICO, ICG, ICF, ICC, JCC, ICOL, IOR0, II, JJ, KK, INBR
TYPE (SCARC_COMPACT_TYPE), POINTER :: SCF, SCC
TYPE (OSCARC_COMPACT_TYPE), POINTER :: OSCF, OSCC

SELECT CASE (NL)

   !!! -----------------------------------------------------------------------------------------
   !!! Take original IJKW from finest level
   !!! -----------------------------------------------------------------------------------------
   CASE (1)

      MESHES_LOOP1: DO NM = NMESHES_MIN, NMESHES_MAX
   
         SCF => SCARC(NM)%COMPACT(NL)
         SCC => SCARC(NM)%COMPACT(NL+1)

!WRITE(SCARC_LU,*) '================== NM=',NM, SCF%NW, SCC%NW
         OTHER_MESHES_LOOP1: DO NOM = NMESHES_MIN, NMESHES_MAX

            IF (NOM == NM) CYCLE OTHER_MESHES_LOOP1

            OSCF => SCARC(NM)%OSCARC(NOM)%COMPACT(NL)                     
            OSCC => SCARC(NM)%OSCARC(NOM)%COMPACT(NL+1)                     

            IWC = 0
            DO IWF = 1, SCF%NW

               IF (SCF%IJKW(9,IWF) /= NOM) CYCLE

               ICF = SCF%ADJACENT_CELL(IWF)
               ICG = SCF%GHOST_CELL(IWF)
!WRITE(SCARC_LU,*) 'adjacent cell:', ICF, ' ghost cell=',ICG

            
!WRITE(SCARC_LU,*) '----->  IWF=',IWF, OSCF%NX, OSCF%NY
!WRITE(SCARC_LU,'(A,15i3)') '----->  IJKW:',OSCF%IJKW(1:15,IWF)
               !!! Determine the neighboring fine cells to ICF
               INBR = 0
               DO KK=OSCF%IJKW(12,IWF),OSCF%IJKW(15,IWF)
                  DO JJ=OSCF%IJKW(11,IWF),OSCF%IJKW(14,IWF)
                     DO II=OSCF%IJKW(10,IWF),OSCF%IJKW(13,IWF)
                        INBR = INBR+1
                !        ICF = (KK-1) * SCF%NX*SCF%NY + (JJ-1)*SCF%NX + II
                     ENDDO
                  ENDDO
               ENDDO
               ICO = (OSCF%IJKW(8,IWF)-1) * OSCF%NX*OSCF%NY + (OSCF%IJKW(7,IWF)-1)*OSCF%NX + OSCF%IJKW(6,IWF)
WRITE(SCARC_LU,*) 'Own IC:', ICF, ' NBR IC=',ICO

               !!! ONLY TEMPORARILY: Only allow same resolution on neighboring grid!
               IF (INBR > 1) THEN
                  WRITE(SCARC_LU,*) 'Error in SCARC_SETUP_WALLCELLS_AMG: different resolutions on neighboring grids! Stop!'
                  STOP
               ENDIF

               ICC = SCF%CELLTYPE(ICF)
               IF (ICC < NSCARC_CELLTYPE_COARSE) CYCLE

               !!! save orientation and neighbor information from 'finer' IJKW 
               IOR0 = SCF%IJKW(4, IWF)

               IWC = IWC + 1
               SCC%IJKW(4, IWC) = IOR0
               SCC%IJKW(9, IWC) = NOM
               SCC%IJKW(1, IWC) = ICC

WRITE(SCARC_LU,*) 'ICC=',ICC

               !!! analyze its neighbors on new coarse grid
               ICPL = 0
               DO ICOL = SCC%A_ROW(ICC)+1, SCC%A_ROW(ICC+1)-1
                  JCC = SCC%A_COL(ICOL)
                  IF (JCC > SCC%NC) THEN
                     SCC%IJKW( 5+ICPL,IWC) = JCC
                     SCC%IJKW(10+ICPL,IWC) = SCF%CELLTYPE(SCC%GHOST_PTR(JCC))
                     ICPL = ICPL + 1
                     IF (ICPL > 4) THEN
                        WRITE(*,*) 'Error in SCARC_SETUP_WALLCELLS_AMG: Too many couplings ',ICPL,' for IC=',ICC,' on level ',NL
                        STOP
                     ENDIF
                  ENDIF
               ENDDO
               SCC%IJKW(2, IWC) = ICPL
                  
WRITE(SCARC_LU,'(a,i3,a,i3,a,15i3)') 'ICC=',ICC,': IJKW(1:15,',IWC,')=',SCC%IJKW(1:15,IWC)
            ENDDO

         ENDDO OTHER_MESHES_LOOP1

      ENDDO MESHES_LOOP1

   !!! -----------------------------------------------------------------------------------------
   !!! Take 'own' IJKW's from coarser level (note that they have another structure)
   !!! -----------------------------------------------------------------------------------------
   CASE DEFAULT

      MESHES_LOOP2: DO NM = NMESHES_MIN, NMESHES_MAX
   

      ENDDO MESHES_LOOP2

END SELECT

   
END SUBROUTINE SCARC_SETUP_WALLCELLS_AMG


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Determine if a cell IC is strongly coupled to another cell JC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE PRINT_MATRIX(A, A_ROW, A_COL, NC1, NC2, NM, NL, CNAME)
REAL(EB), DIMENSION(:), INTENT(IN) :: A
INTEGER , DIMENSION(:), INTENT(IN) :: A_ROW
INTEGER , DIMENSION(:), INTENT(IN) :: A_COL
INTEGER , INTENT(IN) :: NM, NL, NC1, NC2
CHARACTER(1), INTENT(IN):: CNAME
INTEGER :: IC, JC, ICOL, MMATRIX
CHARACTER(60):: CMATRIX
REAL(EB):: MATRIX_LINE(1000)

WRITE (CMATRIX, '(A,A1,A,i2.2,A,i2.2,A)') 'matrix/',CNAME,'_mesh',NM,'_level',NL,'.txt'
MMATRIX=GET_FILE_NUMBER()
OPEN(MMATRIX,FILE=CMATRIX)
      
DO IC = 1, NC1
   MATRIX_LINE=0.0_EB
   DO JC = 1, NC2
      DO  ICOL= A_ROW(IC), A_ROW(IC+1)-1       
         IF (A_COL(ICOL)==JC) MATRIX_LINE(JC)=A(ICOL)
      ENDDO
   ENDDO
   WRITE(MMATRIX,1000) (MATRIX_LINE(JC),JC=1,NC2)
ENDDO

CLOSE(MMATRIX)

1000 FORMAT(F11.6,',',F11.6,',',F11.6,',',F11.6,',',F11.6,',',F11.6,',',F11.6,',',&
            F11.6,',',F11.6,',',F11.6,',',F11.6,',',F11.6,',',F11.6,',', &
            F11.6,',',F11.6,',',F11.6,',',F11.6,',',F11.6,',',F11.6,',', &
            F11.6,',',F11.6,',',F11.6,',',F11.6,',',F11.6,',',F11.6,',', &
            F11.6,',',F11.6,',',F11.6,',',F11.6,',',F11.6,',',F11.6,',', &
            F11.6,',',F11.6,',',F11.6,',',F11.6,',',F11.6,',',F11.6,',', &
            F11.6,',',F11.6,',',F11.6,',',F11.6,',',F11.6,',',F11.6,',', &
            F11.6,',',F11.6,',',F11.6,',',F11.6,',',F11.6,',',F11.6,';')
END SUBROUTINE PRINT_MATRIX

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Determine if a cell IC is strongly coupled to another cell JC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE PRINT_MATRIX2(A, A_ROW, A_COL, NC1, NC2, NM, NL, CNAME)
REAL(EB), DIMENSION(:), INTENT(IN) :: A
INTEGER , DIMENSION(:), INTENT(IN) :: A_ROW
INTEGER , DIMENSION(:), INTENT(IN) :: A_COL
INTEGER , INTENT(IN) :: NM, NL, NC1, NC2
CHARACTER(1), INTENT(IN):: CNAME
INTEGER :: IC, JC, ICOL, MMATRIX
CHARACTER(60):: CMATRIX
REAL(EB):: MATRIX_LINE(1000)

WRITE (CMATRIX, '(A,A1,A,i2.2,A,i2.2,A)') 'matrix2/',CNAME,'_mesh',NM,'_level',NL,'.txt'
MMATRIX=GET_FILE_NUMBER()
OPEN(MMATRIX,FILE=CMATRIX)
      
DO IC = 1, NC1
   MATRIX_LINE=0.0_EB
   DO JC = 1, NC2
      DO  ICOL= A_ROW(IC), A_ROW(IC+1)-1       
         IF (A_COL(ICOL)==JC) MATRIX_LINE(JC)=A(ICOL)
      ENDDO
   ENDDO
   WRITE(MMATRIX,1000) (MATRIX_LINE(JC),JC=1,NC2)
ENDDO

CLOSE(MMATRIX)

1000 FORMAT(40f8.2)
END SUBROUTINE PRINT_MATRIX2


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Determine measure of cells corresponding to requested coarsening type
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Determine if a cell IC is strongly coupled to another cell JC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
LOGICAL FUNCTION STRONGLY_COUPLED(A, A_ROW, IC, ICOL)
REAL(EB), DIMENSION(:), INTENT(IN) :: A
INTEGER , DIMENSION(:), INTENT(IN) :: A_ROW
INTEGER , INTENT(IN) :: IC, ICOL
INTEGER :: JCOL
REAL(EB) :: AMG_TOL, VAL_MAX

AMG_TOL = 0.25_EB
VAL_MAX = 0.00_EB

DO  JCOL= A_ROW(IC)+1, A_ROW(IC+1)-1       
   IF (JCOL /= ICOL) VAL_MAX = MAX(VAL_MAX, A(JCOL))
ENDDO

IF (A(ICOL) >= AMG_TOL * VAL_MAX) THEN
   STRONGLY_COUPLED=.TRUE.
ELSE
   STRONGLY_COUPLED=.FALSE.
ENDIF
RETURN

END FUNCTION STRONGLY_COUPLED


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Determine measure of cells corresponding to requested coarsening type
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SETUP_MEASURES(NL)
INTEGER, PARAMETER  :: NCOUPLINGS = 20
INTEGER, INTENT(IN) :: NL
INTEGER  :: NM, NTYPE, IC, JC, KC, ICOL, JCOL, IS, IW, JW, ICW
INTEGER  :: NX1, NX2, NY1, NY2, NZ1, NZ2, IX, IY, IZ
INTEGER  :: SCOUPLED(NCOUPLINGS), WCOUPLED(NCOUPLINGS)
REAL(EB) :: RAND_NUM, RAND_NUM2
LOGICAL  :: BFIVE, BSFINE, BXFIRST, BZFIRST
TYPE (SCARC_COMPACT_TYPE), POINTER :: SC, SCC

!!! Use special coarsening strategies only on finest level, on coarser levels always use RS3-coarsening
IF (NL> NLEVEL_MIN) THEN
   NTYPE = NSCARC_COARSENING_TEST
ELSE
   NTYPE = TYPE_COARSENING
ENDIF

IF (TYPE_DEBUG >= NSCARC_DEBUG_MEDIUM) WRITE(SCARC_LU,*) '=========== DEFINING MEASURE on level ', NL
!!! Select coarsening strategy
SELECT CASE (NTYPE)

   !!!-------------------------------------------------------------------------------------------------
   !!! TEST
   !!!-------------------------------------------------------------------------------------------------
   CASE (NSCARC_COARSENING_TEST)


      TEST_LOOP1: DO NM = NMESHES_MIN, NMESHES_MAX

         SC => SCARC(NM)%COMPACT(NL)

         TEST_MEASURE_LOOP1: DO IC = 1, SC%NC
            DO ICOL = SC%A_ROW(IC)+1, SC%A_ROW(IC+1)-1
               IF (STRONGLY_COUPLED(SC%A, SC%A_ROW, IC, ICOL)) SC%MEASURE(IC) = SC%MEASURE(IC) + 1.0_EB
            ENDDO
         ENDDO TEST_MEASURE_LOOP1

IF (TYPE_DEBUG >= NSCARC_DEBUG_MEDIUM) THEN
   WRITE(SCARC_LU,*) '1:MEASURES in SETUP_MEASURE_TEST: =====  NM=',NM
   IF (NM==1) THEN
      WRITE(SCARC_LU,'(9f8.2)') (SC%MEASURE(IC), IC=1,SC%NC)
   ELSE
      WRITE(SCARC_LU,'(8f8.2)') (SC%MEASURE(IC), IC=1,SC%NC)
   ENDIF
ENDIF
      ENDDO TEST_LOOP1

      IF (NMESHES > 1) THEN

         TEST_LOOP3: DO NM = NMESHES_MIN, NMESHES_MAX

            SC => SCARC(NM)%COMPACT(NL)
            TEST_MEASURE_LOOP3: DO IW = 1, SC%NW
               IF (SC%GHOST_CELL(IW) == 0) CYCLE TEST_MEASURE_LOOP3
               IC = SC%ADJACENT_CELL(IW)
               SC%MEASURE(IC)=0.0_EB
               DO ICOL = SC%A_ROW(IC)+1, SC%A_ROW(IC+1)-1
                  IF (STRONGLY_COUPLED(SC%A, SC%A_ROW, IC, ICOL)) SC%MEASURE(IC) = SC%MEASURE(IC) + 1.0_EB
               ENDDO
               CALL RANDOM_NUMBER(RAND_NUM)
               SC%MEASURE(IC) = SC%MEASURE(IC) + 0.1_EB + REAL(INT(RAND_NUM*8.9_EB),EB)/10.0_EB

            ENDDO TEST_MEASURE_LOOP3

IF (TYPE_DEBUG >= NSCARC_DEBUG_MEDIUM) THEN
   WRITE(SCARC_LU,*) '2:MEASURES in SETUP_MEASURE_TEST: =====  NM=',NM
   IF (NM==1) THEN
      WRITE(SCARC_LU,'(9f8.2)') (SC%MEASURE(IC), IC=1,SC%NC)
   ELSE
      WRITE(SCARC_LU,'(8f8.2)') (SC%MEASURE(IC), IC=1,SC%NC)
   ENDIF
ENDIF

         ENDDO TEST_LOOP3

      ENDIF


   !!!-------------------------------------------------------------------------------------------------
   !!! RS3-coarsening: 
   !!!      - Original Ruge-Stuben method with parallel postprocessing
   !!!      - Produces good C/F splittings but is inherently serial.  
   !!!      - May produce AMG hierarchies with relatively high operator complexities.
   !!!-------------------------------------------------------------------------------------------------
   CASE (NSCARC_COARSENING_RS3)

      RS3_LOOP1: DO NM = NMESHES_MIN, NMESHES_MAX
      
         SC => SCARC(NM)%COMPACT(NL)            
      
         RS3_MEASURE_LOOP1: DO IC = 1, SC%NC
            DO ICOL = SC%A_ROW(IC)+1, SC%A_ROW(IC+1)-1
               IF (STRONGLY_COUPLED(SC%A, SC%A_ROW, IC, ICOL)) SC%MEASURE(IC) = SC%MEASURE(IC) + 1.0_EB
            ENDDO
         ENDDO RS3_MEASURE_LOOP1
      
IF (TYPE_DEBUG >= NSCARC_DEBUG_MEDIUM) THEN
   WRITE(SCARC_LU,*) '1:MEASURES in SETUP_MEASURE_RS3: =====  NM=',NM
   IF (NM==1) THEN
      WRITE(SCARC_LU,'(9f8.2)') (SC%MEASURE(IC), IC=1,SC%NC)
   ELSE
      WRITE(SCARC_LU,'(8f8.2)') (SC%MEASURE(IC), IC=1,SC%NC)
   ENDIF
ENDIF
      ENDDO RS3_LOOP1

      IF (NMESHES > 1) THEN

         RS3_LOOP3: DO NM = NMESHES_MIN, NMESHES_MAX

            SC => SCARC(NM)%COMPACT(NL)            
            RS3_MEASURE_LOOP3: DO IW = 1, SC%NW

               IF (NL == 1) THEN
                  !IF (SC%GHOST_CELL(IW) == 0) CYCLE RS3_MEASURE_LOOP3
                  IF (SC%IJKW(9,IW) == 0) CYCLE RS3_MEASURE_LOOP3
                  IC = SC%ADJACENT_CELL(IW)
               ELSE
                  IC = SC%IJKW(1,IW)
               ENDIF
               SC%MEASURE(IC)=0.0_EB
               DO ICOL = SC%A_ROW(IC)+1, SC%A_ROW(IC+1)-1
                  IF (STRONGLY_COUPLED(SC%A, SC%A_ROW, IC, ICOL)) SC%MEASURE(IC) = SC%MEASURE(IC) + 1.0_EB
               ENDDO
               CALL RANDOM_NUMBER(RAND_NUM)
               CALL RANDOM_NUMBER(RAND_NUM2)
               RAND_NUM = 0.99_EB* RAND_NUM + 0.099_EB * RAND_NUM2
               !SC%MEASURE(IC) = SC%MEASURE(IC) + 0.1_EB + REAL(INT(RAND_NUM*8.9_EB),EB)/10.0_EB
               SC%MEASURE(IC) = SC%MEASURE(IC) + RAND_NUM
           
            ENDDO RS3_MEASURE_LOOP3

IF (TYPE_DEBUG >= NSCARC_DEBUG_MEDIUM) THEN
   WRITE(SCARC_LU,*) '2:MEASURES in SETUP_MEASURE_RS3: =====  NM=',NM
   IF (NM==1) THEN
      WRITE(SCARC_LU,'(8f8.2)') (SC%MEASURE(IC), IC=1,SC%NC)
   ELSE
      WRITE(SCARC_LU,'(8f8.2)') (SC%MEASURE(IC), IC=1,SC%NC)
   ENDIF
ENDIF
   
         ENDDO RS3_LOOP3
   
      ENDIF

   !!!-------------------------------------------------------------------------------------------------
   !!! FDSA2-coarsening:
   !!!     -  path l=2, length l=2  ==> S_i^(2,2) with parallel postprocessing
   !!!-------------------------------------------------------------------------------------------------
   CASE (NSCARC_COARSENING_FDSA2)

      FDSA2_LOOP: DO NM = NMESHES_MIN, NMESHES_MAX

         SC  => SCARC(NM)%COMPACT(NL)
         SCC => SCARC(NM)%COMPACT(NL+1)

         IF (NMESHES >= 1) THEN

            !!! Set internal measures corresponding to A2-coarsening
            SELECT CASE (TYPE_DIMENSION)
   
               !!! -------------------------------- 2D ------------------------------------------------
               CASE (NSCARC_DIMENSION_TWO)
 
                  SCC%NX=0
                  SCC%NY=1
                  SCC%NZ=0
                  BZFIRST=.TRUE.
                  DO IZ = 1, SC%NZ
                     BXFIRST=.TRUE.
                     DO IX = 1, SC%NX
                        IC = (IZ-1)*SC%NX + IX
                        IF (MOD(IX,2)==1.AND.MOD(IZ,2)==1) THEN
                           SC%MEASURE(IC) = 6.0_EB
                           IF (BZFIRST) SCC%NX=SCC%NX+1
                           IF (BXFIRST) SCC%NZ=SCC%NZ+1
                        ELSE IF (MOD(IX,2)==0.AND.MOD(IZ,2)==0) THEN
                           SC%MEASURE(IC) = 4.0_EB
                        ELSE
                           SC%MEASURE(IC) = 5.0_EB
                        ENDIF
                        BXFIRST=.FALSE.
                     ENDDO
                     BZFIRST=.FALSE.
                  ENDDO
   
IF (TYPE_DEBUG >= NSCARC_DEBUG_MEDIUM) THEN
   WRITE(SCARC_LU,*) 'NX on new coarse level:',SCC%NX
   WRITE(SCARC_LU,*) 'NZ on new coarse level:',SCC%NZ
ENDIF
               !!! -------------------------------- 3D ------------------------------------------------
               CASE (NSCARC_DIMENSION_THREE)
                  DO IZ = 1, SC%NZ
                     DO IY = 1, SC%NY
                        DO IX = 1, SC%NX
                           IC = (IZ-1)*SC%NX*SC%NY + (IY-1)*SC%NX + IX
                           IF (MOD(IZ,2)==1) THEN
                              IF (MOD(IX,2)==1.AND.MOD(IY,2)==1) THEN
                                 SC%MEASURE(IC) = 6.0_EB
                              ELSE
                                 SC%MEASURE(IC) = 5.0_EB
                              ENDIF
                           ELSE
                              IF (MOD(IX,2)==0.AND.MOD(IY,2)==0) THEN
                                 SC%MEASURE(IC) = 6.0_EB
                              ELSE
                                 SC%MEASURE(IC) = 5.0_EB
                              ENDIF
                           ENDIF
                        ENDDO
                     ENDDO
                  ENDDO
            END SELECT
   
         ELSE

            !!! First Initialize internal cells to be weakly coupled fine cells
            !!! Then set internal measures corresponding to A2-coarsening
            SELECT CASE (TYPE_DIMENSION)
   
               !!! -------------------------------- 2D ------------------------------------------------
               CASE (NSCARC_DIMENSION_TWO)
                  DO IZ = 2, SC%NZ-1
                     DO IX = 2, SC%NX-1
                        IC = (IZ-1)*SC%NX + IX
                        CALL RANDOM_NUMBER(RAND_NUM)
                        SC%MEASURE(IC) = 4.0_EB + REAL(INT(RAND_NUM*8.5_EB),EB)/10.0_EB + 0.1_EB
                        SC%MEASURE(IC) = 5.0_EB + RAND_NUM
                     ENDDO
                  ENDDO
                  DO IZ = 3, SC%NZ-2, 2
                     DO IX = 3, SC%NX-2, 2
                        IC = (IZ-1)*SC%NX + IX
                        SC%MEASURE(IC)=7.0_EB
                        DO ICOL = SC%A_ROW(IC)+1, SC%A_ROW(IC+1)-1
                           JC = SC%A_COL(ICOL)
                           IF (STRONGLY_COUPLED(SC%A, SC%A_ROW, IC, ICOL)) SC%MEASURE(JC) = 6.0_EB
                        ENDDO
                     ENDDO
                  ENDDO
   
               !!! -------------------------------- 3D ------------------------------------------------
               CASE (NSCARC_DIMENSION_THREE)
                  DO IZ = 2, SC%NZ-1
                     DO IY = 2, SC%NY-1
                        DO IX = 2, SC%NX-1
                           IC = (IZ-1)*SC%NX*SC%NY + (IY-1)*SC%NX + IX
                           CALL RANDOM_NUMBER(RAND_NUM)
                           SC%MEASURE(IC) = 4.0_EB + REAL(INT(RAND_NUM*8.5_EB),EB)/10.1_EB + 0.1_EB
                           SC%MEASURE(IC) = 4.0_EB + RAND_NUM
                        ENDDO
                     ENDDO
                  ENDDO
                  DO IZ = 3, SC%NZ-2, 2
                     DO IY = 3, SC%NY-2, 2
                        DO IX = 3, SC%NX-2, 2
                           IC = (IZ-1)*SC%NX*SC%NY + (IY-1)*SC%NX + IX
                           SC%MEASURE(IC)=6.0_EB
                           DO ICOL = SC%A_ROW(IC)+1, SC%A_ROW(IC+1)-1
                              JC = SC%A_COL(ICOL)
                              IF (STRONGLY_COUPLED(SC%A, SC%A_ROW, IC, ICOL)) SC%MEASURE(JC) = 5.0_EB
                           ENDDO
                        ENDDO
                     ENDDO
                  ENDDO

            END SELECT

         ENDIF

IF (TYPE_DEBUG >= NSCARC_DEBUG_MEDIUM) THEN
   WRITE(SCARC_LU,*) '1:MEASURES in SETUP_MEASURE_FDSA2: =====  NM=',NM
   IF (NM==1) THEN
      WRITE(SCARC_LU,'(8f8.2)') (SC%MEASURE(IC), IC=1,SC%NC)
   ELSE
      WRITE(SCARC_LU,'(8f8.2)') (SC%MEASURE(IC), IC=1,SC%NC)
   ENDIF
ENDIF
      ENDDO FDSA2_LOOP

      IF (NMESHES > 1) THEN

         FDSA2_LOOP3: DO NM = NMESHES_MIN, NMESHES_MAX

            SC => SCARC(NM)%COMPACT(NL)            
            FDSA2_MEASURE_LOOP3: DO IW = 1, SC%NW
               IF (SC%GHOST_CELL(IW) == 0) CYCLE FDSA2_MEASURE_LOOP3
               IC = SC%ADJACENT_CELL(IW)                  !!! ONLY FOR NL == 1
               SC%MEASURE(IC)=0.0_EB
               DO ICOL = SC%A_ROW(IC)+1, SC%A_ROW(IC+1)-1
                  JC = SC%A_COL(ICOL)
                  IF (STRONGLY_COUPLED(SC%A, SC%A_ROW, IC, ICOL).AND.SC%MEASURE(JC)/=NSCARC_MEASURE_COARSE) THEN
                     SC%MEASURE(IC) = SC%MEASURE(IC) + 1.0_EB
                     !SC%MEASURE(IC) = 0.0_EB
                  ENDIF
               ENDDO
               CALL RANDOM_NUMBER(RAND_NUM)
               SC%MEASURE(IC) = SC%MEASURE(IC) + 0.1_EB + REAL(INT(RAND_NUM*8.9_EB),EB)/10.0_EB
           
            ENDDO FDSA2_MEASURE_LOOP3

IF (TYPE_DEBUG >= NSCARC_DEBUG_MEDIUM) THEN
   WRITE(SCARC_LU,*) '2:MEASURES in SETUP_MEASURE_RS3: =====  NM=',NM
   IF (NM==1) THEN
      WRITE(SCARC_LU,'(8f8.2)') (SC%MEASURE(IC), IC=1,SC%NC)
   ELSE
      WRITE(SCARC_LU,'(8f8.2)') (SC%MEASURE(IC), IC=1,SC%NC)
   ENDIF
ENDIF


   
         ENDDO FDSA2_LOOP3
   
      ENDIF



   !!!-------------------------------------------------------------------------------------------------
   !!! PMIS: Parallel Modified Independent Set 
   !!!     - Very fast construction with low operator complexity.  
   !!!     - Convergence can deteriorate with increasing problem size on structured meshes.  
   !!!     - Uses method similar to Luby's Maximal Independent Set algorithm.
   !!!-------------------------------------------------------------------------------------------------
   CASE (NSCARC_COARSENING_PMIS)

      PMIS_LOOP: DO NM = NMESHES_MIN, NMESHES_MAX

         SC => SCARC(NM)%COMPACT(NL)

         !!! First determine the measure of each cell and add a random number between (0,1)
         PMIS_MEASURE_LOOP: DO IC = 1, SC%NC

            DO ICOL = SC%A_ROW(IC)+1, SC%A_ROW(IC+1)-1
               IF (STRONGLY_COUPLED(SC%A, SC%A_ROW, IC, ICOL)) SC%MEASURE(IC) = SC%MEASURE(IC) + 1.0_EB
            ENDDO

            CALL RANDOM_NUMBER(RAND_NUM)
            CALL RANDOM_NUMBER(RAND_NUM2)
            RAND_NUM = 0.05_EB + 0.9_EB* RAND_NUM + 0.09_EB * RAND_NUM2
            SC%MEASURE(IC) = SC%MEASURE(IC) + RAND_NUM

         ENDDO PMIS_MEASURE_LOOP

      ENDDO PMIS_LOOP


END SELECT

IF (NMESHES > 1) CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_MEASURE, NL)

END SUBROUTINE SCARC_SETUP_MEASURES


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Setup celltypes of mesh corresponding to requested coarsening strategy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SETUP_CELLTYPES(NL)
INTEGER, PARAMETER  :: NCOUPLINGS = 20, NCYC_MAX=1000
INTEGER, INTENT(IN)    :: NL
INTEGER  :: NM, NTYPE
INTEGER  :: IC, JC, KC, ICW, LC, ICOL, JCOL, KCOL, ICP, ICYC, ICE
INTEGER  :: IDIR, JDIR, IS, IW, JW, IX, IY, IZ
INTEGER  :: NX1, NX2, NY1, NY2, NZ1, NZ2
INTEGER  :: SCOUPLED(NCOUPLINGS), WCOUPLED(NCOUPLINGS)
REAL(EB) :: MEASURE_MAX, MEASURE_LOCAL_MAX, MEASURE_GLOBAL_MAX, EPS
LOGICAL  :: BIGGEST
TYPE (SCARC_COMPACT_TYPE), POINTER :: SC

EPS = 1.0E-12
MEASURE_MAX = 0.0_EB

CALL SCARC_DEBUG_QUANTITY(NSCARC_DEBUG_MEASURE, NL, 'SETUP_MEASURES', 'MEASURE INITIAL')

!!! Use special coarsening strategies only on finest level, on coarser levels always use RS3-coarsening
IF (NL> NLEVEL_MIN) THEN
   NTYPE = NSCARC_COARSENING_TEST
ELSE
   NTYPE = TYPE_COARSENING
ENDIF

!!! Define celltypes for corresponding coarsening strategy
SELECT CASE (NTYPE)

   !!!-------------------------------------------------------------------------------------------------
   !!! TEST-coarsening
   !!!-------------------------------------------------------------------------------------------------
   CASE (NSCARC_COARSENING_TEST)

      TEST_LOOP: DO NM = NMESHES_MIN, NMESHES_MAX
      
         SC => SCARC(NM)%COMPACT(NL)
      
         TEST_CYCLE_LOOP: DO
         
            !!! get maximum (remaining) measure for all cells
            MEASURE_MAX = MAXVAL(SC%MEASURE(1:SC%NC))
            IF (MEASURE_MAX <= EPS) EXIT TEST_CYCLE_LOOP
         
            TEST_CELL_LOOP: DO IC = 1, SC%NC
         
               !!! Take first cell with maximum measure as next coarse cell
               IF (MATCH(MEASURE_MAX, SC%MEASURE(IC))) THEN
         
                  SC%MEASURE(IC)  = 0.0_EB
                  SC%CELLTYPE(IC) = NSCARC_CELLTYPE_COARSE
         
IF (TYPE_DEBUG>NSCARC_DEBUG_LESS) THEN
   WRITE(SCARC_LU,*) '========================= SETUP_CELLTYPE for IC=',IC
   WRITE(SCARC_LU,*) 'MEASURE (',IC,')=',SC%MEASURE(IC)
   WRITE(SCARC_LU,*) 'CELLTYPE(',IC,')=',SC%CELLTYPE(IC)
ENDIF
                  !!! Determine set of fine cells 
                  DO ICOL = SC%A_ROW(IC)+1, SC%A_ROW(IC+1)-1
         
                     !!! JC is set to be a fine cell which is no longer measured
                     JC = SC%A_COL(ICOL)
         
                     SC%MEASURE(JC)  = 0.0_EB
                     SC%CELLTYPE(JC) = NSCARC_CELLTYPE_SFINE
         
IF (TYPE_DEBUG>NSCARC_DEBUG_LESS) THEN
   WRITE(SCARC_LU,*) 'MEASURE (',JC,')=',SC%MEASURE(JC)
   WRITE(SCARC_LU,*) 'CELLTYPE(',JC,')=',SC%CELLTYPE(JC)
ENDIF
                     !!!  increase measures of cells KC adjacent to fine cells JC based on strong couplings
                     DO JCOL = SC%A_ROW(JC)+1, SC%A_ROW(JC+1)-1
                        KC = SC%A_COL(JCOL)
                        IF (KC > 0 .AND. KC /= IC .AND. (SC%CELLTYPE(KC)==NSCARC_CELLTYPE_NONE)) THEN
                           SC%MEASURE(KC) = SC%MEASURE(KC) + 1.0_EB
IF (TYPE_DEBUG>NSCARC_DEBUG_LESS) THEN
   WRITE(SCARC_LU,'(i3,a,i3,a,f8.2,8i6)') KC,':MEASURE (',KC,')=',SC%MEASURE(KC), SC%A_COL(SC%A_ROW(JC):SC%A_ROW(JC+1)-1)
ENDIF
                           MEASURE_MAX = MAX(MEASURE_MAX, SC%MEASURE(KC))
                        ENDIF
                     ENDDO 
         
                  ENDDO 
         
CALL SCARC_DEBUG_QUANTITY(NSCARC_DEBUG_MEASURE, NL, 'SETUP_MEASURES', 'MEASURE LOOP')
                  EXIT TEST_CELL_LOOP
               ENDIF
         
            ENDDO TEST_CELL_LOOP
         
         ENDDO TEST_CYCLE_LOOP

         SC%NCW = 0 
      
      ENDDO TEST_LOOP

      IF (NMESHES > 1) THEN

         TEST_CELLTYPE_LOOP2: DO ICYC = 1, NCYC_MAX
            
            MEASURE_GLOBAL_MAX = 0.0_EB
            TEST_LOOP2: DO NM = NMESHES_MIN, NMESHES_MAX
         
               SC => SCARC(NM)%COMPACT(NL)
               MEASURE_LOCAL_MAX  = MAXVAL(SC%MEASURE(1:SC%NC))
               MEASURE_GLOBAL_MAX = MAX(MEASURE_GLOBAL_MAX, MEASURE_LOCAL_MAX)
               IF (MEASURE_LOCAL_MAX <= EPS) CYCLE
            
               LOCAL_MAX_LOOP_TEST: DO IW = 1, SC%NW
            
                  IF (SC%GHOST_CELL(IW) == 0) CYCLE LOCAL_MAX_LOOP_TEST

                  IC = SC%ADJACENT_CELL(IW)
                  IF (SC%CELLTYPE(IC) /= NSCARC_CELLTYPE_NONE) CYCLE LOCAL_MAX_LOOP_TEST

IF (TYPE_DEBUG> NSCARC_DEBUG_LESS) THEN
   WRITE(SCARC_LU,*) 'NM=',NM,' IW=',IW,': IC=',IC
ENDIF
         
                  !!! compare measure of IC with measures of its strongly coupled neighbors
                  BIGGEST = .TRUE.
                  DO ICOL = SC%A_ROW(IC)+1, SC%A_ROW(IC+1)-1       
                     IF (STRONGLY_COUPLED(SC%A, SC%A_ROW, IC, ICOL)) THEN
                        JC = SC%A_COL(ICOL)
                        IF (SC%MEASURE(JC) > SC%MEASURE(IC)) BIGGEST = .FALSE.
                     ENDIF
                  ENDDO 
               
                  !!! if IC has biggest measure set it to be a coarse cell 
                  IF (BIGGEST) THEN
                     SC%CELLTYPE(IC) = NSCARC_CELLTYPE_COARSE
                     SC%NCW = SC%NCW + 1
IF (TYPE_DEBUG >= NSCARC_DEBUG_MEDIUM) WRITE(SCARC_LU,*) 'IC =',IC,': NCW=',SC%NCW
                  ENDIF
   
               ENDDO LOCAL_MAX_LOOP_TEST
            
               !!! set all cells which are strongly coupled to a coarse cell to be a fine cell
               INDEPENDENT_SET_LOOP_TEST: DO IW = 1, SC%NW
                  IF (SC%GHOST_CELL(IW) == 0) CYCLE INDEPENDENT_SET_LOOP_TEST
                  IC = SC%ADJACENT_CELL(IW)

IF (TYPE_DEBUG> NSCARC_DEBUG_LESS) THEN
   WRITE(SCARC_LU,*) '---> NM=',NM,' IW=',IW,': IC=',IC
ENDIF
                  IF (SC%CELLTYPE(IC) == NSCARC_CELLTYPE_COARSE) THEN
                     SC%MEASURE(IC) = 0.0_EB

IF (TYPE_DEBUG> NSCARC_DEBUG_LESS) THEN
   WRITE(SCARC_LU,*) '---> MEASURE(',IC,')=',SC%MEASURE(IC)
ENDIF
!                     IF (IC>SC%NC) CYCLE INDEPENDENT_SET_LOOP_TEST
                     DO ICOL = SC%A_ROW(IC)+1, SC%A_ROW(IC+1)-1       
                        JC = SC%A_COL(ICOL)
                        SELECT CASE (SC%CELLTYPE(JC))
                           CASE (NSCARC_CELLTYPE_NONE)
                              SC%CELLTYPE(JC) = NSCARC_CELLTYPE_SFINE
                           CASE (NSCARC_CELLTYPE_WFINE)
                              SC%CELLTYPE(JC) = NSCARC_CELLTYPE_SFINE
                        END SELECT
                        SC%MEASURE(JC) = 0.0_EB
IF (TYPE_DEBUG> NSCARC_DEBUG_LESS) THEN
   WRITE(SCARC_LU,*) '---> CELLTYPE(',IC,')=',SC%CELLTYPE(IC)
ENDIF
                     ENDDO 
                  ENDIF
               ENDDO INDEPENDENT_SET_LOOP_TEST
   
   CALL SCARC_DEBUG_QUANTITY(NSCARC_DEBUG_MEASURE , NL, 'SETUP_CELLTYPES', 'MEASURE LOOP')
   CALL SCARC_DEBUG_QUANTITY(NSCARC_DEBUG_CELLTYPE, NL, 'SETUP_CELLTYPES', 'CELLTYPE LOOP')
   
            ENDDO TEST_LOOP2
   
            IF (NMESHES > 1) THEN
               CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_MEASURE , NL)
               CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_CELLTYPE, NL)
            ENDIF
            IF (MEASURE_GLOBAL_MAX <= EPS) EXIT TEST_CELLTYPE_LOOP2
   CALL SCARC_DEBUG_QUANTITY(NSCARC_DEBUG_MEASURE , NL, 'SETUP_CELLTYPES', 'MEASURE EXCH')
   CALL SCARC_DEBUG_QUANTITY(NSCARC_DEBUG_CELLTYPE, NL, 'SETUP_CELLTYPES', 'CELLTYPE EXCH')
   
         ENDDO TEST_CELLTYPE_LOOP2
         

      ENDIF
      
   !!!-------------------------------------------------------------------------------------------------
   !!! RS3-coarsening
   !!!-------------------------------------------------------------------------------------------------
   CASE (NSCARC_COARSENING_RS3)

IF (TYPE_DEBUG> NSCARC_DEBUG_LESS) THEN
   WRITE(SCARC_LU,*) ' SETTING UP CELLTYPES RS3'
ENDIF
      RS3_LOOP: DO NM = NMESHES_MIN, NMESHES_MAX
      
         SC => SCARC(NM)%COMPACT(NL)
      
         RS3_CYCLE_LOOP: DO
         
            !!! get maximum (remaining) measure for all cells
            MEASURE_MAX = MAXVAL(SC%MEASURE(1:SC%NC))
            IF (MEASURE_MAX <= EPS) EXIT RS3_CYCLE_LOOP
         
            RS3_CELL_LOOP: DO IC = 1, SC%NC
         
               !!! Take first cell with maximum measure as next coarse cell
               IF (MATCH(MEASURE_MAX, SC%MEASURE(IC))) THEN
         
                  SC%MEASURE(IC)  = 0.0_EB
                  SC%CELLTYPE(IC) = NSCARC_CELLTYPE_COARSE
         
IF (TYPE_DEBUG>NSCARC_DEBUG_LESS) THEN
   WRITE(SCARC_LU,*) '========================= SETUP_CELLTYPE for IC=',IC
   WRITE(SCARC_LU,*) 'MEASURE (',IC,')=',SC%MEASURE(IC)
   WRITE(SCARC_LU,*) 'CELLTYPE(',IC,')=',SC%CELLTYPE(IC)
ENDIF
                  !!! Determine set of fine cells 
                  DO ICOL = SC%A_ROW(IC)+1, SC%A_ROW(IC+1)-1

IF (TYPE_DEBUG>NSCARC_DEBUG_LESS) THEN
   WRITE(SCARC_LU,*) '--------------------------- ICOL =',ICOL
ENDIF
                     IF (STRONGLY_COUPLED(SC%A, SC%A_ROW, IC, ICOL)) THEN

                     !!! JC is set to be a fine cell which is no longer measured
                     JC = SC%A_COL(ICOL)

IF (TYPE_DEBUG>NSCARC_DEBUG_LESS) THEN
   WRITE(SCARC_LU,*) '--------------------------- yes = JC',JC
ENDIF
                     SC%MEASURE(JC)  = 0.0_EB
                     SC%CELLTYPE(JC) = NSCARC_CELLTYPE_FINE
         
IF (TYPE_DEBUG>NSCARC_DEBUG_LESS) THEN
   WRITE(SCARC_LU,*) 'MEASURE (',JC,')=',SC%MEASURE(JC)
   WRITE(SCARC_LU,*) 'CELLTYPE(',JC,')=',SC%CELLTYPE(JC)
ENDIF
                     !!!  increase measures of cells KC adjacent to fine cells JC based on strong couplings
                     DO JCOL = SC%A_ROW(JC)+1, SC%A_ROW(JC+1)-1
                        IF (STRONGLY_COUPLED(SC%A, SC%A_ROW, JC, JCOL)) THEN
                        KC = SC%A_COL(JCOL)
                        IF (KC /= IC .AND. (SC%CELLTYPE(KC)==NSCARC_CELLTYPE_NONE)) THEN
                           SC%MEASURE(KC) = SC%MEASURE(KC) + 1.0_EB
IF (TYPE_DEBUG>NSCARC_DEBUG_LESS) THEN
   WRITE(SCARC_LU,*) 'MEASURE (',KC,')=',SC%MEASURE(KC)
ENDIF
                           MEASURE_MAX = MAX(MEASURE_MAX, SC%MEASURE(KC))
                        ENDIF
                        ENDIF
                     ENDDO 
                     ENDIF
         
                  ENDDO 
         
                  EXIT RS3_CELL_LOOP
               ENDIF
         
            ENDDO RS3_CELL_LOOP
         
         ENDDO RS3_CYCLE_LOOP

         SC%NCW = 0 
      
      ENDDO RS3_LOOP

IF (TYPE_DEBUG> NSCARC_DEBUG_LESS) THEN
   WRITE(SCARC_LU,'(8f8.2)') (SC%MEASURE(IC), IC=1,SC%NC)
   WRITE(SCARC_LU,'(8i8)') (SC%CELLTYPE(IC), IC=1,SC%NC)
ENDIF

      IF (NMESHES > 1) THEN

IF (TYPE_DEBUG> NSCARC_DEBUG_LESS) THEN
   WRITE(SCARC_LU,*) ' BOUNDARY LOOP OF SETUP_CELLTYPES on LEVEL ', NL 
ENDIF
         RS3_CELLTYPE_LOOP2: DO ICYC = 1, NCYC_MAX
            
IF (TYPE_DEBUG> NSCARC_DEBUG_LESS) THEN
   WRITE(SCARC_LU,*) ' ================= ICYC=',ICYC
ENDIF
            MEASURE_GLOBAL_MAX = 0.0_EB
            RS3_LOOP2: DO NM = NMESHES_MIN, NMESHES_MAX
         
               SC => SCARC(NM)%COMPACT(NL)
               MEASURE_LOCAL_MAX  = MAXVAL(SC%MEASURE(1:SC%NC))
               MEASURE_GLOBAL_MAX = MAX(MEASURE_GLOBAL_MAX, MEASURE_LOCAL_MAX)
               IF (MEASURE_LOCAL_MAX <= EPS) CYCLE
            
IF (TYPE_DEBUG> NSCARC_DEBUG_LESS) THEN
   WRITE(SCARC_LU,*) 'MEASURE_LOCAL_MAX=',MEASURE_LOCAL_MAX
   WRITE(SCARC_LU,*) 'MEASURE_GLOBAL_MAX=',MEASURE_GLOBAL_MAX
ENDIF
               LOCAL_MAX_LOOP_RS3: DO IW = 1, SC%NW
            
                  IF (NL == 1) THEN
                     !IF (SC%GHOST_CELL(IW) == 0) CYCLE LOCAL_MAX_LOOP_RS3
                     IF (SC%IJKW(9,(IW)) == 0) CYCLE LOCAL_MAX_LOOP_RS3
                     IC = SC%ADJACENT_CELL(IW)
                  ELSE
                     IC = SC%IJKW(1,IW)    !!!! ????? wirklich richtig ??????????
                  ENDIF

                  IF (SC%CELLTYPE(IC) /= NSCARC_CELLTYPE_NONE) CYCLE LOCAL_MAX_LOOP_RS3

IF (TYPE_DEBUG> NSCARC_DEBUG_LESS) THEN
   WRITE(SCARC_LU,*) 'NM=',NM,' IW=',IW,': IC=',IC
ENDIF
         
                  !!! compare measure of IC with measures of its strongly coupled neighbors
                  BIGGEST = .TRUE.
                  DO ICOL = SC%A_ROW(IC)+1, SC%A_ROW(IC+1)-1       
                     IF (STRONGLY_COUPLED(SC%A, SC%A_ROW, IC, ICOL)) THEN
                        JC = SC%A_COL(ICOL)
                        IF (SC%MEASURE(JC) > SC%MEASURE(IC)) BIGGEST = .FALSE.
                     ENDIF
                  ENDDO 
               
                  !!! if IC has biggest measure set it to be a coarse cell 
                  IF (BIGGEST) THEN
                     SC%CELLTYPE(IC) = NSCARC_CELLTYPE_COARSE
                     SC%NCW = SC%NCW + 1
IF (TYPE_DEBUG >= NSCARC_DEBUG_MEDIUM) WRITE(SCARC_LU,*) 'IC =',IC,': NCW=',SC%NCW
                  ENDIF
   
               ENDDO LOCAL_MAX_LOOP_RS3
            
               !!! set all cells which are strongly coupled to a coarse cell to be a fine cell
               INDEPENDENT_SET_LOOP_RS3: DO IW = 1, SC%NW
                  IF (NL == 1) THEN
                     !IF (SC%GHOST_CELL(IW) == 0) CYCLE INDEPENDENT_SET_LOOP_RS3
                     IF (SC%IJKW(9,(IW)) == 0) CYCLE INDEPENDENT_SET_LOOP_RS3
                     IC = SC%ADJACENT_CELL(IW)
                  ELSE
                     IC = SC%IJKW(1,IW)
                  ENDIF

IF (TYPE_DEBUG> NSCARC_DEBUG_LESS) THEN
   WRITE(SCARC_LU,*) '---> NM=',NM,' IW=',IW,': IC=',IC
ENDIF
                  IF (SC%CELLTYPE(IC) == NSCARC_CELLTYPE_COARSE) THEN
                     SC%MEASURE(IC) = 0.0_EB

IF (TYPE_DEBUG> NSCARC_DEBUG_LESS) THEN
   WRITE(SCARC_LU,*) '---> MEASURE(',IC,')=',SC%MEASURE(IC)
ENDIF
!                     IF (IC>SC%NC) CYCLE INDEPENDENT_SET_LOOP_RS3
                     DO ICOL = SC%A_ROW(IC)+1, SC%A_ROW(IC+1)-1       
                        JC = SC%A_COL(ICOL)
                        SELECT CASE (SC%CELLTYPE(JC))            !!! hier noch mal schauen, sinnvoll ?
                           CASE (NSCARC_CELLTYPE_NONE)
                              SC%CELLTYPE(JC) = NSCARC_CELLTYPE_SFINE
                           CASE (NSCARC_CELLTYPE_WFINE)
                              SC%CELLTYPE(JC) = NSCARC_CELLTYPE_SFINE
                        END SELECT
                        SC%MEASURE(JC) = 0.0_EB
IF (TYPE_DEBUG> NSCARC_DEBUG_LESS) THEN
   WRITE(SCARC_LU,*) '---> CELLTYPE(',IC,')=',SC%CELLTYPE(IC)
ENDIF
                     ENDDO 
                  ENDIF
               ENDDO INDEPENDENT_SET_LOOP_RS3
   
   CALL SCARC_DEBUG_QUANTITY(NSCARC_DEBUG_MEASURE , NL, 'SETUP_CELLTYPES', 'MEASURE LOOP')
   CALL SCARC_DEBUG_QUANTITY(NSCARC_DEBUG_CELLTYPE, NL, 'SETUP_CELLTYPES', 'CELLTYPE LOOP')
   
            ENDDO RS3_LOOP2
   
            IF (NMESHES > 1) THEN
               CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_MEASURE , NL)
               CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_CELLTYPE, NL)
            ENDIF
            IF (MEASURE_GLOBAL_MAX <= EPS) EXIT RS3_CELLTYPE_LOOP2
   CALL SCARC_DEBUG_QUANTITY(NSCARC_DEBUG_MEASURE , NL, 'SETUP_CELLTYPES', 'MEASURE EXCH')
   CALL SCARC_DEBUG_QUANTITY(NSCARC_DEBUG_CELLTYPE, NL, 'SETUP_CELLTYPES', 'CELLTYPE EXCH')
   
         ENDDO RS3_CELLTYPE_LOOP2
         

      ENDIF
      

   !!!-------------------------------------------------------------------------------------------------
   !!! FDSA2-coarsening
   !!!-------------------------------------------------------------------------------------------------
   CASE (NSCARC_COARSENING_FDSA2)

IF (TYPE_DEBUG> NSCARC_DEBUG_LESS) THEN
   WRITE(SCARC_LU,*) ' SETTING UP CELLTYPES FDSA2'
ENDIF
      FDSA2_LOOP1: DO NM = NMESHES_MIN, NMESHES_MAX

         SC => SCARC(NM)%COMPACT(NL)
         FDSA2_CELLTYPE_LOOP: DO IC = 1, SC%NC
            IF (MATCH(NSCARC_MEASURE_COARSE, SC%MEASURE(IC))) THEN
               SC%CELLTYPE(IC) = NSCARC_CELLTYPE_COARSE
               SC%MEASURE(IC)  = NSCARC_MEASURE_NONE
            ELSE IF (MATCH(NSCARC_MEASURE_SFINE, SC%MEASURE(IC))) THEN
               SC%CELLTYPE(IC) = NSCARC_CELLTYPE_SFINE
               SC%MEASURE(IC)  = NSCARC_MEASURE_NONE
            ELSE IF (MATCH(NSCARC_MEASURE_WFINE, SC%MEASURE(IC))) THEN
               SC%CELLTYPE(IC) = NSCARC_CELLTYPE_WFINE
               SC%MEASURE(IC)  = NSCARC_MEASURE_NONE
            ENDIF
         ENDDO FDSA2_CELLTYPE_LOOP

         SC%NCW = 0 

      ENDDO FDSA2_LOOP1

CALL SCARC_DEBUG_QUANTITY(NSCARC_DEBUG_MEASURE, NL, 'SETUP_CELLTYPES', 'CELLTYPES INITIAL1')
CALL SCARC_DEBUG_QUANTITY(NSCARC_DEBUG_CELLTYPE, NL, 'SETUP_CELLTYPES', 'CELLTYPE INITIAL1')

      IF (NMESHES > 1) THEN

         FDSA2_CELLTYPE_LOOP2: DO ICYC = 1, NCYC_MAX
            
IF (TYPE_DEBUG> NSCARC_DEBUG_LESS) THEN
   WRITE(SCARC_LU,*) '============= CyC =',ICYC
ENDIF
            MEASURE_GLOBAL_MAX = 0.0_EB
            FDSA2_LOOP2: DO NM = NMESHES_MIN, NMESHES_MAX
         
               SC => SCARC(NM)%COMPACT(NL)
               MEASURE_LOCAL_MAX  = MAXVAL(SC%MEASURE(1:SC%NC))
               MEASURE_GLOBAL_MAX = MAX(MEASURE_GLOBAL_MAX, MEASURE_LOCAL_MAX)
IF (TYPE_DEBUG> NSCARC_DEBUG_LESS) THEN
   WRITE(SCARC_LU,*) 'MEASURES:' , MEASURE_LOCAL_MAX, MEASURE_GLOBAL_MAX
ENDIF
               IF (MEASURE_LOCAL_MAX <= EPS) CYCLE FDSA2_LOOP2
IF (TYPE_DEBUG> NSCARC_DEBUG_LESS) THEN
   WRITE(SCARC_LU,*) 'MEASURES: after '
ENDIF
            
               LOCAL_MAX_LOOP_FDSA2: DO IW = 1, SC%NW
            
                  IF (SC%GHOST_CELL(IW) == 0) CYCLE LOCAL_MAX_LOOP_FDSA2

                  IC = SC%ADJACENT_CELL(IW)
                  IF (SC%CELLTYPE(IC) /= NSCARC_CELLTYPE_NONE) CYCLE LOCAL_MAX_LOOP_FDSA2

IF (TYPE_DEBUG> NSCARC_DEBUG_LESS) THEN
   WRITE(SCARC_LU,*) 'NM=',NM,' IW=',IW,': IC=',IC
ENDIF
         
                  !!! compare measure of IC with measures of its strongly coupled neighbors
                  BIGGEST = .TRUE.
                  DO ICOL = SC%A_ROW(IC)+1, SC%A_ROW(IC+1)-1       
                     JC = SC%A_COL(ICOL)
                     IF (SC%CELLTYPE(JC) == NSCARC_CELLTYPE_COARSE) CYCLE
                     IF (STRONGLY_COUPLED(SC%A, SC%A_ROW, IC, ICOL)) THEN
                        JC = SC%A_COL(ICOL)
                        IF (SC%MEASURE(JC) > SC%MEASURE(IC)) BIGGEST = .FALSE.
                     ENDIF
                  ENDDO 
               
                  !!! if IC has biggest measure set it to be a coarse cell 
                  IF (BIGGEST) THEN
                     SC%CELLTYPE(IC) = NSCARC_CELLTYPE_COARSE
                     SC%NCW = SC%NCW + 1
IF (TYPE_DEBUG >= NSCARC_DEBUG_MEDIUM) WRITE(SCARC_LU,*) 'IC =',IC,': NCW=',SC%NCW
                  ENDIF
   
               ENDDO LOCAL_MAX_LOOP_FDSA2
            
               !!! set all cells which are strongly coupled to a coarse cell to be a fine cell
               INDEPENDENT_SET_LOOP_FDSA2: DO IW = 1, SC%NW
                  IF (SC%GHOST_CELL(IW) == 0) CYCLE INDEPENDENT_SET_LOOP_FDSA2
                  IC = SC%ADJACENT_CELL(IW)

IF (TYPE_DEBUG> NSCARC_DEBUG_LESS) THEN
   WRITE(SCARC_LU,*) '-------------------'
   WRITE(SCARC_LU,*) '---> NM=',NM,' IW=',IW,': IC=',IC
ENDIF
                  IF (SC%CELLTYPE(IC) == NSCARC_CELLTYPE_COARSE) THEN
                     SC%MEASURE(IC) = 0.0_EB

                     DO ICOL = SC%A_ROW(IC)+1, SC%A_ROW(IC+1)-1
         
                     !!! JC is set to be a fine cell which is no longer measured
                        JC = SC%A_COL(ICOL)
         
                        IF (SC%CELLTYPE(JC) == NSCARC_CELLTYPE_COARSE) THEN
                           WRITE(*,*) 'Wrong coarsening at node ',JC,' on mesh ', NM,' main IC=',IC
                           stop
                        ENDIF

                        SC%MEASURE(JC)  = 0.0_EB
                        SC%CELLTYPE(JC) = NSCARC_CELLTYPE_FINE
         
IF (TYPE_DEBUG>NSCARC_DEBUG_LESS) THEN
   WRITE(SCARC_LU,*) 'MEASURE (',JC,')=',SC%MEASURE(JC)
   WRITE(SCARC_LU,*) 'CELLTYPE(',JC,')=',SC%CELLTYPE(JC)
ENDIF
                     !!!  increase measures of cells KC adjacent to fine cells JC based on strong couplings
                        DO JCOL = SC%A_ROW(JC)+1, SC%A_ROW(JC+1)-1
                           KC = SC%A_COL(JCOL)
                           IF (KC /= IC .AND. (SC%CELLTYPE(KC)==NSCARC_CELLTYPE_NONE)) THEN
                              SC%MEASURE(KC) = SC%MEASURE(KC) + 1.0_EB
IF (TYPE_DEBUG>NSCARC_DEBUG_LESS) THEN
   WRITE(SCARC_LU,*) 'MEASURE (',KC,')=',SC%MEASURE(KC)
ENDIF
                              MEASURE_LOCAL_MAX = MAX(MEASURE_LOCAL_MAX, SC%MEASURE(KC))
                           ENDIF
                        ENDDO 
         
                     ENDDO 
IF (TYPE_DEBUG> NSCARC_DEBUG_LESS) THEN
   WRITE(SCARC_LU,*) 'A: ---> MEASURE (',IC,')=',SC%MEASURE(IC)
ENDIF
                     ENDIF
               ENDDO INDEPENDENT_SET_LOOP_FDSA2
   
   CALL SCARC_DEBUG_QUANTITY(NSCARC_DEBUG_MEASURE , NL, 'SETUP_CELLTYPES', 'MEASURE LOOP')
   CALL SCARC_DEBUG_QUANTITY(NSCARC_DEBUG_CELLTYPE, NL, 'SETUP_CELLTYPES', 'CELLTYPE LOOP')
   
            ENDDO FDSA2_LOOP2
   
WRITE(SCARC_LU,*) 'AFTER FDSA2_LOOP2 A'

            IF (MEASURE_GLOBAL_MAX <= EPS) EXIT FDSA2_CELLTYPE_LOOP2
WRITE(SCARC_LU,*) 'AFTER FDSA2_LOOP2 B'

            IF (NMESHES > 1) THEN
               CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_MEASURE , NL)
               CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_CELLTYPE, NL)
            ENDIF

   CALL SCARC_DEBUG_QUANTITY(NSCARC_DEBUG_MEASURE , NL, 'SETUP_CELLTYPES', 'MEASURE LOOP2')
   CALL SCARC_DEBUG_QUANTITY(NSCARC_DEBUG_CELLTYPE, NL, 'SETUP_CELLTYPES', 'CELLTYPE LOOP2')
   
            DO NM = NMESHES_MIN, NMESHES_MAX
IF (TYPE_DEBUG> NSCARC_DEBUG_LESS) THEN
   WRITE(SCARC_LU,*) '------------------- NM = ',NM,' NW=',SC%NW
ENDIF
               SC => SCARC(NM)%COMPACT(NL)
               DO IW = 1, SC%NW
                  ICW = SC%GHOST_CELL(IW)
                  IF (ICW== 0) CYCLE
                  IC = SC%ADJACENT_CELL(IW)
              
IF (TYPE_DEBUG> NSCARC_DEBUG_LESS) THEN
   WRITE(SCARC_LU,*) 'NM=',NM,' IW=',IW,': IC=',IC
ENDIF
                  IF (SC%CELLTYPE(ICW) == NSCARC_CELLTYPE_COARSE) THEN
                     SC%MEASURE(IC)  = 0.0_EB
                     SC%CELLTYPE(IC) = NSCARC_CELLTYPE_SFINE
                  ENDIF

               ENDDO
            ENDDO

   CALL SCARC_DEBUG_QUANTITY(NSCARC_DEBUG_MEASURE , NL, 'SETUP_CELLTYPES', 'MEASURE EXCH')
   CALL SCARC_DEBUG_QUANTITY(NSCARC_DEBUG_CELLTYPE, NL, 'SETUP_CELLTYPES', 'CELLTYPE EXCH')
   
         ENDDO FDSA2_CELLTYPE_LOOP2
         

      ENDIF
      

   !!!-------------------------------------------------------------------------------------------------
   !!! PMIS-coarsening
   !!!-------------------------------------------------------------------------------------------------
   CASE (NSCARC_COARSENING_PMIS, NSCARC_COARSENING_PMIS2)

      
      PMIS_CELLTYPE_LOOP: DO ICYC = 1, NCYC_MAX
         
         MEASURE_GLOBAL_MAX = 0.0_EB
         PMIS_LOOP: DO NM = NMESHES_MIN, NMESHES_MAX
      
            SC => SCARC(NM)%COMPACT(NL)
            MEASURE_LOCAL_MAX  = MAXVAL(SC%MEASURE(1:SC%NC))
            MEASURE_GLOBAL_MAX = MAX(MEASURE_GLOBAL_MAX, MEASURE_LOCAL_MAX)
            IF (MEASURE_LOCAL_MAX <= EPS) CYCLE
         
            LOCAL_MAX_LOOP_PMIS: DO IC = 1, SC%NC
         
               IF (SC%CELLTYPE(IC) /= NSCARC_CELLTYPE_NONE) CYCLE LOCAL_MAX_LOOP_PMIS
      
               !!! compare measure of IC with measures of its strongly coupled neighbors
               BIGGEST = .TRUE.
               DO ICOL = SC%A_ROW(IC)+1, SC%A_ROW(IC+1)-1       
                  IF (STRONGLY_COUPLED(SC%A, SC%A_ROW, IC, ICOL)) THEN
                     JC = SC%A_COL(ICOL)
                     IF (SC%MEASURE(JC) > SC%MEASURE(IC)) BIGGEST = .FALSE.
                  ENDIF
               ENDDO 
            
               !!! if IC has biggest measure set it to be a coarse cell 
               IF (BIGGEST) SC%CELLTYPE(IC) = NSCARC_CELLTYPE_COARSE

            ENDDO LOCAL_MAX_LOOP_PMIS
         
            !!! set all cells which are strongly coupled to a coarse cell to be a fine cell
            INDEPENDENT_SET_LOOP_PMIS: DO IC = 1, SC%NCE
               IF (SC%CELLTYPE(IC) == NSCARC_CELLTYPE_COARSE) THEN
                  SC%MEASURE(IC) = 0.0_EB
                  IF (IC>SC%NC) CYCLE INDEPENDENT_SET_LOOP_PMIS
                  DO ICOL = SC%A_ROW(IC)+1, SC%A_ROW(IC+1)-1       
                     JC = SC%A_COL(ICOL)
                     IF (SC%CELLTYPE(JC) == NSCARC_CELLTYPE_NONE) THEN
                        SC%CELLTYPE(JC) = NSCARC_CELLTYPE_SFINE
                        SC%MEASURE(JC) = 0.0_EB
                     ENDIF
                  ENDDO 
               ENDIF
            ENDDO INDEPENDENT_SET_LOOP_PMIS

CALL SCARC_DEBUG_QUANTITY(NSCARC_DEBUG_MEASURE , NL, 'SETUP_CELLTYPES', 'MEASURE LOOP')
CALL SCARC_DEBUG_QUANTITY(NSCARC_DEBUG_CELLTYPE, NL, 'SETUP_CELLTYPES', 'CELLTYPE LOOP')

         ENDDO PMIS_LOOP

         IF (NMESHES > 1) THEN
            CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_MEASURE , NL)
            CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_CELLTYPE, NL)
         ENDIF
         IF (MEASURE_GLOBAL_MAX <= EPS) EXIT PMIS_CELLTYPE_LOOP
CALL SCARC_DEBUG_QUANTITY(NSCARC_DEBUG_MEASURE , NL, 'SETUP_CELLTYPES', 'MEASURE EXCH')
CALL SCARC_DEBUG_QUANTITY(NSCARC_DEBUG_CELLTYPE, NL, 'SETUP_CELLTYPES', 'CELLTYPE EXCH')

      ENDDO PMIS_CELLTYPE_LOOP
         
      

END SELECT


!!!---------------------------------------------------------------------------------------------------------
!!! Exchange celltypes along internal boundaries
!!!---------------------------------------------------------------------------------------------------------
CALL SCARC_DEBUG_QUANTITY(NSCARC_DEBUG_CELLTYPE, NL, 'SETUP_CELLTYPES', 'CELLTYPE FINAL0')
IF (NMESHES > 1) CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_CELLTYPE, NL)
CALL SCARC_DEBUG_QUANTITY(NSCARC_DEBUG_CELLTYPE, NL, 'SETUP_CELLTYPES', 'CELLTYPE FINAL1')


!!!---------------------------------------------------------------------------------------------------------
!!! Define sizes for transfer matrices (prolongation/restriction)
!!!---------------------------------------------------------------------------------------------------------
MESHES_LOOP: DO NM = NMESHES_MIN, NMESHES_MAX

   SC => SCARC(NM)%COMPACT(NL)
   
   !!! Determine dimensions of restriction and prolongation matrices
   SC%NCF  = 0
   SC%NCC  = 0
   SC%NP   = 0
   SC%NR   = 0
   DO IC = 1, SC%NC
      SELECT CASE (SC%CELLTYPE(IC))
         CASE (NSCARC_CELLTYPE_COARSE)
            SC%NCC = SC%NCC + 1
            SC%NP  = SC%NP  + 1
         CASE (NSCARC_CELLTYPE_SFINE, NSCARC_CELLTYPE_WFINE)
            SC%NCF = SC%NCF + 1
            SC%NP  = SC%NP  + SC%A_ROW(IC+1)-SC%A_ROW(IC) - 1
      END SELECT
   ENDDO

   !!! Determine number of coarse and fine cells and check correctness of computation
   IF (SC%NCC + SC%NCF /= SC%NC) THEN
      WRITE(*,*) 'Error in AMG standard coarsening, N_CELLS_COARSE + N_CELLS_FINE = ', SC%NCC + SC%NCF, &
                 ' differs from N_CELLS = ', SC%NC, ' on level ', NL
   !   STOP
   ENDIF
   
   !!! Determine new numbering for coarse cells in interior of mesh
   ICP   = 0
   DO IC = 1, SC%NC
      IF (SC%CELLTYPE(IC) == NSCARC_CELLTYPE_COARSE) THEN
         ICP = ICP + 1
         SC%CELLTYPE(IC) = ICP
      ENDIF
   ENDDO

ENDDO MESHES_LOOP

CALL SCARC_DEBUG_QUANTITY(NSCARC_DEBUG_CELLTYPE, NL, 'SETUP_CELLTYPES', 'CELLTYPE FINAL2')
CALL SCARC_DEBUG_QUANTITY(NSCARC_DEBUG_MEASURE , NL, 'SETUP_MEASURES', 'MEASURE FINAL2')

END SUBROUTINE SCARC_SETUP_CELLTYPES


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Set cell types for A1 coarsening
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SET_TYPES_A1(MEASURE, CELLTYPE, JC, OFFSET)
REAL(EB), DIMENSION(:), INTENT(OUT) :: MEASURE
INTEGER , DIMENSION(:), INTENT(OUT) :: CELLTYPE
INTEGER , INTENT(IN) :: JC, OFFSET

MEASURE(JC)         = 0.0_EB
MEASURE(JC-OFFSET)  = 0.0_EB
MEASURE(JC+OFFSET)  = 0.0_EB

CELLTYPE(JC)        = NSCARC_CELLTYPE_SFINE
CELLTYPE(JC-OFFSET) = NSCARC_CELLTYPE_WFINE
CELLTYPE(JC+OFFSET) = NSCARC_CELLTYPE_WFINE
                     
END SUBROUTINE SET_TYPES_A1


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Get maximum measure in internal of mesh (cells along boundaries not regarded)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL(EB) FUNCTION INTERNAL_MAX(MEASURE, NX, NX1, NX2, NY, NY1, NY2, NZ1, NZ2)
REAL(EB), DIMENSION(:), INTENT(IN) :: MEASURE
INTEGER , INTENT(IN) :: NX, NX1, NX2, NY, NY1, NY2, NZ1, NZ2
INTEGER  :: IX, IY, IZ, IC
INTERNAL_MAX = 0.0_EB
DO IZ = NZ1, NZ2
   DO IY = NY1, NY2
      DO IX = NX1, NX2
         IC = (IZ-1) * NX * NY + (IY-1) * NX + IX
         INTERNAL_MAX = MAX(INTERNAL_MAX, MEASURE(IC))
      ENDDO
   ENDDO
ENDDO
RETURN
END FUNCTION INTERNAL_MAX


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Get cell number of internal neighbor for ghost cell IW
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
INTEGER FUNCTION GET_ADJACENT_CELL(IJKW, IW, NX, NY)
INTEGER, DIMENSION(:,:), INTENT(IN) :: IJKW
INTEGER, INTENT(IN) :: IW, NX, NY
INTEGER :: IC

SELECT CASE (TYPE_DIMENSION)
   CASE (NSCARC_DIMENSION_TWO)
      IF (ABS(IJKW(4,IW)) == 2) THEN
         IC = -1
      ELSE
         IC = (IJKW(8,IW)-1)*NX + IJKW(6,IW)
      ENDIF
   CASE (NSCARC_DIMENSION_THREE)
      IC = (IJKW(8,IW)-1)*NX*NY + (IJKW(7,IW)-1)*NX + IJKW(6,IW)
END SELECT

GET_ADJACENT_CELL = IC
RETURN

END FUNCTION GET_ADJACENT_CELL


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Perform interpolation to coarse level
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SETUP_PROLONGATION(NL)
INTEGER , INTENT(IN) :: NL
INTEGER  :: NM, IP, IC, JC, KC, ICOL, JCOL, ICG, ICA
INTEGER  :: IW, IW0, IW2, IDIAG, JDIAG, NOM, NXO, NYO, JCO
REAL(EB) :: SUM_COUPLED, SUM_CPOINTS, SCAL, SUM_COARSE, SUM_DIAG
REAL(EB) :: VALUES(20), WEIGHTS(20)
INTEGER  :: NEIGHBOR(20), NWEIGHTS, NWEIGHTS2
INTEGER  :: COARSE_CELL(20), COARSE_INDEX(20)
TYPE (SCARC_COMPACT_TYPE), POINTER :: SC

SELECT_INTERPOLATION: SELECT CASE (TYPE_INTERPOL)

   !!! -------------------------------------------------------------------------------------------------
   !!! Standard interpolation
   !!! -------------------------------------------------------------------------------------------------
   CASE (NSCARC_INTERPOL_STANDARD)

      MESHES_LOOP_STANDARD: DO NM = NMESHES_MIN, NMESHES_MAX
         WRITE(*,*) 'Standard interpolatiion not yet implemented'
      ENDDO MESHES_LOOP_STANDARD
      stop

   !!! -------------------------------------------------------------------------------------------------
   !!! Classical interpolation
   !!! -------------------------------------------------------------------------------------------------
   CASE (NSCARC_INTERPOL_CLASSICAL)

      MESHES_LOOP_CLASSICAL: DO NM = NMESHES_MIN, NMESHES_MAX

         SC => SCARC(NM)%COMPACT(NL)
         CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_CELLTYPE, NL, 'SETUP_PROLONGATION','CELLTYPE INIT')
         
         IP = 1
         INTERNAL_CELL_LOOP_CLASSICAL: DO IC = 1, SC%NC
         
            VALUES   = 0.0_EB
            WEIGHTS  = 0.0_EB
            NEIGHBOR = 0
         
            !!!
            !!! If IC is a coarse cell, its value is taken
            !!!
            IF (SC%CELLTYPE(IC) >= NSCARC_CELLTYPE_COARSE) THEN
         
               NEIGHBOR(1)= SC%CELLTYPE(IC)
               WEIGHTS(1) = 1.0_EB
               NWEIGHTS   = 1

IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) THEN
WRITE(SCARC_LU,*) '==============================================!!!!! IC=',IC
WRITE(SCARC_LU,*) 'IC=',IC,': COARSE CELL WEIGHTS=',WEIGHTS(1)
ENDIF
         
            !!!
            !!! If IC is a fine cell, a mean value must be computed based on surrounding coarse cells
            !!!
            ELSE 
         
               !!! Get main diagonal entry a_ii for that fine cell
               IDIAG = SC%A_ROW(IC)
               SUM_DIAG = SC%A(IDIAG)
         
               !!! First search for all neighboring coarse grid cells (store them in NEIGHBOR)
               IW = 1
               DO ICOL = SC%A_ROW(IC)+1, SC%A_ROW(IC+1)-1
                  JC = SC%A_COL(ICOL)
                  IF (SC%CELLTYPE(JC) >= NSCARC_CELLTYPE_COARSE) THEN
                     NEIGHBOR(IW) = JC
                     WEIGHTS(IW)  = -SC%A(ICOL)
                     IW = IW + 1
                  ENDIF
               ENDDO
               NWEIGHTS = IW - 1


IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) THEN
WRITE(SCARC_LU,*) '=================== COARSE =============!!!!! IC=',IC
WRITE(SCARC_LU,'(a,i3,a,9i8)') 'IC=',IC,': SFINE CELL has coarse NEIGHBOR=',NEIGHBOR(1:9)
WRITE(SCARC_LU,'(a,i3,a,9f8.2)') 'IC=',IC,': SFINE CELL has coarse WEIGHTS= ',WEIGHTS(1:9)
ENDIF

               !!! Then search for the strongly and weakly coupled fine grid cells
               DO ICOL = SC%A_ROW(IC)+1, SC%A_ROW(IC+1)-1

                  IW2 = 1
                  SUM_COARSE = 0.0_EB

                  JC = SC%A_COL(ICOL)

IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) THEN
WRITE(SCARC_LU,*) '=================== FINE   =============!!!!! IC=',IC
WRITE(SCARC_LU,*) '        ICOL=',ICOL,': JC=',JC
ENDIF
                  SELECT CASE (SC%CELLTYPE(JC))

                     CASE (NSCARC_CELLTYPE_SFINE)

                        !!! search for couplings KC of the strongly coupled JC which belong to the 
                        !!! coarse interpolatory set of IC
                        DO JCOL = SC%A_ROW(JC)+1, SC%A_ROW(JC+1)-1
                           KC = SC%A_COL(JCOL)

IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) THEN
WRITE(SCARC_LU,*) '        JCOL=',JCOL,': KC=',KC
ENDIF
                           NWEIGHTS_SFINE_LOOP: DO IW = 1, NWEIGHTS
                              IF (KC == NEIGHBOR(IW) ) THEN
                                 COARSE_CELL (IW2) = JCOL
                                 COARSE_INDEX(IW2) = IW
                                 SUM_COARSE = SUM_COARSE + SC%A(JCOL)

IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) THEN
WRITE(SCARC_LU,*) '        --->  COARSE_CELL (',IW2,')=',COARSE_CELL(IW2)
WRITE(SCARC_LU,*) '              COARSE_INDEX(',IW2,')=',COARSE_INDEX(IW2)
WRITE(SCARC_LU,*) '              SUM_COARSE=',SUM_COARSE
ENDIF
                                 IW2 = IW2 + 1
                                 EXIT NWEIGHTS_SFINE_LOOP
                              ENDIF
                           ENDDO NWEIGHTS_SFINE_LOOP
                        ENDDO
                        NWEIGHTS2 = IW2 - 1

                        DO IW2 = 1, NWEIGHTS2
                           JCOL = COARSE_CELL(IW2)
                           IW   = COARSE_INDEX(IW2) 
                           WEIGHTS(IW) = WEIGHTS(IW) - SC%A(ICOL)*SC%A(JCOL)/REAL(SUM_COARSE,EB)
IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) THEN
WRITE(SCARC_LU,*) '        --->  JCOL=',JCOL,': IW=',IW
WRITE(SCARC_LU,*) '              WEIGHTS(',IW,')=',WEIGHTS(IW)
ENDIF
                        ENDDO
                        
                     CASE (NSCARC_CELLTYPE_WFINE)
                        
                        SUM_DIAG = SUM_DIAG + SC%A(ICOL) 
                        
                  END SELECT
               ENDDO

               DO IW = 1, NWEIGHTS
                  NEIGHBOR(IW) = SC%CELLTYPE(NEIGHBOR(IW))
                  WEIGHTS(IW)  = WEIGHTS(IW)/SUM_DIAG
IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) THEN
WRITE(SCARC_LU,*) '        --->  final data:'
WRITE(SCARC_LU,*) '              NEIGHBOR(',IW,')=',NEIGHBOR(IW)
WRITE(SCARC_LU,*) '              WEIGHTS (',IW,')=',WEIGHTS(IW)
ENDIF
               ENDDO

            ENDIF

            !!! 
            !!! Define corresponding entry in interpolation matrix P by means of the upper weights
            !!!
WRITE(*,*) 'HIER NOCH MAL CHECKEN, NOCH ALTE VERSION, AN NEUE KOMMUNIKATION ANPASSEN !!!'
WRITE(SCARC_LU,*) 'HIER NOCH MAL CHECKEN, NOCH ALTE VERSION, AN NEUE KOMMUNIKATION ANPASSEN !!!'
            SC%P_ROW(IC) = IP
            DO IW = 1, NWEIGHTS
               IF  (NEIGHBOR(IW) /= -1) THEN
                  SC%P_COL(IP) = NEIGHBOR(IW)
                  SC%P(IP)     = WEIGHTS(IW)
IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) THEN
WRITE(SCARC_LU,*) '        --->  prolongation data:'
WRITE(SCARC_LU,*) '              P_COL(',IP,')=',SC%P_COL(ip)
WRITE(SCARC_LU,*) '              P    (',IP,')=',SC%P(ip)
ENDIF
                  IP = IP +1
               ENDIF
            ENDDO
            SC%P_ROW(SC%NCE+1) = IP

         ENDDO INTERNAL_CELL_LOOP_CLASSICAL
         
      ENDDO MESHES_LOOP_CLASSICAL

   !!! -------------------------------------------------------------------------------------------------
   !!! Direct interpolation
   !!! -------------------------------------------------------------------------------------------------
   CASE (NSCARC_INTERPOL_DIRECT)

      MESHES_LOOP_DIRECT: DO NM = NMESHES_MIN, NMESHES_MAX
      
         SC => SCARC(NM)%COMPACT(NL)
         CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_CELLTYPE, NL, 'SETUP_PROLONGATION','CELLTYPE INIT')
         
         IP = 1
WRITE(SCARC_LU,*) '------------------ COARSE', SC%NC
         INTERNAL_CELL_LOOP: DO IC = 1, SC%NC
         
            VALUES   = 0.0_EB
            WEIGHTS  = 0.0_EB
            NEIGHBOR = 0
         
            !!!
            !!! If IC is a coarse cell, its value is taken
            !!!
            IF (SC%CELLTYPE(IC) >= NSCARC_CELLTYPE_COARSE) THEN
         
               NEIGHBOR(1)= SC%CELLTYPE(IC)
               NEIGHBOR(1)= IC
               WEIGHTS(1) = 1.0_EB
               NWEIGHTS   = 1
         

IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) THEN
WRITE(SCARC_LU,*) '------------------ COARSE'
WRITE(SCARC_LU,*) 'IC=',IC,': WEIGHTS=',WEIGHTS(1)
ENDIF

            !!!
            !!! If IC is a fine cell, a mean value must be computed based on surrounding coarse cells
            !!!
            ELSE 
         
               !!! Get main diagonal entry a_ii for that fine cell
               IDIAG = SC%A_ROW(IC)
         
IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) THEN
WRITE(SCARC_LU,*) 'IDIAG=',IDIAG,' A(IDIAG)=', SC%A(IDIAG)
ENDIF

               !!! Select type of fine cell (weakly/strongly coupled)
               SELECT_FPOINT_TYPE: SELECT CASE(SC%CELLTYPE(IC))

                  !!!
                  !!! Strongly coupled fine cell IC
                  !!! approximate IC by weighted sum of surrounding strongly coupled coarse cells JC
                  !!! Note: N_i: all coupled neighboring cells, C_i: all strongly coupled neighboring cells
                  !!!
                  CASE (NSCARC_CELLTYPE_SFINE)
         
                     !!! Compute scaling factor: SCAL = [sum_(k in N_i) a_ik]/[sum_(l in C_i) a_il]
                     SUM_COUPLED = 0.0_EB
                     SUM_CPOINTS = 0.0_EB
         
IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) THEN
WRITE(SCARC_LU,*) '-------------> SFINE:'
WRITE(SCARC_LU,*) 'SC%A_ROW(',IC,')+1=',SC%A_ROW(IC)+1
WRITE(SCARC_LU,*) 'SC%A_ROW(',IC+1,')-1=',SC%A_ROW(IC+1)-1
ENDIF
                     DO ICOL = SC%A_ROW(IC)+1, SC%A_ROW(IC+1)-1
                        JC = SC%A_COL(ICOL)
                        SUM_COUPLED = SUM_COUPLED + SC%A(ICOL)
                        IF (SC%CELLTYPE(JC) > 0) SUM_CPOINTS = SUM_CPOINTS + SC%A(ICOL)
IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) THEN
WRITE(SCARC_LU,'(a,i3,a,f12.6,a,f12.6)') ' --- JC = ',JC,' SUM_COUPLED=',SUM_COUPLED,' SUM_CPOINTS=',SUM_CPOINTS
ENDIF
                     ENDDO
         
                     SCAL = - SUM_COUPLED/SUM_CPOINTS
IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) THEN
WRITE(SCARC_LU,'(a,i3,a,f12.6,2f12.6)') 'IC=',IC,': SCAL=',SCAL, SUM_COUPLED, SUM_CPOINTS
ENDIF
         
                     !!! for each coupling of IC compute interpolation weight SCAL * a_ij/a_ii
                     IW = 1
                     DO ICOL = SC%A_ROW(IC)+1, SC%A_ROW(IC+1)-1
                        JC = SC%A_COL(ICOL)
                        IF (SC%CELLTYPE(JC) > 0 ) THEN
                           NEIGHBOR(IW) = SC%CELLTYPE(JC)
                           NEIGHBOR(IW) = JC
                           WEIGHTS(IW)  = SCAL * SC%A(ICOL)/SC%A(IDIAG)
IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) THEN
WRITE(SCARC_LU,'(a,i3,a,i3,a,i3,i3,i3,3f12.6)') 'IC=',IC,': JC=',JC,': ',SC%CELLTYPE(JC), IW, NEIGHBOR(IW), WEIGHTS(IW), SC%A(ICOL), SC%A(IDIAG)
ENDIF
                           IW = IW +1
                        ENDIF
                     ENDDO
                     NWEIGHTS = IW - 1
IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) THEN
WRITE(SCARC_LU,'(a,i3,a,8i7)') 'IC=',IC,': NEIGHBORS==',NEIGHBOR(1:IW-1)
WRITE(SCARC_LU,'(a,i3,a,8f7.3)') 'IC=',IC,': WEIGHTS   ==',WEIGHTS(1:IW-1)
ENDIF
         

                  !!!
                  !!! Weakly coupled fine cell IC:
                  !!! Determine strongly coupled fine cells JC surrounding IC and, in turn, replace 
                  !!! each of them by a mean value of their surrounding strongly coupled coarse cells
                  !!!
                  CASE (NSCARC_CELLTYPE_WFINE)
         
                     IW = 1                                               ! weights counter
                     DO ICOL = SC%A_ROW(IC)+1, SC%A_ROW(IC+1)-1           ! loop over couplings of weakly coupled fine cell
         
                        !!! regard subdiagonal matrix entries a_ij of weakly coupled fine cell JC to IC
                        JC = SC%A_COL(ICOL)
         
                        !!! Find all surrounding (coupled) points of JC and compute scaling factor 
                        !!! compute scaling factor SCAL = a_ij/a_ii * [sum_(k in N_j) a_jk]/[sum_(l in C_j) ajl] 
                        SUM_COUPLED = 0.0_EB
                        SUM_CPOINTS = 0.0_EB
         
                        DO JCOL = SC%A_ROW(JC)+1, SC%A_ROW(JC+1)-1
                           KC = SC%A_COL(JCOL)
                           SUM_COUPLED = SUM_COUPLED + SC%A(JCOL)
                           IF (SC%CELLTYPE(KC) > 0 ) SUM_CPOINTS = SUM_CPOINTS + SC%A(JCOL)
                        ENDDO
         
                        IF (SUM_CPOINTS == 0.0_EB) THEN
                           WRITE(*,*) 'Error in SCARC_INTERPOLATION on mesh ',NM,' for cell ',IC,' with neighbor ', JC, &
                                      'on level ', NL,' stopping program!'
                           WRITE(*,*) '  SUM_COUPLED =',SUM_COUPLED
                           WRITE(*,*) '  SUM_CPOINTS =',SUM_CPOINTS
                           STOP
                        ENDIF
                        SCAL =  SC%A(ICOL)/SC%A(IDIAG) * SUM_COUPLED/SUM_CPOINTS
                            
IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) THEN
WRITE(SCARC_LU,*) '------------------ WFINE'
WRITE(SCARC_LU,*) 'IC=',IC,': JC=',JC,': SCAL=',SCAL
ENDIF

                        !!! Get diagonal matrix a_jj for point JC
                        JDIAG = SC%A_ROW(JC)
         
                        !!! Compute interpolation weights for all strong coarse cell couplings KC of JC
                        !!! note that a coarse cell KC may be considered several times for different JC's
                        DO JCOL = SC%A_ROW(JC)+1, SC%A_ROW(JC+1)-1
                           KC = SC%A_COL(JCOL)
                           IF (SC%CELLTYPE(KC) > 0) THEN
                             NEIGHBOR(IW) = SC%CELLTYPE(KC)
                             NEIGHBOR(IW) = KC
                             WEIGHTS(IW)  = SCAL * SC%A(JCOL)/SC%A(JDIAG)
IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) THEN
WRITE(SCARC_LU,'(a,i3,a,i3,a,i3,i3,i3,3f12.6)') 'JC=',JC,': KC=',KC,': ',SC%CELLTYPE(KC), IW, NEIGHBOR(IW), WEIGHTS(IW), SC%A(JCOL), SC%A(JDIAG)
ENDIF
                             IW = IW +1
                           ENDIF
                        ENDDO
                           
                     ENDDO
                     NWEIGHTS = IW - 1
IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) THEN
WRITE(SCARC_LU,'(a,i3,a,10i8)') '1: IC=',IC,': NEIGHBORS==',NEIGHBOR(1:NWEIGHTS)
WRITE(SCARC_LU,'(a,i3,a,10f8.2)') '1: IC=',IC,': WEIGHTS   ==',WEIGHTS(1:NWEIGHTS)
ENDIF
         
                     !!! make weights unique (add weights for multiple coarse cells)
                     DO IW0 = 1, NWEIGHTS
                        DO IW = IW0+1, NWEIGHTS
                           IF  (NEIGHBOR(IW0) /= -1 .AND. NEIGHBOR(IW) == NEIGHBOR(IW0)) THEN
                              WEIGHTS(IW0) = WEIGHTS(IW0) + WEIGHTS(IW)
                              WEIGHTS(IW)  = 0.0_EB
                              NEIGHBOR(IW) = -1
                           ENDIF
                        ENDDO
                     ENDDO
IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) THEN
WRITE(SCARC_LU,'(a,i3,a,10i8)') '2: IC=',IC,': NEIGHBORS==',NEIGHBOR(1:NWEIGHTS)
WRITE(SCARC_LU,'(a,i3,a,10f8.2)') '2: IC=',IC,': WEIGHTS   ==',WEIGHTS(1:NWEIGHTS)
ENDIF
         
               END SELECT SELECT_FPOINT_TYPE
         
            ENDIF
         
         
            !!! 
            !!! Define corresponding entry in interpolation matrix P by means of the upper weights
            !!!
            SC%P_ROW(IC) = IP
            DO IW = 1, NWEIGHTS
               IF  (NEIGHBOR(IW) /= -1) THEN
                  SC%P_COL(IP) = NEIGHBOR(IW)
                  SC%P(IP)     = WEIGHTS(IW)
                  IP = IP +1
               ENDIF
            ENDDO
         
         ENDDO INTERNAL_CELL_LOOP

         !GHOST_CELL_LOOP: DO IC = SC%NC+1, SC%NCE
         !   IF (SC%CELLTYPE(IC) <= NSCARC_CELLTYPE_FINE) THEN
         !      SC%P_ROW(IC) = IP            
         !   ELSE
         !      SC%P(IP)     = 1.0_EB
         !      SC%P_COL(IP) = SC%CELLTYPE(IC)
         !      SC%P_ROW(IC) = IP            
         !   ENDIF
         !   IP = IP +1
         !ENDDO GHOST_CELL_LOOP

         SC%P_ROW(SC%NC+1) = IP
         !SC%NCCE = SC%NCC
         
         
         INTERNAL_CELL_LOOP2: DO IC = 1, SC%NC
            DO ICOL = SC%P_ROW(IC), SC%P_ROW(IC+1)-1
              JC = SC%P_COL(ICOL)
              IF (JC <= SC%NC) THEN                  !!! internal cell : get new coarse number
                 SC%P_COL(ICOL) = SC%CELLTYPE(JC)         
              ENDIF
IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) WRITE(SCARC_LU,'(a,i3,a,i3,a,i3)') 'A:PROL2:IC=',IC,': P_COL(',ICOL,')=',SC%P_COL(ICOL)
            ENDDO
         ENDDO INTERNAL_CELL_LOOP2

         IF (NMESHES > 1) THEN
         IF (NL == 1) THEN
            GHOST_CELL_LOOP: DO IW = 1, SC%NW
IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) WRITE(SCARC_LU,*) '------------> IW=',IW
               ICG = SC%GHOST_CELL(IW)
               IF (ICG == 0) CYCLE GHOST_CELL_LOOP
               ICA = SC%ADJACENT_CELL(IW)
               NOM = SC%IJKW(9,IW)
               NXO = MESHES(NOM)%IBAR
               NYO = MESHES(NOM)%JBAR
IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) WRITE(SCARC_LU,*) '---> IW=',IW,ICG, ICA
               DO ICOL = SC%P_ROW(ICA), SC%P_ROW(ICA+1)-1
                 JC = SC%P_COL(ICOL)
                 IF (JC > SC%NC) THEN
                    JCO = (SC%IJKW(12,IW)-1)*NXO*NYO + (SC%IJKW(11,IW)-1)*NXO + SC%IJKW(10,IW)
                    SC%P_COL(ICOL) = - JCO
IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) WRITE(SCARC_LU,'(a,i3,a,i3,a,i3)') 'B:PROL2:IC=',ICA,': P_COL(',ICOL,')=',JCO
                 ENDIF
               ENDDO
            ENDDO GHOST_CELL_LOOP
         ELSE
           WRITE(SCARC_LU,*) 'STILL MISSING3'
         ENDIF
         ENDIF

         
      ENDDO MESHES_LOOP_DIRECT
   

   !!! -------------------------------------------------------------------------------------------------
   !!! Multipass interpolation
   !!! -------------------------------------------------------------------------------------------------
   CASE (NSCARC_INTERPOL_MULTIPASS)

      MESHES_LOOP_MULTIPASS: DO NM = NMESHES_MIN, NMESHES_MAX
         WRITE(*,*) 'Multipass interpolation not yet implemented'
      ENDDO MESHES_LOOP_MULTIPASS
      stop

END SELECT SELECT_INTERPOLATION

IF (TYPE_DEBUG >= NSCARC_DEBUG_MEDIUM) WRITE(SCARC_LU,*) 'SUSI PROLONG2'
CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_PROLONGATION, NL, 'SETUP_PROLONGATION', 'FINAL PROLONG')
IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) WRITE(SCARC_LU,*) 'BEFORE NSCARC_EXCHANGE_WEIGHTS'

IF (NMESHES > 1) CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_WEIGHTS, NL)

IF (TYPE_DEBUG > NSCARC_DEBUG_LESS) WRITE(SCARC_LU,*) 'AFTER NSCARC_EXCHANGE_WEIGHTS'
CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_PROLONGATION, NL, 'SETUP_PROLONGATION', 'FINAL PROLONG2')
   
END SUBROUTINE SCARC_SETUP_PROLONGATION




SUBROUTINE SCARC_SETUP_PROLONGATION3(NL)
INTEGER , INTENT(IN) :: NL
INTEGER  :: NM, ICP, ICP2, IC, IC2, ICOL
TYPE (SCARC_COMPACT_TYPE), POINTER :: SCF

MESHES_LOOP: DO NM = NMESHES_MIN, NMESHES_MAX

   SCF => SCARC(NM)%COMPACT(NL)

   DO IC = SCF%NC+1, SCF%NCEP
      DO ICOL = SCF%P_ROW(IC), SCF%P_ROW(IC+1)-1
      ENDDO
 !           IF (ABS(NNUM)==1) THEN
!
!               ICF = NINT(RECV_BUF(LL))
!               ICC = NINT(RECV_BUF(LL+1))
!               PW  = RECV_BUF(LL+2)
!
!               SC%P_COL(ICOL) = - ICC
!               SC%P(ICOL) = PW
!
!   IF (TYPE_DEBUG > NSCARC_DEBUG_NONE)  THEN
!      WRITE(SCARC_LU,*) '=========== A =============== '
!      WRITE(SCARC_LU,'(a,i3,a,i3,a,i3,a,i3,a,i3,a,i3,a,f12.3,a,i3,a,i3,a,i3,a,f12.3)') &
!           'IW=',IW,': IC0=',IC0,': ICA=',ICA,': ICG=',ICG,': ICF=',ICF,': ICC=',ICC,&
!          ': PW=',PW,': P_COL(',ICOL,')=',SC%P_COL(ICOL),': P(',ICOL,')=',SC%P(ICOL)
!   ENDIF
!
!               DO ICOL0 = SC%P_ROW(ICA), SC%P_ROW(ICA+1)-1
! WRITE(SCARC_LU,*) 'ICA=',ICA,': ICOL0=',ICOL0,SC%P_COL(ICOL0), ICG
!                 IF (SC%P_COL(ICOL0) == - ICF) THEN
!                    SC%P_COL(ICOL0) = -ICC
!   IF (TYPE_DEBUG > NSCARC_DEBUG_NONE)  THEN
!!      write(SCARC_LU,'(a,i3,a,f12.3)') WRITE(SCARC_LU,*) 'ICA: P_COL(',ICOL0,')=',SC%P_COL(ICOL0)
!   ENDIF
!                 ENDIF
!               ENDDO
!
!               IC = ICG
!               LL = LL + 3
!               ICOL = ICOL + 1

   ENDDO

ENDDO MESHES_LOOP

END SUBROUTINE SCARC_SETUP_PROLONGATION3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Define restriction matrix R (currently transpose of prolongation matrix P)
!!!  - In spite of the additinal need for the storing of R, this is done to save computational time
!!!  - during the later matrix transfer operations 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SETUP_RESTRICTION(NL)
INTEGER , INTENT(IN) :: NL
INTEGER  :: NM, ICP, ICP2, IC, IC2
LOGICAL  :: BFIRST
TYPE (SCARC_COMPACT_TYPE), POINTER :: SC

MESHES_LOOP: DO NM = NMESHES_MIN, NMESHES_MAX

   SC => SCARC(NM)%COMPACT(NL)
   
   IC2 = 1
   DO ICP = 1, SC%NCCE
   
      BFIRST = .TRUE.
      DO IC = 1, SC%NCE
         
         ROW_LOOP: DO ICP2 = SC%P_ROW(IC),SC%P_ROW(IC+1)-1
            IF (SC%P_COL(ICP2) == ICP) THEN
               SC%R(IC2) = SC%P(ICP2)
               IF (BFIRST) SC%R_ROW(ICP) = IC2
               SC%R_COL(IC2) = IC
               IC2 = IC2 + 1
               BFIRST = .FALSE.
               EXIT ROW_LOOP
            ENDIF
         ENDDO ROW_LOOP
   
      ENDDO
   
   ENDDO
   SC%R_ROW(SC%NCCE+1)=IC2
   
ENDDO MESHES_LOOP

CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_RESTRICTION, NL, 'SETUP_RESTRICTION', 'FINAL RESTRICT')
END SUBROUTINE SCARC_SETUP_RESTRICTION


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Compute coarser matrix A_coarse by use of interpolation matrices:
!!!
!!!      A_coarse = I_fine^coarse * A_fine * I_coarse_fine
!!! 
!!! Note the different storage techniques (compact storage technique for 
!!! transfer matrices and coarse matrix, banded for finest matrix)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_TRANSFER_MATRIX(NL)
INTEGER , INTENT(IN) :: NL
INTEGER :: NM, IC, ICO, IERR, ICP1, ICP2, IP, ICOL, IDIAG, III
REAL (EB), ALLOCATABLE, DIMENSION(:) :: VAL
REAL:: AUX1, AUX2, PW
LOGICAL:: BONE, BTWO, BDEBUG
TYPE (SCARC_COMPACT_TYPE), POINTER :: SCF, SCC

IERR = 0

MESHES_LOOP: DO NM = NMESHES_MIN, NMESHES_MAX

   SCF => SCARC(NM)%COMPACT(NL)             ! pointer to fine level
   SCC => SCARC(NM)%COMPACT(NL+1)           ! pointer to coarse level
   
WRITE(SCARC_LU,*) 'TRANSFER_MATRIX: SCF%NCCE=',SCF%NCCE
WRITE(SCARC_LU,*) 'HALLO7: NCCE=',SCF%NCCE

   ALLOCATE (VAL(SCF%NCEP), STAT=IERR)
   CALL CHKMEMERR ('SCARC_TRANSFER_MATRIX', 'VAL', IERR)
   VAL = 0.0_EB
   
WRITE(SCARC_LU,*) 'SC%NC  =',SCF%NC
WRITE(SCARC_LU,*) 'SC%NCC =',SCF%NCC
WRITE(SCARC_LU,*) 'SC%NCCE=',SCF%NCCE

   IP  = 1
   ICP1_LOOP: DO ICP1 = 1, SCF%NCCE
WRITE(SCARC_LU,*) '================================================== ICP1=',ICP1
   !DO ICP1 = 1, SCF%NCC
      ICP2_LOOP1: DO ICP2 = 1, SCF%NCEP

         
WRITE(SCARC_LU,*) '------------- ICP2=',ICP2
      !DO ICP2 = 1, SCF%NCC
         AUX2 = 0.0_EB
         DO IC = 1, SCF%NC
!WRITE(SCARC_LU,*) '      ------> IC  =',IC  
   
BONE = (TYPE_DEBUG>NSCARC_DEBUG_LESS).AND.NMESHES==1.AND.(IC==5.OR.IC==10)
BTWO = (TYPE_DEBUG>NSCARC_DEBUG_LESS).AND.NMESHES==2.AND.(IC==5.OR.IC==10)
BDEBUG = BONE.OR.BTWO

            IDIAG = SCF%A_ROW(IC)
            PW = P_WEIGHT(IC, ICP2, NM, NL)
            AUX1  = SCF%A(IDIAG) * PW
   
    IF (BDEBUG) WRITE(SCARC_LU,'(a,i3,a,i3,a,i3,a,4f9.3, 2i4,XX,6i4)') 'IC=',IC,': IDIAG =',IDIAG, ' ICP2=',ICP2,' AUX1=',AUX1, SCF%A(IDIAG), PW, AUX1, SCF%P_ROW(IC), SCF%P_ROW(IC+1)-1, (SCF%P_COL(III),III=SCF%P_ROW(IC),SCF%P_ROW(IC+1)-1)

            COUPLINGS_LOOP: DO ICOL = SCF%A_ROW(IC)+1, SCF%A_ROW(IC+1)-1
               ICO  = SCF%A_COL(ICOL)

               PW = P_WEIGHT(ICO, ICP2, NM, NL)
               AUX1 = AUX1 + SCF%A(ICOL) * PW
  IF (BDEBUG) WRITE(SCARC_LU,'(a,i3,a,i3,a,i3,a,4f9.3,2i4,XX,6i4)') 'IC=',IC,':   ICO =',ICO, ' ICP2=',ICP2,' AUX1=',AUX1, SCF%A(ICOL), PW, AUX2, SCF%P_ROW(ICO), SCF%P_ROW(ICO+1)-1, (SCF%P_COL(III),III=SCF%P_ROW(ICO),SCF%P_ROW(ICO+1)-1)
            ENDDO COUPLINGS_LOOP
   
            AUX2 = AUX2 + P_WEIGHT(IC, ICP1, NM, NL) * AUX1

  IF (BDEBUG) WRITE(SCARC_LU,'(a,3f12.6)') '     --- > AUX2=',AUX2, P_WEIGHT(IC,ICP1,NM,NL), AUX1
    
         ENDDO
         VAL(ICP2) = AUX2
  !IF (BDEBUG) WRITE(SCARC_LU,'(a,i3,a,f12.6)') '----------> VAL(',ICP2,')=',VAL(ICP2)
      ENDDO ICP2_LOOP1
   
      ICP2_LOOP2: DO ICP2 = SCF%NC+1, SCF%NCEP

         
WRITE(SCARC_LU,*) '------------- ICP2=',ICP2
      !DO ICP2 = 1, SCF%NCC
         AUX2 = 0.0_EB
         DO IC = 1, SCF%NC
!WRITE(SCARC_LU,*) '      ------> IC  =',IC  
   
BONE = (TYPE_DEBUG>NSCARC_DEBUG_LESS).AND.NMESHES==1.AND.(IC==5.OR.IC==10)
BTWO = (TYPE_DEBUG>NSCARC_DEBUG_LESS).AND.NMESHES==2.AND.(IC==5.OR.IC==10)
BDEBUG = BONE.OR.BTWO

            IDIAG = SCF%A_ROW(IC)
            PW = P_WEIGHT(IC, ICP2, NM, NL)
            AUX1  = SCF%A(IDIAG) * PW
   
    IF (BDEBUG) WRITE(SCARC_LU,'(a,i3,a,i3,a,i3,a,4f9.3, 2i4,XX,6i4)') 'IC=',IC,': IDIAG =',IDIAG, ' ICP2=',ICP2,' AUX1=',AUX1, SCF%A(IDIAG), PW, AUX1, SCF%P_ROW(IC), SCF%P_ROW(IC+1)-1, (SCF%P_COL(III),III=SCF%P_ROW(IC),SCF%P_ROW(IC+1)-1)

            COUPLINGS_LOOP: DO ICOL = SCF%A_ROW(IC)+1, SCF%A_ROW(IC+1)-1
               ICO  = SCF%A_COL(ICOL)

               PW = P_WEIGHT(ICO, ICP2, NM, NL)
               AUX1 = AUX1 + SCF%A(ICOL) * PW
  IF (BDEBUG) WRITE(SCARC_LU,'(a,i3,a,i3,a,i3,a,4f9.3,2i4,XX,6i4)') 'IC=',IC,':   ICO =',ICO, ' ICP2=',ICP2,' AUX1=',AUX1, SCF%A(ICOL), PW, AUX2, SCF%P_ROW(ICO), SCF%P_ROW(ICO+1)-1, (SCF%P_COL(III),III=SCF%P_ROW(ICO),SCF%P_ROW(ICO+1)-1)
            ENDDO COUPLINGS_LOOP
   
            AUX2 = AUX2 + P_WEIGHT(IC, ICP1, NM, NL) * AUX1

  IF (BDEBUG) WRITE(SCARC_LU,'(a,3f12.6)') '     --- > AUX2=',AUX2, P_WEIGHT(IC,ICP1,NM,NL), AUX1
    
         ENDDO
         VAL(ICP2) = AUX2
  !IF (BDEBUG) WRITE(SCARC_LU,'(a,i3,a,f12.6)') '----------> VAL(',ICP2,')=',VAL(ICP2)
      ENDDO ICP2_LOOP1
   


      !!! analyze new matrix line and store it corresponding to compact storage technique:
      !!! (diagonal entry first)
      SCC%A(IP)       = VAL(ICP1)
      SCC%A_ROW(ICP1) = IP
      SCC%A_COL(IP)   = ICP1
  IF (BDEBUG) WRITE(SCARC_LU,'(a,i3,a,f12.6,a,i3,a,i3,a,i3,a,i3)') &
          '----------> A(',IP,')=',SCC%A(IP),&
                     ' A_ROW(',ICP1,')=',SCC%A_ROW(ICP1),&
                     ' A_COL(',IP,')=',SCC%A_COL(IP)
   
      IP  = IP + 1
      ICP2_LOOP3: DO ICP2 = 1, SCF%NCCE
      !DO ICP2 = 1, SCF%NCC
         IF (ICP2 /= ICP1 .AND. ABS(VAL(ICP2)) >= 1.0E-12_EB) THEN
            SCC%A(IP)     = VAL(ICP2)
            SCC%A_COL(IP) = ICP2
            IP  = IP + 1
         ENDIF
      ENDDO ICP2_LOOP3
   
WRITE(SCARC_LU,'(20f9.3)') (VAL(ICP2),ICP2=1,SCF%NCCE)
   ENDDO ICP1_LOOP

   SCC%A_ROW(SCC%NC+1) = IP
   SCC%NA = IP - 1
   
   DEALLOCATE(VAL)
   
ENDDO MESHES_LOOP

CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_MATRIX, NL+1, 'TRANSFER_MATRIX', 'Next coarser matrix1')
TYPE_MATRIX = NSCARC_MATRIX_SUBDIAG
IF (NMESHES > 1) CALL SCARC_EXCHANGE(NSCARC_EXCHANGE_MATRIX, NL)

CALL SCARC_DEBUG_QUANTITY (NSCARC_DEBUG_MATRIX, NL+1, 'TRANSFER_MATRIX', 'Next coarser matrix2')

stop

END SUBROUTINE SCARC_TRANSFER_MATRIX


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Determine if corresponding coarse point contributes a non-zero interpolation weight or not
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL(EB) FUNCTION P_WEIGHT(IC, ICP, NM, NL)
INTEGER, INTENT(IN):: IC, ICP, NM, NL
INTEGER :: ICOL
REAL(EB) :: VAL
TYPE (SCARC_COMPACT_TYPE), POINTER :: SCF, SCC

SCF => SCARC(NM)%COMPACT(NL)           
SCC => SCARC(NM)%COMPACT(NL+1)           

VAL = 0.0_EB

IF (ICP <= SCF%NCC) THEN

   P_WEIGHT_LOOP1: DO ICOL = SCF%P_ROW(IC), SCF%P_ROW(IC+1)-1
      IF (ICOL > 0) THEN
        IF (SCF%P_COL(ICOL) /= ICP) CYCLE P_WEIGHT_LOOP1
      ELSE
      ENDIF
      VAL = SCF%P(ICOL)
   ENDDO P_WEIGHT_LOOP1

ELSE
ENDIF

P_WEIGHT = VAL
RETURN 
END FUNCTION P_WEIGHT


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Determine if corresponding coarse point contributes a non-zero interpolation weight or not
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL(EB) FUNCTION R_WEIGHT(IC, ICP, NM, NL)
INTEGER, INTENT(IN):: IC, ICP, NM, NL
INTEGER :: ICOL
TYPE (SCARC_COMPACT_TYPE), POINTER :: SC

SC => SCARC(NM)%COMPACT(NL)           

R_WEIGHT = 0.0_EB
R_WEIGHT_LOOP: DO ICOL = SC%R_ROW(ICP), SC%R_ROW(ICP+1)-1
   IF (SC%R_COL(ICOL) == IC) THEN
      R_WEIGHT =  SC%R(ICOL)
      EXIT R_WEIGHT_LOOP
   ENDIF
ENDDO R_WEIGHT_LOOP

RETURN
END FUNCTION R_WEIGHT


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Interface for the call of ScaRC-solver with requested storage technique
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SOLVER
REAL (EB) :: TNOW_SOLVER

TNOW_SOLVER = SECOND()

SELECT_METHOD: SELECT CASE (TYPE_METHOD)

   !!! Krylov method (CG/BICG)
   CASE (NSCARC_METHOD_KRYLOV)

      SELECT_KRYLOV: SELECT CASE (TYPE_KRYLOV)
         CASE (NSCARC_KRYLOV_CG)                           
            CALL SCARC_METHOD_CG  (NSCARC_SCOPE_MAIN, NSCARC_VECTOR_F)
         CASE (NSCARC_KRYLOV_BICG)                         
            CALL SCARC_METHOD_BICG(NSCARC_SCOPE_MAIN, NSCARC_VECTOR_F)
      END SELECT SELECT_KRYLOV

   !!! Multigrid method
   CASE (NSCARC_METHOD_MULTIGRID)

      CALL SCARC_METHOD_MULTIGRID(NSCARC_SCOPE_MAIN, NSCARC_VECTOR_F)

END SELECT SELECT_METHOD

TUSED_SCARC(NSCARC_TIME_SOLVER  ,:)=TUSED_SCARC(NSCARC_TIME_SOLVER  ,:)+SECOND()-TNOW_SOLVER
TUSED_SCARC(NSCARC_TIME_TOTAL,:)=TUSED_SCARC(NSCARC_TIME_TOTAL,:)+SECOND()-TNOW_SOLVER
END SUBROUTINE SCARC_SOLVER


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Set pointer to chosen vector for banded storage technique
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION POINT_TO_BVECTOR(NVECTOR, NM, NL)
REAL(EB), POINTER, DIMENSION(:,:,:) :: POINT_TO_BVECTOR
INTEGER, INTENT(IN):: NVECTOR, NM, NL

SELECT CASE (NVECTOR)
   CASE (NSCARC_VECTOR_X)
      POINT_TO_BVECTOR => SCARC(NM)%BANDED(NL)%X
   CASE (NSCARC_VECTOR_F)
      POINT_TO_BVECTOR => SCARC(NM)%BANDED(NL)%F
   CASE (NSCARC_VECTOR_Y)
      POINT_TO_BVECTOR => SCARC(NM)%BANDED(NL)%Y
   CASE (NSCARC_VECTOR_G)
      POINT_TO_BVECTOR => SCARC(NM)%BANDED(NL)%G
   CASE (NSCARC_VECTOR_W)
      POINT_TO_BVECTOR => SCARC(NM)%BANDED(NL)%W
   CASE (NSCARC_VECTOR_D)
      POINT_TO_BVECTOR => SCARC(NM)%BANDED(NL)%D
   CASE (NSCARC_VECTOR_Z)
      POINT_TO_BVECTOR => SCARC(NM)%BANDED(NL)%Z
   CASE (NSCARC_VECTOR_X2)
      POINT_TO_BVECTOR => SCARC(NM)%BANDED(NL)%X2
   CASE (NSCARC_VECTOR_D2)
      POINT_TO_BVECTOR => SCARC(NM)%BANDED(NL)%D2
   CASE (NSCARC_VECTOR_W2)
      POINT_TO_BVECTOR => SCARC(NM)%BANDED(NL)%W2
   CASE (NSCARC_VECTOR_Y2)
      POINT_TO_BVECTOR => SCARC(NM)%BANDED(NL)%Y2
   CASE (NSCARC_VECTOR_H)
      POINT_TO_BVECTOR => MESHES(NM)%H
   CASE (NSCARC_VECTOR_HS)
      POINT_TO_BVECTOR => MESHES(NM)%HS
END SELECT

RETURN
END FUNCTION POINT_TO_BVECTOR


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Set pointer to chosen vector for banded storage technique
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION POINT_TO_CVECTOR(NVECTOR, NM, NL)
REAL(EB), POINTER, DIMENSION(:) :: POINT_TO_CVECTOR
INTEGER, INTENT(IN):: NVECTOR, NM, NL

SELECT CASE (NVECTOR)
   CASE (NSCARC_VECTOR_X)
      POINT_TO_CVECTOR => SCARC(NM)%COMPACT(NL)%X
   CASE (NSCARC_VECTOR_F)
      POINT_TO_CVECTOR => SCARC(NM)%COMPACT(NL)%F
   CASE (NSCARC_VECTOR_Y)
      POINT_TO_CVECTOR => SCARC(NM)%COMPACT(NL)%Y
   CASE (NSCARC_VECTOR_G)
      POINT_TO_CVECTOR => SCARC(NM)%COMPACT(NL)%G
   CASE (NSCARC_VECTOR_W)
      POINT_TO_CVECTOR => SCARC(NM)%COMPACT(NL)%W
   CASE (NSCARC_VECTOR_D)
      POINT_TO_CVECTOR => SCARC(NM)%COMPACT(NL)%D
   CASE (NSCARC_VECTOR_Z)
      POINT_TO_CVECTOR => SCARC(NM)%COMPACT(NL)%Z
   CASE (NSCARC_VECTOR_X2)
      POINT_TO_CVECTOR => SCARC(NM)%COMPACT(NL)%X2
   CASE (NSCARC_VECTOR_D2)
      POINT_TO_CVECTOR => SCARC(NM)%COMPACT(NL)%D2
   CASE (NSCARC_VECTOR_W2)
      POINT_TO_CVECTOR => SCARC(NM)%COMPACT(NL)%W2
   CASE (NSCARC_VECTOR_Y2)
      POINT_TO_CVECTOR => SCARC(NM)%COMPACT(NL)%Y2
   CASE (NSCARC_VECTOR_MEASURE)
      POINT_TO_CVECTOR => SCARC(NM)%COMPACT(NL)%MEASURE
END SELECT

RETURN
END FUNCTION POINT_TO_CVECTOR

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Set pointer to chosen vector for banded storage technique
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION POINT_TO_CVECTOR_INT(NVECTOR, NM, NL)
INTEGER, POINTER, DIMENSION(:) :: POINT_TO_CVECTOR_INT
INTEGER, INTENT(IN):: NVECTOR, NM, NL

SELECT CASE (NVECTOR)
   CASE (NSCARC_VECTOR_CELLTYPE)
      POINT_TO_CVECTOR_INT => SCARC(NM)%COMPACT(NL)%CELLTYPE
END SELECT

END FUNCTION POINT_TO_CVECTOR_INT


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Compute global matrix-vector product (including data exchange along internal boundaries)
!!! for banded storage technique
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_MATVEC_PRODUCT(NVECTOR1, NVECTOR2, NL)
INTEGER , INTENT(IN):: NVECTOR1, NVECTOR2, NL
REAL(EB), POINTER, DIMENSION(:,:,:) :: VB1, VB2
REAL(EB), POINTER, DIMENSION(:)     :: VC1, VC2
REAL(EB), POINTER, DIMENSION(:,:)   :: AB
REAL(EB), POINTER, DIMENSION(:)     :: AC
INTEGER , POINTER, DIMENSION(:)     :: AC_ROW, AC_COL
INTEGER , POINTER :: NX, NY, NZ, NC
INTEGER :: NM, I, J, K, IC, JC, ICOL


!!!----------------------------------------------------------------------------------------------------
!!! Exchange internal boundary values of vector1 such that the ghost values contain the corresponding
!!! overlapped values of adjacent neighbor
!!!----------------------------------------------------------------------------------------------------
TYPE_VECTOR = NVECTOR1
IF (NMESHES > 1) CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_VECTOR, NL)


!!!----------------------------------------------------------------------------------------------------
!!! Perform global matrix-vector product:
!!! Note: - matrix already contains subdiagonal values from neighbor along internal boundaries
!!!       - if vector1 contains neighboring values, then correct values of global matvec are achieved
!!!----------------------------------------------------------------------------------------------------
SELECT CASE (TYPE_SYSTEM)

   !!! ---------------------------- Banded storage technique ------------------------------------------
   CASE (NSCARC_SYSTEM_BANDED)

      DO NM = NMESHES_MIN, NMESHES_MAX
      
         VB1 => POINT_TO_BVECTOR (NVECTOR1, NM, NL)
         VB2 => POINT_TO_BVECTOR (NVECTOR2, NM, NL)
         
         NX => SCARC(NM)%BANDED(NL)%NX                                     
         NY => SCARC(NM)%BANDED(NL)%NY
         NZ => SCARC(NM)%BANDED(NL)%NZ

         AB => SCARC(NM)%BANDED(NL)%A                                        ! system matrix
      
         SELECT CASE (TYPE_DIMENSION)
            
            CASE (NSCARC_DIMENSION_TWO)
      
               DO K = 1, NZ
                  DO I = 1, NX
                     IC = (K-1) * NX + I
                     VB2(I, 1, K) =   (  AB(IC, ID ) * VB1(I  , 1, K  )   &         ! diagonal component
                                       + AB(IC, ILZ) * VB1(I  , 1, K-1)   &         ! lower z-component
                                       + AB(IC, ILX) * VB1(I-1, 1, K  )   &         ! lower x-component
                                       + AB(IC, IUX) * VB1(I+1, 1, K  )   &         ! upper x-component
                                       + AB(IC, IUZ) * VB1(I  , 1, K+1) )           ! upper z-component
!  IF (NM==1.and.I==16.and.k==16) &
!         WRITE(SCARC_LU,'(2i3,11e11.3)') I,K,VB2(I,1,K),AB(IC,ID),AB(IC,ILZ),AB(IC,ILX),AB(IC,IUX),AB(IC,IUZ),&
!                                 VB1(I,1,K),VB1(I,1,K-1),VB1(I-1,1,K),VB1(I+1,1,K),VB1(I,1,K+1)


                  ENDDO
               ENDDO
         
            CASE (NSCARC_DIMENSION_THREE)
      
               DO K = 1, NZ
                  DO J = 1, NY
                     DO I = 1, NX
                        IC = (K-1) * NX * NY + (J-1) * NX + I
                        VB2(I, J, K) =  ( AB(IC, ID ) * VB1(I  , J  , K  )   &      ! diagonal component
                                      +   AB(IC, ILZ) * VB1(I  , J  , K-1)   &      ! lower z-component
                                      +   AB(IC, ILY) * VB1(I  , J-1, K  )   &      ! lower y-component
                                      +   AB(IC, ILX) * VB1(I-1, J  , K  )   &      ! lower x-component
                                      +   AB(IC, IUX) * VB1(I+1, J  , K  )   &      ! upper x-component
                                      +   AB(IC, IUY) * VB1(I  , J+1, K  )   &      ! upper y-component
                                      +   AB(IC, IUZ) * VB1(I  , J  , K+1) )        ! upper z-component
                     ENDDO
                  ENDDO
               ENDDO
         
         END SELECT
      
      ENDDO 
   
   !!! ---------------------------- Compact storage technique -----------------------------------------
   CASE (NSCARC_SYSTEM_COMPACT)

      DO NM = NMESHES_MIN, NMESHES_MAX
      
         VC1 => POINT_TO_CVECTOR (NVECTOR1, NM, NL)
         VC2 => POINT_TO_CVECTOR (NVECTOR2, NM, NL)
         
         NC     => SCARC(NM)%COMPACT(NL)%NC                                   ! number of cells

         AC     => SCARC(NM)%COMPACT(NL)%A                                    ! system matrix
         AC_ROW => SCARC(NM)%COMPACT(NL)%A_ROW                                ! row pointer
         AC_COL => SCARC(NM)%COMPACT(NL)%A_COL                                ! column pointer
      
         DO IC = 1, NC
                 
            !!! diagonal entry
            ICOL = AC_ROW(IC)                             
            JC   = AC_COL(ICOL)
      
            VC2 (IC) = AC(ICOL)* VC1(JC)
      
            !!! subdiagonal entries
            DO ICOL = AC_ROW(IC)+1, AC_ROW(IC+1)-1          
               JC = AC_COL(ICOL)
               VC2(IC) =  VC2(IC) + AC(ICOL)* VC1(JC)
            ENDDO
      
         ENDDO
      
      ENDDO 
      
END SELECT

END SUBROUTINE SCARC_MATVEC_PRODUCT


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Compute global scalarproductt (including global data exchange) for banded storage technique
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL(EB) FUNCTION SCARC_SCALAR_PRODUCT(NVECTOR1, NVECTOR2, NL)
INTEGER, INTENT(IN):: NVECTOR1, NVECTOR2, NL
REAL(EB), DIMENSION(:,:,:), POINTER ::  VB1, VB2
REAL(EB), DIMENSION(:)    , POINTER ::  VC1, VC2
INTEGER , POINTER :: NX, NY, NZ, NC
INTEGER  :: NM, IERR, NL0, I, J, K, IC
REAL(EB) :: SP_GLOBAL
! MKL code
REAL(EB) :: DDOT
EXTERNAL :: DDOT

!!!----------------------------------------------------------------------------------------------------
!!! Compute local scalarproduct
!!!----------------------------------------------------------------------------------------------------
SELECT CASE (TYPE_SYSTEM)

   !!! ---------------------------- Banded storage technique ------------------------------------------
   CASE (NSCARC_SYSTEM_BANDED)

      DO NM = NMESHES_MIN, NMESHES_MAX
      
         VB1 => POINT_TO_BVECTOR (NVECTOR1, NM, NL)
         VB2 => POINT_TO_BVECTOR (NVECTOR2, NM, NL)
         
         NX  => SCARC(NM)%BANDED(NL)%NX
         NY  => SCARC(NM)%BANDED(NL)%NY
         NZ  => SCARC(NM)%BANDED(NL)%NZ
      
         SP_LOCAL(NM) = 0.0_EB
         DO K = 1, NZ
            DO J = 1, NY
               DO I = 1, NX
                  SP_LOCAL(NM) = SP_LOCAL(NM) + VB1 (I, J, K) * VB2 (I, J, K)
               ENDDO
            ENDDO
         ENDDO
      
      ENDDO

   !!! ---------------------------- Compact storage technique -----------------------------------------
   CASE (NSCARC_SYSTEM_COMPACT)

      DO NM = NMESHES_MIN, NMESHES_MAX
      
         VC1 => POINT_TO_CVECTOR (NVECTOR1, NM, NL)
         VC2 => POINT_TO_CVECTOR (NVECTOR2, NM, NL)
         
         NC  => SCARC(NM)%COMPACT(NL)%NC
      
         IF (SCARC_MKL) THEN
            SP_LOCAL(NM) = DDOT(NC, VC1, 1, VC2, 1)
         ELSE
            SP_LOCAL(NM) = 0.0_EB
            DO IC = 1, NC
               SP_LOCAL(NM) = SP_LOCAL(NM) + VC1(IC) * VC2(IC)
            ENDDO
         ENDIF
      
      ENDDO

END SELECT

!!!----------------------------------------------------------------------------------------------------
!!! get global scalarproduct by a global summation of the local values
!!!----------------------------------------------------------------------------------------------------
IERR = 0
NL0  = NL
SP_GLOBAL   = 0.0_EB
IF (NMESHES>1 .AND. USE_MPI) THEN
   CALL MPI_ALLREDUCE (SP_LOCAL(MYID+1), SP_GLOBAL, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, IERR)
ELSE
   DO NM=1,NMESHES
      SP_GLOBAL = SP_GLOBAL + SP_LOCAL(NM)
   ENDDO
ENDIF

SCARC_SCALAR_PRODUCT = SP_GLOBAL
RETURN

END FUNCTION SCARC_SCALAR_PRODUCT


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Compute global L2-norm (including global data exchange) for banded storage technique
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL(EB) FUNCTION SCARC_L2NORM(NVECTOR1, NL)
INTEGER, INTENT(IN):: NVECTOR1, NL
REAL(EB), DIMENSION(:,:,:), POINTER ::  VB
REAL(EB), DIMENSION(:)    , POINTER ::  VC
INTEGER , POINTER :: NX, NY, NZ, NC
INTEGER  :: NM, IERR, I, J, K, IC
REAL(EB) :: SP_GLOBAL  
REAL(EB) :: DDOT
EXTERNAL :: DDOT

!!!----------------------------------------------------------------------------------------------------
!!! Compute local scalarproduct
!!!----------------------------------------------------------------------------------------------------
SELECT CASE (TYPE_SYSTEM)

   !!! ---------------------------- Banded storage technique ------------------------------------------
   CASE (NSCARC_SYSTEM_BANDED)

      DO NM = NMESHES_MIN, NMESHES_MAX
      
         VB => POINT_TO_BVECTOR (NVECTOR1, NM, NL)

         NX => SCARC(NM)%BANDED(NL)%NX
         NY => SCARC(NM)%BANDED(NL)%NY
         NZ => SCARC(NM)%BANDED(NL)%NZ
      
         SP_LOCAL(NM) = 0.0_EB
         DO K = 1, NZ
            DO J = 1, NY
               DO I = 1, NX
                  SP_LOCAL(NM) = SP_LOCAL(NM) + VB (I, J, K) * VB (I, J, K)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      
   !!! ---------------------------- Compact storage technique -----------------------------------------
   CASE (NSCARC_SYSTEM_COMPACT)

      DO NM = NMESHES_MIN, NMESHES_MAX
      
         VC => POINT_TO_CVECTOR (NVECTOR1, NM, NL)
      
         NC => SCARC(NM)%COMPACT(NL)%NC
      
         IF (SCARC_MKL) THEN
            SP_LOCAL(NM) = DDOT(NC, VC, 1, VC, 1)
         ELSE
            SP_LOCAL(NM) = 0.0_EB
            DO IC = 1, NC
               SP_LOCAL(NM) = SP_LOCAL(NM) + VC(IC) * VC(IC)
            ENDDO
         ENDIF

      ENDDO

END SELECT

!!!----------------------------------------------------------------------------------------------------
!!! get global scalarproduct by a global summation of the local values
!!!----------------------------------------------------------------------------------------------------
IERR = 0
SP_GLOBAL = 0.0_EB
IF (NMESHES>1 .AND. USE_MPI) THEN
   CALL MPI_ALLREDUCE (SP_LOCAL(MYID+1), SP_GLOBAL, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, IERR)
ELSE
   DO NM=1,NMESHES
      SP_GLOBAL = SP_GLOBAL + SP_LOCAL(NM)
   ENDDO
ENDIF
SP_GLOBAL = SQRT (SP_GLOBAL/REAL(NC_GLOBAL(NL), EB))   

SCARC_L2NORM = SP_GLOBAL
RETURN

END FUNCTION SCARC_L2NORM



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Compute linear combination of two vectors for banded storage technique
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_VECTOR_SUM(NVECTOR1, NVECTOR2, SCAL1, SCAL2, NL)
INTEGER , INTENT(IN):: NVECTOR1, NVECTOR2, NL
REAL(EB), INTENT(IN):: SCAL1, SCAL2
REAL(EB), DIMENSION(:,:,:), POINTER ::  VB1, VB2
REAL(EB), DIMENSION(:)    , POINTER ::  VC1, VC2
INTEGER  :: NM
! MKL code
INTEGER , POINTER :: NC
EXTERNAL :: DAXPBY

SELECT CASE (TYPE_SYSTEM)

   !!! ---------------------------- Banded storage technique ------------------------------------------
   CASE (NSCARC_SYSTEM_BANDED)

      DO NM = NMESHES_MIN, NMESHES_MAX

         VB1 => POINT_TO_BVECTOR(NVECTOR1, NM, NL)
         VB2 => POINT_TO_BVECTOR(NVECTOR2, NM, NL)

         VB2 = SCAL1 * VB1 + SCAL2 * VB2

      ENDDO

   !!! ---------------------------- Compact storage technique ------------------------------------------
   CASE (NSCARC_SYSTEM_COMPACT)

      DO NM = NMESHES_MIN, NMESHES_MAX

         VC1 => POINT_TO_CVECTOR(NVECTOR1, NM, NL)
         VC2 => POINT_TO_CVECTOR(NVECTOR2, NM, NL)

         IF (SCARC_MKL) THEN
            NC  => SCARC(NM)%COMPACT(NL)%NC
            CALL DAXPBY(NC, SCAL1, VC1, 1, SCAL2, VC2, 1)
         ELSE
           VC2 = SCAL1 * VC1 + SCAL2 * VC2
         ENDIF

      ENDDO

END SELECT

END SUBROUTINE SCARC_VECTOR_SUM



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Define vector2 to be a scaled copy of vector 1 for banded storage technique
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_VECTOR_COPY(NVECTOR1, NVECTOR2, SCAL1, NL)
INTEGER , INTENT(IN):: NVECTOR1, NVECTOR2, NL
REAL(EB), INTENT(IN):: SCAL1
REAL(EB), DIMENSION(:,:,:), POINTER ::  VB1, VB2
REAL(EB), DIMENSION(:)    , POINTER ::  VC1, VC2
INTEGER  :: NM
!!! MKL code
INTEGER , POINTER :: NC
EXTERNAL :: DCOPY, DSCAL

SELECT CASE (TYPE_SYSTEM)

   !!! ---------------------------- Banded storage technique ------------------------------------------
   CASE (NSCARC_SYSTEM_BANDED)

      DO NM = NMESHES_MIN, NMESHES_MAX

         VB1 => POINT_TO_BVECTOR(NVECTOR1, NM, NL)
         VB2 => POINT_TO_BVECTOR(NVECTOR2, NM, NL)

         VB2 = SCAL1 * VB1 

      ENDDO

   !!! ---------------------------- Compact storage technique ------------------------------------------
   CASE (NSCARC_SYSTEM_COMPACT)

      DO NM = NMESHES_MIN, NMESHES_MAX

         VC1 => POINT_TO_CVECTOR(NVECTOR1, NM, NL)
         VC2 => POINT_TO_CVECTOR(NVECTOR2, NM, NL)

         IF (SCARC_MKL) THEN
            NC  => SCARC(NM)%COMPACT(NL)%NC
            CALL DCOPY(NC, VC1, 1, VC2, 1)
            CALL DSCAL(NC, SCAL1, VC2, 1)
         ELSE
            VC2 = SCAL1 * VC1
         ENDIF

      ENDDO

END SELECT

END SUBROUTINE SCARC_VECTOR_COPY


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Clear vector
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_VECTOR_CLEAR(NVECTOR, NL)
INTEGER , INTENT(IN):: NVECTOR, NL
REAL(EB), DIMENSION(:,:,:), POINTER ::  VB
REAL(EB), DIMENSION(:)    , POINTER ::  VC
INTEGER  :: NM

SELECT CASE (TYPE_SYSTEM)

   !!! ---------------------------- Banded storage technique ------------------------------------------
   CASE (NSCARC_SYSTEM_BANDED)

      DO NM = NMESHES_MIN, NMESHES_MAX
         VB => POINT_TO_BVECTOR(NVECTOR, NM, NL)
         VB =  0.0_EB
      ENDDO

   !!! ---------------------------- Compact storage technique ------------------------------------------
   CASE (NSCARC_SYSTEM_COMPACT)

      DO NM = NMESHES_MIN, NMESHES_MAX
         VC => POINT_TO_CVECTOR(NVECTOR, NM, NL)
         VC =  0.0_EB
      ENDDO

END SELECT

END SUBROUTINE SCARC_VECTOR_CLEAR

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Define vector to a scalar value
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_VECTOR_DEFINE(NVECTOR, SCAL, NL)
INTEGER , INTENT(IN):: NVECTOR, NL
REAL(EB), INTENT(IN):: SCAL
REAL(EB), DIMENSION(:,:,:), POINTER ::  VB
REAL(EB), DIMENSION(:)    , POINTER ::  VC
INTEGER  :: NM

SELECT CASE (TYPE_SYSTEM)

   !!! ---------------------------- Banded storage technique ------------------------------------------
   CASE (NSCARC_SYSTEM_BANDED)

      DO NM = NMESHES_MIN, NMESHES_MAX
         VB => POINT_TO_BVECTOR(NVECTOR, NM, NL)
         VB =  SCAL
      ENDDO

   !!! ---------------------------- Compact storage technique ------------------------------------------
   CASE (NSCARC_SYSTEM_COMPACT)

      DO NM = NMESHES_MIN, NMESHES_MAX
         VC => POINT_TO_CVECTOR(NVECTOR, NM, NL)
         VC =  SCAL
      ENDDO

END SELECT

END SUBROUTINE SCARC_VECTOR_DEFINE



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Perform preconditioning for banded storage technique
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_PRECONDITIONING (NVECTOR1, NVECTOR2, NL)
USE POIS, ONLY: H2CZSS, H3CZSS
INTEGER, INTENT(IN):: NVECTOR1, NVECTOR2, NL
REAL(EB), DIMENSION(:,:,:), POINTER ::  VB1, VB2, FFT
REAL(EB), DIMENSION(:)    , POINTER ::  VC1, VC2
REAL(EB), DIMENSION(:,:)  , POINTER ::  AB
REAL(EB), DIMENSION(:)    , POINTER ::  AC
INTEGER , DIMENSION(:)    , POINTER ::  AC_ROW, AC_COL
INTEGER , POINTER:: NX, NY, NZ, NC
INTEGER  :: NM, I, J, K, IC, ICOL
REAL(EB) :: AUX, OMEGA=1.5_EB
TYPE (MESH_TYPE), POINTER :: M


SELECT CASE (TYPE_SYSTEM)

   !!! ---------------------------- Banded storage technique ------------------------------------------
   CASE (NSCARC_SYSTEM_BANDED)

      SELECT CASE (TYPE_PRECON)
      
         !!!--------------------------------------------------------------------------------------------
         !!! Multigrid preconditioner
         !!!--------------------------------------------------------------------------------------------
         CASE (NSCARC_PRECON_MULTIGRID)
      
            CALL SCARC_METHOD_MULTIGRID (NSCARC_SCOPE_PRECON, NVECTOR2)
      
         !!!--------------------------------------------------------------------------------------------
         !!! Jacobi preconditioner
         !!!--------------------------------------------------------------------------------------------
         CASE (NSCARC_PRECON_JACOBI)
      
            DO NM = NMESHES_MIN, NMESHES_MAX
      
               VB2 => POINT_TO_BVECTOR(NVECTOR2, NM, NL)
         
               NX => SCARC(NM)%BANDED(NL)%NX
               NY => SCARC(NM)%BANDED(NL)%NY
               NZ => SCARC(NM)%BANDED(NL)%NZ
               AB => SCARC(NM)%BANDED(NL)%A
      
               DO K = 1, NZ
                  DO J = 1, NY
                     DO I = 1, NX
                        IC = (K-1) * NX * NY + (J-1) * NX + I
                        VB2 (I, J, K) = VB2 (I, J, K) / AB(IC, ID)
                     ENDDO
                  ENDDO
               ENDDO
      
            ENDDO
         
         !!!--------------------------------------------------------------------------------------------
         !!! SSOR preconditioner
         !!!--------------------------------------------------------------------------------------------
         CASE (NSCARC_PRECON_SSOR)
         
            DO NM = NMESHES_MIN, NMESHES_MAX
            
               VB2 => POINT_TO_BVECTOR(NVECTOR2, NM, NL)
         
               NX => SCARC(NM)%BANDED(NL)%NX
               NY => SCARC(NM)%BANDED(NL)%NY
               NZ => SCARC(NM)%BANDED(NL)%NZ
               NC => SCARC(NM)%BANDED(NL)%NC
               AB => SCARC(NM)%BANDED(NL)%A
         

               SELECT CASE (TYPE_DIMENSION)
               
                  CASE (NSCARC_DIMENSION_TWO)
               
                     DO K = 1, NZ
                        DO I = 1, NX
                           IC = (K-1) * NX + I
                           AUX =    AB(IC,ILZ) * VB2 (I  , 1, K-1)  &
                                  + AB(IC,ILX) * VB2 (I-1, 1, K  )
                           VB2 (I, 1, K) = (VB2(I, 1, K) - AUX*OMEGA) / AB(IC,ID)
                        ENDDO
                     ENDDO
                     DO K = NZ, 1, - 1
                        DO I = NX, 1, - 1
                           IC = (K-1) * NX + I
                           IF (IC==NC) CYCLE
                           AUX =    AB(IC,IUZ) * VB2 (I  , 1, K+1) &
                                  + AB(IC,IUX) * VB2 (I+1, 1, K  )
                           VB2 (I, 1, K) = VB2 (I, 1, K) - AUX * OMEGA / AB(IC,ID)
                        ENDDO
                     ENDDO
               
                  CASE (NSCARC_DIMENSION_THREE)
               
                     DO K = 1, NZ
                        DO J = 1, NY
                           DO I = 1, NX
                              IC = (K-1) * NX * NY + (J-1) * NX + I
                              AUX =    AB(IC,ILZ) * VB2 (I  , J  , K-1) &
                                     + AB(IC,ILY) * VB2 (I  , J-1, K  ) &
                                     + AB(IC,ILX) * VB2 (I-1, J  , K  )
                             VB2 (I, J, K) = (VB2(I, J, K) - AUX * OMEGA) / AB(IC,ID)
                           ENDDO
                        ENDDO
                     ENDDO
                     DO K = NZ, 1, - 1
                        DO J = NY, 1, - 1
                           DO I = NX, 1, - 1
                              IC = (K-1) * NX * NY + (J-1) * NX + I
                              IF (IC==NC) CYCLE
                              AUX =    AB(IC,IUZ) * VB2 (I  , J  , K+1) &
                                     + AB(IC,IUY) * VB2 (I  , J+1, K  ) &
                                     + AB(IC,IUX) * VB2 (I+1, J  , K  )
                              VB2 (I, J, K) = VB2 (I, J, K) - AUX * OMEGA / AB(IC,ID)
                           ENDDO
                        ENDDO
                     ENDDO
               END SELECT
         
            ENDDO
         
         !!!--------------------------------------------------------------------------------------------
         !!! GSTRIX preconditioner
         !!!--------------------------------------------------------------------------------------------
         !CASE (NSCARC_PRECON_GSTRIX)
         !   CALL SCARC_GSTRIX_BANDED(V2, ...)
            
         !!!--------------------------------------------------------------------------------------------
         !!! FFT preconditioner
         !!!--------------------------------------------------------------------------------------------
         CASE (NSCARC_PRECON_FFT)
            
            DO NM = NMESHES_MIN, NMESHES_MAX
            
               VB1 => POINT_TO_BVECTOR(NVECTOR1, NM, NL)
               VB2 => POINT_TO_BVECTOR(NVECTOR2, NM, NL)
         
               M   => MESHES(NM)
               FFT => SCARC(NM)%PRECON(NLEVEL_MIN)%FFT
               
               FFT(1:M%IBAR, 1:M%JBAR, 1:M%KBAR) = VB1(1:M%IBAR, 1:M%JBAR, 1:M%KBAR)
               SELECT CASE (TYPE_DIMENSION)
                  CASE (NSCARC_DIMENSION_TWO)
                     CALL H2CZSS (M%BXS, M%BXF, M%BZS, M%BZF, &
                                  M%IBAR+1, FFT, M%POIS_PTB, M%SAVE1, M%WORK, M%HX)
                  CASE (NSCARC_DIMENSION_THREE)
                     CALL H3CZSS (M%BXS, M%BXF, M%BYS, M%BYF, M%BZS, M%BZF, &
                                  M%IBAR+1, M%JBAR+1, FFT, M%POIS_PTB, M%SAVE1, M%WORK, M%HX)
               END SELECT
               VB2(1:M%IBAR, 1:M%JBAR, 1:M%KBAR) = FFT(1:M%IBAR, 1:M%JBAR, 1:M%KBAR)
               
            ENDDO
                     
      END SELECT
      

   !!! ---------------------------- Compact storage technique ------------------------------------------
   CASE (NSCARC_SYSTEM_COMPACT)

      SELECT CASE (TYPE_PRECON)
      
         !!!--------------------------------------------------------------------------------------------
         !!! Multigrid preconditioner
         !!!--------------------------------------------------------------------------------------------
         CASE (NSCARC_PRECON_MULTIGRID)
      
            CALL SCARC_METHOD_MULTIGRID (NSCARC_SCOPE_PRECON, NVECTOR2)
      
      
         !!!--------------------------------------------------------------------------------------------
         !!! Jacobi preconditioner
         !!!--------------------------------------------------------------------------------------------
         CASE (NSCARC_PRECON_JACOBI)
      
            DO NM = NMESHES_MIN, NMESHES_MAX
      
               VC2 => POINT_TO_CVECTOR(NVECTOR2, NM, NL)
         
               NC     => SCARC(NM)%COMPACT(NL)%NC
               AC     => SCARC(NM)%COMPACT(NL)%A
               AC_ROW => SCARC(NM)%COMPACT(NL)%A_ROW
      
               DO IC = 1, NC
                  VC2 (IC) = VC2 (IC) / AC (AC_ROW(IC))
!WRITE(SCARC_LU,'(a,i4,a,f12.6,a,i4,a,f12.6)') 'JACOBI: VC2(',IC,')=',VC2(IC),' : A(',IC,')=',AC(AC_ROW(IC))
               ENDDO
      
            ENDDO
         
         !!!--------------------------------------------------------------------------------------------
         !!! SSOR preconditioner
         !!!--------------------------------------------------------------------------------------------
         CASE (NSCARC_PRECON_SSOR)
         
            DO NM = NMESHES_MIN, NMESHES_MAX
            
               VC2 => POINT_TO_CVECTOR(NVECTOR2, NM, NL)
         
               NC     => SCARC(NM)%COMPACT(NL)%NC
               AC     => SCARC(NM)%COMPACT(NL)%A
               AC_ROW => SCARC(NM)%COMPACT(NL)%A_ROW
               AC_COL => SCARC(NM)%COMPACT(NL)%A_COL
         
               !!! use only matrix superdiagonals
               FORWARD_CELL_LOOP: DO IC = 1, NC
               
                  AUX = 0.0_EB
                  LOWER_DIAG_LOOP: DO ICOL = AC_ROW(IC)+1, AC_ROW(IC+1)-1
                     IF (AC_COL(ICOL) >= IC) EXIT LOWER_DIAG_LOOP
                     AUX = AUX + AC(ICOL) * VC2(AC_COL(ICOL))
                  ENDDO LOWER_DIAG_LOOP
                  VC2(IC) = (VC2(IC) - AUX * OMEGA) / AC(AC_ROW(IC))
               
               ENDDO FORWARD_CELL_LOOP
               
               !!! use only matrix subdiagonals
               BACKWARD_CELL_LOOP: DO IC = NC-1, 1, -1
                  
                  AUX = 0.0_EB
                  UPPER_DIAG_LOOP: DO ICOL = AC_ROW(IC)+1, AC_ROW(IC+1)-1
                     IF (AC_COL(ICOL) <= IC) CYCLE
                     AUX = AUX + AC(ICOL) * VC2(AC_COL(ICOL))
                  ENDDO UPPER_DIAG_LOOP
                  VC2(IC) = VC2(IC) - AUX * OMEGA / AC(AC_ROW(IC))
               
               ENDDO BACKWARD_CELL_LOOP
         
            ENDDO
            
         !!!--------------------------------------------------------------------------------------------
         !!! GSTRIX preconditioner
         !!!--------------------------------------------------------------------------------------------
         !CASE (NSCARC_PRECON_GSTRIX)
         !   CALL SCARC_GSTRIX_BANDED(V2, ...)
            
         !!!--------------------------------------------------------------------------------------------
         !!! FFT preconditioner
         !!!--------------------------------------------------------------------------------------------
         CASE (NSCARC_PRECON_FFT)
            
            DO NM = NMESHES_MIN, NMESHES_MAX
            
               VC1 => POINT_TO_CVECTOR(NVECTOR1, NM, NL)
               VC2 => POINT_TO_CVECTOR(NVECTOR2, NM, NL)
         
               M   => MESHES(NM)
               FFT => SCARC(NM)%PRECON(NLEVEL_MIN)%FFT
               
               DO K = 1, M%KBAR
                  DO J = 1, M%JBAR
                     DO I = 1, M%IBAR
                        IC = (K-1) * M%IBAR * M%JBAR + (J-1) * M%IBAR + I
                        FFT(I, J, K) = VC1(IC)
                     ENDDO
                  ENDDO
               ENDDO
               SELECT CASE (TYPE_DIMENSION)
                  CASE (NSCARC_DIMENSION_TWO)
                     CALL H2CZSS (M%BXS, M%BXF, M%BZS, M%BZF,&
                                  M%IBAR+1, FFT, M%POIS_PTB, M%SAVE1, M%WORK, M%HX)
                  CASE (NSCARC_DIMENSION_THREE)
                     CALL H3CZSS (M%BXS, M%BXF, M%BYS, M%BYF, M%BZS, M%BZF, &
                                  M%IBAR+1, M%JBAR+1, FFT, M%POIS_PTB, M%SAVE1, M%WORK, M%HX)
               END SELECT
               DO K = 1, M%KBAR
                  DO J = 1, M%JBAR
                     DO I = 1, M%IBAR
                        IC = (K-1) * M%IBAR * M%JBAR + (J-1) * M%IBAR + I
                        VC2(IC) = FFT(I, J, K)
                     ENDDO
                  ENDDO
               ENDDO
      
            ENDDO
      
      END SELECT

END SELECT

END SUBROUTINE SCARC_PRECONDITIONING 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Perform global CG-method based on global possion-matrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_METHOD_CG(NSCOPE, NVECTOR)
INTEGER, INTENT(IN) :: NSCOPE, NVECTOR
INTEGER   :: NL 
INTEGER   :: ITE, NIT, ISTATE, NPRECON0
REAL (EB) :: SIGMA0, SIGMA1, ALPHA0, GAMMA0
REAL (EB) :: RES, RESIN, EPS
REAL (EB) :: TNOW_KRYLOV
CHARACTER(30) :: CROUTINE = 'null'

!!!----------------------------------------------------------------------------------------------------
!!! Initialization:
!!!   - Set environment variables and define working level
!!!   - Initialize solution, right hand side vector and auxiliary vectors
!!!   - Define iterations parameters
!!!----------------------------------------------------------------------------------------------------
TNOW_KRYLOV = SECOND()
TYPE_SCOPE  = NSCOPE
NPRECON0 = TYPE_PRECON
!IF (TYPE_SCOPE == NSCARC_SCOPE_COARSE) TYPE_PRECON= NSCARC_PRECON_JACOBI

EPS = SCARC_KRYLOV_ACCURACY
NIT = SCARC_KRYLOV_ITERATIONS

CALL SCARC_SETUP_ENVIRONMENT(CROUTINE, NSCOPE, NVECTOR, NL)
CALL SCARC_SETUP_SOLVER(NL)

!CALL SCARC_VECTOR_DEFINE(VEC_X, 0.0_EB, NLEVEL_MAX)
!CALL SCARC_VECTOR_DEFINE(VEC_F, 1.0_EB, NLEVEL_MAX)

!!!----------------------------------------------------------------------------------------------------
!!! Compute initial residual and perform initial preconditioning
!!!----------------------------------------------------------------------------------------------------
CALL SCARC_MATVEC_PRODUCT (VEC_X, VEC_W, NL)                                !  W := A*X
CALL SCARC_VECTOR_SUM     (VEC_F, VEC_W, -1.0_EB, 1.0_EB, NL)               !  W := W - F

CALL SCARC_DEBUG_VECTOR (VEC_X, 'SCARC_METHOD_CG', 'X INIT')
CALL SCARC_DEBUG_VECTOR (VEC_F, 'SCARC_METHOD_CG', 'F INIT')
CALL SCARC_DEBUG_VECTOR (VEC_W, 'SCARC_METHOD_CG', 'W INIT')

RESIN = SCARC_L2NORM (VEC_W, NL)                                            !  RESIN := ||W||

CALL SCARC_CONVERGENCE_INFO(RESIN, 0, NL, CROUTINE)

CALL SCARC_DEBUG_VECTOR (VEC_G, 'SCARC_METHOD_CG', 'G INIT')

CALL SCARC_VECTOR_COPY     (VEC_W, VEC_G, 1.0_EB, NL)                       !  G := W
CALL SCARC_PRECONDITIONING (VEC_W, VEC_G, NL)                               !  G := PRECON(W)

CALL SCARC_DEBUG_VECTOR (VEC_G, 'SCARC_METHOD_CG', 'G PRECON')

SIGMA0 = SCARC_SCALAR_PRODUCT(VEC_W, VEC_G, NL)                             !  SIGMA0 := (W,G)

CALL SCARC_VECTOR_COPY (VEC_G, VEC_D, -1.0_EB, NL)                          !  D := -G

!CALL SCARC_DEBUG_VECTOR (VEC_D, 'SCARC_METHOD_CG', 'D COPY')
!!!----------------------------------------------------------------------------------------------------
!!! Perform conjugate gradient looping
!!!----------------------------------------------------------------------------------------------------
CG_LOOP: DO ITE = 1, NIT
 
   CALL SCARC_MATVEC_PRODUCT (VEC_D, VEC_Y, NL)                             !  Y := A*D

CALL SCARC_DEBUG_VECTOR (VEC_Y, 'SCARC_METHOD_CG', 'Y MATVEC')

   ALPHA0 = SCARC_SCALAR_PRODUCT (VEC_D, VEC_Y, NL)                         !  ALPHA0 := (D,Y)
   ALPHA0 = SIGMA0/ALPHA0                                                 

   CALL SCARC_VECTOR_SUM (VEC_D, VEC_X, ALPHA0, 1.0_EB, NL)                 !  X := ALPHA0*D + X
   CALL SCARC_VECTOR_SUM (VEC_Y, VEC_W, ALPHA0, 1.0_EB, NL)                 !  W := ALPHA0*Y + W

CALL SCARC_DEBUG_VECTOR (VEC_X, 'SCARC_METHOD_CG', 'X VECTORSUM')
CALL SCARC_DEBUG_VECTOR (VEC_W, 'SCARC_METHOD_CG', 'W VECTORSUM')

   RES = SCARC_L2NORM (VEC_W, NL)                                           !  RES := ||W||
 
   ISTATE = SCARC_CONVERGENCE_STATE (RESIN, RES, EPS, ITE, NL, CROUTINE)    !  RES < TOL ??
   IF (ISTATE /= NSCARC_STATE_PROCEED) EXIT CG_LOOP
 
   CALL SCARC_VECTOR_COPY     (VEC_W, VEC_G, 1.0_EB, NL)                    !  G := W
   CALL SCARC_PRECONDITIONING (VEC_W, VEC_G, NL)                            !  G := PRECON(W)

   SIGMA1 = SCARC_SCALAR_PRODUCT (VEC_W, VEC_G, NL)                         !  SIGMA1 := (W,G)
   GAMMA0 = SIGMA1/SIGMA0                                                   
   SIGMA0 = SIGMA1                                                         

   CALL SCARC_VECTOR_SUM (VEC_G, VEC_D, -1.0_EB, GAMMA0, NL)                !  D := -G + GAMMA0*D

ENDDO CG_LOOP
 
!!!----------------------------------------------------------------------------------------------------
!!! Determine convergence rate and print corresponding information
!!! In case of CG as main solver:
!!!   - Transfer ScaRC solution vector X to FDS pressure vector 
!!!   - Set ghost cell values along external boundaries
!!!   - Exchange values along internal boundaries 
!!!----------------------------------------------------------------------------------------------------
CALL SCARC_CONVERGENCE_RATE(RESIN, RES, ITE, ISTATE, CROUTINE)

IF (TYPE_SCOPE == NSCARC_SCOPE_MAIN) THEN
   CALL SCARC_TERMINATE_SOLVER  (NLEVEL_MIN)
   CALL SCARC_UPDATE_GHOSTCELLS (NLEVEL_MIN)
ENDIF

CALL SCARC_DEBUG_VECTOR (VEC_X, 'SCARC_METHOD_CG', 'X FINAL')
TYPE_PRECON = NPRECON0

TUSED_SCARC(NSCARC_TIME_KRYLOV,:)=TUSED_SCARC(NSCARC_TIME_KRYLOV,:)+SECOND()-TNOW_KRYLOV
TUSED_SCARC(NSCARC_TIME_TOTAL ,:)=TUSED_SCARC(NSCARC_TIME_TOTAL ,:)+SECOND()-TNOW_KRYLOV
END SUBROUTINE SCARC_METHOD_CG


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Perform global BICGstab-method based on global possion-matrix - banded storage technique
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_METHOD_BICG(NSCOPE, NVECTOR)
INTEGER, INTENT(IN) :: NSCOPE, NVECTOR
INTEGER   :: NL = NSCARC_LEVEL_NONE
INTEGER   :: ITE, NIT, ISTATE
REAL (EB) :: ALPHA0, ALPHA1, ALPHA2, RHO0, RHO1, DTHETA, DBETA
REAL (EB) :: RES, RESIN, EPS
REAL (EB) :: TNOW_KRYLOV
CHARACTER(30) :: CROUTINE = 'null'

!!!----------------------------------------------------------------------------------------------------
!!! Initialization
!!!   - Set environment variables and define working level
!!!   - Initialize solution, right hand side vector and auxiliary vectors
!!!   - Define iterations parameters
!!!----------------------------------------------------------------------------------------------------
TNOW_KRYLOV = SECOND()
TYPE_SCOPE  = NSCOPE

EPS    = SCARC_KRYLOV_ACCURACY
NIT    = SCARC_KRYLOV_ITERATIONS

ALPHA0 = 1.0_EB
RHO0   = 1.0_EB
RHO1   = 0.0_EB
DBETA  = 0.0_EB
DTHETA = 1.0_EB
 
CALL SCARC_SETUP_ENVIRONMENT(CROUTINE, NSCOPE, NVECTOR, NL)
CALL SCARC_SETUP_SOLVER(NL)
   
!!!----------------------------------------------------------------------------------------------------
!!! Compute initial defect and perform (double) initial preconditioning
!!!----------------------------------------------------------------------------------------------------
!CALL SCARC_VECTOR_COPY    (VEC_F, VEC_W, 1.0_EB, NL)                         !  W := F
!CALL SCARC_PRECONDITIONING (VEC_W, VEC_W, NL)                                !  W := PRECON(W)
CALL SCARC_MATVEC_PRODUCT  (VEC_X, VEC_W, NL)                                 !  W := A*X
CALL SCARC_VECTOR_SUM      (VEC_F, VEC_W, 1.0_EB, -1.0_EB, NL)                !  W := F - W
CALL SCARC_PRECONDITIONING (VEC_W, VEC_W, NL)                                 !  W := PRECON(W)

RESIN = SCARC_L2NORM (VEC_W, NL)                                              !  RESIN := ||W||
CALL SCARC_CONVERGENCE_INFO (RESIN, 0, NL, CROUTINE)

CALL SCARC_VECTOR_COPY (VEC_W, VEC_G, 1.0_EB, NL)                             !  G := W
   
!!!----------------------------------------------------------------------------------------------------
!!! Perform bi-conjugate gradient looping:
!!!----------------------------------------------------------------------------------------------------
BICG_LOOP: DO ITE = 1, NIT

   RHO1  = SCARC_SCALAR_PRODUCT (VEC_G, VEC_W, NL)                            ! RHO1 := (G,W)
   DBETA = (RHO1*DTHETA)/(RHO0*ALPHA0)                                     
   RHO0  = RHO1

   CALL SCARC_VECTOR_SUM      (VEC_W, VEC_Z, 1.0_EB       , DBETA , NL)       ! Z := W + DBETA*Z
   CALL SCARC_VECTOR_SUM      (VEC_Y, VEC_Z, -DBETA*ALPHA0, 1.0_EB, NL)       ! Z := -DBETA*ALPHA0*Y + Z
   CALL SCARC_MATVEC_PRODUCT  (VEC_Z, VEC_Y, NL)                              ! Y := A*Z
   CALL SCARC_PRECONDITIONING (VEC_Y, VEC_Y, NL)                              ! Z := PRECON(Z)

   DTHETA = SCARC_SCALAR_PRODUCT (VEC_G, VEC_Y, NL)                           ! DTHETA := (G,Y)
   DTHETA = RHO1/DTHETA

   CALL SCARC_VECTOR_SUM      (VEC_Y, VEC_W, -DTHETA, 1.0_EB, NL)             ! W := -DTHETA*Y + W
   CALL SCARC_MATVEC_PRODUCT  (VEC_W, VEC_D, NL)                              ! D := A*W
   CALL SCARC_PRECONDITIONING (VEC_D, VEC_D, NL)                              ! D := PRECON(D)
   
   ALPHA1 = SCARC_SCALAR_PRODUCT (VEC_D, VEC_W, NL)                           ! ALPHA1 := (D,W)
   ALPHA2 = SCARC_SCALAR_PRODUCT (VEC_D, VEC_D, NL)                           ! ALPHA2 := (D,D)
   ALPHA0 = ALPHA1/ALPHA2

   CALL SCARC_VECTOR_SUM (VEC_Z, VEC_X,  DTHETA, 1.0_EB, NL)                  ! X :=  DTHETA*Z + X
   CALL SCARC_VECTOR_SUM (VEC_W, VEC_X,  ALPHA0, 1.0_EB, NL)                  ! X :=  ALPHA0*W + X
   CALL SCARC_VECTOR_SUM (VEC_D, VEC_W, -ALPHA0, 1.0_EB, NL)                  ! W := -ALPHA0*D + W

   RES = SCARC_L2NORM (VEC_W, NL)                                             ! RES := ||W||

   ISTATE = SCARC_CONVERGENCE_STATE(RESIN, RES, EPS, ITE, NL, CROUTINE)       ! RES < TOL ???
   IF (ISTATE /= NSCARC_STATE_PROCEED) EXIT BICG_LOOP
 
ENDDO BICG_LOOP

!!!----------------------------------------------------------------------------------------------------
!!! Determine convergence rate and print corresponding information
!!! In case of BICG as main solver:
!!!   - Transfer ScaRC solution vector X to FDS pressure vector 
!!!   - Set ghost cell values along external boundaries
!!!   - Exchange values along internal boundaries 
!!!----------------------------------------------------------------------------------------------------
CALL SCARC_CONVERGENCE_RATE(RESIN, RES, ITE, ISTATE, CROUTINE)

IF (TYPE_SCOPE == NSCARC_SCOPE_MAIN) THEN
   CALL SCARC_TERMINATE_SOLVER(NLEVEL_MIN)
   CALL SCARC_UPDATE_GHOSTCELLS(NLEVEL_MIN)
ENDIF

TUSED_SCARC(NSCARC_TIME_KRYLOV,:)=TUSED_SCARC(NSCARC_TIME_KRYLOV,:)+SECOND()-TNOW_KRYLOV
TUSED_SCARC(NSCARC_TIME_TOTAL ,:)=TUSED_SCARC(NSCARC_TIME_TOTAL ,:)+SECOND()-TNOW_KRYLOV
END SUBROUTINE SCARC_METHOD_BICG


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Perform geometric multigrid method based on global possion-matrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_METHOD_MULTIGRID(NSCOPE, NVECTOR)
INTEGER, INTENT(IN) :: NSCOPE, NVECTOR
INTEGER   :: NL = NSCARC_LEVEL_NONE
INTEGER   :: ITE, NIT, ISTATE, ICYCLE
REAL (EB) :: RES, RESIN, EPS
REAL (EB) :: TNOW_MULTIGRID
CHARACTER(30):: CROUTINE='null'

!!!----------------------------------------------------------------------------------------------------
!!! General initialization:
!!!   - Set environment variables and define working level
!!!   - Initialize solution, right hand side vector 
!!!   - Define iterations parameters (NL is set to finest level)
!!!----------------------------------------------------------------------------------------------------
TNOW_MULTIGRID = SECOND()
TYPE_SCOPE     = NSCOPE

EPS = SCARC_MULTIGRID_ACCURACY
NIT = SCARC_MULTIGRID_ITERATIONS

CALL SCARC_SETUP_ENVIRONMENT(CROUTINE, NSCOPE, NVECTOR, NL)
CALL SCARC_SETUP_SOLVER(NL)

!!!----------------------------------------------------------------------------------------------------
!!! Compute initial defect:  RESIN := || F - A*X ||
!!!   - Initialize cycle counts for MG-iteration
!!!   - Perform initial matrix-vector product on finest level
!!!   - calculate norm of initial residual on finest level
!!!----------------------------------------------------------------------------------------------------
CALL SCARC_MATVEC_PRODUCT (VEC_X, VEC_D, NL)                                           !  D := A*X
CALL SCARC_VECTOR_SUM     (VEC_F, VEC_D, 1.0_EB, -1.0_EB, NL)                          !  D := F - D 

ICYCLE = SCARC_CYCLE_CONTROL(NSCARC_CYCLE_SETUP, NL)
RESIN  = SCARC_L2NORM (VEC_D, NL)                                                      !  RESIN := ||D||
IF (TYPE_DEBUG >=NSCARC_DEBUG_NONE.AND.MYID==0) &
   WRITE(0,'(a,i3,a,e14.5,a,e14.5)') ' MG-Iteration  =',0,': Residuum=',RESIN

CALL SCARC_CONVERGENCE_INFO(RESIN, 0, NL, CROUTINE)

!!!----------------------------------------------------------------------------------------------------
!!! Perform multigrid-looping (start each iteration on finest level) 
!!!----------------------------------------------------------------------------------------------------
MULTIGRID_LOOP: DO ITE = 1, NIT
 
   NL = NLEVEL_MIN
   ICYCLE = SCARC_CYCLE_CONTROL(NSCARC_CYCLE_RESET, NL)

   CYCLE_LOOP: DO WHILE (ICYCLE /= NSCARC_CYCLE_EXIT)

      !!! presmoothing  (smoothing/restriction till coarsest level is reached)
      PRESMOOTHING_LOOP: DO WHILE (NL < NLEVEL_MAX)

         CALL SCARC_SMOOTHING    (VEC_X, VEC_F, VEC_D, NSCARC_CYCLE_PRESMOOTH, NL)     ! D_fine   := smooth(defect)
         CALL SCARC_RESTRICTION  (VEC_D, VEC_F, NL, NL+1)                              ! F_coarse := rest(D_fine)
         CALL SCARC_VECTOR_CLEAR (VEC_X, NL+1)                                         ! X_coarse := 0.0
         NL = NL + 1                                                                   ! set coarser level

      ENDDO PRESMOOTHING_LOOP

      !!! coarse grid solver
      !CALL SCARC_VECTOR_DEFINE(VEC_F, 1.0_EB, NLEVEL_MAX)
CALL SCARC_DEBUG_LEVEL(VEC_F, 'MULTIGRID', 'F before coarse', NLEVEL_MAX)
      SELECT CASE (TYPE_COARSE)
         CASE (NSCARC_COARSE_ITERATIVE)
            CALL SCARC_METHOD_CG (NSCARC_SCOPE_COARSE, VEC_F)                          ! X_coarse := exact_sol(.)
         CASE (NSCARC_COARSE_DIRECT)
            CALL SCARC_METHOD_GE (VEC_X, VEC_F)
      END SELECT
      TYPE_SCOPE = NSCOPE
CALL SCARC_DEBUG_LEVEL(VEC_F, 'MULTIGRID', 'F after coarse', NLEVEL_MAX)
CALL SCARC_DEBUG_LEVEL(VEC_X, 'MULTIGRID', 'X after coarse', NLEVEL_MAX)

      !!! postsmoothing (smoothing/restriction till finest level is reached again)
      POSTSMOOTHING_LOOP: DO WHILE (NL > NLEVEL_MIN) 

         NL=NL-1
         CALL SCARC_PROLONGATION (VEC_X, VEC_D, NL+1, NL)                             ! D_fine := prol(X_coarse)
         CALL SCARC_VECTOR_SUM   (VEC_D, VEC_X, 1.0_EB, 1.0_EB, NL)                   ! X_fine := D_fine + X_fine
         CALL SCARC_SMOOTHING    (VEC_X, VEC_F, VEC_D, NSCARC_CYCLE_POSTSMOOTH, NL)   ! D_fine := smooth(defect)

         ICYCLE = SCARC_CYCLE_CONTROL(NSCARC_CYCLE_PROCEED, NL)                       ! perform requested cycle
         IF (ICYCLE /= NSCARC_CYCLE_POSTSMOOTH) CYCLE CYCLE_LOOP

      ENDDO POSTSMOOTHING_LOOP

   ENDDO CYCLE_LOOP

   IF (NL /= NLEVEL_MIN) THEN
      WRITE(*,*) 'ERROR in SCARC_MULTIGRID, wrong level ', NL
      STOP
   ENDIF

   !!!-------------------------------------------------------------------------------------------------
   !!! Compute norm of new residual on finest level and  leave loop correspondingly
   !!!-------------------------------------------------------------------------------------------------
   CALL SCARC_MATVEC_PRODUCT (VEC_X, VEC_D, NL)                                       ! D := A*X
   CALL SCARC_VECTOR_SUM     (VEC_F, VEC_D, 1.0_EB, -1.0_EB, NL)                      ! D := F - D

   RES = SCARC_L2NORM (VEC_D, NL)                                                     ! RES := ||D||

   ISTATE = SCARC_CONVERGENCE_STATE(RESIN, RES, EPS, ITE, NL, CROUTINE)               ! convergence ?
IF (TYPE_DEBUG >=NSCARC_DEBUG_NONE.AND.MYID==0) &
   WRITE(0,'(a,i3,a,e14.5,a,e14.5)') ' MG-Iteration  =',ITE,': Residuum=',SCARC_RESIDUAL
   IF (ISTATE /= NSCARC_STATE_PROCEED) EXIT MULTIGRID_LOOP
 
ENDDO MULTIGRID_LOOP

!!!----------------------------------------------------------------------------------------------------
!!! Determine convergence rate and print corresponding information
!!! In case of MG as main solver:
!!!   - Transfer ScaRC solution vector X to FDS pressure vector 
!!!   - Set ghost cell values along external boundaries 
!!!   - Exchange values along internal boundaries (consistency!)
!!!----------------------------------------------------------------------------------------------------
CALL SCARC_CONVERGENCE_RATE(RESIN, RES, ITE, ISTATE, CROUTINE)
IF (TYPE_DEBUG >=NSCARC_DEBUG_NONE.AND.MYID==0) &
   WRITE(0,'(a,e14.5)') '                                        ---->  Konvergenzrate=',SCARC_CAPPA

IF (TYPE_SCOPE == NSCARC_SCOPE_MAIN) THEN
   CALL SCARC_TERMINATE_SOLVER(NLEVEL_MIN)
   CALL SCARC_UPDATE_GHOSTCELLS(NLEVEL_MIN)
ENDIF
CALL SCARC_DEBUG_VECTOR (VEC_X, 'SCARC_METHOD_MULTIGRID', 'X FINAL')
 
TUSED_SCARC(NSCARC_TIME_MULTIGRID,:)=TUSED_SCARC(NSCARC_TIME_MULTIGRID,:)+SECOND()-TNOW_MULTIGRID
TUSED_SCARC(NSCARC_TIME_TOTAL    ,:)=TUSED_SCARC(NSCARC_TIME_TOTAL    ,:)+SECOND()-TNOW_MULTIGRID
END SUBROUTINE SCARC_METHOD_MULTIGRID


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Control multigrid cycling (F/V/W)
!!! Note: NLEVEL_MIN corresponds to finest level, NLEVEL_MAX to coarsest level
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
INTEGER FUNCTION SCARC_CYCLE_CONTROL(NSCOPE, NL)
INTEGER, INTENT(IN) :: NSCOPE, NL
INTEGER :: NM, NL0, ICYCLE

SELECT CASE (NSCOPE)

   !!!-------------------------------------------------------------------------------------------------
   !!! initialize cycle counts at beginning of multigrid method
   !!!-------------------------------------------------------------------------------------------------
   CASE (NSCARC_CYCLE_SETUP)

      DO NM = NMESHES_MIN, NMESHES_MAX
         SCARC(NM)%CYCLE_COUNT(2,NL)=1
         DO NL0 = NLEVEL_MIN+1, NLEVEL_MAX - 1
            IF (TYPE_CYCLE==NSCARC_CYCLE_F) THEN
               SCARC(NM)%CYCLE_COUNT(2,NL0)=2
            ELSE
               SCARC(NM)%CYCLE_COUNT(2,NL0)=TYPE_CYCLE
            ENDIF
         ENDDO
      ENDDO

      ICYCLE = NSCARC_CYCLE_NEXT
      
   !!!-------------------------------------------------------------------------------------------------
   !!! reset cycle counts at beginning of each new multigrid iteration
   !!!-------------------------------------------------------------------------------------------------
   CASE (NSCARC_CYCLE_RESET)

      DO NM = NMESHES_MIN, NMESHES_MAX
         DO NL0 = NLEVEL_MIN, NLEVEL_MAX
            SCARC(NM)%CYCLE_COUNT(1,NL0)=SCARC(NM)%CYCLE_COUNT(2,NL0)
         ENDDO
      ENDDO
      ICYCLE = NSCARC_CYCLE_NEXT

   !!!-------------------------------------------------------------------------------------------------
   !!! determine where to proceed with cycling
   !!!-------------------------------------------------------------------------------------------------
   CASE (NSCARC_CYCLE_PROCEED)

      DO NM = NMESHES_MIN, NMESHES_MAX

         SCARC(NM)%CYCLE_COUNT(1,NL)=SCARC(NM)%CYCLE_COUNT(1,NL)-1

         IF (SCARC(NM)%CYCLE_COUNT(1,NL)==0) THEN
            IF (TYPE_CYCLE==NSCARC_CYCLE_F) THEN
               SCARC(NM)%CYCLE_COUNT(1,NL)=1
            ELSE
               SCARC(NM)%CYCLE_COUNT(1,NL)=SCARC(NM)%CYCLE_COUNT(2,NL)
            ENDIF
            IF (NL == NLEVEL_MIN) THEN
               ICYCLE = NSCARC_CYCLE_EXIT
            ELSE
               ICYCLE = NSCARC_CYCLE_POSTSMOOTH
            ENDIF
         ELSE
            IF (NL == NLEVEL_MIN) THEN
               ICYCLE = NSCARC_CYCLE_EXIT
            ELSE
               ICYCLE = NSCARC_CYCLE_NEXT
            ENDIF
         ENDIF

      ENDDO

END SELECT

SCARC_CYCLE_CONTROL = ICYCLE
RETURN

END FUNCTION SCARC_CYCLE_CONTROL


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Perform smoothing - banded storage technique
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SMOOTHING(VEC0_X, VEC0_F, VEC0_D, NTYPE, NL)
INTEGER , INTENT(IN) :: VEC0_X, VEC0_F, VEC0_D, NTYPE, NL
INTEGER :: ITE, NIT, ISTATE
REAL(EB):: RES, RESIN, EPS, OMEGA
REAL(EB):: TNOW_SMOOTH
LOGICAL :: BMATVEC, BL2NORM
CHARACTER(30) :: CROUTINE

TNOW_SMOOTH = SECOND()

!!!----------------------------------------------------------------------------------------------------
!!! Initialization
!!!----------------------------------------------------------------------------------------------------
SELECT CASE (NTYPE)
   CASE (NSCARC_CYCLE_PRESMOOTH)
      CROUTINE = 'SCARC_PRESMOOTHER'
   CASE (NSCARC_CYCLE_POSTSMOOTH)
      CROUTINE = 'SCARC_POSTSMOOTHER'
END SELECT
BL2NORM  = .TRUE.
IF (NTYPE == NSCARC_CYCLE_PRESMOOTH.AND.NL==1) THEN
   BMATVEC = .FALSE.
ELSE
   BMATVEC = .TRUE.
ENDIF
   BMATVEC = .TRUE.

NIT      = SCARC_SMOOTH_ITERATIONS
EPS      = SCARC_SMOOTH_ACCURACY
OMEGA    = SCARC_SMOOTH_OMEGA
ITE      = 0

!!!----------------------------------------------------------------------------------------------------
!!! Calculate initial defect on level NL (only if BMATVEC = .TRUE.)
!!! Because initial vector is set to zero, this defect corresponds to F
!!!----------------------------------------------------------------------------------------------------
IF (BMATVEC) THEN
   CALL SCARC_MATVEC_PRODUCT (VEC0_X, VEC0_D, NL)                             !  D := A*X
   CALL SCARC_VECTOR_SUM     (VEC0_F, VEC0_D, 1.0_EB, -1.0_EB, NL)            !  D := F - D
ENDIF

IF (BL2NORM.AND.BMATVEC) THEN
   RESIN = SCARC_L2NORM (VEC0_D, NL)                                          !  RESIN := ||D||
   CALL SCARC_CONVERGENCE_INFO(RESIN, ITE, NL, CROUTINE)
ELSE
   RESIN = SCARC_RESIDUAL
ENDIF

!!!----------------------------------------------------------------------------------------------------
!!! Smoothing loop
!!!----------------------------------------------------------------------------------------------------
SMOOTH_LOOP: DO ITE=1, NIT
 
   CALL SCARC_PRECONDITIONING (VEC0_D, VEC0_D, NL)                            !  D := PRECON (D)
   CALL SCARC_VECTOR_SUM      (VEC0_D, VEC0_X, OMEGA, 1.0_EB, NL)             !  X := OMEGA*D + X
   CALL SCARC_MATVEC_PRODUCT  (VEC0_X, VEC0_D, NL)                            !  D := A*X
   CALL SCARC_VECTOR_SUM      (VEC0_F, VEC0_D, 1.0_EB, -1.0_EB, NL)           !  D := F - D

   IF (BL2NORM.OR.ITE==NIT) THEN
      RES    = SCARC_L2NORM (VEC0_D, NL)                                      !  RES := ||D||
      ISTATE = SCARC_CONVERGENCE_STATE(RESIN, RES, EPS, ITE, NL, CROUTINE)
      IF (ISTATE /= NSCARC_STATE_PROCEED) EXIT SMOOTH_LOOP                    !  RES < TOL ?
   ENDIF

ENDDO SMOOTH_LOOP

!!!----------------------------------------------------------------------------------------------------
!!! Determine convergence rate and print corresponding information
!!!----------------------------------------------------------------------------------------------------
IF (BL2NORM) CALL SCARC_CONVERGENCE_RATE(RESIN, RES, ITE, ISTATE, CROUTINE)

TUSED_SCARC(NSCARC_TIME_SMOOTH,:)=TUSED_SCARC(NSCARC_TIME_SMOOTH,:)+SECOND()-TNOW_SMOOTH
TUSED_SCARC(NSCARC_TIME_TOTAL ,:)=TUSED_SCARC(NSCARC_TIME_TOTAL ,:)+SECOND()-TNOW_SMOOTH
END SUBROUTINE SCARC_SMOOTHING


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Gaussian elimination for coarse grid solution
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_METHOD_GE(NVECTORX, NVECTORF)
INTEGER, INTENT(IN) :: NVECTORX, NVECTORF
INTEGER :: NM, IC, IX, IY, IZ, IOFFSET, IERR, DISP
REAL(EB), POINTER, DIMENSION(:,:,:) :: VBX, VBF
REAL(EB), POINTER, DIMENSION(:)     :: VCX, VCF
REAL (EB) :: TNOW_COARSE
TYPE (SCARC_TYPE), POINTER :: SM
TYPE (SCARC_BANDED_TYPE),  POINTER :: SB
TYPE (SCARC_COMPACT_TYPE), POINTER :: SC

TNOW_COARSE = SECOND()

IERR=0

SM => SCARC(NMASTER)

IF (TYPE_DEBUG >= NSCARC_DEBUG_LESS) WRITE(SCARC_LU,*) 'CALLING SCARC_METHOD_GE', USE_MPI

!!!----------------------------------------------------------------------------------------------------
!!! Parallel version
!!!----------------------------------------------------------------------------------------------------
IF (USE_MPI) THEN

   DO NM = 1, NMESHES
      SM%COUNTS1(NM-1) = NC_COARSE(NM)
      SM%DISPLS1(NM-1) = SM%OFFSET(NM)
   ENDDO

IF (TYPE_DEBUG >= NSCARC_DEBUG_LESS) WRITE(SCARC_LU,*) 'CALLING SCARC_METHOD_GE: USE_MPI'

   DO NM = NMESHES_MIN, NMESHES_MAX
   
      SELECT CASE (TYPE_SYSTEM)
   
         CASE (NSCARC_SYSTEM_BANDED)
   
            SB  => SCARC(NM)%BANDED(NLEVEL_MAX)
            VBF => POINT_TO_BVECTOR (NVECTORF, NM, NLEVEL_MAX)
   
            IOFFSET = SM%OFFSET(NM)
   
            DISP = SM%DISPLS1(MYID)+1
            SELECT CASE (TYPE_DIMENSION)
   
               CASE (NSCARC_DIMENSION_TWO)
   
                  DO IZ = 1, SB%NZ
                     DO IX = 1, SB%NX
                        IC = (IZ-1)*SB%NX + IX + IOFFSET
                        SM%X_BUF (DISP) = VBF(IX, 1, IZ)
                        DISP = DISP + 1
                     ENDDO
                  ENDDO
   
               CASE (NSCARC_DIMENSION_THREE)
   
                  DO IZ = 1, SB%NZ
                     DO IY = 1, SB%NY
                        DO IX = 1, SB%NX
                           IC = (IZ-1)*SB%NX*SB%NY + (IY-1)*SB%NX + IX + IOFFSET
                           SM%X_BUF (DISP) = VBF(IX, IY, IZ)
                           DISP = DISP + 1
                        ENDDO
                     ENDDO
                  ENDDO
   
            END SELECT

            DISP = SM%DISPLS1(MYID)+1
            CALL MPI_ALLGATHERV(SM%X_BUF(DISP),SM%COUNTS1(MYID),MPI_DOUBLE_PRECISION, &
                                SM%X_COARSE,SM%COUNTS1,SM%DISPLS1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,IERR)

IF (TYPE_DEBUG >= NSCARC_DEBUG_LESS) WRITE(SCARC_LU,*) 'X_COARSE (A): after'
IF (TYPE_DEBUG >= NSCARC_DEBUG_LESS) WRITE(SCARC_LU,'(5f12.5)') (SM%X_COARSE(1:4))

   
         CASE (NSCARC_SYSTEM_COMPACT)
   
            SC => SCARC(NM)%COMPACT(NLEVEL_MAX)
            IOFFSET = SM%OFFSET(NM)
   
            DO IC = 1, SC%NZ
               SM%X_COARSE (IC) = VCF(IC)
            ENDDO
   
            WRITE(*,*) 'SCARC_METHOD_GE, COMPACT, NOT FINISHED YET'
            STOP
      END SELECT
   
   ENDDO
   
   IF (TYPE_DEBUG > NSCARC_DEBUG_MEDIUM) THEN
      WRITE(SCARC_LU,*) 'F_COARSE'
      WRITE(SCARC_LU,'(4e12.4)') (SM%X_COARSE(IX), IX=1, NC_COARSE0)
   ENDIF
   
   !IF (MYID+1 == NMASTER) THEN
      CALL DGETRS('N', NC_COARSE0, 1, SM%A_COARSE, NC_COARSE0, SM%PIVOT, &
                  SM%X_COARSE, NC_COARSE0, IERR)
   !ENDIF
   
   IF (TYPE_DEBUG >= NSCARC_DEBUG_LESS) THEN
      WRITE(SCARC_LU,*) 'X_COARSE'
      WRITE(SCARC_LU,'(4e12.4)') (SM%X_COARSE(IX), IX=1, NC_COARSE0)
   ENDIF
   
   DO NM = NMESHES_MIN, NMESHES_MAX
   
      SELECT CASE (TYPE_SYSTEM)
   
         CASE (NSCARC_SYSTEM_BANDED)
   
            SB => SCARC(NM)%BANDED(NLEVEL_MAX)
            IOFFSET = SM%OFFSET(NM)
   
            VBX => POINT_TO_BVECTOR (NVECTORX, NM, NLEVEL_MAX)
   
            SELECT CASE (TYPE_DIMENSION)
   
               CASE (NSCARC_DIMENSION_TWO)
   
                  DO IZ = 1, SB%NZ
                     DO IX = 1, SB%NX
                        IC = (IZ-1)*SB%NX + IX + IOFFSET
                        VBX(IX, 1, IZ) = SM%X_COARSE (IC) 
                     ENDDO
                  ENDDO
   
               CASE (NSCARC_DIMENSION_THREE)
   
                  DO IZ = 1, SB%NZ
                     DO IY = 1, SB%NY
                        DO IX = 1, SB%NX
                           IC = (IZ-1)*SB%NX*SB%NY + (IY-1)*SB%NX + IX + IOFFSET
                           VBX(IX, IY, IZ) = SM%X_COARSE (IC)  
                        ENDDO
                     ENDDO
                  ENDDO
   
            END SELECT
   
         CASE (NSCARC_SYSTEM_COMPACT)
   
            SC => SCARC(NM)%COMPACT(NLEVEL_MAX)
            IOFFSET = SM%OFFSET(NM)
   
            VBX => POINT_TO_BVECTOR (NVECTORX, NM, NLEVEL_MAX)
   
            DO IC = 1, SC%NC
               VCX(IC) = SM%X_COARSE (IC) 
            ENDDO
   
      END SELECT
   
   ENDDO

!!!----------------------------------------------------------------------------------------------------
!!! Serial version
!!!----------------------------------------------------------------------------------------------------
ELSE

   DO NM = NMESHES_MIN, NMESHES_MAX
   
      SELECT CASE (TYPE_SYSTEM)
   
         CASE (NSCARC_SYSTEM_BANDED)
   
            SB => SCARC(NM)%BANDED(NLEVEL_MAX)
   
            VBF => POINT_TO_BVECTOR (NVECTORF, NM, NLEVEL_MAX)
   
            IOFFSET = SM%OFFSET(NM)
   
            SELECT CASE (TYPE_DIMENSION)
   
               CASE (NSCARC_DIMENSION_TWO)
   
                  DO IZ = 1, SB%NZ
                     DO IX = 1, SB%NX
                           IC = (IZ-1)*SB%NX + IX + IOFFSET
                           SM%X_COARSE (IC) = VBF(IX, 1, IZ)
                     ENDDO
                  ENDDO
   
               CASE (NSCARC_DIMENSION_THREE)
   
                  DO IZ = 1, SB%NZ
                     DO IY = 1, SB%NY
                        DO IX = 1, SB%NX
                           IC = (IZ-1)*SB%NX*SB%NY + (IY-1)*SB%NX + IX + IOFFSET
                           SM%X_COARSE (IC) = VBF(IX, IY, IZ)
                        ENDDO
                     ENDDO
                  ENDDO
   
            END SELECT
   
         CASE (NSCARC_SYSTEM_COMPACT)
   
            SC => SCARC(NM)%COMPACT(NLEVEL_MAX)
            IOFFSET = SM%OFFSET(NM)
   
            DO IC = 1, SC%NZ
               SM%X_COARSE (IC) = VCF(IC)
            ENDDO
   
      END SELECT
   
   ENDDO
   
   IF (TYPE_DEBUG >= NSCARC_DEBUG_MEDIUM) THEN
      WRITE(SCARC_LU,*) 'F_COARSE'
      WRITE(SCARC_LU,'(4e12.4)') (SM%X_COARSE(IX), IX=1, NC_COARSE0)
   ENDIF
   
   IF (MYID+1 == NMASTER) THEN
      CALL DGETRS('N', NC_COARSE0, 1, SM%A_COARSE, NC_COARSE0, SM%PIVOT, &
                  SM%X_COARSE, NC_COARSE0, IERR)
   ENDIF
   
   IF (TYPE_DEBUG >= NSCARC_DEBUG_MEDIUM) THEN
      WRITE(SCARC_LU,*) 'X_COARSE'
      WRITE(SCARC_LU,'(4e12.4)') (SM%X_COARSE(IX), IX=1, NC_COARSE0)
   ENDIF
   
   DO NM = NMESHES_MIN, NMESHES_MAX
   
      SELECT CASE (TYPE_SYSTEM)
   
         CASE (NSCARC_SYSTEM_BANDED)
   
            SB => SCARC(NM)%BANDED(NLEVEL_MAX)
            IOFFSET = SM%OFFSET(NM)
   
            VBX => POINT_TO_BVECTOR (NVECTORX, NM, NLEVEL_MAX)
   
            SELECT CASE (TYPE_DIMENSION)
   
               CASE (NSCARC_DIMENSION_TWO)
   
                  DO IZ = 1, SB%NZ
                     DO IX = 1, SB%NX
                        IC = (IZ-1)*SB%NX + IX + IOFFSET
                        VBX(IX, 1, IZ) = SM%X_COARSE (IC) 
                     ENDDO
                  ENDDO
   
               CASE (NSCARC_DIMENSION_THREE)
   
                  DO IZ = 1, SB%NZ
                     DO IY = 1, SB%NY
                        DO IX = 1, SB%NX
                           IC = (IZ-1)*SB%NX*SB%NY + (IY-1)*SB%NX + IX + IOFFSET
                           VBX(IX, IY, IZ) = SM%X_COARSE (IC)  
                        ENDDO
                     ENDDO
                  ENDDO
   
            END SELECT
   
         CASE (NSCARC_SYSTEM_COMPACT)
   
            SC => SCARC(NM)%COMPACT(NLEVEL_MAX)
            IOFFSET = SM%OFFSET(NM)
   
            VBX => POINT_TO_BVECTOR (NVECTORX, NM, NLEVEL_MAX)
   
            DO IC = 1, SC%NC
               VCX(IC) = SM%X_COARSE (IC) 
            ENDDO
   
      END SELECT
   
   ENDDO

ENDIF

TUSED_SCARC(NSCARC_TIME_COARSE,:)=TUSED_SCARC(NSCARC_TIME_COARSE,:)+SECOND()-TNOW_COARSE
TUSED_SCARC(NSCARC_TIME_TOTAL ,:)=TUSED_SCARC(NSCARC_TIME_TOTAL ,:)+SECOND()-TNOW_COARSE
END SUBROUTINE SCARC_METHOD_GE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Initialize GSTRIX preconditioner/smoother
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SETUP_GSTRIX(NM, NL)
INTEGER, INTENT(IN) :: NM, NL
INTEGER :: IC, IERR
TYPE (SCARC_PRECON_TYPE), POINTER :: SP
TYPE (SCARC_BANDED_TYPE), POINTER :: SB

IERR=0

SB => SCARC(NM)%BANDED(NL)
SP => SCARC(NM)%PRECON(NL)

SELECT CASE (TYPE_DIMENSION)

   !!!----------------------------------------------------------------------------------------------------
   !!! 2D-version
   !!!----------------------------------------------------------------------------------------------------
   CASE (NSCARC_DIMENSION_TWO)

      ALLOCATE (SP%MDX(1:SB%NC), STAT=IERR)
      CALL CHKMEMERR ('SCARC_PRECON_INIT_GSTRIX', 'SP%MDX', IERR)
      
      ALLOCATE (SP%UDX(1:SB%NC), STAT=IERR)
      CALL CHKMEMERR ('SCARC_PRECON_INIT_GSTRIX', 'SP%UDX', IERR)
      
      ALLOCATE (SP%LDX(1:SB%NC), STAT=IERR)
      CALL CHKMEMERR ('SCARC_PRECON_INIT_GSTRIX', 'SP%LDX', IERR)
      
      ALLOCATE (SP%MDZ(1:SB%NC), STAT=IERR)
      CALL CHKMEMERR ('SCARC_PRECON_INIT_GSTRIX', 'SP%MDZ', IERR)
      
      ALLOCATE (SP%DWORK(1:SB%NC), STAT=IERR)
      CALL CHKMEMERR ('SCARC_PRECON_INIT_GSTRIX', 'SP%DWORK', IERR)
      
      !!! save different diagonals of matrix ('D' corresponds to main, 'L' to lower, 'U' to upper)
      IF (TYPE_SYSTEM == NSCARC_SYSTEM_BANDED) THEN
         DO IC = 1, SB%NC
            SP%MDX(IC) = SB%A(IC,ID)         ! main  diagonal in x-direction
            SP%MDZ(IC) = SB%A(IC,ILZ)        ! main  diagonal in z-direction
            SP%LDX(IC) = SB%A(IC,ILX)        ! lower diagonal in x-direction
            SP%UDX(IC) = SB%A(IC,IUX)        ! upper diagonal in x-direction
         ENDDO
      ELSE
         WRITE(*,*) 'HIER NOCH ANPASSEN !!'
         WRITE(SCARC_LU,*) 'HIER NOCH ANPASSEN !!'
         DO IC = 1, SB%NC
   !          SP%MDX(IC) = SC%A(IC,ID)         ! main  diagonal in x-direction
   !         SP%MDZ(IC) = SC%A(IC,ILZ)        ! main  diagonal in z-direction
   !         SP%LDX(IC) = SC%A(IC,ILX)        ! lower diagonal in x-direction
   !         SP%UDX(IC) = SC%A(IC,IUX)        ! upper diagonal in x-direction
         ENDDO
      ENDIF
       
      !!! perform LU-factorization of matrix AG according to banded storage technique
      DO IC = 2, SB%NC
         SP%LDX (IC) = SP%LDX(IC) / SP%MDX(IC-1)
         SP%MDX (IC) = SP%MDX(IC) - SP%LDX(IC) * SP%UDX(IC-1)
      ENDDO
       
      !!! replace diagonal values diag(i) by 1/diag(i) and multiply UDX with 1/diag(i)
      DO IC = 1, SB%NC
         SP%MDX (IC) = 1.0_EB / SP%MDX(IC)
         SP%UDX (IC) = SP%MDX(IC) * SP%UDX(IC)
      ENDDO
      
   !!!----------------------------------------------------------------------------------------------------
   !!! 3D-version
   !!!----------------------------------------------------------------------------------------------------
   CASE (NSCARC_DIMENSION_THREE)
   
      ALLOCATE (SP%MDX(1:SB%NC), STAT=IERR)
      CALL CHKMEMERR ('SCARC_PRECON_INIT_GSTRIX3D', 'SP%MDX', IERR)
      
      ALLOCATE (SP%UDX(1:SB%NC), STAT=IERR)
      CALL CHKMEMERR ('SCARC_PRECON_INIT_GSTRIX3D', 'SP%UDX', IERR)
      
      ALLOCATE (SP%LDX(1:SB%NC), STAT=IERR)
      CALL CHKMEMERR ('SCARC_PRECON_INIT_GSTRIX3D', 'SP%LDX', IERR)
      
      ALLOCATE (SP%MDY(1:SB%NC), STAT=IERR)
      CALL CHKMEMERR ('SCARC_PRECON_INIT_GSTRIX3D', 'SP%MDY', IERR)
      
      ALLOCATE (SP%MDZ(1:SB%NC), STAT=IERR)
      CALL CHKMEMERR ('SCARC_PRECON_INIT_GSTRIX3D', 'SP%MDZ', IERR)
      
      ALLOCATE (SP%DWORK(1:SB%NC), STAT=IERR)
      CALL CHKMEMERR ('SCARC_PRECON_INIT_GSTRIX3D', 'SP%DWORK', IERR)
      
      !!! save different diagonals of matrix ('D' corresponds to main, 'L' to lower, 'U' to upper)
      IF (TYPE_SYSTEM == NSCARC_SYSTEM_BANDED) THEN
         DO IC = 1, SB%NC
            SP%MDX(IC) = SB%A(IC,ID)        ! main  diagonal in x-direction
            SP%MDZ(IC) = SB%A(IC,ILZ)       ! main  diagonal in z-direction
            SP%MDY(IC) = SB%A(IC,ILY)       ! main  diagonal in y-direction
            SP%LDX(IC) = SB%A(IC,ILX)       ! lower diagonal in x-direction
            SP%UDX(IC) = SB%A(IC,IUX)       ! upper diagonal in x-direction
         ENDDO
      ELSE
         WRITE(*,*) 'HIER NOCH ANPASSEN !!'
         WRITE(SCARC_LU,*) 'HIER NOCH ANPASSEN !!'
         DO IC = 1, SB%NC
   !         SP%MDX(IC) = SC%A(IC,ID)        ! main  diagonal in x-direction
   !         SP%MDZ(IC) = SC%A(IC,ILZ)       ! main  diagonal in z-direction
   !         SP%MDY(IC) = SC%A(IC,ILY)       ! main  diagonal in y-direction
   !         SP%LDX(IC) = SC%A(IC,ILX)       ! lower diagonal in x-direction
   !         SP%UDX(IC) = SC%A(IC,IUX)       ! upper diagonal in x-direction
         ENDDO
      ENDIF
       
      !!! perform LU-factorization of matrix AG according to banded storage technique
      DO IC = 2, SB%NC
         SP%LDX (IC) = SP%LDX(IC) / SP%MDX(IC-1)
         SP%MDX (IC) = SP%MDX(IC) - SP%LDX(IC) * SP%UDX(IC-1)
      ENDDO
       
      !!! replace diagonal values diag(i) by 1/diag(i) and multiply UDX with 1/diag(i)
      DO IC = 1, SB%NC
         SP%MDX (IC) = 1.0_EB / SP%MDX(IC)
         SP%UDX (IC) = SP%MDX(IC) * SP%UDX(IC)
      ENDDO
   
END SELECT

END SUBROUTINE SCARC_SETUP_GSTRIX
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! GSTRIX preconditioner/smoother
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_GSTRIX (NVECTOR, NL)
INTEGER, INTENT(IN) :: NVECTOR, NL
REAL(EB), POINTER, DIMENSION(:,:,:) :: VB
REAL(EB), POINTER, DIMENSION(:)     :: VC
INTEGER , POINTER :: NX, NY, NZ
INTEGER :: NM, I, J, K, IC, JC
TYPE (SCARC_PRECON_TYPE) , POINTER :: SP
TYPE (SCARC_BANDED_TYPE) , POINTER :: SB
TYPE (SCARC_COMPACT_TYPE), POINTER :: SC


SELECT_SYSTEM: SELECT CASE (TYPE_SYSTEM)
   
   !!!-------------------------------------------------------------------------------------------------
   !!! Banded system
   !!!-------------------------------------------------------------------------------------------------
   CASE (NSCARC_SYSTEM_BANDED)

      DO NM = NMESHES_MIN, NMESHES_MAX
      
         VB => POINT_TO_BVECTOR (NVECTOR, NM, NL)

         SB => SCARC(NM)%BANDED(NL)
         SP => SCARC(NM)%PRECON(NL)
         
         NX => SB%NX
         NY => SB%NY
         NZ => SB%NZ
         WRITE(SCARC_LU,*) 'ACHTUNG, NOCH ANPASSEN!'
         WRITE(*,*) 'ACHTUNG, NOCH ANPASSEN!'
         
         SELECT_BANDED_DIMENSION: SELECT CASE (TYPE_DIMENSION)
         
            !!!--------------- 2D -----------------
            CASE (NSCARC_DIMENSION_TWO)
         
               ! backward elimination of first NX unkowns (may be solved by tridiagonal system)
               DO I=2,NX
                  VB(I,1,1) = VB(I,1,1)-SP%LDX(I-1)*VB(I-1,1,1)
               ENDDO
               DO I=1,NX
                  VB(I,1,1) = VB(I,1,1)*SP%MDX(I)
               ENDDO
               DO I=NX-1,1,-1
                  VB(I,1,1) = VB(I,1,1)-SP%UDX(I)*VB(I+1,1,1)
               ENDDO
                 
               
               ! backward elimination of following unknowns (here the subdiagonal in z-direction must be taken into account)
               DO K=2,NZ
               
                  IC = (K-1)*NX + I
               
                  ! eliminate subdiagonal in z-direction into right hand side (related solution component already known)
                  DO I=1,NX
                     VB(I,1,K) = VB(I,1,K) - SP%MDZ(IC)*VB(I,1,K-1) 
                  ENDDO
               
                  ! perform elimination of matrix lines corresponding to K
                  DO I=2,NX
                     VB(I,1,K) = VB(I,1,K)-SP%LDX(IC-1)*VB(I-1,1,K)
                  ENDDO
                  DO I=1,NX
                     VB(I,1,K) = VB(I,1,K)*SP%MDX(IC)
                  ENDDO
                  DO I=NX-1,1,-1
                     VB(I,1,K) = VB(I,1,K)-SP%UDX(IC)*VB(I+1,1,K)
                  ENDDO
                 
               ENDDO
               
            !!!--------------- 3D -----------------
            CASE (NSCARC_DIMENSION_THREE)
            
               !!! NOT WORKING YET, has still to be adopted to 3D !!!!!
               DO J=1,NY
               
                  ! backward elimination of first NX unkowns (may be solved by tridiagonal system)
                  DO I=2,NX
                     VB(I,J,1) = VB(I,J,1)-SP%LDX(I-1)*VB(I-1,J,1)
                  ENDDO
                  DO I=1,NX
                     VB(I,J,1) = VB(I,J,1)*SP%MDX(I)
                  ENDDO
                  DO I=NX-1,1,-1
                     VB(I,J,1) = VB(I,J,1)-SP%UDX(I)*VB(I+1,J,1)
                  ENDDO
                    
                  ! backward elimination of following unknowns 
                  ! (here the subdiagonals in y- and z-direction must be taken into account)
                  DO K=2,NZ
                  
                     IC = (K-1)*NX + I
                  
                     ! bring subdiagonal in z-direction to right hand side (corresponding solution component already known)
                     DO I=1,NX
                        VB(I,J,K) = VB(I,J,K) - SP%MDZ(IC)*VB(I,J,K-1) 
                     ENDDO
                  
                     ! perform elimination of matrix lines corresponding to K
                     DO I=2,NX
                        VB(I,J,K) = VB(I,J,K)-SP%LDX(IC-1)*VB(I-1,J,K)
                     ENDDO
                     DO I=1,NX
                        VB(I,J,K) = VB(I,J,K)*SP%MDX(IC)
                     ENDDO
                     DO I=NX-1,1,-1
                        VB(I,J,K) = VB(I,J,K)-SP%UDX(IC)*VB(I+1,J,K)
                     ENDDO
                    
                  ENDDO
               
               ENDDO
               
         END SELECT SELECT_BANDED_DIMENSION

      ENDDO
       

   !!!-------------------------------------------------------------------------------------------------
   !!! Compact system
   !!!-------------------------------------------------------------------------------------------------
   CASE (NSCARC_SYSTEM_COMPACT)

      DO NM = NMESHES_MIN, NMESHES_MAX

         VC => POINT_TO_CVECTOR (NVECTOR, NM, NL)

         SC => SCARC(NM)%COMPACT(NL)
         SP => SCARC(NM)%PRECON(NL)
         
         NX => SC%NX
         NY => SC%NY
         NZ => SC%NZ
         
         SELECT_COMPACT_DIMENSION: SELECT CASE (TYPE_DIMENSION)
         
            !!!--------------- 2D -----------------
            CASE (NSCARC_DIMENSION_TWO)
         
               ! backward elimination of first NX unkowns (may be solved by tridiagonal system)
               DO IC=2,NX
                  VC(IC) = VC(IC)-SP%LDX(I-1)*VC(IC-1)
               ENDDO
               DO I=1,NX
                  VC(IC) = VC(IC)*SP%MDX(I)
               ENDDO
               DO I=NX-1,1,-1
                  VC(IC) = VC(IC)-SP%UDX(I)*VC(IC+1)
               ENDDO
                 
               
               ! backward elimination of following unknowns (here the subdiagonal in z-direction must be taken into account)
               DO K=2,NZ
               
                  IC = (K-1)*NX + I
               
                  ! eliminate subdiagonal in z-direction into right hand side (related solution component already known)
                  DO I=1,NX
                     VC(IC) = VC(IC) - SP%MDZ(IC)*VC(IC-1) 
                  ENDDO
               
                  ! perform elimination of matrix lines corresponding to K
                  DO I=2,NX
                     VC(IC) = VC(IC)-SP%LDX(IC-1)*VC(IC-1)
                  ENDDO
                  DO I=1,NX
                     VC(IC) = VC(IC)*SP%MDX(IC)
                  ENDDO
                  DO I=NX-1,1,-1
                     VC(IC) = VC(IC)-SP%UDX(IC)*VC(IC+1)
                  ENDDO
                 
               ENDDO
               
            !!!--------------- 3D -----------------
            CASE (NSCARC_DIMENSION_THREE)
            
               !!! NOT WORKING YET, has still to be adopted to 3D !!!!!
               DO J=1,NY
               
                  JC = J*NX
            
                  ! backward elimination of first NX unkowns (may be solved by tridiagonal system)
                  DO I=2,NX
                     VC(JC+I) = VC(JC+I)-SP%LDX(I-1)*VC(JC+I-1)
                  ENDDO
                  DO I=1,NX
                     VC(JC+I) = VC(JC+I)*SP%MDX(I)
                  ENDDO
                  DO I=NX-1,1,-1
                     VC(JC+I) = VC(JC+I)-SP%UDX(I)*VC(JC+I+1)
                  ENDDO
                    
                  ! backward elimination of following unknowns 
                  ! (here the subdiagonals in y- and z-direction must be taken into account)
                  DO K=2,NZ
                  
                     IC = (K-1)*NX*NY + (J-1)*NX + I
                  
                     ! bring subdiagonal in z-direction to right hand side (corresponding solution component already known)
                     DO I=1,NX
                        VC(IC) = VC(IC) - SP%MDZ(IC)*VC(IC-1) 
                     ENDDO
                  
                     ! perform elimination of matrix lines corresponding to K
                     DO I=2,NX
                        VC(IC) = VC(IC)-SP%LDX(IC-1)*VC(IC-1)
                     ENDDO
                     DO I=1,NX
                        VC(IC) = VC(IC)*SP%MDX(IC)
                     ENDDO
                     DO I=NX-1,1,-1
                        VC(IC) = VC(IC)-SP%UDX(IC)*VC(IC+1)
                     ENDDO
                    
                  ENDDO
               
               ENDDO
            
         END SELECT SELECT_COMPACT_DIMENSION
          
      ENDDO

END SELECT SELECT_SYSTEM

END SUBROUTINE SCARC_GSTRIX

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Setup environement in every solver call (i.e. set pointers to used vectors)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SETUP_ENVIRONMENT(CROUTINE, NSCOPE, NRHS, NL)
CHARACTER(*), INTENT(OUT) :: CROUTINE
INTEGER, INTENT(IN)  :: NSCOPE, NRHS
INTEGER, INTENT(OUT) :: NL
INTEGER:: NMETHOD


SELECT CASE (NSCOPE)
   CASE (NSCARC_SCOPE_COARSE)
      NMETHOD = NSCARC_METHOD_KRYLOV
   CASE (NSCARC_SCOPE_PRECON)
      NMETHOD = NSCARC_METHOD_MULTIGRID
   CASE DEFAULT
      NMETHOD = TYPE_METHOD
END SELECT

SELECT CASE (NMETHOD)
   
   !!!-------------------------------------------------------------------------------------------------------
   !!! Krylov method
   !!!-------------------------------------------------------------------------------------------------------
   CASE (NSCARC_METHOD_KRYLOV)

      SELECT CASE (TYPE_KRYLOV)
      
         !!! CG-method
         CASE (NSCARC_KRYLOV_CG)
      

            VEC_X = NSCARC_VECTOR_X
            VEC_D = NSCARC_VECTOR_D
            VEC_G = NSCARC_VECTOR_G
            VEC_Y = NSCARC_VECTOR_Y
            VEC_W = NSCARC_VECTOR_W

            VEC_F = NRHS                                             ! set correct right hand side vector
      
            SELECT CASE (NSCOPE)
               CASE (NSCARC_SCOPE_MAIN)
                  CROUTINE = 'SCARC_GLOBAL_CG'
                  NL = NLEVEL_MIN
                  !TYPE_ACCURACY = NSCARC_ACCURACY_RELATIVE
               CASE (NSCARC_SCOPE_COARSE)
                  CROUTINE = 'SCARC_COARSE_CG'
                  NL = NLEVEL_MAX
                  !TYPE_ACCURACY = NSCARC_ACCURACY_ABSOLUTE
            END SELECT
      
         !!! BICG-method
         CASE (NSCARC_KRYLOV_BICG)
      
            VEC_X = NSCARC_VECTOR_X
            VEC_D = NSCARC_VECTOR_D
            VEC_G = NSCARC_VECTOR_G
            VEC_Y = NSCARC_VECTOR_Y
            VEC_W = NSCARC_VECTOR_W
            VEC_Z = NSCARC_VECTOR_Z
      
            VEC_F = NRHS                                             ! set correct right hand side vector

            SELECT CASE (NSCOPE)
               CASE (NSCARC_SCOPE_MAIN)
                  CROUTINE = 'SCARC_GLOBAL_BICG'
                  NL = NLEVEL_MIN
                  !TYPE_ACCURACY = NSCARC_ACCURACY_RELATIVE
               CASE (NSCARC_SCOPE_COARSE)
                  CROUTINE = 'SCARC_COARSE_BICG'
                  NL = NLEVEL_MAX
                  !TYPE_ACCURACY = NSCARC_ACCURACY_ABSOLUTE
            END SELECT
      
      END SELECT
   
   !!!-------------------------------------------------------------------------------------------------------
   !!! Multigrid method
   !!!-------------------------------------------------------------------------------------------------------
   CASE (NSCARC_METHOD_MULTIGRID)

      NL = NLEVEL_MIN

      VEC_F = NRHS                                                   ! set correct right hand side vector

      !!! select scope (multigrid as main solver or preconditioner)
      SELECT CASE (NSCOPE)
         CASE (NSCARC_SCOPE_MAIN)

            CROUTINE = 'SCARC_GLOBAL_MULTIGRID'

            VEC_X = NSCARC_VECTOR_X
            VEC_D = NSCARC_VECTOR_D

            !TYPE_ACCURACY = NSCARC_ACCURACY_RELATIVE

         CASE (NSCARC_SCOPE_PRECON)

            CROUTINE = 'SCARC_PRECON_MULTIGRID'

            VEC_X = NSCARC_VECTOR_X2
            VEC_D = NSCARC_VECTOR_D2

            !TYPE_ACCURACY = NSCARC_ACCURACY_RELATIVE

      END SELECT

END SELECT
   
!WRITE(*,*) 'SETUP_ENVIRONEMENT: NSCOPE=',NSCOPE,', NRHS=',NRHS,', NL=',NL, ', NMETHOD=',NMETHOD
END SUBROUTINE SCARC_SETUP_ENVIRONMENT


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Set initial solution corresponding to boundary data in BXS, BXF, ...
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SETUP_SOLVER(NL)
INTEGER, INTENT(IN) :: NL
INTEGER :: NM, IW, IOR0, I, J, K, IC
TYPE (MESH_TYPE), POINTER ::  M
TYPE (SCARC_BANDED_TYPE) , POINTER :: SB
TYPE (SCARC_COMPACT_TYPE), POINTER :: SC


!WRITE(*,*) 'SETUP_SOLVER: NL=',NL, ', TYPE_SCOPE=',TYPE_SCOPE,': TYPE_PRECON=',TYPE_PRECON

SELECT_SYSTEM: SELECT CASE (TYPE_SYSTEM)

   !!!-------------------------------------------------------------------------------------------------------
   !!! Banded system
   !!!-------------------------------------------------------------------------------------------------------
   CASE (NSCARC_SYSTEM_BANDED)

      SELECT_BANDED_METHOD: SELECT CASE (TYPE_METHOD)

         !!! In case of a Krylov method initialize auxiliary arrays
         CASE (NSCARC_METHOD_KRYLOV)

            DO NM = NMESHES_MIN, NMESHES_MAX
               SCARC(NM)%BANDED(NL)%D = 0.0_EB
               SCARC(NM)%BANDED(NL)%G = 0.0_EB
               SCARC(NM)%BANDED(NL)%Y = 0.0_EB
               SCARC(NM)%BANDED(NL)%W = 0.0_EB
            ENDDO

            IF (TYPE_KRYLOV == NSCARC_KRYLOV_BICG) THEN
               DO NM = NMESHES_MIN, NMESHES_MAX
                  SCARC(NM)%BANDED(NL)%Z = 0.0_EB
               ENDDO
            ENDIF
            
            IF (TYPE_PRECON == NSCARC_PRECON_MULTIGRID) THEN
               
            ENDIF

         !!! In case of a multigrid method with coarse grid solution by CG, clear CG-vectors on max level
         CASE (NSCARC_METHOD_MULTIGRID)
            
            IF (TYPE_COARSE == NSCARC_COARSE_ITERATIVE) THEN
               DO NM = NMESHES_MIN, NMESHES_MAX
                  SCARC(NM)%BANDED(NLEVEL_MAX)%X = 0.0_EB
                  SCARC(NM)%BANDED(NLEVEL_MAX)%D = 0.0_EB
                  SCARC(NM)%BANDED(NLEVEL_MAX)%G = 0.0_EB
                  SCARC(NM)%BANDED(NLEVEL_MAX)%Y = 0.0_EB
                  SCARC(NM)%BANDED(NLEVEL_MAX)%W = 0.0_EB
               ENDDO
            ENDIF

      END SELECT SELECT_BANDED_METHOD

      IF (TYPE_SCOPE == NSCARC_SCOPE_MAIN) THEN
               
!WRITE(*,*) 'IN WALLCELL_LOOP'
         !!! Initialize solution and right hand side vector corresponding to boundary conditions
         SELECT_BANDED_DIMENSION: SELECT CASE (TYPE_DIMENSION)
         
            !!! -------------- 2D ------------------
            CASE (NSCARC_DIMENSION_TWO)

               DO NM = NMESHES_MIN, NMESHES_MAX
         
                  M  => MESHES(NM)
                  SB => SCARC(NM)%BANDED(NL)
               
                  ! get right hand side vector
                  SB%F (1:M%IBAR, 1, 1:M%KBAR) = M%PRHS (1:M%IBAR, 1, 1:M%KBAR)
            
                  ! set ghost values of solution vector to zero and use values of H/HS in interior
                  SB%X(1:M%IBAR, 1, 0       ) = 0.0_EB
                  SB%X(1:M%IBAR, 1, M%KBAR+1) = 0.0_EB
                  SB%X(0       , 1, 1:M%KBAR) = 0.0_EB
                  SB%X(M%IBAR+1, 1, 1:M%KBAR) = 0.0_EB
                  IF (PREDICTOR) THEN
                     SB%X (1:M%IBAR, 1, 1:M%KBAR) = M%H (1:M%IBAR, 1, 1:M%KBAR)
                  ELSE
                     SB%X (1:M%IBAR, 1, 1:M%KBAR) = M%HS(1:M%IBAR, 1, 1:M%KBAR)
                  ENDIF
            
                  ! set boundary conditions at exterior boundaries (corresponding to pois.f90)
                  BANDED_WALLCELL_LOOP2D: DO IW = 1, M%N_EXTERNAL_WALL_CELLS
               
                     I = M%OMESH(NM)%IJKW(6,IW)!****CHECK
                     J = M%OMESH(NM)%IJKW(7,IW)!****CHECK
                     K = M%OMESH(NM)%IJKW(8,IW)!****CHECK
            
                     IF (J /= 1) THEN
                        WRITE(*,*) 'Wrong index for J =',J,' in SCARC_SETUP_SOLVER !!!'
                        STOP
                     ENDIF
                  
                     IOR0 = M%OMESH(NM)%IJKW(4,IW)!****CHECK
                  
                     SELECT CASE (IOR0)
                        CASE (1)
                           IF (M%WALL(IW)%PRESSURE_BC_INDEX==DIRICHLET) THEN
                              SB%F(I,J,K) = SB%F(I,J,K) - 2.0_EB * SB%DXI2 * M%BXS(1,K)         ! Dirichlet
                           ELSE IF (M%WALL(IW)%PRESSURE_BC_INDEX==NEUMANN) THEN
                              SB%F(I,J,K) = SB%F(I,J,K) + SB%DXI * M%BXS(1,K)                   ! Neumann
                           ENDIF
                        CASE (-1)
                           IF (M%WALL(IW)%PRESSURE_BC_INDEX==DIRICHLET) THEN
                              SB%F(I,J,K) = SB%F(I,J,K) - 2.0_EB * SB%DXI2 *M%BXF(1,K)          ! Dirichlet
                           ELSE IF (M%WALL(IW)%PRESSURE_BC_INDEX==NEUMANN) THEN
                              SB%F(I,J,K) = SB%F(I,J,K) - SB%DXI *M%BXF(1,K)                    ! Neumann
                           ENDIF
                        CASE (3)
                           IF (M%WALL(IW)%PRESSURE_BC_INDEX==DIRICHLET) THEN
                              SB%F(I,J,K) = SB%F(I,J,K) - 2.0_EB * SB%DZI2 * M%BZS(I,1)         ! Dirichlet
                           ELSE IF (M%WALL(IW)%PRESSURE_BC_INDEX==NEUMANN) THEN
                              SB%F(I,J,K) = SB%F(I,J,K) + SB%DZI * M%BZS(I,1)                   ! Neumann
                           ENDIF
                        CASE (-3)
                           IF (M%WALL(IW)%PRESSURE_BC_INDEX==DIRICHLET) THEN
                              SB%F(I,J,K) = SB%F(I,J,K) - 2.0_EB * SB%DZI2 * M%BZF(I,1)         ! Dirichlet
                           ELSE IF (M%WALL(IW)%PRESSURE_BC_INDEX==NEUMANN) THEN
                              SB%F(I,J,K) = SB%F(I,J,K) - SB%DZI  * M%BZF(I,1)                  ! Neumann
                           ENDIF
                     END SELECT
                  
                  ENDDO BANDED_WALLCELL_LOOP2D
               
               ENDDO
            
            !!! -------------- 3D ------------------
            CASE (NSCARC_DIMENSION_THREE)
            
               DO NM = NMESHES_MIN, NMESHES_MAX
               
                  M  => MESHES(NM)
                  SB => SCARC(NM)%BANDED(NL)
            
                  ! get right hand side vector
                  SB%F (1:M%IBAR, 1:M%JBAR, 1:M%KBAR) = M%PRHS (1:M%IBAR, 1:M%JBAR, 1:M%KBAR)
            
                  ! set ghost values of solution vector to zero and use values of H/HS in interior
                  SB%X(1:M%IBAR, 1:M%JBAR, 0       ) = 0.0_EB
                  SB%X(1:M%IBAR, 1:M%JBAR, M%KBAR+1) = 0.0_EB
                  SB%X(1:M%IBAR, 0       , 1:M%KBAR) = 0.0_EB
                  SB%X(1:M%IBAR, M%JBAR+1, 1:M%KBAR) = 0.0_EB
                  SB%X(0       , 1:M%JBAR, 1:M%KBAR) = 0.0_EB
                  SB%X(M%IBAR+1, 1:M%JBAR, 1:M%KBAR) = 0.0_EB
                  IF (PREDICTOR) THEN
                     SB%X (1:M%IBAR, 1:M%JBAR, 1:M%KBAR) = M%H (1:M%IBAR, 1:M%JBAR, 1:M%KBAR)
                  ELSE
                     SB%X (1:M%IBAR, 1:M%JBAR, 1:M%KBAR) = M%HS(1:M%IBAR, 1:M%JBAR, 1:M%KBAR)
                  ENDIF
            
                  BANDED_WALLCELL_LOOP3D: DO IW = 1, M%N_EXTERNAL_WALL_CELLS
                  
                     I = M%OMESH(NM)%IJKW(6,IW)!****CHECK
                     J = M%OMESH(NM)%IJKW(7,IW)!****CHECK
                     K = M%OMESH(NM)%IJKW(8,IW)!****CHECK
                  
                     IOR0 = M%OMESH(NM)%IJKW(4,IW)!****CHECK
                  
                     SELECT CASE (IOR0)
                        CASE (1)
                           IF (M%WALL(IW)%PRESSURE_BC_INDEX==DIRICHLET) THEN
                              SB%F(I,J,K) = SB%F(I,J,K) - 2.0_EB * SB%DXI2 * M%BXS(J,K)         ! Dirichlet
                           ELSE IF (M%WALL(IW)%PRESSURE_BC_INDEX==NEUMANN) THEN
                              SB%F(I,J,K) = SB%F(I,J,K) + SB%DXI * M%BXS(J,K)                   ! Neumann
                           ENDIF
                        CASE (-1)
                           IF (M%WALL(IW)%PRESSURE_BC_INDEX==DIRICHLET) THEN
                              SB%F(I,J,K) = SB%F(I,J,K) - 2.0_EB * SB%DXI2 *M%BXF(J,K)          ! Dirichlet
                           ELSE IF (M%WALL(IW)%PRESSURE_BC_INDEX==NEUMANN) THEN
                              SB%F(I,J,K) = SB%F(I,J,K) - SB%DXI *M%BXF(J,K)                    ! Neumann
                           ENDIF
                        CASE (2)
                           IF (M%WALL(IW)%PRESSURE_BC_INDEX==DIRICHLET) THEN
                              SB%F(I,J,K) = SB%F(I,J,K) - 2.0_EB * SB%DYI2 * M%BYS(I,K)         ! Dirichlet
                           ELSE IF (M%WALL(IW)%PRESSURE_BC_INDEX==NEUMANN) THEN
                              SB%F(I,J,K) = SB%F(I,J,K) + SB%DYI * M%BYS(I,K)                   ! Neumann
                           ENDIF
                        CASE (-2)
                           IF (M%WALL(IW)%PRESSURE_BC_INDEX==DIRICHLET) THEN
                              SB%F(I,J,K) = SB%F(I,J,K) - 2.0_EB * SB%DYI2 *M%BYF(I,K)          ! Dirichlet
                           ELSE IF (M%WALL(IW)%PRESSURE_BC_INDEX==NEUMANN) THEN
                              SB%F(I,J,K) = SB%F(I,J,K) - SB%DYI *M%BYF(I,K)                    ! Neumann
                           ENDIF
                        CASE (3)
                           IF (M%WALL(IW)%PRESSURE_BC_INDEX==DIRICHLET) THEN
                              SB%F(I,J,K) = SB%F(I,J,K) - 2.0_EB * SB%DZI2 * M%BZS(I,J)         ! Dirichlet
                           ELSE IF (M%WALL(IW)%PRESSURE_BC_INDEX==NEUMANN) THEN
                              SB%F(I,J,K) = SB%F(I,J,K) + SB%DZI * M%BZS(I,J)                   ! Neumann
                           ENDIF
                        CASE (-3)
                           IF (M%WALL(IW)%PRESSURE_BC_INDEX==DIRICHLET) THEN
                              SB%F(I,J,K) = SB%F(I,J,K) - 2.0_EB * SB%DZI2 * M%BZF(I,J)         ! Dirichlet
                           ELSE IF (M%WALL(IW)%PRESSURE_BC_INDEX==NEUMANN) THEN
                              SB%F(I,J,K) = SB%F(I,J,K) - SB%DZI  * M%BZF(I,J)                  ! Neumann
                           ENDIF
                     END SELECT
                  
                  ENDDO BANDED_WALLCELL_LOOP3D
            
               ENDDO
         
         END SELECT SELECT_BANDED_DIMENSION

   

      ENDIF
         
   !!!-------------------------------------------------------------------------------------------------------
   !!! Compact system
   !!!-------------------------------------------------------------------------------------------------------
   CASE (NSCARC_SYSTEM_COMPACT)

      SELECT_COMPACT_METHOD: SELECT CASE (TYPE_METHOD)

         !!! In case of a Krylov method initialize auxiliary arrays
         CASE (NSCARC_METHOD_KRYLOV)

            DO NM = NMESHES_MIN, NMESHES_MAX
               SCARC(NM)%COMPACT(NL)%D = 0.0_EB
               SCARC(NM)%COMPACT(NL)%G = 0.0_EB
               SCARC(NM)%COMPACT(NL)%Y = 0.0_EB
               SCARC(NM)%COMPACT(NL)%W = 0.0_EB
            ENDDO

            IF (TYPE_KRYLOV == NSCARC_KRYLOV_BICG) THEN
               DO NM = NMESHES_MIN, NMESHES_MAX
                  SCARC(NM)%COMPACT(NL)%Z = 0.0_EB
               ENDDO
            ENDIF
            
         !!! In case of a multigrid method with coarse grid solution by CG, clear CG-vectors on max level
         CASE (NSCARC_METHOD_MULTIGRID)
            
            DO NM = NMESHES_MIN, NMESHES_MAX
               SCARC(NM)%COMPACT(NLEVEL_MAX)%X = 0.0_EB
               SCARC(NM)%COMPACT(NLEVEL_MAX)%D = 0.0_EB
               SCARC(NM)%COMPACT(NLEVEL_MAX)%G = 0.0_EB
               SCARC(NM)%COMPACT(NLEVEL_MAX)%Y = 0.0_EB
               SCARC(NM)%COMPACT(NLEVEL_MAX)%W = 0.0_EB
            ENDDO

      END SELECT SELECT_COMPACT_METHOD
      
      IF (TYPE_SCOPE == NSCARC_SCOPE_MAIN) THEN
         
         !!! Initialize solution and right hand side vector corresponding to boundary conditions
         SELECT_COMPACT_DIMENSION: SELECT CASE (TYPE_DIMENSION)
         
            !!! -------------- 2D ------------------
            CASE (NSCARC_DIMENSION_TWO)
         
               DO NM = NMESHES_MIN, NMESHES_MAX
         
                  M  => MESHES(NM)
                  SC => SCARC(NM)%COMPACT(NL)
               
                  ! get right hand side and initial vector
                  DO K = 1, M%KBAR
                     DO I = 1, M%IBAR
                        IC = (K-1)*M%IBAR + I
                        SC%F(IC) = M%PRHS (I, 1, K)
                        IF (PREDICTOR) THEN
                           SC%X(IC) = M%H (I, 1, K)
                        ELSE
                           SC%X(IC) = M%HS(I, 1, K)
                        ENDIF
                     ENDDO
                  ENDDO
            
            
                  ! set boundary conditions at exterior boundaries (corresponding to pois.f90)
                  COMPACT_WALLCELL_LOOP2D: DO IW = 1, M%N_EXTERNAL_WALL_CELLS
               
                     I = M%OMESH(NM)%IJKW(6,IW)!****CHECK
                     J = M%OMESH(NM)%IJKW(7,IW)!****CHECK
                     K = M%OMESH(NM)%IJKW(8,IW)!****CHECK
            
                     IF (J /= 1) THEN
                        WRITE(*,*) 'Wrong index for J =',J,' in SCARC_SETUP_SOLVER !!!'
                        STOP
                     ENDIF
                  
                     IOR0 = M%OMESH(NM)%IJKW(4,IW)!****CHECK
                     IC = (K-1)*M%IBAR + I
                     
                     SELECT CASE (IOR0)
                        CASE (1)
                           IF (M%WALL(IW)%PRESSURE_BC_INDEX==DIRICHLET) THEN
                              SC%F(IC) = SC%F(IC) - 2.0_EB * SC%DXI2 * M%BXS(1,K)         ! Dirichlet
                           ELSE IF (M%WALL(IW)%PRESSURE_BC_INDEX==NEUMANN) THEN
                              SC%F(IC) = SC%F(IC) + SC%DXI * M%BXS(1,K)                   ! Neumann
                           ENDIF
                        CASE (-1)
                           IF (M%WALL(IW)%PRESSURE_BC_INDEX==DIRICHLET) THEN
                              SC%F(IC) = SC%F(IC) - 2.0_EB * SC%DXI2 *M%BXF(1,K)          ! Dirichlet
                           ELSE IF (M%WALL(IW)%PRESSURE_BC_INDEX==NEUMANN) THEN
                              SC%F(IC) = SC%F(IC) - SC%DXI *M%BXF(1,K)                    ! Neumann
                           ENDIF
                        CASE (3)
                           IF (M%WALL(IW)%PRESSURE_BC_INDEX==DIRICHLET) THEN
                              SC%F(IC) = SC%F(IC) - 2.0_EB * SC%DZI2 * M%BZS(I,1)         ! Dirichlet
                           ELSE IF (M%WALL(IW)%PRESSURE_BC_INDEX==NEUMANN) THEN
                              SC%F(IC) = SC%F(IC) + SC%DZI * M%BZS(I,1)                   ! Neumann
                           ENDIF
                        CASE (-3)
                           IF (M%WALL(IW)%PRESSURE_BC_INDEX==DIRICHLET) THEN
                              SC%F(IC) = SC%F(IC) - 2.0_EB * SC%DZI2 * M%BZF(I,1)         ! Dirichlet
                           ELSE IF (M%WALL(IW)%PRESSURE_BC_INDEX==NEUMANN) THEN
                              SC%F(IC) = SC%F(IC) - SC%DZI  * M%BZF(I,1)                  ! Neumann
                           ENDIF
                     END SELECT
                  
                  ENDDO COMPACT_WALLCELL_LOOP2D
         
               ENDDO
            !!! -------------- 3D ------------------
            CASE (NSCARC_DIMENSION_THREE)
            
               DO NM = NMESHES_MIN, NMESHES_MAX
         
                  M  => MESHES(NM)
                  SC => SCARC(NM)%COMPACT(NL)
            
                  ! get right hand side and initial vector
                  DO K = 1, M%KBAR
                     DO J = 1, M%JBAR
                        DO I = 1, M%IBAR
                           IC = (K-1) * M%IBAR * M%JBAR + (J-1) * M%IBAR + I
                           SC%F(IC) = M%PRHS (I, J, K)
                           IF (PREDICTOR) THEN
                              SC%X(IC) = M%H (I, J, K)
                           ELSE
                              SC%X(IC) = M%HS(I, J, K)
                           ENDIF
                        ENDDO
                     ENDDO
                  ENDDO
            
                  COMPACT_WALLCELL_LOOP3D: DO IW = 1, M%N_EXTERNAL_WALL_CELLS
                  
                     I = M%OMESH(NM)%IJKW(6,IW)!****CHECK
                     J = M%OMESH(NM)%IJKW(7,IW)!****CHECK
                     K = M%OMESH(NM)%IJKW(8,IW)!****CHECK
                  
                     IOR0 = M%OMESH(NM)%IJKW(4,IW)!****CHECK
                     IC = (K-1) * M%IBAR * M%JBAR + (J-1) * M%IBAR + I
                  
                     SELECT CASE (IOR0)
                        CASE (1)
                           IF (M%WALL(IW)%PRESSURE_BC_INDEX==DIRICHLET) THEN
                              SC%F(IC) = SC%F(IC) - 2.0_EB * SC%DXI2 * M%BXS(J,K)         ! Dirichlet
                           ELSE IF (M%WALL(IW)%PRESSURE_BC_INDEX==NEUMANN) THEN
                              SC%F(IC) = SC%F(IC) + SC%DXI * M%BXS(J,K)                   ! Neumann
                           ENDIF
                        CASE (-1)
                           IF (M%WALL(IW)%PRESSURE_BC_INDEX==DIRICHLET) THEN
                              SC%F(IC) = SC%F(IC) - 2.0_EB * SC%DXI2 *M%BXF(J,K)          ! Dirichlet
                           ELSE IF (M%WALL(IW)%PRESSURE_BC_INDEX==NEUMANN) THEN
                              SC%F(IC) = SC%F(IC) - SC%DXI *M%BXF(J,K)                    ! Neumann
                           ENDIF
                        CASE (2)
                           IF (M%WALL(IW)%PRESSURE_BC_INDEX==DIRICHLET) THEN
                              SC%F(IC) = SC%F(IC) - 2.0_EB * SC%DYI2 * M%BYS(I,K)         ! Dirichlet
                           ELSE IF (M%WALL(IW)%PRESSURE_BC_INDEX==NEUMANN) THEN
                              SC%F(IC) = SC%F(IC) + SC%DYI * M%BYS(I,K)                   ! Neumann
                           ENDIF
                        CASE (-2)
                           IF (M%WALL(IW)%PRESSURE_BC_INDEX==DIRICHLET) THEN
                              SC%F(IC) = SC%F(IC) - 2.0_EB * SC%DYI2 *M%BYF(I,K)          ! Dirichlet
                           ELSE IF (M%WALL(IW)%PRESSURE_BC_INDEX==NEUMANN) THEN
                              SC%F(IC) = SC%F(IC) - SC%DYI *M%BYF(I,K)                    ! Neumann
                           ENDIF
                        CASE (3)
                           IF (M%WALL(IW)%PRESSURE_BC_INDEX==DIRICHLET) THEN
                              SC%F(IC) = SC%F(IC) - 2.0_EB * SC%DZI2 * M%BZS(I,J)         ! Dirichlet
                           ELSE IF (M%WALL(IW)%PRESSURE_BC_INDEX==NEUMANN) THEN
                              SC%F(IC) = SC%F(IC) + SC%DZI * M%BZS(I,J)                   ! Neumann
                           ENDIF
                        CASE (-3)
                           IF (M%WALL(IW)%PRESSURE_BC_INDEX==DIRICHLET) THEN
                              SC%F(IC) = SC%F(IC) - 2.0_EB * SC%DZI2 * M%BZF(I,J)         ! Dirichlet
                           ELSE IF (M%WALL(IW)%PRESSURE_BC_INDEX==NEUMANN) THEN
                              SC%F(IC) = SC%F(IC) - SC%DZI  * M%BZF(I,J)                  ! Neumann
                           ENDIF
                     END SELECT
                  
                  ENDDO COMPACT_WALLCELL_LOOP3D
         
               ENDDO
         
         END SELECT SELECT_COMPACT_DIMENSION

      ENDIF

END SELECT SELECT_SYSTEM

END SUBROUTINE SCARC_SETUP_SOLVER


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Print out residual information for loop ITE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_CONVERGENCE_INFO(RES, ITE, NL, CROUTINE)
INTEGER, INTENT(IN) :: ITE, NL
REAL(EB), INTENT(IN) :: RES
CHARACTER(*), INTENT(IN) :: CROUTINE
INTEGER:: NM

DO NM = NMESHES_MIN, NMESHES_MAX
   IF (TYPE_DEBUG > NSCARC_DEBUG_NONE.AND.NM==1) WRITE(SCARC_LU,1000) TRIM(CROUTINE), NM, NL, ITE,  RES
   IF (TYPE_DEBUG > NSCARC_DEBUG_NONE.AND.NM==1) write(*       ,1000) TRIM(CROUTINE), NM, NL, ITE,  RES
ENDDO

1000 FORMAT (5X,A30,': mesh', i4,': level=',i4,': #ite= ',i4,': res =',e14.6)
END SUBROUTINE SCARC_CONVERGENCE_INFO


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Check if solver converges or diverges and print out residual information
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
INTEGER FUNCTION SCARC_CONVERGENCE_STATE(RESIN, RES, EPS, ITE, NL, CROUTINE)

INTEGER, INTENT(IN) :: ITE, NL
INTEGER :: NM, ISTATE
REAL(EB), INTENT(IN) :: RESIN, RES, EPS
CHARACTER(*), INTENT(IN) :: CROUTINE

ISTATE = NSCARC_STATE_PROCEED
SCARC_RESIDUAL = RES

DO NM = NMESHES_MIN, NMESHES_MAX
   IF (TYPE_DEBUG > NSCARC_DEBUG_NONE.AND.NM==1) WRITE(SCARC_LU,1000) TRIM(CROUTINE), NM, NL, ITE, RES
   IF (TYPE_DEBUG > NSCARC_DEBUG_NONE.AND.NM==1) WRITE(*       ,1000) TRIM(CROUTINE), NM, NL, ITE, RES
ENDDO

SELECT CASE (TYPE_ACCURACY)
   CASE (NSCARC_ACCURACY_RELATIVE)
      IF (RES <= RESIN*EPS)  ISTATE = NSCARC_STATE_CONV
      IF (RES <= 1.0E-15)    ISTATE = NSCARC_STATE_CONV
   CASE (NSCARC_ACCURACY_ABSOLUTE)
      IF (RES <= EPS .AND. RES <= RESIN*SCARC_ACCURACY_RELATIVE) ISTATE = NSCARC_STATE_CONV
END SELECT
IF (RES > SCARC_ACCURACY_DIVERGENCE) ISTATE = NSCARC_STATE_DIVG

SCARC_CONVERGENCE_STATE = ISTATE
RETURN

1000 FORMAT (5X,A30,': mesh', i4,': level=',i4,': #ite= ',i4,': res =',e14.6)
END FUNCTION SCARC_CONVERGENCE_STATE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Compute convergence rate and print out resiual information for final loop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_CONVERGENCE_RATE(RESIN, RES, ITE, ISTATE, CROUTINE)

INTEGER, INTENT(IN) :: ITE, ISTATE
INTEGER :: NM
REAL(EB), INTENT(IN) :: RESIN, RES
CHARACTER(*), INTENT(IN) :: CROUTINE

SCARC_RESIDUAL = RES
IF (ISTATE == NSCARC_STATE_DIVG) THEN                       
   SCARC_ITERATIONS = - 1
   SCARC_CAPPA      = 1.0_EB
ELSE 
   IF (ISTATE == NSCARC_STATE_CONV) THEN
     SCARC_ITERATIONS = ITE
   ELSE
     SCARC_ITERATIONS = ITE-1
   ENDIF
   IF (RESIN >= 1.0E-70_EB) THEN
      SCARC_CAPPA = (RES/RESIN) ** (1.0_EB/SCARC_ITERATIONS)
   ELSE 
      SCARC_CAPPA = 0.0_EB
   ENDIF
ENDIF

DO NM = NMESHES_MIN, NMESHES_MAX
   IF (TYPE_DEBUG > NSCARC_DEBUG_NONE.AND.NM==1) WRITE(SCARC_LU,2000) TRIM(CROUTINE), SCARC_CAPPA 
   IF (TYPE_DEBUG > NSCARC_DEBUG_NONE.AND.NM==1) WRITE(*       ,2000) TRIM(CROUTINE), SCARC_CAPPA 
ENDDO

2000 FORMAT (5X,A30,':',10X,'---> convergence rate =',e14.6,/, &
             5X ,'----------------------------------------------------------------------------------------')
END SUBROUTINE SCARC_CONVERGENCE_RATE

 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Perform restriction from finer to coarser grid in multigrid method
!!!    - 'FI' corresponds to finer   grid
!!!    - 'CO' corresponds to coarser grid
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_RESTRICTION (NVECTOR_FI, NVECTOR_CO, NL_FI, NL_CO)
INTEGER, INTENT(IN) :: NVECTOR_FI, NVECTOR_CO, NL_FI, NL_CO
INTEGER :: NM, ICOL, IC
REAL(EB), POINTER, DIMENSION(:,:,:) :: FB_CO, DB_FI
REAL(EB), POINTER, DIMENSION(:)     :: FC_CO, DC_FI, R
INTEGER , POINTER, DIMENSION(:)     :: R_ROW, R_COL
INTEGER , POINTER :: NX_CO, NY_CO, NZ_CO, NC_CO
INTEGER  :: NX_FI, NY_FI, NZ_FI
INTEGER  :: IX_FI, IY_FI, IZ_FI, IC_FI(8)
INTEGER  :: IX_CO, IY_CO, IZ_CO, IC_CO
REAL(EB) :: AUX

TYPE_VECTOR = NVECTOR_FI
IF (TYPE_MULTIGRID == NSCARC_MULTIGRID_ALGEBRAIC) CALL SCARC_EXCHANGE(NSCARC_EXCHANGE_VECTOR, NL_FI)
CALL SCARC_DEBUG_LEVEL (NVECTOR_FI, 'RESTRICTION','D before restrict', NL_FI)

SELECT_SYSTEM: SELECT CASE (TYPE_SYSTEM)

   !!!----------------------------------------------------------------------------------------------------
   !!! Banded system
   !!!----------------------------------------------------------------------------------------------------
   CASE (NSCARC_SYSTEM_BANDED)

      SELECT_BANDED_MULTIGRID: SELECT CASE (TYPE_MULTIGRID)

         !!!------------------------- Geometric multigrid -----------------------------------------------
         CASE (NSCARC_MULTIGRID_GEOMETRIC)

            DO NM = NMESHES_MIN, NMESHES_MAX
            
               NX_CO => SCARC(NM)%BANDED(NL_CO)%NX
               NY_CO => SCARC(NM)%BANDED(NL_CO)%NY
               NZ_CO => SCARC(NM)%BANDED(NL_CO)%NZ
            
               DB_FI => POINT_TO_BVECTOR(NVECTOR_FI, NM, NL_FI)
               FB_CO => POINT_TO_BVECTOR(NVECTOR_CO, NM, NL_CO)
      
               SELECT_BANDED_DIMENSION: SELECT CASE (TYPE_DIMENSION)

                  !!!----------------- 2D -------------------
                  CASE (NSCARC_DIMENSION_TWO)
               
                     DO IZ_CO = 1, NZ_CO
                        DO IX_CO = 1, NX_CO
                 
                           IX_FI = 2*IX_CO
                           IY_FI = 1
                           IZ_FI = 2*IZ_CO
                
                           FB_CO(IX_CO, 1, IZ_CO) = 0.25_EB * (  DB_FI(IX_FI  , IY_FI, IZ_FI-1)  &
                                                               + DB_FI(IX_FI-1, IY_FI, IZ_FI-1)  &
                                                               + DB_FI(IX_FI  , IY_FI, IZ_FI  )  &
                                                               + DB_FI(IX_FI-1, IY_FI, IZ_FI  ) )
                        ENDDO
                     ENDDO
               
                  !!!----------------- 3D -------------------
                  CASE (NSCARC_DIMENSION_THREE)
                     
                     DO IZ_CO = 1, NZ_CO
                        DO IY_CO = 1, NY_CO
                           DO IX_CO = 1, NX_CO
               
                              IX_FI = 2*IX_CO
                              IY_FI = 2*IY_CO
                              IZ_FI = 2*IZ_CO
               
                              FB_CO(IX_CO, IY_CO , IZ_CO) = 0.125_EB * (  DB_FI(IX_FI-1, IY_FI-1, IZ_FI-1)  &
                                                                        + DB_FI(IX_FI  , IY_FI-1, IZ_FI-1)  &
                                                                        + DB_FI(IX_FI-1, IY_FI  , IZ_FI-1)  &
                                                                        + DB_FI(IX_FI  , IY_FI  , IZ_FI-1)  &
                                                                        + DB_FI(IX_FI-1, IY_FI-1, IZ_FI  )  &
                                                                        + DB_FI(IX_FI  , IY_FI-1, IZ_FI  )  &
                                                                        + DB_FI(IX_FI-1, IY_FI  , IZ_FI  )  &
                                                                        + DB_FI(IX_FI  , IY_FI  , IZ_FI  ) )
                           ENDDO
                        ENDDO
                     ENDDO
            
               END SELECT SELECT_BANDED_DIMENSION

            ENDDO
                  
      
         !!!------------------------- Geometric multigrid -----------------------------------------------
         CASE (NSCARC_MULTIGRID_ALGEBRAIC)

            WRITE(*,*) 'No banded restriction for algebraic multigrid available, stopping program!'
            STOP

      END SELECT SELECT_BANDED_MULTIGRID

      
   !!!----------------------------------------------------------------------------------------------------
   !!! Compact system
   !!!----------------------------------------------------------------------------------------------------
   CASE (NSCARC_SYSTEM_COMPACT)

      SELECT_COMPACT_MULTIGRID: SELECT CASE (TYPE_MULTIGRID)

         !!!------------------------- Geometric multigrid -----------------------------------------------
         CASE (NSCARC_MULTIGRID_GEOMETRIC)

            DO NM = NMESHES_MIN, NMESHES_MAX
         
               NX_CO => SCARC(NM)%COMPACT(NL_CO)%NX
               NY_CO => SCARC(NM)%COMPACT(NL_CO)%NY
               NZ_CO => SCARC(NM)%COMPACT(NL_CO)%NZ
         
               NX_FI = 2*NX_CO
               NY_FI = 2*NY_CO
               NZ_FI = 2*NZ_CO
   
               DC_FI  => POINT_TO_CVECTOR(NVECTOR_FI, NM, NL_FI)
               FC_CO  => POINT_TO_CVECTOR(NVECTOR_CO, NM, NL_CO)

               SELECT_COMPACT_DIMENSION: SELECT CASE (TYPE_DIMENSION)
               
                  !!!----------------- 2D -------------------
                  CASE (NSCARC_DIMENSION_TWO)
               
                     DO IZ_CO = 1, NZ_CO
                        DO IX_CO = 1, NX_CO
                  
                           IX_FI = 2*IX_CO
                           IZ_FI = 2*IZ_CO
                 
                           IC_CO     = (IZ_CO-1)*NX_CO + IX_CO
                
                           IC_FI(1) = (IZ_FI-2)*NX_FI + IX_FI - 1
                           IC_FI(2) = (IZ_FI-2)*NX_FI + IX_FI   
                           IC_FI(3) = (IZ_FI-1)*NX_FI + IX_FI - 1
                           IC_FI(4) = (IZ_FI-1)*NX_FI + IX_FI    
               
                           FC_CO(IC_CO) = 0.25_EB * (  DC_FI(IC_FI(1)) &
                                                     + DC_FI(IC_FI(2)) &
                                                     + DC_FI(IC_FI(3)) &
                                                     + DC_FI(IC_FI(4)) )
               
                        ENDDO
                     ENDDO
                  
                  !!!----------------- 3D -------------------
                  CASE (NSCARC_DIMENSION_THREE)
                     
                     DO IZ_CO = 1, NZ_CO
                        DO IY_CO = 1, NY_CO
                           DO IX_CO = 1, NX_CO
                  
                              IX_FI = 2*IX_CO
                              IY_FI = 2*IY_CO
                              IZ_FI = 2*IZ_CO
                 
                              IC_CO    = (IZ_CO-1)*NX_CO*NY_CO + (IY_CO-1)*NX_CO + IX_CO
                
                              IC_FI(1) = (IZ_FI-2)*NX_FI*NY_FI + (IY_FI-2)*NX_FI + IX_FI - 1
                              IC_FI(2) = (IZ_FI-2)*NX_FI*NY_FI + (IY_FI-2)*NX_FI + IX_FI    
                              IC_FI(3) = (IZ_FI-2)*NX_FI*NY_FI + (IY_FI-1)*NX_FI + IX_FI - 1
                              IC_FI(4) = (IZ_FI-2)*NX_FI*NY_FI + (IY_FI-1)*NX_FI + IX_FI    
                              IC_FI(5) = (IZ_FI-1)*NX_FI*NY_FI + (IY_FI-2)*NX_FI + IX_FI - 1
                              IC_FI(6) = (IZ_FI-1)*NX_FI*NY_FI + (IY_FI-2)*NX_FI + IX_FI    
                              IC_FI(7) = (IZ_FI-1)*NX_FI*NY_FI + (IY_FI-1)*NX_FI + IX_FI - 1
                              IC_FI(8) = (IZ_FI-1)*NX_FI*NY_FI + (IY_FI-1)*NX_FI + IX_FI    
               
!WRITE(SCARC_LU,'(12i4)') IX_FI, IY_FI, IZ_FI, IC_CO, IC_FI(1:8)
                              FC_CO(IC_CO) = 0.125_EB * (  DC_FI(IC_FI(1)) &
                                                         + DC_FI(IC_FI(2)) &
                                                         + DC_FI(IC_FI(3)) &
                                                         + DC_FI(IC_FI(4)) &
                                                         + DC_FI(IC_FI(5)) &
                                                         + DC_FI(IC_FI(6)) &
                                                         + DC_FI(IC_FI(7)) &
                                                         + DC_FI(IC_FI(8)) )
                           ENDDO
                        ENDDO
                     ENDDO
               
               END SELECT SELECT_COMPACT_DIMENSION

            ENDDO

         !!!------------------------- Algebraic multigrid -----------------------------------------------
         CASE (NSCARC_MULTIGRID_ALGEBRAIC)

            DO NM = NMESHES_MIN, NMESHES_MAX

               DC_FI => POINT_TO_CVECTOR(NVECTOR_FI, NM, NL_FI)
               FC_CO => POINT_TO_CVECTOR(NVECTOR_CO, NM, NL_CO)

               NC_CO => SCARC(NM)%COMPACT(NL_CO)%NC

               R     => SCARC(NM)%COMPACT(NL_FI)%R
               R_ROW => SCARC(NM)%COMPACT(NL_FI)%R_ROW
               R_COL => SCARC(NM)%COMPACT(NL_FI)%R_COL
   
               DO IC_CO = 1, NC_CO
                  AUX = 0.0_EB
                  DO ICOL = R_ROW(IC_CO), R_ROW(IC_CO+1)-1
                     IC = R_COL(ICOL)
                     AUX = AUX + DC_FI(IC) * R(ICOL)
                  ENDDO
                  FC_CO(IC_CO) = AUX
               ENDDO

            ENDDO

      END SELECT SELECT_COMPACT_MULTIGRID

END SELECT SELECT_SYSTEM

END SUBROUTINE SCARC_RESTRICTION


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Perform prolongation from coarser to finer grid in multigrid method
!!!    - 'CO' corresponds to coarser grid
!!!    - 'FI' corresponds to finer   grid
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_PROLONGATION (NVECTOR_CO, NVECTOR_FI, NL_CO, NL_FI)
INTEGER, INTENT(IN) :: NVECTOR_CO, NVECTOR_FI, NL_CO, NL_FI
INTEGER :: NM, ICOL, IC, I
REAL(EB), POINTER, DIMENSION(:,:,:) :: XB_CO, DB_FI
REAL(EB), POINTER, DIMENSION(:)     :: XC_CO, DC_FI, P
INTEGER , POINTER, DIMENSION(:)     :: P_ROW, P_COL
INTEGER , POINTER :: NX_CO, NY_CO, NZ_CO, NC_FI
INTEGER  :: NX_FI, NY_FI, NZ_FI
INTEGER  :: IX_FI, IY_FI, IZ_FI, IC_FI(8)
INTEGER  :: IX_CO, IY_CO, IZ_CO, IC_CO
REAL(EB) :: AUX

SELECT_SYSTEM: SELECT CASE (TYPE_SYSTEM)

   !!!----------------------------------------------------------------------------------------------------
   !!! Banded system
   !!!----------------------------------------------------------------------------------------------------
   CASE (NSCARC_SYSTEM_BANDED)

      SELECT_BANDED_MULTIGRID: SELECT CASE (TYPE_MULTIGRID)

         !!!------------------------- Geometric multigrid -----------------------------------------------
         CASE (NSCARC_MULTIGRID_GEOMETRIC)

            DO NM = NMESHES_MIN, NMESHES_MAX
            
               NX_CO => SCARC(NM)%BANDED(NL_CO)%NX
               NY_CO => SCARC(NM)%BANDED(NL_CO)%NY
               NZ_CO => SCARC(NM)%BANDED(NL_CO)%NZ
            
               XB_CO  => POINT_TO_BVECTOR(NVECTOR_CO, NM, NL_CO)
               DB_FI  => POINT_TO_BVECTOR(NVECTOR_FI, NM, NL_FI)
      
               SELECT_BANDED_DIMENSION: SELECT CASE (TYPE_DIMENSION)

                  !!!----------------- 2D -------------------
                  CASE (NSCARC_DIMENSION_TWO)
               
                     DO IZ_CO = 1, NZ_CO
                        DO IX_CO = 1, NX_CO
               
                           IX_FI = 2*IX_CO
                           IY_FI = 1
                           IZ_FI = 2*IZ_CO
               
                           DB_FI(IX_FI-1, 1, IZ_FI-1) = XB_CO(IX_CO, 1, IZ_CO)
                           DB_FI(IX_FI  , 1, IZ_FI-1) = XB_CO(IX_CO, 1, IZ_CO)
                           DB_FI(IX_FI-1, 1, IZ_FI  ) = XB_CO(IX_CO, 1, IZ_CO)
                           DB_FI(IX_FI  , 1, IZ_FI  ) = XB_CO(IX_CO, 1, IZ_CO)
               
                        ENDDO
                     ENDDO
               
                  !!!----------------- 3D -------------------
                  CASE (NSCARC_DIMENSION_THREE)
                     
                     DO IZ_CO = 1, NZ_CO
                        DO IY_CO = 1, NY_CO
                           DO IX_CO = 1, NX_CO
               
                              IX_FI = 2*IX_CO
                              IY_FI = 2*IY_CO
                              IZ_FI = 2*IZ_CO
               
                              DB_FI(IX_FI-1, IY_FI-1, IZ_FI-1) = XB_CO(IX_CO, IY_CO, IZ_CO)
                              DB_FI(IX_FI  , IY_FI-1, IZ_FI-1) = XB_CO(IX_CO, IY_CO, IZ_CO)
                              DB_FI(IX_FI-1, IY_FI  , IZ_FI-1) = XB_CO(IX_CO, IY_CO, IZ_CO)
                              DB_FI(IX_FI  , IY_FI  , IZ_FI-1) = XB_CO(IX_CO, IY_CO, IZ_CO)
                              DB_FI(IX_FI-1, IY_FI-1, IZ_FI  ) = XB_CO(IX_CO, IY_CO, IZ_CO)
                              DB_FI(IX_FI  , IY_FI-1, IZ_FI  ) = XB_CO(IX_CO, IY_CO, IZ_CO)
                              DB_FI(IX_FI-1, IY_FI  , IZ_FI  ) = XB_CO(IX_CO, IY_CO, IZ_CO)
                              DB_FI(IX_FI  , IY_FI  , IZ_FI  ) = XB_CO(IX_CO, IY_CO, IZ_CO)
               
                           ENDDO
                        ENDDO
                     ENDDO
                        
               END SELECT SELECT_BANDED_DIMENSION

            ENDDO
      
         !!!------------------------- Geometric multigrid -----------------------------------------------
         CASE (NSCARC_MULTIGRID_ALGEBRAIC)
            WRITE(*,*) 'No banded restriction for algebraic multigrid available, stopping program!'
            STOP

      END SELECT SELECT_BANDED_MULTIGRID

      
   !!!----------------------------------------------------------------------------------------------------
   !!! Compact system
   !!!----------------------------------------------------------------------------------------------------
   CASE (NSCARC_SYSTEM_COMPACT)

      SELECT_COMPACT_MULTIGRID: SELECT CASE (TYPE_MULTIGRID)

         !!!------------------------- Geometric multigrid -----------------------------------------------
         CASE (NSCARC_MULTIGRID_GEOMETRIC)

            DO NM = NMESHES_MIN, NMESHES_MAX
         
               NX_CO => SCARC(NM)%COMPACT(NL_CO)%NX
               NY_CO => SCARC(NM)%COMPACT(NL_CO)%NY
               NZ_CO => SCARC(NM)%COMPACT(NL_CO)%NZ
         
               NX_FI = 2*NX_CO
               NY_FI = 2*NY_CO
               NZ_FI = 2*NZ_CO
   
               XC_CO  => POINT_TO_CVECTOR(NVECTOR_CO, NM, NL_CO)
               DC_FI  => POINT_TO_CVECTOR(NVECTOR_FI, NM, NL_FI)
   
               SELECT_COMPACT_DIMENSION: SELECT CASE (TYPE_DIMENSION)
               
                  !!!----------------- 2D -------------------
                  CASE (NSCARC_DIMENSION_TWO)
               
                     DO IZ_CO = 1, NZ_CO
                        DO IX_CO = 1, NX_CO
               
                           IX_FI = 2*IX_CO
                           IY_FI = 1
                           IZ_FI = 2*IZ_CO
               
                           IC_CO = (IZ_CO-1)*NX_CO + IX_CO
               
                           IC_FI(1) = (IZ_FI-2)*NX_FI + IX_FI - 1
                           IC_FI(2) = (IZ_FI-2)*NX_FI + IX_FI   
                           IC_FI(3) = (IZ_FI-1)*NX_FI + IX_FI - 1
                           IC_FI(4) = (IZ_FI-1)*NX_FI + IX_FI    
               
                           DO I = 1, 4
                              DC_FI(IC_FI(I)) = XC_CO(IC_CO)
                           ENDDO
               
                        ENDDO
                     ENDDO
               
                              
                  !!!----------------- 3D -------------------
                  CASE (NSCARC_DIMENSION_THREE)
                     
                     DO IZ_CO = 1, NZ_CO
                        DO IY_CO = 1, NY_CO
                           DO IX_CO = 1, NX_CO
               
                              IX_FI = 2*IX_CO
                              IY_FI = 2*IY_CO
                              IZ_FI = 2*IZ_CO
               
                              IC_CO    = (IZ_CO-1)*NX_CO*NY_CO + (IY_CO-1)*NX_CO + IX_CO
               
                              IC_FI(1) = (IZ_FI-2)*NX_FI*NY_FI + (IY_FI-2)*NX_FI + IX_FI - 1
                              IC_FI(2) = (IZ_FI-2)*NX_FI*NY_FI + (IY_FI-2)*NX_FI + IX_FI    
                              IC_FI(3) = (IZ_FI-2)*NX_FI*NY_FI + (IY_FI-1)*NX_FI + IX_FI - 1
                              IC_FI(4) = (IZ_FI-2)*NX_FI*NY_FI + (IY_FI-1)*NX_FI + IX_FI    
                              IC_FI(5) = (IZ_FI-1)*NX_FI*NY_FI + (IY_FI-2)*NX_FI + IX_FI - 1
                              IC_FI(6) = (IZ_FI-1)*NX_FI*NY_FI + (IY_FI-2)*NX_FI + IX_FI    
                              IC_FI(7) = (IZ_FI-1)*NX_FI*NY_FI + (IY_FI-1)*NX_FI + IX_FI - 1
                              IC_FI(8) = (IZ_FI-1)*NX_FI*NY_FI + (IY_FI-1)*NX_FI + IX_FI    
               
                              DO I = 1, 8
                                 DC_FI(IC_FI(I)) = XC_CO(IC_CO)
                              ENDDO
               
                           ENDDO
                        ENDDO
                     ENDDO
   
               END SELECT SELECT_COMPACT_DIMENSION

            ENDDO

         !!!------------------------- Algebraic multigrid -----------------------------------------------
         CASE (NSCARC_MULTIGRID_ALGEBRAIC)

            DO NM = NMESHES_MIN, NMESHES_MAX

               XC_CO  => POINT_TO_CVECTOR(NVECTOR_CO, NM, NL_CO)
               DC_FI  => POINT_TO_CVECTOR(NVECTOR_FI, NM, NL_FI)

               NC_FI => SCARC(NM)%COMPACT(NL_FI)%NC

               P     => SCARC(NM)%COMPACT(NL_FI)%P
               P_ROW => SCARC(NM)%COMPACT(NL_FI)%P_ROW
               P_COL => SCARC(NM)%COMPACT(NL_FI)%P_COL
   
               DO IC = 1, NC_FI
                  AUX = 0.0_EB
                  DO ICOL = P_ROW(IC), P_ROW(IC+1)-1
                     IC_CO = P_COL(ICOL)
                     AUX = XC_CO(IC_CO) * P(ICOL)
                  ENDDO
                  DC_FI(IC) = AUX
               ENDDO

            ENDDO
            
      END SELECT SELECT_COMPACT_MULTIGRID

END SELECT SELECT_SYSTEM

END SUBROUTINE SCARC_PROLONGATION


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Finalize data - banded storage technique
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_TERMINATE_SOLVER(NL)
INTEGER, INTENT(IN) :: NL
INTEGER :: NM, I, J, K, IC
REAL(EB), POINTER, DIMENSION(:,:,:) :: HP
TYPE (MESH_TYPE), POINTER :: M
TYPE (SCARC_BANDED_TYPE) , POINTER :: SB
TYPE (SCARC_COMPACT_TYPE), POINTER :: SC

SELECT CASE (TYPE_SYSTEM)
 
   !!!----------------------------------------------------------------------------------------------------
   !!! Banded system
   !!!----------------------------------------------------------------------------------------------------
   CASE (NSCARC_SYSTEM_BANDED)

      DO NM = NMESHES_MIN, NMESHES_MAX
      
         M  => MESHES(NM)
         SB => SCARC(NM)%BANDED(NL)
      
         IF (PREDICTOR) THEN
            HP => M%H
         ELSE
            HP => M%HS
         ENDIF
      
         !!! Overwrite internal values of H or HS by corresponding data of X
         HP(1:M%IBAR, 1:M%JBAR, 1:M%KBAR) = SB%X(1:M%IBAR, 1:M%JBAR, 1:M%KBAR)
        
      ENDDO 
      
   !!!----------------------------------------------------------------------------------------------------
   !!! Compact system
   !!!----------------------------------------------------------------------------------------------------
   CASE (NSCARC_SYSTEM_COMPACT)

      DO NM = NMESHES_MIN, NMESHES_MAX
      
         M  => MESHES(NM)
         SC => SCARC(NM)%COMPACT(NL)
      
         IF (PREDICTOR) THEN
            HP => M%H
         ELSE
            HP => M%HS
         ENDIF

         DO K = 1, M%KBAR
            DO J = 1, M%JBAR
               DO I = 1, M%IBAR
                  IC = (K-1) * M%IBAR * M%JBAR + (J-1) * M%IBAR + I
                  HP(I, J, K) = SC%X(IC)
               ENDDO
            ENDDO
         ENDDO
      
      ENDDO

END SELECT

END SUBROUTINE SCARC_TERMINATE_SOLVER


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Set correct boundary values at external and internal boundaries
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_UPDATE_GHOSTCELLS(NL)
INTEGER, INTENT(IN) :: NL
INTEGER :: NM, IW, IOR0, I0, J0, K0, I1, J1, K1
REAL(EB), POINTER, DIMENSION(:,:,:) :: HP=>NULL()
TYPE (MESH_TYPE), POINTER :: M


!!!---------------------------------------------------------------------------------------------------
!!! Adjust ghost values along external boundaries according to boundary arrays BXS, BXF, ...
!!!---------------------------------------------------------------------------------------------------
DO NM = NMESHES_MIN, NMESHES_MAX

   !!! point to correct pressure vector on mesh 'NM'
   M => MESHES(NM)
   IF (PREDICTOR) THEN
      HP => M%H
   ELSE
      HP => M%HS
   ENDIF
   
   !!! compute ghost cell values
   WALL_CELL_LOOP: DO IW = 1, M%N_EXTERNAL_WALL_CELLS
   
      I0 = M%OMESH(NM)%IJKW(1,IW)!****CHECK
      J0 = M%OMESH(NM)%IJKW(2,IW)!****CHECK
      K0 = M%OMESH(NM)%IJKW(3,IW)!****CHECK
   
      I1 = M%OMESH(NM)%IJKW(6,IW)!****CHECK
      J1 = M%OMESH(NM)%IJKW(7,IW)!****CHECK
      K1 = M%OMESH(NM)%IJKW(8,IW)!****CHECK
   
      IOR0 = M%OMESH(NM)%IJKW(4,IW)!****CHECK
   
      IF (IOR0 == 1) THEN
         IF (M%WALL(IW)%PRESSURE_BC_INDEX==DIRICHLET) THEN
            HP(I0,J0,K0) = -HP(I1,J1,K1) + 2.0_EB * M%BXS(J1,K1)
         ELSE IF (M%WALL(IW)%PRESSURE_BC_INDEX==NEUMANN) THEN
            HP(I0,J0,K0) =  HP(I1,J1,K1) - DXI *M%BXS(J1,K1)
         ENDIF
      ELSE IF (IOR0 == -1) THEN
         IF (M%WALL(IW)%PRESSURE_BC_INDEX==DIRICHLET) THEN
            HP(I0,J0,K0) = -HP(I1,J1,K1) + 2.0_EB * M%BXF(J1,K1)
         ELSE IF (M%WALL(IW)%PRESSURE_BC_INDEX==NEUMANN) THEN
            HP(I0,J0,K0) =  HP(I1,J1,K1) + DXI *M%BXF(J1,K1)
         ENDIF
      ELSE IF (IOR0 == 2) THEN
         IF (M%WALL(IW)%PRESSURE_BC_INDEX==DIRICHLET) THEN
            HP(I0,J0,K0) = -HP(I1,J1,K1) + 2.0_EB * M%BYS(I1,K1)
         ELSE IF (M%WALL(IW)%PRESSURE_BC_INDEX==NEUMANN) THEN
            HP(I0,J0,K0) =  HP(I1,J1,K1) - DETA *M%BYS(I1,K1)
         ENDIF
      ELSE IF (IOR0 == -2) THEN
         IF (M%WALL(IW)%PRESSURE_BC_INDEX==DIRICHLET) THEN
            HP(I0,J0,K0) = -HP(I1,J1,K1) + 2.0_EB * M%BYF(I1,K1)
         ELSE IF (M%WALL(IW)%PRESSURE_BC_INDEX==NEUMANN) THEN
            HP(I0,J0,K0) =  HP(I1,J1,K1) + DETA *M%BYF(I1,K1)
         ENDIF
      ELSE IF (IOR0 == 3) THEN
         IF (M%WALL(IW)%PRESSURE_BC_INDEX==DIRICHLET) THEN
            HP(I0,J0,K0) = -HP(I1,J1,K1) + 2.0_EB * M%BZS(I1,J1)
         ELSE IF (M%WALL(IW)%PRESSURE_BC_INDEX==NEUMANN) THEN
            HP(I0,J0,K0) =  HP(I1,J1,K1) - DZETA *M%BZS(I1,J1)
         ENDIF
      ELSE IF (IOR0 == -3) THEN
         IF (M%WALL(IW)%PRESSURE_BC_INDEX==DIRICHLET) THEN
            HP(I0,J0,K0) = -HP(I1,J1,K1) + 2.0_EB * M%BZF(I1,J1)
         ELSE IF (M%WALL(IW)%PRESSURE_BC_INDEX==NEUMANN) THEN
            HP(I0,J0,K0) =  HP(I1,J1,K1) + DZETA *M%BZF(I1,J1)
         ENDIF
      ENDIF
   ENDDO WALL_CELL_LOOP
   
ENDDO

!!!---------------------------------------------------------------------------------------------------
!!! Perform data exchange to achieve consistency of ghost values along internal boundaries 
!!!---------------------------------------------------------------------------------------------------
IF (NMESHES > 1) CALL SCARC_EXCHANGE(NSCARC_EXCHANGE_BDRY, NL)

END SUBROUTINE SCARC_UPDATE_GHOSTCELLS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  Perform data exchange corresponding to requested exchange type (call receive and send-routines)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_EXCHANGE (NTYPE, NL)
INTEGER, INTENT(IN):: NTYPE, NL

NREQ_SCARC    = 0
TYPE_EXCHANGE = NTYPE

CALL SCARC_RECEIVE(NL)
CALL SCARC_SEND(NL)

END SUBROUTINE SCARC_EXCHANGE
 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  Receive data from neighbors (corresponds to POST_RECEIVES)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_RECEIVE (NL)

INTEGER, INTENT(IN) :: NL
INTEGER :: NM, NOM, IERR, NLEN, KAMG(10)
TYPE (SCARC_TYPE)         , POINTER ::  S
TYPE (OSCARC_TYPE)        , POINTER ::  OS
TYPE (OSCARC_BANDED_TYPE) , POINTER ::  OSB
TYPE (OSCARC_COMPACT_TYPE), POINTER ::  OSC

IERR=0
RECEIVE_MESH_LOOP: DO NM = NMESHES_MIN, NMESHES_MAX
   RECEIVE_OMESH_LOOP: DO NOM = 1, NMESHES
    
      SNODE = PROCESS(NOM)
      IF (PROCESS(NM)==SNODE) CYCLE RECEIVE_OMESH_LOOP

      S   => SCARC(NM)                            ! corresponds to M
      OS  => SCARC(NM)%OSCARC(NOM)                ! corresponds to M3

      IF (OS%NICMAX_S==0 .AND. OS%NICMAX_R==0) CYCLE RECEIVE_OMESH_LOOP

      SELECT_EXCHANGE_TYPE: SELECT CASE (TYPE_EXCHANGE)

         !!!-------------------------------------------------------------------------------------------
         !!! Initialize communication structures for the receiving of data
         !!!-------------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_IJKW)

            SELECT CASE (TYPE_SYSTEM)
               CASE (NSCARC_SYSTEM_BANDED)
                  OSB => OS%BANDED(NL)
                  IF (NL/=NLEVEL_MIN) THEN  ! for lower levels get neighboring IJKW's
                     IF (USE_MPI) THEN
                        NREQ_SCARC = NREQ_SCARC+1
                        CALL MPI_IRECV(OSB%IJKW(1,1),15*OSB%NW,MPI_INTEGER,SNODE, &
                                       TAG_SCARC,MPI_COMM_WORLD,REQ_SCARC(NREQ_SCARC),IERR)
                     ENDIF
                  ELSE                                ! on maximum level neighboring IJKW already exists
                     OSB%IJKW => MESHES(NM)%OMESH(NOM)%IJKW     
                  ENDIF
               CASE (NSCARC_SYSTEM_COMPACT)
                  OSC => OS%COMPACT(NL)
                  IF (NL/=NLEVEL_MIN) THEN  ! for lower levels get neighboring IJKW's
                     IF (USE_MPI) THEN
                        NREQ_SCARC = NREQ_SCARC+1
                        CALL MPI_IRECV(OSC%IJKW(1,1),15*OSC%NW,MPI_INTEGER,SNODE, &
                                       TAG_SCARC,MPI_COMM_WORLD,REQ_SCARC(NREQ_SCARC),IERR)
                     ENDIF
                  ELSE                                ! on maximum level neighboring IJKW already exists
                     OSC%IJKW => MESHES(NM)%OMESH(NOM)%IJKW     
                  ENDIF
            END SELECT
   
            IF (OS%NIC_R > 0) THEN
               !NLEN = (MAX(OS%NIC_R, OS%NIC_S)+2)*2+1
               SELECT CASE (TYPE_DIMENSION)
                  CASE (NSCARC_DIMENSION_TWO)
                     !NLEN = MAX(OS%NIC_R, OS%NIC_S)*5+1
                     NLEN = MAX(OS%NIC_R, OS%NIC_S)*30+1              !!! ONLY TEMPORARILY
                  CASE (NSCARC_DIMENSION_THREE)
                     !NLEN = MAX(OS%NIC_R, OS%NIC_S)*7+1
                     NLEN = MAX(OS%NIC_R, OS%NIC_S)*80+1
               END SELECT
               ALLOCATE (OS%RECV_BUF(NLEN))
               OS%RECV_BUF = 0.0_EB
            ENDIF

         !!!-------------------------------------------------------------------------------------------
         !!! Exchange information about neighboring grid dimensions
         !!!-------------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_GRID)

            OSC => OS%COMPACT(NL)
            IF (USE_MPI) THEN
               NREQ_SCARC = NREQ_SCARC+1
               CALL MPI_IRECV(KAMG(1),5,MPI_INTEGER,SNODE, &
                              TAG_SCARC,MPI_COMM_WORLD,REQ_SCARC(NREQ_SCARC),IERR)
               OSC%NX = KAMG(1)
               OSC%NY = KAMG(2)
               OSC%NZ = KAMG(3)
               OSC%NC = KAMG(4)
               OSC%NW = KAMG(5)
            ELSE
               OSC%NX = SCARC(NM)%OSCARC(NOM)%COMPACT(NL)%NX
               OSC%NY = SCARC(NM)%OSCARC(NOM)%COMPACT(NL)%NY
               OSC%NZ = SCARC(NM)%OSCARC(NOM)%COMPACT(NL)%NZ
               OSC%NC = SCARC(NM)%OSCARC(NOM)%COMPACT(NL)%NC
               OSC%NW = SCARC(NM)%OSCARC(NOM)%COMPACT(NL)%NW
            ENDIF

         !!!-------------------------------------------------------------------------------------------
         !!! Exchange number of neighboring cells for AMG method (compact type only)
         !!!-------------------------------------------------------------------------------------------
         CASE (NSCARC_EXCHANGE_AMG)

            OSC => OS%COMPACT(NL)
            IF (USE_MPI) THEN
               NREQ_SCARC = NREQ_SCARC+1
               CALL MPI_IRECV(KAMG(1),3,MPI_INTEGER,SNODE, &
                              TAG_SCARC,MPI_COMM_WORLD,REQ_SCARC(NREQ_SCARC),IERR)
               OSC%NC  = KAMG(1)
               OSC%NCE = KAMG(2)
               OSC%NW  = KAMG(3)
            ELSE
               OSC%NC  = SCARC(NM)%OSCARC(NOM)%COMPACT(NL)%NC
               OSC%NCE = SCARC(NM)%OSCARC(NOM)%COMPACT(NL)%NCE
               OSC%NW  = SCARC(NM)%OSCARC(NOM)%COMPACT(NL)%NW
            ENDIF


         !!!-------------------------------------------------------------------------------------------
         !!! Perform exchanges for 
         !!!    - internal values for matrix-vector multiplication 
         !!!    - internal boundariy values
         !!!    - internal subdiagonal matrix values
         !!!    - internal subdiagonal or ghost matrix values
         !!!    - internal measure/celltype values
         !!!-------------------------------------------------------------------------------------------
         CASE DEFAULT

            NREQ_SCARC = NREQ_SCARC+1
            IF (USE_MPI) CALL MPI_IRECV(OS%RECV_BUF(1),SIZE(OS%RECV_BUF),MPI_DOUBLE_PRECISION,&
                                        SNODE,TAG_SCARC,MPI_COMM_WORLD,REQ_SCARC(NREQ_SCARC),IERR)
!WRITE(SCARC_LU,*) '============= RECV VECTOR:', SIZE(OS%RECV_BUF)
!WRITE(SCARC_LU,'(10f8.3)') OS%RECV_BUF(1:50)


      END SELECT SELECT_EXCHANGE_TYPE
   ENDDO RECEIVE_OMESH_LOOP
ENDDO RECEIVE_MESH_LOOP

END SUBROUTINE SCARC_RECEIVE
 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Send data to neighbors (corresponds to MESH_EXCHANGE)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SEND (NL)
INTEGER, INTENT(IN) :: NL
INTEGER :: NM, NOM
INTEGER :: IERR
INTEGER :: NLEN, KAMG(5)
INTEGER, POINTER:: NW, NX, NY, NC
REAL(EB), POINTER, DIMENSION(:)    :: BUFFER
INTEGER , POINTER, DIMENSION(:,:)  :: IJKW
TYPE (SCARC_TYPE)         , POINTER :: S 
TYPE (SCARC_BANDED_TYPE)  , POINTER :: SB, SOB
TYPE (SCARC_COMPACT_TYPE) , POINTER :: SC, SOC
TYPE (OSCARC_TYPE)        , POINTER :: OS, OSO
TYPE (OSCARC_BANDED_TYPE) , POINTER :: OSB
TYPE (OSCARC_COMPACT_TYPE), POINTER :: OSC

IERR = 0

EXCHANGE_SEND_LOOP1: DO NM = NMESHES_MIN, NMESHES_MAX

   IF (PROCESS(NM)/=MYID)  CYCLE EXCHANGE_SEND_LOOP1

   EXCHANGE_RECV_LOOP1: DO NOM = 1, NMESHES
    
      SNODE = PROCESS(NOM)
      RNODE = PROCESS(NM)

      S   => SCARC(NM)                           ! corresponds to M
      OS  => SCARC(NM)%OSCARC(NOM)               ! corresponds to M3

      IF (OS%NICMAX_S == 0 .AND. OS%NICMAX_R == 0) CYCLE EXCHANGE_RECV_LOOP1



      !!!-------------------------------------------------------------------------------------------
      !!! Initialize communication structures for the sending of data
      !!!-------------------------------------------------------------------------------------------
      IF (TYPE_EXCHANGE == NSCARC_EXCHANGE_IJKW) THEN

         ! on max level   : take M%IJKW from neighbors (ScaRC pointers already defined, nothing to do ...)
         ! on lover levels: send own IJKW to neighbors
         SELECT CASE (TYPE_SYSTEM)

            CASE (NSCARC_SYSTEM_BANDED)
               SB  =>  S%BANDED(NL)                       ! corresponds to M  for the level 'NL'
               OSB => OS%BANDED(NL)                       ! corresponds to M3 for the level 'NL'
               IF (RNODE /= SNODE) THEN
                  IF (NL /= NLEVEL_MIN) THEN
                     NREQ_SCARC = NREQ_SCARC+1
                     IF (USE_MPI) CALL MPI_ISEND(SB%IJKW(1,1),15*SB%NW,MPI_INTEGER,SNODE, &
                                                 TAG_SCARC,MPI_COMM_WORLD,REQ_SCARC(NREQ_SCARC),IERR)
                  ENDIF
               ELSE
                  SOB      => SCARC(NOM)%BANDED(NL)
                  OSB%IJKW => SOB%IJKW(:,1:SOB%NW)
               ENDIF

            CASE (NSCARC_SYSTEM_COMPACT)
               SC  =>  S%COMPACT(NL)                       ! corresponds to M  for the level 'NL'
               OSC => OS%COMPACT(NL)                       ! corresponds to M3 for the level 'NL'
               IF (RNODE /= SNODE) THEN
                  IF (NL /= NLEVEL_MIN) THEN
                     NREQ_SCARC = NREQ_SCARC+1
                     IF (USE_MPI) CALL MPI_ISEND(SC%IJKW(1,1),15*SC%NW,MPI_INTEGER,SNODE, &
                                                 TAG_SCARC,MPI_COMM_WORLD,REQ_SCARC(NREQ_SCARC),IERR)
                  ENDIF
               ELSE
                  SOC      => SCARC(NOM)%COMPACT(NL)
                  OSC%IJKW => SOC%IJKW(:,1:SOC%NW)
               ENDIF

         END SELECT

         !NLEN=(MAX(OS%NIC_R, OS%NIC_S)+2)*2+1   
         SELECT CASE (TYPE_DIMENSION)
            CASE (NSCARC_DIMENSION_TWO)
               !NLEN = MAX(OS%NIC_R, OS%NIC_S)*5+1
               NLEN = MAX(OS%NIC_R, OS%NIC_S)*30+1
            CASE (NSCARC_DIMENSION_THREE)
               !NLEN = MAX(OS%NIC_R, OS%NIC_S)*7+1
               NLEN = MAX(OS%NIC_R, OS%NIC_S)*80+1
         END SELECT
         ALLOCATE (OS%SEND_BUF(NLEN))
         OS%SEND_BUF = 0.0_EB

      !!!-------------------------------------------------------------------------------------------
      !!! Exchange neighboring grid information
      !!!-------------------------------------------------------------------------------------------
      ELSE IF (TYPE_EXCHANGE == NSCARC_EXCHANGE_GRID) THEN

         SC  =>  S%COMPACT(NL)
         OSC => OS%COMPACT(NL)
         IF (RNODE /= SNODE) THEN
            NREQ_SCARC = NREQ_SCARC+1
            KAMG(1)=SC%NX
            KAMG(2)=SC%NY
            KAMG(3)=SC%NZ
            KAMG(4)=SC%NC
            KAMG(4)=SC%NW
            IF (USE_MPI) CALL MPI_ISEND(KAMG(1),5,MPI_INTEGER,SNODE, &
                                        TAG_SCARC,MPI_COMM_WORLD,REQ_SCARC(NREQ_SCARC),IERR)
         ELSE
            SOC    => SCARC(NOM)%COMPACT(NL)
            OSC%NX =  SOC%NX
            OSC%NY =  SOC%NY
            OSC%NZ =  SOC%NZ
            OSC%NC =  SOC%NC
            OSC%NW =  SOC%NW
!WRITE(SCARC_LU,'(a,i3,a,i3,a,i3,a,i3)') 'SCARC(,',NM,')%OSCARC(,',NOM,')%COMPACT(',NL,')%NC=',OSC%NC
!WRITE(SCARC_LU,'(a,i3,a,i3,a,i3,a,i3)') 'SCARC(,',NM,')%OSCARC(,',NOM,')%COMPACT(',NL,')%NC=',OSC%NCE
!WRITE(SCARC_LU,'(a,i3,a,i3,a,i3,a,i3)') 'SCARC(,',NM,')%OSCARC(,',NOM,')%COMPACT(',NL,')%NC=',OSC%NW
         ENDIF

      !!!-------------------------------------------------------------------------------------------
      !!! Exchange number of neighboring cells for AMG method (compact type only)
      !!!-------------------------------------------------------------------------------------------
      ELSE IF (TYPE_EXCHANGE == NSCARC_EXCHANGE_AMG) THEN

         SC  =>  S%COMPACT(NL)
         OSC => OS%COMPACT(NL)
         IF (RNODE /= SNODE) THEN
            NREQ_SCARC = NREQ_SCARC+1
            KAMG(1)=SC%NC
            KAMG(2)=SC%NCE
            KAMG(3)=SC%NW
            IF (USE_MPI) CALL MPI_ISEND(KAMG(1),3,MPI_INTEGER,SNODE, &
                                        TAG_SCARC,MPI_COMM_WORLD,REQ_SCARC(NREQ_SCARC),IERR)
         ELSE
            SOC    => SCARC(NOM)%COMPACT(NL)
            OSC%NC  =  SOC%NC
            OSC%NCE =  SOC%NCE
            OSC%NW  =  SOC%NW
!WRITE(SCARC_LU,'(a,i3,a,i3,a,i3,a,i3)') 'SCARC(,',NM,')%OSCARC(,',NOM,')%COMPACT(',NL,')%NC=',OSC%NC
!WRITE(SCARC_LU,'(a,i3,a,i3,a,i3,a,i3)') 'SCARC(,',NM,')%OSCARC(,',NOM,')%COMPACT(',NL,')%NC=',OSC%NCE
!WRITE(SCARC_LU,'(a,i3,a,i3,a,i3,a,i3)') 'SCARC(,',NM,')%OSCARC(,',NOM,')%COMPACT(',NL,')%NC=',OSC%NW
         ENDIF


      !!!-------------------------------------------------------------------------------------------
      !!! Exchange data along internal boundaries corresponding to requested exchange type
      !!!-------------------------------------------------------------------------------------------
      ELSE

         NLEN = 0
         BUFFER => OS%SEND_BUF

         SELECT CASE(TYPE_SYSTEM)

            !!! -------------------------- Banded system  ------------------------------------------
            CASE (NSCARC_SYSTEM_BANDED)

               SB  =>  S%BANDED(NL)                       ! corresponds to M  for the level 'NL'
               OSB => OS%BANDED(NL)                       ! corresponds to M3 for the level 'NL'

               IJKW   => OSB%IJKW
               NW     => OSB%NW
               NX     => SB%NX
               NY     => SB%NY
               NC     => SB%NC

               SELECT_EXCHANGE_BANDED: SELECT CASE (TYPE_EXCHANGE)

                  CASE (NSCARC_EXCHANGE_VECTOR) 
                     CALL SCARC_PACK_BVECTOR_REAL(BUFFER, IJKW, TYPE_VECTOR, NLEN, NW, NM, NL)

                  CASE (NSCARC_EXCHANGE_BDRY) 
                     IF (PREDICTOR) THEN
                        CALL SCARC_PACK_BVECTOR_REAL(BUFFER, IJKW, NSCARC_VECTOR_H , NLEN, NW, NM, NL)
                     ELSE
                        CALL SCARC_PACK_BVECTOR_REAL(BUFFER, IJKW, NSCARC_VECTOR_HS, NLEN, NW, NM, NL)
                     ENDIF

                  CASE (NSCARC_EXCHANGE_MATRIX) 
                     CALL SCARC_PACK_BMATRIX(BUFFER, IJKW, TYPE_MATRIX, NLEN, NW, NM, NL)
            
                  CASE DEFAULT
                     WRITE (*,*) 'BANDED: TYPE_EXCHANGE =',TYPE_EXCHANGE,' NOT ALLOWED,/,STOPPING PROGRAM !'
                     STOP

               END SELECT SELECT_EXCHANGE_BANDED

            !!! -------------------------- Compact system  -----------------------------------------
            CASE (NSCARC_SYSTEM_COMPACT)

               SC  =>  S%COMPACT(NL)                       ! corresponds to M  for the level 'NL'
               OSC => OS%COMPACT(NL)                       ! corresponds to M3 for the level 'NL'

               IJKW   => OSC%IJKW
               NW     => OSC%NW
               NX     => SC%NX
               NY     => SC%NY
               NC     => SC%NC

!WRITE(SCARC_LU,*) TYPE_EXCHANGE, NW, NX, NY, NC

               SELECT_EXCHANGE_COMPACT: SELECT CASE (TYPE_EXCHANGE)

                  CASE (NSCARC_EXCHANGE_VECTOR) 
                     CALL SCARC_PACK_CVECTOR_REAL(BUFFER, IJKW, TYPE_VECTOR, NLEN, NW, NM, NL)
            
                  CASE (NSCARC_EXCHANGE_BDRY) 
                     IF (PREDICTOR) THEN
                        CALL SCARC_PACK_BVECTOR_REAL(BUFFER, IJKW, NSCARC_VECTOR_H , NLEN, NW, NM, NL)
                     ELSE
                        CALL SCARC_PACK_BVECTOR_REAL(BUFFER, IJKW, NSCARC_VECTOR_HS, NLEN, NW, NM, NL)
                     ENDIF

                  CASE (NSCARC_EXCHANGE_MATRIX) 
                     CALL SCARC_PACK_CMATRIX(BUFFER, IJKW, TYPE_MATRIX, NLEN, NW, NM, NL)
            
                  CASE (NSCARC_EXCHANGE_MEASURE) 
                     CALL SCARC_PACK_CVECTOR_REAL(BUFFER, IJKW, NSCARC_VECTOR_MEASURE, NLEN, NW, NM, NL)
            
                  CASE (NSCARC_EXCHANGE_CELLTYPE) 
                     CALL SCARC_PACK_CVECTOR_INT (BUFFER, IJKW, NSCARC_VECTOR_CELLTYPE, NLEN, NW, NM, NL)
      
                  CASE (NSCARC_EXCHANGE_WEIGHTS) 
IF (TYPE_DEBUG> NSCARC_DEBUG_NONE) WRITE(SCARC_LU,*) 'BEFORE PACK_CMATRIX NM=',NM,': NOM=',NOM
                     CALL SCARC_PACK_CMATRIX(BUFFER, IJKW, NSCARC_MATRIX_PROLONGATION, NLEN, NW, NM, NL)
IF (TYPE_DEBUG> NSCARC_DEBUG_NONE) WRITE(SCARC_LU,*) 'AFTER PACK_CMATRIX'

                  CASE DEFAULT
                     WRITE (*,*) 'COMPACT: TYPE_EXCHANGE =',TYPE_EXCHANGE,' NOT ALLOWED,/,STOPPING PROGRAM !'
                     STOP
         
               END SELECT SELECT_EXCHANGE_COMPACT

          END SELECT

         !!! Finally exchange send buffer with corresponding neighbors
         IF (RNODE/=SNODE) THEN
            NREQ_SCARC=NREQ_SCARC+1
!WRITE(SCARC_LU,*) '============= SEND VECTOR:', NLEN
!WRITE(SCARC_LU,'(10f8.3)') BUFFER(1:50)
            IF (USE_MPI) CALL MPI_ISEND(BUFFER, NLEN, MPI_DOUBLE_PRECISION, SNODE, &
                                        TAG_SCARC, MPI_COMM_WORLD, REQ_SCARC(NREQ_SCARC), IERR)
         ENDIF

      ENDIF

   ENDDO EXCHANGE_RECV_LOOP1

ENDDO EXCHANGE_SEND_LOOP1


!!!----------------------------------------------------------------------------------------------------
!!! Information from Mesh NM is received by Mesh NOM  (NOM receiver, NM sender)
!!!----------------------------------------------------------------------------------------------------
IF (USE_MPI.AND.NREQ_SCARC/=0) CALL MPI_WAITALL(NREQ_SCARC,REQ_SCARC(1:NREQ_SCARC),MPI_STATUS_IGNORE,IERR)


!!!----------------------------------------------------------------------------------------------------
!!! Extract communication data from corresponding RECEIVE-buffers
!!!----------------------------------------------------------------------------------------------------
EXCHANGE_SEND_LOOP2: DO NOM = NMESHES_MIN, NMESHES_MAX

   EXCHANGE_RECV_LOOP2: DO NM=1,NMESHES
    
      SNODE  = PROCESS(NOM)
      RNODE  = PROCESS(NM)

      OSO => SCARC(NOM)%OSCARC(NM)

      EXCHANGE_RECV_IF: IF (OSO%NICMAX_S/=0 .AND. OSO%NICMAX_R/=0) THEN

         IF (RNODE/=SNODE) THEN
            BUFFER => OSO%RECV_BUF
         ELSE
            OS => SCARC(NM)%OSCARC(NOM)
            BUFFER => OS%SEND_BUF
         ENDIF

!WRITE(SCARC_LU,*) 'EXTRACTING NOM=',NOM,': NM=',NM
         SELECT CASE(TYPE_SYSTEM)

            !!! -------------------------- Banded system  ------------------------------------------
            CASE (NSCARC_SYSTEM_BANDED)

               SOB  => SCARC(NOM)%BANDED(NL)

               IJKW => SOB%IJKW
               NX   => SOB%NX
               NY   => SOB%NY
      
               SELECT_EXCHANGE_BANDED2: SELECT CASE (TYPE_EXCHANGE)
      
                  CASE (NSCARC_EXCHANGE_VECTOR) 
                     CALL SCARC_UNPACK_BVECTOR_REAL(BUFFER, IJKW, TYPE_VECTOR, NOM, NL)
     
                  CASE (NSCARC_EXCHANGE_BDRY)
                     IF (PREDICTOR) THEN
                        CALL SCARC_UNPACK_BVECTOR_REAL(BUFFER, IJKW, NSCARC_VECTOR_H, NOM, NL)
                     ELSE
                        CALL SCARC_UNPACK_BVECTOR_REAL(BUFFER, IJKW, NSCARC_VECTOR_HS, NOM, NL)
                     ENDIF
      
                  CASE (NSCARC_EXCHANGE_MATRIX)
                     CALL SCARC_UNPACK_BMATRIX(BUFFER, IJKW, TYPE_MATRIX, NOM, NL)
      
                END SELECT SELECT_EXCHANGE_BANDED2

            !!! -------------------------- Compact system  -----------------------------------------
            CASE (NSCARC_SYSTEM_COMPACT)

               SOC  => SCARC(NOM)%COMPACT(NL)

               IJKW => SOC%IJKW
               NX   => SOC%NX
               NY   => SOC%NY
      
               SELECT_EXCHANGE_COMPACT2: SELECT CASE (TYPE_EXCHANGE)
      
                  CASE (NSCARC_EXCHANGE_VECTOR) 
                     CALL SCARC_UNPACK_CVECTOR_REAL(BUFFER, IJKW, TYPE_VECTOR, NOM, NL)
      
                  CASE (NSCARC_EXCHANGE_BDRY)
                     IF (PREDICTOR) THEN
                        CALL SCARC_UNPACK_BVECTOR_REAL(BUFFER, IJKW, NSCARC_VECTOR_H, NOM, NL)
                     ELSE
                        CALL SCARC_UNPACK_BVECTOR_REAL(BUFFER, IJKW, NSCARC_VECTOR_HS, NOM, NL)
                     ENDIF
      
                  CASE (NSCARC_EXCHANGE_MATRIX)
!WRITE(SCARC_LU,*) 'NOM=',NOM,': UNPACK CMATRIX : NL=',NL
                     CALL SCARC_UNPACK_CMATRIX(BUFFER, IJKW, TYPE_MATRIX, NOM, NL)
      
                  CASE (NSCARC_EXCHANGE_MEASURE)
                     CALL SCARC_UNPACK_CVECTOR_REAL(BUFFER, IJKW, NSCARC_VECTOR_MEASURE, NOM, NL)
      
                  CASE (NSCARC_EXCHANGE_CELLTYPE)
                     CALL SCARC_UNPACK_CVECTOR_INT (BUFFER, IJKW, NSCARC_VECTOR_CELLTYPE, NOM, NL)
      
                  CASE (NSCARC_EXCHANGE_WEIGHTS)
IF (TYPE_DEBUG> NSCARC_DEBUG_NONE) WRITE(SCARC_LU,*) 'BEFORE UNPACK_CMATRIX'
                     CALL SCARC_UNPACK_CMATRIX(BUFFER, IJKW, NSCARC_MATRIX_PROLONGATION, NOM, NL)
IF (TYPE_DEBUG> NSCARC_DEBUG_NONE) WRITE(SCARC_LU,*) 'AFTER UNPACK_CMATRIX'
      
                END SELECT SELECT_EXCHANGE_COMPACT2

         END SELECT

      ENDIF EXCHANGE_RECV_IF
   ENDDO EXCHANGE_RECV_LOOP2
ENDDO EXCHANGE_SEND_LOOP2

END SUBROUTINE SCARC_SEND
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Pack send buffer for exchange of internal vector in banded storage technique
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_PACK_BVECTOR_REAL (SEND_BUF, IJKW, NTYPE, NLEN, NW, NM, NL)
REAL (EB), DIMENSION (:)  , INTENT(OUT) :: SEND_BUF
INTEGER  , DIMENSION (:,:), INTENT(IN)  :: IJKW
INTEGER, INTENT(IN)  :: NTYPE, NW, NM, NL
INTEGER, INTENT(OUT) :: NLEN
INTEGER ::  IW, IWW, LL, II, JJ, KK
REAL(EB), POINTER, DIMENSION(:,:,:) :: BVECTOR

BVECTOR => POINT_TO_BVECTOR (NTYPE, NM, NL)

LL  = 0
IWW = 0
PACK_SEND: DO IW=1,NW
   IF (IJKW(9,IW)/=NM) CYCLE PACK_SEND
   DO KK=IJKW(12,IW),IJKW(15,IW)
      DO JJ=IJKW(11,IW),IJKW(14,IW)
         DO II=IJKW(10,IW),IJKW(13,IW)
            IWW = IWW + 1
            SEND_BUF(LL+1) = REAL(IW,EB)
            SEND_BUF(LL+2) = BVECTOR(II,JJ,KK)
            LL = LL+2
         ENDDO
      ENDDO
   ENDDO
ENDDO PACK_SEND
NLEN=2*IWW+1

SEND_BUF(NLEN) = -999.0_EB

END SUBROUTINE SCARC_PACK_BVECTOR_REAL
 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Unpack receive buffer for internal update
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_UNPACK_BVECTOR_REAL (RECV_BUF, IJKW, NTYPE, NM, NL)
REAL (EB), DIMENSION (:)    , INTENT(IN)  :: RECV_BUF
INTEGER  , DIMENSION (:,:)  , INTENT(IN)  :: IJKW
INTEGER  , INTENT(IN):: NTYPE, NM, NL
REAL (EB):: ZSUM=0.0_EB
INTEGER IW, LL, ISUM, I, J, K, II, JJ, KK
REAL(EB), POINTER, DIMENSION(:,:,:) :: BVECTOR


BVECTOR => POINT_TO_BVECTOR (NTYPE, NM, NL)

LL = 0
UNPACK_RECV: DO

   IW = NINT(RECV_BUF(LL+1))
   IF (IW==-999) EXIT UNPACK_RECV
   ZSUM=0.0_EB
   DO KK=IJKW(12,IW),IJKW(15,IW)
      DO JJ=IJKW(11,IW),IJKW(14,IW)
         DO II=IJKW(10,IW),IJKW(13,IW)
            ZSUM=ZSUM+RECV_BUF(LL+2)
            LL = LL+2
         ENDDO
      ENDDO
   ENDDO

   ISUM = (IJKW(13,IW)-IJKW(10,IW)+1) * &
          (IJKW(14,IW)-IJKW(11,IW)+1) * &
          (IJKW(15,IW)-IJKW(12,IW)+1)

   I=IJKW(1,IW)
   J=IJKW(2,IW)
   K=IJKW(3,IW)
   BVECTOR(I, J, K) = ZSUM/REAL(ISUM,EB)

!IF (TYPE_DEBUG >= NSCARC_DEBUG_MEDIUM) WRITE(SCARC_LU,*) 'UNPACK_BVECTOR: X(',I,',',J,',',K,')=',BVECTOR(I,J,K)
ENDDO UNPACK_RECV

END SUBROUTINE SCARC_UNPACK_BVECTOR_REAL


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Pack send buffer for exchange of internal vector in compact storage technique
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_PACK_CVECTOR_REAL (SEND_BUF, IJKW, NTYPE, NLEN, NW, NM, NL)
REAL (EB), DIMENSION (:)  , INTENT(OUT) :: SEND_BUF
INTEGER  , DIMENSION (:,:), INTENT(IN)  :: IJKW
INTEGER, INTENT(IN)  :: NTYPE, NW, NM, NL
INTEGER, INTENT(OUT) :: NLEN
INTEGER ::  IW, IWW, LL, II, JJ, KK, IC, ICPL
REAL(EB), POINTER, DIMENSION(:) :: CVECTOR
TYPE (SCARC_COMPACT_TYPE), POINTER:: SC

SC => SCARC(NM)%COMPACT(NL)
CVECTOR => POINT_TO_CVECTOR (NTYPE, NM, NL)

IF (TYPE_MULTIGRID == NSCARC_MULTIGRID_ALGEBRAIC .AND. NL > 1) THEN

   LL  = 0
   IWW = 0
   PACK_SEND1: DO IW=1,NW
      IF (IJKW(9,IW)/=NM) CYCLE PACK_SEND1
      DO ICPL = 0, IJKW(2,IW)
         IWW = IWW + 1
         IC = IJKW(10+ICPL,IW)
         SEND_BUF(LL+1) = REAL(IW,EB)
         SEND_BUF(LL+2) = CVECTOR(IC)
   IF (TYPE_DEBUG > NSCARC_DEBUG_MEDIUM) WRITE(SCARC_LU,*) NM,': PACK_CVECTOR: SEND_BUF(',LL+2,')=',SEND_BUF(LL+2)
         LL = LL+2
      ENDDO
   ENDDO PACK_SEND1
   NLEN=2*IWW+1
   
   SEND_BUF(NLEN) = -999.0_EB

ELSE

   LL  = 0
   IWW = 0
   PACK_SEND2: DO IW=1,NW
      IF (IJKW(9,IW)/=NM) CYCLE PACK_SEND2
      DO KK=IJKW(12,IW),IJKW(15,IW)
         DO JJ=IJKW(11,IW),IJKW(14,IW)
            DO II=IJKW(10,IW),IJKW(13,IW)
               IWW = IWW + 1
               IC = (KK-1)*SC%NX*SC%NY + (JJ-1)*SC%NX + II
               SEND_BUF(LL+1) = REAL(IW,EB)
               SEND_BUF(LL+2) = CVECTOR(IC)
   IF (TYPE_DEBUG > NSCARC_DEBUG_MEDIUM) WRITE(SCARC_LU,*) NM,': PACK_CVECTOR: SEND_BUF(',LL+2,')=',SEND_BUF(LL+2)
               LL = LL+2
            ENDDO
         ENDDO
      ENDDO
   ENDDO PACK_SEND2
   NLEN=2*IWW+1
   
   SEND_BUF(NLEN) = -999.0_EB

ENDIF

END SUBROUTINE SCARC_PACK_CVECTOR_REAL
 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Unpack receive buffer for internal update
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_UNPACK_CVECTOR_REAL (RECV_BUF, IJKW, NTYPE, NM, NL)
REAL (EB), DIMENSION (:)  , INTENT(IN) :: RECV_BUF
INTEGER  , DIMENSION (:,:), INTENT(IN) :: IJKW
INTEGER  , INTENT(IN) :: NTYPE, NM, NL
REAL (EB):: ZSUM
INTEGER  :: IW, LL, ISUM, II, JJ, KK, IC, ICPL
REAL(EB), POINTER, DIMENSION(:) :: CVECTOR
TYPE (SCARC_COMPACT_TYPE), POINTER:: SC

SC => SCARC(NM)%COMPACT(NL)

CVECTOR => POINT_TO_CVECTOR (NTYPE, NM, NL)

IF (TYPE_MULTIGRID == NSCARC_MULTIGRID_ALGEBRAIC .AND. NL > 1) THEN

   LL = 0
   UNPACK_RECV1: DO
   
      IW = NINT(RECV_BUF(LL+1))
      IF (IW==-999) EXIT UNPACK_RECV1

      DO ICPL = 0, IJKW(2, IW)
         IC = IJKW(5+ICPL, IW)
         CVECTOR(IC) = RECV_BUF(LL+2)
   IF (TYPE_DEBUG > NSCARC_DEBUG_MEDIUM) WRITE(SCARC_LU,*) NM,': UNPACK_CVECTOR: X(',IC,')=',CVECTOR(IC), ICPL
         LL = LL + 2
      ENDDO
   
   ENDDO UNPACK_RECV1

ELSE

   LL = 0
   UNPACK_RECV2: DO
   
      IW = NINT(RECV_BUF(LL+1))
      IF (IW==-999) EXIT UNPACK_RECV2
      ZSUM=0.0_EB
      DO KK=IJKW(12,IW),IJKW(15,IW)
         DO JJ=IJKW(11,IW),IJKW(14,IW)
            DO II=IJKW(10,IW),IJKW(13,IW)
               ZSUM=ZSUM+RECV_BUF(LL+2)
               LL = LL+2
            ENDDO
         ENDDO
      ENDDO
   
      ISUM = (IJKW(13,IW)-IJKW(10,IW)+1) * &
             (IJKW(14,IW)-IJKW(11,IW)+1) * &
             (IJKW(15,IW)-IJKW(12,IW)+1)
   
      !IC = SC%ADJACENT_CELL(IW)
      IC = SC%GHOST_CELL(IW)
!WRITE(*,*) 'NM=',NM,' :NL=',NL,': IC=',IC, ' IW=',IW, SIZE(SC%GHOST_CELL)
      CVECTOR(IC) = ZSUM/REAL(ISUM,EB)
   
   IF (TYPE_DEBUG > NSCARC_DEBUG_MEDIUM) WRITE(SCARC_LU,*) NM,': UNPACK_CVECTOR: X(',IC,')=',CVECTOR(IC)
   
   ENDDO UNPACK_RECV2

ENDIF

END SUBROUTINE SCARC_UNPACK_CVECTOR_REAL


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Pack send buffer for exchange of internal vector in compact storage technique
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_PACK_CVECTOR_INT (SEND_BUF, IJKW, NTYPE, NLEN, NW, NM, NL)
REAL (EB), DIMENSION (:)  , INTENT(OUT) :: SEND_BUF
INTEGER  , DIMENSION (:,:), INTENT(IN)  :: IJKW
INTEGER  , INTENT(IN)  :: NTYPE, NW, NM, NL
INTEGER  , INTENT(OUT) :: NLEN
INTEGER  , POINTER, DIMENSION(:) :: CVECTOR
INTEGER ::  IW, IWW, LL, II, JJ, KK, IC
TYPE (SCARC_COMPACT_TYPE), POINTER:: SC

SC => SCARC(NM)%COMPACT(NL)
CVECTOR => POINT_TO_CVECTOR_INT (NTYPE, NM, NL)

LL  = 0
IWW = 0
PACK_SEND: DO IW=1,NW
   IF (IJKW(9,IW)/=NM) CYCLE PACK_SEND
   DO KK=IJKW(12,IW),IJKW(15,IW)
      DO JJ=IJKW(11,IW),IJKW(14,IW)
         DO II=IJKW(10,IW),IJKW(13,IW)
            IWW = IWW + 1
            IC = (KK-1)*SC%NX*SC%NY + (JJ-1)*SC%NX + II
            SEND_BUF(LL+1) = REAL(IW,EB)
            SEND_BUF(LL+2) = REAL(CVECTOR(IC),EB)
   IF (TYPE_DEBUG > NSCARC_DEBUG_MEDIUM) WRITE(SCARC_LU,*) NM,': PACK_CVECTOR: SEND_BUF(',LL+2,')=',SEND_BUF(LL+2)
            LL = LL+2
         ENDDO
      ENDDO
   ENDDO
ENDDO PACK_SEND
NLEN=2*IWW+1

SEND_BUF(NLEN) = -999.0_EB

END SUBROUTINE SCARC_PACK_CVECTOR_INT
 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Unpack receive buffer for internal update
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_UNPACK_CVECTOR_INT (RECV_BUF, IJKW, NTYPE, NM, NL)
REAL (EB), DIMENSION (:)  , INTENT(IN) :: RECV_BUF
INTEGER  , DIMENSION (:,:), INTENT(IN) :: IJKW
INTEGER  , INTENT(IN) :: NTYPE, NM, NL
REAL (EB):: ZSUM
INTEGER  :: IW, LL, ISUM, II, JJ, KK, IC
INTEGER,  POINTER, DIMENSION(:) :: CVECTOR
TYPE (SCARC_COMPACT_TYPE), POINTER:: SC

SC => SCARC(NM)%COMPACT(NL)

CVECTOR => POINT_TO_CVECTOR_INT (NTYPE, NM, NL)

LL = 0
UNPACK_RECV: DO

   IW = NINT(RECV_BUF(LL+1))
   IF (IW==-999) EXIT UNPACK_RECV
   ZSUM=0.0_EB
   DO KK=IJKW(12,IW),IJKW(15,IW)
      DO JJ=IJKW(11,IW),IJKW(14,IW)
         DO II=IJKW(10,IW),IJKW(13,IW)
            ZSUM=ZSUM+RECV_BUF(LL+2)
            LL = LL+2
         ENDDO
      ENDDO
   ENDDO

   ISUM = (IJKW(13,IW)-IJKW(10,IW)+1) * &
          (IJKW(14,IW)-IJKW(11,IW)+1) * &
          (IJKW(15,IW)-IJKW(12,IW)+1)

   IC = SC%GHOST_CELL(IW)
   CVECTOR(IC) = NINT(ZSUM/REAL(ISUM,EB))

ENDDO UNPACK_RECV

END SUBROUTINE SCARC_UNPACK_CVECTOR_INT


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Pack send buffer for exchange of matrix subdiagonals (for both system types)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_PACK_BMATRIX (SEND_BUF, IJKW, NTYPE, NLEN, NW, NM, NL)
REAL (EB), DIMENSION (:)  , INTENT(OUT) :: SEND_BUF
INTEGER  , DIMENSION (:,:), INTENT(IN)  :: IJKW
INTEGER  , INTENT(IN)  :: NTYPE, NW, NM, NL
INTEGER  , INTENT(OUT) :: NLEN
INTEGER ::  IW, IWW, LL, II, JJ, KK, IOR0
TYPE (SCARC_BANDED_TYPE), POINTER:: SB

SB => SCARC(NM)%BANDED(NL)

!IF (NM==1.AND.NL==3) THEN
!WRITE(SCARC_LU,*) '=================================='
!WRITE(SCARC_LU,*) 'PACKING BMATRIX on mesh ', NM, ' level ', NL
!WRITE(SCARC_LU,*) '=================================='
!WRITE(SCARC_LU,*) 'NW(1)=',SCARC(1)%BANDED(3)%NW
!WRITE(SCARC_LU,*) 'NW(2)=',SCARC(2)%BANDED(3)%NW
!WRITE(SCARC_LU,*) 'NW(3)=',SCARC(3)%BANDED(3)%NW
!WRITE(SCARC_LU,*) 'NW(4)=',SCARC(4)%BANDED(3)%NW
!WRITE(SCARC_LU,*) 'NW=',NW
!WRITE(SCARC_LU,*) 'IJKW='
!DO IW=1,NW
!WRITE(SCARC_LU,'(15i3)') (IJKW(II,IW),II=1,15)
!ENDDO
!ENDIF

LL = 0
IWW = 0
SELECT CASE (NTYPE)

   !!! -------------------------------------------------------------------------------------------
   !!! System matrix
   !!! -------------------------------------------------------------------------------------------
   CASE (NSCARC_MATRIX_SUBDIAG)

      PACK_SEND_SYSTEM: DO IW=1,NW
         IF (IJKW(9,IW)/=NM) CYCLE PACK_SEND_SYSTEM
         IOR0 = ABS(IJKW(4,IW))
         DO KK=IJKW(12,IW),IJKW(15,IW)
            DO JJ=IJKW(11,IW),IJKW(14,IW)
               DO II=IJKW(10,IW),IJKW(13,IW)
                  IWW = IWW + 1
                  SEND_BUF(LL+1) = REAL(IW,EB)
                  SEND_BUF(LL+2) = SB%DI2(IOR0)
!IF (NM==1.AND.NL==3) WRITE(SCARC_LU,*) 'SENDING  ', II, KK, IOR0, SB%DI2(IOR0)
                  LL = LL+2
               ENDDO
            ENDDO
         ENDDO
      ENDDO PACK_SEND_SYSTEM

   CASE DEFAULT
   
      WRITE(*,*) 'NTYPE=',NTYPE,' NOT ALLOWED FOR SCARC_PACK_BMATRIX'
      STOP

END SELECT

NLEN=2*IWW+1
SEND_BUF(NLEN) = -999.0_EB
      
END SUBROUTINE SCARC_PACK_BMATRIX
 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Unpack receive buffer for matrix-vector multiplication for compact storage technique
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_UNPACK_BMATRIX (RECV_BUF, IJKW, NTYPE, NM, NL)
REAL (EB), DIMENSION (:)  , INTENT(IN)  :: RECV_BUF
INTEGER  , DIMENSION (:,:), INTENT(IN)  :: IJKW
INTEGER  , INTENT(IN) :: NTYPE, NM, NL
REAL (EB):: ZSUM
INTEGER IW, LL, ISUM, I, J, K, II, JJ, KK, IC, ICPL
TYPE (SCARC_BANDED_TYPE), POINTER:: SB

SB => SCARC(NM)%BANDED(NL)

!IF (NM==1.AND.NL==3) THEN
!WRITE(SCARC_LU,*) '=================================='
!WRITE(SCARC_LU,*) 'UNPACKING BMATRIX on mesh ', NM, ' level ', NL
!WRITE(SCARC_LU,*) '=================================='
!ENDIF

LL = 0
SELECT CASE (NTYPE)

   !!! -------------------------------------------------------------------------------------------
   !!! System matrix
   !!! -------------------------------------------------------------------------------------------
   CASE (NSCARC_MATRIX_SUBDIAG)
      UNPACK_RECV: DO
      
         IW = NINT(RECV_BUF(LL+1))
         IF (IW==-999) EXIT UNPACK_RECV
         ZSUM=0.0_EB
         DO KK=IJKW(12,IW),IJKW(15,IW)
            DO JJ=IJKW(11,IW),IJKW(14,IW)
               DO II=IJKW(10,IW),IJKW(13,IW)
                  ZSUM=ZSUM+RECV_BUF(LL+2)
                  LL = LL+2
               ENDDO
            ENDDO
         ENDDO
      
         ISUM = (IJKW(13,IW)-IJKW(10,IW)+1) * &
                (IJKW(14,IW)-IJKW(11,IW)+1) * &
                (IJKW(15,IW)-IJKW(12,IW)+1)
      
         I=IJKW(6,IW)
         J=IJKW(7,IW)
         K=IJKW(8,IW)
      
         IC = (K-1)*SB%NX*SB%NY + (J-1)*SB%NX + I
         DO ICPL = 2, SB%NCPL
            IF (IW == -SB%A(IC,ICPL)) THEN
               SB%A(IC,ICPL) = ZSUM/REAL(ISUM,EB)
               SB%ADJACENT_CELL (IW) = IC
!IF (NM==1.AND.NL==3) &
!   WRITE(SCARC_LU,'(A,i3,a,i3,a,e14.6,a,i3,a,i3)') 'GETTING  A(',IC,',',ICPL,')=',SB%A(IC,ICPL),': IW=',IW,': ADJACENT=',IC
            ENDIF
         ENDDO
      
      ENDDO UNPACK_RECV

   CASE DEFAULT

      WRITE(*,*) 'NTYPE =',NTYPE,' NOT ALLOWED FOR SCARC_UNPACK_BMATRIX'
      STOP

END SELECT
      
END SUBROUTINE SCARC_UNPACK_BMATRIX


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Pack send buffer for exchange of matrix subdiagonals (for both system types)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_PACK_CMATRIX (SEND_BUF, IJKW, NTYPE, NLEN, NW, NM, NL)
REAL (EB), DIMENSION (:)  , INTENT(OUT) :: SEND_BUF
INTEGER  , DIMENSION (:,:), INTENT(IN)  :: IJKW
INTEGER  , INTENT(IN)  :: NTYPE, NW, NM, NL
INTEGER  , INTENT(OUT) :: NLEN
INTEGER ::  IC, JC, ICOL, ICOL0, IW, IWW, LL, II, JJ, KK, IOR0, NCOL, NNUM
LOGICAL :: BLOWER_NBR, BUPPER_NBR
TYPE (SCARC_COMPACT_TYPE), POINTER:: SC, SCC

SC  => SCARC(NM)%COMPACT(NL)
SCC => SCARC(NM)%COMPACT(NL+1)

SELECT CASE (NTYPE)

   !!! -------------------------------------------------------------------------------------------
   !!! subdiagonal entries of System matrix (from adjacent cells)
   !!! -------------------------------------------------------------------------------------------
   CASE (NSCARC_MATRIX_SUBDIAG)

      LL = 0
      IWW = 0

      PACK_SEND_SUBDIAG: DO IW=1,NW
         IF (IJKW(9,IW)/=NM) CYCLE PACK_SEND_SUBDIAG
         IOR0 = ABS(IJKW(4,IW))
         DO KK=IJKW(12,IW),IJKW(15,IW)
            DO JJ=IJKW(11,IW),IJKW(14,IW)
               DO II=IJKW(10,IW),IJKW(13,IW)
                  IWW = IWW + 1
                  SEND_BUF(LL+1) = REAL(IW,EB)
                  SEND_BUF(LL+2) = SC%DI2(IOR0)
                  LL = LL+2
               ENDDO
            ENDDO
         ENDDO
      ENDDO PACK_SEND_SUBDIAG

      NLEN=2*IWW+1

   !!! -------------------------------------------------------------------------------------------
   !!! subdiagonal entries of System matrix (from adjacent cells)
   !!! -------------------------------------------------------------------------------------------
   CASE (NSCARC_MATRIX_SYSTEM1)

      LL = 1
      SEND_BUF = 0.0_EB

   IF (TYPE_DEBUG > NSCARC_DEBUG_MEDIUM)  WRITE(SCARC_LU,'(a)') 'A_ROW::'
   IF (TYPE_DEBUG > NSCARC_DEBUG_MEDIUM)  WRITE(SCARC_LU,'(10i3)') SC%A_ROW(1:30)
   IF (TYPE_DEBUG > NSCARC_DEBUG_MEDIUM)  WRITE(SCARC_LU,'(a)') 'A_COL::'
   IF (TYPE_DEBUG > NSCARC_DEBUG_MEDIUM)  WRITE(SCARC_LU,'(10i3)') SC%A_COL(1:30)
   IF (TYPE_DEBUG > NSCARC_DEBUG_MEDIUM)  WRITE(SCARC_LU,'(a)') 'A::'
   IF (TYPE_DEBUG > NSCARC_DEBUG_MEDIUM)  WRITE(SCARC_LU,'(10f8.2)') SC%A(1:30)
      PACK_SEND_SYSTEM: DO IW=1,NW
         IF (NL == 1) THEN
            IF (IJKW(9,IW)/=NM) CYCLE PACK_SEND_SYSTEM
   
            DO KK=IJKW(12,IW),IJKW(15,IW)
               DO JJ=IJKW(11,IW),IJKW(14,IW)
                  DO II=IJKW(10,IW),IJKW(13,IW)
   
                     IC = (KK-1)*SC%NX*SC%NY + (JJ-1)*SC%NX + II
      IF (TYPE_DEBUG > NSCARC_DEBUG_MEDIUM) write(scarc_lu,'(a,5i4)') '=================== IC=',IC, IW, II, JJ, KK

                     SEND_BUF(LL) = REAL(IW,EB)

      IF (TYPE_DEBUG > NSCARC_DEBUG_MEDIUM) &
         WRITE(SCARC_LU,'(2i3,a,i3,a,f12.6,i3)') NM,iW,': PACK_WEIGHTS: SEND_BUF(',LL,')=',SEND_BUF(LL), IC
                     LL = LL + 1

                     NCOL = SC%A_ROW(IC+1)-SC%A_ROW(IC)
                     SEND_BUF(LL) = NCOL

      IF (TYPE_DEBUG > NSCARC_DEBUG_MEDIUM) &
         WRITE(SCARC_LU,'(2i3,a,i3,a,f12.6,2i3)') NM,iW,': PACK_WEIGHTS: SEND_BUF(',LL,')=',SEND_BUF(LL), IC, NCOL
                     LL = LL + 1

                     COLUMN_SYSTEM_LOOP: DO ICOL = SC%A_ROW(IC), SC%A_ROW(IC+1)-1
                        JC = SC%A_COL(ICOL)
                        IF (JC >= SC%NC) THEN
                           SEND_BUF(LL) = - REAL(JC,EB)
                        ELSE
                           SEND_BUF(LL) =   REAL(JC,EB)
                        ENDIF
                        SEND_BUF(LL+1) = SC%A(ICOL)
      IF (TYPE_DEBUG > NSCARC_DEBUG_MEDIUM) &
         WRITE(SCARC_LU,'(2i3,a,i3,a,f12.6,i3)') NM,iW,': PACK_WEIGHTS: SEND_BUF(',LL  ,')=',SEND_BUF(LL), IC
      IF (TYPE_DEBUG > NSCARC_DEBUG_MEDIUM) &
         WRITE(SCARC_LU,'(2i3,a,i3,a,f12.6,i3)') NM,iW,': PACK_WEIGHTS: SEND_BUF(',LL+1,')=',SEND_BUF(LL+1), IC
                        LL = LL + 2
                     ENDDO COLUMN_SYSTEM_LOOP

                  ENDDO
               ENDDO
            ENDDO
         ELSE
            WRITE(*,*) 'STILL MISSING1'
         ENDIF
      ENDDO PACK_SEND_SYSTEM
   IF (TYPE_DEBUG > NSCARC_DEBUG_MEDIUM)  WRITE(SCARC_LU,'(a)') '2:A_ROW::'
   IF (TYPE_DEBUG > NSCARC_DEBUG_MEDIUM)  WRITE(SCARC_LU,'(10i3)') SC%A_ROW(1:30)
   IF (TYPE_DEBUG > NSCARC_DEBUG_MEDIUM)  WRITE(SCARC_LU,'(a)') '2:A_COL::'
   IF (TYPE_DEBUG > NSCARC_DEBUG_MEDIUM)  WRITE(SCARC_LU,'(10i3)') SC%A_COL(1:30)
   IF (TYPE_DEBUG > NSCARC_DEBUG_MEDIUM)  WRITE(SCARC_LU,'(a)') '2:A::'
   IF (TYPE_DEBUG > NSCARC_DEBUG_MEDIUM)  WRITE(SCARC_LU,'(10f8.2)') SC%A(1:30)
   IF (TYPE_DEBUG > NSCARC_DEBUG_MEDIUM)  WRITE(SCARC_LU,'(a)') 'SEND_BUF::'
   IF (TYPE_DEBUG > NSCARC_DEBUG_MEDIUM)  WRITE(SCARC_LU,'(10f8.3)') SEND_BUF(1:70)

      NLEN=LL

   !!! -------------------------------------------------------------------------------------------
   !!! subdiagonal entries of System matrix (from adjacent cells)
   !!! -------------------------------------------------------------------------------------------
   CASE (NSCARC_MATRIX_SYSTEM2)

      LL = 1
      SEND_BUF = 0.0_EB

   IF (TYPE_DEBUG > NSCARC_DEBUG_MEDIUM)  WRITE(SCARC_LU,'(a)') 'A_ROW::'
   IF (TYPE_DEBUG > NSCARC_DEBUG_MEDIUM)  WRITE(SCARC_LU,'(10i3)') SC%A_ROW(1:30)
   IF (TYPE_DEBUG > NSCARC_DEBUG_MEDIUM)  WRITE(SCARC_LU,'(a)') 'A_COL::'
   IF (TYPE_DEBUG > NSCARC_DEBUG_MEDIUM)  WRITE(SCARC_LU,'(10i3)') SC%A_COL(1:30)
   IF (TYPE_DEBUG > NSCARC_DEBUG_MEDIUM)  WRITE(SCARC_LU,'(a)') 'A::'
   IF (TYPE_DEBUG > NSCARC_DEBUG_MEDIUM)  WRITE(SCARC_LU,'(10f8.2)') SC%A(1:30)
      PACK_SEND_SYSTEM2: DO IW=1,NW
         IF (NL == 1) THEN
            IF (IJKW(9,IW)/=NM) CYCLE PACK_SEND_SYSTEM2
   
            DO KK=IJKW(12,IW),IJKW(15,IW)
               DO JJ=IJKW(11,IW),IJKW(14,IW)
                  DO II=IJKW(10,IW),IJKW(13,IW)
   
                     IC = (KK-1)*SC%NX*SC%NY + (JJ-1)*SC%NX + II
      IF (TYPE_DEBUG > NSCARC_DEBUG_MEDIUM) write(scarc_lu,'(a,5i4)') '=================== IC=',IC, IW, II, JJ, KK

                     SEND_BUF(LL) = REAL(IW,EB)

      IF (TYPE_DEBUG > NSCARC_DEBUG_MEDIUM) &
         WRITE(SCARC_LU,'(2i3,a,i3,a,f12.6,i3)') NM,iW,': PACK_WEIGHTS: SEND_BUF(',LL,')=',SEND_BUF(LL), IC
                     LL = LL + 1

                     NCOL = SC%A_ROW(IC+1)-SC%A_ROW(IC)
                     SEND_BUF(LL) = NCOL

      IF (TYPE_DEBUG > NSCARC_DEBUG_MEDIUM) &
         WRITE(SCARC_LU,'(2i3,a,i3,a,f12.6,2i3)') NM,iW,': PACK_WEIGHTS: SEND_BUF(',LL,')=',SEND_BUF(LL), IC, NCOL
                     LL = LL + 1

                     COLUMN_SYSTEM2_LOOP: DO ICOL = SC%A_ROW(IC), SC%A_ROW(IC+1)-1
                        JC = SC%A_COL(ICOL)
                        IF (JC >= SC%NC) THEN
                           SEND_BUF(LL) = - REAL(JC,EB)
                        ELSE
                           SEND_BUF(LL) =   REAL(JC,EB)
                        ENDIF
                        SEND_BUF(LL+1) = SC%A(ICOL)
      IF (TYPE_DEBUG > NSCARC_DEBUG_MEDIUM) &
         WRITE(SCARC_LU,'(2i3,a,i3,a,f12.6,i3)') NM,iW,': PACK_WEIGHTS: SEND_BUF(',LL  ,')=',SEND_BUF(LL), IC
      IF (TYPE_DEBUG > NSCARC_DEBUG_MEDIUM) &
         WRITE(SCARC_LU,'(2i3,a,i3,a,f12.6,i3)') NM,iW,': PACK_WEIGHTS: SEND_BUF(',LL+1,')=',SEND_BUF(LL+1), IC
                        LL = LL + 2
                     ENDDO COLUMN_SYSTEM2_LOOP

                  ENDDO
               ENDDO
            ENDDO
         ELSE
            WRITE(*,*) 'STILL MISSING2'
         ENDIF
      ENDDO PACK_SEND_SYSTEM2
   IF (TYPE_DEBUG > NSCARC_DEBUG_MEDIUM)  WRITE(SCARC_LU,'(a)') '2:A_ROW::'
   IF (TYPE_DEBUG > NSCARC_DEBUG_MEDIUM)  WRITE(SCARC_LU,'(10i3)') SC%A_ROW(1:30)
   IF (TYPE_DEBUG > NSCARC_DEBUG_MEDIUM)  WRITE(SCARC_LU,'(a)') '2:A_COL::'
   IF (TYPE_DEBUG > NSCARC_DEBUG_MEDIUM)  WRITE(SCARC_LU,'(10i3)') SC%A_COL(1:30)
   IF (TYPE_DEBUG > NSCARC_DEBUG_MEDIUM)  WRITE(SCARC_LU,'(a)') '2:A::'
   IF (TYPE_DEBUG > NSCARC_DEBUG_MEDIUM)  WRITE(SCARC_LU,'(10f8.2)') SC%A(1:30)
   IF (TYPE_DEBUG > NSCARC_DEBUG_MEDIUM)  WRITE(SCARC_LU,'(a)') 'SEND_BUF::'
   IF (TYPE_DEBUG > NSCARC_DEBUG_MEDIUM)  WRITE(SCARC_LU,'(10f8.3)') SEND_BUF(1:70)

      NLEN=LL


   !!! -------------------------------------------------------------------------------------------
   !!! Prolongation matrix
   !!! -------------------------------------------------------------------------------------------
   CASE (NSCARC_MATRIX_PROLONGATION)

      IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) write(scarc_lu,*) 'SEND MATRIX_PROLONATION'
      IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) write(scarc_lu,'(a,5i4)') '============PACK_CMATRIX=============== NM=',NM, SC%NCC
      LL = 1
      IF (NL == 1) THEN

         PACK_SEND_PROL: DO IW=1,NW

      IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) write(scarc_lu,*) 'IJKW(9,',IW,')=',IJKW(9,IW),NM
            IF (IJKW(9,IW)/=NM) CYCLE PACK_SEND_PROL
   
            !!! -------------------------------------------------------------------------------
            !!! For each cell adajacent to a ghost cell store the following information:
            !!! - number of cells in the stencil
            !!! - corresponding IW (positive, if coarse cell; negative, if fine cell)
            !!! - corresponding fine cell number (number itself, if coarse cell; numbers in stencil, if fine cell)
            !!! - corresponding coarse cell number
            !!! - corresponding prolongation weight
            !!! -------------------------------------------------------------------------------
            DO KK=IJKW(12,IW),IJKW(15,IW)
               DO JJ=IJKW(11,IW),IJKW(14,IW)
                  DO II=IJKW(10,IW),IJKW(13,IW)
   
                     IC = (KK-1)*SC%NX*SC%NY + (JJ-1)*SC%NX + II
      IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) write(scarc_lu,'(a,8i4)') '=================== IC=',IC, IW, II, JJ, KK,SIZE(SC%P_ROW)

                     IF (SC%CELLTYPE(IC) >= NSCARC_CELLTYPE_COARSE) THEN
                        NNUM = 1
                        ICOL = SC%P_ROW(IC)
                        SEND_BUF(LL)   = REAL(NNUM,EB)
                        SEND_BUF(LL+1) = REAL(NM,EB)
                        SEND_BUF(LL+2) = REAL(IC,EB)
                        SEND_BUF(LL+3) = REAL(IW,EB)
                        SEND_BUF(LL+4) = REAL(IC,EB)
                        SEND_BUF(LL+5) = REAL(SC%CELLTYPE(IC),EB)
                        SEND_BUF(LL+6) = SC%P(ICOL)
      IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) THEN
         WRITE(SCARC_LU,*) '-------------SENDING ------------------- COARSE'
         WRITE(SCARC_LU,'(2i3,a,i3,a,f12.6,i3)') NM,iW,': PWA: SEND_BUF(',LL  ,')=',SEND_BUF(LL), IC
         WRITE(SCARC_LU,'(2i3,a,i3,a,f12.6,i3)') NM,iW,': PWA: SEND_BUF(',LL+1,')=',SEND_BUF(LL+1), IC
         WRITE(SCARC_LU,'(2i3,a,i3,a,f12.6,i3)') NM,iW,': PWA: SEND_BUF(',LL+2,')=',SEND_BUF(LL+2), IC
         WRITE(SCARC_LU,'(2i3,a,i3,a,f12.6,i3)') NM,iW,': PWA: SEND_BUF(',LL+3,')=',SEND_BUF(LL+3), IC
         WRITE(SCARC_LU,*) '-----------'
         WRITE(SCARC_LU,'(2i3,a,i3,a,f12.6,i3)') NM,iW,': PWA: SEND_BUF(',LL+4,')=',SEND_BUF(LL+4), IC
         WRITE(SCARC_LU,'(2i3,a,i3,a,f12.6,i3)') NM,iW,': PWA: SEND_BUF(',LL+5,')=',SEND_BUF(LL+5), IC
         WRITE(SCARC_LU,'(2i3,a,i3,a,f12.6,i3)') NM,iW,': PWA: SEND_BUF(',LL+6,')=',SEND_BUF(LL+6), IC
      ENDIF
                        LL = LL + 7
                     ELSE
                        NNUM  = SC%A_ROW(IC+1)-SC%A_ROW(IC)-1
                        ICOL0 = SC%A_ROW(IC)
                        SEND_BUF(LL)   =   REAL(NNUM,EB)
                        SEND_BUF(LL+1) =   REAL(NM,EB)
                        SEND_BUF(LL+2) =   REAL(IC,EB)
                        SEND_BUF(LL+3) = - REAL(IW,EB)
      IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) THEN
         WRITE(SCARC_LU,*) '---------------------SENDING ----------- FINE'
         WRITE(SCARC_LU,'(2i3,a,i3,a,f12.6,i3)') NM,iW,': PWB: SEND_BUF(',LL  ,')=',SEND_BUF(LL), IC
         WRITE(SCARC_LU,'(2i3,a,i3,a,f12.6,i3)') NM,iW,': PWB: SEND_BUF(',LL+1,')=',SEND_BUF(LL+1), IC
         WRITE(SCARC_LU,'(2i3,a,i3,a,f12.6,i3)') NM,iW,': PWB: SEND_BUF(',LL+2,')=',SEND_BUF(LL+2), IC
         WRITE(SCARC_LU,'(2i3,a,i3,a,f12.6,i3)') NM,iW,': PWB: SEND_BUF(',LL+3,')=',SEND_BUF(LL+3), IC
      ENDIF
                        LL = LL + 4
                        COLUMN_PROL_LOOP: DO ICOL = SC%P_ROW(IC), SC%P_ROW(IC+1)-1
                           ICOL0 = ICOL0 + 1
                           JC = SC%A_COL(ICOL0)
                           IF (JC <= SC%NC) THEN
                              SEND_BUF(LL)   = REAL(JC,EB)
                              SEND_BUF(LL+1) = REAL(SC%CELLTYPE(JC),EB)
                           ELSE
                              SEND_BUF(LL)   = - REAL(JC,EB)
                              SEND_BUF(LL+1) = -1.0_EB
                           ENDIF
                           SEND_BUF(LL+2) = SC%P(ICOL)

      IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) THEN
         WRITE(SCARC_LU,*) '-----------'
         WRITE(SCARC_LU,'(2i3,a,i3,a,f12.6,i3)') NM,iW,': PWB: SEND_BUF(',LL  ,')=',SEND_BUF(LL), IC
         WRITE(SCARC_LU,'(2i3,a,i3,a,f12.6,i3)') NM,iW,': PWB: SEND_BUF(',LL+1,')=',SEND_BUF(LL+1), IC
         WRITE(SCARC_LU,'(2i3,a,i3,a,f12.6,i3)') NM,iW,': PWB: SEND_BUF(',LL+2,')=',SEND_BUF(LL+2), IC
      ENDIF
                           LL = LL + 3
                        ENDDO COLUMN_PROL_LOOP

                        DO ICOL = SC%A_ROW(IC)+1, SC%A_ROW(IC+1)-1
                           JC = SC%A_COL(ICOL)
WRITE(SCARC_LU,*) 'ICOL=',ICOL,': JC=',JC,': SLN=',SECOND_LINE_NBR2(IW,IJKW,SC%NX,SC%NY)
                           IF (JC == SECOND_LINE_NBR2(IW,IJKW,SC%NX,SC%NY)) THEN
                              SEND_BUF(LL)   = -1.0_EB
                              SEND_BUF(LL+1) =  REAL(NM,EB)
                              SEND_BUF(LL+2) =  REAL(JC,EB)
                              SEND_BUF(LL+3) = -REAL(IW,EB)
                              SEND_BUF(LL+4) =  REAL(JC,EB)
                              SEND_BUF(LL+5) =  REAL(SC%CELLTYPE(JC),EB)
                              SEND_BUF(LL+6) =  SC%P(SC%P_ROW(JC))
      IF (TYPE_DEBUG > NSCARC_DEBUG_NONE) THEN
         WRITE(SCARC_LU,*) '------------------SENDING -------------- 2. line: JC=', JC
         WRITE(SCARC_LU,*) 'JC=',JC,' SECOND_LINE_NBR(',IW,')=',SECOND_LINE_NBR2(IW, IJKW, SC%NX, SC%NY)
         WRITE(SCARC_LU,'(2i3,a,i3,a,f12.6,i3)') NM,iW,': PWC: SEND_BUF(',LL  ,')=',SEND_BUF(LL), IC
         WRITE(SCARC_LU,'(2i3,a,i3,a,f12.6,i3)') NM,iW,': PWC: SEND_BUF(',LL+1,')=',SEND_BUF(LL+1), IC
         WRITE(SCARC_LU,'(2i3,a,i3,a,f12.6,i3)') NM,iW,': PWC: SEND_BUF(',LL+2,')=',SEND_BUF(LL+2), IC
         WRITE(SCARC_LU,'(2i3,a,i3,a,f12.6,i3)') NM,iW,': PWC: SEND_BUF(',LL+3,')=',SEND_BUF(LL+3), IC
         WRITE(SCARC_LU,*) '-----------'
         WRITE(SCARC_LU,'(2i3,a,i3,a,f12.6,i3)') NM,iW,': PWC: SEND_BUF(',LL+4,')=',SEND_BUF(LL+4), IC
         WRITE(SCARC_LU,'(2i3,a,i3,a,f12.6,i3)') NM,iW,': PWC: SEND_BUF(',LL+5,')=',SEND_BUF(LL+5), IC
         WRITE(SCARC_LU,'(2i3,a,i3,a,f12.6,i3)') NM,iW,': PWC: SEND_BUF(',LL+6,')=',SEND_BUF(LL+6), IC
      ENDIF
                              LL = LL + 7
                           ENDIF
                        ENDDO

                     ENDIF
                  ENDDO
               ENDDO
            ENDDO

         ENDDO PACK_SEND_PROL
   
      NLEN=LL

      ENDIF

END SELECT

SEND_BUF(NLEN) = -999.0_EB
      
END SUBROUTINE SCARC_PACK_CMATRIX
 


INTEGER FUNCTION SECOND_LINE_NBR1(IW, IJKW, NX, NY)
INTEGER, INTENT(IN) :: IW, NX, NY
INTEGER  , DIMENSION (:,:), INTENT(IN)    :: IJKW
INTEGER:: JC

IF (TYPE_DIMENSION == NSCARC_DIMENSION_TWO) THEN
   SELECT CASE(IJKW(4,IW))
      CASE ( 1)
         JC = (IJKW(8,IW)-1)*NX + IJKW(6,IW)+1
      CASE (-1)
         JC = (IJKW(8,IW)-1)*NX + IJKW(6,IW)-1
      CASE ( 3)
         JC = (IJKW(8,IW)  )*NX + IJKW(6,IW)
      CASE (-3)
         JC = (IJKW(8,IW)-2)*NX + IJKW(6,IW)
   END SELECT
ELSE
   SELECT CASE(IJKW(4,IW))
      CASE ( 1)
         JC = (IJKW(8,IW)-1)*NX*NY + (IJKW(7,IW)-1)*NX + IJKW(6,IW)+1
      CASE (-1)
         JC = (IJKW(8,IW)-1)*NX*NY + (IJKW(7,IW)-1)*NX + IJKW(6,IW)-1
      CASE ( 2)
         JC = (IJKW(8,IW)-1)*NX*NY + (IJKW(7,IW)  )*NX + IJKW(6,IW)
      CASE (-2)
         JC = (IJKW(8,IW)-1)*NX*NY + (IJKW(7,IW)-2)*NX + IJKW(6,IW)
      CASE ( 3)
         JC = (IJKW(8,IW)  )*NX*NY + (IJKW(7,IW)-1)*NX + IJKW(6,IW)
      CASE (-3)
         JC = (IJKW(8,IW)-2)*NX*NY + (IJKW(7,IW)-1)*NX + IJKW(6,IW)
   END SELECT
ENDIF

SECOND_LINE_NBR1 = JC
RETURN

END FUNCTION SECOND_LINE_NBR1


INTEGER FUNCTION SECOND_LINE_NBR2(IW, IJKW, NX, NY)
INTEGER, INTENT(IN) :: IW, NX, NY
INTEGER  , DIMENSION (:,:), INTENT(IN)    :: IJKW
INTEGER:: JC

IF (TYPE_DIMENSION == NSCARC_DIMENSION_TWO) THEN
   SELECT CASE(IJKW(4,IW))
      CASE ( 1)
         JC = (IJKW(12,IW)-1)*NX + IJKW(10,IW)-1
      CASE (-1)
         JC = (IJKW(12,IW)-1)*NX + IJKW(10,IW)+1
      CASE ( 3)
         JC = (IJKW(12,IW)-2)*NX + IJKW(10,IW)
      CASE (-3)
         JC = (IJKW(12,IW)  )*NX + IJKW(10,IW)
   END SELECT
ELSE
   SELECT CASE(IJKW(4,IW))
      CASE ( 1)
         JC = (IJKW(12,IW)-1)*NX*NY + (IJKW(11,IW)-1)*NX + IJKW(10,IW)-1
      CASE (-1)
         JC = (IJKW(12,IW)-1)*NX*NY + (IJKW(11,IW)-1)*NX + IJKW(10,IW)+1
      CASE ( 2)
         JC = (IJKW(12,IW)-1)*NX*NY + (IJKW(11,IW)-2)*NX + IJKW(10,IW)
      CASE (-2)
         JC = (IJKW(12,IW)-1)*NX*NY + (IJKW(11,IW)  )*NX + IJKW(10,IW)
      CASE ( 3)
         JC = (IJKW(12,IW)-2)*NX*NY + (IJKW(11,IW)-1)*NX + IJKW(10,IW)
      CASE (-3)
         JC = (IJKW(12,IW)  )*NX*NY + (IJKW(11,IW)-1)*NX + IJKW(10,IW)
   END SELECT
ENDIF

SECOND_LINE_NBR2 = JC
RETURN

END FUNCTION SECOND_LINE_NBR2


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Unpack receive buffer for matrix-vector multiplication for compact storage technique
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_UNPACK_CMATRIX (RECV_BUF, IJKW, NTYPE, NM, NL)
REAL (EB), DIMENSION (:)  , INTENT(IN)    :: RECV_BUF
INTEGER  , DIMENSION (:,:), INTENT(IN)    :: IJKW
INTEGER  , INTENT(IN)    :: NTYPE, NM, NL
REAL (EB):: ZSUM, PW
INTEGER IW, LL, ISUM, I, J, K, II, JJ, KK, IC, IC0, JC, IROW, ICOL, ICOL0, JCOL, NCOL, NROW, INUM, NNUM
INTEGER ICG, ICF, ICC, ICA, ICE, ICE2, NOM
TYPE (SCARC_COMPACT_TYPE), POINTER:: SC

SC => SCARC(NM)%COMPACT(NL)

SELECT CASE (NTYPE)

   !!! -------------------------------------------------------------------------------------------
   !!! Subdiagonal matrix entries along adjacent cells
   !!! -------------------------------------------------------------------------------------------
   CASE (NSCARC_MATRIX_SUBDIAG)

      LL = 0
      UNPACK_RECV_SUBDIAG: DO
      
         IW = NINT(RECV_BUF(LL+1))
         IF (IW==-999) EXIT UNPACK_RECV_SUBDIAG
         ZSUM=0.0_EB
         DO KK=IJKW(12,IW),IJKW(15,IW)
            DO JJ=IJKW(11,IW),IJKW(14,IW)
               DO II=IJKW(10,IW),IJKW(13,IW)
                  ZSUM=ZSUM+RECV_BUF(LL+2)
                  LL = LL+2
               ENDDO
            ENDDO
         ENDDO
      
         ISUM = (IJKW(13,IW)-IJKW(10,IW)+1) * &
                (IJKW(14,IW)-IJKW(11,IW)+1) * &
                (IJKW(15,IW)-IJKW(12,IW)+1)
      
         I=IJKW(6,IW)
         J=IJKW(7,IW)
         K=IJKW(8,IW)
      
         IC   = (K-1)*SC%NX*SC%NY + (J-1)*SC%NX + I
         IROW = SC%A_ROW(IC)
         DO ICOL = SC%A_ROW(IC)+1, SC%A_ROW(IC+1)-1

   IF (TYPE_DEBUG >= NSCARC_DEBUG_LESS.AND.NM==2) THEN
      WRITE(SCARC_LU,'(a,i6,a,7i6)') 'IW=',IW, ': ICOL=', ICOL,I, J, K, IC, IROW, SC%A_COL(ICOL)
   ENDIF
            IF (IW == -SC%A_COL(ICOL)) THEN
               SC%A(ICOL)     = ZSUM/REAL(ISUM,EB)
               SC%NCE         = SC%NCE + 1
               SC%A_COL(ICOL) = SC%NCE
               SC%ADJACENT_CELL (IW) = IC
               SC%GHOST_CELL (IW)    = SC%NCE
   IF (TYPE_DEBUG >= NSCARC_DEBUG_LESS.AND.NM==2) THEN
      WRITE(SCARC_LU,'(a,i2,a,i6)') '==================UNPACK_CMATRIX FOR NM=',NM,':  IC=',IC
      WRITE(SCARC_LU,*) I, J, K, SC%NX, SC%NY, SC%NZ
      WRITE(SCARC_LU,*) 'SC%A(',ICOL,')=',SC%A(ICOL)
      WRITE(SCARC_LU,*) 'SC%A_ICOL(',ICOL,')=',SC%A_COL(ICOL)
      WRITE(SCARC_LU,*) 'SC%ADJACENT_CELL(',IW,')=',SC%ADJACENT_CELL(IW)
      WRITE(SCARC_LU,*) 'SC%GHOST_CELL(',IW,')=',SC%GHOST_CELL(IW)
   ENDIF
            ENDIF
         ENDDO
      
      ENDDO UNPACK_RECV_SUBDIAG
   
   !!! -------------------------------------------------------------------------------------------
   !!! System matrix
   !!! -------------------------------------------------------------------------------------------
   CASE (NSCARC_MATRIX_SYSTEM1)

SC => SCARC(NM)%COMPACT(NL)
      IF (TYPE_DEBUG >= NSCARC_DEBUG_MEDIUM) write(scarc_lu,'(a,5i4)') '===================================== NM=',NM, SC%NC, NL
      LL = 1
   IF (TYPE_DEBUG > NSCARC_DEBUG_MEDIUM)  WRITE(SCARC_LU,*) 'UNPACK_WEIGHTS: NROW=',NROW
   IF (TYPE_DEBUG > NSCARC_DEBUG_MEDIUM)  WRITE(SCARC_LU,'(a)') 'A_ROW::'
   IF (TYPE_DEBUG > NSCARC_DEBUG_MEDIUM)  WRITE(SCARC_LU,'(10i3)') SC%A_ROW(1:30)
   IF (TYPE_DEBUG > NSCARC_DEBUG_MEDIUM)  WRITE(SCARC_LU,'(a)') 'A_COL::'
   IF (TYPE_DEBUG > NSCARC_DEBUG_MEDIUM)  WRITE(SCARC_LU,'(10i3)') SC%A_COL(1:30)
   IF (TYPE_DEBUG > NSCARC_DEBUG_MEDIUM)  WRITE(SCARC_LU,'(a)') 'A::'
   IF (TYPE_DEBUG > NSCARC_DEBUG_MEDIUM)  WRITE(SCARC_LU,'(10f8.2)') SC%A(1:30)
   IF (TYPE_DEBUG > NSCARC_DEBUG_MEDIUM)  WRITE(SCARC_LU,'(a)') 'RECV_BUF::'
   IF (TYPE_DEBUG > NSCARC_DEBUG_MEDIUM)  WRITE(SCARC_LU,'(10f8.3)') RECV_BUF(1:70)

      IF (NL == 1) THEN

         UNPACK_RECV_SYSTEM: DO
      
            IW = NINT(RECV_BUF(LL))
   IF (TYPE_DEBUG > NSCARC_DEBUG_MEDIUM)  WRITE(SCARC_LU,*) '============ IW=',IW, ' READING RECV_BUF(',LL,')=',RECV_BUF(LL)
            IF (IW==-999) EXIT UNPACK_RECV_SYSTEM
            LL = LL + 1
      
            ICG = SC%GHOST_CELL(IW)
            NROW = SC%A_ROW(ICG)
            NCOL = NINT(RECV_BUF(LL))
            
   IF (TYPE_DEBUG > NSCARC_DEBUG_MEDIUM)  WRITE(SCARC_LU,*) 'NCOL=',NCOL,': NROW=',NROW, ICG , LL
            LL = LL + 1

            COLUMN_SYSTEM_LOOP: DO ICOL = NROW, NROW+NCOL-1
               SC%A_COL(ICOL)= NINT(RECV_BUF(LL))
               SC%A(ICOL)    = RECV_BUF(LL+1)
   IF (TYPE_DEBUG > NSCARC_DEBUG_MEDIUM) &
      WRITE(SCARC_LU,'(i3,a,i3,a,f12.6,a,i3,a,3i5)') NM,': UNPACK_WEIGHTS: A(',ICOL,')=',SC%A(ICOL), ' A_COL(',ICOL,')=',&
                                              SC%A_COL(ICOL),ICOL, ICG
               LL = LL + 2
            ENDDO COLUMN_SYSTEM_LOOP
            NROW = NROW + NCOL
            SC%A_ROW(ICG+1)= NROW
   IF (TYPE_DEBUG > NSCARC_DEBUG_MEDIUM) &
      WRITE(SCARC_LU,'(a,i3,a,3i3)') ' A_ROW(',ICG+1,')=',SC%A_ROW(ICG+1), NROW, NCOL

         ENDDO UNPACK_RECV_SYSTEM
         SC%A_ROW(ICG+1)= NROW

      ENDIF
      

   !!! -------------------------------------------------------------------------------------------
   !!! System matrix
   !!! -------------------------------------------------------------------------------------------
   CASE (NSCARC_MATRIX_SYSTEM2)

SC => SCARC(NM)%COMPACT(NL)
      IF (TYPE_DEBUG >= NSCARC_DEBUG_MEDIUM) write(scarc_lu,'(a,5i4)') '===================================== NM=',NM, SC%NC, NL
      LL = 1
   IF (TYPE_DEBUG > NSCARC_DEBUG_MEDIUM)  WRITE(SCARC_LU,*) 'UNPACK_WEIGHTS: NROW=',NROW
   IF (TYPE_DEBUG > NSCARC_DEBUG_MEDIUM)  WRITE(SCARC_LU,'(a)') 'A_ROW::'
   IF (TYPE_DEBUG > NSCARC_DEBUG_MEDIUM)  WRITE(SCARC_LU,'(10i3)') SC%A_ROW(1:30)
   IF (TYPE_DEBUG > NSCARC_DEBUG_MEDIUM)  WRITE(SCARC_LU,'(a)') 'A_COL::'
   IF (TYPE_DEBUG > NSCARC_DEBUG_MEDIUM)  WRITE(SCARC_LU,'(10i3)') SC%A_COL(1:30)
   IF (TYPE_DEBUG > NSCARC_DEBUG_MEDIUM)  WRITE(SCARC_LU,'(a)') 'A::'
   IF (TYPE_DEBUG > NSCARC_DEBUG_MEDIUM)  WRITE(SCARC_LU,'(10f8.2)') SC%A(1:30)
   IF (TYPE_DEBUG > NSCARC_DEBUG_MEDIUM)  WRITE(SCARC_LU,'(a)') 'RECV_BUF::'
   IF (TYPE_DEBUG > NSCARC_DEBUG_MEDIUM)  WRITE(SCARC_LU,'(10f8.3)') RECV_BUF(1:70)

      IF (NL == 1) THEN

         UNPACK_RECV_SYSTEM2: DO
      
            IW = NINT(RECV_BUF(LL))
   IF (TYPE_DEBUG > NSCARC_DEBUG_MEDIUM)  WRITE(SCARC_LU,*) '============ IW=',IW, ' READING RECV_BUF(',LL,')=',RECV_BUF(LL)
            IF (IW==-999) EXIT UNPACK_RECV_SYSTEM2
            LL = LL + 1
      
            ICG = SC%GHOST_CELL(IW)
            NROW = SC%A_ROW(ICG)
            NCOL = NINT(RECV_BUF(LL))
            
   IF (TYPE_DEBUG > NSCARC_DEBUG_MEDIUM)  WRITE(SCARC_LU,*) 'NCOL=',NCOL,': NROW=',NROW, ICG , LL
            LL = LL + 1

            COLUMN_SYSTEM2_LOOP: DO ICOL = NROW, NROW+NCOL-1
               SC%A_COL(ICOL)= NINT(RECV_BUF(LL))
               SC%A(ICOL)    = RECV_BUF(LL+1)
   IF (TYPE_DEBUG > NSCARC_DEBUG_MEDIUM) &
      WRITE(SCARC_LU,'(i3,a,i3,a,f12.6,a,i3,a,3i5)') NM,': UNPACK_WEIGHTS: A(',ICOL,')=',SC%A(ICOL), ' A_COL(',ICOL,')=',&
                                              SC%A_COL(ICOL),ICOL, ICG
               LL = LL + 2
            ENDDO COLUMN_SYSTEM2_LOOP
            NROW = NROW + NCOL
            SC%A_ROW(ICG+1)= NROW
   IF (TYPE_DEBUG > NSCARC_DEBUG_MEDIUM) &
      WRITE(SCARC_LU,'(a,i3,a,3i3)') ' A_ROW(',ICG+1,')=',SC%A_ROW(ICG+1), NROW, NCOL

         ENDDO UNPACK_RECV_SYSTEM2
         SC%A_ROW(ICG+1)= NROW

      ENDIF
      

   !!! -------------------------------------------------------------------------------------------
   !!! Prolongation matrix
   !!! -------------------------------------------------------------------------------------------
   CASE (NSCARC_MATRIX_PROLONGATION)

      IF (TYPE_DEBUG >= NSCARC_DEBUG_NONE) write(scarc_lu,'(a,5i4)') '==========UNPACK_CMATRIX=================== NM=',NM
   IF (TYPE_DEBUG > NSCARC_DEBUG_NONE)  WRITE(SCARC_LU,'(a,i3)') 'NC  =:', SC%NC 
   IF (TYPE_DEBUG > NSCARC_DEBUG_NONE)  WRITE(SCARC_LU,'(a,i3)') 'NCE =:', SC%NCE
   IF (TYPE_DEBUG > NSCARC_DEBUG_NONE)  WRITE(SCARC_LU,'(a,i3)') 'NCC =:', SC%NCC
   IF (TYPE_DEBUG > NSCARC_DEBUG_NONE)  WRITE(SCARC_LU,'(a,i3)') 'NCCE=:', SC%NCCE
   IF (TYPE_DEBUG > NSCARC_DEBUG_NONE)  WRITE(SCARC_LU,'(a,i3)') 'NCEP=:', SC%NCEP
   IF (TYPE_DEBUG > NSCARC_DEBUG_NONE)  WRITE(SCARC_LU,'(a,i3)') 'RECV_BUF::'
   IF (TYPE_DEBUG > NSCARC_DEBUG_NONE)  WRITE(SCARC_LU,'(7f8.3)') RECV_BUF(1:100)
   IF (TYPE_DEBUG > NSCARC_DEBUG_NONE)  WRITE(SCARC_LU,*) 'P_ROW(',SC%NC+1,')=',SC%P_ROW(SC%NC+1), SC%NC

      SC%NCCE = SC%NCC
      SC%NCEP = SC%NC
      IF (NL == 1) THEN

         ICOL  = SC%P_ROW(SC%NC+1)
         ICE   = SC%NC +1
         ICE2  = SC%NCE+1

         !!! ------------------------------------------------------------------------------------
         !!! Treat first-line neighbors only
         !!! ------------------------------------------------------------------------------------
         LL = 1
         UNPACK_RECV_PROL: DO
      
            NNUM = NINT(RECV_BUF(LL))
            IF (NNUM==-999) EXIT UNPACK_RECV_PROL
            IF (NNUM==-  1) THEN
               LL = LL + 7
               CYCLE UNPACK_RECV_PROL
            ENDIF

            NOM = ABS(NINT(RECV_BUF(LL+1)))
            IC0 = ABS(NINT(RECV_BUF(LL+2)))
            IW  = ABS(NINT(RECV_BUF(LL+3)))
            LL  = LL + 4

            ICA = SC%ADJACENT_CELL(IW)
            ICG = SC%GHOST_CELL(IW)

            SC%NCCE = SC%NCCE+1
            SC%NCEP = SC%NCEP+1


            !!!
            !!! if coarse first-cell-neighbor, take new coarse cell number from neighbor in SC%P
            !!!
            IF (NNUM==1) THEN

               ICF = NINT(RECV_BUF(LL))
               ICC = NINT(RECV_BUF(LL+1))
               PW  = RECV_BUF(LL+2)

               SC%P_PTR(ICE,1) = ICC
               SC%P_PTR(ICE,2) = NOM

               SC%P_COL(ICOL) = - ICC
               SC%P(ICOL) = PW

   IF (TYPE_DEBUG > NSCARC_DEBUG_NONE)  THEN
      WRITE(SCARC_LU,*) '=========== B =============== '
      WRITE(SCARC_LU,'(a,2i3,a,i3,a,4i3)') &
                     'NCCE, NCEP=',SC%NCCE, SC%NCEP,': NUM=',NNUM,': NOM=',NOM, NM, IC0, IW
      WRITE(SCARC_LU,'(a,i3,a,i3,a,i3,a,i3,a,i3,a,i3,a,f12.3,a,i3,a,i3,a,i3,a,f12.3,2i3)') &
           'IW=',IW,': IC0=',IC0,': ICA=',ICA,': ICG=',ICG,': ICF=',ICF,': ICC=',ICC,&
          ': PW=',PW,': P_COL(',ICOL,')=',SC%P_COL(ICOL),': P(',ICOL,')=',SC%P(ICOL), &
           SC%P_PTR(ICE,1), SC%P_PTR(ICE,2)
   ENDIF

               DO ICOL0 = SC%P_ROW(ICA), SC%P_ROW(ICA+1)-1
 WRITE(SCARC_LU,*) 'ICA=',ICA,': ICOL0=',ICOL0,SC%P_COL(ICOL0), ICG
                 IF (SC%P_COL(ICOL0) == - ICF) THEN
                    SC%P_COL(ICOL0) = -ICC
   IF (TYPE_DEBUG > NSCARC_DEBUG_NONE)  THEN
!      write(SCARC_LU,'(a,i3,a,f12.3)') WRITE(SCARC_LU,*) 'ICA: P_COL(',ICOL0,')=',SC%P_COL(ICOL0)
   ENDIF
                 ENDIF
               ENDDO

               IC = ICG
               LL = LL + 3
               ICOL = ICOL + 1

               ICE = ICE + 1
               SC%P_ROW(ICE) = ICOL

            !!!
            !!! if fine first-cell-neighbor, get surrounding cells
            !!!
            ELSE

   IF (TYPE_DEBUG > NSCARC_DEBUG_NONE)  THEN
      WRITE(SCARC_LU,*) '=========== C =============== ', NNUM
      WRITE(SCARC_LU,'(a,2i3,a,i3,a,5i3)') &
                     'NCCE, NCEP=',SC%NCCE, SC%NCEP,': NUM=',NNUM,': NOM=',NOM, NM, IC0, IW, LL
   ENDIF
               !!! take SC%P information from neighbor for ghost cells
               DO INUM = 1, ABS(NNUM)
                  ICF = NINT(RECV_BUF(LL))
                  ICC = NINT(RECV_BUF(LL+1))
                  PW  = RECV_BUF(LL+2)
                  IF (ICF < 0) THEN
                     SC%P_COL(ICOL) = SC%CELLTYPE(ICA)
                  ELSE
                     SC%P_COL(ICOL) = - ICC
                  ENDIF
                  SC%P(ICOL) =  PW
   IF (TYPE_DEBUG > NSCARC_DEBUG_NONE)  THEN
      WRITE(SCARC_LU,'(a,i3,a,i3,a,i3,a,i3,a,i3,a,i3,a,f12.3,a,i3,a,i3,a,i3,a,f12.3)') &
           'IW=',IW,': IC0=',IC0,': ICA=',ICA,': ICG=',ICG,': ICF=',ICF,': ICC=',ICC,&
          ': PW=',PW,': P_COL(',ICOL,')=',SC%P_COL(ICOL),': P(',ICOL,')=',SC%P(ICOL)
   ENDIF
                  IC = ICG
                  LL = LL + 3
                  ICOL = ICOL + 1
               ENDDO
               ICE = ICE + 1
               SC%P_ROW(ICE) = ICOL
            ENDIF
   IF (TYPE_DEBUG > NSCARC_DEBUG_NONE)  WRITE(SCARC_LU,'(a,i3,a,i3)') ' P_ROW(',ICE,')=',SC%P_ROW(ICE)

         ENDDO UNPACK_RECV_PROL
   IF (TYPE_DEBUG > NSCARC_DEBUG_NONE)  WRITE(SCARC_LU,'(a,i3,i3)') 'NCCE_final=',SC%NCCE, SC%NCEP

         !!! ------------------------------------------------------------------------------------
         !!! Then treat second-line neighbors 
         !!! ------------------------------------------------------------------------------------
         LL = 1
         UNPACK_RECV_PROL2: DO
      
            NNUM = NINT(RECV_BUF(LL))
            IF (NNUM==-999) EXIT UNPACK_RECV_PROL2
WRITE(SCARC_LU,*) 'SECOND_LINE: NNUM=',NNUM, LL
            IF (NNUM/=-1) THEN
               LL = LL + 4 + NNUM*3
               CYCLE UNPACK_RECV_PROL2
            ENDIF
WRITE(SCARC_LU,*) 'SECOND_LINE: uff'

            NOM = ABS(NINT(RECV_BUF(LL+1)))
            IC0 = ABS(NINT(RECV_BUF(LL+2)))
            IW  = ABS(NINT(RECV_BUF(LL+3)))
            LL  = LL + 4

            ICA = SC%ADJACENT_CELL(IW)
            ICG = SC%GHOST_CELL(IW)

            SC%NCCE = SC%NCCE+1
            SC%NCEP = SC%NCEP+1

   IF (TYPE_DEBUG > NSCARC_DEBUG_NONE)  WRITE(SCARC_LU,'(a,2i3,a,i3,a,4i3)') 'NCCE, NCEP=',SC%NCCE, SC%NCEP,&
                                                                       ': NUM=',NNUM,': NOM=',NOM, NM, LL

            !!! only take coarse second-cell-neighbor, take new coarse cell number from neighbor in SC%P
            ICF = NINT(RECV_BUF(LL))
            ICC = NINT(RECV_BUF(LL+1))
            PW  = RECV_BUF(LL+2)

            SC%P_PTR(ICE,1) = ICC
            SC%P_PTR(ICE,2) = NOM

            SC%P_COL(ICOL) = - ICC
            SC%P(ICOL) = PW

   IF (TYPE_DEBUG > NSCARC_DEBUG_NONE)  THEN
      WRITE(SCARC_LU,*) '=========== A =============== '
      WRITE(SCARC_LU,'(a,i3,a,i3,a,i3,a,i3,a,i3,a,i3,a,f12.3,a,i3,a,i3,a,i3,a,f12.3,2i3)') &
           'IW=',IW,': IC0=',IC0,': ICA=',ICA,': ICG=',ICG,': ICF=',ICF,': ICC=',ICC,&
          ': PW=',PW,': P_COL(',ICOL,')=',SC%P_COL(ICOL),': P(',ICOL,')=',SC%P(ICOL), &
           SC%P_PTR(SC%NCCE,1), SC%P_PTR(SC%NCCE,2)
   ENDIF

            DO ICOL0 = SC%P_ROW(ICA), SC%P_ROW(ICA+1)-1
 WRITE(SCARC_LU,*) 'ICA=',ICA,': ICOL0=',ICOL0,SC%P_COL(ICOL0), ICG
              IF (SC%P_COL(ICOL0) == - ICF) THEN
                 SC%P_COL(ICOL0) = -ICC
   IF (TYPE_DEBUG > NSCARC_DEBUG_NONE)  THEN
!      write(SCARC_LU,'(a,i3,a,f12.3)') WRITE(SCARC_LU,*) 'ICA: P_COL(',ICOL0,')=',SC%P_COL(ICOL0)
   ENDIF
              ENDIF
            ENDDO

            IC = ICG
            LL = LL + 3
            ICOL = ICOL + 1

            ICE = ICE + 1
            SC%P_ROW(ICE) = ICOL

   IF (TYPE_DEBUG > NSCARC_DEBUG_NONE)  WRITE(SCARC_LU,'(a,i3,a,i3)') ' P_ROW(',ICE,')=',SC%P_ROW(ICE)

         ENDDO UNPACK_RECV_PROL2

      ENDIF

END SELECT
      
END SUBROUTINE SCARC_UNPACK_CMATRIX


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Check if difference of two values is less than a given tolerance
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
LOGICAL FUNCTION MATCH (VAL1, VAL2)
REAL (EB), INTENT(IN) :: VAL1, VAL2
REAL (EB) :: TOL
TOL = 1.0E-10_EB
MATCH = .FALSE.
IF (Abs(VAL1-VAL2) <= TOL) MATCH = .TRUE.
RETURN
END FUNCTION MATCH
 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Debug requested quantity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_DEBUG_QUANTITY(NTYPE, NL, CROUTINE, CNAME)
INTEGER, INTENT(IN) :: NTYPE, NL
INTEGER :: NM, IP, IC, IW, I, J, K
CHARACTER (*), INTENT(IN) :: CROUTINE, CNAME
CHARACTER (20) :: LINE
TYPE (SCARC_BANDED_TYPE) , POINTER :: SB
TYPE (SCARC_COMPACT_TYPE), POINTER :: SC, SC1, SC2

IF (TYPE_DEBUG < NSCARC_DEBUG_MEDIUM) RETURN

SELECT CASE (NTYPE)

   !!!----------------------------------------------------------------------------------------------------
   !!! Debug system matrix A (corresponding to system type)
   !!!----------------------------------------------------------------------------------------------------
   CASE (NSCARC_DEBUG_MATRIX)

      DO NM = NMESHES_MIN, NMESHES_MAX
         WRITE(SCARC_LU,1000) CROUTINE, CNAME, NM, NL
         IF (TYPE_SYSTEM == NSCARC_SYSTEM_BANDED) THEN
            SB  => SCARC(NM)%BANDED(NL)
            WRITE(SCARC_LU,*) '----------------- A(.,.):'
            WRITE(SCARC_LU,*) 'NX=',SB%NX
            WRITE(SCARC_LU,*) 'NY=',SB%NY
            WRITE(SCARC_LU,*) 'NZ=',SB%NZ
            WRITE(SCARC_LU,*) 'NA=',SB%NA
            WRITE(SCARC_LU,*) 'NC=',SB%NC
            IF (TYPE_DIMENSION == NSCARC_DIMENSION_TWO) THEN
               DO IC = 1, SB%NC
                  WRITE(SCARC_LU,'(i7,5f13.2)') IC,SB%A(IC,ILZ),SB%A(IC,ILX),SB%A(IC,ID),SB%A(IC,IUX),SB%A(IC,IUZ)
                  IF (MOD(IC,SB%NX)==0) WRITE(SCARC_LU,*) '----------------------------------------------'
               ENDDO
            ELSE
               DO IC = 1, SB%NC
                  WRITE(SCARC_LU,'(i7,7f13.2)') IC,SB%A(IC,ILZ),SB%A(IC,ILY),SB%A(IC,ILX),SB%A(IC,ID),SB%A(IC,IUX),SB%A(IC,IUY),SB%A(IC,IUZ)
                  IF (MOD(IC,SB%NX)==0) WRITE(SCARC_LU,*) '----------------------------------------------'
               ENDDO
            ENDIF
         ELSE
            SC => SCARC(NM)%COMPACT(NL)
            WRITE(SCARC_LU,*) 'NA=',SC%NA
            WRITE(SCARC_LU,*) 'NC=',SC%NC
            WRITE(SCARC_LU,*) '---------------------- A_ROW:', SC%NC
            WRITE(SCARC_LU,'(4i9)') (SC%A_ROW(IC), IC=1,SC%NC+1)
            WRITE(SCARC_LU,*) '---------------------- A_COL:'
            DO IC = 1, SC%NC
               WRITE(SCARC_LU,'(i5,a,20i9)') IC,':',(SC%A_COL(IP),IP=SC%A_ROW(IC),SC%A_ROW(IC+1)-1)
            ENDDO
            WRITE(SCARC_LU,*) '---------------------- A:'
            DO IC = 1, SC%NC
               WRITE(SCARC_LU,'(i5,a,20f9.2)') IC,':',(SC%A(IP),IP=SC%A_ROW(IC),SC%A_ROW(IC+1)-1)
            ENDDO
            WRITE(SCARC_LU,*) 'SIZE(A) =',SIZE(SC%A)
            WRITE(SCARC_LU,*) 'SIZE(A_COL) =',SIZE(SC%A_COL)
            WRITE(SCARC_LU,*) 'SIZE(A_ROW) =',SIZE(SC%A_ROW)
         ENDIF
      ENDDO

   !!!----------------------------------------------------------------------------------------------------
   !!! Debug extended system matrix A (corresponding to system type)
   !!!----------------------------------------------------------------------------------------------------
   CASE (NSCARC_DEBUG_MATRIXE)

      DO NM = NMESHES_MIN, NMESHES_MAX
         WRITE(SCARC_LU,1000) CROUTINE, CNAME, NM, NL
         IF (TYPE_SYSTEM == NSCARC_SYSTEM_BANDED) THEN
            SB  => SCARC(NM)%BANDED(NL)
            WRITE(SCARC_LU,*) '----------------- A(.,.):'
            DO IC = 1, SB%NC
               WRITE(SCARC_LU,'(5f12.3)') SB%A(IC,ILZ),SB%A(IC,ILX),SB%A(IC,ID),SB%A(IC,IUX),SB%A(IC,IUZ)
               IF (MOD(IC,SB%NX)==0) WRITE(SCARC_LU,*) '----------------------------------------------'
            ENDDO
         ELSE
            SC => SCARC(NM)%COMPACT(NL)
            WRITE(SCARC_LU,*) '---------------------- A_ROW:', SC%NC
            WRITE(SCARC_LU,'(4i9)') (SC%A_ROW(IC), IC=1,SC%NCE+1)
            WRITE(SCARC_LU,*) '---------------------- A_COL:'
            DO IC = 1, SC%NCE
               WRITE(SCARC_LU,'(i5,a,20i9)') IC,':',(SC%A_COL(IP),IP=SC%A_ROW(IC),SC%A_ROW(IC+1)-1)
            ENDDO
            WRITE(SCARC_LU,*) '---------------------- A:'
            DO IC = 1, SC%NCE
               WRITE(SCARC_LU,'(i5,a,20f9.2)') IC,':',(SC%A(IP),IP=SC%A_ROW(IC),SC%A_ROW(IC+1)-1)
            ENDDO
            WRITE(SCARC_LU,*) 'SIZE(A) =',SIZE(SC%A)
            WRITE(SCARC_LU,*) 'SIZE(A_COL) =',SIZE(SC%A_COL)
            WRITE(SCARC_LU,*) 'SIZE(A_ROW) =',SIZE(SC%A_ROW)
         ENDIF
      ENDDO

   !!!----------------------------------------------------------------------------------------------------
   !!! Debug prolongation matrix P (only for compact system)
   !!!----------------------------------------------------------------------------------------------------
   CASE (NSCARC_DEBUG_PROLONGATION)

      DO NM = NMESHES_MIN, NMESHES_MAX
         SC => SCARC(NM)%COMPACT(NL)
         WRITE(SCARC_LU,1000) CROUTINE, CNAME, NM, NL
         WRITE(SCARC_LU,*) '============= P_ROW:', NM
         WRITE(SCARC_LU,'(4i9)') (SC%P_ROW(IC), IC=1, SC%NCE)
         WRITE(SCARC_LU,*) '============= P_COL:', NM
         DO IC=1,SC%NCE+2
            WRITE(SCARC_LU,'(i5,a,10i9)') IC,':',(SC%P_COL(IP), IP=SC%P_ROW(IC),SC%P_ROW(IC+1)-1)
         ENDDO
         WRITE(SCARC_LU,*) '============= P:', NM
         DO IC=1,SC%NCE+2
            WRITE(SCARC_LU,'(i5,a,10f9.2)') IC,':',(SC%P(IP), IP=SC%P_ROW(IC),SC%P_ROW(IC+1)-1)
         ENDDO
         WRITE(SCARC_LU,*) '============= P_PTR:', NM
         DO IC=SC%NC+1,SC%NCE+5
            WRITE(SCARC_LU,'(i5,a,2i8)') IC,':',(SC%P_PTR(IC,IP), IP=1,2)
         ENDDO
      ENDDO

   !!!----------------------------------------------------------------------------------------------------
   !!! Debug restriction matrix R (only for compact system)
   !!!----------------------------------------------------------------------------------------------------
   CASE (NSCARC_DEBUG_RESTRICTION)

      DO NM = NMESHES_MIN, NMESHES_MAX
         SC  => SCARC(NM)%COMPACT(NL)
         WRITE(SCARC_LU,1000) CROUTINE, CNAME, NM, NL
         WRITE(SCARC_LU,*) '============= R_ROW:', NM
         WRITE(SCARC_LU,'(4i9)') (SC%R_ROW(IC), IC=1, SC%NCC)
         WRITE(SCARC_LU,*) '============= R_COL:', NM
         DO IC=1,SC%NCC+2
            WRITE(SCARC_LU,'(i5,a,10i9)') IC,':',(SC%R_COL(IP), IP=SC%R_ROW(IC),SC%R_ROW(IC+1)-1)
         ENDDO
         WRITE(SCARC_LU,*) '============= R:', NM
         DO IC=1,SC%NCC+2
            WRITE(SCARC_LU,'(i5,a,10f9.2)') IC,':',(SC%R(IP), IP=SC%R_ROW(IC),SC%R_ROW(IC+1)-1)
         ENDDO
      ENDDO

   !!!----------------------------------------------------------------------------------------------------
   !!! Debug IJKW
   !!!----------------------------------------------------------------------------------------------------
   CASE (NSCARC_DEBUG_IJKW)

      DO NM = NMESHES_MIN, NMESHES_MAX
         WRITE(SCARC_LU,1000) CROUTINE, CNAME, NM, NL
         IF (TYPE_SYSTEM == NSCARC_SYSTEM_BANDED) THEN
            SB => SCARC(NM)%BANDED(NL)
            DO J = 1,SB%NW
               WRITE(SCARC_LU, '(i4,a,15i6)') J,':',(SB%IJKW(I,J), I=1,15)
            ENDDO
         ELSE
            SC => SCARC(NM)%COMPACT(NL)
            DO J = 1,SC%NW
               WRITE(SCARC_LU, '(i4,a,15i6)') J,':',(SC%IJKW(I,J), I=1,15)
            ENDDO
         ENDIF
      ENDDO

   !!!----------------------------------------------------------------------------------------------------
   !!! Debug PRESSURE_BC_INDEX
   !!!----------------------------------------------------------------------------------------------------
   CASE (NSCARC_DEBUG_BCINDEX)

      DO NM = NMESHES_MIN, NMESHES_MAX
         WRITE(SCARC_LU,1000) CROUTINE, CNAME, NM, NL
         IF (TYPE_SYSTEM == NSCARC_SYSTEM_BANDED) THEN
            SB => SCARC(NM)%BANDED(NL)
            WRITE(SCARC_LU, '(8i4)') (SB%BC_INDEX(J), J=1,SB%NW)
         ELSE
            SC => SCARC(NM)%COMPACT(NL)
            WRITE(SCARC_LU, '(8i4)') (SC%BC_INDEX(J), J=1,SC%NW)
         ENDIF
      ENDDO
      
   !!!----------------------------------------------------------------------------------------------------
   !!! Debug ADJACENT_CELL
   !!!----------------------------------------------------------------------------------------------------
   CASE (NSCARC_DEBUG_ACELL)

      DO NM = NMESHES_MIN, NMESHES_MAX
         WRITE(SCARC_LU,1000) CROUTINE, CNAME, NM, NL
         IF (TYPE_SYSTEM == NSCARC_SYSTEM_BANDED) THEN
            SB => SCARC(NM)%BANDED(NL)
            WRITE(SCARC_LU, '(8i8)') (SB%ADJACENT_CELL(IW), IW=1,SB%NW)
         ELSE
            SC  => SCARC(NM)%COMPACT(NL)
            WRITE(SCARC_LU, '(8i8)') (SC%ADJACENT_CELL(IW), IW=1,SC%NW)
         ENDIF
      ENDDO
      
   !!!----------------------------------------------------------------------------------------------------
   !!! Debug ADJACENT_CELL
   !!!----------------------------------------------------------------------------------------------------
   CASE (NSCARC_DEBUG_GCELL)

      DO NM = NMESHES_MIN, NMESHES_MAX
         WRITE(SCARC_LU,1000) CROUTINE, CNAME, NM, NL
         IF (TYPE_SYSTEM == NSCARC_SYSTEM_BANDED) THEN
            SB => SCARC(NM)%BANDED(NL)
            WRITE(SCARC_LU, '(8i8)') (SB%GHOST_CELL(IW), IW=1,SB%NW)
         ELSE
            SC  => SCARC(NM)%COMPACT(NL)
            WRITE(SCARC_LU, '(8i8)') (SC%GHOST_CELL(IW), IW=1,SC%NW)
         ENDIF
      ENDDO
      
   !!!----------------------------------------------------------------------------------------------------
   !!! Debug SUBDIVISION
   !!!----------------------------------------------------------------------------------------------------
   CASE (NSCARC_DEBUG_SUBDIVISION)

      DO NM = NMESHES_MIN, NMESHES_MAX
         WRITE(SCARC_LU,1000) CROUTINE, CNAME, NM, NL
         IF (TYPE_SYSTEM == NSCARC_SYSTEM_BANDED) THEN
            SB  => SCARC(NM)%BANDED(NL)
            WRITE(SCARC_LU,*) 'IJKW '
            WRITE(SCARC_LU,'(16i4)') (J, (SC%IJKW(I, J), I=1,15), J=1,SC%NW)
            WRITE(SCARC_LU,*) 'SUBDIVISION IOR= 1 '
            WRITE(SCARC_LU,'(3i4)') (SC%SUBDIVISION(I, 1), I=1,3)
            WRITE(SCARC_LU,*) 'SUBDIVISION IOR=-1 '
            WRITE(SCARC_LU,'(3i4)') (SC%SUBDIVISION(I,-1), I=1,3)
            WRITE(SCARC_LU,*) 'SUBDIVISION IOR= 2 '
            WRITE(SCARC_LU,'(3i4)') (SC%SUBDIVISION(I, 2), I=1,3)
            WRITE(SCARC_LU,*) 'SUBDIVISION IOR=-2 '
            WRITE(SCARC_LU,'(3i4)') (SC%SUBDIVISION(I,-2), I=1,3)
            WRITE(SCARC_LU,*) 'SUBDIVISION IOR= 3 '
            WRITE(SCARC_LU,'(3i4)') (SC%SUBDIVISION(I, 3), I=1,3)
            WRITE(SCARC_LU,*) 'SUBDIVISION IOR=-3 '
            WRITE(SCARC_LU,'(3i4)') (SC%SUBDIVISION(I,-3), I=1,3)
         ELSE
            SC  => SCARC(NM)%COMPACT(NL)
            WRITE(SCARC_LU,*) 'IJKW '
            WRITE(SCARC_LU,'(16i4)') (J, (SC%IJKW(I, J), I=1,15), J=1,SC%NW)
            WRITE(SCARC_LU,*) 'SUBDIVISION IOR= 1 '
            WRITE(SCARC_LU,'(3i4)') (SC%SUBDIVISION(I, 1), I=1,3)
            WRITE(SCARC_LU,*) 'SUBDIVISION IOR=-1 '
            WRITE(SCARC_LU,'(3i4)') (SC%SUBDIVISION(I,-1), I=1,3)
            WRITE(SCARC_LU,*) 'SUBDIVISION IOR= 2 '
            WRITE(SCARC_LU,'(3i4)') (SC%SUBDIVISION(I, 2), I=1,3)
            WRITE(SCARC_LU,*) 'SUBDIVISION IOR=-2 '
            WRITE(SCARC_LU,'(3i4)') (SC%SUBDIVISION(I,-2), I=1,3)
            WRITE(SCARC_LU,*) 'SUBDIVISION IOR= 3 '
            WRITE(SCARC_LU,'(3i4)') (SC%SUBDIVISION(I, 3), I=1,3)
            WRITE(SCARC_LU,*) 'SUBDIVISION IOR=-3 '
            WRITE(SCARC_LU,'(3i4)') (SC%SUBDIVISION(I,-3), I=1,3)
         ENDIF
      ENDDO
   
   !!!----------------------------------------------------------------------------------------------------
   !!! Debug MEASURE
   !!!----------------------------------------------------------------------------------------------------
   CASE (NSCARC_DEBUG_MEASURE)

      IF (NMESHES == 1.OR.NL>1) THEN                !!! only temporarily
      DO NM = NMESHES_MIN, NMESHES_MAX
         WRITE(SCARC_LU,1000) CROUTINE, CNAME, NM, NL
         SC  => SCARC(NM)%COMPACT(NL)
         IF (NL == 1) THEN
            DO K = SC%NZ, 1, -1
               WRITE(SCARC_LU, '(10f6.2)') (SC%MEASURE((K-1)*SC%NX*SC%NY+I), I=1, SC%NX)
            ENDDO
         ELSE
            WRITE(SCARC_LU, '(4f6.2)') (SC%MEASURE(IC), IC=1, SC%NC)
         ENDIF
      ENDDO
      ELSE IF (NL==1) THEN
      WRITE(SCARC_LU,1000) CROUTINE, CNAME, 12, NL
      SC1 => SCARC(1)%COMPACT(NL)
      SC2 => SCARC(2)%COMPACT(NL)
      IF (SC1%NX == 5) THEN
      WRITE(SCARC_LU,'(6F6.2,6X,5F6.2)') SC1%MEASURE(16:20),SC1%MEASURE(24),SC2%MEASURE(20),SC2%MEASURE(13:16)
      WRITE(SCARC_LU,'(6F6.2,6X,5F6.2)') SC1%MEASURE(11:15),SC1%MEASURE(23),SC2%MEASURE(19),SC2%MEASURE(9:12)
      WRITE(SCARC_LU,'(6F6.2,6X,5F6.2)') SC1%MEASURE(6:10) ,SC1%MEASURE(22),SC2%MEASURE(18),SC2%MEASURE(5:8)
      WRITE(SCARC_LU,'(6F6.2,6X,5F6.2)') SC1%MEASURE(1:5)  ,SC1%MEASURE(21),SC2%MEASURE(17),SC2%MEASURE(1:4)
      ELSE IF (SC1%NX == 8) THEN
      WRITE(SCARC_LU,'(9F6.2,6X,9F6.2)') SC1%MEASURE(57:64),SC1%MEASURE(72),SC2%MEASURE(72),SC2%MEASURE(57:64)
      WRITE(SCARC_LU,'(9F6.2,6X,9F6.2)') SC1%MEASURE(49:56),SC1%MEASURE(71),SC2%MEASURE(71),SC2%MEASURE(49:56)
      WRITE(SCARC_LU,'(9F6.2,6X,9F6.2)') SC1%MEASURE(41:48),SC1%MEASURE(70),SC2%MEASURE(70),SC2%MEASURE(41:48)
      WRITE(SCARC_LU,'(9F6.2,6X,9F6.2)') SC1%MEASURE(33:40),SC1%MEASURE(69),SC2%MEASURE(69),SC2%MEASURE(33:40)
      WRITE(SCARC_LU,'(9F6.2,6X,9F6.2)') SC1%MEASURE(25:32),SC1%MEASURE(68),SC2%MEASURE(68),SC2%MEASURE(25:32)
      WRITE(SCARC_LU,'(9F6.2,6X,9F6.2)') SC1%MEASURE(17:24),SC1%MEASURE(67),SC2%MEASURE(67),SC2%MEASURE(17:24)
      WRITE(SCARC_LU,'(9F6.2,6X,9F6.2)') SC1%MEASURE(9:16) ,SC1%MEASURE(66),SC2%MEASURE(66),SC2%MEASURE(9:16)
      WRITE(SCARC_LU,'(9F6.2,6X,9F6.2)') SC1%MEASURE(1:8)  ,SC1%MEASURE(65),SC2%MEASURE(65),SC2%MEASURE(1:8)
      ELSE IF (SC1%NX == 9) THEN
      WRITE(SCARC_LU,'(10f6.2,6X,9f6.2)') SC1%MEASURE(64:72),SC1%MEASURE(80),SC2%MEASURE(72),SC2%MEASURE(57:64)
      WRITE(SCARC_LU,'(10f6.2,6X,9f6.2)') SC1%MEASURE(55:63),SC1%MEASURE(79),SC2%MEASURE(71),SC2%MEASURE(49:56)
      WRITE(SCARC_LU,'(10f6.2,6X,9f6.2)') SC1%MEASURE(46:54),SC1%MEASURE(78),SC2%MEASURE(70),SC2%MEASURE(41:48)
      WRITE(SCARC_LU,'(10f6.2,6X,9f6.2)') SC1%MEASURE(37:45),SC1%MEASURE(77),SC2%MEASURE(69),SC2%MEASURE(33:40)
      WRITE(SCARC_LU,'(10f6.2,6X,9f6.2)') SC1%MEASURE(28:36),SC1%MEASURE(76),SC2%MEASURE(68),SC2%MEASURE(25:32)
      WRITE(SCARC_LU,'(10f6.2,6X,9f6.2)') SC1%MEASURE(19:27),SC1%MEASURE(75),SC2%MEASURE(67),SC2%MEASURE(17:24)
      WRITE(SCARC_LU,'(10f6.2,6X,9f6.2)') SC1%MEASURE(10:18),SC1%MEASURE(74),SC2%MEASURE(66),SC2%MEASURE(9:16)
      WRITE(SCARC_LU,'(10f6.2,6X,9f6.2)') SC1%MEASURE(1:9)  ,SC1%MEASURE(73),SC2%MEASURE(65),SC2%MEASURE(1:8)
      ENDIF
      ENDIF

   !!!----------------------------------------------------------------------------------------------------
   !!! Debug CELLTYPE
   !!!----------------------------------------------------------------------------------------------------
   CASE (NSCARC_DEBUG_CELLTYPE)

      IF (NMESHES == 1.OR.NL>1) THEN                !!! only temporarily
      DO NM = NMESHES_MIN, NMESHES_MAX
         WRITE(SCARC_LU,1000) CROUTINE, CNAME, NM, NL
         SC  => SCARC(NM)%COMPACT(NL)
         IF (NL == 1) THEN
            DO K = SC%NZ, 1, -1
               WRITE(SCARC_LU, '(10i5)') (SC%CELLTYPE((K-1)*SC%NX*SC%NY+I), I=1, SC%NX)
            ENDDO
         ELSE
            WRITE(SCARC_LU, '(4i5)') (SC%CELLTYPE(IC), IC=1, SC%NC)
         ENDIF
      ENDDO
      ELSE IF (NL==1) THEN
      WRITE(SCARC_LU,1000) CROUTINE, CNAME, 12, NL
      SC1 => SCARC(1)%COMPACT(NL)
      SC2 => SCARC(2)%COMPACT(NL)
      IF (SC1%NX == 4) THEN
      WRITE(SCARC_LU,'(5I6,6X,5I6)') SC1%CELLTYPE(13:16),SC1%CELLTYPE(20),SC2%CELLTYPE(20),SC2%CELLTYPE(13:16)
      WRITE(SCARC_LU,'(5I6,6X,5I6)') SC1%CELLTYPE(9:12),SC1%CELLTYPE(19),SC2%CELLTYPE(19),SC2%CELLTYPE(9:12)
      WRITE(SCARC_LU,'(5I6,6X,5I6)') SC1%CELLTYPE(5:8) ,SC1%CELLTYPE(18),SC2%CELLTYPE(18),SC2%CELLTYPE(5:8)
      WRITE(SCARC_LU,'(5I6,6X,5I6)') SC1%CELLTYPE(1:4)  ,SC1%CELLTYPE(17),SC2%CELLTYPE(17),SC2%CELLTYPE(1:4)
      ELSE IF (SC1%NX == 5) THEN
      WRITE(SCARC_LU,'(6I6,6X,5I6)') SC1%CELLTYPE(16:20),SC1%CELLTYPE(24),SC2%CELLTYPE(20),SC2%CELLTYPE(13:16)
      WRITE(SCARC_LU,'(6I6,6X,5I6)') SC1%CELLTYPE(11:15),SC1%CELLTYPE(23),SC2%CELLTYPE(19),SC2%CELLTYPE(9:12)
      WRITE(SCARC_LU,'(6I6,6X,5I6)') SC1%CELLTYPE(6:10) ,SC1%CELLTYPE(22),SC2%CELLTYPE(18),SC2%CELLTYPE(5:8)
      WRITE(SCARC_LU,'(6I6,6X,5I6)') SC1%CELLTYPE(1:5)  ,SC1%CELLTYPE(21),SC2%CELLTYPE(17),SC2%CELLTYPE(1:4)
      ELSE IF (SC1%NX == 8) THEN
      WRITE(SCARC_LU,'(9I6,6X,9I6)') SC1%CELLTYPE(57:64),SC1%CELLTYPE(72),SC2%CELLTYPE(72),SC2%CELLTYPE(57:64)
      WRITE(SCARC_LU,'(9I6,6X,9I6)') SC1%CELLTYPE(49:56),SC1%CELLTYPE(71),SC2%CELLTYPE(71),SC2%CELLTYPE(49:56)
      WRITE(SCARC_LU,'(9I6,6X,9I6)') SC1%CELLTYPE(41:48),SC1%CELLTYPE(70),SC2%CELLTYPE(70),SC2%CELLTYPE(41:48)
      WRITE(SCARC_LU,'(9I6,6X,9I6)') SC1%CELLTYPE(33:40),SC1%CELLTYPE(69),SC2%CELLTYPE(69),SC2%CELLTYPE(33:40)
      WRITE(SCARC_LU,'(9I6,6X,9I6)') SC1%CELLTYPE(25:32),SC1%CELLTYPE(68),SC2%CELLTYPE(68),SC2%CELLTYPE(25:32)
      WRITE(SCARC_LU,'(9I6,6X,9I6)') SC1%CELLTYPE(17:24),SC1%CELLTYPE(67),SC2%CELLTYPE(67),SC2%CELLTYPE(17:24)
      WRITE(SCARC_LU,'(9I6,6X,9I6)') SC1%CELLTYPE(9:16) ,SC1%CELLTYPE(66),SC2%CELLTYPE(66),SC2%CELLTYPE(9:16)
      WRITE(SCARC_LU,'(9I6,6X,9I6)') SC1%CELLTYPE(1:8)  ,SC1%CELLTYPE(65),SC2%CELLTYPE(65),SC2%CELLTYPE(1:8)
      ELSE IF (SC1%NX == 9) THEN
      WRITE(SCARC_LU,'(10I6,6X,9I6)') SC1%CELLTYPE(64:72),SC1%CELLTYPE(80),SC2%CELLTYPE(72),SC2%CELLTYPE(57:64)
      WRITE(SCARC_LU,'(10I6,6X,9I6)') SC1%CELLTYPE(55:63),SC1%CELLTYPE(79),SC2%CELLTYPE(71),SC2%CELLTYPE(49:56)
      WRITE(SCARC_LU,'(10I6,6X,9I6)') SC1%CELLTYPE(46:54),SC1%CELLTYPE(78),SC2%CELLTYPE(70),SC2%CELLTYPE(41:48)
      WRITE(SCARC_LU,'(10I6,6X,9I6)') SC1%CELLTYPE(37:45),SC1%CELLTYPE(77),SC2%CELLTYPE(69),SC2%CELLTYPE(33:40)
      WRITE(SCARC_LU,'(10I6,6X,9I6)') SC1%CELLTYPE(28:36),SC1%CELLTYPE(76),SC2%CELLTYPE(68),SC2%CELLTYPE(25:32)
      WRITE(SCARC_LU,'(10I6,6X,9I6)') SC1%CELLTYPE(19:27),SC1%CELLTYPE(75),SC2%CELLTYPE(67),SC2%CELLTYPE(17:24)
      WRITE(SCARC_LU,'(10I6,6X,9I6)') SC1%CELLTYPE(10:18),SC1%CELLTYPE(74),SC2%CELLTYPE(66),SC2%CELLTYPE(9:16)
      WRITE(SCARC_LU,'(10I6,6X,9I6)') SC1%CELLTYPE(1:9)  ,SC1%CELLTYPE(73),SC2%CELLTYPE(65),SC2%CELLTYPE(1:8)
      ENDIF
      ENDIF

   !!!----------------------------------------------------------------------------------------------------
   !!! Debug CELLTYPE
   !!!----------------------------------------------------------------------------------------------------
   CASE (NSCARC_DEBUG_COARSE)

      DO NM = NMESHES_MIN, NMESHES_MAX
         WRITE(SCARC_LU,1000) CROUTINE, CNAME, NM, NL
         SC  => SCARC(NM)%COMPACT(NL)
         DO K = SC%NZ,1,-1
            DO I = 1,SC%NX
               IC = (K-1)*SC%NX + I
               IF (SC%CELLTYPE(IC) < 0) THEN
                  LINE(I:I) = 'O'
               ELSE IF (SC%CELLTYPE(IC) > 0) THEN
                  LINE(I:I) = 'X'
               ENDIF
            ENDDO
            WRITE(SCARC_LU,'(12A4)') (LINE(I:I), I=1,SC%NX)
         ENDDO
      ENDDO

END SELECT

1000 FORMAT('======================================================================================',/, &
            '=== ',A20,' : ', A25,' for mesh ',i3,' on level ', i3, /, &
            '======================================================================================')
 
END SUBROUTINE SCARC_DEBUG_QUANTITY

 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Debug full vector information on level NL 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_DEBUG_VECTOR (NVECTOR, CROUTINE, CNAME)
INTEGER, INTENT(IN):: NVECTOR
REAL (EB), POINTER, DIMENSION(:,:,:) :: VB
REAL (EB), POINTER, DIMENSION(:)     :: VC
REAL (EB):: VALUES(10)
INTEGER :: NM, II, JJ, KK, IC, NX8, NY8, NZ8
CHARACTER (*), INTENT(IN) :: CROUTINE, CNAME
TYPE (MESH_TYPE), POINTER :: M
 
IF (TYPE_DEBUG < NSCARC_DEBUG_LESS) RETURN

SELECT_SYSTEM: SELECT CASE (TYPE_SYSTEM)

   !!!----------------------------------------------------------------------------------------------------
   !!! Banded system
   !!!----------------------------------------------------------------------------------------------------
   CASE (NSCARC_SYSTEM_BANDED)

      DO NM = NMESHES_MIN, NMESHES_MAX

         M  => MESHES(NM)
         VB => POINT_TO_BVECTOR (NVECTOR, NM, NLEVEL_MIN)

         NX8=MIN(8,M%IBAR)
         NY8=MIN(8,M%JBAR)
         NZ8=MIN(8,M%KBAR)
         
         WRITE(SCARC_LU,*) '=============================================================================='
         WRITE(SCARC_LU,2000) CROUTINE, CNAME, NM, NLEVEL_MIN
         WRITE(SCARC_LU,*) '=============================================================================='
         SELECT_BANDED_DIMENSION: SELECT CASE (TYPE_DIMENSION)
      
            CASE (NSCARC_DIMENSION_TWO)
      
               DO KK = NZ8, 1, - 1
                  DO II=1,NX8
                     IF (ABS(VB(II,1,KK))<1.0E-14_EB) THEN
                        VALUES(II)=0.0_EB
                     ELSE
                        VALUES(II)=VB(II,1,KK)
                     ENDIF
                  ENDDO
                  WRITE(SCARC_LU, '(a,i3,a,10e13.5)') 'J= ',JJ,' : ',(VALUES(II), II=1, NX8)
               ENDDO
               WRITE(SCARC_LU,*)  '------------------------------------------------',&
                                  '---------------------------------------------------'
      
            CASE (NSCARC_DIMENSION_THREE)
               DO KK = NZ8, 1, - 1
                  WRITE(SCARC_LU,'(a,i3,a)')  '----------------------------------------  K = ', KK,&
                                              '----------------------------------------'
                  DO JJ = NY8, 1, - 1
                     DO II=1,NX8
                        IF (ABS(VB(II,JJ,KK))<1.0E-14_EB) THEN
                           VALUES(II)=0.0_EB
                        ELSE
                           VALUES(II)=VB(II,JJ,KK)
                        ENDIF
                     ENDDO
                     WRITE(SCARC_LU, '(a,i3,a,10e13.5)') 'J= ',JJ,' : ',(VALUES(II), II=1, NX8)
                  ENDDO
               ENDDO
               WRITE(SCARC_LU,*)  '------------------------------------------------',&
                                  '---------------------------------------------------'
         END SELECT SELECT_BANDED_DIMENSION
      
      ENDDO
          
   !!!----------------------------------------------------------------------------------------------------
   !!! Compact system
   !!!----------------------------------------------------------------------------------------------------
   CASE (NSCARC_SYSTEM_COMPACT)

      DO NM = NMESHES_MIN, NMESHES_MAX
   
         M  => MESHES(NM)
         VC => POINT_TO_CVECTOR (NVECTOR, NM, NLEVEL_MIN)
         
         NX8=MIN(8,M%IBAR)
         NY8=MIN(8,M%JBAR)
         NZ8=MIN(8,M%KBAR)
         
         WRITE(SCARC_LU,*) '================================================================================'
         WRITE(SCARC_LU,2000) CROUTINE, CNAME, NM, NLEVEL_MIN
         WRITE(SCARC_LU,*) '================================================================================'
         SELECT_COMPACT_DIMENSION: SELECT CASE (TYPE_DIMENSION)
         
            CASE (NSCARC_DIMENSION_TWO)
               DO KK = NZ8, 1, - 1
                  DO II=1,NX8
                     IC = (KK-1)*M%IBAR + II
                     IF (ABS(VC(IC))<1.0E-14_EB) THEN
                        VALUES(II)=0.0_EB
                     ELSE
                        VALUES(II)=VC(IC)
                     ENDIF
                  ENDDO
                  WRITE(SCARC_LU, '(a,i3,a,8e13.5)') 'J= ',1,' : ',(VALUES(II), II=1, NX8)
               ENDDO
               WRITE(SCARC_LU,*)  '------------------------------------------------',&
                                  '---------------------------------------------------'
         
            CASE (NSCARC_DIMENSION_THREE)
               DO KK = NZ8, 1, - 1
                  WRITE(SCARC_LU,'(a,i3,a)')  '----------------------------------------  K = ', KK,&
                                              '----------------------------------------'
                  DO JJ = NY8, 1, - 1
                     DO II=1,NX8
                        IC = (KK-1)*M%IBAR*M%JBAR + (JJ-1)*M%IBAR + II
                        IF (ABS(VC(IC))<1.0E-14_EB) THEN
                           VALUES(II)=0.0_EB
                        ELSE
                           VALUES(II)=VC(IC)
                        ENDIF
                     ENDDO
                     WRITE(SCARC_LU, '(a,i3,a,8e13.5)') 'J= ',JJ,' : ',(VALUES(II), II=1, NX8)
                  ENDDO
               ENDDO
               WRITE(SCARC_LU,*)  '------------------------------------------------',&
                                  '---------------------------------------------------'
      END SELECT SELECT_COMPACT_DIMENSION

   ENDDO

END SELECT SELECT_SYSTEM
      
2000 FORMAT('=== ',A,' : ',A,' on mesh ',I4,' on level ',I4) 
END SUBROUTINE SCARC_DEBUG_VECTOR


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Only for debugging reasons: print out vector information on level NL 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_DEBUG_LEVEL (NVECTOR, CROUTINE, CNAME, NL)
INTEGER, INTENT(IN):: NVECTOR, NL
REAL (EB), POINTER, DIMENSION(:,:,:) :: VB
REAL (EB), POINTER, DIMENSION(:)     :: VC
REAL (EB):: VALUES(10)
INTEGER :: NM, II, JJ, KK, IC, NX8, NY8, NZ8
CHARACTER (*), INTENT(IN) :: CROUTINE, CNAME
TYPE (SCARC_BANDED_TYPE) , POINTER :: SB
TYPE (SCARC_COMPACT_TYPE), POINTER :: SC
 
IF (TYPE_DEBUG < NSCARC_DEBUG_MEDIUM) RETURN

SELECT_SYSTEM: SELECT CASE (TYPE_SYSTEM)

   !!!----------------------------------------------------------------------------------------------------
   !!! Banded system
   !!!----------------------------------------------------------------------------------------------------
   CASE (NSCARC_SYSTEM_BANDED)

      DO NM = NMESHES_MIN, NMESHES_MAX

         SB => SCARC(NM)%BANDED(NL)
         VB => POINT_TO_BVECTOR (NVECTOR, NM, NL)

         NX8=MIN(8,SB%NX)
         NY8=MIN(8,SB%NY)
         NZ8=MIN(8,SB%NZ)
         
         WRITE(SCARC_LU,*) '=============================================================================='
         WRITE(SCARC_LU,2000) CROUTINE, CNAME, NM, NL, NX8, NY8, NZ8
         WRITE(SCARC_LU,*) '=============================================================================='
         SELECT_BANDED_DIMENSION: SELECT CASE (TYPE_DIMENSION)
         
            CASE (NSCARC_DIMENSION_TWO)
               DO KK = NZ8, 1, - 1
                  DO II=1,NX8
                     IF (ABS(VB(II,1,KK))<1.0E-14_EB) THEN
                        VALUES(II)=0.0_EB
                     ELSE
                        VALUES(II)=VB(II,1,KK)
                     ENDIF
                  ENDDO
                  WRITE(SCARC_LU, '(a,i3,a,10e13.5)') 'J= ',1,' : ',(VALUES(II), II=1, NX8)
               ENDDO
               WRITE(SCARC_LU,*)  '------------------------------------------------',&
                                  '---------------------------------------------------'
      
            CASE (NSCARC_DIMENSION_THREE)
               DO KK = NZ8, 1, - 1
                  WRITE(SCARC_LU,'(a,i3,a)')  '----------------------------------------  K = ', KK,&
                                              '----------------------------------------'
                  DO JJ = NY8, 1, - 1
                     DO II=1,NX8
                        IF (ABS(VB(II,JJ,KK))<1.0E-14_EB) THEN
                           VALUES(II)=0.0_EB
                        ELSE
                           VALUES(II)=VB(II,JJ,KK)
                        ENDIF
                     ENDDO
                     WRITE(SCARC_LU, '(a,i3,a,10e13.5)') 'J= ',JJ,' : ',(VALUES(II), II=1, NX8)
                  ENDDO
               ENDDO
               WRITE(SCARC_LU,*)  '------------------------------------------------',&
                                  '---------------------------------------------------'

         END SELECT SELECT_BANDED_DIMENSION

      ENDDO
         
   !!!----------------------------------------------------------------------------------------------------
   !!! Compact system
   !!!----------------------------------------------------------------------------------------------------
   CASE (NSCARC_SYSTEM_COMPACT)

      DO NM = NMESHES_MIN, NMESHES_MAX

         SC => SCARC(NM)%COMPACT(NL)
         VC => POINT_TO_CVECTOR (NVECTOR, NM, NL)

         WRITE(SCARC_LU,*) '=============================================================================='
         WRITE(SCARC_LU,2001) CROUTINE, CNAME, NM, NL, SC%NC
         WRITE(SCARC_LU,*) '=============================================================================='
         IF (NMESHES == 1) THEN
            IF (NL == 1) THEN
               WRITE(SCARC_LU, '(8e13.5)') (VC(IC), IC=1, SC%NCE)
            ELSE IF (NL == 2) THEN
               WRITE(SCARC_LU, '(4e13.5)') (VC(IC), IC=1, SC%NCE)
            ELSE IF (NL == 3) THEN
               WRITE(SCARC_LU, '(2e13.5)') (VC(IC), IC=1, SC%NCE)
            ENDIF
         ELSE
            IF (NL == 1) THEN
               WRITE(SCARC_LU, '(4e13.5)') (VC(IC), IC=1, SC%NCE)
            ELSE IF (NL == 2) THEN
               WRITE(SCARC_LU, '(2e13.5)') (VC(IC), IC=1, SC%NCE)
            ELSE IF (NL == 3) THEN
               WRITE(SCARC_LU, '(1e13.5)') (VC(IC), IC=1, SC%NCE)
            ENDIF
         ENDIF
         WRITE(SCARC_LU,*)  '------------------------------------------------',&
                            '---------------------------------------------------'
      ENDDO

END SELECT SELECT_SYSTEM

2000 FORMAT('=== ',A,' : ',A,' on mesh ',I4,' on level ',I4, ': NX, NY, NZ=',3i3)
2001 FORMAT('=== ',A,' : ',A,' on mesh ',I4,' on level ',I4, ': NC=',3i3)
END SUBROUTINE SCARC_DEBUG_LEVEL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Produce latex information about coarsening of grid
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_LATEX_INFORMATION(NTYPE, NL)
INTEGER, INTENT(IN) :: NTYPE, NL
INTEGER :: NM, NL0, MLATEX, I, K, IC, IC1, IC2, IC3, IC4, IC5, IC6
CHARACTER (60) :: CLATEX


SELECT CASE(NTYPE)

   !!!-------------------------------------------------------------------------------------------------
   !!! Produce Latex information about grid coarsening with different node sizes for different levels
   !!!-------------------------------------------------------------------------------------------------
   CASE(NSCARC_LATEX_STAGGERED)

      DO NM = NMESHES_MIN, NMESHES_MAX
      
         WRITE (CLATEX, '(A,A,A,i2.2,A,i2.2,A)') 'latex/',TRIM(CHID),'_mesh',NM,'_level',NL+1,'_s.tex'
         MLATEX=GET_FILE_NUMBER()
         OPEN(MLATEX,FILE=CLATEX)
         WRITE(MLATEX,1001)
         WRITE(MLATEX,2001)
      
         NL0 = NLEVEL_MIN
         DO K = SCARC(NM)%COMPACT(NL0)%NZ,1,-1
            DO I = 1,SCARC(NM)%COMPACT(NL0)%NX
               IC = (K-1)*SCARC(NM)%COMPACT(NL0)%NX + I
               WRITE(MLATEX,3002) 'L1','L1',I,K,I,K
               IC1=SCARC(NM)%COMPACT(NL0)%CELLTYPE(IC)
               IF (IC1>0) THEN
                   WRITE(MLATEX,3002) 'L2','L2',I,K,I,K
                   IF (SCARC(NM)%COMPACT(NL0)%NCC>=1) THEN
                   IF (NL==1) CYCLE
                   IC2=SCARC(NM)%COMPACT(NL0+1)%CELLTYPE(IC1)
                   IF (IC2>0) THEN
                      WRITE(MLATEX,3002) 'L3','L3',I,K,I,K
                      IF (SCARC(NM)%COMPACT(NL0+1)%NCC>=1) THEN
                      IF (NL==2) CYCLE
                      IC3=SCARC(NM)%COMPACT(NL0+2)%CELLTYPE(IC2)
                      IF (IC3>0) THEN
                         WRITE(MLATEX,3002) 'L4','L4',I,K,I,K
                         IF (SCARC(NM)%COMPACT(NL0+2)%NCC>=1) THEN
                         IF (NL==3) CYCLE
                         IC4=SCARC(NM)%COMPACT(NL0+3)%CELLTYPE(IC3)
                         IF (IC4>0) THEN
                            WRITE(MLATEX,3002) 'L5','L5',I,K,I,K
                            IF (SCARC(NM)%COMPACT(NL0+3)%NCC>=1) THEN
                            IF (NL==4) CYCLE
                            IC5=SCARC(NM)%COMPACT(NL0+4)%CELLTYPE(IC4)
                            IF (IC5>0) THEN
                               WRITE(MLATEX,3002) 'L6','L6',I,K,I,K
                               IF (SCARC(NM)%COMPACT(NL0+3)%NCC>=1) THEN
                               IF (NL==5) CYCLE
                               IC6=SCARC(NM)%COMPACT(NL0+4)%CELLTYPE(IC4)
                               IF (IC6>0) THEN
                                  WRITE(MLATEX,3002) 'L7','L7',I,K,I,K
                               ENDIF
                               ENDIF
                            ENDIF
                            ENDIF
                         ENDIF
                         ENDIF
                      ENDIF
                      ENDIF
                   ENDIF
                   ENDIF
               ENDIF
            ENDDO
         ENDDO
      
         WRITE(MLATEX,1002)
         CLOSE(MLATEX)
      
      ENDDO
      
      
   !!!-------------------------------------------------------------------------------------------------
   !!! Produce Latex information about grid coarsening with same node sizes for different levels
   !!!-------------------------------------------------------------------------------------------------
   CASE(NSCARC_LATEX_EQUAL)

      DO NM = NMESHES_MIN, NMESHES_MAX
      
         WRITE (CLATEX, '(A,A,A,i2.2,A,i2.2,A)') 'latex/',TRIM(CHID),'_mesh',NM,'_level',NL+1,'_e.tex'
         MLATEX=GET_FILE_NUMBER()
         OPEN(MLATEX,FILE=CLATEX)
         WRITE(MLATEX,1001)
         WRITE(MLATEX,3001)
      
         NL0 = NLEVEL_MIN
         DO K = SCARC(NM)%COMPACT(NL0)%NZ,1,-1
            DO I = 1,SCARC(NM)%COMPACT(NL0)%NX
               IC = (K-1)*SCARC(NM)%COMPACT(NL0)%NX + I
               WRITE(MLATEX,3002) 'L1','L1',I,K,I,K
               IC1=SCARC(NM)%COMPACT(NL0)%CELLTYPE(IC)
               IF (IC1>0) THEN
                   WRITE(MLATEX,3002) 'L2','L2',I,K,I,K
                   IF (SCARC(NM)%COMPACT(NL0)%NCC>=1) THEN
                   IF (NL==1) CYCLE
                   IC2=SCARC(NM)%COMPACT(NL0+1)%CELLTYPE(IC1)
                   IF (IC2>0) THEN
                      WRITE(MLATEX,3002) 'L3','L3',I,K,I,K
                      IF (SCARC(NM)%COMPACT(NL0+1)%NCC>=1) THEN
                      IF (NL==2) CYCLE
                      IC3=SCARC(NM)%COMPACT(NL0+2)%CELLTYPE(IC2)
                      IF (IC3>0) THEN
                         WRITE(MLATEX,3002) 'L4','L4',I,K,I,K
                         IF (SCARC(NM)%COMPACT(NL0+2)%NCC>=1) THEN
                         IF (NL==3) CYCLE
                         IC4=SCARC(NM)%COMPACT(NL0+3)%CELLTYPE(IC3)
                         IF (IC4>0) THEN
                            WRITE(MLATEX,3002) 'L5','L5',I,K,I,K
                            IF (SCARC(NM)%COMPACT(NL0+3)%NCC>=1) THEN
                            IF (NL==4) CYCLE
                            IC5=SCARC(NM)%COMPACT(NL0+4)%CELLTYPE(IC4)
                            IF (IC5>0) THEN
                               WRITE(MLATEX,3002) 'L6','L6',I,K,I,K
                               IF (SCARC(NM)%COMPACT(NL0+3)%NCC>=1) THEN
                               IF (NL==5) CYCLE
                               IC6=SCARC(NM)%COMPACT(NL0+4)%CELLTYPE(IC4)
                               IF (IC6>0) THEN
                                  WRITE(MLATEX,3002) 'L7','L7',I,K,I,K
                               ENDIF
                               ENDIF
                            ENDIF
                            ENDIF
                         ENDIF
                         ENDIF
                      ENDIF
                      ENDIF
                   ENDIF
                   ENDIF
               ENDIF
            ENDDO
         ENDDO
      
         WRITE(MLATEX,1002)
         CLOSE(MLATEX)
      
      ENDDO
      
   !!!-------------------------------------------------------------------------------------------------
   !!! Produce Latex information about grid coarsening with different node sizes for different levels
   !!!-------------------------------------------------------------------------------------------------
   CASE(NSCARC_LATEX_NUMBER)

      DO NM = NMESHES_MIN, NMESHES_MAX
      
         IF (NL == NLEVEL_MIN) THEN
            WRITE (CLATEX, '(A,A,A,i2.2,A,i2.2,A)') 'latex/',TRIM(CHID),'_mesh',NM,'_level',NL,'_n.tex'
            MLATEX=GET_FILE_NUMBER()
            OPEN(MLATEX,FILE=CLATEX)
            WRITE(MLATEX,1001)
            WRITE(MLATEX,4001)
            NL0 = NLEVEL_MIN
            DO K = SCARC(NM)%COMPACT(NL0)%NZ,1,-1
               DO I = 1,SCARC(NM)%COMPACT(NL0)%NX
                  IC = (K-1)*SCARC(NM)%COMPACT(NL0)%NX + I
WRITE(SCARC_LU,*) 'IC=',IC,I,K
                  WRITE(MLATEX,4002) 'L1','L1',I,K,I,K, IC
               ENDDO
            ENDDO
            WRITE(MLATEX,1002)
            CLOSE(MLATEX)
         ENDIF

         WRITE (CLATEX, '(A,A,A,i2.2,A,i2.2,A)') 'latex/',TRIM(CHID),'_mesh',NM,'_level',NL+1,'_n.tex'
         MLATEX=GET_FILE_NUMBER()
         OPEN(MLATEX,FILE=CLATEX)
         WRITE(MLATEX,1001)
         WRITE(MLATEX,4001)
      
         NL0 = NLEVEL_MIN
         DO K = SCARC(NM)%COMPACT(NL0)%NZ,1,-1
            DO I = 1,SCARC(NM)%COMPACT(NL0)%NX
               IC = (K-1)*SCARC(NM)%COMPACT(NL0)%NX + I
               IC1=SCARC(NM)%COMPACT(NL0)%CELLTYPE(IC)
               IF (IC1> 0) THEN
                  WRITE(MLATEX,4002) 'L2','L2',I,K,I,K, IC1
                  IF (SCARC(NM)%COMPACT(NL0)%NCC>=1) THEN
                  IF (NL==1) CYCLE
                  IC2=SCARC(NM)%COMPACT(NL0+1)%CELLTYPE(IC1)
                  IF (IC2>0) THEN
                     WRITE(MLATEX,4002) 'L2','L2',I,K,I,K, IC2
                     IF (SCARC(NM)%COMPACT(NL0+1)%NCC>=1) THEN
                     IF (NL==2) CYCLE
                     IC3=SCARC(NM)%COMPACT(NL0+2)%CELLTYPE(IC2)
                     IF (IC3>0) THEN
                        WRITE(MLATEX,4002) 'L3','L3',I,K,I,K, IC3
                        IF (SCARC(NM)%COMPACT(NL0+2)%NCC>=1) THEN
                        IF (NL==3) CYCLE
                        IC4=SCARC(NM)%COMPACT(NL0+3)%CELLTYPE(IC3)
                        IF (IC4>0) THEN
                           WRITE(MLATEX,4002) 'L4','L4',I,K,I,K, IC4
                           IF (SCARC(NM)%COMPACT(NL0+3)%NCC>=1) THEN
                           IF (NL==4) CYCLE
                           IC5=SCARC(NM)%COMPACT(NL0+4)%CELLTYPE(IC4)
                           IF (IC5>0) THEN
                              WRITE(MLATEX,4002) 'L5','L5',I,K,I,K, IC5
                              IF (SCARC(NM)%COMPACT(NL0+3)%NCC>=1) THEN
                              IF (NL==5) CYCLE
                              IC6=SCARC(NM)%COMPACT(NL0+4)%CELLTYPE(IC4)
                              IF (IC6>0) THEN
                                 WRITE(MLATEX,4002) 'L6','L6',I,K,I,K, IC6
                                 WRITE(MLATEX,4002) 'L7','L7',I,K,I,K, IC6
                              ENDIF
                              ENDIF
                           ENDIF
                           ENDIF
                        ENDIF
                        ENDIF
                     ENDIF
                     ENDIF
                  ENDIF
                  ENDIF
               ELSE
                  WRITE(MLATEX,4002) 'L1','L1',I,K,I,K, IC1
               ENDIF
            ENDDO
         ENDDO
      
         WRITE(MLATEX,1002)
         CLOSE(MLATEX)
      
      ENDDO
END SELECT

1001 FORMAT('\documentclass[11pt]{minimal}' ,/,&
            '\usepackage{tikz}' ,/,&
            '\usetikzlibrary{positioning}' ,/,&
            '\begin{document}' ,/,&
            '\begin{tikzpicture}')
1002 FORMAT('\end{tikzpicture}',/,&
            '\end{document}')

2001 FORMAT ('[inner sep=1mm,',/,&
            ' scale=1.0,',/,&
            ' L7/.style={rectangle, draw=black!40, fill=yellow!100  , thin, minimum size=12mm},',/,&
            ' L6/.style={rectangle, draw=black!40, fill=black!80  , thin, minimum size=10mm},',/,&
            ' L5/.style={rectangle, draw=black!40, fill=blue!80 , thin, minimum size=8mm},',/,&
            ' L4/.style={rectangle, draw=black!40, fill=green!50 , thin, minimum size=6mm},',/,&
            ' L3/.style={rectangle, draw=black!40, fill=red!50   , thin, minimum size=4mm},',/,&
            ' L2/.style={rectangle, draw=black!40, fill=blue!40  , thin, minimum size=2mm},',/,&
            ' L1/.style={rectangle, draw=black!40, fill=black!02, thin, minimum size=10mm}],')
2002 FORMAT('\node [',A,']  (',A2,I2.2,i2.2,')  at (',I2,',',I2,')  {};')

3001 FORMAT ('[inner sep=0.1mm,',/,&
            ' scale=0.5,',/,&
            ' L6/.style={rectangle, draw=black!40, fill=black!80  , thin, minimum size=4mm},',/,&
            ' L5/.style={rectangle, draw=black!40, fill=blue!80 , thin, minimum size=2.4mm},',/,&
            ' L4/.style={rectangle, draw=black!40, fill=green!50 , thin, minimum size=2.6mm},',/,&
            ' L3/.style={rectangle, draw=black!40, fill=red!50   , thin, minimum size=1.8mm},',/,&
            ' L2/.style={rectangle, draw=black!40, fill=blue!40  , thin, minimum size=1mm},',/,&
            ' L1/.style={rectangle, draw=black!40, fill=black!02, thin, minimum size=5mm}],')
3002 FORMAT('\node [',A,']  (',A2,I2.2,i2.2,')  at (',I2,',',I2,')  {};')

4001 FORMAT ('[inner sep=1mm,',/,&
            ' scale=1.0,',/,&
            ' L7/.style={rectangle, draw=black!40, fill=yellow!100,thin, minimum size=10mm},',/,&
            ' L6/.style={rectangle, draw=black!40, fill=red!80  , thin, minimum size=10mm},',/,&
            ' L5/.style={rectangle, draw=black!40, fill=black!60 , thin, minimum size=10mm},',/,&
            ' L4/.style={rectangle, draw=black!40, fill=green!50 , thin, minimum size=10mm},',/,&
            ' L3/.style={rectangle, draw=black!40, fill=red!50   , thin, minimum size=10mm},',/,&
            ' L2/.style={rectangle, draw=black!40, fill=blue!40  , thin, minimum size=10mm},',/,&
            ' L1/.style={rectangle, draw=black!40, fill=black!02, thin, minimum size=10mm}],')
4002 FORMAT('\node [',A,']  (',A2,I2.2,i2.2,')  at (',I2,',',I2,')  {',i2,'};')
 
END SUBROUTINE SCARC_LATEX_INFORMATION


 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Print out timings for ScaRC - not updated at the moment
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_TIMINGS
INTEGER:: NM, IERR, I
INTEGER , ALLOCATABLE, DIMENSION(:)   :: COUNTS_SCARC_TIMERS, DISPLS_SCARC_TIMERS
REAL(EB), ALLOCATABLE, DIMENSION(:,:) :: BUFFER
CHARACTER(40) :: NAME_SCARC(0:N_TIMERS_SCARC)
REAL(EB) :: TPCNT_SCARC(0:N_TIMERS_SCARC)

!!! only temporarily - use routine only in debug mode
IF (TYPE_DEBUG == NSCARC_DEBUG_NONE) RETURN

IERR=0
 
ALLOCATE(COUNTS_SCARC_TIMERS(0:NUMPROCS-1))
ALLOCATE(DISPLS_SCARC_TIMERS(0:NUMPROCS-1))
ALLOCATE(BUFFER(N_TIMERS_SCARC, NMESHES))

COUNTS_SCARC_TIMERS = COUNTS_SCARC*N_TIMERS_SCARC
DISPLS_SCARC_TIMERS = DISPLS_SCARC*N_TIMERS_SCARC

BUFFER = TUSED_SCARC
IF (USE_MPI) CALL MPI_GATHERV(TUSED_SCARC(1,DISPLS_SCARC(MYID)+1),COUNTS_SCARC_TIMERS(MYID),&
                              MPI_DOUBLE_PRECISION,TUSED_SCARC, COUNTS_SCARC_TIMERS,DISPLS_SCARC_TIMERS,&
                              MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)

!!! Printout subroutine timings
!!! outdated version, must be revised (not used at the moment)

IF (MYID==0) THEN
   NAME_SCARC                           = 'null'
   NAME_SCARC(NSCARC_TIME_TOTAL)        = 'ScaRC complete'
   NAME_SCARC(NSCARC_TIME_SETUP)        = 'ScaRC setup'
   NAME_SCARC(NSCARC_TIME_SOLVER)       = 'ScaRC solver'
   NAME_SCARC(NSCARC_TIME_KRYLOV)       = 'ScaRC Krylov method'
   NAME_SCARC(NSCARC_TIME_MULTIGRID)    = 'ScaRC multigrid method'
   NAME_SCARC(NSCARC_TIME_PRECON)       = 'ScaRC preconditioner'
   NAME_SCARC(NSCARC_TIME_SMOOTH)       = 'ScaRC smoother'
   NAME_SCARC(NSCARC_TIME_COARSE)       = 'ScaRC coarse grid solver'
   NAME_SCARC(NSCARC_TIME_MATVEC)       = 'ScaRC matrix-vector product'
   NAME_SCARC(NSCARC_TIME_SCALPROD)     = 'ScaRC scalar product'
   NAME_SCARC(NSCARC_TIME_L2NORM)       = 'ScaRC L2-norm'
   NAME_SCARC(NSCARC_TIME_EXCH_INIT)    = 'ScaRC initialization data exchange'
   NAME_SCARC(NSCARC_TIME_EXCH_VECTOR)  = 'ScaRC exchange of vector values'
   NAME_SCARC(NSCARC_TIME_EXCH_MATRIX)  = 'ScaRC exchange of matrix values'
   NAME_SCARC(NSCARC_TIME_EXCH_MEASURE) = 'ScaRC exchange of measure values'
   
   DO NM=1,NMESHES
      DO I=1,N_TIMERS_SCARC
         TPCNT_SCARC(I) = 100._EB*TUSED_SCARC(I,NM)/TUSED(1,NM)
      ENDDO
      WRITE(LU_OUTPUT,443) NM
      !WRITE(LU_OUTPUT,444) 'MAIN',TUSED(1,NM),100._EB
      WRITE(LU_OUTPUT,444) (NAME_SCARC(I),TUSED_SCARC(I,NM),TPCNT_SCARC(I),I=1,N_TIMERS_SCARC)
   ENDDO
ENDIF

   IF (TYPE_DEBUG /= NSCARC_DEBUG_NONE ) CLOSE(SCARC_LU)

443 FORMAT(//' ScaRC: CPU Time Usage, Mesh ',I3// &
         47X,' CPU (s)        %  '/ &
         7X,' -----------------------------------------------------------------')
444 FORMAT(7X,A40,2F11.2)

END SUBROUTINE SCARC_TIMINGS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Stops the code gracefully after writing a message
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SCARC_SHUTDOWN(CMESSAGE)

CHARACTER(*) CMESSAGE
LOGICAL :: EX

WRITE(LU_ERR,'(/A)') TRIM(CMESSAGE)
INQUIRE(FILE=FN_OUTPUT,EXIST=EX)
IF (.NOT.EX) OPEN(LU_OUTPUT,FILE=TRIM(CHID)//'.out',STATUS='REPLACE',FORM='FORMATTED')
WRITE(LU_OUTPUT,'(/A)') TRIM(CMESSAGE)
STOP

END SUBROUTINE SCARC_SHUTDOWN


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Compute current revision number
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GET_REV_SCRC(MODULE_REV,MODULE_DATE)
INTEGER,INTENT(INOUT) :: MODULE_REV
CHARACTER(255),INTENT(INOUT) :: MODULE_DATE

WRITE(MODULE_DATE,'(A)') SCRCREV(INDEX(SCRCREV,':')+1:LEN_TRIM(SCRCREV)-2)
READ (MODULE_DATE,'(I5)') MODULE_REV
WRITE(MODULE_DATE,'(A)') SCRCDATE

END SUBROUTINE GET_REV_SCRC

END MODULE SCRC
 

