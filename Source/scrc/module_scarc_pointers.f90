! ================================================================================================================
!  MODULE 'SCARC_POINTERS'  
!> \brief Collection of different pointers to specify the different meshes, grid levels, discretizations and matrices
! ================================================================================================================
MODULE SCARC_POINTERS

USE MESH_VARIABLES
USE SCARC_MEMORY_MANAGER
USE SCARC_STACK_ADMINISTRATION
USE SCARC_DISRETIZATION
USE SCARC_METHODS

IMPLICIT NONE

TYPE (MESH_TYPE), POINTER :: M=>NULL()                  !< Pointer to specified mesh (based on MESHES from base code)
TYPE (OMESH_TYPE), POINTER :: OM=>NULL()                !< Pointer to specified neighboring mesh (based on OMESH from base code)
TYPE (WALL_TYPE), POINTER :: MWC=>NULL()                !< Pointer to specified wall cell (based on WALL from base code)

TYPE (EXTERNAL_WALL_TYPE), POINTER :: EWC=>NULL()       !< Pointer to specified external wall cell

TYPE (SCARC_ALLOCATION_TYPE), POINTER :: AL=>NULL()     !< Pointer to allocated structure within ScaRC memory management

TYPE (SCARC_TYPE), POINTER :: S=>NULL()                 !< Pointer to ScaRC-structure on a specified mesh
TYPE (SCARC_OSCARC_TYPE), POINTER :: OS=>NULL()         !< Pointer to ScaRC-structure on a specified neighboring

TYPE (SCARC_LEVEL_TYPE), POINTER :: L=>NULL()           !< Pointer to specified level
TYPE (SCARC_LEVEL_TYPE), POINTER :: LF=>NULL()          !< Pointer to specified fine level
TYPE (SCARC_LEVEL_TYPE), POINTER :: LC=>NULL()          !< Pointer to specified coarse level
TYPE (SCARC_LEVEL_TYPE), POINTER :: OL=>NULL()          !< Pointer to specified level on neighboring mesh
TYPE (SCARC_LEVEL_TYPE), POINTER :: OLF=>NULL()         !< Pointer to specified fine level on neighboring mesh
TYPE (SCARC_LEVEL_TYPE), POINTER :: OLC=>NULL()         !< Pointer to specified coarse level on neighboring mesh

TYPE (SCARC_GRID_TYPE), POINTER :: G=>NULL()            !< Pointer to specified grid discretization 
TYPE (SCARC_GRID_TYPE), POINTER :: GF=>NULL()           !< Pointer to specified fine grid discretization 
TYPE (SCARC_GRID_TYPE), POINTER :: GC=>NULL()           !< Pointer to specified coarse grid discretization 
TYPE (SCARC_GRID_TYPE), POINTER :: OG=>NULL()           !< Pointer to specified grid discretization 
TYPE (SCARC_GRID_TYPE), POINTER :: OGF=>NULL()          !< Pointer to specified fine grid discretization on neighboring mesh
TYPE (SCARC_GRID_TYPE), POINTER :: OGC=>NULL()          !< Pointer to specified coarse grid discretization on neighboring mesh

TYPE (SCARC_FACE_TYPE), POINTER :: F=>NULL()            !< Pointer to specified face of grid
TYPE (SCARC_FACE_TYPE), POINTER :: FF=>NULL()           !< Pointer to specified face of fine grid level
TYPE (SCARC_FACE_TYPE), POINTER :: FC=>NULL()           !< Pointer to specified face of coarse grid level

TYPE (SCARC_OBST_TYPE), POINTER :: OB=>NULL()           !< Pointer to specified obstruction
TYPE (SCARC_OBST_TYPE), POINTER :: OBF=>NULL()          !< Pointer to specified obstruction of fine grid level
TYPE (SCARC_OBST_TYPE), POINTER :: OBC=>NULL()          !< Pointer to specified obstruction of coarse grid level

TYPE (SCARC_WALL_TYPE), POINTER :: GWC=>NULL()          !< Pointer to specified wall cell

TYPE (SCARC_WALL_TYPE), DIMENSION(:), POINTER :: W=>NULL()    !< Pointer to complete wall structure
TYPE (SCARC_WALL_TYPE), DIMENSION(:), POINTER :: WF=>NULL()   !< Pointer to wall structure on fine grid level
TYPE (SCARC_WALL_TYPE), DIMENSION(:), POINTER :: WC=>NULL()   !< Pointer to wall structure on coarse grid level

TYPE (SCARC_SUBDIVISION_TYPE), POINTER :: SUB=>NULL()         !< Pointer to subdivision structure (only shortcut)

TYPE (SCARC_SOLVER_TYPE), POINTER :: SV=>NULL()               !< Pointer to ScaRC solver structure
TYPE (SCARC_SOLVER_TYPE), POINTER :: SVP=>NULL()              !< Pointer to parent ScaRC solver structure

TYPE (SCARC_STAGE_TYPE),  POINTER :: ST=>NULL()               !< Pointer to solver stage structure
TYPE (SCARC_STAGE_TYPE),  POINTER :: STP=>NULL()              !< Pointer to parent solver stage structure

TYPE (SCARC_FFT_TYPE), POINTER :: FFT=>NULL()                 !< Pointer to FFT structure
TYPE (SCARC_MGM_TYPE), POINTER :: MGM=>NULL()                 !< Pointer to McKeeney-Greengard-Mayo structure

TYPE (SCARC_BMATRIX_TYPE), POINTER :: AB=>NULL()       !< Pointer to bandwise matrix structure
TYPE (SCARC_BMATRIX_TYPE), POINTER :: OAB=>NULL()      !< Pointer to neighboring bandwise matrix structure

TYPE (SCARC_CMATRIX_TYPE), POINTER :: A=>NULL()        !< Pointer to compactly stored matrix 
TYPE (SCARC_CMATRIX_TYPE), POINTER :: A1=>NULL()       !< Pointer to compactly stored matrix 
TYPE (SCARC_CMATRIX_TYPE), POINTER :: A2=>NULL()       !< Pointer to compactly stored matrix 
TYPE (SCARC_CMATRIX_TYPE), POINTER :: AC=>NULL()       !< Pointer to compactly stored coarse matrix 
TYPE (SCARC_CMATRIX_TYPE), POINTER :: AF=>NULL()       !< Pointer to compactly stored fine matrix 
TYPE (SCARC_CMATRIX_TYPE), POINTER :: OA=>NULL()       !< Pointer to compactly stored neighboring matrix 
TYPE (SCARC_CMATRIX_TYPE), POINTER :: OAC=>NULL()      !< Pointer to compactly stored coarse neighboring matrix 
TYPE (SCARC_CMATRIX_TYPE), POINTER :: OAF=>NULL()      !< Pointer to compactly stored fine neighboring matrix 

TYPE (SCARC_CMATRIX_TYPE), POINTER :: LM=>NULL()     !< Pointer to compactly stored matrix 
TYPE (SCARC_CMATRIX_TYPE), POINTER :: UM=>NULL()     !< Pointer to compactly stored matrix 

TYPE (SCARC_CMATRIX_TYPE), POINTER :: P=>NULL()        !< Pointer to compactly stored Prolongation matrix 
TYPE (SCARC_CMATRIX_TYPE), POINTER :: PC=>NULL()       !< Pointer to compactly stored coarse Prolongation matrix 
TYPE (SCARC_CMATRIX_TYPE), POINTER :: PF=>NULL()       !< Pointer to compactly stored fine matrix on coarse grid
TYPE (SCARC_CMATRIX_TYPE), POINTER :: OP=>NULL()       !< Pointer to compactly stored neighboring Prolongation matrix 
TYPE (SCARC_CMATRIX_TYPE), POINTER :: OPC=>NULL()      !< Pointer to compactly stored coarse neighboring Prolongation matrix
TYPE (SCARC_CMATRIX_TYPE), POINTER :: OPF=>NULL()      !< Pointer to compactly stored fine neighboring Prolongation matrix 

TYPE (SCARC_CMATRIX_TYPE), POINTER :: R=>NULL()        !< Pointer to compactly stored Restriction matrix 
TYPE (SCARC_CMATRIX_TYPE), POINTER :: RC=>NULL()       !< Pointer to compactly stored fine Restriction matrix 
TYPE (SCARC_CMATRIX_TYPE), POINTER :: RF=>NULL()       !< Pointer to compactly stored coarse Restriction matrix 
TYPE (SCARC_CMATRIX_TYPE), POINTER :: OR=>NULL()       !< Pointer to compactly stored neighboring Restriction matrix 
TYPE (SCARC_CMATRIX_TYPE), POINTER :: ORC=>NULL()      !< Pointer to compactly stored coarse neighboring Restriction matrix
TYPE (SCARC_CMATRIX_TYPE), POINTER :: ORF=>NULL()      !< Pointer to compactly stored fine neighboring Restriction matrix

TYPE (SCARC_CMATRIX_TYPE), POINTER :: C=>NULL()        !< Pointer to compactly stored connection matrix 
TYPE (SCARC_CMATRIX_TYPE), POINTER :: CC=>NULL()       !< Pointer to compactly stored coarse connection matrix 
TYPE (SCARC_CMATRIX_TYPE), POINTER :: CF=>NULL()       !< Pointer to compactly stored fine connection matrix 
TYPE (SCARC_CMATRIX_TYPE), POINTER :: OC=>NULL()       !< Pointer to compactly stored neighboring connection matrix 
TYPE (SCARC_CMATRIX_TYPE), POINTER :: OCC=>NULL()      !< Pointer to compactly stored coarse neighboring connection matrix
TYPE (SCARC_CMATRIX_TYPE), POINTER :: OCF=>NULL()      !< Pointer to compactly stored fine neighboring connection matrix

TYPE (SCARC_CMATRIX_TYPE), POINTER :: Z=>NULL()        !< Pointer to compactly stored Zones matrix 
TYPE (SCARC_CMATRIX_TYPE), POINTER :: ZC=>NULL()       !< Pointer to compactly stored coarse Zones matrix 
TYPE (SCARC_CMATRIX_TYPE), POINTER :: ZF=>NULL()       !< Pointer to compactly stored fine Zones matrix 
TYPE (SCARC_CMATRIX_TYPE), POINTER :: OZ=>NULL()       !< Pointer to compactly stored neighboring Zones matrix 
TYPE (SCARC_CMATRIX_TYPE), POINTER :: OZC=>NULL()      !< Pointer to compactly stored neighboring coarse Zones matrix
TYPE (SCARC_CMATRIX_TYPE), POINTER :: OZF=>NULL()      !< Pointer to compactly stored neighboring fine Zones matrix 

TYPE (SCARC_CMATRIX_TYPE), POINTER :: PP=>NULL()       !< Pointer to compactly stored Poisson-times-Prolongation matrix
TYPE (SCARC_CMATRIX_TYPE), POINTER :: PPF=>NULL()      !< Pointer to compactly stored fine Poisson-times-Prolongation matrix
TYPE (SCARC_CMATRIX_TYPE), POINTER :: PPC=>NULL()      !< Pointer to compactly stored coarse Poisson-times-Prolongation matrix
TYPE (SCARC_CMATRIX_TYPE), POINTER :: OPP=>NULL()      !< Pointer to compactly stored neighboring Poisson-times-Prolongation matrix

TYPE (SCARC_MATRIX_COMPACT_CONDENSED_TYPE),  POINTER :: ACO =>NULL()    !< Pointer to compactly stored condensed Poisson matrix
TYPE (SCARC_MATRIX_BANDWISE_CONDENSED_TYPE), POINTER :: ABCO=>NULL()    !< Pointer to bandwise stored condensed Poisson matrix

TYPE (SCARC_MULTIGRID_TYPE), POINTER :: MG =>NULL()                     !< Pointer to multigrid type

REAL(EB), POINTER, DIMENSION(:) :: XCOR=>NULL()        !< Pointer to vector of node coordinates in x-direction
REAL(EB), POINTER, DIMENSION(:) :: YCOR=>NULL()        !< Pointer to vector of node coordinates in x-direction
REAL(EB), POINTER, DIMENSION(:) :: ZCOR=>NULL()        !< Pointer to vector of node coordinates in x-direction

REAL(EB), POINTER, DIMENSION(:) :: XMID=>NULL()        !< Pointer to vector of cell midpoints in x-direction
REAL(EB), POINTER, DIMENSION(:) :: YMID=>NULL()        !< Pointer to vector of cell midpoints in y-direction
REAL(EB), POINTER, DIMENSION(:) :: ZMID=>NULL()        !< Pointer to vector of cell midpoints in z-direction

REAL(EB), POINTER, DIMENSION(:) :: VC=>NULL()          !< Pointer to vector on coarse grid level
REAL(EB), POINTER, DIMENSION(:) :: VF=>NULL()          !< Pointer to vector on fine grid level
REAL(EB), POINTER, DIMENSION(:) :: V1=>NULL()          !< Pointer to first vector 
REAL(EB), POINTER, DIMENSION(:) :: V2=>NULL()          !< Pointer to second vector 

REAL(EB), POINTER, DIMENSION(:) :: OH1=>NULL()         !< Pointer to other Poisson solution
REAL(EB), POINTER, DIMENSION(:) :: OH2=>NULL()         !< Pointer to other Laplace solution
REAL(EB), POINTER, DIMENSION(:) :: OH3=>NULL()         !< Pointer to other MGM solution

REAL(EB), POINTER, DIMENSION(:,:,:) :: HP=>NULL()      !< Pointer to pressure vector 
REAL(EB), POINTER, DIMENSION(:,:,:) :: H1=>NULL()      !< Pointer to Poisson solution
REAL(EB), POINTER, DIMENSION(:,:,:) :: H2=>NULL()      !< Pointer to Laplace solution
REAL(EB), POINTER, DIMENSION(:,:,:) :: H3=>NULL()      !< Pointer to MGM solution
REAL(EB), POINTER, DIMENSION(:,:,:) :: PRHS=>NULL()    !< Pointer to right hand side vector
REAL(EB), POINTER, DIMENSION(:,:,:) :: UU=>NULL()      !< Pointer to u-velocity vector
REAL(EB), POINTER, DIMENSION(:,:,:) :: VV=>NULL()      !< Pointer to v-velocity vector
REAL(EB), POINTER, DIMENSION(:,:,:) :: WW=>NULL()      !< Pointer to w-velocity vector

REAL(EB), POINTER, DIMENSION(:) ::  RECV_BUFFER_REAL   !< Pointer to double precision receive vector 
INTEGER,  POINTER, DIMENSION(:) ::  RECV_BUFFER_INT    !< Pointer to inter receive vector 

#ifdef WITH_MKL
TYPE (SCARC_MKL_TYPE), POINTER :: MKL=>NULL()          !< Pointer to MKL type
TYPE (SCARC_CMATRIX_TYPE), POINTER :: AS=>NULL()       !< Pointer to symmetric Poisson matrix
TYPE (SCARC_CMATRIX_TYPE), POINTER :: ACS=>NULL()      !< Pointer to coarse symmetric Poisson matrix
TYPE (SCARC_CMATRIX_TYPE), POINTER :: AFS=>NULL()      !< Pointer to fine symmetric Poisson matrix
TYPE (SCARC_CMATRIX_TYPE), POINTER :: OAS=>NULL()      !< Pointer to neighboring symmetric Poisson matrix
TYPE (SCARC_CMATRIX_TYPE), POINTER :: OACS=>NULL()     !< Pointer to coarse neighboring symmetric Poisson matrix
TYPE (SCARC_CMATRIX_TYPE), POINTER :: OAFS=>NULL()     !< Pointer to fine neighboring symmetric Poisson matrix
REAL(FB), DIMENSION(:), POINTER :: V1_FB=>NULL()       !< Pointer to first single precision vector
REAL(FB), DIMENSION(:), POINTER :: V2_FB=>NULL()       !< Pointer to second single precision vector
#endif

#ifdef WITH_SCARC_POSTPROCESSING
TYPE (SCARC_PRESSURE_TYPE), POINTER :: PR=>NULL()                       !< Pointer to pressure type
#endif

END MODULE SCARC_POINTERS


