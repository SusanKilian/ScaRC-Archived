MODULE SCARC_POINTER_ROUTINES
  
USE GLOBAL_CONSTANTS
USE PRECISION_PARAMETERS, ONLY: EB, FB
USE MEMORY_FUNCTIONS, ONLY: CHKMEMERR
USE COMP_FUNCTIONS, ONLY: GET_FILE_NUMBER
USE MESH_VARIABLES, ONLY: MESHES
USE SCARC_CONSTANTS
USE SCARC_VARIABLES
USE SCARC_TYPES

CONTAINS

! -----------------------------------------------------------------------------
!> \brief Point to specified mesh
! -----------------------------------------------------------------------------
SUBROUTINE SCARC_POINT_TO_MESH(NM)
USE SCARC_POINTERS, ONLY: M, S
INTEGER, INTENT(IN) :: NM
M => MESHES(NM)
S => SCARC(NM)
END SUBROUTINE SCARC_POINT_TO_MESH


! -----------------------------------------------------------------------------
!> \brief Point to specified combination of mesh and grid level
! -----------------------------------------------------------------------------
SUBROUTINE SCARC_POINT_TO_LEVEL(NM, NL)
USE SCARC_POINTERS, ONLY: M, S, L
INTEGER, INTENT(IN) :: NM, NL
M => MESHES(NM)
S => SCARC(NM)
L => S%LEVEL(NL)
END SUBROUTINE SCARC_POINT_TO_LEVEL


! ----------------------------------------------------------------------------------------------------
!> \brief Unset ScaRC pointers  
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
!> \brief Point to specified combination of a mesh level and discretization type
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
      G%NW = L%N_WALL_CELLS_EXT                     ! TODO: set it elsewhere
   CASE (NSCARC_GRID_UNSTRUCTURED)
      G => L%UNSTRUCTURED
      G%NW = L%N_WALL_CELLS_EXT + L%N_WALL_CELLS_INT
END SELECT
W => G%WALL

END SUBROUTINE SCARC_POINT_TO_GRID


! ----------------------------------------------------------------------------------------------------
!> \brief Point to specified pairing of mesh levels and discretization types
! ----------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_POINT_TO_MULTIGRID(NM, NL1, NL2)
USE SCARC_POINTERS, ONLY: M, S, LF, LC, GF, GC, WC, WF
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
WC => GC%WALL
WF => GF%WALL

END SUBROUTINE SCARC_POINT_TO_MULTIGRID


! ---------------------------------------------------------------------------------------------------
!> \brief Point to specified combination of a neighboring mesh level and discretization type
! ---------------------------------------------------------------------------------------------------
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


! -----------------------------------------------------------------------------------------------------
!> \brief Point to specified combination of a neighboring mesh level and a discretization type 
! -----------------------------------------------------------------------------------------------------
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
!> \brief Point to specified matrix in compact storage technique 
! ----------------------------------------------------------------------------------------------------
FUNCTION SCARC_POINT_TO_CMATRIX(G, NTYPE)
TYPE(SCARC_CMATRIX_TYPE), POINTER :: SCARC_POINT_TO_CMATRIX
TYPE(SCARC_GRID_TYPE), POINTER, INTENT(IN) :: G
INTEGER, INTENT(IN) :: NTYPE

SELECT CASE(NTYPE)
   CASE (NSCARC_MATRIX_POISSON_PROL)
      SCARC_POINT_TO_CMATRIX => G%POISSON_PROL
   CASE (NSCARC_MATRIX_CONNECTION)
      SCARC_POINT_TO_CMATRIX => G%CONNECTION
   CASE (NSCARC_MATRIX_POISSON)
      SCARC_POINT_TO_CMATRIX => G%POISSON
#ifdef WITH_MKL
   CASE (NSCARC_MATRIX_POISSON_SYM)
      SCARC_POINT_TO_CMATRIX => G%POISSON_SYM
#endif
   CASE (NSCARC_MATRIX_LAPLACE)
      SCARC_POINT_TO_CMATRIX => G%LAPLACE
   CASE (NSCARC_MATRIX_PROLONGATION)
      SCARC_POINT_TO_CMATRIX => G%PROLONGATION
   CASE (NSCARC_MATRIX_RESTRICTION)
      SCARC_POINT_TO_CMATRIX => G%RESTRICTION
   CASE (NSCARC_MATRIX_ZONES)
      SCARC_POINT_TO_CMATRIX => G%ZONES
   CASE (NSCARC_MATRIX_LM)
      SCARC_POINT_TO_CMATRIX => G%LM
   CASE (NSCARC_MATRIX_UM)
      SCARC_POINT_TO_CMATRIX => G%UM
END SELECT

END FUNCTION SCARC_POINT_TO_CMATRIX


! ----------------------------------------------------------------------------------------------------
!> \brief Point to specified matrix in bandwise storage technique 
! ----------------------------------------------------------------------------------------------------
FUNCTION SCARC_POINT_TO_BMATRIX(G, NTYPE)
TYPE(SCARC_BMATRIX_TYPE), POINTER :: SCARC_POINT_TO_BMATRIX
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
!> \brief Point to specified neighboring matrix in compact storage technique 
! ----------------------------------------------------------------------------------------------------
FUNCTION SCARC_POINT_TO_OTHER_CMATRIX(OG, NTYPE)
TYPE(SCARC_CMATRIX_TYPE), POINTER :: SCARC_POINT_TO_OTHER_CMATRIX
TYPE(SCARC_GRID_TYPE), POINTER, INTENT(IN) :: OG
INTEGER, INTENT(IN) :: NTYPE

SELECT CASE(NTYPE)
   CASE (NSCARC_MATRIX_POISSON_PROL)
      SCARC_POINT_TO_OTHER_CMATRIX => OG%POISSON_PROL
   CASE (NSCARC_MATRIX_CONNECTION)
      SCARC_POINT_TO_OTHER_CMATRIX => OG%CONNECTION
   CASE (NSCARC_MATRIX_POISSON)
      SCARC_POINT_TO_OTHER_CMATRIX => OG%POISSON
#ifdef WITH_MKL
   CASE (NSCARC_MATRIX_POISSON_SYM)
      SCARC_POINT_TO_OTHER_CMATRIX => OG%POISSON_SYM
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
!> \brief Point to specified neighboring matrix in bandwise storage technique 
! ----------------------------------------------------------------------------------------------------
FUNCTION SCARC_POINT_TO_OTHER_BMATRIX(OG, NTYPE)
TYPE(SCARC_BMATRIX_TYPE), POINTER :: SCARC_POINT_TO_OTHER_BMATRIX
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
!> \brief Point to specified integer receive buffer for data exchanges
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
!> \brief Point to specified integer receive buffer for data exchanges
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
#ifdef WITH_SCARC_DEBUG2
WRITE(MSG%LU_DEBUG,*) 'POINT_TO BUFFER_REAL:A:'
WRITE(MSG%LU_DEBUG,'(8E14.6)') SCARC_POINT_TO_BUFFER_REAL(1:16)
#endif
      ELSE
         SCARC_POINT_TO_BUFFER_REAL => SCARC(NOM)%OSCARC(NM)%SEND_BUFFER_REAL
#ifdef WITH_SCARC_DEBUG2
WRITE(MSG%LU_DEBUG,*) 'POINT_TO BUFFER_REAL:B:'
WRITE(MSG%LU_DEBUG,'(8E14.6)') SCARC_POINT_TO_BUFFER_REAL(1:16)
#endif
      ENDIF
END SELECT

END FUNCTION SCARC_POINT_TO_BUFFER_REAL


! ------------------------------------------------------------------------------------------------
!> \brief Point to specified vector on a given grid level
! ------------------------------------------------------------------------------------------------
FUNCTION SCARC_POINT_TO_VECTOR (NM, NL, NV)
REAL(EB), POINTER, DIMENSION(:) :: SCARC_POINT_TO_VECTOR
INTEGER, INTENT(IN) :: NM, NL, NV

SELECT CASE (NV)

   ! Stage one vectors (for methods on first hierarchical level)
 
   CASE (NSCARC_VECTOR_ONE_X)
      SCARC_POINT_TO_VECTOR => SCARC(NM)%LEVEL(NL)%STAGE(NSCARC_STAGE_ONE)%X
   CASE (NSCARC_VECTOR_ONE_B)
      SCARC_POINT_TO_VECTOR => SCARC(NM)%LEVEL(NL)%STAGE(NSCARC_STAGE_ONE)%B
   CASE (NSCARC_VECTOR_ONE_D)
      SCARC_POINT_TO_VECTOR => SCARC(NM)%LEVEL(NL)%STAGE(NSCARC_STAGE_ONE)%D
   CASE (NSCARC_VECTOR_ONE_R)
      SCARC_POINT_TO_VECTOR => SCARC(NM)%LEVEL(NL)%STAGE(NSCARC_STAGE_ONE)%R
   CASE (NSCARC_VECTOR_ONE_V)
      SCARC_POINT_TO_VECTOR => SCARC(NM)%LEVEL(NL)%STAGE(NSCARC_STAGE_ONE)%V
   CASE (NSCARC_VECTOR_ONE_Y)
      SCARC_POINT_TO_VECTOR => SCARC(NM)%LEVEL(NL)%STAGE(NSCARC_STAGE_ONE)%Y
   CASE (NSCARC_VECTOR_ONE_Z)
      SCARC_POINT_TO_VECTOR => SCARC(NM)%LEVEL(NL)%STAGE(NSCARC_STAGE_ONE)%Z
#ifdef WITH_SCARC_DEBUG
   CASE (NSCARC_VECTOR_ONE_E)
      SCARC_POINT_TO_VECTOR => SCARC(NM)%LEVEL(NL)%STAGE(NSCARC_STAGE_ONE)%E
#endif

 
   ! Stage two vectors (for methods on second hierarchical level)
 
   CASE (NSCARC_VECTOR_TWO_X)
      SCARC_POINT_TO_VECTOR => SCARC(NM)%LEVEL(NL)%STAGE(NSCARC_STAGE_TWO)%X
   CASE (NSCARC_VECTOR_TWO_B)
      SCARC_POINT_TO_VECTOR => SCARC(NM)%LEVEL(NL)%STAGE(NSCARC_STAGE_TWO)%B
   CASE (NSCARC_VECTOR_TWO_D)
      SCARC_POINT_TO_VECTOR => SCARC(NM)%LEVEL(NL)%STAGE(NSCARC_STAGE_TWO)%D
   CASE (NSCARC_VECTOR_TWO_R)
      SCARC_POINT_TO_VECTOR => SCARC(NM)%LEVEL(NL)%STAGE(NSCARC_STAGE_TWO)%R
   CASE (NSCARC_VECTOR_TWO_V)
      SCARC_POINT_TO_VECTOR => SCARC(NM)%LEVEL(NL)%STAGE(NSCARC_STAGE_TWO)%V
   CASE (NSCARC_VECTOR_TWO_Y)
      SCARC_POINT_TO_VECTOR => SCARC(NM)%LEVEL(NL)%STAGE(NSCARC_STAGE_TWO)%Y
   CASE (NSCARC_VECTOR_TWO_Z)
      SCARC_POINT_TO_VECTOR => SCARC(NM)%LEVEL(NL)%STAGE(NSCARC_STAGE_TWO)%Z
#ifdef WITH_SCARC_DEBUG
   CASE (NSCARC_VECTOR_TWO_E)
      SCARC_POINT_TO_VECTOR => SCARC(NM)%LEVEL(NL)%STAGE(NSCARC_STAGE_TWO)%E
#endif
END SELECT

END FUNCTION SCARC_POINT_TO_VECTOR


! ------------------------------------------------------------------------------------------------
!> \brief Point to specified vector on a given grid level (single precision version)
! ------------------------------------------------------------------------------------------------
FUNCTION SCARC_POINT_TO_VECTOR_FB(NM, NL, NV)
REAL(FB), POINTER, DIMENSION(:) :: SCARC_POINT_TO_VECTOR_FB
INTEGER, INTENT(IN) :: NM, NL, NV


SELECT CASE (NV)

   ! Stage one vectors (for methods on first hierarchical level)
 
   CASE (NSCARC_VECTOR_ONE_X)
      SCARC_POINT_TO_VECTOR_FB => SCARC(NM)%LEVEL(NL)%STAGE(NSCARC_STAGE_ONE)%X_FB
   CASE (NSCARC_VECTOR_ONE_B)
      SCARC_POINT_TO_VECTOR_FB => SCARC(NM)%LEVEL(NL)%STAGE(NSCARC_STAGE_ONE)%B_FB
   CASE (NSCARC_VECTOR_ONE_R)
      SCARC_POINT_TO_VECTOR_FB => SCARC(NM)%LEVEL(NL)%STAGE(NSCARC_STAGE_ONE)%R_FB
   CASE (NSCARC_VECTOR_ONE_V)
      SCARC_POINT_TO_VECTOR_FB => SCARC(NM)%LEVEL(NL)%STAGE(NSCARC_STAGE_ONE)%V_FB

 
   ! Stage two vectors (for methods on second hierarchical level)
 
   CASE (NSCARC_VECTOR_TWO_X)
      SCARC_POINT_TO_VECTOR_FB => SCARC(NM)%LEVEL(NL)%STAGE(NSCARC_STAGE_TWO)%X_FB
   CASE (NSCARC_VECTOR_TWO_B)
      SCARC_POINT_TO_VECTOR_FB => SCARC(NM)%LEVEL(NL)%STAGE(NSCARC_STAGE_TWO)%B_FB
   CASE (NSCARC_VECTOR_TWO_R)
      SCARC_POINT_TO_VECTOR_FB => SCARC(NM)%LEVEL(NL)%STAGE(NSCARC_STAGE_TWO)%R_FB
   CASE (NSCARC_VECTOR_TWO_V)
      SCARC_POINT_TO_VECTOR_FB => SCARC(NM)%LEVEL(NL)%STAGE(NSCARC_STAGE_TWO)%V_FB
END SELECT

END FUNCTION SCARC_POINT_TO_VECTOR_FB


! ------------------------------------------------------------------------------------------------
!> \brief Point to pressure vector in predictor or corrector
! ------------------------------------------------------------------------------------------------
FUNCTION SCARC_POINT_TO_HVECTOR(NM, NV)
REAL(EB), POINTER, DIMENSION(:,:,:) :: SCARC_POINT_TO_HVECTOR
INTEGER, INTENT(IN) :: NM, NV
SELECT CASE (NV)
   CASE (NSCARC_VECTOR_H)
      SCARC_POINT_TO_HVECTOR => MESHES(NM)%H
   CASE (NSCARC_VECTOR_HS)
      SCARC_POINT_TO_HVECTOR => MESHES(NM)%HS
END SELECT
END FUNCTION SCARC_POINT_TO_HVECTOR

END MODULE SCARC_POINTER_ROUTINES
