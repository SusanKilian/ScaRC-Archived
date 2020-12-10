MODULE SCARC_UTILITIES
  
USE GLOBAL_CONSTANTS
USE MESH_VARIABLES, ONLY: MESHES
USE SCARC_CONSTANTS
USE SCARC_VARIABLES

CONTAINS

! -----------------------------------------------------------------------------
!> \brief Check if two meshes are neighbors
! -----------------------------------------------------------------------------
LOGICAL FUNCTION ARE_NEIGHBORS(NM, NOM)
USE SCARC_POINTERS, ONLY: OM
INTEGER, INTENT(IN) :: NM, NOM

ARE_NEIGHBORS = .TRUE.

OM => MESHES(NM)%OMESH(NOM)
IF (OM%NIC_R == 0 .AND. OM%NIC_S == 0) ARE_NEIGHBORS = .FALSE.

END FUNCTION ARE_NEIGHBORS


! ------------------------------------------------------------------------------------------------
!> \brief Get full text information about the data type of the currently processed array
! ------------------------------------------------------------------------------------------------
CHARACTER(10) FUNCTION SCARC_GET_DATA_TYPE(NDATA)
INTEGER, INTENT(IN) :: NDATA

SELECT CASE (NDATA)
   CASE (NSCARC_DATA_INTEGER)
      SCARC_GET_DATA_TYPE = 'INTEGER'
   CASE (NSCARC_DATA_REAL_EB)
      SCARC_GET_DATA_TYPE = 'REAL_EB'
   CASE (NSCARC_DATA_REAL_FB)
      SCARC_GET_DATA_TYPE = 'REAL_FB'
   CASE (NSCARC_DATA_LOGICAL)
      SCARC_GET_DATA_TYPE = 'LOGICAL'
   CASE (NSCARC_DATA_CMATRIX)
      SCARC_GET_DATA_TYPE = 'CMATRIX'
   CASE (NSCARC_DATA_BMATRIX)
      SCARC_GET_DATA_TYPE = 'BMATRIX'
   CASE DEFAULT
      SCARC_GET_DATA_TYPE = ' '
END SELECT

END FUNCTION SCARC_GET_DATA_TYPE


! ------------------------------------------------------------------------------------------------
!> \brief Get full text information about the dimension of the currently processed array
! ------------------------------------------------------------------------------------------------
CHARACTER(10) FUNCTION SCARC_GET_DIMENSION(NDIM)
INTEGER, INTENT(IN) :: NDIM

SELECT CASE (NDIM)
   CASE (1)
      SCARC_GET_DIMENSION = '1'
   CASE (2)
      SCARC_GET_DIMENSION = '2'
   CASE (3)
      SCARC_GET_DIMENSION = '3'
   CASE DEFAULT
      SCARC_GET_DIMENSION = '0'
END SELECT

END FUNCTION SCARC_GET_DIMENSION

! ------------------------------------------------------------------------------------------------
!> \brief Get full text information about the initialization type of the currently processed array
! ------------------------------------------------------------------------------------------------
CHARACTER(10) FUNCTION SCARC_GET_INIT_TYPE(NINIT, NDATA, NSTATE)
INTEGER, INTENT(IN) :: NINIT, NDATA, NSTATE

SELECT CASE (NINIT)
   CASE (NSCARC_INIT_UNDEF)
      SCARC_GET_INIT_TYPE = 'UNDEF'
   CASE (NSCARC_INIT_NONE)
      SCARC_GET_INIT_TYPE = 'NONE'
   CASE (NSCARC_INIT_MINUS)
      SCARC_GET_INIT_TYPE = 'MINUS'
   CASE (NSCARC_INIT_ZERO)
      SCARC_GET_INIT_TYPE = 'ZERO'
   CASE (NSCARC_INIT_ONE)
      SCARC_GET_INIT_TYPE = 'ONE'
   CASE (NSCARC_INIT_TRUE)
      SCARC_GET_INIT_TYPE = 'TRUE'
   CASE (NSCARC_INIT_FALSE)
      SCARC_GET_INIT_TYPE = 'FALSE'
   CASE (NSCARC_INIT_HUGE)
      SCARC_GET_INIT_TYPE = 'HUGE'
   CASE DEFAULT
      SCARC_GET_INIT_TYPE = ' '
END SELECT

IF (NDATA == NSCARC_DATA_CMATRIX .OR. NDATA == NSCARC_DATA_BMATRIX) SCARC_GET_INIT_TYPE = ' '
IF (NSTATE == NSCARC_MEMORY_REMOVE) SCARC_GET_INIT_TYPE = ' '

END FUNCTION SCARC_GET_INIT_TYPE

! ------------------------------------------------------------------------------------------------
!> \brief Get type of matrix storage scheme for specified grid level
! ------------------------------------------------------------------------------------------------
INTEGER FUNCTION SCARC_GET_MATRIX_TYPE(NL)
INTEGER, INTENT(IN) :: NL

IF (NL == NLEVEL_MAX .AND. TYPE_COARSE == NSCARC_COARSE_DIRECT) THEN
   SCARC_GET_MATRIX_TYPE = NSCARC_MATRIX_COMPACT
ELSE
   SCARC_GET_MATRIX_TYPE = TYPE_MATRIX
ENDIF

END FUNCTION SCARC_GET_MATRIX_TYPE

! ------------------------------------------------------------------------------------------------
!> \brief Check if a subdiagonal entry must be computed in a specified coordinate direction
! If a structured discretization is used, then subdiagonals are built in every direction
! Else check the type of the neighboring cell in direction IOR0
! ------------------------------------------------------------------------------------------------
LOGICAL FUNCTION IS_VALID_DIRECTION(IX, IY, IZ, IOR0)
USE SCARC_POINTERS, ONLY: L, G
INTEGER, INTENT(IN)  :: IX, IY, IZ, IOR0
INTEGER :: IC_INDEX, IW_INDEX

IS_VALID_DIRECTION = .FALSE.
IF (TWO_D .AND. ABS(IOR0) == 2) RETURN

SELECT CASE (TYPE_GRID)
   CASE (NSCARC_GRID_STRUCTURED)
      IS_VALID_DIRECTION = .TRUE.                                                ! always build subdiagonals
      RETURN
   CASE (NSCARC_GRID_UNSTRUCTURED)
      IC_INDEX = L%CELL_INDEX_PTR(IX, IY, IZ)                               ! cell index of corresponding cell
      IW_INDEX = 0
      IF (IC_INDEX /= 0) IW_INDEX  = L%WALL_INDEX_PTR(IC_INDEX, -IOR0)      ! check its wall index

      IF (IW_INDEX == 0) THEN                                               ! if zero, build subdiagonal
         IS_VALID_DIRECTION = .TRUE.
         RETURN
      ELSE                                                                  ! if not, only build along interfaces
         IF (TYPE_SCOPE(0) == NSCARC_SCOPE_GLOBAL .AND. G%WALL(IW_INDEX)%BOUNDARY_TYPE== INTERPOLATED_BOUNDARY) THEN
            IS_VALID_DIRECTION = .TRUE.
            RETURN
         ENDIF
      ENDIF
END SELECT
RETURN

END FUNCTION IS_VALID_DIRECTION


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








END MODULE SCARC_UTILITIES
