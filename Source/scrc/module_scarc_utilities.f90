MODULE SCARC_UTILITIES
  
USE MESH_VARIABLES, ONLY: MESHES
USE SCARC_CONSTANTS

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



END MODULE SCARC_UTILITIES
