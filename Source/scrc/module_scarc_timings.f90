MODULE SCARC_TIMINGS
  
USE GLOBAL_CONSTANTS
USE PRECISION_PARAMETERS, ONLY: EB, FB
USE MEMORY_FUNCTIONS, ONLY: CHKMEMERR
USE MPI
USE SCARC_CONSTANTS
USE SCARC_VARIABLES
USE SCARC_MESSAGE_SERVICES

!> \brief Measurement of CPU times

TYPE SCARC_CPU_TYPE
   REAL(EB) :: BUFFER_PACKING   = 0.0_EB                           !< Time for data exchange
   REAL(EB) :: BUFFER_UNPACKING = 0.0_EB                           !< Time for data exchange
   REAL(EB) :: AMG              = 0.0_EB                           !< Time for algebraic multigrid solver
   REAL(EB) :: COARSE           = 0.0_EB                           !< Time for coarse grid solver
   REAL(EB) :: EXCHANGE         = 0.0_EB                           !< Time for data exchange
   REAL(EB) :: ITERATION        = 0.0_EB                           !< Time for Krylov solver
   REAL(EB) :: L2NORM           = 0.0_EB                           !< Time for l2-norm
   REAL(EB) :: MATVEC_PRODUCT   = 0.0_EB                           !< Time for matrix vector multiplication
   REAL(EB) :: OVERALL          = 0.0_EB                           !< Complete time for ScaRC
   REAL(EB) :: RELAXATION       = 0.0_EB                           !< Time for relaxation
   REAL(EB) :: SCALAR_PRODUCT   = 0.0_EB                           !< Time for scalar product
   REAL(EB) :: SETUP            = 0.0_EB                           !< Time for setup of requested ScaRC solver
   REAL(EB) :: SMOOTHER         = 0.0_EB                           !< Time for smoothing
   REAL(EB) :: SOLVER           = 0.0_EB                           !< Time for solver 
   INTEGER  :: N_TIMER          = 13                               !< Total number of timers
END TYPE SCARC_CPU_TYPE

PUBLIC :: DUMP_TIMERS

CONTAINS

! ------------------------------------------------------------------------------------------------
!> \brief Setup time measurements
! ------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_TIMINGS
ALLOCATE (CPU(0:N_MPI_PROCESSES-1), STAT = IERROR)
CALL CHKMEMERR ('SCARC_SETUP_TIMINGS', 'CPU', IERROR)
END SUBROUTINE SCARC_SETUP_TIMINGS


! ----------------------------------------------------------------------------------------------------
!> \brief Dump CPU times of several routines
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

END MODULE SCARC_TIMINGS
