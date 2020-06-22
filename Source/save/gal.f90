! ------------------------------------------------------------------------------------------------------
! Determine coarse grid matrix ACC which is computed as Galerkin operator R*A*P
! ------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_GALERKIN_OPERATOR(NL)
USE SCARC_POINTERS, ONLY: GF, RF, ACC
INTEGER, INTENT(IN) :: NL
INTEGER  :: NM, IPC, IRC, IC, IR, IA, ICOL
REAL(EB) :: DRSUM, DAP, TOL = 1.0E-12_EB

DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   ACC%N_ROW = GF%N_COARSE+1
   ACC%N_VAL = ACC%N_ROW**2             ! only temporarily TODO TOO BIG
   CALL SCARC_ALLOCATE_MATRIX_COMPACT(ACC, NL+1, NSCARC_PRECISION_DOUBLE, NSCARC_INIT_ZERO, 'ACP')
   IF (NL == NLEVEL_MAX - 1 .AND. .NOT.ALLOCATED(ACC%COL_GLOBAL)) &
      CALL SCARC_ALLOCATE_INT1 (ACC%COL_GLOBAL, 1, ACC%N_VAL, NSCARC_INIT_ZERO, 'ACC%COL_GLOBAL')

   IA = 1
   ACC%ROW(1) = 1

   DO IRC = 1, GF%N_COARSE                      ! rows of restriction matrix

      ! First process diagonal entry
      IPC = IRC

      ! Compute entry ACC(ICC,JCC) = SUM( P(:,JCC)*A(IC,:)*R(ICC,IC) )
      DRSUM = 0.0_EB
      DO IR = RF%ROW(IRC), RF%ROW(IRC+1) - 1
         IC = RF%COL(IR)
         DAP = SCARC_VALUE_AP(IC, IPC)
         DRSUM = DRSUM + DAP*RF%VAL(IR)
#ifdef WITH_SCARC_DEBUG2
         WRITE(MSG%LU_DEBUG,'(A,4I4,3E12.4)') 'GAL: IRC, IPC, IR, IC, VAL_AP, RF%VAL, DRSUM:', &
                                               IRC, IPC, IR, IC, DAP, RF%VAL(IR), DRSUM
#endif
      ENDDO
      IF (ABS(DRSUM) > TOL) THEN
         ACC%VAL(IA) = DRSUM
         ACC%COL(IA) = IPC
#ifdef WITH_MKL
IF (IS_MKL_LEVEL(NL+1)) ACC%COL_GLOBAL(IA) = ACC%COL(IA) + GC%NC_OFFSET(NM)
#endif
         IA = IA + 1
#ifdef WITH_SCARC_DEBUG2
         WRITE(MSG%LU_DEBUG,'(A,6I5,E12.4)') 'GALERKIN: IPC, IRC, IR, IC, IA, DRSUM:',&
                                                 IPC,IRC,IR,IC,IA,ACC%COL(IA),ACC%VAL(IA)
#endif
      ENDIF

      ! Then process remaining entries in row
      DO IPC = 1, GF%N_COARSE                   ! columns of prolongation matrix

         IF (IPC == IRC) CYCLE

         ! Compute entry ACC(ICC,JCC) = SUM( P(:,JCC)*A(IC,:)*R(ICC,IC) )
         DRSUM = 0.0_EB
         DO IR = RF%ROW(IRC), RF%ROW(IRC+1) - 1
            IC = RF%COL(IR)
            DAP = SCARC_VALUE_AP(IC, IPC)
            DRSUM = DRSUM + DAP*RF%VAL(IR)
#ifdef WITH_SCARC_DEBUG2
            WRITE(MSG%LU_DEBUG,'(A,4I4,3E12.4)') 'GAL: IRC, IPC, IR, IC, VAL_AP, RF%VAL, DRSUM:', &
                                                  IRC, IPC, IR, IC, DAP, RF%VAL(IR), DRSUM
#endif
         ENDDO
         IF (ABS(DRSUM) > TOL) THEN
            ACC%VAL(IA) = DRSUM
            ACC%COL(IA) = IPC
            IA = IA + 1
#ifdef WITH_SCARC_DEBUG2
            WRITE(MSG%LU_DEBUG,'(A,6I5,E12.4)') 'GALERKIN: IPC, IRC, IR, IC, IA, DRSUM:',&
                                                 IPC,IRC,IR,IC,IA,ACC%COL(IA),ACC%VAL(IA)
#endif
         ENDIF

      ENDDO

      ACC%ROW(IRC+1) = IA
   ENDDO

   ! temporarily
   ACC%COL_GLOBAL = ACC%COL
   ACC%N_STENCIL = ACF%N_STENCIL

   GC%NC =  ACC%N_ROW - 1

#ifdef WITH_SCARC_DEBUG
   WRITE(MSG%LU_DEBUG,*) 'GALERKIN: ACC%ROW'
   WRITE(MSG%LU_DEBUG,*) 'ACC%ROW:', ACC%N_ROW, ACC%N_VAL
   WRITE(MSG%LU_DEBUG,'(8I8)') ACC%ROW
   WRITE(MSG%LU_DEBUG,*) 'ACC%COL:'
   DO IC = 1, ACC%N_ROW-1
      WRITE(MSG%LU_DEBUG,'(16I8)') (ACC%COL(ICOL), ICOL=ACC%ROW(IC), ACC%ROW(IC+1)-1)
   ENDDO
   WRITE(MSG%LU_DEBUG,*) 'ACC%COL_GLOBAL:'
   DO IC = 1, ACC%N_ROW-1
      WRITE(MSG%LU_DEBUG,'(16I8)') (ACC%COL_GLOBAL(ICOL), ICOL=ACC%ROW(IC), ACC%ROW(IC+1)-1)
   ENDDO
   WRITE(MSG%LU_DEBUG,*) 'ACC%VAL:'
   DO IC = 1, ACC%N_ROW-1
      WRITE(MSG%LU_DEBUG,'(16E12.4)') (ACC%VAL(ICOL), ICOL=ACC%ROW(IC), ACC%ROW(IC+1)-1)
   ENDDO
#endif

#ifdef WITH_MKL
   IF (TYPE_MKL_PRECISION == NSCARC_PRECISION_SINGLE) THEN
       IF (TYPE_MKL(NL+1) == NSCARC_MKL_LOCAL .OR. TYPE_MKL(NL+1) == NSCARC_MKL_GLOBAL) &
          CALL SCARC_SETUP_MATRIX_MKL_SINGLE(NM, NL+1)
    ELSE
       IF (TYPE_MKL(NL+1) == NSCARC_MKL_LOCAL .OR. TYPE_MKL(NL+1) == NSCARC_MKL_GLOBAL) &
          CALL SCARC_SETUP_MATRIX_MKL_DOUBLE(NM, NL+1)
    ENDIF
#endif

ENDDO

END SUBROUTINE SCARC_SETUP_GALERKIN_OPERATOR

