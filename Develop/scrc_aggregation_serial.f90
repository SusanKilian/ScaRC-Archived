! ------------------------------------------------------------------------------------------------------
! Multimesh-version
! Compute ZONES for a compact matrix AC
! Returns the number of ZONES (== max(x[:]) + 1 )
! Unaggregated nodes are marked with a -1
! ------------------------------------------------------------------------------------------------------
SUBROUTINE SCARC_SETUP_AGGREGATION_ZONES0(NL)
USE SCARC_POINTERS, ONLY: GF, GC, CF
INTEGER, INTENT(IN) :: NL
INTEGER :: NM, NM2, IC, JC, ICOL, IZONE, JZONE, NMZ
INTEGER :: ACTIVE_MESH, N_ZONES_LOCAL, N_ZONES_PREVIOUS = 0
LOGICAL :: HAS_NEIGHBORS, HAS_AGGREGATED_NEIGHBORS
INTEGER, ALLOCATABLE, DIMENSION(:) :: ZONES_UP_TO_MESH

MESH_INT = 0

CALL SCARC_ALLOCATE_INT1 (ZONES_UP_TO_MESH, 1, NMESHES, NSCARC_INIT_ZERO, 'ZONES_UP_TO_MESH')

MESHES_LOOP: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
   ACTIVE_MESH = 1
   ZONES_LOOP: DO NMZ = 1, NMESHES

#ifdef WITH_SCARC_DEBUG
      WRITE(MSG%LU_DEBUG,*) '================================================================='
      WRITE(MSG%LU_DEBUG,*) '=============== AGGREGATING ON MESH ',NM,' FOR NMZ ', NMZ
      WRITE(MSG%LU_DEBUG,*) '================================================================='
#endif
      ! ---------- on the active mesh perform aggregation of internal cells and broadcast interfaces to neighbors
      ACTIVE_MESH_IF: IF (NM == ACTIVE_MESH) THEN

         N_ZONES_LOCAL = 0
         CALL SCARC_POINT_TO_MULTIGRID(NM, NL, NL+1)                    ! Sets grid pointers GF and GC
         CF => SCARC_POINT_TO_CMATRIX (GF, NSCARC_MATRIX_CONNECTION)

         GF%N_FINE = GF%NCE
         IF (NM > 1) THEN
            N_ZONES_PREVIOUS = ZONES_UP_TO_MESH(NM-1)               ! store previous number of zones
            GF%N_ZONES = N_ZONES_PREVIOUS + 1                       ! note: already been initialized to 1 on mesh 1 
         ENDIF

         !
         ! ----------- Pass 1:  Perform aggregation on internal cells of active mesh
         !
         PASS1_LOOP: DO IC = 1, GF%NC

            IF (GF%ZONES_GLOBAL(IC) /= 0) CYCLE                   ! has cell already been aggregated?

            HAS_NEIGHBORS = .FALSE.
            HAS_AGGREGATED_NEIGHBORS = .FALSE.

            DO ICOL = CF%ROW(IC), CF%ROW(IC+1)-1                  ! are all neighbors free (not already aggregated)?
               JC = CF%COL(ICOL)
               IF (IC /= JC .AND. JC <= GF%NC) THEN               ! only consider internal cells here
                  HAS_NEIGHBORS = .TRUE.
                  IF (GF%ZONES_GLOBAL(JC) /= 0) THEN
                     HAS_AGGREGATED_NEIGHBORS = .TRUE.
                     EXIT
                  ENDIF
               ENDIF
            ENDDO

            IF (.NOT. HAS_NEIGHBORS) THEN                        ! do not aggregate isolated cells
               GF%ZONES_GLOBAL(IC) = NSCARC_HUGE_INT
            ELSE IF (.NOT. HAS_AGGREGATED_NEIGHBORS) THEN        ! build aggregate of this cell and its neighbors
               GF%ZONES_GLOBAL(IC) = GF%N_ZONES
               GF%CPOINTS(GF%N_ZONES) = IC                
               DO ICOL = CF%ROW(IC), CF%ROW(IC+1)-1 
                  JC = CF%COL(ICOL)
                  IF (JC <= GF%NC) GF%ZONES_GLOBAL(CF%COL(ICOL)) = GF%N_ZONES
               ENDDO
               GF%N_ZONES = GF%N_ZONES + 1
            ENDIF

#ifdef WITH_SCARC_DEBUG
            WRITE(MSG%LU_DEBUG,*) 'SETUP_AGGREGATION_ZONES_INTERNAL: =============  IC:',IC
            WRITE(MSG%LU_DEBUG,*) 'ZONES:'
            WRITE(MSG%LU_DEBUG,CFORM_INT) GF%ZONES_GLOBAL(1:GF%NC)
            WRITE(MSG%LU_DEBUG,*) '-------------- overlap:'
            WRITE(MSG%LU_DEBUG,'(10I4)') GF%ZONES_GLOBAL(GF%NC+1: GF%NCE)
            WRITE(MSG%LU_DEBUG,*) 'CPOINTS:'
            WRITE(MSG%LU_DEBUG,CFORM_INT) GF%CPOINTS
#endif

         ENDDO PASS1_LOOP
   
         IF (NMESHES > 1) CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_ZONES, NSCARC_ZONES_INITIAL, NL)

         !
         ! ----------- Pass 2:  Add unaggregated nodes to neighboring aggregate
         !
         PASS2_LOOP: DO IC = 1, GF%NC
            IF (GF%ZONES_GLOBAL(IC) /= 0) CYCLE            
            DO ICOL = CF%ROW(IC), CF%ROW(IC+1)-1
               JC = CF%COL(ICOL)
               JZONE = GF%ZONES_GLOBAL(JC)
               IF (JZONE > 0) THEN
                  IF (JC >= GF%NC .OR. GF%CPOINTS(JZONE)>0) THEN
                     GF%ZONES_GLOBAL(IC) = -JZONE
                     EXIT
                  ENDIF
               ENDIF
            ENDDO
         ENDDO PASS2_LOOP
         GF%N_ZONES = GF%N_ZONES - 1
         !GF%ZONES = ABS(GF%ZONES)
   
#ifdef WITH_SCARC_DEBUG
         WRITE(MSG%LU_DEBUG,*) 'SETUP_AGGREGATION_ZONES_PASS2: =============  ', N_ZONES_LOCAL
         WRITE(MSG%LU_DEBUG,CFORM_INT) GF%ZONES_GLOBAL(1:GF%NC)
         WRITE(MSG%LU_DEBUG,*) '-------------- overlap:'
         WRITE(MSG%LU_DEBUG,'(10I4)') GF%ZONES_GLOBAL(GF%NC+1: GF%NCE)
#endif

         !
         ! ----------- Pass 3  Process remaining nodes which have not been aggregated yet
         !
         PASS3_LOOP: DO IC = 1, GF%NC
   
            IZONE = GF%ZONES_GLOBAL(IC)
   
            ! cell IC has not been aggregated
            IF (IZONE /= 0) THEN
               IF (IZONE > 0) THEN
                  GF%ZONES_GLOBAL(IC) = IZONE 
               ELSE IF (IZONE == NSCARC_HUGE_INT ) THEN
                  GF%ZONES_GLOBAL(IC) = -1
               ELSE
                  GF%ZONES_GLOBAL(IC) = -IZONE 
               ENDIF
               CYCLE PASS3_LOOP
            ENDIF
   
            GF%ZONES_GLOBAL(IC) = GF%N_ZONES
            GF%CPOINTS(GF%N_ZONES) = IC
   
            DO ICOL = CF%ROW(IC), CF%ROW(IC+1)-1
               JC = CF%COL(ICOL)
               IF (JC <= GF%NC .AND. GF%ZONES_GLOBAL(JC) == 0) GF%ZONES_GLOBAL(JC) = GF%N_ZONES
            ENDDO
            GF%N_ZONES = GF%N_ZONES + 1
   
         ENDDO PASS3_LOOP
   
#ifdef WITH_SCARC_DEBUG
         WRITE(MSG%LU_DEBUG,*) 'SETUP_AGGREGATION_ZONES_PASS3: =============  ', N_ZONES_LOCAL
         WRITE(MSG%LU_DEBUG,CFORM_INT) GF%ZONES_GLOBAL(1:GF%NC)
         WRITE(MSG%LU_DEBUG,*) '-------------- overlap:'
         WRITE(MSG%LU_DEBUG,'(10I4)') GF%ZONES_GLOBAL(GF%NC+1: GF%NCE)
#endif
   
         ! Initialize number of coarse and fine grid cells on current SAMG level 
         GF%N_COARSE = GF%N_ZONES
         GC%N_FINE   = GF%N_ZONES
         ZONES_UP_TO_MESH(ACTIVE_MESH) = GF%N_ZONES
      
         GC%NC_LOCAL(NM) = GF%N_ZONES - N_ZONES_PREVIOUS 
         MESH_INT(NM) = GC%NC_LOCAL(NM)
 
         ! only ACTIVE_MESH broadcasts its whole zone information along interfaces
         ! in order to make zones consistent along interfaces 
         IF (NMESHES > 1) CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_ZONES, NSCARC_ZONES_INITIAL, NL)
         IF (N_MPI_PROCESSES>1) &
            CALL MPI_ALLGATHERV(MPI_IN_PLACE, 1, MPI_INTEGER, ZONES_UP_TO_MESH, COUNTS, DISPLS, &
                                MPI_INTEGER, MPI_COMM_WORLD, IERROR)
      
         ACTIVE_MESH = ACTIVE_MESH + 1
      !
      ! ---------- NM is NOT the active mesh
      ! if NM > ACTIVE_MESH, then the complete zone information must be received for both the overlapping 
      ! cells and the adjacent inner cells
      ! if NM < ACTIVE_MESH, then an initial exchange has already taken place and only the overlapping cells
      ! still have to be filled up by the zones of the respective neighbours that have been calculated before
      !
      ELSE ACTIVE_MESH_IF

         IF (NMESHES > 1) THEN
            IF (NM > ACTIVE_MESH) THEN
               CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_ZONES, NSCARC_ZONES_COMPLETE, NL)
            ELSE IF (NM < ACTIVE_MESH) THEN
               CALL SCARC_EXCHANGE (NSCARC_EXCHANGE_ZONES, NSCARC_ZONES_FILLUP, NL)
            ENDIF
         ENDIF
         IF (N_MPI_PROCESSES>1) &
            CALL MPI_ALLGATHERV(MPI_IN_PLACE, 1, MPI_INTEGER, ZONES_UP_TO_MESH, COUNTS, DISPLS, &
                                MPI_INTEGER, MPI_COMM_WORLD, IERROR)

         ACTIVE_MESH = ACTIVE_MESH + 1

      ENDIF ACTIVE_MESH_IF
   ENDDO ZONES_LOOP
ENDDO MESHES_LOOP

DEALLOCATE (ZONES_UP_TO_MESH)

IF (TRIM(CHID) == 'b1') THEN
   CALL SCARC_PRESET_B1_CASE
   MESH_INT(1) = 6
ENDIF

!
! Broadcast number of local mesh cells on level NL to all and build global sum
!
IF (N_MPI_PROCESSES > 1) &
   CALL MPI_ALLGATHERV(MPI_IN_PLACE, 1, MPI_INTEGER, MESH_INT, COUNTS, DISPLS, MPI_INTEGER, MPI_COMM_WORLD, IERROR)

!
! Store information on local and global cells numbers on data structure of corresponding discretization type
!
MESHES_LOOP2: DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX

   CALL SCARC_POINT_TO_MULTIGRID(NM, NL, NL+1)

   GC%NC_LOCAL = MESH_INT
   GC%NC_GLOBAL = SUM(MESH_INT)

   IF (TRIM(CHID) == 'b1') THEN
      GF%NC = 24
      GF%NCE = 24
      GC%NC_LOCAL  = 6
      GC%NC_GLOBAL = 6
      GC%NC = 6
      GC%NCE = 6
   ENDIF
   CALL SCARC_RESIZE_INT1(GF%CPOINTS, GF%N_ZONES, 'GF%CPOINTS')     

#ifdef WITH_SCARC_DEBUG
      WRITE(MSG%LU_DEBUG,*) 'GC%NC_LOCAL=', GC%NC_LOCAL
      WRITE(MSG%LU_DEBUG,*) 'GC%NC_GLOBAL=', GC%NC_GLOBAL
#endif

   DO IC = GF%NC+1, GF%NCE
      GF%ZONES_GLOBAL(IC) = ABS(GF%ZONES_GLOBAL(IC))
   ENDDO
   !IF (MINVAL(GF%ZONES) <= 0) CALL SCARC_SHUTDOWN(NSCARC_ERROR_AGGREGATION, SCARC_NONE, -1)

   IF (NMESHES > 1) THEN
      DO NM2=2,NMESHES
         GC%NC_OFFSET(NM2) = GC%NC_OFFSET(NM2-1) + GC%NC_LOCAL(NM2-1)
      ENDDO
   ENDIF

#ifdef WITH_SCARC_DEBUG
WRITE(MSG%LU_DEBUG,*) 'NC_GLOBAL(',NL+1,')=', NC_GLOBAL(NL+1)
WRITE(MSG%LU_DEBUG,*) 'GC%NC_LOCAL :', GC%NC_LOCAL
WRITE(MSG%LU_DEBUG,*) 'GC%NC_GLOBAL:', GC%NC_GLOBAL
WRITE(MSG%LU_DEBUG,*) 'GC%NC_OFFSET:', GC%NC_OFFSET
#endif

ENDDO MESHES_LOOP2

END SUBROUTINE SCARC_SETUP_AGGREGATION_ZONES0


