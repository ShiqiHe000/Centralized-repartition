

PROGRAM FOUR_PROCS
    USE INPUT, ONLY: MESHFILE, LP_SERIAL
    USE SIZE
    USE MESH
    USE ADAPT
    USE FIELDS
    USE LOAD_PARTITION
    USE LOAD_COAREN
    USE H_COARSENING
    USE MPI
    

    IMPLICIT NONE
    
    !-------------------------------------------------------------------
    INTEGER :: I
    
    !-------------------------------------------------------------------
    
    ! START MPI
    CALL START_MPI
    
    ! READ MESH FILE
    CALL READ_MESH(MESHFILE)
    
    ! INITIALIZE ELEMENT POLYNOMIAL ORDER
    CALL HPADAPT_INIT
    
    ! WRITE DATA TO FILE
    CALL WRITE_FIELDS

    ! PARALLEL LOAD PARTITIONING
    IF( .NOT. LP_SERIAL) THEN

        DO I=1, 45

            CALL DEALLOCATE_GLOBAL
            
            ! H-REFINEMENT
            CALL HADAPT
            
            ! P-REFINEMENT
            CALL PADAPT
            
            ! WRITE DATA TO FILE
            CALL WRITE_FIELDS
            
            CALL DEALLOCATE_GLOBAL
            
            CALL H1_PARTITION_LOCAL
          
            CALL WRITE_FIELDS

        ENDDO
    
    ELSE    ! SERIAL LOAD PARITITONING
        DO I=1,30 
            ! H-REFINEMENT
            CALL HADAPT
            
            ! P-REFINEMENT
            CALL PADAPT
            
            CALL DEALLOCATE_GLOBAL
            
            ! WRITE DATA TO FILE
            CALL WRITE_FIELDS
            
            CALL H1_H2_PARTITION_SERIAL
          
            CALL WRITE_FIELDS
            
        ENDDO
    
    ENDIF
    
    ! START COARSENING
    START_COARSEN = .TRUE.
    
    CALL DEALLOCATE_GLOBAL
    
    DO I=1, 25
    
        CALL P_COARSEN
        
        CALL SET_H_COARSEN_FLAG
    
        CALL WRITE_FIELDS
        
        CALL H_COARSEN
        
        CALL H1_H2_PARTITION_SERIAL
        
        CALL WRITE_FIELDS
        
        CALL DEALLOCATE_GLOBAL
        
    ENDDO
    
    

    
    CALL MPI_FINALIZE(IERROR)

    IF (RANK == 0) THEN
        PRINT *,  "FINISH WRITING DATA TO FILE"
    ENDIF
    

END PROGRAM FOUR_PROCS


