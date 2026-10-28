! =====================================================================
! NAME       : icp_position.f90
! DESCRIPTION: Initial conditions for the position of all ParticleBase
!              solver classes. These initial conditions can be used to
!              initialize positions of particles in other solver
!              classes, using composition of initial condition
!              operators.
!
! Initial conditions avaliable:
!   read_x   : Reads positions from an input file numbered by stat
!   user_x   : Reads positions from a user defined input file
!   random_x : Random initial positions
!
! SPECTER version of GHOST's icp_position.f90. random_x seeds the
! particles only in the PHYSICAL domain (z in [0, nz-Cz-1] grid
! units; the FC-Gram continuation zone gets none), distributing the
! counts across the tasks proportionally to the physical span of
! each slab so that the number density is uniform (tasks whose slab
! lies entirely in the continuation zone get zero particles).
!
! DATE       : 03/21/26 (PDM)
! =====================================================================

module icp_position
  use icpbase_mod
  use equationbase_mod
  use particlebase_mod

  IMPLICIT NONE

  ! ================= Initial conditions supported ====================
  type, extends(icpBase) :: read_pos
    contains
      procedure :: init_GPState => init_readpos
  end type read_pos 
  type, extends(icpBase) :: user_pos
    contains
      procedure :: init_GPState => init_userpos
  end type user_pos 
  type, extends(icpBase) :: rand_pos
    contains
      procedure :: init_GPState => init_randompos
   end type rand_pos

CONTAINS

  ! ===================================================================
  ! Initial conditions
  ! ===================================================================

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Read particle positions from restart files
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_readpos(this,pde,fluidstate,psolver,pstate)
    use gpstate_mod
    use status
    use pstatus
    use filefmt
    use commtypes
    implicit none
    class    (read_pos),         intent   (in) :: this
    class(EquationBase),         intent   (in) :: pde
    type   (GStateComp), target, intent   (in) :: fluidstate(:) 
    class(ParticleBase),         intent(inout) :: psolver
    type  (GPStateComp),         intent(inout) :: pstate(:)

    if ((stat .eq. 0).and.(psolver%myrank_ .eq. 0)) then
      stop 'Cannot read restart files if starting a new run with stat=0'
    endif
    pind = int((stat-1)*lgmult+1)
    WRITE(lgext, lgfmtext) pind
    CALL io_read(psolver,pstate,1,psolver%idir_,trim(psolver%sstate_pos_),lgext)
  end subroutine init_readpos


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Read positions from a user defined input file. Pstates
  !! outside the one used for reading must be resized to the
  !! size of arrays in pstate after calling this routine.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_userpos(this,pde,fluidstate,psolver,pstate)
    use gpstate_mod
    use status
    use pstatus
    use filefmt
    use commtypes
    implicit none
    class    (user_pos),         intent   (in) :: this
    class(EquationBase),         intent   (in) :: pde
    type   (GStateComp), target, intent   (in) :: fluidstate(:) 
    class(ParticleBase),         intent(inout) :: psolver
    type  (GPStateComp),         intent(inout) :: pstate(:)
    real      (kind=GP)                        :: xp,yp,zp
    integer                                    :: nl,nt

    ! Note: each record (line) consists of x y z real positions
    ! within [0,NX-1]x[0,NY-1]x[0,NZ-1] box or the equivalent
    ! in box units depending on wrtunit_ class options.
    OPEN(UNIT=5,FILE=trim(psolver%seedfile_),STATUS='OLD',ACTION='READ',   &
         IOSTAT=psolver%ierr_, IOMSG=psolver%serr_);
    IF ( psolver%ierr_ .NE. 0 ) THEN
      WRITE(*,*)'Init_userpos: file:',trim(psolver%seedfile_),' err: ',    &
         trim(psolver%serr_) 
      STOP
    ENDIF
    READ(5,*,IOSTAT=psolver%ierr_) nt
    IF ( psolver%myrank_.eq.0 .AND. nt.NE.psolver%maxparts_ ) THEN
      WRITE(*,*) 'Init_userpos: Inconsistent seed file: required no. part.=', &
      psolver%maxparts_,' total listed: ',nt,' file:',psolver%seedfile_
      STOP
    ENDIF
    READ(5,*,IOSTAT=psolver%ierr_) xp

    nt = 0  ! global part. record counter
    nl = 0  ! local particle counter
    DO WHILE ( psolver%ierr_.EQ.0 )
      READ(5,*,IOSTAT=psolver%ierr_) xp, y, z
      IF ( psolver%ierr_ .NE. 0 ) THEN
        WRITE(*,*) 'Init_userpos: terminating read; nt=', nt
        EXIT
      ENDIF
      IF ( psolver%wrtunit_ .EQ. 1 ) THEN ! rescale coordinates to grid units
        xp = xp*psolver%invdel_(1)
        yp = yp*psolver%invdel_(2)
        zp = zp*psolver%invdel_(3)
      ENDIF
      IF ( zp.GE.psolver%lxbnds_(3,1) .AND. zp.LT.psolver%lxbnds_(3,2) .AND. &
           yp.GE.psolver%lxbnds_(2,1) .AND. yp.LT.psolver%lxbnds_(2,2) .AND. &
           xp.GE.psolver%lxbnds_(1,1) .AND. xp.LT.psolver%lxbnds_(1,2) ) THEN
        nl = nl + 1
        IF (nl.GT.psolver%partbuff_) THEN
          psolver%partbuff_ = psolver%partbuff_ + psolver%partchunksize_
          CALL psolver%ResizeArrays(psolver%partbuff_,.true.)   ! type-bound: subclasses re-acquire their pool pointers
          call GPState_resize(pstate ,psolver%partbuff_)
        END IF
        psolver%id_(nl)                      = nt
        pstate(psolver%POSITION  )%rcomp(nl) = xp
        pstate(psolver%POSITION+1)%rcomp(nl) = yp
        pstate(psolver%POSITION+2)%rcomp(nl) = zp
      ENDIF
      nt = nt + 1
    ENDDO
    CLOSE(5)

    psolver%nparts_ = nl;
    IF ((psolver%myrank_.NE.0).OR.(psolver%bcollective_.EQ.1)) THEN
      IF (psolver%nparts_.LT.(psolver%partbuff_-psolver%partchunksize_)) THEN
        psolver%partbuff_ = psolver%partbuff_ - ((psolver%partbuff_-psolver%nparts_) &
                        /psolver%partchunksize_-1)*psolver%partchunksize_
        CALL psolver%ResizeArrays(psolver%partbuff_,.false.)  ! type-bound: subclasses re-acquire their pool pointers
        call GPState_resize(pstate ,psolver%partbuff_)
      END IF
    END IF
    CALL MPI_ALLREDUCE(nl,nt,1,MPI_INTEGER,MPI_SUM,psolver%comm_,psolver%ierr_)
    IF ( psolver%myrank_.eq.0 .AND. nt.NE.psolver%maxparts_ ) THEN
      WRITE(*,*) 'Init_userpos: Inconsistent particle count: required no.=',    &
      psolver%maxparts_,' total read: ',nt,' file:',psolver%seedfile_
      STOP
    ENDIF
    IF (psolver%iexchtype_.EQ.GPEXCHTYPE_VDB) THEN
      CALL psolver%gpcomm_%VDBSynch(psolver%vdb_,psolver%maxparts_,psolver%id_, &
           pstate(psolver%POSITION  )%rcomp,                                    &
           pstate(psolver%POSITION+1)%rcomp,                                    &
           pstate(psolver%POSITION+2)%rcomp,                                    &
           psolver%nparts_,psolver%ptmp0_)
    END IF
  end subroutine init_userpos


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Initialize particle positions randomly
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_randompos(this,pde,fluidstate,psolver,pstate)
    use gpstate_mod
    use random
    implicit none
    class    (rand_pos),         intent   (in) :: this
    class(EquationBase),         intent   (in) :: pde
    class(ParticleBase),         intent(inout) :: psolver
    type   (GStateComp), target, intent   (in) :: fluidstate(:)
    type  (GPStateComp),         intent(inout) :: pstate(:)
    real      (kind=GP)                        :: c,r,x1,x2
    integer                                    :: ib,j,nt
    integer                                    :: kphys,klo,khi,nstart,nend

    ! SPECTER: distribute the particle counts across the tasks by the
    ! physical z-COORDINATE SPAN [0, gext_(3)] (grid units; gext_(3) =
    ! nz-Cz-1 is the top-wall coordinate), NOT equally by rank: tasks
    ! whose slab lies entirely in the FC-Gram continuation zone get no
    ! particles, and the number density is uniform in the physical
    ! domain. The integer division telescopes at shared slab
    ! boundaries, so the global total is exactly maxparts and the ids
    ! are unique and contiguous over [0, maxparts-1]. The z positions
    ! are seeded uniformly in [libnds(3,1)-1, min(libnds(3,2),gext_(3))],
    ! the same span the counts are proportional to.
    kphys  = NINT(psolver%gext_(3))                      ! nz-Cz-1
    klo    = MIN(psolver%libnds_(3,1)-1, kphys)          ! span below this slab
    khi    = MIN(psolver%libnds_(3,2),   kphys)          ! span up to slab top
    nstart = INT( (INT(psolver%maxparts_,8)*INT(klo,8)) / INT(kphys,8) )
    nend   = INT( (INT(psolver%maxparts_,8)*INT(khi,8)) / INT(kphys,8) )
    ib     = nstart                                      ! 0-based global id start
    psolver%nparts_ = nend - nstart

    ! Grow the local buffers if this task's physical share exceeds them
    ! (the physical-span distribution concentrates the particles in the
    ! lower slabs)
    IF ( psolver%nparts_.GT.psolver%partbuff_ ) THEN
      psolver%partbuff_ = (1+psolver%nparts_/psolver%partchunksize_)* &
                          psolver%partchunksize_
      CALL psolver%ResizeArrays(psolver%partbuff_,.true.)   ! type-bound: subclasses re-acquire their pool pointers
      call GPState_resize(pstate ,psolver%partbuff_)
    END IF

    DO j = 1, psolver%nparts_
      psolver%id_(j)    = ib + j - 1
      CALL prandom_number(r)
      x1 = real(psolver%libnds_(1,1)-1,kind=GP)
      x2 = real(psolver%libnds_(1,2),kind=GP)
      c = r*(psolver%nd_(1))
      pstate(psolver%POSITION  )%rcomp(j) = min(max(c,x1),x2)
      CALL prandom_number(r)
      x1 = real(psolver%libnds_(2,1)-1,kind=GP)
      x2 = real(psolver%libnds_(2,2),kind=GP)
      c = r*(psolver%nd_(2))
      pstate(psolver%POSITION+1)%rcomp(j) = min(max(c,x1),x2)
      CALL prandom_number(r)
      ! Seed z only within the PHYSICAL domain: the slab upper bound is
      ! clamped to the top wall so no particle lands in the FC-Gram
      ! continuation zone
      x1 = real(psolver%libnds_(3,1)-1,kind=GP)
      x2 = min(real(psolver%libnds_(3,2),kind=GP),psolver%gext_(3))
      pstate(psolver%POSITION+2)%rcomp(j) = min(max(x1+r*(x2-x1),x1),x2)
    ENDDO
    CALL MPI_ALLREDUCE(psolver%nparts_,nt,1,MPI_INTEGER,MPI_SUM,psolver%comm_,  &
                       psolver%ierr_)
    IF ( psolver%myrank_.eq.0 .AND. nt.NE.psolver%maxparts_ ) THEN
      WRITE(*,*) 'init_randompos: Inconsistent particle count: maxparts=',      &
      psolver%maxparts_,' total created: ',nt
      STOP
    ENDIF
    IF ( psolver%iexchtype_.EQ.GPEXCHTYPE_VDB ) THEN
      CALL psolver%gpcomm_%VDBSynch(psolver%vdb_,psolver%maxparts_,psolver%id_, &
                                    pstate(psolver%POSITION  )%rcomp,           &
                                    pstate(psolver%POSITION+1)%rcomp,           &
                                    pstate(psolver%POSITION+2)%rcomp,           &
                                    psolver%nparts_,psolver%ptmp0_)
      CALL GetLocalWrk(psolver,psolver%id_,                                     &
                       pstate(psolver%POSITION  )%rcomp,                        &
                       pstate(psolver%POSITION+1)%rcomp,                        &
                       pstate(psolver%POSITION+2)%rcomp,                        &
                       psolver%nparts_,psolver%vdb_,psolver%maxparts_)  
      IF ( psolver%wrtunit_ .EQ. 1 ) THEN ! rescale coordinates to box units
        psolver%ptmp0_(1,:) = psolver%vdb_(1,:)*psolver%delta_(1)
        psolver%ptmp0_(2,:) = psolver%vdb_(2,:)*psolver%delta_(2)
        psolver%ptmp0_(3,:) = psolver%vdb_(3,:)*psolver%delta_(3)
        CALL ascii_write_lag(psolver, 1, '.',                                   &
                             trim(psolver%sstate_pos_) // 'InitRndSeed',        &
                            '000',0.0_GP,psolver%maxparts_,                     &
                             psolver%ptmp0_(1,:),                               &
                             psolver%ptmp0_(2,:),                               &
                             psolver%ptmp0_(3,:))
      ELSE
        CALL ascii_write_lag(psolver,1,'.',                                     &
                             trim(psolver%sstate_pos_) // 'InitRndSeed',        &
                             '000',0.0_GP,psolver%maxparts_,                    &
                             psolver%vdb_(1,:),                                 &
                             psolver%vdb_(2,:),                                 &
                             psolver%vdb_(3,:))
      ENDIF
    ENDIF
    IF ( .NOT. PartNumConsistent(psolver,psolver%nparts_) ) THEN
      IF ( psolver%myrank_.eq.0 ) THEN
        WRITE(*,*) 'init_randompos: Invalid particle after GetLocalWrk call'
        STOP
      ENDIF
    ENDIF
  end subroutine init_randompos

end module icp_position
