!=================================================================
! FFTP3D v3
! Parallel Fast Fourier Transform in 3D
!
! Performs parallel real-to-complex and complex-to-real FFTs 
! using MPI and the FFTW 3.x library in each node. You should use 
! the FFTPLANS and MPIVARS modules (see the file 'fftp_mod.f90') 
! in each program that calls any of the subroutines in this 
! file. Also, you must create plans for the parallel FFT using 
! the 'fftp3d_create_plan' subroutine, which creates in turn 
! derived data types for message passing using the subroutine 
! 'fftp3d_create_block'.
!
! When dealing with non-periodic data, an FC-Gram scheme
! is utilized to eliminate Gibbs phenomenon.
!
! 2007 Pablo D. Mininni.
!      National Center for Atmospheric Research.
!      e-mail: mininni@ucar.uba.ar 
!
! 16 Feb 2004: Performs complex FFTs in place.
!  8 Jul 2004: itype pointers only used to store datatypes of 
!              blocks in the row and column each processor is.
!  9 Jul 2004: Transposition uses data cache blocking.
! 13 Feb 2007: Transposition uses strip mining (rreddy@psc.edu)
! 25 Aug 2009: Hybrid MPI/OpenMP support (D. Rosenberg & P. Mininni)
! 30 Aug 2009: SINGLE/DOUBLE precision (D. Rosenberg & P. Mininni)
!  3 Jan 2017: Anisotropic boxes (P. Mininni)
! 29 Aug 2019: FC-Gram (M. Fontana)
!
! References:
! Lyon, M.; SIAM J. Sci. Comp., 33(6) (2011). DOI:10.1137/11082436X
! Mininni PD, Gomez DO, Mahajan SM; Astrophys. J. 619, 1019 (2005)
! Gomez DO, Mininni PD, Dmitruk P; Phys. Scripta T116, 123 (2005)
! Gomez DO, Mininni PD, Dmitruk P; Adv. Sp. Res. 35, 899 (2005)
!=================================================================
#include "fftw_wrappers.h"

!*****************************************************************
      SUBROUTINE fftp3d_init_threads(err)
!-----------------------------------------------------------------
!
! Initializes FFTW threads.
!
! Parameters
!     err : if zero, the initialization failed
!-----------------------------------------------------------------

!$    USE threads
      IMPLICIT NONE

      INTEGER, INTENT(INOUT) :: err

!$    CALL GPMANGLE(init_threads)(err)
      IF (err.eq.0) PRINT *,'FFTP threads initialization failed!'

      RETURN
      END SUBROUTINE fftp3d_init_threads

!*****************************************************************
      SUBROUTINE fftp3d_create_plan_rc(plan,n,C,o,tdir,flags)
!-----------------------------------------------------------------
!
! Creates plans for the FFTW in each node.
!
! Parameters
!     plan   : contains the parallel 3D plan [OUT]
!     n      : the size of the dimensions of the input array [IN]
!     C      : number of continuation points in each direction [IN]
!     d      : number of matching points in each direction [IN]
!     tdir   : directory containing the FC-Gram tables [IN]
!     flags  : flags for the FFTW [IN]
!              FFTW_ESTIMATE (sub-optimal but faster)
!              FFTW_MEASURE (optimal but slower to create plans)
!              FFTW_PATIENT AND FFTW_EXHAUSTIVE are also available
!              for extra performance, but may take a long time to
!              create plans (specially when using OpenMP)
!-----------------------------------------------------------------

      USE mpivars
      USE fftplans
!$    USE threads
      USE gtimer
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: n(3),C(3),o(3)
      INTEGER, INTENT(IN) :: flags
      CHARACTER(len=*), INTENT(IN) :: tdir
      TYPE(FFTPLAN), INTENT(OUT) :: plan

      DOUBLE PRECISION, ALLOCATABLE :: A(:,:),Q(:,:)
      INTEGER :: i,j
      CHARACTER(len=5) :: cstr,dstr

      ALLOCATE ( plan%ccarr(n(3),n(2),ista:iend)    )
      ALLOCATE ( plan%carr(n(1)/2+1,n(2),ksta:kend) )
      ALLOCATE ( plan%rarr(n(1),n(2),ksta:kend)     )
!$    CALL GPMANGLE(plan_with_nthreads)(nth)

      CALL GPMANGLE(plan_many_dft_r2c)(plan%planr,2,(/n(1),n(2)/),    &
                         kend-ksta+1,plan%rarr,                       &
                         (/n(1),n(2)*(kend-ksta+1)/),1,n(1)*n(2),     &
                         plan%carr,(/n(1)/2+1,n(2)*(kend-ksta+1)/),1, &
                         (n(1)/2+1)*n(2),flags)
      CALL GPMANGLE(plan_many_dft)(plan%planc,1,n(3),n(2)*(iend-ista+1), &
                         plan%ccarr,(iend-ista+1)*n(2)*n(3),1,n(3),      &
                         plan%ccarr,(iend-ista+1)*n(2)*n(3),1,n(3),      &
                         -1,flags)
      plan%nx = n(1)
      plan%ny = n(2)
      plan%nz = n(3)
      ALLOCATE( plan%itype1(0:nprocs-1) )
      ALLOCATE( plan%itype2(0:nprocs-1) )
      CALL fftp3d_create_block(n,nprocs,myrank,plan%itype1, &
                              plan%itype2)

      plan%Cx = C(1)
      plan%Cy = C(2)
      plan%Cz = C(3)
      plan%ox = o(1)
      plan%oy = o(2)
      plan%oz = o(3)
     
      ! Load continuation matrices. Note that only the transpose of Q
      ! is needed.
      IF ( plan%Cx .ne. 0 ) THEN
         ALLOCATE( A(plan%Cx,plan%ox), Q(plan%ox,plan%ox) )
         WRITE(cstr,'(I5)') plan%Cx
         WRITE(dstr,'(I5)') plan%ox

         OPEN(10, FILE=trim(tdir) // '/A' // trim(adjustl(cstr)) // '-' // &
              trim(adjustl(dstr)) //  '.dat', &
              FORM='unformatted', ACCESS='stream', ACTION='read')
         READ(10) A
         CLOSE(10)
         OPEN(10, FILE=trim(tdir) // '/Q' // trim(adjustl(dstr)) // '.dat',&
              FORM='unformatted', ACCESS='stream', ACTION='read')
         READ(10) Q
         CLOSE(10)
         Q = TRANSPOSE(Q)
         plan%Gxf = MATMUL(A(plan%Cx:1:-1,:),Q)
         plan%Gxl = MATMUL(A,Q)
         DEALLOCATE(A,Q)
      ENDIF
      IF ( plan%Cy .ne. 0 ) THEN
         ALLOCATE( A(plan%Cy,plan%oy), Q(plan%oy,plan%oy) )
         WRITE(cstr,'(I5)') plan%Cy
         WRITE(dstr,'(I5)') plan%oy

         OPEN(10, FILE=trim(tdir) // '/A' // trim(adjustl(cstr)) // '-' // &
              trim(adjustl(dstr)) //  '.dat', &
              FORM='unformatted', ACCESS='stream', ACTION='read')
         READ(10) A
         CLOSE(10)
         OPEN(10, FILE=trim(tdir) // '/Q' // trim(adjustl(dstr)) // '.dat',&
              FORM='unformatted', ACCESS='stream', ACTION='read')
         READ(10) Q
         CLOSE(10)
         Q = TRANSPOSE(Q)
         plan%Gyf = MATMUL(A(plan%Cy:1:-1,:),Q)
         plan%Gyl = MATMUL(A,Q)
         DEALLOCATE(A,Q)
     ENDIF
      IF ( plan%Cz .ne. 0 ) THEN
         ALLOCATE( A(plan%Cz,plan%oz), Q(plan%oz,plan%oz) )
         WRITE(cstr,'(I5)') plan%Cz
         WRITE(dstr,'(I5)') plan%oz

         OPEN(10, FILE=trim(tdir) // '/A' // trim(adjustl(cstr)) // '-' // &
              trim(adjustl(dstr)) //  '.dat', &
              FORM='unformatted', ACCESS='stream', ACTION='read')
         READ(10) A
         CLOSE(10)
         OPEN(10, FILE=trim(tdir) // '/Q' // trim(adjustl(dstr)) // '.dat',&
              FORM='unformatted', ACCESS='stream', ACTION='read')
         READ(10) Q
         CLOSE(10)
         Q = TRANSPOSE(Q)
         plan%Gzf = MATMUL(A(plan%Cz:1:-1,:),Q(:,plan%oz:1:-1))
         plan%Gzl = MATMUL(A,Q)
        DEALLOCATE(A,Q)
    ENDIF

      CALL GTStart(hcom,GT_WTIME)
      CALL GTStart(hfft,GT_WTIME)
      CALL GTStart(htra,GT_WTIME)
      CALL GTStart(hcont,GT_WTIME)
      CALL GTStart(htot,GT_WTIME)


      RETURN
      END SUBROUTINE fftp3d_create_plan_rc

!*****************************************************************
      SUBROUTINE fftp3d_create_plan_cr(plan,n,flags)
!-----------------------------------------------------------------
!
! Creates plans for the FFTW in each node.
!
! Parameters
!     plan   : contains the parallel 3D plan [OUT]
!     n      : the size of the dimensions of the input array [IN]
!     flags  : flags for the FFTW [IN]
!              FFTW_ESTIMATE (sub-optimal but faster)
!              FFTW_MEASURE (optimal but slower to create plans)
!              FFTW_PATIENT AND FFTW_EXHAUSTIVE are also available
!              for extra performance, but may take a long time to
!              create plans (specially when using OpenMP)
!-----------------------------------------------------------------

      USE mpivars
      USE fftplans
!$    USE threads
      USE gtimer
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: n(3)
      INTEGER, INTENT(IN) :: flags
      TYPE(FFTPLAN), INTENT(OUT) :: plan

      ALLOCATE ( plan%ccarr(n(3),n(2),ista:iend)    )
      ALLOCATE ( plan%carr(n(1)/2+1,n(2),ksta:kend) )
      ALLOCATE ( plan%rarr(n(1),n(2),ksta:kend)     )
!$    CALL GPMANGLE(plan_with_nthreads)(nth)

      CALL GPMANGLE(plan_many_dft_c2r)(plan%planr,2,(/n(1),n(2)/),    &
                         kend-ksta+1,plan%carr,                       &
                         (/n(1)/2+1,n(2)*(kend-ksta+1)/),1,           &
                         (n(1)/2+1)*n(2),plan%rarr,                   &
                         (/n(1),n(2)*(kend-ksta+1)/),1,n(1)*n(2),flags)

      CALL GPMANGLE(plan_many_dft)(plan%planc,1,n(3),n(2)*(iend-ista+1), &
                         plan%ccarr,(iend-ista+1)*n(2)*n(3),1,n(3),      &
                         plan%ccarr,(iend-ista+1)*n(2)*n(3),1,n(3),      &
                         1,flags)
      plan%nx = n(1)
      plan%ny = n(2)
      plan%nz = n(3)

      plan%Cx = 0
      plan%Cy = 0
      plan%Cz = 0
      plan%ox = 0
      plan%oy = 0
      plan%oz = 0

      ALLOCATE( plan%itype1(0:nprocs-1) )
      ALLOCATE( plan%itype2(0:nprocs-1) )
   
      CALL fftp3d_create_block(n,nprocs,myrank,plan%itype1, &
                              plan%itype2)

      CALL GTStart(hcom,GT_WTIME)
      CALL GTStart(hfft,GT_WTIME)
      CALL GTStart(htra,GT_WTIME)
      CALL GTStart(hcont,GT_WTIME)
      CALL GTStart(htot,GT_WTIME)


      RETURN
      END SUBROUTINE fftp3d_create_plan_cr


!*****************************************************************
      SUBROUTINE fftp3d_destroy_plan(plan)
!-----------------------------------------------------------------
!
! Destroys FFTW plans in each node.
!
! Parameters
!     plan : the parallel 3D plan [INOUT]
!-----------------------------------------------------------------

      USE fftplans
      USE gtimer
      IMPLICIT NONE

      TYPE(FFTPLAN), INTENT(INOUT) :: plan

      CALL GPMANGLE(destroy_plan)(plan%planr)
      CALL GPMANGLE(destroy_plan)(plan%planc)
      DEALLOCATE( plan%ccarr  )
      DEALLOCATE( plan%carr   )
      DEALLOCATE( plan%rarr   )
      DEALLOCATE( plan%itype1 )
      DEALLOCATE( plan%itype2 )

      IF ( plan%Cx .ne. 0 ) DEALLOCATE( plan%Gxf,plan%Gxl )
      IF ( plan%Cy .ne. 0 ) DEALLOCATE( plan%Gyf,plan%Gyl )
      IF ( plan%Cz .ne. 0 ) DEALLOCATE( plan%Gzf,plan%Gzl )

      CALL GTFree(hcom)
      CALL GTFree(hfft)
      CALL GTFree(htra)
      CALL GTFree(hcont)
      CALL GTFree(htot)

      RETURN
      END SUBROUTINE fftp3d_destroy_plan

!*****************************************************************
      SUBROUTINE fftp3d_create_block(n,nprocs,myrank,itype1,itype2)
!-----------------------------------------------------------------
!
! Defines derived data types for sending and receiving 
! blocks of the 3D matrix between processors. The data 
! types are used to transpose the matrix during the FFT.
!
! Parameters
!     n      : the size of the dimensions of the input array [IN]
!     nprocs : the number of processors [IN]
!     myrank : the rank of the processor [IN]
!     itype1 : contains a derived data type for sending [OUT]
!     itype2 : contains a derived data type for receiving [OUT]
!-----------------------------------------------------------------

      USE commtypes
      IMPLICIT NONE

      INTEGER, INTENT(OUT), DIMENSION(0:nprocs-1) :: itype1,itype2
      INTEGER, INTENT(IN) :: n(3),nprocs
      INTEGER, INTENT(IN) :: myrank

      INTEGER :: ista,iend
      INTEGER :: ksta,kend
      INTEGER :: irank,krank
      INTEGER :: itemp1,itemp2

      CALL range(1,n(3),nprocs,myrank,ksta,kend)
      DO irank = 0,nprocs-1
         CALL range(1,n(1)/2+1,nprocs,irank,ista,iend)
         CALL block3d(1,n(1)/2+1,1,n(2),ksta,ista,iend,1,n(2), &
                     ksta,kend,GC_COMPLEX,itemp1)
         itype1(irank) = itemp1
      END DO
      CALL range(1,n(1)/2+1,nprocs,myrank,ista,iend)
      DO krank = 0,nprocs-1
         CALL range(1,n(3),nprocs,krank,ksta,kend)
         CALL block3d(ista,iend,1,n(2),1,ista,iend,1,n(2),     &
                     ksta,kend,GC_COMPLEX,itemp2)
         itype2(krank) = itemp2
      END DO

      RETURN
      END SUBROUTINE fftp3d_create_block

!*****************************************************************
      SUBROUTINE fftp3d_real_to_complex(plan,in,out,comm)
!-----------------------------------------------------------------
!
! Computes the 3D real-to-complex FFT in parallel. The
! complex output has the same structure than the output
! of the 3D FFTW, but the output is transposed.
!
! Parameters
!     plan : the 3D plan created with fftp3d_create_plan [IN]
!     in   : real input array [IN]
!     out  : complex output array [OUT]
!     comm : the MPI communicator (handle) [IN]
!-----------------------------------------------------------------

      USE commtypes
      USE fprecision
      USE mpivars
      USE fftplans
      USE gtimer
!$    USE threads
      IMPLICIT NONE

      TYPE(FFTPLAN), INTENT(IN) :: plan

      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(plan%nz,plan%ny,ista:iend) :: out
      COMPLEX(KIND=GP), DIMENSION(ista:iend,plan%ny,plan%nz)           :: c1
      REAL(KIND=GP), INTENT(IN), DIMENSION(plan%nx,plan%ny,ksta:kend) :: IN

      DOUBLE PRECISION                    :: t0, t1
 
      INTEGER, DIMENSION(0:nprocs-1)      :: ireq1,ireq2
      INTEGER, DIMENSION(MPI_STATUS_SIZE) :: istatus
      INTEGER, INTENT(IN)                 :: comm


      INTEGER :: i,j,k
      INTEGER :: ii,jj,kk
      INTEGER :: irank
      INTEGER :: isendTo,igetFrom
      INTEGER :: istrip,iproc

      CALL GTStart(htot)
!
! 2D FFT in each node using the FFTW library
!
      CALL GTStart(hfft)
      CALL GPMANGLE(execute_dft_r2c)(plan%planr,in,plan%carr)
      CALL GTStop(hfft); ffttime = ffttime + GTGetTime(hfft); 
   
!
!
! Transposes the result between nodes using 
! strip mining when nstrip>1 (rreddy@psc.edu)
!
      CALL GTStart(hcom)
      do iproc = 0, nprocs-1, nstrip
         do istrip=0, nstrip-1
            irank = iproc + istrip

            isendTo = myrank + irank
            if ( isendTo .ge. nprocs ) isendTo = isendTo - nprocs

            igetFrom = myrank - irank
            if ( igetFrom .lt. 0 ) igetFrom = igetFrom + nprocs
            CALL MPI_IRECV(c1,1,plan%itype2(igetFrom),igetFrom,      & 
                          1,comm,ireq2(irank),ierr)

            CALL MPI_ISEND(plan%carr,1,plan%itype1(isendTo),isendTo, &
                          1,comm,ireq1(irank),ierr)
         enddo

         do istrip=0, nstrip-1
            irank = iproc + istrip
            CALL MPI_WAIT(ireq1(irank),istatus,ierr)
            CALL MPI_WAIT(ireq2(irank),istatus,ierr)
         enddo
      enddo
      CALL GTStop(hcom); comtime = comtime + GTGetTime(hcom)
!
! Cache friendly transposition
!
      CALL GTStart(htra)
!$omp parallel do if ((iend-ista)/csize.ge.nth) private (jj,kk,i,j,k)
      DO ii = ista,iend,csize
!$omp parallel do if ((iend-ista)/csize.lt.nth) private (kk,i,j,k)
         DO jj = 1,plan%ny,csize
            DO kk = 1,plan%nz,csize
               DO i = ii,min(iend,ii+csize-1)
               DO j = jj,min(plan%ny,jj+csize-1)
               DO k = kk,min(plan%nz,kk+csize-1)
                  out(k,j,i) = c1(i,j,k)
               END DO
               END DO
               END DO
            END DO
         END DO
      END DO
      CALL GTStop(htra); tratime = tratime + GTGetTime(htra)

      CALL GTStart(hcont)
! Continuation in Z direction, where the box is nonperiodic
!TODO Explore performance benefits of manual MATMUL 
!$omp parallel do if ((iend-ista)/csize.ge.nth) private (j)
      DO i= ista,iend
!$omp parallel do if ((iend-ista)/csize.lt.nth)
         DO j =1,plan%ny
            ! Continuation using last o points
            out(plan%nz-plan%Cz+1:plan%nz,j,i) = &
               MATMUL(plan%Gzl, out(plan%nz-plan%Cz-plan%oz+1:plan%nz-plan%Cz,j, i))
   
            ! Continuation using first o points
            out(plan%nz-plan%Cz+1:plan%nz,j,i) = out(plan%nz-plan%Cz+1:plan%nz,&
                    j, i) + MATMUL(plan%Gzf, out(1:plan%oz,j,i))
         ENDDO
      ENDDO
      CALL GTStop(hcont); conttime = conttime + GTGetTime(hcont)

!
! 1D FFT in each node using the FFTW library
!
      CALL GTStart(hfft)
      CALL GPMANGLE(execute_dft)(plan%planc,out,out)
      CALL GTStop(hfft); ffttime = ffttime + GTGetTime(hfft)

      CALL GTStop(htot); tottime = tottime + GTGetTime(htot)

      RETURN
      END SUBROUTINE fftp3d_real_to_complex

!*****************************************************************
      SUBROUTINE fftp2d_real_to_complex_xy(plan,in,out,comm)
!-----------------------------------------------------------------
!
! Computes the 3D real-to-complex FFT in parallel. The
! complex output has the same structure than the output
! of the 3D FFTW, but the output is transposed.
!
! Parameters
!     plan : the 3D plan created with fftp3d_create_plan [IN]
!     in   : real input array [IN]
!     out  : complex output array [OUT]
!     comm : the MPI communicator (handle) [IN]
!-----------------------------------------------------------------

      USE commtypes
      USE fprecision
      USE mpivars
      USE fftplans
      USE gtimer
!$    USE threads
      IMPLICIT NONE

      TYPE(FFTPLAN), INTENT(IN) :: plan

      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(plan%nz,plan%ny,ista:iend) :: out
      COMPLEX(KIND=GP), DIMENSION(ista:iend,plan%ny,plan%nz)          :: c1
      REAL(KIND=GP), INTENT(IN), DIMENSION(plan%nx,plan%ny,ksta:kend) :: IN

      INTEGER, DIMENSION(0:nprocs-1)      :: ireq1,ireq2
      INTEGER, DIMENSION(MPI_STATUS_SIZE) :: istatus
      INTEGER, INTENT(IN)                 :: comm
      INTEGER :: i,j,k
      INTEGER :: ii,jj,kk
      INTEGER :: irank
      INTEGER :: isendTo,igetFrom
      INTEGER :: istrip,iproc

      CALL GTStart(htot)
!
! 2D FFT in each node using the FFTW library
!
      CALL GTStart(hfft)
      CALL GPMANGLE(execute_dft_r2c)(plan%planr,in,plan%carr)
      CALL GTStop(hfft); ffttime = ffttime + GTGetTime(hfft); 
   
!
!
! Transposes the result between nodes using 
! strip mining when nstrip>1 (rreddy@psc.edu)
!
      CALL GTStart(hcom)
      do iproc = 0, nprocs-1, nstrip
         do istrip=0, nstrip-1
            irank = iproc + istrip

            isendTo = myrank + irank
            if ( isendTo .ge. nprocs ) isendTo = isendTo - nprocs

            igetFrom = myrank - irank
            if ( igetFrom .lt. 0 ) igetFrom = igetFrom + nprocs
            CALL MPI_IRECV(c1,1,plan%itype2(igetFrom),igetFrom,      & 
                          1,comm,ireq2(irank),ierr)

            CALL MPI_ISEND(plan%carr,1,plan%itype1(isendTo),isendTo, &
                          1,comm,ireq1(irank),ierr)
         enddo

         do istrip=0, nstrip-1
            irank = iproc + istrip
            CALL MPI_WAIT(ireq1(irank),istatus,ierr)
            CALL MPI_WAIT(ireq2(irank),istatus,ierr)
         enddo
      enddo
      CALL GTStop(hcom); comtime = comtime + GTGetTime(hcom)
!
! Cache friendly transposition
!
      CALL GTStart(htra)
!$omp parallel do if ((iend-ista)/csize.ge.nth) private (jj,kk,i,j,k)
      DO ii = ista,iend,csize
!$omp parallel do if ((iend-ista)/csize.lt.nth) private (kk,i,j,k)
         DO jj = 1,plan%ny,csize
            DO kk = 1,plan%nz,csize
               DO i = ii,min(iend,ii+csize-1)
               DO j = jj,min(plan%ny,jj+csize-1)
               DO k = kk,min(plan%nz,kk+csize-1)
                  out(k,j,i) = c1(i,j,k)
               END DO
               END DO
               END DO
            END DO
         END DO
      END DO
      CALL GTStop(htra); tratime = tratime + GTGetTime(htra)
      CALL GTStop(htot); tottime = tottime + GTGetTime(htot)

      RETURN
      END SUBROUTINE fftp2d_real_to_complex_xy

!*****************************************************************
      SUBROUTINE fftp1d_real_to_complex_z(plan,inout,comm)
!-----------------------------------------------------------------
!
! Computes the 3D real-to-complex FFT in parallel. The
! complex output has the same structure than the output
! of the 3D FFTW, but the output is transposed.
!
! Parameters
!     plan : the 3D plan created with fftp3d_create_plan [IN]
!     in   : contains f(z,ky,kx) as input and returns f*(kz,ky,kx) [INOUT]
!     comm : the MPI communicator (handle) [IN]
!-----------------------------------------------------------------

      USE commtypes
      USE fprecision
      USE mpivars
      USE fftplans
      USE gtimer
!$    USE threads
      IMPLICIT NONE

      TYPE(FFTPLAN), INTENT(IN) :: plan
      COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(plan%nz,plan%ny,ista:iend) :: inout
      INTEGER, INTENT(IN)                 :: comm

      INTEGER :: i,j,k

      CALL GTStart(htot)
      CALL GTStart(hcont)
! Continuation in Z direction, where the box is nonperiodic
! Q is already transposed
!$omp paralleldo if ((iend-ista)/csize.ge.nth) private (j)
      DO i= ista,iend
!$omp parallel do if ((iend-ista)/csize.lt.nth)
         DO j =1,plan%ny
            ! Continuation using last o points
            inout(plan%nz-plan%Cz+1:plan%nz,j,i) = &
               MATMUL(plan%Gzl, inout(plan%nz-plan%Cz-plan%oz+1:plan%nz-plan%Cz,j, i))

            ! Continuation using first o points
            inout(plan%nz-plan%Cz+1:plan%nz,j,i) = inout(plan%nz-plan%Cz+1:plan%nz,&
                    j, i) + MATMUL(plan%Gzf, inout(1:plan%oz,j,i))
         ENDDO
      ENDDO
      CALL GTStop(hcont); conttime = conttime + GTGetTime(hcont)

!
! 1D FFT in each node using the FFTW library
!
      CALL GTStart(hfft)
      CALL GPMANGLE(execute_dft)(plan%planc,inout,inout)
      CALL GTStop(hfft); ffttime = ffttime + GTGetTime(hfft)

      CALL GTStop(htot); tottime = tottime + GTGetTime(htot)

      RETURN
      END SUBROUTINE fftp1d_real_to_complex_z

!*****************************************************************
      SUBROUTINE fftp3d_complex_to_real(plan,in,out,comm)
!-----------------------------------------------------------------
!
! Computes the 3D complex-to-real FFT in parallel. The 
! complex input has the same structure than the input 
! of the 3D FFTW, but should be transposed. The real 
! output has the same order than the output of the FFTW.
! The input data is destroyed during the computation.
!
! Parameters
!     plan : the 3D plan created with fftp3d_create_plan [IN]
!     in   : complex input array [IN]
!     out  : real output array [OUT]
!     comm : the MPI communicator (handle) [IN]
!-----------------------------------------------------------------

      USE fprecision
      USE mpivars
      USE commtypes
      USE fftplans
      USE gtimer
!$    USE threads
      IMPLICIT NONE

      TYPE(FFTPLAN), INTENT(IN) :: plan

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(plan%nz,plan%ny,ista:iend) :: in 
      COMPLEX(KIND=GP), DIMENSION(ista:iend,plan%ny,plan%nz)           :: c1
      REAL(KIND=GP), INTENT(OUT), DIMENSION(plan%nx,plan%ny,ksta:kend) :: out

      INTEGER, DIMENSION(0:nprocs-1)      :: ireq1,ireq2
      INTEGER, DIMENSION(MPI_STATUS_SIZE) :: istatus
      INTEGER, INTENT(IN)                 :: comm
      INTEGER :: i,j,k
      INTEGER :: ii,jj,kk
      INTEGER :: irank
      INTEGER :: isendTo, igetFrom
      INTEGER :: istrip,iproc


      CALL GTStart(htot)

!
! 1D FFT in each node using the FFTW library
      CALL GTStart(hfft)
      CALL GPMANGLE(execute_dft)(plan%planc,in,in)
      CALL GTStop(hfft); ffttime = ffttime + GTGetTime(hfft)

      CALL GTStart(htra)
!
! Cache friendly transposition
!
!$omp parallel do if ((iend-ista)/csize.ge.nth) private (jj,kk,i,j,k)
      DO ii = ista,iend,csize
!$omp parallel do if ((iend-ista)/csize.lt.nth) private (kk,i,j,k)
         DO jj = 1,plan%ny,csize
            DO kk = 1,plan%nz,csize
               DO i = ii,min(iend,ii+csize-1)
               DO j = jj,min(plan%ny,jj+csize-1)
               DO k = kk,min(plan%nz,kk+csize-1)
                  c1(i,j,k) = in(k,j,i)
               END DO
               END DO
               END DO
            END DO
         END DO
      END DO
      CALL GTStop(htra); tratime = tratime + GTGetTime(htra)
!
! Transposes the result between nodes using 
! strip mining when nstrip>1 (rreddy@psc.edu)
!
      CALL GTStart(hcom)
      do iproc = 0, nprocs-1, nstrip
         do istrip=0, nstrip-1
            irank = iproc + istrip

            isendTo = myrank + irank
            if ( isendTo .ge. nprocs ) isendTo = isendTo - nprocs

            igetFrom = myrank - irank
            if ( igetFrom .lt. 0 ) igetFrom = igetFrom + nprocs
            CALL MPI_IRECV(plan%carr,1,plan%itype1(igetFrom),igetFrom, & 
                          1,comm,ireq2(irank),ierr)
            CALL MPI_ISEND(c1,1,plan%itype2(isendTo),isendTo, &
                          1,comm,ireq1(irank),ierr)
         enddo

         do istrip=0, nstrip-1
            irank = iproc + istrip
            CALL MPI_WAIT(ireq1(irank),istatus,ierr)
            CALL MPI_WAIT(ireq2(irank),istatus,ierr)
         enddo
      enddo
      CALL GTStop(hcom); comtime = comtime + GTGetTime(hcom)
!
! 2D FFT in each node using the FFTW library
!
      CALL GTStart(hfft)
      CALL GPMANGLE(execute_dft_c2r)(plan%planr,plan%carr,out)
      CALL GTStop(hfft); ffttime = ffttime + GTGetTime(hfft)
      
      CALL GTStop(htot); tottime = tottime + GTGetTime(htot)

      RETURN
      END SUBROUTINE fftp3d_complex_to_real

!*****************************************************************
      SUBROUTINE fftp2d_complex_to_real_xy(plan,in,out,comm)
!-----------------------------------------------------------------
!
! Computes the 3D complex-to-real FFT in parallel. The 
! complex input has the same structure than the input 
! of the 3D FFTW, but should be transposed. The real 
! output has the same order than the output of the FFTW.
! The input data is destroyed during the computation.
!
! Parameters
!     plan : the 3D plan created with fftp3d_create_plan [IN]
!     in   : complex input array f*(z,ky,kx) [IN]
!     out  : real output array f(x,y,z) [OUT]
!     comm : the MPI communicator (handle) [IN]
!-----------------------------------------------------------------

      USE fprecision
      USE mpivars
      USE commtypes
      USE fftplans
      USE gtimer
!$    USE threads
      IMPLICIT NONE

      TYPE(FFTPLAN), INTENT(IN) :: plan

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(plan%nz,plan%ny,ista:iend) :: in 
      COMPLEX(KIND=GP), DIMENSION(ista:iend,plan%ny,plan%nz)           :: c1
      REAL(KIND=GP), INTENT(OUT), DIMENSION(plan%nx,plan%ny,ksta:kend) :: out

      INTEGER, DIMENSION(0:nprocs-1)      :: ireq1,ireq2
      INTEGER, DIMENSION(MPI_STATUS_SIZE) :: istatus
      INTEGER, INTENT(IN)                 :: comm
      INTEGER :: i,j,k
      INTEGER :: ii,jj,kk
      INTEGER :: irank
      INTEGER :: isendTo, igetFrom
      INTEGER :: istrip,iproc

      CALL GTStart(htot)
      CALL GTStart(htra)
!
! Cache friendly transposition
!
!$omp parallel do if ((iend-ista)/csize.ge.nth) private (jj,kk,i,j,k)
      DO ii = ista,iend,csize
!$omp parallel do if ((iend-ista)/csize.lt.nth) private (kk,i,j,k)
         DO jj = 1,plan%ny,csize
            DO kk = 1,plan%nz,csize
               DO i = ii,min(iend,ii+csize-1)
               DO j = jj,min(plan%ny,jj+csize-1)
               DO k = kk,min(plan%nz,kk+csize-1)
                  c1(i,j,k) = in(k,j,i)
               END DO
               END DO
               END DO
            END DO
         END DO
      END DO
      CALL GTStop(htra); tratime = tratime + GTGetTime(htra)
!
! Transposes the result between nodes using 
! strip mining when nstrip>1 (rreddy@psc.edu)
!
      CALL GTStart(hcom)
      do iproc = 0, nprocs-1, nstrip
         do istrip=0, nstrip-1
            irank = iproc + istrip

            isendTo = myrank + irank
            if ( isendTo .ge. nprocs ) isendTo = isendTo - nprocs

            igetFrom = myrank - irank
            if ( igetFrom .lt. 0 ) igetFrom = igetFrom + nprocs
            CALL MPI_IRECV(plan%carr,1,plan%itype1(igetFrom),igetFrom, & 
                          1,comm,ireq2(irank),ierr)
            CALL MPI_ISEND(c1,1,plan%itype2(isendTo),isendTo, &
                          1,comm,ireq1(irank),ierr)
         enddo

         do istrip=0, nstrip-1
            irank = iproc + istrip
            CALL MPI_WAIT(ireq1(irank),istatus,ierr)
            CALL MPI_WAIT(ireq2(irank),istatus,ierr)
         enddo
      enddo
      CALL GTStop(hcom); comtime = comtime + GTGetTime(hcom)
!
! 2D FFT in each node using the FFTW library
!
      CALL GTStart(hfft)
      CALL GPMANGLE(execute_dft_c2r)(plan%planr,plan%carr,out)
      CALL GTStop(hfft); ffttime = ffttime + GTGetTime(hfft)
      
      CALL GTStop(htot); tottime = tottime + GTGetTime(htot)

      RETURN
      END SUBROUTINE fftp2d_complex_to_real_xy

!*****************************************************************
      SUBROUTINE fftp1d_complex_to_real_z(plan,inout,comm)
!-----------------------------------------------------------------
!
! Computes the 3D complex-to-real FFT in parallel. The 
! complex input has the same structure than the input 
! of the 3D FFTW, but should be transposed. The real 
! output has the same order than the output of the FFTW.
! The input data is destroyed during the computation.
!
! Parameters
!     plan : the 3D plan created with fftp3d_create_plan [IN]
!     in   : f*(kz,ky,kx) at input, f(z,ky,kx) at output [INOUT]
!     comm : the MPI communicator (handle) [IN]
!-----------------------------------------------------------------

      USE fprecision
      USE mpivars
      USE commtypes
      USE fftplans
      USE gtimer
!$    USE threads
      IMPLICIT NONE

      TYPE(FFTPLAN), INTENT(IN) :: plan
      COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(plan%nz,plan%ny,ista:iend) :: inout 
      INTEGER, INTENT(IN)                 :: comm

      INTEGER :: i,j,k

!
! 1D FFT in each node using the FFTW library
      CALL GTStart(htot)
      CALL GTStart(hfft)
      CALL GPMANGLE(execute_dft)(plan%planc,inout,inout)
      CALL GTStop(hfft); ffttime = ffttime + GTGetTime(hfft)
      CALL GTStop(htot); tottime = tottime + GTGetTime(htot)

      RETURN
      END SUBROUTINE fftp1d_complex_to_real_z

!*****************************************************************
      SUBROUTINE block3d(imin,imax,jmin,jmax,kmin,ista,iend, &
                        jsta,jend,ksta,kend,ioldtype,inewtype)
!-----------------------------------------------------------------
!
! Soubroutine for defining derived data types in 3D.
!
! Parameters
!     imin : the minimum value in the first dimension [IN]
!     imax : the maximum value in the first dimension [IN]
!     jmin : the minimum value in the second dimension [IN]
!     jmax : the maximum value in the second dimension [IN]
!     kmin : the minimum value in the third dimension [IN]
!     ista : start value of the block in the first dimension [IN]
!     iend : end value of the block in the first dimension [IN]
!     jsta : start value of the block in the second dimension [IN]
!     jend : end value of the block in the second dimension [IN]
!     ksta : start value of the block in the third dimension [IN]
!     kend : end value of the block in the third dimension [IN]
!     ioldtype: data type of the elements in the block [IN]
!     inewtype: the derived data type for the block [OUT]
!-----------------------------------------------------------------

      USE commtypes
      USE fftplans
      IMPLICIT NONE

      INTEGER, DIMENSION (2) :: iblock,idisp,itype

      INTEGER, INTENT(IN)  :: ista,iend
      INTEGER, INTENT(IN)  :: jsta,jend
      INTEGER, INTENT(IN)  :: ksta,kend
      INTEGER, INTENT(IN)  :: imin,imax
      INTEGER, INTENT(IN)  :: jmin,jmax,kmin
      INTEGER, INTENT(IN)  :: ioldtype
      INTEGER, INTENT(OUT) :: inewtype

      INTEGER :: ilen,jlen,klen
      INTEGER :: isize,idist
      INTEGER :: itemp,itemp2
      INTEGER :: ierr

!      INTEGER :: ierr1,nerr=1024
!      CHARACTER*1024 serr
!      
!      nerr = 1024

      CALL MPI_TYPE_EXTENT(ioldtype,isize,ierr)
      ilen = iend-ista+1
      jlen = jend-jsta+1
      klen = kend-ksta+1
      CALL MPI_TYPE_VECTOR(jlen,ilen,imax-imin+1,ioldtype,itemp,ierr)
      idist = (imax-imin+1)*(jmax-jmin+1)*isize
      CALL MPI_TYPE_HVECTOR(klen,1,idist,itemp,itemp2,ierr)
      CALL MPI_TYPE_FREE(itemp,ierr)
      iblock(1) = 1
      iblock(2) = 1
      idisp(1) = 0
      idisp(2) = ((imax-imin+1)*(jmax-jmin+1)*(ksta-kmin) &
                 +(imax-imin+1)*(jsta-jmin)+(ista-imin))*isize
      itype(1) = MPI_LB
      itype(2) = itemp2
      CALL MPI_TYPE_STRUCT(2,iblock,idisp,itype,inewtype,ierr)

      CALL MPI_TYPE_FREE(itemp2,ierr)
      CALL MPI_TYPE_COMMIT(inewtype,ierr)


      RETURN
      END SUBROUTINE block3d
