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
!    Dec 2020: Updated calls to MPI3 and MPI4 practices. (P. Mininni)
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
      SUBROUTINE fftp3d_create_plan_rc(plan,n,flags)
!-----------------------------------------------------------------
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

      ! Create XY -> Z plan
      CALL GPMANGLE(plan_many_dft_r2c)(plan%planrxy,2,(/n(1),n(2)/),          &
                         kend-ksta+1,plan%rarr,                               &
                         (/n(1),n(2)*(kend-ksta+1)/),1,n(1)*n(2),             &
                         plan%carr,(/n(1)/2+1,n(2)*(kend-ksta+1)/),1,         &
                         (n(1)/2+1)*n(2),flags)
      CALL GPMANGLE(plan_many_dft)(plan%plancz,1,n(3),n(2)*(iend-ista+1),     &
                         plan%ccarr,(iend-ista+1)*n(2)*n(3),1,n(3),           &
                         plan%ccarr,(iend-ista+1)*n(2)*n(3),1,n(3),           &
                         -1,flags)

      ! Create X -> YZ plan
      CALL GPMANGLE(plan_many_dft_r2c)(plan%planrx,1,n(1),                    &
                         n(2)*(kend-ksta+1),plan%rarr,n(1)*n(2)*(kend-ksta+1),&
                         1,n(1),plan%carr,(n(1)/2+1)*n(2)*(kend-ksta+1),1,    &
                         (n(1)/2+1),flags)
      CALL GPMANGLE(plan_many_dft)(plan%plancyz,2,(/n(3),n(2)/),(iend-ista+1),&
                         plan%ccarr,(/n(3),n(2)*(iend-ista+1)/),1,n(3)*n(2),  &
                         plan%ccarr,(/n(3),n(2)*(iend-ista+1)/),1,n(3)*n(2),  &
                         -1,flags)
      plan%nx = n(1)
      plan%ny = n(2)
      plan%nz = n(3)
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
      END SUBROUTINE fftp3d_create_plan_rc

!*****************************************************************
      SUBROUTINE fftp3d_create_plan_cr(plan,n,flags)
!-----------------------------------------------------------------
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

      ! Create Z -> YX plan
      CALL GPMANGLE(plan_many_dft_c2r)(plan%planrxy,2,(/n(1),n(2)/),          &
                         kend-ksta+1,plan%carr,                               &
                         (/n(1)/2+1,n(2)*(kend-ksta+1)/),1,                   &
                         (n(1)/2+1)*n(2),plan%rarr,                           &
                         (/n(1),n(2)*(kend-ksta+1)/),1,n(1)*n(2),flags)

      CALL GPMANGLE(plan_many_dft)(plan%plancz,1,n(3),n(2)*(iend-ista+1),     &
                         plan%ccarr,(iend-ista+1)*n(2)*n(3),1,n(3),           &
                         plan%ccarr,(iend-ista+1)*n(2)*n(3),1,n(3),           &
                         1,flags)

      ! Create ZY -> X plan
      CALL GPMANGLE(plan_many_dft_c2r)(plan%planrx,1,n(1),n(2)*(kend-ksta+1), &
                         plan%carr,(n(1)/2+1)*n(2)*(kend-ksta+1),1,n(1)/2+1,  &
                         plan%rarr,n(1)*n(2)*(kend-ksta+1),1,n(1),flags)

      CALL GPMANGLE(plan_many_dft)(plan%plancyz,2,(/n(3),n(2)/),(iend-ista+1),&
                         plan%ccarr,(/n(3),n(2)*(iend-ista+1)/),1,n(3)*n(2),  &
                         plan%ccarr,(/n(3),n(2)*(iend-ista+1)/),1,n(3)*n(2),  &
                         1,flags)


      plan%nx = n(1)
      plan%ny = n(2)
      plan%nz = n(3)

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
! Destroys FFTW plans in each node.
!
! Parameters
!     plan : the parallel 3D plan [INOUT]
!-----------------------------------------------------------------

      USE fftplans
      USE gtimer
      IMPLICIT NONE

      TYPE(FFTPLAN), INTENT(INOUT) :: plan

      CALL GPMANGLE(destroy_plan)(plan%planrxy)
      CALL GPMANGLE(destroy_plan)(plan%planrx)
      CALL GPMANGLE(destroy_plan)(plan%plancyz)
      CALL GPMANGLE(destroy_plan)(plan%plancz)

      DEALLOCATE( plan%ccarr  )
      DEALLOCATE( plan%carr   )
      DEALLOCATE( plan%rarr   )
      DEALLOCATE( plan%itype1 )
      DEALLOCATE( plan%itype2 )

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
! Computes the 3D real-to-complex FFT in parallel. The
! complex output has the same structure than the output
! of the 3D FFTW, but the output is transposed. If the input
! is non-periodic FC-Gram periodic extensions are appropriately
! computed.
!
! Parameters
!     plan : the FCPLAN plan (see fcgram_mod.fpp) [IN]
!     in   : real input array [IN]
!     out  : complex output array [OUT]
!     comm : the MPI communicator (handle) [IN]
!-----------------------------------------------------------------

      USE commtypes
      USE fprecision
      USE mpivars
      USE fcgram
      IMPLICIT NONE

      TYPE(FCPLAN), INTENT(IN) :: plan

      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(plan%z%n,plan%y%n,ista:iend) :: out
      REAL(KIND=GP), INTENT(IN), DIMENSION(plan%x%n,plan%y%n,ksta:kend) :: in

      INTEGER, INTENT(IN)                 :: comm

      IF (plan%y%C .gt. 0) THEN
         CALL fftp1d_real_to_complex_x(plan,in,out,comm)
         CALL fftp2d_real_to_complex_yz(plan,out,comm)
      ELSE
         CALL fftp2d_real_to_complex_xy(plan,in,out,comm)
         CALL fftp1d_real_to_complex_z(plan,out,comm)
      ENDIF

      RETURN
      END SUBROUTINE fftp3d_real_to_complex

!*****************************************************************
      SUBROUTINE fftp2d_real_to_complex_xy(plan,in,out,comm)
!-----------------------------------------------------------------
! Computes a local 2D real-to-complex FFT in the first two indices
! and transposes the result. If the input has shape (nx,ny,nz) the
! output's shape is (nz,ny,nx/2+1).
!
! Parameters
!     plan : the FCPLAN plan (see fcgram_mod.fpp) [IN]
!     in   : real input array [IN]
!     out  : complex output array [OUT]
!     comm : the MPI communicator (handle) [IN]
!-----------------------------------------------------------------

      USE commtypes
      USE fprecision
      USE mpivars
      USE fcgram
      USE gtimer
!$    USE threads
      IMPLICIT NONE

      TYPE(FCPLAN), INTENT(IN) :: plan

      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(plan%z%n,plan%y%n,ista:iend) :: out
      COMPLEX(KIND=GP), DIMENSION(ista:iend,plan%y%n,plan%z%n)          :: c1
      REAL(KIND=GP), INTENT(IN), DIMENSION(plan%x%n,plan%y%n,ksta:kend) :: IN

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
      CALL GPMANGLE(execute_dft_r2c)(plan%planrc%planrxy,in,plan%planrc%carr)
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
            CALL MPI_IRECV(c1,1,plan%planrc%itype2(igetFrom),igetFrom,      & 
                          1,comm,ireq2(irank),ierr)

            CALL MPI_ISEND(plan%planrc%carr,1,plan%planrc%itype1(isendTo),isendTo, &
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
         DO jj = 1,plan%y%n,csize
            DO kk = 1,plan%z%n,csize
               DO i = ii,min(iend,ii+csize-1)
               DO j = jj,min(plan%y%n,jj+csize-1)
               DO k = kk,min(plan%z%n,kk+csize-1)
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
      SUBROUTINE fftp2d_real_to_complex_yz(plan,inout,comm)
!-----------------------------------------------------------------
! Computes the 2D forward FFT along the first two indices of the
! input array. In the case the input array is non-periodic, a 
! suitable FC-Gram periodic extension is computed. The operation
! is performed inplace.
!
! Parameters
!     plan : the FCPLAN plan (see fcgram_mod.fpp) [IN]
!     in   : contains f(z,y,kx) as input and returns f(kz,ky,kx) [INOUT]
!     comm : the MPI communicator (handle) [IN]
!-----------------------------------------------------------------

      USE commtypes
      USE fprecision
      USE mpivars
      USE gtimer
      USE fcgram
!$    USE threads
      IMPLICIT NONE

      TYPE(FCPLAN), INTENT(IN) :: plan
      COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(plan%z%n,plan%y%n,ista:iend) :: inout
      INTEGER, INTENT(IN)                 :: comm

      INTEGER :: ny,Cy,dy,nz,Cz,dz
      INTEGER :: i,j,ii,jj

      ny = plan%y%n
      Cy = plan%y%C
      dy = plan%y%d

      nz = plan%z%n
      Cz = plan%z%C
      dz = plan%z%d

      CALL GTStart(htot)
      CALL GTStart(hcont)

      ! Continuation in Y direction
!$omp paralleldo if ((iend-ista)/csize.ge.nth) private (j,ii,jj)
      DO i= ista,iend
!$omp parallel do if ((iend-ista)/csize.lt.nth) private (ii,jj)
      DO j =1,nz-Cz
         DO ii=1,Cy
            inout(j,ny-Cy+ii,i) = plan%y%dir(ii,1)*inout(j,ny-Cy-dy+1,i) +&
                                plan%y%dir(Cy-ii+1,1)*inout(j,dy,i)
         ENDDO
         DO jj=2,dy
            DO ii=1,Cy
               inout(j,ny-Cy+ii,i) = inout(j,ny-Cy+ii,i) +&
                                   plan%y%dir(ii,jj)*inout(j,ny-Cy-dy+jj,i) +&
                                   plan%y%dir(Cy-ii+1,jj)*inout(j,dy-jj+1,i)
            ENDDO
         ENDDO
      ENDDO
      ENDDO


!$omp paralleldo if ((iend-ista)/csize.ge.nth) private (j,ii,jj)
      DO i= ista,iend
!$omp parallel do if ((iend-ista)/csize.lt.nth) private (ii,jj)
      DO j =1,ny
         DO ii=1,Cz
            inout(nz-Cz+ii,j,i) = plan%z%dir(ii,1)*inout(nz-Cz-dz+1,j,i) +&
                                  plan%z%dir(Cz-ii+1,1)*inout(dz,j,i)
         ENDDO
         DO jj=2,dz
            DO ii=1,Cz
               inout(nz-Cz+ii,j,i) = inout(nz-Cz+ii,j,i) +&
                                   plan%z%dir(ii,jj)*inout(nz-Cz-dz+jj,j,i) +&
                                   plan%z%dir(Cz-ii+1,jj)*inout(dz-jj+1,j,i)
            ENDDO
         ENDDO
      ENDDO
      ENDDO

      CALL GTStop(hcont); conttime = conttime + GTGetTime(hcont)

!
! 2D FFT in each node using the FFTW library
!
      CALL GTStart(hfft)
      CALL GPMANGLE(execute_dft)(plan%planrc%plancyz,inout,inout)
      CALL GTStop(hfft); ffttime = ffttime + GTGetTime(hfft)

      CALL GTStop(htot); tottime = tottime + GTGetTime(htot)

      RETURN
      END SUBROUTINE fftp2d_real_to_complex_yz


!*****************************************************************
      SUBROUTINE fftp1d_real_to_complex_x(plan,in,out,comm)
!-----------------------------------------------------------------
! Computes the 1D real-to-complex FFT along the first index of
! the input array and transposes the result. If the input has 
! shape (nx,ny,nz) the output's shape is (nz,ny,nx/2+1).
!
! Parameters
!     plan : the FCPLAN plan (see fcgram_mod.fpp) [IN]
!     in   : real input array [IN]
!     out  : complex output array [OUT]
!     comm : the MPI communicator (handle) [IN]
!-----------------------------------------------------------------

      USE commtypes
      USE fprecision
      USE mpivars
      USE fcgram
      USE gtimer
!$    USE threads
      IMPLICIT NONE

      TYPE(FCPLAN), INTENT(IN) :: plan

      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(plan%z%n,plan%y%n,ista:iend) :: out
      COMPLEX(KIND=GP), DIMENSION(ista:iend,plan%y%n,plan%z%n)          :: c1
      REAL(KIND=GP), INTENT(IN), DIMENSION(plan%x%n,plan%y%n,ksta:kend) :: IN

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
! 1D FFT in each node using the FFTW library
!
      CALL GTStart(hfft)
      CALL GPMANGLE(execute_dft_r2c)(plan%planrc%planrx,in,plan%planrc%carr)
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
            CALL MPI_IRECV(c1,1,plan%planrc%itype2(igetFrom),igetFrom,      & 
                          1,comm,ireq2(irank),ierr)

            CALL MPI_ISEND(plan%planrc%carr,1,plan%planrc%itype1(isendTo),isendTo, &
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
         DO jj = 1,plan%y%n,csize
            DO kk = 1,plan%z%n,csize
               DO i = ii,min(iend,ii+csize-1)
               DO j = jj,min(plan%y%n,jj+csize-1)
               DO k = kk,min(plan%z%n,kk+csize-1)
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
      END SUBROUTINE fftp1d_real_to_complex_x


!*****************************************************************
      SUBROUTINE fftp1d_real_to_complex_z(plan,inout,comm)
!-----------------------------------------------------------------
! Computes the 1D forward FFT along the first index of the input
! array. In the case the input array is non-periodic, a suitable 
! FC-Gram periodic extension is computed. The operation is
! performed inplace.
!
! Parameters
!     plan : the FCPLAN plan (see fcgram_mod.fpp) [IN]
!     in   : contains f(z,ky,kx) as input and returns f(kz,ky,kx) [INOUT]
!     comm : the MPI communicator (handle) [IN]
!-----------------------------------------------------------------

      USE commtypes
      USE fprecision
      USE mpivars
      USE gtimer
      USE fcgram
!$    USE threads
      IMPLICIT NONE

      TYPE(FCPLAN), INTENT(IN) :: plan
      COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(plan%z%n,plan%y%n,ista:iend) :: inout
      INTEGER, INTENT(IN)                 :: comm

      INTEGER :: n,C,d
      INTEGER :: i,j,ii,jj

      n = plan%z%n
      C = plan%z%C
      d = plan%z%d

      CALL GTStart(htot)
      CALL GTStart(hcont)

      IF (plan%z%C .ge. 0) THEN
!$omp paralleldo if ((iend-ista)/csize.ge.nth) private (j,ii,jj)
      DO i= ista,iend
!$omp parallel do if ((iend-ista)/csize.lt.nth) private (ii,jj)
      DO j =1,plan%y%n
         DO ii=1,C
            inout(n-C+ii,j,i) = plan%z%dir(ii,1)*inout(n-C-d+1,j,i) +&
                                plan%z%dir(C-ii+1,1)*inout(d,j,i)
         ENDDO
         DO jj=2,d
            DO ii=1,C
               inout(n-C+ii,j,i) = inout(n-C+ii,j,i) +&
                                   plan%z%dir(ii,jj)*inout(n-C-d+jj,j,i) +&
                                   plan%z%dir(C-ii+1,jj)*inout(d-jj+1,j,i)
            ENDDO
         ENDDO

! Old approach
!            ! Continuation using last o points
!            inout(n-C+1:n,j,i) = MATMUL(planfc%z%dir, inout(n-C-d+1:n-C,j,i))
!
!            ! Continuation using first o points
!            inout(n-C+1:n,j,i) = inout(n-C+1:n,j, i) + &
!                MATMUL(planfc%z%dir(C:1:-1,:), inout(d:1:-1,j,i))
      ENDDO
      ENDDO
      ENDIF
      CALL GTStop(hcont); conttime = conttime + GTGetTime(hcont)

!
! 1D FFT in each node using the FFTW library
!
      CALL GTStart(hfft)
      CALL GPMANGLE(execute_dft)(plan%planrc%plancz,inout,inout)
      CALL GTStop(hfft); ffttime = ffttime + GTGetTime(hfft)

      CALL GTStop(htot); tottime = tottime + GTGetTime(htot)

      RETURN
      END SUBROUTINE fftp1d_real_to_complex_z

!*****************************************************************
      SUBROUTINE fftp3d_complex_to_real(plan,in,out,comm)
!-----------------------------------------------------------------
! Computes the 3D complex-to-real FFT in parallel. The 
! complex input has the same structure than the input 
! of the 3D FFTW, but should be transposed. The real 
! output has the same order than the output of the FFTW.
! The input data is destroyed during the computation.
!
! Parameters
!     plan : the FCPLAN plan (see fcgram_mod.fpp) [IN]
!     in   : complex input array [IN]
!     out  : real output array [OUT]
!     comm : the MPI communicator (handle) [IN]
!
!-----------------------------------------------------------------

      USE fprecision
      USE mpivars
      USE fcgram
      IMPLICIT NONE

      TYPE(FCPLAN), INTENT(IN) :: plan

      COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(plan%z%n,plan%y%n,ista:iend) :: in 
      REAL(KIND=GP), INTENT(OUT), DIMENSION(plan%x%n,plan%y%n,ksta:kend) :: out

      INTEGER, INTENT(IN)                 :: comm

      CALL fftp1d_complex_to_real_z(plan,in,comm)
      CALL fftp2d_complex_to_real_xy(plan,in,out,comm)
     
      RETURN
      END SUBROUTINE fftp3d_complex_to_real

!*****************************************************************
      SUBROUTINE fftp2d_complex_to_real_xy(plan,in,out,comm)
!-----------------------------------------------------------------
! Computes the 2D real-to-complex FFT along the last couple of
! indices of the input array and transposes the result.
! If the input has shape (nz,ny,nx/2+1) the output's shape is
! (nx,ny,nz). The input data is destroyed during the computation.
!
! Parameters
!     plan : the FCPLAN plan (see fcgram_mod.fpp) [IN]
!     in   : complex input array f*(z,ky,kx) [IN]
!     out  : real output array f(x,y,z) [OUT]
!     comm : the MPI communicator (handle) [IN]
!-----------------------------------------------------------------

      USE fprecision
      USE mpivars
      USE commtypes
      USE fcgram
      USE gtimer
!$    USE threads
      IMPLICIT NONE

      TYPE(FCPLAN), INTENT(IN) :: plan

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(plan%z%n,plan%y%n,ista:iend) :: in 
      COMPLEX(KIND=GP), DIMENSION(ista:iend,plan%y%n,plan%z%n)           :: c1
      REAL(KIND=GP), INTENT(OUT), DIMENSION(plan%x%n,plan%y%n,ksta:kend) :: out

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
         DO jj = 1,plan%y%n,csize
            DO kk = 1,plan%z%n,csize
               DO i = ii,min(iend,ii+csize-1)
               DO j = jj,min(plan%y%n,jj+csize-1)
               DO k = kk,min(plan%z%n,kk+csize-1)
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
            CALL MPI_IRECV(plan%plancr%carr,1,plan%plancr%itype1(igetFrom),igetFrom, & 
                          1,comm,ireq2(irank),ierr)
            CALL MPI_ISEND(c1,1,plan%plancr%itype2(isendTo),isendTo, &
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
      CALL GPMANGLE(execute_dft_c2r)(plan%plancr%planrxy,plan%plancr%carr,out)
      CALL GTStop(hfft); ffttime = ffttime + GTGetTime(hfft)
      
      CALL GTStop(htot); tottime = tottime + GTGetTime(htot)

      RETURN
      END SUBROUTINE fftp2d_complex_to_real_xy

!*****************************************************************
      SUBROUTINE fftp2d_complex_to_real_yz(plan,inout,comm)
!-----------------------------------------------------------------
! Computes the 2D backward FFT along the first couple of indices 
! of the input array. The operation is performed inplace.
!
! Parameters
!     plan : the FCPLAN plan (see fcgram_mod.fpp) [IN]
!     in   : f*(kz,ky,kx) at input, f(z,ky,kx) at output [INOUT]
!     comm : the MPI communicator (handle) [IN]
!-----------------------------------------------------------------

      USE fprecision
      USE mpivars
      USE commtypes
      USE fcgram
      USE gtimer
!$    USE threads
      IMPLICIT NONE

      TYPE(FCPLAN), INTENT(IN) :: plan
      COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(plan%z%n,plan%y%n,ista:iend) :: inout 
      INTEGER, INTENT(IN)                 :: comm

      INTEGER :: i,j,k

!
! 1D FFT in each node using the FFTW library
      CALL GTStart(htot)
      CALL GTStart(hfft)
      CALL GPMANGLE(execute_dft)(plan%plancr%plancyz,inout,inout)
      CALL GTStop(hfft); ffttime = ffttime + GTGetTime(hfft)
      CALL GTStop(htot); tottime = tottime + GTGetTime(htot)

      RETURN
      END SUBROUTINE fftp2d_complex_to_real_yz



!*****************************************************************
      SUBROUTINE fftp1d_complex_to_real_x(plan,in,out,comm)
!-----------------------------------------------------------------
! Computes the 1D real-to-complex FFT along the last index
! of the input array and transposes the result. If the input has
! shape (nz,ny,nx/2+1) the output's shape is  (nx,ny,nz).
! The input data is destroyed during the computation.
!
! Parameters
!     plan : the FCPLAN plan (see fcgram_mod.fpp) [IN]
!     in   : complex input array f(z,y,kx) [IN]
!     out  : real output array f(x,y,z) [OUT]
!     comm : the MPI communicator (handle) [IN]
!-----------------------------------------------------------------

      USE fprecision
      USE mpivars
      USE commtypes
      USE fcgram
      USE gtimer
!$    USE threads
      IMPLICIT NONE

      TYPE(FCPLAN), INTENT(IN) :: plan

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(plan%z%n,plan%y%n,ista:iend) :: in 
      COMPLEX(KIND=GP), DIMENSION(ista:iend,plan%y%n,plan%z%n)           :: c1
      REAL(KIND=GP), INTENT(OUT), DIMENSION(plan%x%n,plan%y%n,ksta:kend) :: out

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
         DO jj = 1,plan%y%n,csize
            DO kk = 1,plan%z%n,csize
               DO i = ii,min(iend,ii+csize-1)
               DO j = jj,min(plan%y%n,jj+csize-1)
               DO k = kk,min(plan%z%n,kk+csize-1)
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
            CALL MPI_IRECV(plan%plancr%carr,1,plan%plancr%itype1(igetFrom),igetFrom, & 
                          1,comm,ireq2(irank),ierr)
            CALL MPI_ISEND(c1,1,plan%plancr%itype2(isendTo),isendTo, &
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
! 1D FFT in each node using the FFTW library
!
      CALL GTStart(hfft)
      CALL GPMANGLE(execute_dft_c2r)(plan%plancr%planrx,plan%plancr%carr,out)
      CALL GTStop(hfft); ffttime = ffttime + GTGetTime(hfft)
      
      CALL GTStop(htot); tottime = tottime + GTGetTime(htot)

      RETURN
      END SUBROUTINE fftp1d_complex_to_real_x


!*****************************************************************
      SUBROUTINE fftp1d_complex_to_real_z(plan,inout,comm)
!-----------------------------------------------------------------
! Computes the 1D backward FFT along the first index of the input
! array. The operation is performed inplace.
!
! Parameters
!     plan : the FCPLAN plan (see fcgram_mod.fpp) [IN]
!     in   : f*(kz,ky,kx) at input, f(z,ky,kx) at output [INOUT]
!     comm : the MPI communicator (handle) [IN]
!-----------------------------------------------------------------

      USE fprecision
      USE mpivars
      USE commtypes
      USE fcgram
      USE gtimer
!$    USE threads
      IMPLICIT NONE

      TYPE(FCPLAN), INTENT(IN) :: plan
      COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(plan%z%n,plan%y%n,ista:iend) :: inout 
      INTEGER, INTENT(IN)                 :: comm

      INTEGER :: i,j,k

!
! 1D FFT in each node using the FFTW library
      CALL GTStart(htot)
      CALL GTStart(hfft)
      CALL GPMANGLE(execute_dft)(plan%plancr%plancz,inout,inout)
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

      INTEGER, INTENT(IN)  :: ista,iend
      INTEGER, INTENT(IN)  :: jsta,jend
      INTEGER, INTENT(IN)  :: ksta,kend
      INTEGER, INTENT(IN)  :: imin,imax
      INTEGER, INTENT(IN)  :: jmin,jmax,kmin
      INTEGER, INTENT(IN)  :: ioldtype
      INTEGER, INTENT(OUT) :: inewtype

      INTEGER(KIND=MPI_ADDRESS_KIND) :: ilb,isize,idist
      INTEGER :: ilen,jlen,klen
      INTEGER :: itemp,itemp2
      INTEGER :: ierr

      CALL MPI_TYPE_GET_EXTENT(ioldtype,ilb,isize,ierr)
      ilen = iend-ista+1
      jlen = jend-jsta+1
      klen = kend-ksta+1
      CALL MPI_TYPE_VECTOR(jlen,ilen,imax-imin+1,ioldtype,itemp,ierr)
      idist = (imax-imin+1)*(jmax-jmin+1)*isize
      CALL MPI_TYPE_CREATE_HVECTOR(klen,1,idist,itemp,itemp2,ierr)
      CALL MPI_TYPE_FREE(itemp,ierr)
      idist = ((imax-imin+1)*(jmax-jmin+1)*(ksta-kmin) &
              +(imax-imin+1)*(jsta-jmin)+(ista-imin))*isize
      CALL MPI_TYPE_CREATE_STRUCT(1,1,idist,itemp2,inewtype,ierr)
      CALL MPI_TYPE_FREE(itemp2,ierr)
      CALL MPI_TYPE_COMMIT(inewtype,ierr)

      RETURN
      END SUBROUTINE block3d

!*****************************************************************
      SUBROUTINE range(n1,n2,nprocs,irank,ista,iend)
!-----------------------------------------------------------------
! Soubroutine for computing the local coordinate range 
! when splitting the original array into the nodes
!
! Parameters
!     n1     : the minimum value in the splitted dimension [IN]
!     n2     : the maximum value in the splitted dimension [IN]
!     nprocs : the number of processors [IN]
!     irank  : the rank of the processor [IN]
!     ista   : start value for the local coordinate [OUT]
!     iend   : end value for the local coordinate [OUT]
!-----------------------------------------------------------------

      IMPLICIT NONE

      INTEGER, INTENT(IN)  :: n1,n2
      INTEGER, INTENT(IN)  :: nprocs,irank
      INTEGER, INTENT(OUT) :: ista,iend

      INTEGER :: myrank
      INTEGER :: iwork1,iwork2

      iwork1 = (n2-n1+1)/nprocs
      iwork2 = MOD(n2-n1+1,nprocs)
      ista = irank*iwork1+n1+MIN(irank,iwork2)
      iend = ista+iwork1-1
      IF (iwork2.gt.irank) iend = iend+1

      RETURN
      END SUBROUTINE range
