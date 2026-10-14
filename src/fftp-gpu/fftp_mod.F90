!=================================================================
! MODULES for FFTP-GPU
! Parallel Fast Fourier Transform in 3D on GPUs
!
! 2007 Pablo D. Mininni.
!      National Center for Atmospheric Research.
!      e-mail: mininni@ucar.uba.ar 
! 2026 GPU version with cuFFT/hipFFT (P.D. Mininni)
!
! Versión de SPECTER del fftp_mod.F90 de fftp-gpu de GHOST: el plan
! es el de src/fftp (los planes de FFTW planrxy y plancz como enteros
! de kind ikind, y el cronómetro de la continuación) más lo que
! agrega GHOST para el dispositivo, con los planes del dispositivo
! llamados dplanrxy y dplancz.
!=================================================================

!=================================================================

  MODULE fftplans
!
! Set the variable ikind to:  4 in 32 bits machines
!                             8 in 64 bits machines
! Set the variable csize to:  8 if L1 cache is <= 64 kb
!                            16 if L1 cache is 128 kb
! The variable nstrip controls strip mining during the 
! transposition. Often set to 1.
!
      USE ISO_C_BINDING
      USE fprecision
      USE mpi_f08, only   :  MPI_Datatype

      IMPLICIT NONE

      INCLUDE 'fftw3.f'
#include "gfft_defs.h"
 
      INTEGER, PARAMETER  :: ikind = IKIND_
      INTEGER, PARAMETER  :: csize = CSIZE_
      INTEGER, PARAMETER  :: nstrip = NSTRIP_
      INTEGER, PARAMETER  :: FFTW_REAL_TO_COMPLEX = FFTW_FORWARD
      INTEGER, PARAMETER  :: FFTW_COMPLEX_TO_REAL = FFTW_BACKWARD
      INTEGER             :: hcom,hfft,hmem,htra,hcont,htot
      DOUBLE PRECISION    :: comtime = 0.0
      DOUBLE PRECISION    :: ffttime = 0.0
      DOUBLE PRECISION    :: tratime = 0.0
      DOUBLE PRECISION    :: conttime = 0.0
      DOUBLE PRECISION    :: tottime = 0.0
      TYPE FFTPLAN
         COMPLEX(KIND=GP), DIMENSION (:,:,:), POINTER  :: ccarr
         COMPLEX(KIND=GP), DIMENSION (:,:,:), POINTER  :: carr
         REAL(KIND=GP),    DIMENSION (:,:,:), POINTER  :: rarr
         TYPE(MPI_Datatype), DIMENSION (:), POINTER :: itype1, itype2
         INTEGER(kind=ikind) :: planrxy,plancz      ! host (FFTW) plans
         GFFT_HANDLE_T       :: dplanrxy,dplancz    ! device plans
         INTEGER :: nx,ny,nz
         INTEGER :: fftdir
         ! Exchange buffer of the host path: the shared gc1 when the
         ! plan has the layout of the first plan created, otherwise a
         ! buffer of its own (host only, see fftp3d_create_plan)
         COMPLEX(KIND=GP), DIMENSION (:,:,:), POINTER :: c1 => NULL()
         LOGICAL :: ownc1 = .FALSE.
      END TYPE FFTPLAN
!
! Buffers of the device path, shared by all plans (the forward and
! the backward plans are never used at the same time) and resident on
! the device: the output of the 2D transform (gcarr), the transposed
! slabs after the exchange (gc1) and the packed messages (gsbuf).
! The host copy of gc1 is also the exchange buffer of the host path
! for the plans with the layout of the first plan created; a plan
! with another grid or task count (e.g. the two grids of BOOTS) gets
! its own host buffer (plan%c1) and can only use the host path.
! The tables hold the kx ranges (gioff/gilen) and z ranges (gkoff/gklen)
! of every task, and the offsets and counts of the packed send
! blocks (gsoff/gscnt) and of the contiguous receive blocks in gc1
! (groff/grcnt). gexchange selects the exchange between tasks:
!   1: device buffers passed to MPI (GPU-aware MPI, default)
!   2: staged through the host copies (for MPI libraries that are
!      not GPU-aware; environment variable GHOST_GPU_EXCHANGE=staged)
      COMPLEX(KIND=GP), DIMENSION (:,:,:), ALLOCATABLE, TARGET :: gcarr,gc1
      COMPLEX(KIND=GP), DIMENSION (:)    , ALLOCATABLE, TARGET :: gsbuf
      INTEGER, DIMENSION (:), ALLOCATABLE :: gioff,gilen,gkoff,gklen
      INTEGER, DIMENSION (:), ALLOCATABLE :: gsoff,gscnt,groff,grcnt
      INTEGER :: gexchange = 1
      LOGICAL :: gbuffers_ready = .FALSE.
      SAVE

  END MODULE fftplans
!=================================================================

  MODULE threads
!
! nth: number of threads used in OpenMP loops and FFTs
      INTEGER :: nth
      SAVE

  END MODULE threads
!=================================================================

  MODULE mpivars
      ! pkend is last k-index corresponding to the
      ! actual physical (non-continuated) domain.
      INTEGER, SAVE :: ista,iend
      INTEGER, SAVE :: jsta,jend
      INTEGER, SAVE :: ksta,kend,pkend
      INTEGER, SAVE :: nprocs,myrank
      INTEGER, SAVE :: provided
      INTEGER, SAVE :: ierr

  END MODULE mpivars
!=================================================================
