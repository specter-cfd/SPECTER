!=================================================================
! MODULES for FFTP v3
! Parallel Fast Fourier Transform in 3D using CUDA
!
! 2011 Duane L. Rosenberg and Pablo D. Mininni.
!      National Center for Atmospheric Research.
!      e-mail: mininni@ucar.uba.ar 
!
! 2021 FC-Gram support (M. Fontana)
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
      USE fprecision
      USE iso_c_binding

      IMPLICIT NONE
 
      INTEGER, PARAMETER  :: ikind = IKIND_
      INTEGER, PARAMETER  :: csize = CSIZE_
      INTEGER, PARAMETER  :: nstrip = NSTRIP_
      INTEGER, PARAMETER  :: nstreams = NSTREAMS_
      INTEGER, PARAMETER  :: FFTCU_REAL_TO_COMPLEX = -1
      INTEGER, PARAMETER  :: FFTCU_COMPLEX_TO_REAL =  1
      INTEGER, PARAMETER  :: FFTW_REAL_TO_COMPLEX = -1
      INTEGER, PARAMETER  :: FFTW_COMPLEX_TO_REAL =  1
      INTEGER, PARAMETER  :: FFTW_MEASURE =  0
      INTEGER, PARAMETER  :: FFTW_PATIENT =  0
      INTEGER, PARAMETER  :: FFTW_ESTIMATE=  0
      INTEGER, SAVE       :: streams_created = 0
      TYPE  (C_PTR)       :: pstream_(nstreams)
      INTEGER, DIMENSION (nstreams) :: issta,issnd
      INTEGER, DIMENSION (nstreams) :: kssta,kssnd
      INTEGER             :: hcom,hfft,hmem,htra,hcont,htot
      DOUBLE PRECISION    :: comtime  = 0.0
      DOUBLE PRECISION    :: ffttime  = 0.0
      DOUBLE PRECISION    :: memtime  = 0.0
      DOUBLE PRECISION    :: tratime  = 0.0
      DOUBLE PRECISION    :: conttime = 0.0
      DOUBLE PRECISION    :: tottime  = 0.0
      TYPE FFTPLAN
         COMPLEX(KIND=GP), DIMENSION (:,:,:), POINTER  :: ccarr
         COMPLEX(KIND=GP), DIMENSION (:,:,:), POINTER  :: carr
         REAL(KIND=GP),    DIMENSION (:,:,:), POINTER  :: rarr
         INTEGER, DIMENSION (:), POINTER               :: itype1, itype2
         INTEGER(kind=ikind)                           :: planrxy,plancz
         INTEGER(kind=ikind)                           :: planrx,plancyz
         INTEGER                                       :: nx,ny,nz

         !Cuda specific bits
         COMPLEX(KIND=GP), POINTER, DIMENSION (:,:,:)  :: ccarrt
         TYPE     (C_PTR)                              :: cu_ccd_,cu_ccd1_
         TYPE     (C_PTR)                              :: cu_cd_,cu_rd_
         TYPE     (C_PTR)                              :: pccarr_,pcarr_
         TYPE     (C_PTR)                              :: prarr_
         INTEGER  (C_INT),        DIMENSION (nstreams) :: icuplanr_
         INTEGER  (C_INT),        DIMENSION (nstreams) :: icuplanc_
         INTEGER(C_SIZE_T)                             :: szccd_,szcd_,szrd_
         INTEGER(C_SIZE_T),       DIMENSION (nstreams) :: str_szccd_
         INTEGER(C_SIZE_T),       DIMENSION (nstreams) :: str_szcd_
         INTEGER(C_SIZE_T),       DIMENSION (nstreams) :: str_szrd_
      END TYPE FFTPLAN
      SAVE

  END MODULE fftplans
!=================================================================

  MODULE cutypes
      USE ISO_C_BINDING
      IMPLICIT NONE
      
   
      INTEGER, PARAMETER :: CUFFT_SUCCESS       =X'0'
      INTEGER, PARAMETER :: CUFFT_INVALID_PLAN  =X'1'
      INTEGER, PARAMETER :: CUFFT_ALLOC_FAILED  =X'2'
      INTEGER, PARAMETER :: CUFFT_INVALID_TYPE  =X'3'
      INTEGER, PARAMETER :: CUFFT_INVALID_VALUE =X'4'
      INTEGER, PARAMETER :: CUFFT_INTERNAL_ERROR=X'5'
      INTEGER, PARAMETER :: CUFFT_EXEC_FAILED   =X'6'
      INTEGER, PARAMETER :: CUFFT_SETUP_FAILED  =X'7'
      INTEGER, PARAMETER :: CUFFT_INVALID_SIZE  =X'8'
      INTEGER, PARAMETER :: CUFFT_UNALIGNED_DATA=X'9'

      ! 
      INTEGER, PARAMETER :: CUFFT_R2C=X'2a'
      INTEGER, PARAMETER :: CUFFT_C2R=X'2c'
      INTEGER, PARAMETER :: CUFFT_C2C=X'29'
      INTEGER, PARAMETER :: CUFFT_D2Z=X'6a'
      INTEGER, PARAMETER :: CUFFT_Z2D=X'6c'
      INTEGER, PARAMETER :: CUFFT_Z2Z=X'69'
      INTEGER, PARAMETER :: cudaHostAllocDefault =0
      INTEGER, PARAMETER :: cudaHostAllocPortable=1
      INTEGER, PARAMETER :: cudaHostAllocMapped  =2
      INTEGER, PARAMETER :: cudaDeviceMapHost    =3
     
      ENUM, BIND(C)
        ENUMERATOR ::                            &
        cudaSuccess                       =0 ,   &
        cudaErrorMissingConfiguration     =1 ,   &
        cudaErrorMemoryAllocation         =2 ,   &
        cudaErrorInitializationError      =3 ,   &
        cudaErrorLaunchFailure            =4 ,   &
        cudaErrorPriorLaunchFailure       =5 ,   &
        cudaErrorLaunchTimeout            =6 ,   &
        cudaErrorLaunchOutOfResources     =7 ,   &
        cudaErrorInvalidDeviceFunction    =8 ,   &
        cudaErrorInvalidConfiguration     =9 ,   &
        cudaErrorInvalidDevice            =10,   &
        cudaErrorInvalidValue             =11,   &
        cudaErrorInvalidPitchValue        =12,   &
        cudaErrorInvalidSymbol            =13,   &
        cudaErrorMapBufferObjectFailed    =14,   &
        cudaErrorUnmapBufferObjectFailed  =15,   &
        cudaErrorInvalidHostPointer       =16,   &
        cudaErrorInvalidDevicePointer     =17,   &
        cudaErrorInvalidTexture           =18,   &
        cudaErrorInvalidTextureBinding    =19,   &
        cudaErrorInvalidChannelDescriptor =20,   &
        cudaErrorInvalidMemcpyDirection   =21,   &
        cudaErrorAddressOfConstant        =22,   &
        cudaErrorTextureFetchFailed       =23,   &
        cudaErrorTextureNotBound          =24,   &
        cudaErrorSynchronizationError     =25,   &
        cudaErrorInvalidFilterSetting     =26,   &
        cudaErrorInvalidNormSetting       =27,   &
        cudaErrorMixedDeviceExecution     =28,   &
        cudaErrorCudartUnloading          =29,   &
        cudaErrorUnknown                  =30,   &
        cudaErrorNotYetImplemented        =31,   &
        cudaErrorMemoryValueTooLarge      =32,   &
        cudaErrorInvalidResourceHandle    =33,   &
        cudaErrorNotReady                 =34,   &
        cudaErrorInsufficientDriver       =35,   &
        cudaErrorSetOnActiveProcess       =36,   &
        cudaErrorNoDevice                 =37,   &
        cudaErrorStartupFailure           =38,   &
        cudaErrorApiFailureBase           =39
      END ENUM

      TYPE, BIND(C) :: cudaDevicePropG
        INTEGER   (C_INT) :: canMapHostMemory
        INTEGER   (C_INT) :: clockRate
        INTEGER   (C_INT) :: computeMode
        INTEGER   (C_INT) :: deviceOverlap
        INTEGER   (C_INT) :: integrated
        INTEGER   (C_INT) :: kernelExecTimeoutEnabled
        INTEGER   (C_INT) :: major
        INTEGER   (C_INT) :: maxGridSize(3)
        INTEGER   (C_INT) :: maxThreadsDim(3)
        INTEGER   (C_INT) :: maxThreadsPerBlock
        INTEGER(C_SIZE_T) :: memPitch
        INTEGER   (C_INT) :: minor
        INTEGER   (C_INT) :: multiProcessorCount
        CHARACTER(C_CHAR) :: name(256) 
        INTEGER   (C_INT) :: regsPerBlock
        INTEGER(C_SIZE_T) :: sharedMemPerBlock
        INTEGER(C_SIZE_T) :: textureAlignment
        INTEGER(C_SIZE_T) :: totalConstMem
        INTEGER(C_SIZE_T) :: totalGlobalMem
        INTEGER   (C_INT) :: warpSize
      END TYPE cudaDevicePropG

  END MODULE cutypes
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
