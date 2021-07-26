!=================================================================
! MODULES for 3D codes
!
! 2009 Duane Rosenberg and Pablo D. Mininni.
!      National Center for Atmospheric Research.
!=================================================================

!=================================================================

  MODULE fprecision
!
! 
#if defined(GDOUBLE_PRECISION)
      INTEGER, PARAMETER :: GP = KIND(0.0D0)
#elif defined(GSINGLE_PRECISION)
      INTEGER, PARAMETER :: GP = KIND(0.0)
#else
#  error 'MODULE FPRECISION: PRECISION must be GDOUBLE_PRECISION or GSINGLE_PRECISION'
#endif
      SAVE

  END MODULE fprecision
!=================================================================

!=================================================================

  MODULE commtypes
      INCLUDE 'mpif.h'
!
! 
#if defined(GDOUBLE_PRECISION)
      INTEGER, SAVE :: GC_REAL    = MPI_DOUBLE_PRECISION
      INTEGER, SAVE :: GC_COMPLEX = MPI_DOUBLE_COMPLEX

      INTEGER, SAVE :: GC_2REAL    = MPI_2DOUBLE_PRECISION
      INTEGER, SAVE :: GC_2COMPLEX = MPI_2DOUBLE_COMPLEX
#elif defined(GSINGLE_PRECISION)
      INTEGER, SAVE :: GC_REAL    = MPI_REAL
      INTEGER, SAVE :: GC_COMPLEX = MPI_COMPLEX

      INTEGER, SAVE :: GC_2REAL    = MPI_2REAL
      INTEGER, SAVE :: GC_2COMPLEX = MPI_2COMPLEX
#else
#  error 'MODULE COMMTYPES: PRECISION must be GDOUBLE_PRECISION or GSINGLE_PRECISION'
#endif

  END MODULE commtypes
!=================================================================
