# 1 "fprecision_mod.fpp"
# 1 "<built-in>"
# 1 "<command-line>"

# 1 "/usr/include/stdc-predef.h" 1 3 4

# 17 "/usr/include/stdc-predef.h" 3 4











































# 1 "<command-line>" 2
# 1 "fprecision_mod.fpp"
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

      INTEGER, PARAMETER :: GP = KIND(0.0D0)





      SAVE

  END MODULE fprecision
!=================================================================

!=================================================================

  MODULE commtypes
      INCLUDE 'mpif.h'
!
! 

      INTEGER, SAVE :: GC_REAL    = MPI_DOUBLE_PRECISION
      INTEGER, SAVE :: GC_COMPLEX = MPI_DOUBLE_COMPLEX

      INTEGER, SAVE :: GC_2REAL    = MPI_2DOUBLE_PRECISION
      INTEGER, SAVE :: GC_2COMPLEX = MPI_2DOUBLE_COMPLEX
# 46 "fprecision_mod.fpp"

  END MODULE commtypes
!=================================================================

