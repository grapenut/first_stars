
Module phch_interface

implicit none

#include "Flash.h"
#include "constants.h"

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Pixel functions
  interface
    subroutine pix2ang(ipix, theta, phi)
      implicit none
      integer,INTENT(IN) :: ipix
      real,INTENT(OUT) :: theta, phi
    end subroutine pix2ang
  end interface
  
  interface
    subroutine ang2pix(theta, phi, ipix)
      implicit none
      real,INTENT(IN) :: theta, phi
      integer,INTENT(OUT) :: ipix
    end subroutine ang2pix
  end interface

  interface
    subroutine vec2pix(vector, ipix)
      implicit none
      REAL,     INTENT(IN), dimension(1:) :: vector
      INTEGER, INTENT(OUT)               :: ipix
    end subroutine vec2pix
  end interface
  
  interface
    subroutine pix2vec  (ipix, vector, vertex)
      implicit none
      INTEGER, INTENT(IN)                             :: ipix
      REAL,     INTENT(OUT),dimension(1:)              :: vector
      REAL,     INTENT(OUT),dimension(1:,1:), optional :: vertex
    end subroutine pix2vec
  end interface

  interface
    subroutine ang2vec (theta, phi, vector)
      implicit none
      real, INTENT(IN) :: theta, phi
      real, INTENT(OUT), dimension(MDIM) :: vector
    end subroutine ang2vec
  end interface
  
  interface
    subroutine vec2ang(x, y, z, r, theta, phi)
      implicit none
      real, INTENT(IN) :: x, y, z
      real, INTENT(OUT) :: r, theta, phi
    end subroutine vec2ang
  end interface

  interface
    subroutine xyz2pix(x, y, z, pix)
      implicit none
      real, INTENT(IN) :: x, y, z
      integer, INTENT(OUT) :: pix
    end subroutine xyz2pix
  end interface

  interface
    subroutine xyz2bin(x, y, z, pix, bin)
      implicit none
      real, INTENT(IN) :: x, y, z
      integer, INTENT(OUT) :: pix, bin
    end subroutine xyz2bin
  end interface

  interface
    subroutine cube2pix(x, y, z, dx, dy, dz, level)
      implicit none
      real, INTENT(IN) :: x, y, z, dx, dy, dz
      integer, INTENT(IN) :: level
    end subroutine cube2pix
  end interface

  interface
    subroutine cell2pix(x, y, z, dx, dy, dz, dens, temp)
      implicit none
      real, INTENT(IN) :: x, y, z, dx, dy, dz, dens, temp
    end subroutine cell2pix
  end interface

  interface
    subroutine phch_mapPixels(numblocks, blocklist)
      implicit none
      integer, INTENT(IN) :: numblocks
      integer, INTENT(IN), dimension(numblocks)  :: blocklist
    end subroutine phch_mapPixels
  end interface

  interface
    subroutine phch_clearPixels(numblocks, blocklist)
      implicit none
      integer, INTENT(IN) :: numblocks
      integer, INTENT(IN), dimension(numblocks)  :: blocklist
    end subroutine phch_clearPixels
  end interface

  interface
    subroutine phch_calculateRate(temp, alpha)
      implicit none
      real, INTENT(IN) :: temp
      real, INTENT(OUT) :: alpha
    end subroutine phch_calculateRate
  end interface

  interface
    subroutine phch_updateSources()
      implicit none
    end subroutine phch_updateSources
  end interface
  
  interface
    subroutine phch_evolveChem(blockCount, blockList, dt)
      implicit none
      integer, INTENT(IN) :: blockCount
      integer, INTENT(IN), dimension(blockCount)  :: blocklist
      real, INTENT(IN) :: dt
    end subroutine phch_evolveChem
  end interface
  
  interface
    subroutine phch_evolvePhoto(y, photoeq, dens, eint, dt)
      implicit none
      real, INTENT(INOUT) :: y(SPECIES_BEGIN:SPECIES_END+1), eint
      real, INTENT(IN) :: photoeq, dens, dt
    end subroutine phch_evolvePhoto
  end interface
  
  interface
    subroutine phch_interpPhoto(y, ix, jx, xi0, metal0)
      implicit none
      real, INTENT(INOUT) :: y(SPECIES_BEGIN:SPECIES_END+1)
      integer, INTENT(IN) :: ix, jx
      real, INTENT(IN) :: xi0, metal0
    end subroutine phch_interpPhoto
  end interface
  
  interface
    subroutine phch_evolveMetal(y, dens, eint, dt)
      implicit none
      real, INTENT(INOUT) :: y(SPECIES_BEGIN:SPECIES_END+1), eint
      real, INTENT(IN) :: dens, dt
    end subroutine phch_evolveMetal
  end interface

  interface
    subroutine phch_interpMetal(y, ix, jx, kx, dens0, metal0, temp0)
      implicit none
      real, INTENT(INOUT) :: y(SPECIES_BEGIN:SPECIES_END+1)
      integer, INTENT(IN) :: ix, jx, kx
      real, INTENT(IN) :: dens0, metal0, temp0
    end subroutine phch_interpMetal
  end interface
  
  interface
    subroutine phch_gasConstant(y, dens, gascon)
      implicit none
      real, dimension(SPECIES_BEGIN:SPECIES_END+1), INTENT(IN) :: y
      real, INTENT(IN) :: dens
      real, INTENT(OUT) :: gascon
    end subroutine phch_gasConstant
  end interface

  interface
    subroutine phch_electronEquilibrium(y)
      implicit none
      real, dimension(SPECIES_BEGIN:SPECIES_END+1), INTENT(IN) :: y
    end subroutine phch_electronEquilibrium
  end interface

  interface
    subroutine phch_makeLogIndices(list, inum, imin, imax)
      implicit none
      integer, INTENT(IN) :: inum
      real, dimension(inum), INTENT(OUT) :: list
      real, INTENT(IN) :: imin, imax
    end subroutine phch_makeLogIndices
  end interface

  interface
    subroutine phch_findNumIndices(list, rows, inum, imin, imax)
      implicit none
      integer, INTENT(IN) :: rows
      real, dimension(rows), INTENT(IN) :: list
      integer, INTENT(OUT) :: inum
      real, INTENT(OUT) :: imin
      real, INTENT(OUT) :: imax
    end subroutine phch_findNumIndices
  end interface      

  interface
    subroutine phch_readFile(myPE, path, arr, rows, cols)
      integer,intent(IN) :: myPE
      character(len=*),intent(IN) :: path
      real,pointer,intent(INOUT) :: arr(:,:)
      integer,intent(OUT) :: rows
      integer,intent(OUT) :: cols
    end subroutine phch_readFile
  end interface

end Module phch_interface

