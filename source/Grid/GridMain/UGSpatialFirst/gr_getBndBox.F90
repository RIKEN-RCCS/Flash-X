!!****if* source/Grid/GridMain/UG/gr_getBndBox
!!
!! NAME
!!  gr_getBndBox
!!
!! SYNOPSIS
!!
!! 
!!  gr_getBndBox(integer(IN)  :: real(OUT) :: boundBox(2, MDIM))
!!  
!! DESCRIPTION 
!!
!!  Gets the physical domain bounding box of the block identified 
!!  by blockId.  For each dimension the left (lower or forward) 
!!  physical coordinate of the block edge and the right (upper or back) 
!!  physical coordinate of the block edge is returned.  See arguments
!!  below for more detail.
!!
!! ARGUMENTS
!!
!!
!!  boundBox - returned array holding the boundBox coordinates in
!!             each dimension
!!
!!            for readability, in constants.h we define IAXIS = 1, JAXIS = 2, KAXIS = 3
!!
!!            boundBox(1,IAXIS) = left edge coordinate of block in x direction
!!            boundBox(2,IAXIS) = right edge coordinate of block in x direction
!!            boundBox(1,JAXIS) = top edge coordinate of block in y direction
!!            boundBox(2,JAXIS) = bottom edge coordinate of block in y direction
!!            boundBox(1,KAXIS) = front edge coordinate of block in z direction
!!            boundBox(2,KAXIS) = back edge coordinate of block in z direction
!!
!! EXAMPLE
!!  
!!   In 2 dimensions, if physical coordinates are ...
!!    
!!     ________________(0.5 1.0)
!!    |                |
!!    |                |
!!    |                |
!!    |                |
!!    |                |
!!    |                |
!!    |                |
!!    |_______________ |
!!  (-0.5, 0.0)
!!
!!
!!
!!     boundBox(1, IAXIS) = -0.5
!!     boundBox(2, IAXIS) = 0.5
!!     boundBox(1, JAXIS) = 0.0
!!     boundBox(2, JAXIS) = 1.0
!!     boundBox(1, KAXIS) = 1 !returned as 1 because only 2 dims
!!     boundBox(1, KAXIS) = 1 !returned as 1 because only 2 dims
!!
!!
!!***


#ifdef DEBUG_ALL
#define DEBUG_GRID
#endif

subroutine gr_getBndBox(boundBox)
  use Grid_data, ONLY :gr_iCoords,gr_jCoords,gr_kCoords
  use Grid_data, ONLY: gr_ilo,gr_ihi,gr_jlo,gr_jhi,gr_klo,gr_khi

  implicit none

#include "constants.h"
#include "Simulation.h"

  real,dimension(2,MDIM),intent(out) :: boundBox


  boundBox = 0.0



  boundBox(1,IAXIS) = gr_iCoords(LEFT_EDGE,gr_ilo,1)
#if NDIM > 1
  boundBox(1,JAXIS) = gr_jCoords(LEFT_EDGE,gr_jlo,1)
#endif
#if NDIM > 2
  boundBox(1,KAXIS) = gr_kCoords(LEFT_EDGE,gr_klo,1)
#endif



  boundBox(2,IAXIS) = gr_iCoords(RIGHT_EDGE,gr_ihi,1)
#if NDIM>1 
  boundBox(2,JAXIS) = gr_jCoords(RIGHT_EDGE,gr_jhi,1)
#endif
#if NDIM>2
  boundBox(2,KAXIS) = gr_kCoords(RIGHT_EDGE,gr_khi,1)
#endif

  return
end subroutine gr_getBndBox














