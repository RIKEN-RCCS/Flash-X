!!****if* source/Grid/GridMain/UG/Grid_getGlobalIndexLimits
!! NOTICE
!!  Copyright 2022 UChicago Argonne, LLC and contributors
!!
!!  Licensed under the Apache License, Version 2.0 (the "License");
!!  you may not use this file except in compliance with the License.
!!
!!  Unless required by applicable law or agreed to in writing, software
!!  distributed under the License is distributed on an "AS IS" BASIS,
!!  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!!  See the License for the specific language governing permissions and
!!  limitations under the License.
!!
!! NAME
!!  Grid_getGlobalIndexLimits
!!
!! SYNOPSIS
!!
!!  call Grid_getGlobalIndexLimits(integer(OUT) :: globalIndexLimits(MDIM))
!!  
!!
!! DESCRIPTION 
!!  Gets the integer index dimensions of the entire grid 
!!  across all processors. Guardcells are not included.
!!
!!  globalIndexLimits(IAXIS) = highest index of grid in i dir  
!!  globalIndexLimits(JAXIS) = highest index of grid in j dir   
!!  globalIndexLimits(KAXIS) = highest index of grid in k dir   
!!
!!  (IAXIS, JAXIS and KAXIS are defined in constants.h
!!  and are set to 1,2 and 3 respectively)
!!
!! ARGUMENTS
!!  globalIndexLimits - returned array
!!
!! EXAMPLE
!!   For a 2d problem with a uniform grid block size of 8 
!!   With 4 blocks laid out in a square (2x2) grid
!! 
!!   globalIndexLimits(IAXIS) = 16
!!   globalIndexLimits(JAXIS) = 16
!!   globalIndexLimits(KAXIS) = 1  !because only 2d
!!
!!
!!***



subroutine Grid_getGlobalIndexLimits(globalIndexLimits)

  use Grid_data, ONLY : gr_gIndexSize

  implicit none

#include "constants.h"

  integer, dimension(MDIM), intent(OUT) :: globalIndexLimits

  globalIndexLimits = gr_gIndexSize

end subroutine Grid_getGlobalIndexLimits
