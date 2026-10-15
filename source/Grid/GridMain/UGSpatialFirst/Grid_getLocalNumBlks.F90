!!****if* source/Grid/GridMain/UG/Grid_getLocalNumBlks
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
!!  Grid_getLocalNumBlks
!!
!! SYNOPSIS
!!
!!  call Grid_getLocalNumBlks(integer(OUT) :: numBlocks)
!!  
!! DESCRIPTION 
!!  Get the number of local blocks on a processor 
!!
!! ARGUMENTS
!!  numBlocks : The number of blocks currently in use on myProcessor
!!
!! HISTORY
!!
!!  CREATED  05/18/04, by AD
!! 
!! 
!! 
!! 
!!
!!***

#ifdef DEBUG_ALL
#define DEBUG_GRID
#endif


subroutine Grid_getLocalNumBlks(numBlocks)
implicit none
  integer,intent(out) :: numBlocks
  numBlocks = 1
  return
end subroutine Grid_getLocalNumBlks
