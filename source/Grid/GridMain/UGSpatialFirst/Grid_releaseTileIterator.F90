!!****if* source/Grid/GridMain/AMR/Grid_releaseTileIterator
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
!!  Grid_releaseTileIterator
!!
!! SYNOPSIS
!!  Grid_releaseTileIterator(Grid_iterator_t(INOUT) :: itor)
!!  
!! DESCRIPTION 
!!  Destroy given block/tile iterator.
!!
!! ARGUMENTS 
!!  itor - the block/tile iterator to destroy.
!!
!! SEE ALSO
!!  Grid_getTileIterator
!!
!!***

#ifdef DEBUG_ALL
#define DEBUG_GRID_ITER
#endif

subroutine Grid_releaseTileIterator(itor)
  use Grid_iterator, ONLY : Grid_iterator_t, destroy_iterator
#ifdef DEBUG_GRID_ITER
  use Grid_data,     ONLY : gr_meshMe
#endif

  implicit none

  type(Grid_iterator_t), intent(INOUT) :: itor

#ifdef DEBUG_GRID_ITER
  integer iterNo

  iterNo = itor%itemNumber()

999 format(A,I3,':',I5,' iterations done.')
  print 999,'Grid_releaseTileIterator @',gr_meshMe,iterNo-1
#endif

  call destroy_iterator(itor)
end subroutine Grid_releaseTileIterator

