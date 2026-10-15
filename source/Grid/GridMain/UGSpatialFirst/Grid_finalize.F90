!!****if* source/Grid/GridMain/UG/Grid_finalize
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
!!
!!  Grid_finalize
!!
!!
!! SYNOPSIS
!!
!!  Grid_finalize()
!!
!!
!! DESCRIPTION
!!
!!  Deallocates memory that has been allocated in the Grid Unit
!!
!!***


subroutine Grid_finalize()

  use gr_bcInterface, ONLY : gr_bcFinalize
  use gr_ptInterface, ONLY : gr_ptFinalize
  use Grid_data, ONLY : gr_iCoords,gr_jCoords,gr_kCoords, gr_gid
  use physicalData, ONLY : unk,facevarx,facevary,facevarz
  use Grid_data, ONLY : scratch,scratch_ctr,&
       &scratch_facevarx,scratch_facevary,scratch_facevarz,gr_flxx, gr_flxy, gr_flxz
  use gr_sbInterface, ONLY : gr_sbFinalize

  implicit none

#include "Simulation.h"

  deallocate(gr_gid)

#ifndef FIXEDBLOCKSIZE  
  if (allocated(gr_iCoords)) deallocate(gr_iCoords)
  if (allocated(gr_jCoords)) deallocate(gr_jCoords)
  if (allocated(gr_kCoords)) deallocate(gr_kCoords)
  if (allocated(unk)) deallocate(unk)
  if (allocated(scratch)) deallocate(scratch)
  if (allocated(facevarx)) deallocate(facevarx)
  if (allocated(facevary)) deallocate(facevary)
  if (allocated(facevarz)) deallocate(facevarz)
  if (allocated(scratch_ctr)) deallocate(scratch_ctr)
  if (allocated(scratch_facevarx)) deallocate(scratch_facevarx)
  if (allocated(scratch_facevary)) deallocate(scratch_facevary)
  if (allocated(scratch_facevarz)) deallocate(scratch_facevarz)
  if (allocated(gr_flxx)) deallocate(gr_flxx)
  if (allocated(gr_flxy)) deallocate(gr_flxy)
  if (allocated(gr_flxz)) deallocate(gr_flxz)
#endif
  call gr_solversFinalize()
  call gr_ptFinalize()
  call gr_bcFinalize()
  call gr_sbFinalize()
end subroutine Grid_finalize
