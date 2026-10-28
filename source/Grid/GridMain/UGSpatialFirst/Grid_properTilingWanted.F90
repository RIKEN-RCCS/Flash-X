logical function Grid_properTilingWanted()
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

#include "Simulation.h"
#include "constants.h"

  use Grid_data, ONLY: gr_enableTiling, gr_tileSize
  use Grid_data, ONLY:           gr_ilo, gr_ihi, &
                                 gr_jlo, gr_jhi, &
                                 gr_klo, gr_khi
  implicit none

  logical :: proper

  Grid_properTilingWanted = .FALSE.

  if (gr_enableTiling) then
     proper = (gr_ilo+gr_tileSize(IAXIS) .LE. gr_ihi)
     if (.NOT. proper .AND. NDIM > 1) proper = (gr_jlo+gr_tileSize(JAXIS) .LE. gr_jhi)
     if (.NOT. proper .AND. NDIM > 2) proper = (gr_klo+gr_tileSize(KAXIS) .LE. gr_khi)
     Grid_properTilingWanted = proper
  end if
end function Grid_properTilingWanted
