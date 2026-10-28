!> @copyright Copyright 2023 UChicago Argonne, LLC and contributors
!!
!! @licenseblock
!! Licensed under the Apache License, Version 2.0 (the "License");
!! you may not use this file except in compliance with the License.
!!
!! Unless required by applicable law or agreed to in writing, software
!! distributed under the License is distributed on an "AS IS" BASIS,
!! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!! See the License for the specific language governing permissions and
!! limitations under the License.
!! @endlicenseblock
!!
!! 
!! 
!! This is nonfixed blocksize mode for Uniform Grid implementation of the Grid unit
!!

subroutine Grid_getMaxcells(maxcells)
  use Grid_data, ONLY : gr_maxCells
  implicit none
  integer, intent(OUT) :: maxCells
  maxcells=gr_maxCells
end subroutine Grid_getMaxcells
