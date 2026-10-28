!!****if* source/Grid/GridMain/UG/gr_specificData
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
!!  Grid_data
!!
!! SYNOPSIS
!!
!!  use gr_specificData
!!
!! DESCRIPTION 
!!  
!!  
!!  Defining data structures for Uniform Grid
!!
!!***

!!REORDER(5):scratch, scratch_ctr, scratch_facevar[xyz]


Module gr_specificData
  implicit none

#include "Simulation.h"
#include "constants.h"


  integer,save :: gr_globalOffset !stores beginning blk offset for a proc
  
end Module gr_specificData
