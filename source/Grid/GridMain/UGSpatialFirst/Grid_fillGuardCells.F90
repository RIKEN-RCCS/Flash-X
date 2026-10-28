!!****if* source/Grid/GridMain/UG/Grid_fillGuardCells
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
!!  Grid_fillGuardCells
!!
!! SYNOPSIS
!!
!!  call Grid_fillGuardCells(integer(IN) :: gridDataStruct,
!!                           integer(IN) :: idir,
!!                  optional,integer(IN) :: minLayers,
!!                  optional,integer(IN) :: eosMode,
!!                  optional,logical(IN) :: doEos,
!!                  optional,integer(IN) :: maskSize,
!!                  optional,logical(IN) :: mask(maskSize),
!!                  optional,logical(IN) :: makeMaskConsistent,
!!                 optional, logical(IN) :: doLogMask,
!!                  optional,integer(IN) :: selectBlockType,
!!                  optional,logical(IN) :: unitReadsMeshDataOnly)
!!
!!
!! DESCRIPTION 
!!  
!!  The argument "gridDataStruct" can take on one of many valid 
!!  values to determine a specific grid data structure on which to apply
!!  the guardcell fill operation. The currently available options are listed with
!!  the arguments. Most users will use CENTER as the option, 
!!  since applications typically use the cell centered grid data, and they want
!!  guardcells to be filled for all the variables.
!!  More specialized applications, such as the unsplit methods, may want to use
!!  other options. 
!!  The user can also choose to fill guard cells either in a single direction,
!!  or all of them. 
!!
!!
!! ARGUMENTS 
!!  
!!
!!  gridDataStruct - integer constant, defined in "constants.h", 
!!                   indicating which grid data structure 
!!                   variable's guardcells to fill.
!!                   UG has 4 data structures for grid variables that
!!                   can have their guardcells filled. 
!!
!!                   unk                all cell centered variables in the grid
!!                   facex,facey,facez  all face centered variables along i,j,k 
!!                                      direction respectively
!!                   
!!                   valid values of gridDataStruct are  
!!                   CENTER             unk only
!!                   WORK               has no meaning in UG
!!                   The face variables are not yet implemented in UG
!!                   FACES              facex,facey, and facez
!!                   FACEX              facex
!!                   FACEY              facey
!!                   FACEZ              facez
!!                   CENTER_FACES     unk,facex,facey,facez
!!
!!  idir - direction of guardcell fill.  User can specify ALLDIR for all (x,y,z)
!!         directions, or if for example the algorithm only does one directional
!!         sweep at a time then time can be saved by filling only the guardcell
!!         direction that is needed.  A user would pass in the constants defined
!!         in constants.h IAXIS, JAXIS or KAXIS to fill guardcells in only one 
!!         direction.        
!!         All layers of guardcells in the given direction(s) are filled.
!!         In the current UG implementation, idir is ignored and a full
!!         guardcell fill in all directions is always performed.
!!
!!         THE REMAINING ARGUMENTS HAVE NO MEANING IN UG
!!
!!  minLayers - number of guardcell layers requested for all directions.
!!
!!   eosMode  - The mode in which eos is to be applied
!!   doEos    - the UG implementation does not act upon this argument
!!   maskSize - the size of the mask array. 
!! 
!!  mask -  It is a one-dimensional logical array 
!!          with indices corresponding to variables in the grid data
!!          structures. If a variable should have its guardcells filled,
!!          the corresponding element in "mask" is true, otherwise it is
!!          false.
!!          The mask is always ignored in this implementation.
!!  
!! makeMaskConsistent - IGNORED in this implementation
!!
!! selectBlockType - IGNORED
!!
!! unitReadsMeshDataOnly - specifies that the unit calling Grid_fillGuardCells
!!                         does not update any internal grid data.  This
!!                         allows us to skip the next guard cell fill because
!!                         the guard cells already contain up to date data.
!!                         IGNORED in this implementation.
!!
!! EXAMPLE
!!
!!   #include "Simulation.h"
!!   #include "constants.h"
!!
!!      call Grid_fillGuardCells( CENTER, IAXIS)
!!
!!     This call will fill all guardcells for all cell-centered 
!!     variables in the x direction.
!!     
!! EXAMPLE 2
!!
!!   #include "Simulation.h"
!!   #include "constants.h"
!!
!!      call Grid_fillGuardCells( CENTER_FACES, ALLDIR)
!!     
!!     This call fills guardcells along all directions in both
!!     cell centered and face centered data structures.
!!
!! NOTES
!!
!!   The masking functionality is not yet included in UG
!!  
!!***


#ifdef DEBUG_ALL
#define DEBUG_GRID
#endif


subroutine Grid_fillGuardCells( gridDataStruct,idir,minLayers,eosMode,doEos&
     &, maskSize, mask, makeMaskConsistent, doLogMask,selectBlockType, &
     unitReadsMeshDataOnly)

  use Grid_data, ONLY : gr_axisComm, gr_exch, gr_gridDataStruct, &
       gr_justExchangedGC,gr_domainBC, &
       gr_offset,gr_allPeriodic,gr_bndOrder, gr_meshMe
  use Grid_interface, ONLY : Grid_getTileIterator, Grid_releaseTileIterator
  use Driver_interface, ONLY : Driver_abort
  use gr_bcInterface, ONLY : gr_bcApplyToAllBlks
  use Grid_tile, ONLY : Grid_tile_t
  use Grid_iterator, ONLY: Grid_iterator_t
  implicit none
#include "constants.h"
#include "Simulation.h"

  type(Grid_tile_t) :: tileDesc
  type(Grid_iterator_t) :: itor
  integer, intent(in) :: gridDataStruct
  integer, intent(in) :: idir
  integer, optional,intent(in) :: minLayers
  integer, optional,intent(in) :: eosMode
  logical, optional,intent(in) :: doEos
  integer, optional, intent(IN) :: maskSize
  logical, optional, dimension(:),intent(IN) :: mask
  logical, optional,intent(IN) :: makeMaskConsistent
  logical, optional,intent(IN) :: doLogMask  ! ignored in this implementation
  integer,optional,intent(in) :: selectBlockType
  logical, optional, intent(in) :: unitReadsMeshDataOnly

  ! integer that hold the boundary condition type (ie, PERIODIC)
  integer,dimension(MDIM) :: guard !DEV: unused - KW 2010-09-22
  integer,dimension(LOW:HIGH,MDIM) :: blkLimits, blkLimitsGC

  integer :: i,j,k,beginDataType,endDataType
  integer :: axis,n,localGridDataStruct
  logical :: isWork=.false.
  
  !These are arrays to hold the sending start point indicies and receiving
  !start point indicies for the MPI data types in the shift data subroutine
  !They are used to specify the starting incidies in the unk array
  ! The four entries are i, j, k, and variable, matching spatial-first storage.
  integer, dimension(MDIM, MDIM+1) :: sendRight, sendLeft, recvRight, recvLeft

  logical,save :: firstTime=.true.

  !if we have just exchanged guardcells then we don't need to do in again.
  !FOR FUTURE: not yet implemented
!!$  if (gr_justExchanged) then
!!$     return
!!$  end if

  

#ifdef DEBUG_GRID

  if((gr_domainBC(LOW,IAXIS) == PERIODIC).and.(gr_domainBC(HIGH,IAXIS) /= PERIODIC)) then
     call Driver_abort("Gaurd Cell fill : one edge is periodic, one is not")
  end if
  
  if((gr_domainBC(LOW,JAXIS) == PERIODIC).and.(gr_domainBC(HIGH,JAXIS) /= PERIODIC)) then
     call Driver_abort("Gaurd Cell fill : one edge is periodic, one is not")
  end if
  
  if((gr_domainBC(LOW,KAXIS) == PERIODIC).and.(gr_domainBC(HIGH,KAXIS) /= PERIODIC)) then
     call Driver_abort("Gaurd Cell fill : one edge is periodic, one is not")
  end if
  
  if((gr_domainBC(LOW,IAXIS) /= PERIODIC).and.(gr_domainBC(HIGH,IAXIS) == PERIODIC)) then
     call Driver_abort("Gaurd Cell fill : one edge is periodic, one is not")
  end if
  
  if((gr_domainBC(LOW,JAXIS) /= PERIODIC).and.(gr_domainBC(HIGH,JAXIS) == PERIODIC)) then
     call Driver_abort("Gaurd Cell fill : one edge is periodic, one is not")
  end if
  
  if((gr_domainBC(LOW,KAXIS) /= PERIODIC).and.(gr_domainBC(HIGH,KAXIS) == PERIODIC)) then
     call Driver_abort("Gaurd Cell fill : one edge is periodic, one is not")
  end if
  
#endif
  
  localGridDataStruct=gridDataStruct
  if((gridDataStruct==CENTER_FACES).and.(NFACE_VARS==0))&
       localGridDataStruct=CENTER
  select case(localGridDataStruct)
  case (CENTER)
     beginDataType=CENTER_DATATYPE
     endDataTYPE=CENTER_DATATYPE
  case (CENTER_FACES)
     beginDataType=CENTER_DATATYPE
     if(NDIM==1)endDataType=FACEX_DATATYPE
     if(NDIM==2)endDataType=FACEY_DATATYPE
     if(NDIM==3)endDataType=FACEZ_DATATYPE
  case(FACES)
     beginDataType=FACEX_DATATYPE
     if(NDIM==1)endDataType=FACEX_DATATYPE
     if(NDIM==2)endDataType=FACEY_DATATYPE
     if(NDIM==3)endDataType=FACEZ_DATATYPE
  case(FACEX)
     beginDataType=FACEX_DATATYPE
     endDataType=FACEX_DATATYPE
  case(FACEY)
     if(NDIM>1) then
        beginDataType=FACEY_DATATYPE
        endDataType=FACEY_DATATYPE
     else
        if(gr_meshMe == MASTER_PE)print*,'warning: trying to fill face along Y for 1D problem'
        return
     end if
  case(FACEZ)

     if(NDIM>2) then
        beginDataType=FACEZ_DATATYPE
        endDataType=FACEZ_DATATYPE
     else
        if(gr_meshMe == MASTER_PE)print*,'warning: trying to fill face along Z for 2D problem'
        return
     end if
  case default
     call Driver_abort('Grid_fillGuardCells: unsupported grid data structure')
  end select

  do i = beginDataType,endDataType
     call Grid_getTileIterator(itor, ALL_BLKS, tiling=.FALSE.)
     do while(itor%isValid())
        call itor%currentTile(tileDesc)
        blkLimits(:,:)=tileDesc%limits
        blkLimitsGC(:,:)=tileDesc%grownLimits
        guard(:)=blkLimits(LOW,:)-blkLimitsGC(LOW,:) 
        recvLeft(:,MDIM+1)=1
        sendLeft(:,MDIM+1)=1
        recvRight(:,MDIM+1)=1
        sendRight(:,MDIM+1)=1
        do k=1,MDIM
           recvLeft(:,k)=blkLimitsGC(LOW,k)
           sendLeft(:,k)=blkLimitsGC(LOW,k)
           recvRight(:,k)=blkLimitsGC(LOW,k)
           sendRight(:,k)=blkLimitsGC(LOW,k)
        end do
        !! do a default initialization of all starting points
        !! and then adjust individual ones as needed
        !! index of interior cell which will be first GC in block to right
        do j = 1,NDIM
           sendRight(j,j) = blkLimits(HIGH,j) - guard(j) + 1
           !index of first interior cell to be sent to be GC on block to left
           sendLeft(j,j) = blkLimits(LOW,j)
           recvLeft(j,j) = blkLimits(HIGH,j)+1 !recv index of GC
           
           call gr_shiftData(gr_axisComm(j), gr_exch(i,j), &
                sendRight(j,:), sendLeft(j,:), &
                recvRight(j,:),recvLeft(j,:),gr_gridDataStruct(i))
        end do
        call itor%next()
     end do
     call Grid_releaseTileIterator(itor)
  end do
  if(.not.gr_allPeriodic) then
     do n = 0,NDIM-1
        axis = gr_bndOrder(NDIM-n)
        call gr_bcApplyToAllBlks(axis,isWork)
     end do
  end if
  gr_justExchangedGC = .true.
  return
end subroutine Grid_fillGuardCells

