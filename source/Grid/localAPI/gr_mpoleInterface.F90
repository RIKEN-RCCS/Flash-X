!!****ih* source/Grid/localAPI/gr_mpoleInterface
!! NOTICE
!!  Copyright 2023 UChicago Argonne, LLC and contributors
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
!!  gr_mpoleInterface
!!
!! SYNOPSIS
!!
!!  use gr_mpoleInterface
!!
!! DESCRIPTION
!!
!!  This is the header file for the multipole solver that defines its
!!  interfaces.
!!***

module gr_mpoleInterface

  implicit none

  interface
     subroutine gr_mpoleCen1Dspherical (idensvar)
       integer, intent (in) :: idensvar
     end subroutine gr_mpoleCen1Dspherical
  end interface

  interface
     subroutine gr_mpoleCen2Dcylindrical (idensvar)
       integer, intent (in) :: idensvar
     end subroutine gr_mpoleCen2Dcylindrical
  end interface

  interface
     subroutine gr_mpoleCen2Dspherical (idensvar)
       integer, intent (in) :: idensvar
     end subroutine gr_mpoleCen2Dspherical
  end interface

  interface
     subroutine gr_mpoleCen3Dcartesian (idensvar)
       integer, intent (in) :: idensvar
     end subroutine gr_mpoleCen3Dcartesian
  end interface

  interface
     subroutine gr_mpoleCen3Dcylindrical (idensvar)
       integer, intent (in) :: idensvar
     end subroutine gr_mpoleCen3Dcylindrical
  end interface

  interface
     subroutine gr_mpoleCen3Dspherical (idensvar)
       integer, intent (in) :: idensvar
     end subroutine gr_mpoleCen3Dspherical
  end interface

  interface
     subroutine gr_mpoleCenterOfExpansion (idensvar)
       integer, intent (in) :: idensvar
     end subroutine gr_mpoleCenterOfExpansion
  end interface

  interface
     subroutine gr_mpoleCollectMoments ()
     end subroutine gr_mpoleCollectMoments
  end interface

  interface
     subroutine gr_mpoleDeallocateRadialArrays ()
     end subroutine gr_mpoleDeallocateRadialArrays
  end interface

  interface
     subroutine gr_mpoleDumpMoments ()
     end subroutine gr_mpoleDumpMoments
  end interface

  interface
     subroutine gr_mpoleFinalize()
     end subroutine gr_mpoleFinalize
  end interface

  interface
     subroutine gr_mpoleHeapsort (nElements,Vector)
       integer, intent (in)    :: nElements
       real,    intent (inout) :: Vector (1:nElements)
     end subroutine gr_mpoleHeapsort
  end interface

  interface
     subroutine gr_mpoleInit()
     end subroutine gr_mpoleInit
  end interface

  interface
     subroutine gr_mpoleMom1Dspherical (idensvar)
       integer, intent (in) :: idensvar
     end subroutine gr_mpoleMom1Dspherical
  end interface

  interface
     subroutine gr_mpoleMom2Dcylindrical (idensvar)
       integer, intent(in) :: idensvar
     end subroutine gr_mpoleMom2Dcylindrical
  end interface

  interface
     subroutine gr_mpoleMom2Dspherical (idensvar)
       integer, intent (in) :: idensvar
     end subroutine gr_mpoleMom2Dspherical
  end interface

  interface
     subroutine gr_mpoleMom3Dcartesian (idensvar)
       integer, intent (in) :: idensvar
     end subroutine gr_mpoleMom3Dcartesian
  end interface

  interface
     subroutine gr_mpoleMom3Dcylindrical (idensvar)
       integer, intent(in) :: idensvar
     end subroutine gr_mpoleMom3Dcylindrical
  end interface

  interface
     subroutine gr_mpoleMom3Dspherical (idensvar)
       integer, intent(in) :: idensvar
     end subroutine gr_mpoleMom3Dspherical
  end interface

  interface
     subroutine gr_mpoleMomBins1Dspherical (maxQtype)
       integer, intent (in) :: maxQtype
     end subroutine gr_mpoleMomBins1Dspherical
  end interface

  interface
     subroutine gr_mpoleMomBins2Dcylindrical (maxQtype)
       integer, intent (in) :: maxQtype
     end subroutine gr_mpoleMomBins2Dcylindrical
  end interface

  interface
     subroutine gr_mpoleMomBins2Dspherical (maxQtype)
       integer, intent (in) :: maxQtype
     end subroutine gr_mpoleMomBins2Dspherical
  end interface

  interface
     subroutine gr_mpoleMomBins3Dcartesian (maxQtype)
       integer, intent (in) :: maxQtype
     end subroutine gr_mpoleMomBins3Dcartesian
  end interface

  interface
     subroutine gr_mpoleMomBins3Dcylindrical (maxQtype)
       integer, intent (in) :: maxQtype
     end subroutine gr_mpoleMomBins3Dcylindrical
  end interface

  interface
     subroutine gr_mpoleMomBins3Dspherical (maxQtype)
       integer, intent (in) :: maxQtype
     end subroutine gr_mpoleMomBins3Dspherical
  end interface

  interface
     subroutine gr_mpoleMoments (idensvar)
       integer, intent (in) :: idensvar
     end subroutine gr_mpoleMoments
  end interface

  interface
     subroutine gr_mpolePot1Dspherical (ipotvar)
       integer, intent (in) :: ipotvar
     end subroutine gr_mpolePot1Dspherical
  end interface

  interface
     subroutine gr_mpolePot2Dcylindrical (ipotvar)
       integer, intent (in) :: ipotvar
     end subroutine gr_mpolePot2Dcylindrical
  end interface

  interface
     subroutine gr_mpolePot2Dspherical (ipotvar)
       integer, intent (in) :: ipotvar
     end subroutine gr_mpolePot2Dspherical
  end interface

  interface
     subroutine gr_mpolePot3Dcartesian (ipotvar)
       integer, intent (in) :: ipotvar
     end subroutine gr_mpolePot3Dcartesian
  end interface

  interface
     subroutine gr_mpolePot3Dcylindrical (ipotvar)
       integer, intent (in) :: ipotvar
     end subroutine gr_mpolePot3Dcylindrical
  end interface

  interface
     subroutine gr_mpolePot3Dspherical (ipotvar)
       integer, intent (in) :: ipotvar
     end subroutine gr_mpolePot3Dspherical
  end interface

  interface
     subroutine gr_mpolePotentials (ipotvar, Poisson_factor,  &
                                    gr_mpoleMomentR,          &
                                    gr_mpoleMomentI,          &
                                    gr_mpoleScratch,          &
                                    gr_mpoleMaxQ,             &
                                    gr_mpoleGravityConstant,  &
                                    gr_mpoleFourPiInv,        &
                                    gr_mpoleGeometry,         &
                                    gr_mpoleRequest)
       integer, intent (in) :: ipotvar
       real,    intent (in) :: Poisson_factor
       real,           intent (inout) :: gr_mpoleMomentR (1:,0:)
       real,           intent (inout) :: gr_mpoleMomentI (1:,1:)
       real,           intent (inout) :: gr_mpoleScratch (1:,1:,1:)
       integer, intent (in) :: gr_mpoleMaxQ
       real,    intent (out) :: gr_mpoleGravityConstant
       real,    intent (in)  :: gr_mpoleFourPiInv
       integer, intent (in)  :: gr_mpoleGeometry
       integer, intent (inout) :: gr_mpoleRequest
     end subroutine gr_mpolePotentials
  end interface

  interface
     subroutine gr_mpolePrintRadialInfo ()
     end subroutine gr_mpolePrintRadialInfo
  end interface

  interface
     subroutine gr_mpoleRad1Dspherical()
     end subroutine gr_mpoleRad1Dspherical
  end interface

  interface
     subroutine gr_mpoleRad2Dcylindrical()
     end subroutine gr_mpoleRad2Dcylindrical
  end interface

  interface
     subroutine gr_mpoleRad2Dspherical()
     end subroutine gr_mpoleRad2Dspherical
  end interface

  interface
     subroutine gr_mpoleRad3Dcartesian()
     end subroutine gr_mpoleRad3Dcartesian
  end interface

  interface
     subroutine gr_mpoleRad3Dcylindrical()
     end subroutine gr_mpoleRad3Dcylindrical
  end interface

  interface
     subroutine gr_mpoleRad3Dspherical()
     end subroutine gr_mpoleRad3Dspherical
  end interface

  interface
     subroutine gr_mpoleRadialSampling()
     end subroutine gr_mpoleRadialSampling
  end interface

  interface
     subroutine gr_mpoleSetInnerZoneGrid (nRlocal,     &
                                          nRinnerZone, &
                                          nPinnerZone, &
                                          RinnerZone   )
       integer, intent (in)    :: nRlocal
       integer, intent (in)    :: nRinnerZone
       integer, intent (in)    :: nPinnerZone
       real,    intent (inout) :: RinnerZone (1:nRinnerZone)
     end subroutine gr_mpoleSetInnerZoneGrid
  end interface

  interface
     subroutine gr_mpoleSetOuterZoneGrid ()
     end subroutine gr_mpoleSetOuterZoneGrid
  end interface

  interface
     subroutine gr_mpoleSetRadialBinData ()
     end subroutine gr_mpoleSetRadialBinData
  end interface

end module gr_mpoleInterface
