!!****if* source/physics/Hydro/HydroMain/Spark/reconstruct
!!
!!  NAME
!!
!!  reconstruct
!!
!!  SYNOPSIS
!!
!!  call reconstruct (  )
!!
!!  DESCRIPTION
!!
!!  ARGUMENTS
!!
!!
!!***
real function minmod(a,b)
  implicit none
  real :: a,b
  minmod=.5 * (sign(1.,a) + sign(1.,b))*min(abs(a),abs(b))
end function minmod


subroutine reconstruct(i1,i2,i3,uPlusArray,uMinusArray,snake,flat)
    implicit none

#include "Simulation.h"
#include "Spark.h"
#define NRECON HY_NUM_VARS+NSPECIES+NMASS_SCALARS
  
    integer, intent(IN) :: i1,i2,i3
    real, dimension(:,:,:,:), intent(INOUT) :: uPlusArray, uMinusArray, snake
    real, dimension(:,:,:), intent(INOUT) :: flat
    real, dimension(NRECON) :: delbar, del_p, del_m, invSumAlpha
    real, dimension(NRECON,3) :: W5p, W5m, betaWeno, Alpha5, omega
    integer :: g, i
    real, parameter :: epsilon = 1.e-36
    real, parameter :: n13o12 = 13./12.
    real, dimension(NRECON) :: mags, betaSum
    !! Set WENO5 coefficients once and for all
    !u_{1,i+1/2}= 2/6*u_{i-2} -7/6*u_{i-1} +11/6*u_{i}
    real, dimension(3), parameter :: coeff1p1(1:3) = (/ 2./6., -7./6., 11./6./)
    !u_{2,i+1/2}=-1/6*u_{i-2} +5/6*u_{i-1} + 2/6*u_{i}
    real, dimension(3), parameter :: coeff1p2(1:3) = (/-1./6.,  5./6.,  2./6./)
    !u_{3,i+1/2}= 2/6*u_{i-2} +5/6*u_{i-1} - 1/6*u_{i}
    real, dimension(3), parameter :: coeff1p3(1:3) = (/ 2./6.,  5./6., -1./6./)
    !(gamma1,gamma2,gamma3)
    real, dimension(3), parameter :: coeff2p(1:3)   = (/0.1, 0.6, 0.3/)
    !u_{1,i-1/2}=-1/6*u_{i-2} +5/6*u_{i-1} + 2/6*u_{i}
    real, dimension(3), parameter :: coeff1m1(1:3) = (/-1./6.,  5./6.,  2./6./)
    !u_{2,i-1/2}= 2/6*u_{i-2} +5/6*u_{i-1} - 1/6*u_{i}
    real, dimension(3), parameter :: coeff1m2(1:3) = (/ 2./6.,  5./6., -1./6./)
    !u_{3,i-1/2}=11/6*u_{i-2} -7/6*u_{i-1} + 2/6*u_{i}
    real, dimension(3), parameter :: coeff1m3(1:3) = (/ 11./6.,-7./6.,  2./6./)
    !(gamma1,gamma2,gamma3)
    real, dimension(3), parameter :: coeff2m(1:3)   = (/0.3, 0.6, 0.1/)
    !$omp declare target
    ! Interpolation stencil for weno
    
    !! Calculate interface values at i+1/2
    W5p(:,1) = coeff1p1(1)*snake(:,i1-2,i2+0,i3+0) &
         + coeff1p1(2)*snake(:,i1-1,i2+0,i3+0) + coeff1p1(3)*snake(:,i1,i2,i3)
    W5p(:,2) = coeff1p2(1)*snake(:,i1-1,i2+0,i3+0) &
         + coeff1p2(2)*snake(:,i1,i2,i3) + coeff1p2(3)*snake(:,i1+1,i2+0,i3+0)
    W5p(:,3) = coeff1p3(1)*snake(:,i1,i2,i3) &
         + coeff1p3(2)*snake(:,i1+1,i2+0,i3+0) + coeff1p3(3)*snake(:,i1+2,i2+0,i3+0)
  
    !! Calculate interface values at i-1/2
    W5m(:,1) = coeff1m1(1)*snake(:,i1-2,i2+0,i3+0) &
         + coeff1m1(2)*snake(:,i1-1,i2+0,i3+0) + coeff1m1(3)*snake(:,i1,i2,i3)
    W5m(:,2) = coeff1m2(1)*snake(:,i1-1,i2+0,i3+0) &
         + coeff1m2(2)*snake(:,i1,i2,i3) + coeff1m2(3)*snake(:,i1+1,i2+0,i3+0)
    W5m(:,3) = coeff1m3(1)*snake(:,i1,i2,i3) &
         + coeff1m3(2)*snake(:,i1+1,i2+0,i3+0) + coeff1m3(3)*snake(:,i1+2,i2+0,i3+0)
  
    !! Calculate smoothness indicators at i+1/2
    betaWeno(:,1) = n13o12*(snake(:,i1-2,i2+0,i3+0) &
         - 2.*snake(:,i1-1,i2+0,i3+0) + snake(:,i1,i2,i3))**2 &
         + 0.25*(snake(:,i1-2,i2+0,i3+0) &
         - 4.*snake(:,i1-1,i2+0,i3+0) + 3.*snake(:,i1,i2,i3))**2
    betaWeno(:,2) = n13o12*(snake(:,i1-1,i2+0,i3+0) &
         - 2.*snake(:,i1,i2,i3) + snake(:,i1+1,i2+0,i3+0))**2 &
         + 0.25*(snake(:,i1-1,i2+0,i3+0) - snake(:,i1+1,i2+0,i3+0))**2
    betaWeno(:,3) = n13o12*(snake(:,i1,i2,i3) &
         - 2.*snake(:,i1+1,i2+0,i3+0) + snake(:,i1+2,i2+0,i3+0))**2 &
         + 0.25*(3.*snake(:,i1,i2,i3) &
         - 4.*snake(:,i1+1,i2+0,i3+0) + snake(:,i1+2,i2+0,i3+0))**2
  
    !! Use problem-adaptive epsilong as in Tchekovskoy+07, A3
    ! This does not seem to work with the WENO-Z indicators of Borges+08
    ! mags(:) = snake(:,ind-2)**2 + snake(:,ind-1)**2 + snake(:,ind)**2 &
    !      + snake(:,ind+1)**2 + snake(:,ind+2)**2
    ! betaWeno(:,1) = betaWeno(:,1) + epsilon*mags(:) + TINY(1.0)
    ! betaWeno(:,2) = betaWeno(:,2) + epsilon*mags(:) + TINY(1.0)
    ! betaWeno(:,3) = betaWeno(:,3) + epsilon*mags(:) + TINY(1.0)
  
    !! This is WENO-Z: this is very similar to weno5 with wenoExp=1
    Alpha5(:,1) = coeff2p(1)*(1.+(abs(betaWeno(:,1)-betaWeno(:,3))/(betaWeno(:,1)+epsilon)))
    Alpha5(:,2) = coeff2p(2)*(1.+(abs(betaWeno(:,1)-betaWeno(:,3))/(betaWeno(:,2)+epsilon)))
    Alpha5(:,3) = coeff2p(3)*(1.+(abs(betaWeno(:,1)-betaWeno(:,3))/(betaWeno(:,3)+epsilon)))
  
    !! Normalize weights at i+1/2
    invSumAlpha = 1./(Alpha5(:,1)+Alpha5(:,2)+Alpha5(:,3))
    omega(:,1)  = Alpha5(:,1)*invSumAlpha
    omega(:,2)  = Alpha5(:,2)*invSumAlpha
    omega(:,3)  = Alpha5(:,3)*invSumAlpha
  
    ! Compute interface value at i+1/2
    uPlusArray(:,i1,i2,i3) = omega(:,1)*W5p(:,1) + omega(:,2)*W5p(:,2) + omega(:,3)*W5p(:,3)
    !! Apply flattening
    uPlusArray(:,i1,i2,i3) = flat(i1,i2,i3)*uPlusArray(:,i1,i2,i3) &
         + (1.-flat(i1,i2,i3))*snake(:,i1,i2,i3)
  
    !! Now move on to i-1/2
    !! This is WENO-Z
    Alpha5(:,1) = coeff2m(1)*(1.+(abs(betaWeno(:,1)-betaWeno(:,3))/(betaWeno(:,1)+epsilon)))
    Alpha5(:,2) = coeff2m(2)*(1.+(abs(betaWeno(:,1)-betaWeno(:,3))/(betaWeno(:,2)+epsilon)))
    Alpha5(:,3) = coeff2m(3)*(1.+(abs(betaWeno(:,1)-betaWeno(:,3))/(betaWeno(:,3)+epsilon)))
  
    !! Normalize weights at i-1/2
    invSumAlpha = 1./(Alpha5(:,1)+Alpha5(:,2)+Alpha5(:,3))
    omega(:,1)  = Alpha5(:,1)*invSumAlpha
    omega(:,2)  = Alpha5(:,2)*invSumAlpha
    omega(:,3)  = Alpha5(:,3)*invSumAlpha
  
    !! Compute interface value at i-1/2
    uMinusArray(:,i1,i2,i3) = omega(:,1)*W5m(:,1) + omega(:,2)*W5m(:,2) + omega(:,3)*W5m(:,3)
    !! Apply flattening
    uMinusArray(:,i1,i2,i3) = flat(i1,i2,i3)*uMinusArray(:,i1,i2,i3) &
         + (1.-flat(i1,i2,i3))*snake(:,i1,i2,i3)
  
    !! Check for monotonicity
    where ( (uPlusArray(:,i1,i2,i3)-snake(:,i1,i2,i3)) &
         *(snake(:,i1,i2,i3)-uMinusArray(:,i1,i2,i3)) <= 0. )
      uPlusArray(:,i1,i2,i3)  = snake(:,i1,i2,i3)
      uMinusArray(:,i1,i2,i3) = snake(:,i1,i2,i3)
    end where
  
    ! print *, "gpu: (",i1,i2,i3,")",NEW_LINE('a'), "P ", uPlusArray(:,i1,i2,i3), NEW_LINE('a'),"M ", uMinusArray(:,i1,i2,i3)
    ! if (i1==17 .AND. i2==17) then
    !   do i=i1-2,i1+2
    !     print *, "Snake(",i,i2,i3,")", snake(7,i,i2,i3), flat(i,i2,i3)
    !   enddo
    ! endif

  end subroutine reconstruct




