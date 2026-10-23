#include "Simulation.h"
! One launch covers all blocks; arithmetic helpers are imported unchanged from NewImpl.
module hy_gpu_amr_kernels
  use hy_gpu_amr_types
  use omp_lib
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  integer, parameter :: rk=kind(1.0)
  integer, parameter :: primitiveMap(7)=[1,2,3,4,5,6,7],acousticMap(6)=[1,2,3,4,5,6]
  integer, parameter :: updateMap(6)=[1,2,3,4,8,7],conservedMap(5)=[1,2,3,4,5],xOffset(3)=[1,0,0]
contains
  subroutine flatten(state,flat,nx,ny,nz,nb,ng,lower,upper,params,error)
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    integer,intent(in):: nx,ny,nz,nb,ng,lower(3),upper(3)
    real(rk),intent(in)::state(1-ng:nx+ng,1-ng:ny+ng,1-(NDIM-2)*ng:nz+(NDIM-2)*ng,9,nb)
    real(rk),intent(inout)::flat(1-ng:nx+ng,1-ng:ny+ng,1-(NDIM-2)*ng:nz+(NDIM-2)*ng,nb)
    type(hydro_parameters),intent(in)::params
    integer,intent(out)::error
    integer::i,j,k,block,status
    real(rk)::f
integer :: macro_scope_1_local_d,macro_scope_1_local_e(3)
real(rk) :: macro_scope_1_local_p,macro_scope_1_local_g,macro_scope_1_local_dp,macro_scope_1_local_b, &
  & macro_scope_1_local_z,macro_scope_1_local_a,macro_scope_1_local_fd,macro_scope_1_local_un,macro_scope_1_local_up

    error=0
!$omp target teams distribute parallel do collapse(4) map(alloc:state,flat) &
!$omp& private(i,j,k) &
!$omp& private(block,status,f) &
!$omp& private(macro_scope_1_local_d,macro_scope_1_local_e,macro_scope_1_local_p) &
!$omp& private(macro_scope_1_local_g,macro_scope_1_local_dp,macro_scope_1_local_b) &
!$omp& private(macro_scope_1_local_z,macro_scope_1_local_a,macro_scope_1_local_fd) &
!$omp& private(macro_scope_1_local_un,macro_scope_1_local_up) &
!$omp& reduction(max:error)
    do block=1,nb
      do k=lower(3)-(NDIM-2),upper(3)+(NDIM-2)
      do j=lower(2)-1,upper(2)+1
      do i=lower(1)-1,upper(1)+1
macro_scope_1: do
status=HY_OK; f=1._rk
if (.not.params%flattening) exit macro_scope_1
macro_scope_1_local_p=state(i,j,k,primitiveMap(QP),block)
macro_scope_1_local_g=state(i,j,k,primitiveMap(QGAMMA),block)
if (.not.ieee_is_finite(macro_scope_1_local_p*macro_scope_1_local_g)) then
  status=HY_BAD_STATE; exit macro_scope_1
end if
if (macro_scope_1_local_p*macro_scope_1_local_g<=0._rk) then
  status=HY_BAD_STATE; exit macro_scope_1
end if
do macro_scope_1_local_d=1,NDIM
  macro_scope_1_local_e=0; macro_scope_1_local_e(macro_scope_1_local_d)=1
  macro_scope_1_local_dp=abs(state(i+macro_scope_1_local_e(1),j+macro_scope_1_local_e(2),k+macro_scope_1_local_e(3), &
  & primitiveMap(QP),block) - &
  state(i-macro_scope_1_local_e(1),j-macro_scope_1_local_e(2),k-macro_scope_1_local_e(3),primitiveMap(QP),block))
  macro_scope_1_local_b=macro_scope_1_local_dp/max(params%flattening_eta, &
  abs(state(i+2*macro_scope_1_local_e(1),j+2*macro_scope_1_local_e(2),k+2*macro_scope_1_local_e(3),primitiveMap(QP), &
  & block) - &
  state(i-2*macro_scope_1_local_e(1),j-2*macro_scope_1_local_e(2),k-2*macro_scope_1_local_e(3),primitiveMap(QP),block)))
  macro_scope_1_local_z=macro_scope_1_local_dp/(macro_scope_1_local_g*macro_scope_1_local_p)
  macro_scope_1_local_a=max(0._rk,min(1._rk,(0.85_rk-macro_scope_1_local_b)/0.10_rk))
  macro_scope_1_local_un=state(i,j,k,primitiveMap(QUX+macro_scope_1_local_d-1),block)
  macro_scope_1_local_up=state(i+macro_scope_1_local_e(1),j+macro_scope_1_local_e(2),k+macro_scope_1_local_e(3), &
  & primitiveMap(QUX+macro_scope_1_local_d-1),block)
  macro_scope_1_local_fd=1._rk
  if (macro_scope_1_local_up<macro_scope_1_local_un) macro_scope_1_local_fd=max(macro_scope_1_local_a,min(1._rk,( &
  & 0.75_rk-macro_scope_1_local_z)/0.50_rk))
  f=min(f,macro_scope_1_local_fd)
end do
if (.not.ieee_is_finite(f)) status=HY_BAD_STATE
exit macro_scope_1
end do macro_scope_1
        error=max(error,status)
        flat(i,j,k,block)=f
      end do
    end do
      end do
      end do
    !$omp end target teams distribute parallel do
  end subroutine

  subroutine reconstruct(state,flat,minus,plus,nx,ny,nz,nb,ng,lower,upper,params,d,error)

    integer,intent(in)::nx,ny,nz,nb,ng,lower(3),upper(3),d
    real(rk),intent(in)::state(1-ng:nx+ng,1-ng:ny+ng,1-(NDIM-2)*ng:nz+(NDIM-2)*ng,9,nb),flat(1-ng:nx+ng,1-ng:ny+ng,1-( &
  & NDIM-2)*ng:nz+(NDIM-2)*ng,nb)
    real(rk),intent(inout)::minus(1-ng:nx+ng,1-ng:ny+ng,1-(NDIM-2)*ng:nz+(NDIM-2)*ng,7,nb),plus(1-ng:nx+ng,1-ng:ny+ng, &
  & 1-(NDIM-2)*ng:nz+(NDIM-2)*ng,7,nb)
    integer::e(3)
    type(hydro_parameters),intent(in)::params
    integer,intent(out)::error
    integer::i,j,k,block,v
    real(rk)::stencil(-2:2),cm(3),cp(3),wm(3),wp(3),qm,qp,f
integer :: macro_scope_10_local_s,macro_scope_10_local_ii,macro_scope_10_local_jj,macro_scope_10_local_kk
real(rk) :: macro_scope_13_local_beta(3),macro_scope_13_local_correction(3),macro_scope_13_local_tau

    e=0; e(d)=1
    error=0
!$omp target teams distribute parallel do collapse(4) map(alloc:state,flat,minus,plus) &
!$omp& private(i,j,k) &
!$omp& private(block,v,stencil) &
!$omp& private(cm,cp,wm) &
!$omp& private(wp,qm,qp) &
!$omp& private(f,macro_scope_10_local_s,macro_scope_10_local_ii) &
!$omp& private(macro_scope_10_local_jj,macro_scope_10_local_kk,macro_scope_13_local_beta) &
!$omp& private(macro_scope_13_local_correction,macro_scope_13_local_tau) &
!$omp& reduction(max:error)
    do block=1,nb
      do k=lower(3)-e(3),upper(3)+e(3)
      do j=lower(2)-e(2),upper(2)+e(2)
      do i=lower(1)-e(1),upper(1)+e(1)
        f=flat(i,j,k,block)
        do v=1,7
macro_scope_10: do
do macro_scope_10_local_s=-2,2
  macro_scope_10_local_ii=i+macro_scope_10_local_s*e(1); macro_scope_10_local_jj=j+macro_scope_10_local_s*e(2); &
  & macro_scope_10_local_kk=k+macro_scope_10_local_s*e(3)
  stencil(macro_scope_10_local_s)=state(macro_scope_10_local_ii,macro_scope_10_local_jj,macro_scope_10_local_kk, &
  & primitiveMap(v),block)
  if (v==QEINT) stencil(macro_scope_10_local_s)=stencil(macro_scope_10_local_s)*state(macro_scope_10_local_ii, &
  & macro_scope_10_local_jj,macro_scope_10_local_kk,primitiveMap(QRHO),block)
end do
exit macro_scope_10
end do macro_scope_10
macro_scope_13: do
cp(1)=(2*stencil(-2)-7*stencil(-1)+11*stencil(0))/6._rk
cp(2)=(-stencil(-1)+5*stencil(0)+2*stencil(1))/6._rk
cp(3)=(2*stencil(0)+5*stencil(1)-stencil(2))/6._rk
cm(1)=(-stencil(-2)+5*stencil(-1)+2*stencil(0))/6._rk
cm(2)=(2*stencil(-1)+5*stencil(0)-stencil(1))/6._rk
cm(3)=(11*stencil(0)-7*stencil(1)+2*stencil(2))/6._rk
macro_scope_13_local_beta(1)=13._rk/12*(stencil(-2)-2*stencil(-1)+stencil(0))**2+0.25_rk*(stencil(-2)-4*stencil(-1)+3 &
  & *stencil(0))**2
macro_scope_13_local_beta(2)=13._rk/12*(stencil(-1)-2*stencil(0)+stencil(1))**2+0.25_rk*(stencil(-1)-stencil(1))**2
macro_scope_13_local_beta(3)=13._rk/12*(stencil(0)-2*stencil(1)+stencil(2))**2+0.25_rk*(3*stencil(0)-4*stencil(1) &
  & +stencil(2))**2
macro_scope_13_local_tau=abs(macro_scope_13_local_beta(1)-macro_scope_13_local_beta(3)); &
  & macro_scope_13_local_correction=1._rk+macro_scope_13_local_tau/(macro_scope_13_local_beta+params%weno_epsilon)
wp=[0.1_rk,0.6_rk,0.3_rk]*macro_scope_13_local_correction; wp=wp/sum(wp)
wm=[0.3_rk,0.6_rk,0.1_rk]*macro_scope_13_local_correction; wm=wm/sum(wm)
exit macro_scope_13
end do macro_scope_13
macro_scope_14: do
qm=f*dot_product(wm,cm)+(1._rk-f)*(stencil(0))
qp=f*dot_product(wp,cp)+(1._rk-f)*(stencil(0))
if ((qp-(stencil(0)))*((stencil(0))-qm)<=0._rk) then
  qm=(stencil(0)); qp=(stencil(0))
end if
exit macro_scope_14
end do macro_scope_14
          minus(i,j,k,v,block)=qm; plus(i,j,k,v,block)=qp
        end do
      end do
    end do
      end do
      end do
    !$omp end target teams distribute parallel do
  end subroutine

  subroutine faces(state,minus,plus,flux,nx,ny,nz,nb,ng,lower,upper,params,d,error)
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    integer,intent(in)::nx,ny,nz,nb,ng,lower(3),upper(3),d
    real(rk),intent(in)::state(1-ng:nx+ng,1-ng:ny+ng,1-(NDIM-2)*ng:nz+(NDIM-2)*ng,9,nb),minus(1-ng:nx+ng,1-ng:ny+ng, &
  & 1-(NDIM-2)*ng:nz+(NDIM-2)*ng,7,nb),plus(1-ng:nx+ng,1-ng:ny+ng,1-(NDIM-2)*ng:nz+(NDIM-2)*ng,7,nb)
    real(rk),intent(inout)::flux(1-ng:nx+ng+1,1-ng:ny+ng+1,1-(NDIM-2)*ng:nz+(NDIM-2)*ng+(NDIM-2),5,NDIM,nb)
    type(hydro_parameters),intent(in)::params
    integer,intent(out)::error
    integer::inputMap(7),outputMap(5),shockVariable,e(3)
    integer::i,block,j,k,status
    real(rk)::ql(7),qr(7),ul(5),ur(5),fl(5),fr(5),usl(5),usr(5),f(5),sl,sr,speed,contact,pstar
    logical::vacuum,inShock
integer :: macro_scope_15_local_v
real(rk) :: macro_scope_20_local_al2,macro_scope_20_local_ar2
real(rk) :: macro_scope_23_local_un
real(rk) :: macro_scope_24_local_un
real(rk) :: macro_scope_29_local_unL,macro_scope_29_local_unR,macro_scope_29_local_leftMass, &
  & macro_scope_29_local_rightMass,macro_scope_29_local_denom
real(rk) :: macro_scope_30_local_numer,macro_scope_30_local_denom,macro_scope_30_local_un
real(rk) :: macro_scope_31_local_numer,macro_scope_31_local_denom,macro_scope_31_local_un
real(rk) :: macro_scope_33_local_ql(7),macro_scope_33_local_qr(7),macro_scope_33_local_ul(5),macro_scope_33_local_ur( &
  & 5),macro_scope_33_local_nu
integer :: macro_scope_34_local_v
integer :: macro_scope_36_local_v
integer :: macro_scope_41_local_v

    e=0; e(d)=1
    inputMap=primitiveMap; outputMap=conservedMap; shockVariable=9
    error=0
!$omp target teams distribute parallel do collapse(4) map(alloc:state,minus,plus,flux) &
!$omp& private(i,block,j) &
!$omp& private(k,status,ql) &
!$omp& private(qr,ul,ur) &
!$omp& private(fl,fr,usl) &
!$omp& private(usr,f,sl) &
!$omp& private(sr,speed,contact) &
!$omp& private(pstar,vacuum,inShock) &
!$omp& private(macro_scope_15_local_v,macro_scope_20_local_al2,macro_scope_20_local_ar2) &
!$omp& private(macro_scope_23_local_un,macro_scope_24_local_un,macro_scope_29_local_unL) &
!$omp& private(macro_scope_29_local_unR,macro_scope_29_local_leftMass,macro_scope_29_local_rightMass) &
!$omp& private(macro_scope_29_local_denom,macro_scope_30_local_numer,macro_scope_30_local_denom) &
!$omp& private(macro_scope_30_local_un,macro_scope_31_local_numer,macro_scope_31_local_denom) &
!$omp& private(macro_scope_31_local_un,macro_scope_33_local_ql,macro_scope_33_local_qr) &
!$omp& private(macro_scope_33_local_ul,macro_scope_33_local_ur,macro_scope_33_local_nu) &
!$omp& private(macro_scope_34_local_v,macro_scope_36_local_v,macro_scope_41_local_v) &
!$omp& reduction(max:error)
    do block=1,nb
      do k=lower(3),upper(3)+e(3)
      do j=lower(2),upper(2)+e(2)
      do i=lower(1),upper(1)+e(1)
        status=HY_OK
        cell_work: do
macro_scope_15: do
do macro_scope_15_local_v=1,7
  ql(macro_scope_15_local_v)=plus(i-e(1),j-e(2),k-e(3),macro_scope_15_local_v,block)
  qr(macro_scope_15_local_v)=minus(i,j,k,macro_scope_15_local_v,block)
end do
exit macro_scope_15
end do macro_scope_15
            if (.not.all(ieee_is_finite(ql)) .or. .not.all(ieee_is_finite(qr))) then
              status=HY_BAD_STATE; exit cell_work
            end if
macro_scope_18: do
ql(QRHO)=max(params%density_floor,ql(QRHO))
ql(QP)=max(params%pressure_floor,ql(QP))
exit macro_scope_18
end do macro_scope_18
macro_scope_19: do
qr(QRHO)=max(params%density_floor,qr(QRHO))
qr(QP)=max(params%pressure_floor,qr(QP))
exit macro_scope_19
end do macro_scope_19
macro_scope_20: do
status=HY_OK; vacuum=.false.; sl=0; sr=0; speed=0
if (.not.all(ieee_is_finite(ql)) .or. .not.all(ieee_is_finite(qr))) then
  status=HY_BAD_STATE; exit macro_scope_20
end if
if (min(ql(QRHO),qr(QRHO),ql(QP),qr(QP))<=0._rk) then
  status=HY_BAD_STATE; exit macro_scope_20
end if
vacuum=min(ql(QRHO),qr(QRHO),ql(QP),qr(QP))<params%tiny_state
if (vacuum) exit macro_scope_20
macro_scope_20_local_al2=ql(QGAMMA)*ql(QP)/ql(QRHO); macro_scope_20_local_ar2=qr(QGAMMA)*qr(QP)/qr(QRHO)
if (.not.ieee_is_finite(macro_scope_20_local_al2) .or. .not.ieee_is_finite(macro_scope_20_local_ar2)) then
  status=HY_BAD_STATE; exit macro_scope_20
end if
if (min(macro_scope_20_local_al2,macro_scope_20_local_ar2)<0._rk) then
  status=HY_BAD_STATE; exit macro_scope_20
end if
sl=min(ql(QUX+d-1)-sqrt(macro_scope_20_local_al2),qr(QUX+d-1)-sqrt(macro_scope_20_local_ar2))
sr=max(ql(QUX+d-1)+sqrt(macro_scope_20_local_al2),qr(QUX+d-1)+sqrt(macro_scope_20_local_ar2))
speed=max(abs(sl),abs(sr))
exit macro_scope_20
end do macro_scope_20
            if (status/=HY_OK) exit cell_work
            
            f=0
            if (.not.vacuum) then
macro_scope_21: do
ul(URHO)=ql(QRHO); ul(UMX:UMZ)=ql(QRHO)*ql(QUX:QUZ)
ul(UENER)=0.5_rk*ql(QRHO)*sum(ql(QUX:QUZ)**2)+ql(QEINT)
exit macro_scope_21
end do macro_scope_21
macro_scope_22: do
ur(URHO)=qr(QRHO); ur(UMX:UMZ)=qr(QRHO)*qr(QUX:QUZ)
ur(UENER)=0.5_rk*qr(QRHO)*sum(qr(QUX:QUZ)**2)+qr(QEINT)
exit macro_scope_22
end do macro_scope_22
macro_scope_23: do
macro_scope_23_local_un=ql(QUX+d-1)
fl=ul*macro_scope_23_local_un
fl(UMX+d-1)=fl(UMX+d-1)+ql(QP)
fl(UENER)=(ul(UENER)+ql(QP))*macro_scope_23_local_un
exit macro_scope_23
end do macro_scope_23
macro_scope_24: do
macro_scope_24_local_un=qr(QUX+d-1)
fr=ur*macro_scope_24_local_un
fr(UMX+d-1)=fr(UMX+d-1)+qr(QP)
fr(UENER)=(ur(UENER)+qr(QP))*macro_scope_24_local_un
exit macro_scope_24
end do macro_scope_24
              inShock=.false.
              if (params%hybrid .and. params%shock_marker_present) then
macro_scope_25: do
inShock=(state(i,j,k,shockVariable,block)/=0._rk .or. &
state(i-e(1),j-e(2),k-e(3),shockVariable,block)/=0._rk)
exit macro_scope_25
end do macro_scope_25
              end if
              if (inShock) then
macro_scope_28: do
status=HY_OK; f=0
if (sl>0._rk) then
  f=fl
else if (sr<0._rk) then
  f=fr
else if (sr==sl) then
  status=HY_DEGENERATE
else
  f=(sr*fl-sl*fr+sr*sl*(ur-ul))/(sr-sl)
end if
exit macro_scope_28
end do macro_scope_28
                if (status/=HY_OK) exit cell_work
              else if (sl>=0._rk) then
                f=fl
              else if (sr<0._rk) then
                f=fr
              else
macro_scope_29: do
status=HY_OK; contact=0; pstar=0
macro_scope_29_local_unL=ql(QUX+d-1); macro_scope_29_local_unR=qr(QUX+d-1)
macro_scope_29_local_leftMass=ql(QRHO)*(sl-macro_scope_29_local_unL); macro_scope_29_local_rightMass=qr(QRHO)*( &
  & sr-macro_scope_29_local_unR)
macro_scope_29_local_denom=macro_scope_29_local_rightMass-macro_scope_29_local_leftMass
if (macro_scope_29_local_denom==0._rk) then
  status=HY_DEGENERATE; exit macro_scope_29
end if
contact=(macro_scope_29_local_rightMass*macro_scope_29_local_unR-macro_scope_29_local_leftMass &
  & *macro_scope_29_local_unL+ql(QP)-qr(QP))/macro_scope_29_local_denom
pstar=macro_scope_29_local_leftMass*(contact-macro_scope_29_local_unL)+ql(QP)
if (.not.ieee_is_finite(contact) .or. .not.ieee_is_finite(pstar)) status=HY_BAD_STATE
exit macro_scope_29
end do macro_scope_29
                if (status/=HY_OK) exit cell_work
macro_scope_30: do
status=HY_OK; usl=0
macro_scope_30_local_un=ql(QUX+d-1); macro_scope_30_local_numer=sl-macro_scope_30_local_un; &
  & macro_scope_30_local_denom=sl-contact
if (macro_scope_30_local_denom==0._rk) then
  status=HY_DEGENERATE; exit macro_scope_30
end if
usl(URHO)=ql(QRHO)*macro_scope_30_local_numer/macro_scope_30_local_denom
usl(UMX:UMZ)=ul(UMX:UMZ)*macro_scope_30_local_numer/macro_scope_30_local_denom
usl(UMX+d-1)=usl(URHO)*contact
usl(UENER)=(ul(UENER)*macro_scope_30_local_numer+pstar*contact-ql(QP)*macro_scope_30_local_un) &
  & /macro_scope_30_local_denom
if (.not.all(ieee_is_finite(usl))) status=HY_BAD_STATE
exit macro_scope_30
end do macro_scope_30
                if (status/=HY_OK) exit cell_work
macro_scope_31: do
status=HY_OK; usr=0
macro_scope_31_local_un=qr(QUX+d-1); macro_scope_31_local_numer=sr-macro_scope_31_local_un; &
  & macro_scope_31_local_denom=sr-contact
if (macro_scope_31_local_denom==0._rk) then
  status=HY_DEGENERATE; exit macro_scope_31
end if
usr(URHO)=qr(QRHO)*macro_scope_31_local_numer/macro_scope_31_local_denom
usr(UMX:UMZ)=ur(UMX:UMZ)*macro_scope_31_local_numer/macro_scope_31_local_denom
usr(UMX+d-1)=usr(URHO)*contact
usr(UENER)=(ur(UENER)*macro_scope_31_local_numer+pstar*contact-qr(QP)*macro_scope_31_local_un) &
  & /macro_scope_31_local_denom
if (.not.all(ieee_is_finite(usr))) status=HY_BAD_STATE
exit macro_scope_31
end do macro_scope_31
                if (status/=HY_OK) exit cell_work
macro_scope_32: do
if (sl>=0._rk) then
  f=fl
else if (contact>=0._rk) then
  f=fl+sl*(usl-ul)
else if (sr>=0._rk) then
  f=fr+sr*(usr-ur)
else
  f=fr
end if
exit macro_scope_32
end do macro_scope_32
              end if
            end if
macro_scope_33: do
if (params%viscosity==0._rk) exit macro_scope_33
macro_scope_34: do
do macro_scope_34_local_v=1,7
  macro_scope_33_local_ql(macro_scope_34_local_v)=state((i-e(1)),(j-e(2)),(k-e(3)),inputMap(macro_scope_34_local_v), &
  & block)
end do
macro_scope_33_local_ql(QEINT)=macro_scope_33_local_ql(QRHO)*macro_scope_33_local_ql(QEINT)
exit macro_scope_34
end do macro_scope_34
macro_scope_36: do
do macro_scope_36_local_v=1,7
  macro_scope_33_local_qr(macro_scope_36_local_v)=state(i,j,k,inputMap(macro_scope_36_local_v),block)
end do
macro_scope_33_local_qr(QEINT)=macro_scope_33_local_qr(QRHO)*macro_scope_33_local_qr(QEINT)
exit macro_scope_36
end do macro_scope_36
macro_scope_38: do
macro_scope_33_local_nu=params%viscosity*max(macro_scope_33_local_ql(QUX+d-1)-macro_scope_33_local_qr(QUX+d-1),0._rk)
exit macro_scope_38
end do macro_scope_38
macro_scope_39: do
macro_scope_33_local_ul(URHO)=macro_scope_33_local_ql(QRHO); macro_scope_33_local_ul(UMX:UMZ)=macro_scope_33_local_ql( &
  & QRHO)*macro_scope_33_local_ql(QUX:QUZ)
macro_scope_33_local_ul(UENER)=0.5_rk*macro_scope_33_local_ql(QRHO)*sum(macro_scope_33_local_ql(QUX:QUZ)**2) &
  & +macro_scope_33_local_ql(QEINT)
exit macro_scope_39
end do macro_scope_39
macro_scope_40: do
macro_scope_33_local_ur(URHO)=macro_scope_33_local_qr(QRHO); macro_scope_33_local_ur(UMX:UMZ)=macro_scope_33_local_qr( &
  & QRHO)*macro_scope_33_local_qr(QUX:QUZ)
macro_scope_33_local_ur(UENER)=0.5_rk*macro_scope_33_local_qr(QRHO)*sum(macro_scope_33_local_qr(QUX:QUZ)**2) &
  & +macro_scope_33_local_qr(QEINT)
exit macro_scope_40
end do macro_scope_40
f=f+macro_scope_33_local_nu*(macro_scope_33_local_ul-macro_scope_33_local_ur)
exit macro_scope_33
end do macro_scope_33
            if (.not.all(ieee_is_finite(f))) then
              status=HY_BAD_STATE; exit cell_work
            end if
macro_scope_41: do
do macro_scope_41_local_v=1,5
  flux(i,j,k,outputMap(macro_scope_41_local_v),d,block)=f(macro_scope_41_local_v)
end do
exit macro_scope_41
end do macro_scope_41
          exit cell_work
        end do cell_work
        error=max(error,status)
      end do
    end do
      end do
      end do
    !$omp end target teams distribute parallel do
  end subroutine

  subroutine update(reference,state,flux,nx,ny,nz,nb,ng,lower,upper,dx,dy,dz,dt,coeff,densityFloor,energyFloor,error)

    integer,intent(in)::nx,ny,nz,nb,ng,lower(3),upper(3)
    real(rk)::spacing(3)
    real(rk),intent(in)::reference(1-ng:nx+ng,1-ng:ny+ng,1-(NDIM-2)*ng:nz+(NDIM-2)*ng,9,nb),flux(1-ng:nx+ng+1, &
  & 1-ng:ny+ng+1,1-(NDIM-2)*ng:nz+(NDIM-2)*ng+(NDIM-2),5,NDIM,nb),dx,dy,dz,dt,coeff(3),densityFloor,energyFloor
    real(rk),intent(inout)::state(1-ng:nx+ng,1-ng:ny+ng,1-(NDIM-2)*ng:nz+(NDIM-2)*ng,9,nb)
    integer,intent(out)::error
    integer::i,j,k,block,v
    real(rk)::u0(5),us(5),source(5),divergence(5),updated(5),rho,velocity(3),totalEnergy,internalEnergy
integer :: macro_scope_43_local_d
real(rk) :: macro_scope_43_local_rho
integer :: macro_scope_47_local_d
real(rk) :: macro_scope_47_local_rho
integer :: macro_scope_51_local_v
integer :: macro_scope_60_local_d

    spacing=[dx,dy,dz]
    error=0
!$omp target teams distribute parallel do collapse(4) map(alloc:reference,state,flux) &
!$omp& private(i,j,k) &
!$omp& private(block,v,u0) &
!$omp& private(us,source,divergence) &
!$omp& private(updated,rho,velocity) &
!$omp& private(totalEnergy,internalEnergy,macro_scope_43_local_d) &
!$omp& private(macro_scope_43_local_rho,macro_scope_47_local_d,macro_scope_47_local_rho) &
!$omp& private(macro_scope_51_local_v,macro_scope_60_local_d) &
!$omp& reduction(max:error)
    do block=1,nb
      do k=lower(3),upper(3)
      do j=lower(2),upper(2)
      do i=lower(1),upper(1)
macro_scope_43: do
macro_scope_43_local_rho=reference(i,j,k,updateMap(1),block)
u0(1)=macro_scope_43_local_rho
do macro_scope_43_local_d=1,3
  u0(1+macro_scope_43_local_d)=macro_scope_43_local_rho*reference(i,j,k,updateMap(1+macro_scope_43_local_d),block)
end do
u0(5)=macro_scope_43_local_rho*reference(i,j,k,updateMap(5),block)
exit macro_scope_43
end do macro_scope_43
macro_scope_47: do
macro_scope_47_local_rho=state(i,j,k,updateMap(1),block)
us(1)=macro_scope_47_local_rho
do macro_scope_47_local_d=1,3
  us(1+macro_scope_47_local_d)=macro_scope_47_local_rho*state(i,j,k,updateMap(1+macro_scope_47_local_d),block)
end do
us(5)=macro_scope_47_local_rho*state(i,j,k,updateMap(5),block)
exit macro_scope_47
end do macro_scope_47
        source=0
macro_scope_51: do
do macro_scope_51_local_v=1,5
  divergence(macro_scope_51_local_v)=(flux(i+1,j,k,conservedMap(macro_scope_51_local_v),1,block) - &
                 flux(i,j,k,conservedMap(macro_scope_51_local_v),1,block))/spacing(1)
  if (NDIM>1) then
    divergence(macro_scope_51_local_v)=divergence(macro_scope_51_local_v)+(flux(i,j+1,k,conservedMap( &
  & macro_scope_51_local_v),2,block) - &
                                flux(i,j,k,conservedMap(macro_scope_51_local_v),2,block))/spacing(2)
  end if
  if (NDIM>2) then
    divergence(macro_scope_51_local_v)=divergence(macro_scope_51_local_v)+(flux(i,j,k+1,conservedMap( &
  & macro_scope_51_local_v),3,block) - &
                                flux(i,j,k,conservedMap(macro_scope_51_local_v),3,block))/spacing(3)
  end if
end do
exit macro_scope_51
end do macro_scope_51
macro_scope_58: do
updated=coeff(1)*u0+coeff(2)*us+(coeff(3)*dt)*(source-divergence)
exit macro_scope_58
end do macro_scope_58
macro_scope_59: do
rho=max(updated(1),densityFloor)
velocity=updated(2:4)/rho
totalEnergy=max(energyFloor,updated(5)/rho)
internalEnergy=max(energyFloor,totalEnergy-0.5_rk*sum(velocity**2))
exit macro_scope_59
end do macro_scope_59
macro_scope_60: do
state(i,j,k,updateMap(1),block)=rho
do macro_scope_60_local_d=1,3
  state(i,j,k,updateMap(1+macro_scope_60_local_d),block)=velocity(macro_scope_60_local_d)
end do
state(i,j,k,updateMap(5),block)=totalEnergy
state(i,j,k,updateMap(6),block)=internalEnergy
exit macro_scope_60
end do macro_scope_60
        if (.not.all(ieee_is_finite(updated))) error=max(error,HY_BAD_STATE)
      end do
    end do
      end do
      end do
    !$omp end target teams distribute parallel do
  end subroutine
end module
