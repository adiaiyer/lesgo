!!
!!  Copyright (C) 2009-2017  Johns Hopkins University
!!
!!  This file is part of lesgo.
!!
!!  lesgo is free software: you can redistribute it and/or modify
!!  it under the terms of the GNU General Public License as published by
!!  the Free Software Foundation, either version 3 of the License, or
!!  (at your option) any later version.
!!
!!  lesgo is distributed in the hope that it will be useful,
!!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!  GNU General Public License for more details.
!!
!!  You should have received a copy of the GNU General Public License
!!  along with lesgo.  If not, see <http://www.gnu.org/licenses/>.
!!

!*******************************************************************************
subroutine wallstress
!*******************************************************************************
!
! This subroutine calculates the wall stress txz, tyz (w-nodes) and dudz,
! dvdz (w-nodes) at the first z-location k = 1. The wall stress is calculated
! depending on lower boundary condition lbc_mom. This subroutine should only
! be called after ensuring coord==0
!
! Options for lbc_mom:
!   0 - stress free
!       txz, tyz, dudz, and dvdz are all 0
!
!   1 - DNS wall boundary conditions
!       calculates wall stress values from the first grid point
!
!   2 - Equilibirum wall model
!       See John D. Albertson's dissertation, eqns (2.46)-(2.52)
!       Also see E. Bou-Zeid, C. Meneveau & M.B. Parlange, "A scale-dependent
!           Lagrangian dynamic model for large eddy simulation of complex
!           turbulent flows" (2005) -- Appendix
!
!   3 - Integral wall model
!       See X.I.A. Yang, J. Sadique, R. Mittal & C. Meneveau, "Integral wall
!           model for large eddy simulations of wall-bounded turbulent flows." (2015)
!
use types, only : rprec
use param, only : lbc_mom
use param, only : ubc_mom, coord, nproc, nz ! these necessary only for upper bc
use messages, only : error
use iwmles, only : iwm_wallstress
use sim_param, only : txz, tyz, dudz, dvdz
implicit none
character(*), parameter :: sub_name = 'wallstress'

! Lower boundary condition
if (coord == 0) then
    select case (lbc_mom)
        ! Stress free
        case (0)
            call ws_free_lbc

        ! DNS wall
        case (1)
            call ws_dns_lbc

        ! Equilibrium wall model
        case (2)
            call ws_equilibrium_lbc

        ! Integral wall model (not implemented for top wall)
        case (3)
            call iwm_wallstress()

        ! Otherwise, invalid
        case default
            call error (sub_name, 'invalid lbc_mom')
    end select
end if

if (coord == nproc-1) then
    select case (ubc_mom)
        ! Stress free
        case (0)
            call ws_free_ubc

        ! DNS wall
        case (1)
            call ws_dns_ubc

        ! Equilibrium wall model
        case (2)
            call ws_equilibrium_ubc

        ! Integral wall model (not implemented for top wall)
        case (3)
            call error(sub_name, 'invalid ubc_mom')

        ! Otherwise, invalid
        case default
            call error(sub_name, 'invalid ubc_mom')
    end select
end if

contains

!*******************************************************************************
subroutine ws_free_lbc
!*******************************************************************************
implicit none

txz(:, :, 1) = 0._rprec
tyz(:, :, 1) = 0._rprec
dudz(:, :, 1) = 0._rprec
dvdz(:, :, 1) = 0._rprec

end subroutine ws_free_lbc

!*******************************************************************************
subroutine ws_free_ubc
!*******************************************************************************
implicit none

txz(:, :,nz) = 0._rprec
tyz(:, :,nz) = 0._rprec
dudz(:,:,nz) = 0._rprec
dvdz(:,:,nz) = 0._rprec

end subroutine ws_free_ubc

!*******************************************************************************
subroutine ws_dns_lbc
!*******************************************************************************
use param, only : nx, ny, nu_molec, z_i, u_star, dz
use param, only : ubot
use sim_param , only : u, v
implicit none
integer :: i, j

do j = 1, ny
    do i = 1, nx
        dudz(i,j,1) = ( u(i,j,1) - ubot ) / (0.5_rprec*dz)
        dvdz(i,j,1) = v(i,j,1) / (0.5_rprec*dz)
        txz(i,j,1) = -nu_molec/(z_i*u_star)*dudz(i,j,1)
        tyz(i,j,1) = -nu_molec/(z_i*u_star)*dvdz(i,j,1)
    end do
end do

end subroutine ws_dns_lbc

!*******************************************************************************
subroutine ws_dns_ubc
!*******************************************************************************
use param, only : nx, ny, nu_molec, z_i, u_star, dz
use param, only : utop
use sim_param , only : u, v
implicit none
integer :: i, j

do j = 1, ny
    do i = 1, nx
        dudz(i,j,nz) = ( utop - u(i,j,nz-1) ) / (0.5_rprec*dz)
        dvdz(i,j,nz) = -v(i,j,nz-1) / (0.5_rprec*dz)
        txz(i,j,nz) = -nu_molec/(z_i*u_star)*dudz(i,j,nz)
        tyz(i,j,nz) = -nu_molec/(z_i*u_star)*dvdz(i,j,nz)
    end do
end do

end subroutine ws_dns_ubc

!*******************************************************************************
subroutine ws_equilibrium_lbc
!*******************************************************************************
use param, only : coord, read_endian, dz, ld, nx, ny, vonk, zo, use_sea_drag_model, jt_total, total_time, path, nsteps_wavy
use param, only : sea_drag_io_flag, sea_drag_io_nstart, sea_drag_io_nend, sea_drag_io_nskip, write_endian,is_swell
use sim_param, only : u, v, ustar_lbc
use sea_surface_drag_model
use string_util
use test_filtermodule
#ifdef PPSCALARS
use scalars, only : obukhov, phi_m
#endif

implicit none

integer :: i, j
character (64) :: fname
real(rprec), dimension(nx, ny) :: denom, u_avg
real(rprec), dimension(ld, ny) :: u1, v1
real(rprec) :: const, time_wavy,factor,zloc,swell_term
logical :: exst

if(use_sea_drag_model) then
        call sea_surface_drag_model_forces()
        !do i=1,nx
        u1 = u_rel(:,:)
        v1 = v_rel(:,:)
        !end do
!! AA BOC
!! We should not use the first grid point for stress calculaktions., Ideally use U10?.
!! For now lets just use the third grid point. Only will change the wave routines
!! z_3 = 5*dz/2    
 zloc  = 2.5_rprec*dz
        denom = log(zloc/zo-eta(1:nx,:)/zo)

        do i=1,nx
        do j=1,ny
           if(isnan(denom(i,j))) then
                   print *,"denom has NaN"
                   stop
           elseif(denom(i,j).eq.0)then
                   print *,"denom has zero:i,j,dz,eta,dz-eta,zo,log(dz-eta)",i,j,zloc,eta(i,j) &
                                                      ,zloc-eta(i,j),zo,log(zloc/zo-eta(i,j)/zo)
           endif
        enddo
        enddo
        if (jt_total .eq. nsteps_wavy) then
              time_wavy = total_time
              if (coord==0) then
                  open(12, file='time_wavy.out', form='unformatted', convert=write_endian)
                  write(12) time_wavy
                  close(12)
                  if (coord==0) print*,"saved time_wavy", coord, jt_total, nsteps_wavy, time_wavy
              endif
        endif
        if (jt_total .gt. nsteps_wavy) then
              inquire (file='time_wavy.out', exist=exst)
              if (exst) then
                  open(13, file='time_wavy.out', form='unformatted', convert=read_endian)
                  read(13) time_wavy
                  close(13)
                  if(coord==0) print*,"time_wavy",jt_total, nsteps_wavy, time_wavy
              end if
        endif
else
        u1 = u(:,:,1)
        v1 = v(:,:,1)
        denom = log(0.5_rprec*dz/zo)
endif

call test_filter(u1)
call test_filter(v1)


u_avg = sqrt(u1(1:nx,1:ny)**2+v1(1:nx,1:ny)**2)
#ifdef PPSCALARS
call obukhov(u_avg,denom)
#else
ustar_lbc = u_avg*vonk/denom
#endif
!ustar_av = sum(ustar_lbc)/(nx*ny)
do j = 1, ny
    do i = 1, nx
        const = -(ustar_lbc(i,j)**2)/u_avg(i,j)
        
       
        txz(i,j,1) = const*u1(i,j)
        tyz(i,j,1) = const*v1(i,j)
        ustar_lbc(i,j) = sqrt(sqrt(txz(i,j,1)**2+tyz(i,j,1)**2))  ! An aproximation of u_* based on tau_t
        if (use_sea_drag_model) then
           if (is_swell .and. (u_rel_c(i,j) .lt. 0.0_rprec) .and. (fd_u(i,j) .gt. 1.0E-8_rprec) )   then
                !! This assumes wave propagating in x-direction. This is ok for now as we tilt the flow and
                !! not the wave.
                swell_term = -0.5_rprec/dz*(24.975_rprec*ustar_lbc(i,j)**2 - 0.964*(c_phase*ustar_lbc(i,j)))*(ak)**2
!                if (coord .eq.0) write(*,*) "i,j,urel,ustar,swell,wdm",i,j,u_rel_c(i,j),ustar_lbc(i,j), swell_term, fd_u(i,j)
                fd_u(i,j) = swell_term 
            endif
           txz(i,j,1) = txz(i,j,1) + fd_u(i,j)*dz*(1-exp(-(total_time-time_wavy)**2))
           tyz(i,j,1) = tyz(i,j,1) + fd_v(i,j)*dz*(1-exp(-(total_time-time_wavy)**2))
           ustar_lbc(i,j) = sqrt(sqrt(txz(i,j,1)**2+tyz(i,j,1)**2))! Recalculate u_* based on total stress
           !this is as in Moeng 84
           
! AA Displacenment height correction to gradientsS
#ifdef PPSCALARS
           dudz(i,j,1) = ustar_lbc(i,j)/((zloc-eta(i,j))*vonK)*u1(i,j)/u_avg(i,j)   &
            * phi_m(i,j) 
           dvdz(i,j,1) = ustar_lbc(i,j)/((zloc-eta(i,j))*vonK)*v1(i,j)/u_avg(i,j)   &
            * phi_m(i,j)
#else
           dudz(i,j,1) = ustar_lbc(i,j)/((zloc-eta(i,j))*vonK)*u1(i,j)/u_avg(i,j)
           dvdz(i,j,1) = ustar_lbc(i,j)/((zloc-eta(i,j))*vonK)*v1(i,j)/u_avg(i,j)
#endif
           dudz(i,j,1) = merge(0._rprec,dudz(i,j,1),u_rel(i,j).eq.0._rprec)
           dvdz(i,j,1) = merge(0._rprec,dvdz(i,j,1),v_rel(i,j).eq.0._rprec)
        else
           !this is as in Moeng 84
#ifdef PPSCALAR
           dudz(i,j,1) = ustar_lbc(i,j)/((0.5_rprec*dz)*vonK)*u(i,j,1)/u_avg(i,j)   &
            * phi_m(i,j)
           dvdz(i,j,1) = ustar_lbc(i,j)/(0.5_rprec*dz*vonK)*v(i,j,1)/u_avg(i,j)   &
            * phi_m(i,j)
#else
           dudz(i,j,1) = ustar_lbc(i,j)/(0.5_rprec*dz*vonK)*u(i,j,1)/u_avg(i,j)
           dvdz(i,j,1) = ustar_lbc(i,j)/(0.5_rprec*dz*vonK)*v(i,j,1)/u_avg(i,j)
#endif
           dudz(i,j,1) = merge(0._rprec,dudz(i,j,1),u(i,j,1).eq.0._rprec)
           dvdz(i,j,1) = merge(0._rprec,dvdz(i,j,1),v(i,j,1).eq.0._rprec)
        endif
    end do
end do

if(use_sea_drag_model) then
if (sea_drag_io_flag) then
    if ((jt_total >= sea_drag_io_nstart).and.(jt_total <= sea_drag_io_nend)) then
        if ( mod(jt_total-sea_drag_io_nstart,sea_drag_io_nskip)==0 ) then
                call string_splice(fname, path //'output/sea_drag_io.', jt_total)
                ! Write binary Output
                call string_concat(fname, '.bin')
                open(unit=13, file=fname, form='unformatted', convert=write_endian,        &
                     access='direct', recl=nx*ny*rprec)
                write(13,rec=1)    u(1:nx,1:ny,1)
                write(13,rec=2)    v(1:nx,1:ny,1)
                write(13,rec=3)   u1(1:nx,1:ny)
                write(13,rec=4)   v1(1:nx,1:ny)
                write(13,rec=5)  txz(1:nx,1:ny,1)
                write(13,rec=6)  tyz(1:nx,1:ny,1)
                write(13,rec=7)  eta(1:nx,1:ny)
                write(13,rec=8) detadx(1:nx,1:ny)
                write(13,rec=9) detady(1:nx,1:ny)
                write(13,rec=10) us_orb(1:nx,1:ny)
                write(13,rec=11) ws_orb(1:nx,1:ny)
                write(13,rec=12) u_rel(1:nx,1:ny)
                write(13,rec=13) v_rel(1:nx,1:ny)
                write(13,rec=14) w_rel(1:nx,1:ny)
                write(13,rec=15) u_rel_c(1:nx,1:ny)
                write(13,rec=16) v_rel_c(1:nx,1:ny)
                write(13,rec=17) n_u(1:nx,1:ny)
                write(13,rec=18) n_v(1:nx,1:ny)
                write(13,rec=19) fd_u(1:nx,1:ny)
                write(13,rec=20) fd_v(1:nx,1:ny) 
                write(13,rec=21) total_time
                close(13)
        endif
    endif
endif
endif

end subroutine ws_equilibrium_lbc

!!*******************************************************************************
!subroutine ws_equilibrium_waves_lbc
!!*******************************************************************************
!use param, only : dz, ld, nx, ny, vonk, zo
!use sim_param, only : u, v, ustar_lbc
!use test_filtermodule
!#ifdef PPSCALARS
!use scalars, only : obukhov, phi_m
!#endif

!implicit none

!integer :: i, j
!real(rprec), dimension(nx, ny) :: denom, u_avg
!real(rprec), dimension(ld, ny) :: u1, v1
!real(rprec) :: const

!u1 = u(:,:,1)
!v1 = v(:,:,1)
!call test_filter(u1)
!call test_filter(v1)
!denom = log(0.5_rprec*dz/zo)
!u_avg = sqrt(u1(1:nx,1:ny)**2+v1(1:nx,1:ny)**2)
!#ifdef PPSCALARS
!call obukhov(u_avg)
!#else
!ustar_lbc = u_avg*vonk/denom
!#endif

!do j = 1, ny
!    do i = 1, nx
!        const = -(ustar_lbc(i,j)**2)/u_avg(i,j)
!        txz(i,j,1) = const*u1(i,j)
!        tyz(i,j,1) = const*v1(i,j)
!        !this is as in Moeng 84
!#ifdef PPSCALARS
!        dudz(i,j,1) = ustar_lbc(i,j)/(0.5_rprec*dz*vonK)*u(i,j,1)/u_avg(i,j)   &
!            * phi_m(i,j)
!        dvdz(i,j,1) = ustar_lbc(i,j)/(0.5_rprec*dz*vonK)*v(i,j,1)/u_avg(i,j)   &
!            * phi_m(i,j)
!#else
!        dudz(i,j,1) = ustar_lbc(i,j)/(0.5_rprec*dz*vonK)*u(i,j,1)/u_avg(i,j)
!        dvdz(i,j,1) = ustar_lbc(i,j)/(0.5_rprec*dz*vonK)*v(i,j,1)/u_avg(i,j)
!#endif
!        dudz(i,j,1) = merge(0._rprec,dudz(i,j,1),u(i,j,1).eq.0._rprec)
!        dvdz(i,j,1) = merge(0._rprec,dvdz(i,j,1),v(i,j,1).eq.0._rprec)
!    end do
!end do
!end subroutine

!*******************************************************************************
!subroutine sea_surface_drag_model() !GN
!*******************************************************************************
!1use param, only : ak, c_by_ustar, wave_angle, u_star, pi, z_i, total_time, dz, nx, ny
!use grid_m, only : grid
!use sim_param, only : u, v, w
!use functions, only : heaviside
!implicit none
!real(rprec), public, dimension(:,:), allocatable :: eta, detadx, detady, us_orb, vs_orb, ws_orb
!real(rprec), public, dimension(:,:), allocatable :: u_rel, v_rel, w_rel, u_rel_c, v_rel_c, w_rel_c 
!real(rprec), public, dimension(:,:), allocatable :: n_u, n_v, n_w, fd_u, fd_v
!real(rprec) :: a_amp, kx_wavno, ky_wavno, k_wavno,  cx_phase, cy_phase, c_phase, omega_freq,  g_accl
!integer :: i
!allocate ( eta   (nx,ny)) ; eta    = 0._rprec
!allocate ( detadx(nx,ny)) ; detadx = 0._rprec
!allocate ( detady(nx,ny)) ; detady = 0._rprec

!allocate ( us_orb(nx,ny)) ; us_orb = 0._rprec
!allocate ( vs_orb(nx,ny)) ; vs_orb = 0._rprec
!allocate ( ws_orb(nx,ny)) ; ws_orb = 0._rprec

!allocate ( u_rel (nx,ny)) ; u_rel  = 0._rprec
!allocate ( v_rel (nx,ny)) ; v_rel  = 0._rprec
!allocate ( w_rel (nx,ny)) ; w_rel  = 0._rprec

!allocate ( u_rel_c (nx,ny)) ; u_rel_c  = 0._rprec
!allocate ( v_rel_c (nx,ny)) ; v_rel_c  = 0._rprec
!allocate ( w_rel_c (nx,ny)) ; w_rel_c  = 0._rprec

!allocate ( n_u (nx,ny)) ; n_u  = 0._rprec
!allocate ( n_v (nx,ny)) ; n_v  = 0._rprec
!allocate ( n_w (nx,ny)) ; n_w  = 0._rprec

!allocate ( fd_u(nx,ny)) ; fd_u = 0._rprec
!allocate ( fd_v(nx,ny)) ; fd_v = 0._rprec


!g_accl         = 9.81_rprec
!c_phase        = c_by_ustar*u_star
!cx_phase       = c_phase*cos(wave_angle)
!cy_phase       = c_phase*u_star*sin(wave_angle)
!k_wavno        = g_accl/c_phase**2
!kx_wavno       = k_wavno*cos(wave_angle)
!ky_wavno       = k_wavno*sin(wave_angle)
!a_amp          = ak*k_wavno
!omega_freq     = c_phase*k_wavno

!g_accl         = g_accl*(z_i/(u_star**2))
!c_phase        = c_phase/u_star
!cx_phase       = cx_phase/u_star
!cy_phase       = cy_phase/u_star  
!kx_wavno       = kx_wavno*z_i
!ky_wavno       = ky_wavno*z_i
!k_wavno        = k_wavno*z_i
!a_amp          = a_amp/z_i
!omega_freq     = omega_freq*(z_i/u_star)

!do i = 1, nx
!  eta(i,:)    =              a_amp*cos( kx_wavno*grid%x(i) + ky_wavno*grid%y(1:ny) - omega_freq*total_time )
!  detadx(i,:) =    -a_amp*kx_wavno*sin( kx_wavno*grid%x(i) + ky_wavno*grid%y(1:ny) - omega_freq*total_time )
!  detady(i,:) =    -a_amp*ky_wavno*sin( kx_wavno*grid%x(i) + ky_wavno*grid%y(1:ny) - omega_freq*total_time )
!  us_orb(i,:) =   a_amp*omega_freq*cos( kx_wavno*grid%x(i) + ky_wavno*grid%y(1:ny) - omega_freq*total_time )
!  ws_orb(i,:) =   a_amp*omega_freq*sin( kx_wavno*grid%x(i) + ky_wavno*grid%y(1:ny) - omega_freq*total_time )

!  u_rel(i,:)  =  u(i,:,1) - us_orb(i,:)
!  v_rel(i,:)  =  v(i,:,1) - vs_orb(i,:)
!  w_rel(i,:)  =  0.5*(w(i,:,1)+w(i,:,2)) - ws_orb(i,:)

!  u_rel_c(i,:)=  u(i,:,1) - cx_phase
!  v_rel_c(i,:)=  v(i,:,1) - cy_phase 
!  w_rel_c(i,:)=  0.5*(w(i,:,1)+w(i,:,2))


!  n_u(i,:)    = u_rel_c(i,:)/SQRT(u_rel_c(i,:)**2+v_rel_c(i,:)**2)
!  n_v(i,:)    = v_rel_c(i,:)/SQRT(u_rel_c(i,:)**2+v_rel_c(i,:)**2)
!  n_w(i,:)    = w_rel_c(i,:)/SQRT(u_rel_c(i,:)**2+v_rel_c(i,:)**2)

!  fd_u(i,:) = -1.2/(1+6*a_amp**2*k_wavno**2)*(a_amp*k_wavno)*u(i,:,1)/dz*SQRT(u_rel_c(i,:)**2+v_rel_c(i,:)**2) &
!               *(n_u(i,:)*detadx(i,:)+n_v(i,:)*detady(i,:))*heaviside(n_u(i,:)*detadx(i,:)+n_v(i,:)*detady(i,:));
  
  !fd_v(i,:) = -1.2/(1+6*a_amp**2*k_wavno**2)*(a_amp*k_wavno)*v(i,:,1)/dz*SQRT(u_rel_c(i,:)**2+v_rel_c(i,:)**2) &
 !              *(n_u(i,:)*detadx(i,:)+n_v(i,:)*detady(i,:))*heaviside(n_u(i,:)*detadx(i,:)+n_v(i,:)*detady(i,:));

!end do


!end subroutine
!*******************************************************************************
subroutine ws_equilibrium_ubc
!*******************************************************************************
use param, only : dz, ld, nx, ny, vonk, zo
use sim_param, only : u, v
use test_filtermodule
implicit none
integer :: i, j
real(rprec), dimension(nx, ny) :: denom, u_avg, ustar
real(rprec), dimension(ld, ny) :: u1, v1
real(rprec) :: const


u1 = u(:,:,nz-1)
v1 = v(:,:,nz-1)
call test_filter(u1)
call test_filter(v1)
denom = log(0.5_rprec*dz/zo)
u_avg = sqrt(u1(1:nx,1:ny)**2+v1(1:nx,1:ny)**2)
ustar = u_avg*vonk/denom

do j = 1, ny
    do i = 1, nx
        const = (ustar(i,j)**2)/u_avg(i,j) ! diff sign for upper b.c.
        txz(i,j,nz) = const*u1(i,j)
        tyz(i,j,nz) = const*v1(i,j)
        !this is as in Moeng 84
        dudz(i,j,nz) = -ustar(i,j)/(0.5_rprec*dz*vonK)*u(i,j,nz-1)/u_avg(i,j)
        dvdz(i,j,nz) = -ustar(i,j)/(0.5_rprec*dz*vonK)*v(i,j,nz-1)/u_avg(i,j)
        dudz(i,j,nz) = merge(0._rprec,dudz(i,j,nz),u(i,j,nz-1).eq.0._rprec)
        dvdz(i,j,nz) = merge(0._rprec,dvdz(i,j,nz),v(i,j,nz-1).eq.0._rprec)
    end do
end do

end subroutine ws_equilibrium_ubc

end subroutine wallstress
