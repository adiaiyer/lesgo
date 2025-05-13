!*******************************************************************************
module sea_surface_drag_model
!*******************************************************************************
use types, only : rprec
use param, only : ld, nx, ny, nz, lbz, u_star, z_i, dz, total_time
use param, only : ak, c_by_ustar, wave_angle, pi

implicit none

save
public

public :: sea_surface_drag_model_init, sea_surface_drag_model_forces

real(rprec), public, dimension(:,:), allocatable :: eta, detadx, detady, us_orb, vs_orb, ws_orb
real(rprec), public, dimension(:,:), allocatable :: u_rel, v_rel, w_rel, u_rel_c, v_rel_c, w_rel_c
real(rprec), public, dimension(:,:), allocatable :: n_u, n_v, n_w, fd_u, fd_v
real(rprec) :: a_amp, kx_wavno, ky_wavno, k_wavno,  cx_phase, cy_phase, c_phase, omega_freq,  g_accl

contains

!*******************************************************************************
subroutine sea_surface_drag_model_init() !GN
!*******************************************************************************
use param, only : ld, nx,ny,z_i,u_star, ak, c_by_ustar, dz, coord
implicit none
allocate ( eta   (ld,ny)) ; eta    = 0._rprec
allocate ( detadx(ld,ny)) ; detadx = 0._rprec
allocate ( detady(ld,ny)) ; detady = 0._rprec

allocate ( us_orb(ld,ny)) ; us_orb = 0._rprec
allocate ( vs_orb(ld,ny)) ; vs_orb = 0._rprec
allocate ( ws_orb(ld,ny)) ; ws_orb = 0._rprec

allocate ( u_rel (ld,ny)) ; u_rel  = 0._rprec
allocate ( v_rel (ld,ny)) ; v_rel  = 0._rprec
allocate ( w_rel (ld,ny)) ; w_rel  = 0._rprec

allocate ( u_rel_c (ld,ny)) ; u_rel_c  = 0._rprec
allocate ( v_rel_c (ld,ny)) ; v_rel_c  = 0._rprec
allocate ( w_rel_c (ld,ny)) ; w_rel_c  = 0._rprec

allocate ( n_u (ld,ny)) ; n_u  = 0._rprec
allocate ( n_v (ld,ny)) ; n_v  = 0._rprec
allocate ( n_w (ld,ny)) ; n_w  = 0._rprec

allocate ( fd_u(ld,ny)) ; fd_u = 0._rprec
allocate ( fd_v(ld,ny)) ; fd_v = 0._rprec


g_accl         = 9.81_rprec
c_phase        = c_by_ustar*u_star
cx_phase       = c_phase*cos(wave_angle)
cy_phase       = c_phase*sin(wave_angle)
k_wavno        = g_accl/c_phase**2
kx_wavno       = k_wavno*cos(wave_angle)
ky_wavno       = k_wavno*sin(wave_angle)
a_amp          = ak/k_wavno
omega_freq     = c_phase*k_wavno

if (coord .eq. 0) then
   print *, 'c,cx,cy,k,kx,ky,a,omega (dimensional)',c_phase,cx_phase,cy_phase,k_wavno,kx_wavno,ky_wavno,a_amp,omega_freq
endif

g_accl         = g_accl*(z_i/(u_star**2))
c_phase        = c_phase/u_star
cx_phase       = cx_phase/u_star
cy_phase       = cy_phase/u_star
kx_wavno       = kx_wavno*z_i
ky_wavno       = ky_wavno*z_i
k_wavno        = k_wavno*z_i
a_amp          = a_amp/z_i
omega_freq     = omega_freq*(z_i/u_star)

if (coord .eq. 0) then
   print *, 'c,cx,cy,k,kx,ky,a,omega (non-dimensional)',c_phase,cx_phase,cy_phase,k_wavno,kx_wavno,ky_wavno,a_amp,omega_freq
endif

if (dz/2<a_amp) then
    print *,"dz < a_amp. Check grid size: dz/2, a_amp", dz*0.5,a_amp
    STOP
endif
end subroutine

!*******************************************************************************
subroutine sea_surface_drag_model_forces() !GN
!*******************************************************************************
use param, only : ld, ak, c_by_ustar, wave_angle, u_star, pi, z_i, total_time, dz, nx, ny, &
                  wave_orientation
use grid_m, only : grid
use sim_param, only : u, v, w
use functions, only : heaviside, heaviside_scalar
implicit none
integer :: i,j

!call sea_surface_drag_model_init()


SELECT CASE (WAVE_ORIENTATION)

  CASE (2)       

    !WIND OPPOSING WAVES        
    do i = 1, nx
    do j = 1, ny
      eta(i,j)    =              a_amp*cos( kx_wavno*grid%x(i) + ky_wavno*grid%y(j) + omega_freq*total_time )
      detadx(i,j) =    -a_amp*kx_wavno*sin( kx_wavno*grid%x(i) + ky_wavno*grid%y(j) + omega_freq*total_time )
      detady(i,j) =    -a_amp*ky_wavno*sin( kx_wavno*grid%x(i) + ky_wavno*grid%y(j) + omega_freq*total_time )
      us_orb(i,j) =   a_amp*omega_freq*cos( kx_wavno*grid%x(i) + ky_wavno*grid%y(j) + omega_freq*total_time )
      ws_orb(i,j) =   a_amp*omega_freq*sin( kx_wavno*grid%x(i) + ky_wavno*grid%y(j) + omega_freq*total_time )
    
      u_rel(i,j)  =  u(i,j,1) + us_orb(i,j)
      v_rel(i,j)  =  v(i,j,1) + vs_orb(i,j)
      w_rel(i,j)  =  w(i,j,1) + ws_orb(i,j)
    
      u_rel_c(i,j)=  u(i,j,1) + cx_phase
      v_rel_c(i,j)=  v(i,j,1) + cy_phase
    
      n_u(i,j)    = u_rel_c(i,j)/SQRT(u_rel_c(i,j)**2+v_rel_c(i,j)**2)
      n_v(i,j)    = v_rel_c(i,j)/SQRT(u_rel_c(i,j)**2+v_rel_c(i,j)**2)
      fd_u(i,j)   = -SIGN(1._rprec,u_rel_c(i,j))*1.2/(1+6*a_amp**2*k_wavno**2)*(a_amp*k_wavno)*abs(u(i,j,1))/dz &
                   *SQRT(u_rel_c(i,j)**2+v_rel_c(i,j)**2) &
                   *(n_u(i,j)*detadx(i,j)+n_v(i,j)*detady(i,j))*heaviside_scalar(n_u(i,j)*detadx(i,j)+n_v(i,j)*detady(i,j)) 
      fd_v(i,j)   = -SIGN(1._rprec,v_rel_c(i,j))*1.2/(1+6*a_amp**2*k_wavno**2)*(a_amp*k_wavno)*abs(v(i,j,1))/dz &
                   *SQRT(u_rel_c(i,j)**2+v_rel_c(i,j)**2) &
                   *(n_u(i,j)*detadx(i,j)+n_v(i,j)*detady(i,j))*heaviside_scalar(n_u(i,j)*detadx(i,j)+n_v(i,j)*detady(i,j))
    
    end do
    end do
  CASE DEFAULT
! WIND FOLLOWING WAVES
    do i = 1, nx
    do j = 1, ny
      eta(i,j)    =              a_amp*cos( kx_wavno*grid%x(i) + ky_wavno*grid%y(j) - omega_freq*total_time )
      detadx(i,j) =    -a_amp*kx_wavno*sin( kx_wavno*grid%x(i) + ky_wavno*grid%y(j) - omega_freq*total_time )
      detady(i,j) =    -a_amp*ky_wavno*sin( kx_wavno*grid%x(i) + ky_wavno*grid%y(j) - omega_freq*total_time )
      us_orb(i,j) =   a_amp*omega_freq*cos( kx_wavno*grid%x(i) + ky_wavno*grid%y(j) - omega_freq*total_time )
      ws_orb(i,j) =   a_amp*omega_freq*sin( kx_wavno*grid%x(i) + ky_wavno*grid%y(j) - omega_freq*total_time )
    
      u_rel(i,j)  =  u(i,j,1) - us_orb(i,j)
      v_rel(i,j)  =  v(i,j,1) - vs_orb(i,j)
      w_rel(i,j)  =  w(i,j,1) - ws_orb(i,j)
   
      u_rel_c(i,j)=  u(i,j,1) - cx_phase
      v_rel_c(i,j)=  v(i,j,1) - cy_phase
    
      n_u(i,j)    = u_rel_c(i,j)/SQRT(u_rel_c(i,j)**2+v_rel_c(i,j)**2)
      n_v(i,j)    = v_rel_c(i,j)/SQRT(u_rel_c(i,j)**2+v_rel_c(i,j)**2)
      fd_u(i,j)   = -SIGN(1._rprec,u_rel_c(i,j))*1.2/(1+6*a_amp**2*k_wavno**2)*(a_amp*k_wavno)*abs(u(i,j,1))/dz &
                   *SQRT(u_rel_c(i,j)**2+v_rel_c(i,j)**2) &
                   *(n_u(i,j)*detadx(i,j)+n_v(i,j)*detady(i,j))*heaviside_scalar(n_u(i,j)*detadx(i,j)+n_v(i,j)*detady(i,j))
    
      fd_v(i,j)   = -SIGN(1._rprec,v_rel_c(i,j))*1.2/(1+6*a_amp**2*k_wavno**2)*(a_amp*k_wavno)*abs(v(i,j,1))/dz &
                   *SQRT(u_rel_c(i,j)**2+v_rel_c(i,j)**2) &
                   *(n_u(i,j)*detadx(i,j)+n_v(i,j)*detady(i,j))*heaviside_scalar(n_u(i,j)*detadx(i,j)+n_v(i,j)*detady(i,j))
    end do
    end do


END SELECT

do i = 1,nx
        do j = 1,ny
             ! Check for NaN in the array
             !print*, "eta,i,j", eta(i,j),i,j
             if (isnan(eta(i,j))) then
                print *, "Array contains NaN"
                stop   ! Stop the program immediately
             end if
        enddo
enddo

end subroutine

!subroutine JONSWAP(EK)



!end subroutine 


end module sea_surface_drag_model
