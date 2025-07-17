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

!*******************************************************************************
module io
!*******************************************************************************
use types, only : rprec
use param, only : ld, nx, ny, nz, nz_tot, path, coord, rank, nproc, jt_total
use param, only : total_time, total_time_dim, lbz, jzmin, jzmax
use param, only : cumulative_time
use sim_param , only : w, dudz, dvdz
use sgs_param , only : Cs_opt2
use string_util
use messages
use time_average
#ifdef PPMPI
use mpi
#endif

#ifdef PPCGNS
use cgns
#ifdef PPMPI
use param, only: ierr
#endif
#endif

implicit none
save
private

public jt_total, openfiles, energy, output_loop, output_final, output_init,    &
    write_tau_wall_bot, write_tau_wall_top, height_ABL

! Where to end with nz index.
integer :: nz_end

! time averaging
type(tavg_t) :: tavg

character(:), allocatable :: fcumulative_time

contains

!*******************************************************************************
subroutine openfiles()
!*******************************************************************************
use param, only : use_cfl_dt, dt, cfl_f, checkpoint_file
implicit none
logical :: exst

! Temporary values used to read time step and CFL from file
real(rprec) :: dt_r, cfl_r

! Create file names
allocate(fcumulative_time, source = path // 'total_time.dat')
allocate(checkpoint_file , source = path // 'vel.out')

if (cumulative_time) then
    inquire (file=fcumulative_time, exist=exst)
    if (exst) then
        open (1, file=fcumulative_time)
        read(1, *) jt_total, total_time, total_time_dim, dt_r, cfl_r
        close (1)
    else
        ! assume this is the first run on cumulative time
        if ( coord == 0 ) then
            write (*, *) '--> Assuming jt_total = 0, total_time = 0.0'
        end if
        jt_total = 0
        total_time = 0._rprec
        total_time_dim = 0._rprec
    end if
end if

! Update dynamic time stepping info if required; otherwise discard.
if ( use_cfl_dt ) then
    dt = dt_r
    cfl_f = cfl_r
end if

end subroutine openfiles

!*******************************************************************************
subroutine energy (ke)
!*******************************************************************************
use types, only : rprec
use param
use sim_param, only : u, v, w
use messages
implicit none
integer :: jx, jy, jz, nan_count
real(rprec)::KE,temp_w
#ifdef PPMPI
real(rprec) :: ke_global
#endif

! Initialize variables
nan_count = 0
ke = 0._rprec

do jz = 1, nz-1
do jy = 1, ny
do jx = 1, nx
    temp_w = 0.5_rprec*(w(jx,jy,jz)+w(jx,jy,jz+1))
    ke = ke + (u(jx,jy,jz)**2+v(jx,jy,jz)**2+temp_w**2)
end do
end do
end do

! Perform spatial averaging
ke = ke*0.5_rprec/(nx*ny*(nz-1))

#ifdef PPMPI
call mpi_reduce (ke, ke_global, 1, MPI_RPREC, MPI_SUM, 0, comm, ierr)
if (rank == 0) then  ! note that it's rank here, not coord
    ke = ke_global/nproc
#endif
    open(2,file=path // 'output/check_ke.dat', status='unknown',               &
        form='formatted', position='append')
    write(2,*) total_time,ke
    close(2)
#ifdef PPMPI
end if
#endif

end subroutine energy

!*******************************************************************************
subroutine height_ABL ()
!*******************************************************************************
use types, only : rprec
use param
use string_util
use sim_param, only : u, v, w, txz, tyz
use messages
use functions, only : binary_search, interp_to_w_grid
use grid_m, only : grid
use coriolis, only : alpha
#ifdef PPSCALARS
use scalars, only : theta, scal_bot, pi_x, pi_y, pi_z, lbc_scal
#endif

implicit none
character (64) :: fname
integer :: jx, jy, jz, k
integer :: displs(nproc), rcounts(nproc)
real(rprec), allocatable, dimension(:,:,:) :: u_w_temp, v_w_temp
real(rprec), dimension(:), allocatable :: u_w_xy_avg, v_w_xy_avg, uw_xy_avg, vw_xy_avg, w_xy_avg
real(rprec), dimension(:), allocatable :: tau_Re_xz, tau_Re_yz, tau_SGS_xz, tau_SGS_yz, tau_total_xz, tau_total_yz, tau_total
real(rprec) :: h_ABL,tau_wall, tau_wall_global, tau_wall_dummy 
real(rprec), dimension(:), allocatable :: tau_total_global,tau_total_xz_global, tau_total_yz_global
real(rprec), dimension(:), allocatable :: tau_Re_xz_global, tau_Re_yz_global, tau_SGS_xz_global, tau_SGS_yz_global
real(rprec), dimension(:), allocatable :: v_w_xy_avg_global, u_w_xy_avg_global
#ifdef PPSCALARS
real(rprec), allocatable, dimension(:,:,:) :: theta_w_temp
real(rprec), dimension(:), allocatable :: theta_w_xy_avg, wtheta_w_xy_avg, wpthetap_w_xy_avg, pi_z_xy_avg, pi_z_total_xy_avg
real(rprec), dimension(:), allocatable :: theta_w_xy_avg_global, wpthetap_w_xy_avg_global, pi_z_xy_avg_global, pi_z_total_xy_avg_global
real(rprec) :: use_scal_bot
#endif
! Initialize variables
allocate ( tau_Re_xz(lbz:nz) )           ; tau_Re_xz    = 0.0_rprec
allocate ( tau_Re_yz(lbz:nz) )           ; tau_Re_yz    = 0.0_rprec
allocate ( tau_SGS_xz(lbz:nz) )          ; tau_SGS_xz   = 0.0_rprec
allocate ( tau_SGS_yz(lbz:nz) )          ; tau_SGS_yz   = 0.0_rprec
allocate ( tau_total_xz(lbz:nz) )        ; tau_total_xz = 0.0_rprec
allocate ( tau_total_yz(lbz:nz) )        ; tau_total_yz = 0.0_rprec
allocate ( tau_total(lbz:nz) )           ; tau_total    = 0.0_rprec 
allocate ( u_w_temp(nx,ny,lbz:nz) )      ; u_w_temp     = 0.0_rprec
allocate ( v_w_temp(nx,ny,lbz:nz) )      ; v_w_temp     = 0.0_rprec
allocate (  w_xy_avg(lbz:nz) )           ; w_xy_avg     = 0.0_rprec
allocate ( uw_xy_avg(lbz:nz) )           ; uw_xy_avg    = 0.0_rprec
allocate ( vw_xy_avg(lbz:nz) )           ; vw_xy_avg    = 0.0_rprec
allocate ( u_w_xy_avg(lbz:nz) )          ; u_w_xy_avg   = 0.0_rprec
allocate ( v_w_xy_avg(lbz:nz) )          ; v_w_xy_avg   = 0.0_rprec
allocate ( tau_Re_xz_global(1:nz_tot))   ; tau_Re_xz_global    = 0.0_rprec
allocate ( tau_Re_yz_global(1:nz_tot))   ; tau_Re_yz_global    = 0.0_rprec
allocate ( tau_SGS_xz_global(1:nz_tot))  ; tau_SGS_xz_global   = 0.0_rprec
allocate ( tau_SGS_yz_global(1:nz_tot))  ; tau_SGS_yz_global   = 0.0_rprec
allocate ( tau_total_global(1:nz_tot))   ; tau_total_global    = 0.0_rprec
allocate ( tau_total_xz_global(1:nz_tot)); tau_total_xz_global = 0.0_rprec
allocate ( tau_total_yz_global(1:nz_tot)); tau_total_yz_global = 0.0_rprec
allocate ( u_w_xy_avg_global(1:nz_tot))  ; u_w_xy_avg_global = 0.0_rprec
allocate ( v_w_xy_avg_global(1:nz_tot))  ; v_w_xy_avg_global = 0.0_rprec
#ifdef PPSCALARS
allocate ( theta_w_temp(nx,ny,lbz:nz) )      ; theta_w_temp         = 0.0_rprec
allocate ( theta_w_xy_avg(lbz:nz)     )      ; theta_w_xy_avg       = 0.0_rprec
allocate ( wtheta_w_xy_avg(lbz:nz)    )      ; wtheta_w_xy_avg      = 0.0_rprec
allocate ( wpthetap_w_xy_avg(lbz:nz)  )      ; wpthetap_w_xy_avg    = 0.0_rprec
allocate ( pi_z_xy_avg(lbz:nz)        )      ; pi_z_xy_avg          = 0.0_rprec
allocate ( pi_z_total_xy_avg(lbz:nz)  )      ; pi_z_total_xy_avg    = 0.0_rprec
allocate ( theta_w_xy_avg_global(1:nz_tot))    ;   theta_w_xy_avg_global    = 0.0_rprec
allocate ( wpthetap_w_xy_avg_global(1:nz_tot)) ; wpthetap_w_xy_avg_global   = 0.0_rprec
allocate ( pi_z_xy_avg_global(1:nz_tot))       ;   pi_z_xy_avg_global       = 0.0_rprec
allocate ( pi_z_total_xy_avg_global(1:nz_tot)) ;   pi_z_total_xy_avg_global = 0.0_rprec
#endif
!Interpolate u,v to w grid
u_w_temp(1:nx,1:ny,lbz:nz) = interp_to_w_grid(u(1:nx,1:ny,lbz:nz), lbz )
v_w_temp(1:nx,1:ny,lbz:nz) = interp_to_w_grid(v(1:nx,1:ny,lbz:nz), lbz )
#ifdef PPSCALARS

#ifdef PPMPI
call mpi_bcast(scal_bot,1,mpi_rprec,0,comm,ierr)
#endif

theta_w_temp(1:nx,1:ny,lbz:nz) = interp_to_w_grid(theta(1:nx,1:ny,lbz:nz), lbz )
if ( rank == 0 ) then
    !print*,"scal_bot in io", rank, total_time, scal_bot
    theta_w_temp(1:nx,1:ny,1)=scal_bot
endif

if (lbc_scal == 2) then
   use_scal_bot = 1.0
else
   use_scal_bot = 0.0
endif 
#endif

do jz = 1, nz
      u_w_xy_avg(jz)   =   u_w_xy_avg(jz) + sum( u_w_temp(1:nx,1:ny,jz) )
      v_w_xy_avg(jz)   =   v_w_xy_avg(jz) + sum( v_w_temp(1:nx,1:ny,jz) )
       uw_xy_avg(jz)   =    uw_xy_avg(jz) + sum( u_w_temp(1:nx,1:ny,jz)*w(1:nx,1:ny,jz) )
       vw_xy_avg(jz)   =    vw_xy_avg(jz) + sum( v_w_temp(1:nx,1:ny,jz)*w(1:nx,1:ny,jz) )
        w_xy_avg(jz)   =     w_xy_avg(jz) + sum(   w(1:nx,1:ny,jz) )
    tau_SGS_xz(jz)     =   tau_SGS_xz(jz) + sum( txz(1:nx,1:ny,jz) )
    tau_SGS_yz(jz)     =   tau_SGS_yz(jz) + sum( tyz(1:nx,1:ny,jz) )
#ifdef PPSCALARS
      theta_w_xy_avg(jz) =  theta_w_xy_avg(jz) + sum( theta_w_temp(1:nx,1:ny,jz) - use_scal_bot*scal_bot)
     wtheta_w_xy_avg(jz) = wtheta_w_xy_avg(jz) + sum( (theta_w_temp(1:nx,1:ny,jz) - use_scal_bot*scal_bot) *w(1:nx,1:ny,jz) )
         pi_z_xy_avg(jz) =     pi_z_xy_avg(jz) + sum( pi_z(1:nx,1:ny,jz) )
#endif
end do

! Perform spatial averaging
u_w_xy_avg = u_w_xy_avg/(nx*ny)
v_w_xy_avg = v_w_xy_avg/(nx*ny)
w_xy_avg   = w_xy_avg/(nx*ny)
uw_xy_avg  = uw_xy_avg/(nx*ny)
vw_xy_avg  = vw_xy_avg/(nx*ny)
tau_SGS_xz = tau_SGS_xz/(nx*ny)
tau_SGS_yz = tau_SGS_yz/(nx*ny)
#ifdef PPSCALARS
 theta_w_xy_avg =  theta_w_xy_avg/(nx*ny)
wtheta_w_xy_avg = wtheta_w_xy_avg/(nx*ny)
    pi_z_xy_avg =     pi_z_xy_avg/(nx*ny)
#endif
!if (coord == 0) then
!   do k=1,nz
!    print*,"coord, t01,t, wt, pi, scal_bot", coord, total_time, theta_w_xy_avg(k), wtheta_w_xy_avg(k), pi_z_xy_avg(k), scal_bot
!   end do
!endif
!if (coord == 1) then
!   do k=1,nz
!    print*,"coord, t01,t, wt, pi, scal_bot", coord, total_time, theta_w_xy_avg(k), wtheta_w_xy_avg(k), pi_z_xy_avg(k), scal_bot
!   end do
!endif
!Compute Re stress, SGS stress and total stress
tau_Re_xz    = uw_xy_avg-u_w_xy_avg*w_xy_avg
tau_Re_yz    = vw_xy_avg-v_w_xy_avg*w_xy_avg
tau_total_xz = -tau_Re_xz-tau_SGS_xz
tau_total_yz = -tau_Re_yz-tau_SGS_yz
tau_total    = sqrt(tau_total_xz**2+tau_total_yz**2)
#ifdef PPSCALARS
wpthetap_w_xy_avg =   wtheta_w_xy_avg - theta_w_xy_avg*w_xy_avg
pi_z_total_xy_avg = - wpthetap_w_xy_avg - pi_z_xy_avg
#endif 
!Get the wall stress
if (rank == 0) then
tau_wall = tau_total(1)
end if

#ifdef PPMPI
!tau_wall locally is 0 in procs with rank>0. Therefore, mpi_allreduce with mpi_sum communicates 
!tau_wall from rank 0 to all other procs and stores the value in tau_wall_dummy
call mpi_allreduce (tau_wall, tau_wall_dummy, 1, MPI_RPREC, MPI_SUM,  comm, ierr)
tau_wall_global=tau_wall_dummy
!Gather total stress and its components into global array using MPI_GATHERV to proc 0 for visualization purposes
do k=1, nproc
   displs(k)=(k-1)*(nz-1)
  rcounts(k)=nz-1
end do
call MPI_Gatherv (   tau_total(1:nz-1), nz-1, MPI_RPREC, tau_total_global   , rcounts, displs, MPI_RPREC, 0, comm, ierr)
call MPI_Gatherv (tau_total_xz(1:nz-1), nz-1, MPI_RPREC, tau_total_xz_global, rcounts, displs, MPI_RPREC, 0, comm, ierr)
call MPI_Gatherv (tau_total_yz(1:nz-1), nz-1, MPI_RPREC, tau_total_yz_global, rcounts, displs, MPI_RPREC, 0, comm, ierr)
call MPI_Gatherv (   tau_Re_xz(1:nz-1), nz-1, MPI_RPREC, tau_Re_xz_global   , rcounts, displs, MPI_RPREC, 0, comm, ierr)
call MPI_Gatherv (   tau_Re_yz(1:nz-1), nz-1, MPI_RPREC, tau_Re_yz_global   , rcounts, displs, MPI_RPREC, 0, comm, ierr)
call MPI_Gatherv (  tau_SGS_xz(1:nz-1), nz-1, MPI_RPREC, tau_SGS_xz_global  , rcounts, displs, MPI_RPREC, 0, comm, ierr)
call MPI_Gatherv (  tau_SGS_yz(1:nz-1), nz-1, MPI_RPREC, tau_SGS_yz_global  , rcounts, displs, MPI_RPREC, 0, comm, ierr)
call MPI_Gatherv (  u_w_xy_avg(1:nz-1), nz-1, MPI_RPREC, u_w_xy_avg_global  , rcounts, displs, MPI_RPREC, 0, comm, ierr)
call MPI_Gatherv (  v_w_xy_avg(1:nz-1), nz-1, MPI_RPREC, v_w_xy_avg_global  , rcounts, displs, MPI_RPREC, 0, comm, ierr)

#ifdef PPSCALARS
call MPI_Gatherv (     theta_w_xy_avg(1:nz-1), nz-1, MPI_RPREC,    theta_w_xy_avg_global  , rcounts, displs, MPI_RPREC, 0, comm, ierr)
call MPI_Gatherv (  wpthetap_w_xy_avg(1:nz-1), nz-1, MPI_RPREC, wpthetap_w_xy_avg_global  , rcounts, displs, MPI_RPREC, 0, comm, ierr)
call MPI_Gatherv (        pi_z_xy_avg(1:nz-1), nz-1, MPI_RPREC,       pi_z_xy_avg_global  , rcounts, displs, MPI_RPREC, 0, comm, ierr)
call MPI_Gatherv (  pi_z_total_xy_avg(1:nz-1), nz-1, MPI_RPREC, pi_z_total_xy_avg_global  , rcounts, displs, MPI_RPREC, 0, comm, ierr)
#endif
!if (rank ==0) then
!do jz=1,nz_tot
!print*, "t2, tau_tot_global, jz", total_time, tau_total_global(jz), jz
!end do
!end if

!if (rank == 0) then
!do jz=1,nz_tot
!   print*,"t2,t ,pi, wt, pit ,jz",total_time, theta_w_xy_avg_global(jz), pi_z_xy_avg_global(jz), wpthetap_w_xy_avg_global(jz), pi_z_total_xy_avg_global(jz) ,jz
!end do
!end if
!if (rank == 1) then
!do k=1,nz-1
!   print*,"t3,t, pi,wt,pit,jz", total_time, theta_w_xy_avg(k), pi_z_xy_avg(k), wpthetap_w_xy_avg(k),pi_z_total_xy_avg(k), k 
!end do
!end if
#endif


!Compute height of ABL from total stress using the criteria
!tau_total(h_ABL)=0.05*tau_wall_global
if (rank==0) then
  jz=1
  DO WHILE ( tau_total_global(jz) .gt. 0.05_rprec*tau_wall_global .and. jz<nz_tot) 
        h_ABL=(jz - 0.5_rprec) * dz - dz/2._rprec
          jz = jz + 1
  END DO

  !Save ABL height to file height_ABL.dat
  open(2,file=path // 'output/height_ABL.dat', status='unknown',               &
        form='formatted', position='append')
  write(2,*) total_time,h_ABL, alpha
  close(2)

  !Save x-y averaged instantaneous stress to file 
  if (tau_xy_avg_calc) then
    if ((jt_total >= tau_xy_avg_nstart).and.(jt_total <= tau_xy_avg_nend)) then
        if ( mod(jt_total-tau_xy_avg_nstart,tau_xy_avg_nskip)==0 ) then
           call string_splice(fname, path //'output/tau_xy.', jt_total)
           ! Write binary Output
           call string_concat(fname, '.bin')
           open(unit=13, file=fname, form='unformatted', convert=write_endian,        &
           access='direct', recl=nz_tot*rprec)
           write(13,rec=1) tau_total_global(1:nz_tot)
           write(13,rec=2) tau_total_xz_global(1:nz_tot)
           write(13,rec=3) tau_total_yz_global(1:nz_tot)
           write(13,rec=4) tau_Re_xz_global(1:nz_tot)
           write(13,rec=5) tau_Re_yz_global(1:nz_tot)
           write(13,rec=6) tau_SGS_xz_global(1:nz_tot)
           write(13,rec=7) tau_SGS_yz_global(1:nz_tot)
           write(13,rec=8) u_w_xy_avg_global(1:nz_tot)
           write(13,rec=9) v_w_xy_avg_global(1:nz_tot)
#ifdef PPSCALARS
           write(13,rec=10) theta_w_xy_avg_global(1:nz_tot)
           write(13,rec=11) wpthetap_w_xy_avg_global(1:nz_tot)
           write(13,rec=12) pi_z_xy_avg_global(1:nz_tot)
           write(13,rec=13) pi_z_total_xy_avg_global(1:nz_tot)
#endif
           close(13)
        end if
    end if
  end if
end if
end subroutine height_ABL
!*******************************************************************************
subroutine write_tau_wall_bot()
!*******************************************************************************
use types ,only: rprec
use param ,only: jt_total, total_time, total_time_dim, dt, dt_dim, wbase
use param ,only: L_x, z_i, u_star
use functions ,only: get_tau_wall_bot
implicit none

real(rprec) :: turnovers

turnovers = total_time_dim / (L_x * z_i / u_star)

open(2,file=path // 'output/tau_wall_bot.dat', status='unknown',               &
    form='formatted', position='append')

!! one time header output
if (jt_total==wbase) write(2,*)                                                &
    'jt_total, total_time, total_time_dim, turnovers, dt, dt_dim, 1.0, tau_wall'

!! continual time-related output
write(2,*) jt_total, total_time, total_time_dim, turnovers, dt, dt_dim,        &
    1.0, get_tau_wall_bot()
close(2)

end subroutine write_tau_wall_bot

!*******************************************************************************
subroutine write_tau_wall_top()
!*******************************************************************************
use types, only : rprec
use param, only : jt_total, total_time, total_time_dim, dt, dt_dim, wbase
use param, only : L_x, z_i, u_star
use functions, only : get_tau_wall_top
implicit none

real(rprec) :: turnovers

turnovers = total_time_dim / (L_x * z_i / u_star)

open(20,file=path // 'output/tau_wall_top.dat', status='unknown',               &
    form='formatted', position='append')

! one time header output
if (jt_total==wbase) write(2,*)                                                &
    'jt_total, total_time, total_time_dim, turnovers, dt, dt_dim, 1.0, tau_wall'

! continual time-related output
write(20,*) jt_total, total_time, total_time_dim, turnovers, dt, dt_dim,        &
    1.0, get_tau_wall_top()
close(20)

end subroutine write_tau_wall_top

#ifdef PPCGNS
#ifdef PPMPI
!*******************************************************************************
subroutine write_parallel_cgns (file_name, nx, ny, nz, nz_tot, start_n_in,     &
    end_n_in, xin, yin, zin, num_fields, fieldNames, input )
!*******************************************************************************
implicit none

integer, intent(in) :: nx, ny, nz, nz_tot, num_fields
! Name of file to be written
character(*), intent(in) :: file_name
! Name of fields we are writing
character(*), intent(in), dimension(:) :: fieldNames
! Data to be written
real(rprec), intent(in), dimension(:) :: input
! Coordinates to write
real(rprec), intent(in), dimension(:) :: xin, yin, zin
! Where the total node counter starts nodes
integer, intent(in) :: start_n_in(3)
! Where the total node counter ends nodes
integer, intent(in) :: end_n_in(3)

integer :: fn=1        ! CGNS file index number
integer :: ier         ! CGNS error status
integer :: base=1      ! base number
integer :: zone=1      ! zone number
integer :: nnodes      ! Number of nodes in this processor
integer :: sol =1      ! solution number
integer :: field       ! section number
integer(cgsize_t) :: sizes(3,3)  ! Sizes

! Convert input to right data type
integer(cgsize_t) :: start_n(3)  ! Where the total node counter starts nodes
integer(cgsize_t) :: end_n(3)  ! Where the total node counter ends nodes

! Building the lcoal mesh
integer :: i,j,k
real(rprec), dimension(nx,ny,nz) :: xyz

!  ! Set the parallel communicator
!  call cgp_mpi_comm_f(cgnsParallelComm, ierr)

! Convert types such that CGNS libraries can handle the input
start_n(1) = int(start_n_in(1), cgsize_t)
start_n(2) = int(start_n_in(2), cgsize_t)
start_n(3) = int(start_n_in(3), cgsize_t)
end_n(1) = int(end_n_in(1), cgsize_t)
end_n(2) = int(end_n_in(2), cgsize_t)
end_n(3) = int(end_n_in(3), cgsize_t)

! The total number of nodes in this processor
nnodes = nx*ny*nz

! Sizes, used to create zone
sizes(:,1) = (/int(nx, cgsize_t),int(ny, cgsize_t),int(nz_tot, cgsize_t)/)
sizes(:,2) = (/int(nx-1, cgsize_t),int(ny-1, cgsize_t),int(nz_tot-1, cgsize_t)/)
sizes(:,3) = (/int(0, cgsize_t) , int(0, cgsize_t), int(0, cgsize_t)/)

! Open CGNS file
call cgp_open_f(file_name, CG_MODE_WRITE, fn, ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

! Write base
call cg_base_write_f(fn, 'Base', 3, 3, base, ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

! Write zone
call cg_zone_write_f(fn, base, 'Zone', sizes, Structured, zone, ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

! Write print info to screen
if (coord .eq. 0) then
    write(*,*) 'Writing, ', file_name
end if

! Create data nodes for coordinates
call cgp_coord_write_f(fn, base, zone, RealDouble, 'CoordinateX', nnodes, ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

call cgp_coord_write_f(fn, base, zone, RealDouble, 'CoordinateY', nnodes, ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

call cgp_coord_write_f(fn, base, zone, RealDouble, 'CoordinateZ', nnodes, ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

! Write the coordinate data in parallel to the queue
!  call cgp_queue_set_f(1, ier)
!  if (ier .ne. CG_OK) call cgp_error_exit_f

! This is done for the 3 dimensions x,y and z
! It writes the coordinates
! Create grid points
do k = 1, nz
do j = 1, ny
do i = 1, nx
    xyz(i,j,k) = xin(i)
end do
end do
end do

call cgp_coord_write_data_f(fn, base, zone, 1,                                 &
    start_n, end_n, xyz(1:nx,1:ny,1:nz), ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

! Write out the queued coordinate data
!  call cgp_queue_flush_f(ier)
!  if (ier .ne. CG_OK) call cgp_error_exit_f
!  call cgp_queue_set_f(0, ier)

! Write the coordinate data in parallel to the queue
!  call cgp_queue_set_f(1, ier)
!  if (ier .ne. CG_OK) call cgp_error_exit_f

do k = 1, nz
do j = 1, ny
do i = 1, nx
    xyz(i,j,k) = yin(j)
end do
end do
end do
call cgp_coord_write_data_f(fn, base, zone, 2,   &
    start_n, end_n, xyz(1:nx,1:ny,1:nz), ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

! Write out the queued coordinate data
!  call cgp_queue_flush_f(ier)
!  if (ier .ne. CG_OK) call cgp_error_exit_f
!  call cgp_queue_set_f(0, ier)

! Write the coordinate data in parallel to the queue
!  call cgp_queue_set_f(1, ier)
!  if (ier .ne. CG_OK) call cgp_error_exit_f

do k = 1, nz
do j = 1, ny
do i = 1, nx
    xyz(i,j,k) = zin(k)
end do
end do
end do
call cgp_coord_write_data_f(fn, base, zone, 3,   &
                            start_n, end_n, xyz(1:nx,1:ny,1:nz), ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

! Write out the queued coordinate data
!  call cgp_queue_flush_f(ier)
!  if (ier .ne. CG_OK) call cgp_error_exit_f
!  call cgp_queue_set_f(0, ier)

! Create a centered solution
call cg_sol_write_f(fn, base, zone, 'Solution', Vertex, sol, ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

! Write the solution
do i=1,num_fields
    call cgp_field_write_f(fn, base, zone, sol, RealDouble, fieldNames(i),     &
        field, ier)
    if (ier .ne. CG_OK) call cgp_error_exit_f

    call cgp_field_write_data_f(fn, base, zone, sol, field, start_n, end_n,    &
        input((i-1)*nnodes+1:(i)*nnodes), ier)
    if (ier .ne. CG_OK) call cgp_error_exit_f

end do

! Close the file
call cgp_close_f(fn, ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

end subroutine write_parallel_cgns

!*******************************************************************************
subroutine write_null_cgns (file_name, nx, ny, nz, nz_tot, start_n_in,         &
    end_n_in, xin, yin, zin, num_fields, fieldNames )
!*******************************************************************************
implicit none

integer, intent(in) :: nx, ny, nz, nz_tot, num_fields
! Name of file to be written
character(*), intent(in) :: file_name
! Name of fields we are writing
character(*), intent(in), dimension(:) :: fieldNames
! Coordinates to write
real(rprec), intent(in), dimension(:) :: xin, yin, zin
! Where the total node counter starts nodes
integer, intent(in) :: start_n_in(3)
! Where the total node counter ends nodes
integer, intent(in) :: end_n_in(3)

integer :: fn=1        ! CGNS file index number
integer :: ier         ! CGNS error status
integer :: base=1      ! base number
integer :: zone=1      ! zone number
integer :: nnodes      ! Number of nodes in this processor
integer :: sol =1      ! solution number
integer :: field       ! section number
integer(cgsize_t) :: sizes(3,3)  ! Sizes

! Convert input to right data type
integer(cgsize_t) :: start_n(3)  ! Where the total node counter starts nodes
integer(cgsize_t) :: end_n(3)  ! Where the total node counter ends nodes

! Building the lcoal mesh
integer :: i,j,k
real(rprec), dimension(nx,ny,nz) :: xyz

!  ! Set the parallel communicator
!  call cgp_mpi_comm_f(cgnsParallelComm, ierr)

! Convert types such that CGNS libraries can handle the input
start_n(1) = int(start_n_in(1), cgsize_t)
start_n(2) = int(start_n_in(2), cgsize_t)
start_n(3) = int(start_n_in(3), cgsize_t)
end_n(1) = int(end_n_in(1), cgsize_t)
end_n(2) = int(end_n_in(2), cgsize_t)
end_n(3) = int(end_n_in(3), cgsize_t)

! The total number of nodes in this processor
nnodes = nx*ny*nz

! Sizes, used to create zone
sizes(:,1) = (/int(nx, cgsize_t),int(ny, cgsize_t),int(nz_tot, cgsize_t)/)
sizes(:,2) = (/int(nx-1, cgsize_t),int(ny-1, cgsize_t),int(nz_tot-1, cgsize_t)/)
sizes(:,3) = (/int(0, cgsize_t) , int(0, cgsize_t), int(0, cgsize_t)/)

! Open CGNS file
call cgp_open_f(file_name, CG_MODE_WRITE, fn, ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

! Write base
call cg_base_write_f(fn, 'Base', 3, 3, base, ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

! Write zone
call cg_zone_write_f(fn, base, 'Zone', sizes, Structured, zone, ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

! Write print info to screen
if (coord .eq. 0) then
    write(*,*) 'Writing, ', file_name
end if

! Create data nodes for coordinates
call cgp_coord_write_f(fn, base, zone, RealDouble, 'CoordinateX', nnodes, ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

call cgp_coord_write_f(fn, base, zone, RealDouble, 'CoordinateY', nnodes, ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

call cgp_coord_write_f(fn, base, zone, RealDouble, 'CoordinateZ', nnodes, ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

! This is done for the 3 dimensions x,y and z
! It writes the coordinates
! Create grid points
do k = 1, nz
do j = 1, ny
do i = 1, nx
    xyz(i,j,k) = xin(i)
end do
end do
end do
write(*,*) "HERE 0.8"

call cgp_coord_write_data_f(fn, base, zone, 1, start_n, end_n, %VAL(0), ier)
write(*,*) "HERE 0.85"
if (ier .ne. CG_OK) call cgp_error_exit_f
write(*,*) "HERE 0.9"

! Write out the queued coordinate data
!  call cgp_queue_flush_f(ier)
!  if (ier .ne. CG_OK) call cgp_error_exit_f
!  call cgp_queue_set_f(0, ier)

! Write the coordinate data in parallel to the queue
!  call cgp_queue_set_f(1, ier)
!  if (ier .ne. CG_OK) call cgp_error_exit_f

do k = 1, nz
do j = 1, ny
do i = 1, nx
    xyz(i,j,k) = yin(j)
end do
end do
end do
call cgp_coord_write_data_f(fn, base, zone, 2, start_n, end_n, %VAL(0), ier)
if (ier .ne. CG_OK) call cgp_error_exit_f
write(*,*) "HERE 1.0"

! Write out the queued coordinate data
!  call cgp_queue_flush_f(ier)
!  if (ier .ne. CG_OK) call cgp_error_exit_f
!  call cgp_queue_set_f(0, ier)

! Write the coordinate data in parallel to the queue
!  call cgp_queue_set_f(1, ier)
!  if (ier .ne. CG_OK) call cgp_error_exit_f

do k = 1, nz
do j = 1, ny
do i = 1, nx
    xyz(i,j,k) = zin(k)
end do
end do
end do
write(*,*) "HERE 1.1"

call cgp_coord_write_data_f(fn, base, zone, 3, start_n, end_n, %VAL(0), ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

! Create a centered solution
call cg_sol_write_f(fn, base, zone, 'Solution', Vertex, sol, ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

! Write the solution
do i = 1, num_fields
    call cgp_field_write_f(fn, base, zone, sol, RealDouble, fieldNames(i),     &
                           field, ier)
    if (ier .ne. CG_OK) call cgp_error_exit_f

    call cgp_field_write_data_f(fn, base, zone, sol, field, start_n, end_n,    &
                                %VAL(0), ier)
    if (ier .ne. CG_OK) call cgp_error_exit_f

end do

! Close the file
call cgp_close_f(fn, ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

write(*,*) "end of write_null_cgns"

end subroutine write_null_cgns
#endif
#endif

!*******************************************************************************
subroutine output_loop()
!*******************************************************************************
!
!  This subroutine is called every time step and acts as a driver for
!  computing statistics and outputing instantaneous data. No actual
!  calculations are performed here.
!
use param, only : jt_total, dt
use param, only : checkpoint_data, checkpoint_nskip
use param, only : tavg_calc, tavg_nstart, tavg_nend, tavg_nskip
use param, only : point_calc, point_nstart, point_nend, point_nskip
use param, only : domain_xz_calc, domain_xz_nstart, domain_xz_nend, domain_xz_nskip
use param, only : domain_calc, domain_nstart, domain_nend, domain_nskip
use param, only : xplane_calc, xplane_nstart, xplane_nend, xplane_nskip
use param, only : yplane_calc, yplane_nstart, yplane_nend, yplane_nskip
use param, only : zplane_calc, zplane_nstart, zplane_nend, zplane_nskip
implicit none

! Determine if we are to checkpoint intermediate times
if( checkpoint_data ) then
    ! Now check if data should be checkpointed this time step
    if ( modulo (jt_total, checkpoint_nskip) == 0) call checkpoint()
end if

!  Determine if time summations are to be calculated
if (tavg_calc) then
    ! Are we between the start and stop timesteps?
    if ((jt_total >= tavg_nstart).and.(jt_total <= tavg_nend)) then
        ! Every timestep (between nstart and nend), add to tavg%dt
        tavg%dt = tavg%dt + dt

        ! Are we at the beginning or a multiple of nstart?
        if ( mod(jt_total-tavg_nstart,tavg_nskip)==0 ) then
            ! Check if we have initialized tavg
            if (.not.tavg%initialized) then
                if (coord == 0) then
                    write(*,*) '-------------------------------'
                    write(*,"(1a,i9,1a,i9)")                                   &
                        'Starting running time summation from ',               &
                        tavg_nstart, ' to ', tavg_nend
                    write(*,*) '-------------------------------'
                end if

                call tavg%init()
            else
                call tavg%compute()
            end if
        end if
    end if
end if

!  Determine if instantaneous point velocities are to be recorded
if(point_calc) then
    if (jt_total >= point_nstart .and. jt_total <= point_nend .and.            &
        ( mod(jt_total-point_nstart,point_nskip)==0) ) then
        if (jt_total == point_nstart) then
            if (coord == 0) then
                write(*,*) '-------------------------------'
                write(*,"(1a,i9,1a,i9)")                                       &
                    'Writing instantaneous point velocities from ',            &
                    point_nstart, ' to ', point_nend
                write(*,"(1a,i9)") 'Iteration skip:', point_nskip
                write(*,*) '-------------------------------'
            end if
        end if
        call inst_write(1)
    end if
end if

!  Determine if instantaneous domain velocities are to be recorded
if(domain_calc) then
    if (jt_total >= domain_nstart .and. jt_total <= domain_nend .and.          &
        ( mod(jt_total-domain_nstart,domain_nskip)==0) ) then
        if (jt_total == domain_nstart) then
            if (coord == 0) then
                write(*,*) '-------------------------------'
                write(*,"(1a,i9,1a,i9)")                                       &
                    'Writing instantaneous domain velocities from ',           &
                    domain_nstart, ' to ', domain_nend
                write(*,"(1a,i9)") 'Iteration skip:', domain_nskip
                write(*,*) '-------------------------------'
            end if

        end if
        call inst_write(2)
    end if
end if

!  Determine if instantaneous x-plane velocities are to be recorded
if(xplane_calc) then
    if (jt_total >= xplane_nstart .and. jt_total <= xplane_nend .and.          &
        ( mod(jt_total-xplane_nstart,xplane_nskip)==0) ) then
    if (jt_total == xplane_nstart) then
        if (coord == 0) then
            write(*,*) '-------------------------------'
            write(*,"(1a,i9,1a,i9)")                                           &
                'Writing instantaneous x-plane velocities from ',              &
                xplane_nstart, ' to ', xplane_nend
            write(*,"(1a,i9)") 'Iteration skip:', xplane_nskip
            write(*,*) '-------------------------------'
            end if
        end if

        call inst_write(3)
    end if
end if

!  Determine if instantaneous y-plane velocities are to be recorded
if(yplane_calc) then
    if (jt_total >= yplane_nstart .and. jt_total <= yplane_nend .and.          &
        ( mod(jt_total-yplane_nstart,yplane_nskip)==0) ) then
        if (jt_total == yplane_nstart) then
            if (coord == 0) then
                write(*,*) '-------------------------------'
                write(*,"(1a,i9,1a,i9)")                                       &
                    'Writing instantaneous y-plane velocities from ',          &
                    yplane_nstart, ' to ', yplane_nend
                write(*,"(1a,i9)") 'Iteration skip:', yplane_nskip
                write(*,*) '-------------------------------'
            end if
        end if

        call inst_write(4)
    end if
end if

!  Determine if instantaneous z-plane velocities are to be recorded
if(zplane_calc) then
    if (jt_total >= zplane_nstart .and. jt_total <= zplane_nend .and.          &
        ( mod(jt_total-zplane_nstart,zplane_nskip)==0) ) then
        if (jt_total == zplane_nstart) then
            if (coord == 0) then
                write(*,*) '-------------------------------'
                write(*,"(1a,i9,1a,i9)")                                       &
                    'Writing instantaneous z-plane velocities from ',          &
                    zplane_nstart, ' to ', zplane_nend
                write(*,"(1a,i9)") 'Iteration skip:', zplane_nskip
                write(*,*) '-------------------------------'
            end if
        end if

        call inst_write(5)
    end if
end if

!  Determine if instantaneous spanwise averaged domain velocities are to be recorded
if(domain_xz_calc) then
    if (jt_total >= domain_xz_nstart .and. jt_total <= domain_xz_nend .and.          &
        ( mod(jt_total-domain_xz_nstart,domain_xz_nskip)==0) ) then
        if (jt_total == domain_xz_nstart) then
            if (coord == 0) then
                write(*,*) '-------------------------------'
                write(*,"(1a,i9,1a,i9)")                                       &
                    'Writing instantaneous y-averaged domain velocities from ',           &
                    domain_xz_nstart, ' to ', domain_xz_nend
                write(*,"(1a,i9)") 'Iteration skip:', domain_xz_nskip
                write(*,*) '-------------------------------'
            end if

        end if
        call inst_write(6)
    end if
end if
end subroutine output_loop

!*******************************************************************************
subroutine inst_write(itype)
!*******************************************************************************
!
! This subroutine is used to write all of the instantaneous data from
! lesgo to file. The types of data written are:
!
!   points   : itype=1
!   domain   : itype=2
!   x-planes : itype=3
!   y-planes : itype=4
!   z-planes : itype=5
!   xz-planes: itype=6
!
! For the points and planar data, this subroutine writes using the
! locations specfied from the param module.
! If additional instantenous values are
! desired to be written, they should be done so using this subroutine.
!
use functions, only : linear_interp, trilinear_interp, interp_to_uv_grid
use param, only : point_nloc, point_loc
use param, only : xplane_nloc, xplane_loc
use param, only : yplane_nloc, yplane_loc
use param, only : zplane_nloc, zplane_loc
use param, only : dx, dy
use param, only : write_endian
use param, only : use_exp_decay, use_sea_drag_model
use sgs_param, only : Cs_opt2, F_LM, F_MM, F_QN, F_NN, Nu_t
use grid_m
use sim_param, only : u, v, w, p
use sim_param, only : txx,tyy,tzz,txy,txz,tyz !GN
use sim_param, only : fxa, fya, fza !AA
use sim_param, only : dwdy, dwdx, dvdx, dudy
use functions, only : interp_to_w_grid
use sea_surface_drag_model, only : fd_u, fd_v
use stat_defs, only : xplane, yplane
#ifdef PPMPI
use stat_defs, only : zplane, point
use param, only : ny, nz, dz
#endif
#ifdef PPLVLSET
use level_set_base, only : phi
use sim_param, only : fx, fy, fz, fxa, fya, fza
use sgs_param, only : delta
#endif
#ifdef PPSCALARS
use scalars, only : theta, pi_x, pi_y, pi_z
#endif

implicit none

integer, intent(in) :: itype
character (64) :: fname
integer :: n, i, j, k
real(rprec), allocatable, dimension(:,:,:) :: ui, vi, wi,w_uv
real(rprec), pointer, dimension(:) :: x, y, z, zw
! Vorticity
real(rprec), dimension (:,:,:), allocatable :: vortx, vorty, vortz

! Pressure
real(rprec), dimension(:,:,:), allocatable :: pres_real
#ifndef PPCGNS
character(64) :: bin_ext

#ifdef PPLVLSET
real(rprec), allocatable, dimension(:,:,:) :: fx_tot, fy_tot, fz_tot
#endif

#ifdef PPMPI
call string_splice(bin_ext, '.c', coord, '.bin')
#else
bin_ext = '.bin'
#endif
#endif

! Nullify pointers
nullify(x,y,z,zw)

! Set grid pointers
x => grid % x
y => grid % y
z => grid % z
zw => grid % zw

!  Allocate space for the interpolated w values
allocate(w_uv(nx,ny,lbz:nz))

!  Make sure w has been interpolated to uv-grid
w_uv = interp_to_uv_grid(w(1:nx,1:ny,lbz:nz), lbz)

!  Instantaneous velocity sampled at point
if(itype==1) then
    do n = 1, point_nloc
        ! Common file name for all output types
        call string_splice(fname, path // 'output/vel.x-', point_loc(n)%xyz(1),&
            '.y-', point_loc(n)%xyz(2), '.z-', point_loc(n)%xyz(3), '.dat')

#ifdef PPMPI
        if(point(n) % coord == coord) then
#endif
            open(unit=13, position="append", file=fname)
            write(13,*) total_time,                                            &
            trilinear_interp(u(1:nx,1:ny,lbz:nz), lbz, point_loc(n)%xyz),      &
            trilinear_interp(v(1:nx,1:ny,lbz:nz), lbz, point_loc(n)%xyz),      &
            trilinear_interp(w_uv(1:nx,1:ny,lbz:nz), lbz, point_loc(n)%xyz)
            close(13)
#ifdef PPMPI
        end if
#endif
    end do

!  Instantaneous write for entire domain
elseif(itype==2) then
    ! Common file name for all output types
    call string_splice(fname, path //'output/vel.', jt_total)

#if defined(PPCGNS) && defined(PPMPI)
    ! Write CGNS Output
    call string_concat(fname, '.cgns')
    call write_parallel_cgns(fname, nx, ny, nz - nz_end, nz_tot,               &
        (/ 1, 1,   (nz-1)*coord + 1 /),                                        &
        (/ nx, ny, (nz-1)*(coord+1) + 1 - nz_end /),                           &
        x(1:nx) , y(1:ny) , z(1:(nz-nz_end) ),                                 &
        3, (/ 'VelocityX', 'VelocityY', 'VelocityZ' /),                        &
        (/ u(1:nx,1:ny,1:(nz-nz_end)), v(1:nx,1:ny,1:(nz-nz_end)),             &
         w_uv(1:nx,1:ny,1:(nz-nz_end)) /) )
#else
    ! Write binary Output
    call string_concat(fname, bin_ext)
    open(unit=13, file=fname, form='unformatted', convert=write_endian,        &
        access='direct', recl=nx*ny*nz*rprec)
    write(13,rec=1) u(:nx,:ny,1:nz)
    write(13,rec=2) v(:nx,:ny,1:nz)
    write(13,rec=3) w_uv(:nx,:ny,1:nz)
    close(13)
#endif

    ! Common file name for all output types
    call string_splice(fname, path //'output/sgs.', jt_total)

#if defined(PPCGNS) && defined(PPMPI)
    ! Write CGNS Output
    call string_concat(fname, '.cgns')
    call write_parallel_cgns(fname, nx, ny, nz - nz_end, nz_tot,               &
        (/ 1, 1,   (nz-1)*coord + 1 /),                                        &
        (/ nx, ny, (nz-1)*(coord+1) + 1 - nz_end /),                           &
        x(1:nx) , y(1:ny) , z(1:(nz-nz_end) ),                                 &
        5, (/ 'Cs_opt2', 'LM', 'MM','QN','NN','Nu_t' /),                        &
        (/ Cs_opt2(1:nx,1:ny,1:(nz-nz_end)), F_LM(1:nx,1:ny,1:(nz-nz_end)),             &
         F_MM(1:nx,1:ny,1:(nz-nz_end)),F_QN(1:nx,1:ny,1:(nz-nz_end)), F_NN(1:nx,1:ny,1:(nz-nz_end)), Nu_t(1:nx,1:ny,1:(nz-nz_end)) /) )
#else
    ! Write binary Output
    call string_concat(fname, bin_ext)
    open(unit=13, file=fname, form='unformatted', convert=write_endian,        &
        access='direct', recl=nx*ny*nz*rprec)
    write(13,rec=1) Cs_opt2(:nx,:ny,1:nz)
    write(13,rec=2) F_LM(:nx,:ny,1:nz)
    write(13,rec=3) F_MM(:nx,:ny,1:nz)
    write(13,rec=4) F_QN(:nx,:ny,1:nz)
    write(13,rec=5) F_NN(:nx,:ny,1:nz)
    write(13,rec=6) Nu_t(:nx,:ny,1:nz)
    close(13)
#endif

    ! Common file name for all output types !GN
    call string_splice(fname, path //'output/tau_sgs.', jt_total)

#if defined(PPCGNS) && defined(PPMPI)
    ! Write CGNS Output
    call string_concat(fname, '.cgns')
    call write_parallel_cgns(fname, nx, ny, nz - nz_end, nz_tot,               &
        (/ 1, 1,   (nz-1)*coord + 1 /),                                        &
        (/ nx, ny, (nz-1)*(coord+1) + 1 - nz_end /),                           &
        x(1:nx) , y(1:ny) , z(1:(nz-nz_end) ),                                 &
        5, (/ 'txx', 'tyy', 'tzz','txy','txz','tyz' /),                        &
        (/ txx(1:nx,1:ny,1:(nz-nz_end)), tyy(1:nx,1:ny,1:(nz-nz_end)),             &
         tzz(1:nx,1:ny,1:(nz-nz_end)),txy(1:nx,1:ny,1:(nz-nz_end)), txz(1:nx,1:ny,1:(nz-nz_end)), tyz(1:nx,1:ny,1:(nz-nz_end)) /) )
#else
    ! Write binary Output
    call string_concat(fname, bin_ext)
    open(unit=13, file=fname, form='unformatted', convert=write_endian,        &
        access='direct', recl=nx*ny*nz*rprec)
    write(13,rec=1) txx(:nx,:ny,1:nz)
    write(13,rec=2) tyy(:nx,:ny,1:nz)
    write(13,rec=3) tzz(:nx,:ny,1:nz)
    write(13,rec=4) txy(:nx,:ny,1:nz)
    write(13,rec=5) txz(:nx,:ny,1:nz)
    write(13,rec=6) tyz(:nx,:ny,1:nz)
    close(13)
#endif

    ! Common file name for all output types !AA
    call string_splice(fname, path //'output/wavestress.', jt_total)

#if defined(PPCGNS) && defined(PPMPI)
    ! Write CGNS Output
    call string_concat(fname, '.cgns')
    call write_parallel_cgns(fname, nx, ny, nz - nz_end, nz_tot,               &
        (/ 1, 1,   (nz-1)*coord + 1 /),                                        &
        (/ nx, ny, (nz-1)*(coord+1) + 1 - nz_end /),                           &
        x(1:nx) , y(1:ny) , z(1:(nz-nz_end) ),                                 &
        5, (/ 'txx', 'tyy', 'tzz','txy','txz','tyz' /),                        &
        (/ fxa(1:nx,1:ny,1:(nz-nz_end)), fya(1:nx,1:ny,1:(nz-nz_end)),             &
         fza(1:nx,1:ny,1:(nz-nz_end)) /) )
#else
    if (use_sea_drag_model .and. use_exp_decay) then
            fxa(:nx,:ny,1) = fd_u(:nx,:ny)*dz 
            fya(:nx,:ny,1) = fd_v(:nx,:ny)*dz
    ! Write binary Output
    call string_concat(fname, bin_ext)
    open(unit=13, file=fname, form='unformatted', convert=write_endian,        &
        access='direct', recl=nx*ny*nz*rprec)
    write(13,rec=1) fxa(:nx,:ny,1:nz)
    write(13,rec=2) fya(:nx,:ny,1:nz)
    write(13,rec=3) fza(:nx,:ny,1:nz)
    close(13)
    endif
#endif

    ! Compute vorticity
    allocate(vortx(nx,ny,lbz:nz), vorty(nx,ny,lbz:nz), vortz(nx,ny,lbz:nz))
    vortx(1:nx,1:ny,lbz:nz) = 0._rprec
    vorty(1:nx,1:ny,lbz:nz) = 0._rprec
    vortz(1:nx,1:ny,lbz:nz) = 0._rprec

    ! Use vorticityx as an intermediate step for performing uv-w interpolation
    ! Vorticity is written in w grid
    vortx(1:nx,1:ny,lbz:nz) = dvdx(1:nx,1:ny,lbz:nz) - dudy(1:nx,1:ny,lbz:nz)
    vortz(1:nx,1:ny,lbz:nz) = interp_to_w_grid( vortx(1:nx,1:ny,lbz:nz), lbz)
    vortx(1:nx,1:ny,lbz:nz) = dwdy(1:nx,1:ny,lbz:nz) - dvdz(1:nx,1:ny,lbz:nz)
    vorty(1:nx,1:ny,lbz:nz) = dudz(1:nx,1:ny,lbz:nz) - dwdx(1:nx,1:ny,lbz:nz)

    if (coord == 0) then
        vortz(1:nx,1:ny, 1) = 0._rprec
    end if

    ! Common file name for all output types
    call string_splice(fname, path //'output/vort.', jt_total)

#if defined(PPCGNS) && defined(PPMPI)
    ! Write CGNS Output
    call string_concat(fname, '.cgns')
    call write_parallel_cgns(fname,nx,ny, nz - nz_end, nz_tot,                 &
        (/ 1, 1,   (nz-1)*coord + 1 /),                                        &
        (/ nx, ny, (nz-1)*(coord+1) + 1 - nz_end /),                           &
        x(1:nx) , y(1:ny) , zw(1:(nz-nz_end) ),                                &
        3, (/ 'VorticityX', 'VorticityY', 'VorticityZ' /),                     &
        (/ vortx(1:nx,1:ny,1:(nz-nz_end)), vorty(1:nx,1:ny,1:(nz-nz_end)),     &
        vortz(1:nx,1:ny,1:(nz-nz_end)) /) )

#else
    ! Write binary Output
    call string_concat(fname, bin_ext)
    open(unit=13, file=fname, form='unformatted', convert=write_endian,        &
        access='direct', recl=nx*ny*nz*rprec)
    write(13,rec=1) vortx(:nx,:ny,1:nz)
    write(13,rec=2) vorty(:nx,:ny,1:nz)
    write(13,rec=3) vortz(:nx,:ny,1:nz)
    close(13)
#endif

    deallocate(vortx, vorty, vortz)

    ! Compute pressure
    allocate(pres_real(nx,ny,lbz:nz))
    pres_real(1:nx,1:ny,lbz:nz) = 0._rprec

    ! Calculate real pressure
    pres_real(1:nx,1:ny,lbz:nz) = p(1:nx,1:ny,lbz:nz)                          &
        - 0.5 * ( u(1:nx,1:ny,lbz:nz)**2                                       &
        + interp_to_uv_grid( w(1:nx,1:ny,lbz:nz), lbz)**2                      &
        + v(1:nx,1:ny,lbz:nz)**2 )

    ! Common file name for all output types
    call string_splice(fname, path //'output/pres.', jt_total)

#if defined(PPCGNS) && defined(PPMPI)
    ! Write CGNS Output
    call string_concat(fname, '.cgns')
    call write_parallel_cgns(fname, nx, ny, nz - nz_end, nz_tot,               &
        (/ 1, 1,   (nz-1)*coord + 1 /),                                        &
        (/ nx, ny, (nz-1)*(coord+1) + 1 - nz_end /),                           &
        x(1:nx) , y(1:ny) , z(1:(nz-nz_end) ),                                 &
        1, (/ 'Pressure' /), (/ pres_real(1:nx,1:ny,1:(nz-nz_end)) /) )

#else
    ! Write binary Output
    call string_concat(fname, bin_ext)
    open(unit=13, file=fname, form='unformatted', convert=write_endian,        &
        access='direct', recl=nx*ny*nz*rprec)
    write(13,rec=1) pres_real(:nx,:ny,1:nz)
    close(13)
#endif

     deallocate(pres_real)

#ifdef PPSCALARS
    ! Common file name for all output types
    call string_splice(fname, path //'output/theta.', jt_total)
#if defined(PPCGNS) && defined(PPMPI)
    ! Write CGNS Output
    call string_concat(fname, '.cgns')
    call write_parallel_cgns(fname, nx, ny, nz - nz_end, nz_tot,               &
     (/ 1, 1,   (nz-1)*coord + 1 /),                                           &
     (/ nx, ny, (nz-1)*(coord+1) + 1 - nz_end /),                              &
     x(1:nx) , y(1:ny) , z(1:(nz-nz_end) ),                                    &
     1, (/ 'Theta' /), (/ theta(1:nx,1:ny,1:(nz-nz_end)) /) )
#else
    ! Write binary Output
    call string_concat(fname, bin_ext)
    open(unit=13, file=fname, form='unformatted', convert=write_endian,        &
     access='direct', recl=nx*ny*nz*rprec)
    write(13,rec=1) theta(:nx,:ny,1:nz)
    write(13,rec=2) pi_x(:nx,:ny,1:nz)!GN
    write(13,rec=3) pi_y(:nx,:ny,1:nz)!GN
    write(13,rec=4) pi_z(:nx,:ny,1:nz)!GN
    close(13)
#endif
#endif

!  Write instantaneous x-plane values
elseif(itype==3) then

    allocate(ui(1,ny,nz), vi(1,ny,nz), wi(1,ny,nz))

    !  Loop over all xplane locations
    do i = 1, xplane_nloc
        do k = 1, nz
            do j = 1, ny
                ui(1,j,k) = linear_interp(u(xplane(i) % istart,j,k),    &
                     u(xplane(i) % istart+1,j,k), dx, xplane(i) % ldiff)
                vi(1,j,k) = linear_interp(v(xplane(i) % istart,j,k),    &
                     v(xplane(i) % istart+1,j,k), dx, xplane(i) % ldiff)
                wi(1,j,k) = linear_interp(w_uv(xplane(i) % istart,j,k), &
                     w_uv(xplane(i) % istart+1,j,k), dx, &
                     xplane(i) % ldiff)
            end do
        end do

        ! Common file name portion for all output types
        call string_splice(fname, path // 'output/vel.x-', xplane_loc(i), '.', jt_total)

#if defined(PPCGNS) && defined(PPMPI)
        ! Write CGNS Output
        call string_concat(fname, '.cgns')
        call write_parallel_cgns (fname,1,ny, nz - nz_end, nz_tot,     &
                        (/ 1, 1,   (nz-1)*coord + 1 /),                &
                        (/ 1, ny, (nz-1)*(coord+1) + 1 - nz_end /),    &
                    xplane_loc(i:i) , y(1:ny) , z(1:(nz-nz_end) ),     &
              3, (/ 'VelocityX', 'VelocityY', 'VelocityZ' /),          &
              (/ ui(1,1:ny,1:(nz-nz_end)), vi(1,1:ny,1:(nz-nz_end)),   &
                 wi(1,1:ny,1:(nz-nz_end)) /) )

#else
        ! Write binary output
        call string_concat(fname, bin_ext)
        open(unit=13,file=fname,form='unformatted',convert=write_endian, access='direct',recl=ny*nz*rprec)
        write(13,rec=1) ui
        write(13,rec=2) vi
        write(13,rec=3) wi
        close(13)
#endif
    end do

    deallocate(ui,vi,wi)

!  Write instantaneous y-plane values
elseif(itype==4) then

    allocate(ui(nx,1,nz), vi(nx,1,nz), wi(nx,1,nz))

    !  Loop over all yplane locations
    do j = 1, yplane_nloc
        do k = 1, nz
            do i = 1, nx

                ui(i,1,k) = linear_interp(u(i,yplane(j) % istart,k),           &
                     u(i,yplane(j) % istart+1,k), dy, yplane(j) % ldiff)
                vi(i,1,k) = linear_interp(v(i,yplane(j) % istart,k),           &
                     v(i,yplane(j) % istart+1,k), dy, yplane(j) % ldiff)
                wi(i,1,k) = linear_interp(w_uv(i,yplane(j) % istart,k),        &
                     w_uv(i,yplane(j) % istart+1,k), dy, yplane(j) % ldiff)
            end do
        end do

        ! Common file name portion for all output types
        call string_splice(fname, path // 'output/vel.y-', yplane_loc(j), '.', &
             jt_total)

#if defined(PPCGNS) && defined(PPMPI)
        call string_concat(fname, '.cgns')
        call write_parallel_cgns (fname,nx,1, nz - nz_end, nz_tot,             &
            (/ 1, 1,   (nz-1)*coord + 1 /),                                    &
            (/ nx, 1, (nz-1)*(coord+1) + 1 - nz_end /),                        &
            x(1:nx) , yplane_loc(j:j) , z(1:(nz-nz_end) ),                     &
            3, (/ 'VelocityX', 'VelocityY', 'VelocityZ' /),                    &
            (/ ui(1:nx,1,1:(nz-nz_end)), vi(1:nx,1,1:(nz-nz_end)),             &
            wi(1:nx,1,1:(nz-nz_end)) /) )
#else
        ! Write binary output
        call string_concat(fname, bin_ext)
        open(unit=13,file=fname,form='unformatted',convert=write_endian, access='direct',recl=nx*nz*rprec)
        write(13,rec=1) ui
        write(13,rec=2) vi
        write(13,rec=3) wi
        close(13)
#endif

    end do

    deallocate(ui,vi,wi)

!  Write instantaneous z-plane values
elseif (itype==5) then

    allocate(ui(nx,ny,1), vi(nx,ny,1), wi(nx,ny,1))

    !  Loop over all zplane locations
    do k = 1, zplane_nloc
        ! Common file name portion for all output types
        call string_splice(fname, path // 'output/vel.z-',                     &
                zplane_loc(k), '.', jt_total)

#ifdef PPCGNS
        call string_concat(fname, '.cgns')
#endif

#ifdef PPMPI
        if(zplane(k) % coord == coord) then
            do j = 1, Ny
                do i = 1, Nx
                    ui(i,j,1) = linear_interp(u(i,j,zplane(k) % istart),       &
                         u(i,j,zplane(k) % istart+1), dz, zplane(k) % ldiff)
                    vi(i,j,1) = linear_interp(v(i,j,zplane(k) % istart),       &
                         v(i,j,zplane(k) % istart+1), dz, zplane(k) % ldiff)
                    wi(i,j,1) = linear_interp(w_uv(i,j,zplane(k) % istart),    &
                         w_uv(i,j,zplane(k) % istart+1), dz, zplane(k) % ldiff)
                end do
            end do

#ifdef PPCGNS
            call warn("inst_write","Z plane writting is currently disabled.")
!            ! Write CGNS Data
!            ! Only the processor with data writes, the other one is written
!            ! using null arguments with 'write_null_cgns'
!            call write_parallel_cgns (fname ,nx, ny, 1, 1,                     &
!                (/ 1, 1,   1 /),                                               &
!                (/ nx, ny, 1 /),                                               &
!                x(1:nx) , y(1:ny) , zplane_loc(k:k), 3,                        &
!                (/ 'VelocityX', 'VelocityY', 'VelocityZ' /),                   &
!                (/ ui(1:nx,1:ny,1), vi(1:nx,1:ny,1), wi(1:nx,1:ny,1) /) )
#else
            call string_concat(fname, bin_ext)
            open(unit=13,file=fname,form='unformatted',convert=write_endian,   &
                            access='direct',recl=nx*ny*1*rprec)
            write(13,rec=1) ui(1:nx,1:ny,1)
            write(13,rec=2) vi(1:nx,1:ny,1)
            write(13,rec=3) wi(1:nx,1:ny,1)
            close(13)
#endif
!
! #ifdef PPMPI
!         else
! #ifdef PPCGNS
!            write(*,*) "At write_null_cgns"
!            call write_null_cgns (fname ,nx, ny, 1, 1,                         &
!            (/ 1, 1,   1 /),                                                   &
!            (/ nx, ny, 1 /),                                                   &
!            x(1:nx) , y(1:ny) , zplane_loc(k:k), 3,                            &
!            (/ 'VelocityX', 'VelocityY', 'VelocityZ' /) )
!#endif
        end if
#endif
    end do
    deallocate(ui,vi,wi)

!  Instantaneous y-averaged  write for entire domain
elseif(itype==6) then
    ! Common file name for all output types
    call string_splice(fname, path //'output/velxz.', jt_total)

#if defined(PPCGNS) && defined(PPMPI)
    ! Write CGNS Output
    call string_concat(fname, '.cgns')
    call write_parallel_cgns(fname, nx, ny, nz - nz_end, nz_tot,               &
        (/ 1, 1,   (nz-1)*coord + 1 /),                                        &
        (/ nx, ny, (nz-1)*(coord+1) + 1 - nz_end /),                           &
        x(1:nx) , y(1:ny) , z(1:(nz-nz_end) ),                                 &
        3, (/ 'VelocityX', 'VelocityY', 'VelocityZ' /),                        &
        (/ u(1:nx,1:ny,1:(nz-nz_end)), v(1:nx,1:ny,1:(nz-nz_end)),             &
         w_uv(1:nx,1:ny,1:(nz-nz_end)) /) )
#else
    ! Write binary Output
    call string_concat(fname, bin_ext)
    open(unit=13, file=fname, form='unformatted', convert=write_endian,        &
        access='direct', recl=nx*nz*rprec)
    write(13,rec=1) SUM(u(:nx,:ny,1:nz),DIM=2)/ny
    write(13,rec=2) SUM(v(:nx,:ny,1:nz),DIM=2)/ny
    write(13,rec=3) SUM(w_uv(:nx,:ny,1:nz),DIM=2)/ny
    close(13)
#endif

    ! Common file name for all output types
    call string_splice(fname, path //'output/sgsxz.', jt_total)

#if defined(PPCGNS) && defined(PPMPI)
    ! Write CGNS Output
    call string_concat(fname, '.cgns')
    call write_parallel_cgns(fname, nx, ny, nz - nz_end, nz_tot,               &
        (/ 1, 1,   (nz-1)*coord + 1 /),                                        &
        (/ nx, ny, (nz-1)*(coord+1) + 1 - nz_end /),                           &
        x(1:nx) , y(1:ny) , z(1:(nz-nz_end) ),                                 &
        5, (/ 'Cs_opt2', 'LM', 'MM','QN','NN','Nu_t' /),                        &
        (/ Cs_opt2(1:nx,1:ny,1:(nz-nz_end)), F_LM(1:nx,1:ny,1:(nz-nz_end)),             &
         F_MM(1:nx,1:ny,1:(nz-nz_end)),F_QN(1:nx,1:ny,1:(nz-nz_end)), F_NN(1:nx,1:ny,1:(nz-nz_end)), Nu_t(1:nx,1:ny,1:(nz-nz_end)) /) )
#else
    ! Write binary Output
    call string_concat(fname, bin_ext)
    open(unit=13, file=fname, form='unformatted', convert=write_endian,        &
        access='direct', recl=nx*nz*rprec)
    write(13,rec=1) SUM( Cs_opt2(:nx,:ny,1:nz),DIM=2)/ny
    write(13,rec=2) SUM( F_LM(:nx,:ny,1:nz), DIM=2 )/ny
    write(13,rec=3) SUM( F_MM(:nx,:ny,1:nz), DIM=2 )/ny
    write(13,rec=4) SUM( F_QN(:nx,:ny,1:nz), DIM=2 )/ny
    write(13,rec=5) SUM( F_NN(:nx,:ny,1:nz), DIM=2 )/ny
    write(13,rec=6) SUM( Nu_t(:nx,:ny,1:nz), DIM=2 )/ny
    close(13)
#endif

    ! Common file name for all output types !GN
    call string_splice(fname, path //'output/tau_sgsxz.', jt_total)

#if defined(PPCGNS) && defined(PPMPI)
    ! Write CGNS Output
    call string_concat(fname, '.cgns')
    call write_parallel_cgns(fname, nx, ny, nz - nz_end, nz_tot,               &
        (/ 1, 1,   (nz-1)*coord + 1 /),                                        &
        (/ nx, ny, (nz-1)*(coord+1) + 1 - nz_end /),                           &
        x(1:nx) , y(1:ny) , z(1:(nz-nz_end) ),                                 &
        5, (/ 'txx', 'tyy', 'tzz','txy','txz','tyz' /),                        &
        (/ txx(1:nx,1:ny,1:(nz-nz_end)), tyy(1:nx,1:ny,1:(nz-nz_end)),             &
         tzz(1:nx,1:ny,1:(nz-nz_end)),txy(1:nx,1:ny,1:(nz-nz_end)), txz(1:nx,1:ny,1:(nz-nz_end)), tyz(1:nx,1:ny,1:(nz-nz_end)) /) )
#else
    ! Write binary Output
    call string_concat(fname, bin_ext)
    open(unit=13, file=fname, form='unformatted', convert=write_endian,        &
        access='direct', recl=nx*nz*rprec)
    write(13,rec=1) SUM(txx(:nx,:ny,1:nz), DIM = 2)/ny
    write(13,rec=2) SUM(tyy(:nx,:ny,1:nz), DIM = 2)/ny
    write(13,rec=3) SUM(tzz(:nx,:ny,1:nz), DIM = 2)/ny
    write(13,rec=4) SUM(txy(:nx,:ny,1:nz), DIM = 2)/ny
    write(13,rec=5) SUM(txz(:nx,:ny,1:nz), DIM = 2)/ny
    write(13,rec=6) SUM(tyz(:nx,:ny,1:nz), DIM = 2)/ny
    close(13)
#endif

    ! Common file name for all output types !AA
    call string_splice(fname, path //'output/wavestressxz.', jt_total)

#if defined(PPCGNS) && defined(PPMPI)
    ! Write CGNS Output
    call string_concat(fname, '.cgns')
    call write_parallel_cgns(fname, nx, ny, nz - nz_end, nz_tot,               &
        (/ 1, 1,   (nz-1)*coord + 1 /),                                        &
        (/ nx, ny, (nz-1)*(coord+1) + 1 - nz_end /),                           &
        x(1:nx) , y(1:ny) , z(1:(nz-nz_end) ),                                 &
        5, (/ 'txx', 'tyy', 'tzz','txy','txz','tyz' /),                        &
        (/ fxa(1:nx,1:ny,1:(nz-nz_end)), fya(1:nx,1:ny,1:(nz-nz_end)),             &
         fza(1:nx,1:ny,1:(nz-nz_end)) /) )
#else
    if (use_sea_drag_model .and. use_exp_decay) then
            fxa(:nx,:ny,1) = fd_u(:nx,:ny)*dz 
            fya(:nx,:ny,1) = fd_v(:nx,:ny)*dz
    ! Write binary Output
    call string_concat(fname, bin_ext)
    open(unit=13, file=fname, form='unformatted', convert=write_endian,        &
        access='direct', recl=nx*nz*rprec)
    write(13,rec=1) SUM(fxa(:nx,:ny,1:nz), DIM = 2)/ny
    write(13,rec=2) SUM(fya(:nx,:ny,1:nz), DIM = 2)/ny
    write(13,rec=3) SUM(fza(:nx,:ny,1:nz), DIM = 2)/ny
    close(13)
    endif
#endif

    ! Compute vorticity
    allocate(vortx(nx,ny,lbz:nz), vorty(nx,ny,lbz:nz), vortz(nx,ny,lbz:nz))
    vortx(1:nx,1:ny,lbz:nz) = 0._rprec
    vorty(1:nx,1:ny,lbz:nz) = 0._rprec
    vortz(1:nx,1:ny,lbz:nz) = 0._rprec

    ! Use vorticityx as an intermediate step for performing uv-w interpolation
    ! Vorticity is written in w grid
    vortx(1:nx,1:ny,lbz:nz) = dvdx(1:nx,1:ny,lbz:nz) - dudy(1:nx,1:ny,lbz:nz)
    vortz(1:nx,1:ny,lbz:nz) = interp_to_w_grid( vortx(1:nx,1:ny,lbz:nz), lbz)
    vortx(1:nx,1:ny,lbz:nz) = dwdy(1:nx,1:ny,lbz:nz) - dvdz(1:nx,1:ny,lbz:nz)
    vorty(1:nx,1:ny,lbz:nz) = dudz(1:nx,1:ny,lbz:nz) - dwdx(1:nx,1:ny,lbz:nz)

    if (coord == 0) then
        vortz(1:nx,1:ny, 1) = 0._rprec
    end if

    ! Common file name for all output types
    call string_splice(fname, path //'output/vortxz.', jt_total)

#if defined(PPCGNS) && defined(PPMPI)
    ! Write CGNS Output
    call string_concat(fname, '.cgns')
    call write_parallel_cgns(fname,nx,ny, nz - nz_end, nz_tot,                 &
        (/ 1, 1,   (nz-1)*coord + 1 /),                                        &
        (/ nx, ny, (nz-1)*(coord+1) + 1 - nz_end /),                           &
        x(1:nx) , y(1:ny) , zw(1:(nz-nz_end) ),                                &
        3, (/ 'VorticityX', 'VorticityY', 'VorticityZ' /),                     &
        (/ vortx(1:nx,1:ny,1:(nz-nz_end)), vorty(1:nx,1:ny,1:(nz-nz_end)),     &
        vortz(1:nx,1:ny,1:(nz-nz_end)) /) )

#else
    ! Write binary Output
    call string_concat(fname, bin_ext)
    open(unit=13, file=fname, form='unformatted', convert=write_endian,        &
        access='direct', recl=nx*nz*rprec)
    write(13,rec=1) SUM(vortx(:nx,:ny,1:nz),DIM=2)/ny
    write(13,rec=2) SUM(vorty(:nx,:ny,1:nz),DIM=2)/ny
    write(13,rec=3) SUM(vortz(:nx,:ny,1:nz),DIM=2)/ny
    close(13)
#endif

    deallocate(vortx, vorty, vortz)

    ! Compute pressure
    allocate(pres_real(nx,ny,lbz:nz))
    pres_real(1:nx,1:ny,lbz:nz) = 0._rprec

    ! Calculate real pressure
    pres_real(1:nx,1:ny,lbz:nz) = p(1:nx,1:ny,lbz:nz)                          &
        - 0.5 * ( u(1:nx,1:ny,lbz:nz)**2                                       &
        + interp_to_uv_grid( w(1:nx,1:ny,lbz:nz), lbz)**2                      &
        + v(1:nx,1:ny,lbz:nz)**2 )

    ! Common file name for all output types
    call string_splice(fname, path //'output/presxz.', jt_total)

#if defined(PPCGNS) && defined(PPMPI)
    ! Write CGNS Output
    call string_concat(fname, '.cgns')
    call write_parallel_cgns(fname, nx, ny, nz - nz_end, nz_tot,               &
        (/ 1, 1,   (nz-1)*coord + 1 /),                                        &
        (/ nx, ny, (nz-1)*(coord+1) + 1 - nz_end /),                           &
        x(1:nx) , y(1:ny) , z(1:(nz-nz_end) ),                                 &
        1, (/ 'Pressure' /), (/ pres_real(1:nx,1:ny,1:(nz-nz_end)) /) )

#else
    ! Write binary Output
    call string_concat(fname, bin_ext)
    open(unit=13, file=fname, form='unformatted', convert=write_endian,        &
        access='direct', recl=nx*nz*rprec)
    write(13,rec=1) SUM(pres_real(:nx,:ny,1:nz), DIM =2)/ny
    close(13)
#endif

     deallocate(pres_real)

#ifdef PPSCALARS
    ! Common file name for all output types
    call string_splice(fname, path //'output/thetaxz.', jt_total)
#if defined(PPCGNS) && defined(PPMPI)
    ! Write CGNS Output
    call string_concat(fname, '.cgns')
    call write_parallel_cgns(fname, nx, ny, nz - nz_end, nz_tot,               &
     (/ 1, 1,   (nz-1)*coord + 1 /),                                           &
     (/ nx, ny, (nz-1)*(coord+1) + 1 - nz_end /),                              &
     x(1:nx) , y(1:ny) , z(1:(nz-nz_end) ),                                    &
     1, (/ 'Theta' /), (/ theta(1:nx,1:ny,1:(nz-nz_end)) /) )
#else
    ! Write binary Output
    call string_concat(fname, bin_ext)
    open(unit=13, file=fname, form='unformatted', convert=write_endian,        &
     access='direct', recl=nx*nz*rprec)
    write(13,rec=1) SUM(theta(:nx,:ny,1:nz),DIM=2)/ny
    write(13,rec=2) SUM(pi_x(:nx,:ny,1:nz),DIM=2)/ny    !GN
    write(13,rec=3) SUM(pi_y(:nx,:ny,1:nz),DIM=2)/ny    !GN
    write(13,rec=4) SUM(pi_z(:nx,:ny,1:nz),DIM=2)/ny    !GN
    close(13)
#endif
#endif
else
    write(*,*) 'Error: itype not specified properly to inst_write!'
    stop
end if

deallocate(w_uv)
nullify(x,y,z,zw)

#ifdef PPLVLSET
contains
!*******************************************************************************
subroutine force_tot()
!*******************************************************************************
#ifdef PPMPI
use mpi_defs, only : mpi_sync_real_array, MPI_SYNC_DOWN
#endif
implicit none

! Zero bogus values
fx(:,:,nz) = 0._rprec
fy(:,:,nz) = 0._rprec
fz(:,:,nz) = 0._rprec

!  Sum both the induced and applied forces
allocate(fx_tot(nx,ny,nz), fy_tot(nx,ny,nz), fz_tot(nx,ny,nz))

#ifdef PPTURBINES
fx_tot = fxa(1:nx,1:ny,1:nz)
fy_tot = fya(1:nx,1:ny,1:nz)
fz_tot = fza(1:nx,1:ny,1:nz)

#elif PPATM
fx_tot = fxa(1:nx,1:ny,1:nz)
fy_tot = fya(1:nx,1:ny,1:nz)
fz_tot = fza(1:nx,1:ny,1:nz)

#elif PPLVLSET
fx_tot = fx(1:nx,1:ny,1:nz)+fxa(1:nx,1:ny,1:nz)
fy_tot = fy(1:nx,1:ny,1:nz)+fya(1:nx,1:ny,1:nz)
fz_tot = fz(1:nx,1:ny,1:nz)+fza(1:nx,1:ny,1:nz)
#else
fx_tot = 0._rprec
fy_tot = 0._rprec
fz_tot = 0._rprec
#endif
!! AA BOC
if (use_sea_drag_model .and. use_exp_decay) then

fx_tot = fxa(1:nx,1:ny,1:nz)
fy_tot = fya(1:nx,1:ny,1:nz)
fz_tot = fza(1:nx,1:ny,1:nz)
endif

!! AA EOC
#ifdef PPMPI
!  Sync forces
call mpi_sync_real_array( fx_tot, 1, MPI_SYNC_DOWN )
call mpi_sync_real_array( fy_tot, 1, MPI_SYNC_DOWN )
call mpi_sync_real_array( fz_tot, 1, MPI_SYNC_DOWN )
#endif

! Put fz_tot on uv-grid
fz_tot(1:nx,1:ny,1:nz) = interp_to_uv_grid( fz_tot(1:nx,1:ny,1:nz), 1 )

return
end subroutine force_tot
#endif

!*******************************************************************************
!subroutine pressure_sync()
!!*******************************************************************************
!use mpi_defs, only : mpi_sync_real_array, MPI_SYNC_DOWN
!use param, only : ld
!implicit none
!
!! Reset bogus values
!p(:,:,nz) = p(:,:,nz-1)
!dpdx(:,:,nz) = dpdx(:,:,nz-1)
!dpdy(:,:,nz) = dpdy(:,:,nz-1)
!dpdz(:,:,nz) = dpdz(:,:,nz-1)
!
!#ifdef PPMPI
!!  Sync pressure
!call mpi_sync_real_array( p, 0 , MPI_SYNC_DOWN )
!call mpi_sync_real_array( dpdx, 1 , MPI_SYNC_DOWN )
!call mpi_sync_real_array( dpdy, 1 , MPI_SYNC_DOWN )
!call mpi_sync_real_array( dpdz, 1 , MPI_SYNC_DOWN )
!#endif
!
!return
!end subroutine pressure_sync
!
!!*******************************************************************************
!subroutine RHS_sync()
!!*******************************************************************************
!use param, only : ld
!use mpi_defs, only : mpi_sync_real_array, MPI_SYNC_DOWN
!implicit none
!
!! Reset bogus values
!RHSx(:,:,nz) = RHSx(:,:,nz-1)
!RHSy(:,:,nz) = RHSy(:,:,nz-1)
!RHSz(:,:,nz) = RHSz(:,:,nz-1)
!
!#ifdef PPMPI
!!  Sync RHS
!call mpi_sync_real_array( RHSx, 0 , MPI_SYNC_DOWN )
!call mpi_sync_real_array( RHSy, 0 , MPI_SYNC_DOWN )
!call mpi_sync_real_array( RHSz, 0 , MPI_SYNC_DOWN )
!#endif
!
!return
!end subroutine RHS_sync

end subroutine inst_write

!*******************************************************************************
subroutine checkpoint ()
!*******************************************************************************
use iwmles
use param, only : nz, checkpoint_file, tavg_calc, lbc_mom, L_x, L_y, L_z, path
#ifdef PPMPI
use param, only : comm, ierr
#endif
use sim_param, only : u, v, w, RHSx, RHSy, RHSz
use sgs_param, only : Cs_opt2, F_LM, F_MM, F_QN, F_NN
use param, only : jt_total, total_time, total_time_dim, dt,                    &
    use_cfl_dt, cfl, write_endian
use cfl_util, only : get_max_cfl
use string_util, only : string_concat
#if PPUSE_TURBINES
use turbines, only : turbines_checkpoint
#endif
#ifdef PPSCALARS
use scalars, only : scalars_checkpoint
#endif
use coriolis

! HIT Inflow
#ifdef PPHIT
use hit_inflow, only : hit_write_restart
#endif

implicit none
character(64) :: fname
real(rprec) :: cfl_w

fname = checkpoint_file
#ifdef PPMPI
call string_concat( fname, '.c', coord )
#endif

!  Open vel.out (lun_default in io) for final output
open(11, file=fname, form='unformatted', convert=write_endian,                 &
    status='unknown', position='rewind')
write (11) u(:, :, 1:nz), v(:, :, 1:nz), w(:, :, 1:nz),                        &
    RHSx(:, :, 1:nz), RHSy(:, :, 1:nz), RHSz(:, :, 1:nz),                      &
    Cs_opt2(:,:,1:nz), F_LM(:,:,1:nz), F_MM(:,:,1:nz),                         &
    F_QN(:,:,1:nz), F_NN(:,:,1:nz)
close(11)

! Open grid.out for final output
if (coord == 0) then
    open(11, file= path // 'grid.out', form='unformatted', convert=write_endian)
    write(11) nproc, Nx, Ny, Nz, L_x, L_y, L_z
    close(11)
end if

#ifdef PPMPI
call mpi_barrier( comm, ierr )
#endif

! Checkpoint time averaging restart data
if ( tavg_calc .and. tavg%initialized ) call tavg%checkpoint()

! Write time and current simulation state
! Set the current cfl to a temporary (write) value based whether CFL is
! specified or must be computed
if( use_cfl_dt ) then
    cfl_w = cfl
else
    cfl_w = get_max_cfl()
end if

!xiang check point for iwm
if(lbc_mom==3)then
    if (coord == 0) call iwm_checkPoint()
end if

#ifdef PPHIT
    if (coord == 0) call hit_write_restart()
#endif

#if PPUSE_TURBINES
call turbines_checkpoint
#endif

#ifdef PPSCALARS
call scalars_checkpoint
#endif

call coriolis_finalize()

!  Update total_time.dat after simulation
if (coord == 0) then
    !--only do this for true final output, not intermediate recording
    open (1, file=fcumulative_time)
    write(1, *) jt_total, total_time, total_time_dim, dt, cfl_w
    close(1)
end if

end subroutine checkpoint

!*******************************************************************************
subroutine output_final()
!*******************************************************************************
use param, only : tavg_calc
implicit none

! Perform final checkpoing
call checkpoint()

!  Check if average quantities are to be recorded
if (tavg_calc .and. tavg%initialized ) call tavg%finalize()

end subroutine output_final

!*******************************************************************************
subroutine output_init ()
!*******************************************************************************
!
!  This subroutine allocates the memory for arrays used for statistical
!  calculations
!
use param, only : dx, dy, dz, nz, lbz
use param, only : point_calc, point_nloc, point_loc
use param, only : xplane_calc, xplane_nloc, xplane_loc
use param, only : yplane_calc, yplane_nloc, yplane_loc
use param, only : zplane_calc, zplane_nloc, zplane_loc
use grid_m
use functions, only : cell_indx
use stat_defs, only : point, xplane, yplane, zplane
implicit none

integer :: i,j,k
real(rprec), pointer, dimension(:) :: x,y,z


#ifdef PPMPI
! This adds one more element to the last processor (which contains an extra one)
! Processor nproc-1 has data from 1:nz
! Rest of processors have data from 1:nz-1
if ( coord == nproc-1 ) then
    nz_end = 0
else
    nz_end = 1
end if
#else
nz_end = 0
#endif

nullify(x,y,z)

x => grid % x
y => grid % y
z => grid % z

! Initialize information for x-planar stats/data
if (xplane_calc) then
    allocate(xplane(xplane_nloc))
    xplane(:) % istart = -1
    xplane(:) % ldiff = 0.

    !  Compute istart and ldiff
    do i = 1, xplane_nloc
        xplane(i) % istart = cell_indx('i', dx, xplane_loc(i))
        xplane(i) % ldiff = xplane_loc(i) - x(xplane(i) % istart)
    end do
end if

! Initialize information for y-planar stats/data
if (yplane_calc) then
    allocate(yplane(yplane_nloc))
    yplane(:) % istart = -1
    yplane(:) % ldiff = 0.

    !  Compute istart and ldiff
    do j = 1, yplane_nloc
        yplane(j) % istart = cell_indx('j', dy, yplane_loc(j))
        yplane(j) % ldiff = yplane_loc(j) - y(yplane(j) % istart)
    end do
end if

! Initialize information for z-planar stats/data
if(zplane_calc) then
    allocate(zplane(zplane_nloc))

    !  Initialize
    zplane(:) % istart = -1
    zplane(:) % ldiff = 0.
    zplane(:) % coord = -1

    !  Compute istart and ldiff
    do k = 1, zplane_nloc

#ifdef PPMPI
        if (zplane_loc(k) >= z(1) .and. zplane_loc(k) < z(nz)) then
            zplane(k) % coord = coord
            zplane(k) % istart = cell_indx('k',dz,zplane_loc(k))
            zplane(k) % ldiff = zplane_loc(k) - z(zplane(k) % istart)
        end if
#else
        zplane(k) % coord = 0
        zplane(k) % istart = cell_indx('k',dz,zplane_loc(k))
        zplane(k) % ldiff = zplane_loc(k) - z(zplane(k) % istart)
#endif
    end do
end if

!  Open files for instantaneous writing
if (point_calc) then
    allocate(point(point_nloc))

    !  Intialize the coord values
    ! (-1 shouldn't be used as coord so initialize to this)
    point % coord=-1
    point % fid = -1

    do i = 1, point_nloc
        !  Find the processor in which this point lives
#ifdef PPMPI
        if (point_loc(i)%xyz(3) >= z(1) .and. point_loc(i)%xyz(3) < z(nz)) then
#endif

            point(i) % coord = coord

            point(i) % istart = cell_indx('i',dx,point_loc(i)%xyz(1))
            point(i) % jstart = cell_indx('j',dy,point_loc(i)%xyz(2))
            point(i) % kstart = cell_indx('k',dz,point_loc(i)%xyz(3))

            point(i) % xdiff = point_loc(i)%xyz(1) - x(point(i) % istart)
            point(i) % ydiff = point_loc(i)%xyz(2) - y(point(i) % jstart)
            point(i) % zdiff = point_loc(i)%xyz(3) - z(point(i) % kstart)

#ifdef PPMPI
        end if
#endif
    end do
end if

nullify(x,y,z)

end subroutine output_init

end module io
