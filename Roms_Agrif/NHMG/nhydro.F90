#include "cppdefs.h"
#ifdef NHMG
module nhydro

  use mg_mpi
  use mg_grids
  use mg_namelist
  use mg_tictoc
  use mg_mpi_exchange
  use mg_netcdf_out
  use mg_compute_rhs
  use mg_correct_uvw
  use mg_solvers
  use mg_compute_barofrc

  implicit none

contains

  !--------------------------------------------------------------
  subroutine nhydro_init(nx,ny,nz,npxg,npyg)
      
    integer(kind=ip), intent(in) :: nx, ny, nz
    integer(kind=ip), intent(in) :: npxg, npyg

    if (myrank==0) write(*,*)' nhydro_init:'

    call mg_mpi_init()

    call read_nhnamelist(vbrank=myrank)

    call define_grids(npxg,npyg,nx,ny,nz)

    call define_neighbours()

    call print_grids()

  end subroutine nhydro_init

  !--------------------------------------------------------------
  subroutine nhydro_matrices(nx,ny,nz,dxa,dya,zetaa,ha,hc,theta_b,theta_s)

    integer(kind=ip), intent(in) :: nx, ny, nz

    real(kind=rp), dimension(0:nx+1,0:ny+1), target, intent(in) :: dxa, dya, zetaa, ha
    real(kind=rp),                           intent(in) :: hc, theta_b, theta_s

    real(kind=rp), dimension(:,:)  , pointer :: dx, dy, zeta, h

    if (myrank==0) write(*,*)' nhydro_matrices:'

    dx   => dxa
    dy   => dya
    zeta => zetaa
    h    => ha

    hlim      = hc
    nhtheta_b = theta_b
    nhtheta_s = theta_s

    call define_matrices(dx,dy,zeta,h)

  end subroutine nhydro_matrices

  !--------------------------------------------------------------
  subroutine nhydro_solve(nx,ny,nz,ua,va,wa,rua,rva)

    integer(kind=ip), intent(in) :: nx, ny, nz

    real(kind=rp), dimension(1:nx+1,0:ny+1,1:nz), target, intent(inout) :: ua
    real(kind=rp), dimension(0:nx+1,1:ny+1,1:nz), target, intent(inout) :: va
    real(kind=rp), dimension(0:nx+1,0:ny+1,0:nz), target, intent(inout) :: wa
    real(kind=rp), dimension(1:nx+1,0:ny+1),      target, intent(out)   :: rua
    real(kind=rp), dimension(0:nx+1,1:ny+1),      target, intent(out)   :: rva

    real(kind=rp), dimension(:,:,:), pointer :: u, v, w
    real(kind=rp), dimension(:,:)  , pointer :: ru, rv

    real(kind=rp)    :: tol
    integer(kind=ip) :: maxite
    integer(kind=ip) :: i,j,k

    integer(kind=ip), save :: iter_solve=0
    iter_solve = iter_solve + 1

    if (myrank==0) write(*,*)' nhydro_solve:'

    call tic(1,'nhydro_solve')

    tol    = solver_prec    ! solver_prec    is defined in the namelist file
    maxite = solver_maxiter ! solver_maxiter is defined in the namelist file

    u => ua
    v => va
    w => wa
    ru => rua
    rv => rva

    if (netcdf_output) then
       call write_netcdf(u,vname='uin',netcdf_file_name='uin.nc',rank=myrank,iter=iter_solve)
       call write_netcdf(v,vname='vin',netcdf_file_name='vin.nc',rank=myrank,iter=iter_solve)
       call write_netcdf(w,vname='win',netcdf_file_name='win.nc',rank=myrank,iter=iter_solve)
    endif

    !- step 1 - 
    call tic(1,'compute_rhs')
    call compute_rhs(u, v, w)
    call toc(1,'compute_rhs')

    if (netcdf_output) then
       call write_netcdf(grid(1)%b,vname='b',netcdf_file_name='b.nc',rank=myrank,iter=iter_solve)
    endif

    !- step 2 -
    call solve_p(tol,maxite)

    if (netcdf_output) then
       call write_netcdf(grid(1)%p,vname='p',netcdf_file_name='p.nc',rank=myrank,iter=iter_solve)
       call write_netcdf(grid(1)%r,vname='r',netcdf_file_name='r.nc',rank=myrank,iter=iter_solve)
    endif

    !- step 3 -
    call correct_uvw(u,v,w)

    if (netcdf_output) then
       call write_netcdf(u,vname='uout',netcdf_file_name='uout.nc',rank=myrank,iter=iter_solve)
       call write_netcdf(v,vname='vout',netcdf_file_name='vout.nc',rank=myrank,iter=iter_solve)
       call write_netcdf(w,vname='wout',netcdf_file_name='wout.nc',rank=myrank,iter=iter_solve)
    endif

    !- step 4 -
    call compute_barofrc(ru,rv)

    if (netcdf_output) then
       call write_netcdf(ru,vname='ru',netcdf_file_name='ru.nc',rank=myrank,iter=iter_solve)
       call write_netcdf(rv,vname='rv',netcdf_file_name='rv.nc',rank=myrank,iter=iter_solve)
    endif

    if (associated(u)) u => null()
    if (associated(v)) v => null()
    if (associated(w)) w => null()

    call toc(1,'nhydro_solve')	

    if (myrank==0) write(*,*)' nhydro_solve end !!!'

    stop
 
  end subroutine nhydro_solve

  !--------------------------------------------------------------
  subroutine nhydro_check_nondivergence(nx,ny,nz,ua,va,wa)

    integer(kind=ip), intent(in) :: nx, ny, nz

    real(kind=rp), dimension(1:nx+1,0:ny+1,1:nz), target, intent(inout) :: ua
    real(kind=rp), dimension(0:nx+1,1:ny+1,1:nz), target, intent(inout) :: va
    real(kind=rp), dimension(0:nx+1,0:ny+1,0:nz), target, intent(inout) :: wa

    real(kind=rp), dimension(:,:,:), pointer :: u, v, w

    integer(kind=ip), save :: iter_check=0
    iter_check = iter_check + 1

    u => ua
    v => va
    w => wa

    if (myrank==0) write(*,*)'- check non-divergence:'

    call compute_rhs(u,v,w)

    if (netcdf_output) then
       call write_netcdf(grid(1)%b,vname='b',netcdf_file_name='check.nc',rank=myrank,iter=iter_check)
    endif

  end subroutine nhydro_check_nondivergence

  !--------------------------------------------------------------
  subroutine nhydro_clean()

    real(kind=rp) :: tstart,tend,perf

    call cpu_time(tstart)

    call grids_dealloc()

    call print_tictoc()

    call cpu_time(tend)

    if (myrank == 0) write(*,*)'nhydro_clean time:',tend-tstart

  end subroutine nhydro_clean

end module nhydro
#else
        module nhydro_empty
        end module
#endif
