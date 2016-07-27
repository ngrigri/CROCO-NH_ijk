#include "cppdefs.h"
#ifdef NHMG
module mg_compute_barofrc

  use mg_mpi
  use mg_tictoc
  use mg_namelist
  use mg_grids
  use mg_mpi_exchange
  use mg_netcdf_out

  implicit none

contains
  !-------------------------------------------------------------------------     
  subroutine compute_barofrc(ru,rv)

    real(kind=rp), dimension(:,:)  , pointer, intent(out) :: ru,rv

    integer(kind=ip):: i, j, k
    integer(kind=ip):: nx, ny, nz

    real(kind=rp), dimension(:,:)  , pointer :: dx,dy
    real(kind=rp), dimension(:,:,:), pointer :: zw
    real(kind=rp), dimension(:,:)  , pointer :: dxu,dyv
    real(kind=rp), dimension(:,:,:), pointer :: dz
    real(kind=rp), dimension(:,:,:), pointer :: p

    !NG comment: constants in a mg_cst.f90 file ?
    real(kind=rp), parameter :: hlf  = 0.5_rp

    nx = grid(1)%nx
    ny = grid(1)%ny
    nz = grid(1)%nz

    dx => grid(1)%dx
    dy => grid(1)%dy
    zw => grid(1)%zw

    if (myrank==0) write(*,*)'- compute barofrc:'

    !! Cell heights
    allocate(dz(0:nx+1,0:ny+1,nz))
    do k = 1,nz
       do j = 0,ny+1
          do i = 0,nx+1
             dz(i,j,k) = zw(i,j,k+1)-zw(i,j,k)
          enddo
       enddo
    enddo

    !! Cell widths
    allocate(dxu(nx+1,0:ny+1))
    do j = 0,ny+1
       do i = 1,nx+1
          dxu(i,j) = hlf * (dx(i,j)+dx(i-1,j))
       enddo
    enddo
    allocate(dyv(0:nx+1,ny+1))
    do j = 1,ny+1
       do i = 0,nx+1
          dyv(i,j) = hlf * (dy(i,j)+dy(i,j-1))
       enddo
    enddo

    !! Compute
    p => grid(1)%p

    ru(:,:) = 0._8
    rv(:,:) = 0._8

    do k = 1,nz
       do j = 0,ny+1 
          do i = 1,nx+1
             ru(i,j) = ru(i,j) - hlf*(dz(i,j,k)+dz(i-1,j,k)) / dxu(i,j) *(p(i,j,k)-p(i-1,j,k))
          enddo

       enddo
    enddo

    do k = 1,nz
       do j = 1,ny+1 
          do i = 0,nx+1
             rv(i,j) = rv(i,j) - hlf*(dz(i,j,k)+dz(i,j-1,k)) / dyv(i,j) *(p(i,j,k)-p(i,j-1,k))
          enddo
       enddo
    enddo

    deallocate(dz)
    deallocate(dxu)
    deallocate(dyv)

  end subroutine compute_barofrc

end module mg_compute_barofrc
#else
        module mg_compute_barofrc_empty
        end module
#endif
