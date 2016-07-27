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

    integer(kind=ip):: k, j, i
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
    allocate(dz(nz,0:ny+1,0:nx+1))
    do i = 0,nx+1
       do j = 0,ny+1
          do k = 1,nz
             dz(k,j,i) = zw(k+1,j,i)-zw(k,j,i)
          enddo
       enddo
    enddo

    !! Cell widths
    allocate(dxu(0:ny+1,nx+1))
    do i = 1,nx+1
       do j = 0,ny+1
          dxu(j,i) = hlf * (dx(j,i)+dx(j,i-1))
       enddo
    enddo
    allocate(dyv(ny+1,0:nx+1))
    do i = 0,nx+1
       do j = 1,ny+1
          dyv(j,i) = hlf * (dy(j,i)+dy(j-1,i))
       enddo
    enddo

    !! Compute
    p => grid(1)%p

    ru(:,:) = 0._8
    rv(:,:) = 0._8

    do i = 1,nx+1
       do j = 0,ny+1 
!       do j = 1,ny

          do k = 1,nz
             ru(i,j) = ru(i,j) - hlf*(dz(k,j,i)+dz(k,j,i-1)) / dxu(j,i) *(p(k,j,i)-p(k,j,i-1))
          enddo

       enddo
    enddo
 
    do i = 0,nx+1
!    do i = 1,nx
       do j = 1,ny+1 

          do k = 1,nz
             rv(i,j) = rv(i,j) - hlf*(dz(k,j,i)+dz(k,j-1,i)) / dyv(j,i) *(p(k,j,i)-p(k,j-1,i))
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
