#include "cppdefs.h"
#ifdef NHMG
module mg_correct_uvw

  use mg_mpi
  use mg_tictoc
  use mg_namelist
  use mg_grids
  use mg_mpi_exchange
  use mg_netcdf_out

  implicit none

contains
  !-------------------------------------------------------------------------     
  subroutine correct_uvw(u,v,w)

    real(kind=rp), dimension(:,:,:), pointer, intent(inout) :: u,v,w

    integer(kind=ip):: i, j, k
    integer(kind=ip):: nx, ny, nz

    real(kind=rp), dimension(:,:)  , pointer :: dx,dy
    real(kind=rp), dimension(:,:,:), pointer :: zr,zw
    real(kind=rp), dimension(:,:,:), pointer :: p
    real(kind=rp) :: dxu,dyv
    real(kind=rp) :: dzw

    !NG comment: constants in a mg_cst.f90 file ?
    real(kind=rp), parameter :: two  = 2._rp
    real(kind=rp), parameter :: one  = 1._rp
    real(kind=rp), parameter :: hlf  = 0.5_rp
    real(kind=rp), parameter :: qrt  = 0.25_rp
    real(kind=rp), parameter :: zero  = 0._rp

    nx = grid(1)%nx
    ny = grid(1)%ny
    nz = grid(1)%nz

    dx => grid(1)%dx
    dy => grid(1)%dy
    zr => grid(1)%zr
    zw => grid(1)%zw

    if (myrank==0) write(*,*)'- correct u,v,w:'

    !! Correct
    p => grid(1)%p

    do k = 1,nz
       do j = 0,ny+1 
          do i = 1,nx+1
             dxu = hlf * (dx(i,j)+dx(i-1,j))
             u(i,j,k) = u(i,j,k) - one / dxu * (p(i,j,k)-p(i-1,j,k))
          enddo
       enddo
    enddo

    do k = 1,nz
       do j = 1,ny+1 
          do i = 0,nx+1

             dyv = hlf * (dy(i,j)+dy(i,j-1))

             v(i,j,k) = v(i,j,k) - one / dyv * (p(i,j,k)-p(i,j-1,k))

          enddo
       enddo
    enddo

    do k = 2,nz !interior and upper levels
       do j = 0,ny+1
          do i = 0,nx+1
             dzw = zr(i,j,k)-zr(i,j,k-1)
             w(i,j,k-1) = w(i,j,k-1) - one / dzw * (p(i,j,k)-p(i,j,k-1))
          enddo

       enddo
    enddo

    k = nz+1 !surface
    do j = 0,ny+1
       do i = 0,nx+1
          dzw = zw(i,j,i,k)-zr(i,j,k-1)
          w(i,j,k-1) = w(i,j,k-1) - one / dzw * (-p(i,j,k-1))
       enddo
    enddo

  end subroutine correct_uvw

end module mg_correct_uvw
#else
        module mg_correct_uvw_empty
        end module
#endif
