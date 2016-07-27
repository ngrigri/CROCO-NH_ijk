#include "cppdefs.h"
#ifdef NHMG
module mg_compute_rhs

  use mg_mpi
  use mg_tictoc
  use mg_namelist
  use mg_grids
  use mg_mpi_exchange
  use mg_netcdf_out

  implicit none

contains
  !-------------------------------------------------------------------------     
  subroutine compute_rhs(u,v,w)

    real(kind=rp), dimension(:,:,:), pointer, intent(in) :: u,v,w

    integer(kind=ip):: i, j, k
    integer(kind=ip):: nx, ny, nz

    real(kind=rp), dimension(:,:)  , pointer :: dx,dy
    real(kind=rp), dimension(:,:,:), pointer :: zr,zw

    real(kind=rp) :: Arx, Ary
    real(kind=rp), dimension(:,:,:), pointer :: dzw
    real(kind=rp), dimension(:,:,:), pointer :: zydx,zxdy
    real(kind=rp), dimension(:,:,:), pointer :: cw
    real(kind=rp), dimension(:,:,:), pointer :: uf,vf,wf,rhs

    !NG comment: constants in a mg_cst.f90 file ?
    real(kind=rp), parameter :: two  = 2._rp
    real(kind=rp), parameter :: one  = 1._rp
    real(kind=rp), parameter :: hlf  = 0.5_rp
    real(kind=rp), parameter :: qrt  = 0.25_rp
    real(kind=rp), parameter :: zero = 0._rp

    integer(kind=ip), save :: iter_rhs=0
    iter_rhs = iter_rhs + 1

    nx = grid(1)%nx
    ny = grid(1)%ny
    nz = grid(1)%nz

    dx => grid(1)%dx
    dy => grid(1)%dy
    zr => grid(1)%zr
    zw => grid(1)%zw

    if (myrank==0) write(*,*)'- compute rhs:'

    !- Cell heights (computed in define matices)
    dzw => grid(1)%dzw

    !- Slopes in x- and y-direction defined at rho-points (computed in define matices)
    zxdy => grid(1)%zxdy
    zydx => grid(1)%zydx

    !- (computed in define matices)
    cw  => grid(1)%cw

    !--------------!
    !- DIVERGENCE -!
    !--------------!
    !- Divergence is stored in rhs array pointer
    !- It is calculated progressively using first uf, 
    !- then using vf and at the e wf
    rhs => grid(1)%b
    rhs(:,:,:) = zero

    !------!
    !- UF -!
    !------!
    !- if (myrank==0) write(*,*)'- compute rhs: uf'

    uf => grid(1)%dummy3Dnz

    k = 1 ! lower level
    do j = 1,ny ! South to North
       do i = 1,nx+1  ! West  to East

          uf(i,j,k) =  &
               qrt                                                         & 
               * ( zw(i,j,k+1) - zw(i,j,k) + zw(i-1,j,k+1) - zw(i-1,j,k) ) &
               * ( dy(i,j) + dy(i-1,j) ) * u(i,j,k)                        &
               - qrt * (                                                   &
               + zxdy(i  ,j,k) * dzw(i  ,j,k+1) * w(i  ,j,k+1-1)   &
               + zxdy(i-1,j,k) * dzw(i-1,j,k+1) * w(i-1,j,k+1-1) ) &
               -( &
               + zxdy(i  ,j,k)*zxdy(  i,j,k)/(cw(i,j,k  )+cw(i  ,j,k+1)) &
               + zxdy(i-1,j,k)*zxdy(i-1,j,k)/(cw(i-1,j,k)+cw(i-1,j,k+1)) ) &
               * (hlf * ( dx(i,j) + dx(i-1,j) )) * u(i,j,k) &
               -( & 
               + zxdy(i,j,k) * zydx(i,j,k)/(cw(i,j,k)+cw(i,j,k+1))   &
               * hlf * ( &
               hlf * ( dy(i,j  ) + dy(i,j-1) ) * v(i,j  ,k) + &
               hlf * ( dy(i,j+1) + dy(i,j  ) ) * v(i,j+1,k))   & 
               + zxdy(i-1,j,k) * zydx(i-1,j,k)/(cw(i-1,j,k)+cw(i-1,j,k+1))   &
               * hlf * ( &
               hlf * ( dy(i-1,j  ) + dy(i-1,j-1) ) * v(i-1,j  ,k) + &
               hlf * ( dy(i-1,j+1) + dy(i-1,j  ) ) * v(i-1,j+1,k)) )

       enddo
    enddo


    do k = 2,nz-1 !interior levels
       do j = 1,ny ! South to North
          do i = 1,nx+1  ! West  to East

             uf(i,j,k) =  qrt                                                  & 
                  * ( zw(i,j,k+1) - zw(i,j,k) + zw(i-1,j,k+1) - zw(i-1,j,k,j) ) &
                  * ( dy(i,j) + dy(i-1,j) ) * u(i,j,k) &
                  - qrt * ( &
                  + zxdy(i  ,j,k) * dzw(i  ,j,k  ) * w(i  ,j,k-1  ) &
                  + zxdy(i  ,j,k) * dzw(i  ,j,k+1) * w(i  ,j,k+1-1) &
                  + zxdy(i-1,j,k) * dzw(i-1,j,k  ) * w(i-1,j,k-1  ) &
                  + zxdy(i-1,j,k) * dzw(i-1,j,k+1) * w(i-1,j,k+1-1) &
                  )  ! umask
          enddo

       enddo
    enddo

    k = nz !upper level

    do j = 1,ny ! South to North
       do i = 1,nx+1  ! West  to East

          uf(i,j,k) = &
               qrt                                                       * & 
               ( zw(i,j,k+1) - zw(i,j,k) + zw(i-1,j,k+1) - zw(i-1,j,k) ) * &
               ( dy(i,j) + dy(i-1,j) ) &
               * u(i,j,k) &
               - qrt * ( &
               + zxdy(i  ,j,k)*       dzw(i  ,j,k  ) * w(i  ,j,k-1  ) &
               + zxdy(i  ,j,k)* two * dzw(i  ,j,k+1) * w(i  ,j,k+1-1) &
               + zxdy(i-1,j,k)*       dzw(i-1,j,k  ) * w(i-1,j,k-1  ) &
               + zxdy(i-1,j,k)* two * dzw(i-1,j,k+1) * w(i-1,j,k+1-1) &
               )  ! umask

       enddo
    enddo

    call fill_halo(1,uf,lbc_null='u')

    if (netcdf_output) then
       call write_netcdf(uf,vname='uf',netcdf_file_name='uf.nc',rank=myrank,iter=iter_rhs)
    endif

    do k = 1,nz
       do j = 1,ny 
          do i = 1,nx
             rhs(i,j,k) = uf(i+1,j,k) - uf(i,j,k) 
          enddo
       enddo
    enddo

    uf => null()

    !------!
    !- VF -!
    !------!
    !- if (myrank==0) write(*,*)'- compute rhs: vf'

    vf => grid(1)%dummy3Dnz

    k = 1 !lower level
    do j = 1,ny+1
       do i = 1,nx
          vf(i,j,k) = &
               qrt                                                         &
               * ( zw(i,j,k+1) - zw(i,j,k) + zw(i,j-1,k+1) - zw(i,j-1,k) ) &
               * ( dx(i,j) + dx(i,j-1) ) * v(i,j,k) &
               - qrt * ( &
               + zydx(i,j  ,k) * dzw(i,j  ,k+1) * w(i,j,k+1-1  ) &
               + zydx(i,j-1,k) * dzw(i,j-1,k+1) * w(i,j-1,k+1-1) &
               ) &  
               -( &
               + zydx(i,j  ,k) * zydx(i,j  ,k)/(cw(i,j  ,k)+cw(i,j  ,k+1)) &
               + zydx(i,j-1,k) * zydx(i,j-1,k)/(cw(i,j-1,k)+cw(i,j-1,k+1)) ) &
               * hlf * ( dy(i,j) + dy(i,j-1) ) * v(i,j,k) &
               - ( &
               + zxdy(i,j,k) * zydx(i,j,k)/(cw(i,j,k)+cw(i,j,k+1))  &
               * hlf * ( &
               hlf * ( dx(i  ,j) + dx(i-1,j) ) * u(i  ,j,k) + &
               hlf * ( dx(i+1,j) + dx(i  ,j) ) * u(i+1,j,k)) &
               + zxdy(i,j-1,k) * zydx(i,j-1,k)/(cw(i,j-1,k)+cw(i,j-1,k+1))   &
               * hlf * ( &
               hlf * ( dx(i  ,j-1) + dx(i-1,j-1) ) * u(i  ,j-1,k) + &
               hlf * ( dx(i+1,j-1) + dx(i  ,j-1) ) * u(i+1,j-1,k)) &
               ) 
       enddo
    enddo

    do k = 2,nz-1 !interior levels
       do j = 1,ny+1
          do i = 1,nx
             vf(i,j,k) = &
                  qrt                                                       * &
                  ( zw(i,j,k+1) - zw(i,j,k) + zw(i,j-1,k+1) - zw(i,j-1,k) ) * &
                  ( dx(i,j) + dx(i,j-1) ) * v(i,j,k) &
                  - qrt * ( &
                  + zydx(i,j  ,k) * dzw(i,j  ,k  ) * w(i,j  ,k-1  ) &
                  + zydx(i,j  ,k) * dzw(i,j  ,k+1) * w(i,j  ,k+1-1) &
                  + zydx(i,j-1,k) * dzw(i,j-1,k  ) * w(i,j-1,k-1  ) &
                  + zydx(i,j-1,k) * dzw(i,j-1,k+1) * w(i,j-1,k+1-1) &
                  )  !* vmask(j,i)
          enddo
       enddo
    enddo

    k = nz !upper level
    do j = 1,ny+1
       do i = 1,nx
          vf(i,j,k) =  &
               qrt                                                       * &
               ( zw(i,j,k+1) - zw(i,j,k) + zw(i,j-1,k+1) - zw(i,j-1,k) ) * &
               ( dx(i,j) + dx(i,j-1) ) &
               * v(i,j,k) &
               - qrt * ( &
               + zydx(i,j  ,k)*       dzw(i,j  ,k  ) * w(i,j  ,k-1  ) &
               + zydx(i,j  ,k)* two * dzw(i,j  ,k+1) * w(i,j  ,k+1-1) &
               + zydx(i,j-1,k)*       dzw(i,j-1,k  ) * w(i,j-1,k-1  ) &
               + zydx(i,j-1,k)* two * dzw(i,j-1,k+1) * w(i,j-1,k+1-1) &
               ) !* vmask(j,i)
       enddo
    enddo

    call fill_halo(1,vf,lbc_null='v')

    if (netcdf_output) then
       call write_netcdf(vf,vname='vf',netcdf_file_name='vf.nc',rank=myrank,iter=iter_rhs)
    endif

    do k = 1,nz
       do j = 1,ny 
          do i = 1,nx  
             rhs(i,j,k) =  rhs(i,j,k) + vf(i,j+1,k) - vf(i,j,k)
          enddo
       enddo
    enddo

    vf => null()

    !------!
    !- WF -!
    !------!
    !- if (myrank==0) write(*,*)'- compute rhs: wf'

    wf => grid(1)%dummy3Dnzp

    k = 1 !bottom
    do j = 1,ny
       do i = 1,nx
          wf(i,j,k) = zero
       enddo
    enddo

    do k = 2,nz !interior levels
       do j = 1,ny
          do i = 1,nx
             wf(i,j,k) = cw(i,j,k) * dzw(i,j,k) * w(i,j,k-1) &
                  - qrt * hlf * ( &
                  + zxdy(i,j,k  ) * ( dx(i  ,j) + dx(i-1,j) ) * u(i  ,j,k  ) &
                  + zxdy(i,j,k  ) * ( dx(i+1,j) + dx(i  ,j) ) * u(i+1,j,k  ) &
                  + zxdy(i,j,k-1) * ( dx(i  ,j) + dx(i-1,j) ) * u(i  ,j,k-1) &
                  + zxdy(i,j,k-1) * ( dx(i+1,j) + dx(i  ,j) ) * u(i+1,j,k-1) ) &
                  - qrt * hlf * ( &
                  + zydx(i,j,k  ) * ( dy(i,j  ) + dy(i,j-1) ) * v(i,j  ,k  ) &
                  + zydx(i,j,k  ) * ( dy(i,j+1) + dy(i,j  ) ) * v(i,j+1,k  ) &
                  + zydx(i,j,k-1) * ( dy(i,j  ) + dy(i,j-1) ) * v(i,j  ,k-1) &
                  + zydx(i,j,k-1) * ( dy(i,j+1) + dy(i,j  ) ) * v(i,j+1,k-1) )
          enddo
       enddo
    enddo

    k = nz+1 !surface
    do j = 1,n
       do i = 1,nxy
          wf(i,j,k) = cw(i,j,k) * dzw(i,j,k) * w(i,j,k-1)&
               - hlf * hlf *( &
               + zxdy(i,j,k-1) * ( dx(i  ,j) + dx(i-1,j) ) * u(i  ,j,k-1)   &
               + zxdy(i,j,k-1) * ( dx(i+1,j) + dx(i  ,j) ) * u(i+1,j,k-1) ) &
               - hlf * hlf * ( &
               + zydx(i,j,k-1) * ( dy(i,j  ) + dy(i,j-1) ) * v(i,j,k-1)     &
               + zydx(i,j,k-1) * ( dy(i,j+1) + dy(i,j  ) ) * v(i,j+1,k-1) )
       enddo
    enddo

    if (netcdf_output) then
       call write_netcdf(wf,vname='wf',netcdf_file_name='wf.nc',rank=myrank,iter=iter_rhs)
    endif

    do k = 1,nz
       do j = 1,ny 
          do i = 1,nx
             rhs(i,j,k) = rhs(i,j,k) + wf(i,j,k+1) - wf(i,j,k)
          enddo
       enddo
    enddo

    wf => null()

    !- if (myrank==0) write(*,*)'- compute rhs (finish)'

  end subroutine compute_rhs

end module mg_compute_rhs
#else
        module mg_compute_rhs_empty
        end module
#endif

