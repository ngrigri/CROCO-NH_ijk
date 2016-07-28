#include "cppdefs.h"
#ifdef NHMG
module mg_define_matrix

  use mg_mpi
  use mg_tictoc
  use mg_namelist
  use mg_grids
  use mg_zr_zw
  use mg_mpi_exchange
  use mg_gather
  use mg_netcdf_out

  implicit none

  interface define_matrices
     module procedure           &
          define_matrices_topo
  end interface define_matrices

  !NG comment: constants in a mg_cst.f90 file ?
  real(kind=rp), parameter :: one  = 1._rp
  real(kind=rp), parameter :: eigh = one/8._rp
  real(kind=rp), parameter :: qrt  = 0.25_rp
  real(kind=rp), parameter :: hlf  = 0.5_rp

contains

  !-------------------------------------------------------------------------  
  subroutine define_matrices_topo(dx, dy, zeta, h)

    real(kind=rp), dimension(:,:), intent(in) :: dx
    real(kind=rp), dimension(:,:), intent(in) :: dy
    real(kind=rp), dimension(:,:), intent(in) :: zeta
    real(kind=rp), dimension(:,:), intent(in) :: h

    integer(kind=ip)::  lev

    real(kind=rp), dimension(:,:), pointer :: dxf
    real(kind=rp), dimension(:,:), pointer :: dxc

    real(kind=rp), dimension(:,:), pointer :: dyf
    real(kind=rp), dimension(:,:), pointer :: dyc

    real(kind=rp), dimension(:,:), pointer :: zetaf
    real(kind=rp), dimension(:,:), pointer :: zetac

    real(kind=rp), dimension(:,:), pointer :: hf
    real(kind=rp), dimension(:,:), pointer :: hc

    real(kind=rp), dimension(:,:,:), pointer :: zrc
    real(kind=rp), dimension(:,:,:), pointer :: zwc

    integer(kind=ip) :: nz,ny,nx
    integer(kind=ip) :: nzf,nyf,nxf
    integer(kind=ip) :: nyc,nxc

    if (myrank==0) write(*,*)'- define matrices from topography (h):'

    do lev = 1, nlevs

       if (myrank==0) write(*,*)'   lev=',lev

       nx=grid(lev)%nx
       ny=grid(lev)%ny
       nz=grid(lev)%nz

       if (lev == 1) then

          grid(lev)%dx  (0:nx+1,0:ny+1) = dx
          grid(lev)%dy  (0:nx+1,0:ny+1) = dy
          grid(lev)%zeta(0:nx+1,0:ny+1) = zeta
          grid(lev)%h   (0:nx+1,0:ny+1) = h

       else
          nxf =grid(lev-1)%nx
          nyf =grid(lev-1)%ny
          nzf =grid(lev-1)%nz

          dxf   => grid(lev-1)%dx
          dyf   => grid(lev-1)%dy
          zetaf => grid(lev-1)%zeta
          hf    => grid(lev-1)%h

          if (grid(lev)%gather == 1) then
             nxc= nx/grid(lev)%ngx
             nyc= ny/grid(lev)%ngy

             allocate(  dxc(0:nxc+1,0:nyc+1))
             allocate(  dyc(0:nxc+1,0:nyc+1))
             allocate(zetac(0:nxc+1,0:nyc+1))
             allocate(   hc(0:nxc+1,0:nyc+1))
          else
             nxc = nx
             nyc = ny
             dxc   => grid(lev)%dx
             dyc   => grid(lev)%dy
             zetac => grid(lev)%zeta
             hc    => grid(lev)%h
          endif

          if ((aggressive).and.(lev==1)) then

             write(*,*) ' define matrices (aggressive).and.(lev==1) not available !'
             STOP

          else

             ! Coarsen dx, dy and h
             dxc(1:nxc,1:nyc) = hlf      * ( &
                  dxf(1:nxf  :2,1:nyf  :2) + &
                  dxf(1:nxf  :2,2:nyf+1:2) + &
                  dxf(2:nxf+1:2,1:nyf  :2) + &
                  dxf(2:nxf+1:2,2:nyf+1:2) )

             dyc(1:nxc,1:nyc) = hlf      * ( &
                  dyf(1:nxf  :2,1:nyf  :2) + &
                  dyf(1:nxf  :2,2:nyf+1:2) + &
                  dyf(2:nxf+1:2,1:nyf  :2) + &
                  dyf(2:nxf+1:2,2:nyf+1:2) )

             zetac(1:nxc,1:nyc)  = qrt      * ( &
                  zetaf(1:nxf  :2,1:nyf  :2)  + &
                  zetaf(1:nxf  :2,2:nyf+1:2)  + &
                  zetaf(2:nxf+1:2,1:nyf  :2)  + &
                  zetaf(2:nxf+1:2,2:nyf+1:2)  )

             hc(1:nxc,1:nyc)  = qrt      * ( &
                  hf(1:nxf  :2,1:nyf  :2)  + &
                  hf(1:nxf  :2,2:nyf+1:2)  + &
                  hf(2:nxf+1:2,1:nyf  :2)  + &
                  hf(2:nxf+1:2,2:nyf+1:2)  )

          endif ! aggressive

          if (grid(lev)%gather == 1) then

             call gather(lev,dxc,grid(lev)%dx)
             call gather(lev,dyc,grid(lev)%dy)
             call gather(lev,zetac,grid(lev)%zeta)
             call gather(lev, hc,grid(lev)%h)

             deallocate(dxc)
             deallocate(dyc)
             deallocate(zetac)
             deallocate( hc)
          endif

       endif ! lev == 1

       call fill_halo(lev, grid(lev)%dx  ,nhalo=1)
       call fill_halo(lev, grid(lev)%dy  ,nhalo=1)
       call fill_halo(lev, grid(lev)%zeta,nhalo=2)
       call fill_halo(lev, grid(lev)%h   ,nhalo=2)

       zrc => grid(lev)%zr
       zwc => grid(lev)%zw

       ! Compute zr and zw
       call setup_zr_zw                (  & 
            hlim,nhtheta_b,nhtheta_s,     &
            grid(lev)%zeta,grid(lev)%h,   &  ! input args
            grid(lev)%zr, grid(lev)%zw,   &  ! output args
            coord_type='new_s_coord'      )  ! optional

       if (netcdf_output) then
          call write_netcdf(grid(lev)%dx  ,vname='dx'  ,netcdf_file_name='dx.nc'  ,rank=myrank,iter=lev)
          call write_netcdf(grid(lev)%dy  ,vname='dy'  ,netcdf_file_name='dy.nc'  ,rank=myrank,iter=lev)
          call write_netcdf(grid(lev)%zeta,vname='zeta',netcdf_file_name='zeta.nc',rank=myrank,iter=lev)
          call write_netcdf(grid(lev)%h   ,vname='h'   ,netcdf_file_name='h.nc'   ,rank=myrank,iter=lev)
          call write_netcdf(grid(lev)%zr  ,vname='zr'  ,netcdf_file_name='zr.nc'  ,rank=myrank,iter=lev)
          call write_netcdf(grid(lev)%zw  ,vname='zw'  ,netcdf_file_name='zw.nc'  ,rank=myrank,iter=lev)
       endif

       ! Define matrix coefficients from dx, dy, zr and zw coarsened
       call define_matrix(lev, grid(lev)%dx, grid(lev)%dy, grid(lev)%zr, grid(lev)%zw)

    enddo ! lev

  end subroutine define_matrices_topo

  !-----------------------------------------------------------------------------------
  subroutine define_matrix(lev, dx, dy, zr, zw)

    integer(kind=ip),intent(in):: lev
    real(kind=rp), dimension(:,:),   pointer, intent(in) :: dx, dy
    real(kind=rp), dimension(:,:,:), pointer, intent(in) :: zr, zw

    ! Define matrix coefficients cA
    ! Coefficients are stored in order of diagonals
    ! cA(1,:,:,:)      -> p(i  ,j  ,k  )
    ! cA(2,:,:,:)      -> p(i  ,j  ,k-1)
    ! cA(3,:,:,:)      -> p(i  ,j-1,k+1)
    ! cA(4,:,:,:)      -> p(i  ,j-1,k  )
    ! cA(5,:,:,:)      -> p(i  ,j-1,k-1)
    ! cA(6,:,:,:)      -> p(i-1,j  ,k+1)
    ! cA(7,:,:,:)      -> p(i-1,j  ,k  )
    ! cA(8,:,:,:)      -> p(i-1,j  ,k-1)

    integer(kind=ip):: i,j, k
    integer(kind=ip):: nx, ny, nz

    real(kind=rp) :: Arz
    real(kind=rp), dimension(:,:,:),   pointer :: dzw
    real(kind=rp), dimension(:,:,:),   pointer :: zydx,zxdy
    real(kind=rp), dimension(:,:,:),   pointer :: cw
    real(kind=rp), dimension(:,:,:,:), pointer :: cA

    nx = grid(lev)%nx
    ny = grid(lev)%ny
    nz = grid(lev)%nz

    cA => grid(lev)%cA 

    !- Used in compute_rhs -!
    if (lev == 1) then
       dzw => grid(lev)%dzw

       dzw(:,:,1) = zr(:,:,1) - zw(:,:,1) !!
       do k = 2,nz
          do j = 0,ny+1
             do i = 0,nx+1
                dzw(i,j,k) = zr(i,j,k) - zr(i,j,k-1) !!  cell height at w-points
             enddo
          enddo
       enddo
       do j = 0,ny+1
          do i = 0,nx+1
             dzw(i,j,nz+1) = zw(i,j,nz+1) - zr(i,j,nz) !!
          enddo
       enddo

       !! Slopes in x- and y-direction defined at rho-points
       zxdy => grid(lev)%zxdy
       zydx => grid(lev)%zydx
       do k = 1, nz
          do j = 0,ny+1
             do i = 0,nx+1
                zydx(i,j,k) = hlf * (( zr(i  ,j+1,k) - zr(i  ,j-1,k) ) / dy(i,j) ) * dx(i,j)
                zxdy(i,j,k) = hlf * (( zr(i+1,j  ,k) - zr(i-1,j  ,k) ) / dx(i,j) ) * dy(i,j)
             enddo
          enddo
       enddo
    endif

    !- Used also in compute_rhs -!
    cw => grid(lev)%cw

    k=1
    do j = 0,ny+1
       do i = 0,nx+1
          Arz = dx(i,j)*dy(i,j)
          cw(i,j,k) = ( Arz / (zr(i,j,k)-zw(i,j,k)) ) * &
               (one + &
               ( hlf * (zw(i+1,j  ,k)-zw(i-1,j  ,k)) / dx(i,j) )**2 + &
               ( hlf * (zw(i  ,j+1,k)-zw(i  ,j-1,k)) / dy(i,j) )**2 )
       enddo
    enddo

    do k = 2,nz
       do j = 0,ny+1
          do i = 0,nx+1
             cw(i,j,k) = ( Arz / (zr(i,j,k)-zr(i,j,k-1)) ) * &
                  (one + &
                  ( hlf * (zw(i+1,j  ,k)-zw(i-1,j  ,k)) / dx(i,j) )**2 + &
                  ( hlf * (zw(i  ,j+1,k)-zw(i  ,j-1,k)) / dy(i,j) )**2 )
          enddo
       enddo
    enddo

    k=nz+1
    do j = 0,ny+1
       do i = 0,nx+1
          cw(i,j,k) = ( Arz / (zw(i,j,k)-zr(i,j,k-1)) ) * &
               (one + &
               ( hlf * (zw(i+1,j  ,k)-zw(i-1,j  ,k)) / dx(i,j) )**2 + &
               ( hlf * (zw(i  ,j+1,k)-zw(i  ,j-1,k)) / dy(i,j) )**2 )
       enddo
    enddo

    !! interaction coeff with neighbours
    !XX
    !---------!
    !- K = 1 -! lower level
    !---------!
    k = 1
    do j = 1,ny+1
       do i = 1,nx

          cA(3,i,j,k) = qrt * ( &
               ( hlf * (zr(i,j+1,k+1)-zr(i,j-1,k+1)) / dy(i,j  ) ) * dx(i,j  ) + &
               ( hlf * (zr(i,j  ,k  )-zr(i,j-2,k  )) / dy(i,j-1) ) * dx(i,j-1) )

          cA(4,i,j,k) =                                                                &
                                ! couples with j-1
               ( qrt                                                                 * &
               ( zw(i,j,k+1) - zw(i,j,k) + zw(i,j-1,k+1) - zw(i,j-1,k) )             * &
               (dx(i,j)+dx(i,j-1)) )/ (hlf * (dy(i,j)+dy(i,j-1)))                      & 
                                ! topo terms 
               - ( (( hlf * (zr(i,j+1,k)-zr(i,j-1,k)) / dy(i,j) ) * dx(i,j))**2 / &
               ( cw(i,j  ,k) + cw(i,j  ,k+1) )     &
               +   (( hlf * (zr(i,j,k)-zr(i,j-2,k)) / dy(i,j-1) ) * dx(i,j-1))**2 / &
               ( cw(i,j-1,k) + cw(i,j-1,k+1) )   ) &
                                ! from j,k cross terms
               - qrt * ( &
               (( hlf * (zr(i,j,k  )-zr(i,j-2,k  )) / dy(i,j-1) ) * dx(i,j-1)) -  &
               ( hlf * (zr(i,j+1,k  )-zr(i,j-1,k  )) / dy(i,j) ) * dx(i,j))  
       enddo
    enddo

    do j = 1,ny
       do i = 1,nx+1

          cA(6,i,j,k) = qrt * ( &
               ( hlf * (zr(i+1,j,k+1)-zr(i-1,j,k+1)) / dx(i,j) ) * dy(i,j) + &
               ( hlf * (zr(i  ,j,k  )-zr(i-2,j,k  )) / dx(i-1,j) ) * dy(i-1,j) )  ! couples with k+1 i-1

          cA(7,i,j,k) =                                                                &
                                ! couples with i-1
               ( qrt                                                                 * &
               ( zw(i,j,k+1) - zw(i,j,k) + zw(i-1,j,k+1) - zw(i-1,j,k) )             * &
               (dy(i,j)+dy(i-1,j)) )                                                 / &
               ( hlf * (dx(i,j)+dx(i-1,j)) )                                           &  
                                ! topo terms
               - ( (( hlf * (zr(i+1,j  ,k)-zr(i-1,j  ,k)) / dx(i,j) ) * dy(i,j))**2 / &
               ( cw(i,j,k  ) + cw(i,j,  k+1) )     &
               +   (( hlf * (zr(i,j  ,k)-zr(i-2,j  ,k)) / dx(i-1,j) ) * dy(i-1,j))**2 / &
               ( cw(i-1,j,k) + cw(i-1,j,k+1) )   ) &
                                ! from i,k cross terms
               - qrt * ( &
               ( hlf * (zr(i  ,j,k) - zr(i-2,j,k)) / dx(i-1,j) ) * dy(i-1,j) - &
               ( hlf * (zr(i+1,j,k) - zr(i-1,j,k)) / dx(i  ,j) ) * dy(i  ,j) )
       enddo
    enddo

    do j = 0,ny 
       do i = 1,nx+1
          ! only for k==1, couples with j+1,i-1
          cA(5,i,j,k) = &
               + hlf * &
               (( hlf * (zr(i+1,j+1,k)-zr(i-1,j+1,k)) / dx(i,j+1) ) * dy(i,j+1)) * &
               (( hlf * (zr(i  ,j+2,k)-zr(i  ,j  ,k)) / dy(i,j+1) ) * dx(i,j+1)) / &
               ( cw(i,j+1,k  ) + cw(i,j+1,k+1))  &
               + hlf * &
               (( hlf * (zr(i,j  ,k)-zr(i-2,j  ,k)) / dx(i-1,j) ) * dy(i-1,j)) * &
               (( hlf * (zr(i-1,j+1,k)-zr(i-1,j-1,k)) / dy(i-1,j) ) * dx(i-1,j)) / &
               ( cw(i-1,j,k) + cw(i-1,j,k+1) )
       enddo
    enddo

    do j = 1,ny+1
       do i = 1,nx+1
          ! only for k==1, couples with j-1,i-1
          cA(8,i,j,k) = &
               - hlf * &
               (( hlf * (zr(i+1,j-1,k)-zr(i-1,j-1,k)) / dx(i,j-1) ) * dy(i,j-1)) * &
               (( hlf * (zr(i  ,j  ,k)-zr(i  ,j-2,k)) / dy(i,j-1) ) * dx(i,j-1)) / &
               (cw(i,j-1,k) + cw(i,j-1,k+1)) & 
               - hlf * &
               (( hlf * (zr(i  ,j  ,k)-zr(i-2,j  ,k)) / dx(i-1,j) ) * dy(i-1,j)) * &
               (( hlf * (zr(i-1,j+1,k)-zr(i-1,j-1,k)) / dy(i-1,j) ) * dx(i-1,j)) / &
               (cw(i-1,j,k) + cw(i-1,j,k+1)) 
       enddo
    enddo

    !XX
    !---------------!
    !- K = 2, nz-1 -! interior levels
    !---------------!

    do k = 2,nz-1 
       do j = 1,ny
          do i = 1,nx
             cA(2,i,j,k) =  cw(i,j,k)                                  ! couples with k-1
          enddo
       enddo
    enddo

    do k = 2,nz-1 
       do j = 1,ny+1
          do i = 1,nx
             cA(3,i,j,k) =  qrt * ( &
                  ( hlf * (zr(i,j+1,k+1)-zr(i,j-1,k+1)) / dy(i,j  ) ) * dx(i,j  ) + &
                  ( hlf * (zr(i,j  ,k  )-zr(i,j-2,k  )) / dy(i,j-1) ) * dx(i,j-1) )     ! couples with k+1 j-1
             cA(4,i,j,k) =  ( qrt * &
                  ( zw(i,j,k+1) - zw(i,j,k) + zw(i,j-1,k+1) - zw(i,j-1,k) ) * &
                  (dx(i,j)+dx(i,j-1)) ) / ( hlf * (dy(i,j)+dy(i,j-1)) ) ! couples with j-1
             cA(5,i,j,k) =- qrt * ( &
                  (( hlf * (zr(i,j+1,k-1)-zr(i,j-1,k-1)) / dy(i,j  ) ) * dx(i,j  )) + &
                  (( hlf * (zr(i,j  ,k  )-zr(i,j-2,k  )) / dy(i,j-1) ) * dx(i,j-1)) )     ! couples with k-1 j-1
          enddo
       enddo
    enddo

    do k = 2,nz-1 
       do j = 1,ny 
          do i = 1,nx+1
             cA(6,i,j,k) =  qrt * ( &
                  (( hlf * (zr(i+1,j  ,k+1)-zr(i-1,j  ,k+1)) / dx(i  ,j) ) * dy(i,j  )) + &
                  (( hlf * (zr(i  ,j  ,k  )-zr(i-2,j  ,k  )) / dx(i-1,j) ) * dy(i-1,j)) )     ! Couples with k+1 i-1
             cA(7,i,j,k) =   (qrt * &
                  ( zw(i,j,k+1) - zw(i,j,k) + zw(i-1,j,k+1) - zw(i-1,j,k) ) * &
                  (dy(i,j)+dy(i-1,j)) ) / ( hlf * (dx(i,j)+dx(i-1,j)) ) ! Couples with i-1
             cA(8,i,j,k) =- qrt * ( &
                  (( hlf * (zr(i+1,j,k-1)-zr(i-1,j,k-1)) / dx(i  ,j) ) * dy(i  ,j)) + &
                  (( hlf * (zr(i  ,j,k  )-zr(i-2,j,k  )) / dx(i-1,j) ) * dy(i-1,j))  )     ! Couples with k-1 i-1
          enddo
       enddo
    enddo

    !XX
    !----------!
    !- K = nz -! upper level
    !----------!
    k = nz

    do j = 1,ny 
       do i = 1,nx    
          cA(2,i,j,k) = cw(i,j,k)                                    ! couples with k-1
       enddo
    enddo

    do j = 1,ny+1
       do i = 1,nx
          cA(4,i,j,k) = ( qrt * &
               ( zw(i,j,k+1) - zw(i,j,k) + zw(i,j-1,k+1) - zw(i,j-1,k) ) * &
               (dx(i,j)+dx(i,j-1)) ) / ( hlf * (dy(i,j)+dy(i,j-1)) ) & ! couples with j-1
               + qrt * ( &
               - (( hlf * (zr(i,j  ,k)-zr(i,j-2,k)) / dy(i,j-1) ) * dx(i,j-1)) &
               + (( hlf * (zr(i,j+1,k)-zr(i,j-1,k)) / dy(i,j  ) ) * dx(i,j  )) )
          cA(5,i,j,k) =- qrt * ( &
               (( hlf * (zr(i,j+1,k-1)-zr(i,j-1,k-1)) / dy(i,j  ) ) * dx(i,j  )) + &
               (( hlf * (zr(i,j  ,k  )-zr(i,j-2,k  )) / dy(i,j-1) ) * dx(i,j-1)) )     ! couples with k-1 j-1
       enddo
    enddo


    do j = 1,ny 
       do i = 1,nx+1
          cA(7,i,j,k) = (qrt * &
               ( zw(i,j,k+1) - zw(i,j,k) + zw(i-1,j,k+1) - zw(i-1,j,k) ) * &
               (dy(i,j)+dy(i-1,j)) ) / ( hlf * (dx(i,j)+dx(i-1,j)) ) & ! Couples with i-1
               + qrt * ( &
               - (( hlf * (zr(i  ,j,k)-zr(i-2,j,k)) / dx(i-1,j)) * dy(i-1,j)) &
               + (( hlf * (zr(i+1,j,k)-zr(i-1,j,k)) / dx(i  ,j)) * dy(i  ,j)) )
          cA(8,i,j,k) =- qrt * ( &
               (( hlf * (zr(i+1,j,k-1)-zr(i-1,j,k-1)) / dx(i  ,j) ) * dy(i  ,j)) + &
               (( hlf * (zr(i  ,j,k  )-zr(i-2,j,k  )) / dx(i-1,j) ) * dy(i-1,j)) )     ! Couples with k-1 i-1
       enddo
    enddo

    call fill_halo(lev,cA)

    !! interaction coeff with itself
    k = 1 !lower level
    do j = 1,ny
       do i = 1,nx
          cA(1,i,j,k) =                         &
               -cA(2,i,j,k+1)                   &
               -cA(4,i,j,k) - cA(4,i  ,j+1,k  ) &
               -cA(7,i,j,k) - cA(7,i+1,j  ,k  ) &
               -cA(6,i,j,k) - cA(8,i+1,j  ,k+1) &
               -cA(3,i,j,k) - cA(5,i  ,j+1,k+1) &
               -cA(5,i,j,k) - cA(5,i+1,j-1,k  ) &
               -cA(8,i,j,k) - cA(8,i+1,j+1,k  )
       enddo
    enddo

    do k = 2,nz-1 !interior levels
       do j = 1,ny
          do i = 1,nx
             cA(1,i,j,k) = &
                  -cA(2,i,j,k) - cA(2,i  ,j  ,k+1) &
                  -cA(4,i,j,k) - cA(4,i  ,j+1,k  ) &
                  -cA(7,i,j,k) - cA(7,i  ,j  ,k+1) &
                  -cA(6,i,j,k) - cA(6,i+1,j  ,k-1) &
                  -cA(8,i,j,k) - cA(8,i+1,j  ,k+1) & 
                  -cA(3,i,j,k) - cA(3,i  ,j+1,k-1) &
                  -cA(5,i,j,k) - cA(5,i  ,j+1,k+1)
          enddo
       enddo
    enddo

    k = nz !upper level
    do j = 1,ny
       do i = 1,nx
          cA(1,i,j,k) =                    &
               - cA(2,i,j,k)               &
               - cw(i,j,k+1)               &
               + hlf * (hlf * (zr(i+2,j  ,k) - zr(i  ,j  ,k)) / dx(i+1,j  )) * dy(i+1,j  ) & 
               - hlf * (hlf * (zr(i  ,j  ,k) - zr(i-2,j  ,k)) / dx(i-1,j  )) * dy(i-1,j  ) & 
               + hlf * (hlf * (zr(i  ,j+2,k) - zr(i  ,j  ,k)) / dy(i  ,j+1)) * dx(i  ,j+1) & 
               - hlf * (hlf * (zr(i  ,j  ,k) - zr(i  ,j-2,k)) / dy(i  ,j-1)) * dx(i  ,j-1) &
               - cA(4,i,j,k) - cA(4,i,j+1,k) &
               - cA(7,i,j,k) - cA(7,i+1,j,k) &
               - cA(6,i+1,j,k-1)           &
               - cA(8,i,j,k)               &
               - cA(3,i,j+1,k-1)           &
               - cA(5,i,j,k)
       enddo ! i
    enddo ! j

    if (netcdf_output) then
       if (myrank==0) write(*,*)'       write cA in a netcdf file'
       call write_netcdf(grid(lev)%cA,vname='ca',netcdf_file_name='cA.nc',rank=myrank,iter=lev)
    endif

  end subroutine define_matrix

end module mg_define_matrix
#else
        module mg_define_matrix_empty
        end module
#endif
