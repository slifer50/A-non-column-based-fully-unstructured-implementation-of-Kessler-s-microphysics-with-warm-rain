!----------------------------------------------------------------------!
!>@brief This module contains basic operators; derivatives (1st and 2nd),
!> curl, divergence and laplacian
!>@author Subroutines written by FX Giraldo.
!>Department of Applied Mathematics
!>Naval Postgraduate School
!>The module interface was created by Sasa Gabersek (NRL-MRY) on 11/2011.
!>@update FX Giraldo cleaned up all subroutines and added CGc, CGd, and DG
!>capabilities. Also added Contravariant forms to split horizontal and vertical
!>directions.
!----------------------------------------------------------------------!
module mod_interface

  use mod_basis, only: ngl, dpsi, &
       nglx, ngly, nglz, npts, &
       dpsix, dpsiy, dpsiz, &
       dpsix_tr, dpsiy_tr, dpsiz_tr, &
       wgl, wglx, wgly, wglz, &
       fx, fy, fz, fx_tr, fy_tr, fz_tr
  use mod_grid, only: coord, intma, intma_1d, npoin, nelem, nz
  use mod_input, only: nelz, geometry_type
  !use mod_initial, only: kvector
  use mod_metrics, only: &
       ksi_x, ksi_y, ksi_z, &
       eta_x, eta_y, eta_z, &
       zeta_x, zeta_y, zeta_z, &
       jac, xjac, massinv, massinv_1d
  use mod_parallel, only : num_send_recv_total

  public :: &
       mass_multiply, &
       compute_gradient, &
       compute_gradient_q, &
       compute_divergence, &
       compute_vertical_derivative, &
       compute_vertical_derivative_kessler_contravariant, &
       compute_vorticity, &
       compute_curl, &  
       compute_derivative_column, &
       compute_local_gradient_v2, &
       compute_local_gradient_v3, &       
       compute_local_gradient_v3_ksi_eta, &
       compute_local_gradient_v3_zeta, &       
       compute_local_gradient_transpose_v3, &
       compute_local_gradient_transpose_v3_ksi_eta, &
       compute_local_gradient_transpose_v3_zeta, &              
       compute_local_gradient_filter_v3, &
       perform_mxm_operation, &
       compute_local_curl, &
       compute_canonical_gradient, &
       mod_interface_compute_strain_derivative_IBP, &
       compute_volume_integral, &
       lcontinuous_grad_q
       
  private

  !----------------------------------------------------------------------!
  logical :: lcontinuous_grad_q = .false.

  !----------------------------------------------------------------------!

  interface compute_volume_integral
     module procedure compute_volume_integral_field, compute_volume_integral_product
  end interface compute_volume_integral

contains

  !----------------------------------------------------------------------!
  !>@brief This subroutine integrates a variable over the domain - locally
  !>@author  Sohail Reddy 09/26/2020
  !----------------------------------------------------------------------!
  subroutine mass_multiply(q_in,ndim)

    implicit none

    !global arrays
    real, intent(inout) :: q_in(ndim,npoin)
    integer, intent(in) :: ndim

    !local arrays
    integer :: inode(npts)
    real    :: wq, jac_e(npts), q_tmp(ndim,npoin), qq(ndim,npts)
    integer :: i, j, k, e, m, ip, ii

    !Initialize
    q_tmp=0.0

    !loop thru the elements
    do e=1,nelem

       !Store Element Variables
       ii=0
       do k=1,nglz
          do j=1,ngly
             do i=1,nglx
                ip=intma(i,j,k,e)
                ii=ii+1
                inode(ii)=ip
                do m=1,ndim
                   qq(m,ii) = q_in(m,ip)
                end do
                jac_e(ii)=jac(i,j,k,e)
             end do !i
          end do !j
       end do !k

       !Do Numerical Integration
       do i=1,npts
          ip=inode(i)
          !Gauss-Lobatto Weight and Jacobian
          wq=jac_e(i)

          !Store integrated q
          do m=1,ndim
             q_tmp(m,ip)=q_tmp(m,ip) + wq*qq(m,i)
          end do !m
       end do !i

    end do !e

    q_in = q_tmp

  end subroutine mass_multiply

  !----------------------------------------------------------------------!
  !>@brief This subroutine constructs the GRADIENT of a 1-vector Q
  !>@author  Francis X. Giraldo on 8/2006
  !>Extended to 3D by Francis X. Giraldo on 11/2009
  !>Modified to use MXM calls by F.X. Giraldo on 12/2013
  !>                  Department of Applied Mathematics
  !>                  Naval Postgraduate School
  !>                  Monterey, CA 93943-5216
  !----------------------------------------------------------------------!
  subroutine compute_gradient(grad_q,q,imass,penalize)

    use mod_mpi_utilities, only: irank0, irank

    implicit none

    !global arrays
    real, intent(out) :: grad_q(3,npoin)
    real, intent(in)  :: q(npoin)
    integer, intent(in)  :: imass
    logical, intent(in), optional :: penalize

    !local arrays
    real, dimension(npts) :: p, p_e, p_n, p_c
    real, dimension(3,num_send_recv_total) :: grad_q_nbh
    integer :: inode(npts)
    real :: wq
    real, dimension(npts) :: e_x, e_y, e_z
    real, dimension(npts) :: n_x, n_y, n_z
    real, dimension(npts) :: c_x, c_y, c_z, jac_e
    real :: p_x, p_y, p_z
    integer :: i, j, k, e, ip, ii
    integer :: ndim
    logical :: lpenalize = .false.

    !print*,'In COMPUTE_GRADIENT: irank lcontinuous = ',irank,lcontinuous_grad_q

    !penalize jump
!    pen = .false.
    if(present(penalize)) lpenalize = penalize

    !Constants for MXM
    ndim=1

    !initialize the global matrix
    grad_q=0.0

    !loop thru the elements
    do e=1,nelem

       !Store Element Variables
       ii=0
       do k=1,nglz
          do j=1,ngly
             do i=1,nglx
                ip=intma(i,j,k,e)
                ii=ii+1
                inode(ii)=ip
                p(ii)  = q(ip)

                !Metrics
                e_x(ii)=ksi_x(i,j,k,e);  e_y(ii)=ksi_y(i,j,k,e);  e_z(ii)=ksi_z(i,j,k,e)
                n_x(ii)=eta_x(i,j,k,e);  n_y(ii)=eta_y(i,j,k,e);  n_z(ii)=eta_z(i,j,k,e)
                c_x(ii)=zeta_x(i,j,k,e); c_y(ii)=zeta_y(i,j,k,e); c_z(ii)=zeta_z(i,j,k,e)
                jac_e(ii)=jac(i,j,k,e)
             end do
          end do
       end do

       !Construct Local Derivatives for Prognostics Variables
       call compute_local_gradient_v3(p_e,p_n,p_c,p,nglx,ngly,nglz,ndim)

       !Loop through I (1-Vector of LGL Points)
       do i=1,npts
          ip=inode(i)

          !Gauss-Lobatto Weight and Jacobian
          wq=jac_e(i)

          !Construct Derivatives in Physical Space (via the Chain Rule)
          p_x=p_e(i)*e_x(i) + p_n(i)*n_x(i) + p_c(i)*c_x(i)
          p_y=p_e(i)*e_y(i) + p_n(i)*n_y(i) + p_c(i)*c_y(i)
          p_z=p_e(i)*e_z(i) + p_n(i)*n_z(i) + p_c(i)*c_z(i)

          !Store RHS
          grad_q(1,ip) = grad_q(1,ip) + wq*p_x
          grad_q(2,ip) = grad_q(2,ip) + wq*p_y
          grad_q(3,ip) = grad_q(3,ip) + wq*p_z
       end do !i
    end do !e

    !compute fluxes
    call compute_gradient_penalty_flux(grad_q,q,lpenalize) !in CREATE_JUMP_PENALTY

    !Apply DSS
    call create_global_rhs(grad_q,grad_q_nbh,3,imass)

  end subroutine compute_gradient

  !----------------------------------------------------------------------!
  subroutine compute_vertical_derivative(dqdz, q, imass,factor)
    
    implicit none

    !global arrays
    real, intent(out) :: dqdz(npoin)
    real, intent(in)  :: q(npoin), factor(npoin)
    integer, intent(in)  :: imass
    !logical, intent(in), optional :: penalize

    !local arrays
    real, dimension(npts) :: p, p_e, p_n, p_c
    integer :: inode(npts)
    real :: wq
    real, dimension(npts) :: e_x, e_y, e_z
    real, dimension(npts) :: n_x, n_y, n_z
    real, dimension(npts) :: c_x, c_y, c_z, jac_e
    real :: p_x, p_y, p_z
    integer :: i, j, k, e, ip, ii
    integer :: ndim
    logical :: pen

    real, dimension(1,num_send_recv_total) :: dqdz_nbh

    !penalize jump, include this back if necessary

    !Constants for MXM
    ndim=1

    !initialize the global matrix
    dqdz=0.0

    !loop thru the elements
    do e=1,nelem

       !Store Element Variables
       ii=0
       do k=1,nglz
          do j=1,ngly
             do i=1,nglx
                ip=intma(i,j,k,e)
                ii=ii+1
                inode(ii)=ip
                p(ii)  = q(ip)

                !Metrics
                e_x(ii)=ksi_x(i,j,k,e)
                e_y(ii)=ksi_y(i,j,k,e)
                e_z(ii)=ksi_z(i,j,k,e)
                n_x(ii)=eta_x(i,j,k,e)
                n_y(ii)=eta_y(i,j,k,e)
                n_z(ii)=eta_z(i,j,k,e)
                c_x(ii)= zeta_x(i,j,k,e)
                c_y(ii)= zeta_y(i,j,k,e)
                c_z(ii)=zeta_z(i,j,k,e)
                jac_e(ii)=jac(i,j,k,e)
             end do
          end do
       end do

       !Construct Local Derivatives for Prognostics Variables
       call compute_local_gradient_v3(p_e,p_n,p_c,p,nglx,ngly,nglz,ndim)

       !Loop through I (1-Vector of LGL Points)
       do i=1,npts
          ip=inode(i)

          !Gauss-Lobatto Weight and Jacobian
          wq=jac_e(i)
          
          !Construct Derivatives in Physical Space (via the Chain Rule)
          p_z=p_c(i)*c_z(i) + p_n(i)*n_z(i) + p_e(i)*e_z(i)
          !Store RHS
          !grad_q(1,ip) = grad_q(1,ip) + wq*p_x
          !grad_q(2,ip) = grad_q(2,ip) + wq*p_y
         dqdz(ip)  = dqdz(ip) + wq*p_z/factor(ip)
       enddo !i
    enddo !e
    call create_global_rhs(dqdz,dqdz_nbh,1,imass)
  end subroutine compute_vertical_derivative
  
  !----------------------------------------------------------------------!
  subroutine compute_vertical_derivative_kessler_contravariant(dqdz, q, v ,imass,factor)

    implicit none

    !global arrays
    real, intent(out) :: dqdz(npoin)
    real, intent(in)  :: q(npoin), factor(npoin), v(npoin)
    integer, intent(in)  :: imass
    !logical, intent(in), optional :: penalize

    !local arrays
    real, dimension(npts) :: p, p_e, p_n, p_c
    real, dimension(npts) :: u, u_c, pu, pu_c, r, r_c, qr, qr_c
    integer :: inode(npts)
    real :: wq, x, y, z, radius
    real, dimension(npts) :: e_x, e_y, e_z
    real, dimension(npts) :: n_x, n_y, n_z, kx, ky, kz
    real, dimension(npts) :: c_x, c_y, c_z, jac_e, xjac_e, ijac_e
    real :: p_x, p_y, p_z, uc
    integer :: i, j, k, e, ip, ii
    integer :: ndim
    logical :: pen

    real, dimension(1,num_send_recv_total) :: dqdz_nbh

    !penalize jump, include this back if necessary

    !Constants for MXM
    ndim=1

    !initialize the global matrix
    dqdz=0.0

    !loop thru the elements
    do e=1,nelem

       !Store Element Variables
       ii=0
       do k=1,nglz
          do j=1,ngly
             do i=1,nglx
                ip=intma(i,j,k,e)
                ii=ii+1
                inode(ii)=ip
                !Metrics
                c_x(ii)=zeta_x(i,j,k,e)
                c_y(ii)=zeta_y(i,j,k,e)
                c_z(ii)=zeta_z(i,j,k,e)
                jac_e(ii)=jac(i,j,k,e)
                xjac_e(ii) = xjac(i,j,k,e)
                ijac_e(ii) = 1.0/xjac_e(ii)
                p(ii) = q(ip) !* v(ip)
                if (geometry_type(1:4) == 'cube') then
                  kx = 0; ky = 0; kz = 1
                else
                  x=coord(1,ip); y=coord(2,ip); z=coord(3,ip)
                  radius=sqrt( dot_product(coord(:,ip),coord(:,ip)) )
                  kx(ii)=x/radius; ky(ii)=y/radius; kz(ii)=z/radius
                endif
                !kx(ii) = kvector(1,ip)
                !ky(ii) = kvector(2,ip)
                !kz(ii) = kvector(3,ip)
                r(ii) = factor(ip)
                qr(ii) = p(ii) / r(ii)

                pu(ii) =  xjac_e(ii)*p(ii)*(c_x(ii) * kx(ii) * v(ip) + c_y(ii) *ky(ii) * v(ip) + c_z(ii) * kz(ii) * v(ip))
                u(ii) = xjac_e(ii)*(c_x(ii) * kx(ii) * v(ip) + c_y(ii) *ky(ii) * v(ip) + c_z(ii) * kz(ii) * v(ip))
             end do
          end do
        end do
               !Construct Local Derivatives for Prognostics Variables
       call compute_local_gradient_v3_zeta(p_c,p,nglx,ngly,nglz,ndim)
       call compute_local_gradient_v3_zeta(u_c,u,nglx,ngly,nglz,ndim)
       call compute_local_gradient_v3_zeta(pu_c,pu,nglx,ngly,nglz,ndim)
       call compute_local_gradient_v3_zeta(r_c,r,nglx,ngly,nglz,ndim)
       call compute_local_gradient_v3_zeta(qr_c,qr,nglx,ngly,nglz,ndim)

       !Loop through I (1-Vector of LGL Points)
       do i=1,npts
          ip=inode(i)

          !Gauss-Lobatto Weight and Jacobian
          wq=jac_e(i)
          uc = ijac_e(i) * u(i)
          !Construct Derivatives in Physical Space (via the Chain Rule)
          p_z = pu_c(i) * ijac_e(i)
          !p_z = uc * (qr(i) * r_c(i) + r(i) * qr_c(i) + p_c(i)) + p(i) * u_c(i) * ijac_e(i) + pu_c(i) * ijac_e(i)
          !p_z= uc * p_c(i) + p(i) * u_c(i) * ijac_e(i) + pu_c(i) * ijac_e(i) !p_e(i)*e_z(i) + p_n(i)*n_z(i) + p_c(i)*c_z(i)

          !Store RHS
          !grad_q(1,ip) = grad_q(1,ip) + wq*p_x
          !grad_q(2,ip) = grad_q(2,ip) + wq*p_y
         dqdz(ip)  = wq*p_z/factor(ip)
       enddo !i
    enddo !e
    call create_global_rhs(dqdz,dqdz_nbh,1,imass)
  end subroutine compute_vertical_derivative_kessler_contravariant

      
  !----------------------------------------------------------------------!
  !>@brief This subroutine constructs the GRADIENT of a ndim-vector Q
  !----------------------------------------------------------------------!
  subroutine compute_gradient_q(grad_q,q,ndim,imass)

    implicit none

    !Global Arrays
    real, intent(out) :: grad_q(3,ndim,npoin)
    real, intent(in)  :: q(ndim,npoin)
    integer, intent(in) :: ndim
    integer, intent(in), optional :: imass

    !Local Arrays
    integer :: inode(npts)
    integer :: e, i, j, k, m, ip, ii, imass_l
    real, dimension(ndim,npts) :: qq, qq_e, qq_n, qq_c
    real, dimension(ndim)      :: qq_x, qq_y, qq_z    
    real, dimension(ndim,num_send_recv_total) :: grad_q_nbh

    !Metric Variables and Constants
    real :: e_x, e_y, e_z
    real :: n_x, n_y, n_z
    real :: c_x, c_y, c_z
    real, dimension(npts) :: ksi_x_e, ksi_y_e, ksi_z_e, &
         eta_x_e, eta_y_e, eta_z_e, &
         zeta_x_e,zeta_y_e,zeta_z_e


    imass_l=0
    if(present(imass)) imass_l=imass

    !Initialize
    grad_q=0.        

    !loop thru the elements
    do e=1,nelem

       !Store Element Variables
       ii=0
       do k=1,nglz
          do j=1,ngly
             do i=1,nglx
                ip=intma(i,j,k,e)
                ii=ii + 1
                inode(ii)=ip
                do m=1,ndim
                   qq(m,ii)=q(m,ip)
                end do

                !Metrics
                ksi_x_e(ii)=ksi_x(i,j,k,e); ksi_y_e(ii)=ksi_y(i,j,k,e); ksi_z_e(ii)=ksi_z(i,j,k,e)
                eta_x_e(ii)=eta_x(i,j,k,e); eta_y_e(ii)=eta_y(i,j,k,e); eta_z_e(ii)=eta_z(i,j,k,e)
                zeta_x_e(ii)=zeta_x(i,j,k,e);zeta_y_e(ii)=zeta_y(i,j,k,e);zeta_z_e(ii)=zeta_z(i,j,k,e)
             end do !i
          end do !j
       end do !k

       !Construct Local Derivatives for Prognostics Variables
       call compute_local_gradient_v3(qq_e,qq_n,qq_c,qq,nglx,ngly,nglz,ndim)

       !Add Metric Terms (they are already pre-multiplied by wgl^3*Jac
       do i=1,npts
          ip=inode(i)          

          !Store Metric Terms
          e_x=ksi_x_e(i);  e_y=ksi_y_e(i);  e_z=ksi_z_e(i)
          n_x=eta_x_e(i);  n_y=eta_y_e(i);  n_z=eta_z_e(i)
          c_x=zeta_x_e(i); c_y=zeta_y_e(i); c_z=zeta_z_e(i)

          !Construct Derivatives in Physical Space (via the Chain Rule)
          !----Perturbation Variables------!
          do m=1,ndim
             qq_x(m)=qq_e(m,i)*e_x + qq_n(m,i)*n_x + qq_c(m,i)*c_x
             qq_y(m)=qq_e(m,i)*e_y + qq_n(m,i)*n_y + qq_c(m,i)*c_y
             qq_z(m)=qq_e(m,i)*e_z + qq_n(m,i)*n_z + qq_c(m,i)*c_z 
          end do !m

          !Store Derivatives and Pass them via MOD_VISCOSITY
          do m=1,ndim
             grad_q(1,m,ip)=qq_x(m)
             grad_q(2,m,ip)=qq_y(m)
             grad_q(3,m,ip)=qq_z(m)
          end do
       end do !i

    end do !e

    !Apply DSS
    do m=1,3
       call create_global_rhs(grad_q(m,:,:),grad_q_nbh,ndim,imass_l)
    end do

  end subroutine compute_gradient_q

  !----------------------------------------------------------------------!
  !>@brief This subroutine constructs the canonical GRADIENT of a 1d-vector Q
  !>@author  Francis X. Giraldo on 8/2006
  !>Extended to 3D by Francis X. Giraldo on 11/2009
  !>Modified to use MXM calls by F.X. Giraldo on 12/2013
  !>Modified from compute_gradient to compute the canonical gradient,
  !>   usefull for contravariant form of equations. F.A.V.B. Alves 03/2020
  !>                  Department of Applied Mathematics
  !>                  Naval Postgraduate School
  !>                  Monterey, CA 93943-5216
  !----------------------------------------------------------------------!
  subroutine compute_canonical_gradient(grad_q_e,grad_q_n,grad_q_c,q,ndim)
    use mod_grid, only: npoin_dg, intma_dg

    implicit none

    !global arrays
    integer, intent(in) :: ndim
    real, dimension(ndim,npoin_dg), intent(out) :: grad_q_e,grad_q_n,grad_q_c
    real, intent(in)  :: q(ndim,npoin)

    !local arrays
    real, dimension(ndim,npts) :: p, p_e, p_n, p_c
    real :: jac_e(npts), wq
    integer :: inode(npts)
    integer :: i, j, k, e, ip, ii, m

    grad_q_e = 0
    grad_q_n = 0
    grad_q_c = 0

    !Constants for MXM
    !loop thru the elements
    do e=1,nelem

       !Store Element Variables
       ii=0
       do k=1,nglz
          do j=1,ngly
             do i=1,nglx
                ip=intma(i,j,k,e)
                ii=ii+1
                inode(ii)=intma_dg(i,j,k,e)
                jac_e(ii)=jac(i,j,k,e)
                do m=1,ndim
                   p(m,ii)  = q(m,ip)
                end do
             end do
          end do
       end do

       !Construct Local Derivatives for Prognostics Variables
       call compute_local_gradient_v3(p_e,p_n,p_c,p,nglx,ngly,nglz,ndim)

       !Loop through I (1-Vector of LGL Points)
       do i=1,npts
          ip=inode(i)
          !Store RHS
          do m=1,ndim
             grad_q_e(m,ip) = p_e(m,i)
             grad_q_n(m,ip) = p_n(m,i)
             grad_q_c(m,ip) = p_c(m,i)
          end do
       end do !i
    end do !e

  end subroutine compute_canonical_gradient

  !----------------------------------------------------------------------!
  !>@brief This subroutine constructs the Divergence of a 3-vector Q
  !>@author  Francis X. Giraldo on 8/2006
  !>Extended to 3D by Francis X. Giraldo on 11/2009
  !>Rewritten to use MXM routines by F.X. Giraldo on 12/2013
  !>                  Department of Applied Mathematics
  !>                  Naval Postgraduate School
  !>                  Monterey, CA 93943-5216
  !----------------------------------------------------------------------!
  subroutine compute_divergence(div_q,q,imass,penalize)

    implicit none

    !global arrays
    real, intent(out) :: div_q(npoin)
    real, intent(in)  :: q(3,npoin)
    integer, intent(in)  :: imass
    logical, intent(in), optional :: penalize

    !local arrays
    real, dimension(3,npts) :: qq, qq_e, qq_n, qq_c
    real, dimension(1,num_send_recv_total) :: div_q_nbh
    integer :: inode(npts)
    real :: wq
    real, dimension(npts) :: e_x, e_y, e_z
    real, dimension(npts) :: n_x, n_y, n_z
    real, dimension(npts) :: c_x, c_y, c_z, jac_e
    real :: u_x, v_y, w_z
    integer :: i, j, k, e, ip, ii
    integer :: ndim
    logical pen

    !penalize jump
    pen = .false.
    if(present(penalize)) pen = penalize

    !Constants for MXM
    ndim=3

    !initialize the global matrix
    div_q=0

    !loop thru the elements
    do e=1,nelem

       !Store Element Variables
       ii=0
       do k=1,nglz
          do j=1,ngly
             do i=1,nglx
                ip=intma(i,j,k,e)
                ii=ii+1
                inode(ii)=ip
                qq(1,ii)= q(1,ip)
                qq(2,ii)= q(2,ip)
                qq(3,ii)= q(3,ip)

                !Metrics
                e_x(ii)=ksi_x(i,j,k,e);  e_y(ii)=ksi_y(i,j,k,e);  e_z(ii)=ksi_z(i,j,k,e)
                n_x(ii)=eta_x(i,j,k,e);  n_y(ii)=eta_y(i,j,k,e);  n_z(ii)=eta_z(i,j,k,e)
                c_x(ii)=zeta_x(i,j,k,e); c_y(ii)=zeta_y(i,j,k,e); c_z(ii)=zeta_z(i,j,k,e)
                jac_e(ii)=jac(i,j,k,e)
             end do
          end do
       end do

       !Construct Local Derivatives for Prognostics Variables
       call compute_local_gradient_v3(qq_e,qq_n,qq_c,qq,nglx,ngly,nglz,ndim)

       !Loop through I (1-Vector of LGL Points)
       do i=1,npts
          ip=inode(i)

          !Gauss-Lobatto Weight and Jacobian
          wq=jac_e(i)

          !Construct Derivatives in Physical Space (via the Chain Rule)
          u_x=qq_e(1,i)*e_x(i) + qq_n(1,i)*n_x(i) + qq_c(1,i)*c_x(i)
          v_y=qq_e(2,i)*e_y(i) + qq_n(2,i)*n_y(i) + qq_c(2,i)*c_y(i)
          w_z=qq_e(3,i)*e_z(i) + qq_n(3,i)*n_z(i) + qq_c(3,i)*c_z(i)

          !Store RHS
          div_q(ip) = div_q(ip) + wq*(u_x + v_y + w_z)
       end do !i

    end do !e

    !compute fluxes
    call compute_divergence_penalty_flux(div_q,q,pen)

    !Apply DSS
    if(imass.ne.2) then !if we need to dss
       call create_global_rhs(div_q,div_q_nbh,1,imass)
    end if

  end subroutine compute_divergence

  !----------------------------------------------------------------------!
  !>@brief This subroutine computes the Curl of the 3-vector Q
  !>@author  Francis X. Giraldo on 9/9/2011
  !>Rewritten to use MXM routines by F.X. Giraldo on 12/2013
  !>           Department of Applied Mathematics
  !>           Naval Postgraduate School
  !>           Monterey, CA 93943-5216
  !----------------------------------------------------------------------!
  subroutine compute_curl(curl_q,q,imass)

    implicit none

    !global arrays
    real, intent(out) :: curl_q(3,npoin)
    real, intent(in)  :: q(3,npoin)
    integer, intent(in)  :: imass

    !local arrays
    real, dimension(3,npts) :: qq, qq_e, qq_n, qq_c
    integer :: inode(npts)
    real :: wq
    real, dimension(npts) :: e_x, e_y, e_z
    real, dimension(npts) :: n_x, n_y, n_z
    real, dimension(npts) :: c_x, c_y, c_z, jac_e
    real :: u_y, u_z
    real :: v_x, v_z
    real :: w_x, w_y
    integer :: i, j, k, m, e, ip, ii
    integer :: m1, m2, m3, ndim

    real, dimension(3,num_send_recv_total) :: curl_q_nbh

    !initialize the global matrix
    curl_q=0
    ndim=3

    !loop thru the elements
    do e=1,nelem

       !Store Element Variables
       ii=0
       do k=1,nglz
          do j=1,ngly
             do i=1,nglx
                ip=intma(i,j,k,e)
                ii=ii+1
                inode(ii)=ip
                qq(1,ii)= q(1,ip)
                qq(2,ii)= q(2,ip)
                qq(3,ii)= q(3,ip)

                !Metrics
                e_x(ii)=ksi_x(i,j,k,e);  e_y(ii)=ksi_y(i,j,k,e);  e_z(ii)=ksi_z(i,j,k,e)
                n_x(ii)=eta_x(i,j,k,e);  n_y(ii)=eta_y(i,j,k,e);  n_z(ii)=eta_z(i,j,k,e)
                c_x(ii)=zeta_x(i,j,k,e); c_y(ii)=zeta_y(i,j,k,e); c_z(ii)=zeta_z(i,j,k,e)
                jac_e(ii)=jac(i,j,k,e)
             end do
          end do
       end do

       !Construct Local Derivatives for Prognostics Variables
       call compute_local_gradient_v3(qq_e,qq_n,qq_c,qq,nglx,ngly,nglz,ndim)

       !Loop through I (and L, LGL Points)
       do i=1,npts
          ip=inode(i)

          !Gauss-Lobatto Weight and Jacobian
          wq=jac_e(i)

          !Construct Derivatives in Physical Space (via the Chain Rule)
          u_y=qq_e(1,i)*e_y(i) + qq_n(1,i)*n_y(i) + qq_c(1,i)*c_y(i)
          u_z=qq_e(1,i)*e_z(i) + qq_n(1,i)*n_z(i) + qq_c(1,i)*c_z(i)
          v_x=qq_e(2,i)*e_x(i) + qq_n(2,i)*n_x(i) + qq_c(2,i)*c_x(i)
          v_z=qq_e(2,i)*e_z(i) + qq_n(2,i)*n_z(i) + qq_c(2,i)*c_z(i)
          w_x=qq_e(3,i)*e_x(i) + qq_n(3,i)*n_x(i) + qq_c(3,i)*c_x(i)
          w_y=qq_e(3,i)*e_y(i) + qq_n(3,i)*n_y(i) + qq_c(3,i)*c_y(i)

          !Store the solution
          curl_q(1,ip) = curl_q(1,ip) + wq*(w_y - v_z)
          curl_q(2,ip) = curl_q(2,ip) + wq*(u_z - w_x)
          curl_q(3,ip) = curl_q(3,ip) + wq*(v_x - u_y)
       end do !i
    end do !e

    !Apply DSS
    call create_global_rhs(curl_q,curl_q_nbh,3,imass)

  end subroutine compute_curl

  !----------------------------------------------------------------------!
  !>@brief This subroutine computes the Vorticity in the Radial Direction of a 3-vector Q
  !>@author  Francis X. Giraldo on 9/9/2011
  !>           Department of Applied Mathematics
  !>           Naval Postgraduate School
  !>           Monterey, CA 93943-5216
  !----------------------------------------------------------------------!
  subroutine compute_vorticity(vorticity_q,q,imass)

    implicit none

    !global arrays
    real, intent(out) :: vorticity_q(3,npoin)
    real, intent(in)  :: q(3,npoin)
    integer, intent(in)  :: imass

    !local arrays
    real :: u(ngl,ngl,ngl)
    real :: v(ngl,ngl,ngl)
    real :: w(ngl,ngl,ngl)
    integer :: inode(ngl,ngl,ngl)
    real :: wq, e_x, e_y, e_z, n_x, n_y, n_z, c_x, c_y, c_z
    real :: u_e, u_n, u_c
    real :: v_e, v_n, v_c
    real :: w_e, w_n, w_c
    real :: h_e, h_n, h_c
    real :: u_x, u_y, u_z
    real :: v_x, v_y, v_z
    real :: w_x, w_y, w_z
    integer :: i, j, k, l, e, ip
    integer :: AllocateStatus

    real, dimension(3,num_send_recv_total) :: vorticity_q_nbh

    !initialize the global matrix
    vorticity_q=0

    !loop thru the elements
    do e=1,nelem

       !Store Element Variables
       do k=1,ngl
          do j=1,ngl
             do i=1,ngl
                ip=intma(i,j,k,e)
                inode(i,j,k)=ip
                ip=inode(i,j,k)
                u(i,j,k)  = q(1,ip)
                v(i,j,k)  = q(2,ip)
                w(i,j,k)  = q(3,ip)
             end do
          end do
       end do

       !Loop through I (and L, LGL Points)
       do k=1,ngl
          do j=1,ngl
             do i=1,ngl
                ip=inode(i,j,k)

                !Gauss-Lobatto Weight and Jacobian
                wq=jac(i,j,k,e)

                !Store Metric Terms
                e_x=ksi_x(i,j,k,e);  e_y=ksi_y(i,j,k,e);  e_z=ksi_z(i,j,k,e)
                n_x=eta_x(i,j,k,e);  n_y=eta_y(i,j,k,e);  n_z=eta_z(i,j,k,e)
                c_x=zeta_x(i,j,k,e); c_y=zeta_y(i,j,k,e); c_z=zeta_z(i,j,k,e)

                !construct Derivatives in Computational Space
                u_e=0; u_n=0; u_c=0
                v_e=0; v_n=0; v_c=0
                w_e=0; w_n=0; w_c=0

                do l = 1,ngl
                   ! Derivatives of Basis functions
                   h_e = dpsi(l,i)
                   h_n = dpsi(l,j)
                   h_c = dpsi(l,k)

                   !KSI Derivatives
                   u_e = u_e + u(l,j,k)*h_e
                   v_e = v_e + v(l,j,k)*h_e
                   w_e = w_e + w(l,j,k)*h_e

                   !ETA Derivatives
                   u_n = u_n + u(i,l,k)*h_n
                   v_n = v_n + v(i,l,k)*h_n
                   w_n = w_n + w(i,l,k)*h_n

                   !ZETA Derivatives
                   u_c = u_c + u(i,j,l)*h_c
                   v_c = v_c + v(i,j,l)*h_c
                   w_c = w_c + w(i,j,l)*h_c
                end do !l

                !Construct Derivatives in Physical Space (via the Chain Rule)
                u_y=u_e*e_y + u_n*n_y + u_c*c_y
                u_z=u_e*e_z + u_n*n_z + u_c*c_z
                v_x=v_e*e_x + v_n*n_x + v_c*c_x
                v_z=v_e*e_z + v_n*n_z + v_c*c_z
                w_x=w_e*e_x + w_n*n_x + w_c*c_x
                w_y=w_e*e_y + w_n*n_y + w_c*c_y

                !Store the solution
                vorticity_q(1,ip) = vorticity_q(1,ip) + wq*(w_y - v_z)
                vorticity_q(2,ip) = vorticity_q(2,ip) + wq*(u_z - w_x)
                vorticity_q(3,ip) = vorticity_q(3,ip) + wq*(v_x - u_y)
             end do !i
          end do !j
       end do !k
    end do !e

    !Apply DSS
    call create_global_rhs(vorticity_q,vorticity_q_nbh,3,imass)

  end subroutine compute_vorticity

  !----------------------------------------------------------------------!
  !>@brief This subroutine constructs the 1D Derivative in the Vertical of a 1-vector Q
  !>@author  Francis X. Giraldo on 8/2006
  !>Extended to 3D by Francis X. Giraldo on 11/2009
  !>                  Department of Applied Mathematics
  !>                  Naval Postgraduate School
  !>                  Monterey, CA 93943-5216
  !----------------------------------------------------------------------!
  subroutine compute_derivative_column(dfdz,f,mvar,ix)

    implicit none

    !global arrays
    real, intent(out) :: dfdz(mvar,nz)
    real, intent(in)  :: f(mvar,nz)
    integer, intent(in) :: mvar, ix

    !local arrays
    real, dimension(nglz,mvar) :: qq, qq_e
    integer :: inode(nglz)
    real :: wq
    integer :: k, m, ip, ez, m1, m2, m3
    logical :: pen

    !initialize the global matrix
    pen=.false.
    dfdz=0

    m1=nglz
    m2=nglz
    m3=mvar

    !Vertical DSS Loop
    do ez=1,nelz

       !Store Element Variables
       do k=1,nglz
          ip=intma_1d(k,ez)
          inode(k)=ip
          do m=1,mvar
             qq(k,m)  = f(m,ip)
          end do !l
       end do !k

       call mxm(dpsiz_tr,m1,qq,m2,qq_e,m3)

       !Loop through I (and L, LGL Points)
       do k=1,nglz
          ip=inode(k)

          !Gauss-Lobatto Weight, Jacobian, and Metric Term
          wq=wglz(k)

          !Store RHS
          do m=1,mvar
             dfdz(m,ip) = dfdz(m,ip) + wq*qq_e(k,m)
          end do
       end do !k
    end do !ez

    !Create penalty flux
    call create_derivative_column_penalty_flux(dfdz,f,mvar,ix,pen)

    !Apply DSS
    call create_global_rhs_1d(dfdz,mvar,ix,1) !DSS with Inverse Mass

  end subroutine compute_derivative_column

  !----------------------------------------------------------------------!
  !>@brief This subroutine computes the canonical derivative of the vector Q as one NxNDIMxNxN system.
  !>That is, it builds Q_e, Q_n, and Q_c.
  !>@author  Francis X. Giraldo on 12/2013
  !>                  Department of Applied Mathematics
  !>                  Naval Postgraduate School
  !>                  Monterey, CA 93943-5216
  !----------------------------------------------------------------------!
  subroutine compute_local_gradient(q_e,q_n,q_c,q,nglx,ngly,nglz,ndim)

    implicit none

    !global arrays
    real, intent(out) :: q_e(nglx,ndim,ngly,nglz)
    real, intent(out) :: q_n(nglx,ndim,ngly,nglz)
    real, intent(out) :: q_c(nglx,ndim,ngly,nglz)
    real, intent(in)  :: q(nglx,ndim,ngly,nglz)
    integer, intent(in) :: nglx, ngly, nglz, ndim

    !local arrays
    integer :: m1, m2, m3, k

    if(nglx > 1) then
       m1=nglx !nglx
       m2=nglx !nglx
       m3=ndim*ngly*nglz !nsize*ngly*nglz
       call mxm(dpsix_tr,m1,q,m2,q_e,m3)
    else
       q_e = 0
    endif

    if(ngly > 1) then
       m1=nglx*ndim !nglx*nsize
       m2=ngly !ngly
       m3=ngly !ngly
       do k=1,nglz !nglz
          call mxm(q(1,1,1,k),m1,dpsiy,m2,q_n(1,1,1,k),m3)
       enddo !k
    else
       q_n = 0
    endif

    if(nglz > 1) then
       m1=nglx*ndim*ngly !nglx*nsize*ngly
       m2=nglz !nglz
       m3=nglz !nglz
       call mxm(q,m1,dpsiz,m2,q_c,m3)
    else
       q_c = 0
    endif

  end subroutine compute_local_gradient

  !----------------------------------------------------------------------!
  !>@brief This subroutine computes the canonical derivative of the vector Q as (N,N,N,NDIM).
  !>That is, it builds Q_e, Q_n, and Q_c.
  !>@author  Francis X. Giraldo on 12/2013
  !>                  Department of Applied Mathematics
  !>                  Naval Postgraduate School
  !>                  Monterey, CA 93943-5216
  !----------------------------------------------------------------------!
  subroutine compute_local_gradient_v2(q_e,q_n,q_c,q,nglx,ngly,nglz,ndim)

    implicit none

    !global arrays
    real, dimension(nglx,ngly,nglz,ndim), intent(out) :: q_e, q_n, q_c
    real, dimension(nglx,ngly,nglz,ndim), intent(in)  :: q
    integer, intent(in) :: nglx, ngly, nglz, ndim

    !local arrays
    integer :: m1, m2, m3, m, k
    real, dimension(nglx,ngly,nglz) :: qt, qt_e, qt_n, qt_c

    !Construct Local Derivatives
    do m=1,ndim
       !Store Input
       qt=q(:,:,:,m)

       !KSI Derivative
       if(nglx > 1) then
          m1=nglx !nglx
          m2=nglx !nglx
          m3=ngly*nglz !nsize*ngly*nglz
          call mxm(dpsix_tr,m1,qt,m2,qt_e,m3)
       else
          qt_e = 0
       endif

       !ETA Derivative
       if(ngly > 1) then
          m1=nglx !nglx
          m2=ngly !ngly
          m3=ngly !ngly
          do k=1,ngl
             call mxm(qt(1,1,k),m1,dpsiy,m2,qt_n(1,1,k),m3)
          enddo !k
       else
          qt_n = 0
       endif

       !Zeta Derivative
       if(nglz > 1) then
          m1=nglx*ngly !nglx*nsize*ngly
          m2=nglz !nglz
          m3=nglz !nglz
          call mxm(qt,m1,dpsiz,m2,qt_c,m3)
       else
          qt_c = 0
       endif

       !Store Output
       q_e(:,:,:,m)=qt_e
       q_n(:,:,:,m)=qt_n
       q_c(:,:,:,m)=qt_c
    end do !m

  end subroutine compute_local_gradient_v2

  !----------------------------------------------------------------------!
  !>@brief This subroutine computes the canonical derivative of the vector Q as (NDIM,N,N,N)
  !>That is, it builds Q_e, Q_n, and Q_c.
  !>This is the workhorse routine in NUMA.
  !>@author  Francis X. Giraldo on 12/2013
  !>                  Department of Applied Mathematics
  !>                  Naval Postgraduate School
  !>                  Monterey, CA 93943-5216
  !----------------------------------------------------------------------!
  subroutine compute_local_gradient_v3(q_e,q_n,q_c,q,nglx,ngly,nglz,ndim)

    implicit none

    !global arrays
    real, dimension(ndim,nglx,ngly,nglz), intent(out) :: q_e, q_n, q_c
    real, dimension(ndim,nglx,ngly,nglz), intent(in)  :: q
    integer, intent(in) :: nglx, ngly, nglz, ndim

    !local arrays
    integer :: m1, m2, m3, m, k
    real, dimension(nglx,ngly,nglz) :: qt, qt_e, qt_n, qt_c

    !Construct Local Derivatives
    do m=1,ndim
       !Store Input
       qt=q(m,:,:,:)

       !KSI Derivative
       if(nglx > 1) then
          m1=nglx !nglx
          m2=nglx !nglx
          m3=ngly*nglz !nsize*ngly*nglz
          call mxm(dpsix_tr,m1,qt,m2,qt_e,m3)
       else
          qt_e = 0
       endif

       !ETA Derivative
       if(ngly > 1) then
          m1=nglx !nglx
          m2=ngly !ngly
          m3=ngly !ngly
          do k=1,nglz
             call mxm(qt(1,1,k),m1,dpsiy,m2,qt_n(1,1,k),m3)
          enddo !k
       else
          qt_n = 0
       endif

       !Zeta Derivative
       if(nglz > 1) then
          m1=nglx*ngly !nglx*nsize*ngly
          m2=nglz !nglz
          m3=nglz !nglz
          call mxm(qt,m1,dpsiz,m2,qt_c,m3)
       else
          qt_c = 0
       endif

       !Store Output
       q_e(m,:,:,:)=qt_e
       q_n(m,:,:,:)=qt_n
       q_c(m,:,:,:)=qt_c
    end do !m

  end subroutine compute_local_gradient_v3

  !----------------------------------------------------------------------!
  !>@brief This subroutine computes the horizontal part of the canonical derivative of the vector Q as (NDIM,N,N,N)
  !>That is, it builds Q_e and Q_n.
  !>@author  Francis X. Giraldo on 12/2013
  !>@author  F.A.V.B. Alves on 03/2020 (modified to perform differentiation along horizontal directions only)
  !>                  Department of Applied Mathematics
  !>                  Naval Postgraduate School
  !>                  Monterey, CA 93943-5216
  !----------------------------------------------------------------------!
  subroutine compute_local_gradient_v3_ksi_eta(q_e,q_n,q,nglx,ngly,nglz,ndim)

    implicit none

    !global arrays
    real, dimension(ndim,nglx,ngly,nglz), intent(out) :: q_e, q_n
    real, dimension(ndim,nglx,ngly,nglz), intent(in)  :: q
    integer, intent(in) :: nglx, ngly, nglz, ndim

    !local arrays
    integer :: m1, m2, m3, m, k
    real, dimension(nglx,ngly,nglz) :: qt, qt_e, qt_n

    !Construct Local Derivatives
    do m=1,ndim
       !Store Input
       qt=q(m,:,:,:)

       !KSI Derivative
       if(nglx > 1) then
          m1=nglx !nglx
          m2=nglx !nglx
          m3=ngly*nglz !nsize*ngly*nglz
          call mxm(dpsix_tr,m1,qt,m2,qt_e,m3)
       else
          qt_e = 0
       endif

       !ETA Derivative
       if(ngly > 1) then
          m1=nglx !nglx
          m2=ngly !ngly
          m3=ngly !ngly
          do k=1,nglz
             call mxm(qt(1,1,k),m1,dpsiy,m2,qt_n(1,1,k),m3)
          enddo !k
       else
          qt_n = 0
       endif

       !Store Output
       q_e(m,:,:,:)=qt_e
       q_n(m,:,:,:)=qt_n
    end do !m

  end subroutine compute_local_gradient_v3_ksi_eta

  !----------------------------------------------------------------------!
  !>@brief This subroutine computes the canonical derivative of the vector Q as (NDIM,N,N,N)
  !>That is, it builds Q_e, Q_n, and Q_c.
  !>@author  Francis X. Giraldo on 3/6/2020
  !>                  Department of Applied Mathematics
  !>                  Naval Postgraduate School
  !>                  Monterey, CA 93943-5216
  !>@date Modified by F.X. Giraldo to only perform differentiation along the Zeta direction
  !----------------------------------------------------------------------!
  subroutine compute_local_gradient_v3_zeta(q_c,q,nglx,ngly,nglz,ndim)

    implicit none

    !global arrays
    real, dimension(ndim,nglx,ngly,nglz), intent(out) :: q_c
    real, dimension(ndim,nglx,ngly,nglz), intent(in)  :: q
    integer, intent(in) :: nglx, ngly, nglz, ndim

    !local arrays
    integer :: m1, m2, m3, m, k
    real, dimension(nglx,ngly,nglz) :: qt, qt_c

    !Construct Local Derivatives
    do m=1,ndim
       !Store Input
       qt=q(m,:,:,:)

       !Zeta Derivative
       if(nglz > 1) then
          m1=nglx*ngly !nglx*nsize*ngly
          m2=nglz !nglz
          m3=nglz !nglz
          call mxm(qt,m1,dpsiz,m2,qt_c,m3)
       else
          qt_c = 0
       endif

       !Store Output
       q_c(m,:,:,:)=qt_c
    end do !m

  end subroutine compute_local_gradient_v3_zeta

  !----------------------------------------------------------------------!
  !>@brief This subroutine computes the inner product of Grad Psi .dot. Grad Q
  !>of the vector Q as (NDIM,N,N,N)
  !>That is, it builds the Laplacian Q_xx, Q_yy, and Q_zz assuming that:
  !>q_e=g1*q_e + g2*q_n + g3_g_c,
  !>q_n=g2*q_e + g4*q_n + g5*q_c,
  !>q_c=g3*q_e + g5*q_n + g6*q_c.
  !>@author  Francis X. Giraldo on 1/2014
  !>                  Department of Applied Mathematics
  !>                  Naval Postgraduate School
  !>                  Monterey, CA 93943-5216
  !----------------------------------------------------------------------!
  subroutine compute_local_gradient_transpose_v3(q,q_e,q_n,q_c,nglx,ngly,nglz,ndim)

    implicit none

    !global arrays
    real, dimension(ndim,nglx,ngly,nglz), intent(out)  :: q
    real, dimension(ndim,nglx,ngly,nglz), intent(in)   :: q_e, q_n, q_c
    integer, intent(in) :: nglx, ngly, nglz, ndim

    !local arrays
    integer :: m1, m2, m3, m, k
    real, dimension(nglx,ngly,nglz) :: qt_e, qt_n, qt_c
    real, dimension(nglx,ngly,nglz) :: q1, q2, q3

    !Construct Local Derivatives
    do m=1,ndim
       !Store Input
       qt_e=q_e(m,:,:,:)
       qt_n=q_n(m,:,:,:)
       qt_c=q_c(m,:,:,:)

       !KSI Derivative
       if(nglx > 1) then
          m1=nglx !nglx
          m2=nglx !nglx
          m3=ngly*nglz !nsize*ngly*nglz
          call mxm(dpsix,m1,qt_e,m2,q1,m3)
       else
          q1 = 0
       endif

       !ETA Derivative
       if(ngly > 1) then
          m1=nglx !nglx
          m2=ngly !ngly
          m3=ngly !ngly
          do k=1,nglz
             call mxm(qt_n(1,1,k),m1,dpsiy_tr,m2,q2(1,1,k),m3)
          enddo !k
       else
          q2 = 0
       endif

       !Zeta Derivative
       if(nglz > 1) then
          m1=nglx*ngly !nglx*nsize*ngly
          m2=nglz !nglz
          m3=nglz !nglz
          call mxm(qt_c,m1,dpsiz_tr,m2,q3,m3)
       else
          q3 = 0
       endif

       !Store Output
       q(m,:,:,:)=q1(:,:,:) + q2(:,:,:) + q3(:,:,:)
    end do !m

  end subroutine compute_local_gradient_transpose_v3

  !----------------------------------------------------------------------!
  !>@brief This subroutine computes the inner product of Grad Psi .dot. Grad Q
  !>of the vector Q as (NDIM,N,N,N)
  !>That is, it builds the Laplacian Q_xx, Q_yy, and Q_zz assuming that:
  !>q_e=g1*q_e + g2*q_n + g3_g_c,
  !>q_n=g2*q_e + g4*q_n + g5*q_c,
  !>q_c=g3*q_e + g5*q_n + g6*q_c.
  !>@author  Francis X. Giraldo on 1/2014
  !>                  Department of Applied Mathematics
  !>                  Naval Postgraduate School
  !>                  Monterey, CA 93943-5216
  !----------------------------------------------------------------------!
  subroutine compute_local_gradient_transpose_v3_ksi_eta(q,q_e,q_n,nglx,ngly,nglz,ndim)

    implicit none

    !global arrays
    real, dimension(ndim,nglx,ngly,nglz), intent(out)  :: q
    real, dimension(ndim,nglx,ngly,nglz), intent(in)   :: q_e, q_n
    integer, intent(in) :: nglx, ngly, nglz, ndim

    !local arrays
    integer :: m1, m2, m3, m, k
    real, dimension(nglx,ngly,nglz) :: qt_e, qt_n
    real, dimension(nglx,ngly,nglz) :: q1, q2, q3

    !Construct Local Derivatives
    do m=1,ndim
       !Store Input
       qt_e=q_e(m,:,:,:)
       qt_n=q_n(m,:,:,:)

       !KSI Derivative
       if(nglx > 1) then
          m1=nglx !nglx
          m2=nglx !nglx
          m3=ngly*nglz !nsize*ngly*nglz
          call mxm(dpsix,m1,qt_e,m2,q1,m3)
       else
          q1 = 0
       endif

       !ETA Derivative
       if(ngly > 1) then
          m1=nglx !nglx
          m2=ngly !ngly
          m3=ngly !ngly
          do k=1,nglz
             call mxm(qt_n(1,1,k),m1,dpsiy_tr,m2,q2(1,1,k),m3)
          enddo !k
       else
          q2 = 0
       endif

       !Store Output
       q(m,:,:,:)=q1(:,:,:) + q2(:,:,:)
    end do !m

  end subroutine compute_local_gradient_transpose_v3_ksi_eta

  !----------------------------------------------------------------------!
  !>@brief This subroutine computes the inner product of Grad Psi .dot. Grad Q
  !>of the vector Q as (NDIM,N,N,N)
  !>That is, it builds the Laplacian Q_xx, Q_yy, and Q_zz assuming that:
  !>q_e=g1*q_e + g2*q_n + g3_g_c,
  !>q_n=g2*q_e + g4*q_n + g5*q_c,
  !>q_c=g3*q_e + g5*q_n + g6*q_c.
  !>@author  Francis X. Giraldo on 1/2014
  !>                  Department of Applied Mathematics
  !>                  Naval Postgraduate School
  !>                  Monterey, CA 93943-5216
  !----------------------------------------------------------------------!
  subroutine compute_local_gradient_transpose_v3_zeta(q,q_c,nglx,ngly,nglz,ndim)

    implicit none

    !global arrays
    real, dimension(ndim,nglx,ngly,nglz), intent(out)  :: q
    real, dimension(ndim,nglx,ngly,nglz), intent(in)   :: q_c
    integer, intent(in) :: nglx, ngly, nglz, ndim

    !local arrays
    integer :: m1, m2, m3, m, k
    real, dimension(nglx,ngly,nglz) :: qt_c
    real, dimension(nglx,ngly,nglz) :: q3

    !Construct Local Derivatives
    do m=1,ndim
       !Store Input
       qt_c=q_c(m,:,:,:)

       !Zeta Derivative
       if(nglz > 1) then
          m1=nglx*ngly !nglx*nsize*ngly
          m2=nglz !nglz
          m3=nglz !nglz
          call mxm(qt_c,m1,dpsiz_tr,m2,q3,m3)
       else
          q3 = 0
       endif

       !Store Output
       q(m,:,:,:)=q3(:,:,:)
    end do !m

  end subroutine compute_local_gradient_transpose_v3_zeta

  !----------------------------------------------------------------------!
  !>@brief This subroutine computes the canonical filter of the vector Q as (NDIM,N,N,N)
  !>That is, it builds Q_e, Q_n, and Q_c.
  !>@author  Francis X. Giraldo on 12/2013
  !>                  Department of Applied Mathematics
  !>                  Naval Postgraduate School
  !>                  Monterey, CA 93943-5216
  !----------------------------------------------------------------------!
  subroutine compute_local_gradient_filter_v3(fqf,q,nglx,ngly,nglz,ndim)

    implicit none

    !global arrays
    real, dimension(ndim,nglx,ngly,nglz), intent(out) :: fqf
    real, dimension(ndim,nglx,ngly,nglz), intent(in)  :: q
    integer, intent(in) :: nglx, ngly, nglz, ndim

    !local arrays
    integer :: m1, m2, m3, m, k
    real, dimension(nglx,ngly,nglz) :: qt, qt_i, qt_ij, qt_ijk

    !Construct Local Derivatives
    do m=1,ndim
       !Store Input
       qt=q(m,:,:,:)

       !KSI Derivative
       if(nglx > 1) then
          m1=nglx !nglx
          m2=nglx !nglx
          m3=ngly*nglz !nsize*ngly*nglz
          call mxm(fx,m1,qt,m2,qt_i,m3)
       else
          qt_i = qt
       endif

       !ETA Derivative
       if(ngly > 1) then
          m1=nglx !nglx
          m2=ngly !ngly
          m3=ngly !ngly
          do k=1,nglz
             call mxm(qt_i(1,1,k),m1,fy_tr,m2,qt_ij(1,1,k),m3)
          enddo !k
       else
          qt_ij = qt_i
       endif

       !Zeta Derivative
       if(nglz > 1) then
          m1=nglx*ngly !nglx*nsize*ngly
          m2=nglz !nglz
          m3=nglz !nglz
          call mxm(qt_ij,m1,fz_tr,m2,qt_ijk,m3)
       else
          qt_ijk = qt_ij
       endif

       !Store Output
       fqf(m,:,:,:)=qt_ijk
    end do !m

  end subroutine compute_local_gradient_filter_v3

  !----------------------------------------------------------------------!
  !>@brief Performs the generic mxm matrix multiplication operation
  !>@author  Sohail Reddy 05/29/2020
  !----------------------------------------------------------------------!
  subroutine perform_mxm_operation(fqf,q,fx_in,fy_in,fz_in, nglx,ngly,nglz,ndim)

    implicit none

    !global arrays
    real, dimension(ndim,nglx,ngly,nglz), intent(out) :: fqf
    real, dimension(ndim,nglx,ngly,nglz), intent(in)  :: q
    real, dimension(nglx,nglx), intent(in)            :: fx_in
    real, dimension(ngly,ngly), intent(in)            :: fy_in
    real, dimension(nglz,nglz), intent(in)            :: fz_in

    integer, intent(in) :: nglx, ngly, nglz, ndim

    !local arrays
    integer :: m1, m2, m3, m, k
    real, dimension(nglx,ngly,nglz) :: qt, qt_i, qt_ij, qt_ijk

    !Construct Local Derivatives
    do m=1,ndim
       !Store Input
       qt=q(m,:,:,:)

       !KSI Derivative
       if(nglx > 1) then
          m1=nglx !nglx
          m2=nglx !nglx
          m3=ngly*nglz !nsize*ngly*nglz
          call mxm(fx_in,m1,qt,m2,qt_i,m3)
       else
          qt_i = qt
       endif

       !ETA Derivative
       if(ngly > 1) then
          m1=nglx !nglx
          m2=ngly !ngly
          m3=ngly !ngly
          do k=1,nglz
             call mxm(qt_i(1,1,k),m1,fy_in,m2,qt_ij(1,1,k),m3)
          enddo !k
       else
          qt_ij = qt_i
       endif

       !Zeta Derivative
       if(nglz > 1) then
          m1=nglx*ngly !nglx*nsize*ngly
          m2=nglz !nglz
          m3=nglz !nglz
          call mxm(qt_ij,m1,fz_in,m2,qt_ijk,m3)
       else
          qt_ijk = qt_ij
       endif

       !Store Output
       fqf(m,:,:,:)=qt_ijk
    end do !m

  end subroutine perform_mxm_operation

  !----------------------------------------------------------------------!
  !>@brief This subroutine computes the local curl
  !>@author  Jeremy E Kozdon on 08/2018
  !>         Department of Applied Mathematics
  !>         Naval Postgraduate School
  !>         Monterey, CA 93943-5216
  !----------------------------------------------------------------------!
  subroutine compute_local_curl(CFx, CFy, CFz, Fx, Fy, Fz, nglx, ngly, nglz)

    implicit none

    !global arrays
    real, dimension(nglx, ngly, nglz), intent(out) :: CFx, CFy, CFz
    real, dimension(nglx, ngly, nglz), intent(in)  :: Fx, Fy, Fz
    integer, intent(in) :: nglx, ngly, nglz

    !local arrays
    integer :: m1, m2, m3, k
    real, dimension(nglx,ngly,nglz) :: Fx_y, Fx_z, Fy_x, Fy_z, Fz_x, Fz_y

    !KSI Derivative
    m1=nglx !nglx
    m2=nglx !nglx
    m3=ngly*nglz !nsize*ngly*nglz
    call mxm(dpsix_tr,m1,Fy,m2,Fy_x,m3)
    call mxm(dpsix_tr,m1,Fz,m2,Fz_x,m3)

    !ETA Derivative
    m1=nglx !nglx
    m2=ngly !ngly
    m3=ngly !ngly
    do k=1,nglz
       call mxm(Fx(1,1,k),m1,dpsiy,m2,Fx_y(1,1,k),m3)
       call mxm(Fz(1,1,k),m1,dpsiy,m2,Fz_y(1,1,k),m3)
    enddo !k

    !Zeta Derivative
    m1=nglx*ngly !nglx*nsize*ngly
    m2=nglz !nglz
    m3=nglz !nglz
    call mxm(Fx,m1,dpsiz,m2,Fx_z,m3)
    call mxm(Fy,m1,dpsiz,m2,Fy_z,m3)

    !Store Output
    CFx(:,:,:) = (Fz_y - Fy_z)
    CFy(:,:,:) = (Fx_z - Fz_x)
    CFz(:,:,:) = (Fy_x - Fx_y)

  end subroutine compute_local_curl

  !> @Brief Compute the volume integral of a field and stores result in result
  !> Input is a q(n,npoin) array. Outputs a r(n) array with the volume integral result of each n
  subroutine compute_volume_integral_field(result,q1)

    use mod_basis, only: nglx, ngly, nglz
    use mod_grid, only: intma, nelem
    use mod_metrics, only: jac
    use mod_types, only: r8
    use mod_numeric_utilities, only: sum_pairwise

    real(r8), intent(in) :: q1(:,:)
    real(r8), intent(out) :: result(size(q1,1))

    integer :: e,i,j,k,ip,l
    real(r8) :: wq, result_el(size(q1,1),nglx*ngly*nglz*nelem)

    !loop thru the elements
    result_el = 0.0
    l=0
    do e=1,nelem
       do k=1,nglz
          do j=1,ngly
             do i=1,nglx
                ip=intma(i,j,k,e)
                wq=jac(i,j,k,e)

                l = l+1
                result_el(:,l) = wq*q1(:,ip)
             enddo !i
          enddo !j
       enddo !k
    enddo !element

    result = sum_pairwise(result_el)

  end subroutine compute_volume_integral_field

  !> @Brief Compute the volume integral of a q1*q2 and stores result in 'result'
  subroutine compute_volume_integral_product(result,q1,q2)

    use mod_basis, only: nglx, ngly, nglz
    use mod_grid, only: intma, nelem
    use mod_metrics, only: jac
    use mod_types, only: r8
    use mod_numeric_utilities, only: sum_pairwise

    real(r8), intent(in) :: q1(:,:), q2(size(q1,1),size(q1,2))
    real(r8), intent(out) :: result(size(q1,1))

    integer :: e,i,j,k,ip,l
    real(r8) :: wq
    real(r8) :: result_i(size(q1,1),nglx*ngly*nglz*nelem)

    !loop thru the elements
    result_i = 0
    l = 0
    do e=1,nelem
       do k=1,nglz
          do j=1,ngly
             do i=1,nglx
                ip=intma(i,j,k,e)
                wq=jac(i,j,k,e)

                l = l+1 
                result_i(:,l) = wq*( q1(:,ip) * q2(:,ip) )
             enddo !i
          enddo !j
       enddo !k
    enddo !element

    result = sum_pairwise(result_i)

  end subroutine compute_volume_integral_product

