!----------------------------------------------------------------------!
!>@brief This subroutine Applies the Physics in a column-based fashion calling 
!> the Kessler_Column routine.
!>@see kessler.f by Sasa Gabersek.
!>
!>@date Reorganized in a column-based fashion by Andreas Mueller.
!> 22 March 2013   
!>
!>@date Simone Marras to solve the production equation with spectral elements.
!> 31 August 2013
!>@date FX Giraldo added options to run various Kessler types based on input parameter kessler_type
!> 04 September 2020 (Happy Birthday Dad!!)
!>@date FX Giraldo restructured to handle both CGc and CGd in a similar fashion for kessler_type==default. The idea is
!>to map the CGd space to CGc space, apply Kessler, and map back to CGd space.
!> 15 September 2020
!New non-column based implementation added by Yassine Tissaoui
!> 20 October 2021  
!----------------------------------------------------------------------!
subroutine kessler(q,press,dt)

  use mod_basis, only: nglz
  use mod_grid, only: coord, npoin, ncol, node_column, nz, nz_cg, intma_1d, intma_1d_cg
  use mod_global_grid, only: ncol_g
  use mod_initial, only: nvar, q_ref
  use mod_input, only:  kessler_type, space_method, nelz, nopz, ztop, zbottom
  use mod_kessler, only: rainnc, rainncv
  use mod_ref, only: press_ref

  implicit none

  !global variables
  real, intent(inout) :: q(nvar,npoin)
  real, intent(in) :: press(npoin)
  real, intent(in) :: dt

  !local variables
  real, dimension (:), allocatable :: rho, t, qv, qc, qr, p, z, rnc, rncv, temp1, temp2
!!$  real :: rho(nz_cg), t(nz_cg), qv(nz_cg), qc(nz_cg), qr(nz_cg), p(nz_cg), z(nz_cg)
!!$  real :: rnc(nz_cg), rncv(nz_cg)
  real :: temp3, temp4
  integer :: i, j, k, e, iz, iz_g, icol
  integer :: ierr, irank
  integer :: ip, ip_g, ivar, nz_kessler
  
  if (kessler_type == 'default' .or. kessler_type == 'fdm') then
  !loop over all columns
      do icol=1,ncol

     !Gather vector along the column (C^0)
        nz_kessler=nz_cg
        if (.not. allocated(rho)) then
           allocate(rho(nz_cg), t(nz_cg), qv(nz_cg), qc(nz_cg), qr(nz_cg), p(nz_cg), z(nz_cg), rnc(nz_cg), rncv(nz_cg))
        end if
        if (space_method(1:2) == 'cg') then
           do e=1,nelz
              do k=1,nglz
                 iz_g = intma_1d_cg(k,e)
                 iz   = intma_1d(k,e)
                 ip   = node_column(icol,iz)
                 !Store CGc values
                 rho(iz_g) = q(1,ip) + q_ref(1,ip)     !rho
                 t(iz_g)   = q(5,ip) + q_ref(5,ip)     !theta
                 qv(iz_g)  = q(6,ip) + q_ref(6,ip)     !qv
                 qc(iz_g)  = q(7,ip)                   !qc
                 qr(iz_g)  = q(8,ip)                   !qr
                 p(iz_g)   = press(ip) + press_ref(ip) !p
                 z(iz_g)   = coord(3,ip)
                 rnc(iz_g) = rainnc(ip)
                 rncv(iz_g)= rainncv(ip)              
              end do
           end do
        else
           print*,' Error in Kessler.F90: this space_method not allowed = ',space_method
           stop
        endif
        call kessler_column_Andreas_SEM(rho,t,qv,qc,qr,p,q,z,temp3,temp4,dt,nz_kessler,icol)
        if (space_method(1:2) == 'cg') then
           do e=1,nelz
              do k=1,nglz
                 iz_g = intma_1d_cg(k,e)
                 iz   = intma_1d(k,e)
                 ip   = node_column(icol,iz)

                 !Store Q values
                 q(1,ip) = rho(iz_g) - q_ref(1,ip) !rho'
                 q(5,ip) = t(iz_g)   - q_ref(5,ip) !theta'
                 q(6,ip) = qv(iz_g)  - q_ref(6,ip) !qv'
                 q(7,ip) = qc(iz_g)
                 q(8,ip) = qr(iz_g)

                 !Store accumulated rain
                 rainnc(ip) = temp3
                 rainncv(ip)= temp4
              end do
           end do
         else
           print*,' Error in Kessler.F90: this space_method not allowed = ',space_method
           stop
         endif
       enddo
     elseif (kessler_type == 'sem') then
       do icol=1,ncol

     !Gather vector along the column (C^0)
        nz_kessler=nz
        if (.not. allocated(rho)) then
           allocate(rho(nz), t(nz), qv(nz), qc(nz), qr(nz), p(nz), z(nz), rnc(nz), rncv(nz))
        end if
           do e=1,nelz
              do k=1,nglz
                 iz   = intma_1d(k,e)
                 ip   = node_column(icol,iz)
                 !Store CGc values
                 rho(iz) = q(1,ip) + q_ref(1,ip)     !rho
                 t(iz)   = q(5,ip) + q_ref(5,ip)     !theta
                 qv(iz)  = q(6,ip) + q_ref(6,ip)     !qv
                 qc(iz)  = q(7,ip)                   !qc
                 qr(iz)  = q(8,ip)                   !qr
                 p(iz)   = press(ip) + press_ref(ip) !p
                 z(iz)   = coord(3,ip)
                 rnc(iz) = rainnc(ip)
                 rncv(iz)= rainncv(ip)
              end do
           end do
        call kessler_column_Andreas_SEM(rho,t,qv,qc,qr,p,q,z,temp3,temp4,dt,nz_kessler,icol)
           do e=1,nelz
              do k=1,nglz
                 iz   = intma_1d(k,e)
                 ip   = node_column(icol,iz)

                 !Store Q values
                 q(1,ip) = rho(iz) - q_ref(1,ip) !rho'
                 q(5,ip) = t(iz)   - q_ref(5,ip) !theta'
                 q(6,ip) = qv(iz)  - q_ref(6,ip) !qv'
                 q(7,ip) = qc(iz)
                 q(8,ip) = qr(iz)

                 !Store accumulated rain
                 rainnc(ip) = temp3
                 rainncv(ip)= temp4
              end do
           end do
       enddo
     elseif (kessler_type == 'sem_no_column') then
        nz_kessler=nz
        if (.not. allocated(rho)) then
           allocate(rho(npoin), t(npoin), qv(npoin), qc(npoin), qr(npoin), p(npoin), z(npoin), rnc(npoin), rncv(npoin), temp1(npoin), temp2(npoin))
        end if
            do ip = 1,npoin
              !Store EBG values
              rho(ip) = q(1,ip) + q_ref(1,ip)     !rho
              t(ip)   = q(5,ip) + q_ref(5,ip)     !theta
              qv(ip)  = q(6,ip) + q_ref(6,ip)     !qv
              qc(ip)  = q(7,ip)                   !qc
              qr(ip)  = q(8,ip)                   !qr
              p(ip)   = press(ip) + press_ref(ip) !p
              z(ip)   = coord(3,ip)
              rnc(ip) = rainnc(ip)
              rncv(ip)= rainncv(ip)              
           end do
        temp1(:)=rnc(:)
        temp2(:)=rainncv(:)
        call kessler_Nocolumn_Yassine(rho,t,qv,qc,qr,p,q,z,temp1,temp2,dt,nz_kessler,q_ref)                

            do ip = 1,npoin
              
              !Store Q values
              q(1,ip) = rho(ip) - q_ref(1,ip) !rho'
              q(5,ip) = t(ip)   - q_ref(5,ip) !theta'
              q(6,ip) = qv(ip)  - q_ref(6,ip) !qv'
              q(7,ip) = qc(ip)
              q(8,ip) = qr(ip)
              
              !Store accumulated rain
              rainnc(ip) = temp1(ip)
              rainncv(ip)= temp2(ip)
           end do
     endif
     

end subroutine kessler

!----------------------------------------------------------------------!
!>@brief Performs Kessler Physics without using column data structures: 
!>@author Yassine Tissaoui
!----------------------------------------------------------------------!
subroutine kessler_Nocolumn_Yassine(rho,t,qv,qc,qr,p,q,z,rainnc,rainncv,dt,nz, qref)
  use mod_basis, only: nglx, ngly, nglz
  use mod_constants, only: cp, rgas, p00
  use mod_grid, only: npoin, ncol, node_column, coord, nelem, intma
  use mod_initial, only: nvar
  use mod_input, only: nelx, nely, nelz, kessler_type, ztop, zbottom, nopz

  implicit none

  !global variables
  real, intent(inout) :: rho(npoin), t(npoin), qv(npoin), qc(npoin), qr(npoin), p(npoin), z(npoin)
  real, intent(inout) :: rainnc(npoin), rainncv(npoin)
  real, intent(in) :: q(nvar,npoin), qref(nvar,npoin)
  real, intent(in) :: dt
  integer, intent(in) :: nz

  !local variables
  real :: xlv,ep2,svp1,svp2,svp3,svpt0,rhowater

  !local variables from the original module_mp_kessler.F
  real :: qrprod, ern, gam, rcgs, rcgsi
  real :: prod(npoin), rhs(npoin)
  real :: vt(npoin), prodk(npoin), vtden(npoin),rdzk(npoin),rhok(npoin)
  real :: rdzw(npoin)
  integer :: nfall, n, nfall_new, icol
  real    :: qrr, pressure, temp, es, qvs, dz, pii
  real    :: f5, dtfall, rdz, product
  real    :: max_heating, max_condense, max_rain, maxqrp
  real    :: vtmax, ernmax, crmax, factorn, time_sediment
  real    :: qcr, factorr, ppt
  real, parameter :: max_cr_sedimentation = 0.75

  !local variables
  real    :: q0(npoin), qp(npoin)
  integer :: i, j, k, ie, ip, ip1, e
  integer :: icount, jcount

  !----------------------------------------------------------------
  real , parameter ::  c1 = .001
  real , parameter ::  c2 = .001
  real , parameter ::  c3 = 2.2
  real , parameter ::  c4 = .875
  real , parameter ::  fudge = 1.0
  real , parameter ::  mxfall = 10.0
  !----------------------------------------------------------------
  xlv =     2500000.
  ep2 =         0.6217504

  ! constants of Magnus-Tetens formula for saturation water pressure:
  ! (see Klemp and Wilhelmson (1978) eq. (2.11),
  ! Emanuel (textbook Atmospheric Convection, 1994) eq. 4.4.14)
  svp1 =        0.6112000
  svp2 =         17.67000
  svp3 =         29.65000
  svpt0 =        273.1500
  rhowater =     1000.000

  !  input arrays
  !
  !  t - potential temperature
  !  qv, qc, qr  - mixing ratio (g/g dry air) of water vapor, cloud water
  !                and rain water
  !  pii         - exner function
  !  dt_in - timestep
  !  z  - height of (t,qv,p,rho) points in meters
  !  dz8w - delta z between t points.
  !
  !  See Klemp and Wilhelmson (1978) Journal of the Atmospheric Sciences
  !  Vol 35, pp 1070-1096 for more details
  !
  !  some more explanations (Andreas Mueller):
  !  es: saturation water pressure according to Magnus-Tetens
  !  qvs: saturation mixing ratio
  !  prod, product: production of qc
  !

  !   f5 = 237.3 * 17.27 * 2.5e6 / cp
  f5 = svp2*(svpt0-svp3)*xlv/cp
  ernmax = 0.
  maxqrp = -100.

  !------------------------------------------------------------------------------
  ! parameters for the time split terminal advection
  !------------------------------------------------------------------------------
  max_heating = 0.
  max_condense = 0.
  max_rain = 0.

  ! are all the variables ready for the microphysics?
  ! start the microphysics
  ! do the sedimentation first
  crmax = 0.

  !------------------------------------------------------------------------------
  ! Terminal velocity calculation and advection, set up coefficients and
  ! compute stable timestep
  !------------------------------------------------------------------------------
  !do e = 1, nelem
   ! do i = 1, nglx
    !  do j = 1, ngly
     !   do k = 1, nglz
      !    ip = intma(i,j,k,e)
       !   if (k .eq. nglz) then
        !    ip1 = intma(i,j,k-1,e)
         !   rdzk(ip) = 1 / (coord(3,ip)-coord(3,ip1))
          !else
           ! ip1 = intma(i,j,k+1,e)
            !rdzk(ip) = 1 / (coord(3,ip1)-coord(3,ip))
          !end if
        !enddo
      !enddo
    !enddo
  !enddo
   do ip = 1,npoin
      rdzk(ip)=1 / ((ztop - zbottom)/max(nelz*nopz,1)) !FIX ME to be more general
      prodk(ip)   = qr(ip)
      rhok(ip) = rho(ip)
      qrr = max(0.,qr(ip)*0.001*rhok(ip))
      vtden(ip) = sqrt(1.123 / rhok(ip))
      vt(ip) = 36.34*(qrr**0.1364) * vtden(ip)
      ! vt: terminal fall velocity (Klemp and Wilhelmson (1978) eq. (2.15))
      crmax = max(vt(ip)*dt*rdzk(ip),crmax)
    enddo
   
  nfall = max(1,nint(0.5+crmax/max_cr_sedimentation)) ! courant number for big timestep.
  dtfall = dt / float(nfall) ! splitting so courant number for sedimentation is stable
  time_sediment = dt
  !------------------------------------------------------------------------------
  ! Terminal velocity calculation and advection
  ! Do a time split loop on this for stability.
  !------------------------------------------------------------------------------
  !Initialize
  q0    = prodk
  column_sedimentation: do while ( nfall > 0 )

     time_sediment = time_sediment - dtfall
     do ip= 1, npoin
       ppt=0.0
       ppt     = rhok(ip)*prodk(ip)*vt(ip)*dtfall/rhowater
       if (coord(3,ip) < 1.0) then 
        rainncv(ip) = ppt*1000.0
        rainnc(ip)  = rainnc(ip) + ppt*1000.0 ! unit = mm
       endif 
     enddo
     !Currently only allows sem kessler.
     call ti_rk35_production_no_column(rhs,prodk,q0,dtfall,rhok,vt)
     do ip=1, npoin
     enddo
      !------------------------------------------------------------------------------
     ! compute new sedimentation velocity, and check/recompute new
     ! sedimentation timestep if this isn't the last split step.
     !------------------------------------------------------------------------------
     if ( nfall > 1 ) then ! this wasn't the last split sedimentation timestep

        nfall = nfall - 1
        crmax = 0.
        do ip = 1, npoin
            qrr = max(0.,prodk(ip)*0.001*rhok(ip))
            vt(ip) = 36.34 * (qrr**0.1346) * vtden(ip)
           ! vt: terminal fall velocity (Klemp and Wilhelmson (1978) eq. (2.15))
           !          vtmax = max(vt(k), vtmax)
            crmax = max(vt(ip)*time_sediment*rdzw(ip),crmax)  !FXG rdzw(k) Is = 0
        enddo

        nfall_new = max(1,nint(0.5+crmax/max_cr_sedimentation))
        if (nfall_new /= nfall ) then
           nfall = nfall_new
           dtfall = time_sediment/nfall
        end if

     else  ! this was the last timestep

        do ip=1,npoin
           prod(ip) = prodk(ip)
        enddo
        nfall = 0  ! exit condition for sedimentation loop
     endif

     !
     ! Update Time and Solution Vector
     !
     q0=prodk

  enddo column_sedimentation

  ! now the conversion processes

  !------------------------------------------------------------------------------
  ! Production of rain and deletion of qc
  ! Production of qc from supersaturation
  ! Evaporation of QR
  !------------------------------------------------------------------------------
  do ip = 1, npoin
     ! autoconversion and accretion: (Klemp and Wilhelmson (1978) eq. (2.13))
     factorn = 1.0 / (1.+ c3 * dt * max(0., qr(ip))**c4)
     qrprod = qc(ip) * (1.0 - factorn) + factorn * c1 * dt * max(qc(ip) - c2,0.)
     rcgs = 0.001 * rho(ip)
     qc(ip) = max(qc(ip) - qrprod,0.)
     qr(ip) = (qr(ip) + prod(ip) - qr(ip))
     qr(ip) = max(qr(ip) + qrprod,0.)
     pii = (p(ip) / p00)**(rgas / cp)
     temp = t(ip) * pii ! t: potential temperature, temp: temperature
     pressure = p(ip)

     gam = 2.5e+06/(cp*pii) ! see Klemp and Wilhelmson (1978) below eq. (2.9d)
     ! 2.5e+6: latent heat of vaporization L
     !      qvs       = 380.*exp(17.27*(temp-273.)/(temp- 36.))/pressure
     es        = 1000.*svp1*exp(svp2*(temp-svpt0)/(temp-svp3))
     ! es: saturation water pressure according to Magnus-Tetens (see Klemp and Wilhelmson (1978) eq. (2.11),
     ! Emanuel (textbook Atmospheric Convection, 1994) eq. 4.4.14)
     qvs       = ep2*es/(pressure-es) ! saturation mixing ratio
     !      prod(k) = (qv(k)-qvs) / (1.+qvs*f5/(temp-36.)**2)
     prod(ip) = (qv(ip)-qvs) / (1.+pressure/(pressure-es)*qvs*f5/ &
          (temp-svp3)**2) ! production of qc
     !       print*,k,rho(k),t(k),qv(k),qc(k),qr(k),p(k),z(k),pressure,qvs
     ern  = min(dt*(((1.6+124.9*(rcgs*qr(ip))**.2046) &
          *(rcgs*qr(ip))**.525)/(2.55e8/(pressure*qvs) &
          +5.4e5))*(dim(qvs,qv(ip))/(rcgs*qvs)), &
          max(-prod(ip)-qc(ip),0.),qr(ip))
     ! ern: evaporation of rain according to Ogura and Takshashi (1971) (see Klemp and Wilhelmson (1978) eq. (2.14))

     ! Update all variables
     product = max(prod(ip),-qc(ip)) ! deletion of qc cannot be larger than qc
     t(ip) = t(ip) + gam*(product - ern)
     qv(ip)= max(qv(ip) - product + ern,0.)
     qc(ip) =       qc(ip) + product
     qr(ip) = qr(ip) - ern

  enddo

end subroutine kessler_Nocolumn_Yassine

subroutine ti_rk35_production_no_column(rhs,qp,q0,dtfall,rhok,vt)

  use mod_grid, only: npoin
  use mod_initial, only: nvar
  use mod_input, only: kessler_type

  implicit none

  !global
  real :: dtfall
  real :: qp(npoin)
  real :: q0(npoin)
  real :: vt(npoin)
  real :: rhok(npoin)
  real :: rhs(npoin)

  !local
  integer :: ip
  real    :: a0, a1, a2, a3, a4, b0, b1, b2, b3, b4, beta, dtt
  real    :: q1(npoin), q2(npoin)
  integer :: ik, i, j, K

  !Initialize
  q1=qp
  q2=0
  a2=0

  ik = 1
  a0=1
  a1=0
  a2=0
  beta=0.377268915331368

  !dtt=dtfall*beta
  dtt=dtfall
  !OK

  !Compute RHS
  call create_rhs_sem_no_column(rhs,q1,rhok,vt)
  qp(:)= q0(:) + dtt*rhs(:)
  !OK up to here. Take it from here for RK35.

  !Apply Boundary Conditions
!!$  qp(1)  = 0.0
!!$  qp(nz) = 0.0
  !call apply_boundary_conditions(qp)

  !Update
  q1=qp

  return
end subroutine ti_rk35_production_no_column

subroutine create_rhs_sem_no_column(rhs,prod,densi,vt)

  use mod_basis, only: nglz, dpsiz, wglz
  use mod_grid, only: intma_1d,  npoin, nz, ncol
  use mod_input, only: nelz, imass
  use mod_metrics, only: jac_1d
  use mod_interface, only: compute_vertical_derivative, compute_vertical_derivative_kessler_contravariant

  implicit none

  !global
  real, dimension(npoin)         :: rhs
  real, dimension(npoin)         :: prod, densi, vt, f

  !local
  integer                   :: ip
  !Loop through Elements
  rhs(:)=0
  do ip = 1,npoin
    f(ip) = prod(ip)*densi(ip)*vt(ip)
  enddo
  call compute_vertical_derivative(rhs,f,imass,densi)
  !call compute_vertical_derivative_kessler_contravariant(rhs,f,vt,imass,densi)

end subroutine create_rhs_sem_no_column

