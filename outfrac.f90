function outfrac(HaloMass, z, pop) 

  ! File: outfrac.f90 
  !  Cre: 2012-06-21
  !  Mod: $Date: 2012/10/23 08:10:11 $ ($Revision: 1.6 $) 

  use constants 
  use storage 
  use interfaces, only : vesc_sq, getsfr, imf, interpolate2, getmet, &
       &bi_interpolate2, sfr_rollinde_pop3, sfr_rollinde_pop2, imf_pop3, getsfr_hot, getsfr_cold
  implicit none 
  real(kind=prec), intent(in) :: HaloMass, z ! [HaloMass] = 10^10 M_solar 
  integer, intent(in) :: pop 
  real(kind=prec) :: outfrac 

  real(kind=prec) :: t, m_d, m_up, m1, m2, int1, int2, sum, tmyr, vesq, Haloradius
  integer :: i, j 

  call interpolate2(tarr, zarr, z, t) ! [t] = yr 
  tmyr = t*1.0e-6_prec ! Myr 
  call interpolate2(stellar_mass, stellar_age, tmyr, m_d) ! [md] = M_solar
  m_up = STMASS_UPLIMIT ! M_solar 

  do i = 1, 4 
     stmet(i) = yields(i,1,1) 
     do j = 1, 10 
        e_kin(i,j) = yields(i,j,4) ! log(erg)
     end do
  end do

  m1 = max(m_d,8.0_prec)
  int1 = integrand(m1) 
  sum = 0.0_prec 
  do 
     m2 = m1*STELLAR_INTEGRAL_MULTIPLIER
     if (m2 > m_up) exit 

     int2 = integrand(m2) 
     sum = sum + 0.5_prec*(m2-m1)*(int1+int2) ! erg 
     int1 = int2 
     m1 = m2 
  end do

  HaloRadius = ((3.0_prec*HaloMass)/&
       &(4.0_prec*pi*rho_critical*omega_nr*delta_virial))**(1.0_prec/3.0_prec) ! Mpc
  vesq = 2.0_prec*newtg*HaloMass*1.0e12_prec*msolkg*cmbympc/HaloRadius ! (m/s)^2

  ! The factor 1.0e-7 converts from erg to joule.
  outfrac = 2.0_prec*sum*outflow_efficiency*1.0e-7_prec&
       &/(vesq*msolkg) ! dimensionless (because we divide by 1 M_solar)

contains 

  function integrand(m) 

    implicit none 
    real(kind=prec), intent(in) :: m 
    real(kind=prec) :: integrand
    real(kind=prec) :: ekinetic, st_age, rs, psi2, zmet, lzmet, st_ageyr, &
         &psi3, integrand2, integrand3 
    logical :: m_overflow, m_underflow, zmet_overflow, zmet_underflow 

    call interpolate2(stellar_age, stellar_mass, m, st_age) ! [st_age] = Myr
    st_ageyr = st_age*1.0e6_prec ! yr 
    call interpolate2(zarr, tarr, t-st_ageyr, rs) 
    zmet = getmet(rs) 

    m_overflow = .false. 
    m_underflow = .false. 
    if (m > data_mmax) then 
       m_overflow = .true. 
    else if (m < data_mmin) then 
       m_underflow = .true. 
    end if

    zmet_overflow = .false. 
    zmet_underflow = .false. 
    if (zmet > data_zmetmax) then 
       zmet_overflow = .true. 
    else if (zmet < data_zmetmin) then 
       zmet_underflow = .true. 
    end if

    if (m_overflow) then 
       if (zmet_overflow) then 
          call bi_interpolate2(stmet, stmass, e_kin, data_zmetmax, data_mmax, ekinetic)
       else if (zmet_underflow) then 
          call bi_interpolate2(stmet, stmass, e_kin, data_zmetmin, data_mmax, ekinetic)
       else 
          call bi_interpolate2(stmet, stmass, e_kin, zmet, data_mmax, ekinetic)
       end if
    else if (m_underflow) then 
       if (zmet_overflow) then 
          call bi_interpolate2(stmet, stmass, e_kin, data_zmetmax, data_mmin, ekinetic)
       else if (zmet_underflow) then 
          call bi_interpolate2(stmet, stmass, e_kin, data_zmetmin, data_mmin, ekinetic)
       else 
          call bi_interpolate2(stmet, stmass, e_kin, zmet, data_mmin, ekinetic)
       end if
    else
       if (zmet_overflow) then 
          call bi_interpolate2(stmet, stmass, e_kin, data_zmetmax, m, ekinetic)
       else if (zmet_underflow) then 
          call bi_interpolate2(stmet, stmass, e_kin, data_zmetmin, m, ekinetic)
       else 
          call bi_interpolate2(stmet, stmass, e_kin, zmet, m, ekinetic)
       end if
    end if

    if (pop==3) then 
       integrand = imf_pop3(m)*ekinetic ! dimensionless
    else if (pop==2) then 
       if (m < POP2_UPLIMIT) then 
          integrand = imf(m)*ekinetic ! dimensionless
       else 
          integrand = 0.0_prec 
       end if
    else
       write (0,*) 'outfrac.f90: pop has an illegal value!'
    end if


!!$    if (HaloMetalAbundance < METALLICITY_POP3TRANS) then 
!!$       integrand = imf_pop3(m)*ekinetic ! dimensionless
!!$    else
!!$       if (m < POP2_UPLIMIT) then 
!!$          integrand = imf(m)*ekinetic ! dimensionless
!!$       else 
!!$          integrand = 0.0_prec 
!!$       end if
!!$    end if

!!$    if (z > REDSHIFT_POP3TRANS) then 
!!$       integrand = imf_pop3(m)*ekinetic ! dimensionless
!!$    else 
!!$       integrand = imf(m)*ekinetic ! dimensionless
!!$    end if


!!$    integrand2 = imf(m)*ekinetic ! erg M_solar^-1  
!!$    integrand3 = imf_pop3(m)*ekinetic ! erg M_solar^-1  
!!$    ! integrand = integrand2 + integrand3 ! erg M_solar^-1  
!!$    integrand = integrand2 ! erg M_solar^-1  

  end function integrand
  
end function outfrac
