function sfr_hs(z)

  ! Cosmic SFR density, calculated by using formula derived by
  ! Springel and Hernquist (2003).  I have picked up the formula from
  ! Faucher-Giguere et al. 2008 (eq. 30 and 31).

  use constants 
  implicit none 
  real(kind=prec), intent(in) :: z ! redshift 
  real(kind=prec) :: sfr_hs ! SFR (Msun yr^-1 Mpc^-3)

  real(kind=prec) :: sfr_z0, a, b, chi 

  sfr_z0 = 0.013_prec ! Msun yr^-1 Mpc^-3 
  a = 0.012_prec ! dimensionless 
  b = 0.041_prec ! dimensionless 
  chi = (omega_nr*(1.0_prec+z)**3+omega_lambda)**(1.0_prec/3.0_prec) ! dimensionless 
  
  sfr_hs = sfr_z0 * chi**2/(1.0_prec + a*((chi-1)**3)*exp(b*(chi**(7.0_prec/4.0_prec)))) 

end function sfr_hs

