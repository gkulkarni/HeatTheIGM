
PROGRAM REION 

  !  Main file of the thermal history code.  This file belongs to
  !  `thermal', NOT to `reion-eq'.

  USE CONSTANTS 
  USE STORAGE 
  USE INTERFACES, ONLY : CLUMPFAC, DTDZ, GAMMA_PH, GAMMA_PI, HUBP, &
       &IGMFPREO, IGMVFRAC, JEANS_LENGTH, PSIG, RTBIS, SIGMAH , &
       &SIGMA_BARYON, SOLDLT, SOLFM, LUMFN, gallum, COUNTER, &
       &SOURCEOV_POP2, SOURCEOV_POP3, GPH_KERNEL_POP2, outflow, &
       &GPI_KERNEL_POP2, GPI_KERNEL_POP3, GPH_KERNEL_POP3, hallum, &
       &ACCRATE, ejrate, sfr_rollinde_pop3, sfr_rollinde_pop2, &
       &interpolate2, getsfr2, getsfr3, getjmc, getjmh, haloyield_nonira, &
       &rtnewt, newtondelta, newtonfm, igmfposto, ejfrac_nonira, outfrac_nonira, &
       &haloyield_species_nonira, counter, ngammafrac, haloyield_species, sfr_hs 
  IMPLICIT NONE 

  REAL(KIND=PREC) :: A, B, BOUND, C, DFM, DGRF, DLT, DLTHI, DLTLO, &
       &DSIG, FMHI, FMLO, FV, F_M, GAMMAPH, GAMMA_HEAT, DNLLDZ, &
       &GAMMA_ION, GAMMA_RC, GAMMA_REC, GAMMA_TOTC, GRF, H0T0, HEATCOOL, &
       &HTT, IGMDCRIT, JNSLN, LMFP, M, NGAMMA, NH, NHI, GAMMAPI,&
       &NHII, NH_PROPER, OLDFM, Q, QE, R, R1, R2, RECLAMBDA, RECRATE, &
       &RHO, RHO_BARYON, SIG, SOURCE, SUM, ERROR, T0SEC, T0YR, T1, T2, &
       &TEMPC, TEMPH, TEMPH_EV, TOLZIN, X2INIT, Z, HEAT, COOL, acc, &
       &ZOFT, FREQ, NDOTM, GPH, GPI, TAU1, TAU2, TAU, JNSM, TOTAL_AGE,&
       &TEMPHVA, X_IIVA, OLDXIIVA, LB_AGE, LB_LAMBDA, LB_LUM, InitialTime,&
       &FinalTime, HaloTotalMass, HaloFormationRedshift, t, HaloLuminosity,&
       &VCIRC, DOFZ, SOURCE_POP2, SOURCE_POP3, GPI_POP3, GPI_POP2, ejc, &
       &GPH_POP2, GPH_POP3, ST_MASS, ST_AGE, YIELD_RECORD(25), ofl, &
       &ejr, LOCAL_HUBBLE_0, metallicity, metmass, m_c, m_o, m_fe, &
       &m_h, fe_abundance, c_abundance, o_abundance, source_dummy, &
       &dtrans, ofl_factor, ejrate_fe, ejrate_o, ejrate_c, febyh, &
       &dfebyh, zfe, psi, totalmstar, tfe, tnow, age_fe, mstar_fe, &
       &r_local, mhalo, mhalo_high, mhalo_low, mdot, mstardot, mmetaldot, &
       &return_fraction, fgas_in, mgasdot, mstardot_insitu, Omega_mz, &
       &fout, HaloVirialRadius, d, Delta_c, DiscSpinParameter, &
       &DiscScaleLength, HaloCircularVelocity, tdyn, zeta, mCdot, &
       &mFedot, mOdot, mH_halos, nu, f, grw, sgm, ms, n, dsgm, mo, &
       &log_multiplier, BinsPerDecade, dm, sm, smc, smh, mminc, mminh, &
       &sm_pop3, sm_pop2, smc_pop3, smc_pop2, smh_pop3, smh_pop2, &
       &st_ngamma, metal_term1, metal_term2, metal_term3, term1, &
       &term2, global_t, fb_struct, mNdot, mSidot, mZndot, mMgdot,&
       &ejrate_n, ejrate_si, ejrate_mg, ejrate_zn, foo, decr, mcoolgas,&
       &sdt, sdl, HaloVirialTemp, GasCoolingRate, eta, p, hb, fb, ft,&
       &mcooldot, halo_lum, zlim, maglim, limsfr  

  INTEGER :: I, J, IER, INF, LAST, NEVAL, NUMBER_OF_LINES, NCALC, &
       &comnum

  CHARACTER(100) :: ZLUMARG, RGNOVDARG, FESCARG, RGNSIZEARG, yfooarg 

  LOGICAL :: PREOVERLAP 

  ! ----------------------------------------

  comnum = command_argument_count()
  if (comnum /= 1) then 
     write (0, '(a)') 'Usage: ./thermal FESC'
     write (0, '(a)') 'Example: ./thermal 0.3'
     stop
  end if

  RGNOVD = 4.0 
  RGNSIZE = 7.746 
  ZLUM = 6.5 

  CALL GET_COMMAND_ARGUMENT(1, FESCARG)
  READ(FESCARG, '(F10.5)') FESC

  ! ----------------------------------------

  ! Set initial conditions. 
  INCLUDE 'readin.inc' 

  ! NGAMMA is calculated from Starburst99 model
  ! `reion-generic' by popsyn/ngtot.f90.
  NGAMMA = 2.82E70_PREC ! (10^10 M_solar)^-1 
  NGAMMA_POP2 = NGAMMA 
  ! NGAMMA_POP3 = 2.82E71_PREC ! (10^10 M_solar)^-1 
  ngamma_pop3 = ngammafrac(z)*1.0e10 ! (10^10 M_solar)^-1

  Z = INITIAL_REDSHIFT 
  Q = 1.0E-8_PREC 

  ! Allocate and initialize arrays that will store values of Jeans
  ! mass and filling factor for redshifts at which calculation is 
  ! done.  Length of these arrays will have to be equal to the 
  ! number of such redshift values. 
  NCALC = (FINAL_REDSHIFT-INITIAL_REDSHIFT)/DZ+1 
  ALLOCATE(JMHARR(NCALC)); ALLOCATE(FFARR(NCALC)); ALLOCATE(SFRARR(NCALC))
  allocate(metarr(ncalc)); allocate(febyharr(ncalc)); allocate(zfearr(ncalc))
  allocate(sfrarr_pop2(ncalc)); allocate(sfrarr_pop3(ncalc)) 
  allocate(nuintegralarr_pop2(ncalc)); allocate(nuintegralarr_pop3(ncalc)) 

  JMHARR = 0.0_PREC; FFARR = 0.0_PREC; SFRARR = 0.0_PREC; zfearr = 0.0_prec
  metarr = 0.0_prec; febyharr = 0.0_prec
  sfrarr_pop2 = 0.0_prec; sfrarr_pop3 = 0.0_prec 
  nuintegralarr_pop2 = 0.0_prec; nuintegralarr_pop3 = 0.0_prec 
  COUNTR = 1 
  FFARR(COUNTR) = Q 
  zfearr(countr) = z

  RHO = OMEGA_NR*RHO_CRITICAL ! 10^10 M_solar / Mpc^3 
  RHO_BARYON = OMEGA_B*RHO_CRITICAL ! 10^10 M_solar / Mpc^3 

  ! Comoving hydrogen _number_ density; 1.1891E57 is the number of
  ! protons in a M_solar.  Note the 10^10 due to our mass units.
  DOFZ = GRFARR(COUNTER(Z))
  NH = RHO_BARYON*(1.0_PREC+RGNOVD*DOFZ)*(1.1891E57_PREC*1.0E10_PREC) ! Mpc^-3 
  NH_PROPER = NH*(1.0_PREC+Z)**3 ! Mpc^-3 

  ! TEMPC is temperature of neutral gas at redshift 
  ! INITIAL_REDSHIFT.  Here I have estimated it by 
  ! assuming that CMB and the neutral gas were coupled
  ! till z = 200 after which the gas cooled ~ (1 + z)^2 
  ! while CMB cooled ~ (1 + z).  See Bharadwaj & Ali (2004).
  TEMPC = CMBTEMP*(1.0_PREC+Z)**2 / 201.0_PREC ! K 
  TEMPH = TEMPC ! K

  TEMPHVA = TEMPH ! K 

  ! JNSLN = JEANS_LENGTH(TEMPH,Z) ! Mpc 
  JNSLN = JEANS_LENGTH(TEMPHVA,Z) ! Mpc 

  JNSM = (4.0_PREC*PI*RHO_BARYON*JNSLN**3)/3.0_PREC ! 10^10 M_solar
  JMHARR(COUNTR) = JNSM ! 10^10 M_solar 

  ! Normalise matter power spectrum. 
  RCARRY_PSIG = 8.0_PREC/SMALLH ! Mpc 
  BOUND = 0.0_PREC
  INF = 1 
  CALL DQAGIE(PSIG,BOUND,INF,ABSERR,ABSREL,MAXINT,SUM,ERROR,&
       &NEVAL,IER,ALIST,BLIST,RLIST,ELIST,IORD,LAST)
  IF (IER > 0) WRITE (0,*) 'PSPECNORM: Error.' 
  PSPECNORM = (SIGMA_EIGHT**2/SUM)

  ! Calculate the rms linear mass fluctuation in baryons.
  SIGMAB = SIGMA_BARYON(Z, TEMPH) 

  IGMDCRIT = IGMDCRIT_PREO 
  R = CLUMPFAC(IGMDCRIT)
  r_local = r 
  F_M = IGMFPREO()
  FV = IGMVFRAC(IGMDCRIT) 
  ! LMFP = q**(1.0_prec/3.0_prec)*LMFP0*JNSLN/((1.0_PREC-q*FV)**(2.0_PREC/3.0_PREC)) ! Mpc 
  LMFP = q**(1.0_prec/3.0_prec)*LMFP0*JNSLN/((1.0_PREC-FV)**(2.0_PREC/3.0_PREC)) ! Mpc 

  X2INIT = 1.2E-5_PREC/(SMALLH*OMEGA_B) ! dimensionless 
  X_II = X2INIT
  OLDXII = X_II 

  X_IIVA = X_II
  OLDXIIVA = OLDXII 

  PREOVERLAP = .TRUE. 
  POSTOCOUNTER = 1 

  GPH_POP2 = GPH_KERNEL_POP2() 
  GPI_POP2 = GPI_KERNEL_POP2()
  GPH_POP3 = GPH_KERNEL_POP3() 
  GPI_POP3 = GPI_KERNEL_POP3()
  TAU = 0.0_PREC 
  TAU1 = TAUINT(FV, Q, X_IIVA, Z) 

  do 
     countr = countr + 1 
     z = z + dz 
     if (z < final_redshift) exit 

     ! Calculate rate of ionizing photons (nphdot). 
     source = sfr_hs(z) ! M_solar yr^-1 Mpc^-3
     nphdot = fesc*source*ngamma*1.0e-10 ! yr^-1 Mpc^-3

     !-------------------------

     ! Calculate case-A HII recombination coefficient (alpha_r).  This
     ! expression is from Hui and Gnedin (1997).
     reclambda = 2.0_prec*thi/temph 
     alpha_r = 1.269e-13_prec * (reclambda**1.503_prec) / &
          (1.0_prec + (reclambda/0.522_prec)**0.470_prec)**1.923_prec ! cm^3 s^-1 
     ! Add case-B recombination coeff. to that. 
     alpha_r = alpha_r + 2.753e-14_prec * (reclambda**1.5_prec) / &
          (1.0_prec + (reclambda/2.74_prec)**0.407_prec)**2.242_prec ! cm^3 s^-1 

     !-------------------------

     ! Calculate hydrogen number density (comoving and proper).  The
     ! conversion factor of 1.18e57 converts M_solar to number of
     ! hydrogen atoms. 1.0e10 reduces mass units to M_solar.
     dofz = grfarr(counter(z))
     nh = rho_baryon*(1.0_prec+rgnovd*dofz)*(1.1891e57_prec*1.0e10_prec) ! Mpc^-3 
     nh_proper = nh*(1.0_prec+z)**3 ! Mpc^-3      

     !-------------------------

     ! Update ionized hydrogen fraction (x_ii) in ionized region. 
     gamma_rec = r*nh_proper*alpha_r*cmbympccb*yrbys ! yr^-1 
     gamma_ion = (gpi_pop2*source*1.0e-10_prec)*&
          &fesc*lmfp*(1.0_prec+z)**3*(cmbympc**2)/q ! yr^-1 

     a = gamma_rec*dz*dtdz(z)
     b = 1.0_prec + dz*dtdz(z)*gamma_ion
     c = -(oldxii + dz*dtdz(z)*gamma_ion)

     qe = -0.5_prec*(b+sign(dsqrt(b**2-4.0_prec*a*c),b))
     r1 = qe/a; r2 = c/qe 

     oldxii = x_ii 
     if (r1<=1.0_prec .and. r1>=0.0_prec) then
        x_ii = r1
     else 
        x_ii = r2 
     end if

     !-------------------------

     ! Calculate case-A HII recombination cooling rate (recrate).
     ! This expression is from Hui and Gnedin (1997).
     reclambda = 2.0_prec*thi/temph 
     recrate = 1.778e-29_prec*temph*reclambda**1.965_prec/&
          &(1.0_prec+(reclambda/0.541_prec)**0.502_prec)&
          &**2.697_prec ! erg cm^3 s^-1
     ! Add case-B HII recombination cooling rate. 
     recrate = recrate + 3.435e-30_prec*temph*reclambda**1.970d0 / &
          (1.0_prec+(reclambda/2.250e0_prec)**0.376_prec)**3.720_prec ! erg cm^3 s^-1

     ! Update temperature of ionized regions (temph) and neutral
     ! regions (tempc).
     nhii = nh*oldxii ! Mpc^-3 
     gamma_rc = (nhii**2)*recrate*erg2j*cmbympccb*yrbys ! J Mpc^-3 yr^-1
     gamma_totc = r_local*gamma_rc*(1.0_prec+z)**6 ! J Mpc^-3 yr^-1

     ! Add Compton cooling 
     gamma_totc = gamma_totc + 5.65e-36_prec*((1.0_prec+z)**7)*&
          &(2.726*(1.0_prec+z)-temph)*nhii*erg2j*yrbys*cmbympccb ! J Mpc^-3 yr^-1 

     nhi = nh*(1.0_prec-oldxii) 
     gammaph = (gph_pop2*source*1.0e-10_prec)*fesc*&
          &lmfp*(1.0_prec+z)**3*(cmbympc**2)  ! J/yr 
     gamma_heat = gammaph*nhi*(1.0_prec+z)**3/q ! J mpc^-3 yr^-1 

     heatcool = 2.0_prec*(gamma_heat-gamma_totc)/&
          &(3.0_prec*nh_proper*(1.0_prec+x_ii)*kboltz) ! K yr^-1 
     temph = (temph + dz*dtdz(z)*heatcool) / &
          &(1.0_prec + 2.0_prec*hubp(z)*dz*dtdz(z)+ &
          &((x_ii-oldxii)/(1.0_prec+oldxii)))
     tempc = tempc+dz*dtempcdz(z) ! K

     mu_MeanMolWt = 1.0_prec / (1.0_prec + x_ii)

     term1 = dtdz(z)*heatcool
     term2 = 2.0_prec*hubp(z)*dz*dtdz(z)+ ((x_ii-oldxii)/(1.0_prec+x_ii))

     !------------------------- 

     ! Minimum circular velocity of haloes that can cool in ionized
     ! regions.
     vcirc = sqrt(2.0_prec*kboltz*temph/mproton)

     ! Calculate corresponding Jeans mass. 
     jnsln = jeans_length(temphva,z) ! Mpc 
     global_t = q*temph+(1.0_prec-q)*tempc
     sigmab = sigma_baryon(z, global_t) 
     jnsm = (4.0_prec*pi*rho_baryon*jnsln**3)/3.0_prec ! 10^10 M_solar
     jmharr(countr) = jnsm

     !------------------------- 

     ! Update filling factor for the pre-overlap phase.
     if (preoverlap) then 
        r = clumpfac(igmdcrit)
        r_local = r 

        oldfm = f_m 
        f_m = igmfpreo()
        fv = igmvfrac(igmdcrit) 
        lmfp = q**(1.0_prec/3.0_prec)*lmfp0*jnsln/((1.0_prec-fv)**(2.0_prec/3.0_prec)) ! Mpc 

        dfm = f_m - oldfm 
        t1 = dz*dtdz(z)*nphdot/(nh*f_m)
        t2 = dz*dtdz(z)*alpha_r*yrbys*cmbympccb*nh*x_ii*r*(1.0_prec+z)**3/f_m
        t2 = t2+dfm/f_m
        q = (q+t1)/(1.0_prec+t2)

        ffarr(countr) = q 
     end if

     !-------------------------

     ! Update filling factor (i.e. F_M) for the post-overlap phase.
     if (q>=1.0_prec) then 
        preoverlap = .false. 
        postocounter = postocounter + 1 

        solfm_z = z 
        solfm_igmdcrit = igmdcrit 
        solfm_ff = q 
        fmlo = f_m 
        fmhi = 0.99999999999_prec 
        tolzin = 1.0e-15_prec 
        oldigmfrac = f_m
        f_m = rtbis(solfm, fmlo, fmhi, tolzin, 'main-fm') 

        q = 1.0_prec 
        ffarr(countr) = q 
        soldlt_z = z
        soldlt_igmfrac = f_m 
        dlthi = 1.0e6_prec !igmdcrit+1.0e30_prec
        dltlo = 0.0_prec !igmdcrit-59.0_prec
        tolzin = 1.0e-8_prec
        dlt = rtbis(soldlt, dltlo, dlthi, tolzin, 'main-dlt')
        igmdcrit = dlt 

        r = clumpfac(igmdcrit)
        fv = igmvfrac(igmdcrit) 
        lmfp = q**(1.0_prec/3.0_prec)*lmfp0*jnsln/((1.0_prec-fv)**(2.0_prec/3.0_prec)) ! mpc 
     end if

     !-------------------------

     ! Calculate *mean* values of TEMPH and X_II.  (Above values are
     ! *mass averaged*, not mean.)
     include 'volav.inc'

     tau2 = tauint(fv, q, x_iiva, z) 
     tau = tau + 0.5_prec*dz*(tau1+tau2) 
     tau1 = tau2 

     dnlldz = speed_of_light*cmbympc*yrbys/&
          &(sqrt(pi)*lmfp*hubp(z)*(1.0_prec+z)) ! dimensionless 
     gammapi = (gpi_pop2*source*1.0e-10_prec)*fesc*lmfp*&
          &(1.0_prec+z)**3*(cmbympc**2)/yrbys ! s^-1
     
     write (*,'(F4.1,4E11.3E2)'), z, q, source, gammapi, nphdot 
     
     !-------------------------

     mminc = getjmc(z)
     mminh = getjmh(z)

  END DO

  TAU = TAU-FINAL_REDSHIFT*TAU1 
  WRITE(0,*) 'TAU=', TAU 

  !-------------------------

CONTAINS 

  FUNCTION DTEMPCDZ(REDSHIFT)

    ! Redshift derivative of temperature of neutral region.
    IMPLICIT NONE
    REAL(KIND = PREC), INTENT(IN) :: REDSHIFT
    REAL(KIND = PREC) :: DTEMPCDZ

    DTEMPCDZ = -2.0_PREC*HUBP(REDSHIFT)*TEMPC*DTDZ(REDSHIFT) ! K 

  END FUNCTION DTEMPCDZ

  FUNCTION TAUINT(FM, F2, X2, Z) 

    IMPLICIT NONE
    REAL(PREC), INTENT(IN) :: FM, F2, X2, Z
    REAL(PREC) :: TAUINT
    REAL(PREC) :: ECHARGE, MELECTRON, THOMSCROSS,&
         &ZDOT, HUB0_SMALLUNITS, NHCM, NE

    ECHARGE = 4.803204E-10_PREC ! esu
    MELECTRON = 9.109382E-28_PREC ! grams

    HUB0_SMALLUNITS = 3.24175E-18_PREC * SMALLH ! s^-1 
    THOMSCROSS = (8.0_PREC*PI/3.0_PREC)*(ECHARGE*ECHARGE&
         &/(MELECTRON*SPEED_OF_LIGHT**2))**2 ! cm^2  
    ZDOT = -HUB0_SMALLUNITS*(1.0_PREC+Z)*SQRT(OMEGA_NR*(1.0_PREC+Z)**3&
         & +OMEGA_LAMBDA) ! s^-1 
    NHCM = NH*CMBYMPCCB ! cm^-3 

    NE = FM*NHCM*X2*F2 + FM*(1.0_PREC-F2)*NHCM*X2INIT + &
         &(1.0_PREC-FM)*NHCM*X2INIT
    TAUINT = NE*SPEED_OF_LIGHT*THOMSCROSS&
         &*(1.0_PREC+Z)**3/ZDOT ! dimensionless 
  END FUNCTION TAUINT

  FUNCTION timedyn(rs)

    REAL(KIND = PREC), INTENT(IN) :: rs
    REAL(KIND = PREC) :: timedyn

    REAL(KIND = PREC) :: HDENS, BKGRHO, RHOCS

    RHOCS = 1.879E-29_PREC * SMALLH ** 2 ! g/cm^3
    BKGRHO = RHOCS * OMEGA_NR * (1.0_PREC + rs) ** 3 ! g/cm^3 
    HDENS = DELTA_VIRIAL * BKGRHO ! g/cm^3 
    timedyn = SQRT(3.0_PREC * PI / (16.0_PREC * NEWTG * &
         &HDENS * 1.0E3_PREC)) / YRBYS ! yr
    ! 1.0E3 converts g to kg and cm to m.  

  END FUNCTION timedyn

END PROGRAM REION

