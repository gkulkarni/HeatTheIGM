
program intspectrum 

  ! File: intspectrum.f90 
  ! Cre: 2010-04-20 
  ! Mod: $Date: 2010/04/21 16:43:53 $; $Revision: 1.2 $ 
  ! 
  ! This code reads in Starburst99 output *.spectrum, and for each
  ! wavelength available, integrates over time and divides by
  ! respective photon energy to give number of photons emitten per
  ! unit frequency per unit stellar mass.

  implicit none 
  double precision :: time1, time2, sum, luminosity, lambda, time,&
       &hplanck, erg2j, c, cmbyang, total_mass, yrbys, photon_energy, nu  
  double precision, dimension(:,:), allocatable :: spectrum 
  double precision, dimension(:), allocatable :: age, wavelength, ngnu 
  integer :: tcount, wcount, i, j 

  hplanck = 6.626069e-34 ! Js
  erg2j = 1.0e-7
  c = 2.998e10 ! cm/s 
  cmbyang = 1.0e8 
  total_mass = 1.0e6 
  yrbys = 3.154e7 

  open(unit=11, file='spectrum', status='old', action='read') 
  time1 = 0.0
  tcount = 0 
  wcount = 0 
  do 
     read(11,*,end=911) time2, lambda, luminosity 
     if (time2 /= time1) tcount = tcount+1 
     if (tcount == 1) wcount = wcount+1 
     time1 = time2 
  end do
911 continue 

  allocate(spectrum(wcount,tcount))
  allocate(age(tcount))
  allocate(wavelength(wcount)) 
  allocate(ngnu(wcount))

  rewind 11 
  do i = 1, tcount 
     do j = 1, wcount 
        read(11,*) time, lambda, luminosity 
        spectrum(j,i) = 10.0**luminosity ! erg/s/angstrom 
        age(i) = time ! yr 
        wavelength(j) = lambda ! angstrom 
     end do
  end do
  close (11) 

  do i = 1, wcount 
     do j = 1, tcount 
        photon_energy = hplanck*c*cmbyang/(erg2j*wavelength(i))
        spectrum(i,j) = spectrum(i,j)/photon_energy ! /s/angstrom 
     end do
  end do

  ! Do time integral.
  do i = 1, wcount 
     sum = age(1)*yrbys*spectrum(i,1) ! /angstrom 
     do j = 2, tcount 
        sum = sum + 0.5*(age(j)-age(j-1))*yrbys*(spectrum(i,j)+&
             &spectrum(i,j-1))
     end do
     sum = sum/total_mass ! /angstrom/m_solar 
     ngnu(i) = sum ! /angstrom/m_solar 
  end do

  do i = 1, wcount 
     nu = c*cmbyang/wavelength(i) ! Hz 
     print *, nu, ngnu(i)*wavelength(i)/nu  ! /Hz/m_solar 
  end do

end program intspectrum

