
program lumtot

  ! File: lumtot.f90 
  ! Cre: 2010-05-10
  ! Mod: $Date: 2010/05/14 04:25:46 $; $Revision: 1.1 $ 

  implicit none 
  double precision :: time1, time2, sum, luminosity, lambda, time,&
       &hplanck, erg2j, c, cmbyang, total_mass, yrbys, photon_energy,&
       &nu, total_age 
  double precision, dimension(:,:), allocatable :: spectrum 
  double precision, dimension(:), allocatable :: age, wavelength, ngnu 
  integer :: tcount, wcount, i, j 

  hplanck = 6.626069e-34 ! Js
  erg2j = 1.0e-7
  c = 2.998e10 ! cm/s 
  cmbyang = 1.0e8 
  total_mass = 1.0e6 ! m_solar 
  total_age = 9.8010e7 ! yr 
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

  ! Do time integral.
  do i = 1, wcount 
     sum = age(1)*yrbys*spectrum(i,1) ! erg/angstrom 
     do j = 2, tcount 
        sum = sum + 0.5*(age(j)-age(j-1))*yrbys*(spectrum(i,j)+&
             &spectrum(i,j-1))
     end do
     sum = sum/total_mass ! erg/angstrom/m_solar 
     sum = sum/(total_age*yrbys) ! erg/angstrom/m_solar/s
     ngnu(i) = sum ! erg/angstrom/m_solar/s
  end do

  do i = 1, wcount 
     nu = c*cmbyang/wavelength(i) ! Hz 
     print *, wavelength(i), nu, ngnu(i)*wavelength(i)/nu  ! erg/Hz/m_solar/s 
  end do

end program lumtot

