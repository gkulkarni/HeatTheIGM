
program ngtot 

  ! File: ngtot.f90 
  ! Cre: 2010-04-20
  ! Mod: $Date: 2010/04/20 15:16:42 $; $Revision: 1.1 $ 
  !
  ! This code calculates total number of ionizing photons per unit
  ! star formed, using Starburst99 output.
 
  implicit none 
  double precision :: age1, age2, ngamma, ngamma_total, ngtot_perunitmass, &
       &total_mass, yrbys 

  total_mass = 1.0e6 
  yrbys = 3.154e7 
  open(unit=11, file='quanta', status='old', action='read')
  read(11,*) age1, ngamma 
  age1 = age1*yrbys 
  ngamma = 10.0**ngamma 
  ngamma_total = age1*ngamma 
  do 
     read(11,*,end=911) age2, ngamma 
     age2 = age2*yrbys 
     ngamma = 10.0**ngamma 
     ngamma_total = ngamma_total + ngamma*(age2-age1)*0.5 
     age1 = age2 
  end do
911 continue 

  ngtot_perunitmass = ngamma_total/total_mass 
  print *, ngtot_perunitmass 

end program ngtot


