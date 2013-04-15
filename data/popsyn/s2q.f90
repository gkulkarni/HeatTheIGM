
program s2q 

  ! File: s2q.f90 
  ! Cre: 2010-04-20
  ! Mod: $Date: 2010/04/20 15:16:43 $; $Revision: 1.1 $ 
  ! 
  ! For Starburst99 output, check is *.quanta1 results if we integrate
  ! *.spectrum1 up to 912 Angstroms. 

  implicit none 
  double precision :: quanta, age, lambda1, luminosity1, lambda2,&
       &luminosity2, photon_energy1, quanta1, photon_energy2, &
       &quanta2, hplanck, erg2j, c, cmbyang  

  hplanck = 6.626069e-34 ! Js
  erg2j = 1.0e-7
  c = 2.998e10 ! cm/s 
  cmbyang = 1.0e8 
  
  open(unit=11, file='xaa', status='old', action='read') 
  quanta = 0.0 
  read(11,*) age, lambda1, luminosity1 
  luminosity1 = 10.0**luminosity1 ! erg/s/ang 
  photon_energy1 = hplanck*c*cmbyang/(erg2j*lambda1) 
  quanta1 = luminosity1/photon_energy1 
  do 
     read(11,*) age, lambda2, luminosity2
     if (lambda2 > 912.0) exit 
     luminosity2 = 10.0**luminosity2 ! erg/s/ang
     photon_energy2 = hplanck*c*cmbyang/(erg2j*lambda2) 
     quanta2 = luminosity2/photon_energy2 
     quanta = quanta + 0.5*(quanta1+quanta2)*(lambda2-lambda1) 
     lambda1 = lambda2 
     quanta1 = quanta2 
  end do

  print *, log10(quanta) 
  
end program s2q

