
function bckgnd_flux(z) 

  implicit none 
  use constants 
  real(kind=prec), intent(in) :: z 
  real(kind=prec) :: bckgnd_flux 

  

contains 

  function optical_depth(z_prime) 

    real(kind_prec), intent(in) :: z_prime 
    real(kind_prec) :: optical_depth 
    

  end function optical_depth

end function bckgnd_flux
