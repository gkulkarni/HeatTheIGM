
; File: ldla_talk.pro 
;  Cre: 2013-02-05
;  Mod: $Date: 2013/01/18 17:10:23 $ ($Revision: 1.1 $) 

set_plot, 'ps'
device, filename='ldla_talk.ps', xsize=7.0, ysize=7.0, $
        /inches, color=1, /HELVETICA, yoffset=1.0
!P.font = 0 

!P.charsize = 1.5
!P.thick = 3

TvLCT, 255, 0, 0, 2 
TvLCT, 0, 127, 255, 3
TvLCT, 255, 255, 0, 4

openr, lun, 'ldla4.dat', /get_lun
ldla_dat = fltarr(2,6)
readf, lun, ldla_dat
rs = ldla_dat[0,*]
ldla_avg = ldla_dat[1,*]
plotsym, 0, 1, /FILL
plot, rs, ldla_avg, xrange=[2,4.5], xtitle='redshift', $
      ytitle='dN/dX', yrange=[0.0,0.2], ythick=3, xthick=3 
close, lun
free_lun, lun 

openr, lun, 'ldla_obs.dat', /get_lun 
ldla_obsdata = fltarr(6,6) 
readf, lun, ldla_obsdata
close, lun 
free_lun, lun 

rs = ldla_obsdata[2,*]
ldla_obs = ldla_obsdata[3,*]
oplot, rs, ldla_obs, psym=8, color=2, thick=3

rserr = fltarr(2,6)
rserr[0,*] = rs[*] - ldla_obsdata[0,*] 
rserr[1,*] = ldla_obsdata[1,*] - rs[*]

ldla_err = fltarr(2,6) 
ldla_err[0,*] = ldla_obsdata[4,*]
ldla_err[1,*] = ldla_obsdata[5,*]

dy1 = ldla_err[1,*]
dy2 = ldla_err[0,*]

dx1 = rserr[0,*]
dx2 = rserr[1,*]

plot_err, rs, ldla_obs, dy1, dy2=dy2, dx1=dx1, dx2=dx2, color=2, thick=3

device, /close_file
set_plot, 'X'

